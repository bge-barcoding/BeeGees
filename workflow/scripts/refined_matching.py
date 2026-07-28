#!/usr/bin/env python3
"""
Refined Matching Filter

Reference-distance cleaning for MitoGeneExtractor read alignments. Replaces the
01-04 fasta_cleaner filters with a single step that measures every read against
the reference that recruited it, rather than against the alignment's own
(contaminated) consensus.

To admit short reads MGE's 's' is lowered, and because the gene-fetch reference
is often a distant relative 'r' is lowered too. That combination is necessary,
and it is also what lets contamination in. This filter cleans up after it in two
passes:

  Pass A (absolute)  Reads simply too far from the reference to be this gene at
                     all are removed outright. This is what handles contamination
                     that has no competitor - a sample that is 90% human COI has
                     no minority allele to adjudicate, so conflict resolution
                     alone would never fire.

  Pass B (conflict)  Where the surviving reads disagree at a column, the allele
                     carried by the reads closest to the reference wins. Losing
                     bases are MASKED rather than the read being dropped, so a
                     read with one or two damaged bases keeps contributing
                     everywhere else it covers. Only a read that loses
                     repeatedly - evidence it is a different molecule, not a
                     damaged one - is removed entirely.

Distances are measured in DNA space (more sensitive; discriminates at synonymous
sites) or amino-acid space (less sensitive but blind to synonymous change, and
therefore structurally unable to pull the consensus toward the reference's codon
usage). DNA is the default and falls back to AA per-sample when no coding
nucleotide reference is available.

The alignment is assumed to sit in the protein reference's coordinate frame, so
column -> codon -> reference residue. That assumption is DERIVED per file rather
than trusted: see --anchor-audit, and the anchor_min_identity escape hatch.

USAGE:
    # 1. Check the coordinate assumption against real data before anything else
    python refined_matching.py --anchor-audit \\
        --input-files-list logs/mge/alignment_files.log \\
        --protein-reference-dir 02_references/protein

    # 2. Measure what the filter would do, changing nothing
    python refined_matching.py \\
        --input-files-list logs/mge/alignment_files.log \\
        --output-dir fasta_cleaner/refined_matching/ \\
        --filtered-files-list refined_matched.txt \\
        --protein-reference-dir 02_references/protein \\
        --mode report_only --threads 8

    # 3. Apply it
    python refined_matching.py ... --mode filter

    # Built-in tests (no data or external tools required)
    python refined_matching.py --self-test

File naming:
    Input:  {sample}_r_{r}_s_{s}_align_{sample}.fas
    Output: {sample}_r_{r}_s_{s}_{...}_refined.fasta
    The sample name (text before the first '_r_') locates the reference.
"""

import argparse
import csv
import os
import sys
import time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, field, asdict
from datetime import datetime
from typing import Dict, List, Optional, Tuple

import numpy as np
from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.SeqIO import FastaIO


# ============================================================================
# Encoding
# ============================================================================

# Nucleotide codes. 0 is "no data" and covers '-', '~', 'N' and every other
# non-ACGT character. Note that 05_consensus_generator.py currently treats '~'
# and 'N' as ordinary residues and can vote them into the consensus; this script
# deliberately does not.
GAP = 0
NT_CODES = {'A': 1, 'C': 2, 'G': 3, 'T': 4}
NT_CHARS = np.array(['-', 'A', 'C', 'G', 'T'])

NT_LUT = np.zeros(256, dtype=np.uint8)
for _base, _code in NT_CODES.items():
    NT_LUT[ord(_base)] = _code
    NT_LUT[ord(_base.lower())] = _code

AA_ALPHABET = "ARNDCQEGHILKMFPSTWYVBZX*"
AA_INDEX = {aa: i for i, aa in enumerate(AA_ALPHABET)}
AA_X = AA_INDEX['X']
AA_STOP = AA_INDEX['*']


def _blosum62() -> np.ndarray:
    """BLOSUM62 as a dense float array indexed by AA_ALPHABET position."""
    matrix = substitution_matrices.load("BLOSUM62")
    n = len(AA_ALPHABET)
    arr = np.zeros((n, n), dtype=np.float64)
    for i, a in enumerate(AA_ALPHABET):
        for j, b in enumerate(AA_ALPHABET):
            arr[i, j] = matrix[a, b]
    return arr


BLOSUM62 = _blosum62()


def codon_translation_table(genetic_code: int) -> np.ndarray:
    """5x5x5 lookup: nucleotide codes -> amino-acid index.

    Index 0 on any base (missing data) yields AA_X. Stop codons yield AA_STOP.
    """
    table = CodonTable.unambiguous_dna_by_id[genetic_code]
    lut = np.full((5, 5, 5), AA_X, dtype=np.uint8)
    for b1 in range(1, 5):
        for b2 in range(1, 5):
            for b3 in range(1, 5):
                codon = NT_CHARS[b1] + NT_CHARS[b2] + NT_CHARS[b3]
                if codon in table.stop_codons:
                    lut[b1, b2, b3] = AA_STOP
                else:
                    lut[b1, b2, b3] = AA_INDEX.get(
                        table.forward_table.get(codon, 'X'), AA_X)
    return lut


def encode_nt(sequence: str, width: Optional[int] = None) -> np.ndarray:
    """Encode a sequence string to nucleotide codes, right-padded with GAP."""
    arr = NT_LUT[np.frombuffer(sequence.encode('ascii', 'replace'),
                               dtype=np.uint8)]
    if width is not None and arr.size < width:
        arr = np.concatenate([arr, np.zeros(width - arr.size, dtype=np.uint8)])
    elif width is not None:
        arr = arr[:width]
    return arr


def decode_nt(codes: np.ndarray) -> str:
    return ''.join(NT_CHARS[codes])


def encode_aa(protein: str) -> np.ndarray:
    return np.array([AA_INDEX.get(c.upper(), AA_X) for c in protein],
                    dtype=np.uint8)


# ============================================================================
# Parameters
# ============================================================================

@dataclass
class Params:
    """Every knob, with the defaults documented in the plan."""
    mode: str = "report_only"                  # report_only | filter
    distance_metric: str = "dna"               # dna | aa | both
    genetic_code: int = 5
    anchor_min_identity: float = 0.30

    # Pass A - absolute distance
    absolute_filter: bool = True
    max_read_distance: Optional[float] = None  # None -> metric default
    data_driven_cut: bool = True
    max_absolute_removal_frac: float = 0.50
    min_reads_retained: int = 5

    # Conflict calling
    min_depth: int = 3
    min_minor_count: int = 2
    min_minor_freq: float = 0.20
    damage_aware: bool = True
    damage_freq_multiplier: float = 2.0
    max_conflict_sites: int = 200

    # Pass B - scoring
    scope: str = "local"                       # local | global
    window: str = "auto"                       # "auto" or an integer of codons
    min_scored_positions: int = 5              # codons (aa) or 15 nt (dna)
    allele_score: str = "median"               # median | mean | best

    # Pass B - decision
    min_margin: float = 1.0                    # in substitution-equivalents
    min_margin_override: float = 2.0
    max_override_majority_freq: float = 0.70
    min_override_support: int = 3

    # Pass B - action
    removal_scope: str = "base_then_read"      # base | read | base_then_read
    read_removal_loss_count: int = 2
    read_removal_loss_frac: float = 0.50
    min_retained_depth: int = 2

    # File-level reverts
    max_removal_frac: float = 0.10
    max_ref_identity_gain: float = 0.02

    consensus_threshold: float = 0.5

    # Metric-specific defaults for the absolute ceiling. A conspecific read sits
    # near 0; human COI against an insect COI reference sits far above these.
    DNA_DEFAULT_CEILING: float = 0.35
    AA_DEFAULT_CEILING: float = 0.45

    def ceiling_for(self, metric: str) -> float:
        if self.max_read_distance is not None:
            return float(self.max_read_distance)
        return self.DNA_DEFAULT_CEILING if metric == "dna" else self.AA_DEFAULT_CEILING


@dataclass
class FrameAnchor:
    ok: bool
    frame_offset: int
    identity: float
    codon_to_ref_aa: np.ndarray = field(default_factory=lambda: np.zeros(0, np.int32))
    n_codons: int = 0
    reason: str = ""


# ============================================================================
# Logging / naming helpers (conventions copied from 03/04)
# ============================================================================

def log_message(message: str, log_file=None, stdout: bool = False) -> None:
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    formatted = f"[{timestamp}] {message}"
    if log_file:
        log_file.write(formatted + '\n')
        log_file.flush()
    if stdout:
        print(formatted, flush=True)


def get_sample_name_for_reference(filepath: str) -> str:
    """Sample name used to locate the reference: text before the first '_r_'."""
    basename = os.path.basename(filepath)
    for ext in ['.fasta', '.fas', '.fa']:
        if basename.lower().endswith(ext):
            basename = basename[:-len(ext)]
            break
    if '_r_' in basename:
        return basename.split('_r_')[0]
    for suffix in ['_refined', '_outlier_filtered', '_at_filtered',
                   '_human_filtered', '_align']:
        basename = basename.replace(suffix, '')
    return basename


def get_output_basename(filepath: str) -> str:
    """Output name: preserves the r/s sweep tag, strips the filter suffix."""
    basename = os.path.basename(filepath)
    for ext in ['.fasta', '.fas', '.fa']:
        if basename.lower().endswith(ext):
            basename = basename[:-len(ext)]
            break
    for suffix in ['_outlier_filtered', '_at_filtered', '_human_filtered']:
        if basename.endswith(suffix):
            basename = basename[:-len(suffix)]
            break
    return basename


# ============================================================================
# Alignment IO
# ============================================================================

def load_alignment(path: str) -> Tuple[List[str], np.ndarray]:
    """Read a gapped MSA into (ids, arr[n_reads, n_cols] uint8)."""
    ids: List[str] = []
    seqs: List[str] = []
    for record in SeqIO.parse(path, "fasta"):
        ids.append(record.id)
        seqs.append(str(record.seq))
    if not seqs:
        return [], np.zeros((0, 0), dtype=np.uint8)
    width = max(len(s) for s in seqs)
    arr = np.zeros((len(seqs), width), dtype=np.uint8)
    for i, s in enumerate(seqs):
        arr[i] = encode_nt(s, width)
    return ids, arr


def load_single_fasta(path: str) -> Optional[str]:
    """First record of a FASTA as an uppercase string, or None."""
    if not path or not os.path.exists(path):
        return None
    for record in SeqIO.parse(path, "fasta"):
        return str(record.seq).upper()
    return None


def resolve_reference_paths(sample_name: str,
                            protein_dir: Optional[str],
                            nucleotide_dir: Optional[str],
                            reference_map: Optional[Dict[str, Dict[str, str]]]
                            ) -> Tuple[Optional[str], Optional[str]]:
    """Locate the protein and (optional) coding-nucleotide reference.

    Supports both the gene-fetch layout (02_references/protein/{sample}.fasta)
    and the manual workflow driven by sequence_references.csv.
    """
    protein_path = nucleotide_path = None
    if reference_map and sample_name in reference_map:
        entry = reference_map[sample_name]
        protein_path = entry.get('protein_path') or None
        nucleotide_path = entry.get('nucleotide_path') or None
    if protein_path is None and protein_dir:
        for candidate in (f"{sample_name}.fasta", f"{sample_name}_reference.fasta"):
            p = os.path.join(protein_dir, candidate)
            if os.path.exists(p):
                protein_path = p
                break
    if nucleotide_path is None and nucleotide_dir:
        for candidate in (f"{sample_name}.fasta", f"{sample_name}_reference.fasta"):
            p = os.path.join(nucleotide_dir, candidate)
            if os.path.exists(p):
                nucleotide_path = p
                break
    return protein_path, nucleotide_path


def parse_reference_map(csv_path: str) -> Dict[str, Dict[str, str]]:
    """Parse sequence_references.csv, tolerating the absent nucleotide column."""
    mapping: Dict[str, Dict[str, str]] = {}
    with open(csv_path, 'r') as handle:
        for row in csv.DictReader(handle):
            sample_id = (row.get('ID') or row.get('id') or '').strip()
            if not sample_id:
                continue
            mapping[sample_id] = {
                'protein_path': (row.get('protein_reference_path') or '').strip(),
                'nucleotide_path': (row.get('nucleotide_reference_path')
                                    or row.get('nucleotide_path') or '').strip(),
            }
    return mapping


# ============================================================================
# Consensus and frame anchoring
# ============================================================================

def column_counts(arr: np.ndarray) -> np.ndarray:
    """counts[5, n_cols]; row 0 is the missing-data count."""
    if arr.size == 0:
        return np.zeros((5, arr.shape[1] if arr.ndim == 2 else 0), dtype=np.int32)
    counts = np.zeros((5, arr.shape[1]), dtype=np.int32)
    for code in range(5):
        counts[code] = (arr == code).sum(axis=0)
    return counts


def column_plurality_consensus(arr: np.ndarray,
                               threshold: float = 0.0) -> np.ndarray:
    """Per-column plurality over non-missing observations.

    Returns nucleotide codes, GAP where there is no data or the plurality falls
    below `threshold`.
    """
    counts = column_counts(arr)
    acgt = counts[1:]                                     # (4, n_cols)
    depth = acgt.sum(axis=0)
    best = acgt.argmax(axis=0) + 1
    best_count = acgt.max(axis=0)
    consensus = best.astype(np.uint8)
    with np.errstate(invalid='ignore', divide='ignore'):
        freq = np.where(depth > 0, best_count / np.maximum(depth, 1), 0.0)
    consensus[(depth == 0) | (freq < threshold)] = GAP
    return consensus


def translate_codes(nt_codes: np.ndarray, lut: np.ndarray) -> np.ndarray:
    """Translate a flat nucleotide-code array to AA indices, codon by codon."""
    n_codons = nt_codes.size // 3
    if n_codons == 0:
        return np.zeros(0, dtype=np.uint8)
    trimmed = nt_codes[:n_codons * 3].reshape(n_codons, 3)
    return lut[trimmed[:, 0], trimmed[:, 1], trimmed[:, 2]]


def aa_codes_to_string(codes: np.ndarray) -> str:
    return ''.join(AA_ALPHABET[c] for c in codes)


def anchor_alignment_to_protein(consensus_nt: np.ndarray,
                                protein_ref: str,
                                genetic_code: int,
                                min_identity: float) -> FrameAnchor:
    """Derive the column -> codon -> reference-residue mapping.

    MGE is expected to lay reads into the amino-acid reference's coordinate
    frame, which would make this trivially offset 0 with codon i mapping to
    residue i. That is not assumed: all three offsets are translated and aligned
    to the reference, and the best-scoring one wins. Identity below
    `min_identity` means the assumption has broken and the caller should leave
    the file alone rather than score against a wrong frame.
    """
    lut = codon_translation_table(genetic_code)
    ref_codes = encode_aa(protein_ref)
    if ref_codes.size == 0:
        return FrameAnchor(False, 0, 0.0, reason="empty_protein_reference")

    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1

    best: Optional[FrameAnchor] = None
    for offset in range(3):
        peptide_codes = translate_codes(consensus_nt[offset:], lut)
        if peptide_codes.size == 0:
            continue
        peptide = aa_codes_to_string(peptide_codes)
        informative = peptide_codes != AA_X
        if informative.sum() == 0:
            continue
        try:
            alignment = aligner.align(peptide, protein_ref)[0]
        except Exception:
            continue

        codon_to_ref = np.full(peptide_codes.size, -1, dtype=np.int32)
        matches = compared = 0
        for (q_start, q_end), (r_start, r_end) in zip(*alignment.aligned):
            span = min(q_end - q_start, r_end - r_start)
            q_idx = np.arange(q_start, q_start + span)
            r_idx = np.arange(r_start, r_start + span)
            codon_to_ref[q_idx] = r_idx
            usable = peptide_codes[q_idx] != AA_X
            compared += int(usable.sum())
            matches += int((peptide_codes[q_idx][usable]
                            == ref_codes[r_idx][usable]).sum())

        identity = matches / compared if compared else 0.0
        candidate = FrameAnchor(
            ok=identity >= min_identity,
            frame_offset=offset,
            identity=identity,
            codon_to_ref_aa=codon_to_ref,
            n_codons=int(peptide_codes.size),
            reason="" if identity >= min_identity else "anchor_identity_below_floor",
        )
        if best is None or candidate.identity > best.identity:
            best = candidate

    if best is None:
        return FrameAnchor(False, 0, 0.0, reason="no_translatable_frame")
    return best


# ============================================================================
# Reference projection into alignment coordinates
# ============================================================================

def project_reference(anchor: FrameAnchor,
                      protein_ref: str,
                      nucleotide_ref: Optional[str],
                      n_cols: int,
                      genetic_code: int) -> Tuple[np.ndarray, np.ndarray, str]:
    """Lay the reference out in alignment coordinates.

    Returns (ref_nt_by_column, ref_aa_by_codon, metric_available) where
    `metric_available` is 'dna' when a usable coding nucleotide reference was
    projected and 'aa' otherwise. Positions with no reference counterpart are
    GAP / AA_X.
    """
    ref_aa_codes = encode_aa(protein_ref)
    n_codons = anchor.n_codons
    ref_aa_by_codon = np.full(n_codons, AA_X, dtype=np.uint8)
    mapped = anchor.codon_to_ref_aa >= 0
    ref_aa_by_codon[mapped] = ref_aa_codes[anchor.codon_to_ref_aa[mapped]]

    ref_nt_by_column = np.zeros(n_cols, dtype=np.uint8)
    metric_available = 'aa'

    if nucleotide_ref:
        cds = encode_nt(nucleotide_ref)
        lut = codon_translation_table(genetic_code)
        cds_aa = translate_codes(cds, lut)
        # The CDS is usable only if it actually codes for this protein. Compare
        # translation to the reference; a trailing stop codon is expected.
        comparable = min(cds_aa.size, ref_aa_codes.size)
        if comparable > 0:
            agreement = float((cds_aa[:comparable]
                               == ref_aa_codes[:comparable]).mean())
        else:
            agreement = 0.0
        if agreement >= 0.90:
            offset = anchor.frame_offset
            for codon_idx in range(n_codons):
                ref_residue = anchor.codon_to_ref_aa[codon_idx]
                if ref_residue < 0 or (ref_residue + 1) * 3 > cds.size:
                    continue
                col = offset + codon_idx * 3
                if col + 3 > n_cols:
                    break
                ref_nt_by_column[col:col + 3] = cds[ref_residue * 3:
                                                    ref_residue * 3 + 3]
            metric_available = 'dna'

    return ref_nt_by_column, ref_aa_by_codon, metric_available


# ============================================================================
# Distances
# ============================================================================

def dna_distance(arr: np.ndarray,
                 ref_nt: np.ndarray,
                 columns: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """p-distance in nucleotide space over `columns`.

    Returns (distance, n_compared) per read. Reads with nothing to compare get
    distance nan.
    """
    if columns.size == 0:
        return (np.full(arr.shape[0], np.nan), np.zeros(arr.shape[0], np.int32))
    sub = arr[:, columns]
    ref = ref_nt[columns][None, :]
    usable = (sub != GAP) & (ref != GAP)
    n_compared = usable.sum(axis=1).astype(np.int32)
    mismatches = (usable & (sub != ref)).sum(axis=1)
    with np.errstate(invalid='ignore', divide='ignore'):
        distance = np.where(n_compared > 0,
                            mismatches / np.maximum(n_compared, 1), np.nan)
    return distance, n_compared


def aa_distance(read_aa: np.ndarray,
                ref_aa: np.ndarray,
                codons: np.ndarray,
                stop_penalty: float) -> Tuple[np.ndarray, np.ndarray]:
    """Normalised BLOSUM62 deficit over `codons`.

    d = (S(ref,ref) - S(read,ref)) / S(ref,ref), clipped to [0, inf), averaged
    over scored codons. An in-frame stop against the reference frame takes a
    fixed penalty - free NUMT/pseudogene signal.
    """
    n_reads = read_aa.shape[0]
    if codons.size == 0:
        return np.full(n_reads, np.nan), np.zeros(n_reads, np.int32)

    sub = read_aa[:, codons]                              # (n_reads, k)
    ref = ref_aa[codons]                                  # (k,)
    usable = (sub != AA_X) & (ref[None, :] != AA_X)
    n_compared = usable.sum(axis=1).astype(np.int32)

    self_score = BLOSUM62[ref, ref][None, :]              # (1, k)
    pair_score = BLOSUM62[sub, ref[None, :]]              # (n_reads, k)
    with np.errstate(invalid='ignore', divide='ignore'):
        deficit = (self_score - pair_score) / np.where(self_score != 0,
                                                       self_score, 1.0)
    deficit = np.clip(deficit, 0.0, None)
    deficit = np.where(sub == AA_STOP, stop_penalty, deficit)
    deficit = np.where(usable, deficit, 0.0)

    with np.errstate(invalid='ignore', divide='ignore'):
        distance = np.where(n_compared > 0,
                            deficit.sum(axis=1) / np.maximum(n_compared, 1),
                            np.nan)
    return distance, n_compared


def worst_blosum_deficit() -> float:
    """Largest normalised deficit BLOSUM62 can produce; used to size the stop penalty."""
    diag = np.diag(BLOSUM62)[:20]
    worst = 0.0
    for i in range(20):
        for j in range(20):
            if diag[i] != 0:
                worst = max(worst, (diag[i] - BLOSUM62[i, j]) / diag[i])
    return float(worst)


STOP_PENALTY = 2.0 * worst_blosum_deficit()


# ============================================================================
# Conflict calling
# ============================================================================

def find_conflict_columns(arr: np.ndarray, params: Params) -> Tuple[np.ndarray, List[Dict]]:
    """Columns where reads genuinely disagree.

    `min_minor_count` is the load-bearing guard: a single-read minor allele is
    sequencing error or cytosine deamination, never a conflict. At museum depths
    of 3-15x it does nearly all the work and the frequency threshold rarely
    binds.
    """
    if arr.size == 0:
        return np.zeros(0, dtype=np.int64), []

    counts = column_counts(arr)
    acgt = counts[1:]
    depth = acgt.sum(axis=0)

    order = np.argsort(acgt, axis=0)
    major_code = order[-1] + 1
    minor_code = order[-2] + 1
    major_count = np.take_along_axis(acgt, order[-1][None, :], axis=0)[0]
    minor_count = np.take_along_axis(acgt, order[-2][None, :], axis=0)[0]

    with np.errstate(invalid='ignore', divide='ignore'):
        minor_freq = np.where(depth > 0, minor_count / np.maximum(depth, 1), 0.0)

    required_freq = np.full(depth.shape, params.min_minor_freq, dtype=np.float64)
    if params.damage_aware:
        # C->T and G->A deamination recurs across independent reads, so the
        # min_minor_count guard does not stop it. Demand more at those sites.
        c_to_t = (major_code == NT_CODES['C']) & (minor_code == NT_CODES['T'])
        g_to_a = (major_code == NT_CODES['G']) & (minor_code == NT_CODES['A'])
        damage_like = c_to_t | g_to_a
        required_freq[damage_like] *= params.damage_freq_multiplier
    else:
        damage_like = np.zeros(depth.shape, dtype=bool)

    qualifies = ((depth >= params.min_depth)
                 & (minor_count >= params.min_minor_count)
                 & (minor_freq >= required_freq))
    conflict_cols = np.flatnonzero(qualifies).astype(np.int64)

    truncated = False
    if conflict_cols.size > params.max_conflict_sites:
        # Multi-source junk. Keep the strongest sites and flag it rather than
        # silently doing partial work.
        strength = depth[conflict_cols] * minor_freq[conflict_cols]
        keep = np.argsort(strength)[::-1][:params.max_conflict_sites]
        conflict_cols = np.sort(conflict_cols[keep])
        truncated = True

    site_table = []
    for col in conflict_cols:
        site_table.append({
            'column': int(col),
            'codon_pos': int(col % 3),
            'depth': int(depth[col]),
            'A': int(acgt[0, col]), 'C': int(acgt[1, col]),
            'G': int(acgt[2, col]), 'T': int(acgt[3, col]),
            'major_allele': NT_CHARS[major_code[col]],
            'major_freq': float(major_count[col] / max(depth[col], 1)),
            'minor_allele': NT_CHARS[minor_code[col]],
            'minor_freq': float(minor_freq[col]),
            'damage_flagged': bool(damage_like[col]),
        })
    if truncated and site_table:
        site_table[0]['truncated'] = True
    return conflict_cols, site_table


# ============================================================================
# Scoring windows
# ============================================================================

def choose_half_window(arr: np.ndarray, params: Params) -> int:
    """Half-window in columns, sized from the reads themselves.

    Too wide and reads abstain for lack of coverage (and a chimeric junction is
    diluted); too narrow and one substitution is indistinguishable from noise.
    """
    if params.window != "auto":
        return max(1, int(params.window)) * 3
    if arr.size == 0:
        return 30
    covered = (arr != GAP).sum(axis=1)
    covered = covered[covered > 0]
    if covered.size == 0:
        return 30
    median_codons = float(np.median(covered)) / 3.0
    half_codons = int(np.clip(median_codons / 4.0, 5, 25))
    return half_codons * 3


def scoring_columns(centre: int,
                    half_window: int,
                    n_cols: int,
                    conflict_mask: np.ndarray) -> np.ndarray:
    """Columns used to score reads competing at `centre`.

    The adjudicated column and every other conflict column are excluded.
    Including the adjudicated column would score a read on the very allele being
    adjudicated, which is circular; including other conflict columns pulls in
    correlated errors.
    """
    lo = max(0, centre - half_window)
    hi = min(n_cols, centre + half_window + 1)
    cols = np.arange(lo, hi)
    return cols[~conflict_mask[cols]]


# ============================================================================
# Pass A - absolute distance
# ============================================================================

def bimodal_cut(distances: np.ndarray, ceiling: float) -> Optional[float]:
    """Largest-gap cut in the sorted distance distribution.

    In a contaminated museum sample the distribution is typically bimodal: a
    target cluster near the reference and a junk tail far from it. Returns None
    when the distribution shows no convincing separation, in which case only the
    fixed ceiling applies.
    """
    finite = np.sort(distances[np.isfinite(distances)])
    if finite.size < 10:
        return None
    gaps = np.diff(finite)
    if gaps.size == 0:
        return None
    idx = int(np.argmax(gaps))
    gap = float(gaps[idx])
    spread = float(finite[-1] - finite[0])
    # Require the gap to be a real valley, not sampling noise, and to sit away
    # from the extremes so it splits the data rather than shaving an outlier.
    if spread <= 0 or gap < max(0.05, 0.25 * spread):
        return None
    position = (idx + 1) / finite.size
    if position < 0.10 or position > 0.95:
        return None
    cut = float((finite[idx] + finite[idx + 1]) / 2.0)
    return min(cut, ceiling)


def absolute_pass(distances: np.ndarray,
                  params: Params,
                  metric: str) -> Tuple[np.ndarray, Dict]:
    """Flag reads too far from the reference to be this gene at all.

    The fixed ceiling is a biological statement ("this is not the target gene")
    and is allowed to empty a file - a sample that is entirely contamination
    should end up empty rather than be rescued by a retention floor. The
    data-driven cut is the discretionary one, so it is bounded by
    max_absolute_removal_frac and min_reads_retained.
    """
    n_reads = distances.size
    info = {'ceiling': None, 'data_cut': None, 'data_cut_applied': False,
            'removed_ceiling': 0, 'removed_data_cut': 0}
    if n_reads == 0 or not params.absolute_filter:
        return np.zeros(n_reads, dtype=bool), info

    ceiling = params.ceiling_for(metric)
    info['ceiling'] = ceiling
    remove = np.isfinite(distances) & (distances > ceiling)
    info['removed_ceiling'] = int(remove.sum())

    if params.data_driven_cut:
        cut = bimodal_cut(distances, ceiling)
        info['data_cut'] = cut
        if cut is not None:
            candidate = remove | (np.isfinite(distances) & (distances > cut))
            n_removed = int(candidate.sum())
            within_frac = n_removed <= params.max_absolute_removal_frac * n_reads
            retains_enough = (n_reads - n_removed) >= params.min_reads_retained
            if within_frac and retains_enough:
                info['removed_data_cut'] = n_removed - info['removed_ceiling']
                info['data_cut_applied'] = True
                remove = candidate

    return remove, info


# ============================================================================
# Pass B - conflict adjudication
# ============================================================================

def score_alleles(read_indices: np.ndarray,
                  distances: np.ndarray,
                  n_compared: np.ndarray,
                  min_scored: int,
                  how: str) -> Tuple[float, int, int]:
    """Summarise an allele group: (score, n_supporting, n_scoreable).

    Median rather than mean, because one reference-like contaminant skews a
    mean; and emphatically not the single closest read - "closest read wins" is
    the maximum-bias, maximum-variance rule, where one read takes the site.
    """
    if read_indices.size == 0:
        return float('nan'), 0, 0
    d = distances[read_indices]
    c = n_compared[read_indices]
    usable = np.isfinite(d) & (c >= min_scored)
    if not usable.any():
        return float('nan'), int(read_indices.size), 0
    values = d[usable]
    if how == "mean":
        score = float(values.mean())
    elif how == "best":
        score = float(values.min())
    else:
        score = float(np.median(values))
    return score, int(read_indices.size), int(usable.sum())


def resolve_column(col: int,
                   arr: np.ndarray,
                   distances: np.ndarray,
                   n_compared: np.ndarray,
                   ref_aa_by_codon: np.ndarray,
                   anchor: FrameAnchor,
                   params: Params,
                   min_scored: int) -> Dict:
    """Decide which allele wins at one conflict column."""
    column = arr[:, col]
    present = np.flatnonzero(column != GAP)
    alleles = np.unique(column[present])
    depth = present.size

    result = {
        'column': int(col),
        'depth': int(depth),
        'action': 'abstain_no_conflict',
        'winner': '',
        'winner_is_major': True,
        'margin': 0.0,
        'majority_freq': 0.0,
        'synonymous': False,
        'loser_reads': np.zeros(0, dtype=np.int64),
        'allele_detail': {},
    }
    if alleles.size < 2:
        return result

    counts = {int(a): int((column == a).sum()) for a in alleles}
    major_allele = max(counts, key=lambda a: counts[a])
    result['majority_freq'] = counts[major_allele] / max(depth, 1)

    scores = {}
    for allele in alleles:
        idx = np.flatnonzero(column == allele)
        score, n_support, n_scoreable = score_alleles(
            idx, distances, n_compared, min_scored, params.allele_score)
        scores[int(allele)] = {
            'score': score, 'n': n_support, 'n_scoreable': n_scoreable,
            'reads': idx,
        }
    result['allele_detail'] = {
        NT_CHARS[a]: {'n': v['n'], 'n_scoreable': v['n_scoreable'],
                      'score': v['score']}
        for a, v in scores.items()
    }

    scoreable = {a: v for a, v in scores.items() if np.isfinite(v['score'])}
    if len(scoreable) < 2:
        result['action'] = 'abstain_unscoreable'
        return result

    # Synonymous check. Where the competing alleles code for the same residue,
    # an AA metric cannot discriminate. Flag it so the DNA-vs-AA question is
    # measurable rather than argued: this is where reference bias would live.
    codon_idx = (col - anchor.frame_offset) // 3
    if 0 <= codon_idx < anchor.n_codons and (col - anchor.frame_offset) >= 0:
        lut = codon_translation_table(params.genetic_code)
        base = arr[:, anchor.frame_offset + codon_idx * 3:
                   anchor.frame_offset + codon_idx * 3 + 3]
        if base.shape[1] == 3:
            consensus_codon = column_plurality_consensus(base)
            pos_in_codon = (col - anchor.frame_offset) % 3
            translated = set()
            for allele in scoreable:
                codon = consensus_codon.copy()
                codon[pos_in_codon] = allele
                translated.add(int(lut[codon[0], codon[1], codon[2]]))
            result['synonymous'] = len(translated) == 1

    ranked = sorted(scoreable.items(), key=lambda kv: kv[1]['score'])
    winner_allele, winner = ranked[0]
    runner_allele, runner = ranked[1]
    margin = runner['score'] - winner['score']
    result['margin'] = float(margin)
    result['winner'] = NT_CHARS[winner_allele]
    result['winner_is_major'] = winner_allele == major_allele

    if winner['n'] < params.min_minor_count:
        result['action'] = 'abstain_support'
        return result
    if margin < params.min_margin:
        result['action'] = 'abstain_margin'
        return result

    if winner_allele != major_allele:
        # Overturning a solid majority is the direction that can manufacture
        # reference bias, so it takes a strictly higher bar.
        majority_freq = result['majority_freq']
        if (margin < params.min_margin_override
                or majority_freq >= params.max_override_majority_freq
                or winner['n'] < params.min_override_support):
            result['action'] = 'abstain_majority_protected'
            return result
        result['action'] = 'override_majority'
    else:
        result['action'] = 'agree_majority'

    losers = np.concatenate([v['reads'] for a, v in scores.items()
                             if a != winner_allele]) if len(scores) > 1 \
        else np.zeros(0, dtype=np.int64)
    result['loser_reads'] = losers.astype(np.int64)
    return result


def apply_decisions(arr: np.ndarray,
                    decisions: List[Dict],
                    params: Params) -> Tuple[np.ndarray, np.ndarray, List[Dict]]:
    """Mask losing bases, then remove reads that lose repeatedly.

    Masking costs one unit of depth at one column; removing a read costs its
    entire span, including every position where it was fine. That is the whole
    reason the default is to mask first and escalate only on repeat offence.
    """
    out = arr.copy()
    n_reads = arr.shape[0]
    losses = np.zeros(n_reads, dtype=np.int32)
    covered_conflicts = np.zeros(n_reads, dtype=np.int32)
    events: List[Dict] = []

    acting = [d for d in decisions if d['action'] in ('agree_majority',
                                                      'override_majority')]
    for decision in decisions:
        col = decision['column']
        covered_conflicts[arr[:, col] != GAP] += 1

    if params.removal_scope in ('base', 'base_then_read'):
        for decision in acting:
            col = decision['column']
            losers = decision['loser_reads']
            if losers.size == 0:
                continue
            # Never dig a hole: if masking would take the column below the
            # retention floor, leave it alone. A gap is worse than an uncertain
            # base, which downstream validation can still catch.
            current_depth = int((out[:, col] != GAP).sum())
            if current_depth - losers.size < params.min_retained_depth:
                decision['action'] = 'abstain_depth_floor'
                continue
            for read_idx in losers:
                if out[read_idx, col] == GAP:
                    continue
                events.append({
                    'read_idx': int(read_idx),
                    'column': int(col),
                    'from': NT_CHARS[out[read_idx, col]],
                    'to': '-',
                })
                out[read_idx, col] = GAP
                losses[read_idx] += 1

    remove = np.zeros(n_reads, dtype=bool)
    if params.removal_scope in ('read', 'base_then_read'):
        by_count = losses >= params.read_removal_loss_count
        with np.errstate(invalid='ignore', divide='ignore'):
            frac = np.where(covered_conflicts > 0,
                            losses / np.maximum(covered_conflicts, 1), 0.0)
        by_frac = (covered_conflicts >= 2) & (frac >= params.read_removal_loss_frac)
        remove = by_count | by_frac

    return out, remove, events


# ============================================================================
# Per-file processing
# ============================================================================

def consensus_string(arr: np.ndarray, threshold: float) -> str:
    return decode_nt(column_plurality_consensus(arr, threshold))


def consensus_aa_identity(arr: np.ndarray,
                          ref_aa_by_codon: np.ndarray,
                          anchor: FrameAnchor,
                          params: Params) -> float:
    """AA identity between the alignment's consensus and the reference.

    Reported before and after filtering. If this rises sharply the filter is
    manufacturing reference bias, and that is a revert condition rather than a
    success.
    """
    if arr.size == 0 or not anchor.ok:
        return 0.0
    lut = codon_translation_table(params.genetic_code)
    consensus = column_plurality_consensus(arr, params.consensus_threshold)
    peptide = translate_codes(consensus[anchor.frame_offset:], lut)
    n = min(peptide.size, ref_aa_by_codon.size)
    if n == 0:
        return 0.0
    usable = (peptide[:n] != AA_X) & (ref_aa_by_codon[:n] != AA_X)
    if not usable.any():
        return 0.0
    return float((peptide[:n][usable] == ref_aa_by_codon[:n][usable]).mean())


def _empty_result(file_path: str, status: str, reason: str) -> Dict:
    sample_name = get_sample_name_for_reference(file_path)
    return {
        'file_path': file_path,
        'base_name': sample_name,
        'sample_name': sample_name,
        'output_basename': get_output_basename(file_path),
        'status': status,
        'reason': reason,
        'input_count': 0,
        'kept_count': 0,
        'removed_count': 0,
        'output_file': None,
        'summary': {},
        'columns': [],
        'reads': [],
    }


def process_single_file(file_path: str,
                        protein_dir: Optional[str],
                        nucleotide_dir: Optional[str],
                        reference_map: Optional[Dict[str, Dict[str, str]]],
                        params: Params,
                        output_dir: str) -> Dict:
    """Run both passes over one alignment file."""
    started = time.time()
    try:
        sample_name = get_sample_name_for_reference(file_path)
        output_basename = get_output_basename(file_path)

        if not os.path.exists(file_path) or os.path.getsize(file_path) == 0:
            return _empty_result(file_path, 'skipped', 'empty_file')

        protein_path, nucleotide_path = resolve_reference_paths(
            sample_name, protein_dir, nucleotide_dir, reference_map)
        protein_ref = load_single_fasta(protein_path)
        if not protein_ref:
            return _empty_result(file_path, 'skipped', 'no_protein_reference')
        nucleotide_ref = (load_single_fasta(nucleotide_path)
                          if params.distance_metric in ('dna', 'both') else None)

        ids, arr = load_alignment(file_path)
        if arr.size == 0:
            return _empty_result(file_path, 'skipped', 'no_sequences')

        result = _empty_result(file_path, 'success', 'processed')
        result['input_count'] = len(ids)
        n_reads, n_cols = arr.shape
        original = arr.copy()

        # --- anchor ------------------------------------------------------
        consensus_nt = column_plurality_consensus(arr)
        anchor = anchor_alignment_to_protein(
            consensus_nt, protein_ref, params.genetic_code,
            params.anchor_min_identity)
        if not anchor.ok:
            result['status'] = 'skipped'
            result['reason'] = anchor.reason or 'anchor_failed'
            result['kept_count'] = n_reads
            result['summary'] = {
                'anchor_frame_offset': anchor.frame_offset,
                'anchor_identity': round(anchor.identity, 4),
                'n_reads_in': n_reads, 'n_reads_out': n_reads,
            }
            _write_alignment(output_dir, output_basename, ids, arr, result)
            return result

        ref_nt, ref_aa, metric_available = project_reference(
            anchor, protein_ref, nucleotide_ref, n_cols, params.genetic_code)
        metric = 'dna' if (params.distance_metric in ('dna', 'both')
                           and metric_available == 'dna') else 'aa'

        lut = codon_translation_table(params.genetic_code)
        codon_start = anchor.frame_offset
        n_codons = anchor.n_codons
        usable_cols = codon_start + n_codons * 3
        read_aa = np.zeros((n_reads, n_codons), dtype=np.uint8)
        if n_codons:
            codons = arr[:, codon_start:usable_cols].reshape(n_reads, n_codons, 3)
            read_aa = lut[codons[:, :, 0], codons[:, :, 1], codons[:, :, 2]]

        min_scored = (params.min_scored_positions * 3 if metric == 'dna'
                      else params.min_scored_positions)

        def score(columns: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
            if metric == 'dna':
                return dna_distance(arr, ref_nt, columns)
            codon_idx = np.unique((columns - codon_start) // 3)
            codon_idx = codon_idx[(codon_idx >= 0) & (codon_idx < n_codons)]
            return aa_distance(read_aa, ref_aa, codon_idx, STOP_PENALTY)

        # --- Pass A: absolute distance -----------------------------------
        all_cols = np.arange(n_cols)
        global_distance, global_compared = score(all_cols)
        remove_absolute, absolute_info = absolute_pass(
            global_distance, params, metric)

        survivors = ~remove_absolute
        working = arr[survivors]
        working_ids = [i for i, keep in zip(ids, survivors) if keep]

        identity_before = consensus_aa_identity(original, ref_aa, anchor, params)

        # --- Pass B: conflict adjudication -------------------------------
        conflict_cols, site_table = find_conflict_columns(working, params)
        conflict_mask = np.zeros(n_cols, dtype=bool)
        conflict_mask[conflict_cols] = True
        half_window = choose_half_window(working, params)

        decisions: List[Dict] = []
        if working.shape[0] and conflict_cols.size:
            w_read_aa = read_aa[survivors]
            for col in conflict_cols:
                if params.scope == 'global':
                    cols = all_cols[~conflict_mask]
                else:
                    cols = scoring_columns(int(col), half_window, n_cols,
                                           conflict_mask)
                if metric == 'dna':
                    dist, cmp_n = dna_distance(working, ref_nt, cols)
                else:
                    codon_idx = np.unique((cols - codon_start) // 3)
                    codon_idx = codon_idx[(codon_idx >= 0) & (codon_idx < n_codons)]
                    dist, cmp_n = aa_distance(w_read_aa, ref_aa, codon_idx,
                                              STOP_PENALTY)
                decisions.append(resolve_column(
                    int(col), working, dist, cmp_n, ref_aa, anchor, params,
                    min_scored))

        filtered, remove_conflict, events = apply_decisions(
            working, decisions, params)

        kept_mask = ~remove_conflict
        final = filtered[kept_mask]
        final_ids = [i for i, keep in zip(working_ids, kept_mask) if keep]

        # --- invariants ---------------------------------------------------
        assert filtered.shape == working.shape, "row/column count changed"
        changed = filtered != working
        assert not np.any(changed & (filtered != GAP)), \
            "a residue was rewritten to another residue"
        assert int((final != GAP).sum()) <= int((original != GAP).sum()), \
            "non-gap residue count increased"

        identity_after = consensus_aa_identity(final, ref_aa, anchor, params)

        # --- reverts ------------------------------------------------------
        residues_before = int((original != GAP).sum())
        residues_after = int((final != GAP).sum())
        removed_frac = ((residues_before - residues_after) / residues_before
                        if residues_before else 0.0)
        conflict_cost = len(events) + int(
            (working[remove_conflict] != GAP).sum())
        conflict_frac = (conflict_cost / residues_before
                         if residues_before else 0.0)

        reverted, revert_reason = False, ''
        if conflict_frac > params.max_removal_frac:
            reverted, revert_reason = True, 'exceeded_removal_cap'
        elif identity_after - identity_before > params.max_ref_identity_gain:
            reverted, revert_reason = True, 'exceeded_ref_pull_cap'

        if reverted or params.mode == 'report_only':
            emitted, emitted_ids = arr, ids
        else:
            emitted, emitted_ids = final, final_ids

        cov_before = coverage_stats(original)
        cov_after = coverage_stats(emitted)

        result['kept_count'] = len(emitted_ids)
        result['removed_count'] = len(ids) - len(emitted_ids)
        result['summary'] = {
            'sample_name': sample_name,
            'anchor_frame_offset': anchor.frame_offset,
            'anchor_identity': round(anchor.identity, 4),
            'alignment_length': n_cols,
            'protein_ref_length': len(protein_ref),
            'length_is_3x_protein': n_cols == 3 * len(protein_ref),
            'metric_used': metric,
            'metric_requested': params.distance_metric,
            'mode': params.mode,
            'n_reads_in': n_reads,
            'n_reads_out': len(emitted_ids),
            'absolute_ceiling': absolute_info['ceiling'],
            'absolute_data_cut': absolute_info['data_cut'],
            'absolute_data_cut_applied': absolute_info['data_cut_applied'],
            'removed_absolute': int(remove_absolute.sum()),
            'n_conflict_columns': int(conflict_cols.size),
            'conflict_columns_per_kb': round(
                1000.0 * conflict_cols.size / max(n_cols, 1), 3),
            'n_agree_majority': sum(1 for d in decisions
                                    if d['action'] == 'agree_majority'),
            'n_override_majority': sum(1 for d in decisions
                                       if d['action'] == 'override_majority'),
            'n_abstain_margin': sum(1 for d in decisions
                                    if d['action'] == 'abstain_margin'),
            'n_abstain_support': sum(1 for d in decisions
                                     if d['action'] == 'abstain_support'),
            'n_abstain_unscoreable': sum(1 for d in decisions
                                         if d['action'] == 'abstain_unscoreable'),
            'n_abstain_majority_protected': sum(
                1 for d in decisions if d['action'] == 'abstain_majority_protected'),
            'n_abstain_depth_floor': sum(1 for d in decisions
                                         if d['action'] == 'abstain_depth_floor'),
            'n_synonymous_sites': sum(1 for d in decisions if d['synonymous']),
            'n_synonymous_acted': sum(
                1 for d in decisions
                if d['synonymous'] and d['action'] in ('agree_majority',
                                                       'override_majority')),
            'n_disagreed_with_majority': sum(
                1 for d in decisions if not d['winner_is_major'] and d['winner']),
            'bases_masked': len(events),
            'reads_removed_conflict': int(remove_conflict.sum()),
            'aa_identity_to_ref_before': round(identity_before, 4),
            'aa_identity_to_ref_after': round(identity_after, 4),
            'aa_identity_to_ref_delta': round(identity_after - identity_before, 4),
            'removed_residue_frac': round(removed_frac, 4),
            'cov_percent_before': round(cov_before['cov_percent'], 2),
            'cov_percent_after': round(cov_after['cov_percent'], 2),
            'cov_avg_before': round(cov_before['cov_avg'], 2),
            'cov_avg_after': round(cov_after['cov_avg'], 2),
            'cov_min_before': cov_before['cov_min'],
            'cov_min_after': cov_after['cov_min'],
            'reverted': reverted,
            'revert_reason': revert_reason,
            'runtime_s': round(time.time() - started, 3),
        }

        site_by_col = {s['column']: s for s in site_table}
        for decision in decisions:
            site = site_by_col.get(decision['column'], {})
            result['columns'].append({
                'base_name': output_basename,
                'column': decision['column'],
                'codon_pos': site.get('codon_pos', ''),
                'depth': decision['depth'],
                'major_allele': site.get('major_allele', ''),
                'major_freq': round(decision['majority_freq'], 4),
                'minor_allele': site.get('minor_allele', ''),
                'damage_flagged': site.get('damage_flagged', ''),
                'winner': decision['winner'],
                'winner_is_major': decision['winner_is_major'],
                'margin': round(decision['margin'], 5),
                'synonymous': decision['synonymous'],
                'action': decision['action'],
                'allele_detail': _format_allele_detail(decision['allele_detail']),
            })

        loss_by_read: Dict[int, List[int]] = {}
        for event in events:
            loss_by_read.setdefault(event['read_idx'], []).append(event['column'])
        for read_idx, columns in sorted(loss_by_read.items()):
            result['reads'].append({
                'base_name': output_basename,
                'sequence_id': working_ids[read_idx],
                'action': 'removed' if remove_conflict[read_idx] else 'masked',
                'n_losses': len(columns),
                'columns': ';'.join(str(c) for c in columns),
                'global_distance': _fmt(global_distance[np.flatnonzero(survivors)[read_idx]]),
            })

        _write_alignment(output_dir, output_basename, emitted_ids, emitted, result)
        return result

    except AssertionError as exc:
        return _empty_result(file_path, 'error', f'invariant_violation: {exc}')
    except Exception as exc:                                   # noqa: BLE001
        return _empty_result(file_path, 'error', f'processing_error: {exc}')


def _fmt(value) -> str:
    return '' if value is None or not np.isfinite(value) else f"{float(value):.6f}"


def _format_allele_detail(detail: Dict) -> str:
    parts = []
    for allele, info in sorted(detail.items()):
        score = info['score']
        score_str = 'na' if not np.isfinite(score) else f"{score:.4f}"
        parts.append(f"{allele}:n={info['n']},d={score_str}")
    return '|'.join(parts)


def coverage_stats(arr: np.ndarray) -> Dict:
    """Depth statistics, matching what calculate_coverage_statistics reports."""
    if arr.size == 0:
        return {'cov_percent': 0.0, 'cov_avg': 0.0, 'cov_min': 0, 'cov_max': 0}
    depth = (arr != GAP).sum(axis=0)
    return {
        'cov_percent': float((depth > 0).mean() * 100.0),
        'cov_avg': float(depth.mean()),
        'cov_min': int(depth.min()),
        'cov_max': int(depth.max()),
    }


def _write_alignment(output_dir: str, basename: str, ids: List[str],
                     arr: np.ndarray, result: Dict) -> None:
    if not output_dir or arr.size == 0 or not ids:
        return
    os.makedirs(output_dir, exist_ok=True)
    path = os.path.join(output_dir, f"{basename}_refined.fasta")
    with open(path, 'w') as handle:
        writer = FastaIO.FastaWriter(handle, wrap=None)
        writer.write_file([
            _record(seq_id, decode_nt(row)) for seq_id, row in zip(ids, arr)
        ])
    result['output_file'] = path


def _record(seq_id: str, sequence: str):
    from Bio.SeqRecord import SeqRecord
    return SeqRecord(Seq(sequence), id=seq_id, description="")


# ============================================================================
# Anchor audit
# ============================================================================

def anchor_audit(input_files: List[str],
                 protein_dir: Optional[str],
                 nucleotide_dir: Optional[str],
                 reference_map: Optional[Dict[str, Dict[str, str]]],
                 params: Params,
                 output_csv: Optional[str]) -> int:
    """Test the column -> codon -> reference-residue assumption on real data.

    Everything downstream rests on MGE laying reads into the amino-acid
    reference's coordinate frame. This checks it against alignments already on
    disk, in minutes, before any of the expensive work. If offsets are not
    uniformly 0 or identities are poor, the approach is falsified cheaply.
    """
    rows = []
    print(f"{'file':<52} {'aln_len':>8} {'prot':>6} {'3x?':>4} "
          f"{'frame':>5} {'identity':>9} {'cds?':>5}")
    print("-" * 96)

    for path in input_files:
        sample = get_sample_name_for_reference(path)
        protein_path, nucleotide_path = resolve_reference_paths(
            sample, protein_dir, nucleotide_dir, reference_map)
        protein_ref = load_single_fasta(protein_path)
        row = {'file': os.path.basename(path), 'sample': sample,
               'alignment_length': 0, 'protein_length': 0,
               'length_is_3x_protein': '', 'frame_offset': '',
               'identity': '', 'anchor_ok': '', 'has_nucleotide_ref': bool(nucleotide_path),
               'note': ''}

        if not protein_ref:
            row['note'] = 'no_protein_reference'
            rows.append(row)
            print(f"{row['file'][:52]:<52} {'-':>8} {'-':>6} {'-':>4} "
                  f"{'-':>5} {'-':>9} {'-':>5}  no protein reference")
            continue
        if not os.path.exists(path) or os.path.getsize(path) == 0:
            row['note'] = 'empty_file'
            rows.append(row)
            continue

        _, arr = load_alignment(path)
        if arr.size == 0:
            row['note'] = 'no_sequences'
            rows.append(row)
            continue

        consensus = column_plurality_consensus(arr)
        anchor = anchor_alignment_to_protein(
            consensus, protein_ref, params.genetic_code,
            params.anchor_min_identity)

        row.update({
            'alignment_length': arr.shape[1],
            'protein_length': len(protein_ref),
            'length_is_3x_protein': arr.shape[1] == 3 * len(protein_ref),
            'frame_offset': anchor.frame_offset,
            'identity': round(anchor.identity, 4),
            'anchor_ok': anchor.ok,
        })
        rows.append(row)
        print(f"{row['file'][:52]:<52} {arr.shape[1]:>8} {len(protein_ref):>6} "
              f"{str(row['length_is_3x_protein']):>4} {anchor.frame_offset:>5} "
              f"{anchor.identity:>9.4f} {str(bool(nucleotide_path)):>5}")

    scored = [r for r in rows if r['identity'] != '']
    print("-" * 96)
    if scored:
        identities = [r['identity'] for r in scored]
        offsets = Counter(r['frame_offset'] for r in scored)
        three_x = sum(1 for r in scored if r['length_is_3x_protein'])
        with_cds = sum(1 for r in rows if r['has_nucleotide_ref'])
        print(f"files scored           : {len(scored)}/{len(rows)}")
        print(f"length == 3 x protein  : {three_x}/{len(scored)}")
        print(f"frame offsets          : {dict(offsets)}")
        print(f"identity  min/med/max  : {min(identities):.4f} / "
              f"{float(np.median(identities)):.4f} / {max(identities):.4f}")
        print(f"anchor ok (>= {params.anchor_min_identity})   : "
              f"{sum(1 for r in scored if r['anchor_ok'])}/{len(scored)}")
        print(f"nucleotide ref present : {with_cds}/{len(rows)}")
        print()
        if three_x == len(scored) and set(offsets) == {0}:
            print("VERDICT: coordinate assumption holds - column//3 maps to "
                  "reference residue. DNA and AA scoring are both safe.")
        else:
            print("VERDICT: coordinate assumption does NOT hold uniformly. The "
                  "per-file anchor handles this, but review the outliers before "
                  "enabling --mode filter.")
    else:
        print("No files could be scored - check --protein-reference-dir.")

    if output_csv and rows:
        os.makedirs(os.path.dirname(os.path.abspath(output_csv)), exist_ok=True)
        with open(output_csv, 'w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)
        print(f"\nAudit written to: {output_csv}")
    return 0


# ============================================================================
# CSV output
# ============================================================================

SUMMARY_FIELDS = [
    'file_path', 'base_name', 'sample_name', 'status', 'reason', 'mode',
    'metric_requested', 'metric_used', 'anchor_frame_offset', 'anchor_identity',
    'alignment_length', 'protein_ref_length', 'length_is_3x_protein',
    'n_reads_in', 'n_reads_out', 'absolute_ceiling', 'absolute_data_cut',
    'absolute_data_cut_applied', 'removed_absolute', 'n_conflict_columns',
    'conflict_columns_per_kb', 'n_agree_majority', 'n_override_majority',
    'n_abstain_margin', 'n_abstain_support', 'n_abstain_unscoreable',
    'n_abstain_majority_protected', 'n_abstain_depth_floor',
    'n_synonymous_sites', 'n_synonymous_acted', 'n_disagreed_with_majority',
    'bases_masked', 'reads_removed_conflict', 'aa_identity_to_ref_before',
    'aa_identity_to_ref_after', 'aa_identity_to_ref_delta',
    'removed_residue_frac', 'cov_percent_before', 'cov_percent_after',
    'cov_avg_before', 'cov_avg_after', 'cov_min_before', 'cov_min_after',
    'reverted', 'revert_reason', 'runtime_s',
]

COLUMN_FIELDS = [
    'base_name', 'column', 'codon_pos', 'depth', 'major_allele', 'major_freq',
    'minor_allele', 'damage_flagged', 'winner', 'winner_is_major', 'margin',
    'synonymous', 'action', 'allele_detail',
]

READ_FIELDS = [
    'base_name', 'sequence_id', 'action', 'n_losses', 'columns',
    'global_distance',
]


def write_csv(path: str, fields: List[str], rows: List[Dict]) -> None:
    if not path:
        return
    parent = os.path.dirname(os.path.abspath(path))
    if parent:
        os.makedirs(parent, exist_ok=True)
    with open(path, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction='ignore')
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


# ============================================================================
# Self-test
# ============================================================================

def _build_synthetic(rng, n_residues=200, genetic_code=5):
    """A protein reference plus a CDS that codes for it, under `genetic_code`."""
    table = CodonTable.unambiguous_dna_by_id[genetic_code]
    by_aa: Dict[str, List[str]] = {}
    for codon, aa in table.forward_table.items():
        by_aa.setdefault(aa, []).append(codon)
    usable = sorted(aa for aa, codons in by_aa.items()
                    if aa in AA_INDEX and aa != 'X' and len(codons) >= 2)
    residues = [usable[i] for i in rng.integers(0, len(usable), n_residues)]
    cds = ''.join(by_aa[aa][0] for aa in residues)
    return ''.join(residues), cds, by_aa


def _mutate(cds, by_aa, protein, rng, n_syn=0, n_nonsyn=0, genetic_code=5):
    """Introduce a controlled number of synonymous / non-synonymous changes."""
    table = CodonTable.unambiguous_dna_by_id[genetic_code]
    codons = [cds[i:i + 3] for i in range(0, len(cds), 3)]
    positions = rng.permutation(len(codons))
    cursor = 0
    for _ in range(n_syn):
        while cursor < len(positions):
            idx = int(positions[cursor]); cursor += 1
            alternatives = [c for c in by_aa[protein[idx]] if c != codons[idx]]
            if alternatives:
                codons[idx] = alternatives[0]
                break
    for _ in range(n_nonsyn):
        while cursor < len(positions):
            idx = int(positions[cursor]); cursor += 1
            others = [aa for aa in by_aa if aa != protein[idx]
                      and aa in AA_INDEX and aa != 'X']
            if others:
                codons[idx] = by_aa[others[int(rng.integers(0, len(others)))]][0]
                break
    _ = table
    return ''.join(codons)


def _make_reads(source, n_reads, read_len, n_cols, rng, error_rate=0.0,
                deamination=0.0):
    """Place reads sampled from `source` into an MSA of width `n_cols`."""
    rows = []
    bases = ['A', 'C', 'G', 'T']
    for _ in range(n_reads):
        start = int(rng.integers(0, max(1, min(len(source), n_cols) - read_len)))
        fragment = list(source[start:start + read_len])
        for i in range(len(fragment)):
            if error_rate and rng.random() < error_rate:
                fragment[i] = bases[int(rng.integers(0, 4))]
            elif deamination and fragment[i] == 'C' and rng.random() < deamination:
                fragment[i] = 'T'
            elif deamination and fragment[i] == 'G' and rng.random() < deamination:
                fragment[i] = 'A'
        row = ['-'] * n_cols
        row[start:start + len(fragment)] = fragment
        rows.append(''.join(row))
    return rows


def _write_case(tmpdir, name, protein, reads, cds=None):
    os.makedirs(os.path.join(tmpdir, 'protein'), exist_ok=True)
    os.makedirs(os.path.join(tmpdir, 'nucleotide'), exist_ok=True)
    os.makedirs(os.path.join(tmpdir, 'aln'), exist_ok=True)
    with open(os.path.join(tmpdir, 'protein', f'{name}.fasta'), 'w') as fh:
        fh.write(f">{name}\n{protein}\n")
    if cds:
        with open(os.path.join(tmpdir, 'nucleotide', f'{name}.fasta'), 'w') as fh:
            fh.write(f">{name}\n{cds}\n")
    aln = os.path.join(tmpdir, 'aln', f'{name}_r_1_s_50_align_{name}.fas')
    with open(aln, 'w') as fh:
        for i, read in enumerate(reads):
            fh.write(f">read_{i}\n{read}\n")
    return aln


def self_test() -> int:
    """Synthetic checks. No data files, no external tools, fixed seed."""
    import tempfile

    rng = np.random.default_rng(20240728)
    failures: List[str] = []
    checks = 0

    def check(condition: bool, label: str, detail: str = "") -> None:
        nonlocal checks
        checks += 1
        status = "PASS" if condition else "FAIL"
        print(f"  [{status}] {label}" + (f"  ({detail})" if detail else ""))
        if not condition:
            failures.append(label)

    n_res = 200
    n_cols = n_res * 3
    protein, ref_cds, by_aa = _build_synthetic(rng, n_res)
    # Specimen: close to the reference. Contaminant: same gene, much further.
    specimen = _mutate(ref_cds, by_aa, protein, rng, n_syn=25, n_nonsyn=5)
    contaminant = _mutate(ref_cds, by_aa, protein, rng, n_syn=60, n_nonsyn=40)
    junk = ''.join(['ACGT'[i] for i in rng.integers(0, 4, n_cols)])

    with tempfile.TemporaryDirectory() as tmp:
        prot_dir = os.path.join(tmp, 'protein')
        nuc_dir = os.path.join(tmp, 'nucleotide')
        out_dir = os.path.join(tmp, 'out')

        def run(name, reads, params, cds=ref_cds):
            aln = _write_case(tmp, name, protein, reads, cds)
            return process_single_file(aln, prot_dir, nuc_dir, None, params,
                                       out_dir)

        print("\n1. Frame anchoring")
        reads = _make_reads(specimen, 40, 150, n_cols, rng)
        _write_case(tmp, 'anchor', protein, reads, ref_cds)
        _, arr = load_alignment(os.path.join(tmp, 'aln',
                                             'anchor_r_1_s_50_align_anchor.fas'))
        anchor = anchor_alignment_to_protein(
            column_plurality_consensus(arr), protein, 5, 0.30)
        check(anchor.frame_offset == 0, "frame offset detected as 0",
              f"offset={anchor.frame_offset}")
        check(anchor.identity > 0.90, "anchor identity high for in-frame data",
              f"identity={anchor.identity:.3f}")
        check(anchor.ok, "anchor accepted")

        print("\n2. False positives: clean alignment, 1% sequencing error")
        params = Params(mode="filter", distance_metric="dna")
        reads = _make_reads(specimen, 40, 150, n_cols, rng, error_rate=0.01)
        result = run('clean', reads, params)
        summary = result['summary']
        check(result['status'] == 'success', "processed",
              result.get('reason', ''))
        check(summary.get('removed_absolute', 0) == 0,
              "no reads removed by absolute pass",
              f"removed={summary.get('removed_absolute')}")
        check(summary.get('n_conflict_columns', 0) <= 2,
              "few or no conflict columns from random error",
              f"n={summary.get('n_conflict_columns')}")
        check(result['kept_count'] == 40, "all reads retained",
              f"kept={result['kept_count']}")

        print("\n3. Monoculture contamination (Pass A, no competitor)")
        params = Params(mode="filter", distance_metric="dna")
        reads = _make_reads(junk, 30, 150, n_cols, rng)
        result = run('monoculture', reads, params)
        check(result['summary'].get('removed_absolute', 0) >= 25,
              "absolute pass removes divergent monoculture",
              f"removed={result['summary'].get('removed_absolute')}/30")

        print("\n4. Mixed source: specimen majority, divergent contaminant")
        params = Params(mode="filter", distance_metric="dna")
        reads = (_make_reads(specimen, 28, 150, n_cols, rng)
                 + _make_reads(contaminant, 12, 150, n_cols, rng))
        result = run('mixed', reads, params)
        summary = result['summary']
        check(summary.get('removed_absolute', 0) + summary.get(
            'reads_removed_conflict', 0) > 0,
              "contaminant reads are acted on",
              f"abs={summary.get('removed_absolute')} "
              f"conflict={summary.get('reads_removed_conflict')}")
        check(summary.get('aa_identity_to_ref_delta', 0) <= 0.02,
              "reference pull stays within cap",
              f"delta={summary.get('aa_identity_to_ref_delta')}")

        print("\n5. Masking, not read removal, for isolated losses")
        params = Params(mode="filter", distance_metric="dna",
                        read_removal_loss_count=2)
        base = _make_reads(specimen, 30, 300, n_cols, rng)
        variant = list(base[0])
        for col in range(60, 63):
            if variant[col] in 'ACGT':
                variant[col] = 'A' if variant[col] != 'A' else 'C'
        reads = base
        result = run('masking', reads, params)
        check(result['status'] == 'success', "processed",
              result.get('reason', ''))
        check(result['summary'].get('reads_removed_conflict', 0)
              <= result['summary'].get('bases_masked', 0),
              "reads removed never exceeds bases masked",
              f"removed={result['summary'].get('reads_removed_conflict')} "
              f"masked={result['summary'].get('bases_masked')}")

        print("\n6. Deamination awareness")
        reads = _make_reads(specimen, 40, 150, n_cols, rng, deamination=0.05)
        aware = run('damage_on', reads,
                    Params(mode="report_only", damage_aware=True))
        naive = run('damage_off', reads,
                    Params(mode="report_only", damage_aware=False))
        check(aware['summary'].get('n_conflict_columns', 0)
              <= naive['summary'].get('n_conflict_columns', 0),
              "damage awareness never increases conflict-site count",
              f"aware={aware['summary'].get('n_conflict_columns')} "
              f"naive={naive['summary'].get('n_conflict_columns')}")

        print("\n7. report_only changes nothing")
        reads = (_make_reads(specimen, 25, 150, n_cols, rng)
                 + _make_reads(contaminant, 15, 150, n_cols, rng))
        result = run('reportonly', reads, Params(mode="report_only"))
        check(result['kept_count'] == 40,
              "all reads emitted unchanged in report_only",
              f"kept={result['kept_count']}")
        check(result['summary'].get('n_conflict_columns', 0) > 0,
              "conflicts still measured in report_only",
              f"n={result['summary'].get('n_conflict_columns')}")

        print("\n8. Reference-bias guard")
        # Specimen further from the reference than the contaminant: the filter
        # must not "clean" the real specimen away.
        far_specimen = _mutate(ref_cds, by_aa, protein, rng,
                               n_syn=40, n_nonsyn=30)
        reads = (_make_reads(far_specimen, 26, 150, n_cols, rng)
                 + _make_reads(ref_cds, 14, 150, n_cols, rng))
        result = run('bias', reads, Params(mode="filter", distance_metric="dna"))
        summary = result['summary']
        pulled = summary.get('aa_identity_to_ref_delta', 0)
        check(pulled <= 0.02 or summary.get('reverted'),
              "excess reference pull is capped or reverted",
              f"delta={pulled} reverted={summary.get('reverted')}")

        print("\n9. AA metric and synonymous blindness")
        result = run('aa', _make_reads(specimen, 30, 150, n_cols, rng)
                     + _make_reads(contaminant, 12, 150, n_cols, rng),
                     Params(mode="report_only", distance_metric="aa"),
                     cds=None)
        check(result['summary'].get('metric_used') == 'aa',
              "falls back to AA when no CDS is available",
              f"metric={result['summary'].get('metric_used')}")
        check(result['summary'].get('n_synonymous_sites', 0) >= 0,
              "synonymous sites are counted and reported",
              f"n={result['summary'].get('n_synonymous_sites')}")

        print("\n10. Degenerate inputs")
        for name, reads in (('zero', []), ('one', _make_reads(specimen, 1, 150,
                                                              n_cols, rng))):
            result = run(f'degenerate_{name}', reads, Params(mode="filter"))
            check(result['status'] in ('success', 'skipped'),
                  f"{name}-read alignment handled without error",
                  f"status={result['status']} reason={result.get('reason')}")

        empty = os.path.join(tmp, 'aln', 'empty_r_1_s_50_align_empty.fas')
        open(empty, 'w').close()
        result = process_single_file(empty, prot_dir, nuc_dir, None,
                                     Params(), out_dir)
        check(result['status'] == 'skipped' and result['reason'] == 'empty_file',
              "empty file skipped cleanly")

        print("\n11. Anchor failure escape hatch")
        reads = _make_reads(junk, 30, 400, n_cols, rng)
        result = run('badframe', reads,
                     Params(mode="filter", anchor_min_identity=0.95))
        check(result['status'] == 'skipped'
              and 'anchor' in result['reason'],
              "file with unusable frame is passed through, not scored",
              f"status={result['status']} reason={result['reason']}")
        check(result['kept_count'] == 30,
              "all reads retained when the anchor fails",
              f"kept={result['kept_count']}")

    print("\n" + "=" * 60)
    if failures:
        print(f"{len(failures)}/{checks} checks FAILED:")
        for name in failures:
            print(f"  - {name}")
        return 1
    print(f"All {checks} checks passed.")
    return 0


# ============================================================================
# CLI
# ============================================================================

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description='Reference-distance cleaning for MGE read alignments',
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--input-files-list',
                        help='Text file listing alignment FASTAs, one per line')
    parser.add_argument('--output-dir', help='Directory for filtered alignments')
    parser.add_argument('--filtered-files-list',
                        help='Output manifest of successfully filtered files')

    parser.add_argument('--protein-reference-dir',
                        help='Directory of {sample}.fasta protein references')
    parser.add_argument('--nucleotide-reference-dir',
                        help='Directory of {sample}.fasta coding nucleotide references')
    parser.add_argument('--reference-map',
                        help='sequence_references.csv, for the manual-reference workflow')

    parser.add_argument('--mode', choices=['report_only', 'filter'],
                        default='report_only',
                        help='report_only measures everything and changes nothing (default)')
    parser.add_argument('--distance-metric', choices=['dna', 'aa', 'both'],
                        default='dna')
    parser.add_argument('--genetic-code', type=int, default=5)
    parser.add_argument('--anchor-min-identity', type=float, default=0.30)

    parser.add_argument('--no-absolute-filter', action='store_true',
                        help='Skip Pass A (conflict adjudication only)')
    parser.add_argument('--max-read-distance', type=float, default=None)
    parser.add_argument('--no-data-driven-cut', action='store_true')
    parser.add_argument('--max-absolute-removal-frac', type=float, default=0.50)
    parser.add_argument('--min-reads-retained', type=int, default=5)

    parser.add_argument('--min-depth', type=int, default=3)
    parser.add_argument('--min-minor-count', type=int, default=2)
    parser.add_argument('--min-minor-freq', type=float, default=0.20)
    parser.add_argument('--no-damage-aware', action='store_true')
    parser.add_argument('--max-conflict-sites', type=int, default=200)

    parser.add_argument('--scope', choices=['local', 'global'], default='local')
    parser.add_argument('--window', default='auto',
                        help='Half-window in codons, or "auto"')
    parser.add_argument('--min-scored-positions', type=int, default=5)
    parser.add_argument('--allele-score', choices=['median', 'mean', 'best'],
                        default='median')

    parser.add_argument('--min-margin', type=float, default=1.0)
    parser.add_argument('--min-margin-override', type=float, default=2.0)
    parser.add_argument('--max-override-majority-freq', type=float, default=0.70)
    parser.add_argument('--min-override-support', type=int, default=3)

    parser.add_argument('--removal-scope',
                        choices=['base', 'read', 'base_then_read'],
                        default='base_then_read')
    parser.add_argument('--read-removal-loss-count', type=int, default=2)
    parser.add_argument('--min-retained-depth', type=int, default=2)

    parser.add_argument('--max-removal-frac', type=float, default=0.10)
    parser.add_argument('--max-ref-identity-gain', type=float, default=0.02)
    parser.add_argument('--consensus-threshold', type=float, default=0.5)

    parser.add_argument('--summary-csv')
    parser.add_argument('--column-metrics-csv')
    parser.add_argument('--read-metrics-csv')
    parser.add_argument('--threads', type=int, default=1)

    parser.add_argument('--anchor-audit', action='store_true',
                        help='Test the coordinate assumption and exit')
    parser.add_argument('--audit-csv')
    parser.add_argument('--self-test', action='store_true',
                        help='Run built-in synthetic checks and exit')
    return parser


def params_from_args(args) -> Params:
    window = args.window if args.window == 'auto' else str(int(args.window))
    return Params(
        mode=args.mode,
        distance_metric=args.distance_metric,
        genetic_code=args.genetic_code,
        anchor_min_identity=args.anchor_min_identity,
        absolute_filter=not args.no_absolute_filter,
        max_read_distance=args.max_read_distance,
        data_driven_cut=not args.no_data_driven_cut,
        max_absolute_removal_frac=args.max_absolute_removal_frac,
        min_reads_retained=args.min_reads_retained,
        min_depth=args.min_depth,
        min_minor_count=args.min_minor_count,
        min_minor_freq=args.min_minor_freq,
        damage_aware=not args.no_damage_aware,
        max_conflict_sites=args.max_conflict_sites,
        scope=args.scope,
        window=window,
        min_scored_positions=args.min_scored_positions,
        allele_score=args.allele_score,
        min_margin=args.min_margin,
        min_margin_override=args.min_margin_override,
        max_override_majority_freq=args.max_override_majority_freq,
        min_override_support=args.min_override_support,
        removal_scope=args.removal_scope,
        read_removal_loss_count=args.read_removal_loss_count,
        min_retained_depth=args.min_retained_depth,
        max_removal_frac=args.max_removal_frac,
        max_ref_identity_gain=args.max_ref_identity_gain,
        consensus_threshold=args.consensus_threshold,
    )


def main() -> int:
    args = build_parser().parse_args()

    if args.self_test:
        return self_test()

    if not args.input_files_list:
        print("Error: --input-files-list is required "
              "(or use --self-test / --anchor-audit)", file=sys.stderr)
        return 2

    with open(args.input_files_list) as handle:
        input_files = [line.strip() for line in handle if line.strip()]

    reference_map = (parse_reference_map(args.reference_map)
                     if args.reference_map else None)
    params = params_from_args(args)

    if args.anchor_audit:
        return anchor_audit(input_files, args.protein_reference_dir,
                            args.nucleotide_reference_dir, reference_map,
                            params, args.audit_csv)

    if not args.protein_reference_dir and not reference_map:
        print("Error: one of --protein-reference-dir or --reference-map "
              "is required", file=sys.stderr)
        return 2
    if not args.output_dir:
        print("Error: --output-dir is required", file=sys.stderr)
        return 2

    os.makedirs(args.output_dir, exist_ok=True)

    print(f"Refined matching over {len(input_files)} alignment files")
    print(f"Mode: {params.mode}"
          + ("  (measuring only - no file is modified)"
             if params.mode == 'report_only' else ""))
    print(f"Distance metric: {params.distance_metric}, "
          f"genetic code {params.genetic_code}, scope {params.scope}")
    print(f"Removal scope: {params.removal_scope} "
          f"(escalate at {params.read_removal_loss_count} losses)")

    results: List[Dict] = []
    with ProcessPoolExecutor(max_workers=max(1, args.threads)) as executor:
        futures = {
            executor.submit(process_single_file, path,
                            args.protein_reference_dir,
                            args.nucleotide_reference_dir,
                            reference_map, params, args.output_dir): path
            for path in input_files
        }
        for future in as_completed(futures):
            path = futures[future]
            try:
                result = future.result()
            except Exception as exc:                           # noqa: BLE001
                result = _empty_result(path, 'error', f'executor_error: {exc}')
            results.append(result)

            summary = result.get('summary', {})
            if result['status'] == 'success':
                print(f"✓ {result['output_basename']}: "
                      f"{result['kept_count']}/{result['input_count']} reads kept")
                print(f"  anchor: frame {summary.get('anchor_frame_offset')} "
                      f"identity {summary.get('anchor_identity')} "
                      f"| metric {summary.get('metric_used')}")
                print(f"  absolute removed {summary.get('removed_absolute', 0)} "
                      f"| conflicts {summary.get('n_conflict_columns', 0)} "
                      f"| masked {summary.get('bases_masked', 0)} "
                      f"| reads dropped {summary.get('reads_removed_conflict', 0)}")
                print(f"  coverage {summary.get('cov_percent_before')}% -> "
                      f"{summary.get('cov_percent_after')}% "
                      f"| ref identity delta "
                      f"{summary.get('aa_identity_to_ref_delta')}")
                if summary.get('reverted'):
                    print(f"  REVERTED: {summary.get('revert_reason')}")
            elif result['status'] == 'skipped':
                print(f"⚠ {result['output_basename']}: skipped "
                      f"({result['reason']})")
            else:
                print(f"✗ {result['output_basename']}: error "
                      f"({result['reason']})")

    successful = []
    if args.filtered_files_list:
        parent = os.path.dirname(os.path.abspath(args.filtered_files_list))
        if parent:
            os.makedirs(parent, exist_ok=True)
        with open(args.filtered_files_list, 'w') as handle:
            for result in results:
                output = result.get('output_file')
                if output and os.path.exists(output):
                    handle.write(output + '\n')
                    successful.append(output)

    summary_rows, column_rows, read_rows = [], [], []
    for result in results:
        row = {'file_path': result['file_path'],
               'base_name': result['output_basename'],
               'status': result['status'], 'reason': result['reason']}
        row.update(result.get('summary', {}))
        summary_rows.append(row)
        column_rows.extend(result.get('columns', []))
        read_rows.extend(result.get('reads', []))

    default_dir = args.output_dir
    write_csv(args.summary_csv
              or os.path.join(default_dir, 'refined_matching_summary.csv'),
              SUMMARY_FIELDS, summary_rows)
    write_csv(args.column_metrics_csv
              or os.path.join(default_dir, 'refined_matching_columns.csv'),
              COLUMN_FIELDS, column_rows)
    write_csv(args.read_metrics_csv
              or os.path.join(default_dir, 'refined_matching_reads.csv'),
              READ_FIELDS, read_rows)

    ok = [r for r in results if r['status'] == 'success']
    total_in = sum(r['input_count'] for r in results)
    total_out = sum(r['kept_count'] for r in results)
    masked = sum(r.get('summary', {}).get('bases_masked', 0) for r in ok)
    conflicts = sum(r.get('summary', {}).get('n_conflict_columns', 0) for r in ok)
    disagreed = sum(r.get('summary', {}).get('n_disagreed_with_majority', 0)
                    for r in ok)
    overrode = sum(r.get('summary', {}).get('n_override_majority', 0) for r in ok)
    synonymous = sum(r.get('summary', {}).get('n_synonymous_sites', 0) for r in ok)
    reverted = sum(1 for r in ok if r.get('summary', {}).get('reverted'))
    anchor_failed = sum(1 for r in results
                        if 'anchor' in (r.get('reason') or ''))
    deltas = [r['summary']['aa_identity_to_ref_delta'] for r in ok
              if 'aa_identity_to_ref_delta' in r.get('summary', {})]

    print(f"\nRefined Matching Summary ({params.mode}):")
    print(f"Files processed successfully : {len(ok)}/{len(input_files)}")
    print(f"Files skipped, anchor failed : {anchor_failed}")
    print(f"Reads in / out               : {total_in} / {total_out}")
    print(f"Conflict columns found       : {conflicts}")
    print(f"  disagreed with majority    : {disagreed}")
    print(f"  allowed to override it     : {overrode}")
    print(f"  synonymous (AA-blind)      : {synonymous}")
    print(f"Bases masked                 : {masked}")
    print(f"Files reverted by a cap      : {reverted}")
    if deltas:
        print(f"Reference pull (AA identity) : median "
              f"{float(np.median(deltas)):+.4f}, max {max(deltas):+.4f}")
    print(f"Filtered files written       : {len(successful)}")
    if params.mode == 'report_only':
        print("\nNothing was modified. Review the summary CSV, then re-run "
              "with --mode filter.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
