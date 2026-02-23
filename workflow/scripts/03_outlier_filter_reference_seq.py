#!/usr/bin/env python3
"""
Reference Sequence Outlier Filter

An alternative to 03_statistical_outlier_filter.py that uses the per-sample protein
reference sequence retrieved by gene fetch (stored in 02_references/protein/{sample}.fasta)
instead of generating a consensus from the input reads themselves.

This approach is better suited to detecting non-target sequences (e.g. fungal reads
in an insect dataset) because the reference is an external, trusted anchor rather than
a consensus that can itself be skewed by contamination.

Algorithm:
1. Load the protein reference for each sample from reference_dir/{sample}.fasta
2. Translate each nucleotide read in all three forward reading frames
3. Select the reading frame with the fewest internal stop codons
4. Calculate amino acid identity between the translated read and the protein reference
   using local alignment (PairwiseAligner, match=1 / mismatch=0 / linear gap)
5. Convert identity to deviation score  (deviation = 1.0 - identity)
6. Apply a percentile-based threshold to the distribution of deviation scores:
   sequences whose deviation exceeds the Nth percentile are flagged as outliers
   (default N = 90, meaning the most-divergent 10 % of reads are candidates)
7. Write retained sequences to *_outlier_filtered.fasta and emit two CSVs:
   - summary CSV   : one row per input file with aggregate statistics
   - metrics CSV   : one row per removed sequence with per-read details

Key differences from 03_statistical_outlier_filter.py:
- No consensus generation from the reads themselves
- Comparison is protein-level (translated read vs amino-acid reference)
- Deviation is measured as distance from the external reference, not from a
  self-derived consensus that could already be contaminated
- Adds --reference-dir and --genetic-code arguments; drops --consensus-threshold
"""

import os
import sys
import csv
import argparse
import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Dict, List, Optional, Tuple
from datetime import datetime


# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

def log_message(message: str, log_file=None, stdout: bool = False) -> None:
    """Timestamp a message and optionally write it to a log file / stdout."""
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    formatted = f"[{timestamp}] {message}"
    if log_file:
        log_file.write(formatted + '\n')
        log_file.flush()
    if stdout:
        print(formatted, flush=True)


# ---------------------------------------------------------------------------
# File / sample name helpers  (mirror logic from 04_reference_filter.py)
# ---------------------------------------------------------------------------

def get_sample_name(filepath: str) -> str:
    """
    Extract the sample identifier from an input FASTA path.

    MGE output files follow the pattern:
        {sample}_r_{r_value}_s_{s_value}_{sample}[_suffix].fasta
    We split on '_r_' to recover the sample ID reliably.
    Fallback: strip known filter suffixes.
    """
    basename = os.path.basename(filepath)
    for ext in ('.fasta', '.fas', '.fa'):
        if basename.lower().endswith(ext):
            basename = basename[: -len(ext)]
            break

    if '_r_' in basename:
        return basename.split('_r_')[0]

    # Fallback for files that don't follow the _r_/_s_ pattern
    for suffix in ('_outlier_filtered', '_at_filtered', '_human_filtered', '_align'):
        basename = basename.replace(suffix, '')

    return basename


def get_output_basename(filepath: str) -> str:
    """
    Return the basename used for the output file, preserving the _r_/_s_ pattern
    but stripping any upstream filter suffix so the chain reads cleanly.
    """
    basename = os.path.basename(filepath)
    for ext in ('.fasta', '.fas', '.fa'):
        if basename.lower().endswith(ext):
            basename = basename[: -len(ext)]
            break

    for suffix in ('_outlier_filtered', '_at_filtered', '_human_filtered'):
        if basename.endswith(suffix):
            basename = basename[: -len(suffix)]
            break

    return basename


# ---------------------------------------------------------------------------
# Reference loading
# ---------------------------------------------------------------------------

def load_protein_reference(reference_dir: str, sample_name: str) -> Optional[str]:
    """
    Load the protein (amino-acid) reference sequence for *sample_name* from
    ``reference_dir/{sample_name}.fasta``.

    If the file contains multiple records the longest is used, since gene fetch
    may occasionally write a multi-record file when several equally good
    references are retrieved.

    Returns None if the reference file is missing or empty.
    """
    ref_path = os.path.join(reference_dir, f"{sample_name}.fasta")
    if not os.path.exists(ref_path):
        return None

    try:
        records = list(SeqIO.parse(ref_path, "fasta"))
    except Exception:
        return None

    if not records:
        return None

    # Use the longest record when multiple are present
    best = max(records, key=lambda r: len(r.seq))
    return str(best.seq).upper()


# ---------------------------------------------------------------------------
# Translation
# ---------------------------------------------------------------------------

def translate_best_frame(nucleotide_seq: str, genetic_code: int = 5) -> Tuple[str, int]:
    """
    Translate *nucleotide_seq* in all three forward reading frames and return
    the translation that has the fewest internal stop codons together with the
    frame index (0, 1, or 2).

    The terminal stop codon (if present) is stripped before counting so that
    genuine coding sequences are not penalised.

    Gaps are removed before translation because the reads may still carry
    alignment gap characters from earlier pipeline steps.
    """
    seq = nucleotide_seq.upper().replace('-', '')

    best_protein: str = ''
    best_frame: int = 0
    min_stops: float = float('inf')

    for frame in range(3):
        subseq = seq[frame:]
        remainder = len(subseq) % 3
        if remainder:
            subseq += 'N' * (3 - remainder)
        if len(subseq) < 3:
            continue
        try:
            protein = str(Seq(subseq).translate(table=genetic_code))
            # Strip terminal stop before counting internal stops
            if protein.endswith('*'):
                protein_body = protein[:-1]
            else:
                protein_body = protein
            internal_stops = protein_body.count('*')
            if internal_stops < min_stops:
                min_stops = internal_stops
                best_protein = protein_body  # keep without terminal stop
                best_frame = frame
        except Exception:
            continue

    return best_protein, best_frame


# ---------------------------------------------------------------------------
# Similarity scoring
# ---------------------------------------------------------------------------

def calculate_protein_identity(query_protein: str, reference_protein: str) -> float:
    """
    Calculate the amino-acid identity between *query_protein* and
    *reference_protein* using local pairwise alignment.

    Scoring scheme (match=1, mismatch=0) means the raw score equals the number
    of identical aligned residues.  Normalising by the shorter sequence length
    gives a fraction in [0, 1] regardless of length differences.

    Stop-codon characters ('*') and fully ambiguous residues ('X') are removed
    before alignment because they carry no discriminatory signal.

    Returns 0.0 when either cleaned sequence is empty.
    """
    query_clean = query_protein.replace('*', '').replace('X', '')
    ref_clean = reference_protein.replace('*', '').replace('X', '')

    if not query_clean or not ref_clean:
        return 0.0

    aligner = PairwiseAligner()
    aligner.mode = 'local'
    aligner.match_score = 1
    aligner.mismatch_score = 0
    aligner.open_gap_score = -1
    aligner.extend_gap_score = -0.1

    try:
        raw_score = aligner.score(query_clean, ref_clean)
        normaliser = min(len(query_clean), len(ref_clean))
        return float(raw_score) / normaliser if normaliser > 0 else 0.0
    except Exception:
        return 0.0


# ---------------------------------------------------------------------------
# Per-file processing
# ---------------------------------------------------------------------------

def process_single_file(
    file_path: str,
    reference_dir: str,
    outlier_percentile: float,
    genetic_code: int,
    output_dir: str,
) -> Dict:
    """
    Core worker function executed (potentially in parallel) for one FASTA file.

    Steps
    -----
    1. Derive sample name and load the corresponding protein reference.
    2. Translate every nucleotide read (best reading frame).
    3. Score each translated read against the reference (amino-acid identity).
    4. Compute deviation = 1 - identity; apply percentile threshold on deviations.
    5. Write kept reads to *_outlier_filtered.fasta.
    6. Return a result dict consumed by main() for CSV writing.
    """
    sample_name = 'unknown'
    output_basename = 'unknown'

    try:
        sample_name = get_sample_name(file_path)
        output_basename = get_output_basename(file_path)

        # ---- guard: empty file ----
        if os.path.getsize(file_path) == 0:
            return _skip_result(file_path, sample_name, output_basename, 'empty_file')

        # ---- load reference ----
        reference_protein = load_protein_reference(reference_dir, sample_name)
        if reference_protein is None:
            return _skip_result(file_path, sample_name, output_basename, 'no_reference_found')

        # ---- parse reads ----
        try:
            records = list(SeqIO.parse(file_path, "fasta"))
        except Exception as exc:
            return _error_result(file_path, sample_name, output_basename,
                                 f'parse_error: {exc}')

        if not records:
            return _skip_result(file_path, sample_name, output_basename, 'no_sequences')

        input_count = len(records)

        # ---- score each read ----
        scored: List[Dict] = []
        for record in records:
            nt_seq = str(record.seq).upper()
            protein, frame = translate_best_frame(nt_seq, genetic_code)
            identity = calculate_protein_identity(protein, reference_protein)
            scored.append({
                'record': record,
                'identity': identity,
                'deviation': 1.0 - identity,
                'reading_frame': frame,
                'protein_length': len(protein),
            })

        # ---- threshold ----
        # Higher deviation = lower identity = more likely non-target.
        # We mirror the original script's percentile logic: sequences whose
        # deviation exceeds the Nth percentile are treated as outliers.
        deviations = [s['deviation'] for s in scored]
        deviation_threshold = (
            float(np.percentile(deviations, outlier_percentile))
            if deviations else float('inf')
        )
        identity_threshold = 1.0 - deviation_threshold

        # ---- filter ----
        kept_records: List[SeqRecord] = []
        removed_sequences: List[Dict] = []

        for s in scored:
            if s['deviation'] > deviation_threshold:
                removed_sequences.append({
                    'sequence_id': s['record'].id,
                    'removal_reason': 'reference_outlier',
                    'identity': s['identity'],
                    'deviation': s['deviation'],
                    'identity_threshold': identity_threshold,
                    'deviation_threshold': deviation_threshold,
                    'reading_frame': s['reading_frame'],
                    'protein_length': s['protein_length'],
                })
            else:
                kept_records.append(s['record'])

        # ---- write output ----
        output_file = None
        if kept_records:
            output_file = os.path.join(output_dir, f"{output_basename}_outlier_filtered.fasta")
            with open(output_file, 'w') as handle:
                SeqIO.write(kept_records, handle, "fasta")

        identities = [s['identity'] for s in scored]
        ref_path = os.path.join(reference_dir, f"{sample_name}.fasta")

        return {
            'file_path': file_path,
            'base_name': sample_name,
            'output_basename': output_basename,
            'status': 'success',
            'reason': 'processed',
            'input_count': input_count,
            'kept_count': len(kept_records),
            'removed_count': len(removed_sequences),
            'removed_sequences': removed_sequences,
            'output_file': output_file,
            'reference_file': ref_path,
            'thresholds': {
                'deviation_threshold': deviation_threshold,
                'identity_threshold': identity_threshold,
                'outlier_percentile': outlier_percentile,
                'mean_identity': float(np.mean(identities)),
                'std_identity': float(np.std(identities)),
                'min_identity': float(np.min(identities)),
                'max_identity': float(np.max(identities)),
            },
        }

    except Exception as exc:
        return _error_result(file_path, sample_name, output_basename,
                             f'processing_error: {exc}')


# ---------------------------------------------------------------------------
# Result dict helpers
# ---------------------------------------------------------------------------

def _base_result(file_path, sample_name, output_basename) -> Dict:
    return {
        'file_path': file_path,
        'base_name': sample_name,
        'output_basename': output_basename,
        'input_count': 0,
        'kept_count': 0,
        'removed_count': 0,
        'removed_sequences': [],
    }


def _skip_result(file_path, sample_name, output_basename, reason) -> Dict:
    d = _base_result(file_path, sample_name, output_basename)
    d.update({'status': 'skipped', 'reason': reason})
    return d


def _error_result(file_path, sample_name, output_basename, reason) -> Dict:
    d = _base_result(file_path, sample_name, output_basename)
    d.update({'status': 'error', 'reason': reason})
    return d


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            'Filter sequences that are outliers relative to the gene fetch '
            'protein reference rather than to a self-generated consensus.'
        )
    )
    parser.add_argument(
        '--input-files-list', required=True,
        help='Text file listing one input FASTA path per line.',
    )
    parser.add_argument(
        '--output-dir', required=True,
        help='Directory where filtered FASTA files are written.',
    )
    parser.add_argument(
        '--filtered-files-list', required=True,
        help='Output text file listing paths of successfully filtered FASTA files.',
    )
    parser.add_argument(
        '--summary-csv', required=True,
        help='Output CSV with one row per input file (aggregate statistics).',
    )
    parser.add_argument(
        '--metrics-csv', required=True,
        help='Output CSV with one row per removed sequence (detailed metrics).',
    )
    parser.add_argument(
        '--reference-dir', required=True,
        help=(
            'Directory containing gene fetch protein reference FASTA files '
            'named {sample}.fasta  (typically 02_references/protein/).'
        ),
    )
    parser.add_argument(
        '--outlier-percentile', type=float, default=90.0,
        help=(
            'Deviation-score percentile used as the outlier cutoff.  Sequences '
            'whose deviation exceeds this percentile are removed.  '
            'Default: 90.0  (removes the most-divergent ~10%% of reads).'
        ),
    )
    parser.add_argument(
        '--genetic-code', type=int, default=5,
        help=(
            'NCBI genetic code table used for translation.  '
            'Default: 5 (invertebrate mitochondrial).  '
            'Use 1 for standard / nuclear, 2 for vertebrate mitochondrial.'
        ),
    )
    parser.add_argument(
        '--threads', type=int, default=1,
        help='Number of parallel worker processes.  Default: 1.',
    )

    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    if not os.path.exists(args.reference_dir):
        print(f"Error: reference directory not found: {args.reference_dir}", file=sys.stderr)
        sys.exit(1)

    with open(args.input_files_list, 'r') as fh:
        input_files = [line.strip() for line in fh if line.strip()]

    print(f"Processing {len(input_files)} file(s)")
    print(f"Reference directory : {args.reference_dir}")
    print(f"Outlier percentile  : {args.outlier_percentile}")
    print(f"Genetic code        : {args.genetic_code}")
    print(f"Threads             : {args.threads}")

    # ---- parallel processing ----
    results: List[Dict] = []
    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        future_map = {
            executor.submit(
                process_single_file,
                fp,
                args.reference_dir,
                args.outlier_percentile,
                args.genetic_code,
                args.output_dir,
            ): fp
            for fp in input_files
        }

        for future in as_completed(future_map):
            fp = future_map[future]
            try:
                result = future.result()
                results.append(result)

                if result['status'] == 'success':
                    thr = result.get('thresholds', {})
                    print(
                        f"  {result['base_name']}: "
                        f"{result['kept_count']}/{result['input_count']} sequences kept"
                    )
                    print(
                        f"    Reference      : {os.path.basename(result.get('reference_file', 'N/A'))}"
                    )
                    print(
                        f"    Mean identity  : {thr.get('mean_identity', 0):.4f} "
                        f"+/- {thr.get('std_identity', 0):.4f}  "
                        f"[{thr.get('min_identity', 0):.4f} – {thr.get('max_identity', 0):.4f}]"
                    )
                    print(
                        f"    Identity cutoff: {thr.get('identity_threshold', 0):.4f}  "
                        f"(deviation > {thr.get('deviation_threshold', 0):.4f} removed)"
                    )
                elif result['status'] == 'skipped':
                    print(f"  {result['base_name']}: skipped ({result['reason']})")
                else:
                    print(f"  {result['base_name']}: ERROR ({result['reason']})")

            except Exception as exc:
                print(f"  {fp}: executor error – {exc}", file=sys.stderr)
                try:
                    sn = get_sample_name(fp)
                    ob = get_output_basename(fp)
                except Exception:
                    sn = ob = os.path.basename(fp)
                results.append(_error_result(fp, sn, ob, f'executor_error: {exc}'))

    # ---- write filtered files list ----
    successful_files: List[str] = []
    with open(args.filtered_files_list, 'w') as fh:
        for r in results:
            if r['status'] == 'success' and r.get('output_file') and \
                    os.path.exists(r['output_file']):
                fh.write(r['output_file'] + '\n')
                successful_files.append(r['output_file'])

    # ---- summary CSV ----
    with open(args.summary_csv, 'w', newline='') as fh:
        writer = csv.writer(fh)
        writer.writerow([
            'file_path', 'base_name', 'reference_file',
            'identity_threshold', 'deviation_threshold',
            'mean_identity', 'std_identity', 'min_identity', 'max_identity',
            'input_count', 'kept_count', 'removed_count',
        ])
        for r in results:
            thr = r.get('thresholds', {})
            writer.writerow([
                r['file_path'],
                r['base_name'],
                r.get('reference_file', ''),
                thr.get('identity_threshold', ''),
                thr.get('deviation_threshold', ''),
                thr.get('mean_identity', ''),
                thr.get('std_identity', ''),
                thr.get('min_identity', ''),
                thr.get('max_identity', ''),
                r['input_count'],
                r['kept_count'],
                r['removed_count'],
            ])

    # ---- metrics CSV (per removed sequence) ----
    with open(args.metrics_csv, 'w', newline='') as fh:
        writer = csv.writer(fh)
        writer.writerow([
            'base_name', 'sequence_id', 'removal_reason',
            'identity', 'deviation',
            'identity_threshold', 'deviation_threshold',
            'reading_frame', 'protein_length',
        ])
        for r in results:
            for seq in r['removed_sequences']:
                writer.writerow([
                    r['base_name'],
                    seq['sequence_id'],
                    seq['removal_reason'],
                    f"{seq['identity']:.6f}",
                    f"{seq['deviation']:.6f}",
                    f"{seq['identity_threshold']:.6f}",
                    f"{seq['deviation_threshold']:.6f}",
                    seq['reading_frame'],
                    seq['protein_length'],
                ])

    # ---- console summary ----
    total_in  = sum(r['input_count']   for r in results)
    total_kept = sum(r['kept_count']   for r in results)
    total_rm  = sum(r['removed_count'] for r in results)
    n_ok      = sum(1 for r in results if r['status'] == 'success')
    n_no_ref  = sum(1 for r in results if r['status'] == 'skipped'
                    and r['reason'] == 'no_reference_found')

    print("\nReference-Based Outlier Filter – Summary")
    print(f"  Files processed successfully : {n_ok}/{len(input_files)}")
    print(f"  Files skipped (no reference) : {n_no_ref}")
    print(f"  Total input sequences        : {total_in}")
    print(f"  Sequences kept               : {total_kept}")
    print(f"  Sequences removed            : {total_rm}")
    print(f"  Filtered files written       : {len(successful_files)}")


if __name__ == "__main__":
    main()
