# `{run_name}_final_metrics.csv` — column definitions

## Overview

**Location:** `{output_dir}/{run_name}_final_metrics.csv`
**Written by:** `beegees/workflow/scripts/val_csv_merger.py` (Snakefile rule `val_barcode_recovery_metrics_combine`)

The file is built by **left-joining** three tables. The key is `fasta_header` with any leading `>` removed, matched to `seq_id`:

1. `03_barcode_recovery/{run_name}_barcode_recovery_metrics.csv`: barcode recovery and cleaning stats. It sets the rows and the first 21 columns.
2. `04_barcode_validation/structural/structural_validation.csv`: barcode structural validation stats. 8 columns.
3. `04_barcode_validation/taxonomic/metrics/02_taxonomic_validation.csv`: barcode taxonomic validation stats. 13 columns. Its `length` column is renamed `aligned_length`.

### What each row represents
One row is one **consensus sequence candidate** for a single sample, where each sample has up to 24 consensus candidates. Every sample × MGE parameter combination (`r`, `s`) × preprocessing mode (`merge`/`concat` for PE, `se` for SE) gives **up to 24 rows**:
- **12 'raw' MGE rows**: the consensus taken straight from MitoGeneExtractor. `mge_params` looks like (e.g.) `r_1.3_s_50_merge` or `r_1.3_s_50_concat`.
- **12 fasta_cleaner rows**: the consensus rebuilt after the cleaning filters. `mge_params` looks like `r_1.3_s_50_fcleaner_merge` or `r_1.3_s_50_flceaner_concat`.

Some columns mean different things depending on the row type. The notes in the tables below say where.

### Missing values
- `null` (literal text) means no value was available for this row and metric type.
- An empty cell means the row had no match in the structural or taxonomic validation table. Usually that's because the sequence didn't reach that validation step (only structurally passing sequences go on to BLAST).

---

## 1. Barcode recovery and cleaning (columns 1–21)
> Source: `compile_barcoding_stats.py`

| # | Header | Definition | Type / values | Notes |
|---|---|---|---|---|
| 1 | `Filename` | Name of the barcode consensus candidate attempt (e.g. `BGSNL096-23_r_1.3_s_50_BGSNL096-23`). | text | |
| 2 | `ID` | Sample ID: The user-defined sample ID, formed of verything before the first `_` in the filename. | text | Groups all candidates for one specimen. |
| 3 | `mge_params` | MGE parameters + mode tag: `r_{r}_s_{s}_{mode}` for raw rows, `r_{r}_s_{s}_fcleaner_{mode}` for cleaner rows. `r` = MGE relative score threshold, `s` = MGE exonerate score threshold. | text | |
| 4 | `fasta_header` | Header of this consensus candidated (pre- or post-cleaning), with `>` prefix, e.g. `>BSNHM004-24_r_1.5_s_50_BSNHM004-24_merge`. This is the join key to the validation tables. | text / `null` | `null` if no consensus with that ID + params is in the FASTA (e.g. MGE produced no consensus). |
| 5 | `sample_taxid` | NCBI taxid gene_fetch used to find the reference protein. | integer / `null` | `null` when `run_gene_fetch: false` and reference sequences were manually supplied. |
| 6 | `ref_accession` | Accession of the protein reference retrieved by gene_fetch and given to MGE (i.e. `protein_accession`). | text / `null` | `null` when `run_gene_fetch: false`, as above. |
| 7 | `ref_rank` | Taxonomic rank at which gene_fetch found the reference (i.e. `matched_rank`), e.g. species, genus or family if it had to go up the lineage. | text / `null` | `null` when `run_gene_fetch: false`, as above. |
| 8 | `n_reads_in` | Reads given to MGE for this sample post-QC ("Number of input sequences" in the MGE `.out` file). | integer / `null` | Raw rows only. |
| 9 | `n_reads_aligned` | **Raw rows:** number of reads aligned to the reference by MGE. **FCleaner rows:** number of reads aligned to the reference by MGE, and which are also kept after all fasta_cleaner filtering steps (. | integer | Means different things by row 'type'. |
| 10 | `n_reads_skipped` | Number of reads 'skipped'/omitted by MGE during the read alignment stage, such as due to frameshifts, low `r` value or multiple mapping locations. | integer / `null` | 'Raw' rows only, but applied to fcleaner rows. |
| 11 | `ref_length` | Length (nucleotides; in bp) of retrieved (protein) reference sequence reads were aligned to. | integer / `null` | Raw rows only, but applied to fcleaner rows. |
| 12 | `cov_min` | Minimum read depth (x) across alignment columns. | number | Raw rows: MGE alignment. Cleaner rows: cleaned alignment. |
| 13 | `cov_max` | Maximum read depth (x) across alignment columns. | number | As above. |
| 14 | `cov_avg` | Mean read depth (x) across alignment columns. | number | As above. |
| 15 | `cov_med` | Median read depth (x) across alignment columns. | number | As above. |
| 16 | `cleaning_removed_human` | Reads removed for mapping to the supplied human reference. | integer / `null` | Cleaner rows only. |
| 17 | `cleaning_removed_at` | Reads removed because their AT content differs from the consensus AT by more than `at_difference`. | integer / `null` | Cleaner rows only. |
| 18 | `cleaning_removed_outlier` | Reads removed as statistical outliers: unweighted or conservation-weighted deviation from consensus above the `outlier_percentile` threshold. | integer / `null` | Cleaner rows only. |
| 19 | `cleaning_removed_reference` | Reads removed by BWA-MEM mapping against the sample's reference (optional step). | integer / `null` | Cleaner rows only. `0` when reference filtering is off. |
| 20 | `cleaning_ambig_bases` | Non-A/C/G/T characters in the cleaned consensus after gaps are removed. | integer / `null` | Cleaner rows only. |
| 21 | `cleaning_cov_percent` | Percentage of alignment columns covered by at least one read in the cleaned alignment. | number (0–100) / `null` | Cleaner rows only. |

---

## 2. Structural validation (columns 22–29)
> Source: `structural_validaiton.csv`, generated by `structural_validation.py`. Applies to every consensus with a matching `fasta_header`.

| # | Header | Definition | Type / values | Notes |
|---|---|---|---|---|
| 22 | `seq_id` | Sequence 'ID' (i.e. consensus candidate ID) in the structural validation output, akin to `Filename`. | text | Empty = no structural validation record. |
| 23 | `barcode_length` | Length of the barcode region taken from the consensus by nhmmer alignment to the marker HMM, rebuilt in HMM coordinates, with leading/trailing Ns trimmed. | integer (bp) | |
| 24 | `barcode_ambiguous_bases_original` | Ns already in the sequence before HMM processing. Shows genuine ambiguity in the consensus, and is a sign of a chimeric barcode. | integer | Must be 0 to pass. |
| 25 | `barcode_ambiguous_bases` | All Ns in the final, HMM-aligned barcode region: original Ns plus gaps turned into Ns. | integer | Must be < 30% of `barcode_base_count` to pass. |
| 26 | `barcode_base_count` | Informative bases in the barcode: `barcode_length − barcode_ambiguous_bases`. | integer | Must be > 300 to pass. |
| 27 | `reading_frame` | Reading frame (0, 1 or 2) with the fewest internal stop codons after translation with the configured genetic code. `-1` = no valid frame / empty barcode. | integer (−1, 0, 1, 2) | Must be ≥ 0 to pass. |
| 28 | `stop_codons` | Internal stop codons in the chosen frame. A terminal stop isn't counted. | integer | Must be 0 to pass. |
| 29 | `barcode_rank` | Quality rank, lower is better. **1:** ≥ 500 bp. **2:** 400–499. **3:** 300–399. **4:** 200–299. **5:** 1–199. Ranks 1–5 need no original Ns, no stop codons and a valid frame. **6:** has original Ns, has stop codons / no valid frame, or no bases. bp here means `barcode_base_count`. | integer 1–6 | |

**Pass criteria for structural validation** (all must hold): `barcode_ambiguous_bases_original == 0`, `stop_codons == 0`, `reading_frame >= 0`, `barcode_base_count > 300`, `barcode_ambiguous_bases < 0.3 × barcode_base_count`. Only passing sequences go on to taxonomic validation.

---

## 3. Taxonomic validation (columns 30–42)
> Source: `tv_local_blast.py` & `tv_blast2taxonomy.py`. Only sequences that passed structural validation have values here.

| # | Header | Definition | Type / values | Notes |
|---|---|---|---|---|
| 30 | `taxval_rank` | Configured validation floor: the coarsest rank a match may be accepted at, set in the pipeline config at run inception. | `species` / `genus` / `family` / `order` | Same for every row in a run (default `family`). |
| 31 | `expected_taxonomy` | Expected taxonomic lineage of the sample from the user's expected-taxonomy file, written as `phylum;class;order;family;genus;species`. Missing ranks stay as empty slots. | text | |
| 32 | `expected_taxonomy_rank` | The validation floor actually used for this sample. It's `taxval_rank` if the sample has a value at that rank; otherwise it's the coarsest rank the sample does have. | text | Differs from `taxval_rank` only when the floor was relaxed. |
| 33 | `obs_taxonomy` | Top 10 BLAST hits that pass `min_pident`/`min_length`, sorted by % identity, length, mismatches and e-value. Format: `sseqid: Taxon (rank), …; …`, listing species/genus/family/order for each hit. | text | Ignores the validation floor. The first entry is **not** necessarily the accepted hit. A repeated accession is a separate HSP. |
| 34 | `match_taxonomy` | Whether any hit (after the `min_pident`/`min_length` filters) matches the expected lineage at or below the floor. Names are compared whole, ignoring case and whitespace. | `YES` / `NO` | |
| 35 | `matched_rank` | Most specific rank at which the accepted BLAST hit matches the expected lineage (i.e. the taxonomic rank at which the expected and observed taxonomy agree). | `species`/`genus`/`family`/`order`, or `null` if no match | |
| 36 | `top_matching_hit` | Subject ID (`sseqid`) of the accepted hit: the first hit, by descending % identity, that matches at or below the floor. | text / `No hits` | If `match_taxonomy = NO`, this is the best overall hit by quality instead. |
| 37 | `pident` | BLAST % identity of the hit in `top_matching_hit`. | number (0–100) | |
| 38 | `aligned_length` | BLAST alignment length of that hit. | integer (bp) | |
| 39 | `mismatch` | BLAST mismatches for that hit. | integer | |
| 40 | `evalue` | BLAST e-value for that hit. | number | |
| 41 | `gaps` | BLAST gap openings (`gapopen`) for that hit. | integer | Not total gap positions.|
| 42 | `selected` | `YES` for the one best sequence per Sample ID among `match_taxonomy = YES` rows; `NO` otherwise. Priority: most specific `matched_rank` → fewest gaps → fewest mismatches → highest pident → lowest e-value → longest alignment → highest MGE `s` → highest MGE `r` → `fcleaner` preferred. | `YES` / `NO` | `YES` rows make up `{run_name}_validated_barcodes.fasta`. A Sample ID with no matches has no `YES` row. |
