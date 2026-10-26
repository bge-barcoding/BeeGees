"""Build config.yaml for BeeGees from Galaxy tool form parameters.

Called by the Galaxy wrapper before `beegees run`. All arguments are passed
by the Cheetah command block as resolved strings/paths.

Usage:
    python build_config_yaml.py --output config.yaml [options]
"""

import argparse
import glob
import os
import sys
import yaml


def err(stage, msg):
    print(f'[build_config_yaml] ERROR ({stage}): {msg}', file=sys.stderr)
    sys.exit(1)


def stage_blast_prefix(prefix, stage_dir='blast_db'):
    """Link a pre-built BLAST database, given as a name prefix, into stage_dir and return a path that exists."""
    files = glob.glob(glob.escape(prefix) + '.*')
    if not files:
        err('taxonomic validation', f'BLAST database not found: no files match {prefix}.*')
    os.makedirs(stage_dir, exist_ok=True)
    for f in files:
        link = os.path.join(stage_dir, os.path.basename(f))
        if not os.path.lexists(link):
            os.symlink(f, link)
    # The pipeline checks os.path.exists(blast_db); BLAST itself ignores this file and reads <name>.n*
    placeholder = os.path.join(stage_dir, os.path.basename(prefix))
    open(placeholder, 'a').close()
    print(f'[build_config_yaml] linked {len(files)} BLAST database file(s) from {prefix} into {stage_dir}/')
    return placeholder


def validate(args):
    # ── Stage 1: samples file ──
    if not os.path.isfile(args.samples_file):
        err('samples file', f'File not found: {args.samples_file}')

    # ── Stage 2: reference source ──
    if args.ref_mode == 'manual':
        if not args.sequence_reference_file:
            err('reference', 'ref_mode is "manual" but --sequence-reference-file was not provided')
        if not os.path.isfile(args.sequence_reference_file):
            err('reference', f'Sequence reference file not found: {args.sequence_reference_file}')
    else:
        if not args.gf_email:
            err('gene-fetch', '--gf-email is required when ref_mode is "gene_fetch"')
        if not args.gf_api_key:
            err('gene-fetch', '--gf-api-key is required when ref_mode is "gene_fetch"')

    # ── Stage 3: fastp adapters — both or neither ──
    if bool(args.fastp_adapter_r1) != bool(args.fastp_adapter_r2):
        err('fastp', 'Provide both --fastp-adapter-r1 and --fastp-adapter-r2, or neither')

    # ── Stage 4: barcode recovery ──
    if not args.br_r:
        err('barcode recovery', '--br-r must contain at least one value')
    if not args.br_s:
        err('barcode recovery', '--br-s must contain at least one value')
    if not 0.0 <= args.br_t <= 1.0:
        err('barcode recovery', f'--br-t must be between 0 and 1, got {args.br_t}')

    # ── Stage 5: fasta_cleaner ──
    if args.fc_reference_dir and not os.path.exists(args.fc_reference_dir):
        err('fasta_cleaner', f'Reference filter path not found: {args.fc_reference_dir}')

    # ── Stage 6: taxonomic validation ──
    if not os.path.exists(args.tv_blast_db):
        # Pre-installed databases (blastdb data table) are a name prefix, not a file
        args.tv_blast_db = stage_blast_prefix(args.tv_blast_db)
    if not os.path.isfile(args.tv_db_taxonomy):
        err('taxonomic validation', f'Taxonomy TSV not found: {args.tv_db_taxonomy}')
    if not 0 <= args.tv_min_pident <= 100:
        err('taxonomic validation', f'--tv-min-pident must be 0-100, got {args.tv_min_pident}')


# Rule resource defaults — not exposed in the Galaxy form.
RULE_DEFAULTS = {
    'gene_fetch':                 {'mem_mb': 8192,  'threads': 4,  'partition': 'long'},
    'fastp_qc':                   {'mem_mb': 16384, 'threads': 4},
    'clean_headers_merge':        {'mem_mb': 8192,  'threads': 1},
    'fastq_concat':               {'mem_mb': 8192,  'threads': 1},
    'quality_trim':               {'mem_mb': 8192,  'threads': 4},
    'downsample':                 {'mem_mb': 8192,  'threads': 4},
    'MitoGeneExtractor':          {'mem_mb': 32768, 'threads': 1,  'partition': 'himem'},
    'rename_and_combine_cons':    {'mem_mb': 8192,  'threads': 4},
    'gzip_merged_clean':          {'mem_mb': 8192,  'threads': 4},
    'human_cox1_filter':          {'mem_mb': 16384, 'threads': 4},
    'at_content_filter':          {'mem_mb': 16384, 'threads': 4},
    'statistical_outlier_filter': {'mem_mb': 16384, 'threads': 4},
    'reference_filter':           {'mem_mb': 4096,  'threads': 1},
    'consensus_generation':       {'mem_mb': 8192,  'threads': 4},
    'extract_stats_to_csv':       {'mem_mb': 8192,  'threads': 1},
    'structural_validation':      {'mem_mb': 32768, 'threads': 1,  'partition': 'himem'},
    'taxonomic_validation':       {'mem_mb': 32768, 'threads': 16, 'partition': 'himem'},
    'blast2taxonomy':             {'mem_mb': 16384, 'threads': 4},
    'download_taxdump':           {'mem_mb': 4096,  'threads': 1},
    'multiqc_plots':              {'mem_mb': 8192,  'threads': 2},
    'multiqc':                    {'mem_mb': 8192,  'threads': 1},
}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', required=True, metavar='PATH')

    # Core
    p.add_argument('--samples-file',  default='./samples.csv')
    p.add_argument('--output-dir',    default='./beegees_output')

    # Reference
    p.add_argument('--ref-mode', choices=['manual', 'gene_fetch'], default='manual')
    p.add_argument('--sequence-reference-file', default='')

    # gene-fetch (only used when --ref-mode=gene_fetch)
    p.add_argument('--gf-email',          default='')
    p.add_argument('--gf-api-key',        default='')
    p.add_argument('--gf-gene',           default='coi')
    p.add_argument('--gf-type',           default='both')
    p.add_argument('--gf-minimum-length', default='500')
    p.add_argument('--gf-input-type',     default='hierarchical')
    p.add_argument('--gf-genbank',        action='store_true')

    # fastp
    p.add_argument('--fastp-adapter-r1',    default='')
    p.add_argument('--fastp-adapter-r2',    default='')
    p.add_argument('--fastp-extra-args',    default='')

    # Downsampling
    p.add_argument('--downsampling-enabled', action='store_true')
    p.add_argument('--downsampling-max-reads', type=int, default=0)

    # Barcode recovery
    p.add_argument('--br-r', nargs='+', type=float, default=[1, 1.3, 1.5])
    p.add_argument('--br-s', nargs='+', type=int,   default=[50, 100])
    p.add_argument('--br-n', type=int,   default=0)
    p.add_argument('--br-C', type=int,   default=5)
    p.add_argument('--br-t', type=float, default=0.5)

    # fasta_cleaner
    p.add_argument('--fc-consensus-threshold', type=float, default=0.5)
    p.add_argument('--fc-human-threshold',     type=float, default=0.95)
    p.add_argument('--fc-at-difference',       type=float, default=0.1)
    p.add_argument('--fc-at-mode',             default='absolute')
    p.add_argument('--fc-outlier-percentile',  type=float, default=90.0)
    p.add_argument('--fc-reference-dir',       default=None)
    p.add_argument('--fc-reference-filter-mode', default='remove_similar')

    # Structural validation
    p.add_argument('--sv-target',  default='cox1')
    p.add_argument('--sv-verbose', action='store_true')

    # Taxonomic validation
    p.add_argument('--tv-blast-db',          required=True)
    p.add_argument('--tv-db-taxonomy',       required=True)
    p.add_argument('--tv-taxval-rank',       default='family')
    p.add_argument('--tv-verbose',           action='store_true')
    p.add_argument('--tv-min-pident',        type=int, default=80)
    p.add_argument('--tv-min-length',        type=int, default=100)

    return p.parse_args()


def build_config(args):
    # Pipeline scripts join relative output paths onto their output dir, so paths must be absolute
    samples_file = os.path.abspath(args.samples_file)
    cfg = {
        'run_name':    'galaxy_run',
        'samples_file': samples_file,
        'output_dir':  os.path.abspath(args.output_dir),
        'sequence_reference_file': (os.path.abspath(args.sequence_reference_file)
                                    if args.sequence_reference_file else ''),
        'run_gene_fetch': args.ref_mode == 'gene_fetch',
        'fastp': {
            'adapter_r1':      args.fastp_adapter_r1,
            'adapter_r2':      args.fastp_adapter_r2,
            'extra_fastp_args': args.fastp_extra_args,
        },
        'downsampling': {
            'enabled':   args.downsampling_enabled,
            'max_reads': args.downsampling_max_reads,
        },
        'barcode_recovery': {
            'r': args.br_r,
            's': args.br_s,
            'n': args.br_n,
            'C': args.br_C,
            't': args.br_t,
        },
        'fasta_cleaner': {
            'consensus_threshold':  args.fc_consensus_threshold,
            'human_threshold':      args.fc_human_threshold,
            'at_difference':        args.fc_at_difference,
            'at_mode':              args.fc_at_mode,
            'outlier_percentile':   args.fc_outlier_percentile,
            'reference_dir':        args.fc_reference_dir,
            'reference_filter_mode': args.fc_reference_filter_mode,
        },
        'structural_validation': {
            'target':  args.sv_target,
            'verbose': args.sv_verbose,
        },
        'taxonomic_validation': {
            'blast_db':           os.path.abspath(args.tv_blast_db),
            'db_taxonomy':        args.tv_db_taxonomy,
            'taxval_rank':        args.tv_taxval_rank,
            'expected_taxonomy':  samples_file,
            'verbose':            args.tv_verbose,
            'min_pident':         args.tv_min_pident,
            'min_length':         args.tv_min_length,
        },
        'rules': RULE_DEFAULTS,
    }

    if args.ref_mode == 'gene_fetch':
        cfg['gene_fetch'] = {
            'email':          args.gf_email,
            'api_key':        args.gf_api_key,
            'gene':           args.gf_gene,
            'type':           args.gf_type,
            'minimum_length': args.gf_minimum_length,
            'input_type':     args.gf_input_type,
            'genbank':        args.gf_genbank,
        }

    return cfg


def main():
    args = parse_args()
    validate(args)
    cfg = build_config(args)
    try:
        with open(args.output, 'w') as f:
            yaml.dump(cfg, f, default_flow_style=False, sort_keys=False)
    except OSError as e:
        err('output', f'Could not write {args.output}: {e}')
    print(f'[build_config_yaml] wrote config to {args.output}')


if __name__ == '__main__':
    main()
