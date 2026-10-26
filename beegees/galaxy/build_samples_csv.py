"""Resolve a user-supplied samples CSV for BeeGees.

The user provides a CSV with columns ID, forward, [reverse], the taxonomy
columns and optionally taxid, where forward/reverse are bare filenames (not
full paths). This script maps those filenames to absolute paths using the
Galaxy collection, then writes the resolved CSV that BeeGees expects. All
other columns are passed through unchanged (gene-fetch reads taxid from it).

Usage:
    python build_samples_csv.py
        --samples-csv  /path/to/user_samples.csv
        --output       samples.csv
        --collection-names  name1 name2 ...
        --collection-paths  /abs/path1 /abs/path2 ...
        [--reads-dir   reads]
        [--require-taxid]
        [--refs-csv    /path/to/user_sequence_references.csv
         --ref-names   name1 ...  --ref-paths  /abs/path1 ...
         --refs-output sequence_references.csv  --refs-dir  references]
"""

import argparse
import csv
import os
import sys

REQUIRED_COLS = {
    'ID':      ['ID'],
    'forward': ['forward'],
}
OPTIONAL_READ_COL = {
    'reverse': ['reverse'],
}
TAX_COLS = ['phylum', 'class', 'order', 'family', 'genus', 'species']


def err(stage, msg):
    print(f'[build_samples_csv] ERROR ({stage}): {msg}', file=sys.stderr)
    sys.exit(1)


def resolve_col(fieldnames, aliases):
    for a in aliases:
        if a in fieldnames:
            return a
    return None


def is_gzipped(path):
    with open(path, 'rb') as f:
        return f.read(2) == b'\x1f\x8b'


def link_read(filename, target, reads_dir):
    """Symlink a Galaxy dataset (.dat) under its real filename so tools can detect gzip by extension."""
    name = os.path.basename(filename)
    if name in ('', '.', '..'):
        err('file resolution', f'Invalid filename: "{filename}"')
    gz = is_gzipped(target)
    if gz and not name.endswith('.gz'):
        name += '.gz'
    elif not gz and name.endswith('.gz'):
        name = name[:-3]
    link = os.path.abspath(os.path.join(reads_dir, name))
    if os.path.lexists(link):
        if os.path.realpath(link) != os.path.realpath(target):
            err('file resolution', f'Two different files map to the same name: {name}')
        return link
    os.symlink(os.path.abspath(target), link)
    return link


FASTQ_EXTS = ('.fastqsanger', '.fastq', '.fq')
FASTA_EXTS = ('.fasta', '.fas', '.faa', '.fa')
REF_COLS = ['ID', 'protein_reference_path']


def strip_ext(name, exts):
    if name.endswith('.gz'):
        name = name[:-3]
    for ext in exts:
        if name.endswith(ext):
            return name[:-len(ext)]
    return name


def build_path_map(names, paths, exts):
    """Map collection element names, and their extension-less stems, to dataset paths."""
    exact = dict(zip(names, paths))
    stems = {}
    for name, path in exact.items():
        stems.setdefault(strip_ext(name, exts), []).append(path)
    return exact, stems


def resolve_file(filename, label, path_map, stem_map, exts, link_dir, collection):
    """Find filename in a collection and return a symlink to it under link_dir."""
    if filename in path_map:
        return link_read(filename, path_map[filename], link_dir)
    # Galaxy collections often drop extensions (test_R1 vs test_R1.fastq.gz)
    candidates = stem_map.get(strip_ext(filename, exts), [])
    if len(candidates) == 1:
        return link_read(filename, candidates[0], link_dir)
    if len(candidates) > 1:
        err('file resolution',
            f'{label}: "{filename}" matches more than one element of the {collection} collection; '
            'make the collection element names unique.')
    err('file resolution',
        f'{label}: "{filename}" not found in the {collection} collection.\n'
        f'  Available: {", ".join(sorted(path_map))}')


def resolve_references(args, sample_ids):
    """Write a sequence references CSV whose protein_reference_path points to collection datasets."""
    if len(args.ref_names) != len(args.ref_paths):
        err('references',
            f'--ref-names ({len(args.ref_names)}) and '
            f'--ref-paths ({len(args.ref_paths)}) have different counts')
    path_map, stem_map = build_path_map(args.ref_names, args.ref_paths, FASTA_EXTS)

    if not os.path.isfile(args.refs_csv):
        err('references', f'File not found: {args.refs_csv}')
    with open(args.refs_csv, newline='') as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames or []
        missing = [c for c in REF_COLS if c not in fieldnames]
        if missing:
            err('references', f'Missing column(s) in sequence references CSV: {", ".join(missing)}'
                f'\n  Found: {", ".join(fieldnames)}')
        rows = {r['ID']: r for r in reader}

    no_ref = [i for i in sample_ids if not rows.get(i, {}).get('protein_reference_path')]
    if no_ref:
        err('references', f'No protein_reference_path given for sample(s): {", ".join(no_ref)}')

    os.makedirs(args.refs_dir, exist_ok=True)
    out_rows = []
    for sid in sample_ids:
        r = dict(rows[sid])
        r['protein_reference_path'] = resolve_file(
            r['protein_reference_path'], f'sample {sid} protein reference',
            path_map, stem_map, FASTA_EXTS, args.refs_dir, 'protein reference')
        out_rows.append(r)

    try:
        with open(args.refs_output, 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=fieldnames)
            w.writeheader()
            w.writerows(out_rows)
    except OSError as e:
        err('references', f'Could not write {args.refs_output}: {e}')

    print(f'[build_samples_csv] resolved {len(out_rows)} protein reference(s) to {args.refs_output}')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--samples-csv',       required=True)
    parser.add_argument('--output',            required=True)
    parser.add_argument('--collection-names',  nargs='*', default=[])
    parser.add_argument('--collection-paths',  nargs='*', default=[])
    parser.add_argument('--reads-dir',         default='reads')
    # gene-fetch taxid mode reads a taxid column from the resolved CSV
    parser.add_argument('--require-taxid',     action='store_true')
    # Manual reference mode: references CSV + protein FASTA collection
    parser.add_argument('--refs-csv',          default=None)
    parser.add_argument('--refs-output',       default='sequence_references.csv')
    parser.add_argument('--ref-names',         nargs='*', default=[])
    parser.add_argument('--ref-paths',         nargs='*', default=[])
    parser.add_argument('--refs-dir',          default='references')
    args = parser.parse_args()

    # ── Stage 1: input file exists ──
    if not os.path.isfile(args.samples_csv):
        err('samples CSV', f'File not found: {args.samples_csv}')

    # ── Stage 2: collection consistency ──
    if len(args.collection_names) != len(args.collection_paths):
        err('collection',
            f'--collection-names ({len(args.collection_names)}) and '
            f'--collection-paths ({len(args.collection_paths)}) have different counts')
    path_map, stem_map = build_path_map(args.collection_names, args.collection_paths, FASTQ_EXTS)

    # ── Stage 3: parse and validate CSV columns ──
    with open(args.samples_csv, newline='') as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames or []

        col = {}
        missing = []
        for key, aliases in REQUIRED_COLS.items():
            found = resolve_col(fieldnames, aliases)
            if not found:
                missing.append(f'{key} (accepted: {", ".join(aliases)})')
            else:
                col[key] = found
        if missing:
            err('columns', 'Missing required column(s):\n  ' + '\n  '.join(missing)
                + f'\n  Found: {", ".join(fieldnames)}')

        rev_col = None
        for aliases in OPTIONAL_READ_COL.values():
            found = resolve_col(fieldnames, aliases)
            if found:
                rev_col = found
                break

        tax_missing = [c for c in TAX_COLS if c not in fieldnames]
        if tax_missing:
            err('columns', 'Missing taxonomy column(s) (used as expected taxonomy): '
                f'{", ".join(tax_missing)}\n  Found: {", ".join(fieldnames)}')
        if args.require_taxid and 'taxid' not in fieldnames:
            err('columns', 'Missing taxid column (required when gene-fetch taxonomy input type is "taxid")'
                f'\n  Found: {", ".join(fieldnames)}')
        rows = list(reader)

    # csv.DictReader stores surplus fields under the key None
    long_rows = [r[col['ID']] for r in rows if None in r]
    if long_rows:
        err('columns', f'Row(s) with more fields than the header: {", ".join(long_rows)}')

    # ── Stage 4: sample ID validation ──
    ids = [r[col['ID']] for r in rows]
    bad_ids = [i for i in ids if '_' in i]
    if bad_ids:
        err('sample IDs',
            'Sample IDs must not contain underscores. '
            f'Offending: {", ".join(bad_ids)}. Use hyphens instead.')
    dupes = [i for i in ids if ids.count(i) > 1]
    if dupes:
        err('sample IDs', f'Duplicate IDs: {", ".join(set(dupes))}')
    if args.require_taxid:
        no_taxid = [r[col['ID']] for r in rows if not (r.get('taxid') or '').strip()]
        if no_taxid:
            err('sample IDs', f'No taxid given for sample(s): {", ".join(no_taxid)}')

    # ── Stage 5: resolve filenames to symlinks of collection datasets ──
    os.makedirs(args.reads_dir, exist_ok=True)

    def resolve(filename, label):
        if not filename:
            return ''
        return resolve_file(filename, label, path_map, stem_map, FASTQ_EXTS, args.reads_dir, 'FASTQ')

    # Keep every uploaded column (e.g. taxid); only the read paths are rewritten
    out_rows = []
    for r in rows:
        out = dict(r)
        out[col['forward']] = resolve(r[col['forward']], f"sample {r[col['ID']]} forward")
        if rev_col:
            out[rev_col] = resolve(r.get(rev_col, ''), f"sample {r[col['ID']]} reverse")
        out_rows.append(out)

    # ── Stage 6: write output ──
    try:
        with open(args.output, 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=fieldnames)
            w.writeheader()
            w.writerows(out_rows)
    except OSError as e:
        err('output', f'Could not write {args.output}: {e}')

    print(f'[build_samples_csv] resolved {len(out_rows)} sample(s) to {args.output}')

    if args.refs_csv:
        resolve_references(args, ids)


if __name__ == '__main__':
    main()
