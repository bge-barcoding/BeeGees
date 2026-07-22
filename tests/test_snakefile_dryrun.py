"""Smoke test: Snakemake dryrun with a minimal config and samples.csv.

Requires snakemake and MitoGeneExtractor to be on PATH (i.e. the BeeGees
conda environment must be active).  The test is automatically skipped when
either is absent, so it is safe to run in bare Python environments.

What this test checks:
- The Snakefile parses without syntax errors.
- Snakemake can resolve the full DAG from a valid config (--dryrun).

It does NOT submit or execute any jobs.
"""
import shutil
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

# Skip the whole module when snakemake is not available.
pytest.importorskip("snakemake", reason="snakemake not installed; skipping dryrun tests")

SNAKEFILE = Path(__file__).parent.parent / "workflow" / "Snakefile"


@pytest.fixture()
def minimal_run(tmp_path):
    """Create the minimum files needed for a Snakemake dryrun."""
    # samples.csv with one fake sample
    samples = tmp_path / "samples.csv"
    samples.write_text("ID,forward,reverse\nSAMPLE01,/data/r1.fastq.gz,/data/r2.fastq.gz\n")

    # Minimal config pointing at the fake samples file
    config_text = textwrap.dedent(f"""\
        run_name: "DRYRUN_TEST"
        samples_file: "{samples}"
        output_dir: "{tmp_path / 'output'}"
        run_gene_fetch: false
        sequence_reference_file: "{samples}"

        fastp:
          adapter_r1: ""
          adapter_r2: ""
          extra_fastp_args: ""

        downsampling:
          enabled: false
          max_reads: 0

        r:
          - 1
          - 1.3
          - 1.5
        s:
          - 50
          - 100
        n: 0
        C: 5
        t: 0.5

        fasta_cleaner:
          consensus_threshold: 0.5
          human_threshold: 0.95
          at_difference: 0.1
          at_mode: "absolute"
          outlier_percentile: 90.0
          reference_dir: null
          reference_filter_mode: "remove_similar"

        run_structural_validation: false
        structural_validation:
          target: "cox1"
          verbose: false
          genetic_code: 5

        run_taxonomic_validation: false
        taxonomic_validation:
          database: ""
          database_taxonomy: ""
          taxval_rank: "family"
          expected_taxonomy: ""
          verbose: false
          min_pident: 80
          min_length: 100

        rules:
          gene_fetch:
            mem_mb: 4096
            threads: 1
            partition: short
          fastp_qc:
            mem_mb: 4096
            threads: 1
          clean_headers_merge:
            mem_mb: 4096
            threads: 1
          fastq_concat:
            mem_mb: 4096
            threads: 1
          quality_trim:
            mem_mb: 4096
            threads: 1
          downsample:
            mem_mb: 4096
            threads: 1
          MitoGeneExtractor:
            mem_mb: 4096
            threads: 1
            partition: short
          rename_and_combine_cons:
            mem_mb: 4096
            threads: 1
          gzip_merged_clean:
            mem_mb: 4096
            threads: 1
          human_cox1_filter:
            mem_mb: 4096
            threads: 1
          at_content_filter:
            mem_mb: 4096
            threads: 1
          statistical_outlier_filter:
            mem_mb: 4096
            threads: 1
          reference_filter:
            mem_mb: 4096
            threads: 1
          consensus_generation:
            mem_mb: 4096
            threads: 1
          extract_stats_to_csv:
            mem_mb: 4096
            threads: 1
          structural_validation:
            mem_mb: 4096
            threads: 1
            partition: short
          taxonomic_validation:
            mem_mb: 4096
            threads: 4
            partition: short
          blast2taxonomy:
            mem_mb: 4096
            threads: 1
          download_taxdump:
            mem_mb: 4096
            threads: 1
          multiqc_plots:
            mem_mb: 4096
            threads: 1
          multiqc:
            mem_mb: 4096
            threads: 1
    """)
    config_file = tmp_path / "config.yaml"
    config_file.write_text(config_text)
    return tmp_path, config_file


class TestSnakefileDryrun:
    def test_snakefile_parses_and_dryrun_succeeds(self, minimal_run):
        """snakemake --dryrun should exit 0 with a valid config."""
        tmp_path, config_file = minimal_run
        result = subprocess.run(
            [
                sys.executable, "-m", "snakemake",
                "--snakefile", str(SNAKEFILE),
                "--configfile", str(config_file),
                "--dryrun",
                "--quiet",
                "--cores", "1",
            ],
            capture_output=True,
            text=True,
            cwd=str(tmp_path),
        )
        assert result.returncode == 0, (
            f"Snakemake dryrun failed.\n"
            f"STDOUT:\n{result.stdout}\n"
            f"STDERR:\n{result.stderr}"
        )
