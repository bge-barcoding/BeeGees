"""Validate that config/config.yaml loads correctly and contains required keys/types.

This test loads the bundled config.yaml with PyYAML and checks structure —
it does not invoke Snakemake and requires no external tools.
"""
from pathlib import Path

import pytest
import yaml


CONFIG_PATH = Path(__file__).parent.parent / "config" / "config.yaml"


@pytest.fixture(scope="module")
def config():
    return yaml.safe_load(CONFIG_PATH.read_text())


class TestConfigYamlStructure:
    def test_file_exists(self):
        assert CONFIG_PATH.exists(), f"config.yaml not found at {CONFIG_PATH}"

    def test_parses_without_error(self, config):
        assert isinstance(config, dict)

    # Top-level required keys
    def test_has_run_name(self, config):
        assert "run_name" in config

    def test_has_samples_file(self, config):
        assert "samples_file" in config

    def test_has_output_dir(self, config):
        assert "output_dir" in config

    def test_has_r_list(self, config):
        assert "r" in config
        assert isinstance(config["r"], list)
        assert len(config["r"]) > 0

    def test_has_s_list(self, config):
        assert "s" in config
        assert isinstance(config["s"], list)
        assert len(config["s"]) > 0

    def test_has_run_gene_fetch_bool(self, config):
        assert "run_gene_fetch" in config
        assert isinstance(config["run_gene_fetch"], bool)

    # fastp section
    def test_fastp_section_present(self, config):
        assert "fastp" in config
        assert isinstance(config["fastp"], dict)

    def test_fastp_adapter_r1_is_string(self, config):
        val = config["fastp"].get("adapter_r1", "")
        assert isinstance(val, str)

    def test_fastp_adapter_r2_is_string(self, config):
        val = config["fastp"].get("adapter_r2", "")
        assert isinstance(val, str)

    # fasta_cleaner section
    def test_fasta_cleaner_section_present(self, config):
        assert "fasta_cleaner" in config
        assert isinstance(config["fasta_cleaner"], dict)

    def test_fasta_cleaner_consensus_threshold_in_range(self, config):
        val = config["fasta_cleaner"]["consensus_threshold"]
        assert 0.0 <= float(val) <= 1.0

    def test_fasta_cleaner_human_threshold_in_range(self, config):
        val = config["fasta_cleaner"]["human_threshold"]
        assert 0.0 <= float(val) <= 1.0

    def test_fasta_cleaner_at_difference_in_range(self, config):
        val = config["fasta_cleaner"]["at_difference"]
        assert 0.0 <= float(val) <= 1.0

    def test_fasta_cleaner_at_mode_valid(self, config):
        assert config["fasta_cleaner"]["at_mode"].lower() in {"absolute", "higher", "lower"}

    def test_fasta_cleaner_outlier_percentile_in_range(self, config):
        val = config["fasta_cleaner"]["outlier_percentile"]
        assert 0.0 <= float(val) <= 100.0

    # structural_validation section
    def test_run_structural_validation_is_bool(self, config):
        assert isinstance(config.get("run_structural_validation"), bool)

    def test_structural_validation_target_valid(self, config):
        target = config.get("structural_validation", {}).get("target", "")
        assert target.lower() in {"cox1", "coi", "rbcl"}

    # taxonomic_validation section
    def test_run_taxonomic_validation_is_bool(self, config):
        assert isinstance(config.get("run_taxonomic_validation"), bool)

    def test_taxonomic_validation_taxval_rank_valid(self, config):
        rank = config.get("taxonomic_validation", {}).get("taxval_rank", "")
        valid_ranks = {"phylum", "class", "order", "family", "genus", "species"}
        assert rank.lower() in valid_ranks

    # rules section
    def test_rules_section_present(self, config):
        assert "rules" in config
        assert isinstance(config["rules"], dict)

    def test_rules_have_mem_mb(self, config):
        for rule_name, rule_cfg in config["rules"].items():
            assert "mem_mb" in rule_cfg, f"Rule '{rule_name}' missing mem_mb"

    def test_rules_have_threads(self, config):
        for rule_name, rule_cfg in config["rules"].items():
            assert "threads" in rule_cfg, f"Rule '{rule_name}' missing threads"

    # MGE parameters
    def test_n_is_numeric(self, config):
        assert isinstance(config.get("n"), (int, float))

    def test_C_is_int(self, config):
        assert isinstance(config.get("C"), int)

    def test_t_is_float_in_range(self, config):
        val = config.get("t")
        assert 0.0 <= float(val) <= 1.0
