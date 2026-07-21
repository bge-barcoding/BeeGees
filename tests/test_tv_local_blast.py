"""Tests for BLASTRunner.sanitize_header (static method) in tv_local_blast.py."""
import importlib.util
from pathlib import Path

_script = Path(__file__).parent.parent / "workflow" / "scripts" / "tv_local_blast.py"
_spec = importlib.util.spec_from_file_location("tv_local_blast", _script)
tv_local_blast = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(tv_local_blast)

sanitize = tv_local_blast.BLASTRunner.sanitize_header


class TestSanitizeHeader:
    def test_strips_gt(self):
        assert not sanitize(">SAMPLE123").startswith(">")

    def test_replaces_pipe(self):
        assert "|" not in sanitize(">SAMPLE|extra")

    def test_replaces_space(self):
        assert " " not in sanitize(">SAMPLE with spaces")

    def test_replaces_colon(self):
        assert ":" not in sanitize(">SAMPLE:colon")

    def test_replaces_slash(self):
        assert "/" not in sanitize(">path/to/seq")

    def test_plain_id_unchanged_except_gt(self):
        assert sanitize(">BSNHM089-24") == "BSNHM089-24"

    def test_long_header_truncated_at_100(self):
        long_header = ">" + "A" * 200
        assert len(sanitize(long_header)) <= 100

    def test_empty_string(self):
        assert sanitize("") == ""

    def test_no_gt_prefix_unchanged(self):
        assert sanitize("BSNHM089-24") == "BSNHM089-24"
