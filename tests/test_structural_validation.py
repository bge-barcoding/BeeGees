"""Tests for structural_validation.py pure functions (no nhmmer required)."""
import importlib.util
from pathlib import Path

import pytest
from beegees.utils.configs import get_package_dir

_script = get_package_dir() / "workflow" / "scripts" / "structural_validation.py"
_spec = importlib.util.spec_from_file_location("struct_val", _script)
struct_val = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(struct_val)


class TestTrimNCharacters:
    def test_no_n_unchanged(self):
        assert struct_val.trim_n_characters("ATCG") == "ATCG"

    def test_leading_n_trimmed(self):
        assert struct_val.trim_n_characters("NNATCG") == "ATCG"

    def test_trailing_n_trimmed(self):
        assert struct_val.trim_n_characters("ATCGNN") == "ATCG"

    def test_both_ends_trimmed(self):
        assert struct_val.trim_n_characters("NNATCGNN") == "ATCG"

    def test_internal_n_preserved(self):
        result = struct_val.trim_n_characters("ATNCG")
        assert "N" in result

    def test_all_n_returns_empty(self):
        result = struct_val.trim_n_characters("NNNN")
        assert result == ""

    def test_empty_string(self):
        assert struct_val.trim_n_characters("") == ""


class TestCalculateBarcodeBaseCount:
    def test_counts_non_gap_non_n(self):
        # ATCG → 4 bases
        assert struct_val.calculate_barcode_base_count("ATCG") == 4

    def test_gaps_not_excluded(self):
        # gaps count toward length; only N characters are subtracted
        assert struct_val.calculate_barcode_base_count("AT-CG") == 5

    def test_n_excluded(self):
        assert struct_val.calculate_barcode_base_count("ATNCG") == 4

    def test_empty_sequence(self):
        assert struct_val.calculate_barcode_base_count("") == 0

    def test_all_gaps(self):
        # no N characters, so length is returned unchanged
        assert struct_val.calculate_barcode_base_count("----") == 4


class TestCalculateBarcodeRank:
    def test_rank_1_long_perfect(self):
        rank = struct_val.calculate_barcode_rank(
            barcode_ambiguous_bases_original=0,
            stop_codons=0,
            reading_frame_valid=True,
            barcode_base_count=550,
        )
        assert rank == 1

    def test_rank_2_400_to_499(self):
        rank = struct_val.calculate_barcode_rank(0, 0, True, 450)
        assert rank == 2

    def test_rank_3_300_to_399(self):
        rank = struct_val.calculate_barcode_rank(0, 0, True, 350)
        assert rank == 3

    def test_rank_4_200_to_299(self):
        rank = struct_val.calculate_barcode_rank(0, 0, True, 250)
        assert rank == 4

    def test_rank_5_1_to_199(self):
        rank = struct_val.calculate_barcode_rank(0, 0, True, 100)
        assert rank == 5

    def test_rank_6_original_n_present(self):
        rank = struct_val.calculate_barcode_rank(1, 0, True, 600)
        assert rank == 6

    def test_rank_6_stop_codon(self):
        rank = struct_val.calculate_barcode_rank(0, 1, True, 600)
        assert rank == 6

    def test_rank_6_invalid_reading_frame(self):
        rank = struct_val.calculate_barcode_rank(0, 0, False, 600)
        assert rank == 6


class TestPassesQualityCriteria:
    def _result(self, **kwargs):
        base = {
            "barcode_ambiguous_bases_original": 0,
            "stop_codons": 0,
            "reading_frame": 0,
            "barcode_base_count": 500,
            "barcode_ambiguous_bases": 0,
        }
        base.update(kwargs)
        return base

    def test_good_sequence_passes(self):
        assert struct_val.passes_quality_criteria(self._result()) is True

    def test_original_n_fails(self):
        assert struct_val.passes_quality_criteria(self._result(barcode_ambiguous_bases_original=1)) is False

    def test_stop_codon_fails(self):
        assert struct_val.passes_quality_criteria(self._result(stop_codons=1)) is False

    def test_invalid_reading_frame_fails(self):
        assert struct_val.passes_quality_criteria(self._result(reading_frame=-1)) is False

    def test_too_short_fails(self):
        assert struct_val.passes_quality_criteria(self._result(barcode_base_count=200)) is False

    def test_high_ambiguity_fails(self):
        # 30% ambiguity exactly at threshold → fails (>= 0.30)
        result = self._result(barcode_base_count=100, barcode_ambiguous_bases=30)
        assert struct_val.passes_quality_criteria(result) is False

    def test_just_below_ambiguity_threshold_passes(self):
        # barcode_base_count must be > 300; use 400 to satisfy that check
        result = self._result(barcode_base_count=400, barcode_ambiguous_bases=29)
        assert struct_val.passes_quality_criteria(result) is True


class TestRunNhmmerOnSequence:
    """The nhmmer wrapper: --cpu wiring and temp-file cleanup.

    nhmmer is not available in CI, so subprocess.run is replaced with a fake that
    writes a tabular output file and records the command it was given.
    """

    @staticmethod
    def _fake_nhmmer(recorder, write_tblout=True):
        import subprocess as _subprocess

        def fake_run(cmd, *args, **kwargs):
            recorder.append(list(cmd))
            if write_tblout:
                tblout = cmd[cmd.index('--tblout') + 1]
                # A no-hit tabular output: comment lines only.
                Path(tblout).write_text("# target name  accession  query name\n#\n")
            return _subprocess.CompletedProcess(cmd, 0, "", "")

        return fake_run

    def test_cpu_defaults_to_one(self, monkeypatch):
        calls = []
        monkeypatch.setattr(struct_val.subprocess, "run", self._fake_nhmmer(calls))

        struct_val.run_nhmmer_on_sequence("ACGT" * 20, "SEQ1", "/fake/COI-5P.hmm")

        cmd = calls[0]
        assert cmd[cmd.index('--cpu') + 1] == "1"

    def test_cpu_reflects_threads_argument(self, monkeypatch):
        calls = []
        monkeypatch.setattr(struct_val.subprocess, "run", self._fake_nhmmer(calls))

        struct_val.run_nhmmer_on_sequence("ACGT" * 20, "SEQ1", "/fake/COI-5P.hmm", threads=8)

        cmd = calls[0]
        assert cmd[cmd.index('--cpu') + 1] == "8"

    def test_no_temp_files_left_behind(self, monkeypatch, tmp_path):
        """Regression: the scratch query/tabular files used to leak, two per sequence."""
        tmpdir = tmp_path / "scratch"
        tmpdir.mkdir()
        monkeypatch.setenv("TMPDIR", str(tmpdir))
        monkeypatch.setattr(struct_val.tempfile, "tempdir", None)
        calls = []
        monkeypatch.setattr(struct_val.subprocess, "run", self._fake_nhmmer(calls))

        for i in range(25):
            struct_val.run_nhmmer_on_sequence("ACGT" * 20, f"SEQ{i}", "/fake/COI-5P.hmm")

        assert len(calls) == 25
        # Guard against a vacuous pass: the scratch files must actually have been
        # created under the directory we are then asserting is empty.
        used = [c[c.index('--tblout') + 1] for c in calls]
        assert all(u.startswith(str(tmpdir)) for u in used), used[:2]
        assert list(tmpdir.iterdir()) == []

    def test_temp_files_cleaned_up_when_nhmmer_fails(self, monkeypatch, tmp_path):
        tmpdir = tmp_path / "scratch"
        tmpdir.mkdir()
        monkeypatch.setenv("TMPDIR", str(tmpdir))
        monkeypatch.setattr(struct_val.tempfile, "tempdir", None)

        import subprocess as _subprocess

        seen = []

        def failing_run(cmd, *args, **kwargs):
            seen.append(cmd[cmd.index('--tblout') + 1])
            return _subprocess.CompletedProcess(cmd, 1, "", "nhmmer: fatal error")

        monkeypatch.setattr(struct_val.subprocess, "run", failing_run)

        assert struct_val.run_nhmmer_on_sequence("ACGT" * 20, "SEQ1", "/fake/hmm") is None
        assert seen and seen[0].startswith(str(tmpdir)), seen
        assert list(tmpdir.iterdir()) == []

    def test_missing_tabular_output_is_not_reported_as_missing_nhmmer(self, monkeypatch, caplog):
        """An absent --tblout must read as 'no hits', not 'nhmmer not found in PATH'."""
        calls = []
        monkeypatch.setattr(struct_val.subprocess, "run",
                            self._fake_nhmmer(calls, write_tblout=False))

        result = struct_val.run_nhmmer_on_sequence("ACGT" * 20, "SEQ1", "/fake/hmm")

        assert result is None
        assert "nhmmer not found in PATH" not in caplog.text


# nhmmer --tblout column order:
#   target acc query acc hmmfrom hmmto alifrom alito envfrom envto sqlen strand E score bias desc
def _tblout_row(seq_id, hmm_from, hmm_to, seq_from, seq_to, evalue, score=100.0, bias=10.0):
    return (f"{seq_id} - COI-5P - {hmm_from} {hmm_to} {seq_from} {seq_to} "
            f"{seq_from} {seq_to} 1587 + {evalue} {score} {bias} -")


def _tblout(seq_id, *rows):
    return "\n".join(["# target name  accession  query name", "#"] + list(rows))


class TestParseNhmmerResultEnvelopes:
    """nhmmer reports one envelope per contiguous match. Keeping only the best
    one truncated barcodes to whichever fragment scored highest - on UK016-H10 a
    299-column envelope at E=6.1e-63 was silently dropped in favour of a
    291-column one at E=3.9e-73, yielding 290 bases from a consensus with a
    900-base unambiguous stretch."""

    SEQ = "UK016-H10_r_1.5_s_50_fcleaner_concat"

    def test_real_h10_case_keeps_both_envelopes(self):
        content = _tblout(
            self.SEQ,
            _tblout_row(self.SEQ, 367, 657, 462, 752, "3.9e-73", 235.2, 38.6),
            _tblout_row(self.SEQ, 2, 300, 97, 395, "6.1e-63", 201.6, 29.2),
            _tblout_row(self.SEQ, 444, 461, 52, 69, "3.5", -7.0, 6.7),
        )
        result = struct_val.parse_nhmmer_result(content, self.SEQ)

        # Sorted by hmm_from, and the E=3.5 row is excluded as insignificant.
        assert [(e["hmm_from"], e["hmm_to"]) for e in result] == [(2, 300), (367, 657)]

    def test_single_envelope_returned_as_list(self):
        content = _tblout(self.SEQ, _tblout_row(self.SEQ, 119, 657, 214, 752, "8.8e-124"))
        result = struct_val.parse_nhmmer_result(content, self.SEQ)

        assert len(result) == 1
        assert (result[0]["hmm_from"], result[0]["hmm_to"]) == (119, 657)

    def test_no_significant_envelope_returns_none(self):
        """None, not [] - run_nhmmer_on_sequence hands this to callers testing `is None`."""
        content = _tblout(self.SEQ, _tblout_row(self.SEQ, 444, 461, 52, 69, "3.5"))
        assert struct_val.parse_nhmmer_result(content, self.SEQ) is None

    def test_empty_content_returns_none(self):
        assert struct_val.parse_nhmmer_result("", self.SEQ) is None

    def test_name_mismatch_excluded(self):
        content = _tblout(self.SEQ, _tblout_row("SOME_OTHER_SEQ", 2, 300, 97, 395, "1e-50"))
        assert struct_val.parse_nhmmer_result(content, self.SEQ) is None

    def test_hmm_overlap_keeps_only_best(self):
        content = _tblout(
            self.SEQ,
            _tblout_row(self.SEQ, 100, 400, 100, 400, "1e-90"),
            _tblout_row(self.SEQ, 350, 657, 800, 1107, "1e-50"),  # HMM overlaps 350-400
        )
        result = struct_val.parse_nhmmer_result(content, self.SEQ)

        assert [(e["hmm_from"], e["hmm_to"]) for e in result] == [(100, 400)]

    def test_sequence_overlap_keeps_only_best(self):
        """Disjoint in HMM space but sharing query bases: merging would place the
        same bases at two model positions."""
        content = _tblout(
            self.SEQ,
            _tblout_row(self.SEQ, 2, 300, 100, 398, "1e-90"),
            _tblout_row(self.SEQ, 367, 657, 350, 640, "1e-50"),  # Seq overlaps 350-398
        )
        result = struct_val.parse_nhmmer_result(content, self.SEQ)

        assert [(e["hmm_from"], e["hmm_to"]) for e in result] == [(2, 300)]

    def test_best_envelope_always_survives(self):
        """Greedy best-first: the envelope that won under the old best-only logic
        is kept whatever else is present, so a barcode can never get shorter."""
        content = _tblout(
            self.SEQ,
            _tblout_row(self.SEQ, 300, 400, 300, 400, "1e-99"),
            _tblout_row(self.SEQ, 250, 350, 250, 350, "1e-98"),
            _tblout_row(self.SEQ, 350, 450, 350, 450, "1e-97"),
        )
        result = struct_val.parse_nhmmer_result(content, self.SEQ)

        assert [(e["hmm_from"], e["hmm_to"]) for e in result] == [(300, 400)]

    def test_three_disjoint_envelopes_all_kept_and_sorted(self):
        content = _tblout(
            self.SEQ,
            _tblout_row(self.SEQ, 400, 500, 400, 500, "1e-70"),
            _tblout_row(self.SEQ, 2, 100, 2, 100, "1e-90"),
            _tblout_row(self.SEQ, 200, 300, 200, 300, "1e-80"),
        )
        result = struct_val.parse_nhmmer_result(content, self.SEQ)

        assert [(e["hmm_from"], e["hmm_to"]) for e in result] == [(2, 100), (200, 300), (400, 500)]

    def test_discarded_overlap_is_logged(self, caplog):
        """The old code dropped significant envelopes with no log line at all."""
        content = _tblout(
            self.SEQ,
            _tblout_row(self.SEQ, 100, 400, 100, 400, "1e-90"),
            _tblout_row(self.SEQ, 350, 657, 800, 1107, "1e-50"),
        )
        with caplog.at_level("INFO"):
            struct_val.parse_nhmmer_result(content, self.SEQ)

        assert "Discarded envelope" in caplog.text
        assert "overlaps accepted envelope" in caplog.text


class TestEnvelopesOverlap:
    @staticmethod
    def _env(hmm_from, hmm_to, seq_from, seq_to):
        return {"hmm_from": hmm_from, "hmm_to": hmm_to,
                "seq_from": seq_from, "seq_to": seq_to}

    def test_disjoint_in_both(self):
        overlaps, _ = struct_val.envelopes_overlap(
            self._env(2, 300, 97, 395), self._env(367, 657, 462, 752))
        assert overlaps is False

    def test_hmm_only(self):
        overlaps, where = struct_val.envelopes_overlap(
            self._env(2, 400, 97, 495), self._env(367, 657, 900, 1190))
        assert overlaps is True
        assert where == "HMM"

    def test_sequence_only(self):
        overlaps, where = struct_val.envelopes_overlap(
            self._env(2, 300, 97, 500), self._env(367, 657, 400, 690))
        assert overlaps is True
        assert where == "sequence"

    def test_adjacent_ranges_do_not_overlap(self):
        """Inclusive coordinates: 1-300 and 301-600 touch but do not share a position."""
        overlaps, _ = struct_val.envelopes_overlap(
            self._env(1, 300, 1, 300), self._env(301, 600, 301, 600))
        assert overlaps is False

    def test_single_shared_position_overlaps(self):
        overlaps, _ = struct_val.envelopes_overlap(
            self._env(1, 300, 1, 300), self._env(300, 600, 900, 1200))
        assert overlaps is True

    def test_minus_strand_coordinates_normalised(self):
        """nhmmer reports minus-strand hits with seq_from > seq_to."""
        overlaps, where = struct_val.envelopes_overlap(
            self._env(2, 300, 500, 100), self._env(367, 657, 400, 690))
        assert overlaps is True
        assert where == "sequence"


class TestConstructHmmSpaceMultipleEnvelopes:
    @staticmethod
    def _env(hmm_from, hmm_to, seq_from, seq_to):
        return {"target_name": "SEQ1", "hmm_from": hmm_from, "hmm_to": hmm_to,
                "seq_from": seq_from, "seq_to": seq_to}

    def test_two_envelopes_both_placed_with_gap_between(self):
        seq = "A" * 10 + "C" * 10          # positions 1-10 A, 11-20 C
        result = struct_val.construct_hmm_space_from_alignment(
            [self._env(1, 10, 1, 10), self._env(21, 30, 11, 20)], seq, 30)

        assert result == "A" * 10 + "-" * 10 + "C" * 10

    def test_single_dict_still_accepted(self):
        result = struct_val.construct_hmm_space_from_alignment(
            self._env(1, 5, 1, 5), "ACGTA", 10)
        assert result == "ACGTA" + "-" * 5

    def test_envelope_order_does_not_matter(self):
        seq = "A" * 10 + "C" * 10
        forward = struct_val.construct_hmm_space_from_alignment(
            [self._env(1, 10, 1, 10), self._env(21, 30, 11, 20)], seq, 30)
        reversed_order = struct_val.construct_hmm_space_from_alignment(
            [self._env(21, 30, 11, 20), self._env(1, 10, 1, 10)], seq, 30)
        assert forward == reversed_order

    def test_h10_barcode_spans_both_envelopes_end_to_end(self):
        """Regression for the 290-base truncation: HMM 2-300 + 367-657 must give
        a 656-position barcode with 590 real bases, not 290."""
        envelopes = [self._env(2, 300, 97, 395), self._env(367, 657, 462, 752)]
        hmm_space = struct_val.construct_hmm_space_from_alignment(envelopes, "A" * 1587, 657)
        barcode = struct_val.trim_sequence_ends(struct_val.replace_gaps_with_n(hmm_space))

        assert len(barcode) == 656
        assert barcode.count("N") == 66                      # HMM 301-366 unfilled
        assert struct_val.calculate_barcode_base_count(barcode) == 590
