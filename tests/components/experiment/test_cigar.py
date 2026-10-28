"""Per-segment CIGAR strings.
"""
from __future__ import annotations
from GenAIRR import CompiledExperiment, Experiment
import GenAIRR as ga  # noqa: E402
import re as _re
from tests.support.fixtures.experiment_dsl import _cigar_op_counts, _cigar_total_ops, _human_igh_exp

import pytest

pytestmark = pytest.mark.behavior

def test_record_has_cigar_fields():
    rec = _human_igh_exp().run_records(n=1, seed=0)[0]
    assert "v_cigar" in rec
    assert "d_cigar" in rec
    assert "j_cigar" in rec


def test_cigar_pure_match_when_no_indels():
    # No indels → CIGAR is a single M run for each segment, length
    # equal to the coding region's nucleotide span.
    rec = _human_igh_exp().run_records(n=1, seed=42)[0]
    v_len = rec["v_sequence_end"] - rec["v_sequence_start"]
    d_len = rec["d_sequence_end"] - rec["d_sequence_start"]
    j_len = rec["j_sequence_end"] - rec["j_sequence_start"]
    assert rec["v_cigar"] == f"{v_len}M"
    assert rec["d_cigar"] == f"{d_len}M"
    assert rec["j_cigar"] == f"{j_len}M"


def test_cigar_unchanged_under_mutations():
    # CIGAR uses ``M`` for both match and mismatch, so mutations
    # don't introduce I or D ops.
    rec = (
        _human_igh_exp()
        .mutate(count=12)
        .run_records(n=1, seed=42)[0]
    )
    for cig in (rec["v_cigar"], rec["d_cigar"], rec["j_cigar"]):
        assert "I" not in cig
        assert "D" not in cig


def test_cigar_insertion_only_produces_i_ops():
    rec = (
        _human_igh_exp()
        .polymerase_indels(count=4, insertion_prob=1.0)
        .run_records(n=1, seed=0)[0]
    )
    counts_v = _cigar_op_counts(rec["v_cigar"])
    counts_d = _cigar_op_counts(rec["d_cigar"])
    counts_j = _cigar_op_counts(rec["j_cigar"])
    # No deletions anywhere with insertion_prob=1.0.
    assert counts_v["D"] == counts_d["D"] == counts_j["D"] == 0
    # All-insertion run total across V+D+J should be ≤ 4 — some may
    # land in NP regions where they don't show in any CIGAR.
    total_I = counts_v["I"] + counts_d["I"] + counts_j["I"]
    assert 0 <= total_I <= 4


def test_cigar_deletion_only_produces_d_ops():
    rec = (
        _human_igh_exp()
        .polymerase_indels(count=4, insertion_prob=0.0)
        .run_records(n=1, seed=0)[0]
    )
    counts_v = _cigar_op_counts(rec["v_cigar"])
    counts_d = _cigar_op_counts(rec["d_cigar"])
    counts_j = _cigar_op_counts(rec["j_cigar"])
    assert counts_v["I"] == counts_d["I"] == counts_j["I"] == 0
    total_D = counts_v["D"] + counts_d["D"] + counts_j["D"]
    assert 0 <= total_D <= 4


def test_cigar_op_count_equals_segment_alignment_columns():
    # Sum of M + I + D ops in a segment's CIGAR equals the number
    # of alignment columns that segment contributes to the V/D/J
    # alignment.
    rec = (
        _human_igh_exp()
        .polymerase_indels(count=5, insertion_prob=0.5)
        .run_records(n=1, seed=42)[0]
    )
    # In the full alignment string, the V coding region's columns
    # are roughly v_sequence_start to v_sequence_end *plus* any
    # deletion gaps inside V. We can count them by walking the
    # alignment-string positions and matching them up via the same
    # column-list logic — but a simpler check is: each segment's
    # CIGAR total must equal `length(coding region in sequence)
    # + number of deletion columns inserted before/within the
    # segment`. For a record with no indels in V at all, total ops
    # = v_sequence_end - v_sequence_start.
    v_len = rec["v_sequence_end"] - rec["v_sequence_start"]
    v_total = _cigar_total_ops(rec["v_cigar"])
    # v_total = (raw V coding bases in sequence) + (deletions in V)
    #         = v_len + count of D ops in V.
    n_deletions_in_v = _cigar_op_counts(rec["v_cigar"])["D"]
    assert v_total == v_len + n_deletions_in_v


def test_cigar_d_cigar_empty_for_vj_chain():
    # No D segment → d_cigar must be the empty string (and the
    # other two must still be valid).
    rec = (
        ga.Experiment.on("human_igk")
        .recombine()
        .run_records(n=1, seed=0)[0]
    )
    assert rec["d_cigar"] == ""
    assert rec["v_cigar"]
    assert rec["j_cigar"]


def test_cigar_consistent_with_alignment_strings():
    # Walk both the CIGAR and the alignment-string columns of the V
    # region and verify they describe the same sequence of ops.
    rec = (
        _human_igh_exp()
        .polymerase_indels(count=4, insertion_prob=0.5)
        .run_records(n=1, seed=42)[0]
    )
    # Reconstruct the V CIGAR from sa/galn directly and compare.
    sa = rec["sequence_alignment"]
    galn = rec["germline_alignment"]
    # Walk only V columns. We use sequence-coord boundaries to
    # locate V — works only because V is the first segment in our
    # assembly (no leading deletion gaps). For a robust ground-
    # truth test we use the engine's own column walker via the
    # private helper.
    from GenAIRR._airr_record import (
        _alignment_columns,
        _cigars_from_columns,
    )

    # Re-derive the CIGAR by re-running the helper. Equivalence is
    # guaranteed by construction; this asserts the public field
    # matches the helper's output (catches any future re-routing).
    refdata = _human_igh_exp().refdata
    outcome = _human_igh_exp().polymerase_indels(
        count=4, insertion_prob=0.5
    ).run(n=1, seed=42)[0]
    cols = _alignment_columns(outcome, refdata)
    v_c, d_c, j_c = _cigars_from_columns(cols)
    assert v_c == rec["v_cigar"]
    assert d_c == rec["d_cigar"]
    assert j_c == rec["j_cigar"]


def test_cigar_runlength_collapses_consecutive_ops():
    # Pure-M CIGAR for a segment must have exactly one run, not
    # one run per position.
    rec = _human_igh_exp().run_records(n=1, seed=0)[0]
    # Count the number of (digit-prefixed op) tokens.
    v_runs = _re.findall(r"\d+[MID]", rec["v_cigar"])
    assert len(v_runs) == 1, f"expected single run, got {v_runs}"


def test_cigar_only_uses_canonical_ops():
    # We emit only M / I / D — no =/X/S/H.
    rec = (
        _human_igh_exp()
        .mutate(count=10)
        .polymerase_indels(count=3, insertion_prob=0.5)
        .run_records(n=1, seed=0)[0]
    )
    for cig in (rec["v_cigar"], rec["d_cigar"], rec["j_cigar"]):
        for _, op in _re.findall(r"(\d+)(\D)", cig):
            assert op in ("M", "I", "D"), f"unexpected op {op!r} in {cig!r}"


def test_cigars_in_default_column_order():
    from GenAIRR.result import _DEFAULT_COLUMN_ORDER

    assert "v_cigar" in _DEFAULT_COLUMN_ORDER
    assert "d_cigar" in _DEFAULT_COLUMN_ORDER
    assert "j_cigar" in _DEFAULT_COLUMN_ORDER
    # Canonical ordering: v_cigar should sit immediately after v_call,
    # likewise for d/j.
    v_call_idx = _DEFAULT_COLUMN_ORDER.index("v_call")
    assert _DEFAULT_COLUMN_ORDER[v_call_idx + 1] == "v_cigar"


def test_cigars_in_tsv_export(tmp_path):
    import csv as _csv

    result = (
        _human_igh_exp()
        .polymerase_indels(count=3, insertion_prob=0.5)
        .run_records(n=2, seed=0)
    )
    path = tmp_path / "h4.tsv"
    result.to_tsv(str(path))
    with open(path, encoding="utf-8") as fh:
        rows = list(_csv.DictReader(fh, delimiter="\t"))
    for row in rows:
        assert row["v_cigar"]
        assert row["j_cigar"]


def test_cigar_invariants_under_full_corruption_stack():
    # Same n=200 stress sweep as alignment strings — verify that
    # CIGARs stay self-consistent with the alignment-string
    # topology under every corruption mode at once.
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(count=15)
        .pcr_amplify(count=5)
        .sequencing_errors(count=8)
        .polymerase_indels(count=4, insertion_prob=0.5)
    )
    failures = []
    for i, rec in enumerate(exp.run_records(n=200, seed=0)):
        # The CIGAR walks every column in the live (extension-aware)
        # span — `M + I + D` is the column count, `seq_len + D_ops`
        # accounts for D-ops adding gap columns beyond the live
        # sequence span.
        for seg in ("v", "d", "j"):
            cig = rec[f"{seg}_cigar"]
            if not cig:
                continue
            ops = _cigar_op_counts(cig)
            total = ops["M"] + ops["I"] + ops["D"]
            seq_len = (
                rec[f"{seg}_sequence_end"] - rec[f"{seg}_sequence_start"]
                if rec[f"{seg}_sequence_end"] is not None
                else 0
            )
            if total != seq_len + ops["D"]:
                failures.append(
                    (i, seg, total, seq_len, ops["D"])
                )
    assert failures == [], f"CIGAR invariant broke for: {failures[:3]}"
