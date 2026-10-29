"""Score, identity and support fields.
"""
from __future__ import annotations
import GenAIRR as ga  # noqa: E402
from tests.support.fixtures.experiment_dsl import _human_igh_exp, _human_igh_records

import pytest

pytestmark = pytest.mark.behavior

def test_record_has_h6_score_fields():
    rec = _human_igh_exp().run_records(n=1, seed=0)[0]
    for seg in ("v", "d", "j"):
        for kind in ("score", "identity", "support"):
            assert f"{seg}_{kind}" in rec


def test_score_and_support_are_stubbed_none():
    rec = _human_igh_exp().run_records(n=1, seed=0)[0]
    # No aligner ran → these are aligner-specific concepts and are
    # always None for a simulator's output.
    for seg in ("v", "d", "j"):
        assert rec[f"{seg}_score"] is None
        assert rec[f"{seg}_support"] is None


def test_identity_is_high_for_pure_recombination():
    # Structural columns match the truth allele exactly, but the live
    # call's elastic boundaries may claim a few NP positions whose
    # random bases sometimes diverge from the projected allele's
    # continuation. Identity therefore stays high (>= 0.85) under pure
    # recombination without being strictly 1.0 in every record.
    for rec in _human_igh_exp().run_records(n=5, seed=0):
        assert rec["v_identity"] >= 0.85
        assert rec["d_identity"] >= 0.85
        assert rec["j_identity"] >= 0.85


def test_identity_drops_with_mutations():
    rec_clean = _human_igh_exp().run_records(n=1, seed=42)[0]
    rec_mutated = (
        _human_igh_exp().mutate(count=20).run_records(n=1, seed=42)[0]
    )
    # The mutated record's V identity should be below the clean
    # record's V identity (which is 1.0).
    assert rec_clean["v_identity"] == 1.0
    assert rec_mutated["v_identity"] < 1.0


def test_identity_in_unit_interval_under_corruption_stack():
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(count=15)
        .pcr_amplify(count=5)
        .polymerase_indels(count=3, insertion_prob=0.5)
    )
    for rec in exp.run_records(n=20, seed=0):
        for seg in ("v", "d", "j"):
            ident = rec[f"{seg}_identity"]
            assert ident is None or 0.0 <= ident <= 1.0


def test_d_identity_none_for_vj_chain():
    rec = (
        ga.Experiment.on("human_igk")
        .recombine()
        .run_records(n=1, seed=0)[0]
    )
    assert rec["d_identity"] is None
    # V and J identity still computed normally — elastic-boundary
    # extensions may claim a few NP positions whose random bases
    # sometimes diverge from the projected allele's continuation, so
    # identity stays high but isn't strictly 1.0.
    assert rec["v_identity"] >= 0.85
    assert rec["j_identity"] >= 0.85


def test_identity_low_for_contaminant_pass():
    # Full contaminant → seq is random bases against the original
    # V/D/J germline → identity collapses toward ~25% (random
    # 4-letter alphabet match rate).
    rec = (
        ga.Experiment.on("human_igh")
        .recombine()
        .contaminate(prob=1.0)
        .run_records(n=1, seed=0)[0]
    )
    assert rec["is_contaminant"] is True
    # Loose bound: we expect somewhere around 0.25, but allow a
    # generous range for sampling variance over a single record.
    assert 0.05 < rec["v_identity"] < 0.5


def test_identity_consistent_with_alignment_string_diff_count():
    # GROUND-TRUTH INVARIANT: the identity numerator must equal the
    # count of V/D/J columns in the alignment strings where sa
    # matches ga (case-insensitive, ignoring gaps and N).
    rec = (
        _human_igh_exp()
        .mutate(count=15)
        .polymerase_indels(count=3, insertion_prob=0.5)
        .run_records(n=1, seed=42)[0]
    )
    sa = rec["sequence_alignment"]
    galn = rec["germline_alignment"]
    for seg in ("v", "d", "j"):
        a_s = rec[f"{seg}_alignment_start"]
        a_e = rec[f"{seg}_alignment_end"]
        if a_s is None:
            continue
        total = a_e - a_s
        matches = sum(
            1
            for sa_c, ga_c in zip(sa[a_s:a_e], galn[a_s:a_e])
            if sa_c not in ("-",)
            and ga_c not in ("-", "N")
            and sa_c.upper() == ga_c.upper()
        )
        expected_identity = matches / total if total > 0 else None
        assert (
            rec[f"{seg}_identity"] == expected_identity
        ), f"{seg}: identity {rec[f'{seg}_identity']} != recomputed {expected_identity}"


def test_identity_decreases_monotonically_with_mutation_count():
    # As mutation count grows, average identity over V should drop.
    means = {}
    for count in (0, 5, 25):
        records = (
            _human_igh_exp().mutate(count=count).run_records(n=20, seed=0)
        )
        means[count] = sum(r["v_identity"] for r in records) / len(records)
    assert means[0] > means[5] > means[25]
    # No mutations → identity stays very close to 1.0; tiny drops from
    # the live call's elastic boundaries are tolerable.
    assert means[0] >= 0.99


def test_h6_fields_in_default_column_order():
    from GenAIRR.result import _DEFAULT_COLUMN_ORDER

    for seg in ("v", "d", "j"):
        for kind in ("score", "identity", "support"):
            assert f"{seg}_{kind}" in _DEFAULT_COLUMN_ORDER


def test_h6_fields_in_tsv_export(tmp_path):
    import csv as _csv

    result = _human_igh_records(n=2)
    path = tmp_path / "h6.tsv"
    result.to_tsv(str(path))
    with open(path, encoding="utf-8") as fh:
        rows = list(_csv.DictReader(fh, delimiter="\t"))
    for row in rows:
        # Stubbed fields export as empty strings.
        assert row["v_score"] == ""
        assert row["v_support"] == ""
        # Identity is a float string. Elastic-boundary extensions can
        # nick a percent off in pure-recombination records, so accept
        # any value at or above 0.85.
        assert float(row["v_identity"]) >= 0.85
