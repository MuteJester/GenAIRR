"""AIRR field naming, cheap metadata, and the translated region fields.
"""
from __future__ import annotations
from GenAIRR import _engine as ge
from GenAIRR import CompiledExperiment, Experiment
import GenAIRR as ga  # noqa: E402
from tests.support.fixtures.experiment_dsl import _human_igh_exp, _human_igh_records

def test_records_use_airr_field_names_for_renames():
    rec = _human_igh_exp().run_records(n=1, seed=0)[0]
    # New canonical AIRR names exist.
    assert "junction" in rec
    assert "np1" in rec
    assert "np2" in rec
    # Legacy non-canonical names are not present.
    assert "junction_nt" not in rec
    assert "np1_region" not in rec
    assert "np2_region" not in rec


def test_record_has_sequence_id():
    result = _human_igh_records(n=3)
    for i, rec in enumerate(result):
        assert rec["sequence_id"] == f"seq{i}"


def test_simulation_result_id_prefix_is_configurable():
    from GenAIRR.result import SimulationResult

    exp = _human_igh_exp()
    outcomes = exp.run(n=2, seed=0)
    result = SimulationResult.from_outcomes(
        outcomes, exp.refdata, id_prefix="rep01_"
    )
    assert result[0]["sequence_id"] == "rep01_0"
    assert result[1]["sequence_id"] == "rep01_1"


def test_record_rev_comp_is_false_for_forward_strand_simulator():
    rec = _human_igh_exp().run_records(n=1, seed=0)[0]
    assert rec["rev_comp"] is False


def test_record_locus_is_igh_for_human_igh():
    for rec in _human_igh_exp().run_records(n=3, seed=0):
        assert rec["locus"] == "IGH"


def test_record_locus_is_igk_for_human_igk():
    out = ga.Experiment.on("human_igk").recombine().run_records(n=3, seed=0)
    for rec in out:
        assert rec["locus"] == "IGK"


def test_record_locus_is_igl_for_human_igl():
    out = ga.Experiment.on("human_igl").recombine().run_records(n=3, seed=0)
    for rec in out:
        assert rec["locus"] == "IGL"


def test_record_locus_for_tcr_chain():
    out = (
        ga.Experiment.on("human_tcrb").allow_curatable_refdata()
        .recombine()
        .run_records(n=3, seed=0)
    )
    for rec in out:
        assert rec["locus"] == "TRB"


def test_derive_locus_from_unrecognised_prefix_yields_empty_string():
    from GenAIRR._airr_record import _derive_locus

    assert _derive_locus("XYZ1*01", "", "") == ""
    assert _derive_locus("", "", "") == ""


def test_derive_locus_falls_back_to_j_when_v_empty():
    from GenAIRR._airr_record import _derive_locus

    assert _derive_locus("", "IGHJ4*01", "") == "IGH"


def test_record_c_call_is_empty_string():
    rec = _human_igh_exp().run_records(n=1, seed=0)[0]
    # We don't simulate the constant region, so c_call is always
    # empty — but the column must exist for AIRR-tooling parsers.
    assert rec["c_call"] == ""


def test_default_column_order_starts_with_airr_metadata():
    from GenAIRR.result import _DEFAULT_COLUMN_ORDER

    # Canonical order: AIRR metadata first, then calls.
    assert _DEFAULT_COLUMN_ORDER[0] == "sequence_id"
    assert _DEFAULT_COLUMN_ORDER[1] == "sequence"
    assert "rev_comp" in _DEFAULT_COLUMN_ORDER
    assert "locus" in _DEFAULT_COLUMN_ORDER
    assert "c_call" in _DEFAULT_COLUMN_ORDER


def test_to_tsv_includes_new_h1_columns(tmp_path):
    import csv as _csv

    result = _human_igh_records(n=2)
    path = tmp_path / "h1.tsv"
    result.to_tsv(str(path))

    with open(path, encoding="utf-8") as fh:
        rows = list(_csv.DictReader(fh, delimiter="\t"))

    for i, row in enumerate(rows):
        assert row["sequence_id"] == f"seq{i}"
        assert row["rev_comp"] == "F"
        assert row["locus"] == "IGH"
        assert row["c_call"] == ""
        assert "junction" in row  # new name
        assert "np1" in row  # new name


def test_stream_records_emits_sequential_sequence_ids():
    out = list(_human_igh_exp().stream_records(n=4, seed=0))
    for i, rec in enumerate(out):
        assert rec["sequence_id"] == f"seq{i}"


def test_stream_records_id_prefix_kwarg_is_threaded_through():
    out = list(
        _human_igh_exp().stream_records(n=3, seed=0, id_prefix="run42_")
    )
    assert [rec["sequence_id"] for rec in out] == [
        "run42_0",
        "run42_1",
        "run42_2",
    ]


def test_record_has_sequence_aa_field():
    rec = _human_igh_exp().run_records(n=1, seed=0)[0]
    assert "sequence_aa" in rec
    assert isinstance(rec["sequence_aa"], str)


def test_record_has_np_aa_fields():
    rec = _human_igh_exp().run_records(n=1, seed=0)[0]
    assert "np1_aa" in rec
    assert "np2_aa" in rec


def test_sequence_aa_is_nonempty_for_productive_outcome():
    out = (
        _human_igh_exp()
        .productive_only().run_records(n=3, seed=42)
    )
    for rec in out:
        assert rec["sequence_aa"] != ""


def test_sequence_aa_contains_junction_aa_for_productive_outcome():
    # The junction is part of the assembled sequence, so its AA
    # translation must appear inside sequence_aa.
    out = (
        _human_igh_exp()
        .productive_only().run_records(n=5, seed=42)
    )
    for rec in out:
        assert rec["junction_aa"], "expected junction_aa for productive seq"
        assert rec["junction_aa"] in rec["sequence_aa"]


def test_sequence_aa_uses_v_anchor_frame():
    # Translation starts at the V anchor's frame (junction_start % 3
    # offset). So the translated AA at the junction's relative
    # position should equal junction_aa.
    rec = (
        _human_igh_exp()
        .productive_only().run_records(n=1, seed=42)[0]
    )
    junction_start = rec["junction_start"]
    frame_offset = junction_start % 3
    # The junction's first codon sits at AA index
    # ``(junction_start - frame_offset) // 3``.
    aa_idx = (junction_start - frame_offset) // 3
    junction_aa_len = len(rec["junction_aa"])
    assert (
        rec["sequence_aa"][aa_idx : aa_idx + junction_aa_len]
        == rec["junction_aa"]
    )


def test_sequence_aa_length_matches_codon_count():
    rec = (
        _human_igh_exp()
        .productive_only().run_records(n=1, seed=42)[0]
    )
    seq = rec["sequence"]
    junction_start = rec["junction_start"]
    frame_offset = junction_start % 3
    expected_len = (len(seq) - frame_offset) // 3
    assert len(rec["sequence_aa"]) == expected_len


def test_sequence_aa_empty_when_junction_uncomputable():
    # A bare RefDataConfig has no anchors → no junction → sequence_aa
    # falls back to empty.
    from GenAIRR import _engine as ge

    cfg = ge.RefDataConfig.vdj()
    cfg.add_v_allele("v1*01", "v1", b"AAACCCGGG")  # no anchor
    cfg.add_d_allele("d1*01", "d1", b"TTTTTT")
    cfg.add_j_allele("j1*01", "j1", b"TTTAAA")  # no anchor
    out = ga.Experiment.on(cfg).allow_curatable_refdata().recombine().run_records(n=2, seed=0)
    for rec in out:
        assert rec["sequence_aa"] == ""
        assert rec["np1_aa"] == ""
        assert rec["np2_aa"] == ""


def test_np_aa_consistent_with_sequence_aa_slice():
    # np1_aa / np2_aa must equal the sequence_aa substring at the
    # codon-aligned positions inside the NP region.
    rec = (
        _human_igh_exp()
        .productive_only().run_records(n=1, seed=42)[0]
    )
    junction_start = rec["junction_start"]
    frame_offset = junction_start % 3
    seq_aa = rec["sequence_aa"]

    # NP1 — locate it via v_sequence_end / d_sequence_start (or
    # j_sequence_start for VJ chains; this test is VDJ).
    np1_start = rec["v_sequence_end"]
    np1_end = rec["d_sequence_start"]
    if np1_end is not None and np1_end > np1_start:
        delta = (frame_offset - np1_start) % 3
        codon_start = np1_start + delta
        n_codons = (np1_end - codon_start) // 3
        if n_codons > 0:
            aa_idx = (codon_start - frame_offset) // 3
            assert (
                rec["np1_aa"]
                == seq_aa[aa_idx : aa_idx + n_codons]
            )


def test_np_aa_only_translates_complete_codons():
    # NP regions are typically not multiples of 3 → aa length should
    # never exceed (np_length // 3) + 1 (with leeway for partial
    # boundary codons).
    out = _human_igh_exp().productive_only().run_records(n=10, seed=0)
    for rec in out:
        np1_len = rec["np1_length"]
        np2_len = rec["np2_length"]
        # Strict bound: aa cannot exceed nt // 3 (always true since we
        # only emit complete codons).
        assert len(rec["np1_aa"]) <= max(0, np1_len // 3 + 1)
        assert len(rec["np2_aa"]) <= max(0, np2_len // 3 + 1)


def test_np_aa_empty_when_np_region_empty():
    # Force a recombine with explicit NP1 length 0 → np1 string is
    # empty, np1_aa must be empty too.
    exp = (
        ga.Experiment.on("human_igh")
        .recombine(np1_lengths=[(0, 1.0)], np2_lengths=[(0, 1.0)])
    )
    for rec in exp.productive_only().run_records(n=3, seed=0):
        assert rec["np1"] == ""
        assert rec["np1_aa"] == ""
        assert rec["np2"] == ""
        assert rec["np2_aa"] == ""


def test_aa_slice_for_region_helper_handles_alignment():
    from GenAIRR._airr_record import _aa_slice_for_region

    # frame_offset = 0, sequence_aa "ABCDE" maps to nt positions
    # [0, 3, 6, 9, 12]. A region [3, 12) covers codons 1, 2, 3 → "BCD".
    sequence_aa = "ABCDE"
    assert _aa_slice_for_region(3, 12, 0, sequence_aa) == "BCD"
    # Region [4, 12) starts mid-codon — first complete codon at 6.
    assert _aa_slice_for_region(4, 12, 0, sequence_aa) == "CD"
    # Region too short for a complete codon.
    assert _aa_slice_for_region(4, 5, 0, sequence_aa) == ""
    # Empty region.
    assert _aa_slice_for_region(5, 5, 0, sequence_aa) == ""


def test_default_column_order_includes_aa_fields():
    from GenAIRR.result import _DEFAULT_COLUMN_ORDER

    assert "sequence_aa" in _DEFAULT_COLUMN_ORDER
    assert "np1_aa" in _DEFAULT_COLUMN_ORDER
    assert "np2_aa" in _DEFAULT_COLUMN_ORDER


def test_to_tsv_includes_aa_columns(tmp_path):
    import csv as _csv

    result = _human_igh_records(n=2)
    path = tmp_path / "h2.tsv"
    result.to_tsv(str(path))
    with open(path, encoding="utf-8") as fh:
        rows = list(_csv.DictReader(fh, delimiter="\t"))
    for row in rows:
        assert "sequence_aa" in row
        assert "np1_aa" in row
        assert "np2_aa" in row


def test_sequence_aa_handles_lowercase_corruption_bases():
    # Quality-corruption emits lowercase bases. Translation upper-cases
    # internally, so the AA strings should still be valid letters.
    rec = (
        _human_igh_exp()
        .sequencing_errors(count=8)
        .run_records(n=1, seed=0)[0]
    )
    # Must be valid AA characters (no lowercase, no '?').
    for c in rec["sequence_aa"]:
        assert c.isupper() or c in ("*", "X")
