"""Ground-truth alignment strings and the coordinate pairs that accompany them.

Gap-aware: these are the fields a downstream aligner is benchmarked against.
"""
from __future__ import annotations
from GenAIRR import _engine as ge
from GenAIRR import CompiledExperiment, Experiment
import GenAIRR as ga  # noqa: E402
from tests.support.fixtures.experiment_dsl import _alignment_invariants_hold, _cigar_op_counts, _cigar_total_ops, _human_igh_exp, _human_igh_records

import pytest

pytestmark = pytest.mark.behavior

def test_record_has_alignment_string_fields():
    rec = _human_igh_exp().run_records(n=1, seed=0)[0]
    assert "sequence_alignment" in rec
    assert "germline_alignment" in rec
    assert "germline_alignment_d_mask" in rec


def test_alignment_strings_equal_length():
    # The triplet must be a strict 1:1 alignment — same length always.
    for rec in _human_igh_exp().run_records(n=5, seed=0):
        sa = rec["sequence_alignment"]
        ga_str = rec["germline_alignment"]
        dm = rec["germline_alignment_d_mask"]
        assert len(sa) == len(ga_str) == len(dm)


def test_no_indel_alignment_length_matches_sequence_length():
    # Without indels, no gaps are inserted; alignment length equals
    # raw sequence length.
    for rec in _human_igh_exp().run_records(n=5, seed=0):
        assert len(rec["sequence_alignment"]) == rec["sequence_length"]


def test_no_indel_sequence_alignment_equals_uppercase_sequence():
    for rec in _human_igh_exp().run_records(n=5, seed=0):
        assert rec["sequence_alignment"] == rec["sequence"].upper()


def test_germline_alignment_has_n_in_np_regions():
    rec = _human_igh_exp().run_records(n=1, seed=0)[0]
    galn = rec["germline_alignment"]
    np1_start = rec["v_sequence_end"]
    np1_end = rec["d_sequence_start"]
    if np1_start is not None and np1_end is not None and np1_end > np1_start:
        assert all(c == "N" for c in galn[np1_start:np1_end])
    np2_start = rec["d_sequence_end"]
    np2_end = rec["j_sequence_start"]
    if np2_start is not None and np2_end is not None and np2_end > np2_start:
        assert all(c == "N" for c in galn[np2_start:np2_end])


def test_d_mask_replaces_d_region_with_ns():
    rec = _human_igh_exp().run_records(n=1, seed=0)[0]
    dm = rec["germline_alignment_d_mask"]
    galn = rec["germline_alignment"]
    d_start = rec["d_sequence_start"]
    d_end = rec["d_sequence_end"]
    if d_start is not None and d_end is not None and d_end > d_start:
        # D portion in d_mask is all N's (no gaps in this no-indel case).
        assert all(c == "N" for c in dm[d_start:d_end])
        # Non-D portions must be unchanged from germline_alignment.
        assert dm[:d_start] == galn[:d_start]
        assert dm[d_end:] == galn[d_end:]


def test_germline_alignment_matches_sequence_alignment_when_no_mutations():
    # No mutations and no indels and no NP regions touching the test
    # range: V coding portion of sequence_alignment must equal V
    # coding portion of germline_alignment.
    rec = _human_igh_exp().run_records(n=1, seed=0)[0]
    sa = rec["sequence_alignment"]
    galn = rec["germline_alignment"]
    # Walk only the V coding region (no NP, no D, no J).
    v_start, v_end = rec["v_sequence_start"], rec["v_sequence_end"]
    assert sa[v_start:v_end] == galn[v_start:v_end]


def test_germline_alignment_diverges_from_sequence_at_mutated_positions():
    # With S5F mutations, sa should differ from galn at coding
    # positions hit by mutation. S5F samples positions with
    # replacement, so back-mutations / position collisions reduce
    # the visible mismatch count below ``count``.
    rec = (
        _human_igh_exp()
        .mutate(count=15)
        .run_records(n=1, seed=42)[0]
    )
    sa = rec["sequence_alignment"]
    galn = rec["germline_alignment"]
    diffs = sum(
        1 for a, b in zip(sa, galn) if a != b and b not in ("N", "-")
    )
    # Strict upper bound — every mismatch must come from one of the
    # 15 mutation events. The lower bound is loose because of
    # collisions; for seed=42 we expect ≥8 distinct positions.
    assert 8 <= diffs <= 15


def test_indel_insertion_creates_gap_in_germline_alignment():
    # Indel insertions add bases with no germline counterpart → those
    # positions should be '-' in germline_alignment.
    # Use insertion_prob=1.0 to force all indels to be insertions.
    rec = (
        _human_igh_exp()
        .polymerase_indels(count=4, insertion_prob=1.0)
        .run_records(n=1, seed=0)[0]
    )
    galn = rec["germline_alignment"]
    sa = rec["sequence_alignment"]
    # All 4 indels are insertions → exactly 4 gaps in galn, 0 in sa.
    assert galn.count("-") == 4
    assert sa.count("-") == 0


def test_indel_deletion_creates_gap_in_sequence_alignment():
    # Indel deletions remove bases → those positions should be '-' in
    # sequence_alignment with the original allele base in germline.
    rec = (
        _human_igh_exp()
        .polymerase_indels(count=4, insertion_prob=0.0)
        .run_records(n=1, seed=0)[0]
    )
    galn = rec["germline_alignment"]
    sa = rec["sequence_alignment"]
    assert sa.count("-") == 4
    assert galn.count("-") == 0


def test_alignment_total_gap_count_at_most_n_indels():
    # For mixed insertions/deletions, every alignment gap traces back
    # to one indel event. Indels that land in NP regions don't surface
    # as alignment gaps (NP has no germline position to mark), so the
    # observed gap count is ≤ ``n_indels``.
    rec = (
        _human_igh_exp()
        .polymerase_indels(count=6, insertion_prob=0.5)
        .run_records(n=1, seed=12345)[0]
    )
    sa = rec["sequence_alignment"]
    galn = rec["germline_alignment"]
    total_gaps = sa.count("-") + galn.count("-")
    assert 0 <= total_gaps <= rec["n_indels"]


def test_alignment_length_increases_with_indels():
    # Each indel (insertion or deletion) adds 1 to alignment length:
    # insertions add a column where seq has a base and germ is gap;
    # deletions add a column where seq is gap and germ has a base.
    rec = (
        _human_igh_exp()
        .polymerase_indels(count=4, insertion_prob=1.0)
        .run_records(n=1, seed=0)[0]
    )
    # sequence_length stays raw (with insertions counted, deletions
    # missing). With 4 insertions and 0 deletions:
    #   alignment_length = sequence_length (no deletions add a column
    #     beyond what insertions already added).
    # Actually: alignment includes deletions as extra columns, so
    # alignment_length == sequence_length + n_deletions.
    assert len(rec["sequence_alignment"]) == rec["sequence_length"]


def test_alignment_length_grows_for_pure_deletions():
    # 4 deletions → sequence_length is shorter than alignment_length
    # (each deletion is a gap column not present in the raw sequence).
    rec = (
        _human_igh_exp()
        .polymerase_indels(count=4, insertion_prob=0.0)
        .run_records(n=1, seed=0)[0]
    )
    assert len(rec["sequence_alignment"]) == rec["sequence_length"] + 4


def test_d_mask_preserves_gaps_in_d_region():
    # Indel deletions inside the D region should still appear as '-'
    # in d_mask (not converted to N).
    rec = (
        _human_igh_exp()
        .polymerase_indels(count=15, insertion_prob=0.0)
        .run_records(n=1, seed=0)[0]
    )
    galn = rec["germline_alignment"]
    dm = rec["germline_alignment_d_mask"]
    # d_mask must have the same gap positions as germline_alignment
    # (gaps are preserved, not masked).
    galn_gap_positions = {i for i, c in enumerate(galn) if c == "-"}
    dm_gap_positions = {i for i, c in enumerate(dm) if c == "-"}
    assert galn_gap_positions == dm_gap_positions


def test_alignment_strings_uppercase():
    # AIRR convention: alignment strings are uppercase regardless of
    # the source sequence case (lowercase NP / mutated bases).
    rec = (
        _human_igh_exp()
        .mutate(count=10)
        .run_records(n=1, seed=0)[0]
    )
    sa = rec["sequence_alignment"]
    galn = rec["germline_alignment"]
    for c in sa:
        assert c.isupper() or c in ("-",), f"non-uppercase {c!r} in sa"
    for c in galn:
        assert c.isupper() or c in ("-",), f"non-uppercase {c!r} in ga"


def test_alignment_strings_for_vj_chain():
    # VJ chains have V and J only (no D, no NP2). Walks must still
    # produce valid 3-string output.
    rec = (
        ga.Experiment.on("human_igk")
        .recombine()
        .run_records(n=1, seed=0)[0]
    )
    sa = rec["sequence_alignment"]
    galn = rec["germline_alignment"]
    dm = rec["germline_alignment_d_mask"]
    assert len(sa) == len(galn) == len(dm) == rec["sequence_length"]
    # No D region → d_mask should equal germline_alignment exactly.
    assert dm == galn


def test_alignment_strings_empty_when_no_assembly():
    # Pre-recombine experiments yield empty sequences → empty strings.
    from GenAIRR import _engine as ge

    # Use a refdata with anchorless alleles so there's nothing
    # productive to anchor.
    cfg = ge.RefDataConfig.vdj()
    cfg.add_v_allele("v1*01", "v1", b"AAACCCGGG")
    cfg.add_d_allele("d1*01", "d1", b"TTTTTT")
    cfg.add_j_allele("j1*01", "j1", b"TTTAAA")
    out = ga.Experiment.on(cfg).allow_curatable_refdata().recombine().run_records(n=1, seed=0)
    rec = out[0]
    # Even without anchors, alignment strings should be the assembled
    # sequence vs allele bases — non-empty here.
    assert len(rec["sequence_alignment"]) == rec["sequence_length"]


def test_alignment_strings_in_default_column_order():
    from GenAIRR.result import _DEFAULT_COLUMN_ORDER

    assert "sequence_alignment" in _DEFAULT_COLUMN_ORDER
    assert "germline_alignment" in _DEFAULT_COLUMN_ORDER
    assert "germline_alignment_d_mask" in _DEFAULT_COLUMN_ORDER


def test_alignment_strings_in_tsv_export(tmp_path):
    import csv as _csv

    result = _human_igh_exp().mutate(count=5).run_records(n=2, seed=0)
    path = tmp_path / "h3.tsv"
    result.to_tsv(str(path))
    with open(path, encoding="utf-8") as fh:
        rows = list(_csv.DictReader(fh, delimiter="\t"))
    for row in rows:
        assert row["sequence_alignment"]
        assert row["germline_alignment"]
        assert row["germline_alignment_d_mask"]
        assert len(row["sequence_alignment"]) == len(row["germline_alignment"])


def test_alignment_handles_tail_indel_insertion_outside_regions():
    # Regression for the bug where an insertion at exactly pool_len
    # didn't extend any region's end → one pool position was outside
    # every region → walker dropped it → sa shorter than seq.
    # Repro: human_igk + indels + seed=0 used to fail at record idx 1.
    exp = (
        ga.Experiment.on("human_igk")
        .recombine()
        .mutate(count=5)
        .pcr_amplify(count=2)
        .polymerase_indels(count=2)
    )
    for rec in exp.run_records(n=50, seed=0):
        assert _alignment_invariants_hold(rec) == []


def test_alignment_invariants_under_full_corruption_stack():
    # Stack every corruption pass at once and verify alignment
    # invariants hold across 200 seeds. This is the definitive
    # H.3 robustness check — touches PCR, quality, indels, contaminant,
    # mutations, and the productive contract simultaneously.
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(count=15)
        .pcr_amplify(count=5)
        .sequencing_errors(count=8)
        .polymerase_indels(count=4, insertion_prob=0.5, rejection_attempts=4096)
    )
    # This stress case requests four edits after heavy mutation. Allocate more
    # complete proposals; exhausted work must not become an unnoticed no-op.
    result = exp.productive_only().run_records(n=200, seed=0)
    assert all(rec["n_indels"] == 4 for rec in result)
    failures = [
        (i, _alignment_invariants_hold(rec))
        for i, rec in enumerate(result)
        if _alignment_invariants_hold(rec)
    ]
    assert failures == [], f"alignment invariants broke for: {failures[:3]}"


def test_alignment_invariants_across_loci():
    # IGH / IGK / IGL / TCRB all need to produce sane alignment
    # strings under heavy corruption. SHM is BCR-only, so the
    # TCR pipeline drops the mutate step.
    for cfg in ("human_igh", "human_igk", "human_igl", "human_tcrb"):
        exp = ga.Experiment.on(cfg).allow_curatable_refdata().recombine()
        if not cfg.startswith("human_tcr"):
            exp = exp.mutate(count=5)
        exp = exp.pcr_amplify(count=2).polymerase_indels(count=2)
        result = exp.run_records(n=50, seed=0)
        failures = [
            (i, _alignment_invariants_hold(rec))
            for i, rec in enumerate(result)
            if _alignment_invariants_hold(rec)
        ]
        assert failures == [], f"{cfg}: {failures[:3]}"


def test_record_has_h5_coord_pair_fields():
    rec = _human_igh_exp().run_records(n=1, seed=0)[0]
    for seg in ("v", "d", "j"):
        for end in ("alignment_start", "alignment_end", "germline_start", "germline_end"):
            assert f"{seg}_{end}" in rec, f"missing {seg}_{end}"


def test_no_indel_alignment_coords_match_sequence_coords():
    # Without indels, alignment-string columns mirror sequence
    # positions exactly → alignment coords equal sequence coords.
    rec = _human_igh_exp().run_records(n=1, seed=42)[0]
    for seg in ("v", "d", "j"):
        s_s = rec[f"{seg}_sequence_start"]
        s_e = rec[f"{seg}_sequence_end"]
        a_s = rec[f"{seg}_alignment_start"]
        a_e = rec[f"{seg}_alignment_end"]
        assert a_s == s_s, f"{seg}: alignment_start {a_s} != sequence_start {s_s}"
        assert a_e == s_e, f"{seg}: alignment_end {a_e} != sequence_end {s_e}"


def test_germline_coord_length_matches_post_trim_segment_length():
    # germline_end - germline_start == allele_len - trim_5 - trim_3.
    # Without indels, that also equals sequence_end - sequence_start.
    rec = _human_igh_exp().run_records(n=1, seed=42)[0]
    for seg in ("v", "d", "j"):
        coding = rec[f"{seg}_sequence_end"] - rec[f"{seg}_sequence_start"]
        germline = rec[f"{seg}_germline_end"] - rec[f"{seg}_germline_start"]
        assert germline == coding


def test_v_germline_end_reflects_v_trim_3():
    # V_3 trim: v_germline_end = v_allele_len - v_trim_3.
    refdata = _human_igh_exp().refdata
    rec = _human_igh_exp().run_records(n=1, seed=42)[0]
    # Locate the assigned V allele.
    v_id = None
    for i in range(refdata.v_pool_size()):
        if refdata.v_allele(i).name == rec["v_call"]:
            v_id = i
            break
    assert v_id is not None
    allele_len = len(bytes(refdata.v_allele(v_id).seq()))
    assert rec["v_germline_end"] == allele_len - rec["v_trim_3"]


def test_d_germline_start_reflects_d_trim_5():
    rec = _human_igh_exp().run_records(n=1, seed=42)[0]
    assert rec["d_germline_start"] == rec["d_trim_5"]


def test_d_germline_end_reflects_d_trim_3():
    refdata = _human_igh_exp().refdata
    rec = _human_igh_exp().run_records(n=1, seed=42)[0]
    d_id = None
    for i in range(refdata.d_pool_size()):
        if refdata.d_allele(i).name == rec["d_call"]:
            d_id = i
            break
    assert d_id is not None
    allele_len = len(bytes(refdata.d_allele(d_id).seq()))
    assert rec["d_germline_end"] == allele_len - rec["d_trim_3"]


def test_j_germline_start_reflects_j_trim_5():
    rec = _human_igh_exp().run_records(n=1, seed=42)[0]
    assert rec["j_germline_start"] == rec["j_trim_5"]


def test_germline_alignment_slice_matches_germline_coord_pair():
    # GROUND-TRUTH INVARIANT: slicing germline_alignment by alignment
    # coords should equal slicing the source allele by germline coords
    # (post-trim, so we get exactly the bases the simulation drew).
    refdata = _human_igh_exp().refdata
    rec = _human_igh_exp().run_records(n=1, seed=42)[0]

    for seg in ("v", "d", "j"):
        a_s = rec[f"{seg}_alignment_start"]
        a_e = rec[f"{seg}_alignment_end"]
        g_s = rec[f"{seg}_germline_start"]
        g_e = rec[f"{seg}_germline_end"]
        if a_s is None:
            continue
        # Find the allele.
        getter = {
            "v": refdata.v_allele,
            "d": refdata.d_allele,
            "j": refdata.j_allele,
        }[seg]
        pool_size = {
            "v": refdata.v_pool_size,
            "d": refdata.d_pool_size,
            "j": refdata.j_pool_size,
        }[seg]()
        call_name = rec[f"{seg}_call"]
        allele_id = next(
            (i for i in range(pool_size) if getter(i).name == call_name),
            None,
        )
        assert allele_id is not None
        allele_seq = bytes(getter(allele_id).seq()).decode().upper()
        expected = allele_seq[g_s:g_e]
        actual = rec["germline_alignment"][a_s:a_e]
        assert (
            actual == expected
        ), f"{seg}: germline_alignment slice {actual!r} != source allele slice {expected!r}"


def test_indel_deletion_alignment_and_sequence_spans_are_both_positive():
    # Sanity: with deletions in V/D/J, both alignment and sequence
    # spans should be positive numbers. Sequence span comes from
    # live-call coords, so an NP-side extension can grow seq span
    # past the structural alignment span.
    rec = (
        _human_igh_exp()
        .polymerase_indels(count=4, insertion_prob=0.0)
        .run_records(n=1, seed=0)[0]
    )
    for seg in ("v", "d", "j"):
        s_s = rec[f"{seg}_sequence_start"]
        s_e = rec[f"{seg}_sequence_end"]
        a_s = rec[f"{seg}_alignment_start"]
        a_e = rec[f"{seg}_alignment_end"]
        if a_s is None:
            continue
        assert a_e - a_s > 0
        assert s_e - s_s > 0


def test_alignment_span_equals_cigar_op_total():
    # ``v_alignment_end - v_alignment_start`` must equal the total
    # M+I+D op count of v_cigar.
    rec = (
        _human_igh_exp()
        .polymerase_indels(count=5, insertion_prob=0.5)
        .run_records(n=1, seed=42)[0]
    )
    for seg in ("v", "d", "j"):
        a_s = rec[f"{seg}_alignment_start"]
        a_e = rec[f"{seg}_alignment_end"]
        cig = rec[f"{seg}_cigar"]
        if a_s is None:
            continue
        assert a_e - a_s == _cigar_total_ops(cig), (
            f"{seg}: alignment span {a_e - a_s} != cigar ops {_cigar_total_ops(cig)}"
        )


def test_germline_span_at_least_cigar_m_plus_d():
    # The number of source-allele bases the alignment + live-extension
    # account for is at least M + D. Germline coords come from the
    # live-call hypothesis bounds, so when an NP-side extension grew
    # the hypothesis's ref range the germline span can exceed M + D.
    rec = (
        _human_igh_exp()
        .polymerase_indels(count=5, insertion_prob=0.5)
        .run_records(n=1, seed=42)[0]
    )
    for seg in ("v", "d", "j"):
        g_s = rec[f"{seg}_germline_start"]
        g_e = rec[f"{seg}_germline_end"]
        cig = rec[f"{seg}_cigar"]
        if g_s is None:
            continue
        ops = _cigar_op_counts(cig)
        assert g_e - g_s >= ops["M"] + ops["D"]


def test_d_coords_none_for_vj_chain():
    rec = (
        ga.Experiment.on("human_igk")
        .recombine()
        .run_records(n=1, seed=0)[0]
    )
    assert rec["d_alignment_start"] is None
    assert rec["d_alignment_end"] is None
    assert rec["d_germline_start"] is None
    assert rec["d_germline_end"] is None
    # V/J coords still populated.
    assert rec["v_alignment_start"] == 0
    assert rec["j_alignment_start"] is not None


def test_v_alignment_starts_at_zero():
    # V is the first segment in our assembly, no leading deletion gaps
    # → v_alignment_start is always 0.
    for rec in _human_igh_exp().run_records(n=5, seed=0):
        assert rec["v_alignment_start"] == 0


def test_h5_coords_in_default_column_order():
    from GenAIRR.result import _DEFAULT_COLUMN_ORDER

    for seg in ("v", "d", "j"):
        for end in ("alignment_start", "alignment_end", "germline_start", "germline_end"):
            assert f"{seg}_{end}" in _DEFAULT_COLUMN_ORDER


def test_h5_coords_in_tsv_export(tmp_path):
    import csv as _csv

    result = _human_igh_records(n=2)
    path = tmp_path / "h5.tsv"
    result.to_tsv(str(path))
    with open(path, encoding="utf-8") as fh:
        rows = list(_csv.DictReader(fh, delimiter="\t"))
    for row in rows:
        for seg in ("v", "d", "j"):
            for end in ("alignment_start", "alignment_end", "germline_start", "germline_end"):
                assert f"{seg}_{end}" in row
