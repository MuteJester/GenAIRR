"""Compatibility with airr.validate_rearrangement, including strict mode.
"""
from __future__ import annotations
import pytest
import GenAIRR as ga  # noqa: E402
import re as _re
from tests.support.fixtures.experiment_dsl import _alignment_invariants_hold, _cigar_op_counts, _human_igh_exp

pytestmark = pytest.mark.behavior

def test_to_tsv_passes_airr_validate_rearrangement(tmp_path):
    airr = pytest.importorskip("airr")
    rec_path = tmp_path / "valid.tsv"
    result = (
        _human_igh_exp().mutate(count=5).run_records(n=5, seed=0)
    )
    result.to_tsv(str(rec_path))
    assert airr.validate_rearrangement(str(rec_path)) is True


def test_to_tsv_passes_airr_validation_under_full_corruption_stack(tmp_path):
    airr = pytest.importorskip("airr")
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(count=15)
        .pcr_amplify(count=5)
        .sequencing_errors(count=8)
        .polymerase_indels(count=4, insertion_prob=0.5)
    )
    path = tmp_path / "stack.tsv"
    exp.run_records(n=20, seed=0).to_tsv(str(path))
    assert airr.validate_rearrangement(str(path)) is True


def test_to_tsv_passes_airr_validation_for_vj_chain(tmp_path):
    airr = pytest.importorskip("airr")
    path = tmp_path / "vj.tsv"
    (
        ga.Experiment.on("human_igk")
        .recombine()
        .run_records(n=5, seed=0)
        .to_tsv(str(path))
    )
    assert airr.validate_rearrangement(str(path)) is True


def test_to_tsv_passes_airr_validation_for_tcrb(tmp_path):
    airr = pytest.importorskip("airr")
    path = tmp_path / "tcrb.tsv"
    # TCR doesn't undergo SHM; use PCR errors for sequencing
    # realism instead.
    (
        ga.Experiment.on("human_tcrb").allow_curatable_refdata()
        .recombine()
        .pcr_amplify(count=3)
        .run_records(n=5, seed=0)
        .to_tsv(str(path))
    )
    assert airr.validate_rearrangement(str(path)) is True


def test_to_tsv_passes_airr_validation_with_contamination(tmp_path):
    airr = pytest.importorskip("airr")
    path = tmp_path / "contam.tsv"
    (
        _human_igh_exp()
        .contaminate(prob=0.5)
        .run_records(n=5, seed=0)
        .to_tsv(str(path))
    )
    assert airr.validate_rearrangement(str(path)) is True


def test_airr_strict_increments_start_fields_only(tmp_path):
    # Explicit legacy mode keeps coords 0-based half-open. ``airr_strict=True``
    # bumps every ``*_start`` field by 1 (giving 1-based-inclusive).
    # The ``*_end`` fields are left unchanged because Python's
    # half-open ``end`` already equals the 1-based-inclusive ``end``.
    import csv as _csv

    result = _human_igh_exp().run_records(n=2, seed=42)
    py_path = tmp_path / "py.tsv"
    strict_path = tmp_path / "strict.tsv"
    result.to_tsv(str(py_path), airr_strict=False)
    result.to_tsv(str(strict_path), airr_strict=True)

    with open(py_path) as fh:
        py_row = next(_csv.DictReader(fh, delimiter="\t"))
    with open(strict_path) as fh:
        strict_row = next(_csv.DictReader(fh, delimiter="\t"))

    # Every *_start coord field should be off by exactly 1.
    for field in (
        "v_sequence_start",
        "d_sequence_start",
        "j_sequence_start",
        "v_alignment_start",
        "d_alignment_start",
        "j_alignment_start",
        "v_germline_start",
        "d_germline_start",
        "j_germline_start",
        "junction_start",
    ):
        if py_row[field] == "" or strict_row[field] == "":
            continue
        assert int(strict_row[field]) == int(py_row[field]) + 1, (
            f"{field}: strict {strict_row[field]} != default {py_row[field]} + 1"
        )

    # Every *_end coord field should be unchanged.
    for field in (
        "v_sequence_end",
        "d_sequence_end",
        "j_sequence_end",
        "v_alignment_end",
        "d_alignment_end",
        "j_alignment_end",
        "v_germline_end",
        "d_germline_end",
        "j_germline_end",
        "junction_end",
    ):
        assert py_row[field] == strict_row[field], (
            f"{field}: strict {strict_row[field]} != default {py_row[field]}"
        )


def test_airr_strict_passes_airr_validation(tmp_path):
    airr = pytest.importorskip("airr")
    path = tmp_path / "strict.tsv"
    result = (
        _human_igh_exp()
        .mutate(count=10)
        .polymerase_indels(count=3, insertion_prob=0.5)
        .run_records(n=5, seed=0)
    )
    result.to_tsv(str(path), airr_strict=True)
    assert airr.validate_rearrangement(str(path)) is True


def test_airr_strict_csv_export(tmp_path):
    import csv as _csv

    result = _human_igh_exp().run_records(n=1, seed=42)
    py_path = tmp_path / "py.csv"
    strict_path = tmp_path / "strict.csv"
    result.to_csv(str(py_path))
    result.to_csv(str(strict_path), airr_strict=True)

    with open(py_path) as fh:
        py_row = next(_csv.DictReader(fh))
    with open(strict_path) as fh:
        strict_row = next(_csv.DictReader(fh))

    assert int(strict_row["v_sequence_start"]) == int(py_row["v_sequence_start"]) + 1


def test_airr_strict_dataframe():
    pd = pytest.importorskip("pandas")
    result = _human_igh_exp().run_records(n=2, seed=42)

    df_default = result.to_dataframe()
    df_strict = result.to_dataframe(airr_strict=True)

    assert isinstance(df_default, pd.DataFrame)
    assert isinstance(df_strict, pd.DataFrame)
    # Every *_start in strict must equal default + 1.
    for col in (
        "v_sequence_start",
        "j_sequence_start",
        "junction_start",
        "v_alignment_start",
        "v_germline_start",
    ):
        for i in range(len(df_default)):
            d = df_default.iloc[i][col]
            s = df_strict.iloc[i][col]
            if d is None or (hasattr(d, "is_integer") and not d == d):
                continue
            assert s == d + 1, f"{col}[{i}]: {s} != {d} + 1"


def test_to_airr_strict_helper_keeps_none_starts():
    # If a *_start is ``None`` (e.g. d_alignment_start on a VJ chain),
    # the helper must leave it as None.
    from GenAIRR.result import _to_airr_strict

    rec = {
        "v_sequence_start": 0,
        "v_sequence_end": 100,
        "d_sequence_start": None,
        "d_sequence_end": None,
        "junction_start": None,
        "j_sequence_start": 200,
        "j_sequence_end": 250,
    }
    converted = _to_airr_strict(rec)
    assert converted["v_sequence_start"] == 1
    assert converted["v_sequence_end"] == 100
    assert converted["d_sequence_start"] is None
    assert converted["junction_start"] is None
    assert converted["j_sequence_start"] == 201


def test_to_airr_strict_helper_does_not_mutate_input():
    from GenAIRR.result import _to_airr_strict

    rec = {"v_sequence_start": 0, "v_sequence_end": 100}
    snapshot = dict(rec)
    _to_airr_strict(rec)
    assert rec == snapshot


def test_airr_strict_tsv_default_is_true():
    # Sanity: ``to_tsv`` without the flag should produce the same
    # output as ``airr_strict=True``.
    import tempfile, os, csv as _csv

    result = _human_igh_exp().run_records(n=1, seed=0)
    with tempfile.TemporaryDirectory() as tmp:
        a = os.path.join(tmp, "a.tsv")
        b = os.path.join(tmp, "b.tsv")
        result.to_tsv(a)
        result.to_tsv(b, airr_strict=True)
        with open(a) as fh:
            ra = list(_csv.DictReader(fh, delimiter="\t"))
        with open(b) as fh:
            rb = list(_csv.DictReader(fh, delimiter="\t"))
    assert ra == rb


def test_identity_helper_from_columns():
    from GenAIRR._airr_record import _alignment_columns, _identity_from_columns

    refdata = _human_igh_exp().refdata
    out = _human_igh_exp().mutate(count=10).run(n=1, seed=42)[0]
    cols = _alignment_columns(out, refdata)
    identities = _identity_from_columns(cols)
    # Same as the public field — guard against silent re-routing.
    rec = _human_igh_exp().mutate(count=10).run_records(n=1, seed=42)[0]
    assert identities["V"] == rec["v_identity"]
    assert identities["D"] == rec["d_identity"]
    assert identities["J"] == rec["j_identity"]


_RC_COMPLEMENT = str.maketrans("ACGTacgtNnUu", "TGCAtgcaNnAa")


def _reverse_complement(s: str) -> str:
    return s.translate(_RC_COMPLEMENT)[::-1]


def test_5prime_loss_drops_leading_bases():
    # Fixed-length 5' loss: every record's sequence is shorter by
    # exactly the requested amount. Reference coverage moves with the read.
    base = ga.Experiment.on("human_igh").recombine()
    lost = ga.Experiment.on("human_igh").recombine().primer_trim_5prime(length=15)
    base_recs = list(base.run_records(n=10, seed=0))
    lost_recs = list(lost.run_records(n=10, seed=0))
    for fwd, lr in zip(base_recs, lost_recs):
        # Sequence shrinks by exactly 15.
        assert lr["sequence_length"] == fwd["sequence_length"] - 15, (
            f"len={lr['sequence_length']} vs {fwd['sequence_length']} - 15"
        )
        assert lr["sequence"] == fwd["sequence"][15:]
        assert lr["v_germline_start"] == fwd["v_germline_start"] + 15
        assert "D" not in lr["v_cigar"]
        assert lr["v_identity"] == 1.0


def test_3prime_loss_drops_trailing_bases():
    base = ga.Experiment.on("human_igh").recombine()
    lost = ga.Experiment.on("human_igh").recombine().primer_trim_3prime(length=10)
    for fwd, lr in zip(base.run_records(n=10, seed=0), lost.run_records(n=10, seed=0)):
        assert lr["sequence_length"] == fwd["sequence_length"] - 10
        assert lr["sequence"] == fwd["sequence"][:-10]
        assert lr["j_germline_end"] == fwd["j_germline_end"] - 10
        assert "D" not in lr["j_cigar"]
        assert lr["j_identity"] == 1.0


def test_loss_zero_is_no_op():
    base = ga.Experiment.on("human_igh").recombine()
    lost_5 = ga.Experiment.on("human_igh").recombine().primer_trim_5prime(length=0)
    lost_3 = ga.Experiment.on("human_igh").recombine().primer_trim_3prime(length=0)
    for fwd, l5, l3 in zip(
        base.run_records(n=5, seed=0),
        lost_5.run_records(n=5, seed=0),
        lost_3.run_records(n=5, seed=0),
    ):
        assert l5["sequence"] == fwd["sequence"]
        assert l3["sequence"] == fwd["sequence"]


def test_loss_preserves_alignment_invariants():
    # Clipping shortens coverage while the alignment still preserves query
    # order and the observed M/I spans.
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(count=8)
        .primer_trim_5prime(length=(0, 25))
        .primer_trim_3prime(length=(0, 20))
        .pcr_amplify(count=3)
    )
    failures = []
    for i, rec in enumerate(exp.run_records(n=100, seed=0)):
        issues = _alignment_invariants_hold(rec)
        if issues:
            failures.append((i, issues))
    assert failures == [], f"alignment invariants broke under 5'/3' loss: {failures[:3]}"


def test_loss_rejects_invalid_length():
    exp = ga.Experiment.on("human_igh").recombine()
    for bad in (-1, "ten", True):
        try:
            exp.primer_trim_5prime(length=bad)
        except (ValueError, TypeError):
            pass
        else:
            raise AssertionError(f"expected error for length={bad!r}")


def test_rev_comp_flips_sequence_and_flag():
    # prob=1.0: every record should have rev_comp=True and the
    # `sequence` should be the reverse-complement of the same record
    # without rev_comp applied.
    base = ga.Experiment.on("human_igh").recombine()
    flipped = ga.Experiment.on("human_igh").recombine().random_strand_orientation(prob=1.0)
    base_recs = list(base.run_records(n=20, seed=0))
    flipped_recs = list(flipped.run_records(n=20, seed=0))
    assert len(base_recs) == len(flipped_recs) == 20
    for fwd, rc in zip(base_recs, flipped_recs):
        assert fwd["rev_comp"] is False
        assert rc["rev_comp"] is True
        assert rc["sequence"] == _reverse_complement(fwd["sequence"])


def test_rev_comp_preserves_coding_oriented_sequence_coords():
    # AIRR coordinates refer to the coding-oriented reverse complement
    # of sequence when rev_comp is true.
    base = ga.Experiment.on("human_igh").recombine()
    flipped = ga.Experiment.on("human_igh").recombine().random_strand_orientation(prob=1.0)
    for fwd, rc in zip(base.run_records(n=20, seed=0), flipped.run_records(n=20, seed=0)):
        for seg in ("v", "d", "j"):
            fs, fe = fwd[f"{seg}_sequence_start"], fwd[f"{seg}_sequence_end"]
            if fs is None:
                continue
            rs, re = rc[f"{seg}_sequence_start"], rc[f"{seg}_sequence_end"]
            assert (rs, re) == (fs, fe)
            assert _reverse_complement(rc["sequence"])[rs:re] == fwd["sequence"][fs:fe]


def test_rev_comp_keeps_alignment_strings_in_forward_orientation():
    # AIRR-spec: `sequence_alignment` and `germline_alignment` stay
    # forward-oriented even when `rev_comp=True`. The CIGAR and the
    # `*_alignment_start/end` and `*_germline_start/end` coords are
    # also forward.
    base = ga.Experiment.on("human_igh").recombine()
    flipped = ga.Experiment.on("human_igh").recombine().random_strand_orientation(prob=1.0)
    for fwd, rc in zip(base.run_records(n=10, seed=0), flipped.run_records(n=10, seed=0)):
        assert rc["sequence_alignment"] == fwd["sequence_alignment"]
        assert rc["germline_alignment"] == fwd["germline_alignment"]
        assert rc["v_cigar"] == fwd["v_cigar"]
        assert rc["d_cigar"] == fwd["d_cigar"]
        assert rc["j_cigar"] == fwd["j_cigar"]
        for seg in ("v", "d", "j"):
            assert rc[f"{seg}_alignment_start"] == fwd[f"{seg}_alignment_start"]
            assert rc[f"{seg}_alignment_end"] == fwd[f"{seg}_alignment_end"]
            assert rc[f"{seg}_germline_start"] == fwd[f"{seg}_germline_start"]
            assert rc[f"{seg}_germline_end"] == fwd[f"{seg}_germline_end"]


def test_rev_comp_zero_prob_is_no_op():
    base = ga.Experiment.on("human_igh").recombine()
    flipped = ga.Experiment.on("human_igh").recombine().random_strand_orientation(prob=0.0)
    for fwd, rc in zip(base.run_records(n=20, seed=0), flipped.run_records(n=20, seed=0)):
        # The coin is flipped (records always present in trace) but
        # never fires when prob=0.
        assert rc["rev_comp"] is False
        assert rc["sequence"] == fwd["sequence"]


def test_rev_comp_rejects_invalid_prob():
    exp = ga.Experiment.on("human_igh").recombine()
    import math
    for bad in (-0.1, 1.5, math.nan, float("inf")):
        try:
            exp.random_strand_orientation(prob=bad)
        except ValueError:
            pass
        else:
            raise AssertionError(f"expected ValueError for prob={bad}")


def test_locus_falls_back_to_refdata_under_heavy_corruption():
    # Under heavy SHM the live-call layer can narrow every V/D/J
    # call to the empty set (no allele supports the corrupted
    # sequence). When that happens `derive_locus` has nothing to
    # parse and would return `""`. Walk the refdata pool's first
    # allele names as a fallback so the AIRR record always carries
    # a meaningful locus.
    expected = {
        "human_igh": "IGH",
        "human_igk": "IGK",
        "human_igl": "IGL",
        "human_tcrb": "TRB",
    }
    for cfg, want in expected.items():
        # SHM is BCR-only; drop mutate for TCR.
        exp = ga.Experiment.on(cfg).allow_curatable_refdata().recombine()
        if not cfg.startswith("human_tcr"):
            exp = exp.mutate(count=15)
        exp = (
            exp.pcr_amplify(count=5)
            .sequencing_errors(count=8)
            .polymerase_indels(count=4, insertion_prob=0.5)
        )
        for i, rec in enumerate(exp.run_records(n=200, seed=0)):
            assert rec["locus"] == want, (
                f"{cfg} rec {i}: locus={rec['locus']!r} expected {want!r} "
                f"(v_call={rec['v_call']!r} j_call={rec['j_call']!r})"
            )


def test_per_segment_ga_slice_matches_called_allele():
    # Audit invariant: when SHM mutations narrow the live call to a
    # different allele than the one originally sampled,
    # `germline_alignment[v_align_start:v_align_end].replace('-','')`
    # must equal the *called* allele's bytes at
    # `[v_germline_start, v_germline_end)`. Same for D and J.
    #
    # The column-walker germline byte is routed through the
    # projected (live-call's first) allele rather than the original
    # `nuc.germline` provenance byte from AssembleSegmentPass — the
    # latter caused a divergence in ~1% of records under heavy
    # corruption.
    pipelines = [
        ("human_igh", "vdj"),
        ("human_igk", "vj"),
        ("human_igl", "vj"),
        ("human_tcrb", "vdj"),
    ]
    failures = []
    for cfg, _kind in pipelines:
        # SHM is BCR-only; drop mutate for TCR.
        exp = ga.Experiment.on(cfg).allow_curatable_refdata().recombine()
        if not cfg.startswith("human_tcr"):
            exp = exp.mutate(count=8)
        exp = (
            exp.pcr_amplify(count=2)
            .sequencing_errors(count=4)
            .polymerase_indels(count=2, insertion_prob=0.5)
        )
        refdata = exp.refdata
        alleles = {"V": {}, "D": {}, "J": {}}
        for seg, n_attr, get in [
            ("V", refdata.v_pool_size, refdata.v_allele),
            ("D", refdata.d_pool_size, refdata.d_allele),
            ("J", refdata.j_pool_size, refdata.j_allele),
        ]:
            for i in range(n_attr()):
                a = get(i)
                alleles[seg][a.name] = a.seq().decode("ascii").upper()
        for i, rec in enumerate(exp.run_records(n=300, seed=0)):
            for seg in ("V", "D", "J"):
                calls = rec[f"{seg.lower()}_call"]
                if not calls:
                    continue
                a_s = rec[f"{seg.lower()}_alignment_start"]
                a_e = rec[f"{seg.lower()}_alignment_end"]
                g_s = rec[f"{seg.lower()}_germline_start"]
                g_e = rec[f"{seg.lower()}_germline_end"]
                if a_s is None:
                    continue
                first_call = calls.split(",")[0]
                ref = alleles[seg].get(first_call)
                if ref is None:
                    continue
                ga_str = rec["germline_alignment"][a_s:a_e].replace("-", "")
                if ga_str.upper() != ref[g_s:g_e].upper():
                    failures.append((cfg, i, seg, first_call, ga_str, ref[g_s:g_e]))
    assert failures == [], f"per-segment ga slice mismatches: {failures[:3]}"


def test_h5_coords_invariants_under_corruption_stack():
    # Same n=200 stress sweep as alignment/CIGAR — verify the coord
    # pairs stay self-consistent under every corruption mode.
    #
    # Strict identities: `*_germline_start/end` are derived from the
    # column walker's `ref_ranges` (the union of ref positions
    # consumed by `M` and `D` ops), which enforces `germline_span ==
    # M + D`. Combined with the alignment-span identity
    # `alignment_span == M + I + D`, the AIRR coord triples are
    # fully self-consistent.
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
        for seg in ("v", "d", "j"):
            a_s = rec[f"{seg}_alignment_start"]
            a_e = rec[f"{seg}_alignment_end"]
            g_s = rec[f"{seg}_germline_start"]
            g_e = rec[f"{seg}_germline_end"]
            cig = rec[f"{seg}_cigar"]
            if a_s is None:
                continue
            ops = _cigar_op_counts(cig)
            if a_e - a_s != ops["M"] + ops["I"] + ops["D"]:
                failures.append(
                    (i, seg, "alignment-span", a_e - a_s, ops)
                )
            if g_e - g_s != ops["M"] + ops["D"]:
                failures.append(
                    (i, seg, "germline-span", g_e - g_s, ops)
                )
    assert failures == [], f"H.5 invariants broke for: {failures[:3]}"


def test_cigar_only_uses_m_when_pure_recombination():
    # Strict: pure recombine + no mutation + no indel = single M
    # block per segment.
    exp = ga.Experiment.on("human_igh").recombine()
    for rec in exp.run_records(n=20, seed=0):
        for seg in ("v_cigar", "d_cigar", "j_cigar"):
            cig = rec[seg]
            if cig:
                assert _re.fullmatch(r"\d+M", cig), f"{seg}={cig!r}"


def test_alignment_consistent_under_contaminant_pass():
    # Even when the contaminant pass overwrites the entire sequence
    # with non-immune bases, the IR's germline_position tracking is
    # preserved — so sequence_alignment shows the contaminant bases
    # and germline_alignment shows the original (unrelated) V/D/J
    # source. The user can detect this via ``is_contaminant``.
    rec = (
        ga.Experiment.on("human_igh")
        .recombine()
        .contaminate(prob=1.0)
        .run_records(n=1, seed=0)[0]
    )
    assert rec["is_contaminant"] is True
    assert _alignment_invariants_hold(rec) == []
    # And the contaminant divergence is visible: most positions in
    # the V coding region differ between sa and ga.
    sa = rec["sequence_alignment"]
    galn = rec["germline_alignment"]
    v_start, v_end = rec["v_sequence_start"], rec["v_sequence_end"]
    diffs = sum(
        1
        for a, b in zip(sa[v_start:v_end], galn[v_start:v_end])
        if a != b and b not in ("N", "-")
    )
    # Random bases should mismatch in roughly 75% of V positions.
    v_len = v_end - v_start
    assert diffs >= v_len // 2


def test_germline_alignment_byte_matches_source_allele():
    # Verify each V position in germline_alignment equals the
    # corresponding source allele base (post-trim, pre-mutation).
    # The "source allele" is the originally-sampled provenance
    # allele, not whatever set ``rec['v_call']`` projects to today —
    # the call is a live-call view that can become a multi-allele
    # string after a mutation pass refreshes it. We pull provenance
    # from the simulation IR directly.
    result = (
        _human_igh_exp()
        .mutate(count=8)
        .run_records(n=1, seed=42)
    )
    rec = result[0]
    refdata = _human_igh_exp().refdata
    galn = rec["germline_alignment"]

    sim = result.outcomes[0].final_simulation()
    v_id = sim.v_allele_id()
    assert v_id is not None, "expected V provenance after recombine + mutate"
    v_seq = bytes(refdata.v_allele(v_id).seq()).decode().upper()
    v_start = rec["v_sequence_start"]
    v_end = rec["v_sequence_end"]
    # No 5' V trim in the DSL → V starts at germline pos 0.
    # No indels in this test → 1:1 mapping.
    expected = v_seq[: v_end - v_start]
    actual = galn[v_start:v_end]
    assert actual == expected, f"V germline mismatch:\nexpected: {expected[:60]}\nactual:   {actual[:60]}"
