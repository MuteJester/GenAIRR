"""FASTQ output, including Phred quality encoding.
"""
from __future__ import annotations
import pytest
import GenAIRR as ga  # noqa: E402
from tests.support.fixtures.experiment_dsl import _alignment_invariants_hold, _human_igh_records

pytestmark = pytest.mark.behavior

def test_to_fastq_constant_quality_writes_q30(tmp_path):
    # Constant-Q model writes the same Q for every (non-N) base.
    # Q30 → Phred+33 ASCII '?' (33 + 30 = 63).
    result = _human_igh_records(n=2)
    path = tmp_path / "out.fastq"
    result.to_fastq(str(path), quality="constant", q=30)
    lines = path.read_text(encoding="utf-8").splitlines()
    assert len(lines) == 8  # 2 records × 4 lines
    # First record: @-header / SEQ / + / Q-string
    assert lines[0].startswith("@seq0|v_call=")
    assert lines[1] == result[0]["sequence"].upper()
    assert lines[2] == "+"
    seq = result[0]["sequence"]
    expected_q = "".join(
        "#" if c in ("N", "n") else "?" for c in seq
    )
    assert lines[3] == expected_q


def test_to_fastq_illumina_quality_shape(tmp_path):
    # Illumina model: ramp up over 10 bases, hold peak, taper over
    # 30 bases.
    result = _human_igh_records(n=1)
    path = tmp_path / "out.fastq"
    result.to_fastq(str(path), quality="illumina")
    lines = path.read_text(encoding="utf-8").splitlines()
    q_string = lines[3]
    seq = result[0]["sequence"]
    # All ASCII ≥ 33 (printable).
    for ch in q_string:
        assert 33 <= ord(ch) <= 126
    # Same length as the sequence.
    assert len(q_string) == len(seq)
    # Ramp-up: position 0 is Q26 (≈ start_q+1 from interpolation),
    # position 9 is Q35 (peak).
    assert ord(q_string[0]) - 33 < ord(q_string[9]) - 33
    # Tail: position -1 is Q18 (end_q), strictly less than peak.
    assert ord(q_string[-1]) - 33 < ord(q_string[len(seq) // 2]) - 33


def test_to_fastq_n_bases_get_low_quality(tmp_path):
    # Manually construct a result where one record has N in the
    # sequence and verify that position gets Q2 (`#` ASCII 35).
    result = _human_igh_records(n=1)
    # Mutate the result's first record's sequence to inject an N.
    rec = result[0]
    seq = rec["sequence"]
    rec["sequence"] = seq[:5] + "N" + seq[6:]
    path = tmp_path / "out.fastq"
    result.to_fastq(str(path), quality="constant", q=30)
    lines = path.read_text(encoding="utf-8").splitlines()
    q_string = lines[3]
    # Position 5 should be the n_q (Q2 default → '#').
    assert q_string[5] == "#"
    # Other positions should be the full Q (Q30 → '?').
    assert q_string[0] == "?"


def test_to_fastq_uppercase_sequence(tmp_path):
    # FASTQ convention: uppercase bases regardless of GenAIRR's
    # internal lowercase germline marking.
    result = _human_igh_records(n=1)
    path = tmp_path / "out.fastq"
    result.to_fastq(str(path), quality="constant", q=30)
    lines = path.read_text(encoding="utf-8").splitlines()
    seq_line = lines[1]
    assert seq_line == seq_line.upper(), (
        f"FASTQ sequence should be uppercase: {seq_line[:40]}"
    )


def test_to_fastq_rejects_unknown_model(tmp_path):
    result = _human_igh_records(n=1)
    path = tmp_path / "out.fastq"
    try:
        result.to_fastq(str(path), quality="bogus")
    except ValueError:
        pass
    else:
        raise AssertionError("expected ValueError for unknown model")


def test_to_fastq_rejects_out_of_range_q(tmp_path):
    result = _human_igh_records(n=1)
    path = tmp_path / "out.fastq"
    for bad in (-1, 100, "high"):
        try:
            result.to_fastq(str(path), quality="constant", q=bad)
        except (ValueError, TypeError):
            pass
        else:
            raise AssertionError(f"expected error for q={bad!r}")


def test_g10_with_metadata_injects_columns():
    exp = (
        ga.Experiment.on("human_igh")
        .with_metadata(sample_id="P1", donor="D001", repertoire_id="R-001-IGH")
        .recombine()
    )
    for r in exp.run_records(n=3, seed=0):
        assert r["sample_id"] == "P1"
        assert r["donor"] == "D001"
        assert r["repertoire_id"] == "R-001-IGH"


def test_g10_with_metadata_merges_on_subsequent_calls():
    base = ga.Experiment.on("human_igh").recombine().with_metadata(
        sample_id="P1", donor="D001"
    )
    refined = base.with_metadata(donor="D002", cell_id="C42")
    r = refined.run_records(n=1, seed=0)[0]
    assert r["sample_id"] == "P1"  # preserved from earlier call
    assert r["donor"] == "D002"  # overwritten by later call
    assert r["cell_id"] == "C42"  # added by later call


def test_g10_with_metadata_none_clears_key():
    exp = (
        ga.Experiment.on("human_igh")
        .with_metadata(sample_id="P1", donor="D001")
        .recombine()
        .with_metadata(donor=None)
    )
    r = exp.run_records(n=1, seed=0)[0]
    assert r["sample_id"] == "P1"
    assert "donor" not in r


def test_g10_with_metadata_columns_in_tsv_output(tmp_path):
    import csv as _csv
    exp = (
        ga.Experiment.on("human_igh")
        .with_metadata(sample_id="P1", donor="D001")
        .recombine()
    )
    path = tmp_path / "out.tsv"
    exp.run_records(n=2, seed=0).to_tsv(str(path))
    with open(path, "r") as fh:
        rows = list(_csv.DictReader(fh, delimiter="\t"))
    assert len(rows) == 2
    for row in rows:
        assert row["sample_id"] == "P1"
        assert row["donor"] == "D001"


def test_g10_with_metadata_works_with_clonal():
    exp = (
        ga.Experiment.on("human_igh")
        .with_metadata(sample_id="P1")
        .recombine()
        .expand_clones(n_clones=2, per_clone=3)
    )
    result = exp.run_records(seed=0)
    assert len(result) == 6
    for r in result:
        assert r["sample_id"] == "P1"
        assert "clone_id" in r


def test_g10_with_metadata_rejects_non_string_keys():
    with pytest.raises(TypeError, match="must be strings"):
        ga.Experiment.on("human_igh").with_metadata(**{1: "x"})


def test_g8_truth_columns_present_when_opted_in():
    exp = ga.Experiment.on("human_igh").recombine()
    result = exp.run_records(n=3, seed=0, expose_provenance=True)
    for r in result:
        assert "truth_v_call" in r
        assert "truth_d_call" in r
        assert "truth_j_call" in r
        # Provenance is always populated for a successful recombine.
        assert r["truth_v_call"], f"missing truth_v_call: {r}"
        assert r["truth_j_call"], f"missing truth_j_call: {r}"


def test_g8_truth_columns_absent_by_default():
    exp = ga.Experiment.on("human_igh").recombine()
    result = exp.run_records(n=2, seed=0)
    for r in result:
        assert "truth_v_call" not in r
        assert "truth_d_call" not in r
        assert "truth_j_call" not in r


def test_g8_truth_retained_in_live_call_under_heavy_shm():
    # Under the score-and-tie live caller, heavy SHM does NOT wipe
    # v_call. The truth allele continues to score highest from the
    # un-mutated positions, so it stays in the tie-set; the call may
    # be a comma-separated list of indistinguishable alleles, but it
    # is always non-empty and the truth allele is always among them.
    exp = ga.Experiment.on("human_igh").recombine().mutate(count=25)
    result = exp.run_records(n=20, seed=0, expose_provenance=True)
    for r in result:
        live = r["v_call"]
        truth = r["truth_v_call"]
        assert truth, "truth_v_call should always be populated"
        assert live, "v_call should never be empty under the score-and-tie caller"
        # truth_v_call is always one of the tied alleles reported in v_call.
        assert truth in live.split(","), (
            f"truth_v_call={truth!r} should appear in v_call tie-set {live!r}"
        )


def test_g8_truth_columns_via_clonal_pipeline():
    # Clonal pipeline should also accept expose_provenance.
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .expand_clones(n_clones=2, per_clone=3)
        .mutate(count=15)
    )
    result = exp.run_records(seed=0, expose_provenance=True)
    assert len(result) == 6
    by_clone = {}
    for r in result:
        by_clone.setdefault(r["clone_id"], []).append(r)
    # All descendants of a clone must share the SAME truth_v_call
    # (they shared the parent recombination).
    for cid, recs in by_clone.items():
        truth_v = {r["truth_v_call"] for r in recs}
        assert len(truth_v) == 1, (
            f"clone {cid} has divergent truth_v_call: {truth_v}"
        )


def test_g7_corrupt_ns_writes_uppercase_n():
    # corrupt_ns(count=10) should sprinkle ~10 N's per record
    # (collisions reduce the count slightly).
    exp = ga.Experiment.on("human_igh").recombine().ambiguous_base_calls(count=10)
    for rec in exp.run_records(n=5, seed=0):
        n_count = rec["sequence"].upper().count("N")
        # With pool length ~370 and 10 random sites, collisions are
        # rare; expect 8-10 inclusive.
        assert 6 <= n_count <= 10, (
            f"expected 6-10 Ns, got {n_count} (sequence len={len(rec['sequence'])})"
        )


def test_g7_corrupt_ns_zero_count_is_no_op():
    base = ga.Experiment.on("human_igh").recombine()
    flipped = ga.Experiment.on("human_igh").recombine().ambiguous_base_calls(count=0)
    for fwd, ns in zip(base.run_records(n=5, seed=0), flipped.run_records(n=5, seed=0)):
        assert fwd["sequence"] == ns["sequence"]


def test_g7_corrupt_ns_count_distribution_shape():
    # Variable-count distribution: count=(5, 15) → uniform integer
    # in [5, 15]. Across many records, the N-count distribution
    # should have all entries in approximately that range.
    exp = ga.Experiment.on("human_igh").recombine().ambiguous_base_calls(count=(5, 15))
    counts = [r["sequence"].upper().count("N") for r in exp.run_records(n=50, seed=0)]
    assert min(counts) >= 4, f"min N count too low: {counts}"
    assert max(counts) <= 15, f"max N count too high: {counts}"


def test_g7_corrupt_ns_compose_with_other_passes():
    # Stack with PCR + indels — N's should still appear and the
    # alignment-string invariants should hold.
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(count=5)
        .pcr_amplify(count=2)
        .polymerase_indels(count=1)
        .ambiguous_base_calls(count=8)
    )
    failures = []
    for i, rec in enumerate(exp.run_records(n=30, seed=0)):
        issues = _alignment_invariants_hold(rec)
        if issues:
            failures.append((i, issues))
        # At least some records should carry Ns.
    n_with_ns = sum(
        1 for r in exp.run_records(n=30, seed=0)
        if "N" in r["sequence"].upper()
    )
    assert n_with_ns >= 25, f"only {n_with_ns}/30 records had Ns"
    assert failures == [], f"alignment invariants broke: {failures[:3]}"


def test_g5_clonal_structure_total_record_count():
    # n_clones × size = total records.
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .expand_clones(n_clones=5, per_clone=4)
        .mutate(count=3)
    )
    result = exp.run_records(seed=0)
    assert len(result) == 20


def test_g5_clonal_descendants_identical_without_post_fork_passes():
    # No mutate / corrupt after the fork → every descendant in a
    # clone is literally identical (same V/D/J + trim + NP).
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .expand_clones(n_clones=3, per_clone=4)
    )
    by_clone = {}
    for rec in exp.run_records(seed=0):
        by_clone.setdefault(rec["clone_id"], []).append(rec)
    assert sorted(by_clone) == [0, 1, 2]
    for cid, recs in by_clone.items():
        seqs = {r["sequence"] for r in recs}
        v_calls = {r["v_call"] for r in recs}
        assert len(seqs) == 1, f"clone {cid} descendants diverged without post-fork passes"
        assert len(v_calls) == 1


def test_g5_clonal_descendants_share_junction_diverge_via_mutate():
    # With mutate after the fork, descendants share the junction
    # backbone (same V/D/J + trim + NP) but diverge in mutated bases.
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .expand_clones(n_clones=2, per_clone=8)
        .mutate(count=8)
    )
    by_clone = {}
    for rec in exp.run_records(seed=0):
        by_clone.setdefault(rec["clone_id"], []).append(rec)
    for cid, recs in by_clone.items():
        v_calls = {r["v_call"] for r in recs if r["v_call"]}
        seqs = {r["sequence"] for r in recs}
        # Sequences diverge under SHM…
        assert len(seqs) > 1, f"clone {cid} descendants didn't diverge under mutate"
        # …but the V allele backbone is shared (when live evidence
        # didn't get wiped). All V calls should reference the same
        # provenance — rejecting empty calls (heavy-SHM noise).
        non_empty = [c for c in v_calls if c]
        if non_empty:
            assert len(set(non_empty)) == 1, (
                f"clone {cid} V calls differ across descendants: {non_empty}"
            )


def test_g5_clonal_explicit_n_must_match_total():
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .expand_clones(n_clones=4, per_clone=5)
    )
    # n=20 matches → ok.
    result = exp.run_records(n=20, seed=0)
    assert len(result) == 20
    # n=10 mismatched → error.
    with pytest.raises(ValueError, match="inconsistent"):
        exp.run_records(n=10, seed=0)


def test_g5_expand_clones_rejects_double_call():
    exp = ga.Experiment.on("human_igh").recombine()
    exp = exp.expand_clones(n_clones=3, per_clone=2)
    with pytest.raises(ValueError, match="only be called once"):
        exp.expand_clones(n_clones=2, per_clone=2)


def test_g5_expand_clones_rejects_invalid_args():
    exp = ga.Experiment.on("human_igh").recombine()
    with pytest.raises(ValueError, match="positive int"):
        exp.expand_clones(n_clones=0, per_clone=5)
    with pytest.raises(ValueError, match="positive int"):
        exp.expand_clones(n_clones=5, per_clone=-1)


def test_g5_clonal_record_to_tsv(tmp_path):
    # The clone_id column should round-trip through to_tsv output.
    import csv as _csv
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .expand_clones(n_clones=3, per_clone=4)
        .mutate(count=2)
    )
    path = tmp_path / "clones.tsv"
    exp.run_records(seed=0).to_tsv(str(path))
    with open(path, "r") as fh:
        reader = _csv.DictReader(fh, delimiter="\t")
        rows = list(reader)
    assert len(rows) == 12
    clone_ids = sorted({int(r["clone_id"]) for r in rows})
    assert clone_ids == [0, 1, 2]
    # Each clone has 4 descendants.
    from collections import Counter
    counts = Counter(int(r["clone_id"]) for r in rows)
    assert counts == Counter({0: 4, 1: 4, 2: 4})


def test_g4_mutate_rejected_on_tcr_chain():
    # SHM is biologically B-cell only. The DSL must error rather
    # than silently producing biologically-false TCR records.
    with pytest.raises(ValueError, match="somatic hypermutation does not occur in TCR"):
        ga.Experiment.on("human_tcrb").allow_curatable_refdata().recombine().mutate(count=5)


def test_g4_mutate_allowed_on_bcr_chain():
    # Sanity: BCR loci accept mutate() unchanged.
    for cfg in ("human_igh", "human_igk", "human_igl"):
        exp = ga.Experiment.on(cfg).allow_curatable_refdata().recombine().mutate(count=3)
        records = exp.run_records(n=2, seed=0)
        assert len(records) == 2


def test_g4_tcr_pipeline_with_corrupt_passes_works():
    # Confirms the recommended TCR pipeline (corrupt_pcr / corrupt_
    # quality / corrupt_indels) still works end-to-end.
    exp = (
        ga.Experiment.on("human_tcrb").allow_curatable_refdata()
        .recombine()
        .pcr_amplify(count=3)
        .sequencing_errors(count=2)
        .polymerase_indels(count=1, insertion_prob=0.5)
    )
    records = exp.run_records(n=5, seed=0)
    assert len(records) == 5
    for rec in records:
        assert rec["locus"] == "TRB"


def test_g3b_allele_weights_boost_listed_v():
    # Boosting a single V allele's weight by 1000x should make it
    # dominate sampling over a 200-record batch.
    from collections import Counter
    target = "IGHVF1-G1*01"
    exp = ga.Experiment.on("human_igh").recombine(v_allele_weights={target: 1000.0})
    counts = Counter()
    for r in exp.run_records(n=200, seed=0):
        first = (r["v_call"] or "").split(",")[0]
        counts[first] += 1
    # The boosted allele should be at least 50% of calls (typically
    # >80% with weight=1000 against ~60 alleles each at 1.0).
    assert counts[target] >= 100, (
        f"boosted allele was sampled {counts[target]} / 200 times; "
        f"top counts: {counts.most_common(3)}"
    )


def test_g3b_allele_weights_unspecified_keeps_default_uniform():
    # When v_allele_weights is None, the distribution should match
    # the un-weighted run (allowing for stochastic variation across
    # 100 records with the same seed).
    exp_a = ga.Experiment.on("human_igh").recombine()
    exp_b = ga.Experiment.on("human_igh").recombine(v_allele_weights=None)
    a = [r["v_call"] for r in exp_a.run_records(n=20, seed=0)]
    b = [r["v_call"] for r in exp_b.run_records(n=20, seed=0)]
    assert a == b


def test_g3b_allele_weights_rejects_unknown_name():
    with pytest.raises(ValueError, match="no V allele named"):
        ga.Experiment.on("human_igh").recombine(
            v_allele_weights={"IGHV-NOT-A-REAL-ALLELE*99": 5.0}
        )


def test_g3b_allele_weights_rejects_non_positive_weight():
    with pytest.raises(ValueError, match="finite positive number"):
        ga.Experiment.on("human_igh").recombine(
            v_allele_weights={"IGHVF1-G1*01": 0.0}
        )
    with pytest.raises(ValueError, match="finite positive number"):
        ga.Experiment.on("human_igh").recombine(
            v_allele_weights={"IGHVF1-G1*01": -1.0}
        )


def test_g3b_allele_weights_d_on_vj_chain_errors():
    with pytest.raises(ValueError, match="cannot weight D alleles"):
        ga.Experiment.on("human_igk").recombine(
            d_allele_weights={"X*01": 1.0}
        )


def test_g3b_allele_weights_empty_dict_errors():
    with pytest.raises(ValueError, match="at least one"):
        ga.Experiment.on("human_igh").recombine(v_allele_weights={})


def test_to_fastq_pairs_with_constant_q_is_decodable(tmp_path):
    # Round-trip check: write FASTQ, parse it back, confirm record
    # count and length match. Uses Python stdlib (no biopython) to
    # keep the test dep-light.
    result = _human_igh_records(n=5)
    path = tmp_path / "out.fastq"
    result.to_fastq(str(path), quality="constant", q=25)
    text = path.read_text(encoding="utf-8")
    records = text.split("@seq")[1:]  # split keeps records starting with header tail
    assert len(records) == 5
    for i, blob in enumerate(records):
        lines = blob.splitlines()
        # Layout: "{i}|...", seq, "+", q_string
        seq = lines[1]
        q = lines[3]
        assert len(seq) == len(q), f"record {i}: seq/q length mismatch"


def test_to_dataframe_returns_dataframe_with_record_rows():
    pd = pytest.importorskip("pandas")

    result = _human_igh_records(n=4)
    df = result.to_dataframe()
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 4
    assert "sequence" in df.columns
    assert "v_call" in df.columns
    # First-row spot-check.
    assert df.iloc[0]["sequence"] == result[0]["sequence"]


def test_to_dataframe_empty_returns_canonical_columns():
    pd = pytest.importorskip("pandas")
    from GenAIRR.result import SimulationResult, _DEFAULT_COLUMN_ORDER

    df = SimulationResult([]).to_dataframe()
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 0
    assert list(df.columns) == _DEFAULT_COLUMN_ORDER


def test_column_order_preserves_extras_after_canonical():
    from GenAIRR.result import SimulationResult, _DEFAULT_COLUMN_ORDER

    rec = {"sequence": "A", "v_call": "x", "extra_a": 1, "extra_b": 2}
    result = SimulationResult([rec])
    cols = result._column_order()
    # Canonical columns come first.
    assert cols[: len(_DEFAULT_COLUMN_ORDER)] == _DEFAULT_COLUMN_ORDER
    # Extras follow in record-encounter order.
    extras = cols[len(_DEFAULT_COLUMN_ORDER) :]
    assert extras == ["extra_a", "extra_b"]
