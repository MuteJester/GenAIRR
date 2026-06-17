"""Cohorts: Experiment.run_cohort + CohortResult."""
import pytest

import GenAIRR as ga
import GenAIRR.data as gdata
from GenAIRR.genotype import Genotype
from GenAIRR.cohort import (
    CohortResult,
    CohortSubjectResult,
    _resolve_counts,
    _resolve_subject_ids,
)


def _cfg():
    return gdata.HUMAN_IGH_OGRDB


def test_cohort_result_core_accessors():
    cfg = _cfg()
    g0 = Genotype.sample(cfg, seed=0, subject_id="A")
    g1 = Genotype.sample(cfg, seed=1, subject_id="B")
    r0 = ga.Experiment.on(cfg).with_genotype(g0).recombine().run_records(n=3, seed=1)
    r1 = ga.Experiment.on(cfg).with_genotype(g1).recombine().run_records(n=2, seed=2)
    subjects = [
        CohortSubjectResult(subject_id="A", genotype=g0, result=r0, refdata=object(), seed=11, count=3),
        CohortSubjectResult(subject_id="B", genotype=g1, result=r1, refdata=object(), seed=22, count=2),
    ]
    c = CohortResult(subjects)
    assert c.subject_ids == ["A", "B"]
    assert c.genotypes == [g0, g1]
    assert c.results == [r0, r1]
    assert c.result_for("B") is r1
    assert c.refdata_for("A") is subjects[0].refdata
    assert len(c) == 5                      # total records
    assert len(c.records) == 5              # fresh concatenated list
    c.records.append({"x": 1})              # mutating the returned list...
    assert len(c.records) == 5              # ...does not affect the cohort
    with pytest.raises(KeyError):
        c.result_for("NOPE")


def test_cohort_result_export_union_and_fasta(tmp_path):
    cfg = _cfg()
    g0 = Genotype.sample(cfg, seed=0, subject_id="A")
    g1 = Genotype.sample(cfg, seed=1, subject_id="B")
    # subject A WITH provenance columns, subject B WITHOUT -> schemas differ
    rA = ga.Experiment.on(cfg).with_genotype(g0).recombine().run_records(
        n=2, seed=1, expose_provenance=True)
    rB = ga.Experiment.on(cfg).with_genotype(g1).recombine().run_records(n=2, seed=2)
    # namespace ids the way run_cohort will (so headers are unique)
    for rec in rA.records:
        rec["sequence_id"] = "A_" + rec["sequence_id"]
    for rec in rB.records:
        rec["sequence_id"] = "B_" + rec["sequence_id"]
    c = CohortResult([
        CohortSubjectResult("A", g0, rA, object(), 1, 2),
        CohortSubjectResult("B", g1, rB, object(), 2, 2),
    ])
    df = c.to_dataframe()
    assert len(df) == 4
    # union of columns: truth_v_call exists (from A) even though B lacks it
    assert "truth_v_call" in df.columns

    fasta = tmp_path / "cohort.fasta"
    c.to_fasta(str(fasta))
    headers = [ln[1:].strip() for ln in fasta.read_text().splitlines() if ln.startswith(">")]
    assert len(headers) == 4
    assert len(set(headers)) == 4                  # globally unique
    assert all(h.startswith("A_") or h.startswith("B_") for h in headers)

    csv = tmp_path / "cohort.csv"
    c.to_csv(str(csv))
    assert csv.read_text().count("\n") >= 5        # header + 4 rows


def test_resolve_subject_ids():
    assert _resolve_subject_ids([None, None, None]) == ["subject_0", "subject_1", "subject_2"]
    assert _resolve_subject_ids(["A", "B"]) == ["A", "B"]
    assert _resolve_subject_ids([1, 2]) == ["1", "2"]          # normalized to str
    with pytest.raises(ValueError, match="duplicate"):
        _resolve_subject_ids(["A", "A"])
    with pytest.raises(ValueError, match="duplicate"):
        _resolve_subject_ids([1, "1"])                          # collide after str()
    with pytest.raises(ValueError, match="some .* subject_id"):
        _resolve_subject_ids(["A", None])                       # mixed


def test_resolve_counts():
    assert _resolve_counts(3, 5, None) == [5, 5, 5]
    assert _resolve_counts(3, 5, [1, 2, 0]) == [1, 2, 0]
    with pytest.raises(ValueError, match="length"):
        _resolve_counts(3, 5, [1, 2])
    for bad in (-1, True, "2", 1.5):
        with pytest.raises(ValueError, match="n_per_subject"):
            _resolve_counts(2, bad, None)
    for bad_list in ([1, -1], [1, True], [1, "2"]):
        with pytest.raises(ValueError, match="counts"):
            _resolve_counts(2, 1, bad_list)


def test_run_cohort_happy_path():
    cfg = _cfg()
    gs = [Genotype.sample(cfg, seed=s, subject_id=f"D{s}") for s in range(3)]
    c = (ga.Experiment.on(cfg).recombine()
         .run_cohort(gs, n_per_subject=5, seed=0, expose_provenance=True))
    assert isinstance(c, CohortResult)
    assert c.subject_ids == ["D0", "D1", "D2"]
    assert len(c.genotypes) == 3
    assert len(c) == 15
    # every record is subject-tagged and its sequence_id is namespaced + unique
    sids = [r["sequence_id"] for r in c.records]
    assert len(set(sids)) == 15
    for r in c.records:
        assert r["sequence_id"].startswith(r["subject_id"] + "_")
    # each subject's result exposes its genotype
    assert c.result_for("D1").genotypes[0].subject_id == "D1"


def test_run_cohort_deterministic_and_independent():
    cfg = _cfg()
    gs = [Genotype.sample(cfg, seed=s, subject_id=f"D{s}") for s in range(2)]
    a = ga.Experiment.on(cfg).recombine().run_cohort(gs, n_per_subject=4, seed=7)
    b = ga.Experiment.on(cfg).recombine().run_cohort(gs, n_per_subject=4, seed=7)
    assert [r["sequence"] for r in a.records] == [r["sequence"] for r in b.records]
    d = ga.Experiment.on(cfg).recombine().run_cohort(gs, n_per_subject=4, seed=8)
    assert [r["sequence"] for r in a.records] != [r["sequence"] for r in d.records]


def test_run_cohort_does_not_mutate_base_experiment():
    cfg = _cfg()
    gs = [Genotype.sample(cfg, seed=s, subject_id=f"D{s}") for s in range(2)]
    exp = ga.Experiment.on(cfg).recombine()
    exp.run_cohort(gs, n_per_subject=2, seed=0)
    assert exp._genotype is None
    # still usable as a plain (no-genotype) experiment afterwards
    res = exp.run_records(n=2, seed=1)
    assert len(res) == 2


def test_run_cohort_auto_subject_ids_and_no_mutation_of_originals():
    cfg = _cfg()
    gs = [Genotype.sample(cfg, seed=s) for s in range(3)]   # no subject_id set
    c = ga.Experiment.on(cfg).recombine().run_cohort(gs, n_per_subject=1, seed=0)
    assert c.subject_ids == ["subject_0", "subject_1", "subject_2"]
    # originals are untouched (resolution happened on snapshots)
    assert all(g.subject_id is None for g in gs)


def test_run_cohort_counts_override_and_zero():
    cfg = _cfg()
    gs = [Genotype.sample(cfg, seed=s, subject_id=f"D{s}") for s in range(3)]
    c = ga.Experiment.on(cfg).recombine().run_cohort(gs, seed=0, counts=[4, 0, 2])
    assert [len(r) for r in c.results] == [4, 0, 2]
    assert len(c) == 6
    # zero-count subject: empty records, but present with stamped genotype + refdata
    zero = c.result_for("D1")
    assert zero.records == []
    assert zero.genotypes[0].subject_id == "D1"
    assert c.refdata_for("D1") is not None


def test_run_cohort_counts_length_mismatch_raises():
    cfg = _cfg()
    gs = [Genotype.sample(cfg, seed=s, subject_id=f"D{s}") for s in range(2)]
    with pytest.raises(ValueError, match="length"):
        ga.Experiment.on(cfg).recombine().run_cohort(gs, counts=[1])


def test_run_cohort_duplicate_subject_ids_raise():
    cfg = _cfg()
    gs = [Genotype.sample(cfg, seed=s, subject_id="DUP") for s in range(2)]
    with pytest.raises(ValueError, match="duplicate"):
        ga.Experiment.on(cfg).recombine().run_cohort(gs, n_per_subject=1)


def test_run_cohort_mutual_exclusions():
    cfg = _cfg()
    g = Genotype.sample(cfg, seed=0, subject_id="A")
    with pytest.raises(ValueError, match="non-empty"):
        ga.Experiment.on(cfg).recombine().run_cohort([])
    with pytest.raises(ValueError, match="with_genotype"):
        ga.Experiment.on(cfg).with_genotype(g).recombine().run_cohort([g])
    with pytest.raises(ValueError, match="restrict_alleles"):
        (ga.Experiment.on(cfg).recombine()
         .restrict_alleles(v=cfg.v_alleles[next(iter(cfg.v_alleles))][0].name)
         .run_cohort([g]))


def test_run_cohort_cartridge_hash_mismatch_raises():
    cfg = _cfg()
    g = Genotype.sample(cfg, seed=0, subject_id="A")
    g._source_hash = "sha256:deadbeef"          # forge a mismatch
    with pytest.raises(ValueError, match="different cartridge"):
        ga.Experiment.on(cfg).recombine().run_cohort([g])


def test_run_cohort_failure_midway_leaves_base_experiment_clean():
    cfg = _cfg()
    good = Genotype.sample(cfg, seed=0, subject_id="A")
    bad = Genotype.sample(cfg, seed=1, subject_id="B")
    import GenAIRR.genotype as _gmod
    exp = ga.Experiment.on(cfg).recombine()
    orig_snapshot = _gmod.Genotype._snapshot

    calls = {"n": 0}

    def boom(self):
        calls["n"] += 1
        if calls["n"] == 2:                      # fail on the 2nd subject
            raise RuntimeError("forced compile-time failure")
        return orig_snapshot(self)

    _gmod.Genotype._snapshot = boom
    try:
        with pytest.raises(RuntimeError, match="forced"):
            exp.run_cohort([good, bad], n_per_subject=1, seed=0)
    finally:
        _gmod.Genotype._snapshot = orig_snapshot
    # base experiment is untouched and still runnable
    assert exp._genotype is None
    assert len(exp.run_records(n=2, seed=1)) == 2
