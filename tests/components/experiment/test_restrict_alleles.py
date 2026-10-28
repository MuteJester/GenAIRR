"""Locking the alleles a run may draw, via restrict_alleles().
"""
from __future__ import annotations
import pytest
from GenAIRR import CompiledExperiment, Experiment
import GenAIRR as ga  # noqa: E402

def _v_name(refdata, outcome):
    return refdata.v_allele(outcome.final_simulation().v_allele_id()).name


def _d_name(refdata, outcome):
    return refdata.d_allele(outcome.final_simulation().d_allele_id()).name


def _j_name(refdata, outcome):
    return refdata.j_allele(outcome.final_simulation().j_allele_id()).name


def test_using_returns_self_for_chaining():
    exp = ga.Experiment.on("human_igh")
    same = exp.restrict_alleles(v="IGHVF1-G1*01")
    assert same is exp


def test_using_locks_single_v_allele():
    exp = (
        ga.Experiment.on("human_igh")
        .restrict_alleles(v="IGHVF1-G1*01")
        .recombine()
    )
    for o in exp.run(n=10, seed=0):
        assert _v_name(exp.refdata, o) == "IGHVF1-G1*01"


def test_using_locks_single_d_allele():
    exp = (
        ga.Experiment.on("human_igh")
        .restrict_alleles(d="IGHD3-3*01")
        .recombine()
    )
    for o in exp.run(n=10, seed=0):
        assert _d_name(exp.refdata, o) == "IGHD3-3*01"


def test_using_locks_single_j_allele():
    exp = (
        ga.Experiment.on("human_igh")
        .restrict_alleles(j="IGHJ4*02")
        .recombine()
    )
    for o in exp.run(n=10, seed=0):
        assert _j_name(exp.refdata, o) == "IGHJ4*02"


def test_using_locks_full_vdj_simultaneously():
    exp = (
        ga.Experiment.on("human_igh")
        .restrict_alleles(v="IGHVF1-G1*01", d="IGHD3-3*01", j="IGHJ4*02")
        .recombine()
    )
    for o in exp.run(n=10, seed=0):
        assert _v_name(exp.refdata, o) == "IGHVF1-G1*01"
        assert _d_name(exp.refdata, o) == "IGHD3-3*01"
        assert _j_name(exp.refdata, o) == "IGHJ4*02"


def test_using_accepts_multi_allele_list():
    allowed = ["IGHVF10-G38*02", "IGHVF10-G38*04"]
    exp = (
        ga.Experiment.on("human_igh")
        .restrict_alleles(v=allowed)
        .recombine()
    )
    seen = {_v_name(exp.refdata, o) for o in exp.run(n=40, seed=0)}
    # Sample stays inside the allowed subset.
    assert seen.issubset(set(allowed))
    # Both members get hit at least once across 40 runs.
    assert seen == set(allowed)


def test_using_accepts_tuple_input():
    allowed = ("IGHVF10-G38*02", "IGHVF10-G38*04")
    exp = (
        ga.Experiment.on("human_igh")
        .restrict_alleles(v=allowed)
        .recombine()
    )
    for o in exp.run(n=10, seed=0):
        assert _v_name(exp.refdata, o) in allowed


def test_using_after_recombine_still_applies():
    # Locks are injected at compile time, so order of fluent calls
    # doesn't matter.
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .restrict_alleles(v="IGHVF1-G1*01")
    )
    for o in exp.run(n=5, seed=0):
        assert _v_name(exp.refdata, o) == "IGHVF1-G1*01"


def test_using_called_twice_overlays_per_segment():
    # First .restrict_alleles() locks V, second locks D — both should hold.
    exp = (
        ga.Experiment.on("human_igh")
        .restrict_alleles(v="IGHVF1-G1*01")
        .restrict_alleles(d="IGHD3-3*01")
        .recombine()
    )
    for o in exp.run(n=5, seed=0):
        assert _v_name(exp.refdata, o) == "IGHVF1-G1*01"
        assert _d_name(exp.refdata, o) == "IGHD3-3*01"


def test_using_called_twice_overwrites_same_segment():
    # Second call's V lock replaces the first one.
    exp = (
        ga.Experiment.on("human_igh")
        .restrict_alleles(v="IGHVF1-G1*01")
        .restrict_alleles(v="IGHVF10-G38*02")
        .recombine()
    )
    for o in exp.run(n=5, seed=0):
        assert _v_name(exp.refdata, o) == "IGHVF10-G38*02"


def test_using_none_clears_prior_lock():
    # Clearing reverts to the full pool.
    exp = (
        ga.Experiment.on("human_igh")
        .restrict_alleles(v="IGHVF1-G1*01")
        .restrict_alleles(v=None)
        .recombine()
    )
    seen = {_v_name(exp.refdata, o) for o in exp.run(n=30, seed=0)}
    assert len(seen) > 1


def test_using_omitted_kwarg_leaves_prior_lock_untouched():
    # Calling .restrict_alleles(d=...) shouldn't clobber an earlier V lock.
    exp = (
        ga.Experiment.on("human_igh")
        .restrict_alleles(v="IGHVF1-G1*01")
        .restrict_alleles(d="IGHD3-3*01")  # only D supplied
        .recombine()
    )
    for o in exp.run(n=5, seed=0):
        assert _v_name(exp.refdata, o) == "IGHVF1-G1*01"
        assert _d_name(exp.refdata, o) == "IGHD3-3*01"


def test_using_unknown_v_name_raises():
    with pytest.raises(ValueError, match="no V allele named"):
        ga.Experiment.on("human_igh").restrict_alleles(v="NOPE*01")


def test_using_unknown_d_name_raises():
    with pytest.raises(ValueError, match="no D allele named"):
        ga.Experiment.on("human_igh").restrict_alleles(d="IGHD-NOPE*01")


def test_using_unknown_j_name_raises():
    with pytest.raises(ValueError, match="no J allele named"):
        ga.Experiment.on("human_igh").restrict_alleles(j="IGHJ-NOPE*01")


def test_using_d_lock_on_vj_chain_raises():
    with pytest.raises(ValueError, match="cannot lock D alleles"):
        ga.Experiment.on("human_igk").restrict_alleles(d="IGHD3-3*01")


def test_using_empty_list_raises():
    with pytest.raises(ValueError, match="must be non-empty"):
        ga.Experiment.on("human_igh").restrict_alleles(v=[])


def test_using_duplicate_in_list_raises():
    with pytest.raises(ValueError, match="duplicate"):
        ga.Experiment.on("human_igh").restrict_alleles(
            v=["IGHVF1-G1*01", "IGHVF1-G1*01"]
        )


def test_using_non_string_entry_raises():
    with pytest.raises(TypeError, match="must all be strings"):
        ga.Experiment.on("human_igh").restrict_alleles(v=["IGHVF1-G1*01", 5])  # type: ignore[list-item]


def test_using_does_not_increment_step_count():
    # ``.restrict_alleles()`` is a configuration knob, not a pipeline step.
    exp = ga.Experiment.on("human_igh").restrict_alleles(v="IGHVF1-G1*01")
    assert exp.step_count == 0


def test_using_is_deterministic_under_same_seed():
    locked = ["IGHVF10-G38*02", "IGHVF10-G38*04"]
    a = (
        ga.Experiment.on("human_igh")
        .restrict_alleles(v=locked)
        .recombine()
        .run(n=10, seed=0xCAFE)
    )
    b = (
        ga.Experiment.on("human_igh")
        .restrict_alleles(v=locked)
        .recombine()
        .run(n=10, seed=0xCAFE)
    )
    for x, y in zip(a, b):
        assert x.final_simulation().bases() == y.final_simulation().bases()


def test_using_composes_with_mutation_and_corruption():
    # Locks should work alongside the rest of the DSL pipeline.
    exp = (
        ga.Experiment.on("human_igh")
        .restrict_alleles(v="IGHVF1-G1*01")
        .recombine()
        .mutate(count=8)
        .pcr_amplify(count=2)
    )
    out = exp.run(n=3, seed=42)
    for o in out:
        assert _v_name(exp.refdata, o) == "IGHVF1-G1*01"
        assert o.trace().find("mutate.s5f.count").value == 8
        assert o.trace().find("corrupt.pcr.count").value == 2


def test_using_works_for_vj_chain_v_and_j():
    # Pick the first V / J allele in the IGK pool so the test stays
    # robust to refdata reshuffles.
    refdata = ga.Experiment.on("human_igk").refdata
    v_name = refdata.v_allele(0).name
    j_name = refdata.j_allele(0).name
    exp = (
        ga.Experiment.on("human_igk")
        .restrict_alleles(v=v_name, j=j_name)
        .recombine()
    )
    for o in exp.run(n=5, seed=0):
        assert _v_name(exp.refdata, o) == v_name
        assert _j_name(exp.refdata, o) == j_name


def test_using_with_productive_contract():
    # Locks + productive() should still produce productive sequences.
    exp = (
        ga.Experiment.on("human_igh")
        .restrict_alleles(v="IGHVF1-G1*01", j="IGHJ4*02")
        .recombine()
    )
    out = exp.productive_only().run_records(n=3, seed=0)
    for rec in out:
        assert rec["v_call"] == "IGHVF1-G1*01"
        assert rec["j_call"] == "IGHJ4*02"
        assert rec["productive"] is True


def test_using_compile_is_idempotent_with_locks():
    exp = (
        ga.Experiment.on("human_igh")
        .restrict_alleles(v="IGHVF1-G1*01")
        .recombine()
    )
    plan_a = exp.compile(allow_curatable_refdata=True)
    plan_b = exp.compile(allow_curatable_refdata=True)
    a = plan_a.run(n=3, seed=0)
    b = plan_b.run(n=3, seed=0)
    for x, y in zip(a, b):
        assert x.final_simulation().bases() == y.final_simulation().bases()
