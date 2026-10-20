"""Capability boundaries must reject unsupported duplex paths before effects."""
import pytest

from GenAIRR import _engine as ge
from GenAIRR.experimental import (
    CoupledDuplexRepair,
    DuplexJoin,
    project_duplex_simulation,
)


def fixture():
    refdata = ge.RefDataConfig.vdj()
    refdata.add_v_allele("V*01", "V", b"AAATGT", anchor=3)
    refdata.add_d_allele("D*01", "D", b"CCC")
    refdata.add_j_allele("J*01", "J", b"TTTAAA", anchor=0)
    script = CoupledDuplexRepair("AAATGT", "CCC", "TTTAAA")
    for name, start in (("vd", 6), ("dj", 3)):
        script = script.join(
            DuplexJoin(name).open("left", offset=0).open("right", offset=0)
            .synapse(right_start=start).ligate()
        )
    seed = script.seed(
        refdata, v_call="V*01", d_call="D*01", j_call="J*01",
        reference_ends_are_coding_ends=True,
    )
    return refdata, seed


@pytest.mark.parametrize("strict", [False, True])
def test_duplex_unknown_p_uses_nullable_record_dispatch_and_refuses_external_seed_trace(strict):
    refdata, seed = fixture()
    compiled = ge.PassPlan().compile(refdata=refdata, strict=strict)
    outcome = compiled.run_from(seed, 7)
    projected = project_duplex_simulation(outcome.final_simulation(), refdata)
    assert projected["record"]["sequence"] == "AAATGTCCCTTTAAA"
    assert projected["record"]["p_v_3_length"] is None
    assert projected["record"]["p_d_5_length"] is None
    record = outcome.airr_record(refdata)
    assert record["p_v_3_length"] is None
    assert record["duplex_repair"]["current_paths"] == projected["paths"]
    assert outcome.validate_record(refdata) == []
    with pytest.raises(ValueError, match="initial repair snapshot"):
        compiled.trace_file_from(outcome, 7)


@pytest.mark.parametrize("strict", [False, True])
def test_duplex_contamination_is_unsupported_even_when_configured_never_to_fire(strict):
    refdata, seed = fixture()
    plan = ge.PassPlan()
    plan.push_corrupt_rev_comp(1.0)
    plan.push_corrupt_contaminant(0.0)
    compiled = plan.compile(refdata=refdata, strict=strict)
    with pytest.raises(ge.StrictSamplingError, match="unsupported for a duplex seed"):
        compiled.run_from(seed, 9)
    assert project_duplex_simulation(seed, refdata)["record"]["rev_comp"] is False


@pytest.mark.parametrize("family_outcome", [False, True])
def test_duplex_family_entrypoints_fail_before_mutation(family_outcome):
    refdata, seed = fixture()
    args = dict(
        mutability=[1.0] * 1024, substitution=[0.25] * 4096, rate=0.0,
        lambda_base=1.0, lambda_mut=0.0, max_generations=1,
        n_max=2, n_sample=1, seed=8,
    )
    if family_outcome:
        compiled = ge.PassPlan().compile(refdata=refdata)
        args["founder"] = compiled.run_from(seed, 7)
        args["refdata"] = refdata
        method = ge.simulate_family_outcomes
    else:
        args["founder"] = seed
        method = ge.simulate_lineage
    with pytest.raises(ValueError, match="duplex family simulation is not supported"):
        method(**args)
