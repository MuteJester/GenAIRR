"""Original-end processing: public contracts, sequence biology, and replay."""

import warnings

import pytest

import GenAIRR as ga
from GenAIRR import _engine

from tests.support.cases.productive import (
    assert_productive_record,
    reverse_complement,
    synthetic_cartridge,
)

pytestmark = pytest.mark.behavior


def experiment(*, v="AAACCCTGT", d=None, np1=((0, 1.0),), np2=((0, 1.0),)):
    cfg = synthetic_cartridge({}, v=v, d=d)
    lengths = {"np1_lengths": np1}
    if d is not None:
        lengths["np2_lengths"] = np2
    return ga.Experiment.on(cfg).recombine(**lengths)


def process(exp, **kwargs):
    return exp.process_coding_ends(reference_ends_are_coding_ends=True, **kwargs)


@pytest.mark.parametrize("acknowledgement", [False, None, 1, "yes"])
def test_reference_end_boundary_requires_explicit_true(acknowledgement):
    with pytest.raises(ValueError, match="reference_ends_are_coding_ends=True"):
        experiment().process_coding_ends(
            v_3=[(-1, 1)], reference_ends_are_coding_ends=acknowledgement
        )


def test_configuration_requires_recombination_and_precedes_later_stages():
    with pytest.raises(ValueError, match="after recombine"):
        process(ga.Experiment.on(synthetic_cartridge({})))
    with pytest.raises(ValueError, match="precede"):
        process(experiment().mutate(count=0))
    with pytest.raises(ValueError, match="precede"):
        process(experiment().end_loss_3prime(length=0))
    with pytest.raises(ValueError, match="precede"):
        process(experiment().clonal_repertoire(n_clones=1))
    with pytest.raises(ValueError, match="already configured"):
        process(process(experiment()))


@pytest.mark.parametrize("enabled", [False, True])
def test_explicit_trim_conflicts_in_both_orders(enabled):
    with pytest.raises(ValueError, match="mutually exclusive"):
        process(experiment().trim(enabled=enabled))
    with pytest.raises(ValueError, match="mutually exclusive"):
        process(experiment()).trim(enabled=enabled)


@pytest.mark.parametrize("end", ["d_5", "d_3"])
def test_vj_rejects_even_zero_explicit_d_outcomes(end):
    with pytest.raises(ValueError, match="only valid for VDJ"):
        process(experiment(), **{end: [(0, 1)]})


@pytest.mark.parametrize(
    "pairs",
    [
        [], [(True, 1)], [(1.5, 1)], [(2**63, 1)],
        [(0, 0)], [(0, -1)], [(0, float("nan"))], [(0, float("inf"))],
        [(0, "invalid")], [(0, 1, 2)],
    ],
)
def test_invalid_outcomes_fail_before_native_compilation(pairs):
    with pytest.raises(ValueError):
        process(experiment(), v_3=pairs)


def test_input_iterables_are_snapshotted_and_description_reports_signed_outcomes():
    pairs = [[-1, 2.0], [0, 1.0]]
    exp = process(experiment(), v_3=iter(pairs))
    pairs[0][0] = -99
    step = exp._steps[0]
    assert step.coding_end_outcomes_v == ((0, -1, 2.0), (0, 0, 1.0))
    assert step.coding_end_outcomes_j == ((0, 0, 1.0),)
    assert "original coding-end processing" in exp.describe()
    assert "V3'=-1,0" in exp.describe()
    assert "reference-end boundary confirmed" in exp.describe()


def test_d_outcomes_are_joint_canonical_cartesian_support():
    exp = process(experiment(d="ACGTAC"), d_5=[(0, 2), (1, 3)], d_3=[(-2, 5), (0, 7)])
    assert exp._steps[0].coding_end_outcomes_d == (
        (0, -2, 10.0), (0, 0, 14.0), (1, -2, 15.0), (1, 0, 21.0),
    )


def test_native_api_requires_boundary_confirmation_and_valid_joint_ends():
    plan = _engine.PassPlan()
    with pytest.raises(ValueError, match="reference_ends_are_coding_ends=True"):
        plan.push_process_coding_ends("V", [(0, -1, 1)])
    with pytest.raises(ValueError):
        plan.push_process_coding_ends("V", [(1, 0, 1)], reference_ends_are_coding_ends=True)
    with pytest.raises(ValueError):
        plan.push_process_coding_ends("J", [(0, 1, 1)], reference_ends_are_coding_ends=True)
    with pytest.raises(ValueError):
        plan.push_emit_coding_end_p("V_5")


@pytest.mark.parametrize("five,three", [(False, -1), (0, True), (0, 1.0), (0, 1.5), (0, "1")])
def test_native_signed_outcomes_reject_bool_and_noninteger_coercion(five, three):
    plan = _engine.PassPlan()
    with pytest.raises((ValueError, TypeError)):
        plan.push_process_coding_ends(
            "V", [(five, three, 1.0)], reference_ends_are_coding_ends=True,
        )


def test_original_v_and_j_tip_extensions_have_exact_bases_and_np_boundaries():
    exp = process(experiment(), v_3=[(-1, 1)], j_5=[(-2, 1)])
    result = exp.run_records(n=3, seed=7)
    for rec in result:
        assert rec["sequence"].upper() == "AAACCCTGT" + "A" + "AA" + "TTTGGG"
        assert rec["np1"].upper() == "AAA"
        assert rec["np1_length"] == 3
        assert rec["p_v_3_length"] == 1
        assert rec["p_j_5_length"] == 2
        assert rec["v_trim_3"] == rec["j_trim_5"] == 0
        assert rec["junction_aa"] == "CKF"
        assert rec["coding_end_processing"]["model"] == "original_coding_end_v1"
    assert result.validate_records(exp.refdata).ok


def test_positive_outcomes_delete_reference_and_do_not_create_new_tip_palindromes():
    exp = process(experiment(v="AAACCCTGTAAC"), v_3=[(2, 1)], j_5=[(1, 1)])
    result = exp.run_records(n=1, seed=3)
    rec = result[0]
    assert rec["sequence"].upper() == "AAACCCTGTA" + "TTGGG"
    assert rec["v_trim_3"] == 2
    assert rec["j_trim_5"] == 1
    assert rec["p_v_3_length"] == rec["p_j_5_length"] == 0
    assert result.validate_records(exp.refdata).ok


@pytest.mark.parametrize("inverted", [False, True])
def test_d_processing_uses_canonical_end_before_orientation(inverted):
    exp = process(experiment(d="ACGTAC"), d_5=[(1, 1)], d_3=[(-2, 1)])
    exp.invert_d(prob=float(inverted))
    result = exp.run_records(n=1, seed=4)
    rec = result[0]
    canonical_d = "CGTAC" + reverse_complement("AC")
    observed_d_with_p = reverse_complement(canonical_d) if inverted else canonical_d
    assert rec["sequence"].upper() == "AAACCCTGT" + observed_d_with_p + "TGGAAA"
    assert rec["d_trim_5"] == 1
    assert rec["d_trim_3"] == 0
    assert rec["p_d_5_length"] == (2 if inverted else 0)
    assert rec["p_d_3_length"] == (0 if inverted else 2)
    assert result.validate_records(exp.refdata).ok


def test_chosen_mode_replaces_implicit_cartridge_trim_and_p_models():
    from dataclasses import replace
    from GenAIRR.reference_models import EmpiricalDistributionSpec

    cfg = synthetic_cartridge({"V_3": [(3, 1)], "J_5": [(1, 1)]})
    cfg.reference_models = replace(cfg.reference_models, trims={
        "V_3": EmpiricalDistributionSpec([(2, 1)]),
        "J_5": EmpiricalDistributionSpec([(1, 1)]),
    })
    exp = ga.Experiment.on(cfg).recombine(np1_lengths=[(0, 1)])
    assert exp._steps[0].trim_v_3 == ((2, 1.0),)
    assert exp._steps[0].p_v_3_lengths == ((3, 1.0),)
    process(exp)
    result = exp.run_records(n=1, seed=0)
    rec = result[0]
    assert rec["sequence"].upper() == "AAACCCTGTTTTGGG"
    assert rec["v_trim_3"] == rec["j_trim_5"] == 0
    assert rec["p_v_3_length"] == rec["p_j_5_length"] == 0
    assert not any(name.startswith(("trim.", "p_addition.")) for name in result.outcomes[0].pass_names())


def test_explicit_processing_silences_raw_reference_trim_noop_warning():
    ref = _engine.RefDataConfig.vj()
    ref.add_v_allele("v*01", "v", b"AAACCCTGT", anchor=6)
    ref.add_j_allele("j*01", "j", b"TTTGGG", anchor=0)
    exp = process(ga.Experiment.on(ref).recombine(np1_lengths=[(0, 1)]))
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        exp.compile()
    assert not any("trim" in str(w.message) for w in caught)


def test_productive_only_preserves_fixed_processing_and_conditions_n_support():
    exp = process(experiment(np1=((0, 1), (1, 1), (2, 1))), v_3=[(-1, 1)]).productive_only()
    result = exp.run_records(n=20, seed=11)
    for rec in result:
        assert rec["p_v_3_length"] == 1
        assert rec["np1_length"] == 3  # one P plus two N nucleotides
        assert_productive_record(rec, exp._dataconfig)
    assert result.validate_records(exp.refdata).ok


def test_unavoidable_original_end_stop_fails_without_substituting_zero_p():
    exp = process(experiment(v="AAACCCTGTT"), v_3=[(-2, 1)]).productive_only()
    with pytest.raises(ValueError):
        exp.run_records(n=1, seed=0)


def test_recompilation_and_trace_replay_preserve_processing_and_reject_changed_weights():
    exp = process(experiment(np1=((2, 1),)), v_3=[(-1, 1), (0, 2)])
    compiled = exp.compile()
    original = compiled.simulator.run(seed=9)
    repeated = exp.compile().simulator.run(seed=9)
    assert original.airr_record(compiled.refdata) == repeated.airr_record(compiled.refdata)
    trace = compiled.simulator.trace_file_from(original, seed=9)
    replayed = compiled.simulator.replay_from_trace_file(trace)
    assert original.airr_record(compiled.refdata) == replayed.airr_record(compiled.refdata)
    changed = process(experiment(np1=((2, 1),)), v_3=[(-1, 2), (0, 1)]).compile()
    with pytest.raises(ValueError, match="signature"):
        changed.simulator.replay_from_trace_file(trace)


def test_clonal_plan_split_preserves_processing_once():
    exp = (
        process(experiment(), v_3=[(-1, 1)], j_5=[(-2, 1)])
        .clonal_repertoire(n_clones=2, max_size=2).end_loss_3prime(length=1)
    )
    result = exp.run_records(seed=2)
    for rec in result:
        assert rec["p_v_3_length"] == 1
        assert rec["p_j_5_length"] == 2
        assert rec["sequence"].upper() == "AAACCCTGTAAATTTGG"
        assert len(rec["coding_end_processing"]["choices"]) == 2
    assert result.validate_records(exp.refdata).ok
