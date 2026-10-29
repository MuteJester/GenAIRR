import json
import math

import pytest

import GenAIRR as ga

pytestmark = pytest.mark.behavior


def truth(v=0, d=(0, 0), j=0):
    return {"coding_end_processing": {
        "model": "original_coding_end_v1",
        "boundary_contract": "provided_allele_coding_ends",
        "recombination_complete": True,
        "choices": [{"segment": s, "original_allele": s+"*01",
                     "five_prime_outcome": a, "three_prime_outcome": b}
                    for s, a, b in (("V", 0, v), ("D", *d), ("J", j, 0))]}}


HASH = "sha256:" + "a" * 64


def fit(rows, **kwargs):
    return ga.CodingEndOutcomeModel.fit_generating_records(
        rows, chain_type="vdj", reference_hash=HASH, source="literal fixture",
        training_sampling_law="unconditioned_generating_choices", independent_generating_events=True, **kwargs)


def test_joint_d_recovery_keeps_correlation_and_serializes():
    rows = [truth(d=(-1, 1))] * 75 + [truth(d=(1, -1))] * 25
    model = fit(rows)
    assert model.d == ((-1, 1, .75), (1, -1, .25))
    assert ga.CodingEndOutcomeModel.from_dict(json.loads(json.dumps(model.to_dict()))) == model
    assert model.provenance["counts"]["D"][0]["count"] == 75
    heldout = model.score_generating_records([truth(d=(-1, 1)), truth(d=(1, -1))])
    assert heldout["log_likelihood"] == pytest.approx(math.log(.75) + math.log(.25))
    copy = model.provenance
    copy["counts"].clear()
    assert model.provenance["counts"]


@pytest.mark.parametrize("bad", [True, 1.0, 2**63, -(2**63)-1])
def test_signed_outcomes_are_typed(bad):
    with pytest.raises(ValueError, match="signed 64-bit"):
        ga.CodingEndOutcomeModel(v=((0, bad, 1),), j=((0, 0, 1),))


@pytest.mark.parametrize("bad", [True, 0, -1, math.nan, math.inf])
def test_model_weights_reject_invalid(bad):
    with pytest.raises(ValueError):
        ga.CodingEndOutcomeModel(v=((0, 0, bad),), j=((0, 0, 1),))


def test_training_rejects_unknown_histories_and_selected_law():
    with pytest.raises(ValueError, match="generating provenance"):
        fit([{"junction": "AACGT", "p_v_3_length": 1}])
    row = truth()
    row["coding_end_processing"]["recombination_complete"] = False
    with pytest.raises(ValueError, match="completed"):
        fit([row])
    with pytest.raises(ValueError, match="ascertainment"):
        ga.CodingEndOutcomeModel.fit_generating_records(
            [truth()], chain_type="vdj", reference_hash=HASH, source="selected",
            training_sampling_law="locally_conditioned_productive")


def test_prior_support_is_explicit_and_unknown_outcomes_are_reported():
    support = ga.CodingEndOutcomeModel(v=((0, 0, 1),), j=((0, 0, 1),),
                                     d=((0, 0, 1), (1, 0, 1)))
    model = fit([truth()], support=support, pseudocount=1)
    assert model.d == ((0, 0, 2/3), (1, 0, 1/3))
    with pytest.raises(ValueError, match="outside declared support"):
        fit([truth(d=(2, 0))], support=support)
    unseen = model.score_generating_records([truth(d=(2, 0))])
    assert unseen["unseen_outcomes"] == 1
    assert unseen["log_likelihood"] == -math.inf


def ref():
    r = ga.RefDataConfig.vdj()
    r.add_v_allele("v*01", "v", b"ATGAAATGTAAA", anchor=6)
    r.add_d_allele("d*01", "d", b"AACCGGTT")
    r.add_j_allele("j*01", "j", b"AAATGGAAA", anchor=3)
    return r


def test_model_flows_through_real_engine_and_replay():
    r = ref()
    model = ga.CodingEndOutcomeModel(v=((0, -1, 1),), j=((0, 0, 1),),
                                    d=((-1, 1, 1), (1, -1, 1)),
                                    reference_hashes=(r.content_hash(),))
    compiled = (ga.Experiment.on(r).recombine(np1_lengths=[(0, 1)], np2_lengths=[(0, 1)])
                .process_coding_ends(model=model, reference_ends_are_coding_ends=True)
                .invert_d(prob=.5).compile())
    generated = []
    for seed in range(20):
        outcome = compiled.simulator.run(seed=seed, strict=True)
        row = outcome.airr_record(compiled.refdata)
        generated.append(row)
        trace = compiled.simulator.trace_file_from(outcome, seed=seed)
        replayed = compiled.simulator.replay_from_trace_file(trace, strict=True)
        assert replayed.airr_record(compiled.refdata) == row
    fitted = ga.CodingEndOutcomeModel.fit_generating_records(
        generated, chain_type="vdj", reference_hash=r.content_hash(),
        source="native synthetic generating histories",
        training_sampling_law="unconditioned_generating_choices", independent_generating_events=True)
    assert {(a, b) for a, b, _ in fitted.d} == {(-1, 1), (1, -1)}
    assert all(entry["marginal_95_interval"][0] <= .5 <= entry["marginal_95_interval"][1]
               for entry in fitted.provenance["counts"]["D"])


def test_model_reference_mismatch_and_mixed_inputs_are_rejected():
    model = ga.CodingEndOutcomeModel(v=((0, 0, 1),), j=((0, 0, 1),),
                                    d=((0, 0, 1),), reference_hashes=(HASH,))
    experiment = ga.Experiment.on(ref()).recombine()
    with pytest.raises(ValueError, match="incompatible"):
        experiment.process_coding_ends(model=model, reference_ends_are_coding_ends=True)
    with pytest.raises(ValueError, match="excludes individual"):
        experiment.process_coding_ends(model=model, v_3=[(0, 1)], reference_ends_are_coding_ends=True)
    with pytest.raises(ValueError, match="excludes separate"):
        experiment.process_coding_ends(d_joint=[(0, 0, 1)], d_5=[(0, 1)], reference_ends_are_coding_ends=True)


def test_d_joint_never_expands_to_cartesian_product():
    experiment = (ga.Experiment.on(ref()).recombine()
                  .process_coding_ends(d_joint=[(-1, 1, 2), (1, -1, 3)],
                                       reference_ends_are_coding_ends=True))
    assert experiment._steps[0].coding_end_outcomes_d == ((-1, 1, 2.0), (1, -1, 3.0))


def test_compatible_history_likelihood_sums_ambiguous_explanations():
    model = ga.CodingEndOutcomeModel(v=((0, 0, 1), (0, 1, 1)), j=((0, 0, 1),))
    def observation(a, b):
        return [{"outcomes": {"V": (0, 0), "J": (0, 0)}, "likelihood": a},
                {"outcomes": {"V": (0, 1), "J": (0, 0)}, "likelihood": b}]
    fitted = model.fit_compatible_histories(
        [observation(.9, .1)] * 75 + [observation(.1, .9)] * 25,
        source="independent identifiable two-component observation fixture",
        training_sampling_law="unconditioned_generating_choices",
        independent_observations=True, complete_scenarios=True)
    assert fitted.v[0][2] == pytest.approx(.8125, abs=1e-4)
    report = fitted.provenance
    assert report["converged"]
    assert report["log_likelihood_trajectory"][-1] == pytest.approx(75*math.log(.75)+25*math.log(.25), abs=1e-6)
    assert all(b >= a-1e-10 for a,b in zip(report["log_likelihood_trajectory"], report["log_likelihood_trajectory"][1:]))
    assert report["mean_history_posterior_entropy"] > 0
    assert report["uncertainty_intervals"] is None


def test_equivalent_histories_do_not_create_identifiability_claims():
    model = ga.CodingEndOutcomeModel(v=((0, 0, 3), (0, 1, 1)), j=((0, 0, 1),))
    scenarios = [{"outcomes": {"V": (0, p), "J": (0, 0)}, "likelihood": 1}
                 for p in (0, 1)]
    fitted = model.fit_compatible_histories(
        [scenarios]*4, source="unidentifiable fixture",
        training_sampling_law="unconditioned_generating_choices",
        independent_observations=True, complete_scenarios=True)
    assert fitted.v[0][2] == pytest.approx(.75)
    assert "not established" in fitted.provenance["identifiability"]
    with pytest.raises(ValueError, match="complete_scenarios"):
        model.fit_compatible_histories([scenarios], source="missing declarations",
                                      training_sampling_law="unconditioned_generating_choices")


def test_independent_events_and_clonal_accounting_are_explicit():
    row = truth()
    row["clone_id"] = 7
    with pytest.raises(ValueError, match="repeat generating"):
        fit([row])
    with pytest.raises(ValueError, match="independent_generating_events"):
        ga.CodingEndOutcomeModel.fit_generating_records(
            [truth()], chain_type="vdj", reference_hash=HASH, source="missing declaration",
            training_sampling_law="unconditioned_generating_choices")


def test_model_hash_is_rechecked_after_reference_changes():
    r = ref()
    model = ga.CodingEndOutcomeModel(v=((0, 0, 1),), j=((0, 0, 1),),
                                    d=((0, 0, 1),), reference_hashes=(r.content_hash(),))
    exp = ga.Experiment.on(r).recombine().process_coding_ends(
        model=model, reference_ends_are_coding_ends=True)
    exp._refdata.add_v_allele("v2*01", "v2", b"ATGAAATGTAAA", anchor=6)
    with pytest.raises(ValueError, match="final compiled reference"):
        exp.compile()
