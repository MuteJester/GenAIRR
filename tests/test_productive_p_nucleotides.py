"""Positive productive-only contracts for deterministic and sampled P additions."""

import pytest

import GenAIRR as ga
from GenAIRR._airr_record import outcome_to_airr_record
from GenAIRR.reference_models import NpBaseModelSpec

from ._p_productive_support import (
    assert_productive_record,
    reverse_complement,
    synthetic_cartridge,
)


@pytest.mark.parametrize("strict", [False, True])
def test_vj_forced_positive_p_requires_compensating_n_length(strict):
    cfg = synthetic_cartridge({"V_3": [(1, 1.0)]})
    exp = (
        ga.Experiment.on(cfg).recombine(np1_lengths=[(n, 1.0) for n in range(3)])
        .trim(enabled=False).productive_only()
    )
    result = exp.run_records(n=24, seed=42, strict=strict)
    for record in result:
        assert record["p_v_3_length"] == 1
        assert record["sequence"].upper() == "AAACCCTGTACCTTTGGG"
        assert record["np1"].upper() == "ACC"
        assert_productive_record(record, cfg)
    assert result.validate_records(exp.refdata).ok


@pytest.mark.parametrize("invert", [False, True])
@pytest.mark.parametrize("strict", [False, True])
def test_vdj_all_four_p_ends_jointly_constrain_both_n_lengths(invert, strict):
    cfg = synthetic_cartridge(
        {end: [(1, 1.0)] for end in ("V_3", "D_5", "D_3", "J_5")}, d="CCC",
    )
    exp = (
        ga.Experiment.on(cfg)
        .recombine(np1_lengths=[(0, 1), (1, 1)], np2_lengths=[(0, 1), (1, 1)])
        .trim(enabled=False)
    )
    if invert:
        exp.invert_d(prob=1.0)
    exp.productive_only()
    result = exp.run_records(n=24, seed=42, strict=strict)
    for record in result:
        for field in ("p_v_3_length", "p_d_5_length", "p_d_3_length", "p_j_5_length"):
            assert record[field] == 1
        assert record["np1"].upper() == ("ACC" if invert else "ACG")
        assert record["np2"].upper() == ("CCA" if invert else "GCA")
        assert record["d_inverted"] is invert
        assert_productive_record(record, cfg)
    assert result.validate_records(exp.refdata).ok


@pytest.mark.parametrize("p_end", ["V_3", "J_5"])
@pytest.mark.parametrize("strict", [False, True])
def test_n_base_support_avoids_stop_completed_across_p_boundary(p_end, strict):
    # V=TGT|T, P='A', N=one base. V-side P creates TA[N];
    # J-side P creates T[N]A. Both require N=C/T to avoid a stop.
    cfg = synthetic_cartridge(
        {p_end: [(1, 1.0)]}, v="TGTT", v_anchor=0, np_base=None,
    )
    exp = (
        ga.Experiment.on(cfg).recombine(np1_lengths=[(1, 1)])
        .trim(enabled=False).productive_only()
    )
    result = exp.run_records(n=64, seed=57, strict=strict)
    observed_n = set()
    for record in result:
        gap = record["np1"].upper()
        assert len(gap) == 2
        assert gap[0 if p_end == "V_3" else 1] == "A"
        n_base = gap[1 if p_end == "V_3" else 0]
        assert n_base in "CT"
        observed_n.add(n_base)
        assert record["p_v_3_length" if p_end == "V_3" else "p_j_5_length"] == 1
        assert_productive_record(record, cfg)
    assert observed_n == {"C", "T"}
    assert result.validate_records(exp.refdata).ok


@pytest.mark.parametrize("strict", [False, True])
def test_markov_n_prefix_accounts_for_forced_suffix_and_future_p(strict):
    # V ends T; N can be AA or CA; J-side P is AAA. Choosing the first A
    # commits a future TAA, despite no complete stop at that first choice.
    cfg = synthetic_cartridge(
        {"J_5": [(3, 1.0)]}, v="TGTT", v_anchor=0, np_base=None,
    )
    cfg.reference_models.np_bases["NP1"] = NpBaseModelSpec(
        kind="markov",
        first_base={"A": 1.0, "C": 1.0, "G": 0.0, "T": 0.0},
        transitions={
            previous: {"A": 1.0, "C": 0.0, "G": 0.0, "T": 0.0}
            for previous in "ACGT"
        },
    )
    exp = (
        ga.Experiment.on(cfg).recombine(np1_lengths=[(2, 1)])
        .trim(enabled=False).productive_only()
    )
    result = exp.run_records(n=32, seed=78, strict=strict)
    for record in result:
        assert record["p_j_5_length"] == 3
        assert record["np1"].upper() == "CAAAA"
        assert record["sequence"].upper() == "TGTTCAAAATTTGGG"
        assert_productive_record(record, cfg)
    assert result.validate_records(exp.refdata).ok


@pytest.mark.parametrize("strict", [False, True])
def test_unavoidable_p_stop_fails_explicitly_without_suppressing_forced_p(strict):
    # TGT|T + AA + TTTGGG has a fixed TAA after the Cys anchor.
    cfg = synthetic_cartridge({"V_3": [(2, 1.0)]}, v="TGTT", v_anchor=0)
    exp = (
        ga.Experiment.on(cfg).recombine(np1_lengths=[(0, 1)])
        .trim(enabled=False).productive_only()
    )
    with pytest.raises((ValueError, ga.StrictSamplingError)):
        exp.run_records(n=1, seed=0, strict=strict)


@pytest.mark.parametrize("strict", [False, True])
def test_impossible_fixed_p_frame_fails_explicitly(strict):
    cfg = synthetic_cartridge({"V_3": [(1, 1.0)]})
    exp = (
        ga.Experiment.on(cfg).recombine(np1_lengths=[(0, 1)])
        .trim(enabled=False).productive_only()
    )
    with pytest.raises((ValueError, ga.StrictSamplingError)):
        exp.run_records(n=1, seed=0, strict=strict)


@pytest.mark.parametrize("strict", [False, True])
def test_forced_p_longer_than_retained_template_is_not_silently_changed(strict):
    cfg = synthetic_cartridge({"V_3": [(4, 1.0)]}, v="TGT", v_anchor=0)
    exp = (
        ga.Experiment.on(cfg).recombine(np1_lengths=[(0, 1)])
        .trim(enabled=False).productive_only()
    )
    # Either omitting P or shortening it to 3 would restore frame but
    # violate the user's singleton P-length support.
    with pytest.raises((ValueError, ga.StrictSamplingError)):
        exp.run_records(n=1, seed=0, strict=strict)


def test_mixed_zero_and_positive_p_preserves_both_feasible_paths():
    cfg = synthetic_cartridge({"V_3": [(0, 1.0), (1, 1.0)]})
    exp = (
        ga.Experiment.on(cfg).recombine(np1_lengths=[(0, 1), (2, 1)])
        .trim(enabled=False).productive_only()
    )
    result = exp.run_records(n=96, seed=51, strict=True)
    seen = set()
    for record in result:
        p_length = record["p_v_3_length"]
        seen.add(p_length)
        assert record["np1_length"] == (0 if p_length == 0 else 3)
        assert_productive_record(record, cfg)
    assert seen == {0, 1}, "productive conditioning must not silently force P=0"
    assert result.validate_records(exp.refdata).ok


def test_p_template_and_feasible_n_length_follow_realized_trim():
    # No trim: V=TGTA, P=TA, N=0. Trim 1: V=TGT, P=AC, N=1.
    cfg = synthetic_cartridge({"V_3": [(2, 1.0)]}, v="TGTA", v_anchor=0)
    exp = (
        ga.Experiment.on(cfg).recombine(np1_lengths=[(0, 1), (1, 1)])
        .trim(v_3=[(0, 1), (1, 1)], j_5=[(0, 1)]).productive_only()
    )
    result = exp.run_records(n=64, seed=27, strict=True)
    seen_trims = set()
    for record in result:
        trim = record["v_trim_3"]
        seen_trims.add(trim)
        effective_v = "TGTA"[:4 - trim]
        expected_p = reverse_complement(effective_v[-2:])
        assert record["p_v_3_length"] == 2
        assert record["np1"].upper() == expected_p + "C" * trim
        assert_productive_record(record, cfg)
    assert seen_trims == {0, 1}
    assert result.validate_records(exp.refdata).ok


def test_productive_p_true_trace_replay_preserves_all_annotations():
    cfg = synthetic_cartridge({"V_3": [(0, 1), (1, 1)], "J_5": [(0, 1), (2, 1)]})
    exp = (
        ga.Experiment.on(cfg).recombine(np1_lengths=[(n, 1) for n in range(3)])
        .trim(enabled=False).productive_only()
    )
    compiled = exp.compile()
    seen_positive = False
    for seed in range(12):
        fresh = compiled.simulator.run(seed=seed, strict=True)
        trace = compiled.simulator.trace_file_from(fresh, seed=seed)
        replay = compiled.simulator.replay_from_trace_file(trace, strict=True)
        first = outcome_to_airr_record(fresh, compiled.refdata, sequence_id="same")
        second = outcome_to_airr_record(replay, compiled.refdata, sequence_id="same")
        assert first == second
        assert_productive_record(first, cfg)
        seen_positive |= first["p_v_3_length"] + first["p_j_5_length"] > 0
    assert seen_positive


def test_productive_positive_p_survives_noop_postfork_projection():
    cfg = synthetic_cartridge({"V_3": [(1, 1)], "J_5": [(2, 1)]})
    exp = (
        ga.Experiment.on(cfg).recombine(np1_lengths=[(0, 1)]).trim(enabled=False)
        .productive_only().clonal_repertoire(n_clones=4, max_size=2)
        .sequencing_errors(count=0)
    )
    result = exp.run_records(seed=42, strict=True)
    for record in result:
        assert record["p_v_3_length"] == 1
        assert record["p_j_5_length"] == 2
        assert record["np1"].upper() == "AAA"
        assert_productive_record(record, cfg)
    assert result.validate_records(exp.refdata).ok
