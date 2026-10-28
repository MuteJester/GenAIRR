"""Unfitted equal-weight controls retain history identity and upstream laws."""
from collections import Counter

import pytest

from GenAIRR import Experiment, _engine
from GenAIRR.experimental import DuplexHistory
from tests.support.cases.duplex import history, model, references, script


def _choices(outcome):
    return [(choice.address, choice.value) for choice in outcome.trace().choices()]


@pytest.mark.parametrize("strict", [False, True])
@pytest.mark.parametrize("productive", [False, True])
def test_omitted_history_weights_match_explicit_ones_through_durable_replay(
    tmp_path, strict, productive
):
    refs = references()
    repairs = [("blunt", script()), ("restored", script(restore=True))]
    implicit = model(refs, [DuplexHistory(name, repair, "V1*01", "D*01", "J*01")
                            for name, repair in repairs])
    explicit = model(refs, [DuplexHistory(name, repair, "V1*01", "D*01", "J*01", 1.0)
                            for name, repair in repairs])
    assert implicit.to_json() == explicit.to_json()

    def compile_model(source):
        experiment = Experiment.on(refs).recombine(duplex_repair=source)
        experiment.mutate(model="uniform", count=1).end_loss_3prime(length=1)
        if productive:
            experiment.productive_only()
        return experiment.compile().simulator

    default_simulator = compile_model(implicit)
    explicit_simulator = compile_model(explicit)
    for seed in (0, 17, 31, 1_000_003):
        default = default_simulator.run(seed, strict=strict)
        authored = explicit_simulator.run(seed, strict=strict)
        assert _choices(default) == _choices(authored)
        expected = default.airr_record(refs)
        assert expected == authored.airr_record(refs)
        assert expected["duplex_repair"]["selection"]["natural_weight"] == 1.0
        default_trace = default_simulator.trace_file_from(default, seed)
        explicit_trace = explicit_simulator.trace_file_from(authored, seed)
        assert default_trace.to_json() == explicit_trace.to_json()

        target = tmp_path / f"unfitted-default-{seed}.json"
        default_trace.write_to(str(target))
        loaded = _engine.TraceFile.read_from(str(target))
        # Each simulator must also accept the other's durable model identity.
        for simulator in (default_simulator, explicit_simulator):
            for replay in (simulator.replay_from_trace_file, simulator.rerun_from_trace_file):
                reconstructed = replay(loaded, strict=strict)
                assert reconstructed.airr_record(refs) == expected
                assert _choices(reconstructed) == _choices(default)
                assert reconstructed.validate_record(refs) == []


@pytest.mark.parametrize("productive", [False, True])
def test_unfitted_uniform_history_ids_keep_separate_mass_for_identical_dna(productive):
    refs = references()
    rows = [history("blunt-a"), history("restored", restore=True), history("blunt-b")]
    assert {row.script.run()["sequence"] for row in rows} == {"AAATGTCCCTTTAAA"}
    experiment = Experiment.on(refs).recombine(duplex_repair=model(refs, rows))
    if productive:
        experiment.productive_only()
    compiled = experiment.compile()
    counts = Counter()
    observed_ids = set()
    for outcome in compiled.stream(n=4096, seed=37, strict=True):
        row_index = outcome.trace().find("duplex_repair.v1.history").value
        counts[row_index] += 1
        if rows[row_index].id not in observed_ids:
            record = outcome.airr_record(refs)
            selection = record["duplex_repair"]["selection"]
            assert selection["history_id"] == rows[row_index].id
            assert selection["natural_weight"] == 1.0
            assert record["sequence"] == "AAATGTCCCTTTAAA"
            observed_ids.add(selection["history_id"])
    assert observed_ids == {"blunt-a", "restored", "blunt-b"}
    assert set(counts) == {0, 1, 2}
    for count in counts.values():
        assert count / sum(counts.values()) == pytest.approx(1 / 3, abs=0.035)


@pytest.mark.parametrize("v_weights,expected_v2", [(None, 0.5), ({"V2*01": 3.0}, 0.75)])
def test_unequal_uniform_history_counts_do_not_reweight_earlier_allele_choices(
    v_weights, expected_v2
):
    refs = references(v_count=2)
    balanced = [history("v1"), history("v2", v_call="V2*01")]
    unequal = [history("v1")] + [
        history(f"v2-{i}", v_call="V2*01", restore=bool(i % 2)) for i in range(5)
    ]

    def compile_rows(rows):
        return (Experiment.on(refs)
                .recombine(duplex_repair=model(refs, rows), v_allele_weights=v_weights)
                .productive_only().compile())

    first, second = compile_rows(balanced), compile_rows(unequal)
    allele_counts = Counter()
    history_counts = Counter()
    for baseline, expanded in zip(first.stream(n=4096, seed=71, strict=True),
                                  second.stream(n=4096, seed=71, strict=True), strict=True):
        # Additional authored support changes neither the chosen generating
        # allele nor the RNG consumed by any earlier source draw.
        for segment in ("v", "d", "j"):
            address = f"sample_allele.{segment}"
            assert baseline.trace().find(address).value == expanded.trace().find(address).value
        v = expanded.trace().find("sample_allele.v").value
        row_index = expanded.trace().find("duplex_repair.v1.history").value
        allele_counts[v] += 1
        history_counts[row_index] += 1
        if v == 0:
            assert row_index == 0
        else:
            assert row_index in range(1, 6)
    assert allele_counts[1] / sum(allele_counts.values()) == pytest.approx(expected_v2, abs=0.035)
    assert set(history_counts) == set(range(6))
