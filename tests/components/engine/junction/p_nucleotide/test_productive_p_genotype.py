"""Productive P completion must preserve the genotype's chromosome linkage."""

import pytest

import GenAIRR as ga
from GenAIRR._airr_record import outcome_to_airr_record
from GenAIRR.alleles import JAllele, VAllele
from GenAIRR.genotype import Genotype

from tests.support.cases.productive import assert_productive_record, synthetic_cartridge
from tests.support.fixtures.phased_anchors import _phased_fixture




@pytest.mark.parametrize("strict", [False, True])
@pytest.mark.parametrize("chromosome_weights", [(0.5, 0.5), (0, 1)])
def test_productive_p_conditions_future_alleles_on_the_drawn_haplotype(
    strict, chromosome_weights,
):
    cfg, genotype, expected_v, expected_j = _phased_fixture()
    genotype.chromosome_weights(*chromosome_weights)
    exp = (
        ga.Experiment.on(cfg).with_genotype(genotype)
        .recombine(np1_lengths=[(0, 1)]).trim(enabled=False).productive_only()
    )
    result = exp.run_records(n=32, seed=43, strict=strict, expose_provenance=True)
    for record in result:
        assert record["haplotype"] == 1
        assert record["truth_v_call"] == expected_v
        assert record["truth_j_call"] == expected_j
        assert record["p_v_3_length"] == 1
        assert record["sequence"].upper() == "TGTACCTTTGGG"
        assert_productive_record(record, cfg)
    assert result.validate_records(exp.refdata).ok


@pytest.mark.parametrize("strict", [False, True])
def test_productive_p_cannot_use_alleles_on_a_zero_weight_chromosome(strict):
    cfg, genotype, _, _ = _phased_fixture()
    genotype.chromosome_weights(1, 0)
    exp = (
        ga.Experiment.on(cfg).with_genotype(genotype)
        .recombine(np1_lengths=[(0, 1)]).trim(enabled=False).productive_only()
    )
    with pytest.raises((ValueError, ga.StrictSamplingError)):
        exp.run_records(n=1, seed=0, strict=strict)


def test_productive_p_phased_true_trace_replay_preserves_provenance():
    cfg, genotype, expected_v, expected_j = _phased_fixture()
    exp = (
        ga.Experiment.on(cfg).with_genotype(genotype)
        .recombine(np1_lengths=[(0, 1)]).trim(enabled=False).productive_only()
    )
    compiled = exp.compile()
    fresh = compiled.simulator.run(seed=43, strict=True)
    trace = compiled.simulator.trace_file_from(fresh, seed=43)
    replay = compiled.simulator.replay_from_trace_file(trace, strict=True)
    first = outcome_to_airr_record(fresh, compiled.refdata, sequence_id="same")
    second = outcome_to_airr_record(replay, compiled.refdata, sequence_id="same")
    assert first == second
    for outcome in (fresh, replay):
        sim = outcome.final_simulation()
        assert outcome.trace().find("sample_haplotype").value == 1
        assert compiled.refdata.v_allele(sim.v_allele_id()).name == expected_v
        assert compiled.refdata.j_allele(sim.j_allele_id()).name == expected_j
    assert_productive_record(first, cfg)
