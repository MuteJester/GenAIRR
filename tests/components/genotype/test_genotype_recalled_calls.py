"""Donor-restricted inference must survive SHM and trace replay."""

import GenAIRR as ga
import pytest
from GenAIRR._airr_record import outcome_to_airr_record
from GenAIRR.genotype import Genotype


def _allowed(cfg, genotype):
    return {
        segment.lower(): set().union(*[
            genotype.carried_alleles(segment, gene)
            for gene in getattr(cfg, f"{segment.lower()}_alleles")
        ])
        for segment in ("V", "D", "J")
    }


def _assert_calls_and_parity(result, refdata, allowed):
    for record, outcome in zip(result.records, result.outcomes):
        for segment, names in allowed.items():
            assert record[f"{segment}_call"] in names
            recalled = set(filter(None, record[f"recalled_{segment}_call"].split(",")))
            assert recalled <= names
        for parity in outcome.check_live_call_cache_parity(refdata):
            assert parity["tie_set_matches"]
            assert parity["hypothesis_bounds_match"] is not False
    assert result.validate_records(refdata).ok


def test_high_shm_cohort_recalled_calls_are_restricted_to_each_donor():
    cfg = ga.HUMAN_IGH_OGRDB
    first = Genotype.from_dataconfig(cfg).complete_from_reference().with_subject("first")
    last = Genotype.from_dataconfig(cfg).complete_from_reference().with_subject("last")
    for segment in ("V", "D", "J"):
        for gene, alleles in getattr(cfg, f"{segment.lower()}_alleles").items():
            last.homozygous(gene, alleles[-1].name, segment=segment)
    cohort = (
        ga.Experiment.on(cfg).recombine().mutate(model="uniform", rate=0.3)
        .polymerase_indels(count=2).primer_trim_5prime(length=4)
        .run_cohort([first, last], n_per_subject=24, seed=11)
    )
    for genotype in (first, last):
        _assert_calls_and_parity(
            cohort.result_for(genotype.subject_id), cohort.refdata_for(genotype.subject_id),
            _allowed(cfg, genotype),
        )


def test_donor_call_support_survives_true_trace_replay():
    cfg = ga.HUMAN_IGH_OGRDB
    genotype = Genotype.from_dataconfig(cfg).complete_from_reference()
    exp = ga.Experiment.on(cfg).with_genotype(genotype).recombine().mutate(rate=0.2)
    compiled = exp.compile()
    fresh = compiled.simulator.run(seed=13)
    trace = compiled.simulator.trace_file_from(fresh, seed=13)
    replay = compiled.simulator.replay_from_trace_file(trace)
    first = outcome_to_airr_record(fresh, compiled.refdata, sequence_id="same")
    second = outcome_to_airr_record(replay, compiled.refdata, sequence_id="same")
    assert first == second
    for parity in replay.check_live_call_cache_parity(compiled.refdata):
        assert parity["tie_set_matches"]
        assert parity["hypothesis_bounds_match"] is not False


def test_genotype_clonal_forks_remain_explicitly_unsupported():
    cfg = ga.HUMAN_IGH_OGRDB
    genotype = Genotype.from_dataconfig(cfg).complete_from_reference()
    clonal = (
        ga.Experiment.on(cfg).with_genotype(genotype).recombine()
        .clonal_repertoire(n_clones=5, max_size=3).mutate(rate=0.15)
    )
    with pytest.raises(ValueError, match=r"with_genotype\(\) is not supported"):
        clonal.compile()


@pytest.mark.parametrize("end", ["5prime", "3prime"])
def test_fully_clipped_reads_keep_generating_calls_but_no_recalled_calls(end):
    cfg = ga.HUMAN_IGH_OGRDB
    genotype = Genotype.from_dataconfig(cfg).complete_from_reference()
    exp = ga.Experiment.on(cfg).with_genotype(genotype).recombine()
    exp = getattr(exp, f"primer_trim_{end}")(length=10000)
    result = exp.run_records(n=1, seed=157, expose_provenance=True)
    record = result[0]
    assert record["sequence"] == ""
    for segment in ("v", "d", "j"):
        assert record[f"{segment}_call"]
        assert record[f"{segment}_call"] == record[f"truth_{segment}_call"]
        assert record[f"recalled_{segment}_call"] == ""
    _assert_calls_and_parity(result, exp.refdata, _allowed(cfg, genotype))
