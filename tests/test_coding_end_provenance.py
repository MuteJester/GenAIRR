"""Original-tip ancestry remains meaningful after later biological/read edits."""

import pytest

import GenAIRR as ga
from GenAIRR.alleles import VAllele

from ._p_productive_support import reverse_complement, synthetic_cartridge


def experiment(*, d=None):
    cfg = synthetic_cartridge({}, d=d)
    exp = ga.Experiment.on(cfg).recombine(np1_lengths=[(0, 1)], **(
        {"np2_lengths": [(0, 1)]} if d is not None else {}
    ))
    return exp.process_coding_ends(
        v_3=[(-1, 1)], j_5=[(-2, 1)], reference_ends_are_coding_ends=True,
    )


def check_observed_origins(record):
    metadata = record["coding_end_processing"]
    coding_sequence = reverse_complement(record["sequence"]) if record["rev_comp"] else record["sequence"]
    assert metadata["coordinate_system"] == "zero_based_coding_oriented"
    counts = dict.fromkeys(("V_3", "D_5", "D_3", "J_5"), 0)
    positions = []
    for origin in metadata["p_origins"]:
        position = origin["query_position"]
        if position is None:
            assert origin["current_base"] is None
        else:
            assert origin["current_base"] == coding_sequence[position]
            counts[origin["assembled_end"]] += 1
            positions.append(position)
    assert len(positions) == len(set(positions))
    assert counts == metadata["observed_p"]
    assert metadata["retained_p_at_join"] == {"V_3": 1, "D_5": 0, "D_3": 0, "J_5": 2}


def test_masking_preserves_original_p_identity_and_original_template():
    exp = experiment().ambiguous_base_calls(count=1000)
    result = exp.run_records(n=2, seed=11)
    for record in result:
        check_observed_origins(record)
        origins = record["coding_end_processing"]["p_origins"]
        assert len(origins) == 3
        assert all(origin["current_base"] == "N" for origin in origins)
        assert [origin["original_base"] for origin in origins] == ["A", "A", "A"]
        assert [origin["template_position"] for origin in origins] == [8, 1, 0]
    assert result.validate_records(exp.refdata).ok


@pytest.mark.parametrize("end", ["5prime", "3prime"])
def test_complete_read_loss_preserves_historical_p_and_clears_current_positions(end):
    exp = getattr(experiment(), f"end_loss_{end}")(length=10000)
    result = exp.run_records(n=1, seed=3)
    record = result[0]
    assert record["sequence"] == ""
    check_observed_origins(record)
    assert all(origin["query_position"] is None for origin in record["coding_end_processing"]["p_origins"])
    assert record["p_v_3_length"] == 1
    assert record["p_j_5_length"] == 2
    assert result.validate_records(exp.refdata).ok


def test_indels_read_loss_and_reverse_orientation_keep_coding_coordinate_ancestry_on_replay():
    compiled = (
        experiment().mutate(model="uniform", rate=0.2)
        .polymerase_indels(count=4).end_loss_5prime(length=5)
        .random_strand_orientation(prob=1).compile()
    )
    for seed in range(8):
        original = compiled.simulator.run(seed=seed)
        trace = compiled.simulator.trace_file_from(original, seed=seed)
        replayed = compiled.simulator.replay_from_trace_file(trace)
        record = original.airr_record(compiled.refdata)
        assert record == replayed.airr_record(compiled.refdata)
        check_observed_origins(record)
        assert original.validate_record(compiled.refdata) == []


def test_receptor_revision_preserves_original_v_tip_allele_identity():
    cfg = synthetic_cartridge({}, d="ACGTAC")
    alternate = VAllele("IGHV-PTEST2*01", "AAACCCTGC", 9, anchor_override=6)
    cfg.v_alleles[alternate.gene] = [alternate]
    exp = (
        ga.Experiment.on(cfg).recombine(np1_lengths=[(0, 1)], np2_lengths=[(0, 1)])
        .process_coding_ends(v_3=[(-1, 1)], j_5=[(-2, 1)], reference_ends_are_coding_ends=True)
        .receptor_revision(prob=1)
    )
    result = exp.run_records(n=8, seed=3)
    for record in result:
        check_observed_origins(record)
        assert record["receptor_revision_applied"]
        original = record["original_v_call"]
        assert original != record["v_call"]
        metadata = record["coding_end_processing"]
        assert next(choice for choice in metadata["choices"] if choice["segment"] == "V")["original_allele"] == original
        assert next(origin for origin in metadata["p_origins"] if origin["source_segment"] == "V")["original_allele"] == original
    assert result.validate_records(exp.refdata).ok


def test_clonal_continuation_inherits_choices_and_lost_p_ancestry():
    exp = experiment().clonal_repertoire(n_clones=2, max_size=2).end_loss_5prime(length=10)
    result = exp.run_records(seed=11)
    for record in result:
        check_observed_origins(record)
        metadata = record["coding_end_processing"]
        assert len(metadata["choices"]) == 2
        assert metadata["observed_p"]["V_3"] == 0
        assert metadata["p_origins"][0]["query_position"] is None
        assert [origin["query_position"] for origin in metadata["p_origins"][1:]] == [0, 1]
    assert result.validate_records(exp.refdata).ok
