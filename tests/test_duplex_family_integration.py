"""Public family boundaries and literal contiguous-coding-span controls."""
import warnings

import pytest

from GenAIRR import Experiment, RefDataConfig
from GenAIRR.experimental import CoupledDuplexRepair, DuplexHistory, DuplexJoin, DuplexRepairModel
from tests.duplex_experiment_cases import history, model, references


def _family(exp, kind):
    if kind == "fixed":
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            return exp.expand_clones(n_clones=2, per_clone=3)
    if kind == "repertoire":
        return exp.clonal_repertoire(n_clones=2, max_size=3)
    return exp.clonal_lineage(n_clones=2, max_generations=2, n_max=6,
                             n_sample=6, rate=0, lambda_base=3)


@pytest.mark.parametrize("kind", ["fixed", "repertoire", "lineage"])
@pytest.mark.parametrize("strict", [False, True])
def test_duplex_source_is_once_per_family_and_nullable_truth_survives(kind, strict):
    refs = references()
    exp = Experiment.on(refs).recombine(duplex_repair=model(refs, [
        history("blunt"), history("restored", restore=True),
    ]))
    result = _family(exp, kind).productive_only().run_records(
        seed=19, strict=strict, validate_records=True,
    )
    assert len(result) > 0
    by_clone = {}
    for rec in result:
        assert rec["sequence"] == "AAATGTCCCTTTAAA"
        assert rec["productive"] is True
        assert rec["p_v_3_length"] is None
        assert rec["n_mutations"] == 0
        selected = rec["duplex_repair"]["selection"]["history_id"]
        by_clone.setdefault(rec["clone_id"], set()).add(selected)
    assert len(by_clone) == 2
    assert all(len(ids) == 1 for ids in by_clone.values())
    assert result.validate_records(refs).ok
    assert result.validate_families().ok


def _crossing_np_fixture():
    refs = references()
    repair = CoupledDuplexRepair("AAATGT", "CCC", "TTTAAA")
    for name, left, right, shift in (("vd", "A", "T", 7), ("dj", "CT", "AG", 5)):
        repair = repair.join(
            DuplexJoin(name).open("left", offset=0).open("right", offset=0)
            .add_n("left", left).add_n("right", right)
            .synapse(right_start=shift).ligate()
        )
    authored = DuplexRepairModel.bind(refs, [
        DuplexHistory("crossing", repair, "V1*01", "D*01", "J*01"),
    ], reference_ends_are_coding_ends=True)
    return refs, authored


def test_duplex_lineage_affinity_translates_across_both_np_boundaries():
    refs, authored = _crossing_np_fixture()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        result = (Experiment.on(refs).recombine(duplex_repair=authored)
                  .clonal_lineage(n_clones=1, max_generations=1, n_max=8, n_sample=8,
                                  rate=0, lambda_base=4, target_aa="KCTPFK")
                  .productive_only().run_records(seed=7, validate_records=True))
    assert len(result) > 0
    for rec in result:
        assert rec["sequence"] == "AAATGTACCCCTTTTAAA"
        assert rec["np1"] == "A" and rec["np2"] == "CT"
        assert rec["lineage_affinity"] == 1.0
    assert all(node.affinity == 1.0 for node in result.lineage_trees[0].nodes())


def test_unresolved_duplex_phase_does_not_prevent_neutral_growth():
    refs = RefDataConfig.vdj()
    refs.add_v_allele("V*01", "V", b"ACG")
    refs.add_d_allele("D*01", "D", b"AAT")
    refs.add_j_allele("J*01", "J", b"TCC")
    repair = CoupledDuplexRepair("ACG", "AAT", "TCC")
    for name in ("vd", "dj"):
        repair = repair.join(DuplexJoin(name).open("left", offset=0).open("right", offset=0)
                             .synapse(right_start=3).ligate())
    authored = DuplexRepairModel.bind(refs, [
        DuplexHistory("no_phase", repair, "V*01", "D*01", "J*01"),
    ], reference_ends_are_coding_ends=True)
    result = (Experiment.on(refs).allow_curatable_refdata().recombine(duplex_repair=authored)
              .clonal_lineage(n_clones=1, max_generations=1, n_max=8, n_sample=8,
                              rate=0, lambda_base=4)
              .run_records(seed=7, validate_records=True))
    assert len(result) > 0
    assert all(rec["productive"] is None and rec["lineage_affinity"] == 0.0 for rec in result)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        scoring = (Experiment.on(refs).allow_curatable_refdata().recombine(duplex_repair=authored)
                   .clonal_lineage(n_clones=1, max_generations=1, target_aa="TIS"))
    with pytest.raises(ValueError, match="(?i)(phase|coding|affinity|span)"):
        scoring.run_records(seed=7)


@pytest.mark.parametrize("kind", ["fixed", "repertoire", "lineage"])
def test_unsupported_observation_adapter_fails_before_family_generation(kind):
    refs = references()
    exp = _family(Experiment.on(refs).recombine(duplex_repair=model(refs)), kind)
    exp.contaminate(prob=0)
    with pytest.raises(ValueError, match="(?i)(duplex|unsupported|contamin)"):
        exp.compile()
