"""Biological metadata regressions found by the October 2026 audit."""
import copy

import pytest

import GenAIRR as ga
from GenAIRR import _engine as ge
from GenAIRR.reference_models import EmpiricalDistributionSpec, ReferenceEmpiricalModels

pytestmark = pytest.mark.behavior


def _revision_fixture(equal_lengths=False):
    ref = ge.RefDataConfig.vdj()
    ref.add_v_allele("v1*01", "v1", b"AAACCCTGT", anchor=6)
    ref.add_v_allele(
        "v2*01", "v2", b"GGGAAATGT" if equal_lengths else b"GGGAAATGTCCC", anchor=6,
    )
    ref.add_d_allele("d1*01", "d1", b"ACGTTA")
    ref.add_j_allele("j1*01", "j1", b"TTTGGG", anchor=0)
    return (
        ga.Experiment.on(ref).allow_curatable_refdata()
        .recombine(np1_lengths=[(0, 1)], np2_lengths=[(0, 1)])
        .trim(enabled=False).receptor_revision(prob=1)
    )


def _revision_seed(applied: bool):
    """A seed whose initial V draw leaves the revision pass able to fire, or not.

    Only `v1*01` has a length-eligible alternate in this fixture, so whether
    the revision applies is decided by which allele the first draw picks. That
    was encoded as "seed 0" and "seed 1", which stopped being true when
    per-pass RNG substreams moved the draw; the premise is the case, so it is
    searched.
    """
    exp = _revision_fixture()
    for seed in range(200):
        result = exp.run_records(n=1, seed=seed, expose_provenance=True)
        if result[0]["receptor_revision_applied"] is applied:
            return exp, result
    raise AssertionError(
        f"no seed in 200 produced a record with receptor_revision_applied="
        f"{applied}, so that half of the pass's behaviour is untested here")


def test_revision_records_current_allele_trim_without_fabricated_deletions():
    exp, result = _revision_seed(applied=True)
    rec = result[0]
    assert rec["original_v_call"] == "v1*01"
    assert rec["v_call"] == "v2*01"
    assert rec["v_trim_3"] == 3
    assert rec["v_germline_end"] == 9
    assert rec["v_cigar"] == "9M"
    assert rec["v_identity"] == 1.0
    assert "-" not in rec["sequence_alignment"]
    assert result.validate_records(exp.refdata).ok


def test_revision_always_selects_a_different_eligible_allele():
    result = _revision_fixture(equal_lengths=True).run_records(
        n=64, seed=0, expose_provenance=True,
    )
    assert all(rec["receptor_revision_applied"] for rec in result)
    assert all(rec["original_v_call"] != rec["truth_v_call"] for rec in result)


def test_revision_without_a_length_eligible_alternate_is_not_applied():
    # The searched seed is one whose initial draw takes the 12-base V; the
    # remaining V has only 9, so there is nothing to revise to.
    _, result = _revision_seed(applied=False)
    rec = result[0]
    assert rec["truth_v_call"] == "v2*01"
    assert rec["receptor_revision_applied"] is False
    assert rec["original_v_call"] == ""


def test_noop_postfork_preserves_founder_trims_and_alignment():
    def run(with_post):
        exp = (
            ga.Experiment.on("human_igh").recombine()
            .trim(v_3=[(2, 1)], d_5=[(1, 1)], d_3=[(2, 1)], j_5=[(1, 1)])
            .clonal_repertoire(n_clones=1, max_size=1)
        )
        if with_post:
            exp.sequencing_errors(count=0)
        return exp.run_records(seed=42)[0]

    before, after = run(False), run(True)
    assert before["sequence"] == after["sequence"]
    expected_trims = {"v_trim_3": 2, "d_trim_5": 1, "d_trim_3": 2, "j_trim_5": 1}
    for field, value in expected_trims.items():
        assert before[field] == after[field] == value
    for field in (
        "v_cigar", "d_cigar", "j_cigar", "v_identity", "d_identity", "j_identity",
        "sequence_alignment", "germline_alignment",
    ):
        assert before[field] == after[field], field


@pytest.mark.parametrize("preset,has_d", [("HUMAN_IGK_OGRDB", False), ("HUMAN_IGH_OGRDB", True)])
def test_airr_np_fields_include_all_p_bases_between_alignments(preset, has_d):
    cfg = copy.deepcopy(getattr(ga, preset))
    lengths = {"V_3": 2, "J_5": 1}
    if has_d:
        lengths.update({"D_5": 1, "D_3": 2})
    cfg.reference_models = ReferenceEmpiricalModels(
        p_nucleotide_lengths={
            end: EmpiricalDistributionSpec([(length, 1)])
            for end, length in lengths.items()
        },
    )
    np_lengths = {"np1_lengths": [(3, 1)]}
    if has_d:
        np_lengths["np2_lengths"] = [(3, 1)]
    exp = ga.Experiment.on(cfg).recombine(**np_lengths).trim(enabled=False)
    result = exp.run_records(n=1, seed=0)
    rec = result[0]
    intervals = [(1, rec["v_sequence_end"], rec["d_sequence_start"] if has_d else rec["j_sequence_start"])]
    if has_d:
        intervals.append((2, rec["d_sequence_end"], rec["j_sequence_start"]))
    for index, start, end in intervals:
        assert rec[f"np{index}"] == rec["sequence"][start:end]
        assert rec[f"np{index}_length"] == end - start == 6
    assert result.validate_records(exp.refdata).ok


@pytest.mark.parametrize("preset,has_d", [("HUMAN_IGK_OGRDB", False), ("HUMAN_IGH_OGRDB", True)])
def test_noop_postfork_preserves_founder_p_addition_provenance(preset, has_d):
    cfg = copy.deepcopy(getattr(ga, preset))
    lengths = {"V_3": 2, "J_5": 1}
    if has_d:
        lengths.update({"D_5": 1, "D_3": 2})
    cfg.reference_models = ReferenceEmpiricalModels(
        p_nucleotide_lengths={
            end: EmpiricalDistributionSpec([(length, 1)])
            for end, length in lengths.items()
        },
    )

    def run(with_post):
        np_lengths = {"np1_lengths": [(0, 1)]}
        if has_d:
            np_lengths["np2_lengths"] = [(0, 1)]
        exp = (
            ga.Experiment.on(cfg).recombine(**np_lengths).trim(enabled=False)
            .clonal_repertoire(n_clones=1, max_size=1)
        )
        if with_post:
            exp.sequencing_errors(count=0)
        result = exp.run_records(seed=42)
        assert result.validate_records(exp.refdata).ok
        return result[0]

    before, after = run(False), run(True)
    assert before["sequence"] == after["sequence"]
    for field, expected in {
        "p_v_3_length": 2,
        "p_d_5_length": 1 if has_d else 0,
        "p_d_3_length": 2 if has_d else 0,
        "p_j_5_length": 1,
    }.items():
        assert before[field] == after[field] == expected
