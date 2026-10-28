"""The corruption steps: PCR error, sequencing quality, indels and contamination.

These model what the protocol did to the molecule, not what the cell did.
"""
from __future__ import annotations
import pytest
from GenAIRR import CompiledExperiment, Experiment
import GenAIRR as ga  # noqa: E402
from tests.support.fixtures.experiment_dsl import _human_igh_exp, _vj_refdata

pytestmark = pytest.mark.behavior

def test_corrupt_pcr_appends_pcr_pass():
    o = _human_igh_exp().pcr_amplify(count=5).run(n=1, seed=0)[0]
    assert o.pass_names()[-1] == "corrupt.pcr"
    assert o.trace().find("corrupt.pcr.count").value == 5


def test_corrupt_pcr_count_zero_is_noop():
    o = _human_igh_exp().pcr_amplify(count=0).run(n=1, seed=0)[0]
    assert o.trace().find("corrupt.pcr.count").value == 0
    assert o.trace().prefix_count("corrupt.pcr.error_site") == 0


def test_corrupt_pcr_records_per_error_addresses():
    o = _human_igh_exp().pcr_amplify(count=3).run(n=1, seed=0)[0]
    assert o.trace().prefix_count("corrupt.pcr.error_site") == 3
    assert o.trace().prefix_count("corrupt.pcr.error_base") == 3


def test_corrupt_pcr_range_count_samples_in_range():
    out = _human_igh_exp().pcr_amplify(count=(1, 10)).run(n=20, seed=0)
    sampled = {o.trace().find("corrupt.pcr.count").value for o in out}
    assert all(1 <= c <= 10 for c in sampled)


def test_corrupt_pcr_rejects_negative_count():
    with pytest.raises(ValueError, match="non-negative"):
        _human_igh_exp().pcr_amplify(count=-1)


def test_corrupt_quality_appends_quality_pass():
    o = _human_igh_exp().sequencing_errors(count=4).run(n=1, seed=0)[0]
    assert o.pass_names()[-1] == "corrupt.quality"
    assert o.trace().find("corrupt.quality.count").value == 4


def test_corrupt_quality_writes_lowercase_bases():
    # Quality errors write lowercase; pick a count high enough that
    # the post-pass sequence contains at least one lowercase base.
    o = _human_igh_exp().sequencing_errors(count=20).run(n=1, seed=0)[0]
    bases = o.final_simulation().bases()
    has_lowercase = any(b in b"acgt" for b in bases)
    assert has_lowercase, "expected lowercase bases after quality errors"


def test_corrupt_quality_records_per_error_addresses():
    o = _human_igh_exp().sequencing_errors(count=2).run(n=1, seed=0)[0]
    assert o.trace().prefix_count("corrupt.quality.error_site") == 2
    assert o.trace().prefix_count("corrupt.quality.error_base") == 2


def test_corrupt_indels_appends_indel_pass():
    o = _human_igh_exp().polymerase_indels(count=2).run(n=1, seed=0)[0]
    assert o.pass_names()[-1] == "corrupt.indel"
    assert o.trace().find("corrupt.indel.count").value == 2


def test_corrupt_indels_default_insertion_prob_is_half():
    # Across 50 seeds, 100 events at p=0.5 should produce a mix of
    # insertions and deletions.
    out = _human_igh_exp().polymerase_indels(count=2).run(n=50, seed=0)
    saw_insertion = False
    saw_deletion = False
    for o in out:
        for i in range(2):
            kind = o.trace().find(f"corrupt.indel.kind[{i}]").value
            if kind:
                saw_insertion = True
            else:
                saw_deletion = True
    assert saw_insertion, "expected at least one insertion at default p=0.5"
    assert saw_deletion, "expected at least one deletion at default p=0.5"


def test_corrupt_indels_insertion_prob_one_is_all_insertions():
    out = _human_igh_exp().polymerase_indels(count=3, insertion_prob=1.0).run(n=5, seed=0)
    for o in out:
        for i in range(3):
            assert o.trace().find(f"corrupt.indel.kind[{i}]").value is True


def test_corrupt_indels_insertion_prob_zero_is_all_deletions():
    out = _human_igh_exp().polymerase_indels(count=3, insertion_prob=0.0).run(n=5, seed=0)
    for o in out:
        for i in range(3):
            assert o.trace().find(f"corrupt.indel.kind[{i}]").value is False


def test_corrupt_indels_changes_pool_length_relative_to_recombine_baseline():
    # Pool length after recombine + indels should differ from the
    # recombine-only baseline at the same seed because indels
    # insert/delete bases.
    base = _human_igh_exp().run(n=1, seed=0)[0].final_simulation()
    with_indels = (
        _human_igh_exp()
        .polymerase_indels(count=10, insertion_prob=1.0)
        .run(n=1, seed=0)[0]
        .final_simulation()
    )
    assert len(with_indels) == len(base) + 10


def test_corrupt_indels_rejects_negative_insertion_prob():
    with pytest.raises(ValueError, match=r"\[0\.0, 1\.0\]"):
        _human_igh_exp().polymerase_indels(count=1, insertion_prob=-0.1)


def test_corrupt_indels_rejects_above_one_insertion_prob():
    with pytest.raises(ValueError, match=r"\[0\.0, 1\.0\]"):
        _human_igh_exp().polymerase_indels(count=1, insertion_prob=1.5)


def test_corrupt_indels_rejects_nan_insertion_prob():
    with pytest.raises(ValueError, match=r"\[0\.0, 1\.0\]"):
        _human_igh_exp().polymerase_indels(count=1, insertion_prob=float("nan"))


def test_corrupt_contaminants_appends_contaminant_pass():
    o = _human_igh_exp().contaminate(prob=0.5).run(n=1, seed=0)[0]
    assert o.pass_names()[-1] == "corrupt.contaminant"
    assert o.trace().find("corrupt.contaminant.applied") is not None


def test_corrupt_contaminants_prob_one_always_applies():
    out = _human_igh_exp().contaminate(prob=1.0).run(n=10, seed=0)
    for o in out:
        assert o.trace().find("corrupt.contaminant.applied").value is True


def test_corrupt_contaminants_prob_zero_never_applies():
    out = _human_igh_exp().contaminate(prob=0.0).run(n=10, seed=0)
    for o in out:
        assert o.trace().find("corrupt.contaminant.applied").value is False


def test_corrupt_contaminants_rejects_negative_prob():
    with pytest.raises(ValueError, match=r"\[0\.0, 1\.0\]"):
        _human_igh_exp().contaminate(prob=-0.1)


def test_corrupt_contaminants_rejects_above_one_prob():
    with pytest.raises(ValueError, match=r"\[0\.0, 1\.0\]"):
        _human_igh_exp().contaminate(prob=1.5)


def test_corrupt_contaminants_rejects_nan_prob():
    with pytest.raises(ValueError, match=r"\[0\.0, 1\.0\]"):
        _human_igh_exp().contaminate(prob=float("nan"))


def test_full_pipeline_chains_all_four_corruption_steps_in_order():
    # Recombine → mutate → PCR → quality → indels → contaminants.
    # Pass names should appear in the order they were chained.
    exp = (
        _human_igh_exp()
        .mutate(count=5)
        .pcr_amplify(count=2)
        .sequencing_errors(count=1)
        .polymerase_indels(count=1)
        .contaminate(prob=0.0)
    )
    o = exp.run(n=1, seed=0)[0]
    pn = o.pass_names()
    assert pn[-5:] == [
        "mutate.s5f",
        "corrupt.pcr",
        "corrupt.quality",
        "corrupt.indel",
        "corrupt.contaminant",
    ]


def test_corruption_steps_are_deterministic_under_same_seed():
    exp_a = _human_igh_exp().pcr_amplify(count=3).polymerase_indels(count=2)
    exp_b = _human_igh_exp().pcr_amplify(count=3).polymerase_indels(count=2)
    a = exp_a.run(n=2, seed=0xCAFE)
    b = exp_b.run(n=2, seed=0xCAFE)
    for x, y in zip(a, b):
        assert x.final_simulation().bases() == y.final_simulation().bases()


def test_corruption_step_count_increments_step_count():
    exp = _human_igh_exp().pcr_amplify(count=1).sequencing_errors(count=1)
    assert exp.step_count == 3  # recombine + 2 corruption steps


def test_corruption_under_productive_contract_runs():
    # The PCR/contaminant base draws are filtered by `productive()`.
    # Quality and indel passes don't use contract-aware filtering
    # (intentionally); they just run and the engine relies on
    # post-hoc verify if anyone checks.
    exp = (
        _human_igh_exp()
        .pcr_amplify(count=3)
        .sequencing_errors(count=2)
        .polymerase_indels(count=1)
        .contaminate(prob=0.1)
    )
    out = exp.productive_only().run(n=3, seed=0)
    assert len(out) == 3


def test_pcr_amplify_rate_runs_end_to_end():
    """rate=1e-4 routes to push_corrupt_pcr_rate and runs to completion."""
    result = (
        ga.Experiment.on("human_igh")
        .recombine()
        .pcr_amplify(rate=1e-4)
        .run_records(n=10, seed=0)
    )
    assert len(result.records) == 10


def test_sequencing_errors_rate_runs_end_to_end():
    result = (
        ga.Experiment.on("human_igh")
        .recombine()
        .sequencing_errors(rate=1e-3)
        .run_records(n=10, seed=0)
    )
    assert len(result.records) == 10


def test_pcr_amplify_rate_zero_produces_no_errors():
    result = (
        ga.Experiment.on("human_igh")
        .recombine()
        .pcr_amplify(rate=0.0)
        .run_records(n=5, seed=0)
    )
    for rec in result.records:
        # n_pcr_errors should be 0 at rate=0
        assert int(rec.get("n_pcr_errors") or 0) == 0


def test_pcr_amplify_higher_rate_yields_more_errors():
    """Sanity scaling check: a 1% rate produces strictly more total
    PCR errors than a 0.001% rate on the same seed/n."""
    def total_errors(rate):
        result = (
            ga.Experiment.on("human_igh")
            .recombine()
            .pcr_amplify(rate=rate)
            .run_records(n=20, seed=0)
        )
        return sum(int(r.get("n_pcr_errors") or 0) for r in result.records)
    low = total_errors(1e-5)
    high = total_errors(0.01)
    assert high > low


def test_pcr_amplify_rate_and_count_mutually_exclusive():
    exp = ga.Experiment.on("human_igh").recombine()
    with pytest.raises(ValueError, match="exactly one of"):
        exp.pcr_amplify(count=3, rate=1e-4)


def test_sequencing_errors_rate_and_count_mutually_exclusive():
    exp = ga.Experiment.on("human_igh").recombine()
    with pytest.raises(ValueError, match="exactly one of"):
        exp.sequencing_errors(count=3, rate=1e-3)


def test_pcr_amplify_requires_count_or_rate():
    exp = ga.Experiment.on("human_igh").recombine()
    with pytest.raises(ValueError, match="exactly one of"):
        exp.pcr_amplify()


def test_sequencing_errors_requires_count_or_rate():
    exp = ga.Experiment.on("human_igh").recombine()
    with pytest.raises(ValueError, match="exactly one of"):
        exp.sequencing_errors()


def test_pcr_amplify_rate_validates_range():
    exp = ga.Experiment.on("human_igh").recombine()
    with pytest.raises(ValueError, match=r"\[0\.0, 1\.0\]"):
        exp.pcr_amplify(rate=1.5)


def test_pcr_amplify_rate_renders_in_describe_as_percent():
    out = (
        ga.Experiment.on("human_igh")
        .recombine()
        .pcr_amplify(rate=1e-4)
        .describe()
    )
    assert "0.01% of bases" in out


def test_sequencing_errors_rate_renders_in_describe_as_percent():
    out = (
        ga.Experiment.on("human_igh")
        .recombine()
        .sequencing_errors(rate=1e-3)
        .describe()
    )
    assert "0.10% of bases" in out


def test_describe_npa_source_tag_empirical_on_dataconfig():
    out = ga.Experiment.on("human_igh").recombine().describe()
    assert "empirical" in out  # NP source tag


def test_describe_npa_source_tag_uniform_fallback_on_raw_refdata():
    out = Experiment.on(_vj_refdata()).allow_curatable_refdata().recombine().trim(enabled=False).describe()
    assert "uniform fallback" in out
