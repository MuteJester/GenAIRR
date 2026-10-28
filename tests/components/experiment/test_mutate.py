"""The mutate() step: S5F and uniform somatic hypermutation.

Includes the per-base rate form, which draws a count per sequence rather than
taking a fixed one.
"""
from __future__ import annotations
import pytest
from GenAIRR import Experiment
import GenAIRR as ga  # noqa: E402
from tests.support.fixtures.experiment_dsl import _human_igh_exp

pytestmark = pytest.mark.behavior

def test_mutate_default_model_is_s5f():
    out = _human_igh_exp().mutate(count=10).run(n=1, seed=42)
    assert "mutate.s5f" in out[0].pass_names()
    assert out[0].trace().find("mutate.s5f.count").value == 10


def test_mutate_uniform_model_uses_uniform_pass():
    out = _human_igh_exp().mutate(model="uniform", count=12).run(n=1, seed=42)
    pass_names = out[0].pass_names()
    assert "mutate.uniform" in pass_names
    assert "mutate.s5f" not in pass_names
    assert out[0].trace().find("mutate.uniform.count").value == 12


def test_mutate_rejects_unknown_model():
    with pytest.raises(ValueError, match="model must be"):
        _human_igh_exp().mutate(model="unknown", count=5)


def test_mutate_count_fixed_int_produces_fixed_count():
    out = _human_igh_exp().mutate(count=7).run(n=5, seed=0)
    for o in out:
        assert o.trace().find("mutate.s5f.count").value == 7


def test_mutate_count_range_samples_in_range():
    # count=(5, 25) → uniform integer in [5, 25] inclusive.
    out = _human_igh_exp().mutate(count=(5, 25)).run(n=20, seed=0)
    sampled = sorted({o.trace().find("mutate.s5f.count").value for o in out})
    # All sampled counts must be in [5, 25].
    assert all(5 <= c <= 25 for c in sampled)
    # 20 seeds should hit at least 5 distinct values.
    assert len(sampled) >= 5


def test_mutate_count_empirical_distribution_respects_support():
    # count=[(0, 1.0), (50, 1.0)] → only 0 or 50 ever sampled.
    out = _human_igh_exp().mutate(count=[(0, 1.0), (50, 1.0)]).run(n=20, seed=0)
    sampled = {o.trace().find("mutate.s5f.count").value for o in out}
    assert sampled.issubset({0, 50})


def test_mutate_count_zero_is_a_noop():
    # count=0 → mutate pass runs but applies zero mutations.
    out = _human_igh_exp().mutate(count=0).run(n=1, seed=0)
    o = out[0]
    assert o.trace().find("mutate.s5f.count").value == 0
    assert o.trace().prefix_count("mutate.s5f.site") == 0


def test_mutate_rejects_negative_count():
    with pytest.raises(ValueError, match="non-negative"):
        _human_igh_exp().mutate(count=-3)


def test_mutate_rejects_inverted_range():
    with pytest.raises(ValueError, match="0 <= low <= high"):
        _human_igh_exp().mutate(count=(10, 5))


def test_mutate_rejects_negative_range_low():
    with pytest.raises(ValueError, match="0 <= low <= high"):
        _human_igh_exp().mutate(count=(-1, 5))


def test_mutate_rejects_empty_empirical_distribution():
    with pytest.raises(ValueError, match="at least one entry"):
        _human_igh_exp().mutate(count=[])


def test_mutate_rejects_bool_count():
    # bool is an int subclass — explicitly rejected so True/False
    # don't silently become count=1/0.
    with pytest.raises(TypeError, match="count: expected"):
        _human_igh_exp().mutate(count=True)


def test_mutate_rejects_float_count():
    with pytest.raises(TypeError, match="count: expected"):
        _human_igh_exp().mutate(count=3.5)  # type: ignore[arg-type]


def test_mutate_can_chain_after_recombine():
    # Pipeline shape: recombination first, then mutation. Pass list
    # ends with mutate.s5f.
    exp = Experiment.on("human_igh").recombine().mutate(count=5)
    assert exp.step_count == 2
    o = exp.run(n=1, seed=0)[0]
    pass_names = o.pass_names()
    assert pass_names[-1] == "mutate.s5f"
    assert pass_names[0].startswith("sample_allele")


def test_mutate_respects_productive_contract():
    # SHM under productive contract: the S5F base draws are filtered
    # so mutations don't introduce stop codons in the junction.
    out = _human_igh_exp().mutate(count=10).productive_only().run(n=3, seed=0)
    assert len(out) == 3


def test_mutate_strict_mode_under_productive_succeeds_for_typical_seeds():
    # Strict mode should succeed for the typical case — there's
    # always *some* admissible base under uniform substitution.
    out = (
        _human_igh_exp()
        .mutate(model="uniform", count=5)
        .productive_only().run(n=1, seed=0, strict=True)
    )
    assert len(out) == 1


def test_mutate_supports_alternate_s5f_models():
    # All four bundled S5F kernels load and run.
    for kernel in ("hh_s5f", "hh_s5f_60", "hh_s5f_opposite", "hkl_s5f"):
        out = _human_igh_exp().mutate(s5f_model=kernel, count=5).run(n=1, seed=0)
        assert out[0].trace().find("mutate.s5f.count").value == 5


def test_mutate_unknown_s5f_model_raises_at_run_time():
    # The error surfaces at compile-time when the step's apply()
    # tries to load the kernel (not at builder time, since we don't
    # eagerly validate).
    exp = _human_igh_exp().mutate(s5f_model="not_a_real_kernel", count=5)
    with pytest.raises(ValueError, match="Unknown S5F model"):
        exp.run(n=1, seed=0)


def test_mutate_is_deterministic_under_same_seed():
    a = _human_igh_exp().mutate(count=10).run(n=3, seed=0xCAFE)
    b = _human_igh_exp().mutate(count=10).run(n=3, seed=0xCAFE)
    for x, y in zip(a, b):
        assert x.final_simulation().bases() == y.final_simulation().bases()


def test_mutate_uniform_count_distribution_through_python_path():
    # Direct test: the (low, high) tuple in Python end up as the right
    # Empirical distribution at the engine layer. count=(3, 4) → only
    # 3 or 4 should ever be sampled.
    out = _human_igh_exp().mutate(model="uniform", count=(3, 4)).run(n=10, seed=0)
    sampled = {o.trace().find("mutate.uniform.count").value for o in out}
    assert sampled.issubset({3, 4})


def test_mutate_rate_runs_end_to_end_with_s5f_model():
    """rate=0.03 routes to the engine's S5F rate-mode constructor and
    runs successfully. The realized mutation count per record is
    drawn from Poisson(rate * pool_len), so we don't assert on
    exact counts — just that records come back and are productive
    when productive_only() is attached."""
    result = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(rate=0.03)
        .productive_only()
        .run_records(n=20, seed=42)
    )
    assert len(result.records) == 20
    for rec in result.records:
        assert rec["productive"] is True


def test_mutate_rate_runs_end_to_end_with_uniform_model():
    result = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(model="uniform", rate=0.05)
        .run_records(n=10, seed=42)
    )
    assert len(result.records) == 10


def test_mutate_rate_count_higher_with_higher_rate():
    """The realized mutation count should scale with rate. Compare
    a 1% rate vs 10% rate on the same seed/n; the 10% batch should
    have strictly more total mutations than the 1% batch."""
    def total_mutations(rate: float) -> int:
        result = (
            ga.Experiment.on("human_igh")
            .recombine()
            .mutate(rate=rate)
            .run_records(n=20, seed=0)
        )
        return sum(int(r.get("n_mutations") or 0) for r in result.records)

    low = total_mutations(0.01)
    high = total_mutations(0.10)
    assert high > low, (
        f"expected more mutations at 10% (got {high}) than 1% (got {low})"
    )


def test_mutate_rate_zero_produces_no_mutations():
    result = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(rate=0.0)
        .run_records(n=5, seed=0)
    )
    for rec in result.records:
        assert int(rec.get("n_mutations") or 0) == 0


def test_mutate_rate_and_count_mutually_exclusive():
    exp = ga.Experiment.on("human_igh").recombine()
    with pytest.raises(ValueError, match="exactly one of"):
        exp.mutate(count=8, rate=0.03)


def test_mutate_requires_count_or_rate():
    exp = ga.Experiment.on("human_igh").recombine()
    with pytest.raises(ValueError, match="exactly one of"):
        exp.mutate()


def test_mutate_rate_validates_range():
    exp = ga.Experiment.on("human_igh").recombine()
    with pytest.raises(ValueError, match=r"\[0\.0, 1\.0\]"):
        exp.mutate(rate=1.5)
    with pytest.raises(ValueError, match=r"\[0\.0, 1\.0\]"):
        exp.mutate(rate=-0.1)


def test_mutate_rate_rejects_non_numeric():
    exp = ga.Experiment.on("human_igh").recombine()
    with pytest.raises(TypeError, match="rate"):
        exp.mutate(rate="three percent")  # type: ignore[arg-type]


def test_mutate_rate_renders_in_describe_as_percent():
    out = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(rate=0.03)
        .describe()
    )
    assert "3.0% of bases" in out
    assert "Poisson" in out


def test_mutate_count_still_renders_in_describe_as_count():
    """Count-based path must keep its existing describe shape."""
    out = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(count=(5, 15))
        .describe()
    )
    assert "5–15 mutations" in out


def test_mutate_rate_is_deterministic_under_same_seed():
    """Two runs at the same seed produce the same mutation counts."""
    def n_muts(seed):
        result = (
            ga.Experiment.on("human_igh")
            .recombine()
            .mutate(rate=0.03)
            .run_records(n=10, seed=seed)
        )
        return [int(r.get("n_mutations") or 0) for r in result.records]
    assert n_muts(7) == n_muts(7)
