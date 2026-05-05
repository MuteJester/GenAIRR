"""Tests for the user-facing :class:`GenAIRR.Experiment` DSL.

Exercises the Python-side fluent builder plus the ``CompiledExperiment``
it produces. Verifies:

- Construction from a Rust ``RefDataConfig`` (including the type
  guard on non-refdata inputs).
- Fluent chaining returns ``self`` so ``Experiment.on(...).recombine()``
  composes.
- ``recombine()`` compiles to the canonical 5-pass VJ shape and the
  canonical 8-pass VDJ shape.
- Custom NP-length distributions propagate through to the final
  trace.
- ``compile()`` is idempotent and the same plan can be re-run with
  different seeds.
- ``run(n, seed)`` produces ``n`` outcomes with deterministic per-run
  seeds (``seed + i``).
- Error paths: bad refdata type, ``n < 1``, empty length distribution.
"""
from __future__ import annotations

import genairr_engine as ge
import pytest

from GenAIRR import CompiledExperiment, Experiment


# ──────────────────────────────────────────────────────────────────
# Fixtures
# ──────────────────────────────────────────────────────────────────


def _vj_refdata() -> ge.RefDataConfig:
    cfg = ge.RefDataConfig.vj()
    cfg.add_v_allele("v1*01", "v1", b"AAACCCGGG", anchor=6)
    cfg.add_j_allele("j1*01", "j1", b"TTTAAA", anchor=0)
    return cfg


def _vdj_refdata() -> ge.RefDataConfig:
    cfg = ge.RefDataConfig.vdj()
    cfg.add_v_allele("v1*01", "v1", b"AAACCCGGG", anchor=6)
    cfg.add_d_allele("d1*01", "d1", b"TTTTTT")
    cfg.add_j_allele("j1*01", "j1", b"TTTAAA", anchor=0)
    return cfg


# ──────────────────────────────────────────────────────────────────
# Construction
# ──────────────────────────────────────────────────────────────────


def test_experiment_on_returns_experiment_instance():
    exp = Experiment.on(_vj_refdata())
    assert isinstance(exp, Experiment)
    assert exp.chain_type == "vj"
    assert exp.step_count == 0


def test_experiment_on_passes_through_vdj_chain_type():
    assert Experiment.on(_vdj_refdata()).chain_type == "vdj"


def test_experiment_on_rejects_unknown_config_name_string():
    # Strings are valid (config-name lookup), but unknown names get a
    # clear ValueError that lists alternatives.
    with pytest.raises(ValueError, match="Unknown config name"):
        Experiment.on("not_a_real_config_name")


def test_experiment_on_rejects_non_refdata_object():
    # Anything that is not a string / DataConfig / RefDataConfig is
    # rejected by the type guard with TypeError.
    with pytest.raises(TypeError, match="config-name string"):
        Experiment.on(42)


def test_experiment_on_rejects_none():
    with pytest.raises(TypeError, match="config-name string"):
        Experiment.on(None)


def test_experiment_on_resolves_known_config_name():
    # "human_igh" routes through the alias table to HUMAN_IGH_OGRDB
    # and produces a working VDJ Experiment with a non-empty refdata.
    exp = Experiment.on("human_igh")
    assert exp.chain_type == "vdj"
    assert exp.refdata.v_pool_size() > 0
    assert exp.refdata.j_pool_size() > 0


def test_experiment_on_accepts_dataconfig_object():
    # A loaded DataConfig (V5-style) is also a valid input.
    import GenAIRR

    cfg = GenAIRR.HUMAN_IGK_OGRDB
    exp = Experiment.on(cfg)
    assert exp.chain_type == "vj"
    assert exp.refdata.v_pool_size() > 0
    assert exp.refdata.j_pool_size() > 0


def test_experiment_repr_includes_chain_type_and_step_count():
    exp = Experiment.on(_vj_refdata()).recombine()
    r = repr(exp)
    assert "Experiment" in r
    assert "chain=vj" in r
    assert "steps=1" in r


# ──────────────────────────────────────────────────────────────────
# Fluent chaining
# ──────────────────────────────────────────────────────────────────


def test_recombine_returns_same_experiment_for_chaining():
    exp = Experiment.on(_vj_refdata())
    same = exp.recombine()
    assert same is exp
    assert exp.step_count == 1


def test_recombine_can_be_called_multiple_times():
    # Two recombine() calls = two recombination shapes in the plan.
    # Not biologically meaningful, but the builder permits it.
    exp = Experiment.on(_vj_refdata()).recombine().recombine()
    assert exp.step_count == 2


# ──────────────────────────────────────────────────────────────────
# compile() and the resulting CompiledExperiment
# ──────────────────────────────────────────────────────────────────


def test_compile_produces_compiled_experiment():
    compiled = Experiment.on(_vj_refdata()).recombine().compile()
    assert isinstance(compiled, CompiledExperiment)


def test_compile_vj_recombine_produces_five_pass_plan():
    compiled = Experiment.on(_vj_refdata()).recombine().compile()
    assert len(compiled.plan) == 5


def test_compile_vdj_recombine_produces_eight_pass_plan():
    compiled = Experiment.on(_vdj_refdata()).recombine().compile()
    assert len(compiled.plan) == 8


def test_compile_is_idempotent_and_returns_distinct_instances():
    exp = Experiment.on(_vj_refdata()).recombine()
    a = exp.compile()
    b = exp.compile()
    assert a is not b
    assert len(a.plan) == len(b.plan) == 5


def test_compiled_experiment_repr_includes_plan_len_and_chain():
    compiled = Experiment.on(_vj_refdata()).recombine().compile()
    r = repr(compiled)
    assert "CompiledExperiment" in r
    assert "plan_len=5" in r
    assert "chain=vj" in r


def test_compiled_experiment_exposes_plan_and_refdata():
    cfg = _vj_refdata()
    compiled = Experiment.on(cfg).recombine().compile()
    assert compiled.refdata is cfg
    assert isinstance(compiled.plan, ge.PassPlan)


# ──────────────────────────────────────────────────────────────────
# Pass-name shapes
# ──────────────────────────────────────────────────────────────────


def test_vj_recombine_emits_canonical_pass_names():
    outcomes = Experiment.on(_vj_refdata()).recombine().run(n=1, seed=0)
    assert outcomes[0].pass_names() == [
        "sample_allele.v",
        "sample_allele.j",
        "assemble.v",
        "generate_np.np1",
        "assemble.j",
    ]


def test_vdj_recombine_emits_canonical_pass_names():
    outcomes = Experiment.on(_vdj_refdata()).recombine().run(n=1, seed=0)
    assert outcomes[0].pass_names() == [
        "sample_allele.v",
        "sample_allele.d",
        "sample_allele.j",
        "assemble.v",
        "generate_np.np1",
        "assemble.d",
        "generate_np.np2",
        "assemble.j",
    ]


# ──────────────────────────────────────────────────────────────────
# Custom NP length distributions
# ──────────────────────────────────────────────────────────────────


def test_recombine_custom_np1_length_distribution_propagates():
    # Force NP1 = exactly 4 every time. The trace's NP1 length and
    # the assembled region's length must both be 4.
    cfg = _vj_refdata()
    outcomes = Experiment.on(cfg).recombine(np1_lengths=[(4, 1.0)]).run(
        n=10, seed=0
    )
    for o in outcomes:
        np1 = o.final_simulation().regions()[1]
        assert len(np1) == 4
        assert o.trace().find("np.np1.length").value == 4


def test_recombine_custom_np2_length_only_used_on_vdj():
    # Pass np2_lengths to a VJ experiment — silently ignored, no error.
    cfg = _vj_refdata()
    outcomes = Experiment.on(cfg).recombine(
        np1_lengths=[(2, 1.0)],
        np2_lengths=[(7, 1.0)],
    ).run(n=1, seed=0)
    np1 = outcomes[0].final_simulation().regions()[1]
    assert len(np1) == 2
    # Only NP1 region — VJ has no NP2.
    assert outcomes[0].final_simulation().region_count() == 3


def test_recombine_custom_np2_length_used_on_vdj():
    # VDJ experiment: NP2 length pinned to 5.
    cfg = _vdj_refdata()
    outcomes = Experiment.on(cfg).recombine(
        np1_lengths=[(2, 1.0)],
        np2_lengths=[(5, 1.0)],
    ).run(n=5, seed=0)
    for o in outcomes:
        regions = o.final_simulation().regions()
        # Regions: V, NP1, D, NP2, J
        np1, np2 = regions[1], regions[3]
        assert len(np1) == 2
        assert len(np2) == 5


def test_recombine_rejects_empty_length_distribution():
    cfg = _vj_refdata()
    with pytest.raises(ValueError, match="length distribution"):
        Experiment.on(cfg).recombine(np1_lengths=[])


# ──────────────────────────────────────────────────────────────────
# run(n, seed) — multi-iteration semantics
# ──────────────────────────────────────────────────────────────────


def test_run_with_n_one_returns_single_outcome():
    outcomes = Experiment.on(_vj_refdata()).recombine().run(n=1, seed=42)
    assert len(outcomes) == 1
    assert outcomes[0].pass_names()[0] == "sample_allele.v"


def test_run_returns_n_outcomes_with_offset_seeds():
    # Per-run seed is `seed + i`. Comparing to direct engine calls
    # at those exact seeds must produce byte-identical outcomes.
    cfg = _vj_refdata()
    outcomes = Experiment.on(cfg).recombine().run(n=4, seed=100)
    assert len(outcomes) == 4

    # Build the equivalent plan manually and compare.
    plan = ge.PassPlan()
    plan.push_sample_allele("V", cfg)
    plan.push_sample_allele("J", cfg)
    plan.push_assemble("V")
    plan.push_generate_np("NP1", [(i, 1.0) for i in range(7)])
    plan.push_assemble("J")

    for i, dsl_outcome in enumerate(outcomes):
        baseline = ge.run(plan, seed=100 + i, refdata=cfg)
        assert (
            dsl_outcome.final_simulation().bases()
            == baseline.final_simulation().bases()
        )


def test_run_default_seed_is_zero():
    a = Experiment.on(_vj_refdata()).recombine().run(n=1)
    b = Experiment.on(_vj_refdata()).recombine().run(n=1, seed=0)
    assert a[0].final_simulation().bases() == b[0].final_simulation().bases()


def test_run_default_n_is_one():
    outcomes = Experiment.on(_vj_refdata()).recombine().run()
    assert len(outcomes) == 1


def test_run_rejects_n_below_one():
    cfg = _vj_refdata()
    exp = Experiment.on(cfg).recombine()
    with pytest.raises(ValueError, match="n must be at least 1"):
        exp.run(n=0)
    with pytest.raises(ValueError, match="n must be at least 1"):
        exp.run(n=-3)


def test_run_is_deterministic_under_same_seed():
    a = Experiment.on(_vj_refdata()).recombine().run(n=10, seed=0xCAFE)
    b = Experiment.on(_vj_refdata()).recombine().run(n=10, seed=0xCAFE)
    for x, y in zip(a, b):
        assert x.final_simulation().bases() == y.final_simulation().bases()


def test_run_produces_distinct_outcomes_across_iterations():
    # 50 iterations × 7 NP-length values × randomized bases →
    # overwhelmingly likely that not all outcomes are byte-identical.
    cfg = _vj_refdata()
    cfg.add_v_allele("v2*01", "v2", b"GGGAAACCC", anchor=6)
    outcomes = Experiment.on(cfg).recombine().run(n=50, seed=0)
    base_strings = [o.final_simulation().bases() for o in outcomes]
    assert len(set(base_strings)) > 1


# ──────────────────────────────────────────────────────────────────
# CompiledExperiment.run can be called multiple times
# ──────────────────────────────────────────────────────────────────


def test_compiled_experiment_can_be_run_multiple_times():
    compiled = Experiment.on(_vj_refdata()).recombine().compile()
    a = compiled.run(n=2, seed=0)
    b = compiled.run(n=2, seed=10)
    # Different seeds → different outcomes. Same plan instance
    # (used by both run() calls) → no leftover state.
    assert a[0].final_simulation().bases() != b[0].final_simulation().bases()
    # And re-running with the original seed reproduces the original.
    again = compiled.run(n=2, seed=0)
    assert again[0].final_simulation().bases() == a[0].final_simulation().bases()


# ──────────────────────────────────────────────────────────────────
# F.6 — Contract bridge: productive() + Experiment.run(respect=...)
# ──────────────────────────────────────────────────────────────────


def test_productive_is_reexported_from_top_level():
    import GenAIRR

    assert hasattr(GenAIRR, "productive")
    assert callable(GenAIRR.productive)
    assert GenAIRR.productive is ge.productive


def test_productive_returns_contract_set_with_four_named_contracts():
    cs = ge.productive()
    assert isinstance(cs, ge.ContractSet)
    assert len(cs) == 4
    assert not cs.is_empty()
    names = cs.names()
    assert names == [
        "productive_junction_frame",
        "no_stop_codon_in_junction",
        "anchor_preserved.v",
        "anchor_preserved.j",
    ]


def test_contract_set_repr_includes_count_and_names():
    r = repr(ge.productive())
    assert "ContractSet" in r
    assert "len=4" in r
    assert "productive_junction_frame" in r


def test_run_with_respect_none_matches_run_without_respect():
    # respect=None is the no-op path — same outcome as no respect arg.
    cfg = _vj_refdata()
    exp = Experiment.on(cfg).recombine()
    a = exp.run(n=3, seed=0)
    b = exp.run(n=3, seed=0, respect=None)
    for x, y in zip(a, b):
        assert x.final_simulation().bases() == y.final_simulation().bases()


def test_run_with_productive_contract_constrains_np1_length_to_in_frame():
    # Synthetic VJ refdata: V_anchor_to_end = 3, J_anchor_to_W3 = 3,
    # so junction = 6 + NP1. Productive frame ⇒ NP1 % 3 == 0.
    # Filter NP1 length distribution to {0..6} uniform; with productive
    # active every NP1 length must be in {0, 3, 6}.
    cfg = _vj_refdata()
    outcomes = Experiment.on(cfg).recombine().run(
        n=20, seed=0, respect=ge.productive()
    )
    for o in outcomes:
        np1 = o.final_simulation().regions()[1]
        assert len(np1) % 3 == 0, (
            f"NP1 length {len(np1)} should be divisible by 3 under productive"
        )


def test_run_unconstrained_produces_some_out_of_frame_np1():
    # Negative control for the test above: without respect, the
    # uniform NP1 distribution produces out-of-frame lengths
    # ~67% of the time. Across 20 seeds we should see at least one.
    cfg = _vj_refdata()
    outcomes = Experiment.on(cfg).recombine().run(n=20, seed=0)
    out_of_frame = [
        o for o in outcomes if len(o.final_simulation().regions()[1]) % 3 != 0
    ]
    assert out_of_frame, "expected at least one out-of-frame NP1 without respect"


def test_run_accepts_respect_as_single_contract_set():
    cfg = _vj_refdata()
    out = Experiment.on(cfg).recombine().run(n=1, seed=0, respect=ge.productive())
    assert len(out) == 1


def test_run_accepts_respect_as_length_one_list():
    # V5 muscle-memory form — same effect as the bare ContractSet.
    cfg = _vj_refdata()
    a = Experiment.on(cfg).recombine().run(n=2, seed=0, respect=ge.productive())
    b = Experiment.on(cfg).recombine().run(n=2, seed=0, respect=[ge.productive()])
    for x, y in zip(a, b):
        assert x.final_simulation().bases() == y.final_simulation().bases()


def test_run_accepts_respect_as_empty_list_no_op():
    # Empty list normalises to None — no contracts active.
    cfg = _vj_refdata()
    a = Experiment.on(cfg).recombine().run(n=2, seed=0)
    b = Experiment.on(cfg).recombine().run(n=2, seed=0, respect=[])
    for x, y in zip(a, b):
        assert x.final_simulation().bases() == y.final_simulation().bases()


def test_run_rejects_non_contract_set_in_respect():
    cfg = _vj_refdata()
    with pytest.raises(TypeError, match="ContractSet"):
        Experiment.on(cfg).recombine().run(n=1, seed=0, respect="not a contract")


def test_run_rejects_non_contract_set_inside_respect_list():
    cfg = _vj_refdata()
    with pytest.raises(TypeError, match="respect\\[0\\]"):
        Experiment.on(cfg).recombine().run(n=1, seed=0, respect=[42])


def test_run_rejects_multi_bundle_respect_list():
    cfg = _vj_refdata()
    with pytest.raises(NotImplementedError, match="2 bundles"):
        Experiment.on(cfg).recombine().run(
            n=1, seed=0, respect=[ge.productive(), ge.productive()]
        )


def test_compiled_experiment_run_supports_respect():
    # The same respect= kwarg flows through CompiledExperiment.run.
    cfg = _vj_refdata()
    compiled = Experiment.on(cfg).recombine().compile()
    a = compiled.run(n=3, seed=0, respect=ge.productive())
    for o in a:
        assert len(o.final_simulation().regions()[1]) % 3 == 0


def test_run_with_respect_is_deterministic_under_same_seed():
    cfg = _vj_refdata()
    exp = Experiment.on(cfg).recombine()
    a = exp.run(n=5, seed=0xCAFE, respect=ge.productive())
    b = exp.run(n=5, seed=0xCAFE, respect=ge.productive())
    for x, y in zip(a, b):
        assert x.final_simulation().bases() == y.final_simulation().bases()


# ──────────────────────────────────────────────────────────────────
# F.7 — Strict-mode bridge: PassError → StrictSamplingError
# ──────────────────────────────────────────────────────────────────


def _vj_refdata_with_restrictive_np1(*, np1_lengths):
    """Build a VJ refdata + PassPlan whose NP1-length distribution is
    fully under our control. Used by the strict-mode tests so we can
    construct deliberately-unsatisfiable scenarios.
    """
    cfg = _vj_refdata()
    plan = ge.PassPlan()
    plan.push_sample_allele("V", cfg)
    plan.push_sample_allele("J", cfg)
    plan.push_assemble("V")
    plan.push_generate_np("NP1", np1_lengths)
    plan.push_assemble("J")
    return cfg, plan


def test_strict_sampling_error_is_exposed_at_top_level():
    import GenAIRR

    assert hasattr(GenAIRR, "StrictSamplingError")
    assert GenAIRR.StrictSamplingError is ge.StrictSamplingError
    assert issubclass(ge.StrictSamplingError, Exception)


def test_strict_with_satisfiable_contract_succeeds():
    cfg = _vj_refdata()
    out = (
        Experiment.on(cfg)
        .recombine()
        .run(n=3, seed=0, respect=ge.productive(), strict=True)
    )
    assert len(out) == 3


def test_strict_with_unsatisfiable_length_distribution_raises():
    # NP1 distribution {1, 2} → no in-frame lengths possible. Strict
    # mode must raise StrictSamplingError.
    cfg, plan = _vj_refdata_with_restrictive_np1(np1_lengths=[(1, 1.0), (2, 1.0)])
    with pytest.raises(ge.StrictSamplingError) as exc_info:
        ge.run(plan, seed=0, refdata=cfg, respect=ge.productive(), strict=True)

    pass_name, address, reason = exc_info.value.args
    assert pass_name == "generate_np.np1"
    assert address == "np.np1.length"
    assert reason == "empty_admissible_support"


def test_permissive_fallback_succeeds_where_strict_raises():
    # Same scenario as the strict raise above, but strict=False
    # accepts the unconstrained sample and lets the run continue.
    cfg, plan = _vj_refdata_with_restrictive_np1(np1_lengths=[(1, 1.0), (2, 1.0)])
    out = ge.run(plan, seed=0, refdata=cfg, respect=ge.productive(), strict=False)
    np1 = out.final_simulation().regions()[1]
    # Length must be one of the only-options from the input distribution.
    assert len(np1) in (1, 2)


def test_strict_mode_passes_through_compiled_experiment_run():
    cfg = _vj_refdata()
    compiled = Experiment.on(cfg).recombine().compile()
    # Satisfiable: succeeds.
    out = compiled.run(n=2, seed=0, respect=ge.productive(), strict=True)
    assert len(out) == 2


def test_strict_mode_default_is_false():
    # When strict is omitted, the run uses permissive sampling and
    # must succeed even on an unsatisfiable plan.
    cfg, plan = _vj_refdata_with_restrictive_np1(np1_lengths=[(1, 1.0)])
    out = ge.run(plan, seed=0, refdata=cfg, respect=ge.productive())
    assert out is not None


def test_strict_without_respect_runs_permissively():
    # strict=True with no contracts is a no-op for failure semantics
    # — there's nothing to fail against. The run proceeds normally.
    cfg = _vj_refdata()
    out = Experiment.on(cfg).recombine().run(n=1, seed=0, strict=True)
    assert len(out) == 1


def test_strict_sampling_error_args_form_three_tuple():
    # Public contract: StrictSamplingError.args is always a 3-tuple
    # so callers can destructure deterministically.
    cfg, plan = _vj_refdata_with_restrictive_np1(np1_lengths=[(1, 1.0)])
    try:
        ge.run(plan, seed=0, refdata=cfg, respect=ge.productive(), strict=True)
    except ge.StrictSamplingError as e:
        assert isinstance(e.args, tuple)
        assert len(e.args) == 3
        assert all(isinstance(a, str) for a in e.args)
    else:
        pytest.fail("expected StrictSamplingError")


def test_strict_mode_is_deterministic_under_same_seed():
    cfg, plan = _vj_refdata_with_restrictive_np1(
        np1_lengths=[(0, 1.0), (1, 1.0), (2, 1.0), (3, 1.0)]
    )
    # Both runs hit the satisfiable subset → both succeed and produce
    # identical outputs.
    a = ge.run(plan, seed=0xC0FFEE, refdata=cfg, respect=ge.productive(), strict=True)
    b = ge.run(plan, seed=0xC0FFEE, refdata=cfg, respect=ge.productive(), strict=True)
    assert a.final_simulation().bases() == b.final_simulation().bases()
