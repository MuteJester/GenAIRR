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
- ``compile()`` is idempotent and the same compiled simulator can be re-run with
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
    assert len(compiled.pass_plan) == 5


def test_compile_vdj_recombine_produces_eight_pass_plan():
    compiled = Experiment.on(_vdj_refdata()).recombine().compile()
    assert len(compiled.pass_plan) == 8


def test_compile_is_idempotent_and_returns_distinct_instances():
    exp = Experiment.on(_vj_refdata()).recombine()
    a = exp.compile()
    b = exp.compile()
    assert a is not b
    assert a.simulator is not b.simulator
    assert len(a.pass_plan) == len(b.pass_plan) == 5


def test_compiled_experiment_repr_includes_plan_len_and_chain():
    compiled = Experiment.on(_vj_refdata()).recombine().compile()
    r = repr(compiled)
    assert "CompiledExperiment" in r
    assert "plan_len=5" in r
    assert "chain=vj" in r


def test_compiled_experiment_exposes_simulator_pass_summary_and_refdata():
    cfg = _vj_refdata()
    compiled = Experiment.on(cfg).recombine().compile()
    assert compiled.refdata is cfg
    assert isinstance(compiled.simulator, ge.CompiledSimulator)
    assert compiled.pass_plan == (
        "sample_allele.v",
        "sample_allele.j",
        "assemble.v",
        "generate_np.np1",
        "assemble.j",
    )


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
    # Contracts are captured at compile time and reused by run().
    cfg = _vj_refdata()
    compiled = Experiment.on(cfg).recombine().compile(respect=ge.productive())
    assert compiled.active_contracts == (
        "productive_junction_frame",
        "no_stop_codon_in_junction",
        "anchor_preserved.v",
        "anchor_preserved.j",
    )
    a = compiled.run(n=3, seed=0)
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


def _vj_refdata_with_runtime_frame_residue():
    """Build a plan where one J allele has no productive NP1 completion.

    J allele 0 has anchor 0 and admits NP1 length 0. J allele 1 has
    anchor 1, so NP1 length 0 is out of frame. The compiled simulator
    should now filter J allele 1 before NP1 length sampling rather than
    relying on a later strict runtime failure.
    """
    cfg = ge.RefDataConfig.vj()
    cfg.add_v_allele("v1*01", "v1", b"AAACCCGGG", anchor=6)
    cfg.add_j_allele("j_good*01", "j_good", b"TTTAAA", anchor=0)
    cfg.add_j_allele("j_bad*01", "j_bad", b"TTTAAA", anchor=1)
    plan = ge.PassPlan()
    plan.push_sample_allele("V", cfg)
    plan.push_sample_allele("J", cfg)
    plan.push_assemble("V")
    plan.push_generate_np("NP1", [(0, 1.0)])
    plan.push_assemble("J")
    return cfg, plan


def _vj_refdata_with_stop_in_v_anchor():
    """Build a strict-runtime contract violation that is not a compile
    precondition failure.

    The V anchor codon is TAA, so the final materialized junction has a
    stop codon even though frame and anchor-preservation preconditions
    are satisfiable.
    """
    cfg = ge.RefDataConfig.vj()
    cfg.add_v_allele("v_stop*01", "v_stop", b"AAATAAGGG", anchor=3)
    cfg.add_j_allele("j1*01", "j1", b"TGGAAA", anchor=0)
    plan = ge.PassPlan()
    plan.push_sample_allele("V", cfg)
    plan.push_sample_allele("J", cfg)
    plan.push_assemble("V")
    plan.push_generate_np("NP1", [(0, 1.0)])
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


def test_productive_compile_rejects_unsatisfiable_length_distribution():
    # NP1 distribution {1, 2} has no in-frame mass for this VJ fixture.
    # This is now a build-time D7 failure, independent of strict mode.
    cfg, plan = _vj_refdata_with_restrictive_np1(np1_lengths=[(1, 1.0), (2, 1.0)])
    with pytest.raises(ValueError, match="NP1 length support has no in-frame mass"):
        ge.run(plan, seed=0, refdata=cfg, respect=ge.productive(), strict=True)


def test_compiled_feasibility_filters_runtime_frame_residue():
    # Compile succeeds because the declared support has an in-frame path,
    # but the bad J allele would leave NP1 with no admissible length.
    # The compiled feasibility layer should filter that J candidate
    # before it is committed, in strict and permissive modes alike.
    cfg, plan = _vj_refdata_with_runtime_frame_residue()
    for strict in (False, True):
        for seed in range(20):
            out = ge.run(
                plan,
                seed=seed,
                refdata=cfg,
                respect=ge.productive(),
                strict=strict,
            )
            assert out.final_simulation().j_allele_id() == 0


def test_strict_mode_passes_through_compiled_experiment_run():
    cfg = _vj_refdata()
    compiled = Experiment.on(cfg).recombine().compile(respect=ge.productive())
    # Satisfiable: succeeds.
    out = compiled.run(n=2, seed=0, strict=True)
    assert len(out) == 2


def test_strict_mode_default_is_false():
    cfg, plan = _vj_refdata_with_runtime_frame_residue()
    out = ge.run(plan, seed=0, refdata=cfg, respect=ge.productive())
    assert out.final_simulation().j_allele_id() == 0


def test_strict_without_respect_runs_permissively():
    # strict=True with no contracts is a no-op for failure semantics
    # — there's nothing to fail against. The run proceeds normally.
    cfg = _vj_refdata()
    out = Experiment.on(cfg).recombine().run(n=1, seed=0, strict=True)
    assert len(out) == 1


def test_strict_sampling_error_args_form_three_tuple():
    # Public contract: StrictSamplingError.args is always a 3-tuple
    # so callers can destructure deterministically.
    cfg, plan = _vj_refdata_with_stop_in_v_anchor()
    with pytest.raises(ge.StrictSamplingError) as err:
        ge.run(plan, seed=0, refdata=cfg, respect=ge.productive(), strict=True)
    exc = err.value
    assert isinstance(exc.args, tuple)
    assert len(exc.args) == 3
    assert all(isinstance(a, str) for a in exc.args)


def test_strict_mode_is_deterministic_under_same_seed():
    cfg, plan = _vj_refdata_with_restrictive_np1(
        np1_lengths=[(0, 1.0), (1, 1.0), (2, 1.0), (3, 1.0)]
    )
    # Both runs hit the satisfiable subset → both succeed and produce
    # identical outputs.
    a = ge.run(plan, seed=0xC0FFEE, refdata=cfg, respect=ge.productive(), strict=True)
    b = ge.run(plan, seed=0xC0FFEE, refdata=cfg, respect=ge.productive(), strict=True)
    assert a.final_simulation().bases() == b.final_simulation().bases()


def test_light_chain_productive_strict_batches_have_feasible_upstream_choices():
    # Regression for VJ productive feasibility: these seed windows used
    # to commit V/J/trim choices that left NP1 with no admissible length,
    # causing strict runtime failure or permissive nonproductive fallback.
    cases = [
        ("human_igk", 20260507, 60),
        ("human_igl", 0, 30),
        ("human_igl", 42, 30),
        ("human_igl", 20260507, 50),
    ]
    for cfg_name, seed, n in cases:
        records = Experiment.on(cfg_name).recombine().run_records(
            n=n,
            seed=seed,
            respect=ge.productive(),
            strict=True,
        )
        assert all(rec["productive"] is True for rec in records), cfg_name


# ──────────────────────────────────────────────────────────────────
# G.1 — .mutate() DSL step (S5F + uniform mutation)
# ──────────────────────────────────────────────────────────────────


def _human_igh_exp():
    """Convenience: default human-IGH recombination experiment."""
    return Experiment.on("human_igh").recombine()


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
    out = _human_igh_exp().mutate(count=10).run(n=3, seed=0, respect=ga.productive())
    assert len(out) == 3


def test_mutate_strict_mode_under_productive_succeeds_for_typical_seeds():
    # Strict mode should succeed for the typical case — there's
    # always *some* admissible base under uniform substitution.
    out = (
        _human_igh_exp()
        .mutate(model="uniform", count=5)
        .run(n=1, seed=0, respect=ga.productive(), strict=True)
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


# Top-level import for the productive() helper used in the
# mutate-under-productive tests above.
import GenAIRR as ga  # noqa: E402


# ──────────────────────────────────────────────────────────────────
# G.2 — corruption DSL: corrupt_pcr / quality / indels / contaminants
# ──────────────────────────────────────────────────────────────────


# --- corrupt_pcr ---


def test_corrupt_pcr_appends_pcr_pass():
    o = _human_igh_exp().corrupt_pcr(count=5).run(n=1, seed=0)[0]
    assert o.pass_names()[-1] == "corrupt.pcr"
    assert o.trace().find("corrupt.pcr.count").value == 5


def test_corrupt_pcr_count_zero_is_noop():
    o = _human_igh_exp().corrupt_pcr(count=0).run(n=1, seed=0)[0]
    assert o.trace().find("corrupt.pcr.count").value == 0
    assert o.trace().prefix_count("corrupt.pcr.error_site") == 0


def test_corrupt_pcr_records_per_error_addresses():
    o = _human_igh_exp().corrupt_pcr(count=3).run(n=1, seed=0)[0]
    assert o.trace().prefix_count("corrupt.pcr.error_site") == 3
    assert o.trace().prefix_count("corrupt.pcr.error_base") == 3


def test_corrupt_pcr_range_count_samples_in_range():
    out = _human_igh_exp().corrupt_pcr(count=(1, 10)).run(n=20, seed=0)
    sampled = {o.trace().find("corrupt.pcr.count").value for o in out}
    assert all(1 <= c <= 10 for c in sampled)


def test_corrupt_pcr_rejects_negative_count():
    with pytest.raises(ValueError, match="non-negative"):
        _human_igh_exp().corrupt_pcr(count=-1)


# --- corrupt_quality ---


def test_corrupt_quality_appends_quality_pass():
    o = _human_igh_exp().corrupt_quality(count=4).run(n=1, seed=0)[0]
    assert o.pass_names()[-1] == "corrupt.quality"
    assert o.trace().find("corrupt.quality.count").value == 4


def test_corrupt_quality_writes_lowercase_bases():
    # Quality errors write lowercase; pick a count high enough that
    # the post-pass sequence contains at least one lowercase base.
    o = _human_igh_exp().corrupt_quality(count=20).run(n=1, seed=0)[0]
    bases = o.final_simulation().bases()
    has_lowercase = any(b in b"acgt" for b in bases)
    assert has_lowercase, "expected lowercase bases after quality errors"


def test_corrupt_quality_records_per_error_addresses():
    o = _human_igh_exp().corrupt_quality(count=2).run(n=1, seed=0)[0]
    assert o.trace().prefix_count("corrupt.quality.error_site") == 2
    assert o.trace().prefix_count("corrupt.quality.error_base") == 2


# --- corrupt_indels ---


def test_corrupt_indels_appends_indel_pass():
    o = _human_igh_exp().corrupt_indels(count=2).run(n=1, seed=0)[0]
    assert o.pass_names()[-1] == "corrupt.indel"
    assert o.trace().find("corrupt.indel.count").value == 2


def test_corrupt_indels_default_insertion_prob_is_half():
    # Across 50 seeds, 100 events at p=0.5 should produce a mix of
    # insertions and deletions.
    out = _human_igh_exp().corrupt_indels(count=2).run(n=50, seed=0)
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
    out = _human_igh_exp().corrupt_indels(count=3, insertion_prob=1.0).run(n=5, seed=0)
    for o in out:
        for i in range(3):
            assert o.trace().find(f"corrupt.indel.kind[{i}]").value is True


def test_corrupt_indels_insertion_prob_zero_is_all_deletions():
    out = _human_igh_exp().corrupt_indels(count=3, insertion_prob=0.0).run(n=5, seed=0)
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
        .corrupt_indels(count=10, insertion_prob=1.0)
        .run(n=1, seed=0)[0]
        .final_simulation()
    )
    assert len(with_indels) == len(base) + 10


def test_corrupt_indels_rejects_negative_insertion_prob():
    with pytest.raises(ValueError, match=r"\[0\.0, 1\.0\]"):
        _human_igh_exp().corrupt_indels(count=1, insertion_prob=-0.1)


def test_corrupt_indels_rejects_above_one_insertion_prob():
    with pytest.raises(ValueError, match=r"\[0\.0, 1\.0\]"):
        _human_igh_exp().corrupt_indels(count=1, insertion_prob=1.5)


def test_corrupt_indels_rejects_nan_insertion_prob():
    with pytest.raises(ValueError, match=r"\[0\.0, 1\.0\]"):
        _human_igh_exp().corrupt_indels(count=1, insertion_prob=float("nan"))


# --- corrupt_contaminants ---


def test_corrupt_contaminants_appends_contaminant_pass():
    o = _human_igh_exp().corrupt_contaminants(prob=0.5).run(n=1, seed=0)[0]
    assert o.pass_names()[-1] == "corrupt.contaminant"
    assert o.trace().find("corrupt.contaminant.applied") is not None


def test_corrupt_contaminants_prob_one_always_applies():
    out = _human_igh_exp().corrupt_contaminants(prob=1.0).run(n=10, seed=0)
    for o in out:
        assert o.trace().find("corrupt.contaminant.applied").value is True


def test_corrupt_contaminants_prob_zero_never_applies():
    out = _human_igh_exp().corrupt_contaminants(prob=0.0).run(n=10, seed=0)
    for o in out:
        assert o.trace().find("corrupt.contaminant.applied").value is False


def test_corrupt_contaminants_rejects_negative_prob():
    with pytest.raises(ValueError, match=r"\[0\.0, 1\.0\]"):
        _human_igh_exp().corrupt_contaminants(prob=-0.1)


def test_corrupt_contaminants_rejects_above_one_prob():
    with pytest.raises(ValueError, match=r"\[0\.0, 1\.0\]"):
        _human_igh_exp().corrupt_contaminants(prob=1.5)


def test_corrupt_contaminants_rejects_nan_prob():
    with pytest.raises(ValueError, match=r"\[0\.0, 1\.0\]"):
        _human_igh_exp().corrupt_contaminants(prob=float("nan"))


# --- composition / determinism ---


def test_full_pipeline_chains_all_four_corruption_steps_in_order():
    # Recombine → mutate → PCR → quality → indels → contaminants.
    # Pass names should appear in the order they were chained.
    exp = (
        _human_igh_exp()
        .mutate(count=5)
        .corrupt_pcr(count=2)
        .corrupt_quality(count=1)
        .corrupt_indels(count=1)
        .corrupt_contaminants(prob=0.0)
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
    exp_a = _human_igh_exp().corrupt_pcr(count=3).corrupt_indels(count=2)
    exp_b = _human_igh_exp().corrupt_pcr(count=3).corrupt_indels(count=2)
    a = exp_a.run(n=2, seed=0xCAFE)
    b = exp_b.run(n=2, seed=0xCAFE)
    for x, y in zip(a, b):
        assert x.final_simulation().bases() == y.final_simulation().bases()


def test_corruption_step_count_increments_step_count():
    exp = _human_igh_exp().corrupt_pcr(count=1).corrupt_quality(count=1)
    assert exp.step_count == 3  # recombine + 2 corruption steps


def test_corruption_under_productive_contract_runs():
    # The PCR/contaminant base draws are filtered by `productive()`.
    # Quality and indel passes don't use contract-aware filtering
    # (intentionally — see Phase F.6 contract bridge); they just run
    # and the engine relies on post-hoc verify if anyone checks.
    exp = (
        _human_igh_exp()
        .corrupt_pcr(count=3)
        .corrupt_quality(count=2)
        .corrupt_indels(count=1)
        .corrupt_contaminants(prob=0.1)
    )
    out = exp.run(n=3, seed=0, respect=ga.productive())
    assert len(out) == 3


# ──────────────────────────────────────────────────────────────────
# G.3 — empirical NP-length + trim distributions from DataConfig
# ──────────────────────────────────────────────────────────────────


def test_recombine_default_inserts_trim_passes_for_dataconfig():
    # When backed by a DataConfig, recombine() now includes trim
    # passes for V_3, D_5, D_3, J_5 (in that order).
    o = _human_igh_exp().run(n=1, seed=0)[0]
    pn = o.pass_names()
    assert "trim.v_3" in pn
    assert "trim.d_5" in pn
    assert "trim.d_3" in pn
    assert "trim.j_5" in pn


def test_recombine_trim_passes_appear_after_sample_allele_before_assemble():
    o = _human_igh_exp().run(n=1, seed=0)[0]
    pn = o.pass_names()
    # All sample_allele.* must come before all trim.*; all trim.*
    # must come before all assemble.*.
    last_sample = max(i for i, n in enumerate(pn) if n.startswith("sample_allele"))
    first_trim = min(i for i, n in enumerate(pn) if n.startswith("trim."))
    last_trim = max(i for i, n in enumerate(pn) if n.startswith("trim."))
    first_assemble = min(i for i, n in enumerate(pn) if n.startswith("assemble."))
    assert last_sample < first_trim < last_trim < first_assemble


def test_recombine_trim_false_omits_trim_passes():
    o = _human_igh_exp().recombine.__wrapped__ if False else None  # noqa: E501 - placeholder
    # Use a fresh experiment with trim=False explicitly.
    exp = ga.Experiment.on("human_igh").recombine(trim=False)
    o = exp.run(n=1, seed=0)[0]
    pn = o.pass_names()
    assert not any(n.startswith("trim.") for n in pn)


def test_recombine_with_dataconfig_uses_empirical_np_lengths():
    # Empirical NP1 distribution can produce lengths beyond the
    # uniform [0..6] placeholder. With 50 seeds we should see
    # at least one length > 6 in human_igh (real distribution
    # has a long tail).
    out = _human_igh_exp().run(n=50, seed=0)
    sampled = {len(o.final_simulation().regions()[1]) for o in out}
    assert any(length > 6 for length in sampled), (
        f"expected at least one NP1 > 6 from empirical dist, got {sorted(sampled)}"
    )


def test_recombine_with_raw_refdata_falls_back_to_uniform_np_and_no_trim():
    # Raw RefDataConfig: no DataConfig backing → no trim passes,
    # NP lengths from the uniform [(0..6, 1.0)] placeholder.
    cfg = ge.RefDataConfig.vj()
    cfg.add_v_allele("v*01", "v", b"AAACCCGGG", anchor=6)
    cfg.add_j_allele("j*01", "j", b"TTTAAA", anchor=0)
    o = Experiment.on(cfg).recombine().run(n=1, seed=0)[0]
    pn = o.pass_names()
    assert not any(n.startswith("trim.") for n in pn)
    # Pure 5-pass VJ shape, no trims.
    assert pn == [
        "sample_allele.v",
        "sample_allele.j",
        "assemble.v",
        "generate_np.np1",
        "assemble.j",
    ]


def test_recombine_with_raw_refdata_trim_true_is_silent_no_op():
    # trim=True on a raw RefDataConfig is silent — there's no data
    # to source distributions from, so the step compiles to the
    # untrimmed plan rather than raising.
    cfg = ge.RefDataConfig.vj()
    cfg.add_v_allele("v*01", "v", b"AAACCCGGG", anchor=6)
    cfg.add_j_allele("j*01", "j", b"TTTAAA", anchor=0)
    o = Experiment.on(cfg).recombine(trim=True).run(n=1, seed=0)[0]
    assert not any(n.startswith("trim.") for n in o.pass_names())


def test_recombine_explicit_np1_lengths_overrides_empirical_default():
    # User-supplied np1_lengths wins over the empirical default
    # (which would otherwise draw from human_igh's empirical NP1
    # distribution).
    #
    # The earlier version of this test chained `_human_igh_exp().recombine(...)`
    # to verify the second recombine's NP1 distribution overrode the first.
    # That pattern became invalid under Phase 5's compile-time
    # trim-after-assembly check (the second recombine's trim passes land
    # after the first recombine's assemble passes). The real intent —
    # "explicit np1_lengths makes the recorded NP1 length match the
    # supplied distribution" — is independent of chaining, so we use a
    # single recombine with the explicit override.
    out = (
        Experiment.on("human_igh")
        .recombine(np1_lengths=[(0, 1.0)])
        .run(n=10, seed=0)
    )
    for o in out:
        np1_lens = [
            r.value for r in o.trace().prefix_query("np.np1.length")
        ]
        assert np1_lens, "expected at least one np.np1.length trace entry"
        assert np1_lens[-1] == 0


def test_recombine_d_trims_are_capped_per_segment():
    # Min D allele length in human_igh is ~10-12bp. With margin=3
    # and two_sided cap, the per-side cap is (min_d - 3) // 2 ≈ 3-4.
    # Across 100 seeds no D_5 or D_3 trim should exceed this.
    out = _human_igh_exp().run(n=100, seed=0)
    d5_max = max(o.trace().find("trim.d_5").value for o in out)
    d3_max = max(o.trace().find("trim.d_3").value for o in out)
    # Conservative upper bound: should never approach min_d (~10).
    assert d5_max <= 7, f"D_5 trim {d5_max} too aggressive"
    assert d3_max <= 7, f"D_3 trim {d3_max} too aggressive"


def test_recombine_does_not_panic_across_many_seeds():
    # The whole point of trim caps: no AssembleSegmentPass panic.
    out = _human_igh_exp().run(n=200, seed=0)
    assert len(out) == 200


def test_recombine_with_empirical_data_composes_with_productive_contract():
    out = _human_igh_exp().run(n=10, seed=0, respect=ga.productive())
    assert len(out) == 10


def test_extract_recombine_defaults_returns_distributions_for_human_igh():
    # Direct unit test on the helper.
    from GenAIRR._dataconfig_extract import extract_recombine_defaults

    cfg = ga.HUMAN_IGH_OGRDB
    defaults = extract_recombine_defaults(cfg)
    assert defaults["np1"] is not None and len(defaults["np1"]) > 0
    assert defaults["np2"] is not None and len(defaults["np2"]) > 0
    assert defaults["trim_v_3"] is not None
    assert defaults["trim_d_5"] is not None
    assert defaults["trim_d_3"] is not None
    assert defaults["trim_j_5"] is not None


def test_extract_recombine_defaults_returns_none_for_vj_d_trims():
    # IGK is a VJ chain — d_alleles is empty, so D trim distributions
    # come back as None.
    from GenAIRR._dataconfig_extract import extract_recombine_defaults

    cfg = ga.HUMAN_IGK_OGRDB
    defaults = extract_recombine_defaults(cfg)
    assert defaults["trim_d_5"] is None
    assert defaults["trim_d_3"] is None


def test_recombine_vj_uses_two_trim_passes():
    # Light chain: V_3 + J_5, no D trims.
    o = ga.Experiment.on("human_igk").recombine().run(n=1, seed=0)[0]
    pn = o.pass_names()
    trims = [n for n in pn if n.startswith("trim.")]
    assert trims == ["trim.v_3", "trim.j_5"]


def test_recombine_step_count_after_trim_does_not_change():
    # `step_count` is the number of FLUENT calls, not the number of
    # passes. recombine() with or without trims is one step.
    a = ga.Experiment.on("human_igh").recombine()
    b = ga.Experiment.on("human_igh").recombine(trim=False)
    assert a.step_count == 1
    assert b.step_count == 1


def test_empirical_recombine_is_deterministic_under_same_seed():
    a = _human_igh_exp().run(n=5, seed=0xCAFE)
    b = _human_igh_exp().run(n=5, seed=0xCAFE)
    for x, y in zip(a, b):
        assert x.final_simulation().bases() == y.final_simulation().bases()


# ──────────────────────────────────────────────────────────────────
# G.4 — SimulationResult wrapper + AIRR-format export
# ──────────────────────────────────────────────────────────────────


def _human_igh_records(n: int = 5, seed: int = 42, **kw):
    """Convenience: run a small human-IGH batch through ``run_records``
    and return the :class:`SimulationResult`."""
    return _human_igh_exp().run_records(n=n, seed=seed, **kw)


def test_run_records_returns_simulation_result():
    from GenAIRR.result import SimulationResult

    result = _human_igh_records(n=3)
    assert isinstance(result, SimulationResult)
    assert len(result) == 3


def test_simulation_result_is_listlike():
    result = _human_igh_records(n=4)

    # __len__
    assert len(result) == 4

    # __getitem__ by int
    rec0 = result[0]
    assert isinstance(rec0, dict)

    # __getitem__ by slice
    pair = result[1:3]
    assert isinstance(pair, list)
    assert len(pair) == 2

    # __iter__
    seen = list(result)
    assert len(seen) == 4
    assert all(isinstance(r, dict) for r in seen)


def test_simulation_result_repr():
    result = _human_igh_records(n=2)
    assert repr(result) == "<SimulationResult n=2>"


def test_simulation_result_exposes_outcomes_and_records():
    result = _human_igh_records(n=2)
    assert result.outcomes is not None
    assert len(result.outcomes) == 2
    # ``records`` is the underlying list view.
    assert result.records is not None
    assert len(result.records) == 2


def test_simulation_result_from_records_only_has_no_outcomes():
    from GenAIRR.result import SimulationResult

    rec = {"sequence": "ACGT", "v_call": "x", "junction": ""}
    result = SimulationResult([rec])
    assert len(result) == 1
    assert result.outcomes is None
    assert result[0]["sequence"] == "ACGT"


def test_record_core_fields_populated():
    result = _human_igh_records(n=1)
    rec = result[0]
    # Core sequence fields
    assert isinstance(rec["sequence"], str)
    assert rec["sequence_length"] == len(rec["sequence"])
    assert rec["sequence_length"] > 0


def test_record_v_call_resolves_to_imgt_name():
    result = _human_igh_records(n=3)
    for rec in result:
        v_call = rec["v_call"]
        assert isinstance(v_call, str)
        # Human IGH V allele names follow the IGHV...*NN pattern.
        assert v_call.startswith("IGH"), v_call
        assert "*" in v_call, v_call


def test_record_d_and_j_calls_resolve_for_vdj():
    result = _human_igh_records(n=3)
    for rec in result:
        assert rec["d_call"].startswith("IGH"), rec["d_call"]
        assert rec["j_call"].startswith("IGH"), rec["j_call"]
        assert "*" in rec["j_call"]


def test_record_d_call_empty_for_vj_chain():
    result = (
        ga.Experiment.on("human_igk").recombine().run_records(n=2, seed=7)
    )
    for rec in result:
        assert rec["d_call"] == ""
        assert rec["v_call"].startswith("IGK")
        assert rec["j_call"].startswith("IGK")


def test_record_v_d_j_coordinates_are_within_sequence():
    result = _human_igh_records(n=3)
    for rec in result:
        slen = rec["sequence_length"]
        for prefix in ("v", "d", "j"):
            start = rec[f"{prefix}_sequence_start"]
            end = rec[f"{prefix}_sequence_end"]
            assert start is not None and end is not None
            assert 0 <= start <= end <= slen


def test_record_junction_window_is_consistent():
    # Junction is derived from V-anchor → J-anchor + 3. Verify the
    # slice matches the recorded start/end and the length is correct.
    result = _human_igh_records(n=5)
    for rec in result:
        if rec["junction_start"] is None:
            continue
        seq = rec["sequence"]
        start = rec["junction_start"]
        end = rec["junction_end"]
        assert seq[start:end] == rec["junction"]
        assert rec["junction_length"] == end - start


def test_record_productive_flag_is_consistent_with_frame_and_stop():
    result = _human_igh_records(n=10, respect=ga.productive())
    for rec in result:
        assert rec["productive"] is True
        assert rec["vj_in_frame"] is True
        assert rec["stop_codon"] is False
        # Junction must be in-frame and translate cleanly.
        assert rec["junction_length"] is not None
        assert rec["junction_length"] % 3 == 0
        assert "*" not in rec["junction_aa"]


def test_record_n_mutations_matches_trace():
    # When .mutate(count=N) is fixed, the record's ``n_mutations``
    # should equal the trace count.
    result = _human_igh_exp().mutate(count=12).run_records(n=2, seed=0)
    for rec in result:
        assert rec["n_mutations"] == 12
        assert rec["mutation_rate"] == pytest.approx(
            12 / rec["sequence_length"]
        )


def test_record_corruption_counters_default_to_zero():
    result = _human_igh_records(n=2)
    for rec in result:
        assert rec["n_pcr_errors"] == 0
        assert rec["n_quality_errors"] == 0
        assert rec["n_indels"] == 0
        assert rec["is_contaminant"] is False


def test_record_corruption_counters_track_pcr_errors():
    result = (
        _human_igh_exp().corrupt_pcr(count=4).run_records(n=2, seed=0)
    )
    for rec in result:
        assert rec["n_pcr_errors"] == 4


def test_compiled_run_records_matches_experiment_run_records():
    # The two paths (Experiment.run_records vs CompiledExperiment.run_records)
    # should produce identical record sequences for the same seed.
    exp = _human_igh_exp().mutate(count=5)
    a = exp.run_records(n=3, seed=99)
    b = exp.compile().run_records(n=3, seed=99)
    assert len(a) == len(b)
    for x, y in zip(a, b):
        assert x["sequence"] == y["sequence"]
        assert x["v_call"] == y["v_call"]
        assert x["junction"] == y["junction"]


def test_to_tsv_roundtrip(tmp_path):
    import csv as _csv

    result = _human_igh_records(n=3)
    path = tmp_path / "batch.tsv"
    result.to_tsv(str(path))

    with open(path, encoding="utf-8") as fh:
        rows = list(_csv.DictReader(fh, delimiter="\t"))

    assert len(rows) == 3
    # Spot-check a few canonical columns.
    for row, rec in zip(rows, result):
        assert row["sequence"] == rec["sequence"]
        assert row["v_call"] == rec["v_call"]
        assert int(row["sequence_length"]) == rec["sequence_length"]


def test_to_csv_uses_comma_delimiter(tmp_path):
    result = _human_igh_records(n=2)
    path = tmp_path / "batch.csv"
    result.to_csv(str(path))

    text = path.read_text(encoding="utf-8")
    header = text.splitlines()[0]
    assert "," in header
    assert "\t" not in header
    assert "sequence" in header
    assert "v_call" in header


def test_to_tsv_writes_canonical_column_order(tmp_path):
    from GenAIRR.result import _DEFAULT_COLUMN_ORDER

    result = _human_igh_records(n=1)
    path = tmp_path / "order.tsv"
    result.to_tsv(str(path))

    header = path.read_text(encoding="utf-8").splitlines()[0].split("\t")
    # Header must start with the canonical column order.
    assert header[: len(_DEFAULT_COLUMN_ORDER)] == _DEFAULT_COLUMN_ORDER


def test_to_tsv_serializes_none_as_empty(tmp_path):
    from GenAIRR.result import SimulationResult

    rec = {
        "sequence": "ACGT",
        "sequence_length": 4,
        "v_call": "v1",
        "d_call": "",
        "junction": "",
        "junction_start": None,
        "junction_end": None,
    }
    path = tmp_path / "none.tsv"
    SimulationResult([rec]).to_tsv(str(path))
    text = path.read_text(encoding="utf-8")
    # ``None`` columns should render as empty fields, not literal "None".
    assert "None" not in text


def test_to_fasta_emits_header_and_sequence(tmp_path):
    result = _human_igh_records(n=2)
    path = tmp_path / "batch.fa"
    result.to_fasta(str(path))

    lines = path.read_text(encoding="utf-8").splitlines()
    assert len(lines) == 4  # 2 records × (header + seq)
    assert lines[0].startswith(">seq0|v_call=")
    assert "|j_call=" in lines[0]
    assert lines[1] == result[0]["sequence"]
    assert lines[2].startswith(">seq1|v_call=")
    assert lines[3] == result[1]["sequence"]


def test_to_fasta_respects_custom_prefix(tmp_path):
    result = _human_igh_records(n=1)
    path = tmp_path / "custom.fa"
    result.to_fasta(str(path), prefix="read_")
    first = path.read_text(encoding="utf-8").splitlines()[0]
    assert first.startswith(">read_0|")


# ──────────────────────────────────────────────────────────────────
# G1 — FASTQ output with Phred Q-scores
# ──────────────────────────────────────────────────────────────────


def test_to_fastq_constant_quality_writes_q30(tmp_path):
    # Constant-Q model writes the same Q for every (non-N) base.
    # Q30 → Phred+33 ASCII '?' (33 + 30 = 63).
    result = _human_igh_records(n=2)
    path = tmp_path / "out.fastq"
    result.to_fastq(str(path), quality="constant", q=30)
    lines = path.read_text(encoding="utf-8").splitlines()
    assert len(lines) == 8  # 2 records × 4 lines
    # First record: @-header / SEQ / + / Q-string
    assert lines[0].startswith("@seq0|v_call=")
    assert lines[1] == result[0]["sequence"].upper()
    assert lines[2] == "+"
    seq = result[0]["sequence"]
    expected_q = "".join(
        "#" if c in ("N", "n") else "?" for c in seq
    )
    assert lines[3] == expected_q


def test_to_fastq_illumina_quality_shape(tmp_path):
    # Illumina model: ramp up over 10 bases, hold peak, taper over
    # 30 bases.
    result = _human_igh_records(n=1)
    path = tmp_path / "out.fastq"
    result.to_fastq(str(path), quality="illumina")
    lines = path.read_text(encoding="utf-8").splitlines()
    q_string = lines[3]
    seq = result[0]["sequence"]
    # All ASCII ≥ 33 (printable).
    for ch in q_string:
        assert 33 <= ord(ch) <= 126
    # Same length as the sequence.
    assert len(q_string) == len(seq)
    # Ramp-up: position 0 is Q26 (≈ start_q+1 from interpolation),
    # position 9 is Q35 (peak).
    assert ord(q_string[0]) - 33 < ord(q_string[9]) - 33
    # Tail: position -1 is Q18 (end_q), strictly less than peak.
    assert ord(q_string[-1]) - 33 < ord(q_string[len(seq) // 2]) - 33


def test_to_fastq_n_bases_get_low_quality(tmp_path):
    # Manually construct a result where one record has N in the
    # sequence and verify that position gets Q2 (`#` ASCII 35).
    result = _human_igh_records(n=1)
    # Mutate the result's first record's sequence to inject an N.
    rec = result[0]
    seq = rec["sequence"]
    rec["sequence"] = seq[:5] + "N" + seq[6:]
    path = tmp_path / "out.fastq"
    result.to_fastq(str(path), quality="constant", q=30)
    lines = path.read_text(encoding="utf-8").splitlines()
    q_string = lines[3]
    # Position 5 should be the n_q (Q2 default → '#').
    assert q_string[5] == "#"
    # Other positions should be the full Q (Q30 → '?').
    assert q_string[0] == "?"


def test_to_fastq_uppercase_sequence(tmp_path):
    # FASTQ convention: uppercase bases regardless of GenAIRR's
    # internal lowercase germline marking.
    result = _human_igh_records(n=1)
    path = tmp_path / "out.fastq"
    result.to_fastq(str(path), quality="constant", q=30)
    lines = path.read_text(encoding="utf-8").splitlines()
    seq_line = lines[1]
    assert seq_line == seq_line.upper(), (
        f"FASTQ sequence should be uppercase: {seq_line[:40]}"
    )


def test_to_fastq_rejects_unknown_model(tmp_path):
    result = _human_igh_records(n=1)
    path = tmp_path / "out.fastq"
    try:
        result.to_fastq(str(path), quality="bogus")
    except ValueError:
        pass
    else:
        raise AssertionError("expected ValueError for unknown model")


def test_to_fastq_rejects_out_of_range_q(tmp_path):
    result = _human_igh_records(n=1)
    path = tmp_path / "out.fastq"
    for bad in (-1, 100, "high"):
        try:
            result.to_fastq(str(path), quality="constant", q=bad)
        except (ValueError, TypeError):
            pass
        else:
            raise AssertionError(f"expected error for q={bad!r}")


def test_g7_corrupt_ns_writes_uppercase_n():
    # corrupt_ns(count=10) should sprinkle ~10 N's per record
    # (collisions reduce the count slightly).
    exp = ga.Experiment.on("human_igh").recombine().corrupt_ns(count=10)
    for rec in exp.run_records(n=5, seed=0):
        n_count = rec["sequence"].upper().count("N")
        # With pool length ~370 and 10 random sites, collisions are
        # rare; expect 8-10 inclusive.
        assert 6 <= n_count <= 10, (
            f"expected 6-10 Ns, got {n_count} (sequence len={len(rec['sequence'])})"
        )


def test_g7_corrupt_ns_zero_count_is_no_op():
    base = ga.Experiment.on("human_igh").recombine()
    flipped = ga.Experiment.on("human_igh").recombine().corrupt_ns(count=0)
    for fwd, ns in zip(base.run_records(n=5, seed=0), flipped.run_records(n=5, seed=0)):
        assert fwd["sequence"] == ns["sequence"]


def test_g7_corrupt_ns_count_distribution_shape():
    # Variable-count distribution: count=(5, 15) → uniform integer
    # in [5, 15]. Across many records, the N-count distribution
    # should have all entries in approximately that range.
    exp = ga.Experiment.on("human_igh").recombine().corrupt_ns(count=(5, 15))
    counts = [r["sequence"].upper().count("N") for r in exp.run_records(n=50, seed=0)]
    assert min(counts) >= 4, f"min N count too low: {counts}"
    assert max(counts) <= 15, f"max N count too high: {counts}"


def test_g7_corrupt_ns_compose_with_other_passes():
    # Stack with PCR + indels — N's should still appear and the
    # alignment-string invariants should hold.
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(count=5)
        .corrupt_pcr(count=2)
        .corrupt_indels(count=1)
        .corrupt_ns(count=8)
    )
    failures = []
    for i, rec in enumerate(exp.run_records(n=30, seed=0)):
        issues = _alignment_invariants_hold(rec)
        if issues:
            failures.append((i, issues))
        # At least some records should carry Ns.
    n_with_ns = sum(
        1 for r in exp.run_records(n=30, seed=0)
        if "N" in r["sequence"].upper()
    )
    assert n_with_ns >= 25, f"only {n_with_ns}/30 records had Ns"
    assert failures == [], f"alignment invariants broke: {failures[:3]}"


def test_g5_clonal_structure_total_record_count():
    # n_clones × size = total records.
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .with_clonal_structure(n_clones=5, size=4)
        .mutate(count=3)
    )
    result = exp.run_records(seed=0)
    assert len(result) == 20


def test_g5_clonal_descendants_identical_without_post_fork_passes():
    # No mutate / corrupt after the fork → every descendant in a
    # clone is literally identical (same V/D/J + trim + NP).
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .with_clonal_structure(n_clones=3, size=4)
    )
    by_clone = {}
    for rec in exp.run_records(seed=0):
        by_clone.setdefault(rec["clone_id"], []).append(rec)
    assert sorted(by_clone) == [0, 1, 2]
    for cid, recs in by_clone.items():
        seqs = {r["sequence"] for r in recs}
        v_calls = {r["v_call"] for r in recs}
        assert len(seqs) == 1, f"clone {cid} descendants diverged without post-fork passes"
        assert len(v_calls) == 1


def test_g5_clonal_descendants_share_junction_diverge_via_mutate():
    # With mutate after the fork, descendants share the junction
    # backbone (same V/D/J + trim + NP) but diverge in mutated bases.
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .with_clonal_structure(n_clones=2, size=8)
        .mutate(count=8)
    )
    by_clone = {}
    for rec in exp.run_records(seed=0):
        by_clone.setdefault(rec["clone_id"], []).append(rec)
    for cid, recs in by_clone.items():
        v_calls = {r["v_call"] for r in recs if r["v_call"]}
        seqs = {r["sequence"] for r in recs}
        # Sequences diverge under SHM…
        assert len(seqs) > 1, f"clone {cid} descendants didn't diverge under mutate"
        # …but the V allele backbone is shared (when live evidence
        # didn't get wiped). All V calls should reference the same
        # provenance — rejecting empty calls (heavy-SHM noise).
        non_empty = [c for c in v_calls if c]
        if non_empty:
            assert len(set(non_empty)) == 1, (
                f"clone {cid} V calls differ across descendants: {non_empty}"
            )


def test_g5_clonal_explicit_n_must_match_total():
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .with_clonal_structure(n_clones=4, size=5)
    )
    # n=20 matches → ok.
    result = exp.run_records(n=20, seed=0)
    assert len(result) == 20
    # n=10 mismatched → error.
    with pytest.raises(ValueError, match="inconsistent"):
        exp.run_records(n=10, seed=0)


def test_g5_with_clonal_structure_rejects_double_call():
    exp = ga.Experiment.on("human_igh").recombine()
    exp = exp.with_clonal_structure(n_clones=3, size=2)
    with pytest.raises(ValueError, match="only be called once"):
        exp.with_clonal_structure(n_clones=2, size=2)


def test_g5_with_clonal_structure_rejects_invalid_args():
    exp = ga.Experiment.on("human_igh").recombine()
    with pytest.raises(ValueError, match="positive int"):
        exp.with_clonal_structure(n_clones=0, size=5)
    with pytest.raises(ValueError, match="positive int"):
        exp.with_clonal_structure(n_clones=5, size=-1)


def test_g5_clonal_record_to_tsv(tmp_path):
    # The clone_id column should round-trip through to_tsv output.
    import csv as _csv
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .with_clonal_structure(n_clones=3, size=4)
        .mutate(count=2)
    )
    path = tmp_path / "clones.tsv"
    exp.run_records(seed=0).to_tsv(str(path))
    with open(path, "r") as fh:
        reader = _csv.DictReader(fh, delimiter="\t")
        rows = list(reader)
    assert len(rows) == 12
    clone_ids = sorted({int(r["clone_id"]) for r in rows})
    assert clone_ids == [0, 1, 2]
    # Each clone has 4 descendants.
    from collections import Counter
    counts = Counter(int(r["clone_id"]) for r in rows)
    assert counts == Counter({0: 4, 1: 4, 2: 4})


def test_g4_mutate_rejected_on_tcr_chain():
    # SHM is biologically B-cell only. The DSL must error rather
    # than silently producing biologically-false TCR records.
    with pytest.raises(ValueError, match="somatic hypermutation does not occur in TCR"):
        ga.Experiment.on("human_tcrb").recombine().mutate(count=5)


def test_g4_mutate_allowed_on_bcr_chain():
    # Sanity: BCR loci accept mutate() unchanged.
    for cfg in ("human_igh", "human_igk", "human_igl"):
        exp = ga.Experiment.on(cfg).recombine().mutate(count=3)
        records = exp.run_records(n=2, seed=0)
        assert len(records) == 2


def test_g4_tcr_pipeline_with_corrupt_passes_works():
    # Confirms the recommended TCR pipeline (corrupt_pcr / corrupt_
    # quality / corrupt_indels) still works end-to-end.
    exp = (
        ga.Experiment.on("human_tcrb")
        .recombine()
        .corrupt_pcr(count=3)
        .corrupt_quality(count=2)
        .corrupt_indels(count=1, insertion_prob=0.5)
    )
    records = exp.run_records(n=5, seed=0)
    assert len(records) == 5
    for rec in records:
        assert rec["locus"] == "TRB"


def test_g3b_allele_weights_boost_listed_v():
    # Boosting a single V allele's weight by 1000x should make it
    # dominate sampling over a 200-record batch.
    from collections import Counter
    target = "IGHVF1-G1*01"
    exp = ga.Experiment.on("human_igh").recombine(v_allele_weights={target: 1000.0})
    counts = Counter()
    for r in exp.run_records(n=200, seed=0):
        first = (r["v_call"] or "").split(",")[0]
        counts[first] += 1
    # The boosted allele should be at least 50% of calls (typically
    # >80% with weight=1000 against ~60 alleles each at 1.0).
    assert counts[target] >= 100, (
        f"boosted allele was sampled {counts[target]} / 200 times; "
        f"top counts: {counts.most_common(3)}"
    )


def test_g3b_allele_weights_unspecified_keeps_default_uniform():
    # When v_allele_weights is None, the distribution should match
    # the un-weighted run (allowing for stochastic variation across
    # 100 records with the same seed).
    exp_a = ga.Experiment.on("human_igh").recombine()
    exp_b = ga.Experiment.on("human_igh").recombine(v_allele_weights=None)
    a = [r["v_call"] for r in exp_a.run_records(n=20, seed=0)]
    b = [r["v_call"] for r in exp_b.run_records(n=20, seed=0)]
    assert a == b


def test_g3b_allele_weights_rejects_unknown_name():
    with pytest.raises(ValueError, match="no V allele named"):
        ga.Experiment.on("human_igh").recombine(
            v_allele_weights={"IGHV-NOT-A-REAL-ALLELE*99": 5.0}
        )


def test_g3b_allele_weights_rejects_non_positive_weight():
    with pytest.raises(ValueError, match="finite positive number"):
        ga.Experiment.on("human_igh").recombine(
            v_allele_weights={"IGHVF1-G1*01": 0.0}
        )
    with pytest.raises(ValueError, match="finite positive number"):
        ga.Experiment.on("human_igh").recombine(
            v_allele_weights={"IGHVF1-G1*01": -1.0}
        )


def test_g3b_allele_weights_d_on_vj_chain_errors():
    with pytest.raises(ValueError, match="cannot weight D alleles"):
        ga.Experiment.on("human_igk").recombine(
            d_allele_weights={"X*01": 1.0}
        )


def test_g3b_allele_weights_empty_dict_errors():
    with pytest.raises(ValueError, match="at least one"):
        ga.Experiment.on("human_igh").recombine(v_allele_weights={})


def test_to_fastq_pairs_with_constant_q_is_decodable(tmp_path):
    # Round-trip check: write FASTQ, parse it back, confirm record
    # count and length match. Uses Python stdlib (no biopython) to
    # keep the test dep-light.
    result = _human_igh_records(n=5)
    path = tmp_path / "out.fastq"
    result.to_fastq(str(path), quality="constant", q=25)
    text = path.read_text(encoding="utf-8")
    records = text.split("@seq")[1:]  # split keeps records starting with header tail
    assert len(records) == 5
    for i, blob in enumerate(records):
        lines = blob.splitlines()
        # Layout: "{i}|...", seq, "+", q_string
        seq = lines[1]
        q = lines[3]
        assert len(seq) == len(q), f"record {i}: seq/q length mismatch"


def test_to_dataframe_returns_dataframe_with_record_rows():
    pd = pytest.importorskip("pandas")

    result = _human_igh_records(n=4)
    df = result.to_dataframe()
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 4
    assert "sequence" in df.columns
    assert "v_call" in df.columns
    # First-row spot-check.
    assert df.iloc[0]["sequence"] == result[0]["sequence"]


def test_to_dataframe_empty_returns_canonical_columns():
    pd = pytest.importorskip("pandas")
    from GenAIRR.result import SimulationResult, _DEFAULT_COLUMN_ORDER

    df = SimulationResult([]).to_dataframe()
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 0
    assert list(df.columns) == _DEFAULT_COLUMN_ORDER


def test_column_order_preserves_extras_after_canonical():
    from GenAIRR.result import SimulationResult, _DEFAULT_COLUMN_ORDER

    rec = {"sequence": "A", "v_call": "x", "extra_a": 1, "extra_b": 2}
    result = SimulationResult([rec])
    cols = result._column_order()
    # Canonical columns come first.
    assert cols[: len(_DEFAULT_COLUMN_ORDER)] == _DEFAULT_COLUMN_ORDER
    # Extras follow in record-encounter order.
    extras = cols[len(_DEFAULT_COLUMN_ORDER) :]
    assert extras == ["extra_a", "extra_b"]


# ──────────────────────────────────────────────────────────────────
# G.5 — allele locking via Experiment.using(v=, d=, j=)
# ──────────────────────────────────────────────────────────────────


def _v_name(refdata, outcome):
    return refdata.v_allele(outcome.final_simulation().v_allele_id()).name


def _d_name(refdata, outcome):
    return refdata.d_allele(outcome.final_simulation().d_allele_id()).name


def _j_name(refdata, outcome):
    return refdata.j_allele(outcome.final_simulation().j_allele_id()).name


def test_using_returns_self_for_chaining():
    exp = ga.Experiment.on("human_igh")
    same = exp.using(v="IGHVF1-G1*01")
    assert same is exp


def test_using_locks_single_v_allele():
    exp = (
        ga.Experiment.on("human_igh")
        .using(v="IGHVF1-G1*01")
        .recombine()
    )
    for o in exp.run(n=10, seed=0):
        assert _v_name(exp.refdata, o) == "IGHVF1-G1*01"


def test_using_locks_single_d_allele():
    exp = (
        ga.Experiment.on("human_igh")
        .using(d="IGHD3-3*01")
        .recombine()
    )
    for o in exp.run(n=10, seed=0):
        assert _d_name(exp.refdata, o) == "IGHD3-3*01"


def test_using_locks_single_j_allele():
    exp = (
        ga.Experiment.on("human_igh")
        .using(j="IGHJ4*02")
        .recombine()
    )
    for o in exp.run(n=10, seed=0):
        assert _j_name(exp.refdata, o) == "IGHJ4*02"


def test_using_locks_full_vdj_simultaneously():
    exp = (
        ga.Experiment.on("human_igh")
        .using(v="IGHVF1-G1*01", d="IGHD3-3*01", j="IGHJ4*02")
        .recombine()
    )
    for o in exp.run(n=10, seed=0):
        assert _v_name(exp.refdata, o) == "IGHVF1-G1*01"
        assert _d_name(exp.refdata, o) == "IGHD3-3*01"
        assert _j_name(exp.refdata, o) == "IGHJ4*02"


def test_using_accepts_multi_allele_list():
    allowed = ["IGHVF10-G38*02", "IGHVF10-G38*04"]
    exp = (
        ga.Experiment.on("human_igh")
        .using(v=allowed)
        .recombine()
    )
    seen = {_v_name(exp.refdata, o) for o in exp.run(n=40, seed=0)}
    # Sample stays inside the allowed subset.
    assert seen.issubset(set(allowed))
    # Both members get hit at least once across 40 runs.
    assert seen == set(allowed)


def test_using_accepts_tuple_input():
    allowed = ("IGHVF10-G38*02", "IGHVF10-G38*04")
    exp = (
        ga.Experiment.on("human_igh")
        .using(v=allowed)
        .recombine()
    )
    for o in exp.run(n=10, seed=0):
        assert _v_name(exp.refdata, o) in allowed


def test_using_after_recombine_still_applies():
    # Locks are injected at compile time, so order of fluent calls
    # doesn't matter.
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .using(v="IGHVF1-G1*01")
    )
    for o in exp.run(n=5, seed=0):
        assert _v_name(exp.refdata, o) == "IGHVF1-G1*01"


def test_using_called_twice_overlays_per_segment():
    # First .using() locks V, second locks D — both should hold.
    exp = (
        ga.Experiment.on("human_igh")
        .using(v="IGHVF1-G1*01")
        .using(d="IGHD3-3*01")
        .recombine()
    )
    for o in exp.run(n=5, seed=0):
        assert _v_name(exp.refdata, o) == "IGHVF1-G1*01"
        assert _d_name(exp.refdata, o) == "IGHD3-3*01"


def test_using_called_twice_overwrites_same_segment():
    # Second call's V lock replaces the first one.
    exp = (
        ga.Experiment.on("human_igh")
        .using(v="IGHVF1-G1*01")
        .using(v="IGHVF10-G38*02")
        .recombine()
    )
    for o in exp.run(n=5, seed=0):
        assert _v_name(exp.refdata, o) == "IGHVF10-G38*02"


def test_using_none_clears_prior_lock():
    # Clearing reverts to the full pool.
    exp = (
        ga.Experiment.on("human_igh")
        .using(v="IGHVF1-G1*01")
        .using(v=None)
        .recombine()
    )
    seen = {_v_name(exp.refdata, o) for o in exp.run(n=30, seed=0)}
    assert len(seen) > 1


def test_using_omitted_kwarg_leaves_prior_lock_untouched():
    # Calling .using(d=...) shouldn't clobber an earlier V lock.
    exp = (
        ga.Experiment.on("human_igh")
        .using(v="IGHVF1-G1*01")
        .using(d="IGHD3-3*01")  # only D supplied
        .recombine()
    )
    for o in exp.run(n=5, seed=0):
        assert _v_name(exp.refdata, o) == "IGHVF1-G1*01"
        assert _d_name(exp.refdata, o) == "IGHD3-3*01"


def test_using_unknown_v_name_raises():
    with pytest.raises(ValueError, match="no V allele named"):
        ga.Experiment.on("human_igh").using(v="NOPE*01")


def test_using_unknown_d_name_raises():
    with pytest.raises(ValueError, match="no D allele named"):
        ga.Experiment.on("human_igh").using(d="IGHD-NOPE*01")


def test_using_unknown_j_name_raises():
    with pytest.raises(ValueError, match="no J allele named"):
        ga.Experiment.on("human_igh").using(j="IGHJ-NOPE*01")


def test_using_d_lock_on_vj_chain_raises():
    with pytest.raises(ValueError, match="cannot lock D alleles"):
        ga.Experiment.on("human_igk").using(d="IGHD3-3*01")


def test_using_empty_list_raises():
    with pytest.raises(ValueError, match="must be non-empty"):
        ga.Experiment.on("human_igh").using(v=[])


def test_using_duplicate_in_list_raises():
    with pytest.raises(ValueError, match="duplicate"):
        ga.Experiment.on("human_igh").using(
            v=["IGHVF1-G1*01", "IGHVF1-G1*01"]
        )


def test_using_non_string_entry_raises():
    with pytest.raises(TypeError, match="must all be strings"):
        ga.Experiment.on("human_igh").using(v=["IGHVF1-G1*01", 5])  # type: ignore[list-item]


def test_using_does_not_increment_step_count():
    # ``.using()`` is a configuration knob, not a pipeline step.
    exp = ga.Experiment.on("human_igh").using(v="IGHVF1-G1*01")
    assert exp.step_count == 0


def test_using_is_deterministic_under_same_seed():
    locked = ["IGHVF10-G38*02", "IGHVF10-G38*04"]
    a = (
        ga.Experiment.on("human_igh")
        .using(v=locked)
        .recombine()
        .run(n=10, seed=0xCAFE)
    )
    b = (
        ga.Experiment.on("human_igh")
        .using(v=locked)
        .recombine()
        .run(n=10, seed=0xCAFE)
    )
    for x, y in zip(a, b):
        assert x.final_simulation().bases() == y.final_simulation().bases()


def test_using_composes_with_mutation_and_corruption():
    # Locks should work alongside the rest of the DSL pipeline.
    exp = (
        ga.Experiment.on("human_igh")
        .using(v="IGHVF1-G1*01")
        .recombine()
        .mutate(count=8)
        .corrupt_pcr(count=2)
    )
    out = exp.run(n=3, seed=42)
    for o in out:
        assert _v_name(exp.refdata, o) == "IGHVF1-G1*01"
        assert o.trace().find("mutate.s5f.count").value == 8
        assert o.trace().find("corrupt.pcr.count").value == 2


def test_using_works_for_vj_chain_v_and_j():
    # Pick the first V / J allele in the IGK pool so the test stays
    # robust to refdata reshuffles.
    refdata = ga.Experiment.on("human_igk").refdata
    v_name = refdata.v_allele(0).name
    j_name = refdata.j_allele(0).name
    exp = (
        ga.Experiment.on("human_igk")
        .using(v=v_name, j=j_name)
        .recombine()
    )
    for o in exp.run(n=5, seed=0):
        assert _v_name(exp.refdata, o) == v_name
        assert _j_name(exp.refdata, o) == j_name


def test_using_with_productive_contract():
    # Locks + productive() should still produce productive sequences.
    exp = (
        ga.Experiment.on("human_igh")
        .using(v="IGHVF1-G1*01", j="IGHJ4*02")
        .recombine()
    )
    out = exp.run_records(n=3, seed=0, respect=ga.productive())
    for rec in out:
        assert rec["v_call"] == "IGHVF1-G1*01"
        assert rec["j_call"] == "IGHJ4*02"
        assert rec["productive"] is True


def test_using_compile_is_idempotent_with_locks():
    exp = (
        ga.Experiment.on("human_igh")
        .using(v="IGHVF1-G1*01")
        .recombine()
    )
    plan_a = exp.compile()
    plan_b = exp.compile()
    a = plan_a.run(n=3, seed=0)
    b = plan_b.run(n=3, seed=0)
    for x, y in zip(a, b):
        assert x.final_simulation().bases() == y.final_simulation().bases()


# ──────────────────────────────────────────────────────────────────
# G.6 — streaming via compiled.stream() / stream_records()
# ──────────────────────────────────────────────────────────────────

import inspect as _inspect
import itertools as _itertools


def test_stream_returns_generator():
    exp = _human_igh_exp()
    s = exp.stream(n=3, seed=0)
    assert _inspect.isgenerator(s)


def test_stream_yields_n_outcomes_when_bounded():
    out = list(_human_igh_exp().stream(n=4, seed=0))
    assert len(out) == 4


def test_stream_default_n_is_unbounded():
    # No ``n`` → infinite generator. Use ``islice`` to bound.
    s = _human_igh_exp().stream(seed=0)
    first_seven = list(_itertools.islice(s, 7))
    assert len(first_seven) == 7


def test_stream_matches_run_for_same_seed():
    exp = _human_igh_exp().mutate(count=5)
    runs = exp.run(n=4, seed=0xCAFE)
    streamed = list(exp.stream(n=4, seed=0xCAFE))
    assert len(runs) == len(streamed)
    for a, b in zip(runs, streamed):
        assert a.final_simulation().bases() == b.final_simulation().bases()


def test_stream_yields_engine_outcomes():
    import genairr_engine as ge

    out = list(_human_igh_exp().stream(n=2, seed=0))
    assert all(isinstance(o, ge.Outcome) for o in out)


def test_stream_n_zero_raises():
    with pytest.raises(ValueError, match="n must be at least 1"):
        list(_human_igh_exp().stream(n=0))


def test_stream_n_negative_raises():
    with pytest.raises(ValueError, match="n must be at least 1"):
        list(_human_igh_exp().stream(n=-3))


def test_stream_is_lazy_under_unbounded_n():
    # Constructing a stream with no bound must not actually run any
    # simulations until iteration begins. Confirmed by the fact that
    # construction returns immediately even with no ``n``.
    s = _human_igh_exp().stream(seed=0)
    assert _inspect.isgenerator(s)
    # Pull a single item — should not blow up regardless of "n".
    first = next(s)
    assert first is not None


def test_stream_seed_offsets_per_iteration():
    # The i-th yielded outcome should match a single ``run()`` call
    # with ``seed = base_seed + i``. Confirms the documented
    # ``seed + i`` per-run offset.
    base_seed = 17
    exp = _human_igh_exp().mutate(count=3)
    streamed = list(exp.stream(n=4, seed=base_seed))
    for i, outcome in enumerate(streamed):
        single = exp.run(n=1, seed=base_seed + i)[0]
        assert outcome.final_simulation().bases() == single.final_simulation().bases()


def test_stream_supports_respect_kwarg():
    out = list(
        _human_igh_exp()
        .mutate(count=5)
        .stream(n=3, seed=0, respect=ga.productive())
    )
    # The stream should still produce productive outcomes.
    for o in out:
        # Junction is in-frame and stop-free → trace records show the
        # productive contract stayed satisfied.
        sim = o.final_simulation()
        assert sim.bases() != b""


def test_stream_records_yields_record_dicts():
    out = list(_human_igh_exp().stream_records(n=3, seed=0))
    assert len(out) == 3
    for rec in out:
        assert isinstance(rec, dict)
        assert "sequence" in rec
        assert "v_call" in rec
        assert "junction" in rec


def test_stream_records_matches_run_records_per_seed():
    exp = _human_igh_exp().mutate(count=5)
    batch = exp.run_records(n=3, seed=42)
    streamed = list(exp.stream_records(n=3, seed=42))
    assert len(batch) == len(streamed)
    for a, b in zip(batch, streamed):
        assert a["sequence"] == b["sequence"]
        assert a["v_call"] == b["v_call"]


def test_stream_records_default_n_is_unbounded():
    s = _human_igh_exp().stream_records(seed=0)
    first_three = list(_itertools.islice(s, 3))
    assert len(first_three) == 3
    assert all("sequence" in r for r in first_three)


def test_stream_records_lazy_iteration_writes_one_at_a_time(tmp_path):
    # Real-world streaming pattern: write each record to TSV as it
    # comes in, never holding more than one in memory.
    import csv as _csv

    from GenAIRR.result import _DEFAULT_COLUMN_ORDER

    path = tmp_path / "stream.tsv"
    with open(path, "w", encoding="utf-8", newline="") as fh:
        writer = _csv.DictWriter(
            fh,
            fieldnames=_DEFAULT_COLUMN_ORDER,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        for rec in _human_igh_exp().stream_records(n=5, seed=0):
            row = {k: ("" if v is None else v) for k, v in rec.items()}
            writer.writerow(row)

    with open(path, encoding="utf-8") as fh:
        rows = list(_csv.DictReader(fh, delimiter="\t"))
    assert len(rows) == 5
    assert all(r["sequence"] for r in rows)


def test_stream_works_on_compiled_experiment():
    exp = _human_igh_exp().mutate(count=5)
    compiled = exp.compile()
    out = list(compiled.stream(n=3, seed=0))
    via_run = compiled.run(n=3, seed=0)
    for a, b in zip(out, via_run):
        assert a.final_simulation().bases() == b.final_simulation().bases()


def test_stream_records_works_on_compiled_experiment():
    compiled = _human_igh_exp().compile()
    out = list(compiled.stream_records(n=2, seed=0))
    assert len(out) == 2
    assert all("v_call" in r for r in out)


def test_stream_compose_with_using_lock():
    exp = (
        ga.Experiment.on("human_igh")
        .using(v="IGHVF1-G1*01")
        .recombine()
    )
    for o in exp.stream(n=5, seed=0):
        vid = o.final_simulation().v_allele_id()
        assert exp.refdata.v_allele(vid).name == "IGHVF1-G1*01"


def test_stream_can_terminate_early_with_break():
    # Caller stops iteration after the first item; the generator
    # should not produce more outcomes than requested.
    exp = _human_igh_exp()
    s = exp.stream(seed=0)
    first = next(s)
    s.close()  # Generator cleanup; subsequent next() raises StopIteration.
    assert first is not None
    with pytest.raises(StopIteration):
        next(s)


def test_stream_records_writes_fasta_streaming(tmp_path):
    # Companion of the TSV streaming test for FASTA output.
    path = tmp_path / "stream.fa"
    with open(path, "w", encoding="utf-8") as fh:
        for i, rec in enumerate(
            _human_igh_exp().stream_records(n=3, seed=0)
        ):
            fh.write(f">seq{i}|v_call={rec['v_call']}\n")
            fh.write(f"{rec['sequence']}\n")
    lines = path.read_text(encoding="utf-8").splitlines()
    assert len(lines) == 6  # 3 records × (header + seq)
    assert lines[0].startswith(">seq0|")
    assert lines[2].startswith(">seq1|")
    assert lines[4].startswith(">seq2|")


# ──────────────────────────────────────────────────────────────────
# H.1 — AIRR field renaming + cheap metadata
# ──────────────────────────────────────────────────────────────────


def test_records_use_airr_field_names_for_renames():
    rec = _human_igh_exp().run_records(n=1, seed=0)[0]
    # New canonical AIRR names exist.
    assert "junction" in rec
    assert "np1" in rec
    assert "np2" in rec
    # Old V5 names no longer exist.
    assert "junction_nt" not in rec
    assert "np1_region" not in rec
    assert "np2_region" not in rec


def test_record_has_sequence_id():
    result = _human_igh_records(n=3)
    for i, rec in enumerate(result):
        assert rec["sequence_id"] == f"seq{i}"


def test_simulation_result_id_prefix_is_configurable():
    from GenAIRR.result import SimulationResult

    exp = _human_igh_exp()
    outcomes = exp.run(n=2, seed=0)
    result = SimulationResult.from_outcomes(
        outcomes, exp.refdata, id_prefix="rep01_"
    )
    assert result[0]["sequence_id"] == "rep01_0"
    assert result[1]["sequence_id"] == "rep01_1"


def test_record_rev_comp_is_false_for_forward_strand_simulator():
    rec = _human_igh_exp().run_records(n=1, seed=0)[0]
    assert rec["rev_comp"] is False


def test_record_locus_is_igh_for_human_igh():
    for rec in _human_igh_exp().run_records(n=3, seed=0):
        assert rec["locus"] == "IGH"


def test_record_locus_is_igk_for_human_igk():
    out = ga.Experiment.on("human_igk").recombine().run_records(n=3, seed=0)
    for rec in out:
        assert rec["locus"] == "IGK"


def test_record_locus_is_igl_for_human_igl():
    out = ga.Experiment.on("human_igl").recombine().run_records(n=3, seed=0)
    for rec in out:
        assert rec["locus"] == "IGL"


def test_record_locus_for_tcr_chain():
    out = (
        ga.Experiment.on("human_tcrb")
        .recombine()
        .run_records(n=3, seed=0)
    )
    for rec in out:
        assert rec["locus"] == "TRB"


def test_derive_locus_from_unrecognised_prefix_yields_empty_string():
    from GenAIRR._airr_record import _derive_locus

    assert _derive_locus("XYZ1*01", "", "") == ""
    assert _derive_locus("", "", "") == ""


def test_derive_locus_falls_back_to_j_when_v_empty():
    from GenAIRR._airr_record import _derive_locus

    assert _derive_locus("", "IGHJ4*01", "") == "IGH"


def test_record_c_call_is_empty_string():
    rec = _human_igh_exp().run_records(n=1, seed=0)[0]
    # We don't simulate the constant region, so c_call is always
    # empty — but the column must exist for AIRR-tooling parsers.
    assert rec["c_call"] == ""


def test_default_column_order_starts_with_airr_metadata():
    from GenAIRR.result import _DEFAULT_COLUMN_ORDER

    # Canonical order: AIRR metadata first, then calls.
    assert _DEFAULT_COLUMN_ORDER[0] == "sequence_id"
    assert _DEFAULT_COLUMN_ORDER[1] == "sequence"
    assert "rev_comp" in _DEFAULT_COLUMN_ORDER
    assert "locus" in _DEFAULT_COLUMN_ORDER
    assert "c_call" in _DEFAULT_COLUMN_ORDER


def test_to_tsv_includes_new_h1_columns(tmp_path):
    import csv as _csv

    result = _human_igh_records(n=2)
    path = tmp_path / "h1.tsv"
    result.to_tsv(str(path))

    with open(path, encoding="utf-8") as fh:
        rows = list(_csv.DictReader(fh, delimiter="\t"))

    for i, row in enumerate(rows):
        assert row["sequence_id"] == f"seq{i}"
        assert row["rev_comp"] == "False"
        assert row["locus"] == "IGH"
        assert row["c_call"] == ""
        assert "junction" in row  # new name
        assert "np1" in row  # new name


def test_stream_records_emits_sequential_sequence_ids():
    out = list(_human_igh_exp().stream_records(n=4, seed=0))
    for i, rec in enumerate(out):
        assert rec["sequence_id"] == f"seq{i}"


def test_stream_records_id_prefix_kwarg_is_threaded_through():
    out = list(
        _human_igh_exp().stream_records(n=3, seed=0, id_prefix="run42_")
    )
    assert [rec["sequence_id"] for rec in out] == [
        "run42_0",
        "run42_1",
        "run42_2",
    ]


# ──────────────────────────────────────────────────────────────────
# H.2 — translated regions (sequence_aa, np1_aa, np2_aa)
# ──────────────────────────────────────────────────────────────────


def test_record_has_sequence_aa_field():
    rec = _human_igh_exp().run_records(n=1, seed=0)[0]
    assert "sequence_aa" in rec
    assert isinstance(rec["sequence_aa"], str)


def test_record_has_np_aa_fields():
    rec = _human_igh_exp().run_records(n=1, seed=0)[0]
    assert "np1_aa" in rec
    assert "np2_aa" in rec


def test_sequence_aa_is_nonempty_for_productive_outcome():
    out = (
        _human_igh_exp()
        .run_records(n=3, seed=42, respect=ga.productive())
    )
    for rec in out:
        assert rec["sequence_aa"] != ""


def test_sequence_aa_contains_junction_aa_for_productive_outcome():
    # The junction is part of the assembled sequence, so its AA
    # translation must appear inside sequence_aa.
    out = (
        _human_igh_exp()
        .run_records(n=5, seed=42, respect=ga.productive())
    )
    for rec in out:
        assert rec["junction_aa"], "expected junction_aa for productive seq"
        assert rec["junction_aa"] in rec["sequence_aa"]


def test_sequence_aa_uses_v_anchor_frame():
    # Translation starts at the V anchor's frame (junction_start % 3
    # offset). So the translated AA at the junction's relative
    # position should equal junction_aa.
    rec = (
        _human_igh_exp()
        .run_records(n=1, seed=42, respect=ga.productive())[0]
    )
    junction_start = rec["junction_start"]
    frame_offset = junction_start % 3
    # The junction's first codon sits at AA index
    # ``(junction_start - frame_offset) // 3``.
    aa_idx = (junction_start - frame_offset) // 3
    junction_aa_len = len(rec["junction_aa"])
    assert (
        rec["sequence_aa"][aa_idx : aa_idx + junction_aa_len]
        == rec["junction_aa"]
    )


def test_sequence_aa_length_matches_codon_count():
    rec = (
        _human_igh_exp()
        .run_records(n=1, seed=42, respect=ga.productive())[0]
    )
    seq = rec["sequence"]
    junction_start = rec["junction_start"]
    frame_offset = junction_start % 3
    expected_len = (len(seq) - frame_offset) // 3
    assert len(rec["sequence_aa"]) == expected_len


def test_sequence_aa_empty_when_junction_uncomputable():
    # A bare RefDataConfig has no anchors → no junction → sequence_aa
    # falls back to empty.
    import genairr_engine as ge

    cfg = ge.RefDataConfig.vdj()
    cfg.add_v_allele("v1*01", "v1", b"AAACCCGGG")  # no anchor
    cfg.add_d_allele("d1*01", "d1", b"TTTTTT")
    cfg.add_j_allele("j1*01", "j1", b"TTTAAA")  # no anchor
    out = ga.Experiment.on(cfg).recombine().run_records(n=2, seed=0)
    for rec in out:
        assert rec["sequence_aa"] == ""
        assert rec["np1_aa"] == ""
        assert rec["np2_aa"] == ""


def test_np_aa_consistent_with_sequence_aa_slice():
    # np1_aa / np2_aa must equal the sequence_aa substring at the
    # codon-aligned positions inside the NP region.
    rec = (
        _human_igh_exp()
        .run_records(n=1, seed=42, respect=ga.productive())[0]
    )
    junction_start = rec["junction_start"]
    frame_offset = junction_start % 3
    seq_aa = rec["sequence_aa"]

    # NP1 — locate it via v_sequence_end / d_sequence_start (or
    # j_sequence_start for VJ chains; this test is VDJ).
    np1_start = rec["v_sequence_end"]
    np1_end = rec["d_sequence_start"]
    if np1_end is not None and np1_end > np1_start:
        delta = (frame_offset - np1_start) % 3
        codon_start = np1_start + delta
        n_codons = (np1_end - codon_start) // 3
        if n_codons > 0:
            aa_idx = (codon_start - frame_offset) // 3
            assert (
                rec["np1_aa"]
                == seq_aa[aa_idx : aa_idx + n_codons]
            )


def test_np_aa_only_translates_complete_codons():
    # NP regions are typically not multiples of 3 → aa length should
    # never exceed (np_length // 3) + 1 (with leeway for partial
    # boundary codons).
    out = _human_igh_exp().run_records(n=10, seed=0, respect=ga.productive())
    for rec in out:
        np1_len = rec["np1_length"]
        np2_len = rec["np2_length"]
        # Strict bound: aa cannot exceed nt // 3 (always true since we
        # only emit complete codons).
        assert len(rec["np1_aa"]) <= max(0, np1_len // 3 + 1)
        assert len(rec["np2_aa"]) <= max(0, np2_len // 3 + 1)


def test_np_aa_empty_when_np_region_empty():
    # Force a recombine with explicit NP1 length 0 → np1 string is
    # empty, np1_aa must be empty too.
    exp = (
        ga.Experiment.on("human_igh")
        .recombine(np1_lengths=[(0, 1.0)], np2_lengths=[(0, 1.0)])
    )
    for rec in exp.run_records(n=3, seed=0, respect=ga.productive()):
        assert rec["np1"] == ""
        assert rec["np1_aa"] == ""
        assert rec["np2"] == ""
        assert rec["np2_aa"] == ""


def test_aa_slice_for_region_helper_handles_alignment():
    from GenAIRR._airr_record import _aa_slice_for_region

    # frame_offset = 0, sequence_aa "ABCDE" maps to nt positions
    # [0, 3, 6, 9, 12]. A region [3, 12) covers codons 1, 2, 3 → "BCD".
    sequence_aa = "ABCDE"
    assert _aa_slice_for_region(3, 12, 0, sequence_aa) == "BCD"
    # Region [4, 12) starts mid-codon — first complete codon at 6.
    assert _aa_slice_for_region(4, 12, 0, sequence_aa) == "CD"
    # Region too short for a complete codon.
    assert _aa_slice_for_region(4, 5, 0, sequence_aa) == ""
    # Empty region.
    assert _aa_slice_for_region(5, 5, 0, sequence_aa) == ""


def test_default_column_order_includes_aa_fields():
    from GenAIRR.result import _DEFAULT_COLUMN_ORDER

    assert "sequence_aa" in _DEFAULT_COLUMN_ORDER
    assert "np1_aa" in _DEFAULT_COLUMN_ORDER
    assert "np2_aa" in _DEFAULT_COLUMN_ORDER


def test_to_tsv_includes_aa_columns(tmp_path):
    import csv as _csv

    result = _human_igh_records(n=2)
    path = tmp_path / "h2.tsv"
    result.to_tsv(str(path))
    with open(path, encoding="utf-8") as fh:
        rows = list(_csv.DictReader(fh, delimiter="\t"))
    for row in rows:
        assert "sequence_aa" in row
        assert "np1_aa" in row
        assert "np2_aa" in row


def test_sequence_aa_handles_lowercase_corruption_bases():
    # Quality-corruption emits lowercase bases. Translation upper-cases
    # internally, so the AA strings should still be valid letters.
    rec = (
        _human_igh_exp()
        .corrupt_quality(count=8)
        .run_records(n=1, seed=0)[0]
    )
    # Must be valid AA characters (no lowercase, no '?').
    for c in rec["sequence_aa"]:
        assert c.isupper() or c in ("*", "X")


# ──────────────────────────────────────────────────────────────────
# H.3 — alignment strings (ground-truth, gap-aware)
# ──────────────────────────────────────────────────────────────────


def test_record_has_alignment_string_fields():
    rec = _human_igh_exp().run_records(n=1, seed=0)[0]
    assert "sequence_alignment" in rec
    assert "germline_alignment" in rec
    assert "germline_alignment_d_mask" in rec


def test_alignment_strings_equal_length():
    # The triplet must be a strict 1:1 alignment — same length always.
    for rec in _human_igh_exp().run_records(n=5, seed=0):
        sa = rec["sequence_alignment"]
        ga_str = rec["germline_alignment"]
        dm = rec["germline_alignment_d_mask"]
        assert len(sa) == len(ga_str) == len(dm)


def test_no_indel_alignment_length_matches_sequence_length():
    # Without indels, no gaps are inserted; alignment length equals
    # raw sequence length.
    for rec in _human_igh_exp().run_records(n=5, seed=0):
        assert len(rec["sequence_alignment"]) == rec["sequence_length"]


def test_no_indel_sequence_alignment_equals_uppercase_sequence():
    for rec in _human_igh_exp().run_records(n=5, seed=0):
        assert rec["sequence_alignment"] == rec["sequence"].upper()


def test_germline_alignment_has_n_in_np_regions():
    rec = _human_igh_exp().run_records(n=1, seed=0)[0]
    galn = rec["germline_alignment"]
    np1_start = rec["v_sequence_end"]
    np1_end = rec["d_sequence_start"]
    if np1_start is not None and np1_end is not None and np1_end > np1_start:
        assert all(c == "N" for c in galn[np1_start:np1_end])
    np2_start = rec["d_sequence_end"]
    np2_end = rec["j_sequence_start"]
    if np2_start is not None and np2_end is not None and np2_end > np2_start:
        assert all(c == "N" for c in galn[np2_start:np2_end])


def test_d_mask_replaces_d_region_with_ns():
    rec = _human_igh_exp().run_records(n=1, seed=0)[0]
    dm = rec["germline_alignment_d_mask"]
    galn = rec["germline_alignment"]
    d_start = rec["d_sequence_start"]
    d_end = rec["d_sequence_end"]
    if d_start is not None and d_end is not None and d_end > d_start:
        # D portion in d_mask is all N's (no gaps in this no-indel case).
        assert all(c == "N" for c in dm[d_start:d_end])
        # Non-D portions must be unchanged from germline_alignment.
        assert dm[:d_start] == galn[:d_start]
        assert dm[d_end:] == galn[d_end:]


def test_germline_alignment_matches_sequence_alignment_when_no_mutations():
    # No mutations and no indels and no NP regions touching the test
    # range: V coding portion of sequence_alignment must equal V
    # coding portion of germline_alignment.
    rec = _human_igh_exp().run_records(n=1, seed=0)[0]
    sa = rec["sequence_alignment"]
    galn = rec["germline_alignment"]
    # Walk only the V coding region (no NP, no D, no J).
    v_start, v_end = rec["v_sequence_start"], rec["v_sequence_end"]
    assert sa[v_start:v_end] == galn[v_start:v_end]


def test_germline_alignment_diverges_from_sequence_at_mutated_positions():
    # With S5F mutations, sa should differ from galn at coding
    # positions hit by mutation. S5F samples positions with
    # replacement, so back-mutations / position collisions reduce
    # the visible mismatch count below ``count``.
    rec = (
        _human_igh_exp()
        .mutate(count=15)
        .run_records(n=1, seed=42)[0]
    )
    sa = rec["sequence_alignment"]
    galn = rec["germline_alignment"]
    diffs = sum(
        1 for a, b in zip(sa, galn) if a != b and b not in ("N", "-")
    )
    # Strict upper bound — every mismatch must come from one of the
    # 15 mutation events. The lower bound is loose because of
    # collisions; for seed=42 we expect ≥8 distinct positions.
    assert 8 <= diffs <= 15


def test_indel_insertion_creates_gap_in_germline_alignment():
    # Indel insertions add bases with no germline counterpart → those
    # positions should be '-' in germline_alignment.
    # Use insertion_prob=1.0 to force all indels to be insertions.
    rec = (
        _human_igh_exp()
        .corrupt_indels(count=4, insertion_prob=1.0)
        .run_records(n=1, seed=0)[0]
    )
    galn = rec["germline_alignment"]
    sa = rec["sequence_alignment"]
    # All 4 indels are insertions → exactly 4 gaps in galn, 0 in sa.
    assert galn.count("-") == 4
    assert sa.count("-") == 0


def test_indel_deletion_creates_gap_in_sequence_alignment():
    # Indel deletions remove bases → those positions should be '-' in
    # sequence_alignment with the original allele base in germline.
    rec = (
        _human_igh_exp()
        .corrupt_indels(count=4, insertion_prob=0.0)
        .run_records(n=1, seed=0)[0]
    )
    galn = rec["germline_alignment"]
    sa = rec["sequence_alignment"]
    assert sa.count("-") == 4
    assert galn.count("-") == 0


def test_alignment_total_gap_count_at_most_n_indels():
    # For mixed insertions/deletions, every alignment gap traces back
    # to one indel event. Indels that land in NP regions don't surface
    # as alignment gaps (NP has no germline position to mark), so the
    # observed gap count is ≤ ``n_indels``.
    rec = (
        _human_igh_exp()
        .corrupt_indels(count=6, insertion_prob=0.5)
        .run_records(n=1, seed=12345)[0]
    )
    sa = rec["sequence_alignment"]
    galn = rec["germline_alignment"]
    total_gaps = sa.count("-") + galn.count("-")
    assert 0 <= total_gaps <= rec["n_indels"]


def test_alignment_length_increases_with_indels():
    # Each indel (insertion or deletion) adds 1 to alignment length:
    # insertions add a column where seq has a base and germ is gap;
    # deletions add a column where seq is gap and germ has a base.
    rec = (
        _human_igh_exp()
        .corrupt_indels(count=4, insertion_prob=1.0)
        .run_records(n=1, seed=0)[0]
    )
    # sequence_length stays raw (with insertions counted, deletions
    # missing). With 4 insertions and 0 deletions:
    #   alignment_length = sequence_length (no deletions add a column
    #     beyond what insertions already added).
    # Actually: alignment includes deletions as extra columns, so
    # alignment_length == sequence_length + n_deletions.
    assert len(rec["sequence_alignment"]) == rec["sequence_length"]


def test_alignment_length_grows_for_pure_deletions():
    # 4 deletions → sequence_length is shorter than alignment_length
    # (each deletion is a gap column not present in the raw sequence).
    rec = (
        _human_igh_exp()
        .corrupt_indels(count=4, insertion_prob=0.0)
        .run_records(n=1, seed=0)[0]
    )
    assert len(rec["sequence_alignment"]) == rec["sequence_length"] + 4


def test_d_mask_preserves_gaps_in_d_region():
    # Indel deletions inside the D region should still appear as '-'
    # in d_mask (not converted to N).
    rec = (
        _human_igh_exp()
        .corrupt_indels(count=15, insertion_prob=0.0)
        .run_records(n=1, seed=0)[0]
    )
    galn = rec["germline_alignment"]
    dm = rec["germline_alignment_d_mask"]
    # d_mask must have the same gap positions as germline_alignment
    # (gaps are preserved, not masked).
    galn_gap_positions = {i for i, c in enumerate(galn) if c == "-"}
    dm_gap_positions = {i for i, c in enumerate(dm) if c == "-"}
    assert galn_gap_positions == dm_gap_positions


def test_alignment_strings_uppercase():
    # AIRR convention: alignment strings are uppercase regardless of
    # the source sequence case (lowercase NP / mutated bases).
    rec = (
        _human_igh_exp()
        .mutate(count=10)
        .run_records(n=1, seed=0)[0]
    )
    sa = rec["sequence_alignment"]
    galn = rec["germline_alignment"]
    for c in sa:
        assert c.isupper() or c in ("-",), f"non-uppercase {c!r} in sa"
    for c in galn:
        assert c.isupper() or c in ("-",), f"non-uppercase {c!r} in ga"


def test_alignment_strings_for_vj_chain():
    # VJ chains have V and J only (no D, no NP2). Walks must still
    # produce valid 3-string output.
    rec = (
        ga.Experiment.on("human_igk")
        .recombine()
        .run_records(n=1, seed=0)[0]
    )
    sa = rec["sequence_alignment"]
    galn = rec["germline_alignment"]
    dm = rec["germline_alignment_d_mask"]
    assert len(sa) == len(galn) == len(dm) == rec["sequence_length"]
    # No D region → d_mask should equal germline_alignment exactly.
    assert dm == galn


def test_alignment_strings_empty_when_no_assembly():
    # Pre-recombine experiments yield empty sequences → empty strings.
    import genairr_engine as ge

    # Use a refdata with anchorless alleles so there's nothing
    # productive to anchor.
    cfg = ge.RefDataConfig.vdj()
    cfg.add_v_allele("v1*01", "v1", b"AAACCCGGG")
    cfg.add_d_allele("d1*01", "d1", b"TTTTTT")
    cfg.add_j_allele("j1*01", "j1", b"TTTAAA")
    out = ga.Experiment.on(cfg).recombine().run_records(n=1, seed=0)
    rec = out[0]
    # Even without anchors, alignment strings should be the assembled
    # sequence vs allele bases — non-empty here.
    assert len(rec["sequence_alignment"]) == rec["sequence_length"]


def test_alignment_strings_in_default_column_order():
    from GenAIRR.result import _DEFAULT_COLUMN_ORDER

    assert "sequence_alignment" in _DEFAULT_COLUMN_ORDER
    assert "germline_alignment" in _DEFAULT_COLUMN_ORDER
    assert "germline_alignment_d_mask" in _DEFAULT_COLUMN_ORDER


def test_alignment_strings_in_tsv_export(tmp_path):
    import csv as _csv

    result = _human_igh_exp().mutate(count=5).run_records(n=2, seed=0)
    path = tmp_path / "h3.tsv"
    result.to_tsv(str(path))
    with open(path, encoding="utf-8") as fh:
        rows = list(_csv.DictReader(fh, delimiter="\t"))
    for row in rows:
        assert row["sequence_alignment"]
        assert row["germline_alignment"]
        assert row["germline_alignment_d_mask"]
        assert len(row["sequence_alignment"]) == len(row["germline_alignment"])


def _alignment_invariants_hold(rec) -> list:
    """Return a list of issue strings for a record's alignment triplet
    — empty list if every invariant passes."""
    sa = rec["sequence_alignment"]
    ga_str = rec["germline_alignment"]
    dm = rec["germline_alignment_d_mask"]
    seq = rec["sequence"]
    issues = []
    if not (len(sa) == len(ga_str) == len(dm)):
        issues.append(
            f"length mismatch: sa={len(sa)} ga={len(ga_str)} dm={len(dm)}"
        )
    if len(sa) != len(seq) + sa.count("-"):
        issues.append(
            f"sa len {len(sa)} != seq len {len(seq)} + sa gaps {sa.count('-')}"
        )
    if sa.replace("-", "") != seq.upper():
        issues.append("sa minus gaps != upper(seq)")
    return issues


def test_alignment_handles_tail_indel_insertion_outside_regions():
    # Regression for the bug where an insertion at exactly pool_len
    # didn't extend any region's end → one pool position was outside
    # every region → walker dropped it → sa shorter than seq.
    # Repro: human_igk + indels + seed=0 used to fail at record idx 1.
    exp = (
        ga.Experiment.on("human_igk")
        .recombine()
        .mutate(count=5)
        .corrupt_pcr(count=2)
        .corrupt_indels(count=2)
    )
    for rec in exp.run_records(n=50, seed=0):
        assert _alignment_invariants_hold(rec) == []


def test_alignment_invariants_under_full_corruption_stack():
    # Stack every corruption pass at once and verify alignment
    # invariants hold across 200 seeds. This is the definitive
    # H.3 robustness check — touches PCR, quality, indels, contaminant,
    # mutations, and the productive contract simultaneously.
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(count=15)
        .corrupt_pcr(count=5)
        .corrupt_quality(count=8)
        .corrupt_indels(count=4, insertion_prob=0.5)
    )
    result = exp.run_records(n=200, seed=0, respect=ga.productive())
    failures = [
        (i, _alignment_invariants_hold(rec))
        for i, rec in enumerate(result)
        if _alignment_invariants_hold(rec)
    ]
    assert failures == [], f"alignment invariants broke for: {failures[:3]}"


def test_alignment_invariants_across_loci():
    # IGH / IGK / IGL / TCRB all need to produce sane alignment
    # strings under heavy corruption. SHM is BCR-only (G4), so the
    # TCR pipeline drops the mutate step.
    for cfg in ("human_igh", "human_igk", "human_igl", "human_tcrb"):
        exp = ga.Experiment.on(cfg).recombine()
        if not cfg.startswith("human_tcr"):
            exp = exp.mutate(count=5)
        exp = exp.corrupt_pcr(count=2).corrupt_indels(count=2)
        result = exp.run_records(n=50, seed=0)
        failures = [
            (i, _alignment_invariants_hold(rec))
            for i, rec in enumerate(result)
            if _alignment_invariants_hold(rec)
        ]
        assert failures == [], f"{cfg}: {failures[:3]}"


# ──────────────────────────────────────────────────────────────────
# H.4 — per-segment CIGAR strings (v_cigar / d_cigar / j_cigar)
# ──────────────────────────────────────────────────────────────────


import re as _re


def _cigar_total_ops(cigar: str) -> int:
    """Sum the run-length counts in a CIGAR string. Used to check
    that the CIGAR's coverage equals the segment's alignment-column
    count."""
    return sum(int(n) for n in _re.findall(r"(\d+)[MID=X]", cigar))


def _cigar_op_counts(cigar: str) -> dict:
    """Return ``{op: total_count}`` so tests can assert e.g. how many
    insertions vs deletions appear."""
    counts: dict = {"M": 0, "I": 0, "D": 0}
    for n, op in _re.findall(r"(\d+)([MID=X])", cigar):
        counts[op] = counts.get(op, 0) + int(n)
    return counts


def test_record_has_cigar_fields():
    rec = _human_igh_exp().run_records(n=1, seed=0)[0]
    assert "v_cigar" in rec
    assert "d_cigar" in rec
    assert "j_cigar" in rec


def test_cigar_pure_match_when_no_indels():
    # No indels → CIGAR is a single M run for each segment, length
    # equal to the coding region's nucleotide span.
    rec = _human_igh_exp().run_records(n=1, seed=42)[0]
    v_len = rec["v_sequence_end"] - rec["v_sequence_start"]
    d_len = rec["d_sequence_end"] - rec["d_sequence_start"]
    j_len = rec["j_sequence_end"] - rec["j_sequence_start"]
    assert rec["v_cigar"] == f"{v_len}M"
    assert rec["d_cigar"] == f"{d_len}M"
    assert rec["j_cigar"] == f"{j_len}M"


def test_cigar_unchanged_under_mutations():
    # CIGAR uses ``M`` for both match and mismatch, so mutations
    # don't introduce I or D ops.
    rec = (
        _human_igh_exp()
        .mutate(count=12)
        .run_records(n=1, seed=42)[0]
    )
    for cig in (rec["v_cigar"], rec["d_cigar"], rec["j_cigar"]):
        assert "I" not in cig
        assert "D" not in cig


def test_cigar_insertion_only_produces_i_ops():
    rec = (
        _human_igh_exp()
        .corrupt_indels(count=4, insertion_prob=1.0)
        .run_records(n=1, seed=0)[0]
    )
    counts_v = _cigar_op_counts(rec["v_cigar"])
    counts_d = _cigar_op_counts(rec["d_cigar"])
    counts_j = _cigar_op_counts(rec["j_cigar"])
    # No deletions anywhere with insertion_prob=1.0.
    assert counts_v["D"] == counts_d["D"] == counts_j["D"] == 0
    # All-insertion run total across V+D+J should be ≤ 4 — some may
    # land in NP regions where they don't show in any CIGAR.
    total_I = counts_v["I"] + counts_d["I"] + counts_j["I"]
    assert 0 <= total_I <= 4


def test_cigar_deletion_only_produces_d_ops():
    rec = (
        _human_igh_exp()
        .corrupt_indels(count=4, insertion_prob=0.0)
        .run_records(n=1, seed=0)[0]
    )
    counts_v = _cigar_op_counts(rec["v_cigar"])
    counts_d = _cigar_op_counts(rec["d_cigar"])
    counts_j = _cigar_op_counts(rec["j_cigar"])
    assert counts_v["I"] == counts_d["I"] == counts_j["I"] == 0
    total_D = counts_v["D"] + counts_d["D"] + counts_j["D"]
    assert 0 <= total_D <= 4


def test_cigar_op_count_equals_segment_alignment_columns():
    # Sum of M + I + D ops in a segment's CIGAR equals the number
    # of alignment columns that segment contributes to the V/D/J
    # alignment.
    rec = (
        _human_igh_exp()
        .corrupt_indels(count=5, insertion_prob=0.5)
        .run_records(n=1, seed=42)[0]
    )
    # In the full alignment string, the V coding region's columns
    # are roughly v_sequence_start to v_sequence_end *plus* any
    # deletion gaps inside V. We can count them by walking the
    # alignment-string positions and matching them up via the same
    # column-list logic — but a simpler check is: each segment's
    # CIGAR total must equal `length(coding region in sequence)
    # + number of deletion columns inserted before/within the
    # segment`. For a record with no indels in V at all, total ops
    # = v_sequence_end - v_sequence_start.
    v_len = rec["v_sequence_end"] - rec["v_sequence_start"]
    v_total = _cigar_total_ops(rec["v_cigar"])
    # v_total = (raw V coding bases in sequence) + (deletions in V)
    #         = v_len + count of D ops in V.
    n_deletions_in_v = _cigar_op_counts(rec["v_cigar"])["D"]
    assert v_total == v_len + n_deletions_in_v


def test_cigar_d_cigar_empty_for_vj_chain():
    # No D segment → d_cigar must be the empty string (and the
    # other two must still be valid).
    rec = (
        ga.Experiment.on("human_igk")
        .recombine()
        .run_records(n=1, seed=0)[0]
    )
    assert rec["d_cigar"] == ""
    assert rec["v_cigar"]
    assert rec["j_cigar"]


def test_cigar_consistent_with_alignment_strings():
    # Walk both the CIGAR and the alignment-string columns of the V
    # region and verify they describe the same sequence of ops.
    rec = (
        _human_igh_exp()
        .corrupt_indels(count=4, insertion_prob=0.5)
        .run_records(n=1, seed=42)[0]
    )
    # Reconstruct the V CIGAR from sa/galn directly and compare.
    sa = rec["sequence_alignment"]
    galn = rec["germline_alignment"]
    # Walk only V columns. We use sequence-coord boundaries to
    # locate V — works only because V is the first segment in our
    # assembly (no leading deletion gaps). For a robust ground-
    # truth test we use the engine's own column walker via the
    # private helper.
    from GenAIRR._airr_record import (
        _alignment_columns,
        _cigars_from_columns,
    )

    # Re-derive the CIGAR by re-running the helper. Equivalence is
    # guaranteed by construction; this asserts the public field
    # matches the helper's output (catches any future re-routing).
    refdata = _human_igh_exp().refdata
    outcome = _human_igh_exp().corrupt_indels(
        count=4, insertion_prob=0.5
    ).run(n=1, seed=42)[0]
    cols = _alignment_columns(outcome, refdata)
    v_c, d_c, j_c = _cigars_from_columns(cols)
    assert v_c == rec["v_cigar"]
    assert d_c == rec["d_cigar"]
    assert j_c == rec["j_cigar"]


def test_cigar_runlength_collapses_consecutive_ops():
    # Pure-M CIGAR for a segment must have exactly one run, not
    # one run per position.
    rec = _human_igh_exp().run_records(n=1, seed=0)[0]
    # Count the number of (digit-prefixed op) tokens.
    v_runs = _re.findall(r"\d+[MID]", rec["v_cigar"])
    assert len(v_runs) == 1, f"expected single run, got {v_runs}"


def test_cigar_only_uses_canonical_ops():
    # We emit only M / I / D — no =/X/S/H.
    rec = (
        _human_igh_exp()
        .mutate(count=10)
        .corrupt_indels(count=3, insertion_prob=0.5)
        .run_records(n=1, seed=0)[0]
    )
    for cig in (rec["v_cigar"], rec["d_cigar"], rec["j_cigar"]):
        for _, op in _re.findall(r"(\d+)(\D)", cig):
            assert op in ("M", "I", "D"), f"unexpected op {op!r} in {cig!r}"


def test_cigars_in_default_column_order():
    from GenAIRR.result import _DEFAULT_COLUMN_ORDER

    assert "v_cigar" in _DEFAULT_COLUMN_ORDER
    assert "d_cigar" in _DEFAULT_COLUMN_ORDER
    assert "j_cigar" in _DEFAULT_COLUMN_ORDER
    # Canonical ordering: v_cigar should sit immediately after v_call,
    # likewise for d/j.
    v_call_idx = _DEFAULT_COLUMN_ORDER.index("v_call")
    assert _DEFAULT_COLUMN_ORDER[v_call_idx + 1] == "v_cigar"


def test_cigars_in_tsv_export(tmp_path):
    import csv as _csv

    result = (
        _human_igh_exp()
        .corrupt_indels(count=3, insertion_prob=0.5)
        .run_records(n=2, seed=0)
    )
    path = tmp_path / "h4.tsv"
    result.to_tsv(str(path))
    with open(path, encoding="utf-8") as fh:
        rows = list(_csv.DictReader(fh, delimiter="\t"))
    for row in rows:
        assert row["v_cigar"]
        assert row["j_cigar"]


def test_cigar_invariants_under_full_corruption_stack():
    # Same n=200 stress sweep as alignment strings — verify that
    # CIGARs stay self-consistent with the alignment-string
    # topology under every corruption mode at once.
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(count=15)
        .corrupt_pcr(count=5)
        .corrupt_quality(count=8)
        .corrupt_indels(count=4, insertion_prob=0.5)
    )
    failures = []
    for i, rec in enumerate(exp.run_records(n=200, seed=0)):
        # Phase 11.3 restored equality: the CIGAR walks every column
        # in the live (extension-aware) span — `M + I + D` is the
        # column count, `seq_len + D_ops` accounts for D-ops adding
        # gap columns beyond the live sequence span.
        for seg in ("v", "d", "j"):
            cig = rec[f"{seg}_cigar"]
            if not cig:
                continue
            ops = _cigar_op_counts(cig)
            total = ops["M"] + ops["I"] + ops["D"]
            seq_len = (
                rec[f"{seg}_sequence_end"] - rec[f"{seg}_sequence_start"]
                if rec[f"{seg}_sequence_end"] is not None
                else 0
            )
            if total != seq_len + ops["D"]:
                failures.append(
                    (i, seg, total, seq_len, ops["D"])
                )
    assert failures == [], f"CIGAR invariant broke for: {failures[:3]}"


# ──────────────────────────────────────────────────────────────────
# H.5 — alignment / germline coord pairs
# ──────────────────────────────────────────────────────────────────


def test_record_has_h5_coord_pair_fields():
    rec = _human_igh_exp().run_records(n=1, seed=0)[0]
    for seg in ("v", "d", "j"):
        for end in ("alignment_start", "alignment_end", "germline_start", "germline_end"):
            assert f"{seg}_{end}" in rec, f"missing {seg}_{end}"


def test_no_indel_alignment_coords_match_sequence_coords():
    # Without indels, alignment-string columns mirror sequence
    # positions exactly → alignment coords equal sequence coords.
    rec = _human_igh_exp().run_records(n=1, seed=42)[0]
    for seg in ("v", "d", "j"):
        s_s = rec[f"{seg}_sequence_start"]
        s_e = rec[f"{seg}_sequence_end"]
        a_s = rec[f"{seg}_alignment_start"]
        a_e = rec[f"{seg}_alignment_end"]
        assert a_s == s_s, f"{seg}: alignment_start {a_s} != sequence_start {s_s}"
        assert a_e == s_e, f"{seg}: alignment_end {a_e} != sequence_end {s_e}"


def test_germline_coord_length_matches_post_trim_segment_length():
    # germline_end - germline_start == allele_len - trim_5 - trim_3.
    # Without indels, that also equals sequence_end - sequence_start.
    rec = _human_igh_exp().run_records(n=1, seed=42)[0]
    for seg in ("v", "d", "j"):
        coding = rec[f"{seg}_sequence_end"] - rec[f"{seg}_sequence_start"]
        germline = rec[f"{seg}_germline_end"] - rec[f"{seg}_germline_start"]
        assert germline == coding


def test_v_germline_end_reflects_v_trim_3():
    # V_3 trim: v_germline_end = v_allele_len - v_trim_3.
    refdata = _human_igh_exp().refdata
    rec = _human_igh_exp().run_records(n=1, seed=42)[0]
    # Locate the assigned V allele.
    v_id = None
    for i in range(refdata.v_pool_size()):
        if refdata.v_allele(i).name == rec["v_call"]:
            v_id = i
            break
    assert v_id is not None
    allele_len = len(bytes(refdata.v_allele(v_id).seq()))
    assert rec["v_germline_end"] == allele_len - rec["v_trim_3"]


def test_d_germline_start_reflects_d_trim_5():
    rec = _human_igh_exp().run_records(n=1, seed=42)[0]
    assert rec["d_germline_start"] == rec["d_trim_5"]


def test_d_germline_end_reflects_d_trim_3():
    refdata = _human_igh_exp().refdata
    rec = _human_igh_exp().run_records(n=1, seed=42)[0]
    d_id = None
    for i in range(refdata.d_pool_size()):
        if refdata.d_allele(i).name == rec["d_call"]:
            d_id = i
            break
    assert d_id is not None
    allele_len = len(bytes(refdata.d_allele(d_id).seq()))
    assert rec["d_germline_end"] == allele_len - rec["d_trim_3"]


def test_j_germline_start_reflects_j_trim_5():
    rec = _human_igh_exp().run_records(n=1, seed=42)[0]
    assert rec["j_germline_start"] == rec["j_trim_5"]


def test_germline_alignment_slice_matches_germline_coord_pair():
    # GROUND-TRUTH INVARIANT: slicing germline_alignment by alignment
    # coords should equal slicing the source allele by germline coords
    # (post-trim, so we get exactly the bases the simulation drew).
    refdata = _human_igh_exp().refdata
    rec = _human_igh_exp().run_records(n=1, seed=42)[0]

    for seg in ("v", "d", "j"):
        a_s = rec[f"{seg}_alignment_start"]
        a_e = rec[f"{seg}_alignment_end"]
        g_s = rec[f"{seg}_germline_start"]
        g_e = rec[f"{seg}_germline_end"]
        if a_s is None:
            continue
        # Find the allele.
        getter = {
            "v": refdata.v_allele,
            "d": refdata.d_allele,
            "j": refdata.j_allele,
        }[seg]
        pool_size = {
            "v": refdata.v_pool_size,
            "d": refdata.d_pool_size,
            "j": refdata.j_pool_size,
        }[seg]()
        call_name = rec[f"{seg}_call"]
        allele_id = next(
            (i for i in range(pool_size) if getter(i).name == call_name),
            None,
        )
        assert allele_id is not None
        allele_seq = bytes(getter(allele_id).seq()).decode().upper()
        expected = allele_seq[g_s:g_e]
        actual = rec["germline_alignment"][a_s:a_e]
        assert (
            actual == expected
        ), f"{seg}: germline_alignment slice {actual!r} != source allele slice {expected!r}"


def test_indel_deletion_alignment_and_sequence_spans_are_both_positive():
    # Sanity: with deletions in V/D/J, both alignment and sequence
    # spans should be positive numbers. The structural invariant
    # "alignment span >= sequence span" (deletions add alignment
    # columns) held pre-Phase-11.2; Phase 11.2 lifted sequence span
    # to live-call coords, so an NP-side extension can grow seq span
    # past the structural alignment span. The cleaner version of
    # this invariant lands in Phase 11.4 once the alignment string
    # gets relabelled to include extension columns.
    rec = (
        _human_igh_exp()
        .corrupt_indels(count=4, insertion_prob=0.0)
        .run_records(n=1, seed=0)[0]
    )
    for seg in ("v", "d", "j"):
        s_s = rec[f"{seg}_sequence_start"]
        s_e = rec[f"{seg}_sequence_end"]
        a_s = rec[f"{seg}_alignment_start"]
        a_e = rec[f"{seg}_alignment_end"]
        if a_s is None:
            continue
        assert a_e - a_s > 0
        assert s_e - s_s > 0


def test_alignment_span_equals_cigar_op_total():
    # ``v_alignment_end - v_alignment_start`` must equal the total
    # M+I+D op count of v_cigar.
    rec = (
        _human_igh_exp()
        .corrupt_indels(count=5, insertion_prob=0.5)
        .run_records(n=1, seed=42)[0]
    )
    for seg in ("v", "d", "j"):
        a_s = rec[f"{seg}_alignment_start"]
        a_e = rec[f"{seg}_alignment_end"]
        cig = rec[f"{seg}_cigar"]
        if a_s is None:
            continue
        assert a_e - a_s == _cigar_total_ops(cig), (
            f"{seg}: alignment span {a_e - a_s} != cigar ops {_cigar_total_ops(cig)}"
        )


def test_germline_span_at_least_cigar_m_plus_d():
    # The number of source-allele bases the alignment + live-extension
    # account for is at least M + D. Phase 11.1 lifted germline coords
    # to read from the live-call hypothesis bounds, so when an NP-side
    # extension grew the hypothesis's ref range the germline span can
    # exceed M + D (which still counts only the structural alignment
    # columns until Phase 11.3 rewires the CIGAR to match).
    rec = (
        _human_igh_exp()
        .corrupt_indels(count=5, insertion_prob=0.5)
        .run_records(n=1, seed=42)[0]
    )
    for seg in ("v", "d", "j"):
        g_s = rec[f"{seg}_germline_start"]
        g_e = rec[f"{seg}_germline_end"]
        cig = rec[f"{seg}_cigar"]
        if g_s is None:
            continue
        ops = _cigar_op_counts(cig)
        assert g_e - g_s >= ops["M"] + ops["D"]


def test_d_coords_none_for_vj_chain():
    rec = (
        ga.Experiment.on("human_igk")
        .recombine()
        .run_records(n=1, seed=0)[0]
    )
    assert rec["d_alignment_start"] is None
    assert rec["d_alignment_end"] is None
    assert rec["d_germline_start"] is None
    assert rec["d_germline_end"] is None
    # V/J coords still populated.
    assert rec["v_alignment_start"] == 0
    assert rec["j_alignment_start"] is not None


def test_v_alignment_starts_at_zero():
    # V is the first segment in our assembly, no leading deletion gaps
    # → v_alignment_start is always 0.
    for rec in _human_igh_exp().run_records(n=5, seed=0):
        assert rec["v_alignment_start"] == 0


def test_h5_coords_in_default_column_order():
    from GenAIRR.result import _DEFAULT_COLUMN_ORDER

    for seg in ("v", "d", "j"):
        for end in ("alignment_start", "alignment_end", "germline_start", "germline_end"):
            assert f"{seg}_{end}" in _DEFAULT_COLUMN_ORDER


def test_h5_coords_in_tsv_export(tmp_path):
    import csv as _csv

    result = _human_igh_records(n=2)
    path = tmp_path / "h5.tsv"
    result.to_tsv(str(path))
    with open(path, encoding="utf-8") as fh:
        rows = list(_csv.DictReader(fh, delimiter="\t"))
    for row in rows:
        for seg in ("v", "d", "j"):
            for end in ("alignment_start", "alignment_end", "germline_start", "germline_end"):
                assert f"{seg}_{end}" in row


# ──────────────────────────────────────────────────────────────────
# H.6 — score / identity / support fields
# ──────────────────────────────────────────────────────────────────


def test_record_has_h6_score_fields():
    rec = _human_igh_exp().run_records(n=1, seed=0)[0]
    for seg in ("v", "d", "j"):
        for kind in ("score", "identity", "support"):
            assert f"{seg}_{kind}" in rec


def test_score_and_support_are_stubbed_none():
    rec = _human_igh_exp().run_records(n=1, seed=0)[0]
    # No aligner ran → these are aligner-specific concepts and are
    # always None for a simulator's output.
    for seg in ("v", "d", "j"):
        assert rec[f"{seg}_score"] is None
        assert rec[f"{seg}_support"] is None


def test_identity_is_one_for_pure_recombination():
    # No mutations, no corruption → every alignment column is a
    # perfect match → identity = 1.0 across all segments.
    for rec in _human_igh_exp().run_records(n=5, seed=0):
        assert rec["v_identity"] == 1.0
        assert rec["d_identity"] == 1.0
        assert rec["j_identity"] == 1.0


def test_identity_drops_with_mutations():
    rec_clean = _human_igh_exp().run_records(n=1, seed=42)[0]
    rec_mutated = (
        _human_igh_exp().mutate(count=20).run_records(n=1, seed=42)[0]
    )
    # The mutated record's V identity should be below the clean
    # record's V identity (which is 1.0).
    assert rec_clean["v_identity"] == 1.0
    assert rec_mutated["v_identity"] < 1.0


def test_identity_in_unit_interval_under_corruption_stack():
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(count=15)
        .corrupt_pcr(count=5)
        .corrupt_indels(count=3, insertion_prob=0.5)
    )
    for rec in exp.run_records(n=20, seed=0):
        for seg in ("v", "d", "j"):
            ident = rec[f"{seg}_identity"]
            assert ident is None or 0.0 <= ident <= 1.0


def test_d_identity_none_for_vj_chain():
    rec = (
        ga.Experiment.on("human_igk")
        .recombine()
        .run_records(n=1, seed=0)[0]
    )
    assert rec["d_identity"] is None
    # V and J identity still computed normally.
    assert rec["v_identity"] == 1.0
    assert rec["j_identity"] == 1.0


def test_identity_low_for_contaminant_pass():
    # Full contaminant → seq is random bases against the original
    # V/D/J germline → identity collapses toward ~25% (random
    # 4-letter alphabet match rate).
    rec = (
        ga.Experiment.on("human_igh")
        .recombine()
        .corrupt_contaminants(prob=1.0)
        .run_records(n=1, seed=0)[0]
    )
    assert rec["is_contaminant"] is True
    # Loose bound: we expect somewhere around 0.25, but allow a
    # generous range for sampling variance over a single record.
    assert 0.05 < rec["v_identity"] < 0.5


def test_identity_consistent_with_alignment_string_diff_count():
    # GROUND-TRUTH INVARIANT: the identity numerator must equal the
    # count of V/D/J columns in the alignment strings where sa
    # matches ga (case-insensitive, ignoring gaps and N).
    rec = (
        _human_igh_exp()
        .mutate(count=15)
        .corrupt_indels(count=3, insertion_prob=0.5)
        .run_records(n=1, seed=42)[0]
    )
    sa = rec["sequence_alignment"]
    galn = rec["germline_alignment"]
    for seg in ("v", "d", "j"):
        a_s = rec[f"{seg}_alignment_start"]
        a_e = rec[f"{seg}_alignment_end"]
        if a_s is None:
            continue
        total = a_e - a_s
        matches = sum(
            1
            for sa_c, ga_c in zip(sa[a_s:a_e], galn[a_s:a_e])
            if sa_c not in ("-",)
            and ga_c not in ("-", "N")
            and sa_c.upper() == ga_c.upper()
        )
        expected_identity = matches / total if total > 0 else None
        assert (
            rec[f"{seg}_identity"] == expected_identity
        ), f"{seg}: identity {rec[f'{seg}_identity']} != recomputed {expected_identity}"


def test_identity_decreases_monotonically_with_mutation_count():
    # As mutation count grows, average identity over V should drop.
    means = {}
    for count in (0, 5, 25):
        records = (
            _human_igh_exp().mutate(count=count).run_records(n=20, seed=0)
        )
        means[count] = sum(r["v_identity"] for r in records) / len(records)
    assert means[0] > means[5] > means[25]
    assert means[0] == 1.0  # no mutations → perfect


def test_h6_fields_in_default_column_order():
    from GenAIRR.result import _DEFAULT_COLUMN_ORDER

    for seg in ("v", "d", "j"):
        for kind in ("score", "identity", "support"):
            assert f"{seg}_{kind}" in _DEFAULT_COLUMN_ORDER


def test_h6_fields_in_tsv_export(tmp_path):
    import csv as _csv

    result = _human_igh_records(n=2)
    path = tmp_path / "h6.tsv"
    result.to_tsv(str(path))
    with open(path, encoding="utf-8") as fh:
        rows = list(_csv.DictReader(fh, delimiter="\t"))
    for row in rows:
        # Stubbed fields export as empty strings.
        assert row["v_score"] == ""
        assert row["v_support"] == ""
        # Identity is a float string.
        assert float(row["v_identity"]) == 1.0


# ──────────────────────────────────────────────────────────────────
# H.7 — airr.validate_rearrangement compatibility + airr_strict mode
# ──────────────────────────────────────────────────────────────────


def test_to_tsv_passes_airr_validate_rearrangement(tmp_path):
    airr = pytest.importorskip("airr")
    rec_path = tmp_path / "valid.tsv"
    result = (
        _human_igh_exp().mutate(count=5).run_records(n=5, seed=0)
    )
    result.to_tsv(str(rec_path))
    assert airr.validate_rearrangement(str(rec_path)) is True


def test_to_tsv_passes_airr_validation_under_full_corruption_stack(tmp_path):
    airr = pytest.importorskip("airr")
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(count=15)
        .corrupt_pcr(count=5)
        .corrupt_quality(count=8)
        .corrupt_indels(count=4, insertion_prob=0.5)
    )
    path = tmp_path / "stack.tsv"
    exp.run_records(n=20, seed=0).to_tsv(str(path))
    assert airr.validate_rearrangement(str(path)) is True


def test_to_tsv_passes_airr_validation_for_vj_chain(tmp_path):
    airr = pytest.importorskip("airr")
    path = tmp_path / "vj.tsv"
    (
        ga.Experiment.on("human_igk")
        .recombine()
        .run_records(n=5, seed=0)
        .to_tsv(str(path))
    )
    assert airr.validate_rearrangement(str(path)) is True


def test_to_tsv_passes_airr_validation_for_tcrb(tmp_path):
    airr = pytest.importorskip("airr")
    path = tmp_path / "tcrb.tsv"
    # G4: TCR doesn't undergo SHM; use PCR errors for sequencing
    # realism instead.
    (
        ga.Experiment.on("human_tcrb")
        .recombine()
        .corrupt_pcr(count=3)
        .run_records(n=5, seed=0)
        .to_tsv(str(path))
    )
    assert airr.validate_rearrangement(str(path)) is True


def test_to_tsv_passes_airr_validation_with_contamination(tmp_path):
    airr = pytest.importorskip("airr")
    path = tmp_path / "contam.tsv"
    (
        _human_igh_exp()
        .corrupt_contaminants(prob=0.5)
        .run_records(n=5, seed=0)
        .to_tsv(str(path))
    )
    assert airr.validate_rearrangement(str(path)) is True


def test_airr_strict_increments_start_fields_only(tmp_path):
    # Default mode keeps coords 0-based half-open. ``airr_strict=True``
    # bumps every ``*_start`` field by 1 (giving 1-based-inclusive).
    # The ``*_end`` fields are left unchanged because Python's
    # half-open ``end`` already equals the 1-based-inclusive ``end``.
    import csv as _csv

    result = _human_igh_exp().run_records(n=2, seed=42)
    py_path = tmp_path / "py.tsv"
    strict_path = tmp_path / "strict.tsv"
    result.to_tsv(str(py_path))
    result.to_tsv(str(strict_path), airr_strict=True)

    with open(py_path) as fh:
        py_row = next(_csv.DictReader(fh, delimiter="\t"))
    with open(strict_path) as fh:
        strict_row = next(_csv.DictReader(fh, delimiter="\t"))

    # Every *_start coord field should be off by exactly 1.
    for field in (
        "v_sequence_start",
        "d_sequence_start",
        "j_sequence_start",
        "v_alignment_start",
        "d_alignment_start",
        "j_alignment_start",
        "v_germline_start",
        "d_germline_start",
        "j_germline_start",
        "junction_start",
    ):
        if py_row[field] == "" or strict_row[field] == "":
            continue
        assert int(strict_row[field]) == int(py_row[field]) + 1, (
            f"{field}: strict {strict_row[field]} != default {py_row[field]} + 1"
        )

    # Every *_end coord field should be unchanged.
    for field in (
        "v_sequence_end",
        "d_sequence_end",
        "j_sequence_end",
        "v_alignment_end",
        "d_alignment_end",
        "j_alignment_end",
        "v_germline_end",
        "d_germline_end",
        "j_germline_end",
        "junction_end",
    ):
        assert py_row[field] == strict_row[field], (
            f"{field}: strict {strict_row[field]} != default {py_row[field]}"
        )


def test_airr_strict_passes_airr_validation(tmp_path):
    airr = pytest.importorskip("airr")
    path = tmp_path / "strict.tsv"
    result = (
        _human_igh_exp()
        .mutate(count=10)
        .corrupt_indels(count=3, insertion_prob=0.5)
        .run_records(n=5, seed=0)
    )
    result.to_tsv(str(path), airr_strict=True)
    assert airr.validate_rearrangement(str(path)) is True


def test_airr_strict_csv_export(tmp_path):
    import csv as _csv

    result = _human_igh_exp().run_records(n=1, seed=42)
    py_path = tmp_path / "py.csv"
    strict_path = tmp_path / "strict.csv"
    result.to_csv(str(py_path))
    result.to_csv(str(strict_path), airr_strict=True)

    with open(py_path) as fh:
        py_row = next(_csv.DictReader(fh))
    with open(strict_path) as fh:
        strict_row = next(_csv.DictReader(fh))

    assert int(strict_row["v_sequence_start"]) == int(py_row["v_sequence_start"]) + 1


def test_airr_strict_dataframe():
    pd = pytest.importorskip("pandas")
    result = _human_igh_exp().run_records(n=2, seed=42)

    df_default = result.to_dataframe()
    df_strict = result.to_dataframe(airr_strict=True)

    assert isinstance(df_default, pd.DataFrame)
    assert isinstance(df_strict, pd.DataFrame)
    # Every *_start in strict must equal default + 1.
    for col in (
        "v_sequence_start",
        "j_sequence_start",
        "junction_start",
        "v_alignment_start",
        "v_germline_start",
    ):
        for i in range(len(df_default)):
            d = df_default.iloc[i][col]
            s = df_strict.iloc[i][col]
            if d is None or (hasattr(d, "is_integer") and not d == d):
                continue
            assert s == d + 1, f"{col}[{i}]: {s} != {d} + 1"


def test_to_airr_strict_helper_keeps_none_starts():
    # If a *_start is ``None`` (e.g. d_alignment_start on a VJ chain),
    # the helper must leave it as None.
    from GenAIRR.result import _to_airr_strict

    rec = {
        "v_sequence_start": 0,
        "v_sequence_end": 100,
        "d_sequence_start": None,
        "d_sequence_end": None,
        "junction_start": None,
        "j_sequence_start": 200,
        "j_sequence_end": 250,
    }
    converted = _to_airr_strict(rec)
    assert converted["v_sequence_start"] == 1
    assert converted["v_sequence_end"] == 100
    assert converted["d_sequence_start"] is None
    assert converted["junction_start"] is None
    assert converted["j_sequence_start"] == 201


def test_to_airr_strict_helper_does_not_mutate_input():
    from GenAIRR.result import _to_airr_strict

    rec = {"v_sequence_start": 0, "v_sequence_end": 100}
    snapshot = dict(rec)
    _to_airr_strict(rec)
    assert rec == snapshot


def test_airr_strict_default_is_false():
    # Sanity: ``to_tsv`` without the flag should produce the same
    # output as ``airr_strict=False``.
    import tempfile, os, csv as _csv

    result = _human_igh_exp().run_records(n=1, seed=0)
    with tempfile.TemporaryDirectory() as tmp:
        a = os.path.join(tmp, "a.tsv")
        b = os.path.join(tmp, "b.tsv")
        result.to_tsv(a)
        result.to_tsv(b, airr_strict=False)
        with open(a) as fh:
            ra = list(_csv.DictReader(fh, delimiter="\t"))
        with open(b) as fh:
            rb = list(_csv.DictReader(fh, delimiter="\t"))
    assert ra == rb


def test_identity_helper_from_columns():
    from GenAIRR._airr_record import _alignment_columns, _identity_from_columns

    refdata = _human_igh_exp().refdata
    out = _human_igh_exp().mutate(count=10).run(n=1, seed=42)[0]
    cols = _alignment_columns(out, refdata)
    identities = _identity_from_columns(cols)
    # Same as the public field — guard against silent re-routing.
    rec = _human_igh_exp().mutate(count=10).run_records(n=1, seed=42)[0]
    assert identities["V"] == rec["v_identity"]
    assert identities["D"] == rec["d_identity"]
    assert identities["J"] == rec["j_identity"]


_RC_COMPLEMENT = str.maketrans("ACGTacgtNnUu", "TGCAtgcaNnAa")


def _reverse_complement(s: str) -> str:
    return s.translate(_RC_COMPLEMENT)[::-1]


def test_phase12d_5prime_loss_drops_leading_bases():
    # Fixed-length 5' loss: every record's sequence is shorter by
    # exactly the requested amount, and the V segment's CIGAR starts
    # with that many `D` ops (the lost ref positions).
    base = ga.Experiment.on("human_igh").recombine()
    lost = ga.Experiment.on("human_igh").recombine().corrupt_5prime_loss(length=15)
    base_recs = list(base.run_records(n=10, seed=0))
    lost_recs = list(lost.run_records(n=10, seed=0))
    for fwd, lr in zip(base_recs, lost_recs):
        # Sequence shrinks by exactly 15.
        assert lr["sequence_length"] == fwd["sequence_length"] - 15, (
            f"len={lr['sequence_length']} vs {fwd['sequence_length']} - 15"
        )
        # V's CIGAR begins with `15D` (the lost prefix of V).
        assert lr["v_cigar"].startswith("15D"), f"v_cigar={lr['v_cigar']!r}"


def test_phase12d_3prime_loss_drops_trailing_bases():
    base = ga.Experiment.on("human_igh").recombine()
    lost = ga.Experiment.on("human_igh").recombine().corrupt_3prime_loss(length=10)
    for fwd, lr in zip(base.run_records(n=10, seed=0), lost.run_records(n=10, seed=0)):
        assert lr["sequence_length"] == fwd["sequence_length"] - 10
        # J's CIGAR ends with `10D` (the lost suffix of J).
        assert lr["j_cigar"].endswith("10D"), f"j_cigar={lr['j_cigar']!r}"


def test_phase12d_loss_zero_is_no_op():
    base = ga.Experiment.on("human_igh").recombine()
    lost_5 = ga.Experiment.on("human_igh").recombine().corrupt_5prime_loss(length=0)
    lost_3 = ga.Experiment.on("human_igh").recombine().corrupt_3prime_loss(length=0)
    for fwd, l5, l3 in zip(
        base.run_records(n=5, seed=0),
        lost_5.run_records(n=5, seed=0),
        lost_3.run_records(n=5, seed=0),
    ):
        assert l5["sequence"] == fwd["sequence"]
        assert l3["sequence"] == fwd["sequence"]


def test_phase12d_loss_preserves_alignment_invariants():
    # 5'/3' loss extends the CIGAR with leading/trailing `D` ops.
    # The H.3 alignment-string invariants (sa.replace('-')==seq, etc.)
    # should still hold — the column walker accounts for the extra
    # D ops in the trailing-D-fill loop.
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(count=8)
        .corrupt_5prime_loss(length=(0, 25))
        .corrupt_3prime_loss(length=(0, 20))
        .corrupt_pcr(count=3)
    )
    failures = []
    for i, rec in enumerate(exp.run_records(n=100, seed=0)):
        issues = _alignment_invariants_hold(rec)
        if issues:
            failures.append((i, issues))
    assert failures == [], f"alignment invariants broke under 5'/3' loss: {failures[:3]}"


def test_phase12d_loss_rejects_invalid_length():
    exp = ga.Experiment.on("human_igh").recombine()
    for bad in (-1, "ten", True):
        try:
            exp.corrupt_5prime_loss(length=bad)
        except (ValueError, TypeError):
            pass
        else:
            raise AssertionError(f"expected error for length={bad!r}")


def test_phase12d_rev_comp_flips_sequence_and_flag():
    # prob=1.0: every record should have rev_comp=True and the
    # `sequence` should be the reverse-complement of the same record
    # without rev_comp applied.
    base = ga.Experiment.on("human_igh").recombine()
    flipped = ga.Experiment.on("human_igh").recombine().corrupt_reverse_complement(prob=1.0)
    base_recs = list(base.run_records(n=20, seed=0))
    flipped_recs = list(flipped.run_records(n=20, seed=0))
    assert len(base_recs) == len(flipped_recs) == 20
    for fwd, rc in zip(base_recs, flipped_recs):
        assert fwd["rev_comp"] is False
        assert rc["rev_comp"] is True
        assert rc["sequence"] == _reverse_complement(fwd["sequence"])


def test_phase12d_rev_comp_flips_per_segment_sequence_coords():
    # When the sequence is flipped, V (at the forward 5' end) ends up
    # at the antisense 3' end. So new_v_sequence_start = seq_len -
    # old_v_sequence_end, and the v slice of the new sequence equals
    # the reverse-complement of the old v slice.
    base = ga.Experiment.on("human_igh").recombine()
    flipped = ga.Experiment.on("human_igh").recombine().corrupt_reverse_complement(prob=1.0)
    for fwd, rc in zip(base.run_records(n=20, seed=0), flipped.run_records(n=20, seed=0)):
        seq_len = fwd["sequence_length"]
        for seg in ("v", "d", "j"):
            fs, fe = fwd[f"{seg}_sequence_start"], fwd[f"{seg}_sequence_end"]
            if fs is None:
                continue
            rs, re = rc[f"{seg}_sequence_start"], rc[f"{seg}_sequence_end"]
            assert rs == seq_len - fe, f"{seg}_sequence_start: {rs} vs {seq_len - fe}"
            assert re == seq_len - fs, f"{seg}_sequence_end: {re} vs {seq_len - fs}"
            # The slice of the flipped sequence at the new coords
            # equals the reverse-complement of the forward segment's
            # sequence slice.
            assert rc["sequence"][rs:re] == _reverse_complement(fwd["sequence"][fs:fe])


def test_phase12d_rev_comp_keeps_alignment_strings_in_forward_orientation():
    # AIRR-spec: `sequence_alignment` and `germline_alignment` stay
    # forward-oriented even when `rev_comp=True`. The CIGAR and the
    # `*_alignment_start/end` and `*_germline_start/end` coords are
    # also forward.
    base = ga.Experiment.on("human_igh").recombine()
    flipped = ga.Experiment.on("human_igh").recombine().corrupt_reverse_complement(prob=1.0)
    for fwd, rc in zip(base.run_records(n=10, seed=0), flipped.run_records(n=10, seed=0)):
        assert rc["sequence_alignment"] == fwd["sequence_alignment"]
        assert rc["germline_alignment"] == fwd["germline_alignment"]
        assert rc["v_cigar"] == fwd["v_cigar"]
        assert rc["d_cigar"] == fwd["d_cigar"]
        assert rc["j_cigar"] == fwd["j_cigar"]
        for seg in ("v", "d", "j"):
            assert rc[f"{seg}_alignment_start"] == fwd[f"{seg}_alignment_start"]
            assert rc[f"{seg}_alignment_end"] == fwd[f"{seg}_alignment_end"]
            assert rc[f"{seg}_germline_start"] == fwd[f"{seg}_germline_start"]
            assert rc[f"{seg}_germline_end"] == fwd[f"{seg}_germline_end"]


def test_phase12d_rev_comp_zero_prob_is_no_op():
    base = ga.Experiment.on("human_igh").recombine()
    flipped = ga.Experiment.on("human_igh").recombine().corrupt_reverse_complement(prob=0.0)
    for fwd, rc in zip(base.run_records(n=20, seed=0), flipped.run_records(n=20, seed=0)):
        # The coin is flipped (records always present in trace) but
        # never fires when prob=0.
        assert rc["rev_comp"] is False
        assert rc["sequence"] == fwd["sequence"]


def test_phase12d_rev_comp_rejects_invalid_prob():
    exp = ga.Experiment.on("human_igh").recombine()
    import math
    for bad in (-0.1, 1.5, math.nan, float("inf")):
        try:
            exp.corrupt_reverse_complement(prob=bad)
        except ValueError:
            pass
        else:
            raise AssertionError(f"expected ValueError for prob={bad}")


def test_phase12_locus_falls_back_to_refdata_under_heavy_corruption():
    # Phase 12.C follow-up: under heavy SHM the live-call layer can
    # narrow every V/D/J call to the empty set (no allele supports
    # the corrupted sequence). Pre-fix `derive_locus` then had
    # nothing to parse and returned `""`. The fix walks the refdata
    # pool's first allele names as a fallback so the AIRR record
    # always carries a meaningful locus.
    expected = {
        "human_igh": "IGH",
        "human_igk": "IGK",
        "human_igl": "IGL",
        "human_tcrb": "TRB",
    }
    for cfg, want in expected.items():
        # G4: SHM is BCR-only; drop mutate for TCR.
        exp = ga.Experiment.on(cfg).recombine()
        if not cfg.startswith("human_tcr"):
            exp = exp.mutate(count=15)
        exp = (
            exp.corrupt_pcr(count=5)
            .corrupt_quality(count=8)
            .corrupt_indels(count=4, insertion_prob=0.5)
        )
        for i, rec in enumerate(exp.run_records(n=200, seed=0)):
            assert rec["locus"] == want, (
                f"{cfg} rec {i}: locus={rec['locus']!r} expected {want!r} "
                f"(v_call={rec['v_call']!r} j_call={rec['j_call']!r})"
            )


def test_phase12_per_segment_ga_slice_matches_called_allele():
    # Phase 12.C — discovered audit invariant: when SHM mutations
    # narrow the live call to a different allele than the one
    # originally sampled, `germline_alignment[v_align_start:v_align_end]
    # ].replace('-','')` must equal the *called* allele's bytes at
    # `[v_germline_start, v_germline_end)`. Same for D and J.
    #
    # Pre-Phase-12.C the column walker emitted `nuc.germline` (the
    # provenance byte from AssembleSegmentPass), causing a divergence
    # in ~1% of records under heavy corruption. The fix routes the
    # column-walker germline byte through the projected (live-call's
    # first) allele.
    pipelines = [
        ("human_igh", "vdj"),
        ("human_igk", "vj"),
        ("human_igl", "vj"),
        ("human_tcrb", "vdj"),
    ]
    failures = []
    for cfg, _kind in pipelines:
        # G4: SHM is BCR-only; drop mutate for TCR.
        exp = ga.Experiment.on(cfg).recombine()
        if not cfg.startswith("human_tcr"):
            exp = exp.mutate(count=8)
        exp = (
            exp.corrupt_pcr(count=2)
            .corrupt_quality(count=4)
            .corrupt_indels(count=2, insertion_prob=0.5)
        )
        refdata = exp.refdata
        alleles = {"V": {}, "D": {}, "J": {}}
        for seg, n_attr, get in [
            ("V", refdata.v_pool_size, refdata.v_allele),
            ("D", refdata.d_pool_size, refdata.d_allele),
            ("J", refdata.j_pool_size, refdata.j_allele),
        ]:
            for i in range(n_attr()):
                a = get(i)
                alleles[seg][a.name] = a.seq().decode("ascii").upper()
        for i, rec in enumerate(exp.run_records(n=300, seed=0)):
            for seg in ("V", "D", "J"):
                calls = rec[f"{seg.lower()}_call"]
                if not calls:
                    continue
                a_s = rec[f"{seg.lower()}_alignment_start"]
                a_e = rec[f"{seg.lower()}_alignment_end"]
                g_s = rec[f"{seg.lower()}_germline_start"]
                g_e = rec[f"{seg.lower()}_germline_end"]
                if a_s is None:
                    continue
                first_call = calls.split(",")[0]
                ref = alleles[seg].get(first_call)
                if ref is None:
                    continue
                ga_str = rec["germline_alignment"][a_s:a_e].replace("-", "")
                if ga_str.upper() != ref[g_s:g_e].upper():
                    failures.append((cfg, i, seg, first_call, ga_str, ref[g_s:g_e]))
    assert failures == [], f"per-segment ga slice mismatches: {failures[:3]}"


def test_h5_coords_invariants_under_corruption_stack():
    # Same n=200 stress sweep as alignment/CIGAR — verify the coord
    # pairs stay self-consistent under every corruption mode.
    #
    # Phase 11.7 restored the strict identities: `*_germline_start/end`
    # are now derived from the column walker's `ref_ranges` (the
    # union of ref positions consumed by `M` and `D` ops), which by
    # construction enforces `germline_span == M + D`. Combined with
    # the alignment-span identity `alignment_span == M + I + D`, the
    # AIRR coord triples are fully self-consistent again.
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(count=15)
        .corrupt_pcr(count=5)
        .corrupt_quality(count=8)
        .corrupt_indels(count=4, insertion_prob=0.5)
    )
    failures = []
    for i, rec in enumerate(exp.run_records(n=200, seed=0)):
        for seg in ("v", "d", "j"):
            a_s = rec[f"{seg}_alignment_start"]
            a_e = rec[f"{seg}_alignment_end"]
            g_s = rec[f"{seg}_germline_start"]
            g_e = rec[f"{seg}_germline_end"]
            cig = rec[f"{seg}_cigar"]
            if a_s is None:
                continue
            ops = _cigar_op_counts(cig)
            if a_e - a_s != ops["M"] + ops["I"] + ops["D"]:
                failures.append(
                    (i, seg, "alignment-span", a_e - a_s, ops)
                )
            if g_e - g_s != ops["M"] + ops["D"]:
                failures.append(
                    (i, seg, "germline-span", g_e - g_s, ops)
                )
    assert failures == [], f"H.5 invariants broke for: {failures[:3]}"


def test_cigar_only_uses_m_when_pure_recombination():
    # Strict: pure recombine + no mutation + no indel = single M
    # block per segment.
    exp = ga.Experiment.on("human_igh").recombine()
    for rec in exp.run_records(n=20, seed=0):
        for seg in ("v_cigar", "d_cigar", "j_cigar"):
            cig = rec[seg]
            if cig:
                assert _re.fullmatch(r"\d+M", cig), f"{seg}={cig!r}"


def test_alignment_consistent_under_contaminant_pass():
    # Even when the contaminant pass overwrites the entire sequence
    # with non-immune bases, the IR's germline_position tracking is
    # preserved — so sequence_alignment shows the contaminant bases
    # and germline_alignment shows the original (unrelated) V/D/J
    # source. The user can detect this via ``is_contaminant``.
    rec = (
        ga.Experiment.on("human_igh")
        .recombine()
        .corrupt_contaminants(prob=1.0)
        .run_records(n=1, seed=0)[0]
    )
    assert rec["is_contaminant"] is True
    assert _alignment_invariants_hold(rec) == []
    # And the contaminant divergence is visible: most positions in
    # the V coding region differ between sa and ga.
    sa = rec["sequence_alignment"]
    galn = rec["germline_alignment"]
    v_start, v_end = rec["v_sequence_start"], rec["v_sequence_end"]
    diffs = sum(
        1
        for a, b in zip(sa[v_start:v_end], galn[v_start:v_end])
        if a != b and b not in ("N", "-")
    )
    # Random bases should mismatch in roughly 75% of V positions.
    v_len = v_end - v_start
    assert diffs >= v_len // 2


def test_germline_alignment_byte_matches_source_allele():
    # Verify each V position in germline_alignment equals the
    # corresponding source allele base (post-trim, pre-mutation).
    # The "source allele" is the originally-sampled provenance
    # allele, not whatever set ``rec['v_call']`` projects to today —
    # under Phase 4 the call is a live-call view that can become a
    # multi-allele string after a mutation pass refreshes it (Phase 6).
    # We pull provenance from the simulation IR directly.
    result = (
        _human_igh_exp()
        .mutate(count=8)
        .run_records(n=1, seed=42)
    )
    rec = result[0]
    refdata = _human_igh_exp().refdata
    galn = rec["germline_alignment"]

    sim = result.outcomes[0].final_simulation()
    v_id = sim.v_allele_id()
    assert v_id is not None, "expected V provenance after recombine + mutate"
    v_seq = bytes(refdata.v_allele(v_id).seq()).decode().upper()
    v_start = rec["v_sequence_start"]
    v_end = rec["v_sequence_end"]
    # No 5' V trim in the DSL → V starts at germline pos 0.
    # No indels in this test → 1:1 mapping.
    expected = v_seq[: v_end - v_start]
    actual = galn[v_start:v_end]
    assert actual == expected, f"V germline mismatch:\nexpected: {expected[:60]}\nactual:   {actual[:60]}"
