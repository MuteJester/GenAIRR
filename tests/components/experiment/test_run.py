"""Running a compiled experiment.

Iteration semantics of run(n, seed), repeated runs of one compiled object, and
the bridge that turns a pass failure into a strict sampling error.
"""
from __future__ import annotations
from GenAIRR import _engine as ge
import pytest
from GenAIRR import CompiledExperiment, Experiment
from tests.support.fixtures.experiment_dsl import _vdj_refdata, _vj_refdata

pytestmark = pytest.mark.behavior

def test_recombine_custom_np1_length_distribution_propagates():
    # Force NP1 = exactly 4 every time. The trace's NP1 length and
    # the assembled region's length must both be 4.
    cfg = _vj_refdata()
    outcomes = Experiment.on(cfg).allow_curatable_refdata().recombine(np1_lengths=[(4, 1.0)]).run(
        n=10, seed=0
    )
    for o in outcomes:
        np1 = o.final_simulation().regions()[1]
        assert len(np1) == 4
        assert o.trace().find("np.np1.length").value == 4


def test_recombine_custom_np2_length_rejected_on_vj():
    # v2.0: passing np2_lengths on a VJ chain is now an explicit
    # ValueError (was silently ignored in v1.x). VJ chains have no
    # NP2 region; the kwarg can only ever indicate a user mistake.
    cfg = _vj_refdata()
    with pytest.raises(ValueError, match="np2_lengths"):
        Experiment.on(cfg).allow_curatable_refdata().recombine(
            np1_lengths=[(2, 1.0)],
            np2_lengths=[(7, 1.0)],
        )


def test_recombine_vj_without_np2_still_produces_single_np():
    # The other half of the old test: VJ chains still produce just
    # V / NP1 / J when np2_lengths is omitted (which is the only
    # supported VJ shape).
    cfg = _vj_refdata()
    outcomes = Experiment.on(cfg).allow_curatable_refdata().recombine(np1_lengths=[(2, 1.0)]).trim(enabled=False).run(n=1, seed=0)
    np1 = outcomes[0].final_simulation().regions()[1]
    assert len(np1) == 2
    # Only NP1 region — VJ has no NP2.
    assert outcomes[0].final_simulation().region_count() == 3


def test_recombine_custom_np2_length_used_on_vdj():
    # VDJ experiment: NP2 length pinned to 5.
    cfg = _vdj_refdata()
    outcomes = Experiment.on(cfg).allow_curatable_refdata().recombine(
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
        Experiment.on(cfg).allow_curatable_refdata().recombine(np1_lengths=[])


def test_run_with_n_one_returns_single_outcome():
    outcomes = Experiment.on(_vj_refdata()).allow_curatable_refdata().recombine().run(n=1, seed=42)
    assert len(outcomes) == 1
    assert outcomes[0].pass_names()[0] == "sample_allele.v"


def test_run_returns_n_outcomes_with_offset_seeds():
    # Per-run seed is `seed + i`. Comparing to direct engine calls
    # at those exact seeds must produce byte-identical outcomes.
    cfg = _vj_refdata()
    outcomes = Experiment.on(cfg).allow_curatable_refdata().recombine().run(n=4, seed=100)
    assert len(outcomes) == 4

    # Build the equivalent plan manually and compare.
    plan = ge.PassPlan()
    plan.push_sample_allele("V", cfg)
    plan.push_sample_allele("J", cfg)
    plan.push_assemble("V")
    plan.push_generate_np("NP1", [(i, 1.0) for i in range(7)])
    plan.push_assemble("J")

    for i, dsl_outcome in enumerate(outcomes):
        baseline = ge.run(plan, seed=100 + i, refdata=cfg, allow_curatable_refdata=True)
        assert (
            dsl_outcome.final_simulation().bases()
            == baseline.final_simulation().bases()
        )


def test_run_default_seed_is_zero():
    a = Experiment.on(_vj_refdata()).allow_curatable_refdata().recombine().run(n=1)
    b = Experiment.on(_vj_refdata()).allow_curatable_refdata().recombine().run(n=1, seed=0)
    assert a[0].final_simulation().bases() == b[0].final_simulation().bases()


def test_run_default_n_is_one():
    outcomes = Experiment.on(_vj_refdata()).allow_curatable_refdata().recombine().run()
    assert len(outcomes) == 1


def test_run_rejects_n_below_one():
    cfg = _vj_refdata()
    exp = Experiment.on(cfg).allow_curatable_refdata().recombine()
    with pytest.raises(ValueError, match="n must be at least 1"):
        exp.run(n=0)
    with pytest.raises(ValueError, match="n must be at least 1"):
        exp.run(n=-3)


def test_run_is_deterministic_under_same_seed():
    a = Experiment.on(_vj_refdata()).allow_curatable_refdata().recombine().run(n=10, seed=0xCAFE)
    b = Experiment.on(_vj_refdata()).allow_curatable_refdata().recombine().run(n=10, seed=0xCAFE)
    for x, y in zip(a, b):
        assert x.final_simulation().bases() == y.final_simulation().bases()


def test_run_produces_distinct_outcomes_across_iterations():
    # 50 iterations × 7 NP-length values × randomized bases →
    # overwhelmingly likely that not all outcomes are byte-identical.
    cfg = _vj_refdata()
    cfg.add_v_allele("v2*01", "v2", b"GGGAAACCC", anchor=6)
    outcomes = Experiment.on(cfg).allow_curatable_refdata().recombine().run(n=50, seed=0)
    base_strings = [o.final_simulation().bases() for o in outcomes]
    assert len(set(base_strings)) > 1


def test_compiled_experiment_can_be_run_multiple_times():
    compiled = Experiment.on(_vj_refdata(canonical=True)).allow_curatable_refdata().recombine().compile(allow_curatable_refdata=True)
    a = compiled.run(n=2, seed=0)
    b = compiled.run(n=2, seed=10)
    # Different seeds → different outcomes. Same plan instance
    # (used by both run() calls) → no leftover state.
    assert a[0].final_simulation().bases() != b[0].final_simulation().bases()
    # And re-running with the original seed reproduces the original.
    again = compiled.run(n=2, seed=0)
    assert again[0].final_simulation().bases() == a[0].final_simulation().bases()


def test_productive_is_reexported_from_top_level():
    import GenAIRR

    assert hasattr(GenAIRR, "productive")
    assert callable(GenAIRR.productive)
    assert GenAIRR.productive is ge.productive


def test_productive_returns_contract_set_with_five_named_contracts():
    cs = ge.productive()
    assert isinstance(cs, ge.ContractSet)
    assert len(cs) == 5
    assert not cs.is_empty()
    names = cs.names()
    assert names == [
        "productive_junction_frame",
        "no_stop_codon_in_junction",
        "canonical_anchor.v1.v",
        "canonical_anchor.v1.j",
        "coding_sequence_integrity",
    ]


def test_contract_set_repr_includes_count_and_names():
    r = repr(ge.productive())
    assert "ContractSet" in r
    assert "len=5" in r
    assert "productive_junction_frame" in r


def test_run_without_productive_only_matches_uncalled_chain():
    # Calling productive_only() then immediately compiling without
    # other constraint methods should match a chain where it was
    # never called *only if* the productive() outcome happens to
    # equal the unconstrained outcome at the chosen seed. We instead
    # verify the basic invariant: an experiment without productive_only
    # produces the same output across two runs at the same seed.
    cfg = _vj_refdata()
    exp = Experiment.on(cfg).allow_curatable_refdata().recombine()
    a = exp.run(n=3, seed=0)
    b = exp.run(n=3, seed=0)
    for x, y in zip(a, b):
        assert x.final_simulation().bases() == y.final_simulation().bases()


def test_productive_only_constrains_np1_length_to_in_frame():
    # Synthetic VJ refdata: V_anchor_to_end = 3, J_anchor_to_W3 = 3,
    # so junction = 6 + NP1. Productive frame ⇒ NP1 % 3 == 0.
    # Filter NP1 length distribution to {0..6} uniform; with productive
    # active every NP1 length must be in {0, 3, 6}.
    cfg = _vj_refdata(canonical=True)
    outcomes = (
        Experiment.on(cfg).allow_curatable_refdata().recombine().productive_only().run(n=20, seed=0)
    )
    for o in outcomes:
        np1 = o.final_simulation().regions()[1]
        assert len(np1) % 3 == 0, (
            f"NP1 length {len(np1)} should be divisible by 3 under productive"
        )


def test_run_unconstrained_produces_some_out_of_frame_np1():
    # Negative control for the test above: without productive_only(),
    # the uniform NP1 distribution produces out-of-frame lengths
    # ~67% of the time. Across 20 seeds we should see at least one.
    cfg = _vj_refdata()
    outcomes = Experiment.on(cfg).allow_curatable_refdata().recombine().run(n=20, seed=0)
    out_of_frame = [
        o for o in outcomes if len(o.final_simulation().regions()[1]) % 3 != 0
    ]
    assert out_of_frame, "expected at least one out-of-frame NP1 without productive_only"


def test_productive_only_is_idempotent():
    # Calling .productive_only() twice is the same as once — the
    # bundle is added at most once to the contract list.
    cfg = _vj_refdata(canonical=True)
    a = (
        Experiment.on(cfg).allow_curatable_refdata().recombine().productive_only().run(n=2, seed=0)
    )
    b = (
        Experiment.on(cfg).allow_curatable_refdata()
        .recombine()
        .productive_only()
        .productive_only()
        .run(n=2, seed=0)
    )
    for x, y in zip(a, b):
        assert x.final_simulation().bases() == y.final_simulation().bases()


def test_productive_only_is_order_independent_in_the_chain():
    # Whether productive_only() is called before or after recombine()
    # / mutate(), the compiled simulator is structurally equal.
    cfg = _vj_refdata(canonical=True)
    a = (
        Experiment.on(cfg).allow_curatable_refdata()
        .recombine()
        .mutate(count=3)
        .productive_only()
        .run(n=2, seed=0)
    )
    b = (
        Experiment.on(cfg).allow_curatable_refdata()
        .productive_only()
        .recombine()
        .mutate(count=3)
        .run(n=2, seed=0)
    )
    for x, y in zip(a, b):
        assert x.final_simulation().bases() == y.final_simulation().bases()


def test_compiled_experiment_carries_productive_bundle():
    # Contracts are captured at compile time and reused by run().
    cfg = _vj_refdata(canonical=True)
    compiled = Experiment.on(cfg).allow_curatable_refdata().recombine().productive_only().compile(allow_curatable_refdata=True)
    assert compiled.active_contracts == (
        "productive_junction_frame",
        "no_stop_codon_in_junction",
        "canonical_anchor.v1.v",
        "canonical_anchor.v1.j",
        "coding_sequence_integrity",
    )
    a = compiled.run(n=3, seed=0)
    for o in a:
        assert len(o.final_simulation().regions()[1]) % 3 == 0


def test_run_with_productive_only_is_deterministic_under_same_seed():
    cfg = _vj_refdata(canonical=True)
    exp = Experiment.on(cfg).allow_curatable_refdata().recombine().productive_only()
    a = exp.run(n=5, seed=0xCAFE)
    b = exp.run(n=5, seed=0xCAFE)
    for x, y in zip(a, b):
        assert x.final_simulation().bases() == y.final_simulation().bases()


def _vj_refdata_with_restrictive_np1(*, np1_lengths):
    """Build a VJ refdata + PassPlan whose NP1-length distribution is
    fully under our control. Used by the strict-mode tests so we can
    construct deliberately-unsatisfiable scenarios.
    """
    cfg = _vj_refdata(canonical=True)
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
    cfg.add_v_allele("v1*01", "v1", b"AAACCCTGT", anchor=6)
    cfg.add_j_allele("j_good*01", "j_good", b"TTTAAA", anchor=0)
    cfg.add_j_allele("j_bad*01", "j_bad", b"ATTTAAA", anchor=1)
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
    cfg = _vj_refdata(canonical=True)
    out = (
        Experiment.on(cfg).allow_curatable_refdata()
        .recombine()
        .productive_only()
        .run(n=3, seed=0, strict=True)
    )
    assert len(out) == 3


def test_productive_compile_rejects_unsatisfiable_length_distribution():
    # NP1 distribution {1, 2} has no in-frame mass for this VJ fixture.
    # This is now a build-time D7 failure, independent of strict mode.
    cfg, plan = _vj_refdata_with_restrictive_np1(np1_lengths=[(1, 1.0), (2, 1.0)])
    with pytest.raises(ValueError, match="supports have no productive completion"):
        ge.run(plan, seed=0, refdata=cfg, respect=ge.productive(), strict=True, allow_curatable_refdata=True)


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
                allow_curatable_refdata=True,
            )
            assert out.final_simulation().j_allele_id() == 0


def test_strict_mode_passes_through_compiled_experiment_run():
    cfg = _vj_refdata(canonical=True)
    compiled = Experiment.on(cfg).allow_curatable_refdata().recombine().productive_only().compile(allow_curatable_refdata=True)
    # Satisfiable: succeeds.
    out = compiled.run(n=2, seed=0, strict=True)
    assert len(out) == 2


def test_strict_mode_default_is_false():
    cfg, plan = _vj_refdata_with_runtime_frame_residue()
    out = ge.run(plan, seed=0, refdata=cfg, respect=ge.productive(), allow_curatable_refdata=True)
    assert out.final_simulation().j_allele_id() == 0


def test_strict_without_respect_runs_permissively():
    # strict=True with no contracts is a no-op for failure semantics
    # — there's nothing to fail against. The run proceeds normally.
    cfg = _vj_refdata()
    out = Experiment.on(cfg).allow_curatable_refdata().recombine().run(n=1, seed=0, strict=True)
    assert len(out) == 1


def test_strict_sampling_error_args_form_three_tuple():
    # Public contract: StrictSamplingError.args is always a 3-tuple
    # so callers can destructure deterministically.
    # References and joining support are valid, so compilation succeeds.
    # One structural deletion cannot preserve the two three-base anchors;
    # the complete event support is empty and fails at runtime.
    cfg = ge.RefDataConfig.vj()
    cfg.add_v_allele("v*01", "v", b"TGT", anchor=0)
    cfg.add_j_allele("j*01", "j", b"TGG", anchor=0)
    plan = ge.PassPlan()
    plan.push_sample_allele("V", cfg)
    plan.push_sample_allele("J", cfg)
    plan.push_assemble("V")
    plan.push_generate_np("NP1", [(0, 1.0)])
    plan.push_assemble("J")
    plan.push_corrupt_indel([(1, 1.0)], insertion_prob=0.0)
    with pytest.raises(ge.StrictSamplingError) as err:
        ge.run(plan, seed=0, refdata=cfg, respect=ge.productive(), strict=True, allow_curatable_refdata=True)
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
    a = ge.run(plan, seed=0xC0FFEE, refdata=cfg, respect=ge.productive(), strict=True, allow_curatable_refdata=True)
    b = ge.run(plan, seed=0xC0FFEE, refdata=cfg, respect=ge.productive(), strict=True, allow_curatable_refdata=True)
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
        records = (
            Experiment.on(cfg_name)
            .recombine()
            .productive_only()
            .run_records(n=n, seed=seed, strict=True)
        )
        assert all(rec["productive"] is True for rec in records), cfg_name
