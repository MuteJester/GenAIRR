"""
Tests for GenAIRR seed system reproducibility.

Verifies that set_seed() produces identical results across runs,
using the Experiment DSL.
"""

import pytest

from GenAIRR import (
    set_seed,
    get_seed,
    reset_seed,
    Experiment,
)


# ============================================================
# Seed management API
# ============================================================

class TestSeedSystem:
    """Test seed management functions."""

    def test_set_seed_stores_value(self):
        set_seed(42)
        assert get_seed() == 42
        set_seed(12345)
        assert get_seed() == 12345

    def test_get_seed_returns_none_initially(self):
        reset_seed()
        assert get_seed() is None

    def test_reset_seed_clears_value(self):
        set_seed(42)
        assert get_seed() == 42
        reset_seed()
        assert get_seed() is None


# ============================================================
# Reproducibility via Experiment DSL
# ============================================================

class TestExperimentReproducibility:
    """Test that same seed produces identical results."""

    def test_same_seed_same_sequences(self):
        r1 = Experiment.on("human_igh").run(n=10, seed=42)
        r2 = Experiment.on("human_igh").run(n=10, seed=42)
        for a, b in zip(r1, r2):
            assert a["sequence"] == b["sequence"]
            assert a["v_call"] == b["v_call"]

    def test_different_seeds_different_results(self):
        r1 = Experiment.on("human_igh").run(n=10, seed=42)
        r2 = Experiment.on("human_igh").run(n=10, seed=999)
        n_same = sum(1 for a, b in zip(r1, r2) if a["sequence"] == b["sequence"])
        assert n_same < 10

    def test_shm_reproducibility(self):
        from GenAIRR.ops import rate
        def _run():
            return (Experiment.on("human_igh")
                    .mutate(rate(0.05, 0.10))
                    .run(n=20, seed=42))
        r1 = _run()
        r2 = _run()
        assert [r["sequence"] for r in r1] == [r["sequence"] for r in r2]
        assert [r["mutation_rate"] for r in r1] == [r["mutation_rate"] for r in r2]

    def test_light_chain_reproducibility(self):
        r1 = Experiment.on("human_igk").run(n=10, seed=123)
        r2 = Experiment.on("human_igk").run(n=10, seed=123)
        for a, b in zip(r1, r2):
            assert a["sequence"] == b["sequence"]

    def test_productive_reproducibility(self):
        r1 = Experiment.on("human_igh").run(n=20, seed=42, productive=True)
        r2 = Experiment.on("human_igh").run(n=20, seed=42, productive=True)
        for a, b in zip(r1, r2):
            assert a["sequence"] == b["sequence"]
            assert a["productive"] == b["productive"]

    def test_corruption_reproducibility(self):
        from GenAIRR.ops import with_5prime_loss
        def _run():
            return (Experiment.on("human_igh")
                    .sequence(with_5prime_loss())
                    .run(n=20, seed=42))
        r1 = _run()
        r2 = _run()
        assert [r["sequence"] for r in r1] == [r["sequence"] for r in r2]

    def test_indels_reproducibility(self):
        from GenAIRR.ops import with_indels
        def _run():
            return (Experiment.on("human_igh")
                    .observe(with_indels())
                    .run(n=20, seed=42))
        r1 = _run()
        r2 = _run()
        assert [r["sequence"] for r in r1] == [r["sequence"] for r in r2]

    def test_compiled_simulator_reproducibility(self):
        """Compiled simulator with same seed produces identical results."""
        from GenAIRR.ops import rate
        sim1 = (Experiment.on("human_igh")
                .mutate(rate(0.01, 0.05))
                .compile(seed=42))
        sim2 = (Experiment.on("human_igh")
                .mutate(rate(0.01, 0.05))
                .compile(seed=42))
        r1 = sim1.simulate(n=10)
        r2 = sim2.simulate(n=10)
        for a, b in zip(r1, r2):
            assert a["sequence"] == b["sequence"]


# ============================================================
# Global set_seed flows into the C engine (T0-8 regression)
# ============================================================
#
# Pre-T0-8, set_seed() updated the Python module-global but
# _compile_to_c never consulted it: every Experiment.run()/compile()
# call without an explicit seed= triggered the C engine's auto-seed
# (time + simulator address). The set_seed() docstring claimed
# reproducibility but delivered nondeterminism.
# Post-T0-8: when the caller passes seed=None, _compile_to_c falls
# back to seed.get_seed(); the C simulator is seeded from the global
# value if one was set.

class TestGlobalSeedFallback:

    def setup_method(self):
        reset_seed()

    def teardown_method(self):
        reset_seed()

    def test_global_set_seed_makes_run_reproducible(self):
        """set_seed(N) followed by Experiment.run() (no explicit seed)
        must yield identical sequences across calls."""
        set_seed(42)
        r1 = Experiment.on("human_igh").run(n=10)

        set_seed(42)
        r2 = Experiment.on("human_igh").run(n=10)

        assert [r["sequence"] for r in r1] == [r["sequence"] for r in r2], (
            "set_seed(42) → run() must reproduce; pre-T0-8 the global "
            "seed was ignored by _compile_to_c and runs diverged."
        )

    def test_explicit_seed_overrides_global(self):
        """An explicit seed= on Experiment.run() takes precedence over
        the global set_seed value."""
        set_seed(42)
        r_global = Experiment.on("human_igh").run(n=10)
        r_explicit = Experiment.on("human_igh").run(n=10, seed=99)

        # Same global, but explicit=99 should override and diverge.
        diffs = sum(1 for a, b in zip(r_global, r_explicit)
                    if a["sequence"] != b["sequence"])
        assert diffs >= 8, (
            f"explicit seed=99 produced {diffs}/10 divergent sequences "
            f"vs the global seed=42; expected ≥8"
        )

    def test_compile_also_honors_global_seed(self):
        """The compile() path takes the same fallback as run()."""
        set_seed(123)
        sim_a = Experiment.on("human_igh").compile()
        sim_b = Experiment.on("human_igh").compile()
        # Two simulators, same global seed, same SimConfig → byte-
        # identical output. Pre-T0-8 these would diverge (auto-seeded
        # by time + each simulator's address).
        a = sim_a.simulate(n=10)
        b = sim_b.simulate(n=10)
        assert [r["sequence"] for r in a] == [r["sequence"] for r in b]

    def test_reset_seed_clears_fallback(self):
        """reset_seed() must clear the global so subsequent runs
        without an explicit seed are NOT silently locked to the
        last set_seed value."""
        set_seed(42)
        r1 = Experiment.on("human_igh").run(n=5)

        reset_seed()
        # Without a global, the C engine auto-seeds. Two consecutive
        # auto-seeded runs may collide on the same wall-clock second
        # AND the same simulator address; we just assert the global
        # is cleared (the auto-seed property itself is covered by
        # TestPerSimulatorRng in test_simulation_integrity.py).
        assert get_seed() is None
        # Smoke: a follow-up run still produces output.
        r2 = Experiment.on("human_igh").run(n=5)
        assert len(r2) == 5
        # And r1 was generated under set_seed(42) so it's still
        # deterministic relative to itself.
        set_seed(42)
        r3 = Experiment.on("human_igh").run(n=5)
        assert [r["sequence"] for r in r1] == [r["sequence"] for r in r3]


# ============================================================
# T1-8: streaming RNG semantics + CompiledSimulator.set_seed
# ============================================================

class TestCompiledSimulatorSetSeed:
    """T1-8: CompiledSimulator owns a single PCG32 stream that
    persists across simulate() calls. set_seed() resets the stream."""

    def test_set_seed_method_exists(self):
        """The high-level CompiledSimulator wrapper must expose
        set_seed (not just CSimulator). Pre-T1-8 users had to reach
        into sim._sim.set_seed which is private."""
        sim = Experiment.on("human_igh").compile(seed=42)
        assert hasattr(sim, "set_seed"), \
            "CompiledSimulator.set_seed missing — users would have to reach into _sim"
        assert callable(sim.set_seed)

    def test_set_seed_resets_stream_deterministically(self):
        """sim.set_seed(42); sim.simulate(N) is byte-identical across
        calls with the same N. This is the explicit reset path."""
        sim = Experiment.on("human_igh").compile(seed=42)

        sim.set_seed(42)
        a = list(sim.simulate(10))
        sim.set_seed(42)
        b = list(sim.simulate(10))

        assert [r["sequence"] for r in a] == [r["sequence"] for r in b], \
            "set_seed(42) must reset the stream to a deterministic state"

    def test_streaming_invariant_two_batches_equal_one_batch(self):
        """Streaming invariant: simulate(N + M) from a fresh seed is
        the same as simulate(N) followed by simulate(M) from the same
        fresh seed. This is the core promise of the streaming RNG."""
        sim1 = Experiment.on("human_igh").compile(seed=999)
        twenty = list(sim1.simulate(20))

        sim2 = Experiment.on("human_igh").compile(seed=999)
        first_ten  = list(sim2.simulate(10))
        second_ten = list(sim2.simulate(10))

        assert [r["sequence"] for r in twenty] == \
               [r["sequence"] for r in first_ten + second_ten], \
            "simulate(20) must equal simulate(10) + simulate(10) on same seed"

    def test_simulate_does_not_silently_reseed(self):
        """Calling simulate() multiple times advances the stream — it
        does NOT silently re-seed. Two back-to-back simulate(10) calls
        on the same simulator must produce DIFFERENT records."""
        sim = Experiment.on("human_igh").compile(seed=42)
        a = list(sim.simulate(10))
        b = list(sim.simulate(10))

        # Sequences are 250+ bp; collisions on random V(D)J output
        # are vanishingly unlikely — if equal, the RNG re-seeded.
        a_seqs = [r["sequence"] for r in a]
        b_seqs = [r["sequence"] for r in b]
        assert a_seqs != b_seqs, \
            "back-to-back simulate(10) returned identical records — " \
            "stream is silently re-seeding (regression)"

    def test_set_seed_zero_is_auto_seed(self):
        """seed=0 is the engine's magic 'auto-seed from time + sim
        address' value. Two consecutive set_seed(0) calls on the same
        simulator are NOT guaranteed deterministic relative to each
        other (auto-seed property), but they produce VALID output."""
        sim = Experiment.on("human_igh").compile(seed=42)
        sim.set_seed(0)
        a = list(sim.simulate(5))
        # Must not crash, must produce 5 records
        assert len(a) == 5
        for r in a:
            assert "sequence" in r and len(r["sequence"]) > 0

    def test_set_seed_after_simulate_resets_mid_life(self):
        """set_seed mid-life must reset the stream — even after
        simulate() has already advanced it."""
        sim = Experiment.on("human_igh").compile(seed=42)
        sim.simulate(10)             # advance stream past initial seed
        sim.set_seed(42)             # reset
        a = list(sim.simulate(10))   # records 0..9 from seed 42

        # Compare against a fresh simulator with the same seed
        fresh = Experiment.on("human_igh").compile(seed=42)
        b = list(fresh.simulate(10))

        assert [r["sequence"] for r in a] == [r["sequence"] for r in b], \
            "set_seed mid-life should produce the same output as a " \
            "freshly compiled simulator with the same seed"
