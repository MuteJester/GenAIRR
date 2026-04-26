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
