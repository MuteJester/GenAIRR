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
