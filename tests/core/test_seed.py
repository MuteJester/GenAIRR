"""
Tests for GenAIRR seed system reproducibility.

Verifies that set_seed() produces identical results across runs.
"""

import unittest

from GenAIRR import (
    set_seed,
    get_seed,
    reset_seed,
    simulate,
    Pipeline,
    steps,
    HUMAN_IGH_OGRDB,
    HUMAN_IGK_OGRDB,
    Uniform,
    S5F,
)


class TestSeedSystem(unittest.TestCase):
    """Test seed management functions."""

    def test_set_seed_stores_value(self):
        """set_seed should store the seed value."""
        set_seed(42)
        self.assertEqual(get_seed(), 42)

        set_seed(12345)
        self.assertEqual(get_seed(), 12345)

    def test_get_seed_returns_none_initially(self):
        """get_seed should return None before set_seed is called."""
        reset_seed()
        self.assertIsNone(get_seed())

    def test_reset_seed_clears_value(self):
        """reset_seed should clear the stored seed."""
        set_seed(42)
        self.assertEqual(get_seed(), 42)
        reset_seed()
        self.assertIsNone(get_seed())


class TestReproducibility(unittest.TestCase):
    """Test that seeding produces reproducible results."""

    def test_simulate_reproducibility_heavy_chain(self):
        """Same seed should produce identical heavy chain sequences."""
        set_seed(42)
        result1 = simulate(HUMAN_IGH_OGRDB, Uniform(0.01, 0.05))

        set_seed(42)
        result2 = simulate(HUMAN_IGH_OGRDB, Uniform(0.01, 0.05))

        self.assertEqual(result1.sequence, result2.sequence)
        self.assertEqual(result1.v_call, result2.v_call)
        self.assertEqual(result1.d_call, result2.d_call)
        self.assertEqual(result1.j_call, result2.j_call)

    def test_simulate_reproducibility_light_chain(self):
        """Same seed should produce identical light chain sequences."""
        set_seed(123)
        result1 = simulate(HUMAN_IGK_OGRDB, Uniform(0.01, 0.05))

        set_seed(123)
        result2 = simulate(HUMAN_IGK_OGRDB, Uniform(0.01, 0.05))

        self.assertEqual(result1.sequence, result2.sequence)
        self.assertEqual(result1.v_call, result2.v_call)
        self.assertEqual(result1.j_call, result2.j_call)

    def test_different_seeds_produce_different_results(self):
        """Different seeds should generally produce different sequences."""
        set_seed(42)
        result1 = simulate(HUMAN_IGH_OGRDB, Uniform(0.01, 0.05))

        set_seed(999)
        result2 = simulate(HUMAN_IGH_OGRDB, Uniform(0.01, 0.05))

        # While theoretically could match, extremely unlikely
        self.assertNotEqual(result1.sequence, result2.sequence)

    def test_multiple_sequences_reproducibility(self):
        """Multiple sequences in a batch should be reproducible."""
        set_seed(42)
        results1 = simulate(HUMAN_IGH_OGRDB, Uniform(0.01, 0.05), n=5)

        set_seed(42)
        results2 = simulate(HUMAN_IGH_OGRDB, Uniform(0.01, 0.05), n=5)

        for r1, r2 in zip(results1, results2):
            self.assertEqual(r1.sequence, r2.sequence)
            self.assertEqual(r1.v_call, r2.v_call)

    def test_pipeline_reproducibility(self):
        """Full pipeline should be reproducible with same seed."""
        set_seed(42)
        pipeline1 = Pipeline(
            config=HUMAN_IGH_OGRDB,
            steps=[
                steps.SimulateSequence(Uniform(0.01, 0.05), productive=True),
                steps.FixVPositionAfterTrimmingIndexAmbiguity(),
            ]
        )
        result1 = pipeline1.execute()

        set_seed(42)
        pipeline2 = Pipeline(
            config=HUMAN_IGH_OGRDB,
            steps=[
                steps.SimulateSequence(Uniform(0.01, 0.05), productive=True),
                steps.FixVPositionAfterTrimmingIndexAmbiguity(),
            ]
        )
        result2 = pipeline2.execute()

        self.assertEqual(result1.sequence, result2.sequence)

    def test_pipeline_with_corruption_reproducibility(self):
        """Pipeline with corruption step should be reproducible."""
        set_seed(42)
        pipeline1 = Pipeline(
            config=HUMAN_IGH_OGRDB,
            steps=[
                steps.SimulateSequence(Uniform(0.0, 0.0), productive=True),
                steps.CorruptSequenceBeginning(probability=0.8),
            ]
        )
        result1 = pipeline1.execute()

        set_seed(42)
        pipeline2 = Pipeline(
            config=HUMAN_IGH_OGRDB,
            steps=[
                steps.SimulateSequence(Uniform(0.0, 0.0), productive=True),
                steps.CorruptSequenceBeginning(probability=0.8),
            ]
        )
        result2 = pipeline2.execute()

        self.assertEqual(result1.sequence, result2.sequence)

    def test_pipeline_with_indels_reproducibility(self):
        """Pipeline with indels step should be reproducible."""
        set_seed(42)
        pipeline1 = Pipeline(
            config=HUMAN_IGH_OGRDB,
            steps=[
                steps.SimulateSequence(Uniform(0.0, 0.0), productive=True),
                steps.InsertIndels(
                    probability=0.5,
                    max_indels=3,
                    insertion_probability=0.5,
                    deletion_probability=0.5,
                ),
            ]
        )
        result1 = pipeline1.execute()

        set_seed(42)
        pipeline2 = Pipeline(
            config=HUMAN_IGH_OGRDB,
            steps=[
                steps.SimulateSequence(Uniform(0.0, 0.0), productive=True),
                steps.InsertIndels(
                    probability=0.5,
                    max_indels=3,
                    insertion_probability=0.5,
                    deletion_probability=0.5,
                ),
            ]
        )
        result2 = pipeline2.execute()

        self.assertEqual(result1.sequence, result2.sequence)


class TestSeedWithDifferentOperations(unittest.TestCase):
    """Test seed affects various random operations consistently."""

    def test_seed_affects_allele_selection(self):
        """Seed should control random allele selection."""
        set_seed(42)
        result1 = simulate(HUMAN_IGH_OGRDB, Uniform(0.0, 0.0))

        set_seed(42)
        result2 = simulate(HUMAN_IGH_OGRDB, Uniform(0.0, 0.0))

        # With 0 mutation rate, sequences differ only by allele choices
        self.assertEqual(result1.v_call, result2.v_call)
        self.assertEqual(result1.d_call, result2.d_call)
        self.assertEqual(result1.j_call, result2.j_call)

    def test_seed_affects_mutation_positions(self):
        """Seed should control which positions get mutated."""
        set_seed(42)
        result1 = simulate(HUMAN_IGH_OGRDB, Uniform(0.05, 0.05))  # Fixed 5% rate

        set_seed(42)
        result2 = simulate(HUMAN_IGH_OGRDB, Uniform(0.05, 0.05))

        # Mutations should be at same positions
        self.assertEqual(result1.mutations, result2.mutations)


class TestS5FReproducibility(unittest.TestCase):
    """Test S5F mutation model reproducibility with seed system."""

    def test_s5f_heavy_chain_reproducibility(self):
        """S5F model should produce identical results with same seed."""
        set_seed(42)
        result1 = simulate(HUMAN_IGH_OGRDB, S5F(0.01, 0.05))

        set_seed(42)
        result2 = simulate(HUMAN_IGH_OGRDB, S5F(0.01, 0.05))

        self.assertEqual(result1.sequence, result2.sequence)
        self.assertEqual(result1.mutations, result2.mutations)

    def test_s5f_light_chain_reproducibility(self):
        """S5F model should produce identical results for light chains."""
        set_seed(123)
        result1 = simulate(HUMAN_IGK_OGRDB, S5F(0.01, 0.05))

        set_seed(123)
        result2 = simulate(HUMAN_IGK_OGRDB, S5F(0.01, 0.05))

        self.assertEqual(result1.sequence, result2.sequence)
        self.assertEqual(result1.mutations, result2.mutations)

    def test_s5f_different_seeds_different_results(self):
        """S5F with different seeds should produce different results."""
        set_seed(42)
        result1 = simulate(HUMAN_IGH_OGRDB, S5F(0.05, 0.10))

        set_seed(999)
        result2 = simulate(HUMAN_IGH_OGRDB, S5F(0.05, 0.10))

        # With different seeds, mutations should differ
        self.assertNotEqual(result1.mutations, result2.mutations)


if __name__ == "__main__":
    unittest.main()
