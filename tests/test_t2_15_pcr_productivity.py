"""T2-15: PCR observation errors must honor the productive contract."""
from __future__ import annotations

from GenAIRR import Experiment, Productivity
from GenAIRR.ops import with_pcr


class TestProductiveGuardedPCR:
    """Observed-stage PCR errors should skip unsafe substitutions under
    PRODUCTIVE_ONLY, while mixed mode still permits the same seed to
    fall non-productive."""

    def test_productive_only_rejects_unsafe_pcr_substitutions_seed2(self):
        kwargs = dict(error_rate=1e-3, cycles=30)

        mixed = (Experiment.on("human_igh")
                 .prepare(with_pcr(**kwargs))
                 .compile(seed=2, productivity=Productivity.PRODUCTIVE_MIXED))
        mixed._sim.set_param("max_productive_attempts", 200)
        mixed_rec, _, mixed_trace = mixed._sim.simulate_one_hooked(
            ["post_pcr", "final"]
        )

        prod = (Experiment.on("human_igh")
                .prepare(with_pcr(**kwargs))
                .compile(seed=2, productivity=Productivity.PRODUCTIVE_ONLY))
        prod._sim.set_param("max_productive_attempts", 200)
        prod_rec, _, prod_trace = prod._sim.simulate_one_hooked(
            ["post_pcr", "final"]
        )

        assert mixed_rec["productive"] is False, (
            "Seed 2 should still show that unguarded mixed-mode PCR can "
            "produce a non-productive observed read"
        )
        assert prod_rec["productive"] is True
        assert prod_rec["n_pcr_errors"] > 0
        assert "[guard] reject substitution stage=observed" not in mixed_trace
        assert "[guard] reject substitution stage=observed" in prod_trace
