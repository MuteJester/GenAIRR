"""T2-21: observed indels must honor the productive contract."""
from __future__ import annotations

from GenAIRR import Experiment, Productivity
from GenAIRR.ops import with_indels


class TestProductiveGuardedObservedIndels:
    """Observed indels should be skipped under PRODUCTIVE_ONLY when
    they would break final productivity."""

    def test_productive_only_rejects_unsafe_observed_indels_seed1(self):
        prob = 0.02

        mixed = (
            Experiment.on("human_igh")
            .observe(with_indels(prob=prob))
            .compile(seed=1, productivity=Productivity.PRODUCTIVE_MIXED)
        )
        mixed._sim.set_param("max_productive_attempts", 200)
        mixed_rec, _, mixed_trace = mixed._sim.simulate_one_hooked(["final"])

        prod = (
            Experiment.on("human_igh")
            .observe(with_indels(prob=prob))
            .compile(seed=1, productivity=Productivity.PRODUCTIVE_ONLY)
        )
        prod._sim.set_param("max_productive_attempts", 200)
        prod_rec, _, prod_trace = prod._sim.simulate_one_hooked(["final"])

        assert mixed_rec["productive"] is False, (
            "Seed 1 should still show that mixed-mode observed indels can "
            "produce a non-productive read"
        )
        assert prod_rec["productive"] is True
        assert (
            prod_rec["n_insertions"] + prod_rec["n_deletions"]
            < mixed_rec["n_insertions"] + mixed_rec["n_deletions"]
        ), "Guarded observed-indel run should skip at least one unsafe indel"
        assert "[guard] reject insertion stage=observed" not in mixed_trace
        assert "[guard] reject deletion stage=observed" not in mixed_trace
        assert (
            "[guard] reject insertion stage=observed" in prod_trace
            or "[guard] reject deletion stage=observed" in prod_trace
        )
