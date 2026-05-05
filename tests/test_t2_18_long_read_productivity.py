"""T2-18: long-read homopolymer indels must honor the productive
contract."""
from __future__ import annotations

from GenAIRR import Experiment, Productivity
from GenAIRR.ops import long_read


class TestProductiveGuardedLongRead:
    """Observed long-read homopolymer indels should be filtered under
    PRODUCTIVE_ONLY when they would break productivity."""

    def test_productive_only_rejects_unsafe_long_read_indels_seed1(self):
        cfg = dict(error_rate=0.10, min_run_length=3, insertion_bias=0.7)

        mixed = (Experiment.on("human_igh")
                 .sequence(long_read(**cfg))
                 .compile(seed=1, productivity=Productivity.PRODUCTIVE_MIXED))
        mixed._sim.set_param("max_productive_attempts", 200)
        mixed_rec, _, mixed_trace = mixed._sim.simulate_one_hooked(["final"])

        prod = (Experiment.on("human_igh")
                .sequence(long_read(**cfg))
                .compile(seed=1, productivity=Productivity.PRODUCTIVE_ONLY))
        prod._sim.set_param("max_productive_attempts", 200)
        prod_rec, _, prod_trace = prod._sim.simulate_one_hooked(["final"])

        assert mixed_rec["productive"] is False, (
            "Seed 1 should still show that unguarded mixed-mode long-read "
            "indels can produce a non-productive observed read"
        )
        assert prod_rec["productive"] is True
        assert prod_rec["n_insertions"] < mixed_rec["n_insertions"], (
            "Guarded long-read run should skip at least one unsafe insertion"
        )
        assert "[guard] reject insertion stage=observed" not in mixed_trace
        assert "[guard] reject insertion stage=observed" in prod_trace
