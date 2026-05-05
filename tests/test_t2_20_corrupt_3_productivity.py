"""T2-20: 3' corruption must honor the productive contract."""
from __future__ import annotations

from GenAIRR import Experiment, Productivity
from GenAIRR.ops import with_3prime_loss


class TestProductiveGuardedCorrupt3:
    """Observed 3' clipping should be skipped under PRODUCTIVE_ONLY
    when it would break final productivity."""

    def test_productive_only_skips_unsafe_3prime_clip_seed1(self):
        cfg = dict(min_remove=120, max_remove=180, min_add=0, max_add=0)

        mixed = (
            Experiment.on("human_igh")
            .sequence(with_3prime_loss(**cfg))
            .compile(seed=1, productivity=Productivity.PRODUCTIVE_MIXED)
        )
        mixed._sim.set_param("max_productive_attempts", 200)
        mixed_rec, _, mixed_trace = mixed._sim.simulate_one_hooked(["final"])

        prod = (
            Experiment.on("human_igh")
            .sequence(with_3prime_loss(**cfg))
            .compile(seed=1, productivity=Productivity.PRODUCTIVE_ONLY)
        )
        prod._sim.set_param("max_productive_attempts", 200)
        prod_rec, _, prod_trace = prod._sim.simulate_one_hooked(["final"])

        assert mixed_rec["productive"] is False, (
            "Seed 1 should still show that mixed-mode 3' clipping can "
            "produce a non-productive observed read"
        )
        assert prod_rec["productive"] is True
        assert len(prod_rec["sequence"]) > len(mixed_rec["sequence"]), (
            "Guarded 3' corruption run should skip the unsafe heavy clip"
        )
        assert "[corrupt_3'] removed" in mixed_trace
        assert "[corrupt_3'] skipped productive-unsafe event" not in mixed_trace
        assert "[corrupt_3'] skipped productive-unsafe event" in prod_trace
        assert "[corrupt_3'] removed" not in prod_trace
