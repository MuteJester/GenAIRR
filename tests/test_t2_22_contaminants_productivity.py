"""T2-22: contaminants must honor the productive contract."""
from __future__ import annotations

from GenAIRR import Experiment, Productivity
from GenAIRR.ops import with_contaminants


class TestProductiveGuardedContaminants:
    """Observed contaminant replacement should be skipped under
    PRODUCTIVE_ONLY because it always destroys receptor productivity."""

    def test_productive_only_skips_contaminant_replacement_seed1(self):
        mixed = (
            Experiment.on("human_igh")
            .observe(with_contaminants(rate=1.0, source="random"))
            .compile(seed=1, productivity=Productivity.PRODUCTIVE_MIXED)
        )
        mixed._sim.set_param("max_productive_attempts", 200)
        mixed_rec, _, mixed_trace = mixed._sim.simulate_one_hooked(["final"])

        prod = (
            Experiment.on("human_igh")
            .observe(with_contaminants(rate=1.0, source="random"))
            .compile(seed=1, productivity=Productivity.PRODUCTIVE_ONLY)
        )
        prod._sim.set_param("max_productive_attempts", 200)
        prod_rec, _, prod_trace = prod._sim.simulate_one_hooked(["final"])

        assert mixed_rec["productive"] is False
        assert mixed_rec["is_contaminant"] is True
        assert "contaminant" in mixed_rec["note"]
        assert prod_rec["productive"] is True
        assert prod_rec["is_contaminant"] is False
        assert "[contaminant] type=random" in mixed_trace
        assert "[contaminant] skipped productive-unsafe contamination" not in mixed_trace
        assert "[contaminant] skipped productive-unsafe contamination" in prod_trace
