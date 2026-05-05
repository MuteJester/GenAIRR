"""T2-16: quality-profile sequencing errors must emit real
substitutions and honor the productive contract."""
from __future__ import annotations

from GenAIRR import Experiment, Productivity
from GenAIRR.ops import with_quality_profile


def _parse_changes(change_str: str):
    out = []
    for entry in change_str.split(","):
        entry = entry.strip()
        if ":" not in entry or ">" not in entry:
            continue
        _, change = entry.split(":", 1)
        src, dst = change.split(">", 1)
        out.append((src, dst))
    return out


class TestQualityErrorsEmitRealMutations:
    """Observed quality errors should not record self-substitutions on
    lowercase real-simulation reads."""

    def test_quality_errors_do_not_emit_self_substitutions(self):
        rec = (Experiment.on("human_igh")
               .sequence(with_quality_profile(base=0.5, peak=1.0))
               .run(n=1, seed=42))[0]

        assert rec["n_sequencing_errors"] > 0
        changes = _parse_changes(rec["sequencing_errors"])
        assert changes
        assert all(src != dst for src, dst in changes), (
            "quality_errors recorded self-substitutions on a lowercase "
            "real-simulation read"
        )


class TestProductiveGuardedQualityErrors:
    """Observed quality-profile substitutions should be filtered under
    PRODUCTIVE_ONLY when they would break productivity."""

    def test_productive_only_rejects_unsafe_quality_substitutions_seed1(self):
        kwargs = dict(base=0.3, peak=0.8)

        mixed = (Experiment.on("human_igh")
                 .sequence(with_quality_profile(**kwargs))
                 .compile(seed=1, productivity=Productivity.PRODUCTIVE_MIXED))
        mixed._sim.set_param("max_productive_attempts", 200)
        mixed_rec, _, mixed_trace = mixed._sim.simulate_one_hooked(
            ["post_quality", "final"]
        )

        prod = (Experiment.on("human_igh")
                .sequence(with_quality_profile(**kwargs))
                .compile(seed=1, productivity=Productivity.PRODUCTIVE_ONLY))
        prod._sim.set_param("max_productive_attempts", 200)
        prod_rec, _, prod_trace = prod._sim.simulate_one_hooked(
            ["post_quality", "final"]
        )

        assert mixed_rec["productive"] is False, (
            "Seed 1 should still show that unguarded mixed-mode quality "
            "errors can produce a non-productive observed read"
        )
        assert prod_rec["productive"] is True
        assert prod_rec["n_sequencing_errors"] > 0
        assert "[guard] reject substitution stage=observed" not in mixed_trace
        assert "[guard] reject substitution stage=observed" in prod_trace
