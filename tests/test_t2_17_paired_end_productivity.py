"""T2-17: paired-end observation modeling must emit real substitutions
and honor the productive contract."""
from __future__ import annotations

from GenAIRR import Experiment, Productivity
from GenAIRR.ops import paired_end


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


class TestPairedEndEmitsRealMutations:
    """Observed paired-end substitutions should not record
    self-substitutions on lowercase real-simulation reads."""

    def test_paired_end_does_not_emit_self_substitutions(self):
        rec = (Experiment.on("human_igh")
               .sequence(paired_end(120))
               .run(n=1, seed=42))[0]

        assert rec["n_sequencing_errors"] > 0
        changes = _parse_changes(rec["sequencing_errors"])
        assert changes
        assert all(src != dst for src, dst in changes), (
            "paired_end recorded self-substitutions on a lowercase "
            "real-simulation read"
        )


class TestProductiveGuardedPairedEnd:
    """Paired-end gap masking and substitutions should be filtered under
    PRODUCTIVE_ONLY when they would break productivity."""

    def test_productive_only_rejects_unsafe_paired_end_edits_seed1(self):
        read_len = 75

        mixed = (Experiment.on("human_igh")
                 .sequence(paired_end(read_len))
                 .compile(seed=1, productivity=Productivity.PRODUCTIVE_MIXED))
        mixed._sim.set_param("max_productive_attempts", 200)
        mixed_rec, _, mixed_trace = mixed._sim.simulate_one_hooked(["final"])

        prod = (Experiment.on("human_igh")
                .sequence(paired_end(read_len))
                .compile(seed=1, productivity=Productivity.PRODUCTIVE_ONLY))
        prod._sim.set_param("max_productive_attempts", 200)
        prod_rec, _, prod_trace = prod._sim.simulate_one_hooked(["final"])

        assert mixed_rec["productive"] is False, (
            "Seed 1 should still show that unguarded mixed-mode paired-end "
            "modeling can produce a non-productive observed read"
        )
        assert prod_rec["productive"] is True
        assert "[guard] reject substitution stage=observed" not in mixed_trace
        assert "[guard] reject substitution stage=observed" in prod_trace
        assert prod_rec["sequence"].count("N") < mixed_rec["sequence"].count("N"), (
            "Guarded paired-end run should skip at least one unsafe gap "
            "masking edit"
        )
