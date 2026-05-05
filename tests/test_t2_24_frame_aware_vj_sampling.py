"""T2-24: frame-aware NP1 length sampling for VJ chains.

V5 step 29 makes `step_assemble` constrain NP1 length on VJ chains
(no D segment) so that the assembled junction is in-frame **by
construction**. The retry loop becomes a near no-op on light chains
under PRODUCTIVE_ONLY: only stop codons within an in-frame junction
can still cause a retry, and those are rare.

These tests pin the behavior:

  1. With max_productive_attempts=1 (no retry budget), human_igk and
     human_igl PRODUCTIVE_ONLY both produce >=90% productive
     sequences on the very first rearrangement attempt. Pre-step-29
     the baseline is ~25% (matches the empirical productive rate
     under random sampling).
  2. PRODUCTIVE_MIXED (no contract, no constraint) is unchanged: the
     productive rate stays at the natural baseline of ~25%.
  3. PRODUCTIVE_ONLY with the default retry budget still produces
     ~100% productive sequences (no regression on the existing
     guarantee).
  4. Heavy chain (VDJ) is intentionally unchanged in step 29 — VDJ
     gets joint NP1+NP2 frame-aware sampling in a follow-up step.
     This test pins the VJ-only scope so a later step that improves
     VDJ doesn't accidentally regress this guarantee.
"""
from __future__ import annotations

from GenAIRR import Experiment, Productivity


def _productive_fraction(species_chain: str, productivity, n: int = 500,
                         max_attempts: int = 1, seed: int = 42) -> float:
    sim = Experiment.on(species_chain).compile(
        seed=seed, productivity=productivity)
    sim._sim.set_param("max_productive_attempts", max_attempts)
    records = sim.simulate(n=n)
    productive = sum(1 for r in records if r["productive"])
    return productive / n


class TestFrameAwareVJ:
    """VJ chains (light chain) under PRODUCTIVE_ONLY assemble in-frame
    by construction — no rearrangement retries needed for frame."""

    def test_kappa_productive_on_first_attempt(self):
        """human_igk: with no retry budget, >=90% must still be productive."""
        frac = _productive_fraction("human_igk", Productivity.PRODUCTIVE_ONLY)
        assert frac >= 0.90, (
            f"human_igk single-attempt productive rate {frac:.2%} < 90%; "
            f"frame-aware NP1 sampling broken or empirical NP distribution "
            f"has no in-frame mass for this allele pool."
        )

    def test_lambda_productive_on_first_attempt(self):
        """human_igl: with no retry budget, >=90% must still be productive."""
        frac = _productive_fraction("human_igl", Productivity.PRODUCTIVE_ONLY)
        assert frac >= 0.90, (
            f"human_igl single-attempt productive rate {frac:.2%} < 90%; "
            f"frame-aware NP1 sampling broken."
        )

    def test_kappa_mixed_baseline_unchanged(self):
        """PRODUCTIVE_MIXED must NOT use the constraint — empirical rate
        stays low, demonstrating the constraint is gated on contract."""
        frac = _productive_fraction("human_igk", Productivity.PRODUCTIVE_MIXED)
        # Empirical productive baseline for human_igk is ~25%. The test
        # gives a wide margin so RNG drift doesn't flake it.
        assert frac < 0.50, (
            f"human_igk MIXED productive rate {frac:.2%} >= 50% — "
            f"frame constraint may be leaking into non-productive contracts."
        )

    def test_kappa_strict_invariant_with_default_budget(self):
        """With the default retry budget (25), PRODUCTIVE_ONLY must hit
        ~100% productive — the residual ~5% stop-codon failures get
        absorbed by a tiny number of retries."""
        frac = _productive_fraction(
            "human_igk", Productivity.PRODUCTIVE_ONLY,
            n=500, max_attempts=25)
        assert frac >= 0.99, (
            f"human_igk PRODUCTIVE_ONLY (budget=25) leaked: "
            f"{frac:.2%} < 99%."
        )

    def test_heavy_chain_vdj_unchanged_in_step29(self):
        """VDJ joint NP1+NP2 frame-aware sampling is the next step.
        For now heavy chain on a single attempt should still match the
        legacy ~25% empirical productive baseline."""
        frac = _productive_fraction("human_igh", Productivity.PRODUCTIVE_ONLY)
        # Loose bounds: anything between 10% and 40% confirms VDJ is on
        # the legacy unconstrained path. >50% would mean a future step
        # changed VDJ behaviour and this test should be updated then.
        assert 0.10 <= frac <= 0.50, (
            f"human_igh single-attempt productive rate {frac:.2%}: VDJ "
            f"frame-aware sampling may have been added — update this "
            f"test (it pins step 29 scope = VJ only)."
        )
