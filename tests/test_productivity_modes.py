"""T1-6: Productivity enum filter modes.

Verifies that under each Productivity flag the simulator returns
only the appropriate kind of sequence — strict invariants, no leaks
allowed under default conditions:

  - Productivity.PRODUCTIVE_ONLY      → every record is productive.
  - Productivity.NON_PRODUCTIVE_ONLY  → every record is non-productive.
  - Productivity.PRODUCTIVE_MIXED     → mix of both classes present.

Backward-compatibility tests for the legacy `productive=True/False`
bool API are also included.

Note on retry exhaustion: the C engine retries up to
``max_productive_attempts`` (default 25). With human_igh's pure
recombination productive baseline ≈ 22.8%, the per-record exhaustion
probability for PRODUCTIVE_ONLY is 0.772^25 ≈ 0.17%. We bump the
retry budget for the strict-invariant tests so leakage is genuinely
impossible (0.772^200 ≈ 5e-22). The default budget is left in place
for the back-compat path tests that already accept ≤1 leak.
"""
from __future__ import annotations

import pytest
from GenAIRR import Experiment, Productivity
from GenAIRR.ops import (
    long_read,
    model,
    paired_end,
    rate,
    with_contaminants,
    with_indels,
    with_3prime_loss,
    with_antigen_selection,
    with_5prime_loss,
    with_d_inversion,
    with_pcr,
    with_quality_profile,
    with_receptor_revision,
)


def _boost_retry_budget(sim, attempts: int = 200) -> None:
    """Raise max_productive_attempts so retry exhaustion can't leak."""
    sim._sim.set_param("max_productive_attempts", attempts)


# ── PRODUCTIVE_ONLY ───────────────────────────────────────────────


class TestProductiveOnly:
    """Productivity.PRODUCTIVE_ONLY: every record must be productive."""

    def test_strict_no_non_productive_leaks(self):
        sim = Experiment.on("human_igh").compile(
            seed=42, productivity=Productivity.PRODUCTIVE_ONLY)
        _boost_retry_budget(sim, 200)
        records = sim.simulate(n=200)
        non_prod = [r for r in records if not r["productive"]]
        assert not non_prod, (
            f"PRODUCTIVE_ONLY leaked {len(non_prod)}/200 non-productive "
            f"records (max_productive_attempts=200, exhaustion ≈ 0)."
        )

    def test_works_via_run_method(self):
        result = (Experiment.on("human_igh")
                  .run(n=20, seed=42, productivity=Productivity.PRODUCTIVE_ONLY))
        # Default retries (25) — accept up to 1 retry-exhaustion leak.
        non_prod = sum(1 for r in result if not r["productive"])
        assert non_prod <= 1, f"too many leaks: {non_prod}/20"

    def test_strict_no_non_productive_leaks_with_s5f_mutation(self):
        sim = (
            Experiment.on("human_igh")
            .mutate(rate(0.02, 0.05), model("s5f"))
            .compile(seed=42, productivity=Productivity.PRODUCTIVE_ONLY)
        )
        _boost_retry_budget(sim, 200)
        records = sim.simulate(n=100)
        non_prod = [r for r in records if not r["productive"]]
        assert not non_prod, (
            f"PRODUCTIVE_ONLY + S5F leaked {len(non_prod)}/100 "
            f"non-productive records after mutation."
        )

    def test_strict_no_non_productive_leaks_with_uniform_mutation(self):
        sim = (
            Experiment.on("human_igh")
            .mutate(rate(0.02, 0.05), model("uniform"))
            .compile(seed=42, productivity=Productivity.PRODUCTIVE_ONLY)
        )
        _boost_retry_budget(sim, 200)
        records = sim.simulate(n=100)
        non_prod = [r for r in records if not r["productive"]]
        assert not non_prod, (
            f"PRODUCTIVE_ONLY + uniform leaked {len(non_prod)}/100 "
            f"non-productive records after mutation."
        )

    def test_strict_no_non_productive_leaks_with_selection(self):
        sim = (
            Experiment.on("human_igh")
            .mutate(
                rate(0.05, 0.10),
                model("s5f"),
                with_antigen_selection(1.0),
            )
            .compile(seed=42, productivity=Productivity.PRODUCTIVE_ONLY)
        )
        _boost_retry_budget(sim, 200)
        records = sim.simulate(n=200)
        non_prod = [r for r in records if not r["productive"]]
        assert not non_prod, (
            f"PRODUCTIVE_ONLY + selection leaked {len(non_prod)}/200 "
            f"non-productive records after post-mutation selection."
        )

    def test_strict_no_non_productive_leaks_with_d_inversion(self):
        sim = (
            Experiment.on("human_igh")
            .recombine(with_d_inversion(1.0))
            .compile(seed=42, productivity=Productivity.PRODUCTIVE_ONLY)
        )
        _boost_retry_budget(sim, 200)
        records = sim.simulate(n=200)
        non_prod = [r for r in records if not r["productive"]]
        assert not non_prod, (
            f"PRODUCTIVE_ONLY + d_inversion leaked {len(non_prod)}/200 "
            f"non-productive records after pre-mutation inversion."
        )

    def test_strict_no_non_productive_leaks_with_receptor_revision(self):
        sim = (
            Experiment.on("human_igh")
            .recombine(with_receptor_revision(prob=1.0, footprint=(1, 1)))
            .compile(seed=42, productivity=Productivity.PRODUCTIVE_ONLY)
        )
        _boost_retry_budget(sim, 200)
        records = sim.simulate(n=200)
        non_prod = [r for r in records if not r["productive"]]
        assert not non_prod, (
            f"PRODUCTIVE_ONLY + receptor_revision leaked {len(non_prod)}/200 "
            f"non-productive records after pre-mutation revision."
        )

    def test_strict_no_non_productive_leaks_with_pcr(self):
        sim = (
            Experiment.on("human_igh")
            .prepare(with_pcr(error_rate=1e-3, cycles=30))
            .compile(seed=42, productivity=Productivity.PRODUCTIVE_ONLY)
        )
        _boost_retry_budget(sim, 200)
        records = sim.simulate(n=200)
        non_prod = [r for r in records if not r["productive"]]
        assert not non_prod, (
            f"PRODUCTIVE_ONLY + pcr leaked {len(non_prod)}/200 "
            f"non-productive records after observed-stage PCR errors."
        )

    def test_strict_no_non_productive_leaks_with_quality_profile(self):
        sim = (
            Experiment.on("human_igh")
            .sequence(with_quality_profile(base=0.3, peak=0.8))
            .compile(seed=42, productivity=Productivity.PRODUCTIVE_ONLY)
        )
        _boost_retry_budget(sim, 200)
        records = sim.simulate(n=200)
        non_prod = [r for r in records if not r["productive"]]
        assert not non_prod, (
            f"PRODUCTIVE_ONLY + quality_profile leaked {len(non_prod)}/200 "
            f"non-productive records after observed-stage sequencing errors."
        )

    def test_strict_no_non_productive_leaks_with_paired_end(self):
        sim = (
            Experiment.on("human_igh")
            .sequence(paired_end(100))
            .compile(seed=42, productivity=Productivity.PRODUCTIVE_ONLY)
        )
        _boost_retry_budget(sim, 200)
        records = sim.simulate(n=200)
        non_prod = [r for r in records if not r["productive"]]
        assert not non_prod, (
            f"PRODUCTIVE_ONLY + paired_end leaked {len(non_prod)}/200 "
            f"non-productive records after observed-stage gap/error modeling."
        )

    def test_strict_no_non_productive_leaks_with_long_read(self):
        sim = (
            Experiment.on("human_igh")
            .sequence(long_read())
            .compile(seed=42, productivity=Productivity.PRODUCTIVE_ONLY)
        )
        _boost_retry_budget(sim, 200)
        records = sim.simulate(n=200)
        non_prod = [r for r in records if not r["productive"]]
        assert not non_prod, (
            f"PRODUCTIVE_ONLY + long_read leaked {len(non_prod)}/200 "
            f"non-productive records after observed-stage homopolymer indels."
        )

    def test_strict_no_non_productive_leaks_with_5prime_loss(self):
        sim = (
            Experiment.on("human_igh")
            .sequence(with_5prime_loss(
                min_remove=250, max_remove=320, min_add=0, max_add=0))
            .compile(seed=42, productivity=Productivity.PRODUCTIVE_ONLY)
        )
        _boost_retry_budget(sim, 200)
        records = sim.simulate(n=200)
        non_prod = [r for r in records if not r["productive"]]
        assert not non_prod, (
            f"PRODUCTIVE_ONLY + 5prime_loss leaked {len(non_prod)}/200 "
            f"non-productive records after observed-stage 5' corruption."
        )

    def test_strict_no_non_productive_leaks_with_3prime_loss(self):
        sim = (
            Experiment.on("human_igh")
            .sequence(with_3prime_loss(
                min_remove=120, max_remove=180, min_add=0, max_add=0))
            .compile(seed=42, productivity=Productivity.PRODUCTIVE_ONLY)
        )
        _boost_retry_budget(sim, 200)
        records = sim.simulate(n=200)
        non_prod = [r for r in records if not r["productive"]]
        assert not non_prod, (
            f"PRODUCTIVE_ONLY + 3prime_loss leaked {len(non_prod)}/200 "
            f"non-productive records after observed-stage 3' corruption."
        )

    def test_strict_no_non_productive_leaks_with_observed_indels(self):
        sim = (
            Experiment.on("human_igh")
            .observe(with_indels(prob=0.15))
            .compile(seed=42, productivity=Productivity.PRODUCTIVE_ONLY)
        )
        _boost_retry_budget(sim, 200)
        records = sim.simulate(n=200)
        non_prod = [r for r in records if not r["productive"]]
        assert not non_prod, (
            f"PRODUCTIVE_ONLY + indels leaked {len(non_prod)}/200 "
            f"non-productive records after observed-stage indels."
        )

    def test_strict_no_non_productive_leaks_with_contaminants(self):
        sim = (
            Experiment.on("human_igh")
            .observe(with_contaminants(rate=1.0, source="random"))
            .compile(seed=42, productivity=Productivity.PRODUCTIVE_ONLY)
        )
        _boost_retry_budget(sim, 200)
        records = sim.simulate(n=200)
        non_prod = [r for r in records if not r["productive"]]
        assert not non_prod, (
            f"PRODUCTIVE_ONLY + contaminants leaked {len(non_prod)}/200 "
            f"non-productive records after contaminant injection."
        )


# ── NON_PRODUCTIVE_ONLY ───────────────────────────────────────────


class TestNonProductiveOnly:
    """Productivity.NON_PRODUCTIVE_ONLY: every record must be non-productive."""

    def test_strict_no_productive_leaks(self):
        sim = Experiment.on("human_igh").compile(
            seed=42, productivity=Productivity.NON_PRODUCTIVE_ONLY)
        _boost_retry_budget(sim, 200)
        records = sim.simulate(n=200)
        prod = [r for r in records if r["productive"]]
        assert not prod, (
            f"NON_PRODUCTIVE_ONLY leaked {len(prod)}/200 productive records."
        )

    def test_works_via_run_method(self):
        result = (Experiment.on("human_igh")
                  .run(n=50, seed=42,
                       productivity=Productivity.NON_PRODUCTIVE_ONLY))
        # Non-productive baseline ≈ 78%; 25 retries → exhaustion ≈ 5e-3
        # per record. For n=50 expect ≤1 leak (very generous).
        prod = sum(1 for r in result if r["productive"])
        assert prod <= 1, f"too many leaks: {prod}/50"


# ── PRODUCTIVE_MIXED ──────────────────────────────────────────────


class TestProductiveMixed:
    """Productivity.PRODUCTIVE_MIXED: unfiltered output — both classes
    must appear in any reasonable sample."""

    def test_mixed_contains_both_classes(self):
        result = (Experiment.on("human_igh")
                  .run(n=200, seed=42,
                       productivity=Productivity.PRODUCTIVE_MIXED))
        records = list(result)
        prod = sum(1 for r in records if r["productive"])
        non_prod = len(records) - prod
        assert prod > 0,     f"MIXED produced no productive sequences (0/{len(records)})"
        assert non_prod > 0, f"MIXED produced no non-productive sequences (0/{len(records)})"

    def test_mixed_ratio_matches_biology(self):
        """Pure recombination on human_igh is ~22.8% productive (T1-15).
        Loose band [0.10, 0.50] absorbs MC noise on n=500."""
        result = (Experiment.on("human_igh")
                  .run(n=500, seed=42,
                       productivity=Productivity.PRODUCTIVE_MIXED))
        records = list(result)
        prod_frac = sum(1 for r in records if r["productive"]) / len(records)
        assert 0.10 < prod_frac < 0.50, (
            f"productive fraction {prod_frac:.3f} outside expected band "
            f"[0.10, 0.50] for unfiltered IGH"
        )

    def test_mixed_is_default_no_args(self):
        """Calling .run() with no productivity arg uses MIXED — same
        ratio as explicit MIXED."""
        result = Experiment.on("human_igh").run(n=200, seed=42)
        records = list(result)
        prod = sum(1 for r in records if r["productive"])
        non_prod = len(records) - prod
        assert prod > 0 and non_prod > 0


# ── Backward compatibility (legacy bool API) ──────────────────────


class TestBackwardCompat:
    """Legacy `productive=True/False` continues to work and maps onto
    the enum exactly."""

    def test_productive_true_maps_to_productive_only(self):
        result = Experiment.on("human_igh").run(n=20, seed=42, productive=True)
        # Same as PRODUCTIVE_ONLY with default retries
        non_prod = sum(1 for r in result if not r["productive"])
        assert non_prod <= 1

    def test_productive_false_maps_to_mixed(self):
        result = Experiment.on("human_igh").run(n=200, seed=42, productive=False)
        records = list(result)
        prod = sum(1 for r in records if r["productive"])
        non_prod = len(records) - prod
        assert prod > 0 and non_prod > 0

    def test_passing_both_raises(self):
        with pytest.raises(ValueError, match="not both"):
            Experiment.on("human_igh").run(
                n=1, seed=42,
                productivity=Productivity.PRODUCTIVE_ONLY,
                productive=True,
            )

    def test_invalid_productivity_type_raises(self):
        with pytest.raises(TypeError, match="Productivity"):
            Experiment.on("human_igh").run(
                n=1, seed=42, productivity="productive_only")  # str, not enum


# ── Native stage-status capture (internal, no DSL surface yet) ───────────


class TestNativeStageStatusCapture:
    """Phase-2 native snapshotting should be visible in the trace
    without changing the public Python DSL."""

    def test_trace_records_rearrangement_status_capture(self):
        sim = Experiment.on("human_igh").compile(seed=42)
        sim._sim.set_trace(True)
        try:
            sim._sim.simulate_one()
            trace = sim._sim.get_trace()
        finally:
            sim._sim.set_trace(False)

        rearr = trace.index("[productivity] captured stage=rearrangement")
        molecule = trace.index("[productivity] captured stage=molecule")
        observed = trace.index("[productivity] captured stage=observed")

        assert rearr < molecule < observed
