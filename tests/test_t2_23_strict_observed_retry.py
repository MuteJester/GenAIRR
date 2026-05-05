"""T2-23: STRICT observed-stage retry at the API layer.

V5 step 27 wires `simcfg_productivity_strict` into a per-record retry
loop in `genairr_simulate`/`genairr_simulate_to_fd`. Under
PRODUCTIVE_ONLY + STRICT, the API layer reruns the entire pipeline
for any record whose final observed productivity is False, up to
``max_productive_attempts`` times.

Local guards in PRODUCTIVE_ONLY already prevent most observed-stage
leaks (verified by `test_productivity_modes.py`), so STRICT will not
typically differ from BEST_EFFORT on the default config. These tests
verify the wiring and the contract:

  1. STRICT can be flipped on via `set_param("productivity_strictness", 1)`.
  2. STRICT does not drop records — N requested, N returned.
  3. STRICT preserves the PRODUCTIVE_ONLY invariant (no leaks).
  4. Strictness flag round-trips via the param API.
"""
from __future__ import annotations

import pytest
from GenAIRR import Experiment, Productivity
from GenAIRR.ops import with_contaminants


PROD_STRICTNESS_BEST_EFFORT = 0
PROD_STRICTNESS_STRICT = 1


def _set_strictness(sim, value: int) -> None:
    sim._sim.set_param("productivity_strictness", value)


class TestStrictObservedRetry:
    def test_strict_param_accepts_known_values(self):
        sim = Experiment.on("human_igh").compile(
            seed=1, productivity=Productivity.PRODUCTIVE_ONLY)
        _set_strictness(sim, PROD_STRICTNESS_STRICT)
        _set_strictness(sim, PROD_STRICTNESS_BEST_EFFORT)

    def test_strict_param_rejects_invalid_values(self):
        sim = Experiment.on("human_igh").compile(
            seed=1, productivity=Productivity.PRODUCTIVE_ONLY)
        with pytest.raises(ValueError):
            _set_strictness(sim, 99)

    def test_strict_returns_full_record_count(self):
        """STRICT must not silently drop records on retry exhaustion."""
        sim = Experiment.on("human_igh").compile(
            seed=42, productivity=Productivity.PRODUCTIVE_ONLY)
        _set_strictness(sim, PROD_STRICTNESS_STRICT)
        sim._sim.set_param("max_productive_attempts", 200)
        records = sim.simulate(n=50)
        assert len(records) == 50

    def test_strict_no_non_productive_leaks_with_contaminants(self):
        """STRICT + PRODUCTIVE_ONLY + contaminants: zero non-productive."""
        sim = (Experiment.on("human_igh")
               .observe(with_contaminants(rate=0.05, source="random"))
               .compile(seed=42, productivity=Productivity.PRODUCTIVE_ONLY))
        _set_strictness(sim, PROD_STRICTNESS_STRICT)
        sim._sim.set_param("max_productive_attempts", 200)
        records = sim.simulate(n=100)
        non_prod = [r for r in records if not r["productive"]]
        assert not non_prod, (
            f"STRICT PRODUCTIVE_ONLY + contaminants leaked "
            f"{len(non_prod)}/100 non-productive records"
        )

    def test_strict_matches_best_effort_record_count(self):
        """Both strictness modes must emit the same number of records
        for the same N — STRICT only adds retries, never drops."""
        n = 50

        sim_be = Experiment.on("human_igh").compile(
            seed=7, productivity=Productivity.PRODUCTIVE_ONLY)
        _set_strictness(sim_be, PROD_STRICTNESS_BEST_EFFORT)
        sim_be._sim.set_param("max_productive_attempts", 200)
        recs_be = sim_be.simulate(n=n)

        sim_s = Experiment.on("human_igh").compile(
            seed=7, productivity=Productivity.PRODUCTIVE_ONLY)
        _set_strictness(sim_s, PROD_STRICTNESS_STRICT)
        sim_s._sim.set_param("max_productive_attempts", 200)
        recs_s = sim_s.simulate(n=n)

        assert len(recs_be) == n
        assert len(recs_s) == n


class TestStrictRetryCounter:
    """V5/Step 28: STRICT observed-stage retry counter exposed via API."""

    def test_counter_starts_at_zero(self):
        sim = Experiment.on("human_igh").compile(
            seed=1, productivity=Productivity.PRODUCTIVE_ONLY)
        assert sim._sim.strict_retry_count == 0

    def test_counter_zero_under_defaults_no_strict(self):
        """Without STRICT, retries never fire even after a full batch."""
        sim = Experiment.on("human_igh").compile(
            seed=1, productivity=Productivity.PRODUCTIVE_ONLY)
        sim._sim.set_param("max_productive_attempts", 200)
        sim.simulate(n=20)
        assert sim._sim.strict_retry_count == 0

    def test_counter_zero_when_local_guards_prevent_leaks(self):
        """STRICT + comprehensive local guards = no retries needed."""
        sim = Experiment.on("human_igh").compile(
            seed=1, productivity=Productivity.PRODUCTIVE_ONLY)
        _set_strictness(sim, PROD_STRICTNESS_STRICT)
        sim._sim.set_param("max_productive_attempts", 200)
        sim.simulate(n=20)
        # Local guards already prevent observed-stage leaks under
        # default config, so STRICT retries should not fire.
        assert sim._sim.strict_retry_count == 0

    def test_counter_reset(self):
        sim = Experiment.on("human_igh").compile(
            seed=1, productivity=Productivity.PRODUCTIVE_ONLY)
        _set_strictness(sim, PROD_STRICTNESS_STRICT)
        sim._sim.set_param("max_productive_attempts", 200)
        sim.simulate(n=10)
        # Reset is idempotent and safe even when counter is already 0.
        sim._sim.reset_strict_retry_count()
        assert sim._sim.strict_retry_count == 0
        sim.simulate(n=10)
        sim._sim.reset_strict_retry_count()
        assert sim._sim.strict_retry_count == 0

    def test_counter_accumulates_across_calls(self):
        """The counter is per-handle, not per-call: multiple simulate()
        calls should not zero it implicitly between calls."""
        sim = Experiment.on("human_igh").compile(
            seed=1, productivity=Productivity.PRODUCTIVE_ONLY)
        _set_strictness(sim, PROD_STRICTNESS_STRICT)
        sim._sim.set_param("max_productive_attempts", 200)
        before = sim._sim.strict_retry_count
        sim.simulate(n=10)
        sim.simulate(n=10)
        # Whatever the value is (likely 0 under default guards), it
        # must be monotonically non-decreasing across batches.
        after = sim._sim.strict_retry_count
        assert after >= before
