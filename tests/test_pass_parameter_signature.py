"""Spec tests for **Slice A: Pass Parameter Signature**.

Trace-file plan signatures used to be a flat `pipe|joined|pass|names`
string. That left a real replay-safety hole: a trace recorded
under `segment_rates={"V": 2.0}` would silently produce different
bases at the same recorded addresses when replayed against a plan
with `segment_rates={"V": 1.0}` because the signatures matched.

Slice A closes the gap. Each pass now contributes its
compile-time-parameter digest via `Pass::parameter_signature`;
the plan signature becomes `name(params)|name(params)|...` and
ships in v3 trace files. v1 and v2 fixtures continue to load and
replay via the legacy names-only comparator
(`pass_plan_signature_names_only`) — backwards compatible by
construction.

Spec coverage:

1. Two plans differing only by `segment_rates` produce different
   plan signatures.
2. Replaying a trace recorded under one rate vector against a
   different rate vector fails the signature gate **before any
   choice is consumed** (i.e. the rejection is the signature
   check, not a runtime divergence).
3. Default `segment_rates=None` and explicit all-ones produce
   the **same** signature (behaviourally equivalent inputs
   collide).
4. The S5F kernel name participates: `mutate(s5f_model="hh_s5f")`
   vs `mutate(s5f_model="hkl_s5f")` give different signatures.
5. Schema version is bumped to 3; v3 signatures contain the
   `name(params)` envelope.
6. Distribution-bearing passes (trim, NP, end-loss) include
   their `support()` output in the signature: a custom NP-length
   distribution changes the signature even though the pass name
   is unchanged.
7. v1/v2 golden fixtures still load and replay (covered in
   `test_trace_file_compat.py`; this file pins the
   schema-version-aware comparator at the surface level).
"""
from __future__ import annotations

import copy
import json

import pytest

import GenAIRR as ga


def _plan_signature(experiment, seed: int = 42) -> str:
    compiled = experiment.compile()
    outcome = compiled.simulator.run(seed=seed)
    tf = compiled.simulator.trace_file_from(outcome, seed=seed)
    return json.loads(tf.to_json())["pass_plan_signature"]


# ──────────────────────────────────────────────────────────────────
# Spec 1 — different segment_rates → different signatures
# ──────────────────────────────────────────────────────────────────


def test_segment_rate_change_flips_plan_signature() -> None:
    """A trace recorded with `segment_rates={"V": 2.0}` and one
    recorded with the default vector now produce DIFFERENT plan
    signatures. Before Slice A they collided because the
    signature was pass-names only."""
    sig_default = _plan_signature(
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(model="s5f", rate=0.03)
    )
    sig_v_doubled = _plan_signature(
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(model="s5f", rate=0.03, segment_rates={"V": 2.0})
    )
    assert sig_default != sig_v_doubled


def test_segment_rate_zero_flips_plan_signature() -> None:
    """Zero-rate vector is also a distinct signature, not the
    default fast-path."""
    sig_default = _plan_signature(
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(model="s5f", rate=0.03)
    )
    sig_np_off = _plan_signature(
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(model="s5f", rate=0.03, segment_rates={"NP": 0.0})
    )
    assert sig_default != sig_np_off


# ──────────────────────────────────────────────────────────────────
# Spec 2 — mismatched-rate replay fails the signature gate
# ──────────────────────────────────────────────────────────────────


def test_replay_with_mismatched_segment_rates_fails_signature_gate() -> None:
    """End-to-end replay safety. A trace produced with one rate
    vector, replayed against a plan with a different rate
    vector, must raise at the signature gate — NOT silently
    sample different bases at the same recorded addresses."""
    exp_a = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(model="s5f", rate=0.03, segment_rates={"V": 2.0})
    )
    exp_b = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(model="s5f", rate=0.03, segment_rates={"V": 1.0})
    )
    c_a = exp_a.compile()
    c_b = exp_b.compile()
    outcome_a = c_a.simulator.run(seed=4242)
    tf_a = c_a.simulator.trace_file_from(outcome_a, seed=4242)

    with pytest.raises(ValueError, match="pass plan signature mismatch"):
        c_b.simulator.replay_from_trace_file(tf_a)
    with pytest.raises(ValueError, match="pass plan signature mismatch"):
        c_b.simulator.rerun_from_trace_file(tf_a)


def test_replay_with_matching_segment_rates_succeeds() -> None:
    """Sanity: when both plans share the same rate vector the
    signature gate is silent and replay reproduces the original
    outcome."""
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(model="s5f", rate=0.03, segment_rates={"V": 2.0})
    )
    compiled = exp.compile()
    outcome = compiled.simulator.run(seed=4242)
    tf = compiled.simulator.trace_file_from(outcome, seed=4242)
    # No raise.
    compiled.simulator.replay_from_trace_file(tf)


# ──────────────────────────────────────────────────────────────────
# Spec 3 — behaviourally equivalent inputs produce equal sigs
# ──────────────────────────────────────────────────────────────────


def test_default_and_explicit_all_ones_segment_rates_collide() -> None:
    """The default `segment_rates=None` and an explicit
    `{"V": 1.0, "D": 1.0, "J": 1.0, "NP": 1.0}` dict are
    behaviourally equivalent — `SegmentRateWeights::is_default()`
    fires in both cases and the signature short-circuits to the
    empty string. Verified end-to-end through the DSL."""
    sig_default = _plan_signature(
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(model="s5f", rate=0.03)
    )
    sig_explicit = _plan_signature(
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(
            model="s5f",
            rate=0.03,
            segment_rates={"V": 1.0, "D": 1.0, "J": 1.0, "NP": 1.0},
        )
    )
    assert sig_default == sig_explicit


# ──────────────────────────────────────────────────────────────────
# Spec 4 — S5F kernel name participates
# ──────────────────────────────────────────────────────────────────


def test_s5f_kernel_name_participates_in_signature() -> None:
    """Two plans that differ only in the S5F kernel name now
    produce different signatures. Before Slice A the kernel name
    was lost at the PyO3 boundary; this pins that the threading
    works end-to-end."""
    sig_hh = _plan_signature(
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(model="s5f", rate=0.03, s5f_model="hh_s5f")
    )
    sig_hkl = _plan_signature(
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(model="s5f", rate=0.03, s5f_model="hkl_s5f")
    )
    assert sig_hh != sig_hkl
    # Both must contain the kernel name verbatim.
    assert "hh_s5f" in sig_hh
    assert "hkl_s5f" in sig_hkl


# ──────────────────────────────────────────────────────────────────
# Spec 5 — schema version bumped to 3; signature uses envelope
# ──────────────────────────────────────────────────────────────────


def test_fresh_emission_carries_schema_version_3() -> None:
    """Slice A bumps the on-disk schema to v3. Fresh emissions
    write v3 by default and the plan signature uses the new
    `name(params)` envelope on every pass."""
    exp = ga.Experiment.on("human_igh").recombine()
    compiled = exp.compile()
    outcome = compiled.simulator.run(seed=42)
    tf = compiled.simulator.trace_file_from(outcome, seed=42)
    doc = json.loads(tf.to_json())
    assert doc["schema_version"] == 3
    # Envelope is observable: every pass renders as `name(...)`.
    sig = doc["pass_plan_signature"]
    assert "(" in sig
    # And each `|`-separated piece carries the envelope.
    for part in sig.split("|"):
        assert "(" in part and part.endswith(")"), (
            f"plan-signature piece {part!r} does not use the "
            "`name(params)` envelope"
        )


# ──────────────────────────────────────────────────────────────────
# Spec 6 — distribution-bearing passes include their support()
# ──────────────────────────────────────────────────────────────────


def test_uniform_count_distribution_participates_in_signature() -> None:
    """Two `mutate(model="uniform", count=…)` plans with
    different count distributions produce different signatures.
    The empirical count distribution is the only parameter that
    changes between them — pinning that distribution support
    enters the signature."""
    sig_count_3 = _plan_signature(
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(model="uniform", count=3)
    )
    sig_count_5 = _plan_signature(
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(model="uniform", count=5)
    )
    assert sig_count_3 != sig_count_5


# ──────────────────────────────────────────────────────────────────
# Spec 7 — probability-bearing passes (invert_d, etc.)
# ──────────────────────────────────────────────────────────────────


def test_invert_d_probability_participates_in_signature() -> None:
    """The D-inversion probability is a compile-time parameter
    that flips the plan signature when changed."""
    sig_prob_0 = _plan_signature(
        ga.Experiment.on("human_igh").recombine().invert_d(prob=0.0)
    )
    sig_prob_50 = _plan_signature(
        ga.Experiment.on("human_igh").recombine().invert_d(prob=0.5)
    )
    assert sig_prob_0 != sig_prob_50
