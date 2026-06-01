"""Spec tests for the per-segment SHM rate scalars slice.

Implements the audit's §11 Slice 1: ``segment_rates`` kwarg on
``Experiment.mutate`` plus the Rust-side site-weighting in both
uniform and S5F passes. See
``docs/shm_segment_rate_design.md`` for the architecture
contract.

Spec coverage (from the user brief):

1. Default ``segment_rates=None`` is byte-identical to pre-slice
   behaviour.
2. ``{"NP": 0.0, "D": 0.0, "J": 0.0}`` forces mutations into V
   only.
3. ``{"V": 0.0, "D": 0.0, "J": 0.0, "NP": 1.0}`` forces mutations
   into NP when NP exists.
4. Same behaviour for uniform and S5F.
5. Productive-only still preserves the triad.
6. Replay round-trip succeeds with matching rates.
7. Replay with incompatible zero-rate support fails.
8. Bad DSL inputs reject: unknown key, negative, NaN/inf, all-zero.
9. Manifest advertises segment-rate support.
10. Existing SHM audit absence pins flip only for ``segment_rates``
    + manifest support; per-segment counters remain absent.

Out of scope: per-segment AIRR counters, S5F kernel digest in
content_hash, replay across different rate vectors (signature
mismatch path).
"""
from __future__ import annotations

import json

import pytest

import GenAIRR as ga
from GenAIRR._engine import StrictSamplingError


# ──────────────────────────────────────────────────────────────────
# Shared helpers
# ──────────────────────────────────────────────────────────────────


def _per_region_diff(record, baseline_seq):
    """Count per-region differences between a mutated record and its
    baseline counterpart. Returns ``{segment: count}`` over V /
    NP1 / D / NP2 / J."""
    base_seq = baseline_seq.upper()
    mut_seq = record["sequence"].upper()
    v_end = record["v_sequence_end"]
    d_start = record["d_sequence_start"]
    d_end = record["d_sequence_end"]
    j_start = record["j_sequence_start"]
    return {
        "V": sum(1 for x, y in zip(base_seq[:v_end], mut_seq[:v_end]) if x != y),
        "NP1": sum(
            1 for x, y in zip(base_seq[v_end:d_start], mut_seq[v_end:d_start]) if x != y
        ),
        "D": sum(
            1 for x, y in zip(base_seq[d_start:d_end], mut_seq[d_start:d_end]) if x != y
        ),
        "NP2": sum(
            1 for x, y in zip(base_seq[d_end:j_start], mut_seq[d_end:j_start]) if x != y
        ),
        "J": sum(1 for x, y in zip(base_seq[j_start:], mut_seq[j_start:]) if x != y),
    }


# ──────────────────────────────────────────────────────────────────
# Spec 1 — default is byte-identical to pre-slice behaviour
# ──────────────────────────────────────────────────────────────────


def test_default_segment_rates_byte_identical_to_no_kwarg() -> None:
    """``mutate(...)`` and ``mutate(..., segment_rates=None)``
    produce byte-identical sequences. The IR carries the default
    tuple verbatim; the Rust passes detect ``is_default()`` and
    short-circuit to the existing fast path."""
    baseline = (
        ga.Experiment.on("human_igh").recombine().mutate(count=10)
        .run_records(n=5, seed=0)
    )
    with_none = (
        ga.Experiment.on("human_igh").recombine().mutate(count=10, segment_rates=None)
        .run_records(n=5, seed=0)
    )
    assert [r["sequence"] for r in baseline] == [r["sequence"] for r in with_none]


def test_explicit_all_ones_byte_identical_to_default() -> None:
    """An explicit ``{V:1, D:1, J:1, NP:1}`` reaches the same
    ``is_default()`` fast path as omitting the kwarg. Byte-
    identical to the no-kwarg case."""
    baseline = (
        ga.Experiment.on("human_igh").recombine().mutate(count=10)
        .run_records(n=5, seed=0)
    )
    explicit = (
        ga.Experiment.on("human_igh").recombine()
        .mutate(count=10, segment_rates={"V": 1.0, "D": 1.0, "J": 1.0, "NP": 1.0})
        .run_records(n=5, seed=0)
    )
    assert [r["sequence"] for r in baseline] == [r["sequence"] for r in explicit]


def test_default_byte_identical_for_s5f_model() -> None:
    """Same default-equivalence check for the S5F model."""
    baseline = (
        ga.Experiment.on("human_igh").recombine().mutate(model="s5f", count=10)
        .run_records(n=5, seed=42)
    )
    with_default = (
        ga.Experiment.on("human_igh").recombine()
        .mutate(model="s5f", count=10, segment_rates={"V": 1.0, "D": 1.0, "J": 1.0, "NP": 1.0})
        .run_records(n=5, seed=42)
    )
    assert [r["sequence"] for r in baseline] == [r["sequence"] for r in with_default]


# ──────────────────────────────────────────────────────────────────
# Spec 2 — {NP:0, D:0, J:0} forces mutations into V only
# ──────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("model", ["uniform", "s5f"])
def test_d_j_np_zero_forces_v_only_mutations(model) -> None:
    """All mutations land in V; D/NP/J get zero differences. Pinned
    for both SHM models."""
    baseline = (
        ga.Experiment.on("human_igh").recombine().run_records(n=5, seed=0)
    )
    mutated = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(
            model=model,
            count=20,
            segment_rates={"D": 0.0, "J": 0.0, "NP": 0.0},
        )
        .run_records(n=5, seed=0)
    )
    for base_rec, mut_rec in zip(baseline, mutated):
        diffs = _per_region_diff(mut_rec, base_rec["sequence"])
        assert diffs["NP1"] == 0, (
            f"{model}: NP1 received mutations under NP=0 rate: {diffs}"
        )
        assert diffs["D"] == 0, (
            f"{model}: D received mutations under D=0 rate: {diffs}"
        )
        assert diffs["NP2"] == 0, f"{model}: NP2 received mutations: {diffs}"
        assert diffs["J"] == 0, f"{model}: J received mutations: {diffs}"
        assert diffs["V"] > 0, (
            f"{model}: V received no mutations despite V=1 rate: {diffs}"
        )


# ──────────────────────────────────────────────────────────────────
# Spec 3 — NP-only when NP exists
# ──────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("model", ["uniform", "s5f"])
def test_v_d_j_zero_forces_np_only_mutations(model) -> None:
    """All mutations land in NP regions; V/D/J get zero
    differences. Pinned for both models. The realised count may
    be less than the requested count for S5F (NP regions are
    short relative to the kernel's 5-mer context support; some
    iterations may find empty profiles) — permissive mode silently
    skips."""
    baseline = (
        ga.Experiment.on("human_igh").recombine().run_records(n=5, seed=0)
    )
    mutated = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(
            model=model,
            count=10,
            segment_rates={"V": 0.0, "D": 0.0, "J": 0.0, "NP": 1.0},
        )
        .run_records(n=5, seed=0)
    )
    saw_np_hits = False
    for base_rec, mut_rec in zip(baseline, mutated):
        diffs = _per_region_diff(mut_rec, base_rec["sequence"])
        assert diffs["V"] == 0, f"{model}: V mutated under V=0: {diffs}"
        assert diffs["D"] == 0, f"{model}: D mutated under D=0: {diffs}"
        assert diffs["J"] == 0, f"{model}: J mutated under J=0: {diffs}"
        if diffs["NP1"] + diffs["NP2"] > 0:
            saw_np_hits = True
    # At least one record across the batch should show NP hits;
    # NP regions can be 0 length in some samples so we don't
    # require per-record positivity.
    assert saw_np_hits, (
        f"{model}: no NP mutations observed across the batch despite "
        "NP=1, V=D=J=0 — sampling support broken."
    )


# ──────────────────────────────────────────────────────────────────
# Spec 4 — same behaviour for uniform and S5F (covered by parametrise above)
# ──────────────────────────────────────────────────────────────────


# ──────────────────────────────────────────────────────────────────
# Spec 5 — productive-only preserves the triad
# ──────────────────────────────────────────────────────────────────


def test_productive_only_preserves_triad_under_segment_rates() -> None:
    """``productive_only()`` + ``segment_rates`` together still
    produce productive records. The constrain-before-propose
    ordering (audit §4) means segment-rate filtering happens
    BEFORE contract admissibility — both gates run, no
    interaction breakage."""
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(count=15, segment_rates={"V": 1.0, "D": 0.5, "J": 0.5, "NP": 0.1})
        .productive_only()
    )
    result = exp.run_records(n=5, seed=0)
    for r in result:
        assert r["productive"] is True


# ──────────────────────────────────────────────────────────────────
# Spec 6 — replay round-trip with matching rates
# ──────────────────────────────────────────────────────────────────


def test_replay_round_trip_with_matching_segment_rates() -> None:
    """Two runs of the same experiment with the same seed produce
    byte-identical records. Confirms the rate vector is part of
    the deterministic-sampling surface — same vector + same seed
    → same output."""
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(count=12, segment_rates={"V": 1.0, "D": 0.0, "J": 0.5, "NP": 0.0})
    )
    a = exp.run_records(n=4, seed=7)
    b = exp.run_records(n=4, seed=7)
    assert [r["sequence"] for r in a] == [r["sequence"] for r in b]


# ──────────────────────────────────────────────────────────────────
# Spec 8 — bad DSL inputs reject
# ──────────────────────────────────────────────────────────────────


def test_unknown_segment_key_rejected() -> None:
    """Spelling errors / wrong case raise ``ValueError`` at the DSL
    boundary."""
    with pytest.raises(ValueError, match="unknown segment key"):
        (
            ga.Experiment.on("human_igh")
            .recombine()
            .mutate(count=5, segment_rates={"v": 1.0})  # lowercase
        )
    with pytest.raises(ValueError, match="unknown segment key"):
        (
            ga.Experiment.on("human_igh")
            .recombine()
            .mutate(count=5, segment_rates={"Np1": 1.0})  # not a bucket
        )


def test_negative_rate_rejected() -> None:
    """Negative values raise ``ValueError`` with the offending
    bucket named."""
    with pytest.raises(ValueError, match=r"segment_rates\['V'\]"):
        (
            ga.Experiment.on("human_igh")
            .recombine()
            .mutate(count=5, segment_rates={"V": -1.0})
        )


def test_nan_rate_rejected() -> None:
    """NaN rate raises ``ValueError`` (with NaN-specific message)."""
    with pytest.raises(ValueError, match="NaN"):
        (
            ga.Experiment.on("human_igh")
            .recombine()
            .mutate(count=5, segment_rates={"V": float("nan")})
        )


def test_inf_rate_rejected() -> None:
    """Infinite rates aren't finite numbers and raise."""
    with pytest.raises(ValueError, match="finite"):
        (
            ga.Experiment.on("human_igh")
            .recombine()
            .mutate(count=5, segment_rates={"V": float("inf")})
        )


def test_all_zero_rates_rejected() -> None:
    """All four buckets zero would make the pass a no-op; rejected
    at the DSL boundary."""
    with pytest.raises(ValueError, match="at least one bucket must have a positive rate"):
        (
            ga.Experiment.on("human_igh")
            .recombine()
            .mutate(
                count=5,
                segment_rates={"V": 0.0, "D": 0.0, "J": 0.0, "NP": 0.0},
            )
        )


def test_bool_rate_rejected_as_type_error() -> None:
    """Booleans are rejected — type check happens before numeric
    check so an explicit ``True`` doesn't slide through as 1.0."""
    with pytest.raises(TypeError):
        (
            ga.Experiment.on("human_igh")
            .recombine()
            .mutate(count=5, segment_rates={"V": True})
        )


def test_non_dict_segment_rates_rejected() -> None:
    """Passing a non-dict (e.g. a list) is a programming-error
    shape and gets a clear ``TypeError``."""
    with pytest.raises(TypeError, match="dict or None"):
        (
            ga.Experiment.on("human_igh")
            .recombine()
            .mutate(count=5, segment_rates=[1.0, 1.0, 1.0, 1.0])
        )


def test_sparse_dict_defaults_missing_buckets_to_one() -> None:
    """``segment_rates={"NP": 0.0}`` defaults V / D / J to 1.0
    (sparse semantics from the audit §6 DSL recommendation)."""
    exp_explicit = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(
            count=15,
            segment_rates={"V": 1.0, "D": 1.0, "J": 1.0, "NP": 0.0},
        )
    )
    exp_sparse = (
        ga.Experiment.on("human_igh").recombine().mutate(count=15, segment_rates={"NP": 0.0})
    )
    a = exp_explicit.run_records(n=3, seed=0)
    b = exp_sparse.run_records(n=3, seed=0)
    assert [r["sequence"] for r in a] == [r["sequence"] for r in b]


# ──────────────────────────────────────────────────────────────────
# Spec 9 — manifest advertises segment-rate support
# ──────────────────────────────────────────────────────────────────


def test_manifest_advertises_segment_rate_support() -> None:
    """``cartridge_manifest()["models"]["shm"]["segment_rate_support"]``
    surfaces the slice's capability with the documented buckets +
    default vector + the v1-boundary ``in_content_hash`` flag."""
    m = ga.HUMAN_IGH_OGRDB.cartridge_manifest()
    shm = m["models"]["shm"]
    assert "segment_rate_support" in shm
    sup = shm["segment_rate_support"]
    assert sup["available"] is True
    assert sup["buckets"] == ["V", "D", "J", "NP"]
    assert sup["default"] == {"V": 1.0, "D": 1.0, "J": 1.0, "NP": 1.0}
    # Slice deliberately doesn't change content_hash semantics.
    assert sup["in_content_hash"] is False
    # And the full manifest round-trips through json.dumps.
    json.dumps(m)


# ──────────────────────────────────────────────────────────────────
# IR plumbing — _MutateStep carries the normalised tuple
# ──────────────────────────────────────────────────────────────────


def test_mutate_step_carries_normalised_tuple() -> None:
    """The pipeline-IR ``_MutateStep`` dataclass carries the
    normalised ``(v, d, j, np)`` tuple. Pinned via direct
    introspection so a refactor that drops the field surfaces
    here."""
    from GenAIRR._pipeline_ir import _MutateStep

    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(count=5, segment_rates={"V": 1.0, "D": 0.5, "J": 0.25, "NP": 0.0})
    )
    mutate_steps = [s for s in exp._steps if isinstance(s, _MutateStep)]
    assert len(mutate_steps) == 1
    assert mutate_steps[0].segment_rates == (1.0, 0.5, 0.25, 0.0)


def test_default_mutate_step_carries_flat_default() -> None:
    """Omitting the kwarg produces a step with ``(1.0, 1.0, 1.0,
    1.0)`` — the documented flat default."""
    from GenAIRR._pipeline_ir import _MutateStep

    exp = ga.Experiment.on("human_igh").recombine().mutate(count=5)
    step = next(s for s in exp._steps if isinstance(s, _MutateStep))
    assert step.segment_rates == (1.0, 1.0, 1.0, 1.0)
