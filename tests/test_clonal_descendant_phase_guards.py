"""Spec tests for the unified clonal descendant-phase ordering guard.

This slice extends the earlier ``invert_d`` / ``receptor_revision``
/ ``paired_end`` ordering-guard work to **every** descendant-phase
DSL step. Pre-fork placement of these steps either silently
misreports the AIRR field (Bugs C / E / F — trace-sourced fields
that don't survive the parent→descendant boundary) or collapses
descendant diversity (the pass runs once on the parent IR and
every clone member inherits an identical effect). Either failure
mode is a clonal-semantics violation; the guard at
:meth:`Experiment.expand_clones` rejects them at the DSL boundary.

Guarded descendant-phase steps:

- ``mutate(...)``
- ``pcr_amplify(...)``
- ``ambiguous_base_calls(...)``
- ``sequencing_errors(...)``
- ``polymerase_indels(...)``
- ``end_loss_5prime(...)``
- ``end_loss_3prime(...)``
- ``random_strand_orientation(...)``
- ``paired_end(...)`` — flipped from the dedicated check to the
  unified table; behaviour identical.

Unguarded (left as a follow-up decision): ``contaminate(...)``.

The recombination-/ancestor-phase guards from the previous slice
(``invert_d``, ``receptor_revision``) stay method-level — they
enforce the opposite direction (pre-fork only) and are tested
separately in ``tests/test_clonal_dsl_ordering_guards.py``.

Spec coverage (from the user brief):

1. Each descendant-phase step before ``expand_clones()`` raises.
2. Each same step after ``expand_clones()`` works.
3. Non-clonal pipelines with those steps still work.
4. Bugs E / F repros now fail at the DSL boundary.
5. Existing paired-end guard still works through the unified guard.
6. Canonical full clonal pipeline works end-to-end.
"""
from __future__ import annotations

import pytest

import GenAIRR as ga


# ──────────────────────────────────────────────────────────────────
# Shared fixture builders
# ──────────────────────────────────────────────────────────────────


def _base():
    return ga.Experiment.on("human_igh").recombine()


# ``builder`` keys must match the DSL method names; each value is a
# callable ``Experiment -> Experiment`` that appends the step in
# question. Used by the parametrised tests below to drive 9
# identical assertion shapes without copying boilerplate.
_DESCENDANT_PHASE_BUILDERS = {
    "mutate": lambda e: e.mutate(count=5),
    "pcr_amplify": lambda e: e.pcr_amplify(count=3),
    "ambiguous_base_calls": lambda e: e.ambiguous_base_calls(count=4),
    "sequencing_errors": lambda e: e.sequencing_errors(count=3),
    "polymerase_indels": lambda e: e.polymerase_indels(count=2),
    "end_loss_5prime": lambda e: e.end_loss_5prime(length=10),
    "end_loss_3prime": lambda e: e.end_loss_3prime(length=10),
    "random_strand_orientation": lambda e: e.random_strand_orientation(prob=1.0),
    "paired_end": lambda e: e.paired_end(r1_length=150, insert_size=300),
}


# ──────────────────────────────────────────────────────────────────
# Spec 1 — each descendant-phase step pre-fork raises
# ──────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("method", sorted(_DESCENDANT_PHASE_BUILDERS))
def test_each_descendant_phase_step_pre_fork_raises_at_expand_clones(method) -> None:
    """The guard fires at the ``expand_clones()`` call site, not at
    the descendant-phase method itself (which must still work in
    non-clonal pipelines). Verify the error message names the
    offending DSL method, the correct ordering, and the canonical
    fix instruction."""
    builder = _DESCENDANT_PHASE_BUILDERS[method]
    exp = builder(_base())
    with pytest.raises(ValueError) as exc_info:
        exp.expand_clones(n_clones=1, per_clone=2)
    msg = str(exc_info.value)
    # The error must name the offending DSL method.
    assert method in msg, (
        f"error message for {method} pre-fork doesn't name the step: {msg!r}"
    )
    # The correct ordering ("after expand_clones").
    assert "after expand_clones" in msg, msg
    # The canonical fix instruction.
    assert "Move" in msg or "move" in msg, msg


def test_mutate_message_uses_shm_specific_phrasing() -> None:
    """``mutate`` gets the more precise message naming SHM as the
    biology — the user spec explicitly asked for this."""
    exp = _base().mutate(count=5)
    with pytest.raises(ValueError) as exc_info:
        exp.expand_clones(n_clones=1, per_clone=2)
    msg = str(exc_info.value).lower()
    assert "shm" in msg, f"mutate message lacks SHM phrasing: {msg!r}"
    assert "descendant-specific" in msg


def test_other_descendant_phase_messages_use_uniform_phrasing() -> None:
    """All non-mutate descendant-phase steps share the uniform
    ``descendant-specific and must be sampled independently for
    each clone member`` phrasing."""
    for method, builder in _DESCENDANT_PHASE_BUILDERS.items():
        if method == "mutate":
            continue
        exp = builder(_base())
        with pytest.raises(ValueError) as exc_info:
            exp.expand_clones(n_clones=1, per_clone=2)
        msg = str(exc_info.value)
        assert "descendant-specific" in msg
        assert "sampled independently for each clone member" in msg, (
            f"{method}: unexpected message phrasing: {msg!r}"
        )


# ──────────────────────────────────────────────────────────────────
# Spec 2 — each step works after expand_clones
# ──────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("method", sorted(_DESCENDANT_PHASE_BUILDERS))
def test_each_descendant_phase_step_post_fork_works(method) -> None:
    """Placed after ``expand_clones()`` (the canonical order), every
    descendant-phase step compiles and runs without error.

    ``random_strand_orientation`` is BCR-only in the sense that we
    use ``human_igh`` here — the call itself is locus-agnostic."""
    builder = _DESCENDANT_PHASE_BUILDERS[method]
    exp = builder(_base().expand_clones(n_clones=1, per_clone=2))
    result = exp.run_records(seed=0)
    assert len(result) == 2


# ──────────────────────────────────────────────────────────────────
# Spec 3 — non-clonal pipelines unaffected
# ──────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("method", sorted(_DESCENDANT_PHASE_BUILDERS))
def test_non_clonal_descendant_phase_steps_unchanged(method) -> None:
    """The guard is at ``expand_clones()``, not at the descendant-
    phase methods themselves. Non-clonal pipelines must continue
    to use these steps freely."""
    builder = _DESCENDANT_PHASE_BUILDERS[method]
    exp = builder(_base())
    result = exp.run_records(n=2, seed=0)
    assert len(result) == 2


# ──────────────────────────────────────────────────────────────────
# Spec 4 — Bugs E and F repros now fail at the DSL boundary
# ──────────────────────────────────────────────────────────────────


def test_bug_e_random_strand_pre_fork_now_rejected() -> None:
    """Bug E: ``random_strand_orientation`` pre-fork used to
    silently produce ``rev_comp=False`` on every descendant (the
    rev_comp trace event landed on the parent only; descendant
    projection couldn't find it). Now the misorder is rejected at
    the DSL boundary."""
    exp = _base().random_strand_orientation(prob=1.0)
    with pytest.raises(ValueError, match="random_strand_orientation"):
        exp.expand_clones(n_clones=1, per_clone=2)


def test_bug_f_end_loss_5prime_pre_fork_now_rejected() -> None:
    """Bug F: ``end_loss_5prime`` pre-fork truncated the parent's
    pool (descendants inherited the truncation) but the AIRR field
    ``end_loss_5_length`` projected as 0 on descendants — the
    trace event landed on the parent. Now the misorder is rejected
    at the DSL boundary."""
    exp = _base().end_loss_5prime(length=20)
    with pytest.raises(ValueError, match="end_loss_5prime"):
        exp.expand_clones(n_clones=1, per_clone=2)


def test_bug_f_end_loss_3prime_pre_fork_now_rejected() -> None:
    """Sibling of the 5' case."""
    exp = _base().end_loss_3prime(length=15)
    with pytest.raises(ValueError, match="end_loss_3prime"):
        exp.expand_clones(n_clones=1, per_clone=2)


# ──────────────────────────────────────────────────────────────────
# Spec 5 — paired-end guard still works through the unified guard
# ──────────────────────────────────────────────────────────────────


def test_existing_paired_end_guard_still_fires_via_unified_table() -> None:
    """The paired-end ordering check moved from a dedicated branch
    in ``expand_clones`` to the unified descendant-phase table.
    Behaviour must be identical: ``recombine().paired_end(...).
    expand_clones(...)`` raises with the paired_end-named
    message."""
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .paired_end(r1_length=150, insert_size=300)
    )
    with pytest.raises(ValueError) as exc_info:
        exp.expand_clones(n_clones=1, per_clone=2)
    msg = str(exc_info.value)
    assert "paired_end" in msg
    assert "after expand_clones" in msg


# ──────────────────────────────────────────────────────────────────
# Spec 6 — canonical full clonal pipeline works
# ──────────────────────────────────────────────────────────────────


def test_canonical_full_clonal_pipeline_works() -> None:
    """The user-spec canonical pipeline assembles ALL descendant-
    phase steps after the fork. End-to-end pin: compile + run +
    correctness of inherited / descendant-specific fields."""
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .invert_d(prob=1.0)
        .receptor_revision(prob=1.0)
        .expand_clones(n_clones=2, per_clone=3)
        .mutate(count=5)
        .end_loss_5prime(length=10)
        .end_loss_3prime(length=10)
        .random_strand_orientation(prob=1.0)
        .paired_end(r1_length=150, insert_size=300)
    )
    result = exp.run_records(seed=0)
    assert len(result) == 6
    for r in result:
        # Ancestor-phase decisions inherited by every descendant.
        assert r["d_inverted"] is True
        assert r["receptor_revision_applied"] is True
        assert r["original_v_call"] != ""
        # Descendant-phase effects actually projected (this is the
        # closing-condition Bug E/F headline pin).
        assert r["rev_comp"] is True
        assert r["end_loss_5_length"] >= 1
        assert r["end_loss_3_length"] >= 1
        assert len(r["r1_sequence"]) == 150
        assert len(r["r2_sequence"]) == 150


# ──────────────────────────────────────────────────────────────────
# Cross-slice — recombination-phase guards still fire (the opposite
# direction). These pins ensure the unified descendant-phase guard
# didn't accidentally weaken the ancestor-phase rules.
# ──────────────────────────────────────────────────────────────────


def test_invert_d_after_expand_clones_still_rejected() -> None:
    """Ancestor-phase ``invert_d`` after ``expand_clones`` stays
    rejected by the method-level guard from the previous slice.
    Pinned cross-slice so a refactor that moves the ancestor-phase
    guards into the unified table can't silently swap directions."""
    exp = _base().expand_clones(n_clones=1, per_clone=2)
    with pytest.raises(ValueError, match="invert_d must be called before"):
        exp.invert_d(prob=1.0)


def test_receptor_revision_after_expand_clones_still_rejected() -> None:
    """Sibling of the invert_d cross-slice pin."""
    exp = _base().expand_clones(n_clones=1, per_clone=2)
    with pytest.raises(
        ValueError, match="receptor_revision must be called before"
    ):
        exp.receptor_revision(prob=1.0)


# ──────────────────────────────────────────────────────────────────
# Guard fires on the FIRST offending step (not the last)
# ──────────────────────────────────────────────────────────────────


def test_guard_reports_first_descendant_phase_step_encountered() -> None:
    """When multiple descendant-phase steps are placed pre-fork,
    the guard names the FIRST one encountered in the step list —
    so the message points the user at the call they need to look
    at first, not an arbitrary later one."""
    exp = (
        _base()
        .mutate(count=5)  # FIRST offender
        .pcr_amplify(count=3)
        .paired_end(r1_length=100, insert_size=200)
    )
    with pytest.raises(ValueError) as exc_info:
        exp.expand_clones(n_clones=1, per_clone=2)
    msg = str(exc_info.value)
    assert "mutate must be called after" in msg, (
        f"guard didn't report the FIRST offender; got: {msg!r}"
    )


# ──────────────────────────────────────────────────────────────────
# Contaminant deliberately NOT guarded — pin the decision
# ──────────────────────────────────────────────────────────────────


def test_contaminate_pre_fork_not_guarded_yet() -> None:
    """``contaminate(prob=...)`` is intentionally left out of the
    descendant-phase table for this slice. A follow-up may add it
    once its placement semantics are decided; pin the current
    state so the omission is deliberate rather than accidental."""
    exp = _base().contaminate(prob=0.5).expand_clones(n_clones=1, per_clone=2)
    # If contaminate gets classified later, this test should flip
    # to expect the guard. For now it must compile without raising.
    assert exp is not None
