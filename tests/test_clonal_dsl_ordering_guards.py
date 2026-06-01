"""Spec tests for the clonal DSL ordering guards (Bugs A / B / C).

The plan-split audit's earlier probe found three silent-drop bugs
where a DSL step landed on the wrong side of ``expand_clones()``:

- **Bug A** — ``invert_d(...)`` placed AFTER ``expand_clones(...)``
  silently drops the inversion (no recombine step in the post-fork
  plan to consume the probability; ``d_inverted=False`` for every
  descendant even at ``prob=1.0``).
- **Bug B** — ``receptor_revision(...)`` placed AFTER
  ``expand_clones(...)`` silently drops the revision (same
  mechanism; ``receptor_revision_applied=False`` and
  ``original_v_call=""`` for every descendant).
- **Bug C** — ``paired_end(...)`` placed BEFORE
  ``expand_clones(...)`` silently produces empty
  ``r1_sequence``/``r2_sequence`` on descendants (the paired-end
  pass lands at the end of the **pre-fork** plan; trace records
  sit on the parent; descendant AIRR projection finds nothing).

This slice adds DSL-boundary ordering guards in
:meth:`Experiment.invert_d`, :meth:`Experiment.receptor_revision`,
and :meth:`Experiment.expand_clones`. The guards reject the
misordered calls with a message that names the correct ordering
plus the biological reason.

Spec coverage (from the user brief):

1. ``.recombine().expand_clones(...).invert_d(...)`` raises.
2. ``.recombine().expand_clones(...).receptor_revision(...)``
   raises.
3. ``.recombine().paired_end(...).expand_clones(...)`` raises.
4. Valid canonical clonal pipeline passes (no exception, real
   biological effect at runtime).
5. Non-clonal ``.recombine().paired_end(...)`` still passes.
6. Error messages name the correct ordering AND the biological
   reason.
7. Contract / audit pins for the three bugs flip from
   "currently broken" to "guarded."

Out of scope here: no compiler-lowering changes, no validator
changes. The audit's "silently produced wrong output" failure mode
is closed at the DSL boundary; the lowering still does what it
did, but the misordered call sites can no longer reach it.
"""
from __future__ import annotations

import re

import pytest

import GenAIRR as ga


# ──────────────────────────────────────────────────────────────────
# Spec 1 — invert_d after expand_clones raises (Bug A guard)
# ──────────────────────────────────────────────────────────────────


def test_invert_d_after_expand_clones_raises() -> None:
    """``.recombine().expand_clones(...).invert_d(...)`` must raise
    a clear ordering error. Before this slice the call accepted
    silently and the inversion probability was dropped — every
    descendant projected ``d_inverted=False`` even at ``prob=1.0``."""
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .expand_clones(n_clones=1, per_clone=2)
    )
    with pytest.raises(ValueError) as exc_info:
        exp.invert_d(prob=1.0)
    msg = str(exc_info.value)
    # Must name BOTH the correct ordering and the biological reason.
    assert "invert_d must be called before expand_clones" in msg, msg
    assert "recombination" in msg.lower(), msg
    assert "inherited" in msg, msg


# ──────────────────────────────────────────────────────────────────
# Spec 2 — receptor_revision after expand_clones raises (Bug B guard)
# ──────────────────────────────────────────────────────────────────


def test_receptor_revision_after_expand_clones_raises() -> None:
    """``.recombine().expand_clones(...).receptor_revision(...)``
    must raise. Pre-slice: silently dropped — ``original_v_call=""``
    on every descendant even at ``prob=1.0``."""
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .expand_clones(n_clones=1, per_clone=2)
    )
    with pytest.raises(ValueError) as exc_info:
        exp.receptor_revision(prob=1.0)
    msg = str(exc_info.value)
    assert (
        "receptor_revision must be called before expand_clones" in msg
    ), msg
    assert "recombination" in msg.lower(), msg
    assert "inherited" in msg, msg


# ──────────────────────────────────────────────────────────────────
# Spec 3 — paired_end before expand_clones raises (Bug C guard)
# ──────────────────────────────────────────────────────────────────


def test_paired_end_before_expand_clones_raises_at_expand_clones() -> None:
    """``.recombine().paired_end(...).expand_clones(...)`` raises at
    the ``expand_clones`` call site (paired_end itself can't know
    whether a later ``expand_clones`` will be appended; the guard
    fires the moment the fork is appended on top of an existing
    paired-end step). Pre-slice: silently produced empty r1/r2
    on descendants."""
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .paired_end(r1_length=150, insert_size=300)
    )
    with pytest.raises(ValueError) as exc_info:
        exp.expand_clones(n_clones=1, per_clone=2)
    msg = str(exc_info.value)
    assert "paired_end must be called after expand_clones" in msg, msg
    assert "descendant-specific" in msg, msg


# ──────────────────────────────────────────────────────────────────
# Spec 4 — Valid canonical clonal pipeline passes
# ──────────────────────────────────────────────────────────────────


def test_canonical_clonal_pipeline_with_all_three_passes() -> None:
    """The canonical order
    ``recombine → invert_d → receptor_revision → expand_clones →
    paired_end`` builds cleanly and the runtime output is
    biologically correct (descendants carry the parent's
    ``d_inverted`` and ``original_v_call`` and have non-empty
    r1/r2 sequences)."""
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .invert_d(prob=1.0)
        .receptor_revision(prob=1.0)
        .expand_clones(n_clones=2, per_clone=3)
        .paired_end(r1_length=150, insert_size=300)
    )
    result = exp.run_records(seed=0)
    assert len(result) == 6
    for r in result:
        # invert_d effect inherited.
        assert r["d_inverted"] is True
        # receptor_revision provenance inherited (Bug D fix).
        assert r["receptor_revision_applied"] is True
        assert r["original_v_call"] != ""
        # paired_end sampled per descendant.
        assert len(r["r1_sequence"]) == 150
        assert len(r["r2_sequence"]) == 150


def test_canonical_clonal_pipeline_minimal() -> None:
    """The smallest valid clonal pipeline
    ``recombine → expand_clones`` continues to compile and run."""
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .expand_clones(n_clones=1, per_clone=2)
    )
    result = exp.run_records(seed=0)
    assert len(result) == 2


# ──────────────────────────────────────────────────────────────────
# Spec 5 — Non-clonal paired_end still passes
# ──────────────────────────────────────────────────────────────────


def test_non_clonal_paired_end_still_passes() -> None:
    """``paired_end()`` itself stays unguarded — non-clonal
    pipelines must continue to use it without an ``expand_clones``
    call. The guard lives on ``expand_clones``, not on
    ``paired_end``."""
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .paired_end(r1_length=150, insert_size=300)
    )
    result = exp.run_records(n=3, seed=0)
    assert len(result) == 3
    for r in result:
        assert len(r["r1_sequence"]) == 150
        assert len(r["r2_sequence"]) == 150


def test_non_clonal_invert_d_and_receptor_revision_unaffected() -> None:
    """``invert_d()`` and ``receptor_revision()`` still work in
    non-clonal pipelines — the guard checks for a clonal fork,
    which is absent here."""
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .invert_d(prob=1.0)
        .receptor_revision(prob=1.0)
    )
    result = exp.run_records(n=2, seed=0)
    assert len(result) == 2
    for r in result:
        assert r["d_inverted"] is True
        assert r["receptor_revision_applied"] is True
        assert r["original_v_call"] != ""


# ──────────────────────────────────────────────────────────────────
# Spec 6 — Error messages name the correct ordering + biological reason
# ──────────────────────────────────────────────────────────────────


def test_error_message_invert_d_names_recombination_time() -> None:
    """The error message for the Bug A guard must explain the
    biology (D inversion is recombination-time) so a confused
    builder can self-diagnose without reading the source. We pin
    a few load-bearing phrases."""
    exp = ga.Experiment.on("human_igh").recombine().expand_clones(
        n_clones=1, per_clone=1
    )
    with pytest.raises(ValueError) as exc_info:
        exp.invert_d(prob=0.5)
    msg = str(exc_info.value)
    # The correct ordering.
    assert "before expand_clones" in msg
    # The biological reason — pin the phrase that names the
    # mechanism.
    assert re.search(r"recombination[- ]?time", msg.lower())
    # The fix instruction.
    assert "Move" in msg or "move" in msg


def test_error_message_receptor_revision_names_recombination_time() -> None:
    """Bug B guard message mirrors Bug A's shape — names the
    correct ordering, the biological reason, and the fix."""
    exp = ga.Experiment.on("human_igh").recombine().expand_clones(
        n_clones=1, per_clone=1
    )
    with pytest.raises(ValueError) as exc_info:
        exp.receptor_revision(prob=0.5)
    msg = str(exc_info.value)
    assert "before expand_clones" in msg
    assert re.search(r"recombination[- ]?time", msg.lower())
    assert "Move" in msg or "move" in msg


def test_error_message_paired_end_names_observation_time() -> None:
    """Bug C guard message names the observation/readout-time
    nature of paired-end."""
    exp = ga.Experiment.on("human_igh").recombine().paired_end(
        r1_length=150, insert_size=300
    )
    with pytest.raises(ValueError) as exc_info:
        exp.expand_clones(n_clones=1, per_clone=1)
    msg = str(exc_info.value)
    assert "after expand_clones" in msg
    assert "descendant-specific" in msg
    assert "Move" in msg or "move" in msg


# ──────────────────────────────────────────────────────────────────
# Spec 7 — Guard pins (no more silent drops) — flipped from the audit
# ──────────────────────────────────────────────────────────────────


def test_pin_guarded_invert_d_post_fork_drop_cannot_reach_runtime() -> None:
    """Before this slice the audit's Bug A probe showed
    ``.expand_clones(...).invert_d(prob=1.0)`` producing
    ``d_inverted=[False, False]`` (silent drop). Now the misorder
    is rejected at the DSL boundary so the silent runtime
    behaviour is structurally unreachable. Pin the rejection
    rather than the broken runtime output."""
    with pytest.raises(ValueError, match="before expand_clones"):
        (
            ga.Experiment.on("human_igh")
            .recombine()
            .expand_clones(n_clones=1, per_clone=2)
            .invert_d(prob=1.0)
        )


def test_pin_guarded_receptor_revision_post_fork_drop_cannot_reach_runtime() -> None:
    """Bug B's silent drop closed at the DSL boundary."""
    with pytest.raises(ValueError, match="before expand_clones"):
        (
            ga.Experiment.on("human_igh")
            .recombine()
            .expand_clones(n_clones=1, per_clone=2)
            .receptor_revision(prob=1.0)
        )


def test_pin_guarded_paired_end_pre_fork_drop_cannot_reach_runtime() -> None:
    """Bug C's silent drop closed at the DSL boundary."""
    with pytest.raises(ValueError, match="after expand_clones"):
        (
            ga.Experiment.on("human_igh")
            .recombine()
            .paired_end(r1_length=150, insert_size=300)
            .expand_clones(n_clones=1, per_clone=2)
        )


# ──────────────────────────────────────────────────────────────────
# Spec 7b — Compiler lowering unchanged (no follow-on need)
# ──────────────────────────────────────────────────────────────────


def test_compiler_lowering_unchanged_recombine_only_still_works() -> None:
    """Pin that the guards live entirely at the DSL boundary — the
    compiler lowering is not touched. A trivial recombine-only
    pipeline must still compile and run identically."""
    exp = ga.Experiment.on("human_igh").recombine()
    result = exp.run_records(n=2, seed=0)
    assert len(result) == 2


def test_guard_error_fires_before_pipeline_compiles() -> None:
    """The guard is at the DSL boundary, not the compile step —
    the error must surface during the fluent builder call, before
    ``.compile()`` is reached. We never reach the runtime."""
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .expand_clones(n_clones=1, per_clone=2)
    )
    # The expand_clones call already returned an Experiment with a
    # fork step. compile() on its current state should still work.
    compiled = exp.compile()
    assert compiled is not None
    # The misordered call is what raises.
    with pytest.raises(ValueError):
        exp.invert_d(prob=1.0)


# ──────────────────────────────────────────────────────────────────
# Extra coverage — `_has_clonal_fork` helper behaves correctly
# ──────────────────────────────────────────────────────────────────


def test_has_clonal_fork_helper_reflects_pipeline_state() -> None:
    """The ``_has_clonal_fork`` helper underlies every guard;
    pin its behaviour directly so a refactor that changes its
    semantics surfaces here."""
    base = ga.Experiment.on("human_igh").recombine()
    assert base._has_clonal_fork() is False
    forked = base.expand_clones(n_clones=1, per_clone=1)
    assert forked._has_clonal_fork() is True
