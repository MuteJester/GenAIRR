"""Spec tests for the Slice 2 parent-outcome accessor.

Closes the audit's §12 Slice 2 (Parent-Outcome Read-Only Surface,
recommended next after Slice 1). The slice adds:

- :attr:`GenAIRR.result.SimulationResult.parents` —
  ``Optional[List[Outcome]]``. ``None`` for non-clonal results;
  one parent per clone for clonal results.
- ``record["parent_id"]`` — integer index into ``result.parents``,
  written onto every clonal AIRR record dict (= ``clone_id``
  today, kept separate by semantic).
- ``CompiledClonalExperiment.run_records`` retains the parent
  ``Outcome`` instead of dropping it after extracting
  ``final_simulation()``.

Spec coverage (from the user brief):

1. Clonal result exposes ``parents`` with one parent per clone.
2. Every clonal record carries ``parent_id``.
3. ``record["parent_id"]`` indexes a real parent outcome.
4. Parent outcome has the pre-branch trace; child trace remains
   descendant-only and is not inflated.
5. Non-clonal result keeps ``parents is None`` and ``parent_id``
   absent from records.
6. Existing ``validate_families()`` still passes and does not
   require parents.
7. Companion contract pins flip (covered in
   ``test_clonal_parent_contract.py``).

Out of scope here (per the spec): no parent-aware family
validator, no Rust-side ``Outcome.parent_id`` field, no
``ClonalFamily`` aggregate, no clonal trace-file format, no
change to descendant trace semantics.
"""
from __future__ import annotations

import pytest

import GenAIRR as ga
from GenAIRR import _engine as ge
from GenAIRR.result import SimulationResult


# ──────────────────────────────────────────────────────────────────
# Shared fixtures
# ──────────────────────────────────────────────────────────────────


def _make_clonal(
    *, n_clones: int = 2, per_clone: int = 3, with_mutate: bool = False
) -> ga.Experiment:
    exp = ga.Experiment.on("human_igh").recombine().expand_clones(
        n_clones=n_clones, per_clone=per_clone
    )
    if with_mutate:
        exp = exp.mutate(count=5)
    return exp


# ──────────────────────────────────────────────────────────────────
# Spec 1 — Clonal result exposes one parent per clone
# ──────────────────────────────────────────────────────────────────


def test_clonal_result_exposes_one_parent_per_clone() -> None:
    """``result.parents`` is a list whose length equals the
    pipeline's ``n_clones`` value. Each entry is a real ``Outcome``
    object exposing the standard accessors."""
    result = _make_clonal(n_clones=4, per_clone=3, with_mutate=True).run_records(
        seed=0
    )
    parents = result.parents
    assert parents is not None
    assert isinstance(parents, list)
    assert len(parents) == 4
    # Every parent is an Outcome — quack-typed via the accessors
    # the audit relies on (final_simulation / trace / events).
    for p in parents:
        assert hasattr(p, "final_simulation")
        assert hasattr(p, "trace")
        assert hasattr(p, "event_count")
        assert hasattr(p, "events")
        assert hasattr(p, "pass_names")
        # The parent's final IR is non-empty — recombine produced
        # something.
        sim = p.final_simulation()
        assert len(sim) > 0
        assert sim.v_allele_id() is not None
        assert sim.j_allele_id() is not None


def test_clonal_parents_length_independent_of_per_clone() -> None:
    """``len(.parents) == n_clones`` regardless of ``per_clone``.
    ``len(.outcomes)`` continues to equal ``n_clones * per_clone``."""
    cases = [(2, 3), (5, 1), (1, 10), (3, 5)]
    for n_clones, per_clone in cases:
        result = _make_clonal(
            n_clones=n_clones, per_clone=per_clone
        ).run_records(seed=0)
        assert result.parents is not None
        assert len(result.parents) == n_clones, (
            f"n_clones={n_clones}, per_clone={per_clone}: "
            f"got {len(result.parents)} parents"
        )
        assert result.outcomes is not None
        assert len(result.outcomes) == n_clones * per_clone


# ──────────────────────────────────────────────────────────────────
# Spec 2 — Every clonal record carries `parent_id`
# ──────────────────────────────────────────────────────────────────


def test_every_clonal_record_carries_parent_id() -> None:
    """``record["parent_id"]`` is present on every record produced
    by a clonal pipeline, alongside the existing ``clone_id``."""
    result = _make_clonal(n_clones=3, per_clone=4).run_records(seed=0)
    for r in result:
        assert "clone_id" in r
        assert "parent_id" in r, (
            f"record missing parent_id: keys={sorted(r.keys())[:5]}…"
        )
        assert isinstance(r["parent_id"], int)


def test_parent_id_matches_clone_id_today() -> None:
    """Slice 2 makes ``parent_id == clone_id`` numerically because
    clones are dense and zero-based. The semantics are distinct
    (clone_id = family identity; parent_id = index into
    ``result.parents``) so they are stamped separately."""
    result = _make_clonal(n_clones=4, per_clone=2).run_records(seed=0)
    for r in result:
        assert r["parent_id"] == r["clone_id"]


# ──────────────────────────────────────────────────────────────────
# Spec 3 — `record["parent_id"]` indexes a real parent
# ──────────────────────────────────────────────────────────────────


def test_record_parent_id_indexes_a_real_parent_outcome() -> None:
    """For every record, ``result.parents[record["parent_id"]]`` is
    a real ``Outcome`` whose ``final_simulation()`` matches the
    descendant's initial IR (``revision(0)``) on the shared
    recombination state (V/D/J allele ids, pool length, regions).
    This is the audit's §2 boundary-handoff invariant — Slice 2
    makes it observable; this test proves it through the new API."""
    result = _make_clonal(n_clones=2, per_clone=3, with_mutate=True).run_records(
        seed=0
    )
    assert result.outcomes is not None
    assert result.parents is not None
    for i, rec in enumerate(result):
        parent = result.parents[rec["parent_id"]]
        descendant = result.outcomes[i]
        parent_sim = parent.final_simulation()
        desc_initial = descendant.revision(0)
        assert desc_initial.v_allele_id() == parent_sim.v_allele_id()
        assert desc_initial.d_allele_id() == parent_sim.d_allele_id()
        assert desc_initial.j_allele_id() == parent_sim.j_allele_id()
        assert len(desc_initial) == len(parent_sim)
        assert desc_initial.region_count() == parent_sim.region_count()


# ──────────────────────────────────────────────────────────────────
# Spec 4 — Parent has pre-branch trace; descendant trace not inflated
# ──────────────────────────────────────────────────────────────────


def test_parent_outcome_carries_pre_branch_trace() -> None:
    """The parent's ``trace()`` and ``events()`` reflect the pre-
    fork plan (typically just ``recombine``): the parent has a
    non-empty trace and a non-empty event ledger. Slice 2's value
    is making this observable; without ``parents`` it was
    structurally inaccessible."""
    result = _make_clonal(n_clones=2, per_clone=3).run_records(seed=0)
    assert result.parents is not None
    for parent in result.parents:
        # Pre-fork plan has recombine — trace + events both non-
        # empty.
        assert len(parent.trace()) > 0, (
            "parent trace is empty; the pre-fork plan didn't sample "
            "anything — either the parent retention is wrong or the "
            "fork-split partitioning regressed."
        )
        assert parent.event_count() > 0
        # The parent's pass_names include the recombine machinery.
        # We check substring rather than exact names so the test
        # survives minor pass-ordering refactors.
        names = " ".join(parent.pass_names())
        for substr in ("sample_allele", "assemble"):
            assert substr in names, (
                f"parent pass_names lacks {substr!r}; pre-fork plan "
                f"shape drifted. Pass names: {parent.pass_names()}"
            )


def test_descendant_trace_remains_descendant_only_not_inflated() -> None:
    """Slice 2 must NOT lift parent trace events onto descendants.
    Each descendant's ``trace()`` and ``event_count()`` continue to
    reflect only post-fork passes — strictly less than a non-clonal
    ``recombine + same post-fork plan`` run from a fresh IR.

    This is the audit's §11 backwards-compat guarantee: descendant
    trace semantics unchanged."""
    clonal_result = _make_clonal(
        n_clones=1, per_clone=2, with_mutate=True
    ).run_records(seed=0)
    assert clonal_result.outcomes is not None
    desc_event_counts = [o.event_count() for o in clonal_result.outcomes]
    desc_trace_lens = [len(o.trace()) for o in clonal_result.outcomes]

    non_clonal = (
        ga.Experiment.on("human_igh").recombine().mutate(count=5).compile()
    )
    non_clonal_outcome = non_clonal.run(n=1, seed=0)[0]
    non_clonal_event_count = non_clonal_outcome.event_count()
    non_clonal_trace_len = len(non_clonal_outcome.trace())

    for de in desc_event_counts:
        assert de < non_clonal_event_count, (
            f"descendant event_count ({de}) >= non-clonal "
            f"({non_clonal_event_count}); parent events leaked onto "
            "descendants — Slice 2 contract violated."
        )
    for dt in desc_trace_lens:
        assert dt < non_clonal_trace_len, (
            f"descendant trace_len ({dt}) >= non-clonal "
            f"({non_clonal_trace_len}); parent trace leaked onto "
            "descendants."
        )


def test_parent_trace_not_duplicated_per_descendant() -> None:
    """One parent ``Outcome`` per clone — not one per descendant.
    With ``per_clone=20`` and ``n_clones=2``, ``len(.parents) == 2``
    (not ``40``). Confirms the audit's §4 "Candidate B"
    no-duplication property: parents are stored once, not per-
    descendant."""
    result = _make_clonal(n_clones=2, per_clone=20).run_records(seed=0)
    assert result.parents is not None
    assert len(result.parents) == 2
    # And the references are stable — looking up the same
    # parent_id twice returns the same object.
    p0a = result.parents[result.records[0]["parent_id"]]
    p0b = result.parents[result.records[1]["parent_id"]]
    assert p0a is p0b, (
        "parent retrieval not idempotent — same parent_id returns "
        "different Outcome objects."
    )


# ──────────────────────────────────────────────────────────────────
# Spec 5 — Non-clonal result keeps `parents is None` and no parent_id
# ──────────────────────────────────────────────────────────────────


def test_non_clonal_result_keeps_parents_none_and_no_parent_id() -> None:
    """Non-clonal pipelines (no ``expand_clones``) return a result
    whose ``.parents is None``. Records have no ``parent_id`` /
    ``clone_id`` columns — symmetric with how ``clone_id`` was
    already omitted in Slice 0."""
    result = ga.Experiment.on("human_igh").recombine().run_records(n=5, seed=0)
    assert result.parents is None
    for r in result:
        assert "parent_id" not in r
        assert "clone_id" not in r


def test_simulationresult_built_from_records_has_no_parents() -> None:
    """A ``SimulationResult`` constructed directly from records
    (no ``outcomes``, no ``parents`` — the records-only constructor
    consumers use after loading from TSV) reports
    ``parents is None``. Slice 2 must not regress the records-only
    constructor."""
    records = [
        {"sequence_id": "s0", "sequence": "ACGT"},
        {"sequence_id": "s1", "sequence": "CGAT"},
    ]
    result = SimulationResult(records)
    assert result.parents is None
    assert result.outcomes is None


# ──────────────────────────────────────────────────────────────────
# Spec 6 — `validate_families()` still works and doesn't require parents
# ──────────────────────────────────────────────────────────────────


def test_validate_families_still_works_without_parents() -> None:
    """The Slice 1 field-only ``validate_families`` continues to
    pass on a clean clonal run and continues to work when
    ``parents is None`` (records-only result). Slice 2 must not
    silently introduce a parent dependency into the field-only
    path."""
    # Clean clonal run with parents available — still ok.
    clonal_result = _make_clonal(
        n_clones=3, per_clone=4, with_mutate=True
    ).run_records(seed=0, expose_provenance=True)
    report = clonal_result.validate_families()
    assert report.ok
    assert report.family_count == 3

    # Records-only result (parents=None) — still ok.
    by_hand_records = [
        {"clone_id": 0, "parent_id": 0, "truth_v_call": "V", "truth_d_call": "D", "truth_j_call": "J"},
        {"clone_id": 0, "parent_id": 0, "truth_v_call": "V", "truth_d_call": "D", "truth_j_call": "J"},
    ]
    by_hand_result = SimulationResult(by_hand_records)
    by_hand_report = by_hand_result.validate_families()
    assert by_hand_report.ok
    assert by_hand_result.parents is None  # not needed


def test_validate_records_true_clonal_still_runs_both_gates() -> None:
    """``validate_records=True`` on a clonal ``run_records`` still
    works after Slice 2 — per-record gate + family gate run in
    order, no new exception type, returns the same
    ``SimulationResult``. Slice 2 is additive on the data side
    only."""
    result = _make_clonal(
        n_clones=2, per_clone=3, with_mutate=True
    ).run_records(seed=0, validate_records=True)
    assert len(result) == 6
    assert result.parents is not None
    assert len(result.parents) == 2


# ──────────────────────────────────────────────────────────────────
# Spec 7 — Parent's IR matches descendant's revision(0) end-to-end
# ──────────────────────────────────────────────────────────────────


def test_parent_observability_unlocks_full_recombination_history() -> None:
    """The parent ``Outcome`` exposes the full ``revisions`` history
    of the pre-fork plan: ``revision(0)`` is the initial empty IR
    and ``final_simulation() == revision(revision_count - 1)``
    is the post-recombination IR. This is the data Slice 3+'s
    parent-aware validator will consume for pre-SHM junction,
    mutation-distance distribution, and truth-allele provenance."""
    result = _make_clonal(n_clones=1, per_clone=1).run_records(seed=0)
    assert result.parents is not None
    parent = result.parents[0]
    # At least one revision per pre-fork pass + the initial.
    assert parent.revision_count() >= 2
    # The first revision is the empty starting IR; len(sim) == 0
    # because nothing has been pooled yet.
    initial = parent.revision(0)
    assert len(initial) == 0
    # The final revision matches final_simulation().
    final_by_index = parent.revision(parent.revision_count() - 1)
    final_by_method = parent.final_simulation()
    assert len(final_by_index) == len(final_by_method)
    assert final_by_index.v_allele_id() == final_by_method.v_allele_id()
