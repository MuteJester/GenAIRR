"""Contract pins for the clonal parent-outcome audit.

Companion to
[`docs/clonal_parent_outcome_design.md`](../docs/clonal_parent_outcome_design.md).
The audit is pre-implementation: it freezes today's parent-outcome
lifecycle (``pin_scaffold_*``) and the addressability gaps the
audit identifies (``pin_absence_*``). No Slice 2 implementation
is proposed here; this file is the baseline a future Slice 2
(parent-outcome read-only surface) will land against.

Sibling to
[`tests/test_clonal_family_contract.py`](test_clonal_family_contract.py) —
that file pins the higher-level family / lineage architecture; this
file zooms in on the **parent-outcome lifecycle** specifically.
Some pins necessarily overlap (e.g. "no parent_id field"); we
re-pin here from the parent-lifecycle angle so the doc ↔ test
lockstep stays unambiguous when Slice 2 lands.

The split:

- ``pin_scaffold_*`` tests freeze the surfaces that exist today:
  the ``Outcome`` / ``Simulation`` shapes the boundary handoff
  uses, the Rust ``run_one_from_with_policy`` entry point, the
  fresh-trace allocation inside ``execute_transactional``, the
  Python orchestration loop's ``pre.run`` → ``final_simulation``
  → ``post.run_from`` shape, and the IR-equality contract at the
  boundary (descendant's ``revision(0)`` agrees with parent's
  ``final_simulation()`` on the shared fields).
- ``pin_absence_*`` tests freeze the remaining gaps Slice 3+
  closes: no ``parent_id`` on the Rust ``Outcome`` type, no
  ``run_family`` / ``parent_outcome_for`` direct helpers (the
  parent is reachable via ``result.parents[record["parent_id"]]``
  after Slice 2), no parent-aware validator, no clonal-aware
  trace-file format, no per-clone SHM event aggregation.

Slice-2 + Slice-3 flip history:

Slice 2 (parent observability) lockstep pins:
- ``pin_present_parents_attr_on_clonal_result`` (was absence)
- ``pin_present_parent_id_on_descendant_records`` (was absence)
- ``pin_absence_no_run_family_or_parent_outcome_helpers``
  (partially flipped from
  ``pin_absence_parent_outcome_unreachable_from_descendant``;
  ``.parents`` is now present, but the direct ``run_family``
  helper still isn't)
- ``pin_scaffold_simulationresult_slots_documented_for_extension``
  now includes ``_parents`` in the slot list.

Slice 3 (parent-aware validator) lockstep pins:
- ``pin_present_parent_aware_validator`` (was absence). Validator
  shipped; deliberately NOT auto-invoked by
  ``validate_records=True`` — that decision is pinned in the
  flipped test too.

Closing condition for the audit-only phase: every test here
passes on the current codebase. A failure in
``pin_scaffold_*`` means the parent-lifecycle infrastructure
drifted; a failure in ``pin_absence_*`` means someone shipped
part of Slice 2 without flipping the companion pin.
"""
from __future__ import annotations

import inspect
import re
from pathlib import Path
from typing import List

import pytest

import GenAIRR as ga
from GenAIRR import _compiled
from GenAIRR.result import SimulationResult


# ──────────────────────────────────────────────────────────────────
# Shared fixtures
# ──────────────────────────────────────────────────────────────────


_REPO_ROOT = Path(__file__).resolve().parent.parent


def _read_rust(relative: str) -> str:
    return (_REPO_ROOT / "engine_rs" / "src" / relative).read_text(
        encoding="utf-8"
    )


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
# 1. Scaffold — `Outcome` struct shape today
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_outcome_struct_shape_today() -> None:
    """``Outcome`` carries exactly the four Rust-side fields:
    ``revisions``, ``pass_names``, ``trace``, ``events``. The PyO3
    bindings expose them as ``revision_count`` / ``revision(i)`` /
    ``final_simulation`` / ``pass_names`` / ``revision_after`` /
    ``trace`` / ``event_count`` / ``events``. No ``parent_id``,
    ``parent_trace``, or ``parent_outcome`` accessors exist today.
    Slice 2 / Slice 3 add the back-pointer; this pin flips then."""
    exp = ga.Experiment.on("human_igh").recombine().mutate(count=3)
    outcome = exp.compile().run(n=1, seed=0)[0]

    # Present accessors (today's contract).
    for present in (
        "revision_count",
        "revision",
        "final_simulation",
        "pass_names",
        "revision_after",
        "trace",
        "event_count",
        "events",
    ):
        assert hasattr(outcome, present), (
            f"PyOutcome.{present} missing; parent-outcome audit relies on it"
        )

    # Absent accessors (the audit's gap pins).
    for forbidden in (
        "parent_id",
        "parent_trace",
        "parent_outcome",
        "parent_simulation",
        "parent_revisions",
        "parent_events",
    ):
        assert not hasattr(outcome, forbidden), (
            f"PyOutcome.{forbidden} now exists; Slice 2 or Slice 3 has "
            "landed — flip the companion `pin_absence_*` pin in lockstep."
        )


# ──────────────────────────────────────────────────────────────────
# 2. Scaffold — `Simulation` struct shape today
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_simulation_struct_shape_today() -> None:
    """The ``Simulation`` IR — the *carrier* across the parent →
    descendant boundary — has six fields: ``pool``, ``sequence``,
    ``assignments``, ``segment_calls``, ``dirty_log``,
    ``mutation_count``. The audit's §2 "what is shared" table maps
    every entry; a new field would silently change the boundary
    handoff and must be reviewed against that table.

    Trace-level invariant: the source struct does not carry
    ``trace`` or ``events``; those live on ``Outcome``. That
    separation is what makes the parent trace structurally
    inaccessible from the descendant."""
    src = _read_rust("ir/simulation.rs")
    # Lift the struct body (between `pub struct Simulation {` and
    # the next `}` at column 0). A loose check on field names is
    # sufficient for the audit-doc lockstep.
    body_match = re.search(
        r"pub struct Simulation \{(.*?)\n\}", src, re.DOTALL
    )
    assert body_match, "could not locate Simulation struct body"
    body = body_match.group(1)
    for field in (
        "pub pool",
        "pub sequence",
        "pub assignments",
        "pub segment_calls",
        "pub dirty_log",
        "pub mutation_count",
    ):
        assert field in body, (
            f"Simulation no longer carries `{field}`; audit §2 shared-"
            "state table is now stale."
        )
    # No trace / events on Simulation — that's the boundary
    # invariant the parent-outcome audit hinges on.
    for forbidden in ("pub trace:", "pub events:"):
        assert forbidden not in body, (
            f"Simulation now carries `{forbidden}`; the parent trace is "
            "no longer structurally inaccessible from descendants — "
            "audit §1.4 conclusion has drifted."
        )


# ──────────────────────────────────────────────────────────────────
# 3. Scaffold — Rust entry point shape
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_run_one_from_with_policy_takes_initial_sim() -> None:
    """``CompiledSimulator::run_one_from_with_policy(initial:
    Simulation, seed, policy) -> Result<Outcome, PassError>`` is the
    entry point through which the parent IR crosses into a
    descendant run. Its signature is the load-bearing contract for
    the Python orchestration's ``post.run_from(parent_sim,
    desc_seed)`` call; a refactor that renamed or changed the
    initial-sim parameter would break the audit's §1.2 description."""
    src = _read_rust("compiled/mod.rs")
    # The signature spans multiple lines. Look for the function
    # header and the `initial: Simulation` parameter.
    sig = re.search(
        r"pub fn run_one_from_with_policy\(\s*"
        r"&self,\s*initial:\s*Simulation,\s*seed:\s*u64,\s*"
        r"policy:\s*ExecutionPolicy",
        src,
    )
    assert sig is not None, (
        "CompiledSimulator::run_one_from_with_policy signature "
        "changed; audit §1.2 entry-point reference drifted."
    )


# ──────────────────────────────────────────────────────────────────
# 4. Scaffold — `execute_transactional` allocates fresh trace
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_execute_transactional_allocates_fresh_trace() -> None:
    """At the top of ``execute_transactional`` the trace is
    allocated fresh: ``let mut trace = Trace::new()``. This is what
    makes the descendant's trace structurally independent of the
    parent's — even though the descendant runs from the parent's
    IR, the choice-trace is reset. The audit's §1.2 conclusion
    ("descendant gets brand-new Trace") hinges on this line."""
    src = _read_rust("compiled/execute.rs")
    assert re.search(
        r"let\s+mut\s+trace\s*:\s*Trace?\s*=\s*Trace::new\(\)|"
        r"let\s+mut\s+trace\s*=\s*Trace::new\(\)",
        src,
    ), (
        "execute_transactional no longer allocates a fresh Trace at "
        "the top; audit §1.2 / §4 conclusion has drifted — the "
        "parent's trace may now be reachable from the descendant, "
        "which Slice 2 would need to account for explicitly."
    )


# ──────────────────────────────────────────────────────────────────
# 5. Scaffold — Python orchestration loop shape
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_orchestration_uses_pre_run_then_post_run_from() -> None:
    """The orchestration loop's ``pre.run`` → ``parent.final_simulation()``
    → ``post.run_from(parent_sim, ...)`` pattern is the audit's §1.1
    canonical shape. A refactor that batched parents, shared parents
    across clone indices, or replaced ``run_from`` with a fresh
    ``run`` would silently break the "one parent per clone" claim
    and break every descendant's truth-allele invariance."""
    src = inspect.getsource(_compiled.CompiledClonalExperiment.run_records)
    assert "self._pre.run(" in src
    assert ".final_simulation()" in src
    assert "self._post.run_from(" in src
    # And the loop nesting is parent-outer / descendant-inner.
    # A swap (descendants outer, parents inner) would re-run the
    # parent K times per clone — silently quadrupling work and
    # breaking determinism.
    pre_idx = src.index("self._pre.run(")
    desc_idx = src.index("self._post.run_from(")
    assert pre_idx < desc_idx, (
        "orchestration loop reordered: descendants are now built "
        "before parents. Audit §1.1 canonical shape drifted."
    )


# ──────────────────────────────────────────────────────────────────
# 6. Scaffold — descendant's `revision(0)` agrees with parent on
#    shared recombination state
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_simulation_carries_shared_recombination_state() -> None:
    """At the parent → descendant boundary the descendant's
    ``revision(0)`` (its initial IR) equals the parent's
    ``final_simulation()`` on the shared-state fields the audit's
    §2 table names: V/D/J allele ids, pool length, region count.
    This is the load-bearing invariant Slice 2+ family validators
    will rely on.

    We drive this manually rather than going through
    ``run_records`` so we can compare ``parent`` and the first
    descendant directly."""
    clonal = _make_clonal(n_clones=1, per_clone=2).compile()
    clone_seed = 0
    parent = clonal._pre.run(seed=clone_seed)
    parent_sim = parent.final_simulation()
    desc = clonal._post.run_from(parent_sim, clone_seed + 1)

    desc_initial = desc.revision(0)
    # Truth alleles: descendant's initial revision must agree with
    # the parent's final IR on the V/D/J truth ids.
    assert desc_initial.v_allele_id() == parent_sim.v_allele_id()
    assert desc_initial.d_allele_id() == parent_sim.d_allele_id()
    assert desc_initial.j_allele_id() == parent_sim.j_allele_id()
    # Pool length: the parent's nucleotide arena is fully carried.
    assert len(desc_initial) == len(parent_sim)
    # Region count: V / NP1 / D / NP2 / J regions present and equal.
    assert desc_initial.region_count() == parent_sim.region_count()


# ──────────────────────────────────────────────────────────────────
# 7. Scaffold — descendant sequences differ under post-fork passes
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_descendant_seqs_differ_under_post_fork_passes() -> None:
    """Audit §3: post-fork stochasticity (SHM, PCR, indels)
    produces per-descendant divergence. Re-pinned here for §3
    traceability; the existing ``test_g5_*`` invariants in
    ``test_experiment.py`` carry the broader behaviour. If this
    breaks, the post-fork stage isn't running stochastically and
    the parent-aware validator's mutation-distance check would
    spuriously pass on degenerate inputs."""
    result = _make_clonal(
        n_clones=1, per_clone=4, with_mutate=True
    ).run_records(seed=0)
    seqs = {r["sequence"] for r in result}
    assert len(seqs) > 1, (
        "single clone with mutate(count=5) produced identical "
        "descendants — post-fork stochasticity regressed."
    )


# ──────────────────────────────────────────────────────────────────
# 8. Scaffold — truth-call stability inside a clone under normal
#    fixtures
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_clonal_truth_calls_stable_within_clone_under_normal_fixtures() -> None:
    """Under the standard IGH fixture (no extreme SHM), every
    descendant of a clone has byte-identical ``truth_v_call`` /
    ``truth_d_call`` / ``truth_j_call``. This is the **audit
    baseline** for Slice 3's parent-aware validator: a clean batch
    today produces no truth divergences, so any future
    ``validate_families_with_parents`` failure is a real bug, not
    fixture noise."""
    result = _make_clonal(n_clones=3, per_clone=4, with_mutate=True).run_records(
        seed=0, expose_provenance=True
    )
    by_clone: dict = {}
    for r in result:
        by_clone.setdefault(r["clone_id"], []).append(r)
    for cid, recs in by_clone.items():
        truths = {(r["truth_v_call"], r["truth_d_call"], r["truth_j_call"]) for r in recs}
        assert len(truths) == 1, (
            f"clone {cid} truth tuple diverges across descendants "
            f"under a normal fixture: {truths}. Audit baseline broken."
        )


# ──────────────────────────────────────────────────────────────────
# 9. Scaffold — `SimulationResult` slots documented for extension
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_simulationresult_slots_documented_for_extension() -> None:
    """``SimulationResult.__slots__`` is the documented attribute
    surface. After Slice 2 it carries
    ``("_records", "_outcomes", "_parents")``. A reviewer adding
    a new slot must update this pin in lockstep so the slot-list
    change shows up as a deliberate, audited diff."""
    assert SimulationResult.__slots__ == ("_records", "_outcomes", "_parents"), (
        f"SimulationResult.__slots__ drifted to {SimulationResult.__slots__}; "
        "expected ('_records', '_outcomes', '_parents'). Either Slice 2 "
        "regressed (parent accessor removed) or a new slot landed without "
        "updating the lockstep pin."
    )


# ──────────────────────────────────────────────────────────────────
# 10. Scaffold — `validate_families` today is field-only
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_validate_families_today_is_field_only() -> None:
    """Slice 1's ``validate_families`` is purely dict-based: it
    iterates record dicts, groups by ``clone_id``, compares truth
    fields. It does **not** reference ``.parents`` /
    ``parent_id`` / parent outcomes — by audit §6 design.
    Slice 3's ``validate_families_with_parents`` is the parent-
    aware counterpart; ``validate_families`` stays as the
    field-only fallback. If a contributor folds parent reads into
    the field-only path, this pin surfaces it."""
    src = inspect.getsource(SimulationResult.validate_families)
    for parent_ref in ("self.parents", ".parents", "parent_id"):
        assert parent_ref not in src, (
            f"validate_families now references `{parent_ref}`; the "
            "field-only / parent-aware split has been lost. Slice 3 "
            "should add a sibling `validate_families_with_parents` "
            "instead of modifying this one."
        )


# ──────────────────────────────────────────────────────────────────
# 11. Absence — no `.parents` attribute on SimulationResult
# ──────────────────────────────────────────────────────────────────


def test_pin_present_parents_attr_on_clonal_result() -> None:
    """Slice 2 shipped: ``SimulationResult.parents`` is
    ``Optional[List[Outcome]]`` — ``None`` for non-clonal,
    one ``Outcome`` per clone for clonal. Flipped from
    ``pin_absence_no_parents_attr_on_simulation_result`` when the
    parent-outcome read-only surface landed.

    See [`tests/test_parent_outcomes.py`](test_parent_outcomes.py)
    for the spec-driven behavioural coverage; this pin is the
    audit-doc lockstep counterpart."""
    assert hasattr(SimulationResult, "parents"), (
        "SimulationResult.parents regressed; Slice 2 backed out."
    )
    # And the property is backed by the ``_parents`` slot.
    assert "_parents" in SimulationResult.__slots__
    # Clonal result populates it; non-clonal leaves it None.
    clonal = _make_clonal(n_clones=3, per_clone=2).run_records(seed=0)
    assert clonal.parents is not None
    assert len(clonal.parents) == 3
    non_clonal = ga.Experiment.on("human_igh").recombine().run_records(
        n=2, seed=0
    )
    assert non_clonal.parents is None


# ──────────────────────────────────────────────────────────────────
# 12. Absence — no `parent_id` field on `Outcome`
# ──────────────────────────────────────────────────────────────────


def test_pin_absence_no_parent_id_field_on_outcome() -> None:
    """``Outcome`` (Python / Rust) has no ``parent_id`` field. The
    audit's §11 recommends keeping this *deferred to Slice 3+*; if
    a contributor adds it in Slice 2 they need to flip this pin
    explicitly.

    Catching the field at the Python level catches both:
    (a) a Rust struct addition that gets exposed via PyO3, and
    (b) a Python-only patch on ``PyOutcome``."""
    result = _make_clonal().run_records(seed=0)
    assert result.outcomes is not None
    sample = result.outcomes[0]
    assert not hasattr(sample, "parent_id"), (
        "Outcome.parent_id now exists; flip pin in lockstep with the "
        "slice that added it."
    )


# ──────────────────────────────────────────────────────────────────
# 13. Absence — no `parent_id` on AIRR records
# ──────────────────────────────────────────────────────────────────


def test_pin_present_parent_id_on_descendant_records() -> None:
    """Slice 2 shipped: every clonal AIRR record carries
    ``parent_id`` alongside ``clone_id``. Today ``parent_id ==
    clone_id`` numerically because clones are dense and zero-
    based, but the semantics are distinct: ``clone_id`` is the
    family identity (Slice 0), ``parent_id`` is the index into
    ``result.parents`` (Slice 2).

    Flipped from ``pin_absence_no_parent_id_on_airr_records``."""
    result = _make_clonal(n_clones=2, per_clone=2).run_records(seed=0)
    for r in result:
        assert "clone_id" in r
        assert "parent_id" in r, (
            "AIRR record stopped carrying `parent_id`; Slice 2 "
            "regressed."
        )
        assert r["parent_id"] == r["clone_id"]
    # Non-clonal records still carry neither field — Slice 2's
    # parent_id is clonal-only, matching the clone_id pattern.
    non_clonal = ga.Experiment.on("human_igh").recombine().run_records(
        n=2, seed=0
    )
    for r in non_clonal:
        assert "parent_id" not in r
        assert "clone_id" not in r


# ──────────────────────────────────────────────────────────────────
# 14. Absence — parent outcome unreachable from descendant
# ──────────────────────────────────────────────────────────────────


def test_pin_absence_no_run_family_or_parent_outcome_helpers() -> None:
    """After Slice 2 the parent is reachable via
    ``result.parents[record["parent_id"]]``. The audit's §7 still
    leaves out a direct ``CompiledClonalExperiment.run_family(seed,
    clone_idx) -> (parent, descendants)`` helper and a
    ``Outcome``-to-``Outcome`` back-pointer; those are Slice 3+
    surfaces. Keep the absence pinned so a future contributor adding
    them does it deliberately.

    Flipped (partially) from
    ``pin_absence_parent_outcome_unreachable_from_descendant`` when
    Slice 2 made ``result.parents`` available."""
    clonal = _make_clonal().compile()
    for forbidden in (
        "parent_outcome_for",
        "parents_for",
        "run_family",
        "family_outcome",
    ):
        assert not hasattr(clonal, forbidden), (
            f"CompiledClonalExperiment.{forbidden} now exists; Slice 3+ "
            "has landed — flip pin in lockstep with the slice."
        )
    # And ``Outcome`` itself doesn't carry a back-pointer (Slice 3+
    # decision per audit §11). The Python-level ``parent_id`` lives
    # on the AIRR record dict, not on the ``Outcome`` object.
    result = clonal.run_records(seed=0)
    assert result.outcomes is not None
    sample = result.outcomes[0]
    for forbidden in ("parent_id", "parent_outcome", "parent_trace"):
        assert not hasattr(sample, forbidden), (
            f"Outcome.{forbidden} now exists; Slice 3+ has landed — "
            "flip pin."
        )


# ──────────────────────────────────────────────────────────────────
# 15. Absence — parent trace + events dropped at fork
# ──────────────────────────────────────────────────────────────────


def test_pin_absence_parent_trace_and_events_dropped_at_fork() -> None:
    """The descendant's ``outcome.event_count()`` and
    ``len(outcome.trace())`` reflect only post-fork passes; the
    parent's events / trace are dropped.

    Empirical pin: a clonal pipeline ``recombine + mutate(count=5)``
    produces descendants whose trace + event counts are STRICTLY
    LESS than a non-clonal ``recombine + mutate(count=5)`` run from
    a fresh IR. The difference equals the recombine-pass events the
    parent ran but didn't pass on.

    Slice 2 doesn't change this — the parent is retained but the
    descendant outcome stays unchanged. A pin failure here means
    someone lifted parent events onto descendant outcomes, which is
    a structural change Slice 2 explicitly avoids."""
    # Clonal: recombine pre-fork, mutate post-fork.
    clonal_result = _make_clonal(
        n_clones=1, per_clone=2, with_mutate=True
    ).run_records(seed=0)
    clonal_outcomes = clonal_result.outcomes
    assert clonal_outcomes is not None
    desc_event_counts: List[int] = [o.event_count() for o in clonal_outcomes]
    desc_trace_lens: List[int] = [len(o.trace()) for o in clonal_outcomes]

    # Non-clonal: same overall passes, from a fresh IR.
    non_clonal = (
        ga.Experiment.on("human_igh").recombine().mutate(count=5).compile()
    )
    non_clonal_outcome = non_clonal.run(n=1, seed=0)[0]
    non_clonal_event_count = non_clonal_outcome.event_count()
    non_clonal_trace_len = len(non_clonal_outcome.trace())

    # The non-clonal run carries BOTH recombine + mutate
    # contributions; the clonal descendant carries ONLY the post-
    # fork mutate contributions.
    for desc_events in desc_event_counts:
        assert desc_events < non_clonal_event_count, (
            f"clonal descendant event_count ({desc_events}) >= "
            f"non-clonal ({non_clonal_event_count}); parent events "
            "may have leaked onto descendant outcomes. If intentional, "
            "flip this audit pin (§4 invariant changed)."
        )
    for desc_trace in desc_trace_lens:
        assert desc_trace < non_clonal_trace_len, (
            f"clonal descendant trace_len ({desc_trace}) >= "
            f"non-clonal ({non_clonal_trace_len}); parent trace may "
            "have leaked onto descendants."
        )


# ──────────────────────────────────────────────────────────────────
# 16. Absence — no parent-aware validator
# ──────────────────────────────────────────────────────────────────


def test_pin_present_parent_aware_validator() -> None:
    """Slice 3 shipped: ``SimulationResult.validate_families_with_parents``
    is the parent-aware companion of the field-only
    ``validate_families``. Flipped from
    ``pin_absence_no_parent_aware_validator`` when the
    parent-aware validator landed.

    Spec-driven behavioural coverage lives in
    [`tests/test_validate_families_with_parents.py`](test_validate_families_with_parents.py);
    this pin is the audit-doc lockstep counterpart.

    The validator is **not wired into ``validate_records=True``**
    by Slice 3 — that gate stays at per-record + field-only
    family. Pin that decision here too so a future contributor
    who wires it surfaces the change deliberately."""
    assert hasattr(SimulationResult, "validate_families_with_parents"), (
        "SimulationResult.validate_families_with_parents regressed; "
        "Slice 3 backed out."
    )
    # The field-only validator stays a sibling, not replaced.
    assert hasattr(SimulationResult, "validate_families")

    # Not wired into validate_records=True: the `_compiled.py` clonal
    # run_records body must reference the field-only family raiser,
    # not a parent-aware one.
    import inspect
    from GenAIRR import _compiled
    src = inspect.getsource(_compiled.CompiledClonalExperiment.run_records)
    assert "_raise_on_family_validation_failure" in src, (
        "Clonal run_records no longer raises on field-only family "
        "validation failure; Slice 1 wiring regressed."
    )
    assert "validate_families_with_parents" not in src, (
        "Clonal run_records now invokes validate_families_with_parents "
        "automatically; Slice 3 deliberately kept this as an explicit "
        "diagnostic — flip the pin if the wiring is intentional."
    )


# ──────────────────────────────────────────────────────────────────
# 17. Absence — no clonal trace-file format
# ──────────────────────────────────────────────────────────────────


def test_pin_absence_no_clonal_trace_file_format() -> None:
    """``trace_file_from`` exists per ``CompiledSimulator`` (Rust)
    and bundles ``(plan, refdata, seed, trace)`` for a single
    outcome. No clonal-aware bundler that emits
    ``(pre_plan, post_plan, clone_seed, desc_seed, parent_trace,
    descendant_trace)`` exists. Slice 5+ adds it; not Slice 2.
    Pin the absence so a future contributor flipping it in Slice 2
    surfaces the scope creep."""
    clonal = _make_clonal().compile()
    for forbidden in (
        "trace_file_from_family",
        "clonal_trace_file_from",
        "family_trace_file",
        "trace_file_from_clone",
    ):
        assert not hasattr(clonal, forbidden), (
            f"CompiledClonalExperiment.{forbidden} now exists; Slice 5+ "
            "has landed early — flip pin and audit scope creep."
        )


# ──────────────────────────────────────────────────────────────────
# 18. Absence — no per-descendant SHM event aggregation
# ──────────────────────────────────────────────────────────────────


def test_pin_absence_no_per_descendant_shm_event_aggregation() -> None:
    """No API surfaces a per-clone or per-family view of SHM events.
    A consumer can read each descendant's ``outcome.events()`` and
    filter for mutate events, but there is no ``family_shm_events`` /
    ``mutations_per_clone`` aggregator. Slice 3+'s mutation-distance
    validator will need this; pin its absence so the surface is
    added deliberately in that slice rather than ad-hoc."""
    result = _make_clonal(with_mutate=True).run_records(seed=0)
    for forbidden in (
        "shm_events_per_clone",
        "mutations_per_clone",
        "family_shm_summary",
        "per_clone_shm",
    ):
        assert not hasattr(result, forbidden), (
            f"SimulationResult.{forbidden} now exists; flip pin."
        )


# ──────────────────────────────────────────────────────────────────
# 19. Doc anchor — audit doc exists + references contract file
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_audit_doc_exists_and_references_contract() -> None:
    """``docs/clonal_parent_outcome_design.md`` must continue to
    exist and reference this contract file. If the doc is deleted,
    the lockstep convention breaks before behaviour-level pins
    start drifting unexplained."""
    doc_path = _REPO_ROOT / "docs" / "clonal_parent_outcome_design.md"
    assert doc_path.exists(), "clonal_parent_outcome_design.md missing"
    doc = doc_path.read_text(encoding="utf-8")
    assert "test_clonal_parent_contract.py" in doc, (
        "Audit doc no longer references the contract file; the doc ↔ "
        "test lockstep convention has drifted."
    )
    # The 14-section structure should be intact.
    for marker in (
        "## 1. Current state",
        "## 4. Q4 — What object should eventually exist?",
        "## 13. Test surface",
        "## 14. Out of scope",
    ):
        assert marker in doc, f"audit doc missing section marker {marker!r}"
