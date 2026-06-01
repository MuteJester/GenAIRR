"""Contract pins for the clonal-family / lineage audit.

Companion to
[`docs/clonal_family_design.md`](../docs/clonal_family_design.md).
The audit is pre-implementation: it fixes today's behaviour
(``pin_scaffold_*``) and the architecture gaps it identifies
(``pin_absence_*``). No family-layer slice is proposed here; this
file is the baseline a future slice will land against.

The split:

- ``pin_scaffold_*`` tests freeze the surfaces today's clonal
  expansion ships: the ``_ClonalForkStep`` marker, the compile-
  time fork split into ``CompiledClonalExperiment``, the parent-IR
  sharing path (``pre.run`` → ``final_simulation`` →
  ``post.run_from``), the ``clone_id`` integer tag, the
  ``clone_seed = seed + clone_idx * 1_000_000`` scheme, and the
  ``CompiledSimulator.run_from`` Rust entry point.
- ``pin_absence_*`` tests freeze the gaps still open after Slice 1
  (the Python-only family validator): no parent-trace preservation,
  no parent ``Outcome`` per clone, no ``FamilyRecord`` /
  ``ClonalFamily`` / ``FamilyOutcome`` type, no
  ``SimulationResult.families`` attribute, no ``Outcome.parent_id``
  back-pointer, no ``run_family(seed, clone_idx)`` public method,
  no pre-SHM junction projection, no per-clone mutation-distance
  aggregation. The family-validator absence pin
  (``pin_absence_no_family_level_validator``) **flipped to
  ``pin_present_family_level_validator``** when Slice 1 landed —
  the rest still await their slices.

When a family-layer slice lands, the ``pin_absence_*`` tests flip
to ``pin_*_present`` / ``pin_*_used`` in lockstep with the slice.

Closing condition for the audit-only phase: every test here passes
on the current codebase. A failure here either means today's
clonal surface has drifted (``pin_scaffold_*`` regression) or
someone has shipped part of the family layer without flipping the
companion absence pin (``pin_absence_*`` regression).
"""
from __future__ import annotations

import inspect
import re
from pathlib import Path

import pytest

import GenAIRR as ga
from GenAIRR import _engine as ge
from GenAIRR import _compiled, _pipeline_ir, experiment as _experiment_module
from GenAIRR.result import SimulationResult, ValidationReport


# ──────────────────────────────────────────────────────────────────
# Shared fixtures
# ──────────────────────────────────────────────────────────────────


def _make_clonal_experiment(
    *, n_clones: int = 2, per_clone: int = 3, with_mutate: bool = False
) -> ga.Experiment:
    exp = ga.Experiment.on("human_igh").recombine().expand_clones(
        n_clones=n_clones, per_clone=per_clone
    )
    if with_mutate:
        exp = exp.mutate(count=5)
    return exp


def _records_grouped_by_clone(result: SimulationResult):
    grouped: dict = {}
    for rec in result:
        grouped.setdefault(rec["clone_id"], []).append(rec)
    return grouped


# ──────────────────────────────────────────────────────────────────
# 1. Scaffold — `_ClonalForkStep` marker shape
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_clonal_fork_step_is_frozen_dataclass() -> None:
    """``_ClonalForkStep`` is the structural marker the compile-time
    fork split keys off. It must be a frozen dataclass with exactly
    ``n_clones`` and ``size`` fields — a future refactor that adds
    runtime mutability or a field with a different name would break
    the partition logic in ``Experiment.compile()``."""
    step = _pipeline_ir._ClonalForkStep(n_clones=3, size=5)
    assert step.n_clones == 3
    assert step.size == 5
    # Frozen — mutation raises.
    with pytest.raises(Exception):  # FrozenInstanceError
        step.n_clones = 7  # type: ignore[misc]
    # Exact field set — defensive against silent additions.
    fields = {f.name for f in step.__dataclass_fields__.values()}  # type: ignore[attr-defined]
    assert fields == {"n_clones", "size"}, (
        f"_ClonalForkStep fields drifted; got {fields}"
    )


def test_pin_scaffold_expand_clones_appends_one_fork_step() -> None:
    """``Experiment.expand_clones`` appends exactly one
    ``_ClonalForkStep`` to ``_steps`` and rejects double-call. The
    audit's "single partition point" assumption rests on this."""
    exp = ga.Experiment.on("human_igh").recombine().expand_clones(
        n_clones=4, per_clone=6
    )
    fork_steps = [s for s in exp._steps if isinstance(s, _pipeline_ir._ClonalForkStep)]
    assert len(fork_steps) == 1
    assert fork_steps[0].n_clones == 4
    assert fork_steps[0].size == 6
    with pytest.raises(ValueError, match="only be called once"):
        exp.expand_clones(n_clones=2, per_clone=2)


# ──────────────────────────────────────────────────────────────────
# 2. Scaffold — compile-time fork split
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_compile_returns_clonal_when_fork_present() -> None:
    """``Experiment.compile()`` returns ``CompiledClonalExperiment``
    iff a ``_ClonalForkStep`` is in the step list, and the plain
    ``CompiledExperiment`` otherwise. The orchestration loop in
    ``CompiledClonalExperiment.run_records`` only runs on the
    clonal type — a refactor that returned ``CompiledExperiment`` for
    a clonal pipeline would silently lose the clone_id tag."""
    plain = ga.Experiment.on("human_igh").recombine().compile()
    clonal = (
        ga.Experiment.on("human_igh")
        .recombine()
        .expand_clones(n_clones=2, per_clone=2)
        .compile()
    )
    assert isinstance(plain, _compiled.CompiledExperiment)
    assert not isinstance(plain, _compiled.CompiledClonalExperiment)
    assert isinstance(clonal, _compiled.CompiledClonalExperiment)


def test_pin_scaffold_compiled_clonal_surface_methods_present() -> None:
    """``CompiledClonalExperiment`` exposes ``n_clones``, ``size``,
    ``total_records``, ``refdata``, ``run``, ``run_records``. The
    audit's API stability claim (§11 Backwards compatibility) requires
    these stay present. A future slice that introduces ``run_family``
    must add it without removing any of these."""
    clonal = _make_clonal_experiment(n_clones=3, per_clone=4).compile()
    assert clonal.n_clones == 3
    assert clonal.size == 4
    assert clonal.total_records == 12
    assert clonal.refdata is not None
    assert callable(getattr(clonal, "run", None))
    assert callable(getattr(clonal, "run_records", None))


# ──────────────────────────────────────────────────────────────────
# 3. Scaffold — parent IR sharing path (`pre.run` → `final_simulation`
#    → `post.run_from`)
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_orchestration_uses_run_from() -> None:
    """The orchestration loop in
    ``CompiledClonalExperiment.run_records`` must call ``run_from``
    on the post-fork simulator with the parent's
    ``final_simulation()``. A refactor that ran the post-fork plan
    via a fresh ``run()`` would silently break the "descendants share
    the same recombination" contract that test_g5 already pins."""
    src = inspect.getsource(_compiled.CompiledClonalExperiment.run_records)
    assert "self._pre.run(" in src
    assert ".final_simulation()" in src
    assert "self._post.run_from(" in src


def test_pin_scaffold_rust_run_from_exists_with_clonal_docstring() -> None:
    """The Rust-side ``CompiledSimulator.run_from`` is the entry
    point the orchestration loop calls; its docstring must still
    name "clonal expansion" so a future contributor doesn't
    remove the method on the assumption it's only used for replay."""
    src = (
        Path(__file__).resolve().parent.parent
        / "engine_rs"
        / "src"
        / "python"
        / "compiled.rs"
    ).read_text(encoding="utf-8")
    assert re.search(r"fn\s+run_from\s*\(", src), (
        "CompiledSimulator.run_from removed; orchestration loop will break"
    )
    assert "clonal expansion" in src, (
        "run_from docstring lost the 'clonal expansion' tag; the audit's "
        "vocabulary for this entry point drifted"
    )


# ──────────────────────────────────────────────────────────────────
# 4. Scaffold — `clone_id` tag + record count
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_clone_id_in_range_and_grouped_by_per_clone() -> None:
    """``clone_id`` is stamped in ``[0, n_clones)`` with exactly
    ``per_clone`` records per group. This is the "thin tag" view
    the audit's §11 promises to preserve after any family-layer
    slice."""
    result = _make_clonal_experiment(n_clones=3, per_clone=5).run_records(seed=0)
    assert len(result) == 15
    grouped = _records_grouped_by_clone(result)
    assert sorted(grouped.keys()) == [0, 1, 2]
    for cid, recs in grouped.items():
        assert len(recs) == 5, (
            f"clone {cid} has {len(recs)} descendants, expected 5"
        )


# ──────────────────────────────────────────────────────────────────
# 5. Scaffold — `clone_seed = seed + clone_idx * 1_000_000` scheme
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_clone_seed_formula_is_stable() -> None:
    """The seed-derivation formula
    ``clone_seed = int(seed) + clone_idx * 1_000_000`` is currently
    internal but visible by behavioural side-effect: two runs with
    the same ``seed`` reproduce byte-for-byte. The audit's §11
    flags this as **internal but stable** — a future change to the
    scheme breaks cross-version byte-reproducibility for clonal
    batches. A future ``run_family(seed, clone_idx)`` slice should
    expose the formula instead of changing it."""
    src = inspect.getsource(_compiled.CompiledClonalExperiment.run_records)
    assert "clone_idx * 1_000_000" in src, (
        "Clonal seed scheme has changed; cross-version byte-"
        "reproducibility for clonal batches is now broken."
    )
    # And the per-descendant offset stays at `clone_seed + 1 + desc_idx`.
    assert "clone_seed + 1 + desc_idx" in src


def test_pin_scaffold_byte_reproducibility_same_seed_same_records() -> None:
    """The "two runs with same seed produce identical records"
    contract is what makes the clonal-batch sha256 comparison work
    in release validation. Pin it here so a refactor that
    accidentally introduced non-determinism (e.g. RNG re-seeding)
    surfaces immediately."""
    a = _make_clonal_experiment(with_mutate=True).run_records(seed=42)
    b = _make_clonal_experiment(with_mutate=True).run_records(seed=42)
    assert [r["sequence"] for r in a] == [r["sequence"] for r in b]
    assert [r["clone_id"] for r in a] == [r["clone_id"] for r in b]
    assert [r["sequence_id"] for r in a] == [r["sequence_id"] for r in b]


# ──────────────────────────────────────────────────────────────────
# 6. Scaffold — parent-IR-driven invariants currently surfaced via
#    expose_provenance
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_truth_alleles_invariant_across_descendants() -> None:
    """``truth_v_call`` / ``truth_d_call`` / ``truth_j_call`` are
    family-invariant under the current architecture because they
    come from the parent IR's assignments. This is the §3 baseline
    a future family validator's "shared truth allele" check will
    enforce automatically; today consumers have to assert it
    themselves."""
    result = (
        _make_clonal_experiment(n_clones=2, per_clone=4, with_mutate=True)
        .run_records(seed=0, expose_provenance=True)
    )
    grouped = _records_grouped_by_clone(result)
    for cid, recs in grouped.items():
        truth_v = {r["truth_v_call"] for r in recs}
        truth_d = {r["truth_d_call"] for r in recs}
        truth_j = {r["truth_j_call"] for r in recs}
        assert len(truth_v) == 1, f"clone {cid} truth_v divergent: {truth_v}"
        assert len(truth_d) == 1, f"clone {cid} truth_d divergent: {truth_d}"
        assert len(truth_j) == 1, f"clone {cid} truth_j divergent: {truth_j}"


# ──────────────────────────────────────────────────────────────────
# 7. Absence — parent trace is not carried on any descendant
# ──────────────────────────────────────────────────────────────────


def test_pin_absence_parent_trace_not_carried_on_descendants() -> None:
    """A descendant's ``outcome.trace`` carries only the post-fork
    sampling sites — the pre-fork (recombination) addressed choices
    are gone after ``parent.final_simulation()``. Two descendants
    of the SAME clone have identically-sized traces (no parent
    history per-descendant); a descendant of a non-clonal recombine-
    only run has a strictly LARGER trace (it carries the pre-fork
    events directly). If this asymmetry disappears, a future
    contributor has either (a) lifted parent events onto
    descendant traces — flip this pin to ``pin_descendants_carry_
    parent_trace``, or (b) regressed the non-clonal trace to drop
    pre-fork events, which is a real bug."""
    # Clonal: recombine pre-fork, mutate post-fork.
    clonal_result = (
        _make_clonal_experiment(n_clones=1, per_clone=2, with_mutate=True)
        .run_records(seed=0)
    )
    clonal_outcomes = clonal_result.outcomes
    assert clonal_outcomes is not None
    desc_trace_lens = [len(o.trace()) for o in clonal_outcomes]
    # Descendants of a single clone all see only post-fork events,
    # so their traces are commensurate (same plan, different seeds).
    # We don't pin equality because mutate may sample different
    # numbers of substitutions per descendant; we pin the structural
    # ceiling: trace length must NOT include the recombination
    # passes' sampling sites.
    non_clonal = (
        ga.Experiment.on("human_igh").recombine().mutate(count=5).compile()
    )
    non_clonal_outcome = non_clonal.run(n=1, seed=0)[0]
    non_clonal_trace_len = len(non_clonal_outcome.trace())
    # Non-clonal trace carries BOTH recombine + mutate sites; clonal
    # descendant trace carries ONLY the post-fork mutate sites. So
    # descendant trace length must be strictly less than non-clonal.
    for desc_len in desc_trace_lens:
        assert desc_len < non_clonal_trace_len, (
            f"clonal descendant trace ({desc_len}) is not strictly "
            f"smaller than non-clonal trace ({non_clonal_trace_len}); "
            "parent recombination events may have leaked onto the "
            "descendant trace. If this is intentional, flip the audit's "
            "§5 'parent trace not preserved' gap pin."
        )


# ──────────────────────────────────────────────────────────────────
# 8. Absence — no parent `Outcome` exposed per clone
# ──────────────────────────────────────────────────────────────────


def test_pin_present_parent_outcome_per_clone() -> None:
    """Slice 2 shipped: ``SimulationResult.parents`` exposes one
    parent ``Outcome`` per clone. The flat ``.outcomes`` list
    still carries only descendant outcomes — parents live
    exclusively on ``.parents``. ``CompiledClonalExperiment``
    still does not expose a direct ``parent_outcome_for(clone_idx)``
    helper; consumers index via ``result.parents[record["parent_id"]]``.

    Flipped from ``pin_absence_no_parent_outcome_per_clone`` when
    the parent-outcome read-only surface landed."""
    result = _make_clonal_experiment(n_clones=3, per_clone=2).run_records(seed=0)
    # Slice 2 surface present.
    assert hasattr(result, "parents"), (
        "SimulationResult.parents regressed; Slice 2 backed out."
    )
    assert result.parents is not None
    assert len(result.parents) == 3
    # ``.outcomes`` continues to be the flat descendant list.
    assert result.outcomes is not None
    assert len(result.outcomes) == 6
    # The direct helper still doesn't ship — that's Slice 3+.
    clonal = _make_clonal_experiment(n_clones=3, per_clone=2).compile()
    assert not hasattr(clonal, "parent_outcome_for"), (
        "CompiledClonalExperiment.parent_outcome_for now exists; "
        "Slice 3+ has landed — flip pin."
    )


def test_pin_absence_no_run_family_public_method() -> None:
    """No public ``run_family(seed, clone_idx) -> (parent, descendants)``
    entry point — consumers who want just one family today have to
    reach through ``._pre`` / ``._post`` and replicate the seed
    formula themselves. The audit's §9 Replay flags this; Slice 2
    flips it."""
    clonal = _make_clonal_experiment().compile()
    assert not hasattr(clonal, "run_family"), (
        "`run_family` is now public; flip this pin to "
        "`pin_run_family_returns_parent_and_descendants`."
    )


# ──────────────────────────────────────────────────────────────────
# 9. Absence — no pre-SHM junction projection
# ──────────────────────────────────────────────────────────────────


def test_pin_absence_pre_shm_junction_not_projected() -> None:
    """No AIRR field exposes the pre-SHM junction. ``junction``
    is the observed (post-SHM) junction and diverges per descendant
    under SHM; the pre-SHM junction would be family-invariant but
    is only reachable via ``parent_sim`` (which is never returned).
    The audit's §3 / §6 family validator can't run the "shared
    pre-SHM junction" check until this gap closes."""
    result = (
        _make_clonal_experiment(n_clones=1, per_clone=3, with_mutate=True)
        .run_records(seed=0)
    )
    sample = result[0]
    # If a future projection adds a pre-SHM junction column, it
    # would surface here with one of these candidate names.
    forbidden = {
        "junction_pre_shm",
        "pre_shm_junction",
        "parent_junction",
        "junction_germline",
        "germline_junction",
        "junction_aa_pre_shm",
    }
    intersection = forbidden & set(sample.keys())
    assert not intersection, (
        f"AIRR record now projects pre-SHM junction column(s) {intersection}; "
        "audit §3 gap pin should flip to `pin_pre_shm_junction_projected`."
    )


# ──────────────────────────────────────────────────────────────────
# 10. Absence — no per-clone mutation-distance aggregation API
# ──────────────────────────────────────────────────────────────────


def test_pin_absence_per_descendant_mutation_distance_not_aggregated() -> None:
    """No API surfaces a per-clone mutation-distance distribution.
    ``SimulationResult`` has ``records`` and ``outcomes`` only; no
    ``mutation_distances_per_clone()`` / ``family_distance_distribution()``.
    Consumers today must group by ``clone_id`` and compute Hamming
    distances themselves. The audit's §6.6 family-validator check
    needs this surface."""
    result = _make_clonal_experiment(with_mutate=True).run_records(seed=0)
    for forbidden in (
        "mutation_distances_per_clone",
        "family_distance_distribution",
        "clone_distance_matrix",
        "per_family_shm_distance",
    ):
        assert not hasattr(result, forbidden), (
            f"SimulationResult exposes `{forbidden}`; audit §6.6 may have "
            "landed — flip this pin."
        )


# ──────────────────────────────────────────────────────────────────
# 11. Absence — `validate_records` ignores `clone_id`; no family
#     validator exists
# ──────────────────────────────────────────────────────────────────


def test_pin_validate_records_stays_per_record_only() -> None:
    """``SimulationResult.validate_records`` iterates
    ``(outcome, record)`` pairs independently and stays blind to
    cross-record family invariants by design — the family-layer
    check lives on the sibling ``validate_families`` (shipped in
    Slice 1). This pin freezes the **separation of concerns**: a
    refactor that folded family checks into ``validate_records``
    would surface here.

    To prove the separation: a record-list with a synthetic
    ``truth_v_call`` divergence within a clone — biologically
    impossible but constructible by mutating returned dicts — is
    NOT caught by ``validate_records``. ``validate_families``
    catches it; that's covered in
    [`tests/test_validate_families.py`](test_validate_families.py)."""
    refdata = ga.Experiment.on("human_igh").compile().refdata
    result = (
        _make_clonal_experiment(n_clones=2, per_clone=3, with_mutate=False)
        .run_records(seed=0, expose_provenance=True)
    )

    # Sanity: untampered batch validates clean.
    report = result.validate_records(refdata)
    assert isinstance(report, ValidationReport)
    assert report.ok, (
        f"baseline clonal batch failed projection validation: "
        f"{report.failures[:3]}"
    )

    # Tamper: force truth_v_call divergence within clone 0.
    grouped = _records_grouped_by_clone(result)
    if grouped[0][0]["truth_v_call"] == "IGHV1-2*02":
        grouped[0][1]["truth_v_call"] = "IGHV1-2*04"  # arbitrary other allele
    else:
        grouped[0][1]["truth_v_call"] = "IGHV1-2*02"

    # Re-run validation. Per-record validator looks at outcome state,
    # not record dicts, so it is blind to a dict-level tamper. The
    # tampered batch STILL validates clean — exactly the gap.
    tampered_report = result.validate_records(refdata)
    assert tampered_report.ok, (
        "validate_records caught a within-clone truth_v_call divergence; "
        "this would indicate a family-layer slice has landed. Flip this "
        "pin to `pin_family_validator_catches_truth_v_divergence`."
    )


def test_pin_present_family_level_validator() -> None:
    """Slice 1 (Python-only family validator) shipped. Pin the
    positive surface so a future refactor that removes it surfaces
    here in lockstep with the contract.

    Flipped from ``pin_absence_no_family_level_validator`` when the
    family-layer slice landed: ``SimulationResult.validate_families``
    is present, and ``_validation`` exports
    ``FamilyValidationFailedError`` + the raise helper. See
    [`tests/test_validate_families.py`](test_validate_families.py)
    for the spec-driven behavioural coverage; this pin is the
    audit-doc lockstep counterpart."""
    assert hasattr(SimulationResult, "validate_families"), (
        "SimulationResult.validate_families regressed; the family-"
        "validator slice has been backed out."
    )
    from GenAIRR import _validation
    assert hasattr(_validation, "FamilyValidationFailedError"), (
        "FamilyValidationFailedError missing from _validation."
    )
    assert hasattr(_validation, "_raise_on_family_validation_failure"), (
        "_raise_on_family_validation_failure helper missing; the "
        "clonal `validate_records=True` integration path is broken."
    )


# ──────────────────────────────────────────────────────────────────
# 12. Absence — no `FamilyRecord` / `ClonalFamily` / `FamilyOutcome`
#     type in the public package
# ──────────────────────────────────────────────────────────────────


def test_pin_absence_no_familyrecord_type() -> None:
    """No first-class family type is exported from the public package.
    The audit's §7 verdict — "we need a real ``FamilyRecord`` /
    ``ClonalOutcome`` concept before adding more biology" — rests on
    this absence. Adding the type is a Slice 3 deliverable; this pin
    flips when it lands."""
    forbidden = {
        "FamilyRecord",
        "ClonalFamily",
        "FamilyOutcome",
        "ClonalOutcome",
        "ClonalFamilyRecord",
    }
    exported = set(getattr(ga, "__all__", []))
    intersection = forbidden & exported
    assert not intersection, (
        f"GenAIRR now exports family type(s) {intersection}; audit §7 "
        "Slice 3 has landed — flip this pin."
    )
    # Also assert they aren't defined anywhere obvious in the public namespace.
    for name in forbidden:
        assert not hasattr(ga, name), (
            f"GenAIRR.{name} now exists at the package level; flip pin."
        )


def test_pin_absence_simulationresult_has_no_families_attr() -> None:
    """``SimulationResult`` does not expose a ``.families`` attribute.
    Today consumers group by ``clone_id`` themselves; Slice 3 adds
    ``.families`` as the family-record list."""
    result = _make_clonal_experiment().run_records(seed=0)
    assert not hasattr(result, "families"), (
        "SimulationResult.families now exists; flip pin."
    )
    # And the class itself doesn't promise it.
    assert "families" not in SimulationResult.__slots__


def test_pin_absence_outcome_has_no_parent_id() -> None:
    """``Outcome`` carries no ``parent_id`` / ``family_id`` /
    ``clone_id`` back-pointer. The ``clone_id`` int that ends up on
    AIRR records is set at projection time in
    ``CompiledClonalExperiment.run_records``, not on the underlying
    ``Outcome``. Slice 3 adds the back-pointer."""
    result = _make_clonal_experiment().run_records(seed=0)
    assert result.outcomes is not None
    sample_outcome = result.outcomes[0]
    for forbidden in ("parent_id", "family_id", "clone_id", "parent_outcome_id"):
        assert not hasattr(sample_outcome, forbidden), (
            f"Outcome.{forbidden} now exists; audit §7 Slice 3 may have "
            "landed — flip pin."
        )


# ──────────────────────────────────────────────────────────────────
# 13. Scaffold — pre-/post-fork plan partitioning is structural,
#     not runtime
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_pre_post_steps_partition_at_fork() -> None:
    """``CompiledClonalExperiment._pre_steps`` and ``_post_steps``
    are exactly the slice of ``Experiment._steps`` before/after the
    ``_ClonalForkStep``. The audit's vocabulary
    ("per-clone passes" = pre-fork, "per-descendant passes" =
    post-fork) hinges on this structural split."""
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .expand_clones(n_clones=2, per_clone=2)
        .mutate(count=5)
    )
    clonal = exp.compile()
    # The pre-fork plan has at least the recombine step.
    assert len(clonal._pre_steps) >= 1
    # The post-fork plan has at least the mutate step.
    assert len(clonal._post_steps) >= 1
    # Neither side contains the fork marker itself (it's consumed at
    # split time, not lowered as a pass).
    assert not any(
        isinstance(s, _pipeline_ir._ClonalForkStep)
        for s in clonal._pre_steps
    )
    assert not any(
        isinstance(s, _pipeline_ir._ClonalForkStep)
        for s in clonal._post_steps
    )


def test_pin_scaffold_post_fork_inherits_assembled_ir() -> None:
    """With no post-fork passes, every descendant of a clone is
    byte-identical because the post-fork plan inherits the parent's
    fully-assembled IR. This is the "identity baseline" the family
    validator (when it lands) must accept as a degenerate but valid
    case (§8 edge case 2)."""
    result = (
        ga.Experiment.on("human_igh")
        .recombine()
        .expand_clones(n_clones=2, per_clone=4)
        .run_records(seed=0)
    )
    grouped = _records_grouped_by_clone(result)
    for cid, recs in grouped.items():
        seqs = {r["sequence"] for r in recs}
        assert len(seqs) == 1, (
            f"clone {cid} descendants diverged without post-fork passes — "
            "the parent-IR-inheritance invariant has broken."
        )


# ──────────────────────────────────────────────────────────────────
# 14. Scaffold — `validate_records=True` already accepted on the
#     clonal kwarg surface (per-record only, today)
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_validate_records_kwarg_accepted_on_clonal() -> None:
    """``CompiledClonalExperiment.run_records`` already accepts
    ``validate_records=True`` (it lands per-record validation on the
    flat record list). The audit's §10 plan is for the family-layer
    pass to slot in behind the same kwarg — no new public surface
    needed. Pin the kwarg here so a refactor that drops it from the
    clonal path is caught before Slice 1 lands."""
    sig = inspect.signature(_compiled.CompiledClonalExperiment.run_records)
    assert "validate_records" in sig.parameters, (
        "CompiledClonalExperiment.run_records no longer accepts "
        "`validate_records`; audit §10 integration plan broken."
    )
    # And it must default to False (back-compat + perf).
    assert sig.parameters["validate_records"].default is False
    # End-to-end: the kwarg actually runs validation without errors on
    # a clean batch.
    result = _make_clonal_experiment(with_mutate=False).run_records(
        seed=0, validate_records=True
    )
    assert isinstance(result, SimulationResult)


# ──────────────────────────────────────────────────────────────────
# 15. Doc anchor — the audit doc continues to exist + references the
#     contract file
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_audit_doc_exists_and_references_contract() -> None:
    """The audit doc ``docs/clonal_family_design.md`` must continue
    to exist and reference this contract file. If a future
    contributor deletes the doc, this surfaces the missing reference
    before behaviour-level pins start drifting unexplained."""
    doc_path = (
        Path(__file__).resolve().parent.parent
        / "docs"
        / "clonal_family_design.md"
    )
    assert doc_path.exists(), "clonal_family_design.md missing"
    doc = doc_path.read_text(encoding="utf-8")
    assert "test_clonal_family_contract.py" in doc, (
        "Audit doc no longer references the contract file; the "
        "doc ↔ test lockstep convention has drifted."
    )
    # Sanity that the 14-section structure is intact.
    for marker in (
        "## 1.",
        "## 7. Q6",
        "## 13. Test surface",
        "## 14. Out of scope",
    ):
        assert marker in doc, f"audit doc missing section marker {marker!r}"
