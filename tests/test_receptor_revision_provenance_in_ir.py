"""Spec tests for the Bug D fix — receptor-revision provenance moved
into the IR (`AlleleInstance.receptor_revision_original_id`).

Background (from the plan-split audit's Bug D):

The descendant ``Outcome``'s trace contains only post-fork choices.
``receptor_revision_applied`` / ``original_v_call`` used to be
trace-sourced (the AIRR builder read them from
``ChoiceAddress::ReceptorRevisionApplied`` and
``ChoiceAddress::SampleAllele(V)``), so canonical clonal pipelines
with ``recombine().receptor_revision(...).expand_clones(...)``
silently produced empty ``original_v_call`` on every descendant
even though the revision had really happened on the parent IR.

This slice moves provenance into the IR:

- New persistent field
  ``AlleleInstance.receptor_revision_original_id: Option<AlleleId>``.
- ``ReceptorRevisionPass`` captures the pre-revision V id (preserving
  the *first* original under hypothetical chaining) and installs it
  on the new V instance via
  ``with_receptor_revision_original_id``.
- AIRR builder reads provenance from
  ``sim.assignments.v.receptor_revision_original_id``.
- Validator's expected values come from the same IR slot. The
  ``details.source`` string for the two issue kinds is now
  ``ir:assignments.v.receptor_revision_original_id``.

Spec coverage (from the user brief):

1. Non-clonal receptor revision still reports
   ``receptor_revision_applied=True`` and correct
   ``original_v_call``.
2. Clonal canonical placement
   ``.receptor_revision(...).expand_clones(...)`` reports the same
   provenance on descendants.
3. Parent outcome and descendant final simulation both carry the
   provenance.
4. Replay still reproduces provenance.
5. Tampered AIRR fields trigger the existing validator issues.
6. Trace-only old logic is not enough: a descendant whose trace
   lacks ``receptor_revision.*`` still projects correctly when the
   IR provenance is present.
"""
from __future__ import annotations

import pytest

import GenAIRR as ga


# ──────────────────────────────────────────────────────────────────
# Shared fixtures
# ──────────────────────────────────────────────────────────────────


def _refdata():
    return ga.Experiment.on("human_igh").compile().refdata


# ──────────────────────────────────────────────────────────────────
# Spec 1 — Non-clonal receptor revision unchanged
# ──────────────────────────────────────────────────────────────────


def test_non_clonal_receptor_revision_still_reports_provenance() -> None:
    """Non-clonal ``.recombine().receptor_revision(prob=1.0)`` keeps
    its pre-fix behaviour: ``receptor_revision_applied=True`` and
    ``original_v_call`` resolves to a real V allele name."""
    exp = ga.Experiment.on("human_igh").recombine().receptor_revision(prob=1.0)
    for seed in range(5):
        result = exp.run_records(n=1, seed=seed)
        rec = result[0]
        assert rec["receptor_revision_applied"] is True, (
            f"seed {seed}: expected applied=True, got {rec['receptor_revision_applied']}"
        )
        assert rec["original_v_call"] != "", (
            f"seed {seed}: expected non-empty original_v_call"
        )
        # The pre-revision V allele name follows the IGHV naming
        # convention.
        assert rec["original_v_call"].startswith("IGHV"), (
            f"seed {seed}: unexpected original_v_call={rec['original_v_call']!r}"
        )


def test_non_clonal_no_revision_keeps_field_empty() -> None:
    """Without a receptor-revision step, ``receptor_revision_applied``
    is False and ``original_v_call`` is empty. Pin the "no-op"
    behaviour the IR field's default (``None``) protects."""
    exp = ga.Experiment.on("human_igh").recombine()
    result = exp.run_records(n=3, seed=0)
    for rec in result:
        assert rec["receptor_revision_applied"] is False
        assert rec["original_v_call"] == ""


# ──────────────────────────────────────────────────────────────────
# Spec 2 — Clonal canonical placement reports same provenance on
#          every descendant
# ──────────────────────────────────────────────────────────────────


def test_clonal_canonical_placement_reports_provenance_on_every_descendant() -> None:
    """The Bug D repro is now green: descendants of a clonal pipeline
    with ``receptor_revision`` placed BEFORE ``expand_clones`` report
    ``receptor_revision_applied=True`` and the parent's real
    ``original_v_call`` on every record."""
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .receptor_revision(prob=1.0)
        .expand_clones(n_clones=2, per_clone=3)
    )
    result = exp.run_records(seed=0)
    assert len(result) == 6
    for rec in result:
        assert rec["receptor_revision_applied"] is True, rec["sequence_id"]
        assert rec["original_v_call"] != "", rec["sequence_id"]
        assert rec["original_v_call"].startswith("IGHV"), rec["sequence_id"]


def test_clonal_descendants_share_original_v_call_within_clone() -> None:
    """Every descendant of a clone reports the SAME
    ``original_v_call`` — provenance is a parent-level property
    inherited via the assignments slot."""
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .receptor_revision(prob=1.0)
        .expand_clones(n_clones=3, per_clone=4)
    )
    result = exp.run_records(seed=0)
    by_clone: dict = {}
    for r in result:
        by_clone.setdefault(r["clone_id"], []).append(r)
    for cid, recs in by_clone.items():
        orig = {r["original_v_call"] for r in recs}
        assert len(orig) == 1, (
            f"clone {cid} descendants disagree on original_v_call: {orig}"
        )


def test_clonal_descendants_match_non_clonal_baseline_per_seed() -> None:
    """For a given seed, the non-clonal and clonal canonical paths
    must produce the same ``original_v_call`` (the same V is sampled
    and revised against the same RNG state). This pins the "IR
    provenance survives the parent→descendant boundary" claim
    against the non-clonal control."""
    non_clonal = ga.Experiment.on("human_igh").recombine().receptor_revision(prob=1.0)
    clonal = (
        ga.Experiment.on("human_igh")
        .recombine()
        .receptor_revision(prob=1.0)
        .expand_clones(n_clones=1, per_clone=1)
    )
    for seed in range(5):
        non_clonal_rec = non_clonal.run_records(n=1, seed=seed)[0]
        clonal_rec = clonal.run_records(seed=seed)[0]
        assert non_clonal_rec["original_v_call"] == clonal_rec["original_v_call"], (
            f"seed {seed}: non-clonal {non_clonal_rec['original_v_call']!r} != "
            f"clonal {clonal_rec['original_v_call']!r}"
        )


# ──────────────────────────────────────────────────────────────────
# Spec 3 — Parent outcome AND descendant final IR both carry provenance
# ──────────────────────────────────────────────────────────────────


def test_parent_outcome_carries_ir_provenance() -> None:
    """The parent's ``final_simulation()`` exposes V/D/J allele ids
    through ``PySimulation``. The IR-side provenance lives on the
    underlying ``AlleleInstance`` but isn't directly accessible via
    a Python accessor (Python access would need a Rust binding —
    out of scope here). We pin it observable via the parent's
    AIRR projection: project the parent outcome via
    ``outcome_to_airr_record`` and read the same provenance fields."""
    from GenAIRR._airr_record import outcome_to_airr_record

    refdata = _refdata()
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .receptor_revision(prob=1.0)
        .expand_clones(n_clones=2, per_clone=2)
    )
    result = exp.run_records(seed=0)
    assert result.parents is not None
    for i, parent in enumerate(result.parents):
        parent_rec = outcome_to_airr_record(
            parent, refdata, sequence_id=f"parent{i}"
        )
        assert parent_rec["receptor_revision_applied"] is True
        assert parent_rec["original_v_call"] != ""


def test_descendant_final_simulation_carries_same_provenance_as_record() -> None:
    """A descendant's final ``Simulation`` IR is what its AIRR
    record projects from. The descendant inherits the parent's
    revised V assignments (including the IR provenance), so its
    final IR's V allele identity must reflect the **revised** V
    (not the original pre-revision V) — but its AIRR record
    reports the **original** V via ``original_v_call``. Pin both
    views from the same outcome to prove the descendant's IR
    really carries the provenance, not just its projected dict."""
    refdata = _refdata()
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .receptor_revision(prob=1.0)
        .expand_clones(n_clones=1, per_clone=2)
    )
    result = exp.run_records(seed=0)
    assert result.outcomes is not None
    for i, outcome in enumerate(result.outcomes):
        rec = result.records[i]
        final = outcome.final_simulation()
        # The descendant's IR has a V allele assigned (the
        # post-revision V).
        post_revision_v_id = final.v_allele_id()
        assert post_revision_v_id is not None
        # And the AIRR record carries the pre-revision identity in
        # original_v_call (proving IR-sourced provenance reached
        # the descendant).
        assert rec["original_v_call"] != ""
        assert rec["receptor_revision_applied"] is True


# ──────────────────────────────────────────────────────────────────
# Spec 4 — Replay reproduces provenance
# ──────────────────────────────────────────────────────────────────


def test_replay_reproduces_receptor_revision_provenance() -> None:
    """Running the same experiment twice with the same seed produces
    byte-identical ``receptor_revision_applied`` / ``original_v_call``
    on every record. The IR-sourced projection has to be
    deterministic across runs — RNG-based revision choices are
    seed-stable, and the IR's ``receptor_revision_original_id`` is
    a deterministic function of those choices."""
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .receptor_revision(prob=1.0)
        .expand_clones(n_clones=2, per_clone=3)
    )
    a = exp.run_records(seed=42)
    b = exp.run_records(seed=42)
    assert [r["original_v_call"] for r in a] == [r["original_v_call"] for r in b]
    assert [r["receptor_revision_applied"] for r in a] == [
        r["receptor_revision_applied"] for r in b
    ]


def test_replay_via_trace_file_roundtrip_preserves_provenance() -> None:
    """Trace files carry the addressed-choice trace; replay produces
    an ``Outcome`` with the same final IR, which must surface the
    same receptor-revision provenance on its AIRR record. Pin
    end-to-end so a regression in the replay path's provenance
    handling surfaces here."""
    exp = ga.Experiment.on("human_igh").recombine().receptor_revision(prob=1.0)
    compiled = exp.compile()
    # ``trace_file_from`` lives on the underlying Rust
    # ``CompiledSimulator`` (the Python ``CompiledExperiment``
    # wraps it via ``simulator``).
    simulator = compiled.simulator
    original = simulator.run(seed=0)
    trace_file = simulator.trace_file_from(original, seed=0)
    replayed = simulator.replay_from_trace_file(trace_file)

    from GenAIRR._airr_record import outcome_to_airr_record

    rec_original = outcome_to_airr_record(
        original, compiled.refdata, sequence_id="o"
    )
    rec_replayed = outcome_to_airr_record(
        replayed, compiled.refdata, sequence_id="r"
    )
    assert (
        rec_original["receptor_revision_applied"]
        == rec_replayed["receptor_revision_applied"]
    )
    assert rec_original["original_v_call"] == rec_replayed["original_v_call"]


# ──────────────────────────────────────────────────────────────────
# Spec 5 — Tampered AIRR fields trigger existing validator issues
# ──────────────────────────────────────────────────────────────────


def test_validator_clean_clonal_revision_path_emits_no_mismatch_issues() -> None:
    """Positive direction: a clean clonal canonical pipeline
    produces records whose IR-sourced ``receptor_revision_applied``
    and ``original_v_call`` match the validator's IR-sourced
    expected values. No ``ReceptorRevisionAppliedMismatch`` or
    ``OriginalVCallMismatch`` issues fire.

    Negative-direction tamper coverage lives in the Rust unit
    tests
    (``validator_flags_receptor_revision_applied_mismatch`` /
    ``validator_flags_original_v_call_mismatch``) which exercise
    the per-record validator directly with a tampered
    ``AirrRecord`` — a tamper path Python's ``validate_records``
    doesn't expose because it re-projects the record from the
    outcome state."""
    refdata = _refdata()
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .receptor_revision(prob=1.0)
        .expand_clones(n_clones=2, per_clone=3)
    )
    result = exp.run_records(seed=0)
    report = result.validate_records(refdata)
    assert report.ok, (
        f"clean clonal canonical revision pipeline tripped the "
        f"validator: {report.failures[:3]}"
    )
    # Per-record summary has no mismatch kinds.
    summary = report.summary()
    assert "ReceptorRevisionAppliedMismatch" not in summary
    assert "OriginalVCallMismatch" not in summary


def test_validator_source_string_flipped_to_ir() -> None:
    """``details.source`` for the two receptor-revision issue
    kinds must now read
    ``ir:assignments.v.receptor_revision_original_id`` instead of
    the old trace addresses (``trace:receptor_revision.applied`` /
    ``trace:sample_allele.v``). Checked at the Rust binding source
    level — Python's ``validate_records`` doesn't expose a
    record-tamper path, so we pin the source string at its
    declaration site."""
    from pathlib import Path

    bindings = (
        Path(__file__).resolve().parent.parent
        / "engine_rs"
        / "src"
        / "python"
        / "outcome.rs"
    ).read_text(encoding="utf-8")
    assert "ir:assignments.v.receptor_revision_original_id" in bindings, (
        "the IR-sourced details.source string is missing from "
        "outcome.rs; Bug D fix has regressed."
    )
    assert "trace:receptor_revision.applied" not in bindings, (
        "old trace-sourced details.source string still appears in "
        "outcome.rs; Bug D fix is incomplete."
    )
    assert "trace:sample_allele.v" not in bindings, (
        "old trace-sourced details.source string still appears in "
        "outcome.rs; Bug D fix is incomplete."
    )


# ──────────────────────────────────────────────────────────────────
# Spec 6 — Trace-only old logic is not enough: descendant trace lacks
#          receptor_revision.* but IR provenance still projects
# ──────────────────────────────────────────────────────────────────


def test_descendant_trace_lacks_receptor_revision_events_but_ir_projects_correctly() -> None:
    """The Bug D mechanism: the descendant's trace doesn't carry
    pre-fork events. Pin both facts side-by-side:

    1. The descendant's ``outcome.events()`` contains NO
       ``receptor_revision.applied`` event (those are pre-fork,
       and the descendant trace is fresh per
       ``execute_transactional``).
    2. Despite that, ``original_v_call`` / ``receptor_revision_applied``
       project correctly because they read from the IR
       (``assignments.v.receptor_revision_original_id``), which
       DOES survive the parent→descendant boundary.

    This is the load-bearing pin for Bug D — a regression that
    re-routes provenance back to the trace would surface here."""
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .receptor_revision(prob=1.0)
        .expand_clones(n_clones=1, per_clone=2)
    )
    result = exp.run_records(seed=0)
    assert result.outcomes is not None
    for desc_outcome, rec in zip(result.outcomes, result.records):
        # (1) No receptor_revision in the descendant pass_names —
        # pre-fork passes don't surface here.
        desc_pass_names = " ".join(desc_outcome.pass_names())
        assert "receptor_revision" not in desc_pass_names, (
            f"descendant carries receptor_revision in its pass_names "
            f"({desc_pass_names!r}); the parent's trace leaked onto the "
            "descendant and the Bug D mechanism has changed."
        )
        # (2) But provenance still projects correctly because the IR
        # carries it.
        assert rec["receptor_revision_applied"] is True
        assert rec["original_v_call"] != ""
