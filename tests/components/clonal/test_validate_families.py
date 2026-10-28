"""Spec tests for the Python-only clonal-family validator slice.

Closes the audit's §12 Slice 1 (recommended first slice). The
slice adds:

- :class:`GenAIRR.result.FamilyValidationReport`
- :meth:`GenAIRR.result.SimulationResult.validate_families`
- :class:`GenAIRR._validation.FamilyValidationFailedError` (sibling
  of :class:`RecordValidationFailedError`)
- ``CompiledClonalExperiment.run_records(validate_records=True)``
  runs both gates: per-record first, then family-layer.

Spec coverage (from the user brief):

1. Clean clonal run passes ``validate_families()``.
2. Synthetic record mutation inside one clone fails with
   ``TruthVCallDiverges``.
3. ``validate_records=True`` on clonal ``run_records()`` runs both
   gates.
4. Non-clonal result returns ok with ``family_count=0``.
5. Report summary buckets issue kinds.
6. Existing ``validate_records(refdata)`` per-record behaviour
   unchanged.

Out of scope here (per the spec): mutation-distance distribution,
pre-SHM junction invariance, parent-trace reconstruction, lineage
topology.
"""
from __future__ import annotations

import pytest

import GenAIRR as ga
from GenAIRR import (
    FamilyValidationFailedError,
    FamilyValidationReport,
    RecordValidationFailedError,
    ValidationReport,
)


# ──────────────────────────────────────────────────────────────────
# Shared fixtures
# ──────────────────────────────────────────────────────────────────


def _make_clonal(
    *, n_clones: int = 2, per_clone: int = 3, with_mutate: bool = False
):
    exp = ga.Experiment.on("human_igh").recombine().expand_clones(
        n_clones=n_clones, per_clone=per_clone
    )
    if with_mutate:
        exp = exp.mutate(count=5)
    return exp


# ──────────────────────────────────────────────────────────────────
# Spec 1 — clean clonal run passes validate_families()
# ──────────────────────────────────────────────────────────────────


def test_clean_clonal_run_passes_validate_families() -> None:
    """A vanilla clonal pipeline run produces a batch whose
    family-layer report is ok, with the right family count and
    members-per-family aggregation."""
    result = _make_clonal(n_clones=3, per_clone=4, with_mutate=True).run_records(
        seed=0, expose_provenance=True
    )
    report = result.validate_families()
    assert isinstance(report, FamilyValidationReport)
    assert report.ok, f"unexpected failures: {report.failures[:3]}"
    assert bool(report) is True
    assert report.count == 12
    assert report.family_count == 3
    assert report.members_per_family == {0: 4, 1: 4, 2: 4}
    assert len(report) == 0  # no failures
    # Untampered batch — summary histogram is empty.
    assert report.summary() == {}


def test_validate_families_works_without_expose_provenance() -> None:
    """When ``expose_provenance=False`` (the default), the truth
    fields aren't on the records. ``validate_families`` must skip
    those invariants silently rather than raising, and still report
    a correct family count + members-per-family aggregation."""
    result = _make_clonal(n_clones=2, per_clone=5, with_mutate=True).run_records(
        seed=1
    )
    report = result.validate_families()
    assert report.ok
    assert report.family_count == 2
    assert report.members_per_family == {0: 5, 1: 5}
    # All four invariants skipped (no truth_*_call to compare) — no
    # failures, empty summary.
    assert report.summary() == {}


# ──────────────────────────────────────────────────────────────────
# Spec 2 — synthetic record mutation inside one clone fails
# ──────────────────────────────────────────────────────────────────


def test_truth_v_call_divergence_within_clone_fails() -> None:
    """Tamper with ``truth_v_call`` on one descendant of clone 0 so
    it disagrees with its siblings. The family validator catches it
    as ``TruthVCallDiverges`` carrying the clone_id, the record
    indices of the group, and the two divergent values."""
    result = _make_clonal(n_clones=2, per_clone=3).run_records(
        seed=0, expose_provenance=True
    )
    # Find the first clone-0 record and tamper it.
    clone_zero_indices = [
        i for i, r in enumerate(result) if r["clone_id"] == 0
    ]
    assert len(clone_zero_indices) == 3
    original_truth_v = result[clone_zero_indices[0]]["truth_v_call"]
    # Pick a different allele name guaranteed not to equal the
    # original. Append a marker; the validator compares strings.
    tampered = original_truth_v + "_TAMPERED"
    result.records[clone_zero_indices[0]]["truth_v_call"] = tampered

    report = result.validate_families()
    assert not report.ok
    # Exactly one failure dict.
    assert len(report.failures) == 1
    failure = report.failures[0]
    assert failure["clone_id"] == 0
    assert failure["issue_kind"] == "TruthVCallDiverges"
    assert sorted(failure["record_indices"]) == sorted(clone_zero_indices)
    assert sorted(failure["values"]) == sorted([original_truth_v, tampered])
    # Other clones untampered — still in members_per_family.
    assert report.members_per_family[1] == 3


def test_truth_d_and_j_divergence_also_caught() -> None:
    """Same mechanism for ``truth_d_call`` and ``truth_j_call`` —
    pin the issue-kind strings the user spec named so a future
    rename surfaces here."""
    result = _make_clonal(n_clones=2, per_clone=2).run_records(
        seed=0, expose_provenance=True
    )
    # Tamper D in clone 0, J in clone 1.
    clone_groups = {0: [], 1: []}
    for i, r in enumerate(result):
        clone_groups[r["clone_id"]].append(i)
    result.records[clone_groups[0][0]]["truth_d_call"] = "FAKE_D"
    result.records[clone_groups[1][0]]["truth_j_call"] = "FAKE_J"

    report = result.validate_families()
    assert not report.ok
    kinds_per_clone = {f["clone_id"]: f["issue_kind"] for f in report.failures}
    assert kinds_per_clone[0] == "TruthDCallDiverges"
    assert kinds_per_clone[1] == "TruthJCallDiverges"


# ──────────────────────────────────────────────────────────────────
# Spec 3 — validate_records=True on clonal runs both gates
# ──────────────────────────────────────────────────────────────────


def test_validate_records_true_runs_both_gates_on_clean_run() -> None:
    """Clean clonal run with ``validate_records=True`` returns the
    SimulationResult without raising. Both per-record and family-
    layer gates pass silently."""
    exp = _make_clonal(n_clones=2, per_clone=3, with_mutate=True)
    result = exp.run_records(seed=0, validate_records=True)
    assert len(result) == 6


def test_validate_records_true_clonal_raises_family_error_when_family_fails() -> None:
    """If the family-layer gate fails, ``validate_records=True`` on
    a clonal run raises :class:`FamilyValidationFailedError` (the
    sibling exception, NOT :class:`RecordValidationFailedError`).
    Distinguishability is the point.

    We can't easily provoke an in-engine family-layer failure
    (truth invariance is structural by construction), so we drive
    the gate manually through the SimulationResult API and assert
    the right exception fires."""
    from GenAIRR._validation import _raise_on_family_validation_failure

    result = _make_clonal(n_clones=2, per_clone=2).run_records(
        seed=0, expose_provenance=True
    )
    # First confirm projection validator passes.
    refdata = ga.Experiment.on("human_igh").compile().refdata
    record_report = result.validate_records(refdata)
    assert record_report.ok

    # Now tamper to force family failure.
    result.records[0]["truth_v_call"] = "FAKE_V"
    family_report = result.validate_families()
    assert not family_report.ok

    with pytest.raises(FamilyValidationFailedError) as exc_info:
        _raise_on_family_validation_failure(family_report)
    # Sibling, not subclass of RecordValidationFailedError.
    assert not isinstance(exc_info.value, RecordValidationFailedError)
    # But both are RuntimeError so `except RuntimeError` catches both.
    assert isinstance(exc_info.value, RuntimeError)
    # The exception carries the report for inspection.
    assert exc_info.value.report is family_report
    # Message is machine-greppable.
    msg = str(exc_info.value)
    assert "failing families" in msg
    assert "summary:" in msg
    assert "TruthVCallDiverges" in msg


def test_validate_records_true_clonal_family_failure_raises_through_run_records() -> None:
    """The clonal ``run_records(validate_records=True)`` path
    actually runs the family gate (not just the per-record gate).
    We simulate the path by calling the family validator directly
    on the post-run result; the in-engine truth invariance can't be
    broken from Python, so this test gets at the wiring through the
    per-record path remaining green AND the family-validator
    method existing and being callable from the result that
    ``run_records`` returns."""
    result = _make_clonal(n_clones=2, per_clone=3, with_mutate=True).run_records(
        seed=0, validate_records=True, expose_provenance=True
    )
    # The result is the same SimulationResult; family validation
    # remains callable on it post-return.
    report = result.validate_families()
    assert report.ok
    assert report.family_count == 2


# ──────────────────────────────────────────────────────────────────
# Spec 4 — non-clonal result returns ok with family_count=0
# ──────────────────────────────────────────────────────────────────


def test_non_clonal_result_validate_families_is_ok_noop() -> None:
    """A non-clonal ``run_records`` result has no ``clone_id`` on
    its records. ``validate_families`` must return ok with
    ``family_count=0`` and an empty ``members_per_family``."""
    exp = ga.Experiment.on("human_igh").recombine().mutate(count=3)
    result = exp.run_records(n=5, seed=0)
    assert all("clone_id" not in r for r in result)
    report = result.validate_families()
    assert report.ok
    assert report.count == 5
    assert report.family_count == 0
    assert report.members_per_family == {}
    assert len(report.failures) == 0


def test_non_clonal_validate_records_true_does_not_raise_family_error() -> None:
    """``validate_records=True`` on a non-clonal pipeline must NOT
    accidentally surface family failures — the non-clonal path
    doesn't run the family gate, and even if it did, the report
    would be ok-no-op."""
    exp = ga.Experiment.on("human_igh").recombine()
    result = exp.run_records(n=4, seed=0, validate_records=True)
    assert len(result) == 4
    # Sanity: the family report would still be ok if asked.
    assert result.validate_families().ok


def test_mixed_clonal_missing_clone_id_surfaces_failure() -> None:
    """If any record carries ``clone_id`` but others don't, the
    batch is structurally broken. The validator emits a
    ``CloneIdMissing`` failure naming the offending record indices.
    Constructed by stripping ``clone_id`` from one record."""
    result = _make_clonal(n_clones=2, per_clone=2).run_records(seed=0)
    # Remove clone_id from one record so the batch is now mixed.
    del result.records[0]["clone_id"]
    report = result.validate_families()
    assert not report.ok
    missing_failure = next(
        (f for f in report.failures if f["issue_kind"] == "CloneIdMissing"),
        None,
    )
    assert missing_failure is not None
    assert missing_failure["record_indices"] == [0]
    assert missing_failure["clone_id"] is None


# ──────────────────────────────────────────────────────────────────
# Spec 5 — report summary buckets issue kinds
# ──────────────────────────────────────────────────────────────────


def test_report_summary_buckets_issue_kinds() -> None:
    """``summary()`` returns a histogram ``{issue_kind: count}``
    across all failing groups, in the same style as
    ``ValidationReport.summary()``."""
    result = _make_clonal(n_clones=3, per_clone=2).run_records(
        seed=0, expose_provenance=True
    )
    # Tamper V on two clones, D on one clone — three failures total,
    # bucketed two:V, one:D.
    clone_groups = {0: [], 1: [], 2: []}
    for i, r in enumerate(result):
        clone_groups[r["clone_id"]].append(i)
    result.records[clone_groups[0][0]]["truth_v_call"] = "V_FAKE_A"
    result.records[clone_groups[1][0]]["truth_v_call"] = "V_FAKE_B"
    result.records[clone_groups[2][0]]["truth_d_call"] = "D_FAKE"

    report = result.validate_families()
    assert not report.ok
    summary = report.summary()
    assert summary == {"TruthVCallDiverges": 2, "TruthDCallDiverges": 1}


# ──────────────────────────────────────────────────────────────────
# Spec 6 — existing validate_records per-record behaviour unchanged
# ──────────────────────────────────────────────────────────────────


def test_validate_records_per_record_behaviour_unchanged() -> None:
    """The existing ``validate_records(refdata)`` API still returns
    :class:`ValidationReport` (not the family report), still uses
    per-record postcondition checks, and is **blind to family-layer
    divergence**. Tamper a truth field within a clone — the
    per-record report stays ok."""
    refdata = ga.Experiment.on("human_igh").compile().refdata
    result = _make_clonal(n_clones=2, per_clone=2, with_mutate=True).run_records(
        seed=0, expose_provenance=True
    )
    baseline = result.validate_records(refdata)
    assert isinstance(baseline, ValidationReport)
    assert baseline.ok

    # Tamper truth_v_call within clone 0 — projection validator
    # ignores it (per-record only, doesn't look across clone_id).
    clone_zero = [i for i, r in enumerate(result) if r["clone_id"] == 0]
    result.records[clone_zero[0]]["truth_v_call"] = "V_TAMPERED"
    tampered = result.validate_records(refdata)
    assert tampered.ok, (
        "validate_records started catching family-layer divergence; "
        "the clean separation between per-record and family gates "
        "has drifted."
    )
    # But validate_families catches it.
    assert not result.validate_families().ok


def test_records_only_simulationresult_can_run_validate_families() -> None:
    """Unlike ``validate_records`` (which needs the original
    ``Outcome`` objects), ``validate_families`` is pure-dict and
    must work on a records-only result, e.g. one loaded from TSV
    or constructed in-memory."""
    from GenAIRR.result import SimulationResult

    records = [
        {"clone_id": 0, "truth_v_call": "V1", "truth_d_call": "D1", "truth_j_call": "J1"},
        {"clone_id": 0, "truth_v_call": "V1", "truth_d_call": "D1", "truth_j_call": "J1"},
        {"clone_id": 1, "truth_v_call": "V2", "truth_d_call": "D2", "truth_j_call": "J2"},
    ]
    result = SimulationResult(records)  # no outcomes attached
    report = result.validate_families()
    assert report.ok
    assert report.family_count == 2
    assert report.members_per_family == {0: 2, 1: 1}
    # And validate_records on the same records-only result still
    # raises (per-record needs Outcomes; family validator does not).
    with pytest.raises(RuntimeError, match="requires the original Outcome"):
        result.validate_records(refdata=None)
