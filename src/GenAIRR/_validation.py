"""Helpers for the opt-in ``validate_records=True`` run kwarg.

The public surface of GenAIRR's batch correctness check lives on
:class:`GenAIRR.result.SimulationResult.validate_records` (returns a
:class:`GenAIRR.result.ValidationReport`). This module wires that
report into the run-and-fail-fast path: when callers pass
``validate_records=True`` to :meth:`Experiment.run_records`, the
result is validated immediately and any failure raises a structured
:class:`RecordValidationFailedError` carrying a machine-greppable
message body.

Splitting this off ``result.py`` keeps the report object itself a
pure-data container (no exception coupling) and gives the run-side
plumbing a single, importable home — see ``_compiled.py`` and
``experiment.py`` for the ``run_records`` integrations.
"""
from __future__ import annotations

from typing import Any, List


class RecordValidationFailedError(RuntimeError):
    """Raised by ``run_records(validate_records=True)`` when one or
    more projected AIRR records fail the postcondition validator.

    Subclass of :class:`RuntimeError` because this is a runtime
    state violation (the engine emitted records inconsistent with
    their own outcomes) — matching the existing
    ``RuntimeError`` convention in
    :meth:`SimulationResult.validate_records` (see ``result.py``).

    The :attr:`report` attribute keeps the original
    :class:`~GenAIRR.result.ValidationReport` so callers that
    ``except RecordValidationFailedError as exc`` can inspect
    ``exc.report.failures`` / ``exc.report.summary()`` without
    re-parsing the message.

    The exception message has three machine-greppable lines, each
    on its own line so simple ``grep`` / regex tools can pull them
    out of CI logs:

    - ``"<N> failing records out of <total>"``
    - ``"summary: {kind: count, ...}"``
    - ``"first failure: record_index=<I>, sequence_id='<SID>', issue kinds=[...]"``
    """

    def __init__(self, report: Any) -> None:
        self.report = report
        super().__init__(_format_failure_message(report))


def _format_failure_message(report: Any) -> str:
    """Build the structured, machine-greppable message body for
    :class:`RecordValidationFailedError`. Kept as a free function
    so tests can call it directly on a synthetic report-like."""
    failures: List[Any] = list(report.failures)
    total = int(report.count)
    n_failing = len(failures)
    summary = report.summary()
    summary_repr = (
        "{"
        + ", ".join(f"{k}: {v}" for k, v in summary.items())
        + "}"
    )
    first = failures[0] if failures else None
    if first is not None:
        first_kinds = [
            issue.get("kind", "Unknown") for issue in first.get("issues", [])
        ]
        first_line = (
            f"first failure: record_index={first.get('record_index')}, "
            f"sequence_id={first.get('sequence_id')!r}, "
            f"issue kinds={first_kinds}"
        )
    else:
        first_line = "first failure: <none>"
    return (
        f"{n_failing} failing records out of {total}\n"
        f"summary: {summary_repr}\n"
        f"{first_line}"
    )


def _raise_on_validation_failure(report: Any) -> None:
    """Raise :class:`RecordValidationFailedError` if ``report`` is
    not ok; otherwise return silently. The single chokepoint the
    ``run_records(validate_records=True)`` code path calls into."""
    if report.ok:
        return
    raise RecordValidationFailedError(report)


class FamilyValidationFailedError(RuntimeError):
    """Raised by ``run_records(validate_records=True)`` on a
    clonal-pipeline result when one or more clone groups fail the
    family-consistency gate.

    Sibling to :class:`RecordValidationFailedError` — caught
    separately so users can distinguish projection bugs
    (per-record AIRR postconditions violated) from family-
    consistency bugs (recombination-time truth fields disagree
    across descendants of a clone). Both inherit
    :class:`RuntimeError` so a single ``except RuntimeError`` still
    catches both.

    The :attr:`report` attribute keeps the original
    :class:`~GenAIRR.result.FamilyValidationReport` so callers that
    ``except FamilyValidationFailedError as exc`` can inspect
    ``exc.report.failures`` / ``exc.report.summary()`` without
    re-parsing the message.

    The exception message has three machine-greppable lines, each
    on its own line so simple ``grep`` / regex tools can pull them
    out of CI logs:

    - ``"<N> failing families out of <family_count> (records=<count>)"``
    - ``"summary: {issue_kind: count, ...}"``
    - ``"first failure: clone_id=<C>, issue_kind=<K>, values=[...]"``
    """

    def __init__(self, report: Any) -> None:
        self.report = report
        super().__init__(_format_family_failure_message(report))


def _format_family_failure_message(report: Any) -> str:
    """Build the structured, machine-greppable message body for
    :class:`FamilyValidationFailedError`. Kept as a free function so
    tests can call it directly on a synthetic report-like."""
    failures: List[Any] = list(report.failures)
    total = int(report.count)
    family_count = int(report.family_count)
    n_failing = len(failures)
    summary = report.summary()
    summary_repr = (
        "{"
        + ", ".join(f"{k}: {v}" for k, v in summary.items())
        + "}"
    )
    first = failures[0] if failures else None
    if first is not None:
        first_line = (
            f"first failure: clone_id={first.get('clone_id')!r}, "
            f"issue_kind={first.get('issue_kind')!r}, "
            f"values={first.get('values')!r}"
        )
    else:
        first_line = "first failure: <none>"
    return (
        f"{n_failing} failing families out of {family_count} "
        f"(records={total})\n"
        f"summary: {summary_repr}\n"
        f"{first_line}"
    )


def _raise_on_family_validation_failure(report: Any) -> None:
    """Raise :class:`FamilyValidationFailedError` if ``report`` is
    not ok; otherwise return silently. The single chokepoint the
    clonal ``run_records(validate_records=True)`` code path calls
    into after the per-record gate has passed."""
    if report.ok:
        return
    raise FamilyValidationFailedError(report)
