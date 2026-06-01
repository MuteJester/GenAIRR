"""Tests for the opt-in ``validate_records=True`` kwarg on
:meth:`Experiment.run_records`.

This slice wires
:meth:`GenAIRR.result.SimulationResult.validate_records` (the public
AIRR-output correctness gate) into the run path so CI can request
"simulate and fail fast if any record is internally inconsistent."

Coverage matrix:

1. Happy path — kwarg-on doesn't raise on a clean batch.
2. Backwards compat — kwarg defaults to False, off-path identical.
3. Signature — default is literally ``False`` (introspectable).
4. Failure path — the chokepoint helper raises a structured,
   machine-greppable :class:`RecordValidationFailedError` on a
   non-ok report.
5. Post-hoc TSV reload — the original "requires outcomes" error
   message on :meth:`SimulationResult.validate_records` is
   preserved unchanged by this slice.
"""
from __future__ import annotations

import inspect

import pytest

import GenAIRR as ga
from GenAIRR import (
    Experiment,
    RecordValidationFailedError,
    SimulationResult,
    ValidationReport,
)
from GenAIRR._validation import _raise_on_validation_failure


def _build_experiment() -> Experiment:
    """The chain shared by the happy-path tests. Productive-only
    so projected records have anchored junctions worth validating."""
    return (
        Experiment.on("human_igh")
        .recombine()
        .productive_only()
        .mutate(rate=0.03)
    )


# 1. Happy path ───────────────────────────────────────────────────────


def test_validate_records_true_passes_on_productive_full_stack():
    """A clean productive-only + SHM batch must validate cleanly
    end-to-end — kwarg on, no exception."""
    exp = _build_experiment()
    # Must not raise.
    result = exp.run_records(n=5, seed=0, validate_records=True)
    # Sanity — we actually got 5 records back, not an empty
    # short-circuit.
    assert len(result) == 5


# 2. Backwards compat ─────────────────────────────────────────────────


def test_validate_records_false_preserves_current_behavior():
    """``validate_records=False`` must be a strict no-op vs. the
    no-kwarg call: same seed → same records."""
    exp = _build_experiment()
    baseline = exp.run_records(n=5, seed=0)
    with_kwarg_off = exp.run_records(n=5, seed=0, validate_records=False)

    assert len(baseline) == len(with_kwarg_off) == 5
    # Compare a couple of representative fields per record. Equality
    # across the whole list would also work but reports better on
    # first mismatch this way.
    for i in range(5):
        for field in ("sequence_id", "sequence", "v_call", "j_call"):
            assert baseline[i].get(field) == with_kwarg_off[i].get(field), (
                f"record {i} field {field!r} differs between no-kwarg and "
                f"validate_records=False: {baseline[i].get(field)!r} vs "
                f"{with_kwarg_off[i].get(field)!r}"
            )


# 3. Signature default ────────────────────────────────────────────────


def test_validate_records_default_is_false():
    """The new kwarg must default to False — this is the
    backwards-compat contract; flipping the default would silently
    add validator cost to every existing caller."""
    sig = inspect.signature(Experiment.run_records)
    assert "validate_records" in sig.parameters, (
        "Experiment.run_records is missing the validate_records kwarg"
    )
    param = sig.parameters["validate_records"]
    assert param.default is False, (
        f"validate_records default should be False, got {param.default!r}"
    )
    # Keyword-only — consistent with sibling kwargs (n, seed, etc.).
    assert param.kind is inspect.Parameter.KEYWORD_ONLY


# 4. Failure path via the chokepoint helper ──────────────────────────


def test_validate_records_true_raises_on_tampered_validation_report():
    """Drive ``_raise_on_validation_failure`` directly with a
    synthetic non-ok report and confirm:

    - it raises :class:`RecordValidationFailedError`;
    - the message has the three documented machine-greppable lines
      (failing-count, summary, first failure).

    Doing it at the helper level rather than corrupting a real
    in-flight result keeps the test fast and focused on the
    raise-and-format contract — the kwarg threading and the helper
    raising are independently tested.
    """
    fake = ValidationReport(
        count=5,
        failures=[
            {
                "record_index": 2,
                "sequence_id": "seq2",
                "issues": [
                    {"kind": "ProductiveMismatch", "details": {}},
                    {"kind": "JunctionStop", "details": {}},
                ],
            },
            {
                "record_index": 4,
                "sequence_id": "seq4",
                "issues": [{"kind": "ProductiveMismatch", "details": {}}],
            },
        ],
    )

    with pytest.raises(RecordValidationFailedError) as exc_info:
        _raise_on_validation_failure(fake)

    msg = str(exc_info.value)
    # Machine-greppable line 1 — failing count over total.
    assert "2 failing records out of 5" in msg, msg
    # Machine-greppable line 2 — histogram of issue kinds.
    assert "summary: " in msg, msg
    assert "ProductiveMismatch" in msg and "JunctionStop" in msg, msg
    # Machine-greppable line 3 — first failure detail.
    assert "first failure: " in msg, msg
    assert "record_index=2" in msg, msg
    assert "sequence_id='seq2'" in msg, msg
    assert "ProductiveMismatch" in msg, msg

    # The report itself is reachable for programmatic inspection.
    assert exc_info.value.report is fake
    assert exc_info.value.report.ok is False
    # Subclass relationship matters for ``except RuntimeError`` callers.
    assert isinstance(exc_info.value, RuntimeError)


# 5. Post-hoc TSV reload still errors out ────────────────────────────


def test_validate_records_post_hoc_on_loaded_tsv_still_raises_requires_outcomes(
    tmp_path,
):
    """Saving to TSV and rehydrating a :class:`SimulationResult` from
    record dicts only must still raise the existing "requires
    outcomes" :class:`RuntimeError`. This slice must not regress
    that surface — TSV-only results never have engine state to
    validate against, by design.

    No native ``from_tsv`` loader is exported, so we read the file
    back with ``csv.DictReader`` and feed records into the public
    constructor with ``outcomes=None`` — the same shape a loader
    would produce.
    """
    import csv

    exp = _build_experiment()
    result = exp.run_records(n=3, seed=0)
    tsv_path = tmp_path / "batch.tsv"
    result.to_tsv(str(tsv_path))

    # Round-trip the TSV back into a records-only SimulationResult.
    with open(tsv_path, "r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        loaded_records = list(reader)
    loaded = SimulationResult(records=loaded_records, outcomes=None)

    with pytest.raises(RuntimeError, match="validate_records requires"):
        loaded.validate_records(exp.refdata)
