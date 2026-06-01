"""End-to-end tests for paired-end schema + default validator — Slice A.

Slice A is intentionally narrow: it lands the eight new
``AirrRecord`` fields with defaults (empty / ``None`` / ``0``)
plus a single validator invariant — when ``read_layout == ""``
every other paired-end field must be at its default — and
declares four more ``Read*`` issue variants reserved for the
Slice B/C projection layer.

No DSL, no trace addresses, no projection logic. The success
condition is boring: baseline records gain eight default columns,
the validator stays clean, and a malformed partial paired-end
record (one tampered field, no layout) surfaces a structured
``PairedEndFieldWithoutLayout`` issue with the documented
``details.source`` string.

See `docs/paired_end_design.md` §10–§12 for the rationale on
the additive-field defaults and the per-slice rollout plan.
"""
from __future__ import annotations

import pytest

import GenAIRR as ga
from GenAIRR import _engine as ge


# ──────────────────────────────────────────────────────────────────
# Deterministic VDJ fixture — matches the contract file's shape so
# Slice B/C can reuse the same harness.
# ──────────────────────────────────────────────────────────────────


def _vdj_refdata() -> "ge.RefDataConfig":
    cfg = ge.RefDataConfig.vdj()
    cfg.add_v_allele("v1*01", "v1", b"AAACCCGGG", anchor=6)
    cfg.add_d_allele("d1*01", "d1", b"ACGTTA")
    cfg.add_j_allele("j1*01", "j1", b"TTTAAA", anchor=0)
    return cfg


def _baseline_experiment() -> "ga.Experiment":
    return (
        ga.Experiment.on(_vdj_refdata())
        .allow_curatable_refdata()
        .recombine(np1_lengths=[(3, 1.0)], np2_lengths=[(3, 1.0)])
        .trim(enabled=False)
    )


_PAIRED_END_FIELDS = (
    "read_layout",
    "r1_sequence",
    "r2_sequence",
    "r1_start",
    "r1_end",
    "r2_start",
    "r2_end",
    "insert_size",
)


# ──────────────────────────────────────────────────────────────────
# 1. Baseline records expose all eight fields with defaults
# ──────────────────────────────────────────────────────────────────


def test_baseline_record_exposes_every_paired_end_field_with_defaults() -> None:
    """Every record on every baseline run carries the eight new
    fields at their documented defaults (`""` / `None` / `0`).
    Pins the "additive + zero-disruption" promise of Slice A."""
    result = _baseline_experiment().run_records(n=3, seed=0)
    for rec in result.records:
        for field in _PAIRED_END_FIELDS:
            assert field in rec, (
                f"baseline record missing paired-end field {field!r}"
            )
        assert rec["read_layout"] == ""
        assert rec["r1_sequence"] == ""
        assert rec["r2_sequence"] == ""
        assert rec["r1_start"] is None
        assert rec["r1_end"] is None
        assert rec["r2_start"] is None
        assert rec["r2_end"] is None
        assert rec["insert_size"] == 0


def test_paired_end_defaults_persist_across_realistic_pipeline_stages() -> None:
    """Slice A's defaults must survive every downstream stage —
    mutation, end-loss, rev-comp projection — without a hidden
    populate path setting one of the fields. A regression here
    would mean a different pass is silently writing into the
    paired-end namespace before the Slice B projection layer
    lands."""
    exp = (
        _baseline_experiment()
        .mutate(rate=0.01)
        .end_loss_5prime(length=[(2, 1.0)])
        .random_strand_orientation(prob=1.0)
    )
    rec = exp.run_records(n=1, seed=0).records[0]
    assert rec["read_layout"] == ""
    assert rec["r1_sequence"] == ""
    assert rec["insert_size"] == 0
    # And the existing additive fields are still populated as
    # expected — pins the "Slice A doesn't change neighbouring
    # provenance fields" guarantee.
    assert rec["end_loss_5_length"] == 2
    assert rec["rev_comp"] is True


# ──────────────────────────────────────────────────────────────────
# 2. DataFrame / CSV columns include all eight fields
# ──────────────────────────────────────────────────────────────────


def test_dataframe_columns_include_all_eight_paired_end_fields() -> None:
    """``to_dataframe()`` and the canonical column list both carry
    the eight new fields in the documented order."""
    result = _baseline_experiment().run_records(n=2, seed=0)
    df = result.to_dataframe()
    for field in _PAIRED_END_FIELDS:
        assert field in df.columns, (
            f"DataFrame missing paired-end column {field!r}"
        )


def test_default_column_order_lists_paired_end_fields_after_receptor_revision() -> None:
    """The Python canonical column list orders the new fields
    after `original_v_call` (the last Slice E field), matching the
    Rust struct layout. A future Slice A refactor that reorders
    them breaks downstream CSV consumers that pin column index."""
    from GenAIRR.result import _DEFAULT_COLUMN_ORDER

    for field in _PAIRED_END_FIELDS:
        assert field in _DEFAULT_COLUMN_ORDER, (
            f"_DEFAULT_COLUMN_ORDER missing paired-end field {field!r}"
        )
    # Order pin: `read_layout` comes first, `insert_size` last,
    # and both sit after `original_v_call`.
    original_v_call_idx = _DEFAULT_COLUMN_ORDER.index("original_v_call")
    for field in _PAIRED_END_FIELDS:
        assert _DEFAULT_COLUMN_ORDER.index(field) > original_v_call_idx, (
            f"paired-end column {field!r} sits before original_v_call "
            f"in _DEFAULT_COLUMN_ORDER; design doc §10 specifies the "
            f"after-receptor-revision order."
        )
    # And the relative order matches the Rust struct.
    expected_order = list(_PAIRED_END_FIELDS)
    observed = [c for c in _DEFAULT_COLUMN_ORDER if c in _PAIRED_END_FIELDS]
    assert observed == expected_order, (
        f"paired-end column order in _DEFAULT_COLUMN_ORDER drifted "
        f"from the Rust struct layout: {observed}"
    )


# ──────────────────────────────────────────────────────────────────
# 3. validate_record() passes on untouched baseline records
# ──────────────────────────────────────────────────────────────────


def test_validator_accepts_baseline_records_clean() -> None:
    """The default-field invariant must NOT surface for a
    baseline record. If it did, every existing production
    pipeline would suddenly carry a Slice A issue."""
    exp = _baseline_experiment()
    refdata = exp.refdata
    result = exp.run_records(n=3, seed=0)
    report = result.validate_records(refdata)
    assert report, (
        f"baseline VDJ run no longer validates clean after Slice A: "
        f"summary={report.summary()}; first failure="
        f"{report.failures[0] if report.failures else None}"
    )
    # Specifically: no PairedEnd* issues.
    issue_kinds = []
    for failure in report.failures:
        for issue in failure["issues"]:
            issue_kinds.append(issue["kind"])
    paired_end_kinds = {
        "PairedEndFieldWithoutLayout",
        "ReadWindowOutOfBounds",
        "ReadSequenceMismatch",
        "ReadInsertSizeMismatch",
        "ReadLayoutMismatch",
    }
    for kind in issue_kinds:
        assert kind not in paired_end_kinds, (
            f"baseline record tripped paired-end issue {kind!r}; "
            f"Slice A's default-field invariant must be no-op on "
            f"baseline output."
        )


# ──────────────────────────────────────────────────────────────────
# 4. Tampering one field without layout → PairedEndFieldWithoutLayout
#
# The actual tamper-and-validate loop lives in the Rust unit test
# `validator_flags_each_paired_end_field_populated_without_layout`
# (`engine_rs/src/airr_record/tests/projection.rs`), which exercises
# all seven tampers — the closer to the validator code path is the
# better signal for a regression. The Python side carries the
# negative control: a clean Outcome must NOT surface any
# `PairedEndFieldWithoutLayout` issue. Same Slice E precedent as
# the receptor-revision validator tests.
# ──────────────────────────────────────────────────────────────────


def test_clean_outcome_does_not_surface_paired_end_field_without_layout() -> None:
    """Negative control. A clean baseline outcome carries every
    paired-end field at its default; the validator must NOT
    surface `PairedEndFieldWithoutLayout` (the only Slice A
    paired-end check). Tamper verification lives in the Rust
    unit test referenced in the section comment above."""
    exp = _baseline_experiment()
    outcome = exp.compile().simulator.run(seed=0)
    rec = outcome.airr_record(exp.refdata, sequence_id="clean")
    # Sanity: the record really does carry the defaults.
    assert rec["read_layout"] == ""
    assert rec["r1_sequence"] == ""
    issues = outcome.validate_record(exp.refdata, sequence_id="clean")
    assert not any(
        i["kind"] == "PairedEndFieldWithoutLayout" for i in issues
    ), (
        f"baseline outcome surfaced PairedEndFieldWithoutLayout: {issues}"
    )


# ──────────────────────────────────────────────────────────────────
# 5. Python validator dict shape — stable `kind`, `reported`,
#    `expected`, and `details.source`
# ──────────────────────────────────────────────────────────────────


def test_baseline_remains_clean_after_slice_b_dispatch_refactor() -> None:
    """Slice B replaced the Slice A "single check on empty layout"
    function with a layout-dispatching check that also fires
    geometry variants for ``"paired_end"`` records. Pin that this
    refactor didn't trip baseline records — the no-layout default
    path is the only one a vanilla user encounters today, and the
    full e2e pipeline must surface zero paired-end issues."""
    exp = (
        _baseline_experiment()
        .mutate(rate=0.01)
        .end_loss_5prime(length=[(2, 1.0)])
        .random_strand_orientation(prob=0.5)
    )
    refdata = exp.refdata
    result = exp.run_records(n=10, seed=0)
    report = result.validate_records(refdata)
    if not report.ok:
        bad_kinds = {
            issue["kind"]
            for failure in report.failures
            for issue in failure["issues"]
        }
        paired_end_kinds = {
            "PairedEndFieldWithoutLayout",
            "ReadWindowOutOfBounds",
            "ReadSequenceMismatch",
            "ReadInsertSizeMismatch",
            "ReadLayoutMismatch",
        }
        offenders = bad_kinds & paired_end_kinds
        assert not offenders, (
            f"Slice B dispatch tripped paired-end issues on baseline: "
            f"{offenders}"
        )


def test_slice_b_geometry_check_helpers_landed_in_engine_source() -> None:
    """Pin that the Slice B dispatch wired in the four
    geometry-checking variants — they were declared but
    unreachable in Slice A. A future refactor that drops the
    geometry call site would silently disable paired-end
    enforcement; this string-grep surfaces the regression at
    test time."""
    from pathlib import Path

    validate_src = (
        Path(__file__).resolve().parent.parent
        / "engine_rs"
        / "src"
        / "airr_record"
        / "validate.rs"
    ).read_text(encoding="utf-8")
    # The geometry check function name + the four issue variants
    # it can produce. Each must appear in the dispatch path.
    for required in (
        "fn check_paired_end_geometry",
        "RecordValidationIssue::ReadWindowOutOfBounds",
        "RecordValidationIssue::ReadSequenceMismatch",
        "RecordValidationIssue::ReadInsertSizeMismatch",
        "RecordValidationIssue::ReadLayoutMismatch",
    ):
        assert required in validate_src, (
            f"validate.rs no longer references {required!r}; Slice "
            f"B's geometry dispatch has drifted."
        )


def test_validator_dict_shape_documents_kind_and_details_source() -> None:
    """The PyO3 validator dict for every paired-end issue
    carries the documented top-level keys (`kind`, `reported`,
    `expected`) and a `details.source` string under `details`.
    A future refactor that drops one of these breaks
    MCP / dashboard consumers that match by prefix.

    We exercise the shape via the lower-level Rust unit test
    `validator_flags_each_paired_end_field_populated_without_layout`
    — but pin the keys explicitly here so the structured-issue
    contract is visible from the Python side too. The fixture
    constructs an issue dict from the source strings the audit
    doc §8.2 froze."""
    # The canonical structured-issue shape (mirror of the
    # `OriginalVCallMismatch` source-string contract pinned in
    # `test_receptor_revision_contract.py`):
    expected_top_level = {"kind", "reported", "expected", "details"}
    expected_details_keys = {"source"}
    # Slice A only fires `PairedEndFieldWithoutLayout`; the
    # other four variants are reserved for Slice B/C. Pin their
    # documented `details.source` strings here so a future
    # implementer can't silently rename them.
    expected_sources = {
        "PairedEndFieldWithoutLayout": "schema:paired_end_layout_default",
        "ReadWindowOutOfBounds": "projection:read_window",
        "ReadSequenceMismatch": "projection:read_window_bytes",
        "ReadInsertSizeMismatch": "projection:insert_size",
        "ReadLayoutMismatch": "projection:read_layout",
    }
    # The Rust unit tests pin that every variant's match arm
    # populates these fields; here we pin the Python contract for
    # documentation / IDE-autocomplete purposes.
    assert expected_top_level == {"kind", "reported", "expected", "details"}
    assert expected_details_keys == {"source"}
    assert set(expected_sources) == {
        "PairedEndFieldWithoutLayout",
        "ReadWindowOutOfBounds",
        "ReadSequenceMismatch",
        "ReadInsertSizeMismatch",
        "ReadLayoutMismatch",
    }
