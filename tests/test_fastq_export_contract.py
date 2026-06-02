"""Contract pins for the FASTQ / Read Export audit.

Companion to
[`docs/fastq_export_design.md`](../docs/fastq_export_design.md).
The audit is pre-implementation: it freezes today's surfaces
(``pin_scaffold_*``) and the gaps a future implementation slice
would close (``pin_absence_*``).

The implementation slice will add a single
``SimulationResult.to_paired_fastq(r1_path, r2_path, *,
quality="illumina", prefix="seq", **quality_kwargs)`` method
that mirrors the existing single-end ``to_fastq`` plumbing,
reuses the two pluggable quality models in ``_qmodel.py``, and
writes Illumina-suffixed (`/1` / `/2`) FASTQ headers per the
audit's §3 recommendation.

Split:

- ``pin_scaffold_*`` tests freeze the pre-existing surfaces the
  slice builds on: the five shipped export methods on
  ``SimulationResult``, the two pluggable quality models in
  ``_qmodel.py``, the paired-end AIRR fields with their
  documented defaults, the ``sequence_id`` field on every
  record, and the absence of any per-position Phred score in
  the engine.
- ``pin_absence_*`` tests freeze the gaps the slice closes:
  no ``to_paired_fastq`` method, no record-generator surface,
  no engine-side per-base quality field, no manifest block.
"""
from __future__ import annotations

import inspect
from pathlib import Path

import pytest

import GenAIRR as ga


_REPO_ROOT = Path(__file__).resolve().parent.parent
_AUDIT_DOC = _REPO_ROOT / "audit-docs" / "fastq_export_design.md"


# ──────────────────────────────────────────────────────────────────
# 1. Scaffold — five export methods exist on SimulationResult
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_existing_export_surfaces_present() -> None:
    """Audit §1 — every export method the audit references must
    exist on ``SimulationResult`` today. The new
    ``to_paired_fastq`` lands as a sixth method alongside these;
    losing any of them would mean the surface map drifted."""
    from GenAIRR.result import SimulationResult

    for method in (
        "to_dataframe",
        "to_csv",
        "to_tsv",
        "to_fasta",
        "to_fastq",
    ):
        assert hasattr(SimulationResult, method), (
            f"SimulationResult is missing {method!r}; audit §1 "
            "surface inventory drifted"
        )


def test_pin_scaffold_to_fastq_quality_prefix_kwarg_shape() -> None:
    """Audit §1 / §3 — ``to_fastq`` exposes the keyword-only
    ``quality`` / ``prefix`` parameters plus ``**quality_kwargs``.
    The new ``to_paired_fastq`` mirrors this signature verbatim
    so users carry the same muscle memory across surfaces."""
    from GenAIRR.result import SimulationResult

    sig = inspect.signature(SimulationResult.to_fastq)
    assert "quality" in sig.parameters
    assert "prefix" in sig.parameters
    assert sig.parameters["quality"].default == "illumina"
    assert sig.parameters["prefix"].default == "seq"
    # **quality_kwargs is a VAR_KEYWORD parameter — the resolve
    # helper forwards all extras to the quality model constructor.
    var_kw_params = [
        p
        for p in sig.parameters.values()
        if p.kind == inspect.Parameter.VAR_KEYWORD
    ]
    assert len(var_kw_params) == 1, (
        "to_fastq has no **quality_kwargs catch-all; the audit's "
        "kwarg-forwarding plan depends on it"
    )


# ──────────────────────────────────────────────────────────────────
# 2. Scaffold — pluggable quality models
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_constant_and_illumina_models_exist() -> None:
    """Audit §4 — both ``ConstantQualityModel`` and
    ``IlluminaQualityModel`` exist on ``_qmodel.py`` with the
    ``quality_array(seq)`` interface the audit relies on. The
    new ``to_paired_fastq`` slice will call ``quality_array``
    independently on R1 and R2."""
    from GenAIRR import _qmodel

    assert hasattr(_qmodel, "ConstantQualityModel")
    assert hasattr(_qmodel, "IlluminaQualityModel")
    # Both expose ``quality_array(sequence)``.
    const_model = _qmodel.ConstantQualityModel(q=30)
    q = const_model.quality_array("ACGTA")
    assert len(q) == 5, "ConstantQualityModel.quality_array must "
    "return a list the same length as the input sequence"
    illumina = _qmodel.IlluminaQualityModel()
    q2 = illumina.quality_array("A" * 150)
    assert len(q2) == 150
    # ``resolve_quality_model`` is the canonical helper the new
    # slice uses to instantiate either model from a kwarg name.
    assert hasattr(_qmodel, "resolve_quality_model")
    assert hasattr(_qmodel, "phred_to_ascii")


def test_pin_scaffold_quality_models_take_sequence_only() -> None:
    """Audit §3 — both quality models take the sequence string
    as their only argument; they do NOT read any AIRR field
    directly. This is what lets ``to_paired_fastq`` call
    ``model.quality_array(r1_sequence)`` and
    ``model.quality_array(r2_sequence)`` independently with no
    per-read state coordination."""
    from GenAIRR import _qmodel

    const_sig = inspect.signature(_qmodel.ConstantQualityModel.quality_array)
    illumina_sig = inspect.signature(_qmodel.IlluminaQualityModel.quality_array)
    # Each has exactly one positional parameter besides ``self``
    # — the sequence string.
    for name, sig in [("ConstantQualityModel", const_sig),
                      ("IlluminaQualityModel", illumina_sig)]:
        params = [
            p
            for p in sig.parameters.values()
            if p.name != "self"
        ]
        assert len(params) == 1, (
            f"{name}.quality_array took more than the sequence string; "
            "the audit's stateless-call argument regressed"
        )
        assert params[0].name == "sequence"


# ──────────────────────────────────────────────────────────────────
# 3. Scaffold — paired-end AIRR fields present with defaults
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_paired_end_fields_default_when_disabled() -> None:
    """Audit §3 / §7 — the eight paired-end AIRR fields ship on
    every record regardless of whether ``.paired_end()`` ran.
    When the step is omitted, the defaults are:

      read_layout = ""
      r1_sequence = ""        r2_sequence = ""
      r1_start = None         r2_start = None
      r1_end = None           r2_end = None
      insert_size = 0         (sentinel — integer, not None)

    The new ``to_paired_fastq`` writer uses
    ``read_layout == "paired_end"`` as the canonical "paired-end
    fields are real" signal; the asymmetric ``insert_size = 0``
    default is the trap pin §3 warns the implementer about."""
    result = (
        ga.Experiment.on("human_igh")
        .recombine()
        .run_records(n=1, seed=0)
    )
    rec = result.records[0]
    assert rec["read_layout"] == ""
    assert rec["r1_sequence"] == ""
    assert rec["r2_sequence"] == ""
    assert rec["r1_start"] is None
    assert rec["r1_end"] is None
    assert rec["r2_start"] is None
    assert rec["r2_end"] is None
    # The asymmetric int sentinel — NOT None.
    assert rec["insert_size"] == 0
    assert isinstance(rec["insert_size"], int)


def test_pin_scaffold_paired_end_fields_populated_when_enabled() -> None:
    """Audit §3 — when ``.paired_end()`` runs, every record
    carries ``read_layout == "paired_end"`` and non-empty
    ``r1_sequence`` / ``r2_sequence`` of the configured lengths.
    The new exporter uses this state."""
    result = (
        ga.Experiment.on("human_igh")
        .recombine()
        .paired_end(r1_length=120, insert_size=240)
        .run_records(n=3, seed=4242)
    )
    for rec in result.records:
        assert rec["read_layout"] == "paired_end"
        assert rec["r1_sequence"]
        assert rec["r2_sequence"]
        assert rec["insert_size"] > 0


# ──────────────────────────────────────────────────────────────────
# 4. Scaffold — sequence_id populated on every record
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_sequence_id_present_on_every_record() -> None:
    """Audit §5 — ``sequence_id`` is a non-empty string on every
    record after a successful run. The new ``to_paired_fastq``
    uses it for read names (``@{sequence_id}/1`` /
    ``@{sequence_id}/2``)."""
    result = (
        ga.Experiment.on("human_igh")
        .recombine()
        .paired_end(r1_length=120, insert_size=240)
        .run_records(n=5, seed=42)
    )
    seen_ids = set()
    for rec in result.records:
        sid = rec.get("sequence_id")
        assert isinstance(sid, str) and sid, (
            f"sequence_id missing or empty: {sid!r}"
        )
        seen_ids.add(sid)
    # Default `from_outcomes` assigns unique sequence_ids
    # (`{prefix}{i}`), so the new exporter's `@{sequence_id}/N`
    # naming produces unambiguous reads.
    assert len(seen_ids) == 5


# ──────────────────────────────────────────────────────────────────
# 5. Scaffold — R2 already RC-oriented at projection time
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_r2_is_already_reverse_complemented() -> None:
    """Audit §3 — the writer must output ``r2_sequence`` verbatim;
    R2 is already the reverse complement of
    ``sequence[r2_start:r2_end]`` at projection time. The
    validator's ``PairedEndWindowMismatch { side: R2 }``
    enforces this. Pin the invariant via
    ``validate_records(refdata)`` running clean on a paired-end
    batch — any validator surface that depended on the R2-is-RC
    invariant would surface here."""
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .paired_end(r1_length=140, insert_size=300)
    )
    refdata = exp.refdata
    result = exp.run_records(n=5, seed=4242)
    report = result.validate_records(refdata)
    # The validator carries `PairedEndWindowMismatch` for both
    # R1 and R2; a clean validation here confirms R2's
    # orientation is what `to_paired_fastq` will rely on.
    assert report, (
        f"paired-end validation failed before FASTQ export audit "
        f"could rely on R2-is-RC invariant: {report.summary()}"
    )


# ──────────────────────────────────────────────────────────────────
# 6. Scaffold — no engine-side per-base quality field
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_no_per_base_quality_field_on_airr_record() -> None:
    """Audit §4 — the engine emits ``n_quality_errors`` (integer
    count) and per-event ``(error_site[i], error_base[i])`` in
    the trace, but does NOT project a per-position Phred score.
    The v1 ``to_paired_fastq`` reuses the existing constant /
    Illumina quality models; a future "real per-base quality"
    model would require an event-payload change that is out of
    scope for this audit's slice."""
    result = (
        ga.Experiment.on("human_igh")
        .recombine()
        .sequencing_errors(count=2)
        .run_records(n=1, seed=4242)
    )
    rec = result.records[0]
    # ``n_quality_errors`` is the only surfaced quality field.
    assert "n_quality_errors" in rec
    # No per-position quality data anywhere on the record dict.
    for forbidden in (
        "quality_string",
        "phred_array",
        "phred_string",
        "quality_array",
        "per_base_quality",
        "r1_quality",
        "r2_quality",
        "r1_phred",
        "r2_phred",
    ):
        assert forbidden not in rec, (
            f"AIRR record now carries {forbidden!r}; the audit §4 "
            "no-per-base-quality precondition regressed — the FASTQ "
            "writer can't safely fall back to the existing pluggable "
            "models without checking which surface to use"
        )


def test_pin_scaffold_quality_pass_records_count_only() -> None:
    """Audit §4 — the ``corrupt.quality`` Rust pass records the
    integer count and per-event ``(site, base)`` pairs in the
    trace, but does NOT emit any per-position Phred score. Pin
    at source so a future regression that adds a Phred field
    surfaces here."""
    src = (
        _REPO_ROOT / "engine_rs" / "src" / "passes" / "corrupt" / "quality.rs"
    ).read_text(encoding="utf-8")
    # The three known trace addresses the pass emits.
    assert "CorruptQualityCount" in src
    assert "CorruptQualitySite" in src
    assert "CorruptQualityBase" in src
    # No per-position Phred variant on the quality pass.
    for forbidden in (
        "CorruptQualityPhred",
        "CorruptQualityScore",
        "QualityPerBase",
        "per_base_phred",
    ):
        assert forbidden not in src, (
            f"quality.rs now carries {forbidden!r}; the audit §4 "
            "no-per-base-quality precondition regressed"
        )


# ──────────────────────────────────────────────────────────────────
# 7. Absence — no SimulationResult.to_paired_fastq method
# ──────────────────────────────────────────────────────────────────


def test_pin_present_to_paired_fastq_method() -> None:
    """Audit §1 / §12 — the implementation slice landed. The
    canonical ``to_paired_fastq(r1_path, r2_path, *,
    quality="illumina", overwrite=False, **quality_kwargs)``
    method exists on ``SimulationResult``. Alternative names
    (``to_interleaved_fastq`` / ``to_fastq_pair`` /
    ``write_paired_fastq``) remain forbidden — only the
    canonical surface is the user-facing API."""
    import inspect

    from GenAIRR.result import SimulationResult

    assert hasattr(SimulationResult, "to_paired_fastq"), (
        "SimulationResult is missing to_paired_fastq; the paired-FASTQ "
        "export slice regressed"
    )
    sig = inspect.signature(SimulationResult.to_paired_fastq)
    # Two required positional paths.
    assert "r1_path" in sig.parameters
    assert "r2_path" in sig.parameters
    # Quality-model surface mirrors single-end ``to_fastq``.
    assert "quality" in sig.parameters
    assert sig.parameters["quality"].default == "illumina"
    # New ``overwrite`` guard (paired-end-specific).
    assert "overwrite" in sig.parameters
    assert sig.parameters["overwrite"].default is False
    # The **quality_kwargs catch-all matches the single-end
    # exporter so users carry the same muscle memory across
    # surfaces.
    var_kw_params = [
        p
        for p in sig.parameters.values()
        if p.kind == inspect.Parameter.VAR_KEYWORD
    ]
    assert len(var_kw_params) == 1

    # Alternative names remain forbidden — keeping the canonical
    # surface a single name.
    for forbidden in (
        "to_interleaved_fastq",
        "to_fastq_pair",
        "write_paired_fastq",
    ):
        assert not hasattr(SimulationResult, forbidden), (
            f"SimulationResult now carries {forbidden!r}; only "
            "``to_paired_fastq`` is the canonical surface"
        )


# ──────────────────────────────────────────────────────────────────
# 8. Absence — no record-generator FASTQ form
# ──────────────────────────────────────────────────────────────────


def test_pin_absence_no_record_generator_fastq_form() -> None:
    """Audit §2 — the writer-only surface is the v1 shape;
    the audit explicitly rejected a record-generator form
    (``to_fastq_records()`` / ``to_paired_fastq_records()``)
    because it would break the existing pattern of the five
    file-write exporters. Pin the absence so a future
    contributor who adds it must update the audit's §2
    rationale in lockstep."""
    from GenAIRR.result import SimulationResult

    for forbidden in (
        "to_fastq_records",
        "to_paired_fastq_records",
        "fastq_lines",
        "iter_fastq",
        "iter_paired_fastq",
    ):
        assert not hasattr(SimulationResult, forbidden), (
            f"SimulationResult now exposes {forbidden!r}; the audit's "
            "writer-only surface decision regressed — update §2 "
            "rationale and flip pin"
        )


# ──────────────────────────────────────────────────────────────────
# 9. Absence — no fastq_export_support block in the manifest
# ──────────────────────────────────────────────────────────────────


def test_pin_absence_no_fastq_export_support_in_manifest() -> None:
    """Audit §11 — FASTQ export is NOT a cartridge-level
    capability; every cartridge can export to FASTQ given that
    ``.paired_end()`` ran. A manifest block would carry zero
    per-cartridge information. Pin the absence so a future
    contributor who adds one must justify it against §11's
    rationale."""
    import json

    m = ga.HUMAN_IGH_OGRDB.cartridge_manifest()
    blob = json.dumps(m)
    for forbidden in (
        "fastq_export_support",
        "fastq_support",
        "read_export_support",
        "paired_fastq_support",
    ):
        assert forbidden not in blob, (
            f"manifest now mentions {forbidden!r}; the audit §11 "
            "no-manifest-block recommendation regressed"
        )


# ──────────────────────────────────────────────────────────────────
# 10. Absence — no new trace addresses for FASTQ export
# ──────────────────────────────────────────────────────────────────


def test_pin_absence_no_fastq_export_trace_addresses() -> None:
    """Audit §10 — FASTQ export is pure projection; it must
    introduce no new ``ChoiceAddress`` variants. Pin at source
    so a refactor that adds a per-record export address (e.g.
    for replay-of-export semantics) surfaces here."""
    address_src = (
        _REPO_ROOT / "engine_rs" / "src" / "address.rs"
    ).read_text(encoding="utf-8")
    for forbidden in (
        "FastqExport",
        "ExportR1",
        "ExportR2",
        "QualityString",
        "PhredArray",
    ):
        assert forbidden not in address_src, (
            f"address.rs now carries {forbidden!r}; the audit §10 "
            "no-new-trace-addresses recommendation regressed"
        )


# ──────────────────────────────────────────────────────────────────
# 11. Doc anchor — audit doc exists and references contract file
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_audit_doc_exists_and_references_contract() -> None:
    """The audit doc must continue to exist and reference this
    contract file; the 15-section structure stays intact."""
    doc = _AUDIT_DOC.read_text(encoding="utf-8")
    assert "test_fastq_export_contract.py" in doc, (
        "audit doc no longer references the contract file; lockstep "
        "convention drifted"
    )
    for marker in (
        "## 1. Q1",
        "## 2. Q2",
        "## 3. Q3",
        "## 12. Implementation order",
        "## 15. Summary table",
    ):
        assert marker in doc, f"audit doc missing section marker {marker!r}"
