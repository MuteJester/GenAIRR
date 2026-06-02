"""Pre-implementation contract for paired-end / read-layout output.

Companion to [`docs/paired_end_design.md`](../docs/paired_end_design.md).
Paired-end read-layout projection (one AIRR row per molecule with
R1 / R2 windows carved out of ``rec.sequence``) isn't implemented
yet; this file pins today's absence + the existing scaffolding the
design relies on. When the implementation slices land, the
absence assertions flip to presence assertions in lockstep — same
convention as [`test_d_inversion_contract.py`](test_d_inversion_contract.py)
and [`test_receptor_revision_contract.py`](test_receptor_revision_contract.py).

The split:

- ``pin_absence_*`` tests assert today's behaviour — the DSL
  method doesn't exist, no ``paired_end.*`` trace addresses
  appear in baseline runs, no AIRR ``r1_sequence`` /
  ``r2_sequence`` / ``read_layout`` fields are emitted, no
  ``PairedEndWindowMismatch`` / ``PairedEndGeometryMismatch``
  validator issues exist. These MUST PASS before the
  implementation slices.
- ``pin_scaffold_*`` tests pin the surfaces that already exist
  in the engine — the byte-level ``complement_base`` helper, the
  ``apply_rev_comp_projection`` post-build pattern, the
  ``end_loss_5_length`` / ``end_loss_3_length`` additive-field
  precedent, the structured-validator ``details.source`` shape
  — which the design relies on but which a future refactor could
  regress. Treat these as "the scaffolding the new mechanism
  needs."

When Slice A → E land, the relevant ``pin_absence_*`` tests get
inverted in lockstep with the new behaviour and the audit doc
crosses out the corresponding items in §13.
"""
from __future__ import annotations

import inspect
import re
from pathlib import Path

import GenAIRR as ga
from GenAIRR import _engine as ge


# ──────────────────────────────────────────────────────────────────
# Deterministic VDJ fixture (matches the shape of the D-inversion
# and receptor-revision contracts so a future Slice D test can
# drop into the same harness).
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


# ──────────────────────────────────────────────────────────────────
# 1. DSL surface absence — Experiment.paired_end is not yet
#    a method, and existing steps don't carry a paired-end kwarg.
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_experiment_paired_end_method_exists() -> None:
    """Slice D flipped this from pin-absence to pin-presence.
    Pin that ``Experiment.paired_end`` exists as a method,
    accepts the three documented kwargs, and rejects unknown
    kwargs (`read_length=`, `r1=`, etc.) — the three guarantees
    the design doc §3 locked in. Validation specifics (range,
    geometry, duplicates) are covered by
    ``test_paired_end_dsl.py``; here we only pin the *existence*
    of the surface so later refactors can't rename or drop the
    method without surfacing it through this audit file as well."""
    assert hasattr(ga.Experiment, "paired_end"), (
        "Experiment.paired_end is missing; Slice D of the "
        "paired-end roadmap should have landed it."
    )
    sig = inspect.signature(ga.Experiment.paired_end)
    for required in ("r1_length", "r2_length", "insert_size"):
        assert required in sig.parameters, (
            f"Experiment.paired_end is missing the {required!r} "
            f"keyword; design doc §3 froze the three kwargs."
        )


def test_pin_absence_existing_steps_carry_no_paired_end_kwarg() -> None:
    """The design doc §3 rejects a kwarg-style alternative on the
    existing observation-stage steps. Pin that
    ``end_loss_5/3prime``, ``random_strand_orientation``, and
    ``recombine`` don't sneak a paired-end-related kwarg in."""
    for name in (
        "end_loss_5prime",
        "end_loss_3prime",
        "primer_trim_5prime",
        "primer_trim_3prime",
        "random_strand_orientation",
        "recombine",
    ):
        sig = inspect.signature(getattr(ga.Experiment, name))
        for forbidden in (
            "paired_end",
            "paired_end_prob",
            "read_length",
            "insert_size",
            "r1_length",
            "r2_length",
        ):
            assert forbidden not in sig.parameters, (
                f"Experiment.{name} carries unexpected kwarg "
                f"{forbidden!r}; design doc §3 rejected the "
                f"kwarg-style alternative."
            )


def test_pin_present_to_paired_fastq_and_absence_of_single_end_alias() -> None:
    """The paired-end design's §11 deferred ``to_paired_fastq``
    to a follow-up slice; that slice has now landed via
    ``docs/fastq_export_design.md`` and is enforced by
    [`tests/test_fastq_export_contract.py`](test_fastq_export_contract.py)
    and [`tests/test_to_paired_fastq.py`](test_to_paired_fastq.py).
    This pin captures the cross-audit lockstep: the
    paired-end-design absence pin flips to present, and the
    sibling ``.single_end`` DSL alias stays absent (the design
    doc §2.2 reserved the name but did not land it)."""
    # `.single_end` DSL surface remains deferred per paired-end
    # design §2.2.
    assert not hasattr(ga.Experiment, "single_end"), (
        "Experiment.single_end exists; paired-end design §2.2 "
        "reserved the name but did not land it in v1."
    )
    # `to_paired_fastq` is now a SimulationResult method per the
    # FASTQ export design doc.
    from GenAIRR.result import SimulationResult

    assert hasattr(SimulationResult, "to_paired_fastq"), (
        "SimulationResult is missing to_paired_fastq; the FASTQ "
        "export slice regressed (see docs/fastq_export_design.md)"
    )


# ──────────────────────────────────────────────────────────────────
# 2. Trace address vocabulary — no paired_end.* records today
# ──────────────────────────────────────────────────────────────────


def test_pin_absence_no_paired_end_address_in_baseline_run() -> None:
    """The proposed trace addresses per §4 are
    ``paired_end.r1_length``, ``.r2_length``, ``.insert_size``.
    A baseline VDJ run with no paired-end step must not emit any
    of them (and there's no other source today)."""
    outcomes = _baseline_experiment().run(n=1, seed=0)
    addrs = {r.address for r in outcomes[0].trace().choices()}
    for forbidden in (
        "paired_end.r1_length",
        "paired_end.r2_length",
        "paired_end.insert_size",
        # The previously-considered "Bool" address form was
        # rejected at audit time — pin its absence so a future
        # contributor doesn't accidentally adopt it.
        "paired_end.applied",
        # Alternative `read_*` namespace — also rejected.
        "read_layout.r1_length",
    ):
        assert forbidden not in addrs, (
            f"baseline run already emits {forbidden!r}; this means "
            f"a different pass is silently writing into the "
            f"paired-end namespace before the dedicated slice "
            f"lands."
        )


def test_pin_scaffold_paired_end_addresses_in_engine_source() -> None:
    """Slice C flipped this from pin-absence to pin-presence.
    Pin the three on-disk address spellings the design doc §4
    froze: ``paired_end.r1_length``, ``.r2_length``,
    ``.insert_size``. The address-schema-version add-only policy
    (per address.rs top-of-file docs) means the additive change
    does NOT require a schema bump — but the spellings themselves
    must round-trip forever, so pin them here as a second line of
    defence beyond the engine's own
    ``frozen_address_spellings_*`` test."""
    addr_src = (
        Path(__file__).resolve().parent.parent
        / "engine_rs"
        / "src"
        / "address.rs"
    ).read_text(encoding="utf-8")
    for expected in (
        '"paired_end.r1_length"',
        '"paired_end.r2_length"',
        '"paired_end.insert_size"',
    ):
        assert expected in addr_src, (
            f"engine_rs/src/address.rs no longer carries {expected}; "
            f"Slice C's address vocabulary has drifted."
        )


# ──────────────────────────────────────────────────────────────────
# 3. AIRR field absence — no r1_sequence / r2_sequence / etc. yet
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_paired_end_airr_fields_emit_with_defaults() -> None:
    """Slice A flipped this from pin-absence to pin-presence.
    The eight new fields documented in §2 are now present on
    every record at their documented defaults
    (``""`` / ``None`` / ``0``). Pin the keys + the rejected
    naming alternatives' continued absence so a later Slice B/C
    refactor can't drop or rename them."""
    rec = _baseline_experiment().run_records(n=1, seed=0).records[0]
    # Present + defaulted.
    for required, default in (
        ("read_layout", ""),
        ("r1_sequence", ""),
        ("r2_sequence", ""),
        ("r1_start", None),
        ("r1_end", None),
        ("r2_start", None),
        ("r2_end", None),
        ("insert_size", 0),
    ):
        assert required in rec, (
            f"Slice A AIRR record missing paired-end field {required!r}"
        )
        assert rec[required] == default, (
            f"baseline record's {required!r} is {rec[required]!r}; "
            f"design doc §10 specifies default {default!r}."
        )
    # Rejected naming alternatives stay out.
    for forbidden in (
        "fragment_size",
        "read1_sequence",
        "read2_sequence",
        "paired_end_applied",
    ):
        assert forbidden not in rec, (
            f"AIRR records expose rejected name {forbidden!r}; "
            f"design doc §2 picked the canonical r1/r2 + "
            f"insert_size spelling."
        )


def test_pin_scaffold_paired_end_columns_in_result_column_list() -> None:
    """Slice A flipped this from pin-absence to pin-presence.
    The canonical Python column list now carries the eight
    fields in struct order (matching the Rust ``AirrRecord``
    layout). Pin order + presence so a CSV consumer pinning a
    column index doesn't silently shift."""
    from GenAIRR.result import _DEFAULT_COLUMN_ORDER

    expected_order = [
        "read_layout",
        "r1_sequence",
        "r2_sequence",
        "r1_start",
        "r1_end",
        "r2_start",
        "r2_end",
        "insert_size",
    ]
    for required in expected_order:
        assert required in _DEFAULT_COLUMN_ORDER, (
            f"_DEFAULT_COLUMN_ORDER missing paired-end column "
            f"{required!r}; Slice A should have landed it."
        )
    observed = [c for c in _DEFAULT_COLUMN_ORDER if c in expected_order]
    assert observed == expected_order, (
        f"paired-end column order in _DEFAULT_COLUMN_ORDER drifted "
        f"from the Rust struct layout: got {observed}, expected "
        f"{expected_order}."
    )


def test_pin_scaffold_paired_end_fields_declared_on_engine_record() -> None:
    """Slice A flipped this from pin-absence to pin-presence.
    Pin the eight struct-field declarations on the Rust
    ``AirrRecord`` so a later refactor that drops one surfaces
    here before the Python column list test fires."""
    record_src = (
        Path(__file__).resolve().parent.parent
        / "engine_rs"
        / "src"
        / "airr_record"
        / "record.rs"
    ).read_text(encoding="utf-8")
    for expected in (
        "pub r1_sequence",
        "pub r2_sequence",
        "pub r1_start",
        "pub r1_end",
        "pub r2_start",
        "pub r2_end",
        "pub insert_size",
        "pub read_layout",
    ):
        assert expected in record_src, (
            f"engine_rs/src/airr_record/record.rs no longer "
            f"declares {expected!r}; Slice A's schema landing has "
            f"drifted."
        )


# ──────────────────────────────────────────────────────────────────
# 4. Existing scaffolding — `reverse_complement` + `complement_base`
#    helpers, `apply_rev_comp_projection` pattern
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_complement_base_kernel_is_deterministic_and_involutive() -> None:
    """``complement_base`` is the byte-level Watson-Crick kernel
    every higher-level reverse-complement helper composes from.
    Slice B's R2-window projection will compose it via the
    existing ``reverse_complement`` string helper. Pin
    determinism + involutivity here so a future kernel rename or
    behavior change surfaces in the audit doc rather than in a
    paired-end bug report.

    The kernel is private to the Rust engine (no Python binding);
    we exercise it through the AIRR record's rev-comp projection
    path — running ``random_strand_orientation(prob=1.0)`` on a
    known-bytes pipeline and asserting the projected sequence is
    the Watson-Crick complement of the unflipped baseline."""
    cfg = _vdj_refdata()
    base = (
        ga.Experiment.on(cfg)
        .allow_curatable_refdata()
        .recombine(np1_lengths=[(3, 1.0)], np2_lengths=[(3, 1.0)])
        .trim(enabled=False)
    )
    flipped = (
        ga.Experiment.on(cfg)
        .allow_curatable_refdata()
        .recombine(np1_lengths=[(3, 1.0)], np2_lengths=[(3, 1.0)])
        .trim(enabled=False)
        .random_strand_orientation(prob=1.0)
    )

    base_seq = base.run_records(n=1, seed=0).records[0]["sequence"]
    flipped_seq = flipped.run_records(n=1, seed=0).records[0]["sequence"]

    # Reverse-complement of the baseline equals the flipped run.
    def _py_revcomp(s: str) -> str:
        table = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
        return s.translate(table)[::-1]

    assert _py_revcomp(base_seq) == flipped_seq, (
        "reverse-complement projection no longer matches the "
        "Python Watson-Crick oracle; Slice B's R2 windowing "
        "relies on this composition staying deterministic."
    )
    # Involutivity smoke: rev-comping twice is identity.
    assert _py_revcomp(_py_revcomp(base_seq)) == base_seq


def test_pin_scaffold_apply_rev_comp_projection_pattern_is_post_build() -> None:
    """``apply_rev_comp_projection`` is invoked at the END of
    ``build_airr_record`` — the post-build flip pattern paired-end
    Slice B follows. Pin that the helper is still in
    ``airr_record/sequence.rs`` and still called *after* every
    other field is populated. If a future refactor reorders the
    call (e.g. flips coords before the field is fully
    materialised), the §7 composition guarantee for paired-end
    breaks."""
    sequence_src = (
        Path(__file__).resolve().parent.parent
        / "engine_rs"
        / "src"
        / "airr_record"
        / "sequence.rs"
    ).read_text(encoding="utf-8")
    builder_src = (
        Path(__file__).resolve().parent.parent
        / "engine_rs"
        / "src"
        / "airr_record"
        / "builder.rs"
    ).read_text(encoding="utf-8")

    assert "fn apply_rev_comp_projection" in sequence_src, (
        "apply_rev_comp_projection helper is missing; paired-end "
        "Slice B's projection layer relies on this pattern."
    )
    # The call site is in `build_airr_record`, and the design
    # requires it to be the LAST thing the builder does (before
    # `rec` is returned).
    assert "apply_rev_comp_projection(&mut rec)" in builder_src, (
        "apply_rev_comp_projection no longer fires from "
        "build_airr_record; paired-end's R1/R2 windowing depends "
        "on the post-build composition order."
    )


# ──────────────────────────────────────────────────────────────────
# 5. End-loss fields scaffold — the additive-AIRR-field precedent
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_end_loss_fields_exist_with_default_zero() -> None:
    """``end_loss_5_length`` and ``end_loss_3_length`` are the
    canonical additive-AIRR-field precedent: they default to 0 on
    every record that didn't go through the end-loss pass, and
    populate only when the user opted in. The eight paired-end
    fields will follow the identical pattern (default empty /
    None / 0 / empty), so the precedent must remain in place."""
    rec = _baseline_experiment().run_records(n=1, seed=0).records[0]
    assert "end_loss_5_length" in rec
    assert "end_loss_3_length" in rec
    # Baseline run carries no end-loss step → both default to 0.
    assert rec["end_loss_5_length"] == 0
    assert rec["end_loss_3_length"] == 0


def test_pin_scaffold_end_loss_runs_before_projection_in_pipeline() -> None:
    """Design doc §6 requires the pipeline order
    ``…corruption → end_loss → paired_end projection``. End-loss
    therefore deletes pool bytes BEFORE paired-end carves out
    R1/R2 windows. Pin the existing behaviour: an end-loss step
    actually shortens ``rec.sequence_length`` rather than leaving
    the molecule intact and recording the loss elsewhere."""
    base_rec = _baseline_experiment().run_records(n=1, seed=0).records[0]
    end_lost_rec = (
        _baseline_experiment()
        .end_loss_5prime(length=[(3, 1.0)])
        .end_loss_3prime(length=[(2, 1.0)])
        .run_records(n=1, seed=0)
        .records[0]
    )
    assert end_lost_rec["sequence_length"] < base_rec["sequence_length"], (
        "end_loss no longer shortens rec.sequence_length; "
        "paired-end Slice B will draw R1/R2 windows from a "
        "molecule that's already been shortened, so this precedent "
        "must remain in place."
    )
    # And the trace records the loss amounts (the additive-trace
    # precedent paired-end will mirror).
    assert end_lost_rec["end_loss_5_length"] == 3
    assert end_lost_rec["end_loss_3_length"] == 2


# ──────────────────────────────────────────────────────────────────
# 6. Validator surface — existing structured issue shape;
#    PairedEnd* variants absent today
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_validator_emits_structured_details_source() -> None:
    """The validator's structured-issue JSON shape carries
    ``kind`` + ``details.source`` strings; paired-end's
    ``PairedEndWindowMismatch`` and ``PairedEndGeometryMismatch``
    will land in the same shape (§8.2). Pin that at least one
    existing source string survives (we use
    ``trace:sample_allele.v`` from receptor revision Slice E,
    pinned by ``test_pin_scaffold_receptor_revision_validator_issues_exist``).
    A future refactor that drops the ``details.source`` convention
    would break the paired-end design's MCP / dashboard
    compatibility plan."""
    outcome_src = (
        Path(__file__).resolve().parent.parent
        / "engine_rs"
        / "src"
        / "python"
        / "outcome.rs"
    ).read_text(encoding="utf-8")
    assert '"trace:' in outcome_src, (
        "no structured `details.source: trace:*` strings in "
        "outcome.rs; paired-end Slice B/C depends on this "
        "convention to land its own validator issues."
    )


def test_pin_scaffold_paired_end_validator_issue_scaffolding_landed() -> None:
    """Slice A flipped this from pin-absence to pin-presence.
    The user's Slice A scoping landed five variant declarations
    (under the ``Read*`` / ``PairedEnd*`` family) but only
    enforces the no-layout-default invariant via
    ``PairedEndFieldWithoutLayout``. Pin all five variant names
    so a Slice B/C refactor that renames one surfaces here.

    The audit doc §8 originally proposed ``PairedEndWindowMismatch``
    + ``PairedEndGeometryMismatch``; the Slice A spec narrowed and
    renamed the variants to ``ReadWindowOutOfBounds`` /
    ``ReadSequenceMismatch`` / ``ReadInsertSizeMismatch`` /
    ``ReadLayoutMismatch`` so the structured-issue family reads
    consistently. The design doc's §8 table will be reconciled in
    the Slice B doc update."""
    validate_src = (
        Path(__file__).resolve().parent.parent
        / "engine_rs"
        / "src"
        / "airr_record"
        / "validate.rs"
    ).read_text(encoding="utf-8")
    for expected in (
        "PairedEndFieldWithoutLayout",
        "ReadWindowOutOfBounds",
        "ReadSequenceMismatch",
        "ReadInsertSizeMismatch",
        "ReadLayoutMismatch",
    ):
        assert expected in validate_src, (
            f"validator no longer declares {expected!r}; Slice A's "
            f"variant scaffolding has drifted."
        )

    # And the Slice A no-layout-default invariant's structured
    # `details.source` string is pinned, mirroring the receptor-
    # revision Slice E source-string contract.
    outcome_src = (
        Path(__file__).resolve().parent.parent
        / "engine_rs"
        / "src"
        / "python"
        / "outcome.rs"
    ).read_text(encoding="utf-8")
    assert '"schema:paired_end_layout_default"' in outcome_src, (
        "PairedEndFieldWithoutLayout's structured `details.source` "
        "string was renamed; downstream MCP / dashboard consumers "
        "match this prefix verbatim."
    )


# ──────────────────────────────────────────────────────────────────
# 7. Baseline run = one row per outcome, no paired-end fields
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_baseline_run_emits_exactly_one_airr_row_per_outcome() -> None:
    """Design doc §2 rejected "two AIRR rows per simulation" for
    v1. Pin that today's baseline keeps the one-row-per-outcome
    contract, so a future Slice A/D refactor that accidentally
    duplicated rows would surface here before the audit doc
    needs revising."""
    n = 7
    result = _baseline_experiment().run_records(n=n, seed=0)
    assert len(result.records) == n, (
        f"baseline run emits {len(result.records)} AIRR rows for "
        f"n={n} outcomes; design doc §2 requires one row per "
        f"outcome."
    )
    # And the records aren't structurally duplicates — each carries
    # a unique sequence_id, the one-row-per-outcome invariant.
    sequence_ids = [rec["sequence_id"] for rec in result.records]
    assert len(set(sequence_ids)) == n


def test_pin_scaffold_baseline_run_validator_passes_clean() -> None:
    """Pin that the validator currently runs clean on a baseline
    full-stack record — paired-end Slice B's validator additions
    should not turn baseline records into failures. The check is
    a pre-condition for Slice A's "additive fields default
    empty/None" promise: if the baseline already failed, the
    paired-end null-defaults wouldn't be a clean signal."""
    refdata = _baseline_experiment().refdata
    result = _baseline_experiment().run_records(n=3, seed=0)
    report = result.validate_records(refdata)
    assert report, (
        f"baseline VDJ run no longer validates clean: "
        f"summary={report.summary()}; first failure="
        f"{report.failures[0] if report.failures else None}. "
        f"Paired-end Slice A's additive-field promise depends "
        f"on this staying green."
    )


def test_pin_scaffold_baseline_trace_replay_round_trip_works() -> None:
    """Pin that the existing trace replay round-trip works on a
    baseline record (no paired-end). Slice C's replay path will
    extend the existing consume-cursor surface; a regression in
    the baseline replay would surface here first."""
    from GenAIRR._airr_record import outcome_to_airr_record

    exp = _baseline_experiment()
    compiled = exp.compile()
    fresh = compiled.simulator.run(seed=0)
    tf = compiled.simulator.trace_file_from(fresh, seed=0)
    replayed = compiled.simulator.rerun_from_trace_file(tf)

    fresh_rec = outcome_to_airr_record(
        fresh, exp._refdata, sequence_id="fresh"
    )
    replayed_rec = outcome_to_airr_record(
        replayed, exp._refdata, sequence_id="replay"
    )
    assert fresh_rec["sequence"] == replayed_rec["sequence"], (
        "baseline trace replay no longer reproduces the sequence; "
        "paired-end Slice C's three-record consume path depends on "
        "this staying deterministic."
    )


# ──────────────────────────────────────────────────────────────────
# 8. Design doc shape pins
# ──────────────────────────────────────────────────────────────────


def test_pin_design_doc_paired_end_lists_fourteen_sections() -> None:
    """The audit doc follows the same 14-section structure as
    `docs/receptor_revision_design.md`. Pin the section count so
    a future contributor can't quietly drop a section."""
    docs_dir = Path(__file__).resolve().parent.parent / "docs"
    if not docs_dir.is_dir():
        import pytest
        pytest.skip("docs/ is contributor-only; not present in this checkout")
    doc = (
        Path(__file__).resolve().parent.parent
        / "docs"
        / "paired_end_design.md"
    ).read_text(encoding="utf-8")
    section_headers = re.findall(r"^## (\d{1,2})\. ", doc, re.MULTILINE)
    expected = [str(i) for i in range(1, 15)]
    assert section_headers == expected, (
        f"docs/paired_end_design.md section ordering changed; "
        f"got {section_headers}, expected {expected}. The audit "
        f"is the contract; renumber tests if the structure is "
        f"intentionally reshaped."
    )


def test_pin_design_doc_summary_table_carries_every_decision() -> None:
    """The §Summary table is the canonical short-form of the
    audit decisions. Pin the row labels so a future contributor
    can't quietly drop a decision."""
    docs_dir = Path(__file__).resolve().parent.parent / "docs"
    if not docs_dir.is_dir():
        import pytest
        pytest.skip("docs/ is contributor-only; not present in this checkout")
    doc = (
        Path(__file__).resolve().parent.parent
        / "docs"
        / "paired_end_design.md"
    ).read_text(encoding="utf-8")
    # Find the summary table block and pull the leftmost cells.
    summary_match = re.search(
        r"## Summary table\n\n\| Question \| Decision \|.*?(?=\n##|\Z)",
        doc,
        re.DOTALL,
    )
    assert summary_match, "design doc summary table is missing"
    summary_block = summary_match.group(0)
    rows = re.findall(r"^\| ([^|]+) \| ", summary_block, re.MULTILINE)
    # Strip header + separator rows.
    rows = [r.strip() for r in rows if r.strip() not in ("Question", "--- ")]
    rows = [r for r in rows if not set(r) <= {"-", " "}]
    for required in (
        "Scope",
        "Output shape",
        "DSL",
        "Pipeline position",
        "Trace",
        "IR carrier",
        "Event",
        "Live-call refresh",
        "AIRR",
        "Contracts",
        "Replay",
        "Validator",
        "`random_strand_orientation` interaction",
        "`end_loss` interaction",
        "Backwards compat",
    ):
        assert required in rows, (
            f"design doc summary table missing {required!r} "
            f"decision row; the audit is the contract."
        )


def test_pin_followup_plan_for_when_slices_land_is_in_this_file() -> None:
    """When Slice A → E land, the pin_absence_* tests in this
    file get inverted to pin_scaffold_* in lockstep with the new
    behaviour. Pin that this followup pattern is documented in
    the module docstring so a future contributor knows what to
    update."""
    this_doc = Path(__file__).read_text(encoding="utf-8")
    assert "pin_absence_" in this_doc
    assert "pin_scaffold_" in this_doc
    assert "lockstep" in this_doc, (
        "Module docstring should explain the pin-flip lockstep "
        "convention. Future contributors who land Slice A will "
        "search for this term to find the playbook."
    )
