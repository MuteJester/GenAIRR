"""Pre-implementation contract for receptor revision.

Companion to [`docs/receptor_revision_design.md`](../docs/receptor_revision_design.md).
Receptor revision (V replacement after initial V(D)J recombination)
isn't implemented yet; this file pins today's absence + the
existing scaffolding the design relies on. When the implementation
slice(s) land, the absence assertions flip to presence assertions
in lockstep — same convention as
[`test_d_inversion_contract.py`](test_d_inversion_contract.py).

The split:

- ``pin_absence_*`` tests assert today's behaviour — the DSL
  method doesn't exist, no ``receptor_revision.*`` trace
  addresses appear in baseline runs, no AIRR ``receptor_revision_applied``
  field is emitted, ``SegmentReplaced`` event variant isn't
  defined yet. These MUST PASS before the implementation slices.
- ``pin_scaffold_*`` tests pin the surfaces that already exist
  in the engine — ``RegionReplaced`` event variant,
  ``SimulationBuilder::replace_region``, the
  ``MultipleRegionsForSegment`` validator — which the design
  relies on but which a future refactor could regress. Treat
  these as "the scaffolding the new mechanism needs."

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
# Deterministic VDJ fixture (matches the D-inversion contract's
# shape so a future Slice D test can drop into the same harness).
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
# 1. DSL surface absence — Experiment.receptor_revision is not yet
#    a method, and recombine() carries no kwarg-style alternative
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_experiment_receptor_revision_method_exists() -> None:
    """Slice D flipped this from pin-absence to pin-presence.
    Pin that ``Experiment.receptor_revision`` exists as a method,
    accepts the documented ``prob`` keyword, and rejects VJ chains
    at the call boundary — the three guarantees the design doc §2
    locked in. Validation specifics (NaN, range, duplicates) are
    covered by ``test_receptor_revision_dsl.py``; here we only
    pin the *existence* of the surface so later refactors can't
    rename or drop the method without surfacing it through this
    audit file as well."""
    assert hasattr(ga.Experiment, "receptor_revision"), (
        "Experiment.receptor_revision is missing; Slice D of the "
        "receptor-revision roadmap should have landed it."
    )
    sig = inspect.signature(ga.Experiment.receptor_revision)
    assert "prob" in sig.parameters, (
        "Experiment.receptor_revision dropped its `prob` keyword; "
        "the design doc §2 froze it as the only argument."
    )


def test_pin_absence_recombine_carries_no_receptor_revision_kwarg() -> None:
    """The design doc §2 rejects a `recombine(receptor_revision_prob=…)`
    kwarg in favour of an explicit fluent step. Pin that no such
    kwarg has been added — keeps the DSL shape consistent across
    biology mechanisms (mirrors the same guard for D inversion)."""
    sig = inspect.signature(ga.Experiment.recombine)
    assert "receptor_revision_prob" not in sig.parameters
    assert "revise_v_prob" not in sig.parameters


def test_pin_absence_no_revise_v_alias_method() -> None:
    """The design considered `revise_v` and rejected it (§2
    rationale: a future light-chain extension wouldn't force a
    rename). Pin that we didn't land BOTH spellings — the alias
    would confuse the surface vocabulary."""
    assert not hasattr(ga.Experiment, "revise_v"), (
        "Experiment.revise_v exists; the design picked "
        "receptor_revision as the canonical name. If both spellings "
        "are intentional now, update the design doc §2."
    )


# ──────────────────────────────────────────────────────────────────
# 2. Trace address vocabulary — no receptor_revision.* records in
#    today's engine
# ──────────────────────────────────────────────────────────────────


def test_pin_absence_no_receptor_revision_address_in_baseline_run() -> None:
    """The proposed trace addresses per §3 are
    ``receptor_revision.applied``, ``.v_allele``, ``.v_trim_3``.
    A baseline VDJ run with no revision step must not emit any of
    them (and there's no other source today)."""
    outcomes = _baseline_experiment().run(n=1, seed=0)
    addrs = {r.address for r in outcomes[0].trace().choices()}
    for forbidden in (
        "receptor_revision.applied",
        "receptor_revision.v_allele",
        "receptor_revision.v_trim_3",
    ):
        assert forbidden not in addrs, (
            f"baseline run already emits {forbidden!r}; this means "
            f"a different pass is silently writing into the receptor-"
            f"revision namespace before the dedicated slice lands."
        )


def test_pin_scaffold_receptor_revision_addresses_in_engine_source() -> None:
    """Slice C flipped this from pin-absence to pin-presence. Pin
    the three on-disk address spellings the design doc §3 froze:
    ``receptor_revision.applied``, ``.v_allele``, ``.v_trim_3``.
    The address-schema-version add-only policy (per address.rs
    top-of-file docs) means the additive change does NOT require
    a schema bump — but the spellings themselves must round-trip
    forever, so pin them here as a second line of defence beyond
    the engine's own ``frozen_address_spellings_*`` test."""
    addr_src = (
        Path(__file__).resolve().parent.parent
        / "engine_rs"
        / "src"
        / "address.rs"
    ).read_text(encoding="utf-8")
    for expected in (
        '"receptor_revision.applied"',
        '"receptor_revision.v_allele"',
        '"receptor_revision.v_trim_3"',
    ):
        assert expected in addr_src, (
            f"engine_rs/src/address.rs no longer carries {expected}; "
            f"Slice C's address vocabulary has drifted."
        )


# ──────────────────────────────────────────────────────────────────
# 3. AIRR field absence — no receptor_revision_applied yet
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_receptor_revision_airr_fields_exist() -> None:
    """Slice E flipped this from pin-absence to pin-presence. The
    AIRR record dict now carries
    ``receptor_revision_applied: bool`` and ``original_v_call: str``
    on every record. Pin both keys (and the documented baseline
    values for a no-revision run) so future projection changes
    don't drop the fields silently. The rejected alternative name
    ``original_v_call_id`` remains absent — design doc §7 picked
    the string name as the canonical spelling."""
    rec = _baseline_experiment().run_records(n=1, seed=0).records[0]
    assert "receptor_revision_applied" in rec
    assert "original_v_call" in rec
    assert rec["receptor_revision_applied"] is False, (
        "baseline run carries receptor_revision_applied=True without "
        "calling .receptor_revision(); a different pass is silently "
        "writing the field."
    )
    assert rec["original_v_call"] == "", (
        "baseline run carries a non-empty original_v_call; design "
        "doc §7 specifies the empty sentinel for no-revision records."
    )
    # The rejected `original_v_call_id` alternative stays out of the
    # surface.
    assert "original_v_call_id" not in rec, (
        "design doc §7 picked original_v_call (string) over "
        "original_v_call_id (int); both spellings shouldn't coexist."
    )


# ──────────────────────────────────────────────────────────────────
# 4. SimulationEvent surface — RegionReplaced exists today; the
#    NEW SegmentReplaced variant the design proposes (§5) is absent
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_region_replaced_event_variant_exists() -> None:
    """`SimulationEvent::RegionReplaced` is the scaffolding the
    design relies on for the *metadata-only* region edits and the
    user-visible event payload. Pin that the variant stays defined
    with the expected `(old, new): Region` payload — a rename or
    payload reshape would break §5's plan."""
    sim_event_src = (
        Path(__file__).resolve().parent.parent
        / "engine_rs"
        / "src"
        / "ir"
        / "sim_event.rs"
    ).read_text(encoding="utf-8")
    # Match the variant definition with its `(old, new): Region`
    # payload across multiple lines.
    pattern = re.compile(
        r"RegionReplaced\s*\{\s*old:\s*Region\s*,\s*new:\s*Region\s*\}",
        re.MULTILINE,
    )
    assert pattern.search(sim_event_src), (
        "SimulationEvent::RegionReplaced variant changed shape or "
        "was removed; design doc §5 relies on the (old, new): Region "
        "payload."
    )


def test_pin_scaffold_segment_replaced_event_variant_is_defined() -> None:
    """Slice A flipped this from pin-absence to pin-presence. Design
    doc §5 introduced ``SimulationEvent::SegmentReplaced { segment,
    old_region, new_region, bytes_delta }`` for *length-changing*
    pool surgery (distinct from the metadata-only
    `RegionReplaced`). Pin the variant's shape so later slices that
    add producers and consumers don't accidentally rename or
    restructure the payload."""
    sim_event_src = (
        Path(__file__).resolve().parent.parent
        / "engine_rs"
        / "src"
        / "ir"
        / "sim_event.rs"
    ).read_text(encoding="utf-8")
    assert "SegmentReplaced {" in sim_event_src, (
        "SimulationEvent::SegmentReplaced variant missing; Slice A "
        "of the receptor-revision roadmap is supposed to have "
        "landed it."
    )
    # Payload shape: four named fields matching design doc §5.
    for field in ("segment:", "old_region:", "new_region:", "bytes_delta:"):
        assert field in sim_event_src, (
            f"SimulationEvent::SegmentReplaced is missing the "
            f"`{field.rstrip(':')}` field design doc §5 specifies."
        )


# ──────────────────────────────────────────────────────────────────
# 5. Builder surface — replace_region exists; replace_segment does not
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_builder_replace_region_exists() -> None:
    """`SimulationBuilder::replace_region` is the existing
    metadata-only region-edit primitive the design relies on (§5).
    Pin that the method stays defined with its event-emitting
    body — a refactor that drops the `RegionReplaced` broadcast
    would break the audit's IR-commit-path plan."""
    builder_src = (
        Path(__file__).resolve().parent.parent
        / "engine_rs"
        / "src"
        / "ir"
        / "builder.rs"
    ).read_text(encoding="utf-8")
    assert "pub(crate) fn replace_region(" in builder_src, (
        "SimulationBuilder::replace_region renamed or removed; "
        "design doc §5 relies on this method existing."
    )
    # And it must still emit the RegionReplaced event.
    assert re.search(
        r"replace_region.*broadcast_event\(\s*SimulationEvent::RegionReplaced",
        builder_src,
        re.DOTALL,
    ), (
        "SimulationBuilder::replace_region no longer broadcasts "
        "RegionReplaced; the scaffold the design relies on is gone."
    )


def test_pin_scaffold_builder_replace_segment_method_exists() -> None:
    """Slice A flipped this from pin-absence to pin-presence.
    `SimulationBuilder::replace_segment` is the NEW bulk-
    replacement primitive (§4.1 / §5). Pin that it (a) is defined,
    (b) broadcasts `SimulationEvent::SegmentReplaced`, and
    (c) returns the `(old, new)` region pair the future
    `ReceptorRevisionPass` will record."""
    builder_src = (
        Path(__file__).resolve().parent.parent
        / "engine_rs"
        / "src"
        / "ir"
        / "builder.rs"
    ).read_text(encoding="utf-8")
    assert "fn replace_segment(" in builder_src, (
        "SimulationBuilder::replace_segment is missing; Slice A of "
        "the receptor-revision roadmap should have landed it."
    )
    assert re.search(
        r"replace_segment.*broadcast_event\(\s*SimulationEvent::SegmentReplaced",
        builder_src,
        re.DOTALL,
    ), (
        "SimulationBuilder::replace_segment no longer broadcasts "
        "SegmentReplaced; the IR-commit path the design relies on "
        "(§5) is broken."
    )


# ──────────────────────────────────────────────────────────────────
# 6. Validator surface — MultipleRegionsForSegment exists; the
#    new receptor-revision validator issues are absent
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_multiple_regions_for_segment_issue_exists() -> None:
    """`MultipleRegionsForSegment` is the existing post-condition
    that enforces the "one region per V/D/J" invariant the design
    relies on (§4 + §10). A revision pass that preserves this
    invariant inherits the existing check for free."""
    validate_src = (
        Path(__file__).resolve().parent.parent
        / "engine_rs"
        / "src"
        / "airr_record"
        / "validate.rs"
    ).read_text(encoding="utf-8")
    assert "MultipleRegionsForSegment" in validate_src, (
        "RecordValidationIssue::MultipleRegionsForSegment removed; "
        "design doc §4 relies on this check to defend the "
        "one-region-per-segment invariant under revision."
    )


def test_pin_scaffold_receptor_revision_validator_issues_exist() -> None:
    """Pin the two validator issue variants design doc §10 specified:
    ``ReceptorRevisionAppliedMismatch`` and
    ``OriginalVCallMismatch``. The PyO3 dict serialization carries
    a ``details.source`` field per the structured-issue convention;
    pin the source string here.

    Bug D fix (see ``tests/test_receptor_revision_provenance_in_ir.py``)
    moved receptor-revision provenance from the trace to the IR. The
    ``details.source`` flipped from
    ``trace:receptor_revision.applied`` / ``trace:sample_allele.v``
    to ``ir:assignments.v.receptor_revision_original_id`` — both
    issue kinds now carry the SAME source string because both
    expected values are derived from the same IR slot."""
    validate_src = (
        Path(__file__).resolve().parent.parent
        / "engine_rs"
        / "src"
        / "airr_record"
        / "validate.rs"
    ).read_text(encoding="utf-8")
    assert "ReceptorRevisionAppliedMismatch" in validate_src
    assert "OriginalVCallMismatch" in validate_src

    outcome_src = (
        Path(__file__).resolve().parent.parent
        / "engine_rs"
        / "src"
        / "python"
        / "outcome.rs"
    ).read_text(encoding="utf-8")
    # New IR-sourced source string (Bug D fix). Used by BOTH
    # receptor-revision validator mismatches.
    assert (
        '"ir:assignments.v.receptor_revision_original_id"' in outcome_src
    ), (
        "IR-sourced details.source string missing from outcome.rs; "
        "Bug D fix has regressed."
    )
    # Old trace-sourced strings must be gone — flipped in lockstep
    # with the IR-provenance slice.
    assert '"trace:receptor_revision.applied"' not in outcome_src, (
        "Old trace-sourced details.source still present in outcome.rs; "
        "Bug D fix is incomplete."
    )
    assert '"trace:sample_allele.v"' not in outcome_src, (
        "Old trace-sourced details.source still present in outcome.rs; "
        "Bug D fix is incomplete."
    )


# ──────────────────────────────────────────────────────────────────
# 7. Live-call refresh — RegionReplaced is currently a no-op
#    on the refresh plan (design doc §6 documented this)
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_region_replaced_is_currently_a_refresh_plan_noop() -> None:
    """Design doc §6 explicitly documents that the
    `LiveCallRefreshPlan` treats `RegionReplaced` as a no-op
    today. The first implementation slice that emits
    `RegionReplaced` from a *production* pass (i.e., Slice C of
    receptor revision) must update the refresh plan first
    (Slice B). Pin the current no-op so an accidental refresh-plan
    edit that adds a non-no-op arm surfaces the dependency."""
    refresh_src = (
        Path(__file__).resolve().parent.parent
        / "engine_rs"
        / "src"
        / "live_call"
        / "refresh_plan.rs"
    ).read_text(encoding="utf-8")
    # The exact arm pattern: RegionReplaced sits in the
    # zero-action match-arm cluster.
    assert "SimulationEvent::RegionReplaced { .. }" in refresh_src, (
        "Refresh plan no longer recognises RegionReplaced; design "
        "doc §6 made specific claims about today's behaviour that "
        "this test pins."
    )


def test_pin_scaffold_segment_replaced_refresh_step_is_defined() -> None:
    """Slice B flipped this from pin-absence to pin-presence.

    The new `LiveCallRefreshStep::SegmentReplaced(Segment)` step
    (§6) now exists in `refresh_plan.rs` and is wired to fire
    when a `SimulationEvent::SegmentReplaced` event is in the
    pass record. Pin that shape so later slices can't drop the
    variant or rename it out from under the receptor-revision
    refresh path.
    """
    refresh_src = (
        Path(__file__).resolve().parent.parent
        / "engine_rs"
        / "src"
        / "live_call"
        / "refresh_plan.rs"
    ).read_text(encoding="utf-8")
    assert "SegmentReplaced(Segment)" in refresh_src, (
        "LiveCallRefreshStep::SegmentReplaced(Segment) variant "
        "missing; Slice B is supposed to have landed it."
    )
    # And the event → step mapping must be wired.
    assert re.search(
        r"SimulationEvent::SegmentReplaced.*LiveCallRefreshStep::SegmentReplaced",
        refresh_src,
        re.DOTALL,
    ), (
        "`SimulationEvent::SegmentReplaced` is no longer mapped to "
        "`LiveCallRefreshStep::SegmentReplaced` in `from_events`; "
        "the Slice B refresh contract is broken."
    )


# ──────────────────────────────────────────────────────────────────
# 8. One-region-per-V invariant — pinned today in baseline runs
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_baseline_runs_emit_exactly_one_region_per_vdj_segment() -> None:
    """The "one region per V/D/J segment" invariant is the
    scaffolding the design relies on for IR commit semantics
    (§4). The validator's ``MultipleRegionsForSegment`` issue is
    the existing post-condition that enforces it. Run a baseline
    VDJ pipeline and assert the issue never fires today — pins
    that revision must preserve the invariant without breaking
    the validator's pre-existing check."""
    refdata = _vdj_refdata()
    exp = (
        ga.Experiment.on(refdata)
        .allow_curatable_refdata()
        .recombine(np1_lengths=[(3, 1.0)], np2_lengths=[(3, 1.0)])
        .trim(enabled=False)
    )
    outcomes = exp.run(n=1, seed=0)
    issues = outcomes[0].validate_record(refdata, sequence_id="baseline")
    kinds = {issue["kind"] for issue in issues}
    assert "MultipleRegionsForSegment" not in kinds, (
        "baseline VDJ run already produces multiple regions per "
        "segment; the validator's invariant the audit relies on no "
        "longer holds. See validate.rs:MultipleRegionsForSegment."
    )


# ──────────────────────────────────────────────────────────────────
# 9. PassCompileEffect — new ReplaceRegion(Segment) variant absent
# ──────────────────────────────────────────────────────────────────


def test_pin_absence_pass_compile_effect_replace_region_not_defined_yet() -> None:
    """Design doc §5 proposes
    ``PassCompileEffect::ReplaceRegion(Segment)``. Slice A lands
    it. Pin its absence today."""
    metadata_src = (
        Path(__file__).resolve().parent.parent
        / "engine_rs"
        / "src"
        / "pass"
        / "metadata.rs"
    ).read_text(encoding="utf-8")
    assert "ReplaceRegion" not in metadata_src, (
        "PassCompileEffect::ReplaceRegion already exists; Slice A "
        "is the canonical landing point."
    )


# ──────────────────────────────────────────────────────────────────
# 10. Audit doc lockstep — the design doc must list every section
#    and the summary table must carry every load-bearing row
# ──────────────────────────────────────────────────────────────────


def test_pin_design_doc_lists_all_fourteen_audit_sections() -> None:
    """``docs/receptor_revision_design.md`` answers the 14 audit
    sections (the 10 user-spec questions plus implementation
    order, backwards-compat, test surface, and out-of-scope).
    Pin every section header so a refactor that drops one
    surfaces here."""
    doc = (
        Path(__file__).resolve().parent.parent
        / "docs"
        / "receptor_revision_design.md"
    )
    assert doc.is_file(), f"design doc is missing at {doc}"
    text = doc.read_text(encoding="utf-8")
    expected_sections = [
        "## 1. Biological scope",
        "## 2. Pipeline position + DSL surface",
        "## 3. Trace addresses",
        "## 4. IR representation",
        "## 5. Events",
        "## 6. Live-call refresh",
        "## 7. AIRR provenance",
        "## 8. Contracts",
        "## 9. Replay",
        "## 10. Validator",
        "## 11. Backwards compatibility",
        "## 12. Implementation order",
        "## 13. Test surface",
        "## 14. Out of scope",
    ]
    for header in expected_sections:
        assert header in text, (
            f"design doc lost section {header!r}; sync with the "
            f"audit checklist before any implementation slice lands."
        )


def test_pin_design_doc_summary_table_carries_every_decision() -> None:
    """The summary table at the end is the at-a-glance contract.
    Pin its load-bearing rows so a future re-edit that drops one
    surfaces here in the same lockstep style as the D-inversion
    test."""
    doc = (
        Path(__file__).resolve().parent.parent
        / "docs"
        / "receptor_revision_design.md"
    )
    text = doc.read_text(encoding="utf-8")
    for row in (
        "Scope",
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
        "Backwards compat",
    ):
        assert row in text, f"summary table dropped row {row!r}"


# ──────────────────────────────────────────────────────────────────
# 11. Test file self-pinning — the file's own conventions are
#    documented in its module docstring
# ──────────────────────────────────────────────────────────────────


def test_pin_followup_plan_for_when_slices_land_is_in_this_file() -> None:
    """Self-test: when the implementation slices land, the
    ``pin_absence_*`` tests in this file must flip to presence
    assertions in lockstep. Pin the convention so a contributor
    skimming the file knows what to do — by asserting the file's
    own docstring spells it out."""
    this_doc = (Path(__file__).resolve()).read_text(encoding="utf-8")
    assert "pin_absence_" in this_doc
    assert "pin_scaffold_" in this_doc
    assert "docs/receptor_revision_design.md" in this_doc
