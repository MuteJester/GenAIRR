"""Pre-implementation contract for D-segment inversion.

Companion to [`docs/d_inversion_design.md`](../docs/d_inversion_design.md).
D inversion isn't implemented yet; this file pins **today's
absence** of every surface the design doc proposes. When the
implementation slice lands, the absence assertions flip to
presence assertions (one-line edits per test, tracked by `# TODO`
markers). That flip is itself a useful signal that the slice
touched the right surface — a future contributor who renames
`d_inverted` to `d_orientation` without updating this file will
break it.

The split:

- `pin_absence_*` tests assert today's behavior — the field
  doesn't exist, the DSL method isn't there, no trace address
  starts with `sample_allele.d.inverted`. These MUST PASS
  before the implementation slice.
- `pin_design_*` tests assert intent that's already encoded in
  the engine (the `NucFlags::INVERTED` slot, the
  `ChoiceValue::Bool` docstring, the `AnchorPreserved::D`
  contract surface). These are stability pins on scaffolding
  that's already in place.

When Slice A lands, the relevant `pin_absence_*` tests get
inverted in lockstep with the new behavior and the audit doc
crosses out the corresponding item in §10.
"""
from __future__ import annotations

import inspect
import re
from pathlib import Path

import GenAIRR as ga
from GenAIRR import _engine as ge


# ──────────────────────────────────────────────────────────────────
# Deterministic VDJ fixture (anchorless D — matches IGH reality).
# ──────────────────────────────────────────────────────────────────


def _vdj_refdata() -> "ge.RefDataConfig":
    cfg = ge.RefDataConfig.vdj()
    cfg.add_v_allele("v1*01", "v1", b"AAACCCGGG", anchor=6)
    cfg.add_d_allele("d1*01", "d1", b"AAATTTGGG")
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
# 1. DSL surface absence — `Experiment.invert_d` does not exist
# ──────────────────────────────────────────────────────────────────


def test_pin_presence_experiment_invert_d_method_lands_in_slice_d() -> None:
    """Slice D landed `Experiment.invert_d(prob=...)` per
    `docs/d_inversion_design.md` §1. Pin that the canonical name
    stays under that spelling — a rename to (say) ``Experiment.d_inversion``
    would silently break every example in the documentation that
    references the method."""
    assert hasattr(ga.Experiment, "invert_d"), (
        "Experiment.invert_d went missing; Slice D is the canonical "
        "home for the DSL surface and must keep this method exposed."
    )
    # The signature also pins the keyword-only `prob` parameter so a
    # future refactor that broadens / renames it surfaces here.
    import inspect

    sig = inspect.signature(ga.Experiment.invert_d)
    assert "prob" in sig.parameters
    # `prob` must be keyword-only — positional `prob` would be a
    # silent API expansion that breaks chained `.invert_d(0.5)`
    # callers reading the parameter as `self`.
    assert sig.parameters["prob"].kind == inspect.Parameter.KEYWORD_ONLY


def test_pin_absence_recombine_does_not_carry_d_inversion_kwarg() -> None:
    """The design doc deliberately picks `.invert_d(prob=…)` over a
    `recombine(d_inversion_prob=…)` kwarg. Pin that the kwarg does
    NOT appear today — if it does, the alternative was already
    chosen and the design doc is stale."""
    sig = inspect.signature(ga.Experiment.recombine)
    assert "d_inversion_prob" not in sig.parameters
    assert "invert_d_prob" not in sig.parameters


# ──────────────────────────────────────────────────────────────────
# 2. AIRR field absence — `d_inverted` not in records today
# ──────────────────────────────────────────────────────────────────


def test_pin_presence_d_inverted_field_lands_in_slice_e() -> None:
    """Slice E landed the `d_inverted: bool` field per
    `docs/d_inversion_design.md` §6.3. Pin that the AIRR record
    exposes the field under that exact name (not `d_orientation`,
    not `inverted_d`) — a rename would silently break every
    downstream consumer that depended on the canonical spelling
    after Slice E shipped."""
    result = _baseline_experiment().run_records(n=1, seed=0)
    rec = result.records[0]
    assert "d_inverted" in rec, (
        "AIRR record no longer exposes `d_inverted`; Slice E is "
        "the canonical home for the provenance flag and must keep "
        "this spelling. See docs/d_inversion_design.md §6.3."
    )
    # The baseline pipeline doesn't call invert_d(), so the field
    # must default to False — captures the "absent step ⇒ field
    # defaults to False" invariant.
    assert rec["d_inverted"] is False


def test_pin_absence_d_orientation_field_not_emitted_today() -> None:
    """Same pin, alternative spelling — the design doc considered
    `d_orientation` and rejected it. Pin that we didn't land it
    under that name either."""
    rec = _baseline_experiment().run_records(n=1, seed=0).records[0]
    assert "d_orientation" not in rec


# ──────────────────────────────────────────────────────────────────
# 3. Trace address vocabulary — `sample_allele.d.inverted` is not
#    in today's trace output
# ──────────────────────────────────────────────────────────────────


def test_pin_absence_no_d_inversion_trace_address_in_baseline_run() -> None:
    """The proposed trace address per §2 of the design doc is
    `sample_allele.d.inverted`. A baseline VDJ run with no
    `invert_d` step must not emit it (and there's no other source
    today)."""
    outcomes = _baseline_experiment().run(n=1, seed=0)
    addrs = {r.address for r in outcomes[0].trace().choices()}
    assert "sample_allele.d.inverted" not in addrs
    # And no neighboring addresses leaked in either — pin the
    # naming surface for cleanliness.
    assert not any(a.startswith("recombine.d.") for a in addrs)
    assert not any(a == "sample_allele.d.orientation" for a in addrs)


def test_pin_presence_d_inversion_address_constant_in_engine_source() -> None:
    """Slice C landed `ChoiceAddress::SampleAlleleDInverted` with
    on-disk spelling `sample_allele.d.inverted`. Pin that the
    constant stays present — a rename here would silently break
    every trace file that recorded the address."""
    addr_src = (
        Path(__file__).resolve().parent.parent
        / "engine_rs"
        / "src"
        / "address.rs"
    ).read_text(encoding="utf-8")
    assert "sample_allele.d.inverted" in addr_src, (
        "engine_rs/src/address.rs no longer references "
        "'sample_allele.d.inverted'; Slice C is the canonical "
        "home for the trace address and must keep this spelling."
    )


# ──────────────────────────────────────────────────────────────────
# 4. AlleleInstance shape — no `orientation` field today
# ──────────────────────────────────────────────────────────────────


def test_pin_absence_d_allele_instance_does_not_expose_orientation_today() -> None:
    """§4 of the design adds `AlleleInstance.orientation`. Until
    that lands, the Python-side allele view (`d_allele(0)`) must
    NOT expose `orientation` — and the bundled D pool must NOT
    expose any per-allele orientation metadata."""
    cfg = _vdj_refdata()
    d = cfg.d_allele(0)
    surface = {a for a in dir(d) if not a.startswith("_")}
    assert "orientation" not in surface, (
        "Allele view already carries `orientation`; coordinate "
        "with §4 of the design doc before the slice lands."
    )


# ──────────────────────────────────────────────────────────────────
# 5. SegmentOrientation enum absence
# ──────────────────────────────────────────────────────────────────


def test_pin_absence_segment_orientation_not_exported_from_engine() -> None:
    """§4 proposes a `SegmentOrientation { Forward, ReverseComplement }`
    enum. Before the slice lands no Python export should exist
    under that name (or its plausible misspellings)."""
    for name in (
        "SegmentOrientation",
        "Orientation",
        "DOrientation",
        "AlleleOrientation",
    ):
        assert not hasattr(ge, name), (
            f"GenAIRR._engine.{name} exists; flip this test in lockstep "
            f"with the design doc §4 implementation slice."
        )


# ──────────────────────────────────────────────────────────────────
# 6. Scaffolding stability — what the engine *already* anticipated
# ──────────────────────────────────────────────────────────────────


def test_pin_design_nuc_flags_inverted_slot_stays_documented() -> None:
    """`engine_rs/src/ir/nucleotide.rs` defines `NucFlags::INVERTED`
    (bit 3) with a docstring that names D inversion as the
    intended use. Pin that the slot stays under that name — a
    rename would break the design doc's §5 proposal."""
    nuc_src = (
        Path(__file__).resolve().parent.parent
        / "engine_rs"
        / "src"
        / "ir"
        / "nucleotide.rs"
    ).read_text(encoding="utf-8")
    assert "pub const INVERTED: NucFlags" in nuc_src, (
        "NucFlags::INVERTED was renamed or removed; update "
        "docs/d_inversion_design.md §5 (assembly emission) "
        "before the implementation slice lands."
    )
    # And the docstring still names D inversion as the use case —
    # otherwise a future reader misses the link.
    assert re.search(r"inverted.*D segment|D segment.*inverted", nuc_src, re.I), (
        "NucFlags::INVERTED docstring no longer cites D-segment "
        "inversion; restore the cross-reference to the design doc."
    )


def test_pin_design_choice_value_bool_still_cites_d_inversion() -> None:
    """`ChoiceValue::Bool` in `engine_rs/src/trace.rs` documents
    "D inversion: yes/no" as a canonical use case (§2). Pin the
    docstring so the design doc and the engine stay aligned."""
    trace_src = (
        Path(__file__).resolve().parent.parent
        / "engine_rs"
        / "src"
        / "trace.rs"
    ).read_text(encoding="utf-8")
    assert re.search(r"D inversion.*Bool|Bool.*D inversion", trace_src, re.I), (
        "ChoiceValue::Bool docstring no longer cites D inversion; "
        "this is the canonical reservation for the slice."
    )


def test_pin_design_anchor_preserved_d_compiles_but_is_not_in_productive() -> None:
    """`AnchorPreserved::D` exists at the engine level; the design
    doc §7 explicitly notes it stays OUT of the default productive
    bundle. Pin that today's productive bundle covers V + J only —
    if D ever joins, the design's "productive_only unaffected by
    inversion" claim needs revisiting."""
    names = list(ge.productive().names())
    assert "anchor_preserved.v" in names
    assert "anchor_preserved.j" in names
    assert "anchor_preserved.d" not in names


# ──────────────────────────────────────────────────────────────────
# 7. Audit doc lockstep — design doc must list every section
# ──────────────────────────────────────────────────────────────────


def test_pin_design_doc_lists_all_eight_audit_questions() -> None:
    """`docs/d_inversion_design.md` answers the 8 audit questions
    the user asked for. Pin every section header so a future
    refactor that drops one surfaces here."""
    doc = (
        Path(__file__).resolve().parent.parent
        / "docs"
        / "d_inversion_design.md"
    )
    assert doc.is_file(), f"design doc is missing at {doc}"
    text = doc.read_text(encoding="utf-8")
    expected_sections = [
        "## 1. Where in the pipeline?",
        "## 2. Trace addresses",
        "## 3. SimulationEvent",
        "## 4. IR representation",
        "## 5. Assembly behavior",
        "## 6. AIRR projection",
        "## 7. Contracts",
        "## 8. Replay",
    ]
    for header in expected_sections:
        assert header in text, (
            f"design doc lost section {header!r}; sync with the audit "
            f"checklist before the implementation slice lands."
        )


def test_pin_design_doc_summary_table_carries_every_decision() -> None:
    """The summary table at the end is the at-a-glance contract.
    Pin its load-bearing rows — these are the decisions the
    implementation must match."""
    doc = (
        Path(__file__).resolve().parent.parent
        / "docs"
        / "d_inversion_design.md"
    )
    text = doc.read_text(encoding="utf-8")
    for row in (
        "Pipeline position",
        "DSL surface",
        "Trace address",
        "SimulationEvent",
        "IR carrier",
        "Assembly",
        "AIRR",
        "Productive bundle",
        "Replay",
        "Backwards compat",
    ):
        assert row in text, f"summary table dropped row {row!r}"


# ──────────────────────────────────────────────────────────────────
# 8. Test file declares its own intent — what each block becomes
#    once the implementation slice lands
# ──────────────────────────────────────────────────────────────────


def test_pin_followup_plan_for_when_slice_lands_is_in_this_file() -> None:
    """Self-test: when the implementation slice lands, the
    `pin_absence_*` tests in this file must flip to presence
    assertions in lockstep. Pin the convention so a contributor
    skimming the file knows what to do — by asserting the file's
    own docstring spells it out."""
    this_doc = (Path(__file__).resolve()).read_text(encoding="utf-8")
    # The module docstring must reference both the flip semantics
    # and the design doc, so the convention is discoverable from
    # the file's first paragraph.
    assert "pin_absence_" in this_doc
    assert "docs/d_inversion_design.md" in this_doc
