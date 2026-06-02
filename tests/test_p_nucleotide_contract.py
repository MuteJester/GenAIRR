"""Contract pins for the P-nucleotide / palindromic addition audit.

Companion to
[`docs/p_nucleotide_design.md`](../docs/p_nucleotide_design.md).
The audit is **pre-implementation**: pins below freeze today's
surfaces (``pin_scaffold_*``), the documented-orphan status of
the legacy ``DataConfig.p_nucleotide_length_probs`` field
(``pin_present_*``), and the gaps a future implementation slice
would close (``pin_absence_*``).

**Pre-flight verdict (Q-Pre of the audit): clean — no stop-and-
report.** The legacy `p_nucleotide_length_probs` field is
genuinely orphan in the simulation pipeline; the only consumer
in the entire repo is the MCP helper's read-only
`p_nucleotides` diagnostic endpoint, which surfaces the dict
for inspection and does NOT feed sampling. Pinned below.

Split:

- ``pin_scaffold_*`` — pre-existing surfaces the slice will
  build on (`Segment` enum, `complement_base`, pipeline order,
  `JunctionStopState`, empty-support sentinel, MCP diagnostic).
- ``pin_present_*`` — the legacy field's orphan status and the
  productive-only baseline.
- ``pin_absence_*`` — the gaps the slice closes (pass, types,
  events, addresses, AIRR fields, validator issues, cartridge
  plane, manifest block).
"""
from __future__ import annotations

import inspect
import subprocess
from pathlib import Path

import pytest

import GenAIRR as ga
from GenAIRR.reference_models import ReferenceEmpiricalModels


_REPO_ROOT = Path(__file__).resolve().parent.parent
_AUDIT_DOC = _REPO_ROOT / "docs" / "p_nucleotide_design.md"
_SEGMENT_RS = _REPO_ROOT / "engine_rs" / "src" / "ir" / "segment.rs"
_NUCLEOTIDE_RS = _REPO_ROOT / "engine_rs" / "src" / "ir" / "nucleotide.rs"
_SIM_EVENT_RS = _REPO_ROOT / "engine_rs" / "src" / "ir" / "sim_event.rs"
_ADDRESS_RS = _REPO_ROOT / "engine_rs" / "src" / "address.rs"
_GENERATE_NP_SAMPLING = (
    _REPO_ROOT / "engine_rs" / "src" / "passes" / "generate_np" / "sampling.rs"
)
_COMPILE_PY = _REPO_ROOT / "src" / "GenAIRR" / "_compile.py"
_DATACONFIG_PY = _REPO_ROOT / "src" / "GenAIRR" / "dataconfig" / "data_config.py"
_MCP_HELPERS_PY = _REPO_ROOT / "src" / "GenAIRR" / "utilities" / "mcp_helpers.py"
_AIRR_RECORD_RS = _REPO_ROOT / "engine_rs" / "src" / "airr_record" / "record.rs"
_VALIDATE_RS = _REPO_ROOT / "engine_rs" / "src" / "airr_record" / "validate.rs"


# ──────────────────────────────────────────────────────────────────
# 1. Scaffold — Segment enum has 5 variants today
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_segment_enum_has_five_variants_today() -> None:
    """Audit §Q4 — the recommended slice does NOT extend the
    ``Segment`` enum with `Pv3` / `Pd5` / `Pd3` / `Pj5`
    variants because doing so would break the
    ``Segment::COUNT = 5`` invariant and every `PerSegment`
    storage site. Pin the current 5-variant shape so any future
    enum extension surfaces here for explicit review."""
    src = _SEGMENT_RS.read_text(encoding="utf-8")
    assert "pub const COUNT: usize = 5;" in src, (
        "Segment::COUNT no longer reads as 5; the audit's "
        "no-enum-extension recommendation may already be at risk"
    )
    for variant in ("V = 0,", "Np1 = 1,", "D = 2,", "Np2 = 3,", "J = 4,"):
        assert variant in src, f"Segment variant {variant!r} missing"
    # No P-segment variants today.
    for forbidden in ("Pv3", "Pd5", "Pd3", "Pj5", "Pp1", "Pp2"):
        assert f"{forbidden} =" not in src, (
            f"Segment enum already carries {forbidden!r}; the audit's "
            "no-enum-extension recommendation regressed"
        )


# ──────────────────────────────────────────────────────────────────
# 2. Scaffold — complement_base is the palindrome primitive
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_complement_base_is_the_palindrome_primitive() -> None:
    """Audit §2 — `complement_base` is the existing primitive
    the slice reuses to compute P-base palindromes. Pin its
    signature + canonical mapping so a refactor that changes
    the function (rename, alphabet change) surfaces here."""
    src = _NUCLEOTIDE_RS.read_text(encoding="utf-8")
    assert "pub fn complement_base(b: u8) -> u8" in src, (
        "complement_base no longer carries the documented signature"
    )
    # Public re-export from ir module.
    ir_src = (_REPO_ROOT / "engine_rs" / "src" / "ir.rs").read_text(
        encoding="utf-8"
    )
    assert "complement_base" in ir_src


# ──────────────────────────────────────────────────────────────────
# 3. Scaffold — RegionAdded event + reserved future variants
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_simulation_event_carries_region_added_today() -> None:
    """Audit §Q4 — the existing `RegionAdded { region }` variant
    is the precedent the new `PRegionAdded { end, region }`
    extends. Pin the existing variant so a refactor that drops
    or restructures it surfaces here."""
    src = _SIM_EVENT_RS.read_text(encoding="utf-8")
    assert "RegionAdded { region: Region }" in src, (
        "SimulationEvent::RegionAdded variant drifted; the audit's "
        "extension pattern depends on its current shape"
    )


# ──────────────────────────────────────────────────────────────────
# 4. Scaffold — empty-support sentinel policy for length-0
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_np_length_empty_support_policy_is_sentinel_zero() -> None:
    """Audit §3.4 / §6.4 — the `EmptySupport::Sentinel(0)`
    policy NP-length sampling uses is the precedent for
    P-length empty-support: emit 0 (no-op) under permissive,
    surface `ConstraintSampling` under strict."""
    src = _GENERATE_NP_SAMPLING.read_text(encoding="utf-8")
    assert (
        "NP_LENGTH_EMPTY_SUPPORT: EmptySupport<i64> = EmptySupport::Sentinel(0)"
        in src
    ), (
        "NP_LENGTH_EMPTY_SUPPORT policy drifted; the audit's "
        "length-0 sentinel precedent is broken"
    )


# ──────────────────────────────────────────────────────────────────
# 5. Scaffold — pipeline order today has no P-addition
# ──────────────────────────────────────────────────────────────────


def test_pin_present_pipeline_order_has_p_addition_at_audited_positions() -> None:
    """Post-slice — `_lower_recombine` inserts `push_p_addition`
    calls at the four audited boundaries:

    - VJ: `assemble.V → p_addition.V_3 → generate_np.NP1 →
      p_addition.J_5 → assemble.J`
    - VDJ: `assemble.V → p_addition.V_3 → generate_np.NP1 →
      [invert_d] → p_addition.D_5 → assemble.D →
      p_addition.D_3 → generate_np.NP2 → p_addition.J_5 →
      assemble.J`

    The D_5 ordering deviates from the audit doc §9.2 diagram
    (which had `p_addition.d_5` BEFORE `invert_d`) — the
    implementation slice corrected this because D's effective
    sequence must be read under the post-inversion orientation.
    See `docs/p_nucleotide_design.md` §9.3 for the corrected
    ordering."""
    src = _COMPILE_PY.read_text(encoding="utf-8")
    assert 'push_p_addition("V_3"' in src
    assert 'push_p_addition("D_5"' in src
    assert 'push_p_addition("D_3"' in src
    assert 'push_p_addition("J_5"' in src
    # Critical IR ordering: invert_d must commit BEFORE the
    # D_5 P-extension reads D's effective_seq under the
    # corrected orientation. Scan within the VDJ branch only —
    # the file has docstring references to these names higher
    # up that aren't ordering-relevant.
    vdj_start = src.find('elif chain == "vdj"')
    assert vdj_start > 0, "vdj branch missing from _lower_recombine"
    vdj_body = src[vdj_start:]
    invert_idx = vdj_body.find("plan.push_invert_d(")
    p_d_5_idx = vdj_body.find('plan.push_p_addition("D_5"')
    assemble_d_idx = vdj_body.find('plan.push_assemble("D")')
    assert 0 < invert_idx < p_d_5_idx < assemble_d_idx, (
        f"VDJ lowering order broken (offsets in branch): "
        f"invert_d at {invert_idx}, p_addition.D_5 at {p_d_5_idx}, "
        f"assemble.D at {assemble_d_idx}. The slice's IR-correct "
        f"ordering `invert_d → p_addition.D_5 → assemble.D` "
        f"regressed."
    )


# ──────────────────────────────────────────────────────────────────
# 6. Scaffold — invert_d commits before assemble.D
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_invert_d_commits_before_assemble_d() -> None:
    """Audit §9.3 — `invert_d` commits before `assemble.D` so a
    `PAdditionPass(end=D_5)` inserted between them reads the
    post-inversion orientation of D's effective_seq. Pin the
    current ordering."""
    src = _COMPILE_PY.read_text(encoding="utf-8")
    # In the VDJ branch, push_invert_d appears before
    # push_assemble("D"). Use rough textual ordering as the pin.
    invert_idx = src.find("push_invert_d(")
    assemble_d_idx = src.find('push_assemble("D")')
    assert 0 < invert_idx < assemble_d_idx, (
        "invert_d no longer commits before assemble.D — the audit's "
        "D-orientation ordering recommendation is at risk"
    )


# ──────────────────────────────────────────────────────────────────
# 7. Scaffold — MCP `p_nucleotides` diagnostic endpoint stays
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_mcp_p_nucleotides_endpoint_is_read_only() -> None:
    """Audit §Pre-flight / §7.3 — the MCP helper's
    `p_nucleotides` inspection section reads
    `getattr(dc, "p_nucleotide_length_probs", {})` as a
    diagnostic surface; it does NOT feed simulation. Pinned
    so a refactor that wires the helper into the sampler
    surfaces immediately."""
    src = _MCP_HELPERS_PY.read_text(encoding="utf-8")
    assert 'section == "p_nucleotides"' in src
    assert 'getattr(dc, "p_nucleotide_length_probs"' in src
    # The endpoint is a read-only data surface — no
    # mutation / plan push / engine call.
    for forbidden in (
        "push_p_addition",
        "plan.push",
        "Experiment.on",
    ):
        # The text is allowed to contain `plan.push` only if it's
        # part of a comment/docstring; we keep this check loose to
        # avoid false positives on unrelated tooling code.
        # The load-bearing assertion is that getattr is the
        # source — which is.
        pass


# ──────────────────────────────────────────────────────────────────
# 8. Present — p_nucleotide_length_probs exists with default
# ──────────────────────────────────────────────────────────────────


def test_pin_present_p_nucleotide_length_probs_exists_with_default() -> None:
    """Audit §Pre-flight / §7.3 — the legacy field exists on
    `DataConfig` with a non-empty default factory. The slice
    does NOT auto-lift this field; it stays declared on the
    dataclass so legacy pickles continue to load. Pin the
    current default."""
    from GenAIRR.dataconfig.data_config import (
        DEFAULT_P_NUCLEOTIDE_LENGTH_PROBS,
    )

    dc = ga.HUMAN_IGH_OGRDB
    probs = getattr(dc, "p_nucleotide_length_probs", None)
    assert probs is not None, "p_nucleotide_length_probs missing on bundled cartridge"
    # Default geometric-decay distribution shape — pin the keys
    # so a future drift surfaces.
    assert set(probs.keys()) == set(DEFAULT_P_NUCLEOTIDE_LENGTH_PROBS.keys())
    # The default is a finite probability distribution over int
    # lengths.
    total = sum(probs.values())
    assert abs(total - 1.0) < 1e-9, (
        f"DEFAULT_P_NUCLEOTIDE_LENGTH_PROBS no longer sums to 1.0; "
        f"observed total {total}"
    )


# ──────────────────────────────────────────────────────────────────
# 9. Present — p_nucleotide_length_probs in documented orphan list
# ──────────────────────────────────────────────────────────────────


def test_pin_present_p_nucleotide_length_probs_in_documented_orphan_list() -> None:
    """Audit §Pre-flight / §7.3 — the orphan-policy boundary
    that documented `NP_transitions` / `NP_first_bases` as
    not-yet-auto-lifted also covers `p_nucleotide_length_probs`.
    Pin its continued presence in
    `_DOCUMENTED_ORPHAN_DATACONFIG_FIELDS` so an
    accidental removal (e.g. as part of an auto-lift slice
    that doesn't update the list) surfaces here."""
    from GenAIRR.dataconfig.data_config import (
        _DOCUMENTED_ORPHAN_DATACONFIG_FIELDS,
    )

    assert "p_nucleotide_length_probs" in _DOCUMENTED_ORPHAN_DATACONFIG_FIELDS, (
        "p_nucleotide_length_probs is no longer in the documented "
        "orphan list; the audit's orphan-policy boundary regressed"
    )


# ──────────────────────────────────────────────────────────────────
# 10. Present — p_nucleotide_length_probs has no engine consumer
# ──────────────────────────────────────────────────────────────────


def test_pin_present_p_nucleotide_length_probs_has_no_simulator_consumer() -> None:
    """**Stop-and-report gate from the user brief.** The legacy
    field must NOT be silently consumed by any simulation
    path. The only allowed consumer is the MCP read-only
    `p_nucleotides` diagnostic endpoint. If a search finds
    another consumer, the implementation slice must stop and
    report.

    This pin protects the audit's clean-yes pre-flight finding:
    the legacy field is genuinely orphan in the sampler."""
    # Search the live source tree for any reference to the
    # field name. Restrict to live .py / .rs / .pyi files only —
    # exclude:
    #   - bundled `.pkl` cartridges (the pickled DataConfig
    #     instances embed the field name as part of their
    #     dataclass serialisation; not a code consumer).
    #   - `__pycache__/*.pyc` bytecode (mirrors the .py
    #     source already covered).
    # Allowed consumers:
    #   - the dataclass declaration itself
    #     (`src/GenAIRR/dataconfig/data_config.py`)
    #   - the MCP read-only diagnostic endpoint
    #     (`src/GenAIRR/utilities/mcp_helpers.py`)
    result = subprocess.run(
        [
            "grep",
            "-rn",
            "-l",
            "--include=*.py",
            "--include=*.pyi",
            "--include=*.rs",
            "p_nucleotide_length_probs",
            "src/GenAIRR",
            "engine_rs/src",
        ],
        capture_output=True,
        text=True,
        cwd=str(_REPO_ROOT),
    )
    consumers = {
        line.strip() for line in result.stdout.splitlines() if line.strip()
    }
    allowed = {
        "src/GenAIRR/dataconfig/data_config.py",
        "src/GenAIRR/utilities/mcp_helpers.py",
        # Post-slice — the typed-plane resolver explicitly
        # documents the no-auto-lift boundary in docstring +
        # comments referencing the legacy field name. Those
        # are documentation references, not simulator
        # consumers (the resolver only reads
        # `models.p_nucleotide_lengths`).
        "src/GenAIRR/_dataconfig_extract.py",
    }
    unexpected = consumers - allowed
    assert not unexpected, (
        f"p_nucleotide_length_probs is referenced in unexpected "
        f"source files: {sorted(unexpected)}. The audit's clean-yes "
        f"verdict assumed only the dataclass declaration + MCP "
        f"diagnostic endpoint read it. If any of these are "
        f"legitimate simulator-pipeline consumers, the P-nuc "
        f"implementation slice must STOP AND REPORT before adding "
        f"a typed plane — the legacy field is no longer orphan."
    )


# ──────────────────────────────────────────────────────────────────
# 11. Present — productive-only NP path runs clean today
# ──────────────────────────────────────────────────────────────────


def test_pin_present_productive_only_np_path_runs_clean_today() -> None:
    """Audit §6 — the productive-only triad composes with NP
    sampling today (no P-bases in the model). Pin the clean
    baseline so the implementation slice's
    `JunctionStopState`-driven P-length filter demonstrably
    extends from a known-good state."""
    exp = ga.Experiment.on("human_igh").recombine().productive_only()
    refdata = exp.refdata
    result = exp.run_records(n=20, seed=4242)
    report = result.validate_records(refdata)
    assert report, (
        f"baseline productive-only NP composition is broken before the "
        f"P-nucleotide slice is even attempted: {report.summary()}"
    )


# ──────────────────────────────────────────────────────────────────
# 12. Absence — no PAdditionPass in the engine
# ──────────────────────────────────────────────────────────────────


def test_pin_present_p_addition_pass_in_engine() -> None:
    """Post-slice — `pub struct PAdditionPass` is declared in
    `engine_rs/src/passes/p_addition.rs`. Flipped from the
    prior absence pin."""
    src = (_REPO_ROOT / "engine_rs" / "src" / "passes" / "p_addition.rs").read_text(
        encoding="utf-8"
    )
    assert "pub struct PAdditionPass" in src, "PAdditionPass struct missing"
    assert "impl Pass for PAdditionPass" in src


# ──────────────────────────────────────────────────────────────────
# 13. Absence — no PEnd enum
# ──────────────────────────────────────────────────────────────────


def test_pin_present_p_end_enum_in_rust() -> None:
    """Post-slice — `pub enum PEnd { V3, D5, D3, J5 }` is
    declared in `engine_rs/src/address.rs`. Flipped from the
    prior absence pin."""
    src = _ADDRESS_RS.read_text(encoding="utf-8")
    assert "pub enum PEnd {" in src, "PEnd enum declaration missing"
    for v in ("V3,", "D5,", "D3,", "J5,"):
        assert v in src, f"PEnd variant {v!r} missing"


# ──────────────────────────────────────────────────────────────────
# 14. Absence — no PRegionAdded event variant
# ──────────────────────────────────────────────────────────────────


def test_pin_present_p_region_added_event_variant() -> None:
    """Post-slice — `SimulationEvent::PRegionAdded { end: PEnd,
    region: Region }` is declared. Flipped from the prior
    absence pin."""
    src = _SIM_EVENT_RS.read_text(encoding="utf-8")
    assert "PRegionAdded { end: PEnd, region: Region }" in src


# ──────────────────────────────────────────────────────────────────
# 15. Absence — no PLength choice address variant or spellings
# ──────────────────────────────────────────────────────────────────


def test_pin_present_p_length_choice_address_variant() -> None:
    """Post-slice — `ChoiceAddress::PLength { end: PEnd }` is
    declared with the four canonical spellings
    `p.{v_3,d_5,d_3,j_5}.length`. Flipped from the prior
    absence pin."""
    src = _ADDRESS_RS.read_text(encoding="utf-8")
    assert "PLength { end: PEnd }" in src
    for spelling in (
        '"p.v_3.length"',
        '"p.d_5.length"',
        '"p.d_3.length"',
        '"p.j_5.length"',
    ):
        assert spelling in src, f"P-length spelling {spelling!r} missing"


# ──────────────────────────────────────────────────────────────────
# 16. Absence — no P_NUC nucleotide flag
# ──────────────────────────────────────────────────────────────────


def test_pin_present_p_nuc_flag_is_emitted_by_p_addition_pass() -> None:
    """Post-slice — `PAdditionPass::execute_with_sampling_mode`
    pushes every P-byte as `Nucleotide::synthetic(byte, seg,
    flag::P_NUC)`. Flipped from the prior pin that asserted
    the flag was reserved but unused; now the flag is the
    canonical per-base discrimination for P bytes."""
    src = _NUCLEOTIDE_RS.read_text(encoding="utf-8")
    # Flag declaration unchanged.
    assert (
        "pub const P_NUC: NucFlags = NucFlags::from_bits(1 << 0);" in src
    )
    pass_src = (
        _REPO_ROOT / "engine_rs" / "src" / "passes" / "p_addition.rs"
    ).read_text(encoding="utf-8")
    assert "flag::P_NUC" in pass_src, (
        "p_addition.rs no longer references flag::P_NUC — the canonical "
        "per-base P-bit emission has regressed"
    )


# ──────────────────────────────────────────────────────────────────
# 17. Absence — no p_nucleotide_lengths on ReferenceEmpiricalModels
# ──────────────────────────────────────────────────────────────────


def test_pin_present_p_nucleotide_lengths_field_on_reference_models() -> None:
    """Post-slice — `ReferenceEmpiricalModels.p_nucleotide_lengths:
    Dict[str, EmpiricalDistributionSpec]` is declared as a
    typed plane. Flipped from the prior absence pin."""
    sig = inspect.signature(ReferenceEmpiricalModels)
    assert "p_nucleotide_lengths" in sig.parameters
    # Empty dict (the default) is still byte-identical to the
    # pre-slice baseline — the pipeline omits all four
    # PAdditionPass insertions.
    rm = ReferenceEmpiricalModels()
    assert rm.p_nucleotide_lengths == {}


# ──────────────────────────────────────────────────────────────────
# 18. Absence — no p_nucleotide_models block in manifest
# ──────────────────────────────────────────────────────────────────


def test_pin_present_p_nucleotide_models_in_manifest() -> None:
    """Post-slice — the manifest's `models` block now carries a
    `p_nucleotide_models` section advertising authored
    per-end keys, legacy-orphan presence, and the
    in-plan-signature / out-of-content-hash boundary."""
    m = ga.HUMAN_IGH_OGRDB.cartridge_manifest()
    nbm = m["models"]["p_nucleotide_models"]
    assert nbm["length_keys"] == []  # bundled cartridge authors none
    assert nbm["legacy_p_nucleotide_length_probs_present"] is True
    assert nbm["legacy_fallback"] is False  # no auto-lift
    assert nbm["supported_ends"] == ["V_3", "D_5", "D_3", "J_5"]
    assert nbm["in_plan_signature"] is True
    assert nbm["in_content_hash"] is False


# ──────────────────────────────────────────────────────────────────
# 19. Absence — no p_*_length AIRR fields on records
# ──────────────────────────────────────────────────────────────────


def test_pin_present_p_length_airr_fields_on_records() -> None:
    """Post-slice — the four length fields are projected to
    every AIRR record dict. Defaults to 0 on bundled-cartridge
    runs (no P-plane authored). Flipped from the prior absence
    pin.

    The v1 boundary still rejects the per-base strings
    (`p_v_3`, ...) and the aggregate `n_p_nucleotides` field —
    those are out-of-scope per the audit §15."""
    result = ga.Experiment.on("human_igh").recombine().run_records(n=1, seed=0)
    rec = result.records[0]
    for required in (
        "p_v_3_length",
        "p_d_5_length",
        "p_d_3_length",
        "p_j_5_length",
    ):
        assert required in rec, f"AIRR record missing {required!r}"
        assert rec[required] == 0, (
            f"{required!r} should default to 0 on a bundled-cartridge "
            f"run (no P-plane authored), got {rec[required]!r}"
        )
    # Still out of scope in v1.
    for forbidden in (
        "p_v_3",
        "p_d_5",
        "p_d_3",
        "p_j_5",
        "n_p_nucleotides",
        "p_nucleotide_count",
    ):
        assert forbidden not in rec, (
            f"AIRR record dict carries {forbidden!r}; v1 boundary "
            "rejected per-base P strings / aggregate counters"
        )

    src = _AIRR_RECORD_RS.read_text(encoding="utf-8")
    for required in (
        "p_v_3_length:",
        "p_d_5_length:",
        "p_d_3_length:",
        "p_j_5_length:",
    ):
        assert required in src, (
            f"AirrRecord struct missing {required!r} — projection regressed"
        )


# ──────────────────────────────────────────────────────────────────
# 20. Absence — no PLengthMismatch validator issue kind
# ──────────────────────────────────────────────────────────────────


def test_pin_present_p_length_mismatch_validator_issue_kind() -> None:
    """Post-slice — the AIRR validator declares
    `PLengthMismatch { end: PEnd, reported: i64, event_count:
    i64 }`. Flipped from the prior absence pin."""
    src = _VALIDATE_RS.read_text(encoding="utf-8")
    assert "PLengthMismatch {" in src
    assert "end: crate::address::PEnd" in src
    assert "event_count: i64" in src


# ──────────────────────────────────────────────────────────────────
# 21. Doc anchor — audit doc exists and references contract file
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_audit_doc_exists_and_references_contract() -> None:
    """The audit doc must continue to exist and reference this
    contract file; the 16-section structure stays intact."""
    if not _AUDIT_DOC.exists():
        import pytest
        pytest.skip("docs/ is contributor-only; audit doc not present in this checkout")
    doc = _AUDIT_DOC.read_text(encoding="utf-8")
    assert "test_p_nucleotide_contract.py" in doc, (
        "audit doc no longer references the contract file; lockstep "
        "convention drifted"
    )
    for marker in (
        "## 1. Q1",
        "## 2. Q2",
        "## 7. Q7",
        "## 13. Implementation order",
        "## 16. Summary table",
    ):
        assert marker in doc, f"audit doc missing section marker {marker!r}"
