"""Contract pins for the D-inversion RC-extension audit, post-
implementation.

Companion to
[`docs/d_inversion_extension_design.md`](../docs/d_inversion_extension_design.md).
The orientation-aware NP-region extension walks landed in the
implementation slice: under `ReverseComplement`, pool-left maps
to allele-right (consumes `trim_3`) and pool-right maps to
allele-left (consumes `trim_5`). The `extension_narrows_tie_set`
gate stays unchanged. Walker and validator oracle both route
through the same scoring kernel, so they remain in lockstep on
every inverted-D record.

The split:

- ``pin_narrowing_*`` tests assert the post-implementation
  behaviour — the RC-disable guards are gone, the trim caps
  swap, the canonical missed-evidence fixture now narrows. These
  flipped from the pre-implementation `pin_v1_boundary_*` tests
  in lockstep with the slice.
- ``pin_scaffold_*`` tests pin the surfaces the implementation
  reused — `matches_observed_with_orientation` exists,
  `extension_narrows_tie_set` remains orientation-agnostic,
  `ExtensionWalkState.orientation` is plumbed through, the
  `compatible_alleles_at_oriented` lookup is wired. These stay
  green to surface any refactor that drops a load-bearing
  primitive.

Closing condition: inverted D's adjacent NP evidence narrows the
tie set; the validator continues to pass without exception
filtering; cache parity stays green under `invert_d(prob=1.0)`.
"""
from __future__ import annotations

import re
from pathlib import Path

import GenAIRR as ga
from GenAIRR import _engine as ge


# ──────────────────────────────────────────────────────────────────
# Deterministic VDJ fixture. The contract uses a tiny VDJ refdata
# where two D alleles share their internal bytes but differ at the
# trimmed-off 5' end. That position would be reached by an NP-region
# extension under Forward; under inverted D, the same position is
# at the OTHER allele end, but the audit shows the symmetric case
# applies (an NP byte adjacent to the structural region in pool
# corresponds to an allele-extension byte the trim cap allows).
#
# The fixture is shaped so a deterministic seed lands the assigned
# D allele = d1*01 with trim_3 = 1 under invert_d(prob=1.0). At
# that trim, an extension step would attempt to match the trimmed-
# off byte at allele position `slice_end` (the audit's
# `*state.ref_end` under RC).
# ──────────────────────────────────────────────────────────────────


# ──────────────────────────────────────────────────────────────────
# 1. v1 boundary — RC extension walks are early-returned at entry
# ──────────────────────────────────────────────────────────────────


def test_pin_narrowing_walker_extensions_no_longer_skip_under_rc() -> None:
    """The implementation slice removed the early-return on RC in
    ``walk_left_extension`` and ``walk_right_extension``. Pin the
    inverse: neither walk carries an RC-disable guard anymore.

    A future contributor who re-introduces the v1 boundary would
    surface here before the missed-evidence fixture even runs."""
    src = (
        Path(__file__).resolve().parent.parent
        / "engine_rs"
        / "src"
        / "live_call"
        / "walker"
        / "extensions.rs"
    ).read_text(encoding="utf-8")
    rc_disable_pattern = re.compile(
        r"matches!\(\s*state\.orientation\s*,\s*"
        r"crate::assignment::SegmentOrientation::ReverseComplement\s*\)\s*\{\s*\n\s*return;",
        re.MULTILINE,
    )
    assert not rc_disable_pattern.search(src), (
        "walker/extensions.rs re-introduced an RC-disable early-"
        "return; the audit's implementation slice removed it."
    )
    # And both walks must consult orientation when computing
    # `candidate_ref_pos` (the audit §2 swap).
    assert "let is_rc = matches!" in src, (
        "walker/extensions.rs no longer branches on orientation when "
        "computing candidate_ref_pos; audit §2 contract drifted."
    )


def test_pin_narrowing_kernel_extensions_no_longer_skip_under_rc() -> None:
    """Mirror of the walker pin. ``score_alleles_with_extensions``
    in ``scoring.rs`` removed its ``skip_extensions`` guard."""
    src = (
        Path(__file__).resolve().parent.parent
        / "engine_rs"
        / "src"
        / "live_call"
        / "scoring.rs"
    ).read_text(encoding="utf-8")
    assert "skip_extensions" not in src, (
        "scoring.rs re-introduced the `skip_extensions` guard; the "
        "audit's implementation slice removed it."
    )
    # And the trim-cap swap helper landed.
    assert "let (left_cap, right_cap) = if is_rc" in src, (
        "scoring.rs no longer swaps the trim caps under RC; audit "
        "§3 contract drifted."
    )


def test_pin_narrowing_rust_unit_test_is_now_the_narrowing_pin() -> None:
    """The Rust-side test name flipped from
    ``rc_extensions_are_disabled_in_v1`` to
    ``rc_extensions_narrow_truth_allele``."""
    src = (
        Path(__file__).resolve().parent.parent
        / "engine_rs"
        / "src"
        / "live_call"
        / "tests"
        / "scoring.rs"
    ).read_text(encoding="utf-8")
    assert "rc_extensions_narrow_truth_allele" in src, (
        "Rust narrowing pin missing; the audit's §13 #1 fixture "
        "(NP byte narrows inverted-D tie set) is no longer pinned."
    )
    assert "rc_extensions_are_disabled_in_v1" not in src, (
        "Old `rc_extensions_are_disabled_in_v1` test still present; "
        "the audit's flip-in-lockstep convention requires removing it."
    )
    # And the trim-cap swap test is pinned too.
    assert "rc_extension_trim_cap_swap_pool_left_consumes_trim_3" in src, (
        "trim-cap swap regression test missing; audit §3 / §13 #4 "
        "requires it."
    )


# ──────────────────────────────────────────────────────────────────
# 2. Inverted-D records still validate structurally
# ──────────────────────────────────────────────────────────────────


def test_pin_inverted_d_records_validate_cleanly_under_v1_boundary() -> None:
    """Even with RC extensions disabled, inverted-D records still
    pass ``validate_records(refdata)`` without exception
    stripping. The release-tier IGH stack established this; we
    pin it on a smaller, deterministic fixture so a refactor that
    breaks the walker ↔ oracle lockstep surfaces here before
    reaching the slow release tier."""
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .invert_d(prob=1.0)
        .mutate(rate=0.01)
    )
    result = exp.run_records(n=5, seed=0)
    report = result.validate_records(exp.refdata)
    bad_kinds = [
        issue["kind"]
        for failure in report.failures
        for issue in failure["issues"]
    ]
    assert "AlleleCallTieSetMismatch" not in bad_kinds, (
        "AlleleCallTieSetMismatch surfaced under v1 RC-disable; "
        "the walker ↔ oracle lockstep was the load-bearing "
        "guarantee that lets the boundary stay deferrable."
    )


def test_pin_d_coordinates_unchanged_under_inverted_d() -> None:
    """Inverted D records continue to carry `d_germline_start`
    and `d_germline_end` in original allele orientation — the
    audit's load-bearing premise. The proposed extension
    implementation must NOT change this contract."""
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .invert_d(prob=1.0)
    )
    rec = exp.run_records(n=1, seed=0).records[0]
    assert rec["d_inverted"] is True
    # Half-open interval invariant: start <= end on every
    # populated D coordinate.
    if rec["d_germline_start"] is not None and rec["d_germline_end"] is not None:
        assert rec["d_germline_start"] <= rec["d_germline_end"]
    # And the sequence coords stay valid pool positions.
    if rec["d_sequence_start"] is not None and rec["d_sequence_end"] is not None:
        assert 0 <= rec["d_sequence_start"] <= rec["d_sequence_end"] <= len(rec["sequence"])


# ──────────────────────────────────────────────────────────────────
# 3. Synthetic fixture demonstrating missing evidence under RC
# ──────────────────────────────────────────────────────────────────


def test_pin_narrowing_synthetic_fixture_demonstrates_narrowing_under_rc() -> None:
    """Audit §6 / §13 #1, post-implementation. The Rust-side
    ``rc_extensions_narrow_truth_allele`` test in
    ``engine_rs/src/live_call/tests/scoring.rs`` constructs the
    canonical fixture (two D alleles tied on the retained
    structural bytes but distinguished by an NP byte adjacent to
    the trimmed-off allele position) and asserts the post-
    extension tie set narrows to the truth allele.

    This Python pin re-anchors the file-search contract at the
    test-name level so a refactor that drops the kernel-level
    fixture would surface here too."""
    scoring_tests = (
        Path(__file__).resolve().parent.parent
        / "engine_rs"
        / "src"
        / "live_call"
        / "tests"
        / "scoring.rs"
    ).read_text(encoding="utf-8")

    assert "rc_extensions_narrow_truth_allele" in scoring_tests, (
        "Kernel narrowing fixture missing; the audit's §13 #1 "
        "regression pin requires it."
    )
    # The narrowing assertion must remain explicit — a refactor
    # that changes the fixture's outcome assertions would surface
    # via this string check.
    assert "narrow to the truth allele" in scoring_tests, (
        "Rust narrowing fixture's outcome string drifted; audit §6 "
        "expected the truth-allele narrowing to be explicit."
    )


# ──────────────────────────────────────────────────────────────────
# 4. Existing scaffolding — primitives the audit's implementation
#    reuses (each must remain present for the slice to work)
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_matches_observed_with_orientation_primitive_exists() -> None:
    """The per-byte comparison primitive landed in the live-call
    cleanup slice. The audit's RC-extension implementation routes
    through it for the orientation-aware match decision. Pin the
    primitive's presence + its `Forward` / `ReverseComplement`
    arms so a refactor that drops one branch surfaces here."""
    src = (
        Path(__file__).resolve().parent.parent
        / "engine_rs"
        / "src"
        / "live_call"
        / "scoring.rs"
    ).read_text(encoding="utf-8")
    assert "fn matches_observed_with_orientation(" in src
    # Both orientation arms must be reachable from the primitive.
    assert "SegmentOrientation::Forward" in src
    assert "SegmentOrientation::ReverseComplement" in src


def test_pin_scaffold_extension_narrows_tie_set_remains_orientation_agnostic() -> None:
    """`extension_narrows_tie_set` is the gate the audit's
    implementation reuses unchanged. It must not gain a hidden
    orientation dependency between now and the implementation
    slice — that would silently change the v1 behaviour too."""
    src = (
        Path(__file__).resolve().parent.parent
        / "engine_rs"
        / "src"
        / "live_call"
        / "scoring.rs"
    ).read_text(encoding="utf-8")
    # Extract the function body and assert it doesn't reference
    # SegmentOrientation. The function is small; a substring grep
    # over its definition window is enough.
    match = re.search(
        r"pub\s+fn\s+extension_narrows_tie_set\([^)]*\)\s*->\s*bool\s*\{(.*?)^\}",
        src,
        re.DOTALL | re.MULTILINE,
    )
    assert match, "extension_narrows_tie_set not found in scoring.rs"
    body = match.group(1)
    assert "SegmentOrientation" not in body, (
        "extension_narrows_tie_set grew an orientation reference; "
        "the audit §4 expected it to stay orientation-agnostic so "
        "the implementation could reuse it unchanged."
    )


def test_pin_scaffold_extension_walk_state_carries_orientation() -> None:
    """`ExtensionWalkState` carries an `orientation` field that
    the implementation slice reads to decide which `ref_pos`
    direction the walk consumes. Pin the field's presence."""
    src = (
        Path(__file__).resolve().parent.parent
        / "engine_rs"
        / "src"
        / "live_call"
        / "walker"
        / "extensions.rs"
    ).read_text(encoding="utf-8")
    assert "pub orientation: crate::assignment::SegmentOrientation" in src, (
        "ExtensionWalkState lost its `orientation` field; the "
        "audit §1 scaffolding the implementation reads has drifted."
    )


def test_pin_scaffold_compatible_alleles_at_oriented_exists() -> None:
    """The walker's inverted-index lookup for orientation-aware
    scoring. The audit's implementation calls into this
    unchanged."""
    src = (
        Path(__file__).resolve().parent.parent
        / "engine_rs"
        / "src"
        / "live_call"
        / "reference_index.rs"
    ).read_text(encoding="utf-8")
    assert "fn compatible_alleles_at_oriented(" in src


# ──────────────────────────────────────────────────────────────────
# 5. Walker ↔ oracle lockstep — the load-bearing guarantee
# ──────────────────────────────────────────────────────────────────


def test_pin_cache_parity_under_inverted_d_with_rc_extensions_enabled() -> None:
    """Audit §13 #7: the incremental walker's cache and a fresh
    `from_existing_region` rebuild must agree under RC after
    extensions are enabled. Cache parity is the engine's other
    integrity gate (alongside the AIRR validator) — pinning it
    here surfaces regressions that bypass the validator but
    desynchronise the cached `SegmentLiveCall` from a from-scratch
    recompute."""
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .invert_d(prob=1.0)
        .mutate(rate=0.01)
    )
    refdata = exp.refdata
    result = exp.run_records(n=10, seed=0)
    assert result.outcomes is not None
    for i, outcome in enumerate(result.outcomes):
        parities = outcome.check_live_call_cache_parity(refdata)
        for p in parities:
            assert p["tie_set_matches"], (
                f"inverted-D record {i} segment {p['segment']}: "
                f"cache parity failed after RC extensions enabled — "
                f"cached={p['cached_tie_set']} fresh={p['fresh_tie_set']}"
            )


def test_pin_walker_and_oracle_agree_on_inverted_d_tie_set() -> None:
    """Audit §13 #3 + the closing release-tier guarantee: the
    walker's incremental tie set agrees with the validator's
    re-derivation on every inverted-D record. The release tests
    already pin this on a 50-record IGH stack; we pin it again on
    a smaller deterministic fixture so a regression surfaces in
    the fast tier."""
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .invert_d(prob=1.0)
        .mutate(rate=0.01)
    )
    result = exp.run_records(n=5, seed=0)
    report = result.validate_records(exp.refdata)
    for failure in report.failures:
        for issue in failure["issues"]:
            assert issue["kind"] != "AlleleCallTieSetMismatch", (
                f"walker ↔ oracle desync on inverted D: "
                f"failure={failure}; issue={issue}"
            )


# ──────────────────────────────────────────────────────────────────
# 6. Design-doc shape pins (mirrors the receptor-revision /
#    paired-end audit conventions)
# ──────────────────────────────────────────────────────────────────


def test_pin_design_doc_lists_fourteen_sections() -> None:
    """The audit doc follows the same 14-section structure as the
    receptor-revision and paired-end designs. Pin the section
    count so a future contributor can't quietly drop one."""
    doc = (
        Path(__file__).resolve().parent.parent
        / "audit-docs"
        / "d_inversion_extension_design.md"
    ).read_text(encoding="utf-8")
    headers = re.findall(r"^## (\d{1,2})\. ", doc, re.MULTILINE)
    assert headers == [str(i) for i in range(1, 15)], (
        f"d_inversion_extension_design.md section ordering changed; "
        f"got {headers}, expected 1..14."
    )


def test_pin_design_doc_summary_table_carries_every_decision() -> None:
    """The §Summary table is the canonical short-form of the audit
    decisions. Pin the row labels so a future contributor can't
    quietly drop a decision."""
    doc = (
        Path(__file__).resolve().parent.parent
        / "audit-docs"
        / "d_inversion_extension_design.md"
    ).read_text(encoding="utf-8")
    summary_match = re.search(
        r"## Summary table\n\n\| Question \| Decision \|.*?(?=\n##|\Z)",
        doc,
        re.DOTALL,
    )
    assert summary_match, "design doc summary table is missing"
    block = summary_match.group(0)
    rows = re.findall(r"^\| ([^|]+) \| ", block, re.MULTILINE)
    rows = [r.strip() for r in rows if r.strip() not in ("Question", "--- ")]
    rows = [r for r in rows if not set(r) <= {"-", " "}]
    for required in (
        "Pool direction ↔ allele direction under RC",
        "`trim_5` / `trim_3` cap mapping under RC",
        "Narrowing gate under RC",
        "Tie set vs. coordinates",
        "Missed-evidence fixtures",
        "Implementation surface",
        "Walker ↔ oracle parity",
        "Backwards compatibility",
        "Trace schema",
        "Validator issues",
    ):
        assert required in rows, (
            f"design doc summary table missing {required!r}; the "
            f"audit is the contract."
        )


def test_pin_lockstep_flip_convention_is_documented_in_this_file() -> None:
    """Audit doc closure: the pin-flip convention that took the
    pre-implementation `pin_v1_boundary_*` tests to the post-
    implementation `pin_narrowing_*` tests stays explained in
    this file so the next contributor who lands a similar audit
    can find the playbook."""
    this_doc = Path(__file__).read_text(encoding="utf-8")
    assert "pin_narrowing_" in this_doc
    assert "pin_scaffold_" in this_doc
    assert "lockstep" in this_doc, (
        "Module docstring no longer explains the pin-flip lockstep "
        "convention; the audit's closure playbook is no longer "
        "discoverable."
    )
