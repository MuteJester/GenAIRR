"""Contract pins for the Mutation Provenance Counters audit.

Companion to
[`docs/mutation_provenance_audit.md`](../docs/mutation_provenance_audit.md).
The audit is pre-implementation: it freezes today's global
counter shape and the per-event provenance the future slice will
aggregate (``pin_scaffold_*``), plus the surfaces the slice will
add (``pin_absence_*``).

Split:

- ``pin_scaffold_*`` tests freeze today's contract: the four
  existing AIRR counter fields, ``BaseChanged`` carrying
  ``segment``, the pass-name filter pattern indels use, the
  trace-sourced vs event-sourced vs IR-sourced counter
  classification, and the segment-rate zero-rate behavioural
  pin that the future ``n_np_mutations`` field would
  satisfy.
- ``pin_absence_*`` tests freeze the gaps the implementation
  slice closes: no per-segment SHM AIRR fields, no
  ``MutationCountMismatch`` validator issue kinds, no
  ``per_segment_counters`` manifest field, no Python typed-event
  exposure on ``EventRecord``.

When the slice lands the relevant ``pin_absence_*`` flips to
``pin_present_*`` in lockstep.
"""
from __future__ import annotations

import re
from pathlib import Path

import pytest

import GenAIRR as ga


_REPO_ROOT = Path(__file__).resolve().parent.parent


# ──────────────────────────────────────────────────────────────────
# 1. Scaffold — BaseChanged carries segment + germline_pos
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_basechanged_carries_segment_field() -> None:
    """Audit §3 pre-flight: ``SimulationEvent::BaseChanged``
    carries ``segment: Segment`` AND ``germline_pos:
    Option<u16>``. The implementation slice keys off the
    ``segment`` field for per-bucket SHM counters. Pinned at the
    source level so a refactor that strips the segment carries
    surface before the slice lands."""
    sim_event_src = (
        _REPO_ROOT / "engine_rs" / "src" / "ir" / "sim_event.rs"
    ).read_text(encoding="utf-8")
    # Find the BaseChanged variant body and verify it carries
    # ``segment`` + ``germline_pos``.
    base_changed_match = re.search(
        r"BaseChanged\s*\{(.*?)\}", sim_event_src, re.DOTALL
    )
    assert base_changed_match is not None, (
        "BaseChanged variant not found in sim_event.rs; audit §3 "
        "pre-flight assumption broken."
    )
    body = base_changed_match.group(1)
    assert "segment: Segment" in body, (
        "BaseChanged no longer carries segment; per-segment SHM "
        "counter slice would need an event-shape change first."
    )
    assert "germline_pos:" in body, (
        "BaseChanged no longer carries germline_pos; orthogonal "
        "downstream consumers may break."
    )


# ──────────────────────────────────────────────────────────────────
# 2. Scaffold — pass-name constants exist with canonical values
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_pass_name_constants_pinned() -> None:
    """Audit §4: the per-segment SHM counter slice filters
    ``outcome.events()`` by ``pass_name in {MUTATE_UNIFORM,
    MUTATE_S5F}``. Pin the constants here so a refactor that
    renames them forces the slice author to update both sides in
    lockstep."""
    address_src = (
        _REPO_ROOT / "engine_rs" / "src" / "address.rs"
    ).read_text(encoding="utf-8")
    expected = {
        "MUTATE_UNIFORM": '"mutate.uniform"',
        "MUTATE_S5F": '"mutate.s5f"',
        "CORRUPT_INDEL": '"corrupt.indel"',
        "CORRUPT_PCR": '"corrupt.pcr"',
        "CORRUPT_QUALITY": '"corrupt.quality"',
    }
    for const, value in expected.items():
        pattern = rf"pub const {const}: &str = {re.escape(value)};"
        assert re.search(pattern, address_src), (
            f"pass-name constant {const}={value} missing; rename "
            "blocks the per-segment counter slice."
        )


# ──────────────────────────────────────────────────────────────────
# 3. Scaffold — global counters today
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_existing_airr_counter_fields_present() -> None:
    """Audit §1: enumerate the global counter fields present on
    every AIRR record today. The proposed per-segment SHM slice
    must preserve all of these (it ONLY adds new fields)."""
    result = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(count=5)
        .polymerase_indels(count=2)
        .pcr_amplify(count=3)
        .sequencing_errors(count=2)
        .run_records(n=1, seed=0)
    )
    rec = result.records[0]
    for required in (
        "n_mutations",
        "mutation_rate",
        "n_pcr_errors",
        "n_quality_errors",
        "n_indels",
        "n_v_indels",
        "n_d_indels",
        "n_j_indels",
        "end_loss_5_length",
        "end_loss_3_length",
    ):
        assert required in rec, (
            f"existing AIRR counter field {required!r} missing; "
            "global counter contract regressed."
        )


def test_pin_scaffold_n_mutations_is_realised_count_for_shm() -> None:
    """Audit §2 / §1: ``n_mutations`` reflects the realised SHM
    count (not the attempted draw). With ``count=15`` and no
    constraints the realised count equals 15."""
    result = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(count=15)
        .run_records(n=3, seed=0)
    )
    for r in result:
        assert r["n_mutations"] == 15


def test_pin_scaffold_n_indels_uses_pass_name_filter() -> None:
    """Audit §3: the indel-counter implementation is the
    reference pattern for the proposed SHM slice. Pin its
    presence at source level so a contributor who refactors the
    pattern propagates the change to the future SHM counters."""
    builder_src = (
        _REPO_ROOT
        / "engine_rs"
        / "src"
        / "airr_record"
        / "builder.rs"
    ).read_text(encoding="utf-8")
    assert "address::CORRUPT_INDEL" in builder_src, (
        "n_indels filter no longer uses address::CORRUPT_INDEL; "
        "the future SHM counter slice must follow the same "
        "pattern."
    )
    # And the per-segment switch is the precedent for the
    # NP-rollup decision: V/D/J get dedicated fields, Np1/Np2
    # roll into the total only.
    assert "Segment::Np1 | Segment::Np2" in builder_src, (
        "n_indels NP-rollup pattern no longer present in builder.rs."
    )


# ──────────────────────────────────────────────────────────────────
# 4. Scaffold — corruption passes do NOT bump n_mutations
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_pcr_and_quality_do_not_increment_n_mutations() -> None:
    """Audit §2: ``BaseChanged`` events emitted by corruption
    passes (PCR, quality, N-corruption) must NOT contribute to
    ``n_mutations`` and (by extension, after the slice) must NOT
    contribute to any per-segment SHM counter. Pin the
    behavioural separation now so the slice's pass-name filter
    has a load-bearing baseline."""
    no_shm = ga.Experiment.on("human_igh").recombine().run_records(n=2, seed=0)
    for r in no_shm:
        assert r["n_mutations"] == 0

    pcr_only = (
        ga.Experiment.on("human_igh")
        .recombine()
        .pcr_amplify(count=10)
        .run_records(n=2, seed=0)
    )
    for r in pcr_only:
        assert r["n_mutations"] == 0, (
            "pcr_amplify bumped n_mutations; biological-vs-artefact "
            "split regressed."
        )

    quality_only = (
        ga.Experiment.on("human_igh")
        .recombine()
        .sequencing_errors(count=8)
        .run_records(n=2, seed=0)
    )
    for r in quality_only:
        assert r["n_mutations"] == 0


def test_pin_scaffold_receptor_revision_and_d_inversion_do_not_bump_n_mutations() -> None:
    """Audit §2: receptor revision rewrites the V slice, D
    inversion changes orientation — neither counts as SHM."""
    revision = (
        ga.Experiment.on("human_igh")
        .recombine()
        .receptor_revision(prob=1.0)
        .run_records(n=3, seed=0)
    )
    for r in revision:
        assert r["n_mutations"] == 0

    inversion = (
        ga.Experiment.on("human_igh")
        .recombine()
        .invert_d(prob=1.0)
        .run_records(n=3, seed=0)
    )
    for r in inversion:
        assert r["n_mutations"] == 0


# ──────────────────────────────────────────────────────────────────
# 5. Scaffold — segment-rate zero exclusion behavioural baseline
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_segment_rate_np_zero_keeps_global_count_intact() -> None:
    """Audit §8 + Slice 1: with ``segment_rates={"NP": 0.0}``,
    the realised SHM count is unchanged (the pass still applies
    ``count`` mutations, just confined to V/D/J). Slice 1 landed
    the ``n_np_mutations`` field; the load-bearing pin now reads
    the value directly and asserts it's ``0`` under NP=0
    targeting."""
    result = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(count=15, segment_rates={"NP": 0.0})
        .run_records(n=3, seed=0)
    )
    for r in result:
        assert r["n_mutations"] == 15
        # The Slice 1 field is now present and pinned to 0 under
        # NP=0 targeting (the audit's load-bearing claim is now
        # directly observable).
        assert r["n_np_mutations"] == 0
        # Sum invariant: V + D + J + NP == n_mutations.
        total = (
            r["n_v_mutations"]
            + r["n_d_mutations"]
            + r["n_j_mutations"]
            + r["n_np_mutations"]
        )
        assert total == 15


def test_pin_scaffold_segment_rate_v_only_keeps_global_count_intact() -> None:
    """Audit §8 sister pin: ``{D:0, J:0, NP:0}`` keeps
    ``n_mutations`` at the requested ``count`` — V-only targeting
    doesn't shrink the realised count."""
    result = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(count=20, segment_rates={"D": 0.0, "J": 0.0, "NP": 0.0})
        .run_records(n=3, seed=0)
    )
    for r in result:
        assert r["n_mutations"] == 20


# ──────────────────────────────────────────────────────────────────
# 6. Scaffold — PCR/quality counters are TRACE-sourced (attempted)
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_pcr_quality_counters_are_trace_sourced() -> None:
    """Audit §1: ``n_pcr_errors`` and ``n_quality_errors`` are
    trace-sourced — they reflect the *attempted* draw, not the
    realised event count. Pin via source inspection so the
    distinction stays explicit when the per-segment SHM slice
    lands (which uses event sourcing instead)."""
    builder_src = (
        _REPO_ROOT
        / "engine_rs"
        / "src"
        / "airr_record"
        / "builder.rs"
    ).read_text(encoding="utf-8")
    assert "trace_int_choice(trace, ChoiceAddress::CorruptPcrCount)" in builder_src, (
        "n_pcr_errors no longer sources from the trace; the future "
        "SHM counter slice's event-sourcing rationale (audit §1) "
        "has drifted — re-check the documentation."
    )
    assert "trace_int_choice(trace, ChoiceAddress::CorruptQualityCount)" in builder_src, (
        "n_quality_errors no longer trace-sourced."
    )


# ──────────────────────────────────────────────────────────────────
# 7. Absence — no per-segment SHM counters today
# ──────────────────────────────────────────────────────────────────


def test_pin_present_per_segment_mutation_counters_in_record() -> None:
    """Slice 1 shipped — the four AIRR fields are present and
    satisfy the sum invariant. Flipped from
    ``pin_absence_no_per_segment_mutation_counters_in_record``.

    Spec-driven behavioural coverage lives in
    [`tests/test_per_segment_mutation_counters.py`](test_per_segment_mutation_counters.py);
    this pin is the audit-doc lockstep counterpart. Per-NP1/NP2
    split stays deferred — the audit §13 out-of-scope pin
    documenting it remains as a separate test."""
    result = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(count=10)
        .run_records(n=1, seed=0)
    )
    rec = result.records[0]
    for field in (
        "n_v_mutations",
        "n_d_mutations",
        "n_j_mutations",
        "n_np_mutations",
    ):
        assert field in rec, (
            f"AIRR record missing {field!r}; Slice 1 backed out."
        )
    # Sum invariant — the load-bearing pin.
    total = (
        rec["n_v_mutations"]
        + rec["n_d_mutations"]
        + rec["n_j_mutations"]
        + rec["n_np_mutations"]
    )
    assert total == rec["n_mutations"], (
        f"sum invariant broken: {total} != {rec['n_mutations']}"
    )


def test_pin_absence_no_per_np_split_counters_in_record() -> None:
    """Per-NP1 / per-NP2 split stays out of scope (audit §13).
    NP1+NP2 events roll into the single ``n_np_mutations``
    bucket; the split-bucket fields don't exist."""
    result = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(count=10)
        .run_records(n=1, seed=0)
    )
    rec = result.records[0]
    for forbidden in ("n_np1_mutations", "n_np2_mutations"):
        assert forbidden not in rec, (
            f"AIRR record now carries {forbidden!r}; per-NP split "
            "slice has landed — flip pin in lockstep."
        )


def test_pin_present_per_segment_mutation_counter_default_columns() -> None:
    """Slice 1 shipped — the four fields appear in
    ``_DEFAULT_COLUMN_ORDER`` so TSV / CSV / DataFrame exports
    surface them deterministically. Flipped from
    ``pin_absence_no_per_segment_mutation_counter_default_columns``."""
    from GenAIRR.result import _DEFAULT_COLUMN_ORDER

    for field in (
        "n_v_mutations",
        "n_d_mutations",
        "n_j_mutations",
        "n_np_mutations",
    ):
        assert field in _DEFAULT_COLUMN_ORDER, (
            f"_DEFAULT_COLUMN_ORDER missing {field!r}; Slice 1 "
            "backed out at the Python column-order layer."
        )


# ──────────────────────────────────────────────────────────────────
# 8. Absence — no validator issue kinds for mutation-count mismatch
# ──────────────────────────────────────────────────────────────────


def test_pin_present_mutation_count_validator_kinds() -> None:
    """Slice 1 shipped — five new validator issue kinds added to
    the Rust validator: ``NVMutationsMismatch``,
    ``NDMutationsMismatch``, ``NJMutationsMismatch``,
    ``NNpMutationsMismatch``, ``MutationCountSumMismatch``.

    Flipped from
    ``pin_absence_no_mutation_count_validator_kinds`` when the
    slice landed. The variants follow the existing
    ``N*IndelsMismatch`` naming convention."""
    validate_src = (
        _REPO_ROOT
        / "engine_rs"
        / "src"
        / "airr_record"
        / "validate.rs"
    ).read_text(encoding="utf-8")
    for kind in (
        "NVMutationsMismatch",
        "NDMutationsMismatch",
        "NJMutationsMismatch",
        "NNpMutationsMismatch",
        "MutationCountSumMismatch",
    ):
        assert kind in validate_src, (
            f"validate.rs missing issue kind {kind!r}; Slice 1 "
            "backed out at the validator layer."
        )


# ──────────────────────────────────────────────────────────────────
# 9. Absence — no manifest field for per-segment counter support
# ──────────────────────────────────────────────────────────────────


def test_pin_absence_no_per_segment_counters_field_in_manifest() -> None:
    """Audit §10: the manifest's ``models.shm`` block doesn't
    advertise the per-segment counter capability yet. Slice 1's
    optional manifest extension would surface here."""
    shm = ga.HUMAN_IGH_OGRDB.cartridge_manifest()["models"]["shm"]
    for forbidden in (
        "per_segment_counters",
        "per_segment_counter_support",
        "shm_counters",
    ):
        assert forbidden not in shm, (
            f"manifest.models.shm now carries {forbidden!r}; flip "
            "pin."
        )


# ──────────────────────────────────────────────────────────────────
# 10. Absence — Python EventRecord lacks typed-event list
# ──────────────────────────────────────────────────────────────────


def test_pin_absence_python_eventrecord_lacks_typed_events() -> None:
    """Audit §4 / §11 out-of-scope item: the Python ``PyEventRecord``
    exposes ``simulation_event_count`` but not the typed event
    list. The proposed counter slice doesn't need it (Rust
    aggregates), but pin the gap so a future "where did each
    mutation land" surface lands deliberately."""
    result = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(count=5)
        .run_records(n=1, seed=0)
    )
    assert result.outcomes is not None
    event_record = result.outcomes[0].events()[0]
    # The count accessor is present.
    assert hasattr(event_record, "simulation_event_count")
    # But the typed event list isn't.
    for forbidden in (
        "simulation_events",
        "events_list",
        "typed_events",
    ):
        assert not hasattr(event_record, forbidden), (
            f"EventRecord.{forbidden} now exists; the Python typed-"
            "event observability gap closed — flip pin."
        )


# ──────────────────────────────────────────────────────────────────
# 11. Doc anchor — audit doc references contract file
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_audit_doc_exists_and_references_contract() -> None:
    """The audit doc must continue to exist and reference this
    contract file; the 14-section structure stays intact."""
    doc_path = _REPO_ROOT / "docs" / "mutation_provenance_audit.md"
    assert doc_path.exists(), "mutation_provenance_audit.md missing"
    doc = doc_path.read_text(encoding="utf-8")
    assert "test_mutation_provenance_contract.py" in doc, (
        "audit doc no longer references the contract file; lockstep "
        "convention drifted."
    )
    for marker in (
        "## 1. Q1",
        "## 5. Q5",
        "## 12. Test surface",
        "## 14. Summary table",
    ):
        assert marker in doc, f"audit doc missing section marker {marker!r}"
