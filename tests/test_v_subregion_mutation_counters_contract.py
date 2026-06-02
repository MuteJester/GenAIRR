"""Contract pins for the V-Subregion Mutation Counters audit.

Companion to
[`docs/v_subregion_mutation_counters_audit.md`](../docs/v_subregion_mutation_counters_audit.md).
The audit is pre-implementation: it freezes today's surfaces
(``pin_scaffold_*``) and the gaps a future implementation slice
would close (``pin_absence_*``).

This audit is the natural follow-up to the V-subregion SHM
**rate** slice — Slice B made per-subregion SHM **targeting**
shippable; the counters slice exposes the **realised**
per-subregion SHM counts on the projected AIRR record.

Six new AIRR fields are proposed:

- ``n_fwr1_mutations`` / ``n_cdr1_mutations`` /
  ``n_fwr2_mutations`` / ``n_cdr2_mutations`` /
  ``n_fwr3_mutations`` — the five canonical IMGT V subregion
  buckets.
- ``n_v_unannotated_mutations`` — V SHM events that can't be
  attributed to a subregion (V-side CDR3 stretch, unannotated V
  alleles in mixed cartridges, or indel-inserted V bases when
  the user inverts the canonical SHM-before-indels order).

Together the six partition the existing ``n_v_mutations``:

    n_fwr1 + n_cdr1 + n_fwr2 + n_cdr2 + n_fwr3
    + n_v_unannotated == n_v_mutations

Split:

- ``pin_scaffold_*`` tests freeze the pre-existing surfaces the
  slice builds on (BaseChanged carries germline_pos; the
  per-segment counter aggregation loop exists; the validator's
  independent-recompute pattern exists; bundled cartridges
  have 100% V-subregion coverage; the per-segment sum
  invariant already holds).
- ``pin_absence_*`` tests freeze the gaps the slice closes
  (no new AIRR fields, no manifest block, no validator issue
  kinds).
"""
from __future__ import annotations

import inspect
import json
from pathlib import Path

import pytest

import GenAIRR as ga


_REPO_ROOT = Path(__file__).resolve().parent.parent
_AUDIT_DOC = _REPO_ROOT / "docs" / "v_subregion_mutation_counters_audit.md"

# The six AIRR fields the implementation slice will add.
_SUBREGION_COUNTER_FIELDS: tuple[str, ...] = (
    "n_fwr1_mutations",
    "n_cdr1_mutations",
    "n_fwr2_mutations",
    "n_cdr2_mutations",
    "n_fwr3_mutations",
    "n_v_unannotated_mutations",
)


# ──────────────────────────────────────────────────────────────────
# 1. Scaffold — BaseChanged carries germline_pos: Option<u16>
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_base_changed_carries_germline_pos() -> None:
    """Audit §2 / §3 — the central feasibility surface. The
    ``BaseChanged`` event variant in ``sim_event.rs`` exposes a
    ``germline_pos: Option<u16>`` field. This is the load-bearing
    payload the counters slice depends on for V→subregion
    attribution. Pinned at source so a refactor that changes the
    field shape (or removes it) surfaces here."""
    sim_event_src = (
        _REPO_ROOT / "engine_rs" / "src" / "ir" / "sim_event.rs"
    ).read_text(encoding="utf-8")
    import re

    # The variant block: `BaseChanged { ... }`.
    match = re.search(
        r"BaseChanged\s*\{(.*?)\}",
        sim_event_src,
        re.DOTALL,
    )
    assert match is not None, (
        "BaseChanged variant not found in sim_event.rs; counters "
        "audit reference drifted"
    )
    body = match.group(1)
    for field in (
        "handle:",
        "old_base:",
        "new_base:",
        "segment:",
        "germline_pos:",
    ):
        assert field in body, (
            f"BaseChanged variant missing {field!r}; counters audit "
            "central feasibility argument regressed"
        )
    # The shape is specifically Option<u16> — the audit relies on the
    # None sentinel for indel-inserted V bases.
    assert "germline_pos: Option<u16>" in body, (
        "BaseChanged.germline_pos shape changed from Option<u16>; "
        "the audit's unannotated-bucket discipline depends on the "
        "None sentinel"
    )


# ──────────────────────────────────────────────────────────────────
# 2. Scaffold — assembly stamps germline_pos as allele-local
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_assembly_stamps_allele_local_germline_pos() -> None:
    """Audit §3 — assembly stamps ``germline_pos = slice_start +
    i`` where ``slice_start = trim_5``. This pre-computes the
    pool-to-allele arithmetic so the counters slice can attribute
    events without rerunning ``v_subregion_at_position``-style
    coordinate translation. Pinned via source-level inspection."""
    src = (
        _REPO_ROOT
        / "engine_rs"
        / "src"
        / "passes"
        / "assemble_segment"
        / "execution.rs"
    ).read_text(encoding="utf-8")
    # The assembly loop computes `slice_start + i` per nucleotide
    # and stamps it via `Nucleotide::germline(...)`. Pin both
    # markers.
    assert "slice_start" in src, (
        "assemble_segment execution no longer uses slice_start; the "
        "assembly-time stamping of germline_pos may have changed"
    )
    assert "Nucleotide::germline" in src, (
        "assembly no longer calls Nucleotide::germline; the "
        "counters audit's allele-local-pos invariant regressed"
    )


# ──────────────────────────────────────────────────────────────────
# 3. Scaffold — per-segment counter aggregation precedent
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_per_segment_counter_aggregation_loop_exists() -> None:
    """Audit §2 — the existing per-segment counter loop in
    ``airr_record/builder.rs`` is the precedent the new counters
    aggregation mirrors. Pin the loop body's load-bearing
    ingredients: the pass-name allowlist filter, the
    ``BaseChanged`` match, and the four bucket assignments."""
    src = (
        _REPO_ROOT / "engine_rs" / "src" / "airr_record" / "builder.rs"
    ).read_text(encoding="utf-8")
    # Pass-name allowlist filter — the counters slice must reuse it.
    assert "MUTATE_UNIFORM" in src
    assert "MUTATE_S5F" in src
    # BaseChanged match in the aggregation loop.
    assert "SimulationEvent::BaseChanged" in src
    # Four buckets — the slice extends Segment::V's arm with
    # subregion dispatch.
    for bucket in (
        "n_v_mutations",
        "n_d_mutations",
        "n_j_mutations",
        "n_np_mutations",
    ):
        assert bucket in src, (
            f"per-segment counter aggregation missing {bucket}; "
            "the counters audit's precedent regressed"
        )


def test_pin_scaffold_pass_name_allowlist_is_two_constants() -> None:
    """Audit §2 / §4 — the pass-name allowlist for biological SHM
    counts is exactly two constants (``MUTATE_UNIFORM`` /
    ``MUTATE_S5F``) defined in ``address.rs``. The counters slice
    must reuse the same filter so PCR / quality / receptor-revision
    base changes don't leak into the new V-subregion buckets."""
    src = (
        _REPO_ROOT / "engine_rs" / "src" / "address.rs"
    ).read_text(encoding="utf-8")
    assert 'pub const MUTATE_UNIFORM: &str = "mutate.uniform"' in src
    assert 'pub const MUTATE_S5F: &str = "mutate.s5f"' in src


# ──────────────────────────────────────────────────────────────────
# 4. Scaffold — validator does an independent recompute
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_validator_recomputes_independently() -> None:
    """Audit §5 — the existing per-segment counter validator
    re-walks ``outcome.events()`` from scratch (NOT from the AIRR
    record's reported field) and applies the same pass-name
    allowlist. This is what gives the sum invariant genuine bite.
    The new V-subregion validator extension inherits the same
    discipline."""
    src = (
        _REPO_ROOT / "engine_rs" / "src" / "airr_record" / "validate.rs"
    ).read_text(encoding="utf-8")
    # The four existing per-segment mismatch issue kinds.
    for kind in (
        "NVMutationsMismatch",
        "NDMutationsMismatch",
        "NJMutationsMismatch",
        "NNpMutationsMismatch",
        "MutationCountSumMismatch",
    ):
        assert kind in src, (
            f"validator missing {kind!r}; the counters audit's "
            "precedent regressed"
        )
    # The validator's recompute loop uses the same allowlist.
    assert "MUTATE_UNIFORM" in src
    assert "MUTATE_S5F" in src


# ──────────────────────────────────────────────────────────────────
# 5. Scaffold — bundled cartridges have 100% V-subregion coverage
# ──────────────────────────────────────────────────────────────────


@pytest.mark.parametrize(
    "name",
    ["HUMAN_IGH_OGRDB", "HUMAN_IGK_OGRDB", "HUMAN_IGL_OGRDB"],
)
def test_pin_scaffold_bundled_v_has_full_subregion_coverage(name: str) -> None:
    """Slice 1 prerequisite. The counters slice promises that
    ``n_v_unannotated_mutations == 0`` on every bundled-cartridge
    record under the canonical pass order; that promise depends
    on 100% subregion annotation coverage on the bundled V pool.
    Pinned here so a cartridge regression would make the future
    counter expected-zero invariant impossible."""
    cfg = getattr(ga, name)
    sup = cfg.cartridge_manifest()["models"]["shm"]["v_subregion_support"]
    assert sup["available"] is True
    assert sup["annotated_v_count"] == sup["total_v_count"] > 0


# ──────────────────────────────────────────────────────────────────
# 6. Scaffold — per-segment sum invariant already holds
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_per_segment_counters_partition_n_mutations() -> None:
    """Audit §7 — the existing per-segment sum invariant
    ``n_v + n_d + n_j + n_np == n_mutations`` is the partition the
    new V-subregion counters live UNDER. The new partition
    discipline is ``n_fwr1 + n_cdr1 + n_fwr2 + n_cdr2 + n_fwr3 +
    n_v_unannotated == n_v_mutations``; both invariants must hold
    simultaneously on every record. Pin the existing invariant so
    the audit's two-layer partition rests on a verified base."""
    result = (
        ga.Experiment.on("human_igh")
        .recombine()
        .productive_only()
        .mutate(model="s5f", rate=0.04)
        .run_records(n=15, seed=4242)
    )
    for r in result:
        total = (
            r["n_v_mutations"]
            + r["n_d_mutations"]
            + r["n_j_mutations"]
            + r["n_np_mutations"]
        )
        assert total == r["n_mutations"], (
            f"per-segment sum invariant violated: "
            f"V={r['n_v_mutations']} + D={r['n_d_mutations']} + "
            f"J={r['n_j_mutations']} + NP={r['n_np_mutations']} = "
            f"{total} != n_mutations={r['n_mutations']}"
        )


# ──────────────────────────────────────────────────────────────────
# 7. Scaffold — V allele identity reachable at projection time
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_v_allele_identity_reachable_at_aggregation() -> None:
    """Audit §2 / §11 — the new aggregation needs to reach the
    assigned V allele's subregion table once per record (hoisted
    outside the event loop). The access path
    ``sim.assignments.get(Segment::V) → AlleleInstance →
    refdata.v_pool.get(allele_id) → allele.subregions`` is the
    same one the rate slice's ``v_subregion_at_position`` uses.
    Pin both the assignment field on Simulation and the
    AlleleInstance shape."""
    sim_src = (
        _REPO_ROOT / "engine_rs" / "src" / "ir" / "simulation.rs"
    ).read_text(encoding="utf-8")
    assignment_src = (
        _REPO_ROOT / "engine_rs" / "src" / "assignment.rs"
    ).read_text(encoding="utf-8")
    # Simulation carries the assignments field (the projection
    # entry point uses outcome.final_simulation()).
    assert "pub assignments:" in sim_src and "AlleleAssignments" in sim_src
    # AlleleInstance exposes allele_id (indexes v_pool) + trim_5.
    assert "pub allele_id: AlleleId" in assignment_src
    assert "pub trim_5: u16" in assignment_src


# ──────────────────────────────────────────────────────────────────
# 8. Scaffold — PCR / quality counters surface separately
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_pcr_and_quality_errors_surface_separately() -> None:
    """Audit §2 — PCR / quality base changes are corruption-stage
    events, NOT biological SHM. They surface via dedicated
    ``n_pcr_errors`` / ``n_quality_errors`` fields and do NOT
    enter the per-segment SHM counters. The new V-subregion
    counters inherit the same filter: PCR / quality substitutions
    in V do NOT contribute to the new buckets. Pin the existence
    of the dedicated fields so a refactor that collapses them
    surfaces."""
    result = (
        ga.Experiment.on("human_igh")
        .recombine()
        .productive_only()
        .mutate(model="s5f", rate=0.03)
        .pcr_amplify(count=5)
        .sequencing_errors(count=3)
        .run_records(n=5, seed=4242)
    )
    rec = result.records[0]
    assert "n_pcr_errors" in rec, (
        "n_pcr_errors field missing; the audit's PCR-isolation "
        "argument depends on this dedicated counter"
    )
    assert "n_quality_errors" in rec, (
        "n_quality_errors field missing; the audit's "
        "quality-isolation argument depends on this dedicated "
        "counter"
    )


# ──────────────────────────────────────────────────────────────────
# 9. Absence — no five-label V subregion counter fields
# ──────────────────────────────────────────────────────────────────


def test_pin_present_v_subregion_counter_fields() -> None:
    """Audit §6 — Slice landed. All six new V-subregion counter
    fields (``n_fwr1_mutations`` / ``n_cdr1_mutations`` /
    ``n_fwr2_mutations`` / ``n_cdr2_mutations`` /
    ``n_fwr3_mutations`` / ``n_v_unannotated_mutations``) are now
    present on every AIRR record. Two-bucket aggregates
    (``n_cdr_mutations`` / ``n_fwr_mutations``) remain
    deliberately absent — per audit §1, derivable downstream."""
    result = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(model="s5f", rate=0.03)
        .run_records(n=1, seed=0)
    )
    rec = result.records[0]
    for required in _SUBREGION_COUNTER_FIELDS:
        assert required in rec, (
            f"AIRR record is missing {required!r}; V-subregion "
            "mutation counters slice regressed"
        )
        assert isinstance(rec[required], int), (
            f"{required!r} must be an int counter; got "
            f"{type(rec[required]).__name__}"
        )
        assert rec[required] >= 0, (
            f"{required!r} must be non-negative; got {rec[required]}"
        )
    # Two-bucket aggregates are deliberately NOT in v1.
    for forbidden in ("n_cdr_mutations", "n_fwr_mutations"):
        assert forbidden not in rec, (
            f"AIRR record now carries the two-bucket aggregate "
            f"{forbidden!r}; audit §1 recommended five-label fields "
            "only — flip this pin separately if a future ergonomics "
            "slice adds aggregates"
        )


# ──────────────────────────────────────────────────────────────────
# 10. Absence — no n_v_unannotated_mutations field
# ──────────────────────────────────────────────────────────────────


def test_pin_present_n_v_unannotated_mutations_field() -> None:
    """Audit §4 — Slice landed. The partition-preserving field
    ``n_v_unannotated_mutations`` exists on every AIRR record.
    The other (alternative-name) variants stay forbidden so the
    canonical surface stays one name."""
    result = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(model="s5f", rate=0.03)
        .run_records(n=1, seed=0)
    )
    rec = result.records[0]
    assert "n_v_unannotated_mutations" in rec, (
        "n_v_unannotated_mutations missing; partition discipline "
        "regressed"
    )
    for forbidden in ("n_v_unmapped_mutations", "n_v_other_mutations"):
        assert forbidden not in rec, (
            f"AIRR record now carries alternate-name {forbidden!r}; "
            "only ``n_v_unannotated_mutations`` is the canonical surface"
        )


# ──────────────────────────────────────────────────────────────────
# 11. Absence — no v_subregion_counter_support in manifest
# ──────────────────────────────────────────────────────────────────


def test_pin_present_v_subregion_counter_support_in_manifest() -> None:
    """Audit §11 — Slice landed. The manifest now carries the
    ``v_subregion_counter_support`` block alongside the
    annotation (Slice 1) and rate (Slice B) blocks. The block
    advertises the six fields, the partition-of relationship to
    ``n_v_mutations``, and the ``in_content_hash=False`` boundary
    (counters are per-record observables, not cartridge identity)."""
    m = ga.HUMAN_IGH_OGRDB.cartridge_manifest()
    shm = m["models"]["shm"]
    # Prerequisite blocks present.
    assert "v_subregion_support" in shm
    assert "v_subregion_rate_support" in shm
    # New counter capability block.
    assert "v_subregion_counter_support" in shm, (
        "v_subregion_counter_support block missing from manifest; "
        "Slice surface regressed"
    )
    sup = shm["v_subregion_counter_support"]
    assert sup["available"] is True
    assert sup["fields"] == [
        "n_fwr1_mutations",
        "n_cdr1_mutations",
        "n_fwr2_mutations",
        "n_cdr2_mutations",
        "n_fwr3_mutations",
        "n_v_unannotated_mutations",
    ]
    assert sup["partition_of"] == "n_v_mutations"
    assert sup["requires_annotations"] is True
    assert sup["unannotated_bucket"] == "n_v_unannotated_mutations"
    assert sup["in_content_hash"] is False
    # Alternate-name variants stay absent.
    blob = json.dumps(shm)
    for forbidden in ("v_subregion_count_support",):
        assert forbidden not in blob, (
            f"manifest now mentions alternate-name {forbidden!r}; "
            "only ``v_subregion_counter_support`` is the canonical "
            "manifest key"
        )


# ──────────────────────────────────────────────────────────────────
# 12. Absence — no V-subregion mismatch validator kinds
# ──────────────────────────────────────────────────────────────────


def test_pin_present_v_subregion_mismatch_validator_kinds() -> None:
    """Audit §5 — Slice landed. The Rust validator now carries
    seven new issue kinds: six per-field
    ``N{Fwr1,Cdr1,Fwr2,Cdr2,Fwr3,VUnannotated}MutationsMismatch``
    plus the cross-field
    ``VSubregionMutationCountSumMismatch`` invariant. Pinned at
    source so a refactor that drops one of them surfaces here."""
    validate_src = (
        _REPO_ROOT / "engine_rs" / "src" / "airr_record" / "validate.rs"
    ).read_text(encoding="utf-8")
    for required in (
        "NFwr1MutationsMismatch",
        "NCdr1MutationsMismatch",
        "NFwr2MutationsMismatch",
        "NCdr2MutationsMismatch",
        "NFwr3MutationsMismatch",
        "NVUnannotatedMutationsMismatch",
        "VSubregionMutationCountSumMismatch",
    ):
        assert required in validate_src, (
            f"validate.rs is missing {required!r}; V-subregion "
            "counter validator regressed"
        )


# ──────────────────────────────────────────────────────────────────
# 13. Absence — no new trace addresses introduced
# ──────────────────────────────────────────────────────────────────


def test_pin_absence_no_new_mutate_trace_addresses() -> None:
    """Audit §10 — the counters slice is pure projection. It
    must NOT introduce any new ``ChoiceAddress`` variants — trace
    files produced before the slice continue to replay
    byte-identically. Pinned at the address vocabulary so a
    refactor that adds a per-subregion address surfaces here."""
    address_src = (
        _REPO_ROOT / "engine_rs" / "src" / "address.rs"
    ).read_text(encoding="utf-8")
    for forbidden in (
        "MutateS5fSubregion",
        "MutateUniformSubregion",
        "MutateSubregionCount",
        "VSubregionMutation",
    ):
        assert forbidden not in address_src, (
            f"address.rs now carries {forbidden!r}; the counters "
            "audit promised no new trace addresses — flip pin and "
            "document the new address"
        )


# ──────────────────────────────────────────────────────────────────
# 14. Doc anchor — audit doc exists and references contract file
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_audit_doc_exists_and_references_contract() -> None:
    """The audit doc must continue to exist and reference this
    contract file; structure intact."""
    if not _AUDIT_DOC.exists():
        import pytest
        pytest.skip("docs/ is contributor-only; audit doc not present in this checkout")
    doc = _AUDIT_DOC.read_text(encoding="utf-8")
    assert "test_v_subregion_mutation_counters_contract.py" in doc, (
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
