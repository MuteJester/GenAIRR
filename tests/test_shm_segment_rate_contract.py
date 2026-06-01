"""Contract pins for the Per-Segment SHM Rate Scalars audit.

Companion to
[`docs/shm_segment_rate_design.md`](../docs/shm_segment_rate_design.md).
The audit is pre-implementation: it pins today's full-pool SHM
targeting and the architectural feasibility of the proposed slice
(`pin_scaffold_*`) plus the surfaces the slice will add
(`pin_absence_*`).

Split:

- ``pin_scaffold_*`` tests freeze today's contract: SHM mutates
  the full assembled pool; mutations land in V / NP1 / D / NP2 /
  J; the site→segment lookup is structurally determinable via
  ``Region.segment``; replay byte-deterministic; productive-only
  preserves triad; ``n_mutations`` is the global biological
  counter.
- ``pin_absence_*`` tests freeze the gaps the implementation
  slice closes: no ``segment_rates`` kwarg on
  ``Experiment.mutate``; no ``_segment_rates`` on the
  pipeline-IR ``_MutateStep``; no per-segment AIRR counter
  fields; no ``segment_rate_support`` in the cartridge manifest.

When the slice lands the relevant ``pin_absence_*`` flips to
``pin_present_*`` in lockstep.
"""
from __future__ import annotations

import inspect
import re
from pathlib import Path

import pytest

import GenAIRR as ga


_REPO_ROOT = Path(__file__).resolve().parent.parent


# ──────────────────────────────────────────────────────────────────
# 1. Scaffold — SHM targets the full assembled pool today
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_shm_mutations_land_in_every_region() -> None:
    """Audit §1: today's SHM has NO segment awareness — at high
    mutation count, substitutions land in V, NP1, D, NP2, J in
    proportion to each region's length. Pin a non-trivial number
    of differences per region across multiple records so the
    fact that "NP and D actually do receive SHM today" stays
    observable. Slice 1 will let users zero this out per
    segment."""
    # 100 records under count=50 — high enough that even short
    # NP/D regions land non-zero hits across the batch.
    base = (
        ga.Experiment.on("human_igh").recombine().run_records(n=20, seed=0)
    )
    mutated = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(count=50)
        .run_records(n=20, seed=0)
    )

    # Per-region differences accumulated across the batch.
    totals = {"V": 0, "NP1": 0, "D": 0, "NP2": 0, "J": 0}
    for b_rec, m_rec in zip(base, mutated):
        b_seq, m_seq = b_rec["sequence"].upper(), m_rec["sequence"].upper()
        v_end = m_rec["v_sequence_end"]
        d_start, d_end = m_rec["d_sequence_start"], m_rec["d_sequence_end"]
        j_start = m_rec["j_sequence_start"]
        bounds = {
            "V": (0, v_end),
            "NP1": (v_end, d_start),
            "D": (d_start, d_end),
            "NP2": (d_end, j_start),
            "J": (j_start, len(m_seq)),
        }
        for seg, (s, e) in bounds.items():
            totals[seg] += sum(
                1 for x, y in zip(b_seq[s:e], m_seq[s:e]) if x != y
            )

    # Every region must accumulate at least some SHM hits over the
    # batch. The fact that NP1/D/NP2/J are non-zero is the load-
    # bearing audit baseline — users today CANNOT prevent SHM from
    # landing here.
    for seg, total in totals.items():
        assert total > 0, (
            f"{seg} received zero SHM hits across the batch; full-pool "
            f"targeting baseline regressed (or fixture too thin). "
            f"totals={totals}"
        )


# ──────────────────────────────────────────────────────────────────
# 2. Scaffold — site→segment lookup is structurally determinable
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_region_struct_carries_segment_field() -> None:
    """Audit §1 pre-flight: ``Region.segment`` is the load-bearing
    field that makes a future per-segment rate slice feasible.
    Pin the source-level presence so a refactor that flattens
    `Region` would surface here before breaking the implementation
    plan."""
    region_src = (
        _REPO_ROOT / "engine_rs" / "src" / "ir" / "region.rs"
    ).read_text(encoding="utf-8")
    assert "pub segment: Segment" in region_src, (
        "Region.segment field missing; per-segment SHM rate slice "
        "would lose its site→segment lookup."
    )
    # And `Sequence.regions` must still hold the Vec.
    seq_src = (
        _REPO_ROOT / "engine_rs" / "src" / "ir" / "sequence.rs"
    ).read_text(encoding="utf-8")
    assert "pub regions: Vec<Region>" in seq_src, (
        "Sequence.regions Vec<Region> missing; the audit's site→"
        "segment lookup mechanism regressed."
    )


def test_pin_scaffold_segment_enum_has_five_variants() -> None:
    """Audit §2: the future segment_rates dict maps onto the
    Segment enum's five variants (V / Np1 / D / Np2 / J).
    NP1 + NP2 collapse into the single "NP" bucket per the
    recommendation. Pin the variant set so a Segment refactor
    that adds C or splits NP forces an explicit slice update."""
    seg_src = (
        _REPO_ROOT / "engine_rs" / "src" / "ir" / "segment.rs"
    ).read_text(encoding="utf-8")
    # Each variant's discriminant assignment.
    for variant in ("V = 0", "Np1 = 1", "D = 2", "Np2 = 3", "J = 4"):
        assert variant in seg_src, (
            f"Segment enum no longer carries {variant!r}; per-segment "
            "rate slice's bucket mapping needs updating."
        )
    assert "COUNT: usize = 5" in seg_src, (
        "Segment::COUNT changed; future per-segment rate surface's "
        "bucket count assumption needs updating."
    )


# ──────────────────────────────────────────────────────────────────
# 3. Scaffold — current mutate() lacks segment_rates kwarg
# ──────────────────────────────────────────────────────────────────


def test_pin_present_mutate_accepts_segment_rates() -> None:
    """Slice 1 shipped — ``Experiment.mutate`` now accepts the
    ``segment_rates: Optional[Dict[str, float]] = None`` kwarg.
    Flipped from ``pin_scaffold_mutate_signature_lacks_segment_rates_today``
    when the implementation slice landed. Spec-driven behavioural
    coverage lives in [`tests/test_segment_rates_implementation.py`](test_segment_rates_implementation.py)."""
    sig = inspect.signature(ga.Experiment.mutate)
    assert "segment_rates" in sig.parameters, (
        "Experiment.mutate.segment_rates regressed; Slice 1 backed out."
    )
    # Default must be ``None`` (sparse — flat 1.0 across buckets).
    assert sig.parameters["segment_rates"].default is None
    # Sanity: the existing kwargs are still there.
    for required_kwarg in ("model", "count", "rate", "s5f_model"):
        assert required_kwarg in sig.parameters


def test_pin_present_segment_rates_field_on_mutate_step() -> None:
    """Slice 1 shipped — ``_MutateStep`` carries the normalised
    4-tuple in canonical order ``(v, d, j, np)``. Flipped from
    ``pin_absence_no_segment_rates_field_on_mutate_step``."""
    from GenAIRR._pipeline_ir import _MutateStep

    assert "segment_rates" in _MutateStep.__dataclass_fields__, (
        "_MutateStep.segment_rates regressed; Slice 1 backed out."
    )
    # Default must be the flat-default tuple.
    default_value = _MutateStep.__dataclass_fields__["segment_rates"].default
    assert default_value == (1.0, 1.0, 1.0, 1.0)


# ──────────────────────────────────────────────────────────────────
# 4. Scaffold — replay byte-deterministic for SHM
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_shm_replay_byte_deterministic_today() -> None:
    """Audit §5: same-seed replay produces byte-identical
    sequences for both SHM models. The slice must preserve this
    under default rates (`{V:1, D:1, J:1, NP:1}`), and under
    any same-rate-vector replay. Re-pinned here from the SHM
    model contract for cross-doc traceability — the slice's
    "preserve under default rates" requirement reduces to this
    pin."""
    for model in ("uniform", "s5f"):
        exp = ga.Experiment.on("human_igh").recombine().mutate(
            model=model, count=10
        )
        a = exp.run_records(n=3, seed=42)
        b = exp.run_records(n=3, seed=42)
        assert [r["sequence"] for r in a] == [r["sequence"] for r in b], (
            f"model {model!r}: same-seed replay diverged."
        )


# ──────────────────────────────────────────────────────────────────
# 5. Scaffold — productive-only preserves triad under today's SHM
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_productive_only_preserves_triad_under_shm() -> None:
    """Audit §4: today's productive_only() + mutate(count=N)
    preserves the productive triad. The slice must preserve this
    under any non-degenerate segment_rates vector. Baseline
    re-pinned here for cross-doc traceability."""
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(count=20)
        .productive_only()
    )
    result = exp.run_records(n=5, seed=0)
    for r in result:
        assert r["productive"] is True


# ──────────────────────────────────────────────────────────────────
# 6. Scaffold — n_mutations is the global counter today
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_n_mutations_is_global_and_partitioned_today() -> None:
    """Audit §7: ``n_mutations`` is the canonical biological SHM
    counter on the AIRR record. After the per-segment counter
    slice landed (`docs/mutation_provenance_audit.md`), the four
    per-segment fields partition it by construction. Pin both
    the global field's continued presence AND the sum invariant
    so a future refactor that splits / rolls / drops a bucket
    surfaces here.

    Per-NP1/NP2 split deliberately stays out of scope (audit
    §13)."""
    result = ga.Experiment.on("human_igh").recombine().mutate(
        count=10
    ).run_records(n=1, seed=0)
    rec = result.records[0]
    assert "n_mutations" in rec
    assert rec["n_mutations"] == 10
    # Per-segment counters now present and satisfy the sum
    # invariant against the global counter.
    total = (
        rec["n_v_mutations"]
        + rec["n_d_mutations"]
        + rec["n_j_mutations"]
        + rec["n_np_mutations"]
    )
    assert total == 10
    # Per-NP1/NP2 split still deferred.
    for forbidden in ("n_np1_mutations", "n_np2_mutations"):
        assert forbidden not in rec, (
            f"AIRR record now carries {forbidden!r}; per-NP split "
            "slice landed — flip pin."
        )


# ──────────────────────────────────────────────────────────────────
# 7. Absence — no per-segment mutation counters
# ──────────────────────────────────────────────────────────────────


def test_pin_absence_no_per_segment_mutation_rate_fields() -> None:
    """Audit §13: per-segment mutation **rate** fields stay out of
    scope — derivable from ``n_<seg>_mutations / <seg>_length``
    and don't warrant their own AIRR columns. The per-segment
    counter slice that landed adds the four count fields
    (``n_v/d/j/np_mutations``) but not the four rate fields.

    Flipped from the original counter-absence pin when the
    mutation-provenance counter slice landed."""
    result = ga.Experiment.on("human_igh").recombine().mutate(
        count=5
    ).run_records(n=1, seed=0)
    rec = result.records[0]
    # Per-segment **count** fields are present (slice landed).
    assert "n_v_mutations" in rec
    # Per-segment **rate** fields stay absent.
    for forbidden in (
        "v_mutation_rate",
        "d_mutation_rate",
        "j_mutation_rate",
        "np_mutation_rate",
        "per_segment_mutation_rates",
    ):
        assert forbidden not in rec, (
            f"AIRR record now carries {forbidden!r}; rate-fields "
            "slice landed — flip pin."
        )


# ──────────────────────────────────────────────────────────────────
# 8. Absence — manifest doesn't advertise segment-rate support yet
# ──────────────────────────────────────────────────────────────────


def test_pin_present_segment_rate_support_in_manifest() -> None:
    """Slice 1 shipped — ``manifest["models"]["shm"]
    ["segment_rate_support"]`` advertises the capability with
    documented buckets + flat default + ``in_content_hash=False``
    boundary. Flipped from
    ``pin_absence_no_segment_rate_support_in_manifest``."""
    m = ga.HUMAN_IGH_OGRDB.cartridge_manifest()
    shm = m["models"]["shm"]
    assert "segment_rate_support" in shm, (
        "manifest.models.shm.segment_rate_support regressed; Slice 1 "
        "backed out."
    )
    sup = shm["segment_rate_support"]
    assert sup["available"] is True
    assert sup["buckets"] == ["V", "D", "J", "NP"]
    assert sup["default"] == {"V": 1.0, "D": 1.0, "J": 1.0, "NP": 1.0}
    # The slice does NOT change content_hash semantics; document
    # the v1 boundary explicitly.
    assert sup["in_content_hash"] is False


# ──────────────────────────────────────────────────────────────────
# 9. Absence — DSL boundary doesn't accept segment_rates yet
# ──────────────────────────────────────────────────────────────────


def test_pin_present_mutate_accepts_segment_rates_kwarg_at_runtime() -> None:
    """Slice 1 shipped — the runtime accepts the kwarg without
    raising on the validation path. Flipped from
    ``pin_absence_mutate_rejects_segment_rates_kwarg``."""
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(count=5, segment_rates={"V": 1.0, "NP": 0.0})
    )
    # And the pipeline compiles + runs end-to-end.
    result = exp.run_records(n=1, seed=0)
    assert len(result) == 1
    assert result[0]["n_mutations"] == 5


# ──────────────────────────────────────────────────────────────────
# 10. Doc anchor — audit doc references contract
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_audit_doc_exists_and_references_contract() -> None:
    """The audit doc must continue to exist and reference this
    contract file; the 14-section structure stays intact."""
    doc_path = _REPO_ROOT / "docs" / "shm_segment_rate_design.md"
    assert doc_path.exists(), "shm_segment_rate_design.md missing"
    doc = doc_path.read_text(encoding="utf-8")
    assert "test_shm_segment_rate_contract.py" in doc, (
        "audit doc no longer references the contract file; lockstep "
        "convention drifted."
    )
    for marker in (
        "## 1. Q1",
        "## 6. Q6",
        "## 12. Test surface",
        "## 14. Summary table",
    ):
        assert marker in doc, f"audit doc missing section marker {marker!r}"
