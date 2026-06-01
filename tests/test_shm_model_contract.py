"""Contract pins for the SHM Model / Context-Mutation audit.

Companion to
[`docs/shm_model_audit.md`](../docs/shm_model_audit.md). The audit
freezes the biological-vs-artefact split between SHM and the
corruption family, the IR-sourced `n_mutations` counter, the S5F
kernel's bundled-but-not-cartridge-owned status, the productive-
only interaction, and the validator + replay posture.

Split:

- ``pin_scaffold_*`` tests freeze the current contract: the two
  SHM passes, the single ``add_to_mutation_count`` source of
  truth, ``n_mutations`` IR-sourcing, full-pool targeting,
  productive-only preservation, truth-call stability, replay
  determinism, S5F kernel inventory.
- ``pin_absence_*`` tests freeze the gaps Slices 1-3 close: no
  SHM model metadata in the manifest, no S5F digest in
  ``content_hash``, no ``shm_model`` AIRR field, no per-segment
  mutation counters, no ``validate_mutation_ledger``, no
  mutation-position provenance AIRR field.

When a future slice closes a gap the relevant ``pin_absence_*``
flips to ``pin_present_*`` in lockstep.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

import GenAIRR as ga
from GenAIRR import _engine as ge


_REPO_ROOT = Path(__file__).resolve().parent.parent


def _refdata():
    return ga.Experiment.on("human_igh").compile().refdata


# ──────────────────────────────────────────────────────────────────
# 1. Scaffold — two SHM passes exist
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_two_shm_pass_implementations_exist() -> None:
    """Audit §1: two biological SHM passes — uniform and S5F.
    Pinned by source-level presence of the two files."""
    uniform_path = (
        _REPO_ROOT / "engine_rs" / "src" / "passes" / "mutate" / "uniform.rs"
    )
    s5f_path = (
        _REPO_ROOT / "engine_rs" / "src" / "passes" / "mutate" / "s5f.rs"
    )
    assert uniform_path.exists(), "UniformMutationPass module missing"
    assert s5f_path.exists(), "S5FMutationPass module missing"


# ──────────────────────────────────────────────────────────────────
# 2. Scaffold — only SHM passes call `add_to_mutation_count`
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_only_mutate_passes_call_add_to_mutation_count() -> None:
    """Audit §2: the single architectural source of truth for "is
    this mutation biological?" is the ``add_to_mutation_count``
    call. Grepping production passes returns exactly two hits:
    ``mutate/uniform.rs`` and ``mutate/s5f/execution.rs``. A third
    hit (e.g. from a new corruption pass) would surface a contract
    violation here."""
    passes_dir = _REPO_ROOT / "engine_rs" / "src" / "passes"
    hits = []
    for path in passes_dir.rglob("*.rs"):
        # Skip test / mutation_transaction infrastructure itself.
        if "test" in path.parts or "mutation_transaction" in path.parts:
            continue
        text = path.read_text(encoding="utf-8")
        # Match call sites only — `tx.add_to_mutation_count(...)`.
        for line in text.splitlines():
            stripped = line.strip()
            if stripped.startswith("//") or stripped.startswith("///"):
                continue
            if "tx.add_to_mutation_count(" in stripped:
                hits.append(path.relative_to(_REPO_ROOT).as_posix())
    # Dedupe per file (a pass may call it from multiple branches).
    hit_set = sorted(set(hits))
    assert hit_set == [
        "engine_rs/src/passes/mutate/s5f/execution.rs",
        "engine_rs/src/passes/mutate/uniform.rs",
    ], (
        f"unexpected add_to_mutation_count call sites: {hit_set}. "
        "A new biological mutation pass needs the call; a new "
        "corruption pass must NOT include it."
    )


# ──────────────────────────────────────────────────────────────────
# 3. Scaffold — corruption passes don't bump n_mutations
# ──────────────────────────────────────────────────────────────────


@pytest.mark.parametrize(
    "build_fn",
    [
        lambda e: e.pcr_amplify(count=10),
        lambda e: e.sequencing_errors(count=8),
        lambda e: e.ambiguous_base_calls(count=5),
        lambda e: e.polymerase_indels(count=2),
        lambda e: e.end_loss_5prime(length=10),
        lambda e: e.end_loss_3prime(length=10),
        lambda e: e.random_strand_orientation(prob=1.0),
    ],
)
def test_pin_scaffold_corruption_does_not_bump_n_mutations(build_fn) -> None:
    """Audit §2: every corruption pass produces ``n_mutations=0``
    when run alone. The corruption family fires ``BaseChanged`` /
    ``IndelInserted`` events but does NOT call
    ``add_to_mutation_count``."""
    exp = build_fn(ga.Experiment.on("human_igh").recombine())
    result = exp.run_records(n=2, seed=0)
    for r in result:
        assert r["n_mutations"] == 0, (
            f"corruption pass bumped n_mutations to {r['n_mutations']}; "
            "the biological/artefact split has been broken."
        )


# ──────────────────────────────────────────────────────────────────
# 4. Scaffold — SHM does bump n_mutations
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_uniform_mutate_bumps_n_mutations() -> None:
    """Audit §2: uniform mutate with ``count=N`` produces
    ``n_mutations == N`` for the realised count (constraint-
    rejected sites in permissive mode don't double-count, but
    without ``productive_only`` none are rejected)."""
    exp = ga.Experiment.on("human_igh").recombine().mutate(
        model="uniform", count=7
    )
    result = exp.run_records(n=3, seed=0)
    for r in result:
        assert r["n_mutations"] == 7


def test_pin_scaffold_s5f_mutate_bumps_n_mutations() -> None:
    """Audit §2: S5F mutate with ``count=N`` produces
    ``n_mutations == N`` (realised count)."""
    exp = ga.Experiment.on("human_igh").recombine().mutate(
        model="s5f", count=8
    )
    result = exp.run_records(n=3, seed=0)
    for r in result:
        assert r["n_mutations"] == 8


# ──────────────────────────────────────────────────────────────────
# 5. Scaffold — combined SHM + corruption tracks only SHM
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_combined_shm_plus_corruption_tracks_only_shm() -> None:
    """Audit §2: a pipeline stacking SHM + every corruption pass
    produces ``n_mutations`` equal to the SHM count only.
    Behavioural pin of the biological / artefact split."""
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(count=5)
        .pcr_amplify(count=10)
        .sequencing_errors(count=8)
        .ambiguous_base_calls(count=5)
    )
    result = exp.run_records(n=3, seed=0)
    for r in result:
        assert r["n_mutations"] == 5


# ──────────────────────────────────────────────────────────────────
# 6. Scaffold — mutation_rate is derived
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_mutation_rate_is_derived_from_n_mutations() -> None:
    """Audit §1 / builder.rs:213: ``mutation_rate == n_mutations /
    sequence_length`` exactly. Not a separate sample."""
    exp = ga.Experiment.on("human_igh").recombine().mutate(count=7)
    result = exp.run_records(n=3, seed=0)
    for r in result:
        seq_len = len(r["sequence"])
        expected = r["n_mutations"] / seq_len if seq_len > 0 else 0.0
        assert abs(r["mutation_rate"] - expected) < 1e-12


# ──────────────────────────────────────────────────────────────────
# 7. Scaffold — SHM targets the full assembled pool
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_shm_targets_full_assembled_pool() -> None:
    """Audit §4: SHM mutates the entire assembled pool (V + NP1 +
    D + NP2 + J), not just V. Behavioural pin: at high SHM count,
    differences relative to the unmutated baseline span the full
    sequence length."""
    base = ga.Experiment.on("human_igh").recombine().run_records(
        n=5, seed=0
    )
    mutated = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(count=30)
        .run_records(n=5, seed=0)
    )
    for i in range(5):
        base_seq = base[i]["sequence"].upper()
        mut_seq = mutated[i]["sequence"].upper()
        # Sequences differ.
        assert base_seq != mut_seq
        # Differences exist and are non-trivial (~ count = 30).
        # We don't pin a strict count — the realised count can
        # vary slightly under sentinel collisions — but the
        # observed count must be in a reasonable range.
        diff_count = sum(1 for a, b in zip(base_seq, mut_seq) if a != b)
        assert 15 <= diff_count <= 30, (
            f"record {i}: expected ~30 SHM differences, got {diff_count}"
        )


# ──────────────────────────────────────────────────────────────────
# 8. Scaffold — productive-only preserves triad under SHM
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_productive_only_preserves_triad_under_shm() -> None:
    """Audit §5: ``productive_only() + mutate(count=N)`` produces
    records with ``productive=True`` for every record even at
    high SHM count. Constraint filtering via the
    ``MutationTransaction`` ensures no substitution introduces a
    stop codon in junction or breaks an anchor."""
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(count=20)
        .productive_only()
    )
    result = exp.run_records(n=5, seed=0)
    for r in result:
        assert r["productive"] is True, (
            f"productive_only failed to preserve productive triad "
            f"at n_mutations={r['n_mutations']}"
        )


# ──────────────────────────────────────────────────────────────────
# 9. Scaffold — truth calls stable under SHM
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_truth_calls_stable_under_shm() -> None:
    """Audit §6: ``truth_v_call`` / ``truth_d_call`` /
    ``truth_j_call`` are IR-sourced and stable under any SHM
    count. Heavy SHM widens the live call but doesn't move the
    truth call."""
    exp = ga.Experiment.on("human_igh").recombine().mutate(count=20)
    result = exp.run_records(n=5, seed=0, expose_provenance=True)
    for r in result:
        # Truth calls always present + non-empty + single allele.
        for field in ("truth_v_call", "truth_d_call", "truth_j_call"):
            assert r[field], f"truth field {field!r} is empty"
            assert "," not in r[field], (
                f"truth field {field!r} is a tie-set string: {r[field]!r}"
            )


def test_pin_scaffold_live_call_truth_invariant_under_heavy_shm() -> None:
    """Audit §6: the live ``v_call`` may widen under heavy SHM,
    but the truth allele MUST appear in the tie set. The
    score-and-tie caller's design guarantees this."""
    exp = ga.Experiment.on("human_igh").recombine().mutate(count=25)
    result = exp.run_records(n=10, seed=0, expose_provenance=True)
    for r in result:
        truth = r["truth_v_call"]
        live = r["v_call"]
        assert truth, "truth_v_call empty"
        assert live, "v_call empty"
        assert truth in live.split(","), (
            f"truth_v_call={truth!r} dropped from v_call tie set {live!r}; "
            "the score-and-tie caller violated its invariant."
        )


# ──────────────────────────────────────────────────────────────────
# 10. Scaffold — validate_records re-derives n_mutations
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_validate_records_redrives_n_mutations() -> None:
    """Audit §7: ``validate_records`` re-derives the record from
    the outcome and compares. A clean SHM batch validates ok;
    dict-level tampering of ``n_mutations`` is silently ignored
    (documented behaviour — the supplied record dict isn't
    audited against the outcome). For real tamper detection,
    callers use the Rust-side ``validate_airr_record`` directly."""
    refdata = _refdata()
    result = (
        ga.Experiment.on("human_igh").recombine().mutate(count=5)
        .run_records(n=3, seed=0)
    )
    # Clean: validates.
    assert result.validate_records(refdata).ok
    # Tamper a record's n_mutations; validator still returns ok
    # because it re-projects from outcome state.
    result.records[0]["n_mutations"] = 9999
    assert result.validate_records(refdata).ok, (
        "validate_records started catching dict-level n_mutations "
        "tampering; the documented re-derive contract has drifted."
    )


# ──────────────────────────────────────────────────────────────────
# 11. Scaffold — replay byte-deterministic for SHM
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_shm_replay_byte_deterministic() -> None:
    """Audit §9: two runs of the same SHM-bearing experiment with
    the same seed produce byte-identical sequences. Pinned per
    SHM model so a refactor that introduces non-deterministic
    behavior surfaces here."""
    for model in ("uniform", "s5f"):
        exp = ga.Experiment.on("human_igh").recombine().mutate(
            model=model, count=10
        )
        a = exp.run_records(n=3, seed=42)
        b = exp.run_records(n=3, seed=42)
        assert [r["sequence"] for r in a] == [r["sequence"] for r in b], (
            f"SHM model {model!r}: same-seed replay produced different "
            "sequences."
        )


# ──────────────────────────────────────────────────────────────────
# 12. Scaffold — S5F bundled kernels documented
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_s5f_kernels_bundled_separately() -> None:
    """Audit §3: four canonical S5F kernels are bundled in
    ``_s5f_loader.py``. Pin the names so a refactor that drops
    one (or renames the loader) surfaces here. They are NOT in
    DataConfig — kernel selection happens at ``Experiment.mutate``
    time, not via the cartridge."""
    from GenAIRR._s5f_loader import _BUILTIN_S5F_MODELS

    expected = {"hh_s5f", "hh_s5f_60", "hh_s5f_opposite", "hkl_s5f"}
    assert set(_BUILTIN_S5F_MODELS.keys()) == expected


# ──────────────────────────────────────────────────────────────────
# 13. Scaffold — performance baseline covers SHM
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_performance_baseline_file_exists() -> None:
    """Audit §8: the existing performance-baseline test file
    covers SHM under realistic count / rate values. A future
    SHM-perf audit would extend it; this pin protects the
    file's existence."""
    perf = _REPO_ROOT / "tests" / "test_performance_budgets.py"
    assert perf.exists(), (
        "performance baseline test file missing; SHM performance "
        "coverage regressed."
    )


# ──────────────────────────────────────────────────────────────────
# 14. Absence — no SHM model metadata in cartridge manifest
# ──────────────────────────────────────────────────────────────────


def test_pin_present_shm_metadata_in_cartridge_manifest() -> None:
    """Slice 1 shipped — ``cartridge_manifest()["models"]["shm"]``
    surfaces the available SHM models, bundled S5F kernel names,
    the DSL default kernel, an optional digest of the default
    kernel's bundled bytes, and the explicit
    ``in_content_hash=False`` v1 boundary.

    Flipped from ``pin_absence_no_shm_model_in_cartridge_manifest``
    when the Slice landed. Spec-driven behavioural coverage lives
    in [`tests/test_cartridge_manifest_shm.py`](test_cartridge_manifest_shm.py);
    this pin is the audit-doc lockstep counterpart.

    The companion ``pin_absence_no_s5f_kernel_in_content_hash``
    stays absent — the manifest documents the v1 boundary but
    doesn't close it. A future slice that folds the kernel digest
    into ``content_hash`` would flip BOTH that absence pin AND
    update the ``in_content_hash`` field documented here."""
    m = ga.HUMAN_IGH_OGRDB.cartridge_manifest()
    shm = m.get("models", {}).get("shm")
    assert shm is not None, (
        "models.shm regressed; Slice 1 backed out."
    )
    # The five documented keys.
    assert "available_models" in shm
    assert "s5f_kernels_available" in shm
    assert "default_s5f_kernel" in shm
    assert "s5f_kernel_digest" in shm
    # The v1-boundary documentation field must stay False — this
    # pin protects the manifest's accurate description of the
    # ``content_hash`` semantics.
    assert shm["in_content_hash"] is False, (
        "manifest now claims SHM kernel is in content_hash; either "
        "the v1 boundary closed (flip the companion absence pin too) "
        "or the manifest is lying about the hash semantics."
    )


def test_pin_absence_no_s5f_kernel_in_content_hash() -> None:
    """Audit §3: swapping the S5F kernel doesn't change
    ``refdata.content_hash()`` because the kernel lives outside
    the bridge. The v1 boundary: same as ``reference_models``."""
    # The refdata content_hash is determined by the bridged
    # cartridge state, not by experiment-level configuration. So
    # building two experiments with the same DataConfig but
    # different S5F kernels still produces the same content hash.
    cfg = ga.HUMAN_IGH_OGRDB
    h1 = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(model="s5f", count=5, s5f_model="hh_s5f")
        .compile()
        .refdata
        .content_hash()
    )
    h2 = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(model="s5f", count=5, s5f_model="hkl_s5f")
        .compile()
        .refdata
        .content_hash()
    )
    assert h1 == h2, (
        "content_hash differs across S5F kernel choices; the v1 "
        "boundary closed — flip this pin in lockstep with the "
        "slice that landed kernel-in-identity."
    )


# ──────────────────────────────────────────────────────────────────
# 15. Absence — no `shm_model` AIRR field
# ──────────────────────────────────────────────────────────────────


def test_pin_absence_no_shm_model_airr_field() -> None:
    """Audit §11 Slice 1 follow-up: AIRR records don't carry a
    field naming which SHM model produced ``n_mutations``. A
    consumer comparing two runs must consult the experiment
    config, not the record dict."""
    result = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(model="s5f", count=5)
        .run_records(n=1, seed=0)
    )
    rec = result.records[0]
    for forbidden in (
        "shm_model",
        "mutation_model",
        "shm_kernel",
        "s5f_kernel_name",
    ):
        assert forbidden not in rec, (
            f"AIRR record now carries {forbidden!r}; flip pin."
        )


# ──────────────────────────────────────────────────────────────────
# 16. Absence — no per-segment mutation counters
# ──────────────────────────────────────────────────────────────────


def test_pin_present_per_segment_mutation_counters_with_canonical_names() -> None:
    """Mutation-provenance Slice 1 shipped (see
    ``docs/mutation_provenance_audit.md``) — the four canonical
    per-segment SHM counter fields landed under the
    ``n_<seg>_mutations`` naming convention.

    Flipped from the original SHM-audit
    ``pin_absence_no_per_segment_mutation_counters`` pin when
    the counter slice landed. Pin both the field presence AND
    the rejected alternative-naming forms so a refactor that
    renames or splits the buckets surfaces here. NP1/NP2 split
    and intra-segment (junction / CDR) counters remain out of
    scope."""
    result = ga.Experiment.on("human_igh").recombine().mutate(
        count=10
    ).run_records(n=1, seed=0)
    rec = result.records[0]
    # Global counter still present.
    assert "n_mutations" in rec
    # The four per-segment buckets landed under the canonical names.
    for field in (
        "n_v_mutations",
        "n_d_mutations",
        "n_j_mutations",
        "n_np_mutations",
    ):
        assert field in rec, (
            f"per-segment counter {field!r} regressed; Slice 1 backed out."
        )
    # Alternative naming forms that the audit rejected remain absent.
    for forbidden in (
        "v_n_mutations",
        "d_n_mutations",
        "j_n_mutations",
        "np1_n_mutations",
        "np2_n_mutations",
        "junction_n_mutations",
    ):
        assert forbidden not in rec, (
            f"AIRR record now carries {forbidden!r}; naming convention "
            "drifted — the slice landed under ``n_<seg>_mutations`` "
            "verbatim."
        )


# ──────────────────────────────────────────────────────────────────
# 17. Absence — no mutation-ledger validator
# ──────────────────────────────────────────────────────────────────


def test_pin_absence_no_mutation_ledger_validator() -> None:
    """Audit §7 / §11 Slice 3: no ``validate_mutation_ledger``
    method on ``Outcome`` or ``SimulationResult``. The per-site
    mutation-event vs pool-diff cross-check doesn't run at
    runtime; the ``MutationTransaction`` boundary lockdown
    provides static protection only."""
    result = ga.Experiment.on("human_igh").recombine().mutate(
        count=5
    ).run_records(n=1, seed=0)
    sample_outcome = result.outcomes[0]
    for forbidden in (
        "validate_mutation_ledger",
        "verify_mutation_positions",
    ):
        assert not hasattr(sample_outcome, forbidden), (
            f"Outcome.{forbidden} now exists; Slice 3 has landed — "
            "flip pin."
        )
    from GenAIRR.result import SimulationResult
    for forbidden in (
        "validate_mutation_ledger",
        "validate_mutation_positions",
    ):
        assert not hasattr(SimulationResult, forbidden), (
            f"SimulationResult.{forbidden} now exists; flip pin."
        )


# ──────────────────────────────────────────────────────────────────
# 18. Absence — no mutation-position provenance AIRR field
# ──────────────────────────────────────────────────────────────────


def test_pin_absence_no_mutation_position_provenance_field() -> None:
    """Audit §11 Slice 3 follow-up: AIRR records don't carry the
    list of mutated positions or per-mutation context info. A
    consumer must read ``outcome.events()`` and filter for
    ``mutate.{uniform,s5f}.site[i]`` records."""
    result = ga.Experiment.on("human_igh").recombine().mutate(
        count=5
    ).run_records(n=1, seed=0)
    rec = result.records[0]
    for forbidden in (
        "mutation_positions",
        "shm_positions",
        "mutation_sites",
        "mutated_positions",
        "mutation_context",
    ):
        assert forbidden not in rec, (
            f"AIRR record now carries {forbidden!r}; flip pin."
        )


# ──────────────────────────────────────────────────────────────────
# 19. Doc anchor — audit doc exists + references contract
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_audit_doc_exists_and_references_contract() -> None:
    """The audit doc must continue to exist and reference this
    contract file; the 14-section structure stays intact. Audit
    convention: a regression here means the change-control
    surface drifted."""
    doc_path = _REPO_ROOT / "docs" / "shm_model_audit.md"
    assert doc_path.exists(), "shm_model_audit.md missing"
    doc = doc_path.read_text(encoding="utf-8")
    assert "test_shm_model_contract.py" in doc, (
        "audit doc no longer references the contract file; lockstep "
        "convention drifted."
    )
    for marker in (
        "## 1. Q1",
        "## 5. Q5",
        "## 7. Q7",
        "## 12. Test surface",
        "## 14. Summary table",
    ):
        assert marker in doc, f"audit doc missing section marker {marker!r}"
