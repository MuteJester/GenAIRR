"""Contract pins for the Markov NP Base Generator slice.

Companion to
[`docs/np_markov_base_generator_design.md`](../docs/np_markov_base_generator_design.md).
The audit's implementation slice has shipped; pins below freeze
the **post-implementation** state:

- ``pin_scaffold_*`` tests still freeze pre-existing surfaces
  the slice built on (still passing) and now also the
  post-slice surfaces (new trait module, ``base_generator``
  field, manifest flip).
- ``pin_present_*`` tests freeze the post-slice behaviour
  (``kind="markov"`` lowers and runs end-to-end).

The original ``pin_absence_*`` set has been flipped to
``pin_present_*`` — each gap the implementation slice closed
is now pinned as **present** so a regression that removes any
of those surfaces surfaces here.
"""
from __future__ import annotations

import inspect
import json
from pathlib import Path

import pytest

import GenAIRR as ga
from GenAIRR.reference_models import (
    NpBaseModelSpec,
    ReferenceEmpiricalModels,
)


_REPO_ROOT = Path(__file__).resolve().parent.parent
_AUDIT_DOC = _REPO_ROOT / "docs" / "np_markov_base_generator_design.md"


# ──────────────────────────────────────────────────────────────────
# 1. Scaffold — execute_with_sampling_mode single entry point
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_execute_with_sampling_mode_single_entry() -> None:
    """Audit §1 — ``GenerateNPPass::execute_with_sampling_mode``
    is the only execution entry point. Branching to
    replay / fast / slow / unconstrained paths happens inside
    ``sample_base``. The Markov refactor modifies this one
    entry, not a separate fast path."""
    src = (
        _REPO_ROOT
        / "engine_rs"
        / "src"
        / "passes"
        / "generate_np"
        / "execution.rs"
    ).read_text(encoding="utf-8")
    assert "fn execute_with_sampling_mode" in src, (
        "execute_with_sampling_mode missing; the audit's single-entry "
        "assumption regressed"
    )
    # No alternate fast-RNG entry point.
    for forbidden in (
        "execute_fresh_rng_fast_path",
        "execute_without_contracts",
        "execute_uniform_fast_path",
    ):
        assert forbidden not in src, (
            f"execution.rs now carries {forbidden!r}; the audit's "
            "single-entry assumption changed — verify the Markov "
            "slice's loop-body change covers all paths"
        )


# ──────────────────────────────────────────────────────────────────
# 2. Scaffold — per-position loop carries no surviving state
# ──────────────────────────────────────────────────────────────────


def test_pin_present_per_position_loop_carries_previous_base() -> None:
    """Post-slice — the per-position loop body in
    ``execute_with_sampling_mode`` now carries
    ``let mut previous: Option<u8> = None;`` ahead of the
    loop and updates it after each base is accepted/recorded.
    Flipped from the prior scaffold pin that asserted the
    loop carried no per-iteration state."""
    src = (
        _REPO_ROOT
        / "engine_rs"
        / "src"
        / "passes"
        / "generate_np"
        / "execution.rs"
    ).read_text(encoding="utf-8")
    assert "for i in 0..length" in src
    assert "let mut previous: Option<u8>" in src, (
        "execution.rs no longer carries the `previous: Option<u8>` "
        "loop-local; the Markov slice's per-position state has "
        "regressed"
    )
    # And it's updated after `record_choice` so a mid-stream
    # replay failure can't leak partial state.
    assert "previous = Some(base);" in src


# ──────────────────────────────────────────────────────────────────
# 3. Scaffold — sample helpers are generic over Distribution
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_sample_helpers_consume_raw_pairs_via_support() -> None:
    """Audit §3 / §9 — the central coupling-clean claim.
    ``sample_base_with_admit_mask`` and
    ``sample_filtered_with_policy`` are generic over
    ``D: Distribution + ?Sized`` AND immediately materialise
    ``Vec<(u8, f64)>`` via ``dist.support()``, never calling
    ``dist.sample(rng)``. This is the clean refactor path the
    audit endorses."""
    src = (
        _REPO_ROOT / "engine_rs" / "src" / "dist" / "filtered.rs"
    ).read_text(encoding="utf-8")
    assert "fn sample_base_with_admit_mask" in src
    assert "fn sample_filtered_with_policy" in src
    # Both functions are generic over D: Distribution + ?Sized —
    # not concrete types.
    assert "D: Distribution<Output = u8> + ?Sized" in src
    assert "D: Distribution<Output = T> + ?Sized" in src


# ──────────────────────────────────────────────────────────────────
# 4. Scaffold — NP trace addresses unchanged under the proposal
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_np_trace_addresses_unchanged_under_markov_proposal() -> None:
    """Audit §2 / §10 — the existing per-position addresses
    carry the full Markov state (replay reconstructs
    previous-base from the prior recorded value at
    ``np.np1.bases[i-1]``). The Markov slice promises NO new
    addresses. Pin the existing address vocabulary."""
    src = (
        _REPO_ROOT / "engine_rs" / "src" / "address.rs"
    ).read_text(encoding="utf-8")
    for required in (
        '"np.np1.length"',
        '"np.np2.length"',
        '"np.np1.bases["',
        '"np.np2.bases["',
    ):
        assert required in src, (
            f"address.rs missing {required!r}; NP address vocabulary drifted"
        )
    # No Markov-specific addresses today.
    for forbidden in (
        "MarkovNpBase",
        "NpMarkovState",
        "NpPreviousBase",
    ):
        assert forbidden not in src, (
            f"address.rs now carries {forbidden!r}; the audit's "
            "no-new-addresses promise regressed"
        )


# ──────────────────────────────────────────────────────────────────
# 5. Scaffold — JunctionStopState single-shot per execute call
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_junction_stop_state_builds_once_per_execute() -> None:
    """Audit §2 / §4 — ``JunctionStopState::build`` is called
    once per ``execute_with_sampling_mode`` invocation; the
    per-position admit mask is then recomputed lazily by the
    observer's ``current_admit_mask`` call. Composes naturally
    with per-position Markov state. Pin the existing
    single-shot + per-position split."""
    src = (
        _REPO_ROOT
        / "engine_rs"
        / "src"
        / "passes"
        / "generate_np"
        / "execution.rs"
    ).read_text(encoding="utf-8")
    assert "JunctionStopState::build" in src
    assert "builder.current_admit_mask" in src


# ──────────────────────────────────────────────────────────────────
# 6. Scaffold — empty-support sentinel is b'N', unchanged
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_np_base_empty_support_sentinel_unchanged() -> None:
    """Audit §4 — empty-support sentinel remains uppercase
    ``b'N'`` under Markov composition (the admit-mask × Markov
    support intersection produces the same ``Sentinel(b'N')``
    when empty). Pin so a refactor that changes the sentinel
    surfaces here."""
    src = (
        _REPO_ROOT
        / "engine_rs"
        / "src"
        / "passes"
        / "generate_np"
        / "sampling.rs"
    ).read_text(encoding="utf-8")
    assert "EmptySupport::Sentinel(b'N')" in src


# ──────────────────────────────────────────────────────────────────
# 7. Scaffold — UniformBase plan-signature canonical substring
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_uniform_plan_signature_substring_is_canonical() -> None:
    """Audit §5 / §7 — the existing ``UniformBase`` plan signature
    renders as the canonical 4-way A/C/G/T weight-1 list. The
    Markov slice's ``UniformNpGenerator`` wrapper MUST produce
    this exact substring so legacy traces replay unchanged. Pin
    the current behaviour."""
    exp = ga.Experiment.on("human_igh").recombine()
    compiled = exp.compile()
    seed = 42
    tf = compiled.simulator.trace_file_from(
        compiled.simulator.run(seed=seed), seed=seed
    )
    sig = json.loads(tf.to_json())["pass_plan_signature"]
    assert "(65:1.0),(67:1.0),(71:1.0),(84:1.0)" in sig, (
        "UniformBase plan signature substring drifted; the Markov "
        "slice's backwards-compatibility constraint regressed"
    )


# ──────────────────────────────────────────────────────────────────
# 8. Scaffold — empirical_first_base plan signature stable
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_empirical_first_base_signature_stable() -> None:
    """Audit §7 — the existing ``CategoricalBase`` plan signature
    is byte-identical for a given first-base distribution. The
    Markov slice's ``CategoricalNpGenerator`` wrapper MUST
    reproduce this byte-equal signature so existing cartridges
    with ``kind="empirical_first_base"`` replay unchanged."""
    import copy
    cfg = copy.deepcopy(ga.HUMAN_IGH_OGRDB)
    cfg.reference_models = ReferenceEmpiricalModels(
        np_bases={
            "NP1": NpBaseModelSpec(
                kind="empirical_first_base",
                first_base={"A": 1.0, "C": 1.0, "G": 1.0, "T": 1.0},
            )
        }
    )
    compiled = ga.Experiment.on(cfg).recombine().compile()
    seed = 42
    tf = compiled.simulator.trace_file_from(
        compiled.simulator.run(seed=seed), seed=seed
    )
    sig = json.loads(tf.to_json())["pass_plan_signature"]
    # The 4-way categorical's support renders identically to
    # UniformBase's because all four weights are equal.
    assert "(65:1.0),(67:1.0),(71:1.0),(84:1.0)" in sig, (
        "uniform-weighted empirical_first_base no longer produces the "
        "canonical 4-way support substring; the wrapper-type byte-"
        "equality constraint can't hold after the Markov slice"
    )


# ──────────────────────────────────────────────────────────────────
# 9. Scaffold — Markov spec already validates at Python layer
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_np_base_model_spec_already_validates_markov() -> None:
    """Audit §6 — ``NpBaseModelSpec(kind="markov")`` already
    validates correctly at the Python spec layer (per-base
    alphabet, weight finiteness, complete from-base matrix).
    The Markov slice changes nothing in the spec layer; pin
    the existing validation so a regression surfaces."""
    # Well-formed Markov spec passes.
    spec = NpBaseModelSpec(
        kind="markov",
        first_base={"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
        transitions={
            "A": {"A": 0.1, "C": 0.4, "G": 0.4, "T": 0.1},
            "C": {"A": 0.1, "C": 0.4, "G": 0.4, "T": 0.1},
            "G": {"A": 0.1, "C": 0.4, "G": 0.4, "T": 0.1},
            "T": {"A": 0.1, "C": 0.4, "G": 0.4, "T": 0.1},
        },
    )
    assert spec.kind == "markov"

    # Partial matrix still rejected at the spec layer.
    with pytest.raises(ValueError, match="missing rows"):
        NpBaseModelSpec(
            kind="markov",
            first_base={"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
            transitions={"A": {"C": 1.0}},
        )


# ──────────────────────────────────────────────────────────────────
# 10. Scaffold — Markov spec lowering CURRENTLY raises
# ──────────────────────────────────────────────────────────────────


def test_pin_present_markov_lowering_now_runs_end_to_end() -> None:
    """Post-slice — ``kind="markov"`` lowers through the new
    ``MarkovBaseGenerator`` and produces records without
    raising. Flipped from the prior ``pin_present`` that
    expected ``NotImplementedError``."""
    import copy
    cfg = copy.deepcopy(ga.HUMAN_IGH_OGRDB)
    cfg.reference_models = ReferenceEmpiricalModels(
        np_bases={
            "NP1": NpBaseModelSpec(
                kind="markov",
                first_base={"A": 1.0, "C": 1.0, "G": 1.0, "T": 1.0},
                transitions={
                    "A": {"A": 1.0, "C": 1.0, "G": 1.0, "T": 1.0},
                    "C": {"A": 1.0, "C": 1.0, "G": 1.0, "T": 1.0},
                    "G": {"A": 1.0, "C": 1.0, "G": 1.0, "T": 1.0},
                    "T": {"A": 1.0, "C": 1.0, "G": 1.0, "T": 1.0},
                },
            )
        }
    )
    recs = ga.Experiment.on(cfg).recombine().run_records(n=3, seed=0)
    assert len(recs) == 3, "Markov lowering produced no records"


# ──────────────────────────────────────────────────────────────────
# 11. Scaffold — manifest reports markov as DEFERRED today
# ──────────────────────────────────────────────────────────────────


def test_pin_present_manifest_now_lists_markov_as_supported() -> None:
    """Post-slice — ``manifest["models"]["np_base_models"]``
    promotes ``markov`` from ``deferred_kinds`` to
    ``supported_kinds`` per audit §11. Flipped from the prior
    scaffold pin."""
    m = ga.HUMAN_IGH_OGRDB.cartridge_manifest()
    nbm = m["models"]["np_base_models"]
    assert nbm["supported_kinds"] == [
        "uniform",
        "empirical_first_base",
        "markov",
    ]
    assert nbm["deferred_kinds"] == []


# ──────────────────────────────────────────────────────────────────
# 12. Absence — no NpBaseGenerator trait
# ──────────────────────────────────────────────────────────────────


def test_pin_present_np_base_generator_trait_in_rust() -> None:
    """Post-slice — ``trait NpBaseGenerator`` is declared in the
    new ``engine_rs/src/passes/generate_np/np_base_generator.rs``
    module per audit §1 / §12. Flipped from the prior absence
    pin."""
    src = (
        _REPO_ROOT
        / "engine_rs"
        / "src"
        / "passes"
        / "generate_np"
        / "np_base_generator.rs"
    ).read_text(encoding="utf-8")
    assert "pub trait NpBaseGenerator" in src, (
        "NpBaseGenerator trait declaration missing from the module"
    )
    # Both required methods present.
    assert "fn support(&self, position: usize, previous: Option<u8>)" in src
    assert "fn signature(&self)" in src


# ──────────────────────────────────────────────────────────────────
# 13. Absence — no MarkovBaseGenerator concrete type
# ──────────────────────────────────────────────────────────────────


def test_pin_present_markov_base_generator_struct() -> None:
    """Post-slice — the concrete ``MarkovBaseGenerator`` struct
    exists with the ``first_base: [f64; 4]`` +
    ``transitions: [[f64; 4]; 4]`` shape per audit §3."""
    src = (
        _REPO_ROOT
        / "engine_rs"
        / "src"
        / "passes"
        / "generate_np"
        / "np_base_generator.rs"
    ).read_text(encoding="utf-8")
    assert "pub struct MarkovBaseGenerator" in src
    assert "first_base: [f64; 4]" in src
    assert "transitions: [[f64; 4]; 4]" in src


# ──────────────────────────────────────────────────────────────────
# 14. Absence — no Uniform/Categorical NP wrapper types yet
# ──────────────────────────────────────────────────────────────────


def test_pin_present_np_generator_wrapper_types() -> None:
    """Post-slice — both wrapper types
    (``UniformNpGenerator`` and ``CategoricalNpGenerator``) are
    declared in ``np_base_generator.rs`` per audit §3 / §7.
    The wrappers are the byte-identical-signature compatibility
    surface that keeps legacy plan signatures stable."""
    src = (
        _REPO_ROOT
        / "engine_rs"
        / "src"
        / "passes"
        / "generate_np"
        / "np_base_generator.rs"
    ).read_text(encoding="utf-8")
    assert "pub struct UniformNpGenerator" in src
    assert "pub struct CategoricalNpGenerator" in src


# ──────────────────────────────────────────────────────────────────
# 15. Absence — no previous-base param on replay validator
# ──────────────────────────────────────────────────────────────────


def test_pin_present_previous_base_param_on_sample_base() -> None:
    """Post-slice — ``sample_base`` now threads
    ``previous: Option<u8>`` per audit §2. The replay validator
    consumes the per-position support pairs computed from the
    generator with this `previous` so Markov replay is
    byte-deterministic against the same matrix."""
    src = (
        _REPO_ROOT
        / "engine_rs"
        / "src"
        / "passes"
        / "generate_np"
        / "sampling.rs"
    ).read_text(encoding="utf-8")
    assert "fn validate_replayed_np_base" in src
    # `previous: Option<u8>` now appears on sample_base
    # (threaded into the generator's per-position support).
    assert "previous: Option<u8>" in src, (
        "sample_base no longer carries `previous: Option<u8>`; the "
        "Markov slice's wiring change has regressed"
    )
    # Generator consultation also pins.
    assert "base_generator.support(" in src


# ──────────────────────────────────────────────────────────────────
# 16. Absence — no markov_transitions kwarg on push_generate_np
# ──────────────────────────────────────────────────────────────────


def test_pin_present_markov_transitions_kwarg_on_push_generate_np() -> None:
    """Post-slice — ``push_generate_np`` exposes the
    ``markov_transitions`` kwarg per audit §6. Flipped from
    the prior absence pin."""
    src = (
        _REPO_ROOT / "engine_rs" / "src" / "python" / "plan.rs"
    ).read_text(encoding="utf-8")
    # Existing kwarg still present.
    assert "base_pairs: Option<Vec<(u8, f64)>>" in src
    # New Markov kwarg now present at the bridge signature.
    assert "markov_transitions: Option<Vec<Vec<(u8, f64)>>>" in src, (
        "push_generate_np no longer carries `markov_transitions` — the "
        "Markov slice's bridge surface has regressed"
    )


# ──────────────────────────────────────────────────────────────────
# 17. Absence — manifest doesn't yet list markov as SUPPORTED
# ──────────────────────────────────────────────────────────────────


def test_pin_present_markov_in_supported_kinds() -> None:
    """Post-slice — the manifest's ``supported_kinds`` list
    now includes ``markov`` per audit §11. Flipped from the
    prior absence pin."""
    m = ga.HUMAN_IGH_OGRDB.cartridge_manifest()
    nbm = m["models"]["np_base_models"]
    assert "markov" in nbm["supported_kinds"], (
        "manifest no longer lists markov in supported_kinds; the "
        "Markov slice's manifest flip has regressed"
    )
    assert "markov" not in nbm["deferred_kinds"], (
        "manifest still marks markov as deferred; the slice's "
        "manifest flip has regressed"
    )


# ──────────────────────────────────────────────────────────────────
# 18. Scaffold — GenerateNPPass carries Distribution-typed base_dist
# ──────────────────────────────────────────────────────────────────


def test_pin_present_generate_np_pass_owns_base_generator_field() -> None:
    """Post-slice — ``GenerateNPPass`` carries
    ``base_generator: Box<dyn NpBaseGenerator>`` (renamed
    from ``base_dist``). Flipped from the prior scaffold pin
    on the legacy field name."""
    src = (
        _REPO_ROOT / "engine_rs" / "src" / "passes" / "generate_np.rs"
    ).read_text(encoding="utf-8")
    assert "base_generator: Box<dyn NpBaseGenerator>" in src, (
        "GenerateNPPass no longer carries the renamed `base_generator` "
        "field; the Markov slice's field migration has regressed"
    )
    # The struct field must be named `base_generator`. The
    # legacy `new(... base_dist: Box<dyn Distribution<Output=u8>>)`
    # constructor can still accept the old wire type (it wraps
    # the distribution in an internal `NpBaseGenerator`
    # adapter), so `base_dist` may still appear as a
    # parameter name — pin only the struct-field shape.
    assert "pub struct GenerateNPPass" in src


# ──────────────────────────────────────────────────────────────────
# 19. Scaffold — legacy NP_transitions still orphan
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_legacy_np_transitions_still_orphan_after_markov_slice() -> None:
    """Audit §7 — the Markov slice does NOT auto-lift legacy
    ``NP_transitions`` / ``NP_first_bases``. The
    stop-and-report from the prior slice on legacy auto-lift
    remains in force. Pin both fields' continued presence in
    ``_DOCUMENTED_ORPHAN_DATACONFIG_FIELDS`` and the manifest's
    ``legacy_fallback=False`` flag."""
    from GenAIRR.dataconfig.data_config import (
        _DOCUMENTED_ORPHAN_DATACONFIG_FIELDS,
    )

    for required in ("NP_transitions", "NP_first_bases"):
        assert required in _DOCUMENTED_ORPHAN_DATACONFIG_FIELDS, (
            f"{required!r} no longer orphan; auto-lift may have "
            "landed — verify the audit doc + pin flips are in lockstep"
        )
    m = ga.HUMAN_IGH_OGRDB.cartridge_manifest()
    nbm = m["models"]["np_base_models"]
    assert nbm["legacy_fallback"] is False, (
        "manifest's legacy_fallback flag flipped to True; the Markov "
        "slice was NOT supposed to enable auto-lift"
    )


# ──────────────────────────────────────────────────────────────────
# 20. Doc anchor — audit doc exists and references contract file
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_audit_doc_exists_and_references_contract() -> None:
    """The audit doc must continue to exist and reference this
    contract file; the 15-section structure stays intact."""
    if not _AUDIT_DOC.exists():
        import pytest
        pytest.skip("docs/ is contributor-only; audit doc not present in this checkout")
    doc = _AUDIT_DOC.read_text(encoding="utf-8")
    assert "test_np_markov_base_generator_contract.py" in doc, (
        "audit doc no longer references the contract file; lockstep "
        "convention drifted"
    )
    for marker in (
        "## 1. Q1",
        "## 2. Q2",
        "## 7. Q7",
        "## 12. Implementation order",
        "## 15. Summary table",
    ):
        assert marker in doc, f"audit doc missing section marker {marker!r}"
