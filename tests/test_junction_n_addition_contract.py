"""Contract pins for the Junction / N-Addition Modeling audit.

Companion to
[`docs/junction_n_addition_audit.md`](../docs/junction_n_addition_audit.md).
The audit is pre-implementation: it freezes today's surfaces
(``pin_scaffold_*``) and the gaps a future implementation slice
would close (``pin_absence_*``).

**Stop-condition was NOT triggered.** The legacy
``DataConfig.NP_transitions`` and ``DataConfig.NP_first_bases``
fields exist on bundled cartridges with populated content, but
they are explicitly listed in
``_DOCUMENTED_ORPHAN_DATACONFIG_FIELDS`` and are never read by
the simulation pipeline. The current NP base distribution is a
hardcoded ``UniformBase`` constructed at the PyO3 bridge. No
replay-safety hole. The audit may proceed and the implementation
slice can lift the legacy Markov fields into a typed
``NpBaseModelSpec`` on ``ReferenceEmpiricalModels``.

Split:

- ``pin_scaffold_*`` tests freeze the pre-existing surfaces the
  slice builds on: the `GenerateNPPass` struct shape,
  `UniformBase` as the hardcoded default, the typed `np_lengths`
  plane on `ReferenceEmpiricalModels`, the `np_lengths` resolver
  fallback chain, the NP trace addresses, the
  `parameter_signature` discipline (Slice A), the productive-only
  narrowing via `JunctionStopState`, the AIRR projection of NP
  fields, the orphan-list inclusion of `NP_transitions` /
  `NP_first_bases`, the bundled-cartridge content of those
  legacy fields.
- ``pin_absence_*`` tests freeze the gaps the slice closes:
  no typed `NpBaseModelSpec`, no `np_bases` field on
  `ReferenceEmpiricalModels`, no `np_base_keys` in the manifest,
  no `recombine(np_base_model=…)` DSL kwarg, no P-nucleotide
  surface (deferred to a separate future audit).
"""
from __future__ import annotations

import inspect
import json
from pathlib import Path

import pytest

import GenAIRR as ga


_REPO_ROOT = Path(__file__).resolve().parent.parent
_AUDIT_DOC = _REPO_ROOT / "audit-docs" / "junction_n_addition_audit.md"


# ──────────────────────────────────────────────────────────────────
# 1. Scaffold — GenerateNPPass struct shape
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_generate_np_pass_struct_shape() -> None:
    """Audit §1 / §7 — `GenerateNPPass` carries three fields.
    The Markov NP base generator slice renamed the third field
    from ``base_dist: Box<dyn Distribution<Output = u8>>`` to
    ``base_generator: Box<dyn NpBaseGenerator>`` (see
    [`docs/np_markov_base_generator_design.md`](../docs/np_markov_base_generator_design.md)
    §1) so the per-position support is conditioned on the
    previously emitted base. The constructor still accepts a
    legacy ``Box<dyn Distribution<Output = u8>>`` via an
    internal adapter for backwards compatibility."""
    src = (
        _REPO_ROOT / "engine_rs" / "src" / "passes" / "generate_np.rs"
    ).read_text(encoding="utf-8")
    import re

    match = re.search(
        r"pub struct GenerateNPPass\s*\{(.*?)\n\}",
        src,
        re.DOTALL,
    )
    assert match is not None, (
        "GenerateNPPass struct not found in generate_np.rs; audit "
        "§1 reference drifted"
    )
    body = match.group(1)
    for field in ("np_segment:", "length_dist:", "base_generator:"):
        assert field in body, (
            f"GenerateNPPass missing field {field!r}; the Markov "
            "slice's field migration regressed"
        )
    assert "Box<dyn NpBaseGenerator>" in body


# ──────────────────────────────────────────────────────────────────
# 2. Scaffold — UniformBase is the hardcoded base distribution today
# ──────────────────────────────────────────────────────────────────


def test_pin_present_push_generate_np_routes_uniform_or_categorical_base() -> None:
    """Audit §1 / §3 — the implementation slice extended
    `push_generate_np` with an optional `base_pairs` kwarg. When
    `None`, the bridge constructs `Box::new(UniformBase)` (byte-
    identical to pre-slice). When `Some(pairs)`, the bridge
    constructs `Box::new(CategoricalBase::from_pairs(pairs))`
    consuming the cartridge-owned typed model. Pin both branches
    at source so a regression that drops either surfaces here."""
    src = (
        _REPO_ROOT / "engine_rs" / "src" / "python" / "plan.rs"
    ).read_text(encoding="utf-8")
    # The function still exists and carries the new kwarg shape.
    assert "fn push_generate_np" in src
    assert "base_pairs: Option<Vec<(u8, f64)>>" in src, (
        "push_generate_np is missing the base_pairs kwarg; the typed "
        "NP base model slice regressed"
    )
    # Both branches of the match are present.
    assert "Box::new(UniformBase)" in src, (
        "push_generate_np no longer falls through to UniformBase when "
        "base_pairs is None; the uniform fast-path regressed"
    )
    assert "CategoricalBase::from_pairs" in src, (
        "push_generate_np no longer constructs CategoricalBase from "
        "base_pairs; the typed model branch regressed"
    )


# ──────────────────────────────────────────────────────────────────
# 3. Scaffold — typed np_lengths exists on ReferenceEmpiricalModels
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_reference_models_carries_typed_np_lengths() -> None:
    """Audit §3 — `ReferenceEmpiricalModels.np_lengths` exists
    as a typed plane (Dict keyed by `"NP1"` / `"NP2"` mapping to
    `EmpiricalDistributionSpec`). The implementation slice adds a
    parallel typed `np_bases` plane; pin the precedent here."""
    from GenAIRR.reference_models import (
        EmpiricalDistributionSpec,
        ReferenceEmpiricalModels,
    )

    sig = inspect.signature(ReferenceEmpiricalModels)
    params = sig.parameters
    assert "np_lengths" in params, (
        "ReferenceEmpiricalModels is missing np_lengths; the audit's "
        "typed-plane precedent regressed"
    )
    # The default factory yields a dict.
    instance = ReferenceEmpiricalModels()
    assert hasattr(instance, "np_lengths")
    assert isinstance(instance.np_lengths, dict)
    # And EmpiricalDistributionSpec is the typed shape np_bases will
    # likely mirror.
    assert EmpiricalDistributionSpec is not None


# ──────────────────────────────────────────────────────────────────
# 4. Scaffold — orphan list carries NP_transitions and NP_first_bases
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_orphan_list_carries_np_transitions_and_first_bases() -> None:
    """Audit §3 — `_DOCUMENTED_ORPHAN_DATACONFIG_FIELDS` lists
    both `NP_transitions` and `NP_first_bases` as documented
    completeness gaps. The implementation slice REMOVES them from
    this list when the typed NP base model consumes them; pinning
    here is what makes that flip auditable."""
    from GenAIRR.dataconfig.data_config import (
        _DOCUMENTED_ORPHAN_DATACONFIG_FIELDS,
    )

    for required in ("NP_transitions", "NP_first_bases"):
        assert required in _DOCUMENTED_ORPHAN_DATACONFIG_FIELDS, (
            f"{required!r} is no longer in the orphan-field list; "
            "the implementation slice may have lifted it into the "
            "typed surface — verify pin flips for the typed NpBaseModel"
            " spec landed in lockstep"
        )


# ──────────────────────────────────────────────────────────────────
# 5. Scaffold — bundled cartridges carry populated NP_transitions
# ──────────────────────────────────────────────────────────────────


@pytest.mark.parametrize(
    "preset",
    ["HUMAN_IGH_OGRDB", "HUMAN_IGK_OGRDB", "HUMAN_IGL_OGRDB"],
)
def test_pin_scaffold_bundled_cartridge_carries_legacy_markov_fields(
    preset: str,
) -> None:
    """Audit §2 / §9 — bundled OGRDB cartridges have populated
    `NP_transitions` and `NP_first_bases` dicts despite the
    bridge ignoring them. The implementation slice consumes these
    fields as the legacy fallback when building the typed
    `NpBaseModelSpec`; pinning the bundled data here confirms the
    fallback path has something to consume."""
    cfg = getattr(ga, preset)
    # `NP_transitions["NP1"]` is a Markov from-base → to-base map.
    assert cfg.NP_transitions, (
        f"{preset}.NP_transitions is empty; the implementation "
        "slice's legacy-Markov fallback has nothing to lift"
    )
    # HUMAN IGK is VJ-chain → only NP1 expected; HUMAN IGH / IGL are
    # heavier so we don't assert NP2 specifically here.
    assert "NP1" in cfg.NP_transitions, (
        f"{preset}.NP_transitions is missing the NP1 key"
    )
    assert cfg.NP_first_bases, (
        f"{preset}.NP_first_bases is empty"
    )


# ──────────────────────────────────────────────────────────────────
# 6. Scaffold — trace addresses for NP exist and replay byte-identical
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_np_trace_addresses_exist() -> None:
    """Audit §4 / §10 — the four NP trace addresses
    (`np.np1.length`, `np.np1.bases[i]`, `np.np2.length`,
    `np.np2.bases[i]`) exist today. The implementation slice does
    NOT introduce any new addresses; the per-base recording stays
    on the existing address namespace."""
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
            f"address.rs missing {required!r}; audit §4 address "
            "vocabulary drifted"
        )


def test_pin_scaffold_np_replay_is_byte_deterministic() -> None:
    """Audit §4 / §10 — same-seed replay through
    `rerun_from_trace_file` reproduces NP1/NP2 bytes exactly. The
    implementation slice's typed base distribution must preserve
    this property — pinning the baseline behaviour here gives the
    slice a regression check."""
    exp = ga.Experiment.on("human_igh").recombine()
    compiled = exp.compile()
    refdata = exp.refdata
    seed = 4242
    fresh = compiled.simulator.run(seed=seed)
    tf = compiled.simulator.trace_file_from(fresh, seed=seed)
    replayed = compiled.simulator.rerun_from_trace_file(tf)
    from GenAIRR._airr_record import outcome_to_airr_record

    fresh_rec = outcome_to_airr_record(fresh, refdata, sequence_id="fresh")
    replayed_rec = outcome_to_airr_record(
        replayed, refdata, sequence_id="replay"
    )
    # Sequence + NP fields round-trip exactly.
    assert fresh_rec["sequence"] == replayed_rec["sequence"]
    assert fresh_rec["np1"] == replayed_rec["np1"]
    assert fresh_rec["np2"] == replayed_rec["np2"]
    assert fresh_rec["np1_length"] == replayed_rec["np1_length"]
    assert fresh_rec["np2_length"] == replayed_rec["np2_length"]


# ──────────────────────────────────────────────────────────────────
# 7. Scaffold — generate_np.parameter_signature folds base_dist
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_generate_np_parameter_signature_folds_base_generator() -> None:
    """Audit §4 / §10 — `GenerateNPPass::parameter_signature`
    folds BOTH the length distribution AND the base generator
    (Slice A discipline). The Markov slice replaced the
    ``fmt_byte_dist(base_dist)`` fold with the generator's
    own ``signature()`` method so the Markov first-base row +
    transition rows all participate in plan identity (see
    [`docs/np_markov_base_generator_design.md`](../docs/np_markov_base_generator_design.md)
    §5).

    Pinning the dual fold at source so a regression that drops
    either fold surfaces."""
    src = (
        _REPO_ROOT / "engine_rs" / "src" / "passes" / "generate_np.rs"
    ).read_text(encoding="utf-8")
    import re

    match = re.search(
        r"fn parameter_signature\(&self\) -> String \{(.*?)\n    \}",
        src,
        re.DOTALL,
    )
    assert match is not None, (
        "parameter_signature not found on GenerateNPPass; the Slice A "
        "discipline regressed"
    )
    body = match.group(1)
    assert "length_dist" in body
    assert "fmt_int_dist" in body
    # Base side now folded via the generator's signature, not
    # via fmt_byte_dist directly. The wrappers
    # (UniformNpGenerator / CategoricalNpGenerator) delegate to
    # fmt_byte_dist internally for byte-identical legacy
    # signatures; the MarkovBaseGenerator returns its own
    # 5-row flat string.
    assert "base_generator.signature()" in body


def test_pin_scaffold_uniform_base_plan_signature_substring() -> None:
    """Audit §4 / §10 — a fresh emission carries the
    `UniformBase` 4-way support in the plan signature (the bytes
    `(65:1.0),(67:1.0),(71:1.0),(84:1.0)` — ASCII for ACGT). When
    the implementation slice replaces `UniformBase` with a
    non-uniform distribution, this substring will change to
    reflect the new weights — pin the current behaviour so the
    flip is visible."""
    exp = ga.Experiment.on("human_igh").recombine()
    compiled = exp.compile()
    seed = 42
    tf = compiled.simulator.trace_file_from(
        compiled.simulator.run(seed=seed), seed=seed
    )
    sig = json.loads(tf.to_json())["pass_plan_signature"]
    assert "generate_np.np1" in sig
    # The UniformBase 4-way support renders as a Phred+33-style
    # weight list; the byte literals for A/C/G/T are 65/67/71/84.
    assert "(65:1.0),(67:1.0),(71:1.0),(84:1.0)" in sig, (
        "plan signature no longer renders UniformBase's 4-way "
        "support — the implementation slice may have changed the "
        "base distribution"
    )


# ──────────────────────────────────────────────────────────────────
# 8. Scaffold — productive_only narrowing + N sentinel
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_productive_only_preserves_triad_under_np() -> None:
    """Audit §5 — `productive_only()` composes correctly with
    NP generation. The implementation slice's typed base
    distribution must preserve the productive triad invariants;
    pin the baseline."""
    exp = ga.Experiment.on("human_igh").recombine().productive_only()
    refdata = exp.refdata
    result = exp.run_records(n=20, seed=4242)
    report = result.validate_records(refdata)
    assert report, (
        f"productive_only + NP baseline validation failed: "
        f"{report.summary()}"
    )


def test_pin_scaffold_np_base_empty_support_sentinel_is_capital_n() -> None:
    """Audit §5 — the empty-support sentinel for NP base
    sampling is uppercase `b'N'`. Pinning at source so a future
    slice that changes the sentinel surfaces here."""
    src = (
        _REPO_ROOT
        / "engine_rs"
        / "src"
        / "passes"
        / "generate_np"
        / "sampling.rs"
    ).read_text(encoding="utf-8")
    assert "EmptySupport::Sentinel(b'N')" in src, (
        "NP base empty-support sentinel changed from b'N'; "
        "audit §5 fact regressed"
    )


# ──────────────────────────────────────────────────────────────────
# 9. Scaffold — AIRR NP fields present on every record
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_airr_np_fields_present() -> None:
    """Audit §6 — every AIRR record carries `np1`, `np1_length`,
    `np2`, `np2_length` as first-class fields. The implementation
    slice doesn't change projection — these fields stay
    authoritative."""
    result = ga.Experiment.on("human_igh").recombine().run_records(n=1, seed=0)
    rec = result.records[0]
    for required in (
        "np1",
        "np1_length",
        "np2",
        "np2_length",
        "np1_aa",
        "np2_aa",
    ):
        assert required in rec, (
            f"AIRR record missing {required!r}; the NP projection "
            "surface regressed"
        )


# ──────────────────────────────────────────────────────────────────
# 10. Scaffold — chain-type behaviour: VJ runs NP1 only, VDJ both
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_vj_chain_runs_only_np1() -> None:
    """Audit §1 — VJ chains (e.g. human IGK / IGL) run
    `generate_np.np1` only because there's no D segment between
    V and J. The plan signature for a VJ recombine confirms NP2
    is not in the plan."""
    exp = ga.Experiment.on("human_igk").recombine()
    compiled = exp.compile()
    seed = 42
    tf = compiled.simulator.trace_file_from(
        compiled.simulator.run(seed=seed), seed=seed
    )
    sig = json.loads(tf.to_json())["pass_plan_signature"]
    assert "generate_np.np1" in sig
    assert "generate_np.np2" not in sig, (
        "VJ chain plan signature now includes generate_np.np2; the "
        "chain-type behaviour audit's pin regressed"
    )


def test_pin_scaffold_vdj_chain_runs_both_np1_and_np2() -> None:
    """Audit §1 — VDJ chains (human IGH) run both NP segments."""
    exp = ga.Experiment.on("human_igh").recombine()
    compiled = exp.compile()
    seed = 42
    tf = compiled.simulator.trace_file_from(
        compiled.simulator.run(seed=seed), seed=seed
    )
    sig = json.loads(tf.to_json())["pass_plan_signature"]
    assert "generate_np.np1" in sig
    assert "generate_np.np2" in sig


# ──────────────────────────────────────────────────────────────────
# 11. Absence — typed NpBaseModelSpec doesn't exist yet
# ──────────────────────────────────────────────────────────────────


def test_pin_present_np_base_model_spec_type() -> None:
    """Audit §3 / §7 — the typed `NpBaseModelSpec` dataclass landed
    on `reference_models.py` with three validated `kind` values.
    Alternative names stay forbidden; the canonical surface is the
    single dataclass shipping today."""
    from GenAIRR import reference_models as ref_module

    assert hasattr(ref_module, "NpBaseModelSpec"), (
        "NpBaseModelSpec missing from reference_models; the typed NP "
        "base model slice regressed"
    )
    # The dataclass is also re-exported at the package root for
    # ergonomic authoring.
    import GenAIRR as _ga
    assert hasattr(_ga, "NpBaseModelSpec")
    # Alternative names stay forbidden — keep the canonical surface
    # narrow.
    for forbidden in (
        "NPBaseModelSpec",
        "NpBaseModel",
        "MarkovBaseSpec",
        "FirstBaseSpec",
    ):
        assert not hasattr(ref_module, forbidden), (
            f"reference_models now exposes alternate-name {forbidden!r}; "
            "only `NpBaseModelSpec` is the canonical surface"
        )


def test_pin_present_np_bases_field_on_reference_models() -> None:
    """Audit §3 / §7 — `ReferenceEmpiricalModels.np_bases` landed
    as the typed plane for cartridge-owned NP base distributions.
    Alternate field names stay forbidden."""
    from GenAIRR.reference_models import ReferenceEmpiricalModels

    sig = inspect.signature(ReferenceEmpiricalModels)
    assert "np_bases" in sig.parameters, (
        "ReferenceEmpiricalModels is missing the np_bases field; "
        "the typed NP base model slice regressed"
    )
    # Default factory produces an empty dict (no model means
    # UniformBase by default — byte-identical to pre-slice).
    instance = ReferenceEmpiricalModels()
    assert hasattr(instance, "np_bases")
    assert instance.np_bases == {}
    # Alternate field names stay forbidden.
    for forbidden in ("np_base_models", "np_transition_models"):
        assert forbidden not in sig.parameters, (
            f"ReferenceEmpiricalModels now exposes alternate-name "
            f"{forbidden!r}; only `np_bases` is the canonical surface"
        )


# ──────────────────────────────────────────────────────────────────
# 12. Absence — manifest doesn't yet advertise NP base capability
# ──────────────────────────────────────────────────────────────────


def test_pin_present_np_base_models_block_in_manifest() -> None:
    """Audit §11 — the cartridge manifest now carries an
    ``np_base_models`` block under ``models`` advertising the
    typed NP base model capability, the per-region kinds (if
    any), and the legacy-fallback decision. Pre-existing
    ``np_length_keys`` stays present."""
    m = ga.HUMAN_IGH_OGRDB.cartridge_manifest()
    models = m["models"]
    # Pre-existing length plane still advertised.
    assert "np_length_keys" in models
    # New block landed.
    assert "np_base_models" in models, (
        "manifest is missing np_base_models; the implementation slice "
        "regressed"
    )
    nbm = models["np_base_models"]
    assert "legacy_fallback" in nbm
    assert nbm["legacy_fallback"] is False, (
        "manifest's legacy_fallback is True; that means the resolver "
        "is auto-lifting NP_transitions / NP_first_bases — which the "
        "audit explicitly deferred. Verify a deliberate slice landed "
        "before flipping this pin"
    )
    assert nbm["legacy_np_transitions_present"] is True
    assert nbm["legacy_np_first_bases_present"] is True
    assert nbm["supported_kinds"] == [
        "uniform",
        "empirical_first_base",
        "markov",
    ]
    assert nbm["deferred_kinds"] == []
    assert nbm["in_plan_signature"] is True
    assert nbm["in_content_hash"] is False


# ──────────────────────────────────────────────────────────────────
# 13. Absence — recombine() has no np_base_model kwarg yet
# ──────────────────────────────────────────────────────────────────


def test_pin_absence_no_np_base_model_kwarg_on_recombine() -> None:
    """The implementation slice ships **cartridge-owned only**
    typed NP base distributions — the model is authored on
    ``cfg.reference_models.np_bases`` and consumed at compile
    time via the resolver cascade. No per-experiment
    ``Experiment.recombine(np_base_model=…)`` kwarg landed —
    matching the user brief's scope discipline. A future
    ergonomics slice could add the kwarg following the existing
    ``np1_lengths`` / ``np2_lengths`` precedent; pin its
    continued absence so that ergonomics slice flips this in
    lockstep with the kwarg landing."""
    sig = inspect.signature(ga.Experiment.recombine)
    # The existing length-override kwargs continue to surface.
    assert "np1_lengths" in sig.parameters
    # The per-experiment NP-base override kwargs are deliberately
    # not added by this slice.
    for forbidden in (
        "np_base_model",
        "np1_base_model",
        "np2_base_model",
        "np_bases",
        "np_first_bases",
        "np_transitions",
    ):
        assert forbidden not in sig.parameters, (
            f"recombine() now accepts {forbidden!r}; the per-experiment "
            "NP-base ergonomics slice has landed — flip pin and add "
            "positive tests for the new kwarg surface"
        )


# ──────────────────────────────────────────────────────────────────
# 14. Absence — no P-nucleotide / palindromic-addition surface
# ──────────────────────────────────────────────────────────────────


def test_pin_present_p_nucleotide_surface_post_v1_slice() -> None:
    """Post-slice — the P-nucleotide v1 implementation landed
    end-to-end. Surfaces present:

    - Typed cartridge plane
      `ReferenceEmpiricalModels.p_nucleotide_lengths` (Python).
    - Engine `PAdditionPass`, `PEnd` enum, `PRegionAdded`
      event variant, `ChoiceAddress::PLength` variant, four
      `p_{v_3,d_5,d_3,j_5}.length` trace addresses.
    - AIRR `p_v_3_length` / `p_d_5_length` / `p_d_3_length` /
      `p_j_5_length` int fields.
    - Validator `PLengthMismatch` issue kind.
    - Manifest `models.p_nucleotide_models` block.

    DSL note: v1 is **cartridge-only** — no
    `.p_nucleotides(...)` method on `Experiment`; authors
    populate `ReferenceEmpiricalModels.p_nucleotide_lengths`
    explicitly (same style as `np_bases`). The per-base
    P-string surfaces (`p_v_3`, ...) and the aggregate
    `n_p_nucleotides` field remain out-of-scope per
    [docs/p_nucleotide_design.md](../docs/p_nucleotide_design.md)
    §15.

    Flipped from the prior absence pin."""
    # DSL stays cartridge-only — no kwarg on recombine, no
    # standalone .p_nucleotides(...) method.
    sig = inspect.signature(ga.Experiment.recombine)
    assert "p_nucleotides" not in sig.parameters
    assert "p_nucleotide_lengths" not in sig.parameters
    assert not hasattr(ga.Experiment, "p_nucleotides")

    # AIRR record gains the four length fields.
    result = ga.Experiment.on("human_igh").recombine().run_records(n=1, seed=0)
    rec = result.records[0]
    for required in (
        "p_v_3_length",
        "p_d_5_length",
        "p_d_3_length",
        "p_j_5_length",
    ):
        assert required in rec
    # v2 surfaces still rejected — per-base strings,
    # aggregate counters.
    for forbidden in (
        "n_p_nucleotides",
        "p_nucleotide_count",
        "n_palindromic",
        "p_v_3",  # per-base string (without `_length`)
        "p_d_5",
        "p_d_3",
        "p_j_5",
    ):
        assert forbidden not in rec, (
            f"AIRR record carries {forbidden!r} — v1 boundary "
            "(lengths-only, no per-base strings or aggregate counters) "
            "regressed"
        )

    # Rust addresses gain the PEnd enum + PLength variant.
    address_src = (
        _REPO_ROOT / "engine_rs" / "src" / "address.rs"
    ).read_text(encoding="utf-8")
    assert "pub enum PEnd" in address_src
    assert "PLength { end: PEnd }" in address_src
    # Pass-name constants are exported.
    assert "P_ADDITION_V_3" in address_src
    assert "P_ADDITION_J_5" in address_src


# ──────────────────────────────────────────────────────────────────
# 15. Doc anchor — audit doc exists and references contract file
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_audit_doc_exists_and_references_contract() -> None:
    """The audit doc must continue to exist and reference this
    contract file; the 15-section structure stays intact."""
    doc = _AUDIT_DOC.read_text(encoding="utf-8")
    assert "test_junction_n_addition_contract.py" in doc, (
        "audit doc no longer references the contract file; lockstep "
        "convention drifted"
    )
    for marker in (
        "## 1. Q1",
        "## 3. Q3",
        "## 7. Q7",
        "## 12. Implementation order",
        "## 15. Summary table",
    ):
        assert marker in doc, f"audit doc missing section marker {marker!r}"
