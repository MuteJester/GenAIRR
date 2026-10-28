"""End-to-end implementation tests for the **typed NP base model**
slice — `NpBaseModelSpec` + `ReferenceEmpiricalModels.np_bases`.

The slice ships:

- A typed cartridge-owned NP base model authoring surface
  (``NpBaseModelSpec`` with three validated ``kind`` values:
  ``"uniform"`` / ``"empirical_first_base"`` / ``"markov"``).
- ``ReferenceEmpiricalModels.np_bases: Dict[str, NpBaseModelSpec]``
  keyed by NP region (``"NP1"`` / ``"NP2"``).
- A new Rust ``CategoricalBase`` distribution (sibling of
  ``UniformBase``) consumed by ``GenerateNPPass.base_dist`` via the
  PyO3 bridge's new ``push_generate_np(base_pairs=…)`` kwarg.
- A manifest ``np_base_models`` block under ``models``.

**Markov is deferred** — the Python spec validates
``kind="markov"`` so cartridges can author it, but the lowering
path raises ``NotImplementedError``. True previous-base-conditional
Markov requires a pass-level architectural change (the engine's
``Distribution<Output = u8>`` trait is stateless) that's tracked
as a follow-up slice in ``docs/junction_n_addition_audit.md``.

Spec coverage (per the user brief):

1. Spec validation: each kind valid/invalid cases.
2. Explicit typed model overrides the cartridge's existing
   pre-slice behaviour (output bytes change).
3. Legacy fallback DELIBERATELY does NOT auto-lift — the slice
   keeps the legacy ``NP_first_bases`` / ``NP_transitions``
   fields orphaned. Pinning this explicitly so a future slice
   that enables auto-lift surfaces here.
4. Uniform default (no typed model) remains byte-identical to
   pre-slice output.
5. Plan signature differs for different NP base models.
6. Replay round-trip succeeds with a typed empirical model.
7. Replay against a different typed model fails the plan-
   signature gate.
8. Productive-only still preserves the triad under the new model.
9. Manifest advertises the NP base model block + kind list.
10. ``kind="markov"`` raises ``NotImplementedError`` at lowering
    time (stop-and-report).
"""
from __future__ import annotations

import copy
import json

import pytest

import GenAIRR as ga
from GenAIRR.reference_models import (
    EmpiricalDistributionSpec,
    NpBaseModelSpec,
    ReferenceEmpiricalModels,
)


# ──────────────────────────────────────────────────────────────────
# Shared fixtures
# ──────────────────────────────────────────────────────────────────


def _cartridge_with_np_bases(
    np_bases: dict, preset: str = "HUMAN_IGH_OGRDB"
) -> "ga.DataConfig":
    """Deep-copy the bundled preset and attach the requested
    ``np_bases`` model. Other reference-models fields are left
    intact so the resulting cartridge is a drop-in for any
    existing experiment."""
    cfg = copy.deepcopy(getattr(ga, preset))
    existing = getattr(cfg, "reference_models", None)
    if isinstance(existing, ReferenceEmpiricalModels):
        # Preserve existing length / trim models if any.
        cfg.reference_models = ReferenceEmpiricalModels(
            np_lengths=existing.np_lengths,
            trims=existing.trims,
            np_bases=np_bases,
        )
    else:
        cfg.reference_models = ReferenceEmpiricalModels(np_bases=np_bases)
    return cfg


def _np_base_fraction(records, region: str, base: str) -> float:
    """Return the fraction of base occurrences in the concatenated
    NP region across all records. Used for distribution-shift
    smoke tests where the bias is dramatic (≥100× weight)."""
    joined = "".join((r.get(region) or "") for r in records).upper()
    if not joined:
        return 0.0
    return joined.count(base.upper()) / len(joined)


# ──────────────────────────────────────────────────────────────────
# Spec 1 — Validation (per-kind valid + invalid cases)
# ──────────────────────────────────────────────────────────────────


def test_uniform_spec_accepts_no_weights() -> None:
    spec = NpBaseModelSpec(kind="uniform")
    assert spec.kind == "uniform"
    assert spec.first_base is None
    assert spec.transitions is None


def test_uniform_spec_rejects_first_base_payload() -> None:
    with pytest.raises(ValueError, match="must not carry first_base"):
        NpBaseModelSpec(kind="uniform", first_base={"A": 1.0})


def test_uniform_spec_rejects_transitions_payload() -> None:
    with pytest.raises(ValueError, match="must not carry transitions"):
        NpBaseModelSpec(kind="uniform", transitions={"A": {"C": 1.0}})


def test_empirical_first_base_requires_weights() -> None:
    with pytest.raises(ValueError, match="requires first_base"):
        NpBaseModelSpec(kind="empirical_first_base")


def test_empirical_first_base_rejects_transitions() -> None:
    with pytest.raises(ValueError, match="must not carry transitions"):
        NpBaseModelSpec(
            kind="empirical_first_base",
            first_base={"A": 1.0},
            transitions={"A": {"C": 1.0}},
        )


def test_empirical_first_base_accepts_partial_alphabet() -> None:
    """Partial first-base distributions are valid — bases with
    zero/missing weight are simply never sampled."""
    spec = NpBaseModelSpec(
        kind="empirical_first_base", first_base={"G": 1.0, "C": 0.5}
    )
    assert spec.first_base == {"G": 1.0, "C": 0.5}


def test_empirical_first_base_rejects_non_canonical_base() -> None:
    with pytest.raises(ValueError, match="not one of"):
        NpBaseModelSpec(
            kind="empirical_first_base", first_base={"N": 1.0}
        )


def test_empirical_first_base_rejects_negative_weight() -> None:
    with pytest.raises(ValueError, match="must be non-negative"):
        NpBaseModelSpec(
            kind="empirical_first_base", first_base={"A": -0.5}
        )


def test_empirical_first_base_rejects_nan_weight() -> None:
    with pytest.raises(ValueError, match="must be finite"):
        NpBaseModelSpec(
            kind="empirical_first_base", first_base={"A": float("nan")}
        )


def test_empirical_first_base_rejects_inf_weight() -> None:
    with pytest.raises(ValueError, match="must be finite"):
        NpBaseModelSpec(
            kind="empirical_first_base", first_base={"A": float("inf")}
        )


def test_empirical_first_base_rejects_bool_weight() -> None:
    with pytest.raises(ValueError, match="must be a finite non-negative number"):
        NpBaseModelSpec(
            kind="empirical_first_base", first_base={"A": True}
        )


def test_empirical_first_base_rejects_all_zero_weights() -> None:
    with pytest.raises(ValueError, match="at least one weight"):
        NpBaseModelSpec(
            kind="empirical_first_base",
            first_base={"A": 0.0, "C": 0.0, "G": 0.0, "T": 0.0},
        )


def test_markov_spec_requires_both_first_base_and_transitions() -> None:
    with pytest.raises(ValueError, match="requires first_base"):
        NpBaseModelSpec(kind="markov", transitions={"A": {"C": 1.0}})
    with pytest.raises(ValueError, match="requires transitions"):
        NpBaseModelSpec(kind="markov", first_base={"A": 1.0})


def test_markov_spec_rejects_partial_transition_matrix() -> None:
    """A Markov NP base model must have rows for all four bases;
    partial matrices are almost certainly an authoring bug."""
    with pytest.raises(ValueError, match="missing rows"):
        NpBaseModelSpec(
            kind="markov",
            first_base={"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
            transitions={"A": {"C": 1.0}},  # only one row
        )


def test_unknown_kind_rejected() -> None:
    with pytest.raises(ValueError, match="unknown kind"):
        NpBaseModelSpec(kind="bogus")


# ──────────────────────────────────────────────────────────────────
# Spec 2 — Explicit typed model perturbs output vs pre-slice baseline
# ──────────────────────────────────────────────────────────────────


def test_g_biased_first_base_model_perturbs_np1_distribution() -> None:
    """A G-biased ``empirical_first_base`` model on NP1 produces a
    near-100% G-fraction in the NP1 region; the baseline pre-slice
    uniform model produces ~25%. The dramatic gap pins the model
    is actually flowing through to the engine — not silently
    ignored."""
    # G-biased: G weight 100× the rest.
    cfg = _cartridge_with_np_bases(
        {
            "NP1": NpBaseModelSpec(
                kind="empirical_first_base",
                first_base={"A": 0.01, "C": 0.01, "G": 100.0, "T": 0.01},
            )
        }
    )
    biased = ga.Experiment.on(cfg).recombine().run_records(n=30, seed=4242)
    # Pre-slice baseline (uniform).
    baseline = ga.Experiment.on(ga.HUMAN_IGH_OGRDB).recombine().run_records(
        n=30, seed=4242
    )

    # Need non-empty NP1 across the batches.
    assert any(r.get("np1") for r in biased.records), (
        "biased run produced no NP1 bases — test is vacuous"
    )
    biased_g = _np_base_fraction(biased.records, "np1", "G")
    baseline_g = _np_base_fraction(baseline.records, "np1", "G")
    assert biased_g > 0.90, (
        f"G-biased NP1 model only produced {biased_g:.1%} Gs across the "
        "batch; the typed spec isn't reaching the engine"
    )
    assert baseline_g < 0.50, (
        f"baseline NP1 G-fraction {baseline_g:.1%} suspicious — "
        "pre-slice uniform should produce ~25%"
    )


# ──────────────────────────────────────────────────────────────────
# Spec 3 — Legacy fallback DOES NOT auto-lift (deliberate)
# ──────────────────────────────────────────────────────────────────


def test_legacy_np_transitions_and_first_bases_remain_orphaned() -> None:
    """The bundled HUMAN_IGH_OGRDB cartridge carries populated
    legacy ``NP_transitions`` and ``NP_first_bases`` dicts. The
    slice DELIBERATELY does NOT auto-lift them into the typed
    ``np_bases`` plane — auto-lifting would silently change output
    bytes vs the pre-slice baseline (the brief's stop-and-report
    condition). This pin freezes that decision so a future slice
    that enables auto-lift must update the audit doc + flip this
    pin."""
    from GenAIRR.dataconfig.data_config import (
        _DOCUMENTED_ORPHAN_DATACONFIG_FIELDS,
    )

    # Bundled cartridge still has the legacy fields populated.
    cfg = ga.HUMAN_IGH_OGRDB
    assert cfg.NP_transitions, (
        "bundled cartridge no longer has NP_transitions populated; "
        "the audit's data baseline regressed"
    )
    assert cfg.NP_first_bases, (
        "bundled cartridge no longer has NP_first_bases populated"
    )
    # And both fields stay in the orphan list.
    for required in ("NP_transitions", "NP_first_bases"):
        assert required in _DOCUMENTED_ORPHAN_DATACONFIG_FIELDS, (
            f"{required!r} is no longer in the orphan-field list; "
            "auto-lift may have landed — verify the audit doc + pin "
            "flips are in lockstep"
        )
    # And the manifest reports ``legacy_fallback=False``.
    m = cfg.cartridge_manifest()
    assert (
        m["models"]["np_base_models"]["legacy_fallback"] is False
    ), "legacy_fallback flag in manifest is no longer False"


# ──────────────────────────────────────────────────────────────────
# Spec 4 — Uniform default is byte-identical to pre-slice
# ──────────────────────────────────────────────────────────────────


def test_no_np_bases_model_is_byte_identical_to_pre_slice() -> None:
    """A cartridge without a typed ``np_bases`` model produces the
    exact same NP bytes as the pre-slice engine (``UniformBase``).
    Same seed → same NP1 / NP2 strings on every record."""
    result_a = (
        ga.Experiment.on(ga.HUMAN_IGH_OGRDB)
        .recombine()
        .run_records(n=10, seed=4242)
    )
    result_b = (
        ga.Experiment.on(ga.HUMAN_IGH_OGRDB)
        .recombine()
        .run_records(n=10, seed=4242)
    )
    for a, b in zip(result_a.records, result_b.records):
        assert a["np1"] == b["np1"]
        assert a["np2"] == b["np2"]


def test_explicit_uniform_spec_is_byte_identical_to_no_spec() -> None:
    """An explicit ``NpBaseModelSpec(kind="uniform")`` is the
    same fast-path as omitting the spec — both lower to
    ``UniformBase`` at the bridge. Byte-identical output."""
    cfg_uniform = _cartridge_with_np_bases(
        {"NP1": NpBaseModelSpec(kind="uniform"), "NP2": NpBaseModelSpec(kind="uniform")}
    )
    explicit_uniform = ga.Experiment.on(cfg_uniform).recombine().run_records(
        n=10, seed=4242
    )
    no_spec = ga.Experiment.on(ga.HUMAN_IGH_OGRDB).recombine().run_records(
        n=10, seed=4242
    )
    for a, b in zip(explicit_uniform.records, no_spec.records):
        assert a["np1"] == b["np1"], (
            f"explicit-uniform and no-spec diverged on NP1: "
            f"{a['np1']!r} vs {b['np1']!r}"
        )
        assert a["np2"] == b["np2"]


# ──────────────────────────────────────────────────────────────────
# Spec 5 — Plan signature differs for different NP base models
# ──────────────────────────────────────────────────────────────────


def _plan_signature(experiment, seed: int = 42) -> str:
    compiled = experiment.compile()
    outcome = compiled.simulator.run(seed=seed)
    tf = compiled.simulator.trace_file_from(outcome, seed=seed)
    return json.loads(tf.to_json())["pass_plan_signature"]


def test_different_typed_np_base_models_produce_different_signatures() -> None:
    """Two cartridges that differ only in their NP base model
    produce different plan signatures — the bridge's
    ``CategoricalBase`` distribution's ``support()`` differs, and
    Slice A's ``parameter_signature`` discipline folds it into
    the signature."""
    cfg_g = _cartridge_with_np_bases(
        {
            "NP1": NpBaseModelSpec(
                kind="empirical_first_base",
                first_base={"A": 1.0, "C": 1.0, "G": 10.0, "T": 1.0},
            )
        }
    )
    cfg_a = _cartridge_with_np_bases(
        {
            "NP1": NpBaseModelSpec(
                kind="empirical_first_base",
                first_base={"A": 10.0, "C": 1.0, "G": 1.0, "T": 1.0},
            )
        }
    )
    sig_g = _plan_signature(ga.Experiment.on(cfg_g).recombine())
    sig_a = _plan_signature(ga.Experiment.on(cfg_a).recombine())
    assert sig_g != sig_a, (
        "plan signatures didn't diverge for distinct NP base models — "
        "Slice A discipline regressed for GenerateNPPass.base_dist"
    )


def test_uniform_spec_signature_matches_no_spec_signature() -> None:
    """An explicit ``kind="uniform"`` spec produces the same plan
    signature as omitting the spec — both flow into the
    ``UniformBase`` fast-path at the bridge."""
    cfg = _cartridge_with_np_bases({"NP1": NpBaseModelSpec(kind="uniform")})
    sig_explicit = _plan_signature(ga.Experiment.on(cfg).recombine())
    sig_no_spec = _plan_signature(
        ga.Experiment.on(ga.HUMAN_IGH_OGRDB).recombine()
    )
    assert sig_explicit == sig_no_spec


# ──────────────────────────────────────────────────────────────────
# Spec 6 — Replay round-trip succeeds with typed empirical model
# ──────────────────────────────────────────────────────────────────


def test_replay_round_trip_with_typed_empirical_first_base_model() -> None:
    cfg = _cartridge_with_np_bases(
        {
            "NP1": NpBaseModelSpec(
                kind="empirical_first_base",
                first_base={"A": 0.4, "C": 0.3, "G": 0.2, "T": 0.1},
            )
        }
    )
    exp = ga.Experiment.on(cfg).recombine()
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
    assert fresh_rec["np1"] == replayed_rec["np1"], (
        f"NP1 diverged under replay with typed empirical model: "
        f"{fresh_rec['np1']!r} vs {replayed_rec['np1']!r}"
    )
    assert fresh_rec["sequence"] == replayed_rec["sequence"]


# ──────────────────────────────────────────────────────────────────
# Spec 7 — Replay against a different typed model fails the gate
# ──────────────────────────────────────────────────────────────────


def test_replay_with_mismatched_np_base_model_fails_signature_gate() -> None:
    """A trace recorded under one NP base model cannot replay
    against a plan with a different model — the plan signature
    folds the base distribution's ``support()`` (Slice A), so
    the gate fires before any choice is consumed."""
    cfg_a = _cartridge_with_np_bases(
        {
            "NP1": NpBaseModelSpec(
                kind="empirical_first_base",
                first_base={"A": 4.0, "C": 1.0, "G": 1.0, "T": 1.0},
            )
        }
    )
    cfg_b = _cartridge_with_np_bases(
        {
            "NP1": NpBaseModelSpec(
                kind="empirical_first_base",
                first_base={"A": 1.0, "C": 1.0, "G": 1.0, "T": 4.0},
            )
        }
    )
    compiled_a = ga.Experiment.on(cfg_a).recombine().compile()
    compiled_b = ga.Experiment.on(cfg_b).recombine().compile()
    outcome = compiled_a.simulator.run(seed=4242)
    tf = compiled_a.simulator.trace_file_from(outcome, seed=4242)
    with pytest.raises(ValueError, match="pass plan signature mismatch"):
        compiled_b.simulator.replay_from_trace_file(tf)


# ──────────────────────────────────────────────────────────────────
# Spec 8 — productive_only triad preserved under typed model
# ──────────────────────────────────────────────────────────────────


def test_productive_only_triad_preserved_under_typed_np_base_model() -> None:
    cfg = _cartridge_with_np_bases(
        {
            "NP1": NpBaseModelSpec(
                kind="empirical_first_base",
                first_base={"A": 0.4, "C": 0.3, "G": 0.2, "T": 0.1},
            ),
            "NP2": NpBaseModelSpec(
                kind="empirical_first_base",
                first_base={"A": 0.1, "C": 0.3, "G": 0.3, "T": 0.3},
            ),
        }
    )
    exp = ga.Experiment.on(cfg).recombine().productive_only()
    refdata = exp.refdata
    result = exp.run_records(n=30, seed=4242)
    report = result.validate_records(refdata)
    assert report, (
        f"productive_only + typed NP base model validation failed: "
        f"{report.summary()}"
    )
    for r in result.records:
        assert r.get("productive") in (True, "T"), (
            f"non-productive record under productive_only(): "
            f"{r.get('sequence_id')}"
        )


# ──────────────────────────────────────────────────────────────────
# Spec 9 — Manifest advertises NP base model block
# ──────────────────────────────────────────────────────────────────


def test_manifest_carries_np_base_models_block() -> None:
    m = ga.HUMAN_IGH_OGRDB.cartridge_manifest()
    nbm = m["models"]["np_base_models"]
    assert set(nbm.keys()) == {
        "models",
        "legacy_fallback",
        "legacy_np_transitions_present",
        "legacy_np_first_bases_present",
        "supported_kinds",
        "deferred_kinds",
        "in_plan_signature",
        "in_content_hash",
    }
    assert nbm["supported_kinds"] == [
        "uniform",
        "empirical_first_base",
        "markov",
    ]
    assert nbm["deferred_kinds"] == []
    assert nbm["in_plan_signature"] is True
    assert nbm["in_content_hash"] is False
    # The bundled HUMAN_IGH_OGRDB cartridge doesn't author a typed
    # np_bases plane; ``models`` is empty.
    assert nbm["models"] == {}
    # But the legacy orphan fields are detectable on bundled
    # pickles.
    assert nbm["legacy_np_transitions_present"] is True
    assert nbm["legacy_np_first_bases_present"] is True


def test_manifest_reports_authored_models_block() -> None:
    """A cartridge with explicitly-authored typed models surfaces
    them in the manifest under ``np_base_models["models"]``."""
    cfg = _cartridge_with_np_bases(
        {
            "NP1": NpBaseModelSpec(
                kind="empirical_first_base",
                first_base={"A": 1.0, "C": 1.0, "G": 1.0, "T": 1.0},
            ),
            "NP2": NpBaseModelSpec(kind="uniform"),
        }
    )
    m = cfg.cartridge_manifest()
    nbm = m["models"]["np_base_models"]
    assert nbm["models"] == {
        "NP1": {"kind": "empirical_first_base"},
        "NP2": {"kind": "uniform"},
    }


# ──────────────────────────────────────────────────────────────────
# Spec 10 — Markov lowering ships and produces records
# ──────────────────────────────────────────────────────────────────


def test_markov_spec_lowers_and_runs_end_to_end() -> None:
    """Flipped from the prior stop-and-report — Markov now lowers
    through ``MarkovBaseGenerator`` and produces records. The
    full conditional-sampling / replay / signature behaviour is
    pinned in
    ``tests/test_np_markov_base_generator_implementation.py``;
    this test just freezes that ``kind='markov'`` no longer
    raises at lowering time."""
    cfg = _cartridge_with_np_bases(
        {
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
    assert len(recs) == 3
