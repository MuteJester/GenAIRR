"""Contract pins for the Trim Distribution Estimation audit.

Companion to
[`docs/trim_distribution_estimation_design.md`](../docs/trim_distribution_estimation_design.md).

Pin set:

- ``pin_scaffold_*`` — live surfaces the new estimator
  reuses verbatim:
  - AirrRecord (Rust) carrying all six AIRR trim fields.
  - The Rust builder hard-zeroing `v_trim_5` / `j_trim_3`
    and sourcing the other four from the trace.
  - `ReferenceEmpiricalModels.trims` typed plane +
    `TRIM_KEYS` / `TRIM_KEYS_VJ` constants + the
    chain-aware validator.
  - `EmpiricalDistributionSpec` validator (estimator output
    container — reused verbatim).
  - `_dataconfig_extract.extract_recombine_defaults` /
    `_trim_from_models` resolver.
  - `Experiment.trim(v_3=..., …)` per-experiment override.
  - Engine `TrimPass` + trace addresses + plan-signature
    fold.
  - `ConfigInfo.has_d` chain-type classifier.
  - Builder stage entry shape + idempotency pattern.
  - `csv.DictReader` AIRR-TSV ingestion.
- ``pin_scaffold_*`` — end-loss vs recombination-trim
  separation across every layer (engine pass module, AIRR
  fields, trace addresses, DSL surface).
- ``pin_present_*`` — stop-and-report verification:
  bundled cartridges populate the four trim fields at
  significant rates; the two hard-zero fields stay zero.
- ``pin_present_*`` — documented surface state: trim
  distributions fold into the plan signature (NO inherited
  soft gap, unlike the allele-usage slice).
- ``pin_absence_*`` — the surfaces the implementation slice
  closes (`estimate_trim_distributions` method,
  `min_count` / `pseudocount` kwarg surface, structured
  `models.trims` manifest block).

**Pre-flight verdict (audit §7): clean-yes.** GenAIRR
populates exactly the four AIRR trim fields the estimator
needs at significant rates; the typed plane + bridge
resolver + engine surface + plan-signature folding are
all live; end-loss is rigorously distinct.
"""
from __future__ import annotations

import copy
import inspect
import json
from pathlib import Path

import pytest

import GenAIRR as ga
from GenAIRR.reference_models import (
    EmpiricalDistributionSpec,
    ReferenceEmpiricalModels,
)


_REPO_ROOT = Path(__file__).resolve().parent.parent
_AUDIT_DOC = _REPO_ROOT / "docs" / "trim_distribution_estimation_design.md"


# ──────────────────────────────────────────────────────────────────
# A. pin_scaffold_* — AIRR trim field provenance (Rust + Python)
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_airr_record_carries_six_trim_fields() -> None:
    """The Rust `AirrRecord` struct exposes all six AIRR trim
    fields (`v_trim_5/3`, `d_trim_5/3`, `j_trim_5/3`). Pinned
    at source so a future schema change drops the unused
    pair only after explicit review."""
    src = (
        _REPO_ROOT / "engine_rs" / "src" / "airr_record" / "record.rs"
    ).read_text(encoding="utf-8")
    for field in ("v_trim_5", "v_trim_3", "d_trim_5", "d_trim_3",
                  "j_trim_5", "j_trim_3"):
        assert f"pub {field}: i64" in src, (
            f"AirrRecord.{field} missing — recombination-trim AIRR "
            f"surface regressed"
        )


def test_pin_scaffold_engine_populates_recombination_trim_fields_only() -> None:
    """The Rust AIRR builder sources `v_trim_3` / `d_trim_5` /
    `d_trim_3` / `j_trim_5` from the trace (`trim.v_3` / etc.)
    and hard-zeros `v_trim_5` and `j_trim_3`. The estimator
    audit's hard-zero finding is pinned at source."""
    src = (
        _REPO_ROOT / "engine_rs" / "src" / "airr_record" / "builder.rs"
    ).read_text(encoding="utf-8")
    assert "rec.v_trim_5 = 0;" in src, (
        "v_trim_5 no longer hard-zeroed — audit §1.2 finding regressed"
    )
    assert "rec.j_trim_3 = 0;" in src, (
        "j_trim_3 no longer hard-zeroed — audit §1.2 finding regressed"
    )
    # The four populated fields read from trace via TrimEnd.
    for marker in (
        'segment: VdjSegment::V',
        'segment: VdjSegment::D',
        'segment: VdjSegment::J',
        'end: TrimEnd::Three',
        'end: TrimEnd::Five',
    ):
        assert marker in src, (
            f"trace-sourced trim field plumbing missing {marker!r}"
        )


def test_pin_scaffold_python_projection_hard_zeros_v_trim_5_with_comment() -> None:
    """The legacy Python AIRR projection at
    `_airr_record.py:585` carries the explanatory comment
    "current Experiment.recombine doesn't trim V_5" — the
    documented engine boundary."""
    src = (
        _REPO_ROOT / "src" / "GenAIRR" / "_airr_record.py"
    ).read_text(encoding="utf-8")
    assert "v_trim_5 = 0" in src
    assert "Experiment.recombine doesn't trim V_5" in src, (
        "documentation comment for the v_trim_5 hard-zero is missing — "
        "the audit's boundary explanation is no longer self-documenting"
    )


# ──────────────────────────────────────────────────────────────────
# B. pin_scaffold_* — typed plane + validator + spec
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_reference_empirical_models_trims_plane_exists() -> None:
    """`ReferenceEmpiricalModels` carries the `trims` plane —
    the existing typed surface the estimator writes into."""
    sig = inspect.signature(ReferenceEmpiricalModels)
    assert "trims" in sig.parameters, (
        "ReferenceEmpiricalModels no longer carries trims plane — "
        "typed-plane discipline regressed"
    )


def test_pin_scaffold_trim_keys_constant_holds() -> None:
    """`TRIM_KEYS = ("V_3", "D_5", "D_3", "J_5")` is the
    canonical trim-key vocabulary. Pinned so a future drift
    that adds V_5 / J_3 keys without the matching engine
    pass surfaces here."""
    from GenAIRR.reference_models import TRIM_KEYS, TRIM_KEYS_VJ

    assert TRIM_KEYS == ("V_3", "D_5", "D_3", "J_5"), (
        f"TRIM_KEYS drifted: {TRIM_KEYS!r}; audit §2.1 documented the "
        f"four-key vocabulary as the canonical recombination-trim surface"
    )
    assert TRIM_KEYS_VJ == ("V_3", "J_5"), (
        f"TRIM_KEYS_VJ drifted: {TRIM_KEYS_VJ!r}; audit §2.1 documented "
        f"the VJ subset"
    )


def test_pin_scaffold_trims_validator_rejects_unknown_keys() -> None:
    """The `trims` plane validator rejects keys outside
    `TRIM_KEYS`. Pinned so a future loosening of the
    key-set check surfaces here."""
    spec = ReferenceEmpiricalModels(
        trims={"V_5": EmpiricalDistributionSpec([(0, 1.0)])}  # V_5 NOT in TRIM_KEYS
    )
    with pytest.raises(ValueError, match=r"trims key .* is not recognised"):
        spec.validate()


def test_pin_scaffold_trims_validator_rejects_d_keys_on_vj() -> None:
    """The `trims` plane validator rejects D-trim keys on VJ
    chains. The estimator inherits this boundary verbatim
    when it constructs an `EmpiricalDistributionSpec` for a
    D key on a VJ cartridge."""
    spec = ReferenceEmpiricalModels(
        trims={"D_5": EmpiricalDistributionSpec([(0, 1.0)])}
    )
    with pytest.raises(ValueError, match=r"VJ chain"):
        spec.validate(chain_type="vj")


def test_pin_scaffold_empirical_distribution_spec_rejects_negative_values() -> None:
    """`EmpiricalDistributionSpec` rejects negative integer
    values. The trim estimator passes raw per-row values
    through the spec constructor; the existing validator
    catches malformed cartridges before they corrupt the
    pipeline."""
    with pytest.raises(ValueError):
        EmpiricalDistributionSpec([(-1, 1.0)]).validate(name="trims[V_3]")


def test_pin_scaffold_empirical_distribution_spec_rejects_non_positive_weights() -> None:
    """`EmpiricalDistributionSpec` rejects zero / negative
    weights. The estimator's normalisation step must
    produce strictly positive per-value weights for the
    spec to survive validation."""
    with pytest.raises(ValueError):
        EmpiricalDistributionSpec([(0, 0.0)]).validate(name="trims[V_3]")
    with pytest.raises(ValueError):
        EmpiricalDistributionSpec([(0, -1.0)]).validate(name="trims[V_3]")


# ──────────────────────────────────────────────────────────────────
# C. pin_scaffold_* — bridge resolver + lowering precedence
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_extract_recombine_defaults_consumes_trims_plane() -> None:
    """`extract_recombine_defaults` returns the four trim
    keys with the typed-plane resolver taking precedence
    over the legacy nested-dict path. The estimator's
    output flows into this resolver verbatim."""
    from GenAIRR._dataconfig_extract import extract_recombine_defaults

    defaults = extract_recombine_defaults(ga.HUMAN_IGH_OGRDB)
    for k in ("trim_v_3", "trim_d_5", "trim_d_3", "trim_j_5"):
        assert k in defaults, (
            f"extract_recombine_defaults no longer returns {k!r} — "
            f"the bridge plumbing the estimator relies on regressed"
        )


def test_pin_scaffold_trim_from_models_resolver_is_typed_plane_adapter() -> None:
    """`_dataconfig_extract._trim_from_models` reads the
    typed plane and returns a flat `[(int, float), ...]`
    list. Pinned at source so a future rename / signature
    change is flagged."""
    src = (
        _REPO_ROOT / "src" / "GenAIRR" / "_dataconfig_extract.py"
    ).read_text(encoding="utf-8")
    assert "def _trim_from_models(" in src, (
        "_trim_from_models helper missing — typed-plane lowering "
        "for trims regressed"
    )
    assert "models.trims.get(key)" in src


def test_pin_scaffold_typed_plane_takes_precedence_over_legacy_trim_dicts() -> None:
    """The bridge's precedence is typed plane > legacy
    nested-dict `cfg.trim_dicts`. Verified via direct
    comparison: a cartridge with both an explicit typed
    spec AND legacy `trim_dicts` produces the spec's
    pairs from `extract_recombine_defaults`."""
    from GenAIRR._dataconfig_extract import extract_recombine_defaults

    cfg = copy.deepcopy(ga.HUMAN_IGH_OGRDB)
    # Sanity: legacy trim_dicts is populated on the bundled cartridge.
    assert cfg.trim_dicts, "bundled cartridge has no legacy trim_dicts to test against"
    # Attach a typed plane with a known-single-value V_3 distribution.
    cfg.reference_models = ReferenceEmpiricalModels(
        trims={"V_3": EmpiricalDistributionSpec([(0, 1.0)])}
    )
    defaults = extract_recombine_defaults(cfg)
    typed = defaults["trim_v_3"]
    assert typed is not None
    # Typed-plane wins: single value at 0 with weight 1.0.
    assert typed == [(0, 1.0)], (
        f"trim_v_3 didn't honour the typed plane override: {typed!r}; "
        f"the precedence the estimator relies on regressed"
    )


def test_pin_scaffold_experiment_trim_dsl_accepts_per_segment_overrides() -> None:
    """`Experiment.trim(v_3=…, d_5=…, d_3=…, j_5=…,
    enabled=…)` is the per-experiment override surface
    chained AFTER `recombine()`. The estimator's typed-plane
    output is at priority 2; the explicit `.trim()`
    surface is at priority 1 — verified in the design doc
    §5.4."""
    sig = inspect.signature(ga.Experiment.trim)
    for kw in ("v_3", "d_5", "d_3", "j_5", "enabled"):
        assert kw in sig.parameters, (
            f"Experiment.trim no longer accepts {kw!r} kwarg — the "
            f"DSL override surface for trim distributions regressed"
        )


def test_pin_scaffold_recombine_does_not_take_per_distribution_trim_kwargs() -> None:
    """`Experiment.recombine` does NOT take per-distribution
    trim kwargs — overrides flow through `.trim()`. Pinned
    so a future shim that adds trim kwargs to `recombine`
    (and creates an ambiguity with `.trim()`) surfaces here
    for explicit review."""
    sig = inspect.signature(ga.Experiment.recombine)
    for forbidden in ("v_3_lengths", "d_5_lengths", "d_3_lengths",
                      "j_5_lengths", "v_3_trim", "d_5_trim"):
        assert forbidden not in sig.parameters, (
            f"Experiment.recombine now accepts {forbidden!r} — "
            f"trim DSL surface has drifted from .trim() to recombine()"
        )


# ──────────────────────────────────────────────────────────────────
# D. pin_scaffold_* — engine surface (TrimPass + paramsig)
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_engine_trim_pass_module_exists() -> None:
    """`engine_rs/src/passes/trim.rs` is the canonical
    recombination-trim engine pass. Pinned at source so a
    rename / module move surfaces here."""
    trim_rs = _REPO_ROOT / "engine_rs" / "src" / "passes" / "trim.rs"
    assert trim_rs.exists(), (
        f"{trim_rs} missing — the engine's recombination-trim pass "
        f"was removed or moved"
    )
    src = trim_rs.read_text(encoding="utf-8")
    assert "TrimPass" in src
    # The TrimPass plan-signature uses fmt_int_dist (audit §7.3 #1).
    assert "fmt_int_dist" in src, (
        "TrimPass no longer folds via fmt_int_dist — trim distributions "
        "may have a soft-gap inheritance now, verify"
    )


def test_pin_scaffold_engine_trace_addresses_use_trim_segment_end_convention() -> None:
    """Engine trace addresses for the four trim passes are
    `trim.v_3` / `trim.d_5` / `trim.d_3` / `trim.j_5`. The
    AIRR builder reads exactly these addresses."""
    src = (
        _REPO_ROOT / "engine_rs" / "src" / "address.rs"
    ).read_text(encoding="utf-8")
    # ChoiceAddress::Trim variant exists.
    assert "Trim {" in src or "Trim{" in src or "Trim {\n" in src, (
        "ChoiceAddress::Trim variant missing — trim trace plumbing "
        "regressed"
    )


def test_pin_scaffold_trim_pass_count_matches_segment_end_pairs() -> None:
    """The recombination pipeline contains exactly four
    `TrimPass` instances at compile time: `(V,Three)`,
    `(D,Five)`, `(D,Three)`, `(J,Five)`. Pinned at source
    so a fifth `TrimPass` (e.g. `(V,Five)`) without the
    matching plane key surfaces here."""
    src = (
        _REPO_ROOT / "src" / "GenAIRR" / "_compile.py"
    ).read_text(encoding="utf-8")
    # The lowering site dispatches via .push_trim("V","3"), etc.
    assert 'push_trim("V", "3"' in src
    assert 'push_trim("D", "5"' in src
    assert 'push_trim("D", "3"' in src
    assert 'push_trim("J", "5"' in src
    # No V_5 / J_3 push.
    assert 'push_trim("V", "5"' not in src
    assert 'push_trim("J", "3"' not in src


# ──────────────────────────────────────────────────────────────────
# E. pin_scaffold_* — end-loss vs recombination-trim separation
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_end_loss_pass_module_distinct_from_trim_pass() -> None:
    """`EndLossPass` lives in `engine_rs/src/passes/corrupt/end_loss.rs`
    — a separate corruption-stage module. Pinned at source so a
    future drift that collapses end-loss into the trim path
    surfaces here."""
    end_loss_rs = (
        _REPO_ROOT / "engine_rs" / "src" / "passes" / "corrupt" / "end_loss.rs"
    )
    assert end_loss_rs.exists(), (
        f"{end_loss_rs} missing — end-loss pass was removed or moved"
    )
    src = end_loss_rs.read_text(encoding="utf-8")
    # Trace addresses live under corrupt.end_loss.5 / .3 — distinct from trim.v_3 etc.
    assert "corrupt.end_loss.5" in src
    assert "corrupt.end_loss.3" in src


def test_pin_scaffold_end_loss_airr_fields_distinct_from_trim_fields() -> None:
    """The Rust AirrRecord exposes `end_loss_5_length` and
    `end_loss_3_length` as SEPARATE fields from `v_trim_5` /
    `j_trim_3`. The estimator MUST NOT consume `end_loss_*`
    fields."""
    src = (
        _REPO_ROOT / "engine_rs" / "src" / "airr_record" / "record.rs"
    ).read_text(encoding="utf-8")
    assert "pub end_loss_5_length: i64" in src
    assert "pub end_loss_3_length: i64" in src
    # Both belong on the record alongside (not in place of) the trim fields.
    for trim_field in ("v_trim_5", "j_trim_3"):
        assert f"pub {trim_field}: i64" in src


def test_pin_scaffold_end_loss_record_docstring_distinguishes_from_v_trim_5() -> None:
    """The `end_loss_5_length` docstring at
    `record.rs:204-207` explicitly distinguishes itself from
    `v_trim_5`, naming the latter as the
    "recombination-stage exonuclease trim". Pinned at source
    so a future edit that removes the self-documenting
    distinction surfaces here."""
    src = (
        _REPO_ROOT / "engine_rs" / "src" / "airr_record" / "record.rs"
    ).read_text(encoding="utf-8")
    assert "recombination-stage exonuclease trim" in src, (
        "end_loss_5_length docstring no longer documents the "
        "recombination-trim vs end-loss boundary — verify the "
        "audit's separation discipline is still explicit"
    )


def test_pin_scaffold_end_loss_dsl_wires_to_end_loss_pass_not_trim() -> None:
    """`Experiment.end_loss_5prime` / its
    `primer_trim_5prime` alias both wire to `EndLossPass`,
    NOT `TrimPass`. The estimator MUST NOT consume DSL
    state from the end-loss surface."""
    src = (
        _REPO_ROOT / "src" / "GenAIRR" / "experiment.py"
    ).read_text(encoding="utf-8")
    # Both DSL methods exist.
    assert "def end_loss_5prime(" in src
    assert "def primer_trim_5prime(" in src
    # Alias relationship explicit in source.
    assert "return self.end_loss_5prime(" in src, (
        "primer_trim_5prime no longer documented as alias for "
        "end_loss_5prime — verify the audit's documented mapping"
    )


# ──────────────────────────────────────────────────────────────────
# F. pin_scaffold_* — chain-type classifier + builder shape
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_config_info_has_d_drives_chain_classification() -> None:
    """`ConfigInfo.has_d` is the authoritative classifier
    the trim estimator dispatches on (VJ vs VDJ). Same
    surface the allele-usage estimator uses."""
    from GenAIRR.dataconfig.enums import ChainType
    from GenAIRR.dataconfig.config_info import ConfigInfo

    assert ChainType.BCR_HEAVY.has_d is True
    assert ChainType.BCR_LIGHT_KAPPA.has_d is False
    assert ChainType.TCR_ALPHA.has_d is False
    assert ChainType.TCR_BETA.has_d is True
    sig = inspect.signature(ConfigInfo)
    assert "has_d" in sig.parameters


def test_pin_scaffold_builder_stage_entries_have_inputs_inferred_warnings() -> None:
    """The canonical `{stage, inputs, inferred, warnings}`
    shape is pinned by the release smoke test (`from_fasta`
    + `build` stages). The `estimate_trim_distributions`
    stage MUST follow the same shape."""
    builder = ga.ReferenceCartridgeBuilder.from_fasta(
        v_fasta=">v1*01\nGAGGTG\n",
        j_fasta=">j1*01\nTGGGGC\n",
        chain_type="BCR_LIGHT_KAPPA",
    )
    for entry in builder.report().stages:
        assert set(entry.keys()) >= {"stage", "inputs", "inferred", "warnings"}, (
            f"stage entry missing canonical keys: {entry}"
        )


def test_pin_scaffold_idempotency_pattern_via_replaced_flag_in_v_subregions() -> None:
    """The `replaced=True` flag pattern for idempotency is
    established by `infer_v_subregions` and `estimate_allele_usage`.
    The trim estimator mirrors this discipline."""
    builder = ga.ReferenceCartridgeBuilder.from_fasta(
        v_fasta=">v1*01\nGAG.GTG\n",
        j_fasta=">j1*01\nTGGGGC\n",
        chain_type="BCR_LIGHT_KAPPA",
    )
    builder.infer_v_subregions()
    builder.infer_v_subregions()
    stages = [s for s in builder.report().stages if s["stage"] == "infer_v_subregions"]
    assert len(stages) == 2
    assert stages[0]["inputs"].get("replaced") is False
    assert stages[1]["inputs"].get("replaced") is True


def test_pin_scaffold_csv_dict_reader_stdlib_handles_airr_tsv() -> None:
    """Python stdlib `csv.DictReader` with `delimiter='\\t'`
    is sufficient for AIRR TSV parsing — pinned via the
    same fixture shape the estimator slice will use."""
    import csv
    import io

    fixture = (
        "v_call\tv_trim_3\td_call\td_trim_5\td_trim_3\tj_call\tj_trim_5\n"
        "IGHV1*01\t2\tIGHD1*01\t0\t3\tIGHJ1*01\t1\n"
        "IGHV2*01\t0\tIGHD2*01\t1\t1\tIGHJ2*01\t4\n"
    )
    rows = list(csv.DictReader(io.StringIO(fixture), delimiter="\t"))
    assert len(rows) == 2
    assert rows[0]["v_trim_3"] == "2"
    assert rows[1]["j_trim_5"] == "4"


# ──────────────────────────────────────────────────────────────────
# G. pin_present_* — stop-and-report verification (GenAIRR populates fields)
# ──────────────────────────────────────────────────────────────────


def test_pin_present_genairr_populates_vdj_trim_fields_reliably() -> None:
    """**Stop-and-report gate.** GenAIRR's bundled VDJ
    cartridge `HUMAN_IGH_OGRDB` populates the four
    recombination-trim AIRR fields at significant rates
    over a 100-record fixed-seed smoke. If this regresses
    (e.g. the engine starts hard-zeroing one of the four),
    the audit's clean-yes finding is invalidated and the
    estimator implementation must STOP AND REPORT."""
    result = (
        ga.Experiment.on("human_igh")
        .recombine()
        .run_records(n=100, seed=42)
    )
    counts = {
        "v_trim_3": sum(1 for r in result.records if r["v_trim_3"] != 0),
        "d_trim_5": sum(1 for r in result.records if r["d_trim_5"] != 0),
        "d_trim_3": sum(1 for r in result.records if r["d_trim_3"] != 0),
        "j_trim_5": sum(1 for r in result.records if r["j_trim_5"] != 0),
    }
    for field, nonzero in counts.items():
        assert nonzero >= 50, (
            f"VDJ cartridge populates {field} on only {nonzero}/100 "
            f"records (expected ≥ 50). Audit's clean-yes verdict "
            f"assumed the four recombination-trim fields land at "
            f"significant rates; this regression invalidates the "
            f"estimator's input premise — STOP AND REPORT before "
            f"implementing."
        )


def test_pin_present_genairr_populates_vj_trim_fields_reliably() -> None:
    """**Stop-and-report gate, VJ side.** The bundled VJ
    cartridge `HUMAN_IGK_OGRDB` populates `v_trim_3` and
    `j_trim_5` at significant rates; the D fields are
    correctly zero (no D segment on a VJ chain)."""
    result = (
        ga.Experiment.on("human_igk")
        .recombine()
        .run_records(n=100, seed=42)
    )
    v_trim_3_nonzero = sum(1 for r in result.records if r["v_trim_3"] != 0)
    j_trim_5_nonzero = sum(1 for r in result.records if r["j_trim_5"] != 0)
    assert v_trim_3_nonzero >= 50, (
        f"VJ cartridge populates v_trim_3 on only {v_trim_3_nonzero}/100 "
        f"records (expected ≥ 50) — STOP AND REPORT"
    )
    assert j_trim_5_nonzero >= 50, (
        f"VJ cartridge populates j_trim_5 on only {j_trim_5_nonzero}/100 "
        f"records (expected ≥ 50) — STOP AND REPORT"
    )
    # D fields stay zero on VJ chain — no D segment exists.
    for r in result.records:
        assert r["d_trim_5"] == 0, r
        assert r["d_trim_3"] == 0, r


def test_pin_present_v_trim_5_and_j_trim_3_are_hard_zero() -> None:
    """`v_trim_5` and `j_trim_3` are hard-zero on every
    record on both VJ and VDJ cartridges — confirms the
    engine's documented "no V_5 / J_3 trim pass" boundary.
    The estimator MUST NOT consume these AIRR columns."""
    for cfg_name in ("human_igh", "human_igk"):
        result = (
            ga.Experiment.on(cfg_name)
            .recombine()
            .run_records(n=100, seed=42)
        )
        for r in result.records:
            assert r["v_trim_5"] == 0, (
                f"{cfg_name}: v_trim_5 nonzero — engine pipeline now "
                f"runs a V_5 trim pass, audit's hard-zero boundary "
                f"regressed"
            )
            assert r["j_trim_3"] == 0, (
                f"{cfg_name}: j_trim_3 nonzero — engine pipeline now "
                f"runs a J_3 trim pass, audit's hard-zero boundary "
                f"regressed"
            )


# ──────────────────────────────────────────────────────────────────
# H. pin_present_* — documented surface state (NO soft-gap inherited)
# ──────────────────────────────────────────────────────────────────


def test_pin_present_trim_distributions_fold_into_plan_signature() -> None:
    """Trim distributions already fold into the plan
    signature via `fmt_int_dist` — unlike the allele-usage
    slice's documented soft gap 1, the trim plane has NO
    inherited gap to manage. The estimator's output benefits
    from cross-cartridge replay protection automatically.

    Pinned here so a future regression that drops trim
    distributions from the signature (e.g. fold optimisation
    that elides constant distributions) is caught."""
    cfg_a = copy.deepcopy(ga.HUMAN_IGH_OGRDB)
    cfg_b = copy.deepcopy(ga.HUMAN_IGH_OGRDB)
    cfg_a.reference_models = ReferenceEmpiricalModels(
        trims={"V_3": EmpiricalDistributionSpec([(0, 1.0)])}
    )
    cfg_b.reference_models = ReferenceEmpiricalModels(
        trims={"V_3": EmpiricalDistributionSpec([(0, 0.5), (1, 0.5)])}
    )
    exp_a = ga.Experiment.on(cfg_a).recombine()
    exp_b = ga.Experiment.on(cfg_b).recombine()
    ca = exp_a.compile()
    cb = exp_b.compile()
    sa = json.loads(
        ca.simulator.trace_file_from(ca.simulator.run(seed=0), seed=0).to_json()
    )["pass_plan_signature"]
    sb = json.loads(
        cb.simulator.trace_file_from(cb.simulator.run(seed=0), seed=0).to_json()
    )["pass_plan_signature"]
    assert sa != sb, (
        "Two cartridges differing only in trims['V_3'] now produce "
        "EQUAL plan signatures — trim distributions no longer fold "
        "into the signature. The estimator slice now inherits a "
        "soft gap; the audit's no-gap finding is invalidated."
    )


# ──────────────────────────────────────────────────────────────────
# I. pin_scaffold_* — manifest already exposes minimal trim block
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_manifest_exposes_trim_keys_today() -> None:
    """The bundled cartridges already expose minimal
    `trim_keys` (list of typed-plane keys present) and
    `legacy_trim_dicts_present` (bool for the legacy nested
    dict). The estimator slice extends this into a
    structured `trims` block per audit §8.2 — pinned here
    so the existing minimal entries remain backwards
    compatible."""
    m = ga.HUMAN_IGH_OGRDB.cartridge_manifest()
    assert "trim_keys" in m["models"]
    assert isinstance(m["models"]["trim_keys"], list)
    assert "legacy_trim_dicts_present" in m["models"]
    assert m["models"]["legacy_trim_dicts_present"] is True  # bundled has legacy dict


# ──────────────────────────────────────────────────────────────────
# J. pin_absence_* — gaps the implementation slice closes
# ──────────────────────────────────────────────────────────────────


def test_pin_present_estimate_trim_distributions_method_on_builder() -> None:
    """Post-slice — `ReferenceCartridgeBuilder.estimate_trim_distributions`
    is callable. The method signature carries the
    `min_count` / `pseudocount` / `replace` kwargs
    documented in the user brief. Flipped from the prior
    absence pin."""
    assert hasattr(ga.ReferenceCartridgeBuilder, "estimate_trim_distributions")
    sig = inspect.signature(ga.ReferenceCartridgeBuilder.estimate_trim_distributions)
    for kw in ("min_count", "pseudocount", "replace"):
        assert kw in sig.parameters, (
            f"estimate_trim_distributions missing {kw!r} kwarg — user "
            f"brief signature regressed"
        )
    # min_count is an int (1), pseudocount a float (0.0), replace a bool (True).
    assert sig.parameters["min_count"].default == 1
    assert sig.parameters["pseudocount"].default == 0.0
    assert sig.parameters["replace"].default is True


def test_pin_present_min_count_pseudocount_kwargs_owned_by_estimators() -> None:
    """Post-slice — `min_count` lives on every estimator
    (allele-usage, trim, NP-length, NP-base-model,
    P-nucleotide-length). `pseudocount` lives on every
    empirical-distribution estimator (trim, NP-length,
    NP-base-model, P-nucleotide-length) — NOT on
    allele-usage. Pinned so a future estimator slice that
    smuggles either kwarg onto a non-estimator sibling
    method surfaces here."""
    method_owners_min_count: list[str] = []
    method_owners_pseudocount: list[str] = []
    for method_name in dir(ga.ReferenceCartridgeBuilder):
        if method_name.startswith("_"):
            continue
        method = getattr(ga.ReferenceCartridgeBuilder, method_name, None)
        if not callable(method):
            continue
        try:
            sig = inspect.signature(method)
        except (TypeError, ValueError):
            continue
        if "min_count" in sig.parameters:
            method_owners_min_count.append(method_name)
        if "pseudocount" in sig.parameters:
            method_owners_pseudocount.append(method_name)
    assert sorted(method_owners_min_count) == [
        "estimate_allele_usage",
        "estimate_np_base_model",
        "estimate_np_length_distributions",
        "estimate_p_nucleotide_lengths",
        "estimate_trim_distributions",
    ], f"min_count kwarg present on unexpected methods: {method_owners_min_count}"
    assert sorted(method_owners_pseudocount) == [
        "estimate_np_base_model",
        "estimate_np_length_distributions",
        "estimate_p_nucleotide_lengths",
        "estimate_trim_distributions",
    ], (
        f"pseudocount kwarg owned by methods other than "
        f"the four empirical-distribution estimators: "
        f"{method_owners_pseudocount}"
    )


def test_pin_present_trim_models_block_in_manifest() -> None:
    """Post-slice — `manifest['models']['trim_models']`
    block exposes the typed plane's status with the
    documented keys: `keys` / `source` /
    `in_plan_signature` / `legacy_trim_dicts_present` /
    `legacy_fallback`. Flipped from the prior absence pin."""
    m = ga.HUMAN_IGH_OGRDB.cartridge_manifest()
    assert "trim_models" in m["models"], (
        "manifest['models']['trim_models'] block missing — user brief "
        "specified this block name; verify the manifest extension landed"
    )
    block = m["models"]["trim_models"]
    # Baseline (no typed-plane spec on bundled cartridge).
    assert block["keys"] == []
    assert block["source"] == "ReferenceEmpiricalModels.trims"
    # NO inherited soft gap — distinguishes from the allele-usage block.
    assert block["in_plan_signature"] is True
    assert block["legacy_trim_dicts_present"] is True  # bundled has legacy
    assert block["legacy_fallback"] is False
    # The minimal top-level entries stay for backwards compatibility.
    assert "trim_keys" in m["models"]
    assert "legacy_trim_dicts_present" in m["models"]


def test_pin_absence_no_estimate_trim_distributions_module_or_class() -> None:
    """No sibling module / class with this name was
    accidentally introduced. The estimator method lives on
    `ReferenceCartridgeBuilder`, not as a free function."""
    import GenAIRR.cartridge_builder as cb

    for forbidden in ("TrimDistributionEstimator", "estimate_trim_distributions",
                      "TrimEstimator"):
        assert not hasattr(cb, forbidden), (
            f"cartridge_builder.{forbidden} now exists — verify the "
            f"estimator surface is owned by the builder method, not a "
            f"parallel free-function / dataclass surface"
        )


# ──────────────────────────────────────────────────────────────────
# K. Doc anchor
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_audit_doc_exists_and_references_contract() -> None:
    """The audit doc exists and references the contract
    file by name; section structure intact."""
    if not _AUDIT_DOC.exists():
        import pytest
        pytest.skip("docs/ is contributor-only; audit doc not present in this checkout")
    doc = _AUDIT_DOC.read_text(encoding="utf-8")
    assert "test_trim_distribution_estimation_contract.py" in doc, (
        "audit doc no longer references the contract file"
    )
    for marker in (
        "## 1. Q1",
        "## 2. Q2",
        "## 7. Q7",
        "## 9. Implementation order",
        "## 12. Summary table",
    ):
        assert marker in doc, f"audit doc missing section marker {marker!r}"
