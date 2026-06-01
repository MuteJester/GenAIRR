"""Refdata-validation compile gate — Python surface.

The compile-time gate enforces structural correctness of reference
data before any pass runs. Failures must surface as ``ValueError``
from ``Experiment.compile()`` (and from the lower-level
``CompiledSimulator.compile_from_plan`` path) with the full
aggregated issue list embedded in the message.

These tests complement the cargo-side gate tests (`compiled::tests::
compile_runtime::refdata_validation_gate`); together they prove the
gate fires uniformly across the Rust API and the Python entry point.
"""
from __future__ import annotations

import pytest

import GenAIRR as ga


# ──────────────────────────────────────────────────────────────────
# Bundled presets — the gate must NOT block valid production data.
# ──────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("preset", ["human_igh", "human_igk", "human_igl"])
def test_bundled_preset_compiles_through_gate(preset: str) -> None:
    # Minimal pipeline: just recombine. If the gate is overzealous,
    # this is where it surfaces.
    exp = ga.Experiment.on(preset).recombine()
    compiled = exp.compile()
    assert compiled is not None


# ──────────────────────────────────────────────────────────────────
# Synthetic refdata — gate must fire on malformed input.
# ──────────────────────────────────────────────────────────────────


def _custom_vdj(seq_v: bytes = b"TGTAAACCC",
                anchor_v: int | None = 0,
                v_name: str = "IGHV1-1*01",
                add_dup_j: bool = False,
                extra_v_byte: bytes | None = None):
    """Build a deliberately-customisable VDJ refdata for failure tests."""
    from GenAIRR._engine import RefDataConfig

    cfg = RefDataConfig.vdj()
    cfg.add_v_allele(v_name, "IGHV1-1", seq_v, anchor=anchor_v)
    if extra_v_byte is not None:
        cfg.add_v_allele("IGHV1-2*01", "IGHV1-2", extra_v_byte, anchor=0)
    cfg.add_d_allele("IGHD1-1*01", "IGHD1-1", b"GGGCCCAAA")
    cfg.add_j_allele("IGHJ1*01", "IGHJ1", b"TGGAAACCC", anchor=0)
    if add_dup_j:
        cfg.add_j_allele("IGHJ1*01", "IGHJ1", b"TGGAAACCC", anchor=0)
    return cfg


def test_experiment_compile_rejects_v_anchor_non_cys() -> None:
    # GGG (Gly) at anchor; not Cys.
    cfg = _custom_vdj(seq_v=b"GGGAAACCC", anchor_v=0)
    exp = ga.Experiment.on(cfg).recombine()
    with pytest.raises(ValueError) as exc:
        exp.compile()
    msg = str(exc.value)
    assert "refdata validation issue" in msg
    assert "V allele" in msg
    # Spec point 4: error mentions the kind ("expected C").
    assert "expected C" in msg


def test_experiment_compile_rejects_invalid_byte_in_v_seq() -> None:
    # '.' is rejected; gaps must be resolved before refdata.
    cfg = _custom_vdj(seq_v=b"TGT.ACCCC", anchor_v=0)
    exp = ga.Experiment.on(cfg).recombine()
    with pytest.raises(ValueError) as exc:
        exp.compile()
    msg = str(exc.value)
    assert "refdata validation issue" in msg
    # Mentions byte / position info.
    assert "byte at position" in msg


def test_experiment_compile_rejects_duplicate_allele_name() -> None:
    cfg = _custom_vdj(add_dup_j=True)
    exp = ga.Experiment.on(cfg).recombine()
    with pytest.raises(ValueError) as exc:
        exp.compile()
    msg = str(exc.value)
    assert "refdata validation issue" in msg
    assert "duplicate allele name" in msg
    assert "IGHJ1*01" in msg


def test_experiment_compile_rejects_missing_v_anchor() -> None:
    cfg = _custom_vdj(anchor_v=None)
    exp = ga.Experiment.on(cfg).recombine()
    with pytest.raises(ValueError) as exc:
        exp.compile()
    msg = str(exc.value)
    assert "refdata validation issue" in msg
    assert "has no anchor" in msg


def test_experiment_compile_rejects_anchor_out_of_bounds() -> None:
    # Anchor at position 7 of a 9-base sequence → codon needs [7..10)
    # → out of bounds.
    cfg = _custom_vdj(seq_v=b"AAACCCTGT", anchor_v=7)
    exp = ga.Experiment.on(cfg).recombine()
    with pytest.raises(ValueError) as exc:
        exp.compile()
    assert "anchor at 7" in str(exc.value)


def test_experiment_compile_error_mentions_issue_count() -> None:
    # Two problems at once: Gly V anchor + duplicate J name. The
    # gate must aggregate both into one error.
    cfg = _custom_vdj(seq_v=b"GGGAAACCC", anchor_v=0, add_dup_j=True)
    exp = ga.Experiment.on(cfg).recombine()
    with pytest.raises(ValueError) as exc:
        exp.compile()
    msg = str(exc.value)
    assert "issue(s)" in msg
    # Aggregated form: both kinds visible in one error.
    assert "expected C" in msg
    assert "duplicate allele name" in msg


# ──────────────────────────────────────────────────────────────────
# Compile-time error vs runtime error — the productive failure-mode
# tests must still distinguish ValueError (compile gate) from
# StrictSamplingError (runtime). The gate operates on REFDATA
# structure, not on sampling-time constraint mass.
# ──────────────────────────────────────────────────────────────────


# ──────────────────────────────────────────────────────────────────
# Curation policy — Strict (default) vs AllowCuratable.
# ──────────────────────────────────────────────────────────────────


def _curatable_only_vdj():
    """VDJ refdata with only curatable issues:
    - V anchor codon is GGG (Gly), not Cys.
    - J anchor missing entirely.
    No structural Fatal issues.
    """
    from GenAIRR._engine import RefDataConfig

    cfg = RefDataConfig.vdj()
    cfg.add_v_allele("IGHV1-1*01", "IGHV1-1", b"AAACCCGGG", anchor=6)
    cfg.add_d_allele("IGHD1-1*01", "IGHD1-1", b"GGGCCCAAA")
    cfg.add_j_allele("IGHJ-orphan*01", "IGHJ-orphan", b"TGGAAACCC", anchor=None)
    return cfg


def test_strict_mode_rejects_curatable_only_refdata() -> None:
    """Default `compile()` is strict — pseudogene-shape issues block
    even when there are no structural Fatal problems."""
    cfg = _curatable_only_vdj()
    exp = ga.Experiment.on(cfg).recombine()
    with pytest.raises(ValueError) as exc:
        exp.compile()
    msg = str(exc.value)
    # Severity-tagged.
    assert "[curatable]" in msg
    # Remediation hint surfaces because EVERY issue is curatable.
    assert "allow_curatable_refdata" in msg
    assert "pseudogene/ORF" in msg


def test_allow_curatable_accepts_curatable_only_refdata() -> None:
    """``Experiment.allow_curatable_refdata()`` opts the gate into
    the lenient mode — curatable issues now pass."""
    cfg = _curatable_only_vdj()
    compiled = (
        ga.Experiment.on(cfg).allow_curatable_refdata().recombine().compile()
    )
    assert compiled is not None


def test_allow_curatable_still_rejects_invalid_byte() -> None:
    """Spec point 5: invalid byte stays fatal even under
    ``allow_curatable_refdata``."""
    cfg = _curatable_only_vdj()
    cfg.add_v_allele("IGHV-bad*01", "IGHV-bad", b"TGT.AAACC", anchor=0)
    exp = ga.Experiment.on(cfg).allow_curatable_refdata().recombine()
    with pytest.raises(ValueError) as exc:
        exp.compile()
    msg = str(exc.value)
    assert "[fatal]" in msg
    assert "byte at position" in msg


def test_allow_curatable_still_rejects_duplicate_allele_name() -> None:
    cfg = _curatable_only_vdj()
    cfg.add_v_allele("IGHV1-1*01", "IGHV1-1", b"AAACCCGGG", anchor=6)  # duplicate
    exp = ga.Experiment.on(cfg).allow_curatable_refdata().recombine()
    with pytest.raises(ValueError) as exc:
        exp.compile()
    msg = str(exc.value)
    assert "[fatal]" in msg
    assert "duplicate allele name" in msg


def test_allow_curatable_passes_via_compile_kwarg() -> None:
    """Per-call kwarg overrides the instance default."""
    cfg = _curatable_only_vdj()
    exp = ga.Experiment.on(cfg).recombine()
    # No fluent opt-in, but per-call kwarg flips the mode.
    compiled = exp.compile(allow_curatable_refdata=True)
    assert compiled is not None


# ──────────────────────────────────────────────────────────────────
# Bundled mouse / TCR — include pseudogenes; gate behaviour test.
# ──────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("preset", ["mouse_igh", "human_tcrb"])
def test_bundled_pseudogene_preset_fails_strict_with_clear_message(preset: str) -> None:
    """Bundled mouse_igh and human_tcrb include pseudogene/ORF alleles
    whose anchor codons don't translate to Cys / W / F. Strict mode
    must reject and the error must mention the curation remediation."""
    exp = ga.Experiment.on(preset).recombine()
    with pytest.raises(ValueError) as exc:
        exp.compile()
    msg = str(exc.value)
    assert "[curatable]" in msg
    assert "allow_curatable_refdata" in msg


@pytest.mark.parametrize("preset", ["mouse_igh", "human_tcrb"])
def test_bundled_pseudogene_preset_compiles_under_allow_curatable(preset: str) -> None:
    """The same bundled presets compile cleanly under the opt-in."""
    compiled = (
        ga.Experiment.on(preset).allow_curatable_refdata().recombine().compile()
    )
    assert compiled is not None


# ──────────────────────────────────────────────────────────────────
# RefDataConfig.validate() — issue dicts carry severity tag.
# ──────────────────────────────────────────────────────────────────


def test_validate_dict_carries_severity_tag() -> None:
    """Each issue dict has a ``severity`` key with value
    ``'fatal'`` or ``'curatable'``."""
    cfg = _curatable_only_vdj()
    cfg.add_v_allele("IGHV-bad*01", "IGHV-bad", b"TGT.AAACC", anchor=0)
    issues = cfg.validate()
    severities = {i["severity"] for i in issues}
    assert severities == {"fatal", "curatable"}
    # Invalid byte → fatal.
    byte_issues = [i for i in issues if i["kind"] == "InvalidAlleleByte"]
    assert byte_issues
    assert byte_issues[0]["severity"] == "fatal"
    # V anchor not Cys → curatable.
    anchor_issues = [i for i in issues if i["kind"] == "VAnchorNotCys"]
    assert anchor_issues
    assert anchor_issues[0]["severity"] == "curatable"


def test_validate_with_mode_python_binding() -> None:
    """``RefDataConfig.validate_with_mode`` dispatches to the same
    underlying gate as the compile path."""
    cfg = _curatable_only_vdj()
    # Strict — raises.
    with pytest.raises(ValueError):
        cfg.validate_with_mode(mode="strict")
    # AllowCuratable — passes.
    cfg.validate_with_mode(mode="allow_curatable")
    # Unknown mode — argument error.
    with pytest.raises(ValueError, match="unknown validation mode"):
        cfg.validate_with_mode(mode="bogus")


def test_compile_gate_does_not_swallow_strict_runtime_errors() -> None:
    """Bundled refdata + a strict productive contract that exhausts
    sampling mass should still surface as a runtime
    ``StrictSamplingError``, NOT as a compile-time ``ValueError``.
    The gate must only catch refdata structure problems.
    """
    from GenAIRR._engine import StrictSamplingError

    # Valid refdata + an impossible-at-runtime productive setup is
    # exactly what the bundled productive_only pipeline runs, so
    # compile passes. Use bundled IGH (clean refdata) + productive_only
    # — if the runtime path can sample, it does; if it can't, it
    # raises StrictSamplingError. Either way, compile() must NOT
    # raise.
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .productive_only()
    )
    compiled = exp.compile()  # must not raise ValueError
    # Either runtime succeeds or it raises the runtime error class —
    # never ValueError from the compile gate.
    try:
        compiled.run(n=1, seed=42)
    except StrictSamplingError:
        pass  # runtime-side error is fine; the gate didn't swallow it
    except ValueError as e:
        pytest.fail(f"runtime failure should not be ValueError: {e}")
