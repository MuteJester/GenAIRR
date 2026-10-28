"""Release-tier smoke test for **Allele Usage Estimation v1**.

Proves the end-to-end allele-usage slice composes cleanly
under a release-tier workflow:

1. Build a tiny inline-FASTA VDJ cartridge via the v1 facade.
2. Estimate allele usage from a small AIRR-like list where
   one V allele dominates.
3. Run recombination on the resulting cartridge and observe
   the bias.
4. Assert the manifest exposes the typed plane.
5. Assert the build report carries the estimator stage.
6. Assert an explicit ``recombine(v_allele_weights=...)``
   kwarg still overrides the cartridge plane.

Companion to
[`docs/allele_usage_estimation_design.md`](../docs/allele_usage_estimation_design.md).
"""
from __future__ import annotations

import GenAIRR as ga


# ──────────────────────────────────────────────────────────────────
# Inline VDJ FASTA fixtures — kept small enough to run fast but
# realistic enough to compile through ``Experiment.on(cfg).recombine()``.
# Synthetic anchors are absent; the test compiles under
# ``allow_curatable_refdata()`` per the documented v1 workflow.
# ──────────────────────────────────────────────────────────────────

_V_FASTA = (
    ">IGHV1-MOCK*01\n"
    "GAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTACAGCCTGGGGGGTCCCTGAGACTC"
    "TCCTGTGCAGCCTCTGGATTCACCTTCAGTAGCTATGCCATGAGCTGGGTCCGCCAGGCT\n"
    ">IGHV2-MOCK*01\n"
    "CAGGTCAACTTAAGGGAGTCTGGTCCTGCGCTGGTGAAACCCACACAGACCCTCACACTG"
    "ACCTGCACCTTCTCTGGGTTCTCACTCAGCACTAGTGGAATGTGTGTGAGCTGGATCCGT\n"
)

_D_FASTA = (
    ">IGHD1-MOCK*01\n"
    "GGGTATAGCAGCAGCTGGTAC\n"
    ">IGHD2-MOCK*01\n"
    "AGGATATTGTAGTGGTGGTAGCTGCTACTCC\n"
)

_J_FASTA = (
    ">IGHJ1-MOCK*01\n"
    "TACTACTACGGTATGGACGTCTGGGGCCAAGGGACCACGGTCACCGTCTCCTCAG\n"
    ">IGHJ2-MOCK*01\n"
    "ACTACTGGTACTTCGATCTCTGGGGCCGTGGCACCCTGGTCACTGTCTCCTCAG\n"
)


# Dominant / minor split — 9:1 — large enough to observe at n=100
# without flaking on a single seed (binomial mass at p=0.9 over 100
# trials gives a comfortable margin).
_DOMINANT_V = "IGHV1-MOCK*01"
_MINOR_V = "IGHV2-MOCK*01"
_D_ALLELE = "IGHD1-MOCK*01"
_J_ALLELE = "IGHJ1-MOCK*01"


def _rearrangements_dominant_v1() -> list[dict]:
    """Synthetic AIRR records — 9× dominant V + 1× minor V."""
    return (
        [
            {"v_call": _DOMINANT_V, "d_call": _D_ALLELE, "j_call": _J_ALLELE}
        ] * 9
        + [
            {"v_call": _MINOR_V, "d_call": _D_ALLELE, "j_call": _J_ALLELE}
        ]
    )


def _build_cartridge_with_estimator() -> "ga.DataConfig":
    """Stand the cartridge up through the v1 facade + estimator."""
    return (
        ga.ReferenceCartridgeBuilder
        .from_fasta(
            v_fasta=_V_FASTA,
            d_fasta=_D_FASTA,
            j_fasta=_J_FASTA,
            chain_type="BCR_HEAVY",
        )
        .infer_identity(
            species="HUMAN",
            locus="IGH",
            reference_set="MOCK_RELEASE_AU",
            name="MOCK_IGH_AU",
            source="ReferenceCartridgeBuilder",
        )
        .estimate_allele_usage(
            _rearrangements_dominant_v1(),
            ambiguous="fractional",
            min_count=0.5,
        )
        .build()
    )


# ──────────────────────────────────────────────────────────────────
# 1. End-to-end smoke: estimator → manifest → report → bias visible
# ──────────────────────────────────────────────────────────────────


def test_release_smoke_estimated_allele_usage_biases_recombination_default() -> None:
    """Full v1 happy path:

    - Build the cartridge through the v1 facade.
    - Estimator wrote the dominant V allele into the typed plane.
    - ``manifest["models"]["allele_usage"]["available"]`` is
      ``True`` and ``"V"`` is in ``nonempty_segments``.
    - ``build_report`` carries the ``estimate_allele_usage`` stage.
    - ``Experiment.on(cfg).recombine()`` uses the estimated weights
      by default — the dominant V allele appears strictly more
      often than the minor V allele over n=100 records.
    """
    cfg = _build_cartridge_with_estimator()
    assert isinstance(cfg, ga.DataConfig)

    # Manifest exposes the typed plane.
    au_block = cfg.cartridge_manifest()["models"]["allele_usage"]
    assert au_block["available"] is True, au_block
    assert "V" in au_block["nonempty_segments"], au_block
    # Soft-gap pin still holds at release time.
    assert au_block["in_plan_signature"] is False, au_block

    # Build report carries the estimator stage.
    report = cfg.build_report
    assert report is not None
    stage_names = [s["stage"] for s in report.stages]
    assert "estimate_allele_usage" in stage_names, stage_names

    # Recombination on the cartridge sees the bias by default.
    result = (
        ga.Experiment.on(cfg)
        .allow_curatable_refdata()
        .recombine()
        .run_records(n=100, seed=4242)
    )
    assert len(result.records) == 100

    dominant_hits = sum(
        1 for r in result.records
        if _DOMINANT_V in (r.get("v_call") or "")
    )
    minor_hits = sum(
        1 for r in result.records
        if _MINOR_V in (r.get("v_call") or "")
    )
    # 9:1 input distribution → dominant > minor with very high
    # probability at n=100. Allow plenty of margin so the test
    # doesn't flake on tie-set ambiguity.
    assert dominant_hits > minor_hits, (
        f"estimator-biased recombination didn't favour the dominant "
        f"V allele: dominant_hits={dominant_hits}, minor_hits="
        f"{minor_hits}. Verify the cartridge plane is reaching the "
        f"engine through the precedence chain."
    )


# ──────────────────────────────────────────────────────────────────
# 2. Explicit recombine kwarg overrides the cartridge plane
# ──────────────────────────────────────────────────────────────────


def test_release_smoke_explicit_recombine_kwarg_overrides_cartridge_plane() -> None:
    """Control case: even with the cartridge plane biased toward
    ``IGHV1-MOCK*01``, an explicit
    ``recombine(v_allele_weights={IGHV2-MOCK*01: 100.0})`` flips
    the dominance — proves the precedence kwarg > plane > uniform
    still holds at release time."""
    cfg = _build_cartridge_with_estimator()

    result = (
        ga.Experiment.on(cfg)
        .allow_curatable_refdata()
        .recombine(v_allele_weights={_MINOR_V: 100.0})
        .run_records(n=100, seed=4242)
    )

    minor_hits = sum(
        1 for r in result.records
        if _MINOR_V in (r.get("v_call") or "")
    )
    dominant_hits = sum(
        1 for r in result.records
        if _DOMINANT_V in (r.get("v_call") or "")
    )
    assert minor_hits > dominant_hits, (
        f"explicit recombine kwarg didn't override the cartridge "
        f"plane: minor_hits={minor_hits}, dominant_hits="
        f"{dominant_hits}. Precedence kwarg > plane > uniform is "
        f"broken — verify Experiment.recombine's resolution order."
    )
