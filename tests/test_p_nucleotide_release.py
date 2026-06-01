"""Release-tier composition tests for the P-nucleotide v1 slice.

Confirms the slice composes with the rest of the release stack
(productive IGH, D inversion, receptor revision, SHM,
corruption, AIRR projection, validator, cache parity) and
preserves three load-bearing invariants:

- **Records validate clean under the full stack.** Productive
  IGH with all four P-ends configured + ``invert_d`` + optional
  receptor revision + SHM + corruption produces records whose
  AIRR projection passes
  :py:meth:`SimulationResult.validate_records` AND whose
  cached live-call state matches a from-scratch parity
  recompute.

- **Replay round-trip is byte-identical with ``invert_d``
  fired.** The IR-correct lowering order `invert_d →
  p_addition.D_5 → assemble.D` is exercised end-to-end by
  forcing ``invert_d(prob=1.0)``; same-seed re-execution
  reproduces ``sequence`` and all four ``p_*_length`` fields.

- **Palindrome derivation is byte-exact against the allele.**
  A fixed P-length fixture restricted to a single V allele
  produces P bytes that match the reverse-complement of the
  V's post-trim 3' coding flank — the deterministic property
  the audit guaranteed.

Companion to [`docs/p_nucleotide_design.md`](../docs/p_nucleotide_design.md).
"""
from __future__ import annotations

import copy
import json

import pytest

import GenAIRR as ga
from GenAIRR._airr_record import outcome_to_airr_record
from GenAIRR.reference_models import (
    EmpiricalDistributionSpec,
    ReferenceEmpiricalModels,
)


# ──────────────────────────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────────────────────────


def _complement(byte: str) -> str:
    return {"A": "T", "T": "A", "C": "G", "G": "C"}.get(byte.upper(), byte)


def _expected_p_v_3_bytes(v_seq: str, v_trim_3: int, length: int) -> str:
    """The audit's palindrome rule for V_3 in pool order:
    P[i] = complement(effective_v[len-1-i]) for i in [0, length)
    where effective_v is the post-trim V coding sequence.
    """
    effective = v_seq[: len(v_seq) - v_trim_3]
    return "".join(_complement(effective[-1 - i]) for i in range(length))


def _cartridge_with_p_lengths(p_lengths: dict, preset: str = "HUMAN_IGH_OGRDB") -> "ga.DataConfig":
    """Deep-copy a bundled preset and attach the requested typed
    P-plane while preserving any existing NP / trim / NP-base
    models on the source cartridge."""
    cfg = copy.deepcopy(getattr(ga, preset))
    existing = getattr(cfg, "reference_models", None)
    if isinstance(existing, ReferenceEmpiricalModels):
        cfg.reference_models = ReferenceEmpiricalModels(
            np_lengths=existing.np_lengths,
            trims=existing.trims,
            np_bases=existing.np_bases,
            p_nucleotide_lengths=p_lengths,
        )
    else:
        cfg.reference_models = ReferenceEmpiricalModels(
            p_nucleotide_lengths=p_lengths,
        )
    return cfg


_FOUR_END_RELEASE_LENGTHS = {
    "V_3": EmpiricalDistributionSpec([(0, 0.4), (1, 0.3), (2, 0.2), (3, 0.1)]),
    "D_5": EmpiricalDistributionSpec([(0, 0.5), (1, 0.3), (2, 0.2)]),
    "D_3": EmpiricalDistributionSpec([(0, 0.5), (1, 0.3), (2, 0.2)]),
    "J_5": EmpiricalDistributionSpec([(0, 0.4), (1, 0.3), (2, 0.2), (3, 0.1)]),
}


# ──────────────────────────────────────────────────────────────────
# 1. Productive IGH full stack: validate_records + cache parity
# ──────────────────────────────────────────────────────────────────


def test_productive_igh_full_stack_with_p_addition_validates_and_parity_holds() -> None:
    """Two-layer integrity check under the full release stack:

    - Layer 1 (downstream contract): ``validate_records``
      inspects the projected AIRR records. The four new
      ``p_*_length`` fields participate in the validator's
      event-ledger recompute via ``PLengthMismatch``.
    - Layer 2 (engine-side guard):
      ``check_live_call_cache_parity`` compares cached
      live-call state to a from-scratch recompute. P-bytes
      sit between assembled segments and carry the
      ``P_NUC`` flag; they do not shift segment regions and
      must NOT invalidate cached calls.

    Pipeline: productive IGH with all four P-ends + D
    inversion + receptor revision + SHM + polymerase indels +
    primer trim. At least one record must carry a non-zero
    P field across the batch — otherwise the test isn't
    actually exercising P-addition.
    """
    cfg = _cartridge_with_p_lengths(_FOUR_END_RELEASE_LENGTHS)
    exp = (
        ga.Experiment.on(cfg)
        .recombine()
        .productive_only()
        .invert_d(prob=0.3)
        .receptor_revision(prob=0.2)
        .mutate(rate=0.01)
        .polymerase_indels(count=1)
        .primer_trim_3prime(length=(1, 2))
    )
    refdata = exp.refdata
    result = exp.run_records(n=40, seed=4242)

    report = result.validate_records(refdata)
    assert report, (
        f"projection layer (validator) failed under productive IGH + "
        f"P-addition full stack: {report.summary()}"
    )

    assert result.outcomes is not None
    for i, outcome in enumerate(result.outcomes):
        for p in outcome.check_live_call_cache_parity(refdata):
            assert p["tie_set_matches"], (
                f"engine layer (cache parity) failed at record {i} "
                f"segment {p['segment']}: cached={p['cached_tie_set']} "
                f"fresh={p['fresh_tie_set']}"
            )

    # Sanity: P-addition is genuinely exercised across the batch.
    any_p_emitted = any(
        rec["p_v_3_length"] + rec["p_d_5_length"]
        + rec["p_d_3_length"] + rec["p_j_5_length"] > 0
        for rec in result.records
    )
    assert any_p_emitted, (
        "no P bytes emitted across 40 records — the typed P-plane "
        "isn't reaching the lowering layer, or every length sample "
        "hit zero by accident. Verify the cartridge attachment."
    )


# ──────────────────────────────────────────────────────────────────
# 2. Replay round-trip with invert_d(prob=1.0)
# ──────────────────────────────────────────────────────────────────


def test_replay_round_trip_under_forced_d_inversion_preserves_sequence_and_p_lengths() -> None:
    """Forces ``invert_d(prob=1.0)`` so every record exercises
    the IR-critical ``invert_d → p_addition.D_5 → assemble.D``
    lowering order. Replay round-trip via
    ``rerun_from_trace_file`` must reproduce ``sequence`` and
    all four ``p_*_length`` fields byte-for-byte across
    multiple seeds.

    The D_5 palindrome derives from D's post-inversion
    effective_seq; replay reconstructs the inversion + trim +
    length state from the trace and re-derives the same P
    bytes. Any IR-ordering regression that drops the
    `invert_d` commit before `p_addition.D_5` would produce
    different sequences on replay and fail this gate."""
    cfg = _cartridge_with_p_lengths(_FOUR_END_RELEASE_LENGTHS)
    exp = ga.Experiment.on(cfg).recombine().productive_only().invert_d(prob=1.0)
    refdata = exp.refdata
    compiled = exp.compile()
    seen_d_5_nonzero = False
    for seed in range(8):
        fresh = compiled.simulator.run(seed=seed)
        tf = compiled.simulator.trace_file_from(fresh, seed=seed)
        replayed = compiled.simulator.rerun_from_trace_file(tf)

        fr = outcome_to_airr_record(fresh, refdata, sequence_id=f"fresh-{seed}")
        rr = outcome_to_airr_record(replayed, refdata, sequence_id=f"replay-{seed}")
        for field in (
            "sequence",
            "p_v_3_length",
            "p_d_5_length",
            "p_d_3_length",
            "p_j_5_length",
        ):
            assert fr[field] == rr[field], (
                f"seed {seed}: P-addition round-trip desynced on "
                f"field {field!r} under forced D inversion: "
                f"fresh={fr[field]!r}, replay={rr[field]!r}"
            )
        # Every record is D-inverted by construction.
        assert fr["d_inverted"] is True, (
            f"seed {seed}: invert_d(prob=1.0) should produce d_inverted=True"
        )
        if fr["p_d_5_length"] > 0:
            seen_d_5_nonzero = True
    # Across 8 seeds at least one D_5 P-extension fires.
    assert seen_d_5_nonzero, (
        "no seed in [0..8) produced a non-zero p_d_5_length — the "
        "IR-critical post-inversion D_5 path isn't actually exercised"
    )


# ──────────────────────────────────────────────────────────────────
# 3. Deterministic palindrome derivation on a single V allele
# ──────────────────────────────────────────────────────────────────


def test_deterministic_palindrome_at_v_3_matches_allele_reverse_complement() -> None:
    """Restrict the V allele pool to a single named allele and
    force ``p_nucleotide_lengths["V_3"] = [(3, 1.0)]`` so every
    record's P_V3 region is exactly 3 bytes long. Then verify
    that the P bytes in the assembled ``sequence``, sitting
    immediately after ``v_sequence_end``, equal the
    reverse-complement of the V allele's last 3 post-trim
    coding bytes — the audit's deterministic palindrome
    invariant.

    VJ chain (HUMAN_IGK) so there's no D inversion to thread
    through. NP1 length is sampled from the cartridge's
    typed plane; we read the V's allele sequence + trim from
    the AIRR record and compute the expected palindrome
    independently."""
    cfg = _cartridge_with_p_lengths(
        {"V_3": EmpiricalDistributionSpec([(3, 1.0)])},
        preset="HUMAN_IGK_OGRDB",
    )
    # Pick the first allele of the first V gene to lock the
    # test to. `cfg.v_alleles` is keyed by gene name; the value
    # is a list of `Allele` instances.
    first_gene = sorted(cfg.v_alleles)[0]
    v_allele = cfg.v_alleles[first_gene][0]
    v_allele_name = v_allele.name
    v_seq = v_allele.ungapped_seq.upper()
    exp = (
        ga.Experiment.on(cfg)
        .restrict_alleles(v=v_allele_name)
        .recombine()
    )
    result = exp.run_records(n=10, seed=4242)
    for rec in result.records:
        # `v_call` may carry a tie set (comma-separated allele
        # names) when the live-call walker can't disambiguate
        # the assigned allele from a sibling — the cartridge
        # may have alleles with identical post-trim 3' flanks.
        # The restricted sample is guaranteed to be our allele
        # though, so the palindrome check against `v_seq` is
        # the load-bearing assertion. We just confirm our
        # allele is among the call set.
        assert v_allele_name in rec["v_call"].split(","), rec["v_call"]
        assert rec["p_v_3_length"] == 3, rec
        v_end = rec["v_sequence_end"]
        v_trim_3 = rec["v_trim_3"]
        observed = rec["sequence"][v_end : v_end + 3].upper()
        expected = _expected_p_v_3_bytes(v_seq, v_trim_3, 3)
        assert observed == expected, (
            f"V_3 palindrome mismatch on {rec['sequence_id']}: "
            f"expected {expected!r} (reverse-complement of "
            f"v[-3:] after trim_3={v_trim_3}), got {observed!r}. "
            f"V allele {v_allele_name} length {len(v_seq)}, "
            f"post-trim flank {v_seq[len(v_seq)-v_trim_3-3:len(v_seq)-v_trim_3]!r}"
        )


# ──────────────────────────────────────────────────────────────────
# 4. Fixed P-length fixture produces exact field values
# ──────────────────────────────────────────────────────────────────


def test_fixed_p_length_fixture_produces_exact_aggregate_field_values() -> None:
    """Singleton P-length distributions per end produce
    records whose four `p_*_length` AIRR fields exactly match
    the fixture values on every record — same boundary the
    cartridge's `EmpiricalDistributionSpec` enforces at the
    sampler. Pinned across many seeds to lock in the
    "single-value distribution = constant output" property."""
    cfg = _cartridge_with_p_lengths(
        {
            "V_3": EmpiricalDistributionSpec([(2, 1.0)]),
            "D_5": EmpiricalDistributionSpec([(3, 1.0)]),
            "D_3": EmpiricalDistributionSpec([(1, 1.0)]),
            "J_5": EmpiricalDistributionSpec([(4, 1.0)]),
        }
    )
    result = ga.Experiment.on(cfg).recombine().run_records(n=12, seed=0)
    for rec in result.records:
        assert rec["p_v_3_length"] == 2, rec
        assert rec["p_d_5_length"] == 3, rec
        assert rec["p_d_3_length"] == 1, rec
        assert rec["p_j_5_length"] == 4, rec
