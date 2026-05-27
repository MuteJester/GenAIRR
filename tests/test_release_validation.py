"""Release-tier validation test.

Runs realistic full-stack simulations against bundled species
presets and runs two complementary equivalence checks per record:

  1. `result.validate_records(refdata)` — the user-facing AIRR
     postcondition validator. Asks "does each projected AIRR
     record agree with an independent re-derivation from the
     final Outcome?" Empty failure list = the projected output
     is internally consistent.

  2. `outcome.check_live_call_cache_parity(refdata)` — the
     engine-integrity parity harness. Asks "does the cached
     SegmentLiveCall on the final Simulation equal a from-
     scratch recompute over the same sim + refdata?" Empty
     mismatches = the runtime cache that *feeds* projection is
     itself consistent.

Both gates ship green on every release. They're deliberately kept
separate: the validator is the downstream contract for users; the
parity harness is the engine-side guard. A bug in the live-call
refresh path can land in either gate first (or both); having them
side by side localises the failure quickly.

History: the C4 allele-call oracle mirrors the walker's
NP-extension scoring (`live_call/scoring.rs`); two boundary bugs
in the live-call refresh path were closed (`segment_region_overlaps_dirty`
strict-`<` + `on_indel_inserted` `<=`) after this harness surfaced
them. All four locus configurations validate at 100% across the
seed counts run here.
"""
from __future__ import annotations

import GenAIRR as ga


def test_full_stack_vdj_validates_all_airr_records():
    """Productive human IGH full-stack pipeline (the canonical
    realistic config). Two gates run side by side:

      AIRR validator — every projected record must agree with the
        engine's truth oracle (no filtering).
      Cache parity   — every cached SegmentLiveCall must equal a
        from-scratch recompute over the same sim + refdata.

    Validator covers projected output; parity covers the cached
    state feeding it. If a future refresh-path bug regresses, one
    or both gates surface the divergence — start with whichever
    fails first."""
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .productive_only()
        .mutate(rate=0.03)
        .pcr_amplify(rate=1e-4)
        .polymerase_indels(count=2)
        .primer_trim_5prime(length=(0, 3))
        .primer_trim_3prime(length=(0, 3))
    )
    refdata = exp.refdata

    result = exp.run_records(n=100, seed=4242)

    # Gate 1: projected AIRR records.
    report = result.validate_records(refdata)
    assert report, (
        f"validation failed on {len(report.failures)}/{report.count} "
        f"records — summary={report.summary()}; "
        f"first failure={report.failures[0] if report.failures else None}"
    )

    # Gate 2: cached live calls vs from-scratch recompute. Iterate
    # outcomes (they're attached to the result for non-clonal runs)
    # and assert per-segment parity on every record.
    assert result.outcomes is not None, "result must carry outcomes for parity check"
    for i, outcome in enumerate(result.outcomes):
        for p in outcome.check_live_call_cache_parity(refdata):
            assert p["tie_set_matches"], (
                f"record {i} segment {p['segment']}: cached live call "
                f"diverges from fresh recompute — "
                f"cached={p['cached_tie_set']} fresh={p['fresh_tie_set']}"
            )
            if p["hypothesis_bounds_match"] is not None:
                assert p["hypothesis_bounds_match"], (
                    f"record {i} segment {p['segment']}: hypothesis bounds "
                    f"diverge cached={p['cached_hypothesis']} "
                    f"fresh={p['fresh_hypothesis']}"
                )


def test_full_stack_vdj_non_productive_validates():
    """Same VDJ pipeline minus productive_only — ~30% of records
    have junction stops or OOF junctions, but every AIRR record
    must still agree with the validator's re-derivation."""
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(rate=0.03)
        .pcr_amplify(rate=1e-4)
        .polymerase_indels(count=2)
    )
    refdata = exp.refdata
    result = exp.run_records(n=100, seed=4242)
    report = result.validate_records(refdata)
    assert report, (
        f"non-productive VDJ validation failed: {report.summary()}; "
        f"first failure={report.failures[0] if report.failures else None}"
    )


def test_full_stack_validates_igk_with_no_filtering():
    """Productive IGK full stack — SAME shape as IGH but on a VJ
    chain. Previously this fixture surfaced a ~5% J-segment
    tie-set divergence under end-loss + productive_only;
    the boundary fix in segment_region_overlaps_dirty closed it.

    This test is the regression guard: if the fix gets reverted or
    a similar refresh-skip bug appears at the 3' end-loss boundary,
    this fails first."""
    exp = (
        ga.Experiment.on("human_igk")
        .recombine()
        .productive_only()
        .mutate(rate=0.03)
        .pcr_amplify(rate=1e-4)
        .polymerase_indels(count=2)
        .primer_trim_5prime(length=(0, 3))
        .primer_trim_3prime(length=(0, 3))
    )
    refdata = exp.refdata
    result = exp.run_records(n=100, seed=4242)
    report = result.validate_records(refdata)
    assert report, (
        f"IGK full-stack validation failed: {report.summary()}; "
        f"first failure={report.failures[0] if report.failures else None}"
    )


def test_full_stack_validates_lambda_with_no_filtering():
    """IGL (lambda) productive full stack — the second VJ locus."""
    exp = (
        ga.Experiment.on("human_igl")
        .recombine()
        .productive_only()
        .mutate(rate=0.03)
        .pcr_amplify(rate=1e-4)
    )
    refdata = exp.refdata
    result = exp.run_records(n=100, seed=4242)
    report = result.validate_records(refdata)
    assert report, (
        f"IGL validation failed: {report.summary()}; "
        f"first failure={report.failures[0] if report.failures else None}"
    )


def test_end_loss_three_prime_does_not_strand_stale_live_call():
    """Regression test for the segment_region_overlaps_dirty
    boundary bug: a 3' primer-trim of length≥1 with productive_only
    on a productive IGK simulation must produce a J live-call that
    matches the from-scratch recompute oracle on every seed.

    The bug: an IndelDeleted event at exactly pool_len-1 (the
    deleted byte) leaves a dirty window at [pool_len-1, pool_len).
    Post-deletion J.region.end = pool_len-1. Strict-`<` overlap
    check skipped the segment, leaving the stale pre-deletion live
    call committed. Fixed by making the overlap upper bound
    inclusive: `w.start <= region_end`.
    """
    exp = (
        ga.Experiment.on("human_igk")
        .recombine()
        .productive_only()
        .primer_trim_3prime(length=(1, 3))  # always at least 1 byte trimmed
    )
    refdata = exp.refdata
    result = exp.run_records(n=200, seed=4242)
    report = result.validate_records(refdata)
    assert report, (
        f"3'-end-loss boundary regression: {report.summary()}; "
        f"first failure={report.failures[0] if report.failures else None}"
    )


# ──────────────────────────────────────────────────────────────────
# Two-layer integrity model — the canonical example of running
# both gates together on the same outcomes.
# ──────────────────────────────────────────────────────────────────


def test_validator_and_parity_are_independent_layers():
    """Demonstrates the two-layer model end to end:

      Layer 1 (downstream contract): validate_records inspects the
        *projected* AIRR record. This is what users see when they
        consume engine output.

      Layer 2 (engine-side guard): check_live_call_cache_parity
        inspects the *cached* live-call state on the final
        Simulation. This is the state projection reads from.

    Both gates run on every outcome in this test. If a future
    refresh-path bug regresses, the failing layer tells you where
    to look — projection (rare; usually downstream of a cache
    bug) or refresh-path (more common, caught at the source by
    parity). Both green means projection AND cache are consistent.

    The two layers are independent by design: a projection bug
    could pass parity but fail validation; a cache bug usually
    fails parity first and *may* leak into the validator. Keeping
    them separate lets the failure point at the root cause."""
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .productive_only()
        .mutate(rate=0.03)
        .polymerase_indels(count=2)
        .primer_trim_3prime(length=(1, 3))
    )
    refdata = exp.refdata
    result = exp.run_records(n=50, seed=4242)

    report = result.validate_records(refdata)
    assert report, f"projection layer (validator) failed: {report.summary()}"

    assert result.outcomes is not None
    for i, outcome in enumerate(result.outcomes):
        for p in outcome.check_live_call_cache_parity(refdata):
            assert p["tie_set_matches"], (
                f"engine layer (cache parity) failed at record {i} "
                f"segment {p['segment']}: cached={p['cached_tie_set']} "
                f"fresh={p['fresh_tie_set']}"
            )
