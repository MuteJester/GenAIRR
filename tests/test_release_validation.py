"""Release-tier validation test.

Runs realistic full-stack simulations against bundled species
presets and validates every projected AIRR record through the
engine's postcondition oracle. The load-bearing release gate.

If these tests fail, the failure summary names the issue kinds and
record indices to inspect — start there before profiling or rolling
back.

The C4 allele-call oracle mirrors the walker's NP-extension
scoring (see `live_call/scoring.rs::score_alleles_with_extensions`),
so non-zero trim_5/trim_3 doesn't cause spurious tie-set mismatches.
The IGK J residual under end-loss + productive was traced to a
boundary case in `segment_region_overlaps_dirty` (strict `<` on
the upper bound discarded refresh signals when end-loss 3'
shrunk a region's end to equal the dirty window's site) and fixed
in `ir/builder.rs`. All four locus configurations now validate
at 100% across 500 seeds.
"""
from __future__ import annotations

import GenAIRR as ga


def test_full_stack_vdj_validates_all_airr_records():
    """Productive human IGH full-stack pipeline (the canonical
    realistic config). Every projected AIRR record must validate
    with NO filtering — the extension-aware oracle agrees with the
    walker on V, D, and J calls under arbitrary trim."""
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
    report = result.validate_records(refdata)

    assert report, (
        f"validation failed on {len(report.failures)}/{report.count} "
        f"records — summary={report.summary()}; "
        f"first failure={report.failures[0] if report.failures else None}"
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
