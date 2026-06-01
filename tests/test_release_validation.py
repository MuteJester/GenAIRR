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


# ──────────────────────────────────────────────────────────────────
# D inversion — release-tier inclusion. Pin that the new mechanism
# survives both the AIRR validator and the cache-parity oracle under
# realistic pipeline loads, on a heavy-chain catalogue with every
# corruption stage active.
#
# Locks down the Slice A → E arc against regressions: any bug that
# corrupts the orientation flag, lets `d_inverted` desync from the
# IR, or breaks live-call refresh under inverted D would surface
# here before reaching downstream consumers.
#
# **D-tie oracle under inversion — RESOLVED.** Earlier slices of
# the D-inversion arc carried a documented limitation: the walker
# scored assembled D bytes against the forward-orientation index,
# so the post-Slice-E inverted-D records surfaced spurious
# `AlleleCallTieSetMismatch{segment: D}` issues. The release
# tests below previously routed through a `_strip_inverted_d_tie_set_issues`
# helper that filtered those.
#
# The D-inversion live-call / allele-call cleanup slice replaced
# the boundary primitive with the orientation-aware
# `matches_observed_with_orientation` (and routed the walker, the
# extension walks, the `from_existing_region` rebuild, and the
# validator oracle through it). The filter is no longer needed —
# the validator now agrees with the walker on inverted-D tie sets
# directly. The previous helper has been removed; the release
# tests below assert validator cleanliness without exceptions.
# ──────────────────────────────────────────────────────────────────


def test_productive_igh_full_stack_with_d_inversion_validates_all_records():
    """Productive IGH + every corruption stage + ``invert_d(prob=1.0)``.

    With prob=1, every record's D is committed in reverse-complement
    orientation, so every record exercises the inverted-D code path
    through assembly, walker scoring, AIRR projection, and the
    validator's `d_inverted` consistency check. Both two-layer
    gates must run clean.
    """
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .invert_d(prob=1.0)
        .productive_only()
        .mutate(rate=0.03)
        .pcr_amplify(rate=1e-4)
        .polymerase_indels(count=2)
        .primer_trim_5prime(length=(0, 3))
        .primer_trim_3prime(length=(0, 3))
    )
    refdata = exp.refdata
    result = exp.run_records(n=50, seed=4242)

    # Gate 1: AIRR validator. Every postcondition runs clean now
    # that the orientation-aware walker + validator oracle
    # converge on the inverted-D tie set. No exceptions stripped.
    report = result.validate_records(refdata)
    assert report, (
        f"inverted-D productive IGH validation failed: "
        f"summary={report.summary()}; first failure="
        f"{report.failures[0] if report.failures else None}"
    )

    # Gate 2: live-call cache parity on every outcome.
    assert result.outcomes is not None
    for i, outcome in enumerate(result.outcomes):
        for p in outcome.check_live_call_cache_parity(refdata):
            assert p["tie_set_matches"], (
                f"inverted-D record {i} segment {p['segment']}: "
                f"cache parity failed under inversion — "
                f"cached={p['cached_tie_set']} fresh={p['fresh_tie_set']}"
            )

    # Every record's d_inverted must be True (prob=1 commits RC on
    # every seed). Acts as a smoke check that the AIRR field
    # plumbing fires for the realistic stack.
    assert all(rec["d_inverted"] is True for rec in result.records), (
        "prob=1.0 must commit ReverseComplement on every record"
    )


def test_non_productive_igh_full_stack_with_d_inversion_validates_all_records():
    """Same as the productive variant minus the constraint bundle.

    Without `productive_only`, junctions can carry in-frame stops or
    out-of-frame lengths. The inverted-D bytes propagate into the
    junction-scan path and the validator must still agree with the
    engine's re-derivation on every record — pins that inversion
    doesn't break the validator's structural / counter checks even
    when the record is non-productive."""
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .invert_d(prob=1.0)
        .mutate(rate=0.03)
        .pcr_amplify(rate=1e-4)
        .polymerase_indels(count=2)
    )
    refdata = exp.refdata
    result = exp.run_records(n=50, seed=4242)
    report = result.validate_records(refdata)
    assert report, (
        f"inverted-D non-productive IGH validation failed: "
        f"summary={report.summary()}; first failure="
        f"{report.failures[0] if report.failures else None}"
    )
    assert all(rec["d_inverted"] is True for rec in result.records)


def test_invert_d_trace_replay_round_trip_passes_validator():
    """Replay determinism with inversion in the chain.

    For each fresh outcome: build a TraceFile, rerun it via
    `rerun_from_trace_file`, build the AIRR record from the replayed
    outcome, and assert it (a) carries the same `d_inverted` value
    and (b) passes the validator.

    Uses ``prob=0.5`` so the seed sweep exercises both the True and
    False branches of the inversion decision; the replay must
    reproduce whichever branch fired in the fresh run."""
    from GenAIRR._airr_record import outcome_to_airr_record

    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .invert_d(prob=0.5)
        .productive_only()
        .mutate(rate=0.01)
    )
    refdata = exp.refdata
    compiled = exp.compile()
    seen_true = False
    seen_false = False
    for seed in range(8):
        fresh = compiled.simulator.run(seed=seed)
        tf = compiled.simulator.trace_file_from(fresh, seed=seed)
        replayed = compiled.simulator.rerun_from_trace_file(tf)

        fresh_rec = outcome_to_airr_record(fresh, refdata, sequence_id=f"fresh-{seed}")
        replayed_rec = outcome_to_airr_record(
            replayed, refdata, sequence_id=f"replay-{seed}"
        )

        # Orientation must round-trip.
        assert fresh_rec["d_inverted"] == replayed_rec["d_inverted"], (
            f"seed {seed}: d_inverted desynced through replay "
            f"({fresh_rec['d_inverted']} vs {replayed_rec['d_inverted']})"
        )
        if fresh_rec["d_inverted"]:
            seen_true = True
        else:
            seen_false = True

        # The full AIRR record's sequence bytes must round-trip too.
        assert fresh_rec["sequence"] == replayed_rec["sequence"], (
            f"seed {seed}: assembled sequence diverged under replay"
        )

        # And the replayed record must pass the validator unchanged.
        replayed_issues = replayed.validate_record(refdata, sequence_id=f"replay-{seed}")
        kinds = {issue["kind"] for issue in replayed_issues}
        assert "DInvertedMismatch" not in kinds, (
            f"seed {seed}: replayed record tripped DInvertedMismatch"
        )

    # Both branches must fire across the 8 seeds at prob=0.5 — if
    # not, the test isn't actually exercising the True/False split.
    assert seen_true, "prob=0.5 over 8 seeds: no True branch ever fired"
    assert seen_false, "prob=0.5 over 8 seeds: no False branch ever fired"


# ──────────────────────────────────────────────────────────────────
# Receptor revision (audit-first biology mechanism #2). Same shape
# as the D-inversion release tests above: validator + cache parity
# under a realistic productive IGH stack, plus a replay round-trip
# that pins both AIRR provenance fields. The distribution invariant
# (`prob=0.25`, ±5σ) lives in `test_distribution_invariants.py`
# alongside the equivalent `invert_d` Bernoulli draw.
# ──────────────────────────────────────────────────────────────────


def test_productive_igh_full_stack_with_receptor_revision_validates_all_records():
    """Productive IGH + every corruption stage + ``receptor_revision(prob=1.0)``.

    Every record's V segment is rewritten by the receptor-revision
    pass, so every record exercises the post-recombine V-replacement
    code path through `SegmentReplaced` event emission, the
    AllStructural-equivalent live-call refresh, AIRR projection of
    the new `receptor_revision_applied` / `original_v_call` fields,
    and the validator's two corresponding consistency checks. Both
    two-layer gates must run clean.
    """
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .receptor_revision(prob=1.0)
        .productive_only()
        .mutate(rate=0.03)
        .polymerase_indels(count=2)
    )
    refdata = exp.refdata
    result = exp.run_records(n=50, seed=4242)

    # Gate 1: AIRR validator. Every postcondition, including the
    # new `ReceptorRevisionAppliedMismatch` and `OriginalVCallMismatch`
    # checks from Slice E.
    report = result.validate_records(refdata)
    assert report, (
        f"revised-V productive IGH validation failed: "
        f"summary={report.summary()}; first failure={report.failures[0] if report.failures else None}"
    )

    # Gate 2: live-call cache parity. After a `SegmentReplaced(V)`
    # event the refresh hook runs an AllStructural-equivalent
    # re-walk; the cached `SegmentLiveCall` must agree with a fresh
    # from-scratch recompute on every outcome.
    assert result.outcomes is not None
    for i, outcome in enumerate(result.outcomes):
        for p in outcome.check_live_call_cache_parity(refdata):
            assert p["tie_set_matches"], (
                f"revised-V record {i} segment {p['segment']}: "
                f"cache parity failed after receptor revision — "
                f"cached={p['cached_tie_set']} fresh={p['fresh_tie_set']}"
            )

    # Every record's receptor_revision_applied must be True (prob=1
    # commits revision on every seed). Smoke check that the AIRR
    # field plumbing fires for the realistic stack.
    assert all(
        rec["receptor_revision_applied"] is True for rec in result.records
    ), "prob=1.0 must commit a revision on every record"
    # And every record carries a non-empty original_v_call when
    # applied — the trace-sourced pre-revision V name.
    assert all(rec["original_v_call"] for rec in result.records), (
        "applied=True records must carry a non-empty original_v_call"
    )


def test_receptor_revision_trace_replay_round_trip_passes_validator():
    """Replay determinism with receptor revision in the chain.

    For each fresh outcome: build a TraceFile, rerun it via
    `rerun_from_trace_file`, build the AIRR record from the replayed
    outcome, and assert it (a) carries the same
    `receptor_revision_applied` + `original_v_call` values and (b)
    passes the validator with no `ReceptorRevisionAppliedMismatch`
    or `OriginalVCallMismatch`.

    Uses ``prob=0.5`` over multiple seeds so the sweep exercises
    both branches of the revision decision; the replay must
    reproduce whichever branch fired in the fresh run.
    """
    from GenAIRR._airr_record import outcome_to_airr_record

    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .receptor_revision(prob=0.5)
        .productive_only()
        .mutate(rate=0.01)
    )
    refdata = exp.refdata
    compiled = exp.compile()
    seen_true = False
    seen_false = False
    for seed in range(8):
        fresh = compiled.simulator.run(seed=seed)
        tf = compiled.simulator.trace_file_from(fresh, seed=seed)
        replayed = compiled.simulator.rerun_from_trace_file(tf)

        fresh_rec = outcome_to_airr_record(
            fresh, refdata, sequence_id=f"fresh-{seed}"
        )
        replayed_rec = outcome_to_airr_record(
            replayed, refdata, sequence_id=f"replay-{seed}"
        )

        # Both provenance fields must round-trip.
        assert (
            fresh_rec["receptor_revision_applied"]
            == replayed_rec["receptor_revision_applied"]
        ), (
            f"seed {seed}: receptor_revision_applied desynced through replay "
            f"({fresh_rec['receptor_revision_applied']} vs "
            f"{replayed_rec['receptor_revision_applied']})"
        )
        assert (
            fresh_rec["original_v_call"] == replayed_rec["original_v_call"]
        ), (
            f"seed {seed}: original_v_call desynced through replay "
            f"({fresh_rec['original_v_call']!r} vs "
            f"{replayed_rec['original_v_call']!r})"
        )

        if fresh_rec["receptor_revision_applied"]:
            seen_true = True
        else:
            seen_false = True

        # Full AIRR record's sequence bytes must round-trip too.
        assert fresh_rec["sequence"] == replayed_rec["sequence"], (
            f"seed {seed}: assembled sequence diverged under replay"
        )

        # And the replayed record must pass the two new validator
        # checks with no mismatches.
        replayed_issues = replayed.validate_record(
            refdata, sequence_id=f"replay-{seed}"
        )
        kinds = {issue["kind"] for issue in replayed_issues}
        assert "ReceptorRevisionAppliedMismatch" not in kinds, (
            f"seed {seed}: replayed record tripped "
            f"ReceptorRevisionAppliedMismatch"
        )
        assert "OriginalVCallMismatch" not in kinds, (
            f"seed {seed}: replayed record tripped OriginalVCallMismatch"
        )

    # Both branches must fire across the 8 seeds at prob=0.5; if
    # they don't, the test isn't exercising both halves of the
    # revision decision and a regression in one branch would slip
    # through unnoticed.
    assert seen_true, "prob=0.5 over 8 seeds: no True branch ever fired"
    assert seen_false, "prob=0.5 over 8 seeds: no False branch ever fired"


# ──────────────────────────────────────────────────────────────────
# Paired-end / read layout (audit-first biology mechanism #3 — by
# arc count; biologically it's a sequencing-stage projection, not
# a biology mechanism). Same shape as the receptor-revision tests:
# validator + cache parity under a realistic productive IGH stack
# composed with every other audit mechanism, plus a replay round-
# trip that pins all eight paired-end fields. The distribution
# invariant (`insert_size=(low, high)`, ±5σ) lives in
# `test_distribution_invariants.py` alongside the equivalent
# Bernoulli draws.
# ──────────────────────────────────────────────────────────────────


def test_productive_igh_full_stack_with_paired_end_validates_all_records():
    """Productive IGH composed with **every** audit mechanism —
    receptor revision + D inversion + mutation + end-loss (both
    sides) + random strand orientation + paired-end.

    This is the high-value composition check: paired-end must
    remain a projection layer after every IR-mutating + every
    observation-stage mechanism. Both two-layer gates run clean,
    and every record carries the eight populated paired-end fields.

    The D-tie oracle is still known-incomplete under inversion
    (see `_strip_inverted_d_tie_set_issues`); the remaining issue
    set still includes every paired-end check
    (`PairedEndFieldWithoutLayout`, `ReadWindowOutOfBounds`,
    `ReadSequenceMismatch`, `ReadInsertSizeMismatch`,
    `ReadLayoutMismatch`).
    """
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .receptor_revision(prob=0.5)
        .invert_d(prob=0.5)
        .productive_only()
        .mutate(rate=0.03)
        .end_loss_5prime(length=[(2, 1.0)])
        .end_loss_3prime(length=[(2, 1.0)])
        .random_strand_orientation(prob=0.5)
        .paired_end(r1_length=80, insert_size=200)
    )
    refdata = exp.refdata
    result = exp.run_records(n=50, seed=4242)

    # Gate 1: AIRR validator. No exceptions stripped — the
    # orientation-aware walker + validator oracle converge on the
    # inverted-D tie set directly.
    report = result.validate_records(refdata)
    assert report, (
        f"paired-end full-stack validation failed: "
        f"summary={report.summary()}; first failure="
        f"{report.failures[0] if report.failures else None}"
    )

    # Gate 2: live-call cache parity on every outcome. Paired-end
    # is projection-only and must not invalidate the live-call
    # layer — pin that cache parity stays clean even under the
    # full stack of preceding mechanisms.
    assert result.outcomes is not None
    for i, outcome in enumerate(result.outcomes):
        for p in outcome.check_live_call_cache_parity(refdata):
            assert p["tie_set_matches"], (
                f"paired-end record {i} segment {p['segment']}: "
                f"cache parity failed — cached={p['cached_tie_set']} "
                f"fresh={p['fresh_tie_set']}"
            )

    # Every record carries the eight paired-end fields populated.
    # Smoke check that the AIRR field plumbing fires for the
    # realistic stack.
    assert all(
        rec["read_layout"] == "paired_end" for rec in result.records
    ), "paired-end must populate every record under the full stack"
    assert all(
        len(rec["r1_sequence"]) == 80 for rec in result.records
    ), "r1_sequence length must match the requested r1_length"
    assert all(
        rec["insert_size"] == 200 for rec in result.records
    ), "insert_size must match the requested value"


def test_paired_end_trace_replay_round_trip_passes_validator():
    """Replay determinism with paired-end in the chain over
    multiple seeds with **variable** insert sizes (a
    `(low, high)` uniform-int distribution).

    For each fresh outcome: build a TraceFile, rerun it via
    ``rerun_from_trace_file``, build the AIRR record from the
    replayed outcome, and assert all eight paired-end fields +
    the full ``sequence`` round-trip bit-for-bit. The variable
    insert size pins that the per-seed sampled values flow
    through the trace into the AIRR projection identically on
    replay.
    """
    from GenAIRR._airr_record import outcome_to_airr_record

    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .productive_only()
        .mutate(rate=0.01)
        .paired_end(r1_length=80, insert_size=(150, 280))
    )
    refdata = exp.refdata
    compiled = exp.compile()
    seen_insert_sizes = set()
    for seed in range(8):
        fresh = compiled.simulator.run(seed=seed)
        tf = compiled.simulator.trace_file_from(fresh, seed=seed)
        replayed = compiled.simulator.rerun_from_trace_file(tf)

        fresh_rec = outcome_to_airr_record(
            fresh, refdata, sequence_id=f"fresh-{seed}"
        )
        replayed_rec = outcome_to_airr_record(
            replayed, refdata, sequence_id=f"replay-{seed}"
        )

        for field in (
            "read_layout",
            "r1_sequence",
            "r2_sequence",
            "r1_start",
            "r1_end",
            "r2_start",
            "r2_end",
            "insert_size",
            "sequence",
        ):
            assert fresh_rec[field] == replayed_rec[field], (
                f"seed {seed}: paired-end field {field!r} desynced under "
                f"replay ({fresh_rec[field]!r} vs {replayed_rec[field]!r})"
            )
        seen_insert_sizes.add(fresh_rec["insert_size"])

        # And the replayed record must pass the validator unchanged
        # — none of the five paired-end issue variants surface.
        replayed_issues = replayed.validate_record(
            refdata, sequence_id=f"replay-{seed}"
        )
        kinds = {issue["kind"] for issue in replayed_issues}
        for forbidden in (
            "PairedEndFieldWithoutLayout",
            "ReadWindowOutOfBounds",
            "ReadSequenceMismatch",
            "ReadInsertSizeMismatch",
            "ReadLayoutMismatch",
        ):
            assert forbidden not in kinds, (
                f"seed {seed}: replayed record tripped {forbidden}"
            )

    # Sanity: across 8 seeds we must observe at least two distinct
    # insert sizes — otherwise the test isn't exercising the
    # variable-insert path and a regression that froze the
    # distribution would slip through.
    assert len(seen_insert_sizes) >= 2, (
        f"variable-insert test only saw one distinct insert size "
        f"across 8 seeds: {seen_insert_sizes}"
    )


# ──────────────────────────────────────────────────────────────────
# Targeted SHM (per-segment rate scalars) — release-tier consolidation
#
# Closes the per-segment SHM rates slice the same way D inversion /
# receptor revision / paired-end closed: a productive IGH full-stack
# run with realistic kwargs + replay round-trip + the zero-rate
# exclusion invariant. See ``docs/shm_segment_rate_design.md`` (audit)
# and ``tests/test_segment_rates_implementation.py`` (slice spec
# tests).
# ──────────────────────────────────────────────────────────────────


def test_productive_igh_full_stack_with_segment_rates_validates_all_records():
    """Productive IGH + every corruption stage + ``invert_d`` +
    ``receptor_revision`` + ``paired_end`` with a realistic
    non-default ``segment_rates`` vector.

    Exercises the entire targeted-SHM code path under contracts
    that are most likely to interact with segment-rate filtering
    (productive_only + heavy SHM rate). Asserts both the AIRR
    validator and the live-call cache parity layer run clean —
    same two-gate posture the d_inversion / receptor_revision /
    paired_end stacks already use."""
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .invert_d(prob=0.3)
        .receptor_revision(prob=0.3)
        .productive_only()
        .mutate(
            model="s5f",
            rate=0.03,
            segment_rates={"V": 1.0, "D": 0.2, "J": 0.5, "NP": 0.0},
        )
        .pcr_amplify(rate=1e-4)
        .polymerase_indels(count=2)
        .primer_trim_5prime(length=(0, 3))
        .primer_trim_3prime(length=(0, 3))
        .paired_end(r1_length=150, insert_size=300)
    )
    refdata = exp.refdata
    result = exp.run_records(n=50, seed=4242)

    # Gate 1: AIRR validator. Every postcondition runs clean.
    report = result.validate_records(refdata)
    assert report, (
        f"targeted-SHM productive IGH validation failed: "
        f"summary={report.summary()}; first failure="
        f"{report.failures[0] if report.failures else None}"
    )

    # Gate 2: live-call cache parity on every outcome.
    assert result.outcomes is not None
    for i, outcome in enumerate(result.outcomes):
        for p in outcome.check_live_call_cache_parity(refdata):
            assert p["tie_set_matches"], (
                f"targeted-SHM record {i} segment {p['segment']}: "
                f"cache parity failed — cached={p['cached_tie_set']} "
                f"fresh={p['fresh_tie_set']}"
            )

    # Paired-end fields populated on every record (smoke check that
    # the full stack ran through to the last pass). The downstream
    # corruption passes (polymerase indels, end-loss) can flip
    # ``productive`` to False on some records even when
    # ``productive_only()`` constrained the sampling — that's the
    # documented behaviour shared by the d_inversion and
    # receptor_revision release tests. The two-gate validator is
    # the real correctness check.
    for rec in result.records:
        assert rec["r1_sequence"], "paired-end r1_sequence empty"
        assert rec["r2_sequence"], "paired-end r2_sequence empty"


def test_segment_rates_trace_replay_round_trip_passes_validator():
    """Replay determinism with non-default ``segment_rates`` in the
    chain.

    For each fresh outcome: build a TraceFile, rerun via
    ``rerun_from_trace_file``, project to AIRR and assert the
    replayed record (a) reproduces the assembled sequence and
    ``n_mutations`` exactly, (b) reproduces the live calls and
    junction fields, and (c) passes the per-record validator.

    Uses a mixed-rate vector (V/D/J non-zero, NP zero) so the
    replay path's zero-rate validation actually has work to do —
    a recorded site that fell in V must still validate under the
    same rate vector at replay time.
    """
    from GenAIRR._airr_record import outcome_to_airr_record

    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .productive_only()
        .mutate(
            model="s5f",
            rate=0.03,
            segment_rates={"V": 1.0, "D": 0.5, "J": 0.5, "NP": 0.0},
        )
    )
    refdata = exp.refdata
    compiled = exp.compile()
    seen_any_mutations = False
    for seed in range(8):
        fresh = compiled.simulator.run(seed=seed)
        tf = compiled.simulator.trace_file_from(fresh, seed=seed)
        replayed = compiled.simulator.rerun_from_trace_file(tf)

        fresh_rec = outcome_to_airr_record(
            fresh, refdata, sequence_id=f"fresh-{seed}"
        )
        replayed_rec = outcome_to_airr_record(
            replayed, refdata, sequence_id=f"replay-{seed}"
        )

        # Sequence + n_mutations round-trip.
        assert fresh_rec["sequence"] == replayed_rec["sequence"], (
            f"seed {seed}: assembled sequence diverged under replay"
        )
        assert fresh_rec["n_mutations"] == replayed_rec["n_mutations"], (
            f"seed {seed}: n_mutations desynced "
            f"({fresh_rec['n_mutations']} vs {replayed_rec['n_mutations']})"
        )
        # Live calls + junction also round-trip — the segment-rate
        # replay path doesn't disturb projection.
        for field in ("v_call", "d_call", "j_call", "junction", "junction_aa"):
            assert fresh_rec[field] == replayed_rec[field], (
                f"seed {seed} field {field!r}: replay diverged "
                f"({fresh_rec[field]!r} vs {replayed_rec[field]!r})"
            )

        if fresh_rec["n_mutations"] > 0:
            seen_any_mutations = True

        # Replayed record passes the per-record validator.
        replayed_issues = replayed.validate_record(
            refdata, sequence_id=f"replay-{seed}"
        )
        assert not replayed_issues, (
            f"seed {seed}: replayed record carries unexpected issues "
            f"{[i.get('kind') for i in replayed_issues]}"
        )

    # The sweep must actually exercise SHM (otherwise the replay
    # path's segment-rate validation never runs).
    assert seen_any_mutations, (
        "no SHM mutations observed across 8 seeds — replay sweep "
        "isn't exercising the segment-rate validation path."
    )


def test_segment_rates_zero_rate_exclusion_invariant_full_stack():
    """**Zero-rate exclusion invariant** — the load-bearing
    distribution check for the segment-rates slice. Across N
    records of a productive-IGH full stack with
    ``segment_rates={"D": 0, "J": 0, "NP": 0}``, every mutated
    site must lie within the V region. Any single mutated base
    outside V is a support-leak bug; the test fails immediately
    on the first offender so the diagnostic points at the
    specific record + segment that leaked.

    Uses the AIRR record's V/D/J/J coordinate fields to assign
    each mutated position to a segment. Comparing the
    no-mutation baseline (same seed, same recombination) to the
    mutated batch lets us identify which positions changed; the
    segment classification flags zero-rate leaks."""
    n = 50
    seed = 8181
    base_exp = (
        ga.Experiment.on("human_igh").recombine().productive_only()
    )
    targeted_exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .productive_only()
        .mutate(
            model="s5f",
            count=25,
            segment_rates={"V": 1.0, "D": 0.0, "J": 0.0, "NP": 0.0},
        )
    )
    base = base_exp.run_records(n=n, seed=seed)
    targeted = targeted_exp.run_records(n=n, seed=seed)

    saw_any_mutation = False
    for i, (base_rec, mut_rec) in enumerate(zip(base, targeted)):
        b_seq = base_rec["sequence"].upper()
        m_seq = mut_rec["sequence"].upper()
        if b_seq == m_seq:
            # Same seed gave zero realised mutations on this record
            # (rare but possible with constraint filtering); nothing
            # to classify.
            continue
        # Per-record segment boundaries.
        v_end = mut_rec["v_sequence_end"]
        d_start = mut_rec["d_sequence_start"]
        d_end = mut_rec["d_sequence_end"]
        j_start = mut_rec["j_sequence_start"]
        # Sequences must be the same length — segment_rates doesn't
        # change pool length under SHM (substitutions only).
        assert len(b_seq) == len(m_seq), (
            f"record {i}: sequence length changed under SHM "
            f"(base={len(b_seq)}, mutated={len(m_seq)}); SHM should "
            "only substitute, not insert / delete."
        )
        for pos, (bb, mb) in enumerate(zip(b_seq, m_seq)):
            if bb == mb:
                continue
            saw_any_mutation = True
            if pos < v_end:
                segment_label = "V"
            elif pos < d_start:
                segment_label = "NP1"
            elif pos < d_end:
                segment_label = "D"
            elif pos < j_start:
                segment_label = "NP2"
            else:
                segment_label = "J"
            assert segment_label == "V", (
                f"record {i} pos {pos} (segment {segment_label}): "
                "mutation landed outside V despite "
                "segment_rates={V:1, D:0, J:0, NP:0}. "
                "Zero-rate exclusion broke."
            )

    # The invariant only has bite when the sweep actually exercises
    # SHM on at least one record.
    assert saw_any_mutation, (
        "targeted SHM sweep produced zero mutations across 50 records; "
        "test isn't exercising the path."
    )


# ──────────────────────────────────────────────────────────────────
# Mutation provenance counters — release-tier consolidation
#
# Closes the per-segment SHM counter slice with the same closure
# standard as targeted SHM: full-stack IGH + non-default
# segment_rates + sum-invariant cross-check + corruption isolation.
# See ``docs/mutation_provenance_audit.md`` for the architecture
# contract and ``tests/test_per_segment_mutation_counters.py`` for
# the slice spec tests.
# ──────────────────────────────────────────────────────────────────


def test_per_segment_mutation_counters_full_stack_validates_and_partitions():
    """Productive IGH full stack with non-default ``segment_rates``
    + heavy corruption. Three load-bearing assertions:

    1. ``validate_records(refdata)`` runs clean — the validator's
       re-derived per-segment counts agree with the engine's by
       construction.
    2. **Sum invariant** ``n_v + n_d + n_j + n_np == n_mutations``
       holds for every record. The audit's headline claim.
    3. **Corruption isolation**: ``n_pcr_errors`` /
       ``n_quality_errors`` are populated (heavy corruption ran)
       while the four per-segment SHM counters reflect biological
       SHM only — the pass-name filter excludes PCR / quality
       ``BaseChanged`` events from the per-segment buckets.
    """
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .invert_d(prob=0.3)
        .receptor_revision(prob=0.3)
        .productive_only()
        .mutate(
            model="s5f",
            count=20,
            segment_rates={"V": 1.0, "D": 0.4, "J": 0.6, "NP": 0.2},
        )
        .pcr_amplify(count=10)
        .sequencing_errors(count=8)
        .polymerase_indels(count=2)
        .primer_trim_5prime(length=(0, 3))
        .primer_trim_3prime(length=(0, 3))
    )
    refdata = exp.refdata
    result = exp.run_records(n=50, seed=4242)

    # Gate 1: AIRR validator. Every record's projected counters
    # agree with the engine's re-derivation; the new five issue
    # kinds (`N{V,D,J,Np}MutationsMismatch` + sum-invariant
    # cross-check) all fire zero.
    report = result.validate_records(refdata)
    assert report, (
        f"per-segment counter full stack tripped validator: "
        f"summary={report.summary()}; first failure="
        f"{report.failures[0] if report.failures else None}"
    )

    # Gate 2: live-call cache parity (the standard release
    # posture; per-segment counters don't touch live-call state
    # but the cross-layer sanity check stays).
    assert result.outcomes is not None
    for i, outcome in enumerate(result.outcomes):
        for p in outcome.check_live_call_cache_parity(refdata):
            assert p["tie_set_matches"], (
                f"counter full-stack record {i} segment {p['segment']}: "
                f"cache parity failed — cached={p['cached_tie_set']} "
                f"fresh={p['fresh_tie_set']}"
            )

    # Gate 3: sum invariant on every record. The audit's
    # load-bearing pin promoted to release-tier; any tampering /
    # leakage / off-by-one in the aggregation loop surfaces here
    # before reaching downstream consumers.
    for r in result:
        total = (
            r["n_v_mutations"]
            + r["n_d_mutations"]
            + r["n_j_mutations"]
            + r["n_np_mutations"]
        )
        assert total == r["n_mutations"], (
            f"sum invariant violated: "
            f"V={r['n_v_mutations']} + D={r['n_d_mutations']} + "
            f"J={r['n_j_mutations']} + NP={r['n_np_mutations']} = "
            f"{total} != n_mutations={r['n_mutations']}"
        )

    # Gate 4: corruption isolation. Across the batch we must
    # observe non-zero PCR + quality counters (corruption is
    # configured at non-trivial counts; some records must accrue
    # at least one of each) AND the per-segment SHM counters must
    # reflect biological SHM only. Verify by reconciling the
    # per-segment sum against the segment-rate-filtered
    # ``n_mutations`` — they were verified equal record-by-record
    # above; here we confirm corruption ran independently.
    saw_pcr = any(r["n_pcr_errors"] > 0 for r in result)
    saw_quality = any(r["n_quality_errors"] > 0 for r in result)
    assert saw_pcr, (
        "no PCR errors observed across 50 records despite "
        "pcr_amplify(count=10) — heavy corruption isn't running."
    )
    assert saw_quality, (
        "no quality errors observed despite sequencing_errors(count=8)."
    )

    # Sanity: with NP at 0.2 and other segments higher, every
    # bucket should accrue at least one mutation across the batch.
    # This pins the segment-rate ↔ counter interaction at the
    # release-tier (a future bug that silently routed every NP
    # mutation to another bucket would surface here).
    batch_v = sum(r["n_v_mutations"] for r in result)
    batch_d = sum(r["n_d_mutations"] for r in result)
    batch_j = sum(r["n_j_mutations"] for r in result)
    batch_np = sum(r["n_np_mutations"] for r in result)
    batch_global = sum(r["n_mutations"] for r in result)
    assert batch_v + batch_d + batch_j + batch_np == batch_global, (
        f"batch-level sum invariant violated: "
        f"V={batch_v} + D={batch_d} + J={batch_j} + NP={batch_np} = "
        f"{batch_v + batch_d + batch_j + batch_np} != "
        f"global={batch_global}"
    )
    assert batch_v > 0, "V received zero SHM hits across the batch"
    assert batch_j > 0, "J received zero SHM hits across the batch"


def test_v_subregion_annotation_release_sanity():
    """Release-tier pin for the **V-Subregion Cartridge Annotation
    Surface** (Slice 1).

    The slice promised: V-region substructure annotations are a
    first-class cartridge property — derived from IMGT-gapped
    sequences, inspectable through the manifest, folded into
    ``refdata_content_hash`` — but they DO NOT yet perturb
    simulation behaviour. The Slice 2 sampling kwarg
    (``v_subregion_rates``) and Slice 3 per-region counters
    (``n_cdr1_mutations`` etc.) are still deferred.

    This test pins four release-level invariants together:

      1. Bundled human IGH OGRDB reports 100% V-subregion
         coverage in the cartridge manifest.
      2. A normal simulation (productive_only + S5F SHM)
         validates clean — adding the annotation surface did
         NOT regress any AIRR-record contract.
      3. Editing a single subregion interval flips
         ``refdata_content_hash`` — the slice's
         ``in_content_hash=True`` claim holds at the release
         layer, not just in unit tests.
      4. The deferred surfaces stay absent: no
         ``v_subregion_rates`` kwarg on
         ``Experiment.mutate``; no ``n_cdr*_mutations`` /
         ``n_fwr*_mutations`` fields on AIRR records.

    If any of these fail, either Slice 2 / Slice 3 has
    silently landed (good — update this pin) or a
    cartridge-side regression broke the annotation plumbing
    (bad — diagnose and fix).
    """
    import copy
    import inspect

    cfg = ga.HUMAN_IGH_OGRDB

    # ── (1) Manifest reports 100% V-subregion coverage ──────────
    manifest = cfg.cartridge_manifest()
    sup = manifest["models"]["shm"]["v_subregion_support"]
    assert sup["available"] is True, (
        "v_subregion_support reports available=False on a bundled "
        "cartridge with IMGT gapped_seq coverage; bridge-time "
        "derivation regressed"
    )
    assert sup["annotated_v_count"] == sup["total_v_count"], (
        f"V-subregion coverage incomplete: "
        f"{sup['annotated_v_count']}/{sup['total_v_count']} V alleles "
        "carry subregions on a bundled cartridge with full gapped_seq "
        "coverage; bridge-time derivation regressed"
    )
    assert sup["total_v_count"] > 0, "bundled HUMAN_IGH has zero V alleles"
    assert sup["labels"] == ["FWR1", "CDR1", "FWR2", "CDR2", "FWR3"]
    assert sup["derivation"] == "bridge_imgt_gapped_seq"
    assert sup["in_content_hash"] is True

    # ── (2) Normal simulation validates clean ──────────────────
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .productive_only()
        .mutate(model="s5f", rate=0.03)
    )
    refdata = exp.refdata
    result = exp.run_records(n=50, seed=4242)
    report = result.validate_records(refdata)
    assert report, (
        f"V-subregion annotation surface regressed baseline simulation: "
        f"{len(report.failures)}/{report.count} records failed validation; "
        f"summary={report.summary()}; "
        f"first failure={report.failures[0] if report.failures else None}"
    )

    # ── (3) Editing a subregion interval flips content_hash ────
    from GenAIRR._refdata_resolver import dataconfig_to_refdata

    base_hash = dataconfig_to_refdata(cfg).content_hash()

    edited = copy.deepcopy(cfg)
    target = next(iter(edited.v_alleles.values()))[0]
    seq_len = len(target.ungapped_seq)
    # Choose an interval set that is well-formed but unlikely to
    # collide with the IMGT-derived one for any allele.
    target.subregions = {
        "FWR1": (0, min(60, seq_len // 4)),
        "CDR1": (min(60, seq_len // 4), min(78, seq_len // 3)),
        "FWR2": (min(78, seq_len // 3), min(120, seq_len // 2)),
    }
    edited_hash = dataconfig_to_refdata(edited).content_hash()
    assert edited_hash != base_hash, (
        "editing a subregion interval did NOT flip refdata_content_hash; "
        "the slice's in_content_hash=True claim regressed"
    )

    # ── (4) v_subregion_rates kwarg is now PRESENT (Slice B). ──
    sig = inspect.signature(ga.Experiment.mutate)
    assert "v_subregion_rates" in sig.parameters, (
        "Experiment.mutate is missing the v_subregion_rates kwarg; "
        "Slice B surface regressed"
    )
    # And the alternate / legacy names stay absent — only the
    # canonical spelling is the user-facing surface.
    for forbidden in ("subregion_rates", "v_region_rates"):
        assert forbidden not in sig.parameters, (
            f"Experiment.mutate now accepts {forbidden!r}; only "
            "``v_subregion_rates`` is the canonical surface."
        )

    # ── (5) V-subregion mutation counters are now PRESENT. ────
    # The counters slice landed after Slice B; six AIRR fields
    # partition n_v_mutations cleanly. Two-bucket aggregates stay
    # deliberately absent per the audit's §1 recommendation.
    rec = result.records[0]
    for required in (
        "n_cdr1_mutations",
        "n_cdr2_mutations",
        "n_fwr1_mutations",
        "n_fwr2_mutations",
        "n_fwr3_mutations",
        "n_v_unannotated_mutations",
    ):
        assert required in rec, (
            f"AIRR record missing {required!r}; V-subregion counters "
            "slice surface regressed"
        )
    for forbidden in ("n_cdr_mutations", "n_fwr_mutations"):
        assert forbidden not in rec, (
            f"AIRR record now carries two-bucket aggregate {forbidden!r}; "
            "counters audit §1 recommended five-label fields only"
        )


# ──────────────────────────────────────────────────────────────────
# V-subregion SHM targeting — release-tier consolidation (Slice B)
#
# Closes the V-subregion SHM rate slice with the same closure
# standard as targeted SHM: full-stack IGH + non-default
# (segment_rates × v_subregion_rates) + AIRR validator + cache
# parity + per-segment sum invariant + replay round-trip.
# See ``docs/v_subregion_shm_rate_design.md`` for the architecture
# contract and ``tests/test_v_subregion_rates_implementation.py``
# for the per-spec implementation tests.
# ──────────────────────────────────────────────────────────────────


def test_productive_igh_full_stack_with_v_subregion_rates_validates_all_records():
    """Productive IGH + every corruption stage + ``invert_d`` +
    ``receptor_revision`` + ``paired_end`` with a realistic
    non-default ``(segment_rates, v_subregion_rates)`` pair.

    Exercises the unified ``combined_site_factor`` weighting path
    (segment × V-subregion) under contracts that interact most
    aggressively with the rate filters (productive_only + heavy
    SHM rate + CDR-targeted subregion vector). Five gates:

    1. AIRR validator passes on every record.
    2. Live-call cache parity holds on every outcome (the
       refresh-path bug class the segment-rates release test
       guards against — re-checked here under the combined
       factor).
    3. Per-segment SHM sum invariant
       ``n_v + n_d + n_j + n_np == n_mutations`` holds — Slice B
       doesn't introduce per-region counters, but the existing
       per-segment counters must still partition cleanly.
    4. V receives mutations under the targeted CDR/FWR rates —
       the slice isn't silently dropping all V sites.
    5. Paired-end fields populate on every record — the slice
       didn't break the full plan execution."""
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .invert_d(prob=0.3)
        .receptor_revision(prob=0.3)
        .productive_only()
        .mutate(
            model="s5f",
            rate=0.04,
            segment_rates={"V": 1.0, "D": 0.3, "J": 0.5, "NP": 0.2},
            v_subregion_rates={"CDR": 3.0, "FWR": 0.5},
        )
        .pcr_amplify(rate=1e-4)
        .polymerase_indels(count=2)
        .primer_trim_5prime(length=(0, 3))
        .primer_trim_3prime(length=(0, 3))
        .paired_end(r1_length=150, insert_size=300)
    )
    refdata = exp.refdata
    result = exp.run_records(n=50, seed=4242)

    # Gate 1: AIRR validator.
    report = result.validate_records(refdata)
    assert report, (
        f"V-subregion-targeted SHM productive IGH validation failed: "
        f"summary={report.summary()}; first failure="
        f"{report.failures[0] if report.failures else None}"
    )

    # Gate 2: live-call cache parity.
    assert result.outcomes is not None
    for i, outcome in enumerate(result.outcomes):
        for p in outcome.check_live_call_cache_parity(refdata):
            assert p["tie_set_matches"], (
                f"V-subregion-targeted record {i} segment {p['segment']}: "
                f"cache parity failed — cached={p['cached_tie_set']} "
                f"fresh={p['fresh_tie_set']}"
            )

    # Gate 3: per-segment SHM sum invariant per record. Slice B
    # composes with the per-segment counters slice — every
    # mutation must still be attributed to exactly one segment.
    for r in result:
        total = (
            r["n_v_mutations"]
            + r["n_d_mutations"]
            + r["n_j_mutations"]
            + r["n_np_mutations"]
        )
        assert total == r["n_mutations"], (
            f"per-segment sum invariant violated under v_subregion_rates: "
            f"V={r['n_v_mutations']} + D={r['n_d_mutations']} + "
            f"J={r['n_j_mutations']} + NP={r['n_np_mutations']} = "
            f"{total} != n_mutations={r['n_mutations']}"
        )

    # Gate 4: V hits accrue across the batch — the CDR-targeted /
    # FWR-suppressed vector hasn't dropped the entire V segment.
    batch_v = sum(r["n_v_mutations"] for r in result)
    assert batch_v > 0, (
        "V received zero SHM hits across 50 records under "
        "v_subregion_rates={CDR: 3.0, FWR: 0.5} — the targeted "
        "vector is silently dropping all V sites"
    )

    # Gate 5: paired-end fields populated on every record.
    for rec in result.records:
        assert rec["r1_sequence"], "paired-end r1_sequence empty"
        assert rec["r2_sequence"], "paired-end r2_sequence empty"


def test_v_subregion_rates_trace_replay_round_trip_passes_validator():
    """Replay determinism with non-default ``v_subregion_rates`` —
    a recorded V site that fell in a non-zero subregion must
    still validate under the same rate vector at replay time, and
    a sequence-level + counter-level round-trip must hold.

    Mirrors the ``segment_rates`` replay round-trip discipline
    (eight seeds, AIRR-projected fields compared, per-record
    validator clean) and additionally pins that the V-CDR rate
    factor 0 doesn't desync between fresh sampling and replay.
    """
    from GenAIRR._airr_record import outcome_to_airr_record

    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .productive_only()
        .mutate(
            model="s5f",
            rate=0.03,
            v_subregion_rates={"CDR": 2.5, "FWR": 0.5},
        )
    )
    refdata = exp.refdata
    compiled = exp.compile()
    seen_any_mutations = False
    for seed in range(8):
        fresh = compiled.simulator.run(seed=seed)
        tf = compiled.simulator.trace_file_from(fresh, seed=seed)
        replayed = compiled.simulator.rerun_from_trace_file(tf)

        fresh_rec = outcome_to_airr_record(
            fresh, refdata, sequence_id=f"fresh-{seed}"
        )
        replayed_rec = outcome_to_airr_record(
            replayed, refdata, sequence_id=f"replay-{seed}"
        )

        # Sequence + n_mutations round-trip.
        assert fresh_rec["sequence"] == replayed_rec["sequence"], (
            f"seed {seed}: assembled sequence diverged under "
            "v_subregion_rates replay"
        )
        assert fresh_rec["n_mutations"] == replayed_rec["n_mutations"], (
            f"seed {seed}: n_mutations desynced "
            f"({fresh_rec['n_mutations']} vs {replayed_rec['n_mutations']})"
        )
        # Live calls + junction also round-trip.
        for field in ("v_call", "d_call", "j_call", "junction", "junction_aa"):
            assert fresh_rec[field] == replayed_rec[field], (
                f"seed {seed} field {field!r}: replay diverged "
                f"({fresh_rec[field]!r} vs {replayed_rec[field]!r})"
            )
        # Per-segment counters round-trip.
        for field in (
            "n_v_mutations",
            "n_d_mutations",
            "n_j_mutations",
            "n_np_mutations",
            # Per-V-subregion partition counters (V-subregion mutation
            # counters slice). All six round-trip under replay because
            # they're derived from the deterministic event ledger.
            "n_fwr1_mutations",
            "n_cdr1_mutations",
            "n_fwr2_mutations",
            "n_cdr2_mutations",
            "n_fwr3_mutations",
            "n_v_unannotated_mutations",
        ):
            assert fresh_rec[field] == replayed_rec[field], (
                f"seed {seed} counter {field!r} desynced under replay "
                f"({fresh_rec[field]} vs {replayed_rec[field]})"
            )

        if fresh_rec["n_mutations"] > 0:
            seen_any_mutations = True

        replayed_issues = replayed.validate_record(
            refdata, sequence_id=f"replay-{seed}"
        )
        assert not replayed_issues, (
            f"seed {seed}: replayed record carries unexpected issues "
            f"{[i.get('kind') for i in replayed_issues]}"
        )

    assert seen_any_mutations, (
        "no SHM mutations observed across 8 seeds — replay sweep "
        "isn't exercising the v_subregion_rates path."
    )


def test_v_subregion_rates_zero_rate_exclusion_invariant_full_stack():
    """**V-subregion zero-rate exclusion invariant** — the
    distribution check that pins the slice's headline claim:
    a zero-rate V subregion drops every site in that subregion
    out of proposal support BEFORE contract admissibility.

    Two complementary configurations:

      (a) ``v_subregion_rates={"CDR": 0.0, "FWR": 1.0}`` — no
          CDR1 / CDR2 mutations should land on annotated V
          alleles.
      (b) ``v_subregion_rates={"FWR": 0.0, "CDR": 1.0}`` — no
          FWR1 / FWR2 / FWR3 mutations on annotated V alleles.

    For each configuration, classify every observed V-segment
    mutation by IMGT subregion using the cartridge's subregion
    intervals + the per-record assigned V allele + the
    germline-start offset. The forbidden bucket must be empty
    across N=50 records; the *complementary* bucket must be
    non-empty (otherwise the test is vacuously passing).
    """
    from GenAIRR._refdata_resolver import dataconfig_to_refdata
    import copy

    cfg = ga.HUMAN_IGH_OGRDB

    def _v_classify(mut_rec, base_rec, refdata):
        """Return ``(cdr_hits, fwr_hits)`` over the V region of
        one record, using subregion intervals on the assigned V
        allele to bucket each differing position. Same shape as
        the baseline-diff classifier in
        ``test_v_subregion_rates_implementation.py``."""
        base_seq = base_rec["sequence"].upper()
        mut_seq = mut_rec["sequence"].upper()
        if len(base_seq) != len(mut_seq):
            return 0, 0
        v_end = mut_rec.get("v_sequence_end")
        if v_end is None:
            return 0, 0
        germ_start = mut_rec.get("v_germline_start")
        if germ_start is None:
            return 0, 0
        trim_5 = int(germ_start)
        v_call = (mut_rec.get("v_call") or "").split(",")[0]
        sub_intervals = []
        for v_id in range(refdata.v_pool_size()):
            a = refdata.v_allele(v_id)
            if a.name == v_call:
                sub_intervals = list(a.subregions)
                break
        cdr_hits = 0
        fwr_hits = 0
        for i in range(int(v_end)):
            if base_seq[i] == mut_seq[i]:
                continue
            allele_local = i + trim_5
            for label, s, e in sub_intervals:
                if s <= allele_local < e:
                    if label in ("CDR1", "CDR2"):
                        cdr_hits += 1
                    elif label in ("FWR1", "FWR2", "FWR3"):
                        fwr_hits += 1
                    break
        return cdr_hits, fwr_hits

    seed = 9090
    n = 50
    base_exp = (
        ga.Experiment.on("human_igh").recombine().productive_only()
    )
    base = base_exp.run_records(n=n, seed=seed)

    # ── Configuration (a): CDR zeroed ─────────────────────────
    cdr_off_exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .productive_only()
        .mutate(
            model="s5f",
            rate=0.05,
            v_subregion_rates={"CDR": 0.0, "FWR": 1.0},
        )
    )
    cdr_off_refdata = cdr_off_exp.refdata
    cdr_off = cdr_off_exp.run_records(n=n, seed=seed)
    cdr_total = 0
    fwr_total = 0
    for base_rec, mut_rec in zip(base, cdr_off):
        c, f = _v_classify(mut_rec, base_rec, cdr_off_refdata)
        cdr_total += c
        fwr_total += f
    assert cdr_total == 0, (
        f"CDR-zeroed full stack produced {cdr_total} CDR1/CDR2 "
        "mutations across 50 records; expected 0 (zero-rate V "
        "subregion must drop sites from support before contract "
        "admissibility)"
    )
    assert fwr_total > 0, (
        "CDR-zeroed full stack produced zero FWR mutations either; "
        "the V-subregion rate path isn't exercising the FWR-only "
        "configuration"
    )

    # ── Configuration (b): FWR zeroed (symmetric) ─────────────
    fwr_off_exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .productive_only()
        .mutate(
            model="s5f",
            rate=0.05,
            v_subregion_rates={"FWR": 0.0, "CDR": 1.0},
        )
    )
    fwr_off_refdata = fwr_off_exp.refdata
    fwr_off = fwr_off_exp.run_records(n=n, seed=seed)
    cdr_total = 0
    fwr_total = 0
    for base_rec, mut_rec in zip(base, fwr_off):
        c, f = _v_classify(mut_rec, base_rec, fwr_off_refdata)
        cdr_total += c
        fwr_total += f
    assert fwr_total == 0, (
        f"FWR-zeroed full stack produced {fwr_total} FWR1/FWR2/FWR3 "
        "mutations across 50 records; expected 0"
    )
    assert cdr_total > 0, (
        "FWR-zeroed full stack produced zero CDR mutations either; "
        "the V-subregion rate path isn't exercising the CDR-only "
        "configuration"
    )


# ──────────────────────────────────────────────────────────────────
# V-subregion mutation counters — release-tier consolidation
#
# Closes the V-subregion mutation-counters slice with the same
# closure standard as the rate / per-segment counter slices:
# productive IGH full stack + non-default (segment_rates ×
# v_subregion_rates) + AIRR validator + cache parity +
# two-layer partition invariants. See
# ``docs/v_subregion_mutation_counters_audit.md`` for the
# architecture contract and
# ``tests/test_v_subregion_mutation_counters_implementation.py``
# for the per-spec implementation tests.
# ──────────────────────────────────────────────────────────────────


def test_v_subregion_mutation_counters_full_stack_partition_and_validate():
    """Productive IGH full stack with non-default
    ``segment_rates`` × ``v_subregion_rates`` + heavy corruption +
    paired-end. Four gates:

    1. ``validate_records(refdata)`` runs clean — none of the
       seven new V-subregion validator issue kinds fire on engine-
       projected records.
    2. Live-call cache parity holds on every outcome.
    3. **Two-layer partition invariant** on every record:
         - ``n_v + n_d + n_j + n_np == n_mutations``
           (per-segment partition — pre-existing).
         - ``n_fwr1 + n_cdr1 + n_fwr2 + n_cdr2 + n_fwr3
            + n_v_unannotated == n_v_mutations``
           (V-subregion partition — the counters slice's
           headline claim).
    4. **Non-vacuous distribution**: across the batch, at least
       one FWR bucket AND at least one CDR bucket accrue a non-
       zero count. Catches a hypothetical regression where the
       attribution silently routed every event to the unannotated
       bucket.
    """
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .invert_d(prob=0.3)
        .receptor_revision(prob=0.3)
        .productive_only()
        .mutate(
            model="s5f",
            rate=0.05,
            segment_rates={"V": 1.0, "D": 0.3, "J": 0.5, "NP": 0.2},
            v_subregion_rates={"CDR": 2.5, "FWR": 0.5},
        )
        .pcr_amplify(rate=1e-4)
        .sequencing_errors(count=4)
        .polymerase_indels(count=2)
        .primer_trim_5prime(length=(0, 3))
        .primer_trim_3prime(length=(0, 3))
        .paired_end(r1_length=150, insert_size=300)
    )
    refdata = exp.refdata
    result = exp.run_records(n=50, seed=4242)

    # Gate 1: AIRR validator. The six per-field
    # `N<Region>MutationsMismatch` issue kinds and the cross-field
    # `VSubregionMutationCountSumMismatch` invariant must all be
    # silent on engine-projected records.
    report = result.validate_records(refdata)
    assert report, (
        f"V-subregion-counter full stack tripped validator: "
        f"summary={report.summary()}; first failure="
        f"{report.failures[0] if report.failures else None}"
    )

    # Gate 2: live-call cache parity. Counters are pure
    # projection and don't touch live-call state, but the
    # cross-layer sanity check stays.
    assert result.outcomes is not None
    for i, outcome in enumerate(result.outcomes):
        for p in outcome.check_live_call_cache_parity(refdata):
            assert p["tie_set_matches"], (
                f"counter full-stack record {i} segment {p['segment']}: "
                f"cache parity failed — cached={p['cached_tie_set']} "
                f"fresh={p['fresh_tie_set']}"
            )

    # Gate 3: two-layer partition invariant on every record.
    for r in result:
        per_segment_total = (
            r["n_v_mutations"]
            + r["n_d_mutations"]
            + r["n_j_mutations"]
            + r["n_np_mutations"]
        )
        assert per_segment_total == r["n_mutations"], (
            f"per-segment sum invariant violated: "
            f"V={r['n_v_mutations']} + D={r['n_d_mutations']} + "
            f"J={r['n_j_mutations']} + NP={r['n_np_mutations']} = "
            f"{per_segment_total} != n_mutations={r['n_mutations']}"
        )
        v_subregion_total = (
            r["n_fwr1_mutations"]
            + r["n_cdr1_mutations"]
            + r["n_fwr2_mutations"]
            + r["n_cdr2_mutations"]
            + r["n_fwr3_mutations"]
            + r["n_v_unannotated_mutations"]
        )
        assert v_subregion_total == r["n_v_mutations"], (
            f"V-subregion sum invariant violated: "
            f"FWR1={r['n_fwr1_mutations']} + CDR1={r['n_cdr1_mutations']} + "
            f"FWR2={r['n_fwr2_mutations']} + CDR2={r['n_cdr2_mutations']} + "
            f"FWR3={r['n_fwr3_mutations']} + "
            f"unannotated={r['n_v_unannotated_mutations']} = "
            f"{v_subregion_total} != n_v_mutations={r['n_v_mutations']}"
        )

    # Gate 4: non-vacuous distribution — at least one FWR AND at
    # least one CDR bucket gets a hit across the batch. Catches a
    # hypothetical regression where the attribution silently
    # routed every event to one bucket (e.g. unannotated, or
    # FWR3 catching everything because of an off-by-one).
    batch_fwr_total = (
        sum(r["n_fwr1_mutations"] for r in result)
        + sum(r["n_fwr2_mutations"] for r in result)
        + sum(r["n_fwr3_mutations"] for r in result)
    )
    batch_cdr_total = (
        sum(r["n_cdr1_mutations"] for r in result)
        + sum(r["n_cdr2_mutations"] for r in result)
    )
    assert batch_fwr_total > 0, (
        "V-subregion counter full stack accrued zero FWR mutations "
        "across 50 records — distribution test is vacuous, or the "
        "FWR attribution path silently broke"
    )
    assert batch_cdr_total > 0, (
        "V-subregion counter full stack accrued zero CDR mutations "
        "across 50 records — distribution test is vacuous, or the "
        "CDR attribution path silently broke"
    )

    # Paired-end fields populated on every record (smoke check
    # that the full stack ran through to the last pass).
    for rec in result.records:
        assert rec["r1_sequence"], "paired-end r1_sequence empty"
        assert rec["r2_sequence"], "paired-end r2_sequence empty"


# ──────────────────────────────────────────────────────────────────
# Paired FASTQ Export — release-tier smoke
#
# Closes the paired-end FASTQ export slice with a release-grade
# end-to-end check: a productive IGH paired-end pipeline →
# `to_paired_fastq` → readback two files → assert structural
# invariants. The slice's spec tests live in
# `tests/test_to_paired_fastq.py`; this is the
# "did it survive the full stack" smoke. See
# `docs/fastq_export_design.md`.
# ──────────────────────────────────────────────────────────────────


def test_paired_fastq_export_full_stack_smoke(tmp_path):
    """Productive IGH + heavy stack + paired-end → write two
    FASTQ files → read them back and pin four invariants:

    1. **File record count == result record count.** Each AIRR
       record produces exactly four lines per file.
    2. **Headers match `@{sequence_id}/1` and `@{sequence_id}/2`**
       record-by-record.
    3. **R1 / R2 sequence body lines match the AIRR
       `r1_sequence` / `r2_sequence` fields verbatim** (R2 is
       already RC at projection time; the writer outputs it as-is).
    4. **Quality string lengths match the corresponding sequence
       lengths** on every read (R1 and R2 scored independently).

    Smoke-tier only — the spec tests in
    `tests/test_to_paired_fastq.py` cover the
    error-path / overwrite / quality-model variants. This test
    confirms the export survives the same full-stack pipeline the
    other release-tier tests run.
    """
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .productive_only()
        .mutate(model="s5f", rate=0.03)
        .pcr_amplify(rate=1e-4)
        .polymerase_indels(count=2)
        .primer_trim_5prime(length=(0, 3))
        .primer_trim_3prime(length=(0, 3))
        .paired_end(r1_length=150, insert_size=300)
    )
    result = exp.run_records(n=20, seed=4242)
    r1_path = tmp_path / "reads_R1.fastq"
    r2_path = tmp_path / "reads_R2.fastq"
    result.to_paired_fastq(str(r1_path), str(r2_path))

    r1_lines = r1_path.read_text(encoding="utf-8").splitlines()
    r2_lines = r2_path.read_text(encoding="utf-8").splitlines()
    n = len(result.records)

    # Gate 1: file record count == result record count.
    assert len(r1_lines) == n * 4, (
        f"R1 file has {len(r1_lines)} lines for {n} records; "
        "expected n*4 (4 lines per FASTQ record)"
    )
    assert len(r2_lines) == n * 4, (
        f"R2 file has {len(r2_lines)} lines for {n} records; "
        "expected n*4 (4 lines per FASTQ record)"
    )

    for i, rec in enumerate(result.records):
        sequence_id = rec["sequence_id"]
        # Gate 2: header format.
        assert r1_lines[i * 4] == f"@{sequence_id}/1", (
            f"R1 header mismatch at record {i}: "
            f"got {r1_lines[i * 4]!r}, expected @{sequence_id}/1"
        )
        assert r2_lines[i * 4] == f"@{sequence_id}/2", (
            f"R2 header mismatch at record {i}: "
            f"got {r2_lines[i * 4]!r}, expected @{sequence_id}/2"
        )
        # Gate 3: sequence body equality (R2 verbatim, no second RC).
        r1_seq_line = r1_lines[i * 4 + 1]
        r2_seq_line = r2_lines[i * 4 + 1]
        assert r1_seq_line == rec["r1_sequence"].upper(), (
            f"R1 body mismatch at record {i}: {r1_seq_line[:40]!r}"
        )
        assert r2_seq_line == rec["r2_sequence"].upper(), (
            f"R2 body mismatch at record {i}: {r2_seq_line[:40]!r}"
        )
        # FASTQ separator line.
        assert r1_lines[i * 4 + 2] == "+"
        assert r2_lines[i * 4 + 2] == "+"
        # Gate 4: quality string length parity.
        r1_q = r1_lines[i * 4 + 3]
        r2_q = r2_lines[i * 4 + 3]
        assert len(r1_q) == len(r1_seq_line), (
            f"R1 record {i}: quality len {len(r1_q)} != "
            f"sequence len {len(r1_seq_line)}"
        )
        assert len(r2_q) == len(r2_seq_line), (
            f"R2 record {i}: quality len {len(r2_q)} != "
            f"sequence len {len(r2_seq_line)}"
        )
