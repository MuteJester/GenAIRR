"""End-to-end tests for `Experiment.receptor_revision(prob=...)` — Slice D.

Pins the user-facing DSL surface that exposes V-segment receptor
revision as a fluent step:

- VJ chains reject the method at call time.
- prob=0 records a `Bool(false)` and leaves V untouched.
- prob=1 records the three-record set (`applied`, `v_allele`,
  `v_trim_3`) and rewrites the V slice with the replacement allele's
  retained prefix.
- Replay round-trip preserves the revision decision bit-for-bit.
- The trace's `receptor_revision.applied` record sits AFTER initial
  recombination choices and BEFORE any mutation/corruption records.
- Calling `.receptor_revision()` twice raises `ValueError`.
- prob validation (NaN, out-of-range) raises at the DSL boundary.

No `receptor_revision_applied` / `original_v_call` AIRR field
assertions here — that surface lands in Slice E.
"""
from __future__ import annotations

import pytest

import GenAIRR as ga
from GenAIRR import _engine as ge
from GenAIRR import _refdata_resolver


# ──────────────────────────────────────────────────────────────────
# Deterministic VDJ fixture
# ──────────────────────────────────────────────────────────────────


def _vdj_refdata() -> "ge.RefDataConfig":
    """Tiny VDJ cartridge with two distinguishable V alleles. Both
    are length 9 so the same-length-retained constraint is trivially
    satisfied at ``trim_3 = 0`` regardless of which candidate
    receptor revision picks. The originally-assembled V
    (``AAAAAAAAA``) and the replacement (``CCCCCCCCC``) differ at
    every byte position, so the post-revision pool is unambiguously
    distinguishable from the no-revision baseline."""
    cfg = ge.RefDataConfig.vdj()
    cfg.add_v_allele("v1*01", "v1", b"AAAAAAAAA", anchor=6)
    cfg.add_v_allele("v2*01", "v2", b"CCCCCCCCC", anchor=6)
    cfg.add_d_allele("d1*01", "d1", b"GGGGGG")
    cfg.add_j_allele("j1*01", "j1", b"TTTAAA", anchor=0)
    return cfg


def _baseline(*, revise: bool, prob: float = 1.0) -> "ga.Experiment":
    """VDJ pipeline with fixed-length NP1/NP2 + no trim. When
    ``revise`` is true the experiment carries
    ``receptor_revision(prob=prob)`` so the V slot is rewritten with
    one of the two same-length V candidates."""
    exp = (
        ga.Experiment.on(_vdj_refdata())
        .allow_curatable_refdata()
        .recombine(np1_lengths=[(3, 1.0)], np2_lengths=[(3, 1.0)])
        .trim(enabled=False)
    )
    if revise:
        exp = exp.receptor_revision(prob=prob)
    return exp


# ──────────────────────────────────────────────────────────────────
# 1. VJ chain rejects receptor_revision at call time
# ──────────────────────────────────────────────────────────────────


def test_receptor_revision_on_vj_chain_raises_value_error() -> None:
    """``receptor_revision`` is heavy-chain v1; light-chain VJ
    biology and pipeline don't carry the same revision mechanism.
    The DSL surface rejects the combination as soon as the user
    calls the method, not later at compile time."""
    cfg = ge.RefDataConfig.vj()
    cfg.add_v_allele("v1*01", "v1", b"TGTAAACCC", anchor=0)
    cfg.add_j_allele("j1*01", "j1", b"TTTAAA", anchor=0)

    with pytest.raises(ValueError, match="only valid for VDJ chains"):
        ga.Experiment.on(cfg).allow_curatable_refdata().receptor_revision()


def test_receptor_revision_on_bundled_human_igk_raises_value_error() -> None:
    """The same guard fires for bundled light-chain catalogues —
    ``human_igk`` is a VJ cartridge."""
    with pytest.raises(ValueError, match="only valid for VDJ chains"):
        ga.Experiment.on("human_igk").receptor_revision()


# ──────────────────────────────────────────────────────────────────
# 2. prob=0 records false, V stays intact
# ──────────────────────────────────────────────────────────────────


def test_receptor_revision_prob_zero_records_false_and_no_allele_or_trim() -> None:
    """The Bool is always recorded, even when the revision doesn't
    fire — replay-determinism needs the outcome on the trace
    regardless. The two conditional records must NOT appear."""
    exp = _baseline(revise=True, prob=0.0)
    outcomes = exp.run(n=1, seed=0)
    addrs = {r.address for r in outcomes[0].trace().choices()}

    assert "receptor_revision.applied" in addrs
    applied_rec = next(
        r for r in outcomes[0].trace().choices()
        if r.address == "receptor_revision.applied"
    )
    assert applied_rec.value is False
    assert "receptor_revision.v_allele" not in addrs
    assert "receptor_revision.v_trim_3" not in addrs


def test_receptor_revision_prob_zero_matches_no_revision_baseline_at_v_region() -> None:
    """prob=0 must commit no V mutation. The V region bytes between
    the no-revision baseline and the prob=0 revision experiment
    are byte-identical at the same seed (the revision Bool draws an
    RNG slot but a False outcome doesn't dip into the allele/trim
    streams). Pins that the no-op revision path is genuinely
    pool-neutral."""
    baseline_rec = _baseline(revise=False).run_records(n=1, seed=0).records[0]
    revised_rec = _baseline(revise=True, prob=0.0).run_records(n=1, seed=0).records[0]
    base_v = baseline_rec["sequence"][
        baseline_rec["v_sequence_start"]:baseline_rec["v_sequence_end"]
    ]
    rev_v = revised_rec["sequence"][
        revised_rec["v_sequence_start"]:revised_rec["v_sequence_end"]
    ]
    assert base_v == rev_v


# ──────────────────────────────────────────────────────────────────
# 3. prob=1 records all three choices and rewrites V
# ──────────────────────────────────────────────────────────────────


def test_receptor_revision_prob_one_records_three_choices_and_rewrites_v() -> None:
    """When the revision fires:

    1. ``receptor_revision.applied`` = True
    2. ``receptor_revision.v_allele`` carries the replacement
       allele's id.
    3. ``receptor_revision.v_trim_3`` carries the derived 3' trim
       (zero in this fixture since both V alleles are length 9 and
       the old V region is length 9).
    4. The pool's V region bytes equal the replacement allele's
       retained prefix.

    The downstream D/J regions remain byte-identical to the
    no-revision baseline because the same-length constraint shifts
    nothing.
    """
    revised = _baseline(revise=True).run_records(n=1, seed=0).records[0]
    trace = _baseline(revise=True).run(n=1, seed=0)[0].trace()

    applied = next(r for r in trace.choices() if r.address == "receptor_revision.applied")
    v_allele_rec = next(r for r in trace.choices() if r.address == "receptor_revision.v_allele")
    v_trim_rec = next(r for r in trace.choices() if r.address == "receptor_revision.v_trim_3")

    assert applied.value is True
    # The derived trim equals (allele.len - old_v_len). Both V
    # alleles in the fixture are length 9 and the old V region
    # is length 9; both are length-eligible candidates so the
    # derived trim is zero either way.
    assert v_trim_rec.value == 0
    # The V region bytes equal the recorded allele's full 9 bytes.
    v_start = revised["v_sequence_start"]
    v_end = revised["v_sequence_end"]
    v_bytes = revised["sequence"][v_start:v_end]
    assert v_end - v_start == 9
    # Recorded allele id is one of the two V candidates; whichever
    # one fires, the bytes match its sequence.
    expected_by_id = {0: "AAAAAAAAA", 1: "CCCCCCCCC"}
    assert v_bytes == expected_by_id[v_allele_rec.value]


def test_receptor_revision_does_not_disturb_d_or_j_region_bytes() -> None:
    """Same-length receptor revision must rewrite V only — D and J
    bytes are byte-identical to the no-revision baseline at the
    same seed."""
    baseline_rec = _baseline(revise=False).run_records(n=1, seed=0).records[0]
    revised_rec = _baseline(revise=True).run_records(n=1, seed=0).records[0]

    base_d = baseline_rec["sequence"][
        baseline_rec["d_sequence_start"]:baseline_rec["d_sequence_end"]
    ]
    rev_d = revised_rec["sequence"][
        revised_rec["d_sequence_start"]:revised_rec["d_sequence_end"]
    ]
    assert base_d == rev_d, f"D bytes diverged: base={base_d!r} revised={rev_d!r}"

    base_j = baseline_rec["sequence"][
        baseline_rec["j_sequence_start"]:baseline_rec["j_sequence_end"]
    ]
    rev_j = revised_rec["sequence"][
        revised_rec["j_sequence_start"]:revised_rec["j_sequence_end"]
    ]
    assert base_j == rev_j, f"J bytes diverged: base={base_j!r} revised={rev_j!r}"


# ──────────────────────────────────────────────────────────────────
# 4. Replay round-trip
# ──────────────────────────────────────────────────────────────────


def test_receptor_revision_trace_replays_bit_for_bit() -> None:
    """A trace recorded by a ``receptor_revision(prob=1.0)`` run
    replays through ``rerun_from_trace_file`` and reproduces the
    same trace and assembled sequence."""
    exp = _baseline(revise=True)
    compiled = exp.compile()
    fresh_outcome = compiled.simulator.run(seed=0)
    fresh_trace = fresh_outcome.trace()

    tf = compiled.simulator.trace_file_from(fresh_outcome, seed=0)
    replayed = compiled.simulator.rerun_from_trace_file(tf)

    for address in (
        "receptor_revision.applied",
        "receptor_revision.v_allele",
        "receptor_revision.v_trim_3",
    ):
        fresh_rec = next((r for r in fresh_trace.choices() if r.address == address), None)
        replayed_rec = next(
            (r for r in replayed.trace().choices() if r.address == address),
            None,
        )
        if fresh_rec is None:
            assert replayed_rec is None, (
                f"replay records {address} that fresh run didn't"
            )
        else:
            assert replayed_rec is not None, (
                f"replay drops the {address} record"
            )
            assert fresh_rec.value == replayed_rec.value

    from GenAIRR._airr_record import outcome_to_airr_record
    fresh_rec = outcome_to_airr_record(fresh_outcome, exp._refdata, sequence_id="fresh")
    replayed_rec = outcome_to_airr_record(replayed, exp._refdata, sequence_id="replay")
    assert fresh_rec["sequence"] == replayed_rec["sequence"]


# ──────────────────────────────────────────────────────────────────
# 5. Schedule ordering: receptor_revision after recombine, before mutate
# ──────────────────────────────────────────────────────────────────


def test_receptor_revision_sits_after_recombine_and_before_mutate() -> None:
    """The design doc §2 requires
    ``recombine → receptor_revision → mutate/corrupt``. The trace
    address order is the observable proxy: ``sample_allele.v`` /
    ``np.np2.bases[N]`` must precede ``receptor_revision.applied``,
    which in turn must precede ``mutate.s5f.count``."""
    exp = (
        _baseline(revise=True)
        .mutate(rate=0.01)
    )
    addrs = [r.address for r in exp.run(n=1, seed=0)[0].trace().choices()]

    rev_idx = addrs.index("receptor_revision.applied")
    sample_v_idx = addrs.index("sample_allele.v")
    assert sample_v_idx < rev_idx, (
        f"sample_allele.v at {sample_v_idx} should precede "
        f"receptor_revision.applied at {rev_idx}; trace: {addrs}"
    )
    # The last NP2 base (proxy for "recombination is complete") sits
    # before the revision.
    np2_indices = [i for i, a in enumerate(addrs) if a.startswith("np.np2.bases[")]
    if np2_indices:
        assert max(np2_indices) < rev_idx
    # And the mutation count is recorded *after* the revision.
    mut_idx = addrs.index("mutate.s5f.count")
    assert rev_idx < mut_idx, (
        f"receptor_revision.applied at {rev_idx} should precede "
        f"mutate.s5f.count at {mut_idx}; trace: {addrs}"
    )


def test_baseline_run_without_receptor_revision_omits_the_addresses() -> None:
    """Negative control: when the user does NOT call
    ``receptor_revision``, no ``receptor_revision.*`` records appear
    in the trace. Slice D shouldn't change the baseline-pipeline
    trace shape."""
    exp = _baseline(revise=False)
    addrs = {r.address for r in exp.run(n=1, seed=0)[0].trace().choices()}
    for absent in (
        "receptor_revision.applied",
        "receptor_revision.v_allele",
        "receptor_revision.v_trim_3",
    ):
        assert absent not in addrs


# ──────────────────────────────────────────────────────────────────
# 6. Duplicate / validation guards
# ──────────────────────────────────────────────────────────────────


def test_receptor_revision_called_twice_raises_value_error() -> None:
    """v1 picks a single revision probability per pipeline; the DSL
    rejects the second call so an over-eager builder gets a loud
    error instead of silent last-wins semantics."""
    exp = _baseline(revise=False)
    exp.receptor_revision(prob=0.5)
    with pytest.raises(ValueError, match="already configured"):
        exp.receptor_revision(prob=0.05)


def test_receptor_revision_prob_out_of_range_raises_value_error() -> None:
    """Mirror of the ReceptorRevisionPass Rust-side guard. The DSL
    must surface the bad probability before reaching the engine."""
    exp = _baseline(revise=False)
    with pytest.raises(ValueError, match=r"\[0.0, 1.0\]"):
        exp.receptor_revision(prob=1.5)


def test_receptor_revision_negative_prob_raises_value_error() -> None:
    exp = _baseline(revise=False)
    with pytest.raises(ValueError, match=r"\[0.0, 1.0\]"):
        exp.receptor_revision(prob=-0.1)


def test_receptor_revision_nan_prob_raises_value_error() -> None:
    exp = _baseline(revise=False)
    with pytest.raises(ValueError, match="NaN"):
        exp.receptor_revision(prob=float("nan"))


def test_receptor_revision_non_numeric_prob_raises_value_error() -> None:
    exp = _baseline(revise=False)
    with pytest.raises(ValueError, match="must be a number"):
        exp.receptor_revision(prob="0.5")  # type: ignore[arg-type]


# ──────────────────────────────────────────────────────────────────
# 7. Productive pipeline composes
# ──────────────────────────────────────────────────────────────────


def test_receptor_revision_composes_with_productive_only_at_prob_one() -> None:
    """Productive contracts (``productive_only=True``) must compose
    with receptor revision without surfacing a contract violation —
    the same-length retained slice preserves V anchor positions and
    the downstream D/J geometry is untouched, so the active
    contract bundle's ``admits_post_event`` accepts the post-event
    state.

    This is a behavioral smoke test, not a contract algebra
    statement: if a future contract change ever rejects revision
    state, this test surfaces the regression before AIRR fields
    rely on it."""
    rec = (
        _baseline(revise=True)
        .productive_only()
        .run_records(n=1, seed=0)
        .records[0]
    )
    # Productive records carry a non-empty sequence and a resolved
    # V region — both confirm the pass survived contract arbitration.
    assert rec["sequence"]
    assert rec["v_sequence_start"] is not None


# ──────────────────────────────────────────────────────────────────
# 8. Resolver-bridge sanity: bundled human_igh works
# ──────────────────────────────────────────────────────────────────


def test_receptor_revision_on_bundled_human_igh_records_addresses_in_trace() -> None:
    """Smoke test against the real bundled human_igh cartridge —
    confirms the DSL composes with the rest of the recombine
    pipeline and at least the always-on Bool record fires."""
    cfg = _refdata_resolver._resolve_config_name("human_igh")
    exp = (
        ga.Experiment.on(cfg)
        .allow_curatable_refdata()
        .recombine()
        .receptor_revision(prob=1.0)
    )
    addrs = {r.address for r in exp.run(n=1, seed=42)[0].trace().choices()}
    assert "receptor_revision.applied" in addrs


# ──────────────────────────────────────────────────────────────────
# Slice E: AIRR fields + validator
# ──────────────────────────────────────────────────────────────────


def test_no_revision_step_yields_false_and_empty_airr_fields() -> None:
    """Negative control: the no-revision baseline carries the two
    new AIRR keys with default values (false / empty)."""
    rec = _baseline(revise=False).run_records(n=1, seed=0).records[0]
    assert rec["receptor_revision_applied"] is False
    assert rec["original_v_call"] == ""


def test_revision_prob_zero_yields_false_and_empty_airr_fields() -> None:
    """The pass ran but rolled applied=false. Both fields stay
    at their no-revision defaults — the empty sentinel reads as
    "no revision happened" without forcing a Bool cross-check."""
    rec = _baseline(revise=True, prob=0.0).run_records(n=1, seed=0).records[0]
    assert rec["receptor_revision_applied"] is False
    assert rec["original_v_call"] == ""


def test_revision_prob_one_populates_airr_fields_with_original_v_name() -> None:
    """When the revision fires, `receptor_revision_applied=True` and
    `original_v_call` carries the V allele name the recombine pass
    originally committed (resolved from the trace's first
    `sample_allele.v` record). `v_call` continues to report the
    post-revision identity, so the two strings may differ."""
    rec = _baseline(revise=True, prob=1.0).run_records(n=1, seed=0).records[0]
    assert rec["receptor_revision_applied"] is True
    assert rec["original_v_call"], (
        "applied=True must populate original_v_call from the trace"
    )
    # Both V candidates in the fixture have distinguishable bytes,
    # so the post-revision v_call may differ from the original.
    assert rec["original_v_call"] in ("v1*01", "v2*01")
    assert rec["v_call"] in ("v1*01", "v2*01")


def test_revision_can_produce_v_call_distinct_from_original_v_call() -> None:
    """Sample several seeds: at least one must produce
    `v_call != original_v_call`. The two-V fixture is symmetric so
    a particular seed could land on the same allele either way; the
    multi-seed sweep guarantees we observe the cross-allele case
    that the field's semantics rely on."""
    saw_distinct = False
    for seed in range(8):
        rec = _baseline(revise=True, prob=1.0).run_records(n=1, seed=seed).records[0]
        if rec["v_call"] != rec["original_v_call"]:
            saw_distinct = True
            break
    assert saw_distinct, (
        "across 8 seeds we never saw v_call differ from "
        "original_v_call; the two-V fixture should produce that "
        "case at least once."
    )


def test_revision_airr_fields_round_trip_through_replay() -> None:
    """Replay must preserve both Slice E fields bit-for-bit. A trace
    recorded by a `prob=1.0` run replays through ``rerun_from_trace_file``;
    the post-replay AIRR record carries the same
    `receptor_revision_applied` and `original_v_call`."""
    from GenAIRR._airr_record import outcome_to_airr_record
    exp = _baseline(revise=True, prob=1.0)
    compiled = exp.compile()
    fresh_outcome = compiled.simulator.run(seed=0)
    tf = compiled.simulator.trace_file_from(fresh_outcome, seed=0)
    replayed = compiled.simulator.rerun_from_trace_file(tf)

    fresh_rec = outcome_to_airr_record(fresh_outcome, exp._refdata, sequence_id="fresh")
    replayed_rec = outcome_to_airr_record(replayed, exp._refdata, sequence_id="replay")

    assert fresh_rec["receptor_revision_applied"] == replayed_rec["receptor_revision_applied"]
    assert fresh_rec["original_v_call"] == replayed_rec["original_v_call"]


def test_revision_airr_fields_appear_in_dataframe_columns() -> None:
    """The two new fields are wired into the canonical AIRR column
    list, so dataframe / CSV exports include them."""
    df = _baseline(revise=True, prob=1.0).run_records(n=2, seed=0).to_dataframe()
    assert "receptor_revision_applied" in df.columns
    assert "original_v_call" in df.columns


def test_validator_accepts_clean_revised_record_end_to_end() -> None:
    """Run the full validator over a clean post-revision record and
    confirm no `ReceptorRevisionAppliedMismatch` or
    `OriginalVCallMismatch` surfaces. This pins the end-to-end
    builder/validator agreement against a real run, not just unit-
    test trace fixtures."""
    exp = _baseline(revise=True, prob=1.0)
    compiled = exp.compile()
    outcome = compiled.simulator.run(seed=0)
    issues = outcome.validate_record(exp._refdata, sequence_id="clean")
    kinds = [issue["kind"] for issue in issues]
    assert "ReceptorRevisionAppliedMismatch" not in kinds, (
        f"validator flagged a clean revised record: {issues}"
    )
    assert "OriginalVCallMismatch" not in kinds, (
        f"validator flagged a clean revised record: {issues}"
    )
