"""End-to-end tests for `Experiment.invert_d(prob=...)` — Slice D.

Pins the user-facing DSL surface that exposes V(D)J inversion as a
fluent step:

- VJ chains reject the method at call time.
- prob=0 records a `Bool(false)` and leaves the D segment forward.
- prob=1 records a `Bool(true)` and the pool carries the
  reverse-complement of the D allele at the D region.
- Replay round-trip preserves the inversion decision bit-for-bit.
- The trace's `sample_allele.d.inverted` record sits in the correct
  schedule position (between NP1 generation and NP2 generation —
  i.e., after the D sample, before the D assembly emits bytes).
- Calling `.invert_d()` twice raises `ValueError`.

No `d_inverted` AIRR field assertions here — that surface lands in
Slice E.
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
    """Tiny VDJ cartridge. D allele = ``ACGTTA`` so the reverse-
    complement ``TAACGT`` is distinguishable from the forward bytes
    at every position — every byte changes under inversion."""
    cfg = ge.RefDataConfig.vdj()
    cfg.add_v_allele("v1*01", "v1", b"AAACCCGGG", anchor=6)
    cfg.add_d_allele("d1*01", "d1", b"ACGTTA")
    cfg.add_j_allele("j1*01", "j1", b"TTTAAA", anchor=0)
    return cfg


def _baseline(invert: bool) -> "ga.Experiment":
    """VDJ pipeline with fixed-length NP1/NP2 + no trim. When
    ``invert`` is true the experiment carries `invert_d(prob=1.0)`
    so the D assignment is always reverse-complemented."""
    exp = (
        ga.Experiment.on(_vdj_refdata())
        .allow_curatable_refdata()
        .recombine(np1_lengths=[(3, 1.0)], np2_lengths=[(3, 1.0)])
        .trim(enabled=False)
    )
    if invert:
        exp = exp.invert_d(prob=1.0)
    return exp


# ──────────────────────────────────────────────────────────────────
# 1. VJ chain rejects invert_d at call time
# ──────────────────────────────────────────────────────────────────


def test_invert_d_on_vj_chain_raises_value_error() -> None:
    """``invert_d`` is meaningless on VJ chains — there is no D
    pool to invert. The DSL surface rejects the combination as
    soon as the user calls the method, not later at compile time."""
    cfg = ge.RefDataConfig.vj()
    cfg.add_v_allele("v1*01", "v1", b"TGTAAACCC", anchor=0)
    cfg.add_j_allele("j1*01", "j1", b"TTTAAA", anchor=0)

    with pytest.raises(ValueError, match="only valid for VDJ chains"):
        ga.Experiment.on(cfg).allow_curatable_refdata().invert_d(prob=0.05)


def test_invert_d_on_bundled_human_igk_raises_value_error() -> None:
    """The same guard fires for bundled light-chain catalogues —
    `human_igk` is a VJ cartridge."""
    with pytest.raises(ValueError, match="only valid for VDJ chains"):
        ga.Experiment.on("human_igk").invert_d()


# ──────────────────────────────────────────────────────────────────
# 2. prob=0 records false, D assembles forward
# ──────────────────────────────────────────────────────────────────


def test_invert_d_prob_zero_records_false_and_keeps_d_forward() -> None:
    exp = _baseline(invert=False).invert_d(prob=0.0)
    outcomes = exp.run(n=1, seed=0)
    addrs = {r.address for r in outcomes[0].trace().choices()}
    assert "sample_allele.d.inverted" in addrs
    inv_rec = next(
        r for r in outcomes[0].trace().choices()
        if r.address == "sample_allele.d.inverted"
    )
    assert inv_rec.value is False

    # The D bytes in the assembled sequence equal the forward
    # allele bytes — `ACGTTA` at the D region.
    rec = exp.run_records(n=1, seed=0).records[0]
    d_start = rec["d_sequence_start"]
    d_end = rec["d_sequence_end"]
    assert d_start is not None and d_end is not None
    assert rec["sequence"][d_start:d_end] == "ACGTTA"


# ──────────────────────────────────────────────────────────────────
# 3. prob=1 records true, D assembles reverse-complemented
# ──────────────────────────────────────────────────────────────────


def test_invert_d_prob_one_records_true_and_emits_reverse_complement_d() -> None:
    """Per the design doc: when the inversion decision is true, the
    bytes the pool carries at the D region are the
    reverse-complement of the D allele. D allele = ``ACGTTA`` →
    expected pool D = ``TAACGT``."""
    exp = _baseline(invert=True)
    outcomes = exp.run(n=1, seed=0)
    inv_rec = next(
        r for r in outcomes[0].trace().choices()
        if r.address == "sample_allele.d.inverted"
    )
    assert inv_rec.value is True

    rec = exp.run_records(n=1, seed=0).records[0]
    d_start = rec["d_sequence_start"]
    d_end = rec["d_sequence_end"]
    assert d_start is not None and d_end is not None
    assert rec["sequence"][d_start:d_end] == "TAACGT", (
        f"D region should be the reverse-complement of ACGTTA = TAACGT; "
        f"got {rec['sequence'][d_start:d_end]!r}"
    )


def test_invert_d_only_flips_d_region_bytes() -> None:
    """V and J bytes must be byte-identical between the forward and
    inverted runs (same seed). Only the D region's bytes differ."""
    fwd = _baseline(invert=False).run_records(n=1, seed=0).records[0]
    rc = _baseline(invert=True).run_records(n=1, seed=0).records[0]

    # V byte-for-byte identical.
    v_a = fwd["sequence"][fwd["v_sequence_start"]:fwd["v_sequence_end"]]
    v_b = rc["sequence"][rc["v_sequence_start"]:rc["v_sequence_end"]]
    assert v_a == v_b

    # J byte-for-byte identical (same seed, same NP1/NP2 RNG path
    # before inversion fires).
    j_a = fwd["sequence"][fwd["j_sequence_start"]:fwd["j_sequence_end"]]
    j_b = rc["sequence"][rc["j_sequence_start"]:rc["j_sequence_end"]]
    assert j_a == j_b

    # D differs (forward vs RC).
    d_a = fwd["sequence"][fwd["d_sequence_start"]:fwd["d_sequence_end"]]
    d_b = rc["sequence"][rc["d_sequence_start"]:rc["d_sequence_end"]]
    assert d_a == "ACGTTA"
    assert d_b == "TAACGT"


# ──────────────────────────────────────────────────────────────────
# 4. Replay round-trip
# ──────────────────────────────────────────────────────────────────


def test_invert_d_trace_replays_bit_for_bit() -> None:
    """A trace recorded by an `invert_d(prob=1.0)` run replays
    through `rerun_from_trace_file` and reproduces the same trace
    and the same assembled sequence."""
    exp = _baseline(invert=True)
    compiled = exp.compile()
    fresh_outcome = compiled.simulator.run(seed=0)
    fresh_trace = fresh_outcome.trace()

    tf = compiled.simulator.trace_file_from(fresh_outcome, seed=0)
    replayed = compiled.simulator.rerun_from_trace_file(tf)

    # Trace orientation choice equals the fresh value.
    fresh_inv = next(
        r for r in fresh_trace.choices()
        if r.address == "sample_allele.d.inverted"
    )
    replayed_inv = next(
        r for r in replayed.trace().choices()
        if r.address == "sample_allele.d.inverted"
    )
    assert fresh_inv.value == replayed_inv.value

    # And the assembled-sequence bytes round-trip.
    from GenAIRR._airr_record import outcome_to_airr_record
    fresh_rec = outcome_to_airr_record(
        fresh_outcome, exp._refdata, sequence_id="fresh"
    )
    replayed_rec = outcome_to_airr_record(
        replayed, exp._refdata, sequence_id="replayed"
    )
    assert fresh_rec["sequence"] == replayed_rec["sequence"]


# ──────────────────────────────────────────────────────────────────
# 5. Schedule ordering: invert_d sits between NP1 and assemble.D
# ──────────────────────────────────────────────────────────────────


def test_invert_d_trace_record_sits_after_np1_and_before_np2() -> None:
    """The schedule analyser must place InvertDPass between
    ``generate.np1`` and ``assemble.d``. Since assembly passes
    don't emit trace records, the observable proxy is the
    relative position of ``sample_allele.d.inverted`` versus the
    NP1 / NP2 records — invert_d sits AFTER np.np1.bases[N] and
    BEFORE np.np2.length. This pins the pool layout staying
    V-NP1-D-NP2-J under inversion."""
    exp = _baseline(invert=True)
    outcomes = exp.run(n=1, seed=0)
    addrs = [r.address for r in outcomes[0].trace().choices()]
    inv_idx = addrs.index("sample_allele.d.inverted")
    np1_last_idx = max(
        i for i, a in enumerate(addrs) if a.startswith("np.np1.bases[")
    )
    np2_first_idx = addrs.index("np.np2.length")
    assert np1_last_idx < inv_idx < np2_first_idx, (
        "sample_allele.d.inverted must commit AFTER NP1 generation "
        "and BEFORE NP2 generation so the pool layout stays "
        "V-NP1-D-NP2-J. Got trace order: " + ", ".join(addrs)
    )


def test_baseline_run_without_invert_d_omits_the_address() -> None:
    """Negative control: when the user does NOT call ``invert_d``,
    no `sample_allele.d.inverted` record appears in the trace.
    Slice D shouldn't change the baseline-pipeline trace shape."""
    exp = _baseline(invert=False)
    outcomes = exp.run(n=1, seed=0)
    addrs = {r.address for r in outcomes[0].trace().choices()}
    assert "sample_allele.d.inverted" not in addrs


# ──────────────────────────────────────────────────────────────────
# 6. Duplicate `.invert_d()` calls rejected
# ──────────────────────────────────────────────────────────────────


def test_invert_d_called_twice_raises_value_error() -> None:
    """v1 picks a single inversion probability per pipeline; the
    DSL rejects the second call so an over-eager builder gets a
    loud error instead of silent last-wins semantics."""
    exp = _baseline(invert=False)
    exp.invert_d(prob=0.5)
    with pytest.raises(ValueError, match="already configured"):
        exp.invert_d(prob=0.05)


def test_invert_d_prob_out_of_range_raises_value_error() -> None:
    """Mirror of the InvertDPass Rust-side guard. The DSL must
    surface the bad probability before reaching the engine."""
    exp = _baseline(invert=False)
    with pytest.raises(ValueError, match=r"\[0.0, 1.0\]"):
        exp.invert_d(prob=1.5)


def test_invert_d_nan_prob_raises_value_error() -> None:
    exp = _baseline(invert=False)
    with pytest.raises(ValueError, match="NaN"):
        exp.invert_d(prob=float("nan"))


# ──────────────────────────────────────────────────────────────────
# 7. Resolver-bridge sanity: the DSL works against bundled VDJ data
# ──────────────────────────────────────────────────────────────────


def test_invert_d_on_bundled_human_igh_records_address_in_trace() -> None:
    """Smoke test against the real bundled human_igh cartridge —
    confirms the DSL composes with the rest of the recombine
    pipeline and the trace address fires on a non-synthetic
    catalogue."""
    cfg = _refdata_resolver._resolve_config_name("human_igh")
    exp = (
        ga.Experiment.on(cfg)
        .allow_curatable_refdata()
        .recombine()
        .invert_d(prob=1.0)
    )
    outcomes = exp.run(n=1, seed=42)
    addrs = {r.address for r in outcomes[0].trace().choices()}
    assert "sample_allele.d.inverted" in addrs


# ──────────────────────────────────────────────────────────────────
# Slice E: AIRR `d_inverted` field
# ──────────────────────────────────────────────────────────────────


def test_d_inverted_airr_field_true_under_prob_one() -> None:
    """The canonical Slice E surface: `d_inverted` in the AIRR
    record dict reflects whether the D allele was committed in
    reverse-complement orientation. prob=1.0 → True."""
    rec = _baseline(invert=True).run_records(n=1, seed=0).records[0]
    assert rec["d_inverted"] is True


def test_d_inverted_false_under_prob_zero() -> None:
    """Symmetric: a Bool(false) commit (or no commit) means the D
    allele stayed Forward, and `d_inverted` is False."""
    exp = _baseline(invert=False).invert_d(prob=0.0)
    rec = exp.run_records(n=1, seed=0).records[0]
    assert rec["d_inverted"] is False


def test_d_inverted_false_without_invert_d_step() -> None:
    """Baseline behaviour: a VDJ experiment that NEVER calls
    `.invert_d()` has no inversion step at all; `d_inverted` is
    False everywhere. Pin this so a future refactor that
    accidentally defaults the field to True breaks loudly."""
    rec = _baseline(invert=False).run_records(n=1, seed=0).records[0]
    assert rec["d_inverted"] is False


def test_d_inverted_false_for_vj_chain() -> None:
    """VJ chains have no D pool. The builder's `unwrap_or(false)`
    on the missing D assignment surfaces here as a False value in
    every record."""
    cfg = ge.RefDataConfig.vj()
    cfg.add_v_allele("v1*01", "v1", b"TGTAAACCC", anchor=0)
    cfg.add_j_allele("j1*01", "j1", b"TTTAAA", anchor=0)
    rec = (
        ga.Experiment.on(cfg)
        .allow_curatable_refdata()
        .recombine(np1_lengths=[(3, 1.0)])
        .trim(enabled=False)
        .run_records(n=1, seed=0)
        .records[0]
    )
    assert rec["d_inverted"] is False


def test_d_inverted_replay_round_trip_preserves_value() -> None:
    """A fresh inversion decision survives a trace-replay round
    trip. The field's source is the final IR (not the trace), so
    the replayed sim reproduces the same orientation and the
    AIRR record reflects it."""
    from GenAIRR._airr_record import outcome_to_airr_record

    exp = _baseline(invert=True)
    compiled = exp.compile()
    fresh_outcome = compiled.simulator.run(seed=0)
    tf = compiled.simulator.trace_file_from(fresh_outcome, seed=0)
    replayed = compiled.simulator.rerun_from_trace_file(tf)

    fresh_rec = outcome_to_airr_record(
        fresh_outcome, exp._refdata, sequence_id="fresh"
    )
    replayed_rec = outcome_to_airr_record(
        replayed, exp._refdata, sequence_id="replayed"
    )
    assert fresh_rec["d_inverted"] is True
    assert replayed_rec["d_inverted"] == fresh_rec["d_inverted"]


def test_d_inverted_appears_in_simulation_result_dataframe_columns() -> None:
    """The user-facing CSV / DataFrame exports must include
    `d_inverted`. `SimulationResult` flattens the per-record
    dicts into columns; pin that the new column lands on the
    canonical list (next to `is_contaminant`)."""
    result = _baseline(invert=True).run_records(n=1, seed=0)
    # Quickly verify the column flows through to the column list
    # used by `to_csv` / `to_dataframe`.
    df = result.to_dataframe()
    assert "d_inverted" in df.columns
    assert df["d_inverted"].iloc[0] is True or df["d_inverted"].iloc[0]


def test_validator_check_passes_for_clean_inverted_record() -> None:
    """The AIRR record-validator runs the engine's `d_inverted`
    consistency check among its other postconditions. Pin that a
    clean record (built by the canonical path) never trips the
    `DInvertedMismatch` issue under prob=1."""
    exp = _baseline(invert=True)
    outcomes = exp.run(n=1, seed=0)
    issues = outcomes[0].validate_record(exp._refdata, sequence_id="clean")
    kinds = {issue["kind"] for issue in issues}
    assert "DInvertedMismatch" not in kinds


def test_validator_check_passes_for_forward_baseline_record() -> None:
    """Symmetric: the baseline forward path also passes the
    `d_inverted` consistency check. Together with the
    `prob=1`-clean test, the two pin that builder + validator
    agree on expected values across the full orientation surface."""
    exp = _baseline(invert=False)
    outcomes = exp.run(n=1, seed=0)
    issues = outcomes[0].validate_record(exp._refdata, sequence_id="clean")
    kinds = {issue["kind"] for issue in issues}
    assert "DInvertedMismatch" not in kinds
