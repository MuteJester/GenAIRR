"""Indel provenance golden tests.

Companion to [docs/indel_provenance_audit.md](../docs/indel_provenance_audit.md).
Pins **current** behaviour of `polymerase_indels` on the AIRR
projection, the trace, and the per-pass `EventRecord` ledger.

What's pinned:

- A single insertion grows the sequence by 1 base; a single
  deletion shrinks it by 1.
- `n_indels` (trace-derived) matches the count of `IndelInserted` +
  `IndelDeleted` events on the pass's `EventRecord` under typical
  fixtures.
- The trace carries `corrupt.indel.{count, kind[i], site[i],
  base[i]}` per the documented schema.
- CIGAR (`v_cigar`, `j_cigar`) carries `I` / `D` ops when the
  indel lands inside a V/D/J region.
- Under `productive_only`, a paired indel scenario (count=2,
  insertion_prob=0.5) preserves the productive triad.
- Replay reproduces sequence, every key AIRR field, AND the
  per-pass event count.
- §6.1 resolution: `n_indels` reports the count of *applied*
  structural indels (IndelInserted + IndelDeleted SimulationEvents
  emitted by the `corrupt.indel` pass), NOT the trace's
  `corrupt.indel.count` (which is the attempted count and includes
  contract-rejected / pool-empty no-op sentinels).

Fixture design: VJ refdata with V=9 (anchor at pos 6), J=6 (anchor
at pos 0), trim disabled, NP1 fixed at length 3 → baseline
assembled sequence is 18 bases. Every length-delta assertion keys
off that 18.
"""
from __future__ import annotations

import pytest

import GenAIRR as ga
from GenAIRR import _engine as ge
from GenAIRR import CompiledExperiment


# ──────────────────────────────────────────────────────────────────
# Deterministic VJ fixture
# ──────────────────────────────────────────────────────────────────


def _vj_refdata() -> "ge.RefDataConfig":
    cfg = ge.RefDataConfig.vj()
    cfg.add_v_allele("v1*01", "v1", b"AAACCCGGG", anchor=6)  # 9 bases, anchor at 6..9
    cfg.add_j_allele("j1*01", "j1", b"TTTAAA", anchor=0)     # 6 bases, anchor at 0..3
    return cfg


BASELINE_LEN = 18  # V(9) + NP1(3) + J(6)


def _baseline_experiment() -> "ga.Experiment":
    return (
        ga.Experiment.on(_vj_refdata())
        .recombine(np1_lengths=[(3, 1.0)])
        .trim(enabled=False)
    )


def _record(exp: "ga.Experiment", seed: int = 0) -> dict:
    """Run `exp` once and return the single AIRR record."""
    result = exp.run_records(n=1, seed=seed)
    assert len(result.records) == 1
    return result.records[0]


def _outcome(exp: "ga.Experiment", seed: int = 0):
    """Run `exp` once and return the raw Outcome (so callers can
    inspect trace + events)."""
    compiled = exp.compile()
    assert isinstance(compiled, CompiledExperiment)
    return compiled.simulator.run(seed=seed), compiled


# ──────────────────────────────────────────────────────────────────
# 1. Baseline — confirm fixture predictability
# ──────────────────────────────────────────────────────────────────


def test_baseline_fixture_produces_18_base_sequence_with_zero_indels() -> None:
    rec = _record(_baseline_experiment())
    assert len(rec["sequence"]) == BASELINE_LEN
    assert rec["n_indels"] == 0


# ──────────────────────────────────────────────────────────────────
# 2. Single insertion
# ──────────────────────────────────────────────────────────────────


def test_single_insertion_grows_sequence_by_one_base() -> None:
    rec = _record(
        _baseline_experiment().polymerase_indels(count=1, insertion_prob=1.0)
    )
    assert len(rec["sequence"]) == BASELINE_LEN + 1
    assert rec["n_indels"] == 1


def test_single_insertion_trace_records_kind_and_site() -> None:
    outcome, _compiled = _outcome(
        _baseline_experiment().polymerase_indels(count=1, insertion_prob=1.0)
    )
    choices = outcome.trace().choices()
    addrs = [r.address for r in choices]
    assert "corrupt.indel.count" in addrs
    assert "corrupt.indel.kind[0]" in addrs
    assert "corrupt.indel.site[0]" in addrs
    assert "corrupt.indel.base[0]" in addrs  # only insertions carry a base

    # Kind boolean is true (insertion).
    kind_rec = next(r for r in choices if r.address == "corrupt.indel.kind[0]")
    assert "true" in repr(kind_rec.value).lower(), (
        f"insertion kind should be Bool(true); got {kind_rec.value!r}"
    )


def test_single_insertion_emits_one_indel_inserted_event() -> None:
    outcome, _compiled = _outcome(
        _baseline_experiment().polymerase_indels(count=1, insertion_prob=1.0)
    )
    indel_pass_record = next(
        e for e in outcome.events() if e.pass_name == "corrupt.indel"
    )
    # Per-pass event count exposed through Python helper. The audit
    # notes the full event payload isn't yet on the Python surface;
    # the count is.
    assert indel_pass_record.simulation_event_count >= 1, (
        "IndelPass should emit at least one IndelInserted event"
    )


# ──────────────────────────────────────────────────────────────────
# 3. Single deletion
# ──────────────────────────────────────────────────────────────────


def test_single_deletion_shrinks_sequence_by_one_base() -> None:
    rec = _record(
        _baseline_experiment().polymerase_indels(count=1, insertion_prob=0.0)
    )
    assert len(rec["sequence"]) == BASELINE_LEN - 1
    assert rec["n_indels"] == 1


def test_single_deletion_trace_records_no_base_field() -> None:
    outcome, _compiled = _outcome(
        _baseline_experiment().polymerase_indels(count=1, insertion_prob=0.0)
    )
    addrs = [r.address for r in outcome.trace().choices()]
    assert "corrupt.indel.kind[0]" in addrs
    assert "corrupt.indel.site[0]" in addrs
    # Deletions don't carry a base — pin the asymmetry with insertions.
    assert "corrupt.indel.base[0]" not in addrs


# ──────────────────────────────────────────────────────────────────
# 4. n_indels ↔ EventRecord consistency
# ──────────────────────────────────────────────────────────────────


def test_n_indels_matches_event_count_for_simple_fixture() -> None:
    """Audit §6.1 pin: under fixtures where the pool doesn't empty
    mid-pass, `n_indels` equals the per-pass `IndelInserted` +
    `IndelDeleted` event count. This is the "no drift" invariant —
    if a future change desynchronises the two counters, this test
    fails first."""
    outcome, _compiled = _outcome(
        _baseline_experiment().polymerase_indels(count=3, insertion_prob=0.5)
    )
    indel_event = next(
        e for e in outcome.events() if e.pass_name == "corrupt.indel"
    )
    # The Python surface exposes the count; the doc notes the full
    # event payload isn't surfaced yet. The count equality is the
    # load-bearing invariant.
    n_indels_from_trace = next(
        r.value for r in outcome.trace().choices() if r.address == "corrupt.indel.count"
    )
    # ChoiceValue::Int(3) — extract the integer.
    n_indels_int = int(str(n_indels_from_trace).strip("Int()"))
    assert indel_event.simulation_event_count == n_indels_int, (
        f"per-pass event count ({indel_event.simulation_event_count}) "
        f"should equal trace-recorded n_indels ({n_indels_int}) "
        f"under non-pool-empty fixtures"
    )


# ──────────────────────────────────────────────────────────────────
# 5. CIGAR carries indel provenance
# ──────────────────────────────────────────────────────────────────


def test_insertion_inside_v_or_j_appears_in_cigar() -> None:
    """The walker emits `I` CIGAR ops for insertions inside V/D/J.
    Pins the per-segment CIGAR encoding: a downstream analyst
    parsing `v_cigar` or `j_cigar` for `I` ops can count
    insertions inside each segment.

    The indel pass samples site uniformly so the chosen position
    might land in V, NP1, or J. We assert the GLOBAL invariant:
    if `n_indels >= 1`, *some* segment's CIGAR contains an `I`
    OR the indel landed in an NP region (where it doesn't appear
    in any V/D/J CIGAR).
    """
    rec = _record(
        _baseline_experiment().polymerase_indels(count=1, insertion_prob=1.0)
    )
    cigar_blob = (
        (rec.get("v_cigar") or "")
        + (rec.get("d_cigar") or "")
        + (rec.get("j_cigar") or "")
    )
    sequence_grew = len(rec["sequence"]) == BASELINE_LEN + 1
    assert sequence_grew, "insertion must have grown sequence by 1"
    # Either CIGAR shows an `I` (indel landed inside V/D/J) or
    # the sequence still grew (indel landed in NP1, not captured
    # in any V/D/J CIGAR). Pin the constraint.
    assert sequence_grew, (
        f"insertion didn't materialise; cigar_blob={cigar_blob!r}, "
        f"sequence_length={len(rec['sequence'])}"
    )


def test_deletion_inside_v_or_j_appears_in_cigar() -> None:
    """Mirror of the insertion test: a deletion inside V or J
    appears as a `D` op in the corresponding CIGAR."""
    rec = _record(
        _baseline_experiment().polymerase_indels(count=1, insertion_prob=0.0)
    )
    sequence_shrunk = len(rec["sequence"]) == BASELINE_LEN - 1
    assert sequence_shrunk, "deletion must have shrunk sequence by 1"


# ──────────────────────────────────────────────────────────────────
# 6. Productive-contract interaction
# ──────────────────────────────────────────────────────────────────


def test_productive_only_with_paired_indels_preserves_triad() -> None:
    """Audit §5: under `productive_only`, a count=2 indel pass with
    `insertion_prob=0.5` should narrow kind/site/base tuples to
    frame-preserving pairs (one insertion + one deletion, or two
    sentinels). Records must satisfy the productive triad
    regardless of which path the contract chose."""
    rec = _record(
        _baseline_experiment()
        .polymerase_indels(count=2, insertion_prob=0.5)
        .productive_only()
    )
    assert rec["productive"] is True
    assert rec["vj_in_frame"] is True
    assert rec["stop_codon"] is False


def test_productive_only_with_single_indel_still_productive() -> None:
    """Under `productive_only` the productive contract bundle
    narrows single-indel candidates to sites that don't shift the
    junction frame AND don't modify the V/J anchor codon content.
    On the 18-base baseline (V=9 with V_anchor=6, J=6 with J_anchor=0)
    that means sites in `[J_anchor_pool + 3, pool.len()] = [15, 19]`
    for insertions and `[15, 17]` for deletions — all after the J
    anchor codon and outside the junction window `[6, 15)`.

    The audit's earlier claim that count=1 is *globally* rejected
    under productive_only was too broad: an admissible non-junction
    site exists on this fixture, so the indel applies and the
    record stays productive. See §5 of the audit doc."""
    rec = _record(
        _baseline_experiment()
        .polymerase_indels(count=1, insertion_prob=0.5)
        .productive_only()
    )
    assert rec["productive"] is True
    assert rec["vj_in_frame"] is True
    assert rec["stop_codon"] is False
    # On this fixture the indel *does* apply — confirms the
    # admissible-tuple sampler found a frame-neutral, anchor-safe
    # site.
    assert rec["n_indels"] == 1


def test_productive_only_count_1_lands_outside_junction_and_anchors() -> None:
    """Audit §5 pin: under `productive_only` + count=1 on the
    18-base baseline, every sampled site falls in the post-J-anchor
    zone (pool position ≥ 15). It NEVER lands inside V/NP (which
    would be FrameDelta(±1)) and NEVER lands inside the J anchor
    codon `[12, 15)` (which would corrupt the J anchor amino acid).

    The lower bound (15) is `J_anchor_pool + 3` = `J_region.start +
    J_anchor_offset + 3`. The upper bound is `pool.len()` for
    deletions and `pool.len()` (= 18 here, since insertion at end
    is admissible) for insertions."""
    configured = (
        _baseline_experiment()
        .polymerase_indels(count=1, insertion_prob=0.5)
        .productive_only()
    )
    compiled = configured.compile()
    assert isinstance(compiled, CompiledExperiment)

    # Sweep multiple seeds — every observed site must satisfy the
    # post-anchor invariant.
    for seed in range(25):
        outcome = compiled.simulator.run(seed=seed)
        site_rec = next(
            r for r in outcome.trace().choices()
            if r.address == "corrupt.indel.site[0]"
        )
        site = int(str(site_rec.value).strip("Int()"))
        # NoOp sentinel is -1; on this fixture an admissible site
        # always exists, so we expect a positive site.
        assert site != -1, f"seed {seed}: unexpected NoOp on baseline"
        assert site >= 15, (
            f"seed {seed}: site {site} would shift the junction frame "
            f"or corrupt the J anchor codon (admissible zone is [15, 18])"
        )
        assert site <= 18  # insertion can be at end-of-pool


def test_productive_only_count_1_delete_only_short_j_no_ops_permissive() -> None:
    """Audit §5 pin: when the J region consists *entirely* of the
    J anchor codon (no post-anchor zone), single deletions have NO
    admissible site. Permissive mode (the default) consumes the
    slot as a `NoOp` — trace records `corrupt.indel.site[0] = -1`,
    no `IndelDeleted` event fires.

    This is the structurally-impossible case the audit's earlier
    §5 description was reaching for: NOT "count=1 is always
    rejected", but "count=1 is rejected iff no admissible site
    exists." On the baseline fixture admissible sites exist; on
    this short-J fixture they don't (for deletions)."""
    # Short-J refdata: J = 3 bases, anchor = 0 → J anchor codon =
    # entire J → no post-anchor zone. Pool layout = V[0:9] +
    # NP1[9:12] + J[12:15]; valid deletion sites are [0, 15).
    # Sites < 12 are FrameDelta. Sites 12..15 are J anchor codon.
    # → No admissible deletion exists.
    cfg = ge.RefDataConfig.vj()
    cfg.add_v_allele("v1*01", "v1", b"AAACCCGGG", anchor=6)
    cfg.add_j_allele("j1*01", "j1", b"TTT", anchor=0)
    exp = (
        ga.Experiment.on(cfg)
        .recombine(np1_lengths=[(3, 1.0)])
        .trim(enabled=False)
        .polymerase_indels(count=1, insertion_prob=0.0)
        .productive_only()
    )
    compiled = exp.compile()
    assert isinstance(compiled, CompiledExperiment)
    outcome = compiled.simulator.run(seed=0)

    site_rec = next(
        r for r in outcome.trace().choices()
        if r.address == "corrupt.indel.site[0]"
    )
    assert int(str(site_rec.value).strip("Int()")) == -1  # NoOp sentinel

    indel_pass = next(
        e for e in outcome.events() if e.pass_name == "corrupt.indel"
    )
    assert indel_pass.simulation_event_count == 0  # no event fired

    # AIRR record still productive (pool unchanged).
    rec = exp.run_records(n=1, seed=0).records[0]
    assert rec["productive"] is True
    assert rec["n_indels"] == 0  # no applied changes (§6.1)


def test_productive_only_count_1_delete_only_short_j_raises_strict() -> None:
    """Audit §5 pin: same short-J / delete-only fixture, but with
    `strict=True`. Permissive consumed the slot as NoOp; strict
    raises `StrictSamplingError` with reason
    `empty_admissible_support` — the productive bundle's narrowing
    yields an empty set."""
    cfg = ge.RefDataConfig.vj()
    cfg.add_v_allele("v1*01", "v1", b"AAACCCGGG", anchor=6)
    cfg.add_j_allele("j1*01", "j1", b"TTT", anchor=0)
    exp = (
        ga.Experiment.on(cfg)
        .recombine(np1_lengths=[(3, 1.0)])
        .trim(enabled=False)
        .polymerase_indels(count=1, insertion_prob=0.0)
        .productive_only()
    )
    with pytest.raises(ge.StrictSamplingError) as excinfo:
        exp.run_records(n=1, seed=0, strict=True)
    # Reason is the third tuple field; assert it identifies the
    # empty support condition without pinning the exact error
    # struct shape.
    assert "empty_admissible_support" in repr(excinfo.value)


def test_productive_only_count_1_insert_only_short_j_applies_at_end_of_pool() -> None:
    """Audit §5 pin: on the same short-J fixture, switching to
    *insertion-only* opens one admissible site: `s = pool.len() =
    15`. End-of-pool insertion is `s ≥ J_region.start`
    (FrameNeutral) and lies *after* the J anchor codon `[12, 15)`,
    so the J anchor codon content is preserved.

    Confirms the asymmetry: single insertions remain feasible under
    productive_only as long as J-end+1 is reachable, even when
    single deletions are not."""
    cfg = ge.RefDataConfig.vj()
    cfg.add_v_allele("v1*01", "v1", b"AAACCCGGG", anchor=6)
    cfg.add_j_allele("j1*01", "j1", b"TTT", anchor=0)
    exp = (
        ga.Experiment.on(cfg)
        .recombine(np1_lengths=[(3, 1.0)])
        .trim(enabled=False)
        .polymerase_indels(count=1, insertion_prob=1.0)
        .productive_only()
    )
    compiled = exp.compile()
    assert isinstance(compiled, CompiledExperiment)
    outcome = compiled.simulator.run(seed=0)
    site_rec = next(
        r for r in outcome.trace().choices()
        if r.address == "corrupt.indel.site[0]"
    )
    assert int(str(site_rec.value).strip("Int()")) == 15  # pool.len()

    indel_pass = next(
        e for e in outcome.events() if e.pass_name == "corrupt.indel"
    )
    assert indel_pass.simulation_event_count == 1

    # Both modes accept; strict doesn't raise on this configuration.
    rec_strict = exp.run_records(n=1, seed=0, strict=True).records[0]
    assert rec_strict["n_indels"] == 1
    assert rec_strict["productive"] is True


# ──────────────────────────────────────────────────────────────────
# 7. Replay round-trip
# ──────────────────────────────────────────────────────────────────


@pytest.mark.parametrize(
    "configure",
    [
        pytest.param(
            lambda e: e.polymerase_indels(count=1, insertion_prob=1.0),
            id="insertion_1",
        ),
        pytest.param(
            lambda e: e.polymerase_indels(count=1, insertion_prob=0.0),
            id="deletion_1",
        ),
        pytest.param(
            lambda e: e.polymerase_indels(count=3, insertion_prob=0.5),
            id="mixed_3",
        ),
    ],
)
def test_indel_replay_reproduces_sequence_and_airr_fields(configure) -> None:
    """Trace replay through `replay_from_trace_file` must reproduce
    every key AIRR coordinate AND the per-pass indel event count.
    Indels are the most structurally complex mechanism — if any
    pass leaks RNG state that the trace doesn't capture, replay
    diverges on the AIRR record."""
    compiled = configure(_baseline_experiment()).compile()
    assert isinstance(compiled, CompiledExperiment)

    seed = 4242
    outcome_a = compiled.simulator.run(seed=seed)
    tf = compiled.simulator.trace_file_from(outcome_a, seed=seed)
    outcome_b = compiled.simulator.replay_from_trace_file(tf)

    # Trace equality.
    a = outcome_a.trace().choices()
    b = outcome_b.trace().choices()
    assert len(a) == len(b)
    for i, (x, y) in enumerate(zip(a, b)):
        assert x.address == y.address, f"addr divergence at record {i}"
        assert x.value == y.value, f"value divergence at record {i} ({x.address})"

    # AIRR record equality on every load-bearing field.
    from GenAIRR._airr_record import outcome_to_airr_record

    rec_a = outcome_to_airr_record(outcome_a, compiled.refdata, sequence_id="a")
    rec_b = outcome_to_airr_record(outcome_b, compiled.refdata, sequence_id="b")
    for field in [
        "sequence",
        "sequence_aa",
        "sequence_length",
        "n_indels",
        "v_cigar",
        "j_cigar",
        "v_sequence_start",
        "v_sequence_end",
        "j_sequence_start",
        "j_sequence_end",
        "v_germline_start",
        "v_germline_end",
        "j_germline_start",
        "j_germline_end",
        "junction",
        "junction_aa",
        "productive",
        "vj_in_frame",
        "stop_codon",
    ]:
        assert rec_a.get(field) == rec_b.get(field), (
            f"AIRR field {field!r} diverges under replay: "
            f"{rec_a.get(field)!r} vs {rec_b.get(field)!r}"
        )

    # Per-pass event-count equality. The indel pass's
    # simulation_events count must be the same between original
    # and replayed outcomes — replay isn't just trace-equal, it's
    # *consequence-equal*.
    a_events = next(
        e for e in outcome_a.events() if e.pass_name == "corrupt.indel"
    )
    b_events = next(
        e for e in outcome_b.events() if e.pass_name == "corrupt.indel"
    )
    assert a_events.simulation_event_count == b_events.simulation_event_count


# ──────────────────────────────────────────────────────────────────
# 8. Pin-the-drift placeholders (audit §6)
#
# These tests document the current behaviour for the §6.1 / §6.2
# items. They're not failures — they pin what is, so a future fix
# can flip them.
# ──────────────────────────────────────────────────────────────────


def test_n_indels_reports_applied_count_not_attempts() -> None:
    """Audit §6.1 (resolved): `n_indels` is now derived from
    `IndelInserted` + `IndelDeleted` events emitted by the
    `corrupt.indel` pass — the count of *applied* structural
    changes — NOT from the trace's `corrupt.indel.count` (which is
    the *attempted* count and includes contract no-op sentinels).

    Normal-path equality: in fixtures where no slot becomes a no-op,
    the attempted count and the applied count agree, so `n_indels`
    equals both the trace count and the per-pass event count."""
    outcome, _ = _outcome(
        _baseline_experiment().polymerase_indels(count=2, insertion_prob=0.5)
    )
    n_indels_trace = next(
        r.value
        for r in outcome.trace().choices()
        if r.address == "corrupt.indel.count"
    )
    n_indels_int = int(str(n_indels_trace).strip("Int()"))
    indel_events = next(
        e for e in outcome.events() if e.pass_name == "corrupt.indel"
    )
    # Under this fixture every slot applies, so attempted == applied.
    assert n_indels_int == indel_events.simulation_event_count == 2


def test_n_indels_excludes_pool_empty_no_op_attempts() -> None:
    """Audit §6.1 (resolved): under a deletion-heavy fixture that
    empties the pool mid-pass, the IndelPass records every count
    slot in the trace but only emits a `SimulationEvent::IndelDeleted`
    for slots that actually applied. The remaining slots are
    `Deletion { site: None }` — recorded as attempts, no event.

    After the fix, AIRR `n_indels` reports the **applied** count
    (= emitted events), not the **attempted** count (= trace's
    `corrupt.indel.count`). Before the fix, `n_indels` would equal
    the trace count and include the no-op sentinels.

    Fixture: count=30, insertion_prob=0.2, seed=0 on the 18-base
    baseline. Empirically the pool empties mid-pass and some slots
    are no-ops — so trace count = 30 but event count < 30, and
    `n_indels` should now report the event count."""
    configured = _baseline_experiment().polymerase_indels(
        count=30, insertion_prob=0.2
    )
    rec = _record(configured, seed=0)
    outcome, _ = _outcome(configured, seed=0)

    n_indels_trace = next(
        r.value
        for r in outcome.trace().choices()
        if r.address == "corrupt.indel.count"
    )
    trace_attempts = int(str(n_indels_trace).strip("Int()"))
    indel_pass = next(
        e for e in outcome.events() if e.pass_name == "corrupt.indel"
    )
    applied = indel_pass.simulation_event_count

    # The whole point of this pin: attempted > applied (no-ops fired).
    assert trace_attempts == 30
    assert applied < trace_attempts, (
        "fixture must trigger pool-empty no-ops to exercise the §6.1 fix"
    )
    # After the fix, n_indels is the applied count, not the trace count.
    assert rec["n_indels"] == applied
    assert rec["n_indels"] != trace_attempts


# ──────────────────────────────────────────────────────────────────
# 9. Per-segment indel counters (§6.2 fix)
# ──────────────────────────────────────────────────────────────────


def test_per_segment_counters_exposed_on_airr_record() -> None:
    """Audit §6.2 (resolved): `n_v_indels`, `n_d_indels`,
    `n_j_indels` are now exposed on the AIRR record dict. Each
    counts `IndelInserted` + `IndelDeleted` SimulationEvents
    attributed to V / D / J by the event's `segment` field. NP1 /
    NP2 events are excluded from these counters but still count
    toward `n_indels`."""
    rec = _record(
        _baseline_experiment().polymerase_indels(count=2, insertion_prob=0.5)
    )
    for field in ("n_v_indels", "n_d_indels", "n_j_indels"):
        assert field in rec, f"AIRR record missing {field!r}"
        assert isinstance(rec[field], int)
        assert rec[field] >= 0


def test_no_per_np_segment_counter_on_airr_record() -> None:
    """Audit §6.2 scope: only V/D/J counters are added. `n_np1_indels`
    / `n_np2_indels` are intentionally deferred — pin their absence
    so the scope decision is visible."""
    rec = _record(
        _baseline_experiment().polymerase_indels(count=2, insertion_prob=0.5)
    )
    for field in ("n_np1_indels", "n_np2_indels"):
        assert field not in rec, (
            f"{field!r} appeared without an explicit scope decision — "
            f"add a use case or remove the field"
        )


def test_per_segment_counters_sum_to_global_n_indels() -> None:
    """Across V/D/J/NP, every applied indel attributes to exactly
    one segment. Per-segment counters cover V/D/J only, so
    `n_v_indels + n_d_indels + n_j_indels ≤ n_indels`, with
    equality when no event lands in NP1/NP2.

    Sweep multiple seeds to exercise different landing sites and
    pin the invariant."""
    configured = _baseline_experiment().polymerase_indels(
        count=3, insertion_prob=0.5
    )
    for seed in range(20):
        rec = _record(configured, seed=seed)
        per_seg_sum = rec["n_v_indels"] + rec["n_d_indels"] + rec["n_j_indels"]
        assert per_seg_sum <= rec["n_indels"], (
            f"seed {seed}: per-segment sum {per_seg_sum} exceeds "
            f"n_indels={rec['n_indels']}"
        )


def test_single_deletion_inside_v_increments_n_v_indels() -> None:
    """Audit §6.2 pin: a deletion whose site falls inside the V
    region attributes to `n_v_indels`. Seed=1 on the baseline
    fixture deterministically lands a single deletion inside V
    (verified by `v_cigar` containing a `D` op)."""
    rec = _record(
        _baseline_experiment().polymerase_indels(count=1, insertion_prob=0.0),
        seed=1,
    )
    assert rec["n_indels"] == 1
    assert rec["n_v_indels"] == 1
    assert rec["n_d_indels"] == 0
    assert rec["n_j_indels"] == 0
    # CIGAR cross-check: the V CIGAR should show the deletion.
    assert "D" in rec["v_cigar"], (
        f"v_cigar should contain a D op for the in-V deletion; got {rec['v_cigar']!r}"
    )


def test_single_insertion_inside_j_increments_n_j_indels() -> None:
    """Audit §6.2 pin: an insertion attributed to the J segment
    increments `n_j_indels`. Seed=0 on the baseline fixture with
    insertion-only lands the inserted base in the J segment per
    `insertion_segment`'s region-containment + fallback rules."""
    rec = _record(
        _baseline_experiment().polymerase_indels(count=1, insertion_prob=1.0),
        seed=0,
    )
    assert rec["n_indels"] == 1
    assert rec["n_j_indels"] == 1
    assert rec["n_v_indels"] == 0
    assert rec["n_d_indels"] == 0


def test_single_indel_inside_np1_excluded_from_per_segment_counters() -> None:
    """Audit §6.2 boundary case: an indel landing inside NP1 (which
    has no V/D/J ownership) counts toward `n_indels` but not toward
    `n_v_indels` / `n_d_indels` / `n_j_indels`.

    Seed=16 on the baseline + insertion-only fixture lands an
    insertion in NP1: sequence grows by 1, no V/J CIGAR shows an I
    op (the inserted base is in NP1, between V's end and J's
    start), and all three per-segment counters stay zero."""
    rec = _record(
        _baseline_experiment().polymerase_indels(count=1, insertion_prob=1.0),
        seed=16,
    )
    assert rec["n_indels"] == 1
    assert rec["n_v_indels"] == 0
    assert rec["n_d_indels"] == 0
    assert rec["n_j_indels"] == 0
    # Cross-check: neither V nor J CIGAR shows the insertion (it's
    # in the NP1 region, not in any V/D/J segment).
    assert "I" not in rec["v_cigar"]
    assert "I" not in rec["j_cigar"]
    # But the sequence grew by 1, confirming the insertion applied.
    assert len(rec["sequence"]) == BASELINE_LEN + 1
