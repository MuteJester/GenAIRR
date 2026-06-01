"""Primer-trim / end-loss provenance golden tests.

Companion to [docs/primer_trim_end_loss_audit.md](../docs/primer_trim_end_loss_audit.md).
This file pins **current** behaviour — including the
provenance-gap drift items the audit catalogues. When a future
slice fixes the drift, the corresponding assertions here flip from
"documents the gap" to "asserts the fix." That flip is itself a
useful signal that the fix landed on the right surface.

What's pinned:

- 5'-loss / 3'-loss shorten the sequence by exactly the requested
  length under deterministic NP1 / trim fixtures.
- Recombination-trim AIRR fields (`v_trim_5/3`, `j_trim_5/3`) are
  untouched by end-loss — the engine preserves biological
  provenance at the trace/event layer.
- The trace carries `corrupt.end_loss.5` / `corrupt.end_loss.3`
  records with the realized count.
- `n_indels` is unaffected — end-loss isn't polymerase-indel
  biology.
- Replay through `replay_from_trace_file` reproduces every key
  AIRR coordinate.
- Pin-the-drift: there is no `end_loss_5_length` /
  `end_loss_3_length` / `n_end_loss_*` field on the AIRR record
  today. When one is added, the relevant test below will fail and
  signal the fix.

Fixture design: VJ refdata with V=9 (anchor at pos 6), J=6 (anchor
at pos 0), trim disabled, NP1 fixed at length 3 → final assembled
sequence is exactly 9+3+6=18 bases. Every length-math assertion
keys off that 18.
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
    cfg.add_v_allele("v1*01", "v1", b"AAACCCGGG", anchor=6)  # 9 bases, anchor codon at 6..9
    cfg.add_j_allele("j1*01", "j1", b"TTTAAA", anchor=0)     # 6 bases, anchor codon at 0..3
    return cfg


# Total assembled sequence length when trim is disabled and NP1
# is pinned at 3: V(9) + NP1(3) + J(6) = 18.
BASELINE_LEN = 18


def _baseline_experiment() -> "ga.Experiment":
    """Build the deterministic recombine-only fixture: no trim,
    NP1 fixed at 3 bases. Every test reuses this as its starting
    point so length math is fully predictable."""
    return (
        ga.Experiment.on(_vj_refdata()).allow_curatable_refdata()
        .recombine(np1_lengths=[(3, 1.0)])
        .trim(enabled=False)
    )


def _record(exp: "ga.Experiment", seed: int = 0) -> dict:
    """Run `exp` once and return the single AIRR record."""
    result = exp.run_records(n=1, seed=seed)
    assert len(result.records) == 1
    return result.records[0]


# ──────────────────────────────────────────────────────────────────
# 1. Baseline — confirm fixture predictability
# ──────────────────────────────────────────────────────────────────


def test_baseline_fixture_produces_18_base_sequence_with_zero_trims() -> None:
    rec = _record(_baseline_experiment())
    assert len(rec["sequence"]) == BASELINE_LEN
    # Recombination trim disabled → all trim fields are 0.
    assert rec["v_trim_5"] == 0
    assert rec["v_trim_3"] == 0
    assert rec["j_trim_5"] == 0
    assert rec["j_trim_3"] == 0
    # No indel pass in the plan → n_indels = 0.
    assert rec["n_indels"] == 0


# ──────────────────────────────────────────────────────────────────
# 2. 5' end-loss / primer-trim
# ──────────────────────────────────────────────────────────────────


def test_5prime_loss_shortens_sequence_by_requested_length() -> None:
    rec = _record(_baseline_experiment().primer_trim_5prime(length=2))
    assert len(rec["sequence"]) == BASELINE_LEN - 2


def test_5prime_loss_does_not_update_v_trim_5() -> None:
    # The audit's load-bearing semantic: end-loss is observation-
    # stage and does NOT rewrite the recombination-stage allele
    # trim metadata. The two provenances stay distinct in the
    # AIRR record's `v_trim_5` field.
    rec = _record(_baseline_experiment().primer_trim_5prime(length=2))
    assert rec["v_trim_5"] == 0, (
        f"v_trim_5 should reflect recombination trim only, not end-loss. "
        f"Got v_trim_5={rec['v_trim_5']} after primer_trim_5prime(length=2)."
    )


def test_5prime_loss_records_corrupt_end_loss_5_in_trace() -> None:
    # The trace carries the realized end-loss count. This is the
    # provenance signal at the trace layer — the audit notes it
    # never reaches the projected AIRR record (drift item §6.1).
    compiled = (
        _baseline_experiment().primer_trim_5prime(length=2).compile(allow_curatable_refdata=True)
    )
    assert isinstance(compiled, CompiledExperiment)
    outcome = compiled.simulator.run(seed=0)
    end_loss_records = [
        r for r in outcome.trace().choices() if r.address == "corrupt.end_loss.5"
    ]
    assert len(end_loss_records) == 1
    # The recorded value is the *realized* count (clamped to pool
    # length, etc.). With length=2 and an 18-base pool, no
    # clamping happens.
    rec = end_loss_records[0]
    # ChoiceValue::Int(2) — the realized 5'-loss count.
    assert "2" in repr(rec.value)


# ──────────────────────────────────────────────────────────────────
# 3. 3' end-loss / primer-trim — symmetric assertions
# ──────────────────────────────────────────────────────────────────


def test_3prime_loss_shortens_sequence_by_requested_length() -> None:
    rec = _record(_baseline_experiment().primer_trim_3prime(length=2))
    assert len(rec["sequence"]) == BASELINE_LEN - 2


def test_3prime_loss_does_not_update_j_trim_3() -> None:
    rec = _record(_baseline_experiment().primer_trim_3prime(length=2))
    assert rec["j_trim_3"] == 0, (
        f"j_trim_3 should reflect recombination trim only, not end-loss. "
        f"Got j_trim_3={rec['j_trim_3']} after primer_trim_3prime(length=2)."
    )


def test_3prime_loss_records_corrupt_end_loss_3_in_trace() -> None:
    compiled = (
        _baseline_experiment().primer_trim_3prime(length=2).compile(allow_curatable_refdata=True)
    )
    assert isinstance(compiled, CompiledExperiment)
    outcome = compiled.simulator.run(seed=0)
    end_loss_records = [
        r for r in outcome.trace().choices() if r.address == "corrupt.end_loss.3"
    ]
    assert len(end_loss_records) == 1
    assert "2" in repr(end_loss_records[0].value)


# ──────────────────────────────────────────────────────────────────
# 4. Combined 5' + 3' loss
# ──────────────────────────────────────────────────────────────────


def test_combined_5prime_and_3prime_loss_shortens_by_sum() -> None:
    rec = _record(
        _baseline_experiment()
        .primer_trim_5prime(length=2)
        .primer_trim_3prime(length=1)
    )
    assert len(rec["sequence"]) == BASELINE_LEN - 3
    # Neither side touches recombination trim.
    assert rec["v_trim_5"] == 0
    assert rec["j_trim_3"] == 0


# ──────────────────────────────────────────────────────────────────
# 5. `n_indels` is a separate provenance (audit §6.4)
# ──────────────────────────────────────────────────────────────────


def test_end_loss_does_not_increment_n_indels() -> None:
    # `n_indels` counts polymerase-indel events (`IndelPass`),
    # not end-loss / primer-trim losses. The audit recommends a
    # dedicated counter for end-loss (drift §6.4); until that
    # lands, `n_indels` must remain end-loss-blind.
    rec = _record(
        _baseline_experiment()
        .primer_trim_5prime(length=3)
        .primer_trim_3prime(length=2)
    )
    assert rec["n_indels"] == 0, (
        f"n_indels should not include end-loss / primer-trim. "
        f"Got n_indels={rec['n_indels']} after 5+2 bases of primer trim."
    )


# ──────────────────────────────────────────────────────────────────
# 6. Replay round-trip
# ──────────────────────────────────────────────────────────────────


@pytest.mark.parametrize(
    "configure",
    [
        pytest.param(lambda e: e.primer_trim_5prime(length=2), id="5prime"),
        pytest.param(lambda e: e.primer_trim_3prime(length=2), id="3prime"),
        pytest.param(
            lambda e: e.primer_trim_5prime(length=2).primer_trim_3prime(length=1),
            id="both",
        ),
    ],
)
def test_primer_trim_replay_reproduces_sequence_and_coords(configure) -> None:
    """Trace replay through `replay_from_trace_file` must reproduce
    every key AIRR coordinate for an end-loss experiment.
    Critical because end-loss is a structural mechanism that shifts
    every downstream region's pool coordinate."""
    compiled = configure(_baseline_experiment()).compile(allow_curatable_refdata=True)
    assert isinstance(compiled, CompiledExperiment)

    seed = 4242
    outcome_a = compiled.simulator.run(seed=seed)
    tf = compiled.simulator.trace_file_from(outcome_a, seed=seed)
    outcome_b = compiled.simulator.replay_from_trace_file(tf)

    # Trace records match record-for-record.
    a_records = outcome_a.trace().choices()
    b_records = outcome_b.trace().choices()
    assert len(a_records) == len(b_records)
    for i, (a, b) in enumerate(zip(a_records, b_records)):
        assert a.address == b.address, f"addr divergence at record {i}"
        assert a.value == b.value, f"value divergence at record {i} ({a.address})"

    # AIRR records match on every load-bearing coordinate.
    from GenAIRR._airr_record import outcome_to_airr_record

    record_a = outcome_to_airr_record(outcome_a, compiled.refdata, sequence_id="a")
    record_b = outcome_to_airr_record(outcome_b, compiled.refdata, sequence_id="b")
    for field in [
        "sequence",
        "sequence_aa",
        "v_sequence_start",
        "v_sequence_end",
        "j_sequence_start",
        "j_sequence_end",
        "v_germline_start",
        "v_germline_end",
        "j_germline_start",
        "j_germline_end",
        "v_trim_5",
        "v_trim_3",
        "j_trim_5",
        "j_trim_3",
        "junction",
        "junction_aa",
        "productive",
        "vj_in_frame",
        "stop_codon",
    ]:
        assert record_a.get(field) == record_b.get(field), (
            f"AIRR field {field!r} diverges under replay: "
            f"{record_a.get(field)!r} vs {record_b.get(field)!r}"
        )


# ──────────────────────────────────────────────────────────────────
# 7. Productive_only interaction (audit §5 — already correct)
# ──────────────────────────────────────────────────────────────────


def test_productive_only_keeps_records_productive_after_5prime_loss() -> None:
    # 5'-loss of 6 bytes removes V[0..6] but stops at the V-anchor
    # codon (positions 6..9). Under productive_only the
    # AnchorPreserved contract ensures the loss never bites into
    # the anchor — the productive triad must still hold.
    rec = _record(
        _baseline_experiment()
        .primer_trim_5prime(length=6)
        .productive_only()
    )
    assert rec["productive"] is True
    assert rec["vj_in_frame"] is True
    assert rec["stop_codon"] is False


def test_productive_only_keeps_records_productive_after_3prime_loss() -> None:
    # 3'-loss biting into J would corrupt the J anchor codon
    # (positions 0..3). productive_only must keep the loss at or
    # below the safe bound. Test asserts the end-to-end productive
    # invariant holds regardless of clamping behaviour.
    rec = _record(
        _baseline_experiment()
        .primer_trim_3prime(length=3)
        .productive_only()
    )
    assert rec["productive"] is True
    assert rec["vj_in_frame"] is True
    assert rec["stop_codon"] is False


# ──────────────────────────────────────────────────────────────────
# 8. End-loss provenance exposed on the AIRR record
#
# Closes audit drift §6.1 / §6.4: the realized end-loss byte count
# is now a first-class field on the projected AIRR record, mirroring
# how `n_pcr_errors` / `n_quality_errors` / `n_indels` surface their
# respective per-pass counts. The two fields are populated from the
# `corrupt.end_loss.{5,3}` trace records via `_trace_int`, defaulting
# to 0 when no end-loss pass ran.
# ──────────────────────────────────────────────────────────────────


def test_end_loss_5_length_field_equals_trace_value() -> None:
    """Audit §6.1 fix: the realized 5'-loss count flows into the
    AIRR record as `end_loss_5_length`. Pin: field exists, equals
    the requested length when no clamping fires, and recombination
    trim metadata is still independent provenance."""
    rec = _record(_baseline_experiment().primer_trim_5prime(length=2))
    assert "end_loss_5_length" in rec
    assert rec["end_loss_5_length"] == 2
    # Recombination trim provenance still distinct.
    assert rec["v_trim_5"] == 0
    # 3'-only counter stays at default.
    assert rec["end_loss_3_length"] == 0


def test_end_loss_3_length_field_equals_trace_value() -> None:
    """Audit §6.1 fix mirror for the 3' end."""
    rec = _record(_baseline_experiment().primer_trim_3prime(length=2))
    assert "end_loss_3_length" in rec
    assert rec["end_loss_3_length"] == 2
    assert rec["j_trim_3"] == 0
    assert rec["end_loss_5_length"] == 0


def test_end_loss_fields_default_to_zero_when_no_end_loss_pass_ran() -> None:
    """Baseline experiment with no `primer_trim_*` step: both
    `end_loss_*_length` fields default to 0 (mirrors how
    `n_indels` defaults to 0 when no indel pass ran)."""
    rec = _record(_baseline_experiment())
    assert rec["end_loss_5_length"] == 0
    assert rec["end_loss_3_length"] == 0


def test_combined_end_loss_fields_match_per_end_trace_values() -> None:
    """Both ends populated independently — composes additively in
    sequence length, but the two counters stay separate."""
    rec = _record(
        _baseline_experiment()
        .primer_trim_5prime(length=2)
        .primer_trim_3prime(length=1)
    )
    assert rec["end_loss_5_length"] == 2
    assert rec["end_loss_3_length"] == 1
    # Sequence-length math still holds.
    assert len(rec["sequence"]) == BASELINE_LEN - 3
