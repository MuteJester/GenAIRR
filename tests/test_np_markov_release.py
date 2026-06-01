"""Release-tier composition tests for the Markov NP base generator.

Confirms that the Markov NP base sampling slice composes with
the rest of the release stack (productive IGH, SHM, contracts,
trace replay, AIRR projection + validation) and preserves the
two load-bearing invariants from the audit:

- **Records validate clean.** Productive IGH with a Markov NP1
  cartridge produces records whose AIRR projection passes
  :py:meth:`SimulationResult.validate_records` end-to-end.
- **Replay round-trip is byte-identical.** The recorded NP1
  bytes, the full ``sequence``, and the ``junction`` reproduce
  exactly under ``rerun_from_trace_file``; no new addresses,
  no Markov-specific replay quirks.
- **Dependency invariant.** A deterministic Markov fixture
  (first base = G; G→A, A→T, T→C, C→G) produces the expected
  walk pattern across positions, which a position-independent
  implementation cannot reproduce.

Companion to
[`docs/np_markov_base_generator_design.md`](../docs/np_markov_base_generator_design.md)
+ [`docs/junction_n_addition_audit.md`](../docs/junction_n_addition_audit.md).
"""
from __future__ import annotations

import copy
import json

import pytest

import GenAIRR as ga
from GenAIRR._airr_record import outcome_to_airr_record
from GenAIRR.reference_models import (
    NpBaseModelSpec,
    ReferenceEmpiricalModels,
)


# ──────────────────────────────────────────────────────────────────
# Shared fixtures
# ──────────────────────────────────────────────────────────────────


# Deterministic transition graph from the user brief — first base
# G, G→A, A→T, T→C, C→G — so a non-zero NP1 produces the cycle
# starting `G`, then `G A T C G A T C …`.
_RELEASE_MARKOV = NpBaseModelSpec(
    kind="markov",
    first_base={"G": 1.0},
    transitions={
        "A": {"T": 1.0},
        "C": {"G": 1.0},
        "G": {"A": 1.0},
        "T": {"C": 1.0},
    },
)


def _cartridge_with_markov_np1() -> "ga.DataConfig":
    cfg = copy.deepcopy(ga.HUMAN_IGH_OGRDB)
    existing = getattr(cfg, "reference_models", None)
    if isinstance(existing, ReferenceEmpiricalModels):
        cfg.reference_models = ReferenceEmpiricalModels(
            np_lengths=existing.np_lengths,
            trims=existing.trims,
            np_bases={"NP1": _RELEASE_MARKOV},
        )
    else:
        cfg.reference_models = ReferenceEmpiricalModels(
            np_bases={"NP1": _RELEASE_MARKOV}
        )
    return cfg


def _np1_bytes_from_trace(compiled, *, seed: int) -> str:
    outcome = compiled.simulator.run(seed=seed)
    tf = compiled.simulator.trace_file_from(outcome, seed=seed)
    data = json.loads(tf.to_json())
    return "".join(
        ev["value"]["Base"]
        for ev in data.get("trace", {}).get("choices", [])
        if ev.get("address", "").startswith("np.np1.bases[")
    )


# ──────────────────────────────────────────────────────────────────
# 1. Release composition — productive IGH + Markov NP1 validates
# ──────────────────────────────────────────────────────────────────


def test_productive_igh_full_stack_with_markov_np1_validates_all_records() -> None:
    """Full release composition — productive IGH with a Markov
    NP1 base model, SHM, and the standard end-of-pipeline
    ``validate_records`` postcondition gate. Demonstrates that:

    - The bridge dispatch (`push_generate_np(... markov_transitions=...)`)
      composes through ``Experiment.recombine().productive_only().mutate()``.
    - The admit-mask × Markov-support intersection preserves the
      productive triad (no stop codons in junction, junction
      length divisible by 3, anchor codons intact).
    - The AIRR projection + validator agree on the engine's
      own truth oracle under the new generator surface.
    """
    cfg = _cartridge_with_markov_np1()
    exp = (
        ga.Experiment.on(cfg)
        .recombine()
        .productive_only()
        .mutate(rate=0.005)
    )
    refdata = exp.refdata
    result = exp.run_records(n=50, seed=4242)
    report = result.validate_records(refdata)
    assert report, (
        f"productive IGH + Markov NP1 full-stack validation failed: "
        f"{report.summary()}"
    )
    # And every record actually carries an NP1 region — the
    # cartridge's NP1 length distribution is non-degenerate.
    np1_lengths = [len(r.get("np1") or "") for r in result.records]
    assert any(l > 0 for l in np1_lengths), (
        "no record carried an NP1 segment; the Markov fixture isn't "
        "exercising the per-position sampler"
    )


# ──────────────────────────────────────────────────────────────────
# 2. Replay round-trip — sequence + junction + NP bases preserved
# ──────────────────────────────────────────────────────────────────


def test_markov_release_replay_round_trip_preserves_sequence_junction_np_bases() -> None:
    """Replay determinism across several seeds — for each fresh
    outcome build the AIRR record, replay via
    ``rerun_from_trace_file``, rebuild the AIRR record, and
    assert the load-bearing fields (``sequence``, ``junction``,
    ``np1``) match byte-for-byte. NP1 bytes are sampled by the
    Markov generator with `previous`-base conditioning; the
    replay validator reconstructs `previous` from the prior
    recorded base, so byte-identical replay is the central
    invariant pinned here."""
    cfg = _cartridge_with_markov_np1()
    exp = ga.Experiment.on(cfg).recombine().productive_only()
    refdata = exp.refdata
    compiled = exp.compile()
    matched_any_np1 = False
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
        for field in ("sequence", "junction", "np1"):
            assert fresh_rec[field] == replayed_rec[field], (
                f"seed {seed}: Markov NP1 round-trip desynced on "
                f"field {field!r}: fresh={fresh_rec[field]!r}, "
                f"replay={replayed_rec[field]!r}"
            )
        if fresh_rec["np1"]:
            matched_any_np1 = True
    assert matched_any_np1, (
        "no seed in [0..8) produced a non-empty NP1; the test isn't "
        "actually exercising Markov sampling"
    )


# ──────────────────────────────────────────────────────────────────
# 3. Dependency invariant — deterministic Markov walk
# ──────────────────────────────────────────────────────────────────


def test_deterministic_markov_transitions_produce_expected_walk() -> None:
    """The deterministic fixture (first=G, G→A, A→T, T→C, C→G)
    produces the cyclic byte sequence `G A T C G A T C …`
    starting at position 0. A position-independent generator
    that ignored `previous` could not produce this pattern; a
    Markov generator that mis-keyed the row order (e.g. by
    transposing from/to) would produce a different cycle.
    Either failure mode trips here.
    """
    cfg = _cartridge_with_markov_np1()
    compiled = ga.Experiment.on(cfg).recombine().compile()
    expected_cycle = "GATC"
    # Sweep a handful of seeds until we hit a non-zero NP1.
    np1 = ""
    for seed in (42, 7, 99, 12345, 4242, 1, 11, 113):
        candidate = _np1_bytes_from_trace(compiled, seed=seed)
        if candidate:
            np1 = candidate
            break
    assert np1, "no seed produced a non-zero NP1 — fixture isn't exercising sampler"
    for i, base in enumerate(np1):
        expected = expected_cycle[i % 4]
        assert base == expected, (
            f"Markov walk diverged at position {i}: got {base!r}, "
            f"expected {expected!r} (cycle G→A→T→C; full walk: {np1!r})"
        )
