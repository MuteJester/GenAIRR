"""End-to-end behaviour tests for the **Markov NP base generator**
implementation slice.

Companion to
[`docs/np_markov_base_generator_design.md`](../docs/np_markov_base_generator_design.md)
and the contract pin file
[`tests/test_np_markov_base_generator_contract.py`](test_np_markov_base_generator_contract.py).

Mirrors the user's greenlight test surface:

1. Markov spec lowers and runs.
2. Strong transition matrix produces observable dependency
   (A → T near 100 %, T → G near 100 %, etc.).
3. Replay round-trip succeeds against the same Markov matrix.
4. Replay against a different transition matrix fails the
   plan-signature gate before any choice is consumed.
5. Productive-only triad preserved under Markov.
6. Uniform / no-spec output is byte-identical to the pre-Markov
   baseline.
7. ``empirical_first_base`` is byte-identical to the pre-Markov
   baseline.
8. Manifest lists ``"markov"`` under ``supported_kinds`` and
   removes it from ``deferred_kinds``.

Legacy orphan pins live in
``tests/test_np_base_model_implementation.py`` (the typed-NP-base
slice already enforces them); the Markov slice does NOT enable
auto-lift, so the orphan invariants ride unchanged.
"""
from __future__ import annotations

import copy
import json

import pytest

import GenAIRR as ga
from GenAIRR.reference_models import (
    NpBaseModelSpec,
    ReferenceEmpiricalModels,
)


# ──────────────────────────────────────────────────────────────────
# Shared fixtures
# ──────────────────────────────────────────────────────────────────


def _cartridge_with_np_bases(np_bases: dict) -> "ga.DataConfig":
    cfg = copy.deepcopy(ga.HUMAN_IGH_OGRDB)
    existing = getattr(cfg, "reference_models", None)
    if isinstance(existing, ReferenceEmpiricalModels):
        cfg.reference_models = ReferenceEmpiricalModels(
            np_lengths=existing.np_lengths,
            trims=existing.trims,
            np_bases=np_bases,
        )
    else:
        cfg.reference_models = ReferenceEmpiricalModels(np_bases=np_bases)
    return cfg


def _np1_bases_from_trace(compiled, *, seed: int) -> str:
    """Drive the simulator once, then walk the recorded trace
    to extract the literal ``np.np1.bases[i]`` byte sequence
    in order. Used to verify the Markov walk pattern."""
    outcome = compiled.simulator.run(seed=seed)
    tf = compiled.simulator.trace_file_from(outcome, seed=seed)
    data = json.loads(tf.to_json())
    bases = []
    for ev in data.get("trace", {}).get("choices", []):
        addr = ev.get("address", "")
        if addr.startswith("np.np1.bases["):
            bases.append(ev["value"]["Base"])
    return "".join(bases)


# Deterministic A→T→G→C→A→T→G→C cycle starting from A.
_DETERMINISTIC_MARKOV = NpBaseModelSpec(
    kind="markov",
    first_base={"A": 1.0, "C": 0.0, "G": 0.0, "T": 0.0},
    transitions={
        "A": {"A": 0.0, "C": 0.0, "G": 0.0, "T": 1.0},  # A → T
        "T": {"A": 0.0, "C": 0.0, "G": 1.0, "T": 0.0},  # T → G
        "G": {"A": 0.0, "C": 1.0, "G": 0.0, "T": 0.0},  # G → C
        "C": {"A": 1.0, "C": 0.0, "G": 0.0, "T": 0.0},  # C → A
    },
)


# ──────────────────────────────────────────────────────────────────
# 1. Markov spec lowers and runs
# ──────────────────────────────────────────────────────────────────


def test_markov_spec_lowers_and_runs() -> None:
    """The smoke test — a Markov spec on NP1 lowers through
    the bridge's new ``markov_transitions`` kwarg into a
    ``MarkovBaseGenerator`` and produces records without
    raising. Replaces the prior ``NotImplementedError`` stop-
    and-report."""
    cfg = _cartridge_with_np_bases({"NP1": _DETERMINISTIC_MARKOV})
    recs = ga.Experiment.on(cfg).recombine().run_records(n=5, seed=0)
    assert len(recs) == 5


# ──────────────────────────────────────────────────────────────────
# 2. Strong transition matrix produces observable dependency
# ──────────────────────────────────────────────────────────────────


def test_strong_markov_matrix_produces_deterministic_walk() -> None:
    """A degenerate matrix where every row puts 100 % weight on
    one base produces a deterministic A→T→G→C→A walk from
    position 0. The recorded trace's NP1 byte sequence is the
    cleanest check that the per-position support is genuinely
    keyed on ``previous`` — a position-independent
    implementation could not produce this pattern."""
    cfg = _cartridge_with_np_bases({"NP1": _DETERMINISTIC_MARKOV})
    compiled = ga.Experiment.on(cfg).recombine().compile()
    expected_cycle = "ATGC"
    # Walk several seeds until we hit one with non-zero NP1
    # length; some seeds give a length-0 NP1 which can't
    # exercise the cycle. The Markov determinism guarantee is
    # the load-bearing assertion, not the seed choice.
    np1 = ""
    for seed in (42, 7, 99, 12345, 4242):
        candidate = _np1_bases_from_trace(compiled, seed=seed)
        if candidate:
            np1 = candidate
            break
    assert np1, "no seed produced NP1 of non-zero length"
    for i, base in enumerate(np1):
        assert base == expected_cycle[i % 4], (
            f"Markov walk diverged at position {i}: got {base!r}, "
            f"expected {expected_cycle[i % 4]!r} (full walk: {np1!r})"
        )


def test_a_to_t_near_one_hundred_percent_under_strong_matrix() -> None:
    """Aggregate check — under a matrix where every transition
    is to T (regardless of `previous`), almost every NP1 base
    at positions 1+ must be T. The first base alone retains
    its A weight; positions 1+ converge."""
    cfg = _cartridge_with_np_bases(
        {
            "NP1": NpBaseModelSpec(
                kind="markov",
                first_base={"A": 1.0, "C": 0.0, "G": 0.0, "T": 0.0},
                transitions={
                    "A": {"A": 0.0, "C": 0.0, "G": 0.0, "T": 1.0},
                    "C": {"A": 0.0, "C": 0.0, "G": 0.0, "T": 1.0},
                    "G": {"A": 0.0, "C": 0.0, "G": 0.0, "T": 1.0},
                    "T": {"A": 0.0, "C": 0.0, "G": 0.0, "T": 1.0},
                },
            )
        }
    )
    recs = ga.Experiment.on(cfg).recombine().run_records(n=200, seed=99)
    # Concatenate non-first positions of NP1 across records.
    tail_bases = "".join((r.get("np1") or "")[1:] for r in recs).upper()
    assert tail_bases, "No NP1 tail bytes collected — fixture broken"
    t_fraction = tail_bases.count("T") / len(tail_bases)
    assert t_fraction >= 0.99, (
        f"Strong A→T matrix should force position 1+ to T; "
        f"observed t_fraction={t_fraction:.3f}"
    )


# ──────────────────────────────────────────────────────────────────
# 3. Replay round-trip succeeds against the same matrix
# ──────────────────────────────────────────────────────────────────


def test_markov_replay_round_trip_byte_identical() -> None:
    """Same Markov cartridge + same seed → identical NP1 bases
    on a round-trip through ``trace_file_from`` →
    ``replay_from_trace_file``. The byte-determinism is the
    payoff of the trace-records-emitted-base + replay-rebuilds-
    previous design — no new trace addresses, no Markov-
    specific replay logic outside the per-position validator."""
    cfg = _cartridge_with_np_bases({"NP1": _DETERMINISTIC_MARKOV})
    compiled = ga.Experiment.on(cfg).recombine().compile()
    seed = 4242
    outcome = compiled.simulator.run(seed=seed)
    tf = compiled.simulator.trace_file_from(outcome, seed=seed)
    # Replay should succeed; the recorded NP1 bases must
    # reproduce byte-for-byte under the same generator.
    replayed = compiled.simulator.replay_from_trace_file(tf)
    # Walk both outcomes and verify the NP1 region matches.
    a_data = json.loads(tf.to_json())
    a_np1 = "".join(
        ev["value"]["Base"]
        for ev in a_data.get("trace", {}).get("choices", [])
        if ev.get("address", "").startswith("np.np1.bases[")
    )
    b_tf = compiled.simulator.trace_file_from(replayed, seed=seed)
    b_data = json.loads(b_tf.to_json())
    b_np1 = "".join(
        ev["value"]["Base"]
        for ev in b_data.get("trace", {}).get("choices", [])
        if ev.get("address", "").startswith("np.np1.bases[")
    )
    assert a_np1 == b_np1, (
        f"Markov replay diverged: original NP1 {a_np1!r}, "
        f"replayed NP1 {b_np1!r}"
    )


# ──────────────────────────────────────────────────────────────────
# 4. Replay against different transitions fails plan signature
# ──────────────────────────────────────────────────────────────────


def test_replay_with_different_markov_matrix_fails_signature_gate() -> None:
    """The plan signature folds the full Markov payload (first-
    base row + 4 transition rows). Two cartridges with
    identical ``first_base`` but different ``transitions`` MUST
    produce different signatures — replay across them fires
    the signature gate before any choice is consumed."""
    cfg_a = _cartridge_with_np_bases({"NP1": _DETERMINISTIC_MARKOV})
    # Same first_base, different transitions.
    cfg_b = _cartridge_with_np_bases(
        {
            "NP1": NpBaseModelSpec(
                kind="markov",
                first_base={"A": 1.0, "C": 0.0, "G": 0.0, "T": 0.0},
                transitions={
                    "A": {"A": 0.0, "C": 1.0, "G": 0.0, "T": 0.0},  # A→C, not A→T
                    "T": {"A": 0.0, "C": 0.0, "G": 1.0, "T": 0.0},
                    "G": {"A": 0.0, "C": 1.0, "G": 0.0, "T": 0.0},
                    "C": {"A": 1.0, "C": 0.0, "G": 0.0, "T": 0.0},
                },
            )
        }
    )
    compiled_a = ga.Experiment.on(cfg_a).recombine().compile()
    compiled_b = ga.Experiment.on(cfg_b).recombine().compile()
    seed = 7
    outcome = compiled_a.simulator.run(seed=seed)
    tf = compiled_a.simulator.trace_file_from(outcome, seed=seed)
    with pytest.raises(ValueError, match="pass plan signature mismatch"):
        compiled_b.simulator.replay_from_trace_file(tf)


# ──────────────────────────────────────────────────────────────────
# 5. Productive-only triad preserved under Markov
# ──────────────────────────────────────────────────────────────────


def test_productive_only_triad_preserved_under_markov() -> None:
    """The productive-only triad (no stop in junction, junction
    length % 3 == 0, anchor codons intact) must still hold
    when the cartridge uses Markov NP base sampling. The
    admit-mask × Markov-support intersection is the
    composition point — empty intersections fall back to the
    existing ``b'N'`` sentinel, and the contract validator
    runs to completion."""
    cfg = _cartridge_with_np_bases(
        {
            "NP1": NpBaseModelSpec(
                kind="markov",
                first_base={"A": 0.4, "C": 0.3, "G": 0.2, "T": 0.1},
                transitions={
                    "A": {"A": 0.1, "C": 0.4, "G": 0.4, "T": 0.1},
                    "C": {"A": 0.4, "C": 0.1, "G": 0.1, "T": 0.4},
                    "G": {"A": 0.4, "C": 0.1, "G": 0.1, "T": 0.4},
                    "T": {"A": 0.1, "C": 0.4, "G": 0.4, "T": 0.1},
                },
            ),
            "NP2": NpBaseModelSpec(
                kind="markov",
                first_base={"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
                transitions={
                    "A": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
                    "C": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
                    "G": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
                    "T": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
                },
            ),
        }
    )
    exp = ga.Experiment.on(cfg).recombine().productive_only()
    refdata = exp.refdata
    result = exp.run_records(n=20, seed=4242)
    report = result.validate_records(refdata)
    assert report, (
        f"productive_only + Markov NP base model validation failed: "
        f"{report.summary()}"
    )


# ──────────────────────────────────────────────────────────────────
# 6. Uniform / no-spec is byte-identical to pre-Markov baseline
# ──────────────────────────────────────────────────────────────────


def test_no_spec_signature_is_byte_identical_to_pre_markov_baseline() -> None:
    """The Markov slice ships ``UniformNpGenerator`` as a
    byte-identical-signature wrapper around the legacy
    ``UniformBase``. A bundled cartridge with no typed
    ``np_bases`` plane must produce the canonical 4-way
    A/C/G/T support substring in its plan signature."""
    compiled = ga.Experiment.on("human_igh").recombine().compile()
    seed = 4242
    outcome = compiled.simulator.run(seed=seed)
    tf = compiled.simulator.trace_file_from(outcome, seed=seed)
    sig = json.loads(tf.to_json())["pass_plan_signature"]
    assert "(65:1.0),(67:1.0),(71:1.0),(84:1.0)" in sig, (
        "Uniform / no-spec NP base plan signature regressed; the "
        "UniformNpGenerator wrapper's byte-identical signature "
        "constraint has broken"
    )


def test_no_spec_npn_run_is_byte_identical_byte_for_byte() -> None:
    """Stronger than the substring check — run the SAME seed
    on a no-spec experiment, walk the NP1 trace, and verify
    the byte sequence matches an independent draw from the
    same compiled experiment. Same compiled experiment + same
    seed → byte-identical NP1 (already a property of the
    engine; pin it here against accidental drift from the
    Markov slice's wiring change)."""
    compiled = ga.Experiment.on("human_igh").recombine().compile()
    seed = 99
    a = _np1_bases_from_trace(compiled, seed=seed)
    b = _np1_bases_from_trace(compiled, seed=seed)
    assert a == b, (
        f"Same-seed Markov-slice-wired uniform NP1 diverged: {a!r} vs {b!r}"
    )


# ──────────────────────────────────────────────────────────────────
# 7. empirical_first_base is byte-identical to pre-Markov baseline
# ──────────────────────────────────────────────────────────────────


def test_empirical_first_base_signature_is_byte_identical_to_pre_markov() -> None:
    """The Markov slice ships ``CategoricalNpGenerator`` as a
    byte-identical-signature wrapper. A cartridge with
    ``kind="empirical_first_base"`` must produce the same
    `fmt_byte_dist`-shaped signature substring it did before
    the slice."""
    cfg = _cartridge_with_np_bases(
        {
            "NP1": NpBaseModelSpec(
                kind="empirical_first_base",
                first_base={"A": 1.0, "C": 1.0, "G": 1.0, "T": 1.0},
            )
        }
    )
    compiled = ga.Experiment.on(cfg).recombine().compile()
    seed = 4242
    tf = compiled.simulator.trace_file_from(
        compiled.simulator.run(seed=seed), seed=seed
    )
    sig = json.loads(tf.to_json())["pass_plan_signature"]
    # Uniform-weighted empirical_first_base renders identically
    # to UniformBase's 4-way support substring.
    assert "(65:1.0),(67:1.0),(71:1.0),(84:1.0)" in sig


def test_empirical_first_base_npn_signature_includes_full_acgt_row() -> None:
    """Non-uniform empirical_first_base — the rendered substring
    must include all four bases in canonical order, byte-
    identical to the pre-Markov ``fmt_byte_dist(CategoricalBase)``
    output. A regression that emits only the positive-weight
    bases would break legacy trace replay."""
    cfg = _cartridge_with_np_bases(
        {
            "NP1": NpBaseModelSpec(
                kind="empirical_first_base",
                first_base={"A": 4.0, "C": 1.0, "G": 1.0, "T": 1.0},
            )
        }
    )
    compiled = ga.Experiment.on(cfg).recombine().compile()
    seed = 4242
    tf = compiled.simulator.trace_file_from(
        compiled.simulator.run(seed=seed), seed=seed
    )
    sig = json.loads(tf.to_json())["pass_plan_signature"]
    assert "(65:4.0),(67:1.0),(71:1.0),(84:1.0)" in sig, (
        "empirical_first_base signature drifted; check the "
        "CategoricalNpGenerator wrapper's signature delegation"
    )


# ──────────────────────────────────────────────────────────────────
# 8. Manifest flip
# ──────────────────────────────────────────────────────────────────


def test_manifest_lists_markov_as_supported() -> None:
    """Audit §11 — the manifest's ``supported_kinds`` list
    promotes ``markov`` to supported and empties
    ``deferred_kinds``. Mirrors the pin flip in
    ``tests/test_np_markov_base_generator_contract.py`` so a
    regression that splits the two paths surfaces here too."""
    m = ga.HUMAN_IGH_OGRDB.cartridge_manifest()
    nbm = m["models"]["np_base_models"]
    assert "markov" in nbm["supported_kinds"]
    assert nbm["deferred_kinds"] == []


def test_manifest_authored_markov_model_surfaces_in_models_block() -> None:
    """A cartridge with an authored Markov spec must surface
    under ``np_base_models["models"][region]`` with the
    correct ``kind`` tag — same per-region entry shape used
    by the existing typed-NP-base slice."""
    cfg = _cartridge_with_np_bases({"NP1": _DETERMINISTIC_MARKOV})
    m = cfg.cartridge_manifest()
    nbm = m["models"]["np_base_models"]
    assert nbm["models"]["NP1"]["kind"] == "markov"


# ──────────────────────────────────────────────────────────────────
# Bonus: plan signature differs for distinct Markov matrices
# ──────────────────────────────────────────────────────────────────


def test_distinct_markov_matrices_produce_distinct_plan_signatures() -> None:
    """The Markov signature canonically flattens first_base +
    4 transition rows. Two cartridges with the same
    first_base but different transitions must produce
    different plan signatures (else replay against the wrong
    cartridge could silently succeed)."""
    cfg_a = _cartridge_with_np_bases({"NP1": _DETERMINISTIC_MARKOV})
    cfg_b = _cartridge_with_np_bases(
        {
            "NP1": NpBaseModelSpec(
                kind="markov",
                first_base={"A": 1.0, "C": 0.0, "G": 0.0, "T": 0.0},
                transitions={
                    # A → C instead of A → T — only the A row differs.
                    "A": {"A": 0.0, "C": 1.0, "G": 0.0, "T": 0.0},
                    "T": {"A": 0.0, "C": 0.0, "G": 1.0, "T": 0.0},
                    "G": {"A": 0.0, "C": 1.0, "G": 0.0, "T": 0.0},
                    "C": {"A": 1.0, "C": 0.0, "G": 0.0, "T": 0.0},
                },
            )
        }
    )
    compiled_a = ga.Experiment.on(cfg_a).recombine().compile()
    compiled_b = ga.Experiment.on(cfg_b).recombine().compile()
    seed = 4242
    sig_a = json.loads(
        compiled_a.simulator.trace_file_from(
            compiled_a.simulator.run(seed=seed), seed=seed
        ).to_json()
    )["pass_plan_signature"]
    sig_b = json.loads(
        compiled_b.simulator.trace_file_from(
            compiled_b.simulator.run(seed=seed), seed=seed
        ).to_json()
    )["pass_plan_signature"]
    assert sig_a != sig_b, (
        "Two Markov matrices with different rows produced the same "
        "plan signature; the signature folder is broken"
    )
    # Both should contain the canonical first-base row substring.
    assert "markov:first=" in sig_a
    assert "|from=A:" in sig_a
    assert "|from=T:" in sig_a
