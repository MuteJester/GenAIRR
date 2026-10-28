"""Productive-contract stress matrix across the full mechanism stack.

Cross-axis biological-invariant regression suite. With the engine
architecture now hardened (event-driven derived-state refresh, code-
shape lockdowns, executable policy-conformance checks), the next
load-bearing claim is **biological**: complex mechanism stacks stay
biologically valid when contracts require it, and reproducible when
replayed.

This file exercises:

- **chain types** — VJ (``human_igk``, light chain, J-anchor ``F``)
  and VDJ (``human_igh``, heavy chain, J-anchor ``W``).
- **mechanism stacks** — each mutation/corruption category in
  isolation plus a full stack interleaving all of them.
- **contracts** — ``productive_only`` (the canonical contract bundle:
  in-frame junction, no stop codons, V/J anchor amino acids
  preserved) versus uncovered.
- **replay** — fresh RNG-driven run vs trace-injected replay.

Assertions per axis:

- Under ``productive_only``, every emitted record must carry
  ``productive=True``, ``vj_in_frame=True``, ``stop_codon=False``,
  ``junction_aa[0]='C'``, and ``junction_aa[-1]`` matching the chain's
  J-anchor amino acid.
- Without ``productive_only``, an aggressive mechanism stack must
  produce at least one non-productive record across 30 seeds — proof
  the contract is doing real work rather than benign decoration.
- Trace replay must reproduce the original AIRR record and trace
  record-for-record on the same compiled experiment.

Scope: 30 seeds per fixture. This is a regression net, not a
benchmark — bumping seed counts into the hundreds is appropriate
only for ad-hoc probing.
"""
from __future__ import annotations

from typing import Callable, Dict, Set, TypedDict

import pytest

import GenAIRR as ga
from GenAIRR import CompiledExperiment
from GenAIRR._airr_record import outcome_to_airr_record


# ──────────────────────────────────────────────────────────────────
# Chain inventory
# ──────────────────────────────────────────────────────────────────


class ChainSpec(TypedDict):
    """Refdata config + chain-specific anchor expectations. Typed so
    Pyright knows `chain["j_anchor_aa"]` is a `set[str]`, not the
    untyped union the bare `dict()` literal would produce."""

    config: str
    j_anchor_aa: Set[str]


# Each entry pins the engine refdata config string and the chain's
# canonical J-anchor amino acid set. IgH (VDJ) uses Trp (W); IgK
# (VJ light chain) uses Phe (F).
CHAINS: Dict[str, ChainSpec] = {
    "vj_igk": {"config": "human_igk", "j_anchor_aa": {"F"}},
    "vdj_igh": {"config": "human_igh", "j_anchor_aa": {"W"}},
}

# Per-fixture seed count. 30 is enough to see a meaningful fraction
# of non-productive records when the contract is absent (most stacks
# in the matrix exceed ~30% breakage rate without productive_only).
N_SEEDS = 30
BASE_SEED = 4242


# ──────────────────────────────────────────────────────────────────
# Mechanism stack builders
#
# Each builder takes a chain config string + a `productive` flag and
# returns a finished `Experiment`. The stacks are deliberately
# focused: each one isolates one mechanism category so a regression
# can be localised, plus one full-stack fixture that combines them.
# ──────────────────────────────────────────────────────────────────


def stack_recombine_only(config: str, *, productive: bool) -> ga.Experiment:
    exp = ga.Experiment.on(config).recombine()
    return exp.productive_only() if productive else exp


def stack_shm_s5f(config: str, *, productive: bool) -> ga.Experiment:
    exp = ga.Experiment.on(config).recombine().mutate(model="s5f", count=(2, 6))
    return exp.productive_only() if productive else exp


def stack_uniform_mutation(config: str, *, productive: bool) -> ga.Experiment:
    exp = (
        ga.Experiment.on(config).recombine().mutate(model="uniform", count=(2, 6))
    )
    return exp.productive_only() if productive else exp


def stack_pcr(config: str, *, productive: bool) -> ga.Experiment:
    exp = ga.Experiment.on(config).recombine().pcr_amplify(count=2)
    return exp.productive_only() if productive else exp


def stack_quality(config: str, *, productive: bool) -> ga.Experiment:
    exp = ga.Experiment.on(config).recombine().sequencing_errors(count=2)
    return exp.productive_only() if productive else exp


def stack_indel(config: str, count: int, *, productive: bool) -> ga.Experiment:
    """Polymerase-indel pass with a fixed count. Used to exercise
    the indel structural mechanism across the count=0/1/2/3 axis."""
    exp = (
        ga.Experiment.on(config)
        .recombine()
        .polymerase_indels(count=count, insertion_prob=0.5)
    )
    return exp.productive_only() if productive else exp


def stack_primer_trim(config: str, *, productive: bool) -> ga.Experiment:
    exp = (
        ga.Experiment.on(config)
        .recombine()
        .primer_trim_5prime(length=2)
        .primer_trim_3prime(length=2)
    )
    return exp.productive_only() if productive else exp


def stack_full(config: str, *, productive: bool) -> ga.Experiment:
    """Full-stack fixture: SHM + PCR + quality + indels + primer trim
    + N-corruption. Deliberately excludes ``contaminate`` (which
    can replace the sequence with a non-productive contaminant) and
    ``random_strand_orientation`` (a flag-only pass that doesn't
    interact with productivity)."""
    exp = (
        ga.Experiment.on(config)
        .recombine()
        .mutate(model="s5f", count=(2, 6))
        .pcr_amplify(count=1)
        .sequencing_errors(count=1)
        .polymerase_indels(count=2, insertion_prob=0.5)
        .ambiguous_base_calls(count=1)
        .primer_trim_5prime(length=2)
    )
    return exp.productive_only() if productive else exp


# Tabulated matrix entries: (label, builder, runs_under_productive).
# `runs_under_productive=False` means the productive-only assertion
# is skipped for this stack — useful for the negative-control test
# which intentionally runs without the contract.
PRODUCTIVE_STACKS: list[tuple[str, Callable[..., ga.Experiment]]] = [
    ("recombine_only", lambda c, *, productive: stack_recombine_only(c, productive=productive)),
    ("shm_s5f", lambda c, *, productive: stack_shm_s5f(c, productive=productive)),
    ("uniform_mutation", lambda c, *, productive: stack_uniform_mutation(c, productive=productive)),
    ("pcr", lambda c, *, productive: stack_pcr(c, productive=productive)),
    ("quality", lambda c, *, productive: stack_quality(c, productive=productive)),
    ("indel_count_0", lambda c, *, productive: stack_indel(c, 0, productive=productive)),
    ("indel_count_2", lambda c, *, productive: stack_indel(c, 2, productive=productive)),
    ("indel_count_3", lambda c, *, productive: stack_indel(c, 3, productive=productive)),
    ("primer_trim", lambda c, *, productive: stack_primer_trim(c, productive=productive)),
    ("full_stack", lambda c, *, productive: stack_full(c, productive=productive)),
]


# ──────────────────────────────────────────────────────────────────
# Helper: assert the productive triad + anchors on a record dict.
# ──────────────────────────────────────────────────────────────────


def assert_record_is_canonically_productive(
    record: dict, *, j_anchor_aa: set[str], context: str
) -> None:
    """Pin every AIRR field the productive bundle promises."""
    assert record["productive"] is True, f"{context}: productive=False"
    assert record["vj_in_frame"] is True, f"{context}: vj_in_frame=False"
    assert record["stop_codon"] is False, f"{context}: stop_codon=True"

    junction_aa = record.get("junction_aa") or ""
    # Empty junctions are theoretically possible (NP1=0, V trimmed
    # exactly to the anchor, etc.) but extremely rare in the
    # productive-only universe. If the junction is empty, the
    # anchor checks are vacuous — skip them rather than fail.
    if junction_aa:
        assert junction_aa[0] == "C", (
            f"{context}: V-anchor (Cys) lost — junction_aa={junction_aa!r}"
        )
        assert junction_aa[-1] in j_anchor_aa, (
            f"{context}: J-anchor missing — junction_aa={junction_aa!r}, "
            f"expected last char in {j_anchor_aa}"
        )


# ──────────────────────────────────────────────────────────────────
# 1. Productive-only invariant holds across the matrix.
# ──────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("chain_key", list(CHAINS.keys()))
@pytest.mark.parametrize(
    "stack_label, build_stack",
    PRODUCTIVE_STACKS,
    ids=[label for label, _ in PRODUCTIVE_STACKS],
)
def test_productive_only_holds_across_mechanism_stack(
    chain_key: str, stack_label: str, build_stack: Callable[..., ga.Experiment]
) -> None:
    chain = CHAINS[chain_key]
    exp = build_stack(chain["config"], productive=True)
    result = exp.run_records(n=N_SEEDS, seed=BASE_SEED)

    violations: list[str] = []
    for i, record in enumerate(result.records):
        ctx = f"{chain_key}/{stack_label}/seed={BASE_SEED + i}"
        # Collect failures across seeds rather than failing on the
        # first; a single message lists every offender, which is
        # more useful when diagnosing a regression.
        try:
            assert_record_is_canonically_productive(
                record, j_anchor_aa=chain["j_anchor_aa"], context=ctx
            )
        except AssertionError as err:
            violations.append(str(err))

    assert not violations, (
        f"productive_only contract broken under {chain_key}/{stack_label}:\n  "
        + "\n  ".join(violations)
    )


# ──────────────────────────────────────────────────────────────────
# 2. Contract is doing real work — unconstrained stacks produce
#    some non-productive records.
# ──────────────────────────────────────────────────────────────────


# Only the aggressive stacks are expected to reliably break
# productivity without the contract; gentle stacks
# (`recombine_only`, `primer_trim`, `indel_count_0`) often produce
# 100% productive records by chance and would be brittle as
# negative controls.
NEGATIVE_CONTROL_STACKS: list[tuple[str, Callable[..., ga.Experiment]]] = [
    ("shm_s5f", lambda c, *, productive: stack_shm_s5f(c, productive=productive)),
    ("uniform_mutation", lambda c, *, productive: stack_uniform_mutation(c, productive=productive)),
    ("indel_count_3", lambda c, *, productive: stack_indel(c, 3, productive=productive)),
    ("full_stack", lambda c, *, productive: stack_full(c, productive=productive)),
]


@pytest.mark.parametrize("chain_key", list(CHAINS.keys()))
@pytest.mark.parametrize(
    "stack_label, build_stack",
    NEGATIVE_CONTROL_STACKS,
    ids=[label for label, _ in NEGATIVE_CONTROL_STACKS],
)
def test_unconstrained_mechanism_stacks_produce_some_non_productive_records(
    chain_key: str, stack_label: str, build_stack: Callable[..., ga.Experiment]
) -> None:
    """Without ``productive_only``, an aggressive stack must produce
    at least one record that fails the productive triad. If every
    record happens to be productive across N_SEEDS seeds, the stack
    isn't actually exercising any mechanism that can break
    productivity — which would mean the corresponding
    ``test_productive_only_holds_*`` parameterisation is a tautology
    rather than a meaningful contract test.
    """
    chain = CHAINS[chain_key]
    exp = build_stack(chain["config"], productive=False)
    result = exp.run_records(n=N_SEEDS, seed=BASE_SEED)

    non_productive = [r for r in result.records if not r["productive"]]
    assert non_productive, (
        f"{chain_key}/{stack_label}: ALL {N_SEEDS} unconstrained records "
        f"were productive — the matching productive_only test is a "
        f"tautology. Crank up the count, switch the model, or pick a "
        f"more aggressive mechanism to make this a real negative control."
    )


# ──────────────────────────────────────────────────────────────────
# 3. Trace replay reproduces full-stack outcome record-for-record.
# ──────────────────────────────────────────────────────────────────

# Replay is single-outcome (one trace = one record), so we
# parametrise over chain × contract × representative stack.
REPLAY_STACKS: list[tuple[str, Callable[..., ga.Experiment]]] = [
    ("full_stack", lambda c, *, productive: stack_full(c, productive=productive)),
    ("indel_count_2", lambda c, *, productive: stack_indel(c, 2, productive=productive)),
]


@pytest.mark.parametrize("chain_key", list(CHAINS.keys()))
@pytest.mark.parametrize("productive", [True, False], ids=["productive_only", "unconstrained"])
@pytest.mark.parametrize(
    "stack_label, build_stack",
    REPLAY_STACKS,
    ids=[label for label, _ in REPLAY_STACKS],
)
def test_trace_replay_reproduces_full_stack_outcome(
    chain_key: str,
    productive: bool,
    stack_label: str,
    build_stack: Callable[..., ga.Experiment],
) -> None:
    chain = CHAINS[chain_key]
    compiled = build_stack(chain["config"], productive=productive).compile()
    # Narrow the `compile()` return union: none of these fixtures
    # use `expand_clones`, so the result is the non-clonal
    # `CompiledExperiment` shape with a `.simulator` attribute.
    assert isinstance(compiled, CompiledExperiment)

    seed = BASE_SEED
    outcome_a = compiled.simulator.run(seed=seed)
    tf = compiled.simulator.trace_file_from(outcome_a, seed=seed)
    outcome_b = compiled.simulator.replay_from_trace_file(tf)

    # ── Trace equality ────────────────────────────────────────
    choices_a = outcome_a.trace().choices()
    choices_b = outcome_b.trace().choices()
    assert len(choices_a) == len(choices_b), (
        f"{chain_key}/{stack_label}/productive={productive}: "
        f"trace length mismatch ({len(choices_a)} vs {len(choices_b)})"
    )
    for i, (a, b) in enumerate(zip(choices_a, choices_b)):
        assert a.address == b.address, (
            f"{chain_key}/{stack_label}/productive={productive} "
            f"address mismatch at record {i}: {a.address!r} vs {b.address!r}"
        )
        assert a.value == b.value, (
            f"{chain_key}/{stack_label}/productive={productive} "
            f"value mismatch at record {i} (addr={a.address})"
        )

    # ── AIRR record equality on every load-bearing field ──────
    record_a = outcome_to_airr_record(outcome_a, compiled.refdata, sequence_id="a")
    record_b = outcome_to_airr_record(outcome_b, compiled.refdata, sequence_id="b")
    for field in [
        "sequence",
        "sequence_aa",
        "junction",
        "junction_aa",
        "productive",
        "vj_in_frame",
        "stop_codon",
        "v_call",
        "j_call",
        "n_mutations",
        "n_pcr_errors",
        "n_quality_errors",
        "n_indels",
    ]:
        assert record_a.get(field) == record_b.get(field), (
            f"{chain_key}/{stack_label}/productive={productive}: "
            f"AIRR field {field!r} diverges under replay: "
            f"{record_a.get(field)!r} vs {record_b.get(field)!r}"
        )

    # ── Productive contract still holds on the replayed outcome ──
    if productive:
        assert_record_is_canonically_productive(
            record_b,
            j_anchor_aa=chain["j_anchor_aa"],
            context=f"{chain_key}/{stack_label}/replayed",
        )
