"""Productive-contract failure-mode golden tests.

Companion to [docs/productive_failure_mode_audit.md](../docs/productive_failure_mode_audit.md).
Pins the failure matrix from audit §4 — for each detectable
failure mode of the productive bundle (ProductiveJunctionFrame,
NoStopCodonInJunction, AnchorPreserved.V, AnchorPreserved.J), the
test asserts the strict-mode error, the permissive-mode sentinel,
and the diagnostic surface.

Counterpart to ``test_productive_stress_matrix.py`` which proves
"valid stacks stay productive." This file proves the dual:
"invalid choices fail predictably."
"""
from __future__ import annotations

import pytest

import GenAIRR as ga
from GenAIRR import _engine as ge


# ──────────────────────────────────────────────────────────────────
# Refdata factories
# ──────────────────────────────────────────────────────────────────


def _vj_basic() -> "ge.RefDataConfig":
    """V[0:9] anchor=6; J[0:6] anchor=0. V_anchor_to_end=3,
    J_anchor_offset=0 → junction = V[6:9] + NP1 + J[0:3], length
    = NP1+6. In-frame requires NP1 ≡ 0 mod 3."""
    cfg = ge.RefDataConfig.vj()
    cfg.add_v_allele("v1*01", "v1", b"AAACCCTGT", anchor=6)
    cfg.add_j_allele("j1*01", "j1", b"TTTAAA", anchor=0)
    return cfg


def _vj_short_j() -> "ge.RefDataConfig":
    """V[0:9] anchor=6; J[0:3] anchor=0. J is its own anchor codon
    with no post-anchor zone — used to construct the indel
    no-admissible-tuple fixture."""
    cfg = ge.RefDataConfig.vj()
    cfg.add_v_allele("v1*01", "v1", b"AAACCCTGT", anchor=6)
    cfg.add_j_allele("j1*01", "j1", b"TTT", anchor=0)
    return cfg


def _vj_stop_completing_v() -> "ge.RefDataConfig":
    """Canonical Cys is followed by an incomplete TA codon at the V end.
    One NP base completes that codon; A/G would create a stop."""
    cfg = ge.RefDataConfig.vj()
    cfg.add_v_allele("v1*01", "v1", b"AAACCCTGTTA", anchor=6)
    cfg.add_j_allele("j1*01", "j1", b"TTTAAA", anchor=0)
    return cfg


# ──────────────────────────────────────────────────────────────────
# 1. Compile-time precondition failures (F1, F2, F3)
# ──────────────────────────────────────────────────────────────────


def test_all_invalid_np_lengths_fails_at_compile_time() -> None:
    """Audit F1: when every NP1 length in the distribution
    violates the productive frame constraint (junction not
    divisible by 3), `Experiment.compile(allow_curatable_refdata=True)` raises `ValueError`
    with a precondition-failure message. This is fail-fast — the
    error surfaces before any record runs, and is independent of
    `strict=True/False`."""
    # NP1 lengths 1, 2, 4, 5 all produce junction lengths
    # 7, 8, 10, 11 — none divisible by 3.
    exp = (
        ga.Experiment.on(_vj_basic()).allow_curatable_refdata()
        .recombine(np1_lengths=[(1, 1.0), (2, 1.0), (4, 1.0), (5, 1.0)])
        .trim(enabled=False)
        .productive_only()
    )
    with pytest.raises(ValueError) as excinfo:
        exp.run_records(n=1, seed=0)
    msg = str(excinfo.value)
    assert "productive_junction_frame" in msg
    assert "precondition failed" in msg
    assert "no productive completion" in msg


def test_all_invalid_v_trim_3_fails_at_compile_time() -> None:
    """Audit F2: the terminal Cys occupies the last three V bases.
    Any positive V3 trim removes part of that codon. A
    distribution of {4, 5} eliminates every valid anchor."""
    exp = (
        ga.Experiment.on(_vj_basic()).allow_curatable_refdata()
        .recombine(np1_lengths=[(3, 1.0)])
        .trim(v_3=[(4, 1.0), (5, 1.0)])
        .productive_only()
    )
    with pytest.raises(ValueError) as excinfo:
        exp.run_records(n=1, seed=0)
    msg = str(excinfo.value)
    assert "productive_junction_frame" in msg
    assert "no productive completion" in msg


def test_all_invalid_j_trim_5_fails_at_compile_time() -> None:
    """Audit F3: J_anchor=0 means the J anchor codon starts at
    J[0]. Any j_trim_5 > 0 removes part of the anchor codon. A
    distribution of {1, 2, 3} eliminates every valid anchor."""
    exp = (
        ga.Experiment.on(_vj_basic()).allow_curatable_refdata()
        .recombine(np1_lengths=[(3, 1.0)])
        .trim(j_5=[(1, 1.0), (2, 1.0), (3, 1.0)])
        .productive_only()
    )
    with pytest.raises(ValueError) as excinfo:
        exp.run_records(n=1, seed=0)
    msg = str(excinfo.value)
    assert "productive_junction_frame" in msg
    assert "no productive completion" in msg


def test_compile_time_failure_raises_same_value_error_in_strict_mode() -> None:
    """Compile-time errors fire from `compile()`, before the
    runtime `strict` flag is consulted. The flag doesn't change
    the exception type for these failures."""
    exp = (
        ga.Experiment.on(_vj_basic()).allow_curatable_refdata()
        .recombine(np1_lengths=[(1, 1.0)])  # always out-of-frame
        .trim(enabled=False)
        .productive_only()
    )
    with pytest.raises(ValueError):
        exp.run_records(n=1, seed=0, strict=False)
    with pytest.raises(ValueError):
        exp.run_records(n=1, seed=0, strict=True)


# ──────────────────────────────────────────────────────────────────
# 2. Runtime narrowing on mixed supports (F4, F5)
# ──────────────────────────────────────────────────────────────────


def test_mixed_np_length_support_narrows_at_sample_time() -> None:
    """Audit F4: when the NP length distribution mixes valid and
    invalid candidates, the static precondition passes (at least
    one in-contract candidate exists), and the runtime sampler
    narrows per-call. Across many seeds, only in-contract lengths
    appear in the trace."""
    exp = (
        ga.Experiment.on(_vj_basic()).allow_curatable_refdata()
        .recombine(np1_lengths=[(1, 1.0), (3, 1.0), (6, 1.0)])  # 1 invalid, 3 & 6 valid
        .trim(enabled=False)
        .productive_only()
    )
    observed = set()
    for seed in range(50):
        rec = exp.run_records(n=1, seed=seed).records[0]
        observed.add(rec["np1_length"])
    assert 1 not in observed, (
        f"out-of-frame length 1 leaked through; observed={observed}"
    )
    assert observed.issubset({3, 6})


def test_mixed_v_trim_support_narrows_at_sample_time() -> None:
    """Audit F5: the complete terminal Cys codon must survive.
    Only zero V3 trimming retains all three landmark bases."""
    exp = (
        ga.Experiment.on(_vj_basic()).allow_curatable_refdata()
        .recombine(np1_lengths=[(3, 1.0)])
        .trim(v_3=[(0, 1.0), (1, 1.0), (2, 1.0), (3, 1.0), (5, 1.0)])
        .productive_only()
    )
    observed = set()
    for seed in range(30):
        rec = exp.run_records(n=1, seed=seed).records[0]
        observed.add(rec["v_trim_3"])
    assert observed == {0}, (
        f"anchor-violating trim values leaked through; observed={observed}"
    )


# ──────────────────────────────────────────────────────────────────
# 3. NP base mask vs stop codons (F6)
# ──────────────────────────────────────────────────────────────────


def test_np_base_mask_excludes_stop_completing_bases() -> None:
    """A partial TA codon after canonical Cys must complete as TAC/TAT.

    This exercises a real codon spanning the V/NP boundary. The conserved
    Cys remains intact and the one-base NP makes the joint in-frame.
    """
    exp = (
        ga.Experiment.on(_vj_stop_completing_v())
        .recombine(np1_lengths=[(1, 1.0)])
        .trim(enabled=False)
        .productive_only()
    )
    observed = set()
    for seed in range(50):
        rec = exp.run_records(n=1, seed=seed).records[0]
        np1 = rec["np1"]
        observed.add(np1)
        assert np1 in {"C", "T"}, f"seed {seed}: stop-completing NP base {np1!r}"
        assert rec["junction"].startswith("TGTTA" + np1)
        assert rec["productive"] is True
    assert observed == {"C", "T"}


# ──────────────────────────────────────────────────────────────────
# 4. Indel no-admissible-tuple (F7)
# ──────────────────────────────────────────────────────────────────


def test_indel_no_admissible_tuple_raises_strict_emits_no_op_permissive() -> None:
    """Audit F7: short-J fixture (J = anchor codon only) +
    delete-only count=1 + productive_only → no admissible
    deletion site exists. In permissive mode the slot becomes a
    NoOp (trace site = -1, no event fires). In strict mode a
    `StrictSamplingError` is raised at the indel pass."""
    exp = (
        ga.Experiment.on(_vj_short_j()).allow_curatable_refdata()
        .recombine(np1_lengths=[(3, 1.0)])
        .trim(enabled=False)
        .polymerase_indels(count=1, insertion_prob=0.0)
        .productive_only()
    )

    # Permissive: NoOp sentinel, no event.
    compiled = exp.compile(allow_curatable_refdata=True)
    outcome = compiled.simulator.run(seed=0)
    site_rec = next(
        r for r in outcome.trace().choices()
        if r.address == "corrupt.indel.site[0]"
    )
    assert int(str(site_rec.value).strip("Int()")) == -1

    indel_pass = next(
        e for e in outcome.events() if e.pass_name == "corrupt.indel"
    )
    assert indel_pass.simulation_event_count == 0

    # Strict: raises with structured diagnosis.
    with pytest.raises(ge.StrictSamplingError) as excinfo:
        exp.run_records(n=1, seed=0, strict=True)
    args = excinfo.value.args
    assert args[0] == "corrupt.indel"
    # The finite continuation solver detects empty support for the entire
    # requested count before sampling a concrete edit site.
    assert args[1] == "corrupt.indel.count"
    assert args[2] == "empty_admissible_support"


# ──────────────────────────────────────────────────────────────────
# 5. SHM under productive — high rate doesn't fail strict (F8)
# ──────────────────────────────────────────────────────────────────


def test_shm_high_rate_under_productive_never_fails_strict() -> None:
    """High-rate SHM preserves productivity across synonymous constraints.

    Positions with no valid substitution may remain unchanged; this does not
    claim that every requested mutation can alter an invariant landmark base.
    """
    exp = (
        ga.Experiment.on(_vj_basic()).allow_curatable_refdata()
        .recombine(np1_lengths=[(3, 1.0)])
        .trim(enabled=False)
        .mutate(rate=0.99)
        .productive_only()
    )
    # Should not raise across seeds. Verify records are productive.
    for seed in range(10):
        rec = exp.run_records(n=1, seed=seed, strict=True).records[0]
        assert rec["productive"] is True


# ──────────────────────────────────────────────────────────────────
# 6. Diagnostic quality
# ──────────────────────────────────────────────────────────────────


def test_strict_error_carries_pass_address_reason() -> None:
    """`StrictSamplingError` exposes a 3-tuple of
    `(pass_name, choice_address, reason_string)` via `args`. The
    reason string is one of the canonical
    `FilteredSampleError` variants: `support_unavailable`,
    `empty_admissible_support`, or `invalid_filtered_support`."""
    exp = (
        ga.Experiment.on(_vj_short_j()).allow_curatable_refdata()
        .recombine(np1_lengths=[(3, 1.0)])
        .trim(enabled=False)
        .polymerase_indels(count=1, insertion_prob=0.0)
        .productive_only()
    )
    with pytest.raises(ge.StrictSamplingError) as excinfo:
        exp.run_records(n=1, seed=0, strict=True)
    err = excinfo.value
    assert isinstance(err.args, tuple)
    assert len(err.args) == 3
    pass_name, address, reason = err.args
    assert isinstance(pass_name, str) and pass_name
    assert isinstance(address, str) and address
    assert reason in (
        "support_unavailable",
        "empty_admissible_support",
        "invalid_filtered_support",
    )


# ──────────────────────────────────────────────────────────────────
# 7. Replay drift — sentinel trace under strict mode (§6.2)
# ──────────────────────────────────────────────────────────────────


def test_replay_of_permissive_sentinel_trace_passes_under_strict() -> None:
    """Audit §6.2: a trace recorded in permissive mode with a
    NoOp sentinel value (site = -1 on the indel slot) replays
    under `strict=True` without raising. Replay consumes the
    recorded value verbatim and reproduces the original outcome.

    This is the documented behaviour (see docstrings on
    ``CompiledExperiment.run`` and
    ``CompiledExperiment.replay_from_trace_file``): ``strict=True``
    means "fail when a *fresh* sample would have empty admissible
    support," not "fail when a recorded sample would have failed
    under fresh strict sampling." Users wanting the latter should
    call ``simulator.run(seed=<trace.seed>, strict=True)`` instead
    of ``replay_from_trace_file(tf, strict=True)``.

    This test pins the documented contract end-to-end."""
    exp = (
        ga.Experiment.on(_vj_short_j()).allow_curatable_refdata()
        .recombine(np1_lengths=[(3, 1.0)])
        .trim(enabled=False)
        .polymerase_indels(count=1, insertion_prob=0.0)
        .productive_only()
    )
    compiled = exp.compile(allow_curatable_refdata=True)

    # 1. Permissive run produces a sentinel trace.
    seed = 0
    out_perm = compiled.simulator.run(seed=seed)
    tf = compiled.simulator.trace_file_from(out_perm, seed=seed)
    # Confirm the recorded value IS the sentinel.
    site_rec = next(
        r for r in out_perm.trace().choices()
        if r.address == "corrupt.indel.site[0]"
    )
    assert int(str(site_rec.value).strip("Int()")) == -1

    # 2. Strict replay of that trace doesn't raise — the documented
    #    "replay consumes sentinels verbatim" guarantee.
    out_replay = compiled.simulator.replay_from_trace_file(tf, strict=True)
    indel_pass_replay = next(
        e for e in out_replay.events() if e.pass_name == "corrupt.indel"
    )
    assert indel_pass_replay.simulation_event_count == 0  # same NoOp outcome

    # 3. The trace's recorded sentinel survives unchanged under
    #    strict replay — replay didn't reinterpret it.
    site_rec_replay = next(
        r for r in out_replay.trace().choices()
        if r.address == "corrupt.indel.site[0]"
    )
    assert int(str(site_rec_replay.value).strip("Int()")) == -1

    # 4. The documented escape hatch: a *fresh* strict run on the
    #    same seed raises StrictSamplingError. This is the call
    #    a user should make if they want strict-fresh semantics
    #    on a recorded seed.
    with pytest.raises(ge.StrictSamplingError) as exc_fresh:
        compiled.simulator.run(seed=seed, strict=True)
    assert exc_fresh.value.args[0] == "corrupt.indel"

    # 5. Confirm the divergence is fresh-vs-replay, not a strict-mode-
    #    wide relaxation: rerun_from_trace_file (which re-runs the
    #    sampler against the seed rather than consuming the trace
    #    verbatim) raises under strict, like a fresh run.
    with pytest.raises(ge.StrictSamplingError) as exc_rerun:
        compiled.simulator.rerun_from_trace_file(tf, strict=True)
    assert exc_rerun.value.args[0] == "corrupt.indel"


def test_strict_docstring_documents_fresh_vs_replay_distinction() -> None:
    """Audit §6.2 follow-through: the user-facing docstrings for
    ``CompiledExperiment.run`` and
    ``CompiledExperiment.replay_from_trace_file`` must explicitly
    mention the fresh-vs-replay distinction so a user reading
    ``help(simulator.run)`` learns of it without consulting the
    audit doc.

    If a future refactor strips the documentation, this test
    flags the regression — the documented behaviour is a
    contract, and the docstring is the API surface for it."""
    # Compile a tiny experiment so we can introspect the bound methods.
    exp = (
        ga.Experiment.on(_vj_basic()).allow_curatable_refdata()
        .recombine(np1_lengths=[(3, 1.0)])
        .trim(enabled=False)
    )
    compiled = exp.compile(allow_curatable_refdata=True)
    sim = compiled.simulator

    run_doc = sim.run.__doc__ or ""
    replay_doc = sim.replay_from_trace_file.__doc__ or ""

    # The high-level user-visible distinction: "fresh" appears in run,
    # and the cross-link to replay appears too.
    assert "fresh" in run_doc.lower(), (
        "simulator.run docstring should mention 'fresh' sampling"
    )
    assert "replay" in run_doc.lower(), (
        "simulator.run docstring should cross-reference replay semantics"
    )
    # Replay validates proposals without replacing them. P-aware productive
    # replay also checks natural support/completion under either policy.
    assert "strict" in replay_doc.lower()
    assert "recorded choices are proposals" in replay_doc.lower()
    assert "positive natural support" in replay_doc.lower()
    assert "regardless of" in replay_doc.lower()
    assert "without resampling" in replay_doc.lower()


# ──────────────────────────────────────────────────────────────────
# 8. Pin-the-drift (audit §6)
# ──────────────────────────────────────────────────────────────────


def test_drift_pinned_compile_and_runtime_failures_use_different_exception_types() -> None:
    """Audit §6.1: compile-time empty support raises Python
    `ValueError`; runtime empty support raises
    `StrictSamplingError` (which extends `Exception` directly,
    NOT `ValueError`). There's no shared base class — a user
    catching one won't catch the other. Pin both the class
    hierarchy and the divergence."""
    # Class hierarchy: StrictSamplingError is NOT a ValueError
    # subclass — `except ValueError:` won't catch a runtime
    # contract failure.
    assert not issubclass(ge.StrictSamplingError, ValueError)

    # Compile-time path: NP length all-invalid → ValueError that
    # is NOT a StrictSamplingError.
    exp_compile = (
        ga.Experiment.on(_vj_basic()).allow_curatable_refdata()
        .recombine(np1_lengths=[(1, 1.0)])
        .trim(enabled=False)
        .productive_only()
    )
    with pytest.raises(ValueError) as exc_compile:
        exp_compile.run_records(n=1, seed=0, strict=True)
    assert not isinstance(exc_compile.value, ge.StrictSamplingError)

    # Runtime path: indel no-admissible-tuple → StrictSamplingError
    # that is NOT a ValueError.
    exp_runtime = (
        ga.Experiment.on(_vj_short_j()).allow_curatable_refdata()
        .recombine(np1_lengths=[(3, 1.0)])
        .trim(enabled=False)
        .polymerase_indels(count=1, insertion_prob=0.0)
        .productive_only()
    )
    with pytest.raises(ge.StrictSamplingError) as exc_runtime:
        exp_runtime.run_records(n=1, seed=0, strict=True)
    assert not isinstance(exc_runtime.value, ValueError)


def test_drift_pinned_np_base_sentinel_path_is_structurally_unreachable() -> None:
    """Audit §6.3: finite productive completion prevents an NP dead end.
    A candidate base must retain a complete productive continuation;
    the sampler cannot repair an impossible prefix with a sentinel N.

    Pin this by asserting that across all our failure-mode
    fixtures, no `np1.base[i]` slot is ever sentinel-`N` in the
    trace."""
    exp = (
        ga.Experiment.on(_vj_stop_completing_v()).allow_curatable_refdata()
        .recombine(np1_lengths=[(1, 1.0)])
        .trim(enabled=False)
        .productive_only()
    )
    compiled = exp.compile(allow_curatable_refdata=True)
    for seed in range(30):
        out = compiled.simulator.run(seed=seed)
        for rec in out.trace().choices():
            if rec.address.startswith("np1.base["):
                # The NP-base trace value is a `ChoiceValue::Base(u8)`.
                # If the sentinel ever fires we'd see `N` here.
                val = str(rec.value)
                assert "N" not in val.upper().split("'"), (
                    f"seed {seed}: NP base trace {rec.address}={val} — "
                    f"sentinel N path is no longer unreachable"
                )
