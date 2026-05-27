"""Live-call cache equivalence tests.

Every simulation's `sim.segment_calls` cache must equal a from-scratch
recomputation of `assembled_segment_live_call` against the final IR.
This is a stronger invariant than what the AIRR-record validator checks:
the validator compares the *projected* AIRR record to the engine's
truth oracle (catching cache divergence indirectly, after projection);
this harness compares the cache directly to a fresh recompute, so a
stale-cache bug fails closer to its source.

The two bugs the AIRR validator already surfaced and we fixed
(`segment_region_overlaps_dirty` strict-`<` + NP-claim over-extension)
would each fail one or more of these tests; if those regress, this
file is the first to scream.

Each test sweeps a mechanism stack across many seeds, runs the
parity harness, and asserts `tie_set_matches` and (when both sides
have hypotheses) `hypothesis_bounds_match` for every V/D/J segment
of every record.
"""
from __future__ import annotations

import GenAIRR as ga


def _assert_parity(parity_list, *, label: str, seed: int):
    """Assert tie-set + hypothesis-bound equality across every segment."""
    for p in parity_list:
        assert p["tie_set_matches"], (
            f"[{label}] seed={seed} segment={p['segment']}: "
            f"cached tie {p['cached_tie_set']} != fresh tie {p['fresh_tie_set']}; "
            f"cached_hyp={p['cached_hypothesis']} fresh_hyp={p['fresh_hypothesis']}"
        )
        if p["hypothesis_bounds_match"] is not None:
            assert p["hypothesis_bounds_match"], (
                f"[{label}] seed={seed} segment={p['segment']}: "
                f"hypothesis bounds diverge "
                f"cached={p['cached_hypothesis']} fresh={p['fresh_hypothesis']}"
            )


def _sweep(label: str, builder, *, n_seeds: int = 50, configs=None):
    if configs is None:
        configs = ["human_igh", "human_igk", "human_igl"]
    for config in configs:
        exp = builder(config)
        refdata = exp.refdata
        compiled = exp.compile()
        for seed in range(n_seeds):
            outcome = compiled.simulator.run(seed=seed)
            parity = outcome.check_live_call_cache_parity(refdata)
            _assert_parity(parity, label=f"{label}[{config}]", seed=seed)


# ──────────────────────────────────────────────────────────────────
# Mechanism stacks
# ──────────────────────────────────────────────────────────────────


def test_cache_parity_holds_for_recombine_only():
    """Bare recombine — no SHM, no PCR, no indels, no end-loss.
    Cache should equal from-scratch trivially since no post-assembly
    events fire."""
    _sweep("recombine", lambda cfg: ga.Experiment.on(cfg).recombine())


def test_cache_parity_holds_under_productive_constraint():
    """Productive-only narrows trim and NP samples but otherwise
    leaves the post-assembly state untouched."""
    _sweep(
        "recombine + productive",
        lambda cfg: ga.Experiment.on(cfg).recombine().productive_only(),
    )


def test_cache_parity_holds_under_shm():
    """SHM is the dominant `BaseChanged` source. Each mutation
    triggers an `EditedSegments` refresh; cache must stay in sync
    with the fresh recompute. BCR-only — TCRs reject mutate()."""
    _sweep(
        "+ SHM",
        lambda cfg: ga.Experiment.on(cfg).recombine().mutate(rate=0.05),
        configs=["human_igh", "human_igk", "human_igl"],
    )


def test_cache_parity_holds_under_indels():
    """Indels emit `IndelInserted` / `IndelDeleted` events that
    trigger `AllStructural` refresh. The walker observer marks
    `needs_rebuild` for interior deletions and rebuilds at seal
    time. Both paths must produce the same tie set as a fresh
    from-scratch walk."""
    _sweep(
        "+ indels",
        lambda cfg: ga.Experiment.on(cfg)
        .recombine()
        .polymerase_indels(count=2, insertion_prob=0.5),
    )


def test_cache_parity_holds_under_end_loss():
    """End-loss removes pool bytes at the 5'/3' boundaries. The
    3'-boundary deletion was where the `segment_region_overlaps_dirty`
    strict-`<` bug lived — if it regresses, this test fails first.

    `productive_only` is included specifically because the bug only
    surfaced under productive contracts (when trim/end-loss
    distributions intersected the problematic case)."""
    _sweep(
        "+ end-loss + productive",
        lambda cfg: ga.Experiment.on(cfg)
        .recombine()
        .productive_only()
        .primer_trim_5prime(length=(0, 3))
        .primer_trim_3prime(length=(0, 3)),
    )


def test_cache_parity_holds_under_pcr_and_sequencing_errors():
    """PCR and sequencing errors are also `BaseChanged` sources but
    with different distributions (uppercase A/C/G/T for PCR,
    lowercase for sequencing). Both refresh via `EditedSegments`."""
    _sweep(
        "+ PCR + seq errors",
        lambda cfg: ga.Experiment.on(cfg)
        .recombine()
        .pcr_amplify(rate=1e-4)
        .sequencing_errors(rate=1e-4),
    )


def test_cache_parity_holds_under_full_stack():
    """All mechanisms combined — the workload the release-tier
    validator runs on. If any mechanism's refresh path drifts
    against a from-scratch walk, this catches it."""

    def builder(cfg):
        exp = ga.Experiment.on(cfg).recombine().productive_only()
        if not cfg.startswith("human_tcr"):
            exp = exp.mutate(rate=0.03)
        return (
            exp.pcr_amplify(rate=1e-4)
            .polymerase_indels(count=2, insertion_prob=0.5)
            .primer_trim_5prime(length=(0, 3))
            .primer_trim_3prime(length=(0, 3))
        )

    _sweep("full stack", builder, n_seeds=50)


# ──────────────────────────────────────────────────────────────────
# Direct regression pin for the 3' end-loss boundary bug
# ──────────────────────────────────────────────────────────────────


def test_cache_parity_pins_three_prime_end_loss_boundary_case():
    """Direct regression: the case that exposed the
    `segment_region_overlaps_dirty` strict-`<` boundary bug.
    `primer_trim_3prime(length=(1, 3))` guarantees at least one
    3' byte is deleted per record, which hits the boundary case
    on the J region's tail. Cache parity must hold."""
    exp = (
        ga.Experiment.on("human_igk")
        .recombine()
        .productive_only()
        .primer_trim_3prime(length=(1, 3))
    )
    refdata = exp.refdata
    compiled = exp.compile()
    for seed in range(200):
        outcome = compiled.simulator.run(seed=seed)
        parity = outcome.check_live_call_cache_parity(refdata)
        _assert_parity(parity, label="3'-end-loss-boundary", seed=seed)


# ──────────────────────────────────────────────────────────────────
# Harness shape contract
# ──────────────────────────────────────────────────────────────────


def test_parity_harness_returns_expected_dict_shape():
    """Pin the parity-dict shape so downstream tooling has a
    stable contract."""
    exp = ga.Experiment.on("human_igh").recombine()
    compiled = exp.compile()
    outcome = compiled.simulator.run(seed=0)
    parity = outcome.check_live_call_cache_parity(exp.refdata)
    assert parity, "expected at least one segment with a parity result"
    for p in parity:
        assert set(p.keys()) >= {
            "segment",
            "tie_set_matches",
            "cached_present",
            "fresh_present",
            "cached_tie_set",
            "fresh_tie_set",
            "hypothesis_bounds_match",
            "cached_hypothesis",
            "fresh_hypothesis",
        }
        assert p["segment"] in {"V", "D", "J"}
        assert isinstance(p["tie_set_matches"], bool)
        assert isinstance(p["cached_tie_set"], list)
        assert isinstance(p["fresh_tie_set"], list)
        if p["hypothesis_bounds_match"] is not None:
            assert isinstance(p["hypothesis_bounds_match"], bool)
        if p["cached_hypothesis"] is not None:
            assert {"seq_start", "seq_end", "ref_start", "ref_end"} <= set(
                p["cached_hypothesis"].keys()
            )
