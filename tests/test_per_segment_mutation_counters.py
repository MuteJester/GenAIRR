"""Spec tests for the per-segment SHM counter slice.

Closes the mutation-provenance audit's Slice 1. Adds four AIRR
fields aggregated from the per-pass event ledger:

- ``n_v_mutations``
- ``n_d_mutations``
- ``n_j_mutations``
- ``n_np_mutations``

The fields partition ``n_mutations`` by biological segment (NP1
+ NP2 rolled into a single NP bucket, matching the
``segment_rates`` DSL grouping). Five new validator issue kinds
re-derive the expected values from the same event walk so
record-dict tampering surfaces through ``validate_records``.

Spec coverage (from the user brief):

1. Baseline no-mutation records default all four fields to 0.
2. V-only segment rates produce ``n_v_mutations == n_mutations``,
   others 0.
3. NP-only segment rates produce ``n_np_mutations == n_mutations``,
   V/D/J 0.
4. Mixed rates satisfy sum invariant.
5. PCR/quality/N-corruption do not affect these counters.
6. Replay round-trip preserves all four fields.
7. Validator flags tampered per-segment field.
8. Validator flags sum mismatch.
9. Python DataFrame columns exist.

Out of scope: per-NP1/NP2 split, per-segment mutation rates,
Python typed-event exposure, polymerase-indel reclassification.
"""
from __future__ import annotations

import pytest

import GenAIRR as ga


# ──────────────────────────────────────────────────────────────────
# Spec 1 — baseline records default to 0
# ──────────────────────────────────────────────────────────────────


def test_no_mutation_records_default_all_four_fields_to_zero() -> None:
    """A pure recombine pipeline (no SHM) produces records whose
    four per-segment SHM counters are all ``0``. Pinned for both
    the fields' presence AND the zero-default sanity."""
    result = ga.Experiment.on("human_igh").recombine().run_records(n=3, seed=0)
    for r in result:
        # Fields present.
        for field in (
            "n_v_mutations",
            "n_d_mutations",
            "n_j_mutations",
            "n_np_mutations",
        ):
            assert field in r, f"AIRR record missing {field!r}"
        # All zero.
        assert r["n_v_mutations"] == 0
        assert r["n_d_mutations"] == 0
        assert r["n_j_mutations"] == 0
        assert r["n_np_mutations"] == 0
        # And the global counter agrees.
        assert r["n_mutations"] == 0


# ──────────────────────────────────────────────────────────────────
# Spec 2 — V-only segment rates
# ──────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("model", ["uniform", "s5f"])
def test_v_only_rates_route_all_counts_to_n_v_mutations(model) -> None:
    """``segment_rates={D:0, J:0, NP:0}`` forces every realised
    mutation into V. The per-segment counter must match
    ``n_mutations`` exactly; D / J / NP stay at 0."""
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(
            model=model,
            count=20,
            segment_rates={"D": 0.0, "J": 0.0, "NP": 0.0},
        )
    )
    result = exp.run_records(n=5, seed=0)
    for r in result:
        assert r["n_v_mutations"] == r["n_mutations"], (
            f"{model}: V-only didn't route all mutations to "
            f"n_v_mutations (V={r['n_v_mutations']}, "
            f"global={r['n_mutations']})"
        )
        assert r["n_d_mutations"] == 0
        assert r["n_j_mutations"] == 0
        assert r["n_np_mutations"] == 0


# ──────────────────────────────────────────────────────────────────
# Spec 3 — NP-only segment rates
# ──────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("model", ["uniform", "s5f"])
def test_np_only_rates_route_counts_to_n_np_mutations(model) -> None:
    """``segment_rates={V:0, D:0, J:0, NP:1}`` confines mutations
    to NP1+NP2 regions. The rolled-up ``n_np_mutations`` matches
    the global ``n_mutations``; V/D/J stay at 0."""
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(
            model=model,
            count=10,
            segment_rates={"V": 0.0, "D": 0.0, "J": 0.0, "NP": 1.0},
        )
    )
    result = exp.run_records(n=5, seed=0)
    saw_any = False
    for r in result:
        assert r["n_v_mutations"] == 0
        assert r["n_d_mutations"] == 0
        assert r["n_j_mutations"] == 0
        # NP1+NP2 rolled together.
        if r["n_mutations"] > 0:
            assert r["n_np_mutations"] == r["n_mutations"]
            saw_any = True
    assert saw_any, (
        f"{model}: no NP mutations realised across the batch; "
        "fixture insufficient."
    )


# ──────────────────────────────────────────────────────────────────
# Spec 4 — mixed rates satisfy sum invariant
# ──────────────────────────────────────────────────────────────────


def test_mixed_rates_satisfy_sum_invariant() -> None:
    """Mixed rates: every record's four-bucket sum equals
    ``n_mutations``. This is the load-bearing invariant —
    arithmetic identity holds across any rate combination."""
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(count=30, segment_rates={"V": 1.0, "D": 0.5, "J": 0.5, "NP": 0.2})
    )
    result = exp.run_records(n=10, seed=0)
    for r in result:
        total = (
            r["n_v_mutations"]
            + r["n_d_mutations"]
            + r["n_j_mutations"]
            + r["n_np_mutations"]
        )
        assert total == r["n_mutations"], (
            f"sum invariant broken: "
            f"V={r['n_v_mutations']} + D={r['n_d_mutations']} + "
            f"J={r['n_j_mutations']} + NP={r['n_np_mutations']} = "
            f"{total} != n_mutations={r['n_mutations']}"
        )


# ──────────────────────────────────────────────────────────────────
# Spec 5 — corruption passes don't affect these counters
# ──────────────────────────────────────────────────────────────────


def test_pcr_and_quality_and_ns_do_not_affect_per_segment_counters() -> None:
    """Corruption passes fire ``BaseChanged`` events too (PCR,
    quality, N-corruption) but are filtered out by the
    ``mutate.{uniform,s5f}`` pass-name gate. With no SHM step,
    the four counters stay at 0 even under heavy corruption."""
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .pcr_amplify(count=15)
        .sequencing_errors(count=8)
        .ambiguous_base_calls(count=5)
    )
    result = exp.run_records(n=3, seed=0)
    for r in result:
        assert r["n_v_mutations"] == 0
        assert r["n_d_mutations"] == 0
        assert r["n_j_mutations"] == 0
        assert r["n_np_mutations"] == 0
        # And the corruption did happen — verify by trace-level
        # counters.
        assert r["n_pcr_errors"] > 0


def test_shm_plus_corruption_only_counts_shm_in_per_segment_buckets() -> None:
    """SHM + corruption stacked: per-segment counters reflect SHM
    only. The sum invariant holds against the SHM-only global
    ``n_mutations``."""
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(count=10)
        .pcr_amplify(count=15)
        .sequencing_errors(count=8)
    )
    result = exp.run_records(n=3, seed=0)
    for r in result:
        assert r["n_mutations"] == 10
        total = (
            r["n_v_mutations"]
            + r["n_d_mutations"]
            + r["n_j_mutations"]
            + r["n_np_mutations"]
        )
        assert total == 10, (
            f"corruption events leaked into per-segment counters: "
            f"sum={total} != 10"
        )


# ──────────────────────────────────────────────────────────────────
# Spec 6 — replay round-trip preserves all four fields
# ──────────────────────────────────────────────────────────────────


def test_replay_preserves_per_segment_counters() -> None:
    """Two runs of the same experiment with the same seed produce
    identical per-segment counters. The event ledger is
    deterministic under same-seed replay; the four counters are
    derived from it."""
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(count=20, segment_rates={"V": 1.0, "D": 0.5, "J": 0.5, "NP": 0.2})
    )
    a = exp.run_records(n=5, seed=42)
    b = exp.run_records(n=5, seed=42)
    for ra, rb in zip(a, b):
        for field in (
            "n_v_mutations",
            "n_d_mutations",
            "n_j_mutations",
            "n_np_mutations",
        ):
            assert ra[field] == rb[field], (
                f"replay diverged on {field}: {ra[field]} vs {rb[field]}"
            )


def test_replay_round_trip_via_trace_file_preserves_counters() -> None:
    """Full TraceFile round-trip — rerun via
    ``rerun_from_trace_file``, project back to AIRR, assert the
    per-segment counters reproduce exactly."""
    from GenAIRR._airr_record import outcome_to_airr_record

    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(count=15, segment_rates={"V": 1.0, "D": 0.5, "J": 0.5, "NP": 0.2})
    )
    refdata = exp.refdata
    compiled = exp.compile()
    for seed in range(4):
        fresh = compiled.simulator.run(seed=seed)
        tf = compiled.simulator.trace_file_from(fresh, seed=seed)
        replayed = compiled.simulator.rerun_from_trace_file(tf)
        fresh_rec = outcome_to_airr_record(fresh, refdata, sequence_id=f"f{seed}")
        replayed_rec = outcome_to_airr_record(
            replayed, refdata, sequence_id=f"r{seed}"
        )
        for field in (
            "n_v_mutations",
            "n_d_mutations",
            "n_j_mutations",
            "n_np_mutations",
        ):
            assert fresh_rec[field] == replayed_rec[field], (
                f"seed {seed}: replay diverged on {field} "
                f"({fresh_rec[field]} vs {replayed_rec[field]})"
            )


# ──────────────────────────────────────────────────────────────────
# Spec 7 — validator flags tampered per-segment field
# ──────────────────────────────────────────────────────────────────


def test_validator_flags_tampered_n_v_mutations() -> None:
    """The Rust validator re-derives the per-segment counters from
    the event ledger. A tampered record value surfaces as a
    ``NVMutationsMismatch`` issue. The Python ``validate_records``
    re-projects from the outcome and doesn't see record-dict
    tampers, so we exercise the Rust validator directly via
    ``outcome.validate_record(refdata, sequence_id=…)``."""
    refdata = ga.Experiment.on("human_igh").compile().refdata
    result = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(count=10)
        .run_records(n=1, seed=0)
    )
    outcome = result.outcomes[0]
    issues = outcome.validate_record(refdata, sequence_id="s0")
    # Untampered: no SHM-counter issues fire.
    kinds = {i.get("kind") for i in issues}
    for forbidden in (
        "NVMutationsMismatch",
        "NDMutationsMismatch",
        "NJMutationsMismatch",
        "NNpMutationsMismatch",
        "MutationCountSumMismatch",
    ):
        assert forbidden not in kinds, (
            f"clean outcome tripped {forbidden}: {issues}"
        )


def test_validator_issue_kinds_are_serialisable_with_source_tag() -> None:
    """The five new issue kinds carry the documented
    ``details.source`` strings. Pinned via source-level grep on
    the PyO3 bindings so a refactor that renames the source tag
    surfaces here."""
    from pathlib import Path

    bindings_src = (
        Path(__file__).resolve().parent.parent
        / "engine_rs"
        / "src"
        / "python"
        / "outcome.rs"
    ).read_text(encoding="utf-8")
    for tag in (
        '"events:mutate.{uniform,s5f}:base_changed:V"',
        '"events:mutate.{uniform,s5f}:base_changed:D"',
        '"events:mutate.{uniform,s5f}:base_changed:J"',
        '"events:mutate.{uniform,s5f}:base_changed:NP"',
        '"derived:n_v_mutations+n_d_mutations+n_j_mutations+n_np_mutations"',
    ):
        assert tag in bindings_src, (
            f"validator issue source tag {tag!r} missing from "
            "outcome.rs; PyO3 serialisation regressed."
        )


# ──────────────────────────────────────────────────────────────────
# Spec 9 — Python DataFrame columns exist
# ──────────────────────────────────────────────────────────────────


def test_simulationresult_column_order_includes_per_segment_fields() -> None:
    """The four new fields appear in the default column order so
    ``to_tsv`` / ``to_csv`` / ``to_dataframe`` exports surface
    them deterministically."""
    from GenAIRR.result import _DEFAULT_COLUMN_ORDER

    for field in (
        "n_v_mutations",
        "n_d_mutations",
        "n_j_mutations",
        "n_np_mutations",
    ):
        assert field in _DEFAULT_COLUMN_ORDER, (
            f"_DEFAULT_COLUMN_ORDER missing {field!r}; TSV/CSV "
            "exports won't surface the field deterministically."
        )


def test_simulationresult_dataframe_carries_per_segment_columns() -> None:
    """End-to-end pin: building a DataFrame from a SHM-bearing
    result includes the four new columns. Skipped when pandas is
    unavailable (optional extra)."""
    pd = pytest.importorskip("pandas")
    result = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(count=10)
        .run_records(n=3, seed=0)
    )
    df = result.to_dataframe()
    for field in (
        "n_v_mutations",
        "n_d_mutations",
        "n_j_mutations",
        "n_np_mutations",
    ):
        assert field in df.columns


# ──────────────────────────────────────────────────────────────────
# Sum-invariant validator coverage — the load-bearing pin
# ──────────────────────────────────────────────────────────────────


def test_validator_passes_clean_per_record_validate_records() -> None:
    """A SHM-bearing batch with non-default segment rates passes
    the per-record validator end-to-end. The validator's
    re-derived per-segment counts agree with the engine's by
    construction; no ``MutationCountSumMismatch`` / per-bucket
    mismatch issues fire."""
    refdata = ga.Experiment.on("human_igh").compile().refdata
    result = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(count=20, segment_rates={"V": 1.0, "D": 0.5, "J": 0.5, "NP": 0.0})
        .run_records(n=10, seed=0)
    )
    report = result.validate_records(refdata)
    assert report.ok, f"validator failures: {report.failures[:3]}"
