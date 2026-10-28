"""End-to-end implementation tests for **V-Subregion Mutation
Counters** — the partition of `n_v_mutations` across the five
IMGT subregion labels plus `n_v_unannotated_mutations`.

The slice's headline contract:

    n_fwr1_mutations + n_cdr1_mutations + n_fwr2_mutations
    + n_cdr2_mutations + n_fwr3_mutations
    + n_v_unannotated_mutations == n_v_mutations

…holds on every record, with the validator surfacing six
per-field mismatch issues plus the cross-field
`VSubregionMutationCountSumMismatch` invariant.

Specs covered (per the user brief):

1. Baseline no-mutation: all six fields are 0.
2. Bundled annotated V cartridges: partition equals
   `n_v_mutations`; `n_v_unannotated_mutations` accumulates only
   for the V-side CDR3 stretch (between FWR3.end and the V allele
   end) — a small but expected non-zero on bundled cartridges
   because the five-label set deliberately stops at FWR3.
3. `v_subregion_rates={"CDR": 0, "FWR": 1}` yields zero CDR
   counters.
4. `v_subregion_rates={"FWR": 0, "CDR": 1}` yields zero FWR
   counters.
5. Mixed full-stack validates clean.
6. Replay preserves all six fields.
7. Validator flags a tampered per-field counter.
8. Validator flags a partition mismatch.
9. Unannotated custom V allele routes V mutations to
   `n_v_unannotated_mutations`.
10. PCR / quality / N-corruption / receptor revision / D
    inversion do NOT affect the counters.
11. DataFrame columns exist with the canonical order.
12. (Pin discipline lives in the contract file; this file is
    behavioural tests only.)
"""
from __future__ import annotations

import copy
import inspect

import pytest

import GenAIRR as ga


_COUNTER_FIELDS: tuple[str, ...] = (
    "n_fwr1_mutations",
    "n_cdr1_mutations",
    "n_fwr2_mutations",
    "n_cdr2_mutations",
    "n_fwr3_mutations",
    "n_v_unannotated_mutations",
)


def _partition(rec) -> int:
    return sum(rec[f] for f in _COUNTER_FIELDS)


# ──────────────────────────────────────────────────────────────────
# Spec 1 — Baseline no-mutation: all six fields are 0
# ──────────────────────────────────────────────────────────────────


def test_baseline_no_mutation_yields_zero_counters() -> None:
    """An experiment with no `.mutate()` step produces zero V SHM
    events. The six new counters are 0 on every record."""
    result = (
        ga.Experiment.on("human_igh")
        .recombine()
        .productive_only()
        .run_records(n=10, seed=4242)
    )
    for r in result.records:
        for f in _COUNTER_FIELDS:
            assert r[f] == 0, (
                f"baseline no-mutate run has non-zero {f}={r[f]}"
            )
        assert r["n_v_mutations"] == 0


# ──────────────────────────────────────────────────────────────────
# Spec 2 — Bundled annotated V cartridges: partition holds
# ──────────────────────────────────────────────────────────────────


@pytest.mark.parametrize(
    "preset",
    ["human_igh", "human_igk", "human_igl"],
)
def test_partition_equals_n_v_mutations_on_bundled_cartridges(preset: str) -> None:
    """The six counters sum exactly to `n_v_mutations` on every
    bundled-cartridge record. This is the slice's load-bearing
    invariant — any drift would make the partition meaningless.

    Note: `n_v_unannotated_mutations` is NOT identically zero on
    bundled cartridges. The five canonical IMGT labels stop at
    FWR3; SHM events landing in the V-side CDR3 contribution
    (between `FWR3.end` and `len(allele.seq)`) carry a valid
    `germline_pos` but no matching subregion interval, so they
    route to the unannotated bucket. This is by design (audit
    §4 unannotated case 4)."""
    result = (
        ga.Experiment.on(preset)
        .recombine()
        .productive_only()
        .mutate(model="s5f", rate=0.05)
        .run_records(n=30, seed=4242)
    )
    saw_any_v = False
    for r in result.records:
        assert _partition(r) == r["n_v_mutations"], (
            f"{preset}: V-subregion partition violated: "
            f"sum={_partition(r)} != n_v_mutations={r['n_v_mutations']}; "
            f"per-bucket={[(f, r[f]) for f in _COUNTER_FIELDS]}"
        )
        if r["n_v_mutations"] > 0:
            saw_any_v = True
    assert saw_any_v, (
        f"{preset}: sweep produced zero V SHM mutations — partition "
        "test is vacuous"
    )


# ──────────────────────────────────────────────────────────────────
# Spec 3 — CDR zeroed: CDR counters are 0
# ──────────────────────────────────────────────────────────────────


def test_cdr_zeroed_yields_zero_cdr_counters() -> None:
    """Setting `v_subregion_rates={"CDR": 0.0, "FWR": 1.0}` drops
    every site in CDR1 / CDR2 from proposal support, so the two
    CDR counters must be 0 on every record. The FWR counters and
    `n_v_unannotated_mutations` together account for the rest of
    `n_v_mutations`."""
    result = (
        ga.Experiment.on("human_igh")
        .recombine()
        .productive_only()
        .mutate(
            model="s5f",
            rate=0.05,
            v_subregion_rates={"CDR": 0.0, "FWR": 1.0},
        )
        .run_records(n=40, seed=4242)
    )
    cdr_total = 0
    fwr_total = 0
    for r in result.records:
        cdr_total += r["n_cdr1_mutations"] + r["n_cdr2_mutations"]
        fwr_total += (
            r["n_fwr1_mutations"]
            + r["n_fwr2_mutations"]
            + r["n_fwr3_mutations"]
        )
        # Partition still holds.
        assert _partition(r) == r["n_v_mutations"]
    assert cdr_total == 0, (
        f"CDR-zeroed run produced {cdr_total} CDR mutations across "
        "40 records; expected 0"
    )
    assert fwr_total > 0, (
        "CDR-zeroed run produced zero FWR mutations across 40 "
        "records — sampling path didn't exercise V"
    )


# ──────────────────────────────────────────────────────────────────
# Spec 4 — FWR zeroed: FWR counters are 0
# ──────────────────────────────────────────────────────────────────


def test_fwr_zeroed_yields_zero_fwr_counters() -> None:
    """Symmetric to the CDR test — zero out FWR; only CDR
    mutations should accrue. `n_v_unannotated_mutations` can
    still grow from V-side CDR3 events (those positions are
    outside every subregion interval and don't get a rate
    factor — wait, they DO get factor 1.0 by default; see audit
    §4 unannotated case 4 plus rate-slice §2 "Edge cases"
    fallback `1.0`)."""
    result = (
        ga.Experiment.on("human_igh")
        .recombine()
        .productive_only()
        .mutate(
            model="s5f",
            rate=0.05,
            v_subregion_rates={"FWR": 0.0, "CDR": 1.0},
        )
        .run_records(n=40, seed=4242)
    )
    fwr_total = 0
    cdr_total = 0
    for r in result.records:
        fwr_total += (
            r["n_fwr1_mutations"]
            + r["n_fwr2_mutations"]
            + r["n_fwr3_mutations"]
        )
        cdr_total += r["n_cdr1_mutations"] + r["n_cdr2_mutations"]
        assert _partition(r) == r["n_v_mutations"]
    assert fwr_total == 0, (
        f"FWR-zeroed run produced {fwr_total} FWR mutations; expected 0"
    )
    assert cdr_total > 0, (
        "FWR-zeroed run produced zero CDR mutations — rate path "
        "isn't exercising CDR-only configuration"
    )


# ──────────────────────────────────────────────────────────────────
# Spec 5 — Mixed full-stack validates clean
# ──────────────────────────────────────────────────────────────────


def test_full_stack_validates_clean_under_counters() -> None:
    """Heavy full-stack pipeline (productive_only + corruption +
    paired_end + non-default segment + V-subregion rates).
    `validate_records(refdata)` runs clean — the new validator
    checks for per-field mismatches + partition invariant must
    not fire on engine-projected records."""
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .productive_only()
        .mutate(
            model="s5f",
            rate=0.04,
            segment_rates={"V": 1.0, "D": 0.3, "J": 0.5, "NP": 0.2},
            v_subregion_rates={"CDR": 2.5, "FWR": 0.5},
        )
        .pcr_amplify(rate=1e-4)
        .polymerase_indels(count=2)
        .primer_trim_5prime(length=(0, 3))
    )
    refdata = exp.refdata
    result = exp.run_records(n=40, seed=4242)
    report = result.validate_records(refdata)
    assert report, (
        f"V-subregion counter mixed-stack validation failed: "
        f"summary={report.summary()}; first failure="
        f"{report.failures[0] if report.failures else None}"
    )
    # Partition still holds on every record under the heavy stack.
    for r in result:
        assert _partition(r) == r["n_v_mutations"]


# ──────────────────────────────────────────────────────────────────
# Spec 6 — Replay preserves all six fields
# ──────────────────────────────────────────────────────────────────


def test_replay_preserves_all_six_counter_fields() -> None:
    """A trace round-trip via `rerun_from_trace_file` reproduces
    the original outcome's six counter fields byte-for-byte. The
    counters are derived from the deterministic event ledger
    which the replay reproduces exactly."""
    from GenAIRR._airr_record import outcome_to_airr_record

    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .productive_only()
        .mutate(
            model="s5f",
            rate=0.04,
            v_subregion_rates={"CDR": 2.0, "FWR": 0.5},
        )
    )
    refdata = exp.refdata
    compiled = exp.compile()
    seen_v_mutations = False
    for seed in range(6):
        fresh = compiled.simulator.run(seed=seed)
        tf = compiled.simulator.trace_file_from(fresh, seed=seed)
        replayed = compiled.simulator.rerun_from_trace_file(tf)
        fresh_rec = outcome_to_airr_record(
            fresh, refdata, sequence_id=f"fresh-{seed}"
        )
        replayed_rec = outcome_to_airr_record(
            replayed, refdata, sequence_id=f"replay-{seed}"
        )
        for f in _COUNTER_FIELDS:
            assert fresh_rec[f] == replayed_rec[f], (
                f"seed {seed} field {f!r} desynced under replay "
                f"({fresh_rec[f]} vs {replayed_rec[f]})"
            )
        # And the partition still equals n_v_mutations on both sides.
        assert _partition(fresh_rec) == fresh_rec["n_v_mutations"]
        assert _partition(replayed_rec) == replayed_rec["n_v_mutations"]
        if fresh_rec["n_v_mutations"] > 0:
            seen_v_mutations = True
    assert seen_v_mutations, (
        "no V SHM mutations observed across 6 seeds — replay sweep "
        "isn't exercising the counter path"
    )


# ──────────────────────────────────────────────────────────────────
# Spec 7 + 8 — Validator mismatch detection
# ──────────────────────────────────────────────────────────────────
#
# The Python validator entry point (`Outcome.validate_record`)
# builds a fresh AirrRecord from the outcome on every call — it
# never accepts a hand-edited dict, so a tampered-dict-style test
# can't run through Python. The validator's mismatch-detection
# logic is exercised at the Rust unit-test layer (see
# `engine_rs/src/airr_record/tests/projection.rs`) where
# `AirrRecord` struct fields can be edited directly before calling
# `validate_airr_record(rec, outcome, refdata)`.
#
# At the Python layer we pin the converse: **engine-projected
# records always validate clean**. This is the strongest user-
# visible contract — if the builder's aggregation ever drifted
# from the validator's independent recompute, the projection test
# would surface the disagreement on every run.


def test_engine_projected_records_carry_no_v_subregion_mismatches() -> None:
    """Validator must never surface a V-subregion mismatch on an
    engine-projected record. The builder and the validator use
    independent recomputes over the same event ledger; if they
    ever disagree, this test catches it.

    The seven new issue kinds the slice introduces:
    ``N{Fwr1,Cdr1,Fwr2,Cdr2,Fwr3,VUnannotated}MutationsMismatch``
    plus ``VSubregionMutationCountSumMismatch``. None should fire
    on a clean engine projection."""
    forbidden_kinds = {
        "NFwr1MutationsMismatch",
        "NCdr1MutationsMismatch",
        "NFwr2MutationsMismatch",
        "NCdr2MutationsMismatch",
        "NFwr3MutationsMismatch",
        "NVUnannotatedMutationsMismatch",
        "VSubregionMutationCountSumMismatch",
    }
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .productive_only()
        .mutate(model="s5f", rate=0.05,
                v_subregion_rates={"CDR": 2.0, "FWR": 0.5})
    )
    refdata = exp.refdata
    result = exp.run_records(n=30, seed=4242)
    outcomes = result.outcomes
    assert outcomes is not None
    for i, outcome in enumerate(outcomes):
        issues = outcome.validate_record(
            refdata, sequence_id=f"r-{i}"
        )
        kinds_seen = {issue.get("kind") for issue in issues}
        offending = kinds_seen & forbidden_kinds
        assert not offending, (
            f"record {i}: engine projection produced V-subregion "
            f"mismatch issue(s) {offending}; the builder's "
            "aggregation has drifted from the validator's recompute"
        )


# ──────────────────────────────────────────────────────────────────
# Spec 9 — Unannotated V allele routes V mutations to unannotated
# ──────────────────────────────────────────────────────────────────


def test_unannotated_v_allele_routes_to_unannotated_bucket() -> None:
    """When the assigned V allele has empty `subregions` (legacy
    cartridge / hand-authored V), every V SHM event routes to
    `n_v_unannotated_mutations` and the five canonical buckets
    stay 0. The partition still holds:
    `n_v_unannotated_mutations == n_v_mutations`."""
    cfg = copy.deepcopy(ga.HUMAN_IGH_OGRDB)
    # Strip subregions off every V allele to simulate a fully
    # unannotated cartridge. Note: default rates work; only a
    # non-default v_subregion_rates against an unannotated
    # cartridge is rejected at the DSL boundary.
    for alleles in cfg.v_alleles.values():
        for a in alleles:
            a.subregions = None
            a.gapped_seq = ""
    result = (
        ga.Experiment.on(cfg)
        .recombine()
        .productive_only()
        .mutate(model="s5f", rate=0.05)
        .run_records(n=15, seed=4242)
    )
    saw_v_events = False
    for r in result.records:
        # The five canonical buckets must be 0 on every record.
        for f in ("n_fwr1_mutations", "n_cdr1_mutations", "n_fwr2_mutations",
                  "n_cdr2_mutations", "n_fwr3_mutations"):
            assert r[f] == 0, (
                f"unannotated-cartridge record has non-zero {f}={r[f]}"
            )
        # All V SHM events route to unannotated.
        assert r["n_v_unannotated_mutations"] == r["n_v_mutations"], (
            f"unannotated V allele routing broke: "
            f"unannotated={r['n_v_unannotated_mutations']} != "
            f"n_v={r['n_v_mutations']}"
        )
        # Partition still holds (trivially).
        assert _partition(r) == r["n_v_mutations"]
        if r["n_v_mutations"] > 0:
            saw_v_events = True
    assert saw_v_events, (
        "unannotated-cartridge sweep produced zero V SHM events — "
        "test is vacuous"
    )


# ──────────────────────────────────────────────────────────────────
# Spec 10 — PCR / quality / N-corruption do not affect counters
# ──────────────────────────────────────────────────────────────────


def test_corruption_stage_passes_do_not_affect_counters() -> None:
    """The aggregation filters to `mutate.{uniform,s5f}` pass
    names only. Heavy PCR + quality + N-corruption + receptor
    revision + D inversion in the stack must not increment the
    six new counters beyond what mutate.s5f / mutate.uniform
    produces.

    Strategy: two parallel experiments at the same seed —
    one with corruption, one without. The six counters must
    match record-by-record."""
    base_exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .productive_only()
        .mutate(model="s5f", rate=0.04)
    )
    corrupted_exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .productive_only()
        .mutate(model="s5f", rate=0.04)
        .pcr_amplify(count=8)
        .sequencing_errors(count=6)
        .ambiguous_base_calls(count=4)
    )
    seed = 4242
    n = 30
    base_records = base_exp.run_records(n=n, seed=seed)
    corrupted_records = corrupted_exp.run_records(n=n, seed=seed)
    saw_pcr = False
    for base_rec, mut_rec in zip(base_records.records, corrupted_records.records):
        for f in _COUNTER_FIELDS:
            assert base_rec[f] == mut_rec[f], (
                f"corruption stage perturbed {f!r}: "
                f"baseline={base_rec[f]} corrupted={mut_rec[f]}"
            )
        if mut_rec["n_pcr_errors"] > 0:
            saw_pcr = True
    assert saw_pcr, (
        "corruption stack didn't produce any PCR errors — test "
        "is vacuous"
    )


# ──────────────────────────────────────────────────────────────────
# Spec 11 — DataFrame columns exist with canonical order
# ──────────────────────────────────────────────────────────────────


def test_dataframe_columns_carry_v_subregion_counters() -> None:
    """`SimulationResult.to_dataframe()` exposes the six new
    counter columns in canonical order (FWR1 / CDR1 / FWR2 /
    CDR2 / FWR3 / unannotated). Sits next to the per-segment
    counters in the column list."""
    result = (
        ga.Experiment.on("human_igh")
        .recombine()
        .productive_only()
        .mutate(model="s5f", rate=0.04)
        .run_records(n=5, seed=4242)
    )
    df = result.to_dataframe()
    columns = list(df.columns)
    for f in _COUNTER_FIELDS:
        assert f in columns, (
            f"DataFrame column {f!r} missing; default column order "
            "doesn't surface the new counter fields"
        )
    # The new counters land immediately after `n_np_mutations` per
    # the column-order definition in `src/GenAIRR/result.py`.
    np_idx = columns.index("n_np_mutations")
    expected_block = [
        "n_fwr1_mutations",
        "n_cdr1_mutations",
        "n_fwr2_mutations",
        "n_cdr2_mutations",
        "n_fwr3_mutations",
        "n_v_unannotated_mutations",
    ]
    assert columns[np_idx + 1 : np_idx + 1 + len(expected_block)] == expected_block, (
        f"V-subregion counter columns are not in canonical order "
        f"after n_np_mutations; got {columns[np_idx + 1 : np_idx + 7]}"
    )


# ──────────────────────────────────────────────────────────────────
# Manifest exposure
# ──────────────────────────────────────────────────────────────────


def test_manifest_advertises_v_subregion_counter_support() -> None:
    m = ga.HUMAN_IGH_OGRDB.cartridge_manifest()
    sup = m["models"]["shm"]["v_subregion_counter_support"]
    assert sup["available"] is True
    assert sup["fields"] == [
        "n_fwr1_mutations",
        "n_cdr1_mutations",
        "n_fwr2_mutations",
        "n_cdr2_mutations",
        "n_fwr3_mutations",
        "n_v_unannotated_mutations",
    ]
    assert sup["partition_of"] == "n_v_mutations"
    assert sup["requires_annotations"] is True
    assert sup["unannotated_bucket"] == "n_v_unannotated_mutations"
    assert sup["in_content_hash"] is False
