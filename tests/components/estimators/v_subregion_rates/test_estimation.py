"""End-to-end implementation tests for **Slice B:
`v_subregion_rates`** on `Experiment.mutate(...)`.

Slice B layers per-V-subregion SHM rate scalars on top of the
existing per-segment rates. The site's final weight is

```
weight = base_model_weight × segment_rate(segment) × v_subregion_rate(subregion)
```

for V sites (when the assigned V allele carries IMGT subregion
annotations); non-V sites and unannotated V alleles receive
factor `1.0` for the subregion component so the kwarg composes
cleanly with mixed cartridges.

This file pins the user-visible biology surface end-to-end:

- DSL validation (unknown key, type errors, NaN/inf/bool/negative,
  all-zero, unsatisfiable against unannotated cartridge).
- Alias expansion (`FWR` / `CDR`) with explicit-label override.
- `FWR`-only and `CDR`-only configurations work.
- Zero-rate subregion drops V-CDR1/CDR2 mutations on annotated V
  alleles (the canonical "turn off SHM in CDRs" test).
- Non-V sites are unaffected by the kwarg.
- Both `model="uniform"` and `model="s5f"` accept the kwarg.
- `productive_only()` still passes triad invariants under
  non-default rates.
- Replay with matching rates succeeds; replay with mismatched
  rates fails at the plan-signature gate (Slice A surface).
- Manifest advertises the new rate-support block.
- CDR/FR mutation counter AIRR fields stay absent (Slice C
  territory).
"""
from __future__ import annotations

import copy
import inspect
import json
import math

import pytest

import GenAIRR as ga
from GenAIRR.experiment import _validate_v_subregion_rates


# ──────────────────────────────────────────────────────────────────
# Spec 1 — DSL validation
# ──────────────────────────────────────────────────────────────────


def test_dsl_rejects_unknown_label() -> None:
    with pytest.raises(ValueError, match="unknown label"):
        ga.Experiment.on("human_igh").recombine().mutate(
            model="s5f", rate=0.03, v_subregion_rates={"CDR3": 1.0}
        )


def test_dsl_rejects_negative_value() -> None:
    with pytest.raises(ValueError, match="non-negative"):
        ga.Experiment.on("human_igh").recombine().mutate(
            model="s5f", rate=0.03, v_subregion_rates={"CDR1": -0.5}
        )


def test_dsl_rejects_bool_value() -> None:
    # ``bool`` is a subclass of int in Python; the validator must
    # reject it explicitly so a user typo like {"CDR1": True} doesn't
    # silently become "{CDR1: 1.0}".
    with pytest.raises(TypeError, match="bool"):
        ga.Experiment.on("human_igh").recombine().mutate(
            model="s5f", rate=0.03, v_subregion_rates={"CDR1": True}
        )


def test_dsl_rejects_nan_value() -> None:
    with pytest.raises(ValueError, match="NaN"):
        ga.Experiment.on("human_igh").recombine().mutate(
            model="s5f",
            rate=0.03,
            v_subregion_rates={"CDR1": float("nan")},
        )


def test_dsl_rejects_infinity_value() -> None:
    with pytest.raises(ValueError, match="finite"):
        ga.Experiment.on("human_igh").recombine().mutate(
            model="s5f",
            rate=0.03,
            v_subregion_rates={"CDR1": float("inf")},
        )


def test_dsl_rejects_all_zero_after_expansion() -> None:
    with pytest.raises(ValueError, match="at least one label must have a positive"):
        ga.Experiment.on("human_igh").recombine().mutate(
            model="s5f",
            rate=0.03,
            v_subregion_rates={"FWR": 0.0, "CDR": 0.0},
        )


def test_dsl_rejects_unsatisfiable_against_unannotated_cartridge() -> None:
    """If the cartridge has zero annotated V alleles the
    user's non-default rate vector would be a deterministic no-op.
    Reject at the DSL boundary."""
    import GenAIRR as ga
    cfg = copy.deepcopy(ga.HUMAN_IGH_OGRDB)
    # Strip subregions off every V allele.
    for alleles in cfg.v_alleles.values():
        for a in alleles:
            a.subregions = None
            a.gapped_seq = ""
    exp = ga.Experiment.on(cfg).recombine()
    with pytest.raises(ValueError, match="no V-subregion annotations"):
        exp.mutate(
            model="s5f",
            rate=0.03,
            v_subregion_rates={"CDR1": 2.0},
        )


# ──────────────────────────────────────────────────────────────────
# Spec 2 — Alias expansion + explicit override
# ──────────────────────────────────────────────────────────────────


def test_fwr_alias_expands_to_three_framework_labels() -> None:
    out = _validate_v_subregion_rates({"FWR": 0.5})
    # Order: FWR1, CDR1, FWR2, CDR2, FWR3.
    assert out == (0.5, 1.0, 0.5, 1.0, 0.5)


def test_cdr_alias_expands_to_two_cdr_labels() -> None:
    out = _validate_v_subregion_rates({"CDR": 3.0})
    assert out == (1.0, 3.0, 1.0, 3.0, 1.0)


def test_explicit_label_overrides_alias() -> None:
    out = _validate_v_subregion_rates({"FWR": 0.5, "FWR2": 2.0})
    assert out == (0.5, 1.0, 2.0, 1.0, 0.5)


def test_alias_plus_other_explicit_keys_compose() -> None:
    out = _validate_v_subregion_rates({"CDR": 3.0, "FWR1": 0.1})
    assert out == (0.1, 3.0, 1.0, 3.0, 1.0)


def test_default_and_omitted_kwarg_are_equivalent() -> None:
    """`None` and `{}` both resolve to the flat default."""
    assert _validate_v_subregion_rates(None) == (1.0, 1.0, 1.0, 1.0, 1.0)
    assert _validate_v_subregion_rates({}) == (1.0, 1.0, 1.0, 1.0, 1.0)


def test_explicit_all_ones_equals_default() -> None:
    """The behavioural-equivalence guarantee at the DSL layer:
    the explicit-all-ones dict and the default produce equal
    tuples, which is what makes their plan signatures collide
    via the `is_default()` short-circuit in
    `VSubregionRateWeights`."""
    explicit = _validate_v_subregion_rates(
        {"FWR1": 1.0, "CDR1": 1.0, "FWR2": 1.0, "CDR2": 1.0, "FWR3": 1.0}
    )
    assert explicit == _validate_v_subregion_rates(None)


# ──────────────────────────────────────────────────────────────────
# Spec 3 — FWR-only / CDR-only configurations work
# ──────────────────────────────────────────────────────────────────


def test_fwr_only_configuration_runs_and_validates() -> None:
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .productive_only()
        .mutate(model="s5f", rate=0.03, v_subregion_rates={"FWR": 0.5})
    )
    refdata = exp.refdata
    result = exp.run_records(n=20, seed=4242)
    report = result.validate_records(refdata)
    assert report, f"validation failed: {report.summary()}"


def test_cdr_only_configuration_runs_and_validates() -> None:
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .productive_only()
        .mutate(model="s5f", rate=0.03, v_subregion_rates={"CDR": 3.0})
    )
    refdata = exp.refdata
    result = exp.run_records(n=20, seed=4242)
    report = result.validate_records(refdata)
    assert report, f"validation failed: {report.summary()}"


# ──────────────────────────────────────────────────────────────────
# Spec 4 — Zero-rate subregion prevents mutations in that region
# ──────────────────────────────────────────────────────────────────


def _classify_v_region_mutations(mut_rec, baseline_rec, refdata) -> dict:
    """Diff a mutated record's sequence against a (same-seed)
    baseline record's sequence; classify each differing position
    by V subregion when it lands in the V region.

    Same baseline-diff strategy as the segment_rates implementation
    tests: run two parallel experiments on identical seeds with
    and without the mutate pass; positions where the sequences
    diverge are SHM hits. (S5F's RNG-consumption shape under
    permissive empty-support handling makes this approximation
    not perfectly clean, but for the rate-vector zeroing tests
    we only care about COUNTS in CDR vs FWR — exact byte
    positions are unnecessary.)

    Returns a dict with keys ``FWR1`` / ``CDR1`` / ``FWR2`` /
    ``CDR2`` / ``FWR3`` (annotated V hits) plus ``unannotated_v``
    (V hits outside the annotation span — e.g. V-side CDR3) and
    ``non_v`` (D / NP / J hits).
    """
    counts = {
        "FWR1": 0,
        "CDR1": 0,
        "FWR2": 0,
        "CDR2": 0,
        "FWR3": 0,
        "unannotated_v": 0,
        "non_v": 0,
    }
    base_seq = baseline_rec["sequence"].upper()
    mut_seq = mut_rec["sequence"].upper()
    # Lengths can differ under indel passes, but neither test
    # configures indels — so we expect equal lengths.
    if len(base_seq) != len(mut_seq):
        return counts
    v_end = mut_rec.get("v_sequence_end")
    if v_end is None:
        return counts
    v_end = int(v_end)

    # Subregion intervals for the assigned V allele (use the
    # mutated record's v_call — recombination is shared with the
    # baseline so the assignment matches).
    sub_intervals: list = []
    v_call = (mut_rec.get("v_call") or "").split(",")[0]
    if v_call:
        for v_id in range(refdata.v_pool_size()):
            a = refdata.v_allele(v_id)
            if a.name == v_call:
                sub_intervals = list(a.subregions)
                break

    # Trim_5: number of bases trimmed off the 5' end of the V
    # allele before assembly. AIRR ``v_germline_start`` is the
    # 1-based germline position the V region starts at; the
    # GenAIRR engine writes it 0-based here, but the field can
    # be either — accept either by treating it as the smaller of
    # the two. We use the explicit ``v_alignment_start`` /
    # ``v_germline_start`` field if present (zero-based).
    germ_start = mut_rec.get("v_germline_start")
    if germ_start is None:
        # Defensive — without it we can't map allele-local.
        # Conservatively classify every V-region diff as
        # ``unannotated_v``.
        for i in range(v_end):
            if base_seq[i] != mut_seq[i]:
                counts["unannotated_v"] += 1
        for i in range(v_end, len(mut_seq)):
            if base_seq[i] != mut_seq[i]:
                counts["non_v"] += 1
        return counts
    trim_5 = int(germ_start)

    for i in range(v_end):
        if base_seq[i] == mut_seq[i]:
            continue
        allele_local = i + trim_5
        landed = None
        for label, s, e in sub_intervals:
            if s <= allele_local < e:
                landed = label
                break
        if landed is None:
            counts["unannotated_v"] += 1
        else:
            counts[landed] += 1
    for i in range(v_end, len(mut_seq)):
        if base_seq[i] != mut_seq[i]:
            counts["non_v"] += 1
    return counts


def _run_baseline_and_mutated(rate_kwargs: dict, n: int = 60, seed: int = 7777):
    """Helper: run a no-mutation baseline and a mutated experiment
    on the same seed; return ``(baseline_records, mutated_records,
    refdata)`` aligned by index."""
    baseline = (
        ga.Experiment.on("human_igh")
        .recombine()
        .productive_only()
        .run_records(n=n, seed=seed)
    )
    mut_exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .productive_only()
        .mutate(model="s5f", rate=0.05, **rate_kwargs)
    )
    refdata = mut_exp.refdata
    mutated = mut_exp.run_records(n=n, seed=seed)
    return baseline, mutated, refdata


def test_zero_rate_cdr_prevents_cdr_mutations() -> None:
    """Setting ``CDR: 0`` drops every site inside CDR1 / CDR2 on
    annotated V alleles from proposal support; the empirical CDR
    mutation count across many records must be 0.

    We diff each mutated record against a no-mutation baseline at
    the same seed to recover the per-position SHM hits, then
    classify each V-region hit by IMGT subregion."""
    seed = 7777
    n_records = 60
    # Baseline diff for default rates — confirms the test isn't
    # trivially passing because the baseline has no CDR mutations
    # at all.
    base_default, mut_default, refdata = _run_baseline_and_mutated(
        rate_kwargs={}, n=n_records, seed=seed
    )
    baseline_cdr = 0
    for base_rec, mut_rec in zip(base_default, mut_default):
        bucket = _classify_v_region_mutations(mut_rec, base_rec, refdata)
        baseline_cdr += bucket["CDR1"] + bucket["CDR2"]
    assert baseline_cdr > 0, (
        "default-rate baseline produced no CDR1/CDR2 mutations; "
        "the zero-rate test would be vacuously passing"
    )

    # CDR-zeroed.
    base_off, mut_off, refdata = _run_baseline_and_mutated(
        rate_kwargs={"v_subregion_rates": {"CDR": 0.0}},
        n=n_records,
        seed=seed,
    )
    cdr_hits = 0
    fwr_hits = 0
    for base_rec, mut_rec in zip(base_off, mut_off):
        bucket = _classify_v_region_mutations(mut_rec, base_rec, refdata)
        cdr_hits += bucket["CDR1"] + bucket["CDR2"]
        fwr_hits += bucket["FWR1"] + bucket["FWR2"] + bucket["FWR3"]
    assert cdr_hits == 0, (
        f"CDR-zeroed run produced {cdr_hits} CDR mutations; "
        "expected 0 (zero-rate subregion drops sites from support "
        "before contract admissibility)"
    )
    assert fwr_hits > 0, (
        "zero-CDR run produced no FWR mutations either; "
        "rate-vector path is broken"
    )


def test_zero_rate_fwr_prevents_fwr_mutations() -> None:
    """Symmetric to the CDR test — zero out FWR; only CDR
    mutations should remain on annotated V alleles."""
    base_off, mut_off, refdata = _run_baseline_and_mutated(
        rate_kwargs={"v_subregion_rates": {"FWR": 0.0}},
        n=60,
        seed=7777,
    )
    fwr_hits = 0
    cdr_hits = 0
    for base_rec, mut_rec in zip(base_off, mut_off):
        bucket = _classify_v_region_mutations(mut_rec, base_rec, refdata)
        fwr_hits += bucket["FWR1"] + bucket["FWR2"] + bucket["FWR3"]
        cdr_hits += bucket["CDR1"] + bucket["CDR2"]
    assert fwr_hits == 0, (
        f"FWR-zeroed run produced {fwr_hits} FWR mutations; expected 0"
    )
    assert cdr_hits > 0, (
        "zero-FWR run produced no CDR mutations either; rate-vector "
        "path is broken"
    )


# ──────────────────────────────────────────────────────────────────
# Spec 5 — Non-V sites ignore the v_subregion_rates kwarg
# ──────────────────────────────────────────────────────────────────


def test_non_v_sites_ignore_v_subregion_rates() -> None:
    """A non-V SHM site (D, J, NP1, NP2) sees no V-subregion
    factor. We set every V subregion to zero EXCEPT CDR1 (which
    keeps the kwarg satisfiable — the all-zero-after-expansion
    rule rejects at validation time). The resulting V-region
    mutations must all land in CDR1; non-V mutations must still
    survive unchanged."""
    base_off, mut_off, refdata = _run_baseline_and_mutated(
        rate_kwargs={
            "v_subregion_rates": {
                "FWR": 0.0,
                "CDR1": 1.0,
                "CDR2": 0.0,
            }
        },
        n=40,
        seed=1234,
    )
    cdr1_hits = 0
    other_v_hits = 0
    non_v_hits = 0
    for base_rec, mut_rec in zip(base_off, mut_off):
        bucket = _classify_v_region_mutations(mut_rec, base_rec, refdata)
        cdr1_hits += bucket["CDR1"]
        other_v_hits += (
            bucket["FWR1"] + bucket["FWR2"] + bucket["FWR3"] + bucket["CDR2"]
        )
        non_v_hits += bucket["non_v"]
    assert other_v_hits == 0, (
        f"V-subregion zeroing leaked {other_v_hits} mutations into "
        "FWR / CDR2; expected 0"
    )
    # Non-V SHM should still be present — the kwarg didn't touch
    # D/J/NP, and there's no reason for those mutations to vanish.
    assert non_v_hits > 0, (
        "expected non-V (D/J/NP) SHM mutations to survive the V-only "
        "kwarg, but the records carry none — the kwarg is leaking "
        "into non-V sites"
    )


# ──────────────────────────────────────────────────────────────────
# Spec 6 — Works for both uniform and S5F
# ──────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("model", ["uniform", "s5f"])
def test_v_subregion_rates_works_for_both_models(model: str) -> None:
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .productive_only()
        .mutate(
            model=model,
            rate=0.03,
            v_subregion_rates={"CDR": 2.0, "FWR": 0.5},
        )
    )
    refdata = exp.refdata
    result = exp.run_records(n=20, seed=4242)
    report = result.validate_records(refdata)
    assert report, (
        f"validation failed for model={model!r}: {report.summary()}"
    )


# ──────────────────────────────────────────────────────────────────
# Spec 7 — productive_only triad invariants preserved
# ──────────────────────────────────────────────────────────────────


def test_productive_only_triad_preserved_under_v_subregion_rates() -> None:
    """productive_only() pins three invariants: productive
    junction frame, no stop codon in junction, anchor preserved.
    All must hold under non-default v_subregion_rates."""
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .productive_only()
        .mutate(
            model="s5f",
            rate=0.05,
            v_subregion_rates={"CDR": 4.0, "FWR": 0.25},
        )
    )
    refdata = exp.refdata
    result = exp.run_records(n=40, seed=0xBADC0DE)
    report = result.validate_records(refdata)
    assert report, (
        f"productive triad regressed under v_subregion_rates: "
        f"{report.summary()}"
    )
    # And every record carries productive=T per the bundle.
    for r in result.records:
        assert r.get("productive") in (True, "T"), (
            f"non-productive record under productive_only(): {r.get('sequence_id')}"
        )


# ──────────────────────────────────────────────────────────────────
# Spec 8 — Matching-rates replay succeeds
# ──────────────────────────────────────────────────────────────────


def test_replay_with_matching_v_subregion_rates_succeeds() -> None:
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(
            model="s5f",
            rate=0.03,
            v_subregion_rates={"CDR1": 2.0, "FWR2": 0.5},
        )
    )
    compiled = exp.compile()
    outcome = compiled.simulator.run(seed=4242)
    tf = compiled.simulator.trace_file_from(outcome, seed=4242)
    # Identical compile reproduces the trace exactly.
    compiled.simulator.replay_from_trace_file(tf)


# ──────────────────────────────────────────────────────────────────
# Spec 9 — Mismatched-rates replay fails at the plan-signature gate
# ──────────────────────────────────────────────────────────────────


def test_replay_with_mismatched_v_subregion_rates_fails() -> None:
    """Trace recorded under one V-subregion vector cannot replay
    against a plan with a different V-subregion vector — the
    plan-signature gate (Slice A) catches it before any choice is
    consumed."""
    exp_a = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(
            model="s5f",
            rate=0.03,
            v_subregion_rates={"CDR": 2.0},
        )
    )
    exp_b = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(
            model="s5f",
            rate=0.03,
            v_subregion_rates={"CDR": 3.0},
        )
    )
    c_a = exp_a.compile()
    c_b = exp_b.compile()
    outcome = c_a.simulator.run(seed=4242)
    tf = c_a.simulator.trace_file_from(outcome, seed=4242)
    with pytest.raises(ValueError, match="pass plan signature mismatch"):
        c_b.simulator.replay_from_trace_file(tf)


def test_default_v_subregion_rates_replay_byte_identical_to_pre_slice() -> None:
    """Plans with `v_subregion_rates=None` (or omitted) must
    produce the SAME plan signature as plans without the kwarg —
    no regression for legacy pipelines."""
    exp_a = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(model="s5f", rate=0.03)
    )
    exp_b = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(model="s5f", rate=0.03, v_subregion_rates=None)
    )
    c_a, c_b = exp_a.compile(), exp_b.compile()
    tf_a = c_a.simulator.trace_file_from(
        c_a.simulator.run(seed=42), seed=42
    )
    tf_b = c_b.simulator.trace_file_from(
        c_b.simulator.run(seed=42), seed=42
    )
    sig_a = json.loads(tf_a.to_json())["pass_plan_signature"]
    sig_b = json.loads(tf_b.to_json())["pass_plan_signature"]
    assert sig_a == sig_b


# ──────────────────────────────────────────────────────────────────
# Spec 10 — Manifest advertises support
# ──────────────────────────────────────────────────────────────────


def test_manifest_advertises_v_subregion_rate_support() -> None:
    m = ga.HUMAN_IGH_OGRDB.cartridge_manifest()
    sup = m["models"]["shm"]["v_subregion_rate_support"]
    assert sup["available"] is True
    assert sup["labels"] == ["FWR1", "CDR1", "FWR2", "CDR2", "FWR3"]
    assert sup["aliases"] == {
        "FWR": ["FWR1", "FWR2", "FWR3"],
        "CDR": ["CDR1", "CDR2"],
    }
    assert sup["default"] == {
        "FWR1": 1.0,
        "CDR1": 1.0,
        "FWR2": 1.0,
        "CDR2": 1.0,
        "FWR3": 1.0,
    }
    assert sup["in_plan_signature"] is True
    assert sup["in_content_hash"] is False


# ──────────────────────────────────────────────────────────────────
# Spec 11 — CDR/FR mutation counter fields are now present
# (the counters slice shipped after Slice B; full per-counter
# behaviour lives in ``test_v_subregion_mutation_counters_implementation.py``;
# here we just confirm the fields surface and the partition
# invariant composes with non-default rates).
# ──────────────────────────────────────────────────────────────────


def test_cdr_fwr_mutation_counter_fields_present_and_partition_holds() -> None:
    """The counters slice landed after Slice B. The six new
    V-subregion fields are present on every AIRR record and
    partition ``n_v_mutations`` cleanly — including when a
    non-default ``v_subregion_rates`` vector is in play. The
    two-bucket aggregates (``n_cdr_mutations`` /
    ``n_fwr_mutations``) stay deliberately absent."""
    result = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(
            model="s5f",
            rate=0.03,
            v_subregion_rates={"CDR": 2.0},
        )
        .run_records(n=5, seed=0)
    )
    six_fields = (
        "n_fwr1_mutations",
        "n_cdr1_mutations",
        "n_fwr2_mutations",
        "n_cdr2_mutations",
        "n_fwr3_mutations",
        "n_v_unannotated_mutations",
    )
    for r in result.records:
        for f in six_fields:
            assert f in r, (
                f"AIRR record missing {f!r}; counters slice regressed"
            )
        # Partition holds.
        assert sum(r[f] for f in six_fields) == r["n_v_mutations"]
    # Two-bucket aggregates explicitly not in v1.
    rec = result.records[0]
    for forbidden in ("n_cdr_mutations", "n_fwr_mutations"):
        assert forbidden not in rec, (
            f"AIRR record now carries the two-bucket aggregate "
            f"{forbidden!r}; the audit recommended five-label fields only"
        )


def test_mutate_signature_carries_v_subregion_rates_kwarg() -> None:
    """Source-level pin: the kwarg exists with the documented
    default."""
    sig = inspect.signature(ga.Experiment.mutate)
    assert "v_subregion_rates" in sig.parameters
    param = sig.parameters["v_subregion_rates"]
    assert param.default is None
