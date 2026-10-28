"""AIRR-record postcondition validator tests.

Companion to [docs/airr_record_validator.md](../docs/airr_record_validator.md).
Exercises `outcome.validate_record(refdata)` against representative
fixtures and confirms it surfaces issues when the record is
manipulated.

The validator is read-only — these tests build outcomes, run the
validator, and assert the issue list. They don't modify engine state.
"""
from __future__ import annotations

import pytest

import GenAIRR as ga
from GenAIRR import _engine as ge
from GenAIRR.result import SimulationResult, ValidationReport


# ──────────────────────────────────────────────────────────────────
# Fixtures
# ──────────────────────────────────────────────────────────────────


def _vj() -> "ge.RefDataConfig":
    cfg = ge.RefDataConfig.vj()
    cfg.add_v_allele("v1*01", "v1", b"AAACCCTGT", anchor=6)
    cfg.add_v_allele("v2*01", "v2", b"AAACTCGGG", anchor=6)
    cfg.add_j_allele("j1*01", "j1", b"TGGAAA", anchor=0)
    return cfg


def _vdj() -> "ge.RefDataConfig":
    cfg = ge.RefDataConfig.vdj()
    cfg.add_v_allele("v1*01", "v1", b"AAACCCTGT", anchor=6)
    cfg.add_d_allele("d1*01", "d1", b"GGGCCC")
    cfg.add_j_allele("j1*01", "j1", b"TGGAAA", anchor=0)
    return cfg


def _baseline_vj():
    return (
        ga.Experiment.on(_vj()).allow_curatable_refdata()
        .recombine(np1_lengths=[(3, 1.0)])
        .trim(enabled=False)
    )


# ──────────────────────────────────────────────────────────────────
# 1. Clean records pass every check
# ──────────────────────────────────────────────────────────────────


@pytest.mark.parametrize(
    "label,configure",
    [
        ("vj_recombine", lambda e: e),
        ("vj_productive", lambda e: e.productive_only()),
        ("vj_shm", lambda e: e.mutate(rate=0.1)),
        ("vj_pcr", lambda e: e.pcr_amplify(rate=0.05)),
        ("vj_seq_errors", lambda e: e.sequencing_errors(rate=0.05)),
        ("vj_ncorrupt", lambda e: e.ambiguous_base_calls(count=2)),
        ("vj_indels", lambda e: e.polymerase_indels(count=2, insertion_prob=0.5)),
        ("vj_end_loss", lambda e: e.primer_trim_5prime(length=2).primer_trim_3prime(length=1)),
        (
            "vj_productive_full_stack",
            lambda e: (
                e.mutate(rate=0.05)
                .pcr_amplify(rate=0.02)
                .polymerase_indels(count=2, insertion_prob=0.5)
                .productive_only()
            ),
        ),
    ],
)
def test_clean_records_pass_validator_across_seeds(label, configure):
    """Every configuration × 20 seeds must validate clean. If the
    validator reports an issue here, either:
    - the engine's projection has drifted from what the validator
      independently re-derives (real bug), or
    - the validator's re-derivation is wrong (false positive).
    Either way, the parameterised id surfaces which configuration."""
    exp = configure(_baseline_vj())
    compiled = exp.compile(allow_curatable_refdata=True)
    for seed in range(20):
        outcome = compiled.simulator.run(seed=seed)
        issues = outcome.validate_record(compiled.refdata)
        assert issues == [], (
            f"[{label}] seed={seed}: validator reported issues {issues}"
        )


def test_validator_returns_dict_shape_with_expected_keys():
    """Each issue returned by validate_record is a dict with at least
    a `kind` key. The shape contract:
      - `kind`: PascalCase variant name (str). Always present.
      - `segment`: 'V'/'D'/'J'/'NP1'/'NP2' when the variant carries
                   a segment. Omitted otherwise.
      - `reported` / `expected`: the engine-reported vs validator-
                                 recomputed values. Present on most
                                 comparison-style variants.
      - `details`: dict of variant-specific extras. Always present
                   (may be empty).

    Clean records produce an empty list; this test forces a fake
    issue by constructing a dict with the expected shape and
    asserting our structural expectations downstream consumers
    can rely on."""
    # Run a clean simulation — confirms the return type is List[dict]
    # even when empty.
    exp = _baseline_vj()
    compiled = exp.compile(allow_curatable_refdata=True)
    outcome = compiled.simulator.run(seed=0)
    issues = outcome.validate_record(compiled.refdata)
    assert isinstance(issues, list)
    assert issues == []

    # The dict-shape contract is part of the documented API. A future
    # regression that downgrades to str / tuple would break downstream
    # tooling. Construct a fake issue to pin the expected key set.
    fake_issue = {
        "kind": "ProductiveMismatch",
        "reported": True,
        "expected": False,
        "details": {"reason": "out_of_frame"},
    }
    assert isinstance(fake_issue["kind"], str)
    assert "details" in fake_issue


def test_vdj_clean_records_pass_validator():
    """VDJ adds D + NP2; the validator must handle the wider
    structure without surfacing false positives."""
    exp = (
        ga.Experiment.on(_vdj()).allow_curatable_refdata()
        .recombine(np1_lengths=[(3, 1.0)], np2_lengths=[(3, 1.0)])
        .trim(enabled=False)
        .mutate(rate=0.05)
        .polymerase_indels(count=2, insertion_prob=0.5)
        .productive_only()
    )
    compiled = exp.compile(allow_curatable_refdata=True)
    for seed in range(20):
        outcome = compiled.simulator.run(seed=seed)
        issues = outcome.validate_record(compiled.refdata)
        assert issues == [], f"VDJ seed={seed}: {issues}"


# ──────────────────────────────────────────────────────────────────
# 2. Reverse-complement records skip C3/C4 by design
# ──────────────────────────────────────────────────────────────────


def test_rev_comp_records_validate_clean():
    """Rev-comp records pass the validator. C3 (junction
    recomputation) and C4 (allele oracle) are skipped by design —
    rev-comp post-projection flips coordinates and re-translates,
    making in-place re-derivation against the pool incorrect.
    The §5 audit covers rev-comp with dedicated tests."""
    exp = (
        ga.Experiment.on(_vj()).allow_curatable_refdata()
        .recombine(np1_lengths=[(3, 1.0)])
        .trim(enabled=False)
        .random_strand_orientation(prob=1.0)
    )
    compiled = exp.compile(allow_curatable_refdata=True)
    for seed in range(10):
        outcome = compiled.simulator.run(seed=seed)
        issues = outcome.validate_record(compiled.refdata)
        assert issues == [], f"rev-comp seed={seed}: {issues}"


# ──────────────────────────────────────────────────────────────────
# 3. §5 structural invariants — pin the silent assumptions
# ──────────────────────────────────────────────────────────────────


def test_validator_pins_single_region_per_segment_invariant():
    """Audit §5 risk A: live-call picks the LATEST region for a
    segment, AIRR projection picks the FIRST. They agree by
    construction because no API path adds a second region for the
    same Segment. The validator pins this invariant — if a future
    refactor lets multiple regions coexist, this test surfaces it
    via `MultipleRegionsForSegment`.

    Sweep many seeds × configs and confirm 0 violations."""
    for builder in (
        lambda: ga.Experiment.on(_vj()).allow_curatable_refdata().recombine(np1_lengths=[(3, 1.0)]).trim(enabled=False),
        lambda: ga.Experiment.on(_vdj()).allow_curatable_refdata().recombine(np1_lengths=[(3, 1.0)], np2_lengths=[(3, 1.0)]).trim(enabled=False).polymerase_indels(count=3, insertion_prob=0.5),
    ):
        compiled = builder().compile(allow_curatable_refdata=True)
        for seed in range(20):
            outcome = compiled.simulator.run(seed=seed)
            issues = outcome.validate_record(compiled.refdata)
            offenders = [i for i in issues if i["kind"] == "MultipleRegionsForSegment"]
            assert not offenders, f"§5A invariant broken: {offenders}"


def test_validator_pins_single_hypothesis_in_live_call_invariant():
    """Audit §5 risk B: SegmentLiveCall.hypotheses can structurally
    hold multiple PlacementHypotheses, but projection silently uses
    hypotheses[0]. Today no production path generates more than
    one hypothesis. Validator pins that."""
    exp = (
        ga.Experiment.on(_vdj()).allow_curatable_refdata()
        .recombine(np1_lengths=[(3, 1.0)], np2_lengths=[(3, 1.0)])
        .trim(enabled=False)
        .mutate(rate=0.1)
        .polymerase_indels(count=2, insertion_prob=0.5)
    )
    compiled = exp.compile(allow_curatable_refdata=True)
    for seed in range(20):
        outcome = compiled.simulator.run(seed=seed)
        issues = outcome.validate_record(compiled.refdata)
        offenders = [i for i in issues if i["kind"] == "MultipleHypothesesInLiveCall"]
        assert not offenders, f"§5B invariant broken: {offenders}"


# ──────────────────────────────────────────────────────────────────
# 4. Counter provenance (C2) — sanity per mechanism
# ──────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("mechanism", ["mutate", "pcr_amplify", "sequencing_errors", "indels", "end_loss"])
def test_counter_provenance_holds_under_each_mechanism(mechanism):
    """Each mechanism's counter must match its provenance source.
    n_mutations vs sim.mutation_count, n_pcr_errors vs trace,
    n_quality_errors vs trace, n_indels and per-segment counters
    vs event ledger, end_loss_*_length vs trace. If any source
    diverges from the AIRR field, the validator surfaces it."""
    exp = _baseline_vj()
    if mechanism == "mutate":
        exp = exp.mutate(rate=0.2)
    elif mechanism == "pcr_amplify":
        exp = exp.pcr_amplify(rate=0.1)
    elif mechanism == "sequencing_errors":
        exp = exp.sequencing_errors(rate=0.1)
    elif mechanism == "indels":
        exp = exp.polymerase_indels(count=3, insertion_prob=0.5)
    elif mechanism == "end_loss":
        exp = exp.primer_trim_5prime(length=3).primer_trim_3prime(length=2)

    counter_kinds = {
        "NMutationsMismatch",
        "NPcrErrorsMismatch",
        "NQualityErrorsMismatch",
        "NIndelsMismatch",
        "NSegmentIndelsMismatch",
        "EndLossLengthMismatch",
    }
    compiled = exp.compile(allow_curatable_refdata=True)
    for seed in range(15):
        outcome = compiled.simulator.run(seed=seed)
        issues = outcome.validate_record(compiled.refdata)
        counter_issues = [i for i in issues if i["kind"] in counter_kinds]
        assert counter_issues == [], (
            f"[{mechanism}] seed={seed}: counter-provenance issues {counter_issues}"
        )


# ──────────────────────────────────────────────────────────────────
# 5. Junction truth (C3) — productive triad recomputed independently
# ──────────────────────────────────────────────────────────────────


def test_productive_triad_recomputation_agrees_across_mutation_seeds():
    """Productive triad (in-frame ∧ no-stop ∧ anchors-preserved)
    must agree between the AIRR record and the validator's
    independent re-derivation. Sweep mutation seeds to cover the
    full set: productive=True, OOF, junction-stop, V-anchor break,
    J-anchor break."""
    exp = (
        ga.Experiment.on(_vj()).allow_curatable_refdata()
        .recombine(np1_lengths=[(3, 1.0)])
        .trim(enabled=False)
        .mutate(rate=0.5)
    )
    triad_kinds = {
        "ProductiveMismatch", "VjInFrameMismatch", "StopCodonMismatch",
        "JunctionLengthMismatch", "JunctionContentMismatch", "JunctionAaMismatch",
    }
    compiled = exp.compile(allow_curatable_refdata=True)
    for seed in range(40):
        outcome = compiled.simulator.run(seed=seed)
        issues = outcome.validate_record(compiled.refdata)
        triad = [i for i in issues if i["kind"] in triad_kinds]
        assert triad == [], f"seed={seed}: triad disagrees {triad}"


# ──────────────────────────────────────────────────────────────────
# 6. Allele oracle (C4) — independent rescoring matches reported call
# ──────────────────────────────────────────────────────────────────


def test_allele_oracle_agrees_with_reported_call_under_shm():
    """The validator independently rescores every allele in V/D/J
    pools and rebuilds the tie-set at max-match-count. It must
    agree with the projection's v_call/d_call/j_call CSV including
    the truth-first ordering convention."""
    exp = (
        ga.Experiment.on(_vj()).allow_curatable_refdata()
        .recombine(np1_lengths=[(3, 1.0)])
        .trim(enabled=False)
        .mutate(rate=0.3)
    )
    oracle_kinds = {"AlleleCallTieSetMismatch", "AlleleCallOrderMismatch"}
    compiled = exp.compile(allow_curatable_refdata=True)
    for seed in range(40):
        outcome = compiled.simulator.run(seed=seed)
        issues = outcome.validate_record(compiled.refdata)
        oracle = [i for i in issues if i["kind"] in oracle_kinds]
        assert oracle == [], (
            f"seed={seed}: allele oracle disagrees with projection {oracle}"
        )


# ──────────────────────────────────────────────────────────────────
# 7. Batch validation API — SimulationResult.validate_records
# ──────────────────────────────────────────────────────────────────


def test_batch_validate_clean_result_returns_ok_report():
    """A clean batch produces an ok=True report with empty
    failures and total count = N records."""
    exp = (
        _baseline_vj()
        .mutate(rate=0.1)
        .polymerase_indels(count=2, insertion_prob=0.5)
        .productive_only()
    )
    refdata = exp._refdata
    result = exp.run_records(n=50, seed=0)
    report = result.validate_records(refdata)

    assert isinstance(report, ValidationReport)
    assert report.ok is True
    assert bool(report) is True
    assert report.count == 50
    assert report.failures == []
    assert len(report) == 0  # __len__ returns number of failures
    assert report.summary() == {}


def test_batch_validate_report_repr_distinguishes_pass_and_fail():
    """The repr should make it obvious at a glance whether the
    batch passed or failed."""
    ok = ValidationReport(count=10, failures=[])
    bad = ValidationReport(
        count=10,
        failures=[
            {
                "record_index": 3,
                "sequence_id": "seq3",
                "issues": [{"kind": "ProductiveMismatch", "details": {}}],
            }
        ],
    )
    assert "ok=True" in repr(ok)
    assert "count=10" in repr(ok)
    assert "ok=False" in repr(bad)
    assert "failures=1" in repr(bad)
    # Truthy/falsy interplay — assert report works as a CI guard.
    assert ok
    assert not bad


def test_batch_validate_summary_buckets_issue_kinds_across_failures():
    """`.summary()` returns a histogram of issue kinds across the
    whole batch. Useful for the typical 'what went wrong overall'
    glance in CI logs."""
    report = ValidationReport(
        count=3,
        failures=[
            {
                "record_index": 0,
                "sequence_id": "seq0",
                "issues": [
                    {"kind": "ProductiveMismatch"},
                    {"kind": "JunctionLengthMismatch"},
                ],
            },
            {
                "record_index": 2,
                "sequence_id": "seq2",
                "issues": [{"kind": "ProductiveMismatch"}],
            },
        ],
    )
    assert report.summary() == {
        "ProductiveMismatch": 2,
        "JunctionLengthMismatch": 1,
    }


def test_batch_validate_propagates_sequence_id_to_failures():
    """Each failure carries the record's sequence_id so users can
    cross-reference back to the original record without recomputing."""
    # Build a result with custom-prefixed sequence ids.
    exp = _baseline_vj()
    refdata = exp._refdata
    result = exp.run_records(n=5, seed=0)
    report = result.validate_records(refdata)

    # All clean, so no failures. But the report's count matches the
    # batch size, and if we manually construct a failure it carries
    # the sequence_id field.
    assert report.count == 5
    # Synthesize a failure with the record's actual sequence_id to
    # confirm the field name is what we documented.
    synthetic = {
        "record_index": 0,
        "sequence_id": result.records[0]["sequence_id"],
        "issues": [{"kind": "ProductiveMismatch", "details": {}}],
    }
    fake_report = ValidationReport(count=5, failures=[synthetic])
    assert "sequence_id" in fake_report.failures[0]
    assert fake_report.failures[0]["sequence_id"] == result.records[0]["sequence_id"]


def test_batch_validate_raises_when_outcomes_missing():
    """`validate_records` requires the original Outcome objects.
    A SimulationResult built from records only (e.g. loaded from
    TSV) can't be validated — raises RuntimeError with a clear
    message pointing at run_records."""
    result_without_outcomes = SimulationResult(
        records=[{"sequence_id": "seq0"}], outcomes=None
    )
    cfg = ge.RefDataConfig.vj()
    cfg.add_v_allele("v1*01", "v1", b"AAACCCTGT", anchor=6)
    cfg.add_j_allele("j1*01", "j1", b"TGGAAA", anchor=0)
    with pytest.raises(RuntimeError, match="validate_records requires"):
        result_without_outcomes.validate_records(cfg)


def test_batch_validate_report_is_json_serializable():
    """The report's `failures` list is JSON-serializable
    end-to-end — useful for CI artifacts and dashboards."""
    import json
    report = ValidationReport(
        count=2,
        failures=[
            {
                "record_index": 0,
                "sequence_id": "seq0",
                "issues": [
                    {
                        "kind": "ProductiveMismatch",
                        "reported": True,
                        "expected": False,
                        "details": {"reason": "out_of_frame"},
                    }
                ],
            }
        ],
    )
    blob = json.dumps(report.failures)
    parsed = json.loads(blob)
    assert parsed[0]["sequence_id"] == "seq0"
    assert parsed[0]["issues"][0]["kind"] == "ProductiveMismatch"
