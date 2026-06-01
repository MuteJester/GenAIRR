"""Trace-injected replay (Option B), first slice.

Covers:
- Existing seed-based rerun (`rerun_from_trace_file`) still works.
- Trace-injected replay (`replay_from_trace_file`) reproduces the
  outcome and the migrated-pass consumes from the cursor.
- Error surface: address mismatch and value-kind mismatch produce
  readable `StrictSamplingError` with stable `replay.<kind>` tags.
- Exhausted-cursor (truncated input trace) is reported as a
  structured replay error.

Only `SampleAllelePass` is wired in this slice; the remaining
sampling sites still draw from the RNG. Tests use a `recombine()`
plan so the migrated-pass path is exercised.
"""
from __future__ import annotations

import json

import pytest

import GenAIRR as ga
from GenAIRR import _engine as eng


# ──────────────────────────────────────────────────────────────────
# Fixture: a fast VJ experiment.
# ──────────────────────────────────────────────────────────────────


@pytest.fixture(scope="module")
def vj_compiled():
    return ga.Experiment.on("human_igk").recombine().compile()


# ──────────────────────────────────────────────────────────────────
# rerun_from_trace_file (Option A) — sanity check still works.
# ──────────────────────────────────────────────────────────────────


def test_rerun_from_trace_file_still_works(vj_compiled):
    outcome = vj_compiled.simulator.run(seed=11)
    tf = vj_compiled.simulator.trace_file_from(outcome, seed=11)
    rerun = vj_compiled.simulator.rerun_from_trace_file(tf)
    # Both paths must give the same trace.
    for a, b in zip(outcome.trace().choices(), rerun.trace().choices()):
        assert a.address == b.address
        assert a.value == b.value


# ──────────────────────────────────────────────────────────────────
# replay_from_trace_file (Option B) — trace-injected reproduction.
# ──────────────────────────────────────────────────────────────────


# All sampling sites consume from the cursor in replay mode. Full
# migration was reached when Tier 3 (NP bases, contaminant, end-loss,
# indel tuples) landed; the engine is now deterministic in trace
# rather than in seed.
#
# The set below is retained for forward-compat / future migration
# work — if a new pass is added that records a trace, its prefix
# goes here. The strict-positional cursor requires every record in
# the input trace to be consumed by some pass during replay; that
# invariant is enforced via `assert_drained` in
# `execute_transactional`.
MIGRATED_ADDRESS_PREFIXES = (
    "sample_allele.",
    "trim.",
    "np.",
    "corrupt.",
    "mutate.",
)


def _all_records_have_migrated_prefix(outcome) -> bool:
    """Sanity check used by the full-stack golden test: every trace
    record must fall under one of the migrated prefixes. If a future
    pass emits a non-migrated record, this assertion catches the
    drift before the strict-positional cursor breaks at replay."""
    return all(
        rec.address.startswith(MIGRATED_ADDRESS_PREFIXES)
        for rec in outcome.trace().choices()
    )


def test_replay_from_trace_file_reproduces_full_outcome(vj_compiled):
    # Now that every sampling site is migrated, replay reproduces
    # the ORIGINAL trace byte-for-byte — not just a migrated subset.
    # The engine is deterministic in trace rather than in seed: the
    # cursor injects every value the sampler would have drawn, and
    # the RNG is bypassed end-to-end.
    outcome_a = vj_compiled.simulator.run(seed=4242)
    tf = vj_compiled.simulator.trace_file_from(outcome_a, seed=4242)
    outcome_b = vj_compiled.simulator.replay_from_trace_file(tf)

    assert outcome_a.pass_names() == outcome_b.pass_names()
    a_choices = outcome_a.trace().choices()
    b_choices = outcome_b.trace().choices()
    assert len(a_choices) == len(b_choices)
    for a, b in zip(a_choices, b_choices):
        assert a.address == b.address
        assert a.value == b.value


def test_replay_from_trace_file_injects_allele_value_from_cursor(vj_compiled):
    outcome_a = vj_compiled.simulator.run(seed=99)
    tf = vj_compiled.simulator.trace_file_from(outcome_a, seed=99)

    # Read the original V allele id from the trace.
    original_v = outcome_a.trace().find("sample_allele.v").value

    # Replay → migrated SampleAllelePass consumes that exact id from
    # the cursor (not from the RNG).
    outcome_b = vj_compiled.simulator.replay_from_trace_file(tf)
    replayed_v = outcome_b.trace().find("sample_allele.v").value
    assert original_v == replayed_v


def test_replay_against_disk_round_tripped_trace_file(tmp_path, vj_compiled):
    outcome_a = vj_compiled.simulator.run(seed=7)
    tf = vj_compiled.simulator.trace_file_from(outcome_a, seed=7)
    path = tmp_path / "trace.json"
    tf.write_to(str(path))

    loaded = eng.TraceFile.read_from(str(path))
    outcome_b = vj_compiled.simulator.replay_from_trace_file(loaded)
    a_choices = outcome_a.trace().choices()
    b_choices = outcome_b.trace().choices()
    assert len(a_choices) == len(b_choices)
    for a, b in zip(a_choices, b_choices):
        assert a.address == b.address
        assert a.value == b.value


# ──────────────────────────────────────────────────────────────────
# Signature checks (same shape as rerun_from_trace_file).
# ──────────────────────────────────────────────────────────────────


def test_replay_rejects_plan_mismatch(vj_compiled):
    outcome = vj_compiled.simulator.run(seed=0)
    tf = vj_compiled.simulator.trace_file_from(outcome, seed=0)
    different = (
        ga.Experiment.on("human_igk").recombine().pcr_amplify(rate=0.001).compile()
    )
    with pytest.raises(ValueError, match="pass plan signature mismatch"):
        different.simulator.replay_from_trace_file(tf)


def test_replay_rejects_refdata_mismatch(vj_compiled):
    outcome = vj_compiled.simulator.run(seed=0)
    tf = vj_compiled.simulator.trace_file_from(outcome, seed=0)
    # IGL vs IGK changes both the refdata identity and the
    # cartridge-derived trim / NP distributions baked into the
    # plan parameters. Under Slice A's parameterized plan
    # signature, the plan-parameter gate fires before the refdata
    # gate; accept either rejection — the point of the test is
    # that the replay is refused.
    different_refdata = ga.Experiment.on("human_igl").recombine().compile()
    with pytest.raises(
        ValueError,
        match="(refdata signature mismatch|pass plan signature mismatch)",
    ):
        different_refdata.simulator.replay_from_trace_file(tf)


# ──────────────────────────────────────────────────────────────────
# Replay-error surface: address mismatch, value-kind mismatch,
# exhausted cursor. Each produces a structured StrictSamplingError
# with a stable `replay.<kind>` tag in the reason payload.
# ──────────────────────────────────────────────────────────────────


def _mutate_trace_file_json(tf: "eng.TraceFile", mutate) -> "eng.TraceFile":
    """Round-trip a TraceFile through JSON, mutating the dict in
    between, so tests can inject malformed records without going
    through the engine's recording path."""
    data = json.loads(tf.to_json())
    mutate(data)
    return eng.TraceFile.from_json(json.dumps(data))


def test_replay_address_mismatch_surfaces_structured_error(vj_compiled):
    outcome = vj_compiled.simulator.run(seed=12)
    tf = vj_compiled.simulator.trace_file_from(outcome, seed=12)

    # Corrupt the first record's address: swap V → J. The migrated
    # SampleAllelePass.V will see address "sample_allele.j" instead
    # of "sample_allele.v" and raise an AddressMismatch.
    def swap_first_address(d):
        d["trace"]["choices"][0]["address"] = "sample_allele.j"

    corrupted = _mutate_trace_file_json(tf, swap_first_address)

    with pytest.raises(eng.StrictSamplingError) as exc_info:
        vj_compiled.simulator.replay_from_trace_file(corrupted)
    # StrictSamplingError args = (pass_name, address, reason).
    pass_name, _, reason = exc_info.value.args
    assert pass_name == "sample_allele.v"
    assert "replay.address_mismatch" in reason
    assert "sample_allele.v" in reason  # expected
    assert "sample_allele.j" in reason  # got


def test_replay_value_kind_mismatch_surfaces_structured_error(vj_compiled):
    outcome = vj_compiled.simulator.run(seed=21)
    tf = vj_compiled.simulator.trace_file_from(outcome, seed=21)

    # Corrupt the first record's value: replace AlleleId with Int.
    def swap_first_value_kind(d):
        d["trace"]["choices"][0]["value"] = {"Int": 0}

    corrupted = _mutate_trace_file_json(tf, swap_first_value_kind)

    with pytest.raises(eng.StrictSamplingError) as exc_info:
        vj_compiled.simulator.replay_from_trace_file(corrupted)
    pass_name, _, reason = exc_info.value.args
    assert pass_name == "sample_allele.v"
    assert "replay.value_kind_mismatch" in reason
    assert "AlleleId" in reason
    assert "Int" in reason


def test_replay_exhausted_cursor_surfaces_structured_error(vj_compiled):
    outcome = vj_compiled.simulator.run(seed=33)
    tf = vj_compiled.simulator.trace_file_from(outcome, seed=33)

    # Truncate the trace to zero records: the first migrated pass
    # asks for sample_allele.v, the cursor is empty, ExhaustedTrace.
    def truncate_all(d):
        d["trace"]["choices"] = []

    corrupted = _mutate_trace_file_json(tf, truncate_all)

    with pytest.raises(eng.StrictSamplingError) as exc_info:
        vj_compiled.simulator.replay_from_trace_file(corrupted)
    pass_name, _, reason = exc_info.value.args
    assert pass_name == "sample_allele.v"
    assert "replay.exhausted" in reason
    assert "sample_allele.v" in reason


# ──────────────────────────────────────────────────────────────────
# Tier-1 per-pass replay coverage
#
# Each migrated pass uses a distinct typed `expect_*` helper:
#   - SampleAllelePass → expect_allele_id (covered above)
#   - TrimPass         → expect_int
#   - GenerateNPPass   → expect_int (length only — bases still on RNG)
#   - RevCompPass      → expect_bool
#
# Tests below pin that:
#   (a) the migrated address consumes the cursor's value, and
#   (b) the downstream application/validation path runs unchanged,
#       so the resulting IR matches what fresh sampling would have
#       produced at the same recorded values.
# ──────────────────────────────────────────────────────────────────


@pytest.fixture(scope="module")
def vj_compiled_with_rev_comp():
    """VJ plan that also exercises RevCompPass so the Bool consumer
    path is reachable. `prob=1.0` so the flag deterministically lands
    as True under any seed."""
    return (
        ga.Experiment.on("human_igk")
        .recombine()
        .random_strand_orientation(prob=1.0)
        .compile()
    )


def test_replay_consumes_trim_length_from_cursor(vj_compiled):
    # The IGK recombine plan trims V_3 + J_5. Both must consume from
    # the cursor under replay.
    outcome_a = vj_compiled.simulator.run(seed=55)
    tf = vj_compiled.simulator.trace_file_from(outcome_a, seed=55)
    outcome_b = vj_compiled.simulator.replay_from_trace_file(tf)

    for addr in ("trim.v_3", "trim.j_5"):
        original = outcome_a.trace().find(addr)
        replayed = outcome_b.trace().find(addr)
        assert original is not None, f"plan should record {addr}"
        assert replayed is not None
        assert original.value == replayed.value, (
            f"trim {addr} must consume from cursor, "
            f"got original={original.value} replayed={replayed.value}"
        )


def test_replay_consumes_np_length_from_cursor(vj_compiled):
    outcome_a = vj_compiled.simulator.run(seed=88)
    tf = vj_compiled.simulator.trace_file_from(outcome_a, seed=88)
    outcome_b = vj_compiled.simulator.replay_from_trace_file(tf)

    original = outcome_a.trace().find("np.np1.length")
    replayed = outcome_b.trace().find("np.np1.length")
    assert original is not None
    assert replayed is not None
    assert original.value == replayed.value


# Note: previously gated by the Tier-3 NP-base migration; unblocked
# now that every sampling site (allele/trim/NP-length/NP-base/SHM/
# corruption/end-loss/indel/contaminant) consumes from the cursor.
def test_replay_consumes_rev_comp_bool_from_cursor(vj_compiled_with_rev_comp):
    outcome_a = vj_compiled_with_rev_comp.simulator.run(seed=3)
    tf = vj_compiled_with_rev_comp.simulator.trace_file_from(outcome_a, seed=3)
    outcome_b = vj_compiled_with_rev_comp.simulator.replay_from_trace_file(tf)

    original = outcome_a.trace().find("corrupt.rev_comp.applied")
    replayed = outcome_b.trace().find("corrupt.rev_comp.applied")
    assert original is not None
    assert replayed is not None
    assert original.value == replayed.value


# Wrong-kind / wrong-address errors at each migrated address. The
# error chain (StrictSamplingError carrying a `replay.<kind>` tag) is
# the same across passes — these tests confirm each pass's site
# actually triggers the chain.


def test_replay_trim_value_kind_mismatch(vj_compiled):
    outcome = vj_compiled.simulator.run(seed=0)
    tf = vj_compiled.simulator.trace_file_from(outcome, seed=0)

    # Corrupt the first trim record (trim.v_3 sits after the two
    # sample_allele records → index 2).
    def swap_trim_kind(d):
        choices = d["trace"]["choices"]
        # Locate the first "trim." record.
        idx = next(i for i, c in enumerate(choices) if c["address"].startswith("trim."))
        # Replace Int payload with Base — wrong kind.
        choices[idx]["value"] = {"Base": "A"}

    corrupted = _mutate_trace_file_json(tf, swap_trim_kind)

    with pytest.raises(eng.StrictSamplingError) as exc_info:
        vj_compiled.simulator.replay_from_trace_file(corrupted)
    _, _, reason = exc_info.value.args
    assert "replay.value_kind_mismatch" in reason
    assert "Int" in reason
    assert "Base" in reason


def test_replay_np_length_value_kind_mismatch(vj_compiled):
    outcome = vj_compiled.simulator.run(seed=0)
    tf = vj_compiled.simulator.trace_file_from(outcome, seed=0)

    def swap_np_length_kind(d):
        choices = d["trace"]["choices"]
        idx = next(i for i, c in enumerate(choices) if c["address"] == "np.np1.length")
        choices[idx]["value"] = {"Bool": True}

    corrupted = _mutate_trace_file_json(tf, swap_np_length_kind)

    with pytest.raises(eng.StrictSamplingError) as exc_info:
        vj_compiled.simulator.replay_from_trace_file(corrupted)
    _, _, reason = exc_info.value.args
    assert "replay.value_kind_mismatch" in reason


def test_replay_rev_comp_value_kind_mismatch(vj_compiled_with_rev_comp):
    outcome = vj_compiled_with_rev_comp.simulator.run(seed=0)
    tf = vj_compiled_with_rev_comp.simulator.trace_file_from(outcome, seed=0)

    def swap_rev_comp_kind(d):
        choices = d["trace"]["choices"]
        idx = next(
            i for i, c in enumerate(choices) if c["address"] == "corrupt.rev_comp.applied"
        )
        choices[idx]["value"] = {"Int": 0}

    corrupted = _mutate_trace_file_json(tf, swap_rev_comp_kind)

    with pytest.raises(eng.StrictSamplingError) as exc_info:
        vj_compiled_with_rev_comp.simulator.replay_from_trace_file(corrupted)
    _, _, reason = exc_info.value.args
    assert "replay.value_kind_mismatch" in reason


def test_replay_trim_runs_validation_path_on_bad_recorded_value(vj_compiled):
    # The invariant: replay consumes the value, then runs the same
    # validation/application path as fresh sampling. A trim length
    # out of range (negative or > u16::MAX) must still surface as
    # `InvalidDistributionOutput` — not silently corrupt the IR.
    outcome = vj_compiled.simulator.run(seed=0)
    tf = vj_compiled.simulator.trace_file_from(outcome, seed=0)

    def bad_trim(d):
        choices = d["trace"]["choices"]
        idx = next(i for i, c in enumerate(choices) if c["address"].startswith("trim."))
        choices[idx]["value"] = {"Int": -7}

    corrupted = _mutate_trace_file_json(tf, bad_trim)

    with pytest.raises(eng.StrictSamplingError) as exc_info:
        vj_compiled.simulator.replay_from_trace_file(corrupted)
    pass_name, _, reason = exc_info.value.args
    # The TrimPass validation path produced the error, not the
    # replay layer — proving the migration runs through the shared
    # validation/application path rather than a special-case branch.
    assert pass_name.startswith("trim.")
    assert "negative_trim" in reason


# ──────────────────────────────────────────────────────────────────
# Golden full-stack replay-from-trace-file
#
# Builds a plan exercising EVERY migrated sampling site — allele,
# trim, NP length+bases, SHM count+site+base, PCR count+site+base,
# rev-comp coin, end-loss, contaminant, indel tuple — then commits
# the produced trace as a golden artifact. Replay-from-trace-file
# against the golden must reproduce the trace byte-for-byte under
# the strict-positional cursor with the drained-cursor assertion
# enabled.
#
# A failure here means either:
#   1. A new sampling site was added but its replay path wasn't
#      migrated, OR
#   2. An existing pass's record shape changed in a way that
#      diverges from the committed golden.
#
# Either case warrants explicit investigation before regenerating
# the golden.
# ──────────────────────────────────────────────────────────────────


from pathlib import Path  # noqa: E402

FULL_STACK_GOLDEN = (
    Path(__file__).parent / "golden" / "trace_files" / "full_stack_seed4242.trace.json"
)


@pytest.fixture(scope="module")
def full_stack_compiled():
    """A plan that exercises every migrated sampling site."""
    return (
        ga.Experiment.on("human_igk")
        .recombine()
        .mutate(model="uniform", count=2)
        .pcr_amplify(count=1)
        .sequencing_errors(count=1)
        .polymerase_indels(count=2, insertion_prob=0.5)
        .ambiguous_base_calls(count=1)
        .primer_trim_5prime(length=2)
        .contaminate(prob=0.0)  # disabled to keep the plan compact
        .random_strand_orientation(prob=0.5)
        .compile()
    )


def test_full_stack_replay_reproduces_committed_golden_trace(full_stack_compiled):
    """Emit-once / compare-thereafter golden test for the full
    trace-injected replay surface.

    On first run with no committed golden, the test emits the file
    and skips with a note. On subsequent runs, it verifies the live
    engine reproduces the trace and that replaying the committed
    golden through `replay_from_trace_file` also reproduces it.
    """
    seed = 4242
    outcome = full_stack_compiled.simulator.run(seed=seed)
    tf = full_stack_compiled.simulator.trace_file_from(outcome, seed=seed)

    # Every recorded address must fall under a migrated prefix —
    # otherwise the strict-positional drained-cursor assertion
    # would reject the trace at replay time.
    assert _all_records_have_migrated_prefix(outcome), (
        "Some sampling site emitted a non-migrated trace address. "
        "Update MIGRATED_ADDRESS_PREFIXES or migrate the site to "
        "the consume-trace path before regenerating the golden."
    )

    # Read or emit the golden artifact.
    FULL_STACK_GOLDEN.parent.mkdir(parents=True, exist_ok=True)
    if not FULL_STACK_GOLDEN.exists():
        tf.write_to(str(FULL_STACK_GOLDEN))
        pytest.skip(f"Emitted fresh full-stack golden at {FULL_STACK_GOLDEN}")

    loaded = eng.TraceFile.read_from(str(FULL_STACK_GOLDEN))
    # The fixture is preserved as a v1 backwards-compat anchor
    # (see test_trace_file_compat.py). When its schema is older
    # than the live one, the plan-signature field cannot be equal
    # by construction — the schema-version-aware comparator
    # inside `replay_from_trace_file` handles cross-schema replay
    # transparently. The determinism check is the per-choice
    # address+value equality below, which is what this test is
    # actually pinning.
    assert loaded.refdata_signature == tf.refdata_signature
    if loaded.schema_version == tf.schema_version:
        assert loaded.pass_plan_signature == tf.pass_plan_signature

    # Replay via the trace-injected path — RNG bypassed, every
    # value comes from the cursor, drained-cursor assertion runs
    # at the end of the plan.
    outcome_b = full_stack_compiled.simulator.replay_from_trace_file(loaded)

    # Byte-for-byte trace equality between the live run and the
    # cursor-driven replay.
    a_choices = outcome.trace().choices()
    b_choices = outcome_b.trace().choices()
    assert len(a_choices) == len(b_choices), (
        f"trace length mismatch: live={len(a_choices)} replayed={len(b_choices)}"
    )
    for i, (a, b) in enumerate(zip(a_choices, b_choices)):
        assert a.address == b.address, f"address mismatch at index {i}: {a.address} vs {b.address}"
        assert a.value == b.value, f"value mismatch at index {i} (addr={a.address})"
