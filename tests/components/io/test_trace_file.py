"""Durable trace-file artifact: round-trip + rerun reproducibility.

Covers the Python surface of the Option-A slice:

- `Outcome.trace()` records the run's choices in chronological order
  (already covered by the engine smoke tests; here we use it as a
  building block for the trace-file flow).
- `CompiledSimulator.trace_file_from(outcome, seed)` bundles a trace
  plus the simulator's plan / refdata signatures into a `TraceFile`.
- `TraceFile.write_to(path)` + `TraceFile.read_from(path)` round-trip
  the bundle through disk-resident JSON.
- `CompiledSimulator.rerun_from_trace_file(tf)` reproduces the
  original outcome byte-for-byte (engine determinism).
- Signature mismatches (plan or refdata) raise `ValueError` at rerun
  time rather than silently producing a different outcome.

The fixture trace files live in `tests/golden/trace_files/` so a
future Option-B (trace-injected replay) can consume them as
reference inputs.
"""
from __future__ import annotations

import json

import pytest

import GenAIRR as ga
from GenAIRR import _engine as genairr_engine

pytestmark = pytest.mark.behavior


# ──────────────────────────────────────────────────────────────────
# Fixtures: a tiny VJ experiment that's cheap to compile + run.
# Using an in-tree DataConfig keeps the test runtime under a second.
# ──────────────────────────────────────────────────────────────────


@pytest.fixture(scope="module")
def vj_experiment() -> ga.Experiment:
    return ga.Experiment.on("human_igk").recombine()


@pytest.fixture(scope="module")
def vj_compiled(vj_experiment) -> "ga.CompiledExperiment":
    return vj_experiment.compile()


# ──────────────────────────────────────────────────────────────────
# Module surface
# ──────────────────────────────────────────────────────────────────


def test_trace_file_class_is_exported():
    assert hasattr(genairr_engine, "TraceFile")


def test_trace_file_schema_version_class_attr():
    # The schema-version constant should be a positive int that
    # matches the on-disk shape we produce. Bumping it is an
    # explicit incompatibility marker.
    assert isinstance(genairr_engine.TraceFile.SCHEMA_VERSION, int)
    assert genairr_engine.TraceFile.SCHEMA_VERSION >= 1


# ──────────────────────────────────────────────────────────────────
# Bundle: CompiledSimulator.trace_file_from(outcome, seed)
# ──────────────────────────────────────────────────────────────────


def test_trace_file_from_outcome_carries_metadata(vj_compiled):
    outcome = vj_compiled.simulator.run(seed=42)
    tf = vj_compiled.simulator.trace_file_from(outcome, seed=42)

    assert tf.schema_version == genairr_engine.TraceFile.SCHEMA_VERSION
    assert tf.engine_version  # non-empty string
    assert tf.seed == 42
    # Plan signature is a pipe-joined pass-name string.
    assert "|" in tf.pass_plan_signature
    # Refdata signature is the structural form `chain=...;v=...;d=...;j=...;c=...`.
    assert tf.refdata_signature.startswith("chain=")
    # Trace records survive the bundle.
    assert len(tf.trace) == len(outcome.trace())


def test_trace_file_from_outcome_trace_matches_source(vj_compiled):
    outcome = vj_compiled.simulator.run(seed=7)
    tf = vj_compiled.simulator.trace_file_from(outcome, seed=7)

    source_choices = outcome.trace().choices()
    bundled_choices = tf.trace.choices()
    assert len(source_choices) == len(bundled_choices)
    for a, b in zip(source_choices, bundled_choices):
        assert a.address == b.address
        assert a.value == b.value


# ──────────────────────────────────────────────────────────────────
# Round-trip: JSON string + disk
# ──────────────────────────────────────────────────────────────────


def test_trace_file_round_trips_through_json_string(vj_compiled):
    outcome = vj_compiled.simulator.run(seed=11)
    tf = vj_compiled.simulator.trace_file_from(outcome, seed=11)

    s = tf.to_json()
    back = genairr_engine.TraceFile.from_json(s)

    assert back.schema_version == tf.schema_version
    assert back.engine_version == tf.engine_version
    assert back.seed == tf.seed
    assert back.pass_plan_signature == tf.pass_plan_signature
    assert back.refdata_signature == tf.refdata_signature
    assert len(back.trace) == len(tf.trace)


def test_trace_file_round_trips_through_disk(tmp_path, vj_compiled):
    outcome = vj_compiled.simulator.run(seed=99)
    tf = vj_compiled.simulator.trace_file_from(outcome, seed=99)

    path = tmp_path / "trace.json"
    tf.write_to(str(path))
    back = genairr_engine.TraceFile.read_from(str(path))

    assert back.seed == tf.seed
    assert back.pass_plan_signature == tf.pass_plan_signature
    assert back.refdata_signature == tf.refdata_signature
    assert len(back.trace) == len(tf.trace)


def test_trace_file_json_is_pretty_printed_and_loads_as_dict(tmp_path, vj_compiled):
    outcome = vj_compiled.simulator.run(seed=0)
    tf = vj_compiled.simulator.trace_file_from(outcome, seed=0)
    path = tmp_path / "trace.json"
    tf.write_to(str(path))

    raw = path.read_text()
    # Pretty-printing produces multi-line JSON.
    assert "\n" in raw
    parsed = json.loads(raw)
    assert parsed["schema_version"] == genairr_engine.TraceFile.SCHEMA_VERSION
    assert parsed["seed"] == 0
    assert "trace" in parsed
    assert "choices" in parsed["trace"]


# ──────────────────────────────────────────────────────────────────
# Rerun: signature-checked reproduction
# ──────────────────────────────────────────────────────────────────


def test_rerun_from_trace_file_reproduces_outcome(vj_compiled):
    outcome_a = vj_compiled.simulator.run(seed=1234)
    tf = vj_compiled.simulator.trace_file_from(outcome_a, seed=1234)

    outcome_b = vj_compiled.simulator.rerun_from_trace_file(tf)

    # Same pass count.
    assert outcome_b.pass_names() == outcome_a.pass_names()
    # Same trace, address-for-address.
    a_choices = outcome_a.trace().choices()
    b_choices = outcome_b.trace().choices()
    assert len(a_choices) == len(b_choices)
    for a, b in zip(a_choices, b_choices):
        assert a.address == b.address
        assert a.value == b.value


def test_rerun_from_trace_file_after_disk_round_trip_is_deterministic(
    tmp_path, vj_compiled,
):
    outcome_a = vj_compiled.simulator.run(seed=5678)
    tf = vj_compiled.simulator.trace_file_from(outcome_a, seed=5678)
    path = tmp_path / "trace.json"
    tf.write_to(str(path))

    loaded = genairr_engine.TraceFile.read_from(str(path))
    outcome_b = vj_compiled.simulator.rerun_from_trace_file(loaded)

    a_choices = outcome_a.trace().choices()
    b_choices = outcome_b.trace().choices()
    assert len(a_choices) == len(b_choices)
    for a, b in zip(a_choices, b_choices):
        assert a.address == b.address
        assert a.value == b.value


# ──────────────────────────────────────────────────────────────────
# Signature mismatches must raise
# ──────────────────────────────────────────────────────────────────


def test_rerun_from_trace_file_rejects_plan_mismatch(vj_compiled):
    outcome = vj_compiled.simulator.run(seed=0)
    tf = vj_compiled.simulator.trace_file_from(outcome, seed=0)

    # Compile a *different* plan against the same refdata: add a
    # pcr_amplify step so the pass-plan signature changes.
    different = (
        ga.Experiment.on("human_igk").recombine().pcr_amplify(rate=0.001).compile()
    )

    with pytest.raises(ValueError, match="pass plan signature mismatch"):
        different.simulator.rerun_from_trace_file(tf)


def test_rerun_from_trace_file_rejects_refdata_mismatch(vj_compiled):
    outcome = vj_compiled.simulator.run(seed=0)
    tf = vj_compiled.simulator.trace_file_from(outcome, seed=0)

    # Compile the same plan shape against a *different* refdata
    # (IGL vs IGK). The two compiled plans differ in TWO ways:
    # (a) cartridge identity (refdata content hash), and (b) the
    # cartridge-derived empirical distributions baked into the
    # trim / NP-length pass parameters. Either one is enough for
    # `rerun_from_trace_file` to reject; under Slice A's
    # parameterized plan signature, the plan-parameter mismatch
    # now fires before the refdata gate. The point of this test
    # is that the rerun is rejected — accept any of the relevant
    # mismatch messages.
    different_refdata = ga.Experiment.on("human_igl").recombine().compile()

    with pytest.raises(
        ValueError,
        match="(refdata signature mismatch|pass plan signature mismatch)",
    ):
        different_refdata.simulator.rerun_from_trace_file(tf)


# ──────────────────────────────────────────────────────────────────
# Schema-version guard
# ──────────────────────────────────────────────────────────────────


def test_trace_file_from_json_rejects_future_schema_version():
    fake = {
        "schema_version": 9999,
        "engine_version": "0.0.0",
        "seed": 0,
        "pass_plan_signature": "",
        "refdata_signature": "",
        "trace": {"choices": []},
    }
    with pytest.raises(ValueError, match="schema version"):
        genairr_engine.TraceFile.from_json(json.dumps(fake))


# ──────────────────────────────────────────────────────────────────
# Address validation
# ──────────────────────────────────────────────────────────────────


def test_validate_addresses_accepts_built_in_run(vj_compiled):
    outcome = vj_compiled.simulator.run(seed=0)
    tf = vj_compiled.simulator.trace_file_from(outcome, seed=0)
    # Every recorded address must parse to a built-in ChoiceAddress.
    tf.validate_addresses()


# ──────────────────────────────────────────────────────────────────
# The emitted artifact, round-tripped
# ──────────────────────────────────────────────────────────────────


def test_an_emitted_trace_file_reruns_to_the_same_choices(tmp_path, vj_compiled):
    """A trace file written to disk reruns to the choices it recorded.

    This used to compare a fresh run against
    `tests/golden/trace_files/human_igk_recombine_seed4242.trace.json`,
    emitting it if absent, and say that a failure meant either a schema bump
    or a determinism regression. There is a third case, and it is the one that
    arrived: a deliberate determinism change. Per-pass RNG substreams move
    every recorded value, so seed 4242 no longer draws what that file holds.

    Regenerating the file was the wrong answer. It is the canonical v1
    backwards-compat anchor -- `test_trace_file_compat.py::
    test_v1_fixtures_exist_in_committed_set` exists to keep a fixture written
    by an older build, which is the only way to prove v1 still loads at the
    next schema bump. Overwriting it would have destroyed the one thing it is
    for.

    So the two claims are separated. **Durability** of the on-disk format
    against old files stays with `test_trace_file_compat.py`, over the frozen
    fixtures, which load and replay there. **Round-tripping** is this test,
    over a file it writes itself, which can never go stale.

    What a determinism change costs, and what nothing here catches: a trace
    file records a seed and a plan signature but not which sampling model drew
    it, so `rerun_from_trace_file` on a file written by an older engine
    re-draws a different molecule and reports no mismatch. The record column
    `determinism_version` says it for a dataset; the trace file has no such
    field yet.
    """
    seed = 4242
    outcome = vj_compiled.simulator.run(seed=seed)
    emitted = vj_compiled.simulator.trace_file_from(outcome, seed=seed)

    path = tmp_path / f"human_igk_recombine_seed{seed}.trace.json"
    emitted.write_to(str(path))
    loaded = genairr_engine.TraceFile.read_from(str(path))

    assert loaded.seed == seed
    assert loaded.schema_version == emitted.schema_version
    assert loaded.refdata_signature == emitted.refdata_signature
    assert loaded.pass_plan_signature == emitted.pass_plan_signature
    assert len(loaded.trace) == len(emitted.trace)
    assert loaded.trace.choices(), "the run recorded no choices to round-trip"
    for from_disk, in_memory in zip(loaded.trace.choices(), emitted.trace.choices()):
        assert from_disk.address == in_memory.address
        assert from_disk.value == in_memory.value

    rerun = vj_compiled.simulator.rerun_from_trace_file(loaded)
    for redrawn, recorded in zip(rerun.trace().choices(), emitted.trace.choices()):
        assert redrawn.address == recorded.address
        assert redrawn.value == recorded.value
