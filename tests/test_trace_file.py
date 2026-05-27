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
from pathlib import Path

import pytest

import GenAIRR as ga
from GenAIRR import _engine as genairr_engine


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

    # Compile the same plan against a *different* refdata (IGL vs IGK).
    different_refdata = ga.Experiment.on("human_igl").recombine().compile()

    with pytest.raises(ValueError, match="refdata signature mismatch"):
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
# Golden fixture for future Option B (trace-injected replay)
# ──────────────────────────────────────────────────────────────────


GOLDEN_DIR = Path(__file__).parent / "golden" / "trace_files"


def test_emit_golden_trace_file_for_future_b(vj_compiled):
    """Emit (and check the round-trip of) a golden trace file the
    follow-up trace-injected-replay slice will consume as a fixture.

    The file is committed to source control; this test confirms that
    a fresh build against the same plan + refdata still reproduces
    its trace byte-for-byte. If this test fails after an engine
    change, either the schema bumped (intentional) or determinism
    regressed (a bug); inspect the diff before regenerating.
    """
    seed = 4242
    outcome = vj_compiled.simulator.run(seed=seed)
    tf = vj_compiled.simulator.trace_file_from(outcome, seed=seed)

    # Read or emit the golden artifact.
    GOLDEN_DIR.mkdir(parents=True, exist_ok=True)
    golden_path = GOLDEN_DIR / f"human_igk_recombine_seed{seed}.trace.json"
    if not golden_path.exists():
        tf.write_to(str(golden_path))
        pytest.skip(f"Emitted fresh golden trace file at {golden_path}")

    # Compare engine output against the on-disk golden.
    loaded = genairr_engine.TraceFile.read_from(str(golden_path))
    assert loaded.seed == seed
    assert loaded.pass_plan_signature == tf.pass_plan_signature
    assert loaded.refdata_signature == tf.refdata_signature
    assert len(loaded.trace) == len(tf.trace)
    for a, b in zip(loaded.trace.choices(), tf.trace.choices()):
        assert a.address == b.address
        assert a.value == b.value

    # And: re-running from the golden trace file matches both.
    rerun = vj_compiled.simulator.rerun_from_trace_file(loaded)
    for a, b in zip(rerun.trace().choices(), tf.trace.choices()):
        assert a.address == b.address
        assert a.value == b.value
