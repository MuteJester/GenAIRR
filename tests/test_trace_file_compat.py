"""Versioned trace-file compatibility guardrails.

This is the **"we did not break persisted science artifacts"** test
suite. Every `.trace.json` file committed under
`tests/golden/trace_files/` is treated as a frozen fixture of the
on-disk contract at the time it was emitted. Each fixture must:

1. **Load** through the current `TraceFile.read_from` — no schema
   version mismatch, no JSON parse error, no missing-field panic.
2. **Validate** every recorded address parses to a built-in
   `ChoiceAddress` (the address vocabulary the trace was written
   against must still be parseable today).
3. **Replay** through `replay_from_trace_file` against a
   reconstructed plan + refdata — the full cursor-driven path,
   including the strict-positional drained-cursor assertion.

When the engine bumps `TRACE_FILE_SCHEMA_VERSION` or
`ADDRESS_SCHEMA_VERSION`, this test is the canary that proves old
fixtures still survive the change. A failure here means either:
  * a fixture must be migrated (with explicit version bump notes), or
  * the engine change was a covert breaking change and should be
    reconsidered.

Two fixtures committed today, both v1:
  - `human_igk_recombine_seed4242.trace.json` — recombine-only plan
  - `full_stack_seed4242.trace.json` — every migrated sampling site
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable

import pytest

import GenAIRR as ga
from GenAIRR import _engine as eng


GOLDEN_DIR = Path(__file__).parent / "golden" / "trace_files"


# Each entry: (fixture filename, callable that reconstructs the
# CompiledExperiment that produced the fixture).
#
# When you add a new committed fixture, register it here so the
# compat suite picks it up. The reconstruction callable encodes the
# plan + refdata identity the fixture was emitted against — a
# fixture that can't be reproduced is a fixture that can't be
# replayed.


def _human_igk_recombine_only():
    return ga.Experiment.on("human_igk").recombine().compile()


def _full_stack_seed4242():
    return (
        ga.Experiment.on("human_igk")
        .recombine()
        .mutate(model="uniform", count=2)
        .pcr_amplify(count=1)
        .sequencing_errors(count=1)
        .polymerase_indels(count=2, insertion_prob=0.5)
        .ambiguous_base_calls(count=1)
        .primer_trim_5prime(length=2)
        .contaminate(prob=0.0)
        .random_strand_orientation(prob=0.5)
        .compile()
    )


GOLDEN_FIXTURES = {
    "human_igk_recombine_seed4242.trace.json": _human_igk_recombine_only,
    "full_stack_seed4242.trace.json": _full_stack_seed4242,
}


def _discover_committed_fixtures() -> Iterable[str]:
    """Yield every `.trace.json` filename committed under the
    golden directory. Each must be registered in GOLDEN_FIXTURES."""
    if not GOLDEN_DIR.exists():
        return []
    return sorted(p.name for p in GOLDEN_DIR.glob("*.trace.json"))


def test_every_committed_fixture_is_registered():
    """A fixture sitting in the golden directory without a
    registered reconstruction callable can't be replayed → can't
    enforce the durability contract. This guards against
    drive-by-committed JSON blobs that no one knows how to
    reproduce."""
    discovered = set(_discover_committed_fixtures())
    registered = set(GOLDEN_FIXTURES.keys())
    unregistered = discovered - registered
    assert not unregistered, (
        f"committed golden fixtures with no reconstruction recipe: {unregistered}; "
        f"register them in GOLDEN_FIXTURES so the compat suite can replay them"
    )


@pytest.mark.parametrize("fixture", sorted(GOLDEN_FIXTURES.keys()))
def test_golden_fixture_loads(fixture):
    """The fixture parses through the current `TraceFile.read_from`
    surface — no schema-version rejection, no missing-field error,
    no JSON parse error."""
    path = GOLDEN_DIR / fixture
    if not path.exists():
        pytest.skip(f"fixture {fixture} not yet committed; emit via test_trace_file or test_trace_replay")
    tf = eng.TraceFile.read_from(str(path))
    # Sanity: the fixture's schema_version is one the engine knows.
    assert tf.schema_version >= 1


@pytest.mark.parametrize("fixture", sorted(GOLDEN_FIXTURES.keys()))
def test_golden_fixture_addresses_all_parse(fixture):
    """Every recorded address in the fixture must parse to a
    built-in `ChoiceAddress`. If a future engine drops a variant
    or renames a spelling without bumping `ADDRESS_SCHEMA_VERSION`,
    this test catches the regression at the file-format boundary
    rather than at replay time."""
    path = GOLDEN_DIR / fixture
    if not path.exists():
        pytest.skip(f"fixture {fixture} not yet committed")
    tf = eng.TraceFile.read_from(str(path))
    # The Rust-side `validate_addresses(allow_custom=False)` walks
    # every record and rejects on the first unknown spelling.
    tf.validate_addresses()


@pytest.mark.parametrize("fixture", sorted(GOLDEN_FIXTURES.keys()))
def test_golden_fixture_replays_through_reconstructed_compiled(fixture):
    """The fixture's trace, re-applied through
    `replay_from_trace_file` against the reconstruction recipe,
    produces an outcome whose own trace matches the fixture's
    record-for-record. This is the durability contract end to end."""
    path = GOLDEN_DIR / fixture
    if not path.exists():
        pytest.skip(f"fixture {fixture} not yet committed")
    tf = eng.TraceFile.read_from(str(path))

    factory = GOLDEN_FIXTURES[fixture]
    compiled = factory()

    outcome = compiled.simulator.replay_from_trace_file(tf)

    # Trace round-trip.
    loaded_choices = tf.trace.choices()
    replayed_choices = outcome.trace().choices()
    assert len(loaded_choices) == len(replayed_choices), (
        f"trace length mismatch for {fixture}: "
        f"fixture={len(loaded_choices)} replayed={len(replayed_choices)}"
    )
    for i, (a, b) in enumerate(zip(loaded_choices, replayed_choices)):
        assert a.address == b.address, (
            f"address divergence in {fixture} at index {i}: "
            f"fixture={a.address!r} replayed={b.address!r}"
        )
        assert a.value == b.value, (
            f"value divergence in {fixture} at index {i} (addr={a.address})"
        )


# ──────────────────────────────────────────────────────────────────
# Schema policy regression: v1 fixtures must remain loadable
# ──────────────────────────────────────────────────────────────────


def test_v1_fixtures_exist_in_committed_set():
    """At least one v1 fixture must be committed at all times. v1
    is the historical baseline of the on-disk contract; losing it
    means we can no longer prove v1→current compatibility."""
    found_v1 = False
    for fixture in _discover_committed_fixtures():
        path = GOLDEN_DIR / fixture
        data = json.loads(path.read_text())
        if data.get("schema_version") == 1:
            found_v1 = True
            break
    assert found_v1, (
        "no v1 trace fixture found; commit a v1 fixture so the "
        "schema-evolution compat test has a baseline to verify "
        "against. Existing fixtures may have been regenerated under "
        "a newer schema version — check git history before deleting them."
    )


def test_v1_fixture_has_no_v2_fields_serialised():
    """Sanity-check the on-disk shape of v1 fixtures: they must
    NOT carry the v2-introduced `refdata_content_hash` or
    `address_schema_version` keys. If a v1 fixture ever ends up
    with v2 fields on disk, something regenerated it incorrectly."""
    for fixture in _discover_committed_fixtures():
        path = GOLDEN_DIR / fixture
        data = json.loads(path.read_text())
        if data.get("schema_version") != 1:
            continue
        assert "refdata_content_hash" not in data, (
            f"{fixture} is tagged schema_version=1 but carries the "
            f"v2 field `refdata_content_hash` — fixture has been "
            f"corrupted by a regeneration under a newer engine"
        )
        assert "address_schema_version" not in data, (
            f"{fixture} is tagged schema_version=1 but carries the "
            f"v2 field `address_schema_version`"
        )
