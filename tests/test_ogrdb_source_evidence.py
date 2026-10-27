"""OGRDB evidence classification: claims, mappings and findings from a source inventory.

``evidence.interpret`` turns the uninterpreted inventory of a verified lock into a new
inventory that adds claims, mappings and findings and leaves the records, candidates,
mentioned identities, observations, support and earlier findings untouched.  Claims carry
honest support, derivation and status: coordinate claims are accepted only inside the
profile's validated scope; boundaries, functionality and J phase independence stay
unresolved; provider endpoints are never repaired; all exports of one provider set are
one evidence lineage.

Expected values come from the independently authored fixture files
(``expectations.json``, ``interpretation_evidence.json``,
``human_igh_release10_candidate_selection.json``), from hand-checked literals, or from
reading raw payload fields with the standard library; never from the code under test.
"""
from __future__ import annotations

import pytest

from GenAIRR._artifact_io import canonical_json
from GenAIRR._reference_sources import evidence
from GenAIRR._reference_sources import model as m
from GenAIRR._reference_sources import ogrdb
from GenAIRR._reference_sources import profiles
from tests.test_ogrdb_source_coordinates import interpreted
from tests.test_ogrdb_source_inventory import (
    IGL, SETS, _key, build, error, inventory_of, pinned,
)

UNTOUCHED = ("records", "candidates", "mentioned_identities", "observations", "support")


# ── the entry point ──────────────────────────────────────────────────


def test_interpret_rejects_arguments_of_the_wrong_type():
    lock, inventory = inventory_of(IGL)
    blobs = pinned(IGL)
    for args in ((lock.to_dict(), blobs, inventory), (lock, blobs, inventory.to_dict()), (lock, [], inventory),
                 (lock, {"airr": "text"}, inventory)):
        with pytest.raises(TypeError):
            evidence.interpret(*args)


def test_the_inventory_must_belong_to_the_lock():
    lock, _ = inventory_of(IGL)
    _, other = inventory_of(IGL, ("airr", "airr_ex"))
    found = error("digest_mismatch", evidence.interpret, lock, pinned(IGL), other)
    assert found.path == "/source_lock_hash"


def test_the_inventory_must_be_the_uninterpreted_inventory_of_these_bytes():
    lock, inventory = inventory_of(IGL)
    blobs = pinned(IGL)
    fields = {name: getattr(inventory, name) for name in ("source_lock_hash", "adapter") + UNTOUCHED}
    trimmed = m.SourceInventory.assemble(**fields, claims=(), mappings=(), findings=inventory.findings[1:])
    assert error("identity_mismatch", evidence.interpret, lock, blobs, trimmed).path == "/inventory"
    _, done = interpreted(IGL)
    assert error("identity_mismatch", evidence.interpret, lock, blobs, done).path == "/inventory"


def test_the_bytes_must_be_the_locked_bytes():
    lock, inventory = inventory_of(IGL)
    blobs = pinned(IGL)
    blobs["gapped"] = blobs["gapped"] + b"\n"
    assert error("digest_mismatch", evidence.interpret, lock, blobs, inventory).path == "/assets/gapped"
    blobs = pinned(IGL)
    del blobs["gapped_ex"]
    assert error("missing_asset", evidence.interpret, lock, blobs, inventory).path == "/assets/gapped_ex"


def test_the_profile_is_the_one_the_lock_names_and_is_never_chosen_from_data(monkeypatch):
    lock, inventory = inventory_of(IGL)
    monkeypatch.setattr(profiles, "PROFILES", {})
    error("unsupported_profile", evidence.interpret, lock, pinned(IGL), inventory)


@pytest.mark.parametrize("source", SETS, ids=_key)
def test_interpretation_adds_claims_mappings_and_findings_and_keeps_the_inventory_stage(source):
    lock, before = inventory_of(source)
    _, after = interpreted(source)
    for name in UNTOUCHED:
        assert getattr(after, name) == getattr(before, name), name
    assert set(before.findings) <= set(after.findings)
    assert before.claims == () and before.mappings == ()
    assert after.claims and after.mappings
    assert ogrdb.account(lock, after) == ogrdb.account(lock, before)
    assert ogrdb.comparison(lock, after) == ogrdb.comparison(lock, before)
    m.check_inventory(after, lock)
    for item in after.claims + after.mappings + after.findings:
        assert item.source_lock_hash == m.lock_hash(lock)


def test_interpretation_is_deterministic_and_round_trips():
    lock, inventory = build(IGL, pinned(IGL))
    first = evidence.interpret(lock, pinned(IGL), inventory)
    second = evidence.interpret(lock, pinned(IGL), inventory)
    assert first == second
    assert m.inventory_hash(first) == m.inventory_hash(second)
    assert m.SourceInventory.from_json(canonical_json(first.to_dict())) == first
    assert first == interpreted(IGL)[1]
