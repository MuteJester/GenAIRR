"""The recipe document: what it may say and what it must refuse.
"""
from __future__ import annotations
import copy
import json
import pytest
from GenAIRR._artifact_io import canonical_json, strict_loads, typed_digest
from GenAIRR._reference_sources import model as m
from tests.support.fixtures.reference_source_models import ADAPTER, CLAIM_ID, FINDING_ID, IGH, ROLE, ROOT, SELECTED, SOURCES, _RECORDS, _find, _load_json, _pointer, _raises, _records, _reseal, _source_set, decision, lock_dict, recipe_dict, selection_entry

def test_recipe_keys_are_frozen():
    recipe = m.Recipe.from_dict(recipe_dict())
    assert tuple(sorted(recipe.to_dict())) == tuple(sorted((
        "artifact", "schema_version", "source_lock_hash", "adapter", "projection", "reference", "selection",
        "decisions", "requires", "defaults_profile")))
    assert tuple(recipe.to_dict()["adapter"]) == ("id", "version", "coordinate_profile")
    assert tuple(recipe.to_dict()["projection"]) == ("id", "version")
    assert sorted(recipe.to_dict()["selection"][0]) == ["label", "reason", "record", "segment", "sequence_sha256"]


def test_recipe_round_trips_and_is_byte_stable():
    d = recipe_dict()
    recipe = m.Recipe.from_dict(d)
    assert recipe.to_dict() == d
    text = canonical_json(d)
    assert canonical_json(m.Recipe.from_json(text).to_dict()) == text
    assert m.Recipe.from_dict(recipe.to_dict()) == recipe


def test_recipe_reference_matches_the_annotation_sidecar_reference_shape():
    from GenAIRR._reference.model import ReferenceHeader
    recipe = m.Recipe.from_dict(recipe_dict())
    header = ReferenceHeader.from_dict(recipe.reference.to_dict())
    assert (header.locus, header.chain_type, header.taxon.taxonomy_id) == ("IGH", "vdj", 9606)


def test_recipe_digest_is_computed_over_the_complete_recipe_and_never_stored_in_it():
    d = recipe_dict()
    recipe = m.Recipe.from_dict(d)
    assert m.recipe_hash(recipe) == typed_digest(d, "genairr_reference_source_recipe_v1")
    assert m.recipe_hash(recipe) not in canonical_json(recipe.to_dict())
    for key in ("commitment", "recipe_hash", "digest"):
        bad = copy.deepcopy(d)
        bad[key] = m.recipe_hash(recipe)
        with _raises("malformed", ""):
            m.Recipe.from_dict(bad)


def test_recipe_selection_order_is_preserved_as_authored():
    forward = m.Recipe.from_dict(recipe_dict())
    backward = recipe_dict()
    backward["selection"].reverse()
    backward = m.Recipe.from_dict(backward)
    assert [e.label for e in backward.selection] == list(reversed(SELECTED))
    assert m.recipe_hash(forward) != m.recipe_hash(backward)


def test_recipe_selects_distinct_identities_with_identical_dna():
    recipe = m.Recipe.from_dict(recipe_dict())
    by_label = {e.label: e for e in recipe.selection}
    assert by_label["IGHV1-69*01"].sequence_sha256 == by_label["IGHV1-69D*01"].sequence_sha256
    assert by_label["IGHV1-69*01"].record != by_label["IGHV1-69D*01"].record


def test_recipe_rejects_two_alternatives_for_one_executable_label():
    d = recipe_dict()
    d["selection"].append(selection_entry("IGHV1-58*03", fmt="airr_ex"))
    with _raises("unresolved_selection", "/selection"):
        m.Recipe.from_dict(d)
    d = recipe_dict()
    d["selection"].append(dict(d["selection"][0], label="IGHV1-58*99"))
    with _raises("unresolved_selection", "/selection"):
        m.Recipe.from_dict(d)


@pytest.mark.parametrize("edit,kind,path", [
    (lambda e: e.update(segment="C"), "malformed", "/selection/0/segment"),
    (lambda e: e.update(label="IGHV1-58 *03"), "malformed", "/selection/0/label"),
    (lambda e: e.update(label=""), "malformed", "/selection/0/label"),
    (lambda e: e["record"].update(pointer="GermlineSet/0"), "malformed", "/selection/0/record/pointer"),
    (lambda e: e["record"].update(pointer="/GermlineSet/~2"), "malformed", "/selection/0/record/pointer"),
    (lambda e: e["record"].update(pointer="/" + "a" * m.MAX_POINTER), "resource_limit", "/selection/0/record/pointer"),
    (lambda e: e["record"].update(asset="../airr"), "malformed", "/selection/0/record/asset"),
    (lambda e: e["record"].update(asset="IGH_VDJ.json"), "malformed", "/selection/0/record/asset"),
    (lambda e: e.update(sequence_sha256="sha256:abc"), "malformed", "/selection/0/sequence_sha256"),
    (lambda e: e.update(reason=" "), "malformed", "/selection/0/reason"),
    (lambda e: e.update(reason=None), "malformed", "/selection/0/reason"),
    (lambda e: e.update(reason="x" * (m.MAX_TEXT + 1)), "resource_limit", "/selection/0/reason"),
    (lambda e: e.update(predicate="all_functional"), "malformed", "/selection/0"),
    (lambda e: e.pop("reason"), "malformed", "/selection/0"),
])
def test_recipe_selection_entries_are_strict(edit, kind, path):
    d = recipe_dict()
    edit(d["selection"][0])
    with _raises(kind, path):
        m.Recipe.from_dict(d)


def test_recipe_selection_segments_must_fit_the_chain_type():
    lock = lock_dict("9606.IGKappa_VJ")
    d = recipe_dict(lock, selected=())
    d["selection"] = [dict(selection_entry("IGHD1-1*01"), record={"asset": "airr", "pointer": _pointer(0)})]
    with _raises("identity_mismatch", "/selection/0/segment"):
        m.Recipe.from_dict(d)


def test_recipe_null_and_missing_values_are_distinct():
    d = recipe_dict()
    d["reference"]["build_date"] = None
    d["reference"]["version"] = None
    recipe = m.Recipe.from_dict(d)
    assert recipe.reference.build_date is None and recipe.reference.version is None
    for key in ("build_date", "version"):
        d = recipe_dict()
        del d["reference"][key]
        with _raises("malformed", "/reference"):
            m.Recipe.from_dict(d)
    for bad_date in ("today", "2026-02-30", "2026-10-22T00:00:00", 20261022):
        d = recipe_dict()
        d["reference"]["build_date"] = bad_date
        with _raises("malformed", "/reference/build_date"):
            m.Recipe.from_dict(d)
    d = recipe_dict()
    d["reference"]["name"] = None
    with _raises("malformed", "/reference/name"):
        m.Recipe.from_dict(d)


@pytest.mark.parametrize("edit,kind,path", [
    (lambda d: d.update(schema_version=2), "unsupported_schema", "/schema_version"),
    (lambda d: d.update(schema_version=True), "malformed", "/schema_version"),
    (lambda d: d.update(artifact="genairr_reference_source_lock"), "malformed", "/artifact"),
    (lambda d: d["adapter"].update(version=2), "unsupported_profile", "/adapter"),
    (lambda d: d["projection"].update(version=2), "unsupported_profile", "/projection"),
    (lambda d: d["projection"].update(id="full_span"), "unsupported_profile", "/projection"),
    (lambda d: d.update(defaults_profile="genairr.unfitted_defaults.v2"), "unsupported_profile", "/defaults_profile"),
    (lambda d: d.update(source_lock_hash=None), "malformed", "/source_lock_hash"),
    (lambda d: d["reference"].update(locus="IGH", chain_type="vj"), "malformed", "/reference/chain_type"),
    (lambda d: d["reference"]["taxon"].update(taxonomy_id="9606"), "malformed", "/reference/taxon/taxonomy_id"),
])
def test_recipe_rejects_unknown_schema_and_profiles(edit, kind, path):
    d = recipe_dict()
    edit(d)
    with _raises(kind, path):
        m.Recipe.from_dict(d)


@pytest.mark.parametrize("edit,kind,path", [
    (lambda r: r.append({"capability": "all_functional", "evidence": "source"}), "malformed", "/requires/1/capability"),
    (lambda r: r.append({"capability": "canonical_productive", "evidence": "independent"}), "malformed",
     "/requires/1/evidence"),
    (lambda r: r.append({"capability": "observed_span_simulation", "evidence": "source"}), "malformed", "/requires"),
    (lambda r: r[0].update(level="source"), "malformed", "/requires/0"),
])
def test_recipe_requirements_use_the_runtime_vocabulary(edit, kind, path):
    d = recipe_dict()
    edit(d["requires"])
    with _raises(kind, path):
        m.Recipe.from_dict(d)


def test_recipe_requirements_keep_authored_order():
    d = recipe_dict()
    d["requires"] = [{"capability": "source_annotated_junction", "evidence": "source"},
                     {"capability": "observed_span_simulation", "evidence": "assertion"}]
    recipe = m.Recipe.from_dict(d)
    assert [r.capability for r in recipe.requires] == ["source_annotated_junction", "observed_span_simulation"]


def test_recipe_decisions_support_only_the_three_actions():
    d = recipe_dict()
    d["decisions"] = [decision(), decision("withhold_claim"),
                      decision("accept_claim", claims=["claim:" + "c" * 64])]
    recipe = m.Recipe.from_dict(d)
    assert [x.action for x in recipe.decisions] == ["select_alternative", "withhold_claim", "accept_claim"]
    for action in ("edit_bases", "extend_terminal", "ignore_warnings", "guess_coordinates", "exec"):
        d["decisions"] = [decision(action=action)]
        with _raises("malformed", "/decisions/0/action"):
            m.Recipe.from_dict(d)


@pytest.mark.parametrize("entry,path", [
    (decision(findings=[]), "/decisions/0/findings"),
    (decision("withhold_claim", claims=[]), "/decisions/0/claims"),
    (decision("accept_claim", evidence=[]), "/decisions/0/evidence"),
    (decision("select_alternative", claims=[CLAIM_ID]), "/decisions/0/claims"),
    (decision(findings=["finding:" + "B" * 64]), "/decisions/0/findings/0"),
    (decision(findings=[FINDING_ID, FINDING_ID]), "/decisions/0/findings"),
    (decision("withhold_claim", claims=["IGHV1-58*03/anchor"]), "/decisions/0/claims/0"),
    (decision(sequence_sha256=None), "/decisions/0/sequence_sha256"),
    (decision(reason=""), "/decisions/0/reason"),
    (decision(evidence=[{"asset": "airr", "pointer": "/x"}, {"asset": "airr", "pointer": "/x"}]),
     "/decisions/0/evidence"),
    (dict(decision(), ignore_warnings=True), "/decisions/0"),
])
def test_recipe_decision_entries_are_strict(entry, path):
    d = recipe_dict()
    d["decisions"] = [entry]
    with _raises("malformed", path):
        m.Recipe.from_dict(d)


def test_recipe_one_claim_is_decided_once():
    d = recipe_dict()
    d["decisions"] = [decision("withhold_claim"), decision("accept_claim")]
    with _raises("malformed", "/decisions"):
        m.Recipe.from_dict(d)


def test_check_recipe_accepts_its_own_lock():
    lock = m.SourceLock.from_dict(lock_dict())
    m.check_recipe(m.Recipe.from_dict(recipe_dict()), lock)


def test_rehashed_recipe_for_another_lock_is_a_digest_mismatch():
    lock = m.SourceLock.from_dict(lock_dict())
    other = lock_dict(origin="local_file")
    with _raises("digest_mismatch", "/source_lock_hash"):
        m.check_recipe(m.Recipe.from_dict(recipe_dict(other)), lock)


@pytest.mark.parametrize("edit,path", [
    (lambda r: r["taxon"].update(scientific_name="Mus musculus"), "/reference/taxon"),
    (lambda r: r["taxon"].update(taxonomy_id=10090), "/reference/taxon"),
    (lambda r: r.update(locus="TRB"), "/reference/locus"),
    (lambda r: r.update(locus="IGK", chain_type="vj"), "/reference/locus"),
])
def test_recipe_cannot_relabel_the_source_identity(edit, path):
    lock = m.SourceLock.from_dict(lock_dict())
    d = recipe_dict()
    edit(d["reference"])
    if d["reference"]["chain_type"] == "vj":
        d["selection"] = [e for e in d["selection"] if e["segment"] != "D"]
    with _raises("identity_mismatch", path):
        m.check_recipe(m.Recipe.from_dict(d), lock)


def test_recipe_adapter_must_equal_the_lock_adapter(monkeypatch):
    from GenAIRR._reference_sources import common
    monkeypatch.setitem(common.SUPPORTED_ADAPTERS, ("ogrdb_germline_set", 2),
                        ("ogrdb", (("ogrdb_api_v2_human_ig_observed_v1", 1),)))
    lock = m.SourceLock.from_dict(lock_dict())
    d = recipe_dict()
    d["adapter"]["version"] = 2
    with _raises("identity_mismatch", "/adapter"):
        m.check_recipe(m.Recipe.from_dict(d), lock)


def test_recipe_records_must_name_locked_record_exports():
    lock = m.SourceLock.from_dict(lock_dict(formats=("airr", "gapped", "gapped_ex")))
    d = recipe_dict(lock.to_dict(), selected=SELECTED[:-1])
    d["selection"].append(selection_entry("IGHJ4*02", fmt="airr_ex"))
    with _raises("identity_mismatch", "/selection/4/record/asset"):
        m.check_recipe(m.Recipe.from_dict(d), lock)
    d = recipe_dict(lock.to_dict())
    d["selection"][0]["record"]["asset"] = "gapped"
    with _raises("identity_mismatch", "/selection/0/record/asset"):
        m.check_recipe(m.Recipe.from_dict(d), lock)
    d = recipe_dict(lock.to_dict())
    d["decisions"] = [decision(evidence=[{"asset": "airr_ex", "pointer": "/x"}])]
    with _raises("identity_mismatch", "/decisions/0/evidence/0/asset"):
        m.check_recipe(m.Recipe.from_dict(d), lock)


def test_changed_recipe_text_changes_the_recipe_identity_even_with_identical_dna():
    first = m.Recipe.from_dict(recipe_dict())
    d = recipe_dict()
    d["selection"][0]["reason"] = "Reviewed again; same supplied sequence"
    second = m.Recipe.from_dict(d)
    assert [e.sequence_sha256 for e in first.selection] == [e.sequence_sha256 for e in second.selection]
    assert m.recipe_hash(first) != m.recipe_hash(second)
