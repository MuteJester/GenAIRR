"""The source inventory built from pinned export bytes.
"""
from __future__ import annotations
import copy
import pytest
from GenAIRR._artifact_io import canonical_json, strict_loads, typed_digest
from GenAIRR._reference_sources import model as m
from tests.support.fixtures.reference_source_models import ADAPTER, CLAIM_ID, FINDING_ID, IGH, INV_RECORDS, Inv, L, ROLE, ROOT, SELECTED, SOURCES, _RECORDS, _blocks, _claim, _fasta, _find, _gapped_index, _inventory_dict, _load_json, _nested, _pointer, _raises, _recipe, _records, _reseal, _sorted, _source_set, decision, inv, inventory, lock_dict, recipe_dict, selection_entry

pytestmark = pytest.mark.behavior

def test_inventory_round_trips_and_is_byte_stable(inventory):
    d = inventory.to_dict()
    again = m.SourceInventory.from_dict(d)
    assert again == inventory
    text = canonical_json(d)
    assert canonical_json(m.SourceInventory.from_json(text).to_dict()) == text
    assert d["artifact"] == "genairr_reference_source_inventory" and d["schema_version"] == 1
    assert m.inventory_hash(inventory) == typed_digest(d, "genairr_reference_source_inventory_v1")


def test_inventory_assemble_puts_every_list_in_canonical_order(inv, inventory):
    shuffled = inv.build(records=list(reversed(inv.records)), claims=list(reversed(inv.claims)),
                         observations=list(reversed(inv.observations)), findings=list(reversed(inv.findings)))
    assert shuffled == inventory
    with _raises("malformed", "/records"):
        m.SourceInventory(**{**_fields_of(inventory), "records": tuple(reversed(inventory.records))})


def _fields_of(inventory):
    return {k: getattr(inventory, k) for k in ("source_lock_hash", "adapter", "records", "candidates",
                                                "mentioned_identities", "observations", "support", "claims",
                                                "mappings", "findings")}


def test_inventory_records_are_occurrences_keyed_by_asset_and_pointer(inventory):
    ids = [r.id for r in inventory.records]
    assert len(ids) == len(set(ids)) == len(INV_RECORDS)
    by_provider = {}
    for record in inventory.records:
        by_provider.setdefault(record.provider_id, []).append(record.locator.asset)
    # The same provider id in both exports stays two record occurrences.
    assert sorted(by_provider["OGRDB:A00014"]) == ["airr", "airr_ex"]


def test_inventory_counts_keep_occurrences_candidates_and_mentions_distinct(inventory):
    counts = inventory.export_counts("airr_ex")
    assert counts == {"record_occurrences": {"V": 2, "D": 1, "J": 1, "unassigned": 0},
                      "sequence_candidates": {"V": 2, "D": 1, "J": 1, "unassigned": 0},
                      "mentioned_identities": {"V": 1, "D": 0, "J": 0, "unassigned": 0}}
    assert inventory.export_counts("airr")["mentioned_identities"]["V"] == 0
    assert inventory.export_counts("airr")["record_occurrences"]["V"] == 3


def test_identical_dna_identities_cannot_cite_each_others_evidence(inventory, inv):
    own, twin = inv.cand["airr", "IGHV1-69D*01"], inv.cand["airr", "IGHV1-69*01"]
    assert inventory.candidate(own).sha256 == inventory.candidate(twin).sha256
    assert inventory.candidate(own).record != inventory.candidate(twin).record
    own_cdr3 = inv.obs_at(inv.rec["airr", "IGHV1-69D*01"], "v_gene_delineations/0/cdr3_start")
    twin_cdr3 = inv.obs_at(inv.rec["airr", "IGHV1-69*01"], "v_gene_delineations/0/cdr3_start")
    inv.build(claims=inv.claims + [inv.claim(own, "anchor", 285, [own_cdr3], system="coding_sequence")])
    with _raises("malformed") as caught:
        inv.build(claims=inv.claims + [inv.claim(own, "anchor", 285, [twin_cdr3], system="coding_sequence")])
    assert caught.error.path.startswith("/claims/") and caught.error.path.endswith("/evidence/0")
    assert "IGHV1-69" in caught.error.message or inv.rec["airr", "IGHV1-69D*01"].id in caught.error.message


def test_evidence_roots_belong_to_the_claim_subject(inv):
    target = inv.cand["airr", "IGHV1-69*01"]
    foreign = inv.obs_at(inv.rec["airr", "IGHV1-58*03"], "v_gene_delineations/0/cdr3_start")
    with _raises("malformed") as caught:
        inv.build(claims=inv.claims + [inv.claim(target, "anchor", inv.anchor, [foreign], system="coding_sequence")])
    assert caught.error.path.endswith("/evidence/0")
    # A set-level field is not evidence about one allele.
    with _raises("malformed"):
        inv.build(claims=inv.claims + [inv.claim(target, "anchor", inv.anchor, [L("airr", "/Info/version")],
                                                 system="coding_sequence")])
    # The same-label record of the other export is citable only by a variant comparison.
    ex, cand = inv.cand["airr_ex", "IGHV1-58*03"], inv.cand["airr", "IGHV1-58*03"]
    with _raises("malformed"):
        inv.build(claims=inv.claims + [inv.claim(ex, "boundary_5", "complete", [L("airr", cand.pointer)])])
    inv.build(claims=inv.claims + [inv.claim(ex, "boundary_3", "partial", [L("airr", cand.pointer)],
                                             derivation="variant_comparison")])
    other_label = inv.cand["airr", "IGHV1-69*01"]
    with _raises("malformed"):
        inv.build(claims=inv.claims + [inv.claim(ex, "boundary_3", "partial", [L("airr", other_label.pointer)],
                                                 derivation="variant_comparison")])


def test_claims_cite_only_mappings_of_their_own_candidate(inv):
    cand = inv.cand["airr", "IGHV1-58*03"]
    ex_mapping = [x for x in inv.mappings if x.input.asset == "airr_ex"][0]
    claims = [c for c in inv.claims if not (c.subject == cand and c.derivation == "imgt_alignment")]
    crossed = inv.claim(cand, "anchor", inv.anchor, [ex_mapping.output], system="coding_sequence",
                        derivation="imgt_alignment", mappings=[ex_mapping.id])
    with _raises("malformed") as caught:
        inv.build(claims=claims + [crossed])
    assert caught.error.path.startswith("/claims/") and caught.error.path.endswith("/mappings/0")
    record = inv.rec["airr", "IGHV1-58*03"]
    functional = inv.claim(record, "functionality", "functional", [inv.obs_at(record, "functional")],
                           mappings=[inv.mappings[0].id])
    with _raises("malformed") as caught:
        inv.build(claims=[c for c in inv.claims if c.subject != record] + [functional])
    assert caught.error.path.endswith("/mappings/0")


def test_accept_claim_evidence_is_bound_to_the_claim_subject(inventory, inv):
    jcand, jrec = inv.cand["airr", "IGHJ4*02"], inv.rec["airr", "IGHJ4*02"]
    phase = _claim(inventory, "coding_phase", jcand)
    sha = inventory.candidate(jcand).sha256
    for evidence in ({"asset": "airr", "pointer": "/Info/version"},
                     inv.obs_at(inv.rec["airr", "IGHV1-58*03"], "functional").to_dict()):
        entry = decision("accept_claim", claims=[phase.id], sequence_sha256=sha, evidence=[evidence])
        with _raises("identity_mismatch", "/decisions/0/evidence/0"):
            m.check_decisions(_recipe(inventory, decisions=[entry]), inventory)
    entry = decision("accept_claim", claims=[phase.id], sequence_sha256=sha,
                     evidence=[inv.obs_at(jrec, "j_codon_frame").to_dict()])
    m.check_decisions(_recipe(inventory, decisions=[entry]), inventory)


def test_phase_value_and_phase_independence_are_separate_claims(inventory, inv):
    jcand = inv.cand["airr", "IGHJ4*02"]
    phase = [c for c in inventory.claims if c.subject == jcand and c.predicate == "coding_phase"]
    independence = [c for c in inventory.claims if c.subject == jcand and c.predicate == "phase_independence"]
    assert [(c.value.get(), c.status, c.derivation) for c in phase] == [(2, "accepted", "stated")]
    assert [(c.value.get(), c.status, c.support) for c in independence] == [(None, "unresolved", "none")]
    vcand = inv.cand["airr", "IGHV1-58*03"]
    (v_independence,) = [c for c in inventory.claims if c.subject == vcand and c.predicate == "phase_independence"]
    assert (v_independence.value.get(), v_independence.status) == ("anchor_derived", "accepted")


def test_claim_and_finding_ids_are_content_derived_and_verified(inventory):
    d = _inventory_dict(inventory)
    claim = d["claims"][0]
    body = {k: v for k, v in claim.items() if k != "id"}
    assert claim["id"] == "claim:" + typed_digest(body, "genairr_reference_source_claim_id_v1").split(":")[1]
    d["claims"][0]["id"] = "claim:" + "0" * 64
    with _raises("digest_mismatch", "/claims/0/id"):
        m.SourceInventory.from_dict(d)
    d = _inventory_dict(inventory)
    (index,) = [i for i, f in enumerate(d["findings"]) if f["code"] == "paralog_only_identity"]
    d["findings"][index]["observations"]["label"] = "IGHV3-23D*01"
    with _raises("digest_mismatch", f"/findings/{index}/id"):
        m.SourceInventory.from_dict(d)
    d = _inventory_dict(inventory)
    d["mappings"][0]["blocks"][0]["output"][0] += 0
    d["mappings"][0]["transform"] = "remove_imgt_gaps"
    assert m.SourceInventory.from_dict(d) == inventory


@pytest.mark.parametrize("edit,kind,path_prefix", [
    (lambda d: d["records"].append(copy.deepcopy(d["records"][0])), "malformed", "/records"),
    (lambda d: d["records"][0].update(label=None), "malformed", "/records/0/label"),
    (lambda d: d["records"][0].update(segment=None), "malformed", "/records/0/segment"),
    (lambda d: d["records"][0].update(status="malformed"), "malformed", "/records/0/problem"),
    (lambda d: d["records"][0].update(problem="unreadable"), "malformed", "/records/0/problem"),
    (lambda d: d["records"][0].update(taxonomy_id=True), "malformed", "/records/0/taxonomy_id"),
    (lambda d: d["records"][0]["locator"].update(asset="gapped"), "malformed", "/records/0/locator"),
    (lambda d: d["candidates"][0].update(sha256="sha256:" + "0" * 64), "digest_mismatch", "/candidates/0/sha256"),
    (lambda d: d["candidates"][0].update(sequence=d["candidates"][0]["sequence"].lower(),
                                         sha256=m.sequence_digest(d["candidates"][0]["sequence"].lower())),
     "malformed", "/candidates/0/sequence"),
    (lambda d: d["candidates"][0]["record"].update(pointer="/GermlineSet/0/allele_descriptions/1"), "malformed",
     "/candidates/0"),
    (lambda d: d["observations"][0].update(presence="absent"), "malformed", "/observations/0/value"),
    (lambda d: d["observations"][0].update(interpretation="trusted"), "malformed",
     "/observations/0/interpretation"),
    (lambda d: d["observations"][0]["record"].update(asset="airr_ex"), "malformed", "/observations/0"),
    (lambda d: d["support"][0].update(status="verified"), "malformed", "/support/0/status"),
    (lambda d: d["mentioned_identities"][0].update(label="IGHV1-69*01"), "malformed", "/mentioned_identities/0"),
    (lambda d: d["mentioned_identities"][0].update(segment="D"), "malformed", "/mentioned_identities/0"),
    (lambda d: d["mentioned_identities"][0].update(mentions=[]), "malformed", "/mentioned_identities/0/mentions"),
    (lambda d: d["mentioned_identities"][0]["mentions"][0].update(relationship="ortholog"), "malformed",
     "/mentioned_identities/0/mentions/0/relationship"),
])
def test_inventory_items_are_strict(inventory, edit, kind, path_prefix):
    d = _inventory_dict(inventory)
    edit(d)
    with _raises(kind) as caught:
        m.SourceInventory.from_dict(d)
    assert caught.error.path.startswith(path_prefix), caught.error


def test_inventory_keeps_malformed_and_unsupported_records(inv):
    loc = L("airr", _pointer(0))
    odd = m.SourceRecord(locator=loc, status="unsupported", problem="sequence_type 'C' is not a V, D or J segment",
                         label="IGHC*01 (constant)", provider_id=None, provider_ref=None, segment=None,
                         taxonomy_id=None, locus="IGH")
    inventory = inv.build(records=inv.records + [odd])
    assert inventory.record(loc).status == "unsupported"
    assert inventory.export_counts("airr")["record_occurrences"]["unassigned"] == 1
    with _raises("malformed", "/label"):
        m.SourceRecord(locator=loc, status="parsed", problem=None, label="IGHC*01 (constant)", provider_id=None,
                       provider_ref=None, segment="V", taxonomy_id=None, locus=None)


def _rebuild_claim(inv, index, **changes):
    claims = list(inv.claims)
    fields = {k: getattr(claims[index], k) for k in ("source_lock_hash", "subject", "predicate", "value", "coordinate_system", "support",
                                                      "derivation", "status", "evidence", "mappings")}
    fields["value"] = fields["value"].get()
    fields.update(changes)
    claims[index] = m.Claim(**fields)
    return inv.build(claims=claims)


def _claim_index(inv, predicate, label="IGHV1-58*03", fmt="airr", **match):
    subjects = {inv.cand[fmt, label], inv.rec[fmt, label]}
    for i, claim in enumerate(inv.claims):
        if claim.predicate == predicate and claim.subject in subjects and all(
                getattr(claim, k) == v for k, v in match.items()):
            return i
    raise AssertionError(predicate)


@pytest.mark.parametrize("predicate,label,changes,path", [
    ("phase_independence", "IGHJ4*02", {"value": "independent"}, "/value"),
    ("coding_phase", "IGHJ4*02", {"value": None}, "/value"),
    ("coding_phase", "IGHJ4*02", {"value": 3}, "/value"),
    ("coding_phase", "IGHJ4*02", {"value": True}, "/value"),
    ("coding_phase", "IGHJ4*02", {"support": "none"}, "/support"),
    ("coding_phase", "IGHJ4*02", {"coordinate_system": None}, "/coordinate_system"),
    ("coding_phase", "IGHJ4*02", {"coordinate_system": "imgt_gapped"}, "/coordinate_system"),
    ("coding_phase", "IGHJ4*02", {"status": "candidate"}, "/support"),
    ("coding_phase", "IGHJ4*02", {"evidence": ()}, "/evidence"),
    ("coding_phase", "IGHJ4*02", {"derivation": "guessed"}, "/derivation"),
    ("anchor", "IGHJ4*02", {"value": -1}, "/value"),
    ("anchor", "IGHJ4*02", {"derivation": "anchor_derived"}, "/derivation"),
    ("functionality", "IGHJ4*02", {"value": "unknown"}, "/value"),
    ("functionality", "IGHJ4*02", {"value": False}, "/value"),
    ("functionality", "IGHJ4*02", {"coordinate_system": "coding_sequence"}, "/coordinate_system"),
    ("extension", "IGHV1-58*03", {"value": {"side": "3", "length": 0, "matches_record_context": True}},
     "/value/length"),
    ("extension", "IGHV1-58*03", {"value": {"side": "3", "length": 1}}, "/value"),
])
def test_claim_values_and_vocabularies_are_strict(inv, predicate, label, changes, path):
    fmt = "airr_ex" if predicate == "extension" else "airr"
    index = _claim_index(inv, predicate, label, fmt)
    with _raises("malformed") as caught:
        _rebuild_claim(inv, index, **changes)
    assert caught.error.path.endswith(path), caught.error


def test_claims_are_checked_against_their_subjects_and_evidence(inv):
    jcand = inv.cand["airr", "IGHJ4*02"]
    dcand = inv.cand["airr", "IGHD1-1*01"]
    cdr3 = inv.obs_at(inv.rec["airr", "IGHV1-58*03"], "v_gene_delineations/0/cdr3_start")
    cases = [
        # No universal D coding phase, and no anchors on D.
        (inv.claim(dcand, "coding_phase", 0, [L("airr", dcand.pointer)], system="coding_sequence"), "malformed"),
        (inv.claim(dcand, "anchor", 3, [L("airr", dcand.pointer)], system="coding_sequence"), "malformed"),
        # v_subregions only on V; functionality only on a record.
        (inv.claim(jcand, "v_subregions", {"FWR1": [0, 10]}, [L("airr", jcand.pointer)], system="coding_sequence",
                   derivation="imgt_alignment", mappings=[inv.mappings[0].id]), "malformed"),
        (inv.claim(jcand, "functionality", "functional", [L("airr", jcand.pointer)]), "malformed"),
        # Anchor beyond the candidate's coding sequence.
        (inv.claim(jcand, "anchor", 46, [L("airr", jcand.pointer)], system="coding_sequence"), "malformed"),
        # Evidence roots are observed source fields.
        (inv.claim(jcand, "anchor", 14, [L("airr", jcand.pointer + "/nowhere")], system="coding_sequence"),
         "malformed"),
        # Mappings must exist; alignment-derived claims need one.
        (inv.claim(jcand, "anchor", 14, [L("airr", jcand.pointer)], system="coding_sequence",
                   mappings=["mapping:" + "0" * 64]), "malformed"),
        # Unknown subject.
        (inv.claim(L("airr", "/GermlineSet/0/allele_descriptions/0/coding_sequence"), "anchor", 1, [cdr3],
                   system="coding_sequence"), "malformed"),
    ]
    for claim, kind in cases:
        with _raises(kind):
            inv.build(claims=inv.claims + [claim])
    with _raises("malformed", "/mappings"):
        inv.claim(inv.cand["airr", "IGHV1-58*03"], "anchor", inv.anchor, [cdr3], system="coding_sequence",
                  derivation="imgt_alignment")


def test_claim_statuses_express_agreement_and_conflict(inv):
    vcand = inv.cand["airr", "IGHV1-58*03"]
    cdr3 = inv.obs_at(inv.rec["airr", "IGHV1-58*03"], "v_gene_delineations/0/cdr3_start")
    disagreeing = inv.claim(vcand, "anchor", inv.anchor + 3, [L("airr", vcand.pointer)], system="coding_sequence",
                            derivation="unique_motif")
    with _raises("malformed"):
        inv.build(claims=inv.claims + [disagreeing])
    lone_conflict = inv.claim(vcand, "boundary_5", "complete", [cdr3], status="conflicting")
    with _raises("malformed"):
        inv.build(claims=inv.claims + [lone_conflict])
    pair = [inv.claim(vcand, "boundary_5", "complete", [cdr3], status="conflicting"),
            inv.claim(vcand, "boundary_5", "partial", [L("airr", vcand.pointer)], status="conflicting")]
    assert sum(c.status == "conflicting" for c in inv.build(claims=inv.claims + pair).claims) == 2
    motif = inv.claim(vcand, "anchor", inv.anchor + 3, [L("airr", vcand.pointer)], system="coding_sequence",
                      derivation="unique_motif", status="candidate", support="none")
    inv.build(claims=inv.claims + [motif])
    independent = inv.claim(vcand, "phase_independence", "independent", [cdr3])
    with _raises("malformed"):
        inv.build(claims=[c for c in inv.claims if not (c.subject == vcand and c.predicate == "phase_independence")]
                  + [independent])


def test_mappings_bind_one_record_occurrence_to_one_gapped_record_of_its_own_variant(inv):
    good = inv.mappings[0]
    fields = {k: getattr(good, k) for k in ("source_lock_hash", "kind", "input", "output", "input_sha256", "output_sha256", "input_space",
                                             "output_space", "blocks", "transform", "validation")}
    ex_gapped = inv.mappings[1].output
    with _raises("malformed", "/output"):
        m.Mapping(**{**fields, "output": ex_gapped})
    with _raises("malformed", "/output"):
        m.Mapping(**{**fields, "output": L("gapped", "/GermlineSet/0")})
    with _raises("malformed", "/blocks"):
        m.Mapping(**{**fields, "blocks": ()})
    with _raises("malformed", "/blocks"):
        m.Mapping(**{**fields, "blocks": (m.Block((0, 10), (0, 10)), m.Block((5, 20), (30, 45)))})
    with _raises("malformed", "/output"):
        m.Block((0, 10), (0, 11))
    with _raises("malformed", "/input"):
        m.Block((5, 5), (0, 0))
    with _raises("malformed", "/input"):
        m.Block((0, True), (0, 1))
    with _raises("malformed", "/transform"):
        m.Mapping(**{**fields, "transform": "align_anything"})
    mismatch = m.Mapping(**{**fields, "blocks": (), "validation": "mismatch"})
    assert mismatch.validation == "mismatch"
    with _raises("digest_mismatch"):
        inv.build(mappings=[m.Mapping(**{**fields, "input_sha256": "sha256:" + "1" * 64}), inv.mappings[1]])
    with _raises("malformed"):
        inv.build(mappings=inv.mappings + [mismatch])
    other = inv.cand["airr", "IGHV1-69*01"]
    with _raises("malformed"):
        inv.build(mappings=inv.mappings + [m.Mapping(**{**fields, "input": other, "blocks": (), "validation": "mismatch",
                                                         "input_sha256": inv.build().candidate(other).sha256})])
    short = fields["blocks"][:-1]
    with _raises("malformed"):
        inv.build(mappings=[m.Mapping(**{**fields, "blocks": short}), inv.mappings[1]],
                  claims=[c for c in inv.claims if not c.mappings])


def test_coding_interval_mapping_converts_record_sequence_coordinates(inv):
    index, raw = _find("IGHJ4*02")
    rec, cand = inv.rec["airr", "IGHJ4*02"], inv.cand["airr", "IGHJ4*02"]
    start, end = raw["gene_start"] - 1, raw["gene_end"]
    assert raw["sequence"][start:end] == raw["coding_sequence"]
    seq_loc = L("airr", rec.pointer + "/sequence")
    observation = m.Observation(locator=seq_loc, record=rec, presence="present", value=raw["sequence"],
                                interpretation="interpreted")
    mapping = m.Mapping(source_lock_hash=inv.lock.commitment, kind="coding_interval", input=seq_loc, output=cand,
                        input_sha256=m.text_digest(raw["sequence"]),
                        output_sha256=m.sequence_digest(raw["coding_sequence"]), input_space="record_sequence",
                        output_space="coding_sequence", blocks=(m.Block((start, end), (0, end - start)),),
                        transform="extract_interval", validation="exact")
    inventory = inv.build(observations=inv.observations + [observation], mappings=inv.mappings + [mapping])
    assert mapping in inventory.mappings
    interval = inv.claim(cand, "coding_interval", [start, end], [seq_loc], system="record_sequence",
                         mappings=[mapping.id])
    inv.build(observations=inv.observations + [observation], mappings=inv.mappings + [mapping],
              claims=inv.claims + [interval])
    with _raises("malformed"):
        m.Mapping(source_lock_hash=inv.lock.commitment, kind="coding_interval", input=seq_loc, output=cand,
                  input_sha256=mapping.input_sha256,
                  output_sha256=mapping.output_sha256, input_space="record_sequence", output_space="coding_sequence",
                  blocks=(m.Block((start, end), (1, end - start + 1)),), transform="extract_interval",
                  validation="exact")


def test_findings_are_strict(inv):
    base = inv.findings[0]
    fields = {k: getattr(base, k) for k in ("source_lock_hash", "code", "severity", "category", "subjects", "locators", "observations",
                                             "capabilities")}
    fields["observations"] = fields["observations"].get()
    for change, path in (({"code": "Variant differs"}, "/code"), ({"severity": "fatal"}, "/severity"),
                         ({"category": "biology"}, "/category"), ({"observations": ["x"]}, "/observations"),
                         ({"capabilities": ("canonical_productive", "observed_span_simulation")}, "/capabilities"),
                         ({"capabilities": ("all_functional",)}, "/capabilities/0"),
                         ({"subjects": (), "locators": ()}, "/subjects"),
                         ({"subjects": tuple(reversed(base.subjects))}, "/subjects")):
        with _raises("malformed", path):
            m.Finding(**{**fields, **change})
    unknown = m.Finding(**{**fields, "subjects": (L("airr", _pointer(0)),)})
    with _raises("malformed"):
        inv.build(findings=inv.findings + [unknown])


def test_preserved_values_are_finite_bounded_json(inv):
    loc = L("airr", "/Info/title")
    for value, kind in ((float("nan"), "malformed"), (float("inf"), "malformed"), ({1: "x"}, "malformed"),
                        ((1, 2), "malformed"), ({"\ud800": 1}, "malformed"),
                        (_nested(m.MAX_EMBEDDED_VALUE_DEPTH + 1), "resource_limit"),
                        (["x"] * (m.MAX_VALUE_NODES + 1), "resource_limit")):
        with _raises(kind):
            m.Observation(locator=loc, record=None, presence="present", value=value, interpretation="preserved")
    a = m.Observation(locator=loc, record=None, presence="present", value={"b": 1, "a": [1.5, None]},
                      interpretation="preserved")
    b = m.Observation(locator=loc, record=None, presence="present", value={"a": [1.5, None], "b": 1},
                      interpretation="preserved")
    assert a == b and a.value.get() == {"a": [1.5, None], "b": 1}
    copy_ = a.value.get()
    copy_["a"].append(2)
    assert a.value.get() == {"a": [1.5, None], "b": 1}
    absent = m.Observation(locator=L("airr", "/Info/missing"), record=None, presence="absent", value=None,
                           interpretation="preserved")
    null = m.Observation(locator=L("airr", "/Info/missing"), record=None, presence="present", value=None,
                         interpretation="preserved")
    assert absent != null and absent.to_dict()["presence"] == "absent" and null.to_dict()["presence"] == "present"


def test_check_inventory_against_its_lock(inventory):
    m.check_inventory(inventory, m.SourceLock.from_dict(lock_dict()))
    with _raises("digest_mismatch", "/source_lock_hash"):
        m.check_inventory(inventory, m.SourceLock.from_dict(lock_dict(origin="local_file")))
    narrow = Inv(m.SourceLock.from_dict(lock_dict(formats=("airr", "gapped", "gapped_ex")))).build()
    with _raises("identity_mismatch") as caught:
        m.check_inventory(narrow, m.SourceLock.from_dict(lock_dict(formats=("airr", "gapped", "gapped_ex"))))
    assert "airr_ex" in caught.error.message
