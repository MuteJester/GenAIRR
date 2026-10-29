"""Selecting against an inventory, and the report a build produces.
"""
from __future__ import annotations
import copy
import pytest
from GenAIRR._artifact_io import canonical_json, typed_digest
from GenAIRR._reference_sources import model as m
from tests.support.fixtures.reference_source_models import Inv, L, SELECTED, _claim, _find, _finding, _pointer, _raises, _recipe, decision, lock_dict, report_dict, selection_entry

pytestmark = pytest.mark.behavior

def test_check_selection_binds_exact_record_and_digest(inventory):
    m.check_selection(_recipe(inventory), inventory)


@pytest.mark.parametrize("edit,kind,path", [
    (lambda d: d["selection"][0].update(label="IGHV1-58*04"), "identity_mismatch", "/selection/0/label"),
    (lambda d: d["selection"][3].update(segment="J"), "identity_mismatch", "/selection/3/segment"),
    (lambda d: d["selection"][0]["record"].update(pointer=_pointer(0)), "unresolved_selection",
     "/selection/0/record"),
    (lambda d: d["selection"][0].update(sequence_sha256=selection_entry("IGHV1-58*03", "airr_ex")["sequence_sha256"]),
     "unresolved_selection", "/selection/0/sequence_sha256"),
])
def test_rehashed_selection_that_disagrees_with_its_record_is_rejected(inventory, edit, kind, path):
    with _raises(kind, path):
        m.check_selection(_recipe(inventory, edit=edit), inventory)


def test_selected_records_must_agree_with_the_source_identity(inv):
    records = list(inv.records)
    index = [i for i, r in enumerate(records) if r.locator == inv.rec["airr", "IGHJ4*02"]][0]
    fields = {k: getattr(records[index], k) for k in ("locator", "status", "problem", "label", "provider_id",
                                                       "provider_ref", "segment", "taxonomy_id", "locus")}
    for change in ({"taxonomy_id": 10090}, {"locus": "IGK"}):
        records[index] = m.SourceRecord(**{**fields, **change})
        with _raises("identity_mismatch", "/selection/4/record"):
            m.check_selection(_recipe(None), inv.build(records=records))
    records[index] = m.SourceRecord(**{**fields, "status": "malformed", "problem": "coding_sequence is not DNA"})
    with _raises("unresolved_selection", "/selection/4/record"):
        m.check_selection(_recipe(None), inv.build(records=records))


def test_check_selection_needs_the_same_lock(inventory):
    other = Inv(m.SourceLock.from_dict(lock_dict(origin="local_file"))).build()
    with _raises("digest_mismatch", "/source_lock_hash"):
        m.check_selection(_recipe(inventory), other)


def test_check_decisions_accepts_current_ids_and_digests(inventory, inv):
    finding = _finding(inventory, "variant_coding_sequence_differs")
    jcand = inv.cand["airr", "IGHJ4*02"]
    j_sha = inventory.candidate(jcand).sha256
    dcand = inv.cand["airr", "IGHD1-1*01"]
    boundary = _claim(inventory, "boundary_3", dcand)
    phase = _claim(inventory, "coding_phase", jcand)
    decisions = [decision(findings=[finding.id]),
                 decision("withhold_claim", claims=[boundary.id], sequence_sha256=inventory.candidate(dcand).sha256),
                 decision("accept_claim", claims=[phase.id], sequence_sha256=j_sha,
                          evidence=[inv.obs_at(inv.rec["airr", "IGHJ4*02"], "j_codon_frame").to_dict()])]
    m.check_decisions(_recipe(inventory, decisions=decisions), inventory)


def test_stale_decisions_are_rejected(inventory, inv):
    finding = _finding(inventory, "variant_coding_sequence_differs")
    jcand = inv.cand["airr", "IGHJ4*02"]
    j_sha = inventory.candidate(jcand).sha256
    phase = _claim(inventory, "coding_phase", jcand)
    cases = [
        (decision(findings=["finding:" + "f" * 64]), "stale_decision", "/decisions/0/findings/0"),
        (decision("withhold_claim", claims=["claim:" + "c" * 64], sequence_sha256=j_sha), "stale_decision",
         "/decisions/0/claims/0"),
        (decision("withhold_claim", claims=[phase.id]), "stale_decision", "/decisions/0/sequence_sha256"),
        (decision(findings=[finding.id], evidence=[{"asset": "airr", "pointer": "/GermlineSet/0/nowhere"}]),
         "stale_decision", "/decisions/0/evidence/0"),
        (decision(findings=[_finding(inventory, "paralog_only_identity").id]), "stale_decision",
         "/decisions/0/sequence_sha256"),
        (decision(findings=[finding.id], sequence_sha256=selection_entry("IGHV1-58*03", "airr_ex")["sequence_sha256"]),
         "stale_decision", "/decisions/0/sequence_sha256"),
    ]
    for entry, kind, path in cases:
        with _raises(kind, path):
            m.check_decisions(_recipe(inventory, decisions=[entry]), inventory)


def test_accepted_claims_must_be_representable_by_the_runtime_reference(inventory, inv):
    ex = inv.cand["airr_ex", "IGHV1-58*03"]
    extension = _claim(inventory, "extension", ex)
    jcand = inv.cand["airr", "IGHJ4*02"]
    independence = _claim(inventory, "phase_independence", jcand)
    for claim, sha in ((extension, inventory.candidate(ex).sha256), (independence, inventory.candidate(jcand).sha256)):
        entry = decision("accept_claim", claims=[claim.id], sequence_sha256=sha,
                         evidence=[claim.evidence[0].to_dict()])
        with _raises("projection", "/decisions/0/claims/0"):
            m.check_decisions(_recipe(inventory, decisions=[entry]), inventory)


def test_decisions_on_unselected_sequences_are_stale(inventory, inv):
    ex = inv.cand["airr_ex", "IGHJ4*02"]
    record = inv.rec["airr_ex", "IGHJ4*02"]
    functionality = _claim(inventory, "functionality", record)
    entry = decision("withhold_claim", claims=[functionality.id], sequence_sha256=inventory.candidate(ex).sha256)
    selected = tuple(label for label in SELECTED if label != "IGHJ4*02")
    with _raises("stale_decision", "/decisions/0/claims/0"):
        m.check_decisions(_recipe(inventory, selected=selected, decisions=[entry]), inventory)


def test_report_round_trips_and_is_byte_stable(report, report_d):
    assert report.to_dict() == report_d
    text = canonical_json(report_d)
    assert canonical_json(m.BuildReport.from_json(text).to_dict()) == text
    assert m.report_hash(report) == typed_digest(report_d, "genairr_reference_source_build_report_v1")


def test_report_exposes_every_record_claim_and_distinct_count(report):
    d = report.to_dict()
    assert len(d["records"]) == len(d["inventory"]["records"])
    assert len(d["claims"]) == len(d["inventory"]["claims"])
    airr_ex = d["counts"]["exports"][1]
    assert airr_ex["record_occurrences"]["V"] == 2 and airr_ex["mentioned_identities"]["V"] == 1
    assert d["counts"]["selected_alleles"] == {"V": 3, "D": 1, "J": 1, "unassigned": 0}
    assert d["counts"]["record_outcomes"] == {"selected": 5, "excluded": 4, "rejected": 0}
    outcomes = {c["outcome"] for c in d["claims"]}
    assert outcomes == {"projected", "withheld", "not_representable", "not_selected"}
    assert {r["exclusion"] for r in d["records"] if r["outcome"] == "excluded"} == {
        "alternative_for_selected_identity", "duplicate_occurrence"}


def _edited_report(report_d, edit):
    d = copy.deepcopy(report_d)
    edit(d)
    return d


def _record_index(d, label, asset):
    for i, r in enumerate(d["records"]):
        rec = [x for x in d["inventory"]["records"] if x["locator"] == r["record"]][0]
        if rec["label"] == label and r["record"]["asset"] == asset:
            return i
    raise AssertionError(label)


def _claim_outcome_index(d, outcome):
    return [i for i, c in enumerate(d["claims"]) if c["outcome"] == outcome][0]


@pytest.mark.parametrize("edit,kind,path", [
    (lambda d: d["hashes"].update(build_inputs="sha256:" + "9" * 64), "digest_mismatch", "/hashes/build_inputs"),
    (lambda d: d["hashes"].update(inventory="sha256:" + "9" * 64), "digest_mismatch", "/hashes/inventory"),
    (lambda d: d["recipe"].update(source_lock_hash="sha256:" + "9" * 64), "digest_mismatch", "/recipe/source_lock_hash"),
    (lambda d: d["recipe"]["reference"].update(locus="TRB"), "identity_mismatch", "/recipe/reference/locus"),
    (lambda d: d["runtime"].update(execution_hash="sha256:XYZ"), "malformed", "/runtime/execution_hash"),
    (lambda d: d["generator"].update(name="NotGenAIRR"), "malformed", "/generator/name"),
    (lambda d: d.update(schema_version=2), "unsupported_schema", "/schema_version"),
    (lambda d: d.update(report_hash="sha256:" + "0" * 64), "malformed", ""),
    (lambda d: d["comparison"].update(complete=False), "malformed", "/comparison/complete"),
    (lambda d: d["comparison"].update(findings=["finding:" + "0" * 64]), "malformed", "/comparison/findings/0"),
    (lambda d: d["counts"]["exports"][0]["record_occurrences"].update(V=209), "malformed", "/counts/exports"),
    (lambda d: d["counts"]["selected_alleles"].update(V=190), "malformed", "/counts/selected_alleles"),
    (lambda d: d["counts"]["record_outcomes"].update(rejected=1), "malformed", "/counts/record_outcomes"),
    (lambda d: d["counts"]["selected_alleles"].update(V=True), "malformed", "/counts/selected_alleles/V"),
    (lambda d: d["counts"]["selected_alleles"].update(V=-1), "malformed", "/counts/selected_alleles/V"),
    (lambda d: d["selection"].pop(), "malformed", "/selection"),
    (lambda d: d["selection"][0].update(allele="V:IGHV1-58*04"), "malformed", "/selection/0/allele"),
    (lambda d: d["selection"][0].update(candidate=d["selection"][1]["candidate"]), "malformed",
     "/selection/0/candidate"),
    (lambda d: d["records"].reverse(), "malformed", "/records/0/record"),
    (lambda d: d["records"][_record_index(d, "IGHV1-58*03", "airr")].update(
        outcome="excluded", allele=None, exclusion="outside_explicit_selection"), "malformed", "/records/"),
    (lambda d: d["records"][_record_index(d, "IGHV1-58*03", "airr_ex")].update(exclusion="duplicate_occurrence"),
     "malformed", "/records/"),
    (lambda d: d["records"][_record_index(d, "IGHV1-69*01", "airr_ex")].update(
        exclusion="alternative_for_selected_identity"), "malformed", "/records/"),
    (lambda d: d["records"][_record_index(d, "IGHV1-69*01", "airr_ex")].update(related=None), "malformed",
     "/records/"),
    (lambda d: d["records"][_record_index(d, "IGHV1-69*01", "airr_ex")].update(exclusion="looked_odd"), "malformed",
     "/records/"),
    (lambda d: d["records"][_record_index(d, "IGHV1-69*01", "airr_ex")].update(outcome="rejected", exclusion=None,
                                                                              related=None), "malformed", "/records/"),
    (lambda d: d["claims"][_claim_outcome_index(d, "not_representable")].update(outcome="projected"), "malformed",
     "/claims/"),
    (lambda d: d["claims"][_claim_outcome_index(d, "not_selected")].update(outcome="projected"), "malformed",
     "/claims/"),
    (lambda d: d["claims"][_claim_outcome_index(d, "projected")].update(outcome="not_selected"), "malformed",
     "/claims/"),
    (lambda d: d["claims"][_claim_outcome_index(d, "withheld")].update(reason=None), "malformed", "/claims/"),
    (lambda d: d["claims"][_claim_outcome_index(d, "withheld")].update(outcome="projected", reason=None), "malformed",
     "/claims/"),
    (lambda d: d["claims"][_claim_outcome_index(d, "projected")].update(decision=0), "malformed", "/claims/"),
    (lambda d: d["claims"][_claim_outcome_index(d, "projected")].update(decision=7), "malformed", "/claims/"),
    (lambda d: d["claims"][_claim_outcome_index(d, "projected")].update(decision=True), "malformed", "/claims/"),
    (lambda d: d["claims"].reverse(), "malformed", "/claims/0/claim"),
    (lambda d: d["capabilities"]["statuses"].reverse(), "malformed", "/capabilities/statuses"),
    (lambda d: d["capabilities"]["statuses"][4].update(available=True), "malformed", "/capabilities/statuses/4"),
    (lambda d: d["capabilities"]["statuses"][6].update(evidence="source"), "malformed",
     "/capabilities/statuses/6/evidence"),
    (lambda d: d["capabilities"].update(phase_support="guessed"), "malformed", "/capabilities/phase_support"),
    (lambda d: d["requirements"][0].update(satisfied=False), "malformed", "/requirements/0/satisfied"),
    (lambda d: d["requirements"][0].update(evidence="source"), "malformed", "/requirements/0"),
    (lambda d: d["requirements"].append({"capability": "original_coding_ends", "evidence": "assertion",
                                         "satisfied": False}), "malformed", "/requirements"),
])
def test_report_rejects_inconsistent_content(report_d, edit, kind, path):
    with _raises(kind) as caught:
        m.BuildReport.from_dict(_edited_report(report_d, edit))
    assert caught.error.path.startswith(path), caught.error


def test_report_withheld_and_accepted_claims_follow_the_decision_ledger(report_d):
    d = report_d
    decided = {c["claim"]: c for c in d["claims"] if c["decision"] is not None}
    actions = [x["action"] for x in d["recipe"]["decisions"]]
    assert sorted((actions[c["decision"]], c["outcome"]) for c in decided.values()) == [
        ("accept_claim", "projected"), ("withhold_claim", "withheld")]
    bad = copy.deepcopy(d)
    (index,) = [i for i, c in enumerate(bad["claims"]) if c["decision"] == 1]
    bad["claims"][index].update(decision=None, reason="unresolved")
    with _raises("malformed", "/claims"):
        m.BuildReport.from_dict(bad)


def test_report_records_rejected_source_records_apart_from_exclusions(inventory, inv):
    loc = L("airr", _pointer(0))
    odd = m.SourceRecord(locator=loc, status="malformed", problem="coding_sequence is not a string", label=None,
                         provider_id=None, provider_ref=None, segment=None, taxonomy_id=None, locus=None)
    _, raw = _find("IGHV1-2*02")
    other = m.SourceRecord(locator=L("airr", _pointer(_find("IGHV1-2*02")[0])), status="parsed", problem=None,
                           label=raw["label"], provider_id=raw["allele_description_id"],
                           provider_ref=raw["allele_description_ref"], segment="V", taxonomy_id=9606, locus="IGH")
    d = report_dict(inv.build(records=inv.records + [odd, other]), inv, decisions=[])
    report = m.BuildReport.from_dict(d)
    by_record = {r.record: r for r in report.records}
    assert by_record[loc].outcome == "rejected" and by_record[loc].detail == "coding_sequence is not a string"
    assert by_record[other.locator].exclusion == "outside_explicit_selection"
    assert report.counts.record_outcomes.to_dict() == {"selected": 5, "excluded": 5, "rejected": 1}
    (index,) = [i for i, r in enumerate(d["records"]) if r["record"] == loc.to_dict()]
    d["records"][index].update(outcome="excluded", exclusion="outside_explicit_selection")
    d["counts"]["record_outcomes"] = {"selected": 5, "excluded": 6, "rejected": 0}
    with _raises("malformed", f"/records/{index}/outcome"):
        m.BuildReport.from_dict(d)


def test_report_incomplete_comparison_lists_the_missing_export():
    lock = m.SourceLock.from_dict(lock_dict(formats=("airr", "gapped")))
    comparison = m.Comparison(complete=False, exports=("airr",), missing=("airr_ex",), findings=())
    comparison.check(lock, None)
    with _raises("malformed", "/complete"):
        m.Comparison(complete=True, exports=("airr",), missing=("airr_ex",), findings=())
