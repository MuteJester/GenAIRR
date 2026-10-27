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

import dataclasses
import types

import pytest

from GenAIRR._artifact_io import canonical_json
from GenAIRR._reference_sources import evidence
from GenAIRR._reference_sources import model as m
from GenAIRR._reference_sources import ogrdb
from GenAIRR._reference_sources import profiles
from tests.test_ogrdb_source_coordinates import (
    GAPPED, NOTE_LABELS, assert_status, at, candidate_of, claim_of, findings, interpret,
    interpreted, raw_gapped, raw_record, record_of, small,
)
from tests.test_ogrdb_source_inventory import (
    EXPECT, IGH, IGK, IGL, SELECTION, SETS, _key, build, error, inventory_of, pinned, raw_records,
)

EXTENSION_CODES = ("coding_extension_matches_record_context", "coding_extension_without_record_context",
                   "coding_extension_disagrees_with_record_context")
CODING_END_ROWS = {row["label"]: row for row in EXPECT["igl4_airr_ex_coding_end_observations"]}

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


# ── provider coding extensions ───────────────────────────────────────


def _suffixes(source) -> dict:
    return EXPECT["sets"][f"{source[0]}/{source[1]}"]["airr_ex_coding_suffix_beyond_airr"]


def _context_oracle(source, label, extension) -> tuple:
    """The record bases after gene_end, read from the raw airr_ex record (independent oracle)."""
    record = raw_record(source, "airr_ex", label)
    bases = record["sequence"][record["gene_end"]:record["gene_end"] + len(extension)]
    if bases == extension:
        return "matches", bases
    if extension.startswith(bases):
        return "without_record_context", bases
    return "disagrees", bases


def extension_findings(inventory) -> dict:
    found = {}
    for code in EXTENSION_CODES:
        for finding in findings(inventory, code):
            label = finding.observations.get()["label"]
            assert label not in found
            found[label] = (code, finding)
    return found


@pytest.mark.parametrize("source", SETS, ids=_key)
def test_each_provider_extension_is_classified_against_its_own_record_context(source):
    lock, inventory = interpreted(source)
    found = extension_findings(inventory)
    assert set(found) == set(_suffixes(source))
    for label, extension in _suffixes(source).items():
        code, finding = found[label]
        status, bases = _context_oracle(source, label, extension)
        assert code == {"matches": "coding_extension_matches_record_context",
                        "without_record_context": "coding_extension_without_record_context",
                        "disagrees": "coding_extension_disagrees_with_record_context"}[status]
        observations = finding.observations.get()
        assert (observations["extending_export"], observations["side"], observations["extension"],
                observations["record_context"]) == ("airr_ex", "3", extension, bases)
        assert finding.capabilities == ("original_coding_ends",)
        records = {s.asset: s for s in finding.subjects if s.id in {r.id for r in inventory.records}}
        assert set(records) == {"airr", "airr_ex"}
        long = candidate_of(inventory, "airr_ex", label)
        claim = claim_of(inventory, long, "extension")
        assert_status(claim, "accepted", "source", "variant_comparison",
                      {"side": "3", "length": len(extension), "matches_record_context": status == "matches"})
        assert claim.coordinate_system is None
        short = candidate_of(inventory, "airr", label)
        assert long.sequence == short.sequence + extension
        assert short.locator in claim.evidence and long.locator in claim.evidence
        assert [x for x in inventory.claims if x.subject == short.locator and x.predicate == "extension"] == []


def test_igh_extensions_split_15_with_matching_context_8_without_and_none_disagreeing():
    _, inventory = interpreted(IGH)
    found = extension_findings(inventory)
    by_code = {code: {label for label, (c, _) in found.items() if c == code} for code in EXTENSION_CODES}
    assert [len(by_code[code]) for code in EXTENSION_CODES] == [15, 8, 0]
    choices = {c["label"]: c for c in SELECTION["record_choices"]}
    assert set(choices) <= by_code["coding_extension_matches_record_context"]
    for label, choice in choices.items():
        observations = found[label][1].observations.get()
        assert choice["extension_matches_record_context_after_gene_end"] is True
        assert observations["record_context"] == choice["record_context_after_gene_end"]
        assert observations["extension"] == choice["alternative"]["extension_3prime"]
    omitted_by_choice = {o["label"] for o in SELECTION["omitted"]
                         if any(r["rule"] == "explicit_record_choice" for r in o["failed_rules"])}
    assert by_code["coding_extension_without_record_context"] < omitted_by_choice
    assert len(omitted_by_choice) == 9
    for label in by_code["coding_extension_without_record_context"]:
        assert found[label][1].observations.get()["record_context"] == ""


def test_igl_extensions_two_base_class_and_the_separate_disagreeing_extension():
    _, inventory = interpreted(IGL)
    found = extension_findings(inventory)
    assert set(found) == set(CODING_END_ROWS)
    for label, row in CODING_END_ROWS.items():
        code, finding = found[label]
        observations = finding.observations.get()
        assert (observations["extension"], observations["record_context"], observations["gene_end"],
                observations["v_rs_start"]) == (row["coding_extension"], row["sequence_after_gene_end"],
                                                row["gene_end"], row["v_rs_start"])
        if row["finding"] == "coding_extension_disagrees_with_record_context":
            assert code == row["finding"]
            assert observations["extension_end_vs_v_rss"] == "overlaps_v_rss"
        else:
            assert code == "coding_extension_matches_record_context"
            assert observations["extension_end_vs_v_rss"] == "reaches_v_rss_minus_1"
    assert found["IGLV3-9*01"] == ("coding_extension_disagrees_with_record_context", found["IGLV3-9*01"][1])
    assert (found["IGLV3-9*01"][1].observations.get()["extension"],
            found["IGLV3-9*01"][1].observations.get()["record_context"]) == ("CACCC", "CACAG")


def test_the_igk_extension():
    _, inventory = interpreted(IGK)
    found = extension_findings(inventory)
    assert {label: (code, f.observations.get()["extension"]) for label, (code, f) in found.items()} == {
        "IGKV1D-8*04": ("coding_extension_matches_record_context", "TCC")}


def test_an_extension_whose_position_is_not_proven_has_no_record_context():
    def edit(fmt, record):
        if fmt == "airr_ex" and record["label"] == "IGLV10-54*04":
            record["gene_end"] -= 3                  # the declared interval no longer gives the airr sequence
    lock, inventory = interpret(small(IGL, {"IGLV10-54*04", "IGLJ1*01"}, edit=edit))
    ((code, finding),) = extension_findings(inventory).values()
    assert code == "coding_extension_without_record_context"
    observations = finding.observations.get()
    assert observations["record_context"] is None and "does not reproduce" in observations["reason"]
    claim = claim_of(inventory, candidate_of(inventory, "airr_ex", "IGLV10-54*04"), "extension")
    assert_status(claim, "candidate", "none", "variant_comparison",
                  {"side": "3", "length": 2, "matches_record_context": False})


# ── coding ends before a declared V RSS ──────────────────────────────


def coding_end_findings(inventory) -> dict:
    return {f.observations.get()["label"]: f for f in findings(inventory, "coding_end_before_declared_v_rss")}


def test_the_five_igl_coding_ends_before_a_declared_v_rss():
    _, inventory = interpreted(IGL)
    found = coding_end_findings(inventory)
    expected = {label: row for label, row in CODING_END_ROWS.items()
                if row["finding"] == "coding_end_before_declared_v_rss"}
    assert set(found) == set(expected) == set(NOTE_LABELS)
    notes = {n.label: n for n in profiles.PROFILE.historical_notes}
    for label, row in expected.items():
        finding = found[label]
        observations = finding.observations.get()
        assert (observations["gene_end"], observations["v_rs_start"], observations["bases_through_rss_minus_1"]) == (
            row["gene_end"], row["v_rs_start"], row["sequence_after_gene_end"])
        assert observations["gap"] == 2
        assert observations["extension"] == {"export": "airr_ex", "bases": row["coding_extension"],
                                             "matches_record_context": True, "reaches_rss_minus_1": True}
        assert observations["historical_note"] == {"field": "curation", "quotes": list(notes[label].quotes),
                                                   "affirmed_span": "unresolved"}
        raw_support = raw_record(IGL, "airr", label)["unrearranged_support"]
        assert raw_support
        assert observations["genomic_support"] == [
            {"repository": s["repository_name"], "reference": s["repository_ref"], "status": "embedded_unretrieved"}
            for s in raw_support]
        assert {s.asset for s in finding.subjects} == {"airr", "airr_ex"}
        assert finding.capabilities == ("original_coding_ends",)
        for fmt in GAPPED:   # the shorter declared end is kept; nothing is extended or corrected
            interval = claim_of(inventory, candidate_of(inventory, fmt, label), "coding_interval")
            raw_rec = raw_record(IGL, fmt, label)
            if fmt == "airr":
                assert interval.value.get() == [raw_rec["gene_start"] - 1, row["gene_end"]]
            else:
                assert interval.status == "unresolved"
    assert "IGLV3-9*01" not in found


def test_the_igh_two_base_analogue_is_ighv3_7_04():
    _, inventory = interpreted(IGH)
    found = coding_end_findings(inventory)
    assert set(found) == {"IGHV3-7*04"}
    (choice,) = [c for c in SELECTION["record_choices"] if c["label"] == "IGHV3-7*04"]
    observations = found["IGHV3-7*04"].observations.get()
    assert observations["v_rs_start"] == choice["v_rs_start"] == 457
    assert observations["bases_through_rss_minus_1"] == choice["alternative"]["extension_3prime"] == "A"
    assert observations["extension"] == {"export": "airr_ex", "bases": "A", "matches_record_context": True,
                                         "reaches_rss_minus_1": choice["extension_reaches_rss_minus_1"]}
    assert observations["historical_note"] is None
    _, igk = interpreted(IGK)
    assert coding_end_findings(igk) == {}


# ── declared V RSS positions ─────────────────────────────────────────


def _gap_counts(observations) -> dict:
    counts = {"0": observations["adjacent"]} if observations["adjacent"] else {}
    counts.update({gap: len(labels) for gap, labels in observations["gaps"].items()})
    if observations["no_v_rs_start"]:
        counts["no_v_rs_start"] = len(observations["no_v_rs_start"])
    return counts


def _raw_gap_counts(source, fmt) -> dict:
    counts: dict = {}
    for record in raw_records(source, fmt):
        if record["sequence_type"] == "V":
            key = "no_v_rs_start" if record["v_rs_start"] is None else str(record["v_rs_start"] - record["gene_end"] - 1)
            counts[key] = counts.get(key, 0) + 1
    return counts


@pytest.mark.parametrize("source", SETS, ids=_key)
def test_v_rss_spacing_is_counted_per_export(source):
    _, inventory = interpreted(source)
    summaries = {f.observations.get()["asset"]: f for f in findings(inventory, "v_rss_spacing")}
    assert set(summaries) == set(GAPPED)
    expected = EXPECT["sets"][f"{source[0]}/{source[1]}"]["airr_v_rss_gap_counts"]
    assert _gap_counts(summaries["airr"].observations.get()) == expected
    assert _gap_counts(summaries["airr_ex"].observations.get()) == _raw_gap_counts(source, "airr_ex")
    for asset, finding in summaries.items():
        observations = finding.observations.get()
        assert observations["v_records"] == EXPECT["sets"][f"{source[0]}/{source[1]}"][asset]["V"]
        assert observations["unusable"] == []
        assert finding.locators == (m.Locator(asset, ogrdb.EXPORT_ANCHOR),)


def test_igl_v_rss_spacing_names_the_two_base_intervals_and_the_missing_rss():
    _, inventory = interpreted(IGL)
    for finding in findings(inventory, "v_rss_spacing"):
        observations = finding.observations.get()
        assert observations["adjacent"] == 82
        assert sorted(observations["gaps"]["2"]) == sorted(NOTE_LABELS)
        assert len(observations["no_v_rs_start"]) == 1


def test_a_v_rss_start_inside_the_coding_interval_is_recorded_and_nothing_moves():
    _, inventory = interpreted(IGH)
    found = {f.observations.get()["label"]: f for f in findings(inventory, "v_rss_start_inside_coding_interval")}
    assert set(found) == {"IGHV1-58*01", "IGHV3-38*01"}
    (row,) = [r for r in SELECTION["selected_record_observations"] if r["label"] == "IGHV1-58*01"]
    observations = found["IGHV1-58*01"].observations.get()
    assert (observations["gene_start"], observations["gene_end"], observations["v_rs_start"],
            observations["record_bases_at_v_rs_start"], observations["record_bases_after_gene_end"]) == (
        row["gene_start"], row["gene_end"], row["v_rs_start"], row["sequence_at_v_rs_start"],
        row["sequence_after_gene_end"])
    assert observations["gap"] == -115
    # Interpretation evidence: IGHV3-38*01 declares v_rs_start 452 before gene_end 460.
    assert (found["IGHV3-38*01"].observations.get()["v_rs_start"],
            found["IGHV3-38*01"].observations.get()["gene_end"]) == (452, 460)
    for fmt in GAPPED:
        candidate = candidate_of(inventory, fmt, "IGHV1-58*01")
        assert_status(claim_of(inventory, candidate, "coding_interval"), "accepted", "source", "stated",
                      [row["gene_start"] - 1, row["gene_end"]])
        assert claim_of(inventory, candidate, "anchor").status == "accepted"
    for source in (IGK, IGL):
        assert findings(interpreted(source)[1], "v_rss_start_inside_coding_interval") == []
        assert findings(interpreted(source)[1], "v_rss_start_before_coding_interval") == []


def test_a_v_rss_start_before_the_coding_interval_has_its_own_finding():
    def edit(fmt, record):
        if record["label"] == "IGLV10-54*04":
            record["v_rs_start"] = 100
    lock, inventory = interpret(small(IGL, {"IGLV10-54*04", "IGLJ1*01"}, edit=edit))
    (finding,) = findings(inventory, "v_rss_start_before_coding_interval")
    assert (finding.observations.get()["v_rs_start"], finding.observations.get()["gap"]) == (100, 100 - 463 - 1)
    assert findings(inventory, "coding_end_before_declared_v_rss") == []


# ── historically affirmed spans ──────────────────────────────────────


def test_historical_affirmation_notes_stay_unresolved():
    _, inventory = interpreted(IGL)
    found = {f.observations.get()["label"]: f for f in findings(inventory, "affirmed_span_mapping_unresolved")}
    assert set(found) == set(NOTE_LABELS)
    notes = {n.label: n for n in profiles.PROFILE.historical_notes}
    for label, finding in found.items():
        observations = finding.observations.get()
        assert observations["quotes"] == list(notes[label].quotes)
        assert (observations["field"], observations["numbering"], observations["text_matches"],
                observations["sequence_bound"], observations["gapped_mapping"]) == (
            "curation", "imgt_gapped", True, True, "exact")
        assert "no reviewed interpretation" in observations["reason"]
        assert {s.asset for s in finding.subjects} == {"airr", "airr_ex"}
        for fmt in GAPPED:
            assert at(record_of(inventory, fmt, label), "curation") in finding.locators
    assert findings(inventory, "affirmed_span_mapped") == []
    for source in (IGH, IGK):
        assert findings(interpreted(source)[1], "affirmed_span_mapping_unresolved") == []


def test_a_reviewed_note_would_give_an_affirmed_span_through_the_exact_alignment(monkeypatch):
    reviewed = tuple(dataclasses.replace(n, interpretation=profiles.NoteInterpretation(339))
                     if n.label == "IGLV3-21*04" else n for n in profiles.PROFILE.historical_notes)
    profile = dataclasses.replace(profiles.PROFILE, historical_notes=reviewed)
    monkeypatch.setattr(profiles, "PROFILES", types.MappingProxyType({(profile.id, profile.version): profile}))
    lock, inventory = inventory_of(IGL)
    inventory = evidence.interpret(lock, pinned(IGL), inventory)
    (finding,) = findings(inventory, "affirmed_span_mapped")
    expected = {sum(1 for char in raw_gapped(IGL, fmt, "IGLV3-21*04")[:339] if char != ".")
                for fmt in ("gapped", "gapped_ex")}
    assert len(expected) == 1
    observations = finding.observations.get()
    assert observations["label"] == "IGLV3-21*04"
    assert observations["affirmed_interval"] == [0, expected.pop()]
    assert observations["last_affirmed_column"] == 339
    assert {f.observations.get()["label"] for f in findings(inventory, "affirmed_span_mapping_unresolved")} == \
        set(NOTE_LABELS) - {"IGLV3-21*04"}
    assert coding_end_findings(inventory)["IGLV3-21*04"].observations.get()["historical_note"][
        "affirmed_span"] == "mapped"


def test_a_note_whose_text_changed_is_not_read():
    def edit(fmt, record):
        if record["label"] == "IGLV3-21*04":
            record["curation"] = record["curation"] + " Edited."
    lock, inventory = interpret(small(IGL, {"IGLV3-21*04", "IGLJ1*01"}, edit=edit))
    (finding,) = findings(inventory, "affirmed_span_mapping_unresolved")
    observations = finding.observations.get()
    assert observations["text_matches"] is False and "not the note text" in observations["reason"]
    assert coding_end_findings(inventory)["IGLV3-21*04"].observations.get()["historical_note"] is None


def test_only_genomic_support_entries_are_listed_as_genomic_support():
    def edit(fmt, record):
        if record["label"] == "IGLV10-54*04":
            record["rearranged_support"] = [{"sequence_id": "r1", "repository_name": "VDJbase",
                                             "repository_ref": "P1"}]
    lock, inventory = interpret(small(IGL, {"IGLV10-54*04", "IGLJ1*01"}, edit=edit))
    (finding,) = findings(inventory, "coding_end_before_declared_v_rss")
    listed = finding.observations.get()["genomic_support"]
    raw_support = raw_record(IGL, "airr", "IGLV10-54*04")["unrearranged_support"]
    assert [(s["repository"], s["reference"]) for s in listed] == [
        (s["repository_name"], s["repository_ref"]) for s in raw_support]
    assert ("VDJbase", "P1") not in [(s["repository"], s["reference"]) for s in listed]
