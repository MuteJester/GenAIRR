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
import json
import types

import pytest

from GenAIRR._artifact_io import canonical_json
from GenAIRR._reference_sources import evidence
from GenAIRR._reference_sources import model as m
from GenAIRR._reference_sources import ogrdb
from GenAIRR._reference_sources import profiles
from tests.test_ogrdb_source_coordinates import (
    EVIDENCE, GAPPED, NOTE_LABELS, assert_status, at, candidate_of, claim_of, claims, findings, interpret,
    interpret_in_scope, interpreted, raw_gapped, raw_record, record_of, small, span_limiting_notes,
)
from tests.test_ogrdb_source_inventory import (
    EXPECT, IGH, IGK, IGL, SELECTION, SETS, _key, build, error, inventory_of, pinned, raw_records, seal,
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


@pytest.mark.parametrize("field_name, value", [
    ("j_phase_independence", "independent"), ("v_phase_derivation", "stated"),
    ("functionality_classes", types.MappingProxyType({"false": "pseudogene"})), ("landmark_columns", (309, 311)),
    ("missing_delineation_value", -1), ("rules", types.MappingProxyType({"j_anchor_coding0": "j_cdr3_end - gene_start"})),
    ("v_delineation_scheme", "Kabat"),
])
def test_a_profile_declaring_conventions_the_code_does_not_implement_is_refused(monkeypatch, field_name, value):
    profile = dataclasses.replace(profiles.PROFILE, **{field_name: value})
    monkeypatch.setattr(profiles, "PROFILES", types.MappingProxyType({(profile.id, profile.version): profile}))
    lock, inventory = inventory_of(IGL)
    found = error("unsupported_profile", evidence.interpret, lock, pinned(IGL), inventory)
    assert field_name in found.message


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
        assert_status(claim, "accepted", "assertion", "variant_comparison",
                      {"side": "3", "length": len(extension), "matches_record_context": status == "matches"})
        assert claim.coordinate_system is None
        short = candidate_of(inventory, "airr", label)
        assert long.sequence == short.sequence + extension
        assert short.locator in claim.evidence and long.locator in claim.evidence
        assert [x for x in inventory.claims if x.subject == short.locator and x.predicate == "extension"] == []


@pytest.mark.parametrize("source", SETS, ids=_key)
def test_a_matching_context_that_is_only_the_records_own_end_is_marked(source):
    _, inventory = interpreted(source)
    found = extension_findings(inventory)
    for label, extension in _suffixes(source).items():
        record = raw_record(source, "airr_ex", label)
        status, _ = _context_oracle(source, label, extension)
        expected = (None if status == "without_record_context"
                    else record["gene_end"] + len(extension) == len(record["sequence"]))
        assert found[label][1].observations.get()["context_is_sequence_terminus"] == expected, label


def test_only_ighv3_7_04_and_the_igl_matches_have_bases_beyond_the_extension():
    terminus = {}
    for source in SETS:
        _, inventory = interpreted(source)
        for label, (code, finding) in extension_findings(inventory).items():
            if code == "coding_extension_matches_record_context":
                terminus[label] = finding.observations.get()["context_is_sequence_terminus"]
    assert sorted(label for label, value in terminus.items() if not value) == sorted(
        ["IGHV3-7*04"] + list(NOTE_LABELS))
    assert sum(1 for label, value in terminus.items() if value and label.startswith("IGHV")) == 14
    assert terminus["IGKV1D-8*04"] is True


@pytest.mark.parametrize("prefix_length, terminus", [(28, True), (3, False)])
def test_a_5prime_extension_context_is_a_terminus_only_where_the_record_starts(monkeypatch, prefix_length,
                                                                                terminus):
    label = "IGLJ1*01"
    record = raw_record(IGL, "airr_ex", label)
    assert record["gene_start"] == 29                 # 28 record bases before the coding interval
    prefix = record["sequence"][28 - prefix_length:28]

    def edit(fmt, rec):
        if fmt == "airr_ex" and rec["label"] == label:
            rec["coding_sequence"] = prefix + rec["coding_sequence"]
    lock, inventory = interpret_in_scope(monkeypatch, small(IGL, {label}, edit=edit))
    ((code, finding),) = extension_findings(inventory).values()
    observations = finding.observations.get()
    assert code == "coding_extension_matches_record_context"
    assert (observations["side"], observations["extension"], observations["record_context"]) == ("5", prefix, prefix)
    assert observations["context_is_sequence_terminus"] is terminus


def test_a_partial_disagreeing_context_that_reaches_the_record_end_is_a_terminus():
    label = "IGHV4-4*09"                              # airr_ex adds GA; the record ends at gene_end

    def edit(fmt, record):
        if record["label"] == label:
            record["sequence"] = record["sequence"] + "T"
    lock, inventory = interpret(small(IGH, {label, "IGHJ1*01"}, edit=edit), IGH)
    ((code, finding),) = extension_findings(inventory).values()
    observations = finding.observations.get()
    assert (code, observations["extension"], observations["record_context"]) == (
        "coding_extension_disagrees_with_record_context", "GA", "T")
    assert observations["context_is_sequence_terminus"] is True


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
    (note,) = [n for n in profiles.PROFILE.historical_notes if n.label == "IGHV3-7*04"]
    assert observations["historical_note"] == {"field": "curation", "quotes": list(note.quotes),
                                               "affirmed_span": "unresolved"}
    assert note.stated_column == 319
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
def test_the_gap_between_coding_end_and_v_rss_start_is_counted_per_export(source):
    _, inventory = interpreted(source)
    summaries = {f.observations.get()["asset"]: f for f in findings(inventory, "v_rss_gap_after_coding_end")}
    assert set(summaries) == set(GAPPED)
    expected = EXPECT["sets"][f"{source[0]}/{source[1]}"]["airr_v_rss_gap_counts"]
    assert _gap_counts(summaries["airr"].observations.get()) == expected
    assert _gap_counts(summaries["airr_ex"].observations.get()) == _raw_gap_counts(source, "airr_ex")
    for asset, finding in summaries.items():
        observations = finding.observations.get()
        assert observations["v_records"] == EXPECT["sets"][f"{source[0]}/{source[1]}"][asset]["V"]
        assert observations["unusable"] == []
        assert finding.locators == (m.Locator(asset, ogrdb.EXPORT_ANCHOR),)


def test_igl_v_rss_gaps_name_the_two_base_intervals_and_the_missing_rss():
    _, inventory = interpreted(IGL)
    for finding in findings(inventory, "v_rss_gap_after_coding_end"):
        observations = finding.observations.get()
        assert observations["adjacent"] == 82
        assert sorted(observations["gaps"]["2"]) == sorted(NOTE_LABELS)
        assert len(observations["no_v_rs_start"]) == 1


def _raw_v_rss(source, check) -> dict:
    """label -> raw values of V records (both exports agree) for which ``check(record)`` holds."""
    found = {}
    for fmt in GAPPED:
        for record in raw_records(source, fmt):
            if record["sequence_type"] == "V" and record["v_rs_start"] is not None and check(record):
                found[record["label"]] = record
    return found


def _without_cac(record) -> bool:
    start = record["v_rs_start"]
    return start <= len(record["sequence"]) and record["sequence"][start - 1:start + 2] != "CAC"


@pytest.mark.parametrize("source, count", [(IGH, 9), (IGK, 2), (IGL, 4)], ids=["IGH", "IGK", "IGL"])
def test_a_declared_v_rss_start_without_a_cac_core_is_recorded(source, count):
    _, inventory = interpreted(source)
    found = {f.observations.get()["label"]: f.observations.get()
             for f in findings(inventory, "declared_v_rss_start_without_cac")}
    oracle = _raw_v_rss(source, _without_cac)
    assert set(found) == set(oracle) and len(found) == count
    for label, observations in found.items():
        record = oracle[label]
        start = record["v_rs_start"]
        assert (observations["v_rs_start"], observations["record_bases_at_v_rs_start"]) == (
            start, record["sequence"][start - 1:start + 6])
    for summary in findings(inventory, "v_rss_gap_after_coding_end"):
        assert summary.observations.get()["starts_without_cac"] == sorted(oracle)
    if source == IGH:
        for label in ("IGHV3-16*02", "IGHV3-16*03", "IGHV3-16*04"):
            assert found[label]["record_bases_at_v_rs_start"] == "TCCTGTG" and found[label]["gap"] == 0


def test_a_declared_v_rss_end_beyond_the_record_sequence_is_recorded():
    _, inventory = interpreted(IGH)
    found = {f.observations.get()["label"]: f.observations.get()
             for f in findings(inventory, "v_rss_coordinate_beyond_record_sequence")}
    oracle = _raw_v_rss(IGH, lambda r: r["v_rs_end"] is not None and r["v_rs_end"] > len(r["sequence"]))
    assert set(found) == set(oracle) and len(found) == 11
    assert {k: found["IGHV3-69-1*02"][k] for k in ("field", "value", "sequence_length")} == {
        "field": "v_rs_end", "value": 350, "sequence_length": 339}
    for source in (IGK, IGL):
        assert findings(interpreted(source)[1], "v_rss_coordinate_beyond_record_sequence") == []


def test_a_declared_v_rss_start_beyond_the_record_sequence_gives_no_position():
    def edit(fmt, record):
        if record["label"] == "IGLV10-54*04":
            record["v_rs_start"] = len(record["sequence"]) + 5
    lock, inventory = interpret(small(IGL, {"IGLV10-54*04", "IGLJ1*01"}, edit=edit))
    (finding,) = findings(inventory, "v_rss_coordinate_beyond_record_sequence")
    observations = finding.observations.get()
    assert (observations["field"], observations["value"], observations["sequence_length"]) == ("v_rs_start", 509, 504)
    for summary in findings(inventory, "v_rss_gap_after_coding_end"):
        assert summary.observations.get()["unusable"] == ["IGLV10-54*04"] and summary.observations.get()["gaps"] == {}
    assert findings(inventory, "coding_end_before_declared_v_rss") == []
    assert findings(inventory, "declared_v_rss_start_without_cac") == []


def test_a_declared_v_rss_start_too_close_to_the_record_end_is_not_judged_for_a_cac_core():
    def edit(fmt, record):
        if record["label"] == "IGLV10-54*04":
            record["v_rs_start"] = len(record["sequence"]) - 1          # two record bases from the start
    lock, inventory = interpret(small(IGL, {"IGLV10-54*04", "IGLJ1*01"}, edit=edit))
    assert findings(inventory, "declared_v_rss_start_without_cac") == []
    for summary in findings(inventory, "v_rss_gap_after_coding_end"):
        observations = summary.observations.get()
        assert observations["starts_without_cac"] == []
        assert observations["starts_too_close_to_sequence_end"] == ["IGLV10-54*04"]


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


def note_findings(inventory, code="affirmed_span_mapping_unresolved") -> dict:
    found = {}
    for finding in findings(inventory, code):
        label = finding.observations.get()["label"]
        assert label not in found
        found[label] = finding
    return found


def _span_oracle(source, gapped_fmt, label, column) -> dict:
    """The supplied span's last column against a stated column, read from the raw gapped record."""
    gapped = raw_gapped(source, gapped_fmt, label)
    last = len(gapped.rstrip("."))
    beyond = sum(1 for char in gapped[column:] if char != ".")
    relation = "exceeds" if last > column else "equals" if last == column else "short_of"
    return {"supplied_last_column": last, "span_vs_stated_column": relation, "bases_beyond_stated_column": beyond,
            "columns_short_of_stated_column": max(0, column - last)}


@pytest.mark.parametrize("source", SETS, ids=_key)
def test_every_span_limiting_note_is_reported_unresolved(source):
    _, inventory = interpreted(source)
    found = note_findings(inventory)
    oracle = span_limiting_notes(source)
    assert set(found) == set(oracle)
    assert findings(inventory, "affirmed_span_mapped") == []
    notes = {n.label: n for n in profiles.PROFILE.historical_notes if (n.set_id, n.release) == source}
    for label, finding in found.items():
        observations = finding.observations.get()
        column = min(oracle[label]) if oracle[label] else None
        assert (observations["field"], observations["numbering"], observations["stated_column"]) == (
            "curation", "imgt_gapped", column)
        assert observations["quotes"] == list(notes[label].quotes)
        assert finding.capabilities == ("original_coding_ends",)
        assert [r["asset"] for r in observations["records"]] == ["airr", "airr_ex"]
        for entry in observations["records"]:
            assert (entry["text_matches"], entry["sequence_bound"], entry["gapped_mapping"]) == (True, True, "exact")
            assert "no reviewed interpretation" in entry["reason"] and entry["affirmed_interval"] is None
            expected = (_span_oracle(source, GAPPED[entry["asset"]], label, column) if column is not None else
                        dict.fromkeys(("supplied_last_column", "span_vs_stated_column", "bases_beyond_stated_column",
                                       "columns_short_of_stated_column")))
            assert {k: entry[k] for k in expected} == expected
        assert {s.asset for s in finding.subjects} == {"airr", "airr_ex"}
        for fmt in GAPPED:
            assert at(record_of(inventory, fmt, label), "curation") in finding.locators


@pytest.mark.parametrize("source, label, beyond", [(IGL, "IGLV2-14*04", 2), (IGH, "IGHV2-70*20", 3)])
def test_supplied_coding_spans_that_run_past_the_affirmed_column(source, label, beyond):
    _, inventory = interpreted(source)
    records = note_findings(inventory)[label].observations.get()["records"]
    assert [(r["span_vs_stated_column"], r["bases_beyond_stated_column"]) for r in records] == [
        ("exceeds", beyond), ("exceeds", beyond)]


def test_a_note_without_a_stated_column_is_still_reported():
    _, inventory = interpreted(IGH)
    observations = note_findings(inventory)["IGHV1-3*04"].observations.get()
    assert observations["stated_column"] is None
    assert all(r["span_vs_stated_column"] is None for r in observations["records"])


def test_a_reviewed_note_would_give_an_affirmed_span_through_the_exact_alignment(monkeypatch):
    reviewed = tuple(dataclasses.replace(n, interpretation=profiles.NoteInterpretation(1, 339))
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
    assert (observations["first_affirmed_column"], observations["last_affirmed_column"]) == (1, 339)
    assert [r["affirmed_interval"] for r in observations["records"]] == [[0, expected.pop()]] * 2
    assert set(note_findings(inventory)) == set(span_limiting_notes(IGL)) - {"IGLV3-21*04"}
    assert coding_end_findings(inventory)["IGLV3-21*04"].observations.get()["historical_note"][
        "affirmed_span"] == "mapped"


SPAN_KEYS = ("supplied_last_column", "span_vs_stated_column", "bases_beyond_stated_column",
             "columns_short_of_stated_column")


def test_a_note_whose_text_changed_is_not_read():
    def edit(fmt, record):
        if record["label"] == "IGLV3-21*04":
            record["curation"] = record["curation"] + " Edited."
    lock, inventory = interpret(small(IGL, {"IGLV3-21*04", "IGLJ1*01"}, edit=edit))
    (finding,) = findings(inventory, "affirmed_span_mapping_unresolved")
    observations = finding.observations.get()
    # The profile's reading of the old text is not presented as the record's current note.
    assert (observations["text_matches"], observations["quotes"], observations["stated_column"]) == (False, [], None)
    assert observations["reviewed_text_sha256"] == [n for n in profiles.PROFILE.historical_notes
                                                     if n.label == "IGLV3-21*04"][0].text_sha256
    for entry in observations["records"]:
        assert entry["text_matches"] is False and "not the note text" in entry["reason"]
        assert "not applied" in entry["reason"]
        assert all(entry[k] is None for k in SPAN_KEYS)
    assert coding_end_findings(inventory)["IGLV3-21*04"].observations.get()["historical_note"] is None


def test_an_edited_note_stating_another_column_is_not_compared_with_the_profiles_column():
    def edit(fmt, record):
        if record["label"] == "IGLV2-14*04":
            record["curation"] = "Affirmed up to and including nucleotide 339 after review."
    lock, inventory = interpret(small(IGL, {"IGLV2-14*04", "IGLJ1*01"}, edit=edit))
    (finding,) = findings(inventory, "affirmed_span_mapping_unresolved")
    observations = finding.observations.get()
    assert observations["stated_column"] is None and observations["quotes"] == []
    assert all(entry[k] is None for entry in observations["records"] for k in SPAN_KEYS)


def test_a_note_read_against_another_sequence_is_not_compared_with_this_sequence():
    def edit(fmt, record):
        if record["label"] == "IGLV2-14*04":
            record["coding_sequence"] = record["coding_sequence"][:-1]
            record["gene_end"] -= 1

    def gapped_edit(fmt, header, sequence):
        return header, (sequence[:-1] if header == "IGLV2-14*04" else sequence)

    lock, inventory = interpret(small(IGL, {"IGLV2-14*04", "IGLJ1*01"}, edit=edit, gapped_edit=gapped_edit))
    (finding,) = findings(inventory, "affirmed_span_mapping_unresolved")
    observations = finding.observations.get()
    (note,) = [n for n in profiles.PROFILE.historical_notes if n.label == "IGLV2-14*04"]
    assert (observations["text_matches"], observations["stated_column"]) == (True, 337)   # the text is the note's
    assert observations["quotes"] == list(note.quotes)
    for entry in observations["records"]:
        assert (entry["text_matches"], entry["sequence_bound"], entry["gapped_mapping"]) == (True, False, "exact")
        assert "not the sequence" in entry["reason"] and "not applied" in entry["reason"]
        assert all(entry[k] is None for k in SPAN_KEYS)


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


# ── boundaries, functionality, lineage and support ───────────────────


BOUNDARY_FIELDS = {
    ("V", "boundary_5"): ("gene_start", "curational_tags"),
    ("V", "boundary_3"): ("gene_end", "v_rs_start", "v_rs_end", "curational_tags"),
    ("D", "boundary_5"): ("gene_start", "d_rs_5_prime_start", "d_rs_5_prime_end", "curational_tags"),
    ("D", "boundary_3"): ("gene_end", "d_rs_3_prime_start", "d_rs_3_prime_end", "curational_tags"),
    ("J", "boundary_5"): ("gene_start", "j_rs_start", "j_rs_end", "curational_tags"),
    ("J", "boundary_3"): ("gene_end", "j_donor_splice", "curational_tags"),
}


@pytest.mark.parametrize("source", SETS, ids=_key)
def test_no_coding_boundary_is_ever_claimed_complete(source):
    _, inventory = interpreted(source)
    parsed = [r for r in inventory.records if r.status == "parsed"]
    for record in parsed:
        (candidate,) = inventory.candidates_of(record.locator)
        for predicate in ("boundary_5", "boundary_3"):
            claim = claim_of(inventory, candidate, predicate)
            assert_status(claim, "unresolved", "none", "stated")
            assert claim.value.get() is None and claim.mappings == ()
            assert set(claim.evidence) == {at(record, name) for name in BOUNDARY_FIELDS[(record.segment, predicate)]}
    assert len(claims(inventory, "boundary_5")) == len(claims(inventory, "boundary_3")) == len(parsed)


def test_tags_functionality_rss_adjacency_and_exact_slices_set_no_boundary():
    _, inventory = interpreted(IGL)
    for label in NOTE_LABELS:       # likely_full_length, functional, adjacent context, exact airr slice
        raw_rec = raw_record(IGL, "airr", label)
        assert raw_rec["curational_tags"] == ["likely_full_length"] and raw_rec["functional"] is True
        candidate = candidate_of(inventory, "airr", label)
        assert claim_of(inventory, candidate, "coding_interval").status == "accepted"
        assert claim_of(inventory, candidate, "boundary_3").status == "unresolved"


@pytest.mark.parametrize("source", SETS, ids=_key)
def test_provider_functionality_stays_raw_and_unclassified(source):
    _, inventory = interpreted(source)
    for record in [r for r in inventory.records if r.status == "parsed"]:
        claim = claim_of(inventory, record, "functionality")
        assert_status(claim, "unresolved", "none", "stated")
        assert claim.evidence == (at(record, "functional"),)
    assert not [x for x in inventory.claims if x.predicate == "functionality" and x.value.get() is not None]
    summaries = {f.observations.get()["asset"]: f.observations.get()
                 for f in findings(inventory, "provider_functionality_unclassified")}
    assert set(summaries) == set(GAPPED)
    for fmt, observations in summaries.items():
        values = [r.get("functional", "absent") for r in raw_records(source, fmt)]
        assert observations["values"] == {"true": values.count(True), "false": values.count(False),
                                          "null": values.count(None), "absent": values.count("absent"),
                                          "other": len(values) - sum(values.count(v) for v in (True, False, None,
                                                                                               "absent"))}
        assert observations["classified"] == 0


def test_the_meaning_of_the_functional_flag_is_recorded_as_unresolved():
    detail = {u["id"]: u["detail"] for u in EVIDENCE["unresolved"]}["functional_false_meaning"]
    # The producer sets functional = 'F' in functionality, and 'F' also occurs in 'ORF'.
    assert "'F' in functionality" in detail and "'ORF'" in detail
    assert "true conflates" in detail and "false means only" in detail
    _, inventory = interpreted(IGL)
    for finding in findings(inventory, "provider_functionality_unclassified"):
        assert finding.observations.get()["reason"] == detail


def test_a_false_provider_functional_is_not_a_pseudogene():
    _, inventory = interpreted(IGK)
    record = record_of(inventory, "airr", "IGKV2-29*01")
    assert raw_record(IGK, "airr", "IGKV2-29*01")["functional"] is False
    assert inventory.observation(at(record, "functional")).value.get() is False
    assert claim_of(inventory, record, "functionality").status == "unresolved"


@pytest.mark.parametrize("source", SETS, ids=_key)
def test_all_exports_of_one_provider_set_are_one_evidence_lineage(source):
    lock, inventory = interpreted(source)
    lineage = f"ogrdb:{source[0]}/{source[1]}"
    (finding,) = findings(inventory, "source_evidence_lineage")
    observations = finding.observations.get()
    assert observations["lineage"] == lineage
    assert observations["assets"] == ["airr", "airr_ex", "gapped", "gapped_ex"]
    assert observations["independent_lineages"] == 1
    assert observations["support_status"] == "embedded_unretrieved"
    assert observations["support_references"] == len(inventory.support)
    for claim in inventory.claims:
        roots = list(claim.evidence) + [inventory.mapping(x).output for x in claim.mappings]
        assert {evidence.evidence_lineage(lock, root) for root in roots} == {lineage}
    anchors = claims(inventory, "anchor", "airr", "V")
    assert anchors and all({root.asset for root in a.evidence} == {"airr", "gapped"} for a in anchors)
    for code in ("coding_extension_matches_record_context", "coding_end_before_declared_v_rss"):
        for item in findings(inventory, code):
            assert item.observations.get()["lineage"] == lineage
    error("identity_mismatch", evidence.evidence_lineage, interpreted(IGL, ("airr", "airr_ex"))[0],
          m.Locator("gapped", "/0"))


def test_embedded_support_is_never_retrieved_or_upgraded():
    _, inventory = interpreted(IGL)
    assert {s.status for s in inventory.support} == {"embedded_unretrieved"}
    record = record_of(inventory, "airr", "IGLV6-57*04")
    found = {(s.repository, s.reference) for s in inventory.support if s.record == record.locator}
    assert ("GenBank", "Gibson:HG01106") in found              # a repository label with an opaque study id
    support = {s.locator for s in inventory.support}
    assert not [x for x in inventory.claims if set(x.evidence) & support]


@pytest.mark.parametrize("source", SETS, ids=_key)
def test_support_derivation_and_status_are_honest_for_every_claim(source):
    _, inventory = interpreted(source)
    rules = {   # predicate -> {(segment, status, support, derivation)}
        "anchor": {("V", "accepted", "source", "stated"), ("J", "accepted", "source", "stated")},
        "coding_phase": {("V", "accepted", "assertion", "anchor_derived"), ("J", "accepted", "source", "stated")},
        "phase_independence": {("V", "accepted", "assertion", "anchor_derived"), ("J", "unresolved", "none", "stated")},
        "coding_interval": {(s, "accepted", "source", "stated") for s in "VDJ"} | {
            (s, "unresolved", "none", "stated") for s in "VDJ"},
        "extension": {(s, "accepted", "assertion", "variant_comparison") for s in "VDJ"},
        "boundary_5": {(s, "unresolved", "none", "stated") for s in "VDJ"},
        "boundary_3": {(s, "unresolved", "none", "stated") for s in "VDJ"},
        "functionality": {(s, "unresolved", "none", "stated") for s in "VDJ"},
    }
    seen = set()
    for claim in inventory.claims:
        segment = inventory.subject_record(claim).segment
        key = (segment, claim.status, claim.support, claim.derivation)
        assert key in rules[claim.predicate], (claim.predicate, key)
        seen.add(claim.predicate)
        if claim.predicate == "phase_independence":
            assert claim.value.get() != "independent"
    assert seen == set(rules)
    assert not [x for x in inventory.claims if inventory.subject_record(x).segment == "D"
                and x.predicate in ("anchor", "coding_phase", "phase_independence")]


def test_outside_the_validated_scope_no_claim_is_accepted():
    blobs = pinned(IGL)
    blobs["airr_ex"] = blobs["airr_ex"] + b"\n"           # drifted bytes of a validated release
    lock, inventory = interpret(blobs)
    assert {x.status for x in inventory.claims} == {"candidate", "unresolved"}
    extension = claim_of(inventory, candidate_of(inventory, "airr_ex", "IGLV10-54*04"), "extension")
    assert_status(extension, "candidate", "none", "variant_comparison")


def test_identical_dna_identities_cite_only_their_own_records():
    _, inventory = interpreted(IGL)
    j2, j3 = record_of(inventory, "airr", "IGLJ2*01"), record_of(inventory, "airr", "IGLJ3*01")
    assert candidate_of(inventory, "airr", "IGLJ2*01").sha256 == candidate_of(inventory, "airr", "IGLJ3*01").sha256
    for record, other in ((j2, j3), (j3, j2)):
        mine = [x for x in inventory.claims if inventory.subject_record(x) == record]
        assert mine
        for claim in mine:
            assert all(root.is_under(record.locator) for root in claim.evidence)
            assert not any(root.is_under(other.locator) for root in claim.evidence)
    assert [r for r in inventory.records if r.locator.asset == "airr_ex" and r.label == "IGLJ3*01"] == []
    assert not [x for x in inventory.claims if x.subject.asset == "airr_ex"
                and inventory.subject_record(x).label == "IGLJ3*01"]


def test_interpret_validates_its_result_against_the_lock(monkeypatch):
    calls = []
    original = evidence.consistency.check_inventory

    def recording(inventory, lock):
        calls.append((inventory, lock))
        return original(inventory, lock)

    monkeypatch.setattr(evidence.consistency, "check_inventory", recording)
    lock, inventory = build(IGL, small(IGL, {"IGLJ1*01"}))
    result = evidence.interpret(lock, small(IGL, {"IGLJ1*01"}), inventory)
    assert calls == [(result, lock)]


def test_profile_notes_apply_only_to_their_own_set_and_release():
    blobs = small(IGL, {"IGLV3-21*04", "IGLJ1*01"})
    for fmt in GAPPED:
        document = json.loads(blobs[fmt].decode("utf-8"))
        germline_set = document["GermlineSet"][0]
        germline_set["release_version"] = 5.0
        germline_set["germline_set_ref"] = "OGRDB:IGLambda_VJ.5"
        blobs[fmt] = json.dumps(document).encode("utf-8")
    lock = seal(IGL, blobs, release=5)
    inventory = evidence.interpret(lock, blobs, ogrdb.build_inventory(lock, blobs))
    assert raw_record(IGL, "airr", "IGLV3-21*04")["curation"] in [
        o.value.get() for o in inventory.observations if o.locator.pointer.endswith("/curation")]
    assert findings(inventory, "affirmed_span_mapping_unresolved") == []
    (finding,) = findings(inventory, "coding_end_before_declared_v_rss")
    assert finding.observations.get()["historical_note"] is None
