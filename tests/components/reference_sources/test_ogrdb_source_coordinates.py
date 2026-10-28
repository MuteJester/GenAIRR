"""OGRDB coordinate profile: J and V landmarks, coding intervals and gapped mappings.

The interpretation profile ``ogrdb_api_v2_human_ig_observed_v1`` reads the pinned OGRDB
producer's conventions:

* J anchor = ``j_cdr3_end - 1`` in the coding sequence; reported J phase value =
  ``j_codon_frame - 1``.  The phase value is source-stated; whether the producer derived
  it independently of the landmark is unresolved.
* V anchor = ``cdr3_start - 4`` in the coding sequence, accepted only when the record's own
  label-identical gapped record of the same export variant ungaps to its coding sequence,
  carries a gap-free codon at columns 310-312 and gives the same ``cdr3_start``.  A V phase
  is at most anchor-derived.
* The coding interval ``sequence[gene_start - 1 : gene_end]`` is applied only where it
  reproduces the supplied coding sequence.

Expected values come from the independently authored fixture files
(``expectations.json``, ``interpretation_evidence.json``,
``human_igh_release10_candidate_selection.json``), from hand-checked literals, or from
reading raw payload fields with the standard library; never from the code under test.
The AIRR-prose J reading exists only in this file, as a test oracle that shows the
difference; it is not a production path.
"""
from __future__ import annotations

import dataclasses
import hashlib
import json

import pytest

from GenAIRR._reference_sources import coordinates as c
from GenAIRR._reference_sources import evidence
from GenAIRR._reference_sources import model as m
from GenAIRR._reference_sources import profiles
from tests.support.fixtures.ogrdb_inventory import (
    EXPECT, IGH, IGK, IGL, SETS, _key, encode, error, pinned, raw_records,
)
from tests.support.fixtures.ogrdb_coordinates import EVIDENCE, GAPPED, V_LABEL, assert_status, at, candidate_of, claim_of, claims, findings, interpret, interpret_in_scope, interpreted, raw_gapped, raw_record, record_of, small, span_limiting_notes

pytestmark = pytest.mark.behavior

PROFILE_REF = m.CoordinateProfileRef("ogrdb_api_v2_human_ig_observed_v1", 1)
SET_KEYS = {"9606.IGH_VDJ/10": IGH, "9606.IGKappa_VJ/5": IGK, "9606.IGLambda_VJ/4": IGL}
# Test oracle for span-limiting curation notes (never used by the code under test): an affirmed
# column ("up to and including nucleotide N" / "... base N") or a trailing-dot marker.


# ── helpers ──────────────────────────────────────────────────────────








def imgt_cdr3_start(record) -> object:
    (delineation,) = [d for d in record["v_gene_delineations"] if d["delineation_scheme"] == "IMGT"]
    return delineation["cdr3_start"]


def sha(text: str) -> str:
    return "sha256:" + hashlib.sha256(text.encode("utf-8")).hexdigest()


def airr_prose_j(record) -> tuple:
    """The AIRR 2.0 prose reading of the J fields (TEST ORACLE ONLY, never a production path):
    ``j_cdr3_end`` as a coordinate in the record's sequence field and ``j_codon_frame`` as the
    codon position of the first coding nucleotide."""
    return record["j_cdr3_end"] - record["gene_start"], (1 - record["j_codon_frame"]) % 3


# ── the interpretation profile ───────────────────────────────────────


def test_exactly_one_interpretation_profile_is_implemented():
    assert set(profiles.PROFILES) == {("ogrdb_api_v2_human_ig_observed_v1", 1)}
    profile = profiles.profile_for(PROFILE_REF)
    assert profile is profiles.PROFILE
    assert profile.ref == PROFILE_REF


def test_profile_lookup_is_explicit_and_never_guessed():
    error("unsupported_profile", profiles.profile_for, m.CoordinateProfileRef("ogrdb_api_v2_human_ig_observed_v2", 2))
    with pytest.raises(TypeError):
        profiles.profile_for("ogrdb_api_v2_human_ig_observed_v1")
    with pytest.raises(TypeError):
        profiles.profile_for(None)


def test_profile_scope_is_the_frozen_validated_registry_entry():
    assert profiles.PROFILE.scope == m.PROFILE_SCOPES[("ogrdb_api_v2_human_ig_observed_v1", 1)]


def test_profile_rules_are_the_producer_conventions_and_not_the_airr_prose():
    rules = EXPECT["rules_under_test"]
    producer = ("j_anchor_coding0", "j_reported_phase", "v_anchor_coding0", "v_cdr3_start_from_gapped",
                "coding_interval")
    assert profiles.PROFILE.rules == {name: rules[name] for name in producer}
    assert rules["j_airr_prose_anchor_coding0"] not in profiles.PROFILE.rules.values()
    assert rules["j_airr_prose_phase"] not in profiles.PROFILE.rules.values()


def test_profile_cites_the_pinned_producer_code_and_documents_as_data():
    producer = EVIDENCE["producer_code"]
    evidence = profiles.PROFILE.evidence
    assert (evidence.repository, evidence.commit) == (producer["repository"], producer["commit"])
    assert {(f.path, f.git_blob, f.sha256, tuple(sorted(r for r, _ in f.lines))) for f in evidence.files} == {
        (f["path"], f["git_blob"], "sha256:" + f["sha256"], tuple(sorted(f["lines"]))) for f in producer["files"]}
    tool = producer["related_tool"]
    related = evidence.related_tool
    assert (related.repository, related.commit, related.path, related.git_blob, related.sha256,
            tuple(sorted(r for r, _ in related.lines))) == (
        tool["repository"], tool["commit"], tool["path"], tool["git_blob"], "sha256:" + tool["sha256"],
        tuple(sorted(tool["lines"])))
    assert {(d.id, d.url, d.sha256, d.bytes) for d in evidence.documents} == {
        (d["id"], d.get("pinned_url", d["url"]), "sha256:" + d["sha256"], d["bytes"]) for d in EVIDENCE["documents"]}


def test_every_rule_names_the_pinned_lines_it_rests_on():
    cited = {(f["path"], lines) for f in EVIDENCE["producer_code"]["files"] for lines in f["lines"]}
    tool = EVIDENCE["producer_code"]["related_tool"]
    cited |= {(tool["path"], lines) for lines in tool["lines"]}
    basis = profiles.PROFILE.basis
    assert set(basis) == set(profiles.PROFILE.rules) | {"gapped_alignment", "missing_delineation_value"}
    for rule, citations in basis.items():
        assert citations, rule
        assert set(citations) <= cited, rule


def test_profile_keeps_the_unresolved_conventions_unresolved():
    assert profiles.PROFILE.unresolved == tuple((u["id"], u["detail"]) for u in EVIDENCE["unresolved"])
    assert profiles.PROFILE.j_phase_independence == "unresolved"
    assert profiles.PROFILE.v_phase_derivation == "anchor_derived"
    assert profiles.PROFILE.functionality_classes == {}      # no provider value is classified


def test_the_profile_declares_the_constants_the_rules_use():
    assert profiles.PROFILE.rules == c.RULES
    assert profiles.PROFILE.landmark_columns == (c.LANDMARK_FIRST_COLUMN, c.LANDMARK_LAST_COLUMN) == (310, 312)
    assert profiles.PROFILE.missing_delineation_value == c.MISSING_DELINEATION_VALUE == 0
    for name, implemented in evidence.IMPLEMENTED_CONVENTIONS.items():
        assert getattr(profiles.PROFILE, name) == implemented, name


def test_profile_landmark_residues():
    assert profiles.PROFILE.landmark_residues == {"V": {"IGH": ("C",), "IGK": ("C",), "IGL": ("C",)},
                                                  "J": {"IGH": ("W",), "IGK": ("F",), "IGL": ("F",)}}




def test_the_note_oracle_finds_the_settled_set_of_span_limiting_notes():
    counts = {key: len(span_limiting_notes(source, fmt)) for source, key in ((IGH, "IGH"), (IGK, "IGK"), (IGL, "IGL"))
              for fmt in ("airr", "airr_ex")}
    assert counts == {"IGH": 34, "IGK": 1, "IGL": 6}
    for source in (IGH, IGK, IGL):
        assert span_limiting_notes(source, "airr") == span_limiting_notes(source, "airr_ex")
    assert [label for label, columns in span_limiting_notes(IGH).items() if not columns] == ["IGHV1-3*04"]


def test_every_span_limiting_note_is_a_profile_note_bound_to_its_text_and_sequences():
    notes = profiles.PROFILE.historical_notes
    expected = {(source[0], source[1], label): columns for source in (IGH, IGK, IGL)
                for label, columns in span_limiting_notes(source).items()}
    assert {(n.set_id, n.release, n.label) for n in notes} == set(expected)
    assert len(notes) == len(expected) == 41
    for note in notes:
        source = (note.set_id, note.release)
        assert (note.field, note.numbering) == ("curation", "imgt_gapped")
        assert note.interpretation is None           # no reviewed interpretation exists
        texts = {raw_record(source, fmt, note.label)["curation"] for fmt in ("airr", "airr_ex")}
        (text,) = texts
        assert note.text_sha256 == sha(text)
        assert note.quotes and all(quote in text for quote in note.quotes)
        columns = expected[(note.set_id, note.release, note.label)]
        assert note.stated_column == (min(columns) if columns else None)
        if note.stated_column is not None:
            assert any(f"including nucleotide {note.stated_column}" in q or f"including base {note.stated_column}" in q
                       for q in note.quotes)
        assert dict(note.sequences) == {fmt: sha(raw_record(source, fmt, note.label)["coding_sequence"])
                                        for fmt in ("airr", "airr_ex")}


def test_note_quotes_keep_the_providers_own_confirmations():
    notes = {n.label: n for n in profiles.PROFILE.historical_notes}
    assert any("MK587524" in q and "confirms the inference" in q for q in notes["IGLV2-14*03"].quotes)
    assert any("bases 338-339" in q and "are T and C" in q for q in notes["IGLV2-14*03"].quotes)
    assert any("identical genomic sequence has been accepted by IMGT as the full length sequence" in q
               for q in notes["IGLV3-25*03"].quotes)
    assert any("final three nucleotides of the sequence - TAC" in q for q in notes["IGHV2-70*20"].quotes)
    for label in ("IGHV3-13*06", "IGHV3-43D*05"):
        assert any("missing terminal nucleotide of the sequence can be confirmed" in q for q in notes[label].quotes)


# ── stated positions ─────────────────────────────────────────────────


@pytest.mark.parametrize("value, problem", [
    (None, "null"), (c.ABSENT, "absent"), (True, "Boolean"), (False, "Boolean"), (8.0, "not an integer"),
    ("8", "not an integer"), ([8], "not an integer"), (-1, "below 1"),
])
def test_stated_positions_are_exact_integers(value, problem):
    found = c.stated_position(value)
    assert found.value is None and problem in found.problem


def test_zero_is_a_missing_delineation_coordinate_only_where_the_profile_says_so():
    assert c.stated_position(0, missing_value=0).value is None
    assert "missing" in c.stated_position(0, missing_value=0).problem
    assert c.stated_position(0).value is None and "below 1" in c.stated_position(0).problem
    assert c.stated_position(8) == c.Position(8, None)


# ── J landmark and phase ─────────────────────────────────────────────


def test_j_examples_cover_all_three_phase_classes_and_nonzero_context_offsets():
    assert {e["j_codon_frame"] for e in EXPECT["j_examples"]} == {1, 2, 3}
    assert all(e["gene_start"] > 1 for e in EXPECT["j_examples"])


@pytest.mark.parametrize("example", EXPECT["j_examples"], ids=lambda e: e["label"])
def test_j_landmark_follows_the_producer_convention(example):
    record = raw_record(SET_KEYS[example["set"]], example["export"], example["label"])
    assert (len(record["sequence"]), len(record["coding_sequence"]), record["gene_start"], record["gene_end"],
            record["j_cdr3_end"], record["j_codon_frame"]) == (
        example["sequence_length"], example["coding_length"], example["gene_start"], example["gene_end"],
        example["j_cdr3_end"], example["j_codon_frame"])
    found = c.j_landmark(record["coding_sequence"], record["j_cdr3_end"], record["j_codon_frame"])
    assert (found.anchor, found.codon, found.phase) == (example["anchor_coding0"], example["anchor_codon"],
                                                       example["phase"])
    assert found.anchor_problem is None and found.phase_problem is None and found.consistent is True


def test_hand_checked_j_literals():
    j7 = raw_record(IGL, "airr_ex", "IGLJ7*02")
    found = c.j_landmark(j7["coding_sequence"], j7["j_cdr3_end"], j7["j_codon_frame"])
    assert (found.anchor, found.codon, found.phase) == (7, "TTC", 1)
    j4 = raw_record(IGH, "airr", "IGHJ4*02")
    found = c.j_landmark(j4["coding_sequence"], j4["j_cdr3_end"], j4["j_codon_frame"])
    assert (found.anchor, found.codon, found.phase) == (14, "TGG", 2)


@pytest.mark.parametrize("example", EXPECT["j_examples"], ids=lambda e: e["label"])
def test_the_airr_prose_reading_is_a_different_convention(example):
    """The oracle puts the landmark before the coding sequence, and the producer's number
    read in the record context names another codon; production never yields either."""
    record = raw_record(SET_KEYS[example["set"]], example["export"], example["label"])
    anchor, phase = airr_prose_j(record)
    assert anchor == example["airr_prose_anchor_coding0"] < 0
    offset = record["j_cdr3_end"] - 1
    assert record["sequence"][offset:offset + 3] == example["same_offset_in_sequence"] != example["anchor_codon"]
    found = c.j_landmark(record["coding_sequence"], record["j_cdr3_end"], record["j_codon_frame"])
    assert found.anchor != anchor
    if record["j_codon_frame"] != 1:
        assert found.phase != phase


def test_the_airr_prose_reading_is_never_a_fallback():
    j7 = raw_record(IGL, "airr_ex", "IGLJ7*02")
    # j_cdr3_end 40 puts the producer landmark at 39, beyond the 38 coding bases.  The prose
    # reading would give 40 - 29 = 11, a valid index: it must not be tried.
    assert airr_prose_j(dict(j7, j_cdr3_end=40))[0] == 11
    found = c.j_landmark(j7["coding_sequence"], 40, j7["j_codon_frame"])
    assert found.anchor is None and found.codon is None and "beyond" in found.anchor_problem


@pytest.mark.parametrize("value, problem", [
    (None, "null"), (c.ABSENT, "absent"), (True, "Boolean"), (8.0, "not an integer"), ("8", "not an integer"),
    (0, "below 1"), (-3, "below 1"), (37, "beyond"), (10 ** 12, "beyond"),
])
def test_j_cdr3_end_values_that_give_no_anchor(value, problem):
    j7 = raw_record(IGL, "airr_ex", "IGLJ7*02")
    found = c.j_landmark(j7["coding_sequence"], value, 2)
    assert found.anchor is None and problem in found.anchor_problem
    assert found.phase == 1 and found.consistent is None


@pytest.mark.parametrize("value, problem", [
    (None, "null"), (c.ABSENT, "absent"), (True, "Boolean"), (2.0, "not an integer"), ("2", "not an integer"),
    (0, "1, 2 or 3"), (4, "1, 2 or 3"), (-1, "1, 2 or 3"),
])
def test_j_codon_frame_values_that_give_no_phase(value, problem):
    j7 = raw_record(IGL, "airr_ex", "IGLJ7*02")
    found = c.j_landmark(j7["coding_sequence"], 8, value)
    assert found.phase is None and problem in found.phase_problem
    assert found.anchor == 7 and found.consistent is None


def test_a_stated_phase_that_disagrees_with_the_landmark_is_not_shifted():
    j7 = raw_record(IGL, "airr_ex", "IGLJ7*02")
    found = c.j_landmark(j7["coding_sequence"], 8, 3)
    assert (found.anchor, found.phase, found.consistent) == (7, 2, False)


# ── V landmark ───────────────────────────────────────────────────────


@pytest.mark.parametrize("example", EXPECT["v_examples"], ids=lambda e: f"{e['label']}-{e['export']}")
def test_v_landmark_is_bound_to_its_own_gapped_record(example):
    source = SET_KEYS[example["set"]]
    record = raw_record(source, example["export"], example["label"])
    gapped = raw_gapped(source, example["gapped_export"], example["label"])
    assert (len(record["coding_sequence"]), imgt_cdr3_start(record), record["functional"]) == (
        example["coding_length"], example["cdr3_start"], example["functional"])
    assert gapped[309:312] == example["gapped_columns_310_312"]
    assert c.gapped_cdr3_start(gapped) == example["non_gap_in_columns_1_312"] + 1
    assert c.gapped_to_coding(gapped, 309, 312) == (example["non_gap_in_columns_1_309"],
                                                    example["non_gap_in_columns_1_312"])
    found = c.v_landmark(record["coding_sequence"], imgt_cdr3_start(record), gapped)
    assert (found.anchor, found.codon, found.problem, found.kind) == (example["anchor_coding0"],
                                                                     example["anchor_codon"], None, None)
    offset = example["anchor_coding0"]
    assert record["sequence"][offset:offset + 3] == example["same_offset_in_sequence"]


@pytest.mark.parametrize("label, anchor, codon", [
    ("IGHV4-34*01", 282, "TGT"), ("IGHV5-51*01", 285, "TGT"),
])
def test_hand_checked_v_literals(label, anchor, codon):
    record = raw_record(IGH, "airr", label)
    found = c.v_landmark(record["coding_sequence"], imgt_cdr3_start(record), raw_gapped(IGH, "gapped", label))
    assert (found.anchor, found.codon, found.problem) == (anchor, codon, None)


def test_a_noncanonical_v_landmark_is_never_moved_to_a_nearby_cysteine():
    record = raw_record(IGL, "airr_ex", "IGLV3-2*03")
    coding = record["coding_sequence"]
    found = c.v_landmark(coding, imgt_cdr3_start(record), raw_gapped(IGL, "gapped_ex", "IGLV3-2*03"))
    assert (found.anchor, found.codon) == (258, "GCT")
    # Out-of-frame cysteine codons nearby are never used (hand-checked positions).
    assert coding[254:257] == "TGT" and coding[272:275] == "TGC"
    assert found.anchor not in (254, 272)


@pytest.mark.parametrize("value, problem", [
    (0, "missing"), (None, "null"), (c.ABSENT, "absent"), (True, "Boolean"), (289.0, "not an integer"),
    (3, "below 4"), (298, "beyond"),
])
def test_cdr3_start_values_that_give_no_v_anchor(value, problem):
    record = raw_record(IGH, "airr", "IGHV1-2*02")
    found = c.v_landmark(record["coding_sequence"], value, raw_gapped(IGH, "gapped", "IGHV1-2*02"))
    assert found.anchor is None and found.kind == "field" and problem in found.problem


def _v_case(label="IGHV1-2*02"):
    record = raw_record(IGH, "airr", label)
    return record["coding_sequence"], imgt_cdr3_start(record), raw_gapped(IGH, "gapped", label)


def test_a_v_anchor_needs_a_gapped_record():
    coding, cdr3_start, _ = _v_case()
    found = c.v_landmark(coding, cdr3_start, None)
    assert found.anchor is None and found.kind == "alignment"


def test_a_cdr3_start_that_disagrees_with_the_gapped_record_is_not_bound():
    coding, cdr3_start, gapped = _v_case()
    found = c.v_landmark(coding, cdr3_start + 3, gapped)
    assert found.anchor is None and found.kind == "binding" and "289" in found.problem
    assert found.gapped_cdr3_start == 289


def test_a_gapped_record_that_does_not_ungap_to_the_coding_sequence_is_not_bound():
    coding, cdr3_start, gapped = _v_case()
    position = gapped.index("A")
    changed = gapped[:position] + "C" + gapped[position + 1:]
    found = c.v_landmark(coding, cdr3_start, changed)
    assert found.anchor is None and found.kind == "binding" and "ungap" in found.problem


def test_a_gap_in_the_landmark_codon_columns_is_not_bound():
    coding, cdr3_start, gapped = _v_case()
    shifted = gapped[:310] + "." + gapped[310:]          # a gap at column 311
    assert shifted.replace(".", "") == coding
    found = c.v_landmark(coding, cdr3_start, shifted)
    assert found.anchor is None and found.kind == "binding" and "310-312" in found.problem
    assert c.gapped_to_coding(shifted, 309, 312) is None


def test_a_gapped_record_shorter_than_312_columns_is_not_bound():
    gapped = raw_gapped(IGH, "gapped", "IGHV1-2*02")[:300]
    found = c.v_landmark(gapped.replace(".", ""), 100, gapped)
    assert found.anchor is None and found.kind == "binding" and "312" in found.problem
    assert c.gapped_cdr3_start(gapped) is None


# ── coding interval ──────────────────────────────────────────────────


@pytest.mark.parametrize("example", EXPECT["j_examples"], ids=lambda e: e["label"])
def test_the_coding_interval_is_applied_where_it_reproduces_the_coding_sequence(example):
    record = raw_record(SET_KEYS[example["set"]], example["export"], example["label"])
    found = c.coding_interval(record["sequence"], record["gene_start"], record["gene_end"], record["coding_sequence"])
    assert (found.start, found.end, found.problem, found.reproduces) == (
        example["gene_start"] - 1, example["gene_end"], None, True)


@pytest.mark.parametrize("row", EXPECT["igl4_airr_ex_coding_end_observations"], ids=lambda r: r["label"])
def test_an_extended_coding_sequence_is_not_reproduced_by_its_declared_interval(row):
    record = raw_record(IGL, "airr_ex", row["label"])
    found = c.coding_interval(record["sequence"], record["gene_start"], record["gene_end"], record["coding_sequence"])
    assert (found.start, found.end, found.reproduces) == (record["gene_start"] - 1, row["gene_end"], False)
    assert (found.relation, found.extra) == ("candidate_extends_interval_3prime", row["coding_extension"])
    unextended = raw_record(IGL, "airr", row["label"])
    found = c.coding_interval(unextended["sequence"], unextended["gene_start"], unextended["gene_end"],
                              unextended["coding_sequence"])
    assert found.reproduces is True


@pytest.mark.parametrize("gene_start, gene_end, problem", [
    (None, 66, "gene_start"), (True, 66, "gene_start"), (0, 66, "gene_start"), (29.0, 66, "gene_start"),
    (29, None, "gene_end"), (29, False, "gene_end"), (40, 30, "after gene_end"), (29, 67, "beyond"),
])
def test_coding_interval_values_that_give_no_interval(gene_start, gene_end, problem):
    j7 = raw_record(IGL, "airr_ex", "IGLJ7*02")
    found = c.coding_interval(j7["sequence"], gene_start, gene_end, j7["coding_sequence"])
    assert found.start is None and found.reproduces is None and problem in found.problem


@pytest.mark.parametrize("sequence", [None, c.ABSENT, 5, ["ACGT"]])
def test_a_record_without_a_sequence_text_gives_no_interval(sequence):
    found = c.coding_interval(sequence, 1, 3, "ACG")
    assert found.start is None and "sequence" in found.problem


# ── record context of an extension ───────────────────────────────────


@pytest.mark.parametrize("row", EXPECT["igl4_airr_ex_coding_end_observations"], ids=lambda r: r["label"])
def test_extension_context_is_read_at_the_declared_end(row):
    record = raw_record(IGL, "airr_ex", row["label"])
    found = c.extension_context(record["sequence"], (record["gene_start"] - 1, record["gene_end"]), "3",
                                row["coding_extension"])
    assert found.bases == row["sequence_after_gene_end"]
    expected = "matches" if row["finding"] == "coding_end_before_declared_v_rss" else "disagrees"
    assert found.status == expected


def test_extension_context_without_or_with_partial_record_bases():
    assert c.extension_context("AAACCC", (0, 6), "3", "GA") == c.ExtensionContext("without_record_context", "")
    assert c.extension_context("AAACCCG", (0, 6), "3", "GA") == c.ExtensionContext("without_record_context", "G")
    assert c.extension_context("AAACCCT", (0, 6), "3", "GA") == c.ExtensionContext("disagrees", "T")
    assert c.extension_context("TTAAACCC", (2, 8), "5", "TT") == c.ExtensionContext("matches", "TT")
    assert c.extension_context("AAACCC", (0, 6), "5", "TT") == c.ExtensionContext("without_record_context", "")


# ── declared V RSS position ──────────────────────────────────────────


def test_declared_v_rss_positions():
    # IGHV1-58*01 (interpretation evidence): gene_start 143, gene_end 438, v_rs_start 324.
    assert c.v_rss_position(143, 438, 324) == c.RssPosition(-115, "inside_coding_interval", None)
    assert c.v_rss_position(170, 463, 466) == c.RssPosition(2, "after_coding_end", None)
    assert c.v_rss_position(170, 463, 464) == c.RssPosition(0, "after_coding_end", None)
    assert c.v_rss_position(170, 463, 100) == c.RssPosition(-364, "before_coding_interval", None)
    for missing in (None, c.ABSENT):
        found = c.v_rss_position(170, 463, missing)
        assert found.gap is None and found.position is None
    assert "Boolean" in c.v_rss_position(170, 463, True).problem


# ── historically affirmed spans ──────────────────────────────────────


def _note(label="IGLV3-21*04"):
    (note,) = [n for n in profiles.PROFILE.historical_notes if n.label == label]
    return note


def test_an_affirmed_span_needs_a_reviewed_interpretation():
    record = raw_record(IGL, "airr", "IGLV3-21*04")
    gapped = raw_gapped(IGL, "gapped", "IGLV3-21*04")
    found = c.affirmed_span(_note(), "airr", record["coding_sequence"], gapped)
    assert found.interval is None and "no reviewed interpretation" in found.reason


@pytest.mark.parametrize("export, gapped_export", [("airr", "gapped"), ("airr_ex", "gapped_ex")])
def test_a_reviewed_note_maps_through_an_exact_gapped_alignment(export, gapped_export):
    reviewed = dataclasses.replace(_note(), interpretation=profiles.NoteInterpretation(1, 339))
    record = raw_record(IGL, export, "IGLV3-21*04")
    gapped = raw_gapped(IGL, gapped_export, "IGLV3-21*04")
    expected = sum(1 for char in gapped[:339] if char != ".")        # counted independently
    found = c.affirmed_span(reviewed, export, record["coding_sequence"], gapped)
    assert (found.interval, found.reason) == ((0, expected), None)


def test_a_reviewed_note_stays_unresolved_without_a_matching_alignment_or_sequence():
    reviewed = dataclasses.replace(_note(), interpretation=profiles.NoteInterpretation(1, 339))
    coding = raw_record(IGL, "airr", "IGLV3-21*04")["coding_sequence"]
    gapped = raw_gapped(IGL, "gapped", "IGLV3-21*04")
    assert "no gapped record" in c.affirmed_span(reviewed, "airr", coding, None).reason
    changed = gapped.replace("A", "C", 1)
    assert "does not ungap" in c.affirmed_span(reviewed, "airr", coding, changed).reason
    other = "C" + coding[1:] if coding[0] != "C" else "A" + coding[1:]
    assert "not the sequence" in c.affirmed_span(reviewed, "airr", other, other).reason
    assert "not the sequence" in c.affirmed_span(reviewed, "gapped", coding, gapped).reason
    beyond = dataclasses.replace(reviewed, interpretation=profiles.NoteInterpretation(1, 400))
    assert "beyond" in c.affirmed_span(beyond, "airr", coding, gapped).reason
    for first, last in ((0, 339), (340, 339)):
        empty = dataclasses.replace(reviewed, interpretation=profiles.NoteInterpretation(first, last))
        assert "not an interval" in c.affirmed_span(empty, "airr", coding, gapped).reason


def test_a_reviewed_span_need_not_start_at_column_1():
    gapped = raw_gapped(IGL, "gapped", "IGLV3-21*04")
    coding = gapped.replace(".", "")
    for first, last in ((2, 339), (310, 312), (100, 200)):
        reviewed = dataclasses.replace(_note(), interpretation=profiles.NoteInterpretation(first, last))
        expected = (sum(1 for char in gapped[:first - 1] if char != "."),
                    sum(1 for char in gapped[:last] if char != "."))
        assert c.affirmed_span(reviewed, "airr", coding, gapped).interval == expected


# ── supplied span against a stated column ────────────────────────────


@pytest.mark.parametrize("source, fmt, gapped_fmt, label, column, relation, beyond, short", [
    (IGL, "airr", "gapped", "IGLV2-14*04", 337, "exceeds", 2, 0),
    (IGL, "airr_ex", "gapped_ex", "IGLV2-14*04", 337, "exceeds", 2, 0),
    (IGL, "airr", "gapped", "IGLV2-14*03", 337, "equals", 0, 0),
    (IGL, "airr_ex", "gapped_ex", "IGLV2-14*03", 337, "exceeds", 2, 0),
    (IGH, "airr", "gapped", "IGHV2-70*20", 319, "exceeds", 3, 0),
    (IGH, "airr", "gapped", "IGHV4-4*09", 318, "equals", 0, 0),
    (IGH, "airr", "gapped", "IGHV4-4*09", 330, "short_of", 0, 12),
])
def test_the_supplied_span_is_compared_with_a_stated_column(source, fmt, gapped_fmt, label, column, relation,
                                                            beyond, short):
    gapped = raw_gapped(source, gapped_fmt, label)
    blocks = m.alignment_blocks(gapped)
    found = c.span_against_column(blocks, column)
    assert (found.last_column, found.relation, found.bases_beyond, found.columns_short) == (
        len(gapped.rstrip(".")), relation, beyond, short)
    assert found.bases_beyond == sum(1 for char in gapped[column:] if char != ".")


# ── codons ───────────────────────────────────────────────────────────


@pytest.mark.parametrize("codon, residue", [("TGT", "C"), ("TGC", "C"), ("TGG", "W"), ("TTT", "F"), ("TTC", "F"),
                                            ("GCT", "A"), ("TAT", "Y"), ("TGA", "*"), ("GGT", "G")])
def test_landmark_codons_translate_with_the_standard_code(codon, residue):
    assert c.residue(codon) == residue


# ── interpreted inventories ──────────────────────────────────────────

J_CODONS = {"9606.IGH_VDJ": {"TGG"}, "9606.IGKappa_VJ": {"TTT", "TTC"}, "9606.IGLambda_VJ": {"TTT", "TTC"}}
















def labels_in(inventory, finding) -> set:
    return {inventory.record(s).label for s in finding.subjects if s.id in {r.id for r in inventory.records}}












# ── J in the pinned exports ──────────────────────────────────────────


@pytest.mark.parametrize("source", SETS, ids=_key)
def test_every_j_occurrence_gets_a_source_stated_anchor_and_phase_value(source):
    lock, inventory = interpreted(source)
    key = f"{source[0]}/{source[1]}"
    for fmt in GAPPED:
        anchors = claims(inventory, "anchor", fmt, "J")
        assert len(anchors) == EXPECT["sets"][key][fmt]["J"]
        for anchor in anchors:
            record = inventory.subject_record(anchor)
            raw_j = raw_record(source, fmt, record.label)
            assert_status(anchor, "accepted", "source", "stated", raw_j["j_cdr3_end"] - 1)
            assert anchor.coordinate_system == "coding_sequence"
            assert anchor.evidence == (at(record, "j_cdr3_end"),) and anchor.mappings == ()
            candidate = inventory.candidate(anchor.subject)
            assert candidate.sequence[anchor.value.get():anchor.value.get() + 3] in J_CODONS[source[0]]
            phase = claim_of(inventory, candidate, "coding_phase")
            assert_status(phase, "accepted", "source", "stated", raw_j["j_codon_frame"] - 1)
            assert phase.evidence == (at(record, "j_codon_frame"),)
            independence = claim_of(inventory, candidate, "phase_independence")
            assert_status(independence, "unresolved", "none", "stated")
            assert independence.value.get() is None
            assert set(independence.evidence) == {at(record, "j_cdr3_end"), at(record, "j_codon_frame")}


@pytest.mark.parametrize("example", EXPECT["j_examples"], ids=lambda e: e["label"])
def test_j_examples_in_the_interpreted_inventory(example):
    _, inventory = interpreted(SET_KEYS[example["set"]])
    candidate = candidate_of(inventory, example["export"], example["label"])
    assert claim_of(inventory, candidate, "anchor").value.get() == example["anchor_coding0"]
    assert claim_of(inventory, candidate, "coding_phase").value.get() == example["phase"]


def test_no_j_phase_is_ever_claimed_independent():
    for source in SETS:
        _, inventory = interpreted(source)
        for claim in claims(inventory, "phase_independence", segment="J"):
            assert claim.status == "unresolved" and claim.value.get() is None


# ── V in the pinned exports ──────────────────────────────────────────


@pytest.mark.parametrize("source", SETS, ids=_key)
def test_every_v_occurrence_is_bound_to_its_own_gapped_record_and_accepted(source):
    lock, inventory = interpreted(source)
    key = f"{source[0]}/{source[1]}"
    mappings = {x.id: x for x in inventory.mappings}
    for fmt in GAPPED:
        anchors = claims(inventory, "anchor", fmt, "V")
        assert len(anchors) == EXPECT["sets"][key][fmt]["V"]
        for anchor in anchors:
            record = inventory.subject_record(anchor)
            cdr3_start = imgt_cdr3_start(raw_record(source, fmt, record.label))
            assert_status(anchor, "accepted", "source", "stated", cdr3_start - 4)
            (mapping_id,) = anchor.mappings
            mapping = mappings[mapping_id]
            assert (mapping.kind, mapping.validation, mapping.input, mapping.output.asset) == (
                "gapped_alignment", "exact", anchor.subject, GAPPED[fmt])
            gapped = inventory.observation(mapping.output).value.get()
            assert gapped["header"] == record.label
            assert gapped["sequence"] == raw_gapped(source, GAPPED[fmt], record.label)
            (delineation,) = [j for j, d in enumerate(raw_record(source, fmt, record.label)["v_gene_delineations"])
                              if d["delineation_scheme"] == "IMGT"]
            assert set(anchor.evidence) == {at(record, "v_gene_delineations", delineation, "cdr3_start"),
                                            at(record, "v_gene_delineations", delineation, "delineation_scheme"),
                                            mapping.output}


@pytest.mark.parametrize("example", EXPECT["v_examples"], ids=lambda e: f"{e['label']}-{e['export']}")
def test_v_examples_in_the_interpreted_inventory(example):
    _, inventory = interpreted(SET_KEYS[example["set"]])
    candidate = candidate_of(inventory, example["export"], example["label"])
    anchor = claim_of(inventory, candidate, "anchor")
    assert anchor.value.get() == example["anchor_coding0"]
    assert candidate.sequence[example["anchor_coding0"]:example["anchor_coding0"] + 3] == example["anchor_codon"]


@pytest.mark.parametrize("source", SETS, ids=_key)
def test_a_v_phase_is_only_ever_anchor_derived(source):
    _, inventory = interpreted(source)
    for anchor in claims(inventory, "anchor", segment="V"):
        candidate = inventory.candidate(anchor.subject)
        phase = claim_of(inventory, candidate, "coding_phase")
        # Every V anchor in these exports is a multiple of 3 (interpretation evidence), so the
        # anchor-derived phase is 0; it is GenAIRR's derivation, not a source-stated phase.
        assert_status(phase, "accepted", "assertion", "anchor_derived", 0)
        assert anchor.value.get() % 3 == 0
        independence = claim_of(inventory, candidate, "phase_independence")
        assert_status(independence, "accepted", "assertion", "anchor_derived", "anchor_derived")
    assert not [x for x in inventory.claims if x.predicate == "phase_independence"
                and x.value.get() == "independent"]


@pytest.mark.parametrize("source", SETS, ids=_key)
def test_noncanonical_landmark_codons_are_findings_and_never_relocated(source):
    _, inventory = interpreted(source)
    expected = EVIDENCE["conventions"]["v_landmark"]["residues_other_than_cys_at_imgt_104"][
        f"{source[0].split('.', 1)[1]}.{source[1]}"]
    found = {}
    for finding in findings(inventory, "noncanonical_landmark_codon"):
        observations = finding.observations.get()
        assert finding.capabilities == ("canonical_productive",)
        assert observations["segment"] == "V" and observations["expected_residues"] == ["C"]
        found[observations["label"]] = observations["codon"]
        for subject in finding.subjects:
            if subject.id in {c.id for c in inventory.candidates}:
                anchor = claim_of(inventory, inventory.candidate(subject), "anchor")
                assert anchor.status == "accepted" and anchor.value.get() == observations["anchor"]
                assert inventory.candidate(subject).sequence[observations["anchor"]:observations["anchor"] + 3] == \
                    observations["codon"]
        carried_by = {fmt for fmt in GAPPED if any(r["label"] == observations["label"]
                                                   for r in raw_records(source, fmt))}
        assert {s.asset for s in finding.subjects} == carried_by
    assert found == expected


def test_iglv3_2_03_keeps_its_landmark_at_codon_104():
    _, inventory = interpreted(IGL)
    for fmt in GAPPED:
        candidate = candidate_of(inventory, fmt, "IGLV3-2*03")
        anchor = claim_of(inventory, candidate, "anchor")
        assert anchor.status == "accepted" and anchor.value.get() == 258
        assert candidate.sequence[258:261] == "GCT"
    (finding,) = [f for f in findings(inventory, "noncanonical_landmark_codon")
                  if f.observations.get()["label"] == "IGLV3-2*03"]
    assert finding.observations.get() == {"label": "IGLV3-2*03", "segment": "V", "anchor": 258, "codon": "GCT",
                                          "residue": "A", "expected_residues": ["C"],
                                          "coordinate_system": "coding_sequence"}


def test_identical_dna_identities_keep_their_own_gapped_records():
    _, inventory = interpreted(IGH)
    outputs = {}
    for label in ("IGHV1-69*01", "IGHV1-69D*01"):
        anchor = claim_of(inventory, candidate_of(inventory, "airr", label), "anchor")
        (mapping_id,) = anchor.mappings
        output = inventory.mapping(mapping_id).output
        assert inventory.observation(output).value.get()["header"] == label
        outputs[label] = output
        for root in anchor.evidence:
            assert root.is_under(record_of(inventory, "airr", label).locator) or root == output
    assert outputs["IGHV1-69*01"] != outputs["IGHV1-69D*01"]
    assert candidate_of(inventory, "airr", "IGHV1-69*01").sha256 == candidate_of(inventory, "airr",
                                                                                 "IGHV1-69D*01").sha256


def test_gapped_alignments_never_cross_export_variants():
    for source in SETS:
        _, inventory = interpreted(source)
        for mapping in inventory.mappings:
            if mapping.kind == "gapped_alignment":
                assert mapping.output.asset == GAPPED[mapping.input.asset]


# ── coding intervals in the pinned exports ───────────────────────────


@pytest.mark.parametrize("source", SETS, ids=_key)
def test_the_coding_interval_is_accepted_only_where_it_reproduces_the_coding_sequence(source):
    _, inventory = interpreted(source)
    extended = EXPECT["sets"][f"{source[0]}/{source[1]}"]["airr_ex_coding_suffix_beyond_airr"]
    mappings = {x.id: x for x in inventory.mappings}
    for fmt in GAPPED:
        for claim in claims(inventory, "coding_interval", fmt):
            record = inventory.subject_record(claim)
            raw_rec = raw_record(source, fmt, record.label)
            (mapping_id,) = claim.mappings
            mapping = mappings[mapping_id]
            assert (mapping.kind, mapping.input, mapping.output) == ("coding_interval", at(record, "sequence"),
                                                                    claim.subject)
            assert set(claim.evidence) == {at(record, "gene_start"), at(record, "gene_end"), at(record, "sequence")}
            if fmt == "airr_ex" and record.label in extended:
                assert_status(claim, "unresolved", "none", "stated")
                assert mapping.validation == "mismatch"
            else:
                assert_status(claim, "accepted", "source", "stated", [raw_rec["gene_start"] - 1, raw_rec["gene_end"]])
                assert claim.coordinate_system == "record_sequence"
                assert mapping.validation == "exact"
                assert mapping.blocks == (m.Block((raw_rec["gene_start"] - 1, raw_rec["gene_end"]),
                                                  (0, len(raw_rec["coding_sequence"]))),)
    not_reproduced = findings(inventory, "coding_interval_not_reproduced")
    assert {f.observations.get()["label"]: f.observations.get()["extra"] for f in not_reproduced} == extended
    for finding in not_reproduced:
        assert {s.asset for s in finding.subjects} == {"airr_ex"}
        assert finding.observations.get()["relation"] == "candidate_extends_interval_3prime"


# ── scope ────────────────────────────────────────────────────────────


def test_a_lock_outside_the_profile_scope_gets_no_accepted_coordinate_claim():
    blobs = pinned(IGL)
    document = json.loads(blobs["airr"].decode("utf-8"))
    blobs["airr"] = json.dumps(document, indent=1).encode("utf-8")       # same content, drifted bytes
    lock, inventory = interpret(blobs)
    assert m.profile_scope_problem(lock) is not None
    m.check_inventory(inventory, lock)
    coordinate = [x for x in inventory.claims if x.predicate in m.COORDINATE_PREDICATES]
    assert coordinate and {x.status for x in coordinate} == {"candidate", "unresolved"}
    j7 = candidate_of(inventory, "airr_ex", "IGLJ7*02")
    assert_status(claim_of(inventory, j7, "anchor"), "candidate", "none", "stated", 7)
    v32 = candidate_of(inventory, "airr", "IGLV3-2*03")
    assert_status(claim_of(inventory, v32, "anchor"), "candidate", "none", "stated", 258)
    (finding,) = findings(inventory, "coordinate_profile_out_of_scope")
    assert finding.observations.get()["problem"] == m.profile_scope_problem(lock)
    assert finding.capabilities == ("source_annotated_junction", "canonical_productive")


def test_without_gapped_exports_no_v_anchor_is_accepted():
    lock, inventory = interpreted(IGL, ("airr", "airr_ex"))
    assert m.profile_scope_problem(lock) is None
    for anchor in claims(inventory, "anchor", segment="V"):
        assert_status(anchor, "unresolved", "none", "stated")
        assert anchor.mappings == ()
    for anchor in claims(inventory, "anchor", segment="J"):
        assert anchor.status == "accepted"
    unavailable = findings(inventory, "v_landmark_alignment_unavailable")
    assert len({f.observations.get()["label"] for f in unavailable}) == 88
    assert {f.observations.get()["reason"] for f in unavailable} == {"the gapped export is not locked",
                                                                     "the gapped_ex export is not locked"}


# ── V failures through the whole interpretation ──────────────────────




def _set_cdr3(value):
    def edit(fmt, record):
        if record["label"] == V_LABEL:
            record["v_gene_delineations"][0]["cdr3_start"] = value
    return edit


@pytest.mark.parametrize("value, problem", [
    (0, "missing"), (None, "null"), (True, "Boolean"), (268.0, "not an integer"), (3, "below 4"),
])
def test_an_unusable_cdr3_start_leaves_the_v_anchor_unresolved_even_in_scope(monkeypatch, value, problem):
    lock, inventory = interpret_in_scope(monkeypatch, small(IGL, {V_LABEL}, edit=_set_cdr3(value)))
    for fmt in GAPPED:
        candidate = candidate_of(inventory, fmt, V_LABEL)
        assert_status(claim_of(inventory, candidate, "anchor"), "unresolved", "none", "stated")
        assert_status(claim_of(inventory, candidate, "coding_phase"), "unresolved", "none", "anchor_derived")
        assert claim_of(inventory, candidate, "phase_independence").status == "unresolved"
    (finding,) = findings(inventory, "v_landmark_field_unusable")
    assert problem in finding.observations.get()["reason"]
    assert labels_in(inventory, finding) == {V_LABEL}


def test_a_cdr3_start_the_gapped_record_does_not_bind_is_unresolved_even_in_scope(monkeypatch):
    raw_v = raw_record(IGL, "airr", V_LABEL)
    lock, inventory = interpret_in_scope(monkeypatch,
                                         small(IGL, {V_LABEL}, edit=_set_cdr3(imgt_cdr3_start(raw_v) + 3)))
    candidate = candidate_of(inventory, "airr", V_LABEL)
    anchor = claim_of(inventory, candidate, "anchor")
    assert_status(anchor, "unresolved", "none", "stated")
    (mapping_id,) = anchor.mappings
    assert inventory.mapping(mapping_id).validation == "exact"         # the alignment itself is fine
    (finding,) = findings(inventory, "v_landmark_binding_failed")
    assert "differs" in finding.observations.get()["reason"]


def test_a_gapped_record_that_differs_from_the_coding_sequence_is_a_mismatch_mapping(monkeypatch):
    def change(fmt, header, sequence):
        if fmt == "gapped" and header == V_LABEL:
            position = sequence.index("A")
            sequence = sequence[:position] + "C" + sequence[position + 1:]
        return header, sequence
    lock, inventory = interpret_in_scope(monkeypatch, small(IGL, {V_LABEL}, gapped_edit=change))
    anchor = claim_of(inventory, candidate_of(inventory, "airr", V_LABEL), "anchor")
    assert_status(anchor, "unresolved", "none", "stated")
    (mapping_id,) = anchor.mappings
    mapping = inventory.mapping(mapping_id)
    assert (mapping.validation, mapping.blocks) == ("mismatch", ())
    assert claim_of(inventory, candidate_of(inventory, "airr_ex", V_LABEL), "anchor").status == "accepted"


@pytest.mark.parametrize("gapped_edit, reason", [
    (lambda fmt, h, s: None if fmt == "gapped" else (h, s), "no gapped record is labelled exactly"),
    (lambda fmt, h, s: (h + "_x", s) if fmt == "gapped" else (h, s), "no gapped record is labelled exactly"),
    (lambda fmt, h, s: (h.lower(), s) if fmt == "gapped" else (h, s), "no gapped record is labelled exactly"),
])
def test_only_an_exactly_labelled_gapped_record_of_the_same_variant_binds(monkeypatch, gapped_edit, reason):
    def edit_v_only(fmt, header, sequence):
        return gapped_edit(fmt, header, sequence) if header == V_LABEL else (header, sequence)
    blobs = small(IGL, {V_LABEL, "IGLJ1*01"}, gapped_edit=edit_v_only)
    lock, inventory = interpret_in_scope(monkeypatch, blobs)
    anchor = claim_of(inventory, candidate_of(inventory, "airr", V_LABEL), "anchor")
    assert_status(anchor, "unresolved", "none", "stated")
    assert anchor.mappings == ()
    (finding,) = findings(inventory, "v_landmark_alignment_unavailable")
    assert reason in finding.observations.get()["reason"]
    assert claim_of(inventory, candidate_of(inventory, "airr_ex", V_LABEL), "anchor").status == "accepted"


def test_an_ambiguous_gapped_label_binds_nothing(monkeypatch):
    blobs = small(IGL, {V_LABEL})
    blobs["gapped"] = blobs["gapped"] + blobs["gapped"]
    lock, inventory = interpret_in_scope(monkeypatch, blobs)
    anchor = claim_of(inventory, candidate_of(inventory, "airr", V_LABEL), "anchor")
    assert anchor.status == "unresolved" and anchor.mappings == ()
    (finding,) = findings(inventory, "v_landmark_alignment_unavailable")
    assert "2 gapped records are labelled" in finding.observations.get()["reason"]


def test_a_label_carried_by_two_records_of_one_export_binds_no_gapped_record(monkeypatch):
    blobs = small(IGL, {V_LABEL})
    document = json.loads(blobs["airr"].decode("utf-8"))
    records = document["GermlineSet"][0]["allele_descriptions"]
    records.append(dict(records[0], allele_description_id="X1"))
    blobs["airr"] = encode(document)
    lock, inventory = interpret_in_scope(monkeypatch, blobs)
    for record in [r for r in inventory.records if r.locator.asset == "airr"]:
        (candidate,) = inventory.candidates_of(record.locator)
        assert claim_of(inventory, candidate, "anchor").status == "unresolved"
    (finding,) = findings(inventory, "v_landmark_alignment_unavailable")
    assert "2 records of airr carry the label" in finding.observations.get()["reason"]


def test_a_landmark_codon_with_a_gap_in_the_gapped_record_is_unresolved(monkeypatch):
    def shift(fmt, header, sequence):
        return header, (sequence[:310] + "." + sequence[310:] if fmt == "gapped" else sequence)
    lock, inventory = interpret_in_scope(monkeypatch, small(IGL, {V_LABEL}, gapped_edit=shift))
    assert claim_of(inventory, candidate_of(inventory, "airr", V_LABEL), "anchor").status == "unresolved"
    (finding,) = findings(inventory, "v_landmark_binding_failed")
    assert "310-312" in finding.observations.get()["reason"]


@pytest.mark.parametrize("value, problem", [(None, "no IMGT delineation"), ([], "no IMGT delineation"),
                                            ("x", "no IMGT delineation")])
def test_a_v_record_without_a_single_imgt_delineation_is_unresolved(monkeypatch, value, problem):
    def edit(fmt, record):
        if record["label"] == V_LABEL:
            record["v_gene_delineations"] = value
    lock, inventory = interpret_in_scope(monkeypatch, small(IGL, {V_LABEL}, edit=edit))
    anchor = claim_of(inventory, candidate_of(inventory, "airr", V_LABEL), "anchor")
    assert anchor.status == "unresolved"
    assert anchor.evidence[0] == at(record_of(inventory, "airr", V_LABEL), "v_gene_delineations")
    (finding,) = findings(inventory, "v_landmark_field_unusable")
    assert problem in finding.observations.get()["reason"]


def test_two_imgt_delineations_are_not_chosen_between(monkeypatch):
    def edit(fmt, record):
        if record["label"] == V_LABEL:
            record["v_gene_delineations"] = record["v_gene_delineations"] * 2
    lock, inventory = interpret_in_scope(monkeypatch, small(IGL, {V_LABEL}, edit=edit))
    assert claim_of(inventory, candidate_of(inventory, "airr", V_LABEL), "anchor").status == "unresolved"
    (finding,) = findings(inventory, "v_landmark_field_unusable")
    assert "2 IMGT delineations" in finding.observations.get()["reason"]


# ── J failures through the whole interpretation ──────────────────────


J_LABEL = "IGLJ7*02"


def _set_j(field_name, value):
    def edit(fmt, record):
        if record["label"] == J_LABEL:
            record[field_name] = value
    return edit


def test_a_stated_phase_that_disagrees_with_the_landmark_is_unsupported_even_in_scope(monkeypatch):
    lock, inventory = interpret_in_scope(monkeypatch, small(IGL, {J_LABEL}, edit=_set_j("j_codon_frame", 3)))
    for fmt in GAPPED:
        candidate = candidate_of(inventory, fmt, J_LABEL)
        assert_status(claim_of(inventory, candidate, "anchor"), "candidate", "none", "stated", 7)
        assert_status(claim_of(inventory, candidate, "coding_phase"), "candidate", "none", "stated", 2)
    (finding,) = findings(inventory, "phase_landmark_projection_unsupported")
    assert finding.observations.get() == {"label": J_LABEL, "segment": "J", "anchor": 7, "phase": 2,
                                          "anchor_modulo_3": 1}
    assert finding.capabilities == ("source_annotated_junction", "canonical_productive")


@pytest.mark.parametrize("field_name, value, problem", [
    ("j_cdr3_end", True, "Boolean"), ("j_cdr3_end", None, "null"), ("j_cdr3_end", 8.0, "not an integer"),
    ("j_cdr3_end", 40, "beyond"), ("j_codon_frame", False, "Boolean"), ("j_codon_frame", 4, "1, 2 or 3"),
])
def test_unusable_j_fields_are_unresolved_and_never_reinterpreted(monkeypatch, field_name, value, problem):
    lock, inventory = interpret_in_scope(monkeypatch, small(IGL, {J_LABEL}, edit=_set_j(field_name, value)))
    candidate = candidate_of(inventory, "airr_ex", J_LABEL)
    unresolved, other = ("anchor", "coding_phase") if field_name == "j_cdr3_end" else ("coding_phase", "anchor")
    assert_status(claim_of(inventory, candidate, unresolved), "unresolved", "none", "stated")
    other_claim = claim_of(inventory, candidate, other)
    if other == "anchor":
        assert_status(other_claim, "accepted", "source", "stated", 7)     # the landmark stands on its own
    else:
        assert_status(other_claim, "candidate", "none", "stated", 1)     # a phase needs a consistent landmark
    (finding,) = findings(inventory, "j_landmark_field_unusable")
    observations = finding.observations.get()
    assert observations["field"] == field_name and problem in observations["reason"]


def test_a_noncanonical_j_codon_is_reported_where_the_source_puts_it(monkeypatch):
    lock, inventory = interpret_in_scope(monkeypatch, small(IGL, {J_LABEL}, edit=_set_j("j_cdr3_end", 11)))
    candidate = candidate_of(inventory, "airr_ex", J_LABEL)
    anchor = claim_of(inventory, candidate, "anchor")
    assert_status(anchor, "accepted", "source", "stated", 10)
    (finding,) = findings(inventory, "noncanonical_landmark_codon")
    observations = finding.observations.get()
    assert (observations["anchor"], observations["codon"], observations["expected_residues"]) == (
        10, candidate.sequence[10:13], ["F"])


# ── boundaries of the coordinate rules ───────────────────────────────


def test_a_j_landmark_codon_may_end_exactly_at_the_coding_end():
    coding = "A" * 35 + "TTC"
    found = c.j_landmark(coding, 36, 3)
    assert (found.anchor, found.codon, found.anchor_problem) == (35, "TTC", None)
    found = c.j_landmark(coding, 37, 3)
    assert found.anchor is None and "beyond" in found.anchor_problem


def test_a_v_landmark_codon_may_end_exactly_at_the_coding_end():
    coding = "A" * 309 + "TGT"             # 312 bases; the alignment has no gaps
    found = c.v_landmark(coding, 313, coding)
    assert (found.anchor, found.codon, found.problem) == (309, "TGT", None)
    found = c.v_landmark(coding, 314, coding)
    assert found.anchor is None and found.kind == "field" and "beyond" in found.problem


def test_a_v_rss_starting_at_gene_start_is_inside_the_coding_interval():
    assert c.v_rss_position(5, 10, 5) == c.RssPosition(-6, "inside_coding_interval", None)
    assert c.v_rss_position(5, 10, 4) == c.RssPosition(-7, "before_coding_interval", None)


def test_a_v_phase_is_the_anchor_modulo_3_and_not_a_constant(monkeypatch):
    raw_v = raw_record(IGL, "airr", V_LABEL)
    cdr3_start = imgt_cdr3_start(raw_v)

    def drop_first_base(fmt, record):        # the record starts one base later; the frame shifts
        if record["label"] == V_LABEL:
            record["coding_sequence"] = record["coding_sequence"][1:]
            record["gene_start"] += 1
            record["v_gene_delineations"][0]["cdr3_start"] -= 1

    def gap_first_column(fmt, header, sequence):
        return header, ("." + sequence[1:] if header == V_LABEL else sequence)

    blobs = small(IGL, {V_LABEL}, edit=drop_first_base, gapped_edit=gap_first_column)
    lock, inventory = interpret_in_scope(monkeypatch, blobs)
    for fmt in GAPPED:
        candidate = candidate_of(inventory, fmt, V_LABEL)
        anchor = claim_of(inventory, candidate, "anchor")
        assert_status(anchor, "accepted", "source", "stated", cdr3_start - 5)
        assert (cdr3_start - 5) % 3 == 2
        assert_status(claim_of(inventory, candidate, "coding_phase"), "accepted", "assertion", "anchor_derived", 2)
        assert claim_of(inventory, candidate, "coding_interval").status == "accepted"
