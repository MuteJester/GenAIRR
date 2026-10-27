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
from GenAIRR._reference_sources import model as m
from GenAIRR._reference_sources import profiles
from tests.test_ogrdb_source_inventory import (
    EXPECT, IGH, IGK, IGL, ROOT, error, raw, raw_records,
)

EVIDENCE = json.loads((ROOT / "interpretation_evidence.json").read_bytes().decode("utf-8"))
PROFILE_REF = m.CoordinateProfileRef("ogrdb_api_v2_human_ig_observed_v1", 1)
SET_KEYS = {"9606.IGH_VDJ/10": IGH, "9606.IGKappa_VJ/5": IGK, "9606.IGLambda_VJ/4": IGL}
NOTE_LABELS = tuple(row["label"] for row in EXPECT["igl4_airr_ex_coding_end_observations"]
                    if row["finding"] == "coding_end_before_declared_v_rss")


# ── helpers ──────────────────────────────────────────────────────────


def raw_record(source, fmt, label) -> dict:
    """One record of a raw export, read with the standard library JSON reader."""
    (record,) = [r for r in raw_records(source, fmt) if r["label"] == label]
    return record


def raw_gapped(source, fmt, label) -> str:
    """One gapped FASTA record's sequence, read independently: header line, then the joined lines."""
    text = raw(source, fmt).decode("utf-8")
    found = []
    for chunk in text.split(">")[1:]:
        header, _, body = chunk.partition("\n")
        if header.rstrip("\r") == label:
            found.append("".join(line.rstrip("\r") for line in body.split("\n")))
    (sequence,) = found
    return sequence


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
    assert [u for u, _ in profiles.PROFILE.unresolved] == [u["id"] for u in EVIDENCE["unresolved"]]
    assert profiles.PROFILE.j_phase_independence == "unresolved"
    assert profiles.PROFILE.v_phase_derivation == "anchor_derived"
    assert profiles.PROFILE.functionality_classes == {}      # no provider value is classified


def test_profile_landmark_residues():
    assert profiles.PROFILE.landmark_residues == {"V": {"IGH": ("C",), "IGK": ("C",), "IGL": ("C",)},
                                                  "J": {"IGH": ("W",), "IGK": ("F",), "IGL": ("F",)}}


def test_historical_notes_are_bound_to_their_exact_text_and_sequences():
    notes = profiles.PROFILE.historical_notes
    assert sorted(n.label for n in notes) == sorted(NOTE_LABELS)
    for note in notes:
        assert (note.set_id, note.release, note.field, note.numbering) == ("9606.IGLambda_VJ", 4, "curation",
                                                                          "imgt_gapped")
        assert note.interpretation is None           # no reviewed interpretation exists
        texts = {raw_record(IGL, fmt, note.label)["curation"] for fmt in ("airr", "airr_ex")}
        (text,) = texts
        assert note.text_sha256 == sha(text)
        assert note.quotes and all(quote in text for quote in note.quotes)
        assert dict(note.sequences) == {fmt: sha(raw_record(IGL, fmt, note.label)["coding_sequence"])
                                        for fmt in ("airr", "airr_ex")}


# ── stated positions ─────────────────────────────────────────────────


@pytest.mark.parametrize("value, problem", [
    (None, "null"), (c.ABSENT, "absent"), (True, "Boolean"), (False, "Boolean"), (8.0, "not an integer"),
    ("8", "not an integer"), ([8], "not an integer"), (-1, "below 1"),
])
def test_stated_positions_are_exact_integers(value, problem):
    found = c.stated_position(value)
    assert found.value is None and problem in found.problem


def test_zero_is_a_missing_delineation_coordinate_only_where_the_profile_says_so():
    assert c.stated_position(0, zero_is_missing=True).value is None
    assert "missing" in c.stated_position(0, zero_is_missing=True).problem
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
    reviewed = dataclasses.replace(_note(), interpretation=profiles.NoteInterpretation(last_affirmed_column=339))
    record = raw_record(IGL, export, "IGLV3-21*04")
    gapped = raw_gapped(IGL, gapped_export, "IGLV3-21*04")
    expected = sum(1 for char in gapped[:339] if char != ".")        # counted independently
    found = c.affirmed_span(reviewed, export, record["coding_sequence"], gapped)
    assert (found.interval, found.reason) == ((0, expected), None)


def test_a_reviewed_note_stays_unresolved_without_a_matching_alignment_or_sequence():
    reviewed = dataclasses.replace(_note(), interpretation=profiles.NoteInterpretation(last_affirmed_column=339))
    coding = raw_record(IGL, "airr", "IGLV3-21*04")["coding_sequence"]
    gapped = raw_gapped(IGL, "gapped", "IGLV3-21*04")
    assert "no gapped record" in c.affirmed_span(reviewed, "airr", coding, None).reason
    changed = gapped.replace("A", "C", 1)
    assert "does not ungap" in c.affirmed_span(reviewed, "airr", coding, changed).reason
    other = "C" + coding[1:] if coding[0] != "C" else "A" + coding[1:]
    assert "not the sequence" in c.affirmed_span(reviewed, "airr", other, other).reason
    assert "not the sequence" in c.affirmed_span(reviewed, "gapped", coding, gapped).reason
    beyond = dataclasses.replace(reviewed, interpretation=profiles.NoteInterpretation(last_affirmed_column=400))
    assert "beyond" in c.affirmed_span(beyond, "airr", coding, gapped).reason


# ── codons ───────────────────────────────────────────────────────────


@pytest.mark.parametrize("codon, residue", [("TGT", "C"), ("TGC", "C"), ("TGG", "W"), ("TTT", "F"), ("TTC", "F"),
                                            ("GCT", "A"), ("TAT", "Y"), ("TGA", "*"), ("GGT", "G")])
def test_landmark_codons_translate_with_the_standard_code(codon, residue):
    assert c.residue(codon) == residue
