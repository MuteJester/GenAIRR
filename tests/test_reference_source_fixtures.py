"""Frozen OGRDB source fixtures: exact bytes, provenance and hand-checked values.

``tests/data/reference_sources`` holds complete OGRDB API v2 responses (both AIRR JSON
variants and both IMGT-gapped FASTA variants for human IGH release 10, IGK release 5 and
IGL release 4), a provenance manifest, the interpretation evidence, hand-checked
expectations and a candidate IGH selection.  These tests keep those files from drifting:
every fixture must match its recorded digest and size, and the hand-checked values must
still read off the raw bytes.  Values are recomputed here by plain slicing, independently
of any GenAIRR import code.
"""
import hashlib
import json
import re
from collections import Counter
from pathlib import Path

import pytest

from GenAIRR._reference.capabilities import CAPABILITIES

ROOT = Path(__file__).resolve().parent / "data" / "reference_sources"
SOURCES = json.loads((ROOT / "sources.json").read_bytes())
EVIDENCE = json.loads((ROOT / "interpretation_evidence.json").read_bytes())
EXPECT = json.loads((ROOT / "expectations.json").read_bytes())
SELECTION = json.loads((ROOT / "human_igh_release10_candidate_selection.json").read_bytes())
GENERATED_TEXT = ("sources.json", "interpretation_evidence.json", "expectations.json",
                  "human_igh_release10_candidate_selection.json")
FIXTURE_FILES = {"airr": "airr.json", "airr_ex": "airr_ex.json",
                 "gapped": "gapped.fasta", "gapped_ex": "gapped_ex.fasta"}
CYS = {"TGT", "TGC"}


def _committed_assets():
    for source_set in SOURCES["sets"]:
        for asset in source_set["assets"]:
            if asset["fixture"] is not None:
                yield source_set, asset


def _set_dir(key):
    set_id, release = key.split("/")
    return ROOT / "ogrdb" / set_id / release


def _records(key, export):
    document = json.loads((_set_dir(key) / FIXTURE_FILES[export]).read_bytes())
    (germline_set,) = document["GermlineSet"]
    return germline_set["allele_descriptions"]


def _by_label(key, export):
    records = _records(key, export)
    by_label = {record["label"]: record for record in records}
    assert len(by_label) == len(records), "duplicate labels"
    return by_label


def _fasta(key, export):
    text = (_set_dir(key) / FIXTURE_FILES[export]).read_bytes().decode("ascii")
    records, name = {}, None
    for line in text.split("\n"):
        if line.startswith(">"):
            name = line[1:]
            assert name not in records, f"duplicate FASTA header {name}"
            records[name] = []
        elif line:
            records[name].append(line)
    return {label: "".join(lines) for label, lines in records.items()}


def _natural(label):
    return [int(part) if part.isdigit() else part for part in re.split(r"(\d+)", label)]


# --- exact bytes and provenance -------------------------------------------------------

COMMITTED = list(_committed_assets())


@pytest.mark.parametrize("source_set,asset", COMMITTED, ids=[asset["fixture"] for _, asset in COMMITTED])
def test_committed_fixture_matches_recorded_digest_and_size(source_set, asset):
    data = (ROOT / asset["fixture"]).read_bytes()
    assert len(data) == asset["bytes"]
    assert hashlib.sha256(data).hexdigest() == asset["sha256"]
    assert asset["fixture_content"] == "exact_full_response_bytes"
    expected_url = (f"{SOURCES['provider']['base_url']}/germline/set/{source_set['set_id']}/"
                    f"{source_set['release']}/{asset['format']}")
    assert asset["url"] == expected_url and asset["final_url"] == expected_url
    assert asset["http_status"] == 200


def test_directory_holds_only_listed_files():
    listed = {asset["fixture"] for _, asset in COMMITTED} | set(GENERATED_TEXT)
    present = {path.relative_to(ROOT).as_posix() for path in ROOT.rglob("*") if path.is_file()}
    assert present == listed


@pytest.mark.parametrize("name", GENERATED_TEXT)
def test_generated_text_is_utf8_with_lf_endings(name):
    data = (ROOT / name).read_bytes()
    data.decode("utf-8")
    assert b"\r" not in data
    assert data.endswith(b"\n") and not data.endswith(b"\n\n")


def test_manifest_covers_the_supported_releases_and_both_record_exports():
    releases = {(s["set_id"], s["release"]) for s in SOURCES["sets"]}
    assert releases == {("9606.IGH_VDJ", 10), ("9606.IGKappa_VJ", 5), ("9606.IGLambda_VJ", 4)}
    for source_set in SOURCES["sets"]:
        committed = {a["format"] for a in source_set["assets"] if a["fixture"] is not None}
        assert committed == {"airr", "airr_ex", "gapped", "gapped_ex"}
        by_format = {a["format"]: a for a in source_set["assets"]}
        # The unsuffixed release endpoint serves the '_ex' JSON for these human sets.
        assert by_format[None]["sha256"] == by_format["airr_ex"]["sha256"]
        document = json.loads((_set_dir(f"{source_set['set_id']}/{source_set['release']}")
                               / "airr.json").read_bytes())
        (germline_set,) = document["GermlineSet"]
        assert germline_set["germline_set_ref"] == source_set["germline_set_ref"]
        assert germline_set["release_version"] == source_set["release"]
        # The JSON names the API software licence; the data licence is recorded separately.
        assert document["Info"]["license"]["name"] == "European Union Public License version 1.2"


def test_data_licence_is_recorded_separately_from_software_licence():
    assert SOURCES["data_licence"]["licence"] == "CC0 1.0 Universal"
    assert "source code" in SOURCES["software_licence"]["statement"]
    assert SOURCES["data_licence"]["licence"] != SOURCES["software_licence"]["licence"]
    assert SOURCES["citation"]["doi"] == "10.1093/nar/gkaf1094"


# --- set-level observations -------------------------------------------------------------

@pytest.mark.parametrize("key", sorted(EXPECT["sets"]))
def test_set_level_counts_and_variant_differences(key):
    expected = EXPECT["sets"][key]
    labels = {}
    for export in ("airr", "airr_ex"):
        records = _records(key, export)
        by_label = {r["label"]: r for r in records}
        labels[export] = by_label
        counts = Counter(r["sequence_type"] for r in records)
        want = expected[export]
        assert len(records) == want["records"]
        assert (counts["V"], counts["D"], counts["J"]) == (want["V"], want["D"], want["J"])
        assert sum(r["allele_description_ref"].endswith("_ext") for r in records) == want["refs_ending_ext"]
        mentioned = {p for r in records for p in (r["paralogs"] or [])}
        assert sorted(mentioned - by_label.keys()) == want["paralog_only_identities"]
    assert sorted(labels["airr"].keys() - labels["airr_ex"].keys()) == expected["labels_only_in_airr"]
    assert sorted(labels["airr_ex"].keys() - labels["airr"].keys()) == expected["labels_only_in_airr_ex"]

    suffixes = {}
    for label in labels["airr"].keys() & labels["airr_ex"].keys():
        short = labels["airr"][label]["coding_sequence"]
        long = labels["airr_ex"][label]["coding_sequence"]
        if short != long:
            assert long.startswith(short)
            suffixes[label] = long[len(short):]
    assert suffixes == expected["airr_ex_coding_suffix_beyond_airr"]

    groups = {}
    for record in labels["airr"].values():
        groups.setdefault((record["sequence_type"], record["coding_sequence"]), []).append(record["label"])
    shared = sorted(sorted(group) for group in groups.values() if len(group) > 1)
    assert shared == expected["airr_identical_coding_groups"]

    rss = Counter()
    for record in labels["airr"].values():
        if record["sequence_type"] == "V":
            start = record["v_rs_start"]
            rss["no_v_rs_start" if start is None else str(start - 1 - record["gene_end"])] += 1
    assert dict(rss) == expected["airr_v_rss_gap_counts"]

    for export in ("airr", "airr_ex"):
        sequences = labels[export]
        for record in sequences.values():
            start, end = record["gene_start"], record["gene_end"]
            if export == "airr" or record["label"] not in expected["airr_ex_coding_suffix_beyond_airr"]:
                assert record["sequence"][start - 1:end] == record["coding_sequence"], record["label"]


@pytest.mark.parametrize("key", sorted(EXPECT["sets"]))
@pytest.mark.parametrize("export,gapped_export", [("airr", "gapped"), ("airr_ex", "gapped_ex")])
def test_gapped_records_match_json_records_and_cdr3_start(key, export, gapped_export):
    records = _by_label(key, export)
    gapped = _fasta(key, gapped_export)
    assert len(gapped) == EXPECT["sets"][key]["gapped_records"][gapped_export]
    assert gapped.keys() == records.keys()
    for label, record in records.items():
        aligned = gapped[label]
        assert set(aligned) <= set("ACGT.")
        assert aligned.replace(".", "") == record["coding_sequence"], label
        if record["sequence_type"] != "V":
            assert "." not in aligned
            continue
        (delineation,) = record["v_gene_delineations"]
        assert delineation["delineation_scheme"] == "IMGT"
        assert delineation["aligned_sequence"] is None and delineation["unaligned_sequence"] is None
        assert not aligned.startswith(".")
        assert "." not in aligned[309:312]
        assert delineation["cdr3_start"] == 1 + sum(ch != "." for ch in aligned[:312]), label
        assert sum(ch != "." for ch in aligned[:309]) == delineation["cdr3_start"] - 4


@pytest.mark.parametrize("key", sorted(EXPECT["sets"]))
def test_non_cys_imgt_104_residues_are_the_recorded_ones(key):
    short = key.replace("9606.", "").replace("/", ".")
    recorded = EVIDENCE["conventions"]["v_landmark"]["residues_other_than_cys_at_imgt_104"][short]
    for export in ("airr", "airr_ex"):
        found = {}
        for record in _records(key, export):
            if record["sequence_type"] == "V":
                start = record["v_gene_delineations"][0]["cdr3_start"] - 4
                codon = record["coding_sequence"][start:start + 3]
                if codon not in CYS:
                    found[record["label"]] = codon
        assert found == recorded


def test_igl_coding_end_observations():
    records = _by_label("9606.IGLambda_VJ/4", "airr_ex")
    for row in EXPECT["igl4_airr_ex_coding_end_observations"]:
        record = records[row["label"]]
        sequence, coding = record["sequence"], record["coding_sequence"]
        start, end = record["gene_start"], record["gene_end"]
        assert (end, record["v_rs_start"]) == (row["gene_end"], row["v_rs_start"])
        declared = sequence[start - 1:end]
        assert coding.startswith(declared)
        assert coding[len(declared):] == row["coding_extension"]
        assert sequence[end:end + len(row["coding_extension"])] == row["sequence_after_gene_end"]
        if row["finding"] == "coding_end_before_declared_v_rss":
            assert sequence[end:record["v_rs_start"] - 1] == row["coding_extension"]


# --- hand-checked landmark examples -------------------------------------------------------

@pytest.mark.parametrize("row", EXPECT["j_examples"], ids=lambda row: row["label"])
def test_j_examples(row):
    record = _by_label(row["set"], row["export"])[row["label"]]
    sequence, coding = record["sequence"], record["coding_sequence"]
    assert (len(sequence), len(coding)) == (row["sequence_length"], row["coding_length"])
    assert (record["gene_start"], record["gene_end"]) == (row["gene_start"], row["gene_end"])
    assert (record["j_cdr3_end"], record["j_codon_frame"]) == (row["j_cdr3_end"], row["j_codon_frame"])
    anchor = row["anchor_coding0"]
    assert coding[anchor:anchor + 3] == row["anchor_codon"]
    assert (anchor - row["phase"]) % 3 == 0
    assert sequence[anchor:anchor + 3] == row["same_offset_in_sequence"]
    assert row["airr_prose_anchor_coding0"] == record["j_cdr3_end"] - record["gene_start"] < 0


@pytest.mark.parametrize("row", EXPECT["v_examples"], ids=lambda row: row["label"])
def test_v_examples(row):
    record = _by_label(row["set"], row["export"])[row["label"]]
    aligned = _fasta(row["set"], row["gapped_export"])[row["label"]]
    coding = record["coding_sequence"]
    assert record["functional"] is row["functional"]
    assert len(coding) == row["coding_length"]
    assert record["v_gene_delineations"][0]["cdr3_start"] == row["cdr3_start"]
    anchor = row["anchor_coding0"]
    assert coding[anchor:anchor + 3] == row["anchor_codon"]
    assert aligned[309:312] == row["gapped_columns_310_312"]
    assert sum(ch != "." for ch in aligned[:309]) == row["non_gap_in_columns_1_309"]
    assert sum(ch != "." for ch in aligned[:312]) == row["non_gap_in_columns_1_312"]
    assert record["sequence"][anchor:anchor + 3] == row["same_offset_in_sequence"]


# --- candidate IGH selection --------------------------------------------------------------

CHOSEN_AIRR_RECORDS = {
    # label: 3' extension carried only by the 'airr_ex' record (the alternative not adopted)
    "IGHV1-58*03": "A", "IGHV1-69*20": "A", "IGHV3-7*04": "A", "IGHV3-9*04": "ATA", "IGHV3-21*05": "A",
    "IGHV3-30*20": "A", "IGHV3-30*i02": "A", "IGHV3-53*05": "A", "IGHV4-31*11": "A", "IGHV4-39*08": "A",
    "IGHV4-39*09": "A", "IGHV4-59*12": "A", "IGHV4-59*i03": "A", "IGHV4-61*11": "A",
}
EXTENDED_AND_TRUNCATED = (
    "IGHV1-69*18", "IGHV4-4*09", "IGHV4-4*10", "IGHV4-30-4*i02", "IGHV4-39*i03", "IGHV4-59*11",
    "IGHV4-61*10", "IGHV4-61*i04", "IGHV5-51*06",
)
PARALOG_ONLY_IN_AIRR_EX = (
    "IGHD4-11*01", "IGHD5-18*01", "IGHV1-69D*01", "IGHV2-70*04", "IGHV3-23D*01", "IGHV3-30-3*03",
    "IGHV3-30-5*01", "IGHV3-30-5*02", "IGHV4-30-4*09",
)
FUNCTIONAL_FALSE = ("IGHV3-19*01", "IGHV3-25*02", "IGHV3-25*03", "IGHV3-25*05", "IGHV3-62*03", "IGHV6-1*03",
                    "IGHV8-51-1*02", "IGHV8-51-1*03")


def _differing(record, other):
    return {k for k in record.keys() | other.keys() if record.get(k) != other.get(k)} - {"allele_description_ref"}


def test_igh_candidate_selection_accounts_for_every_record():
    key = "9606.IGH_VDJ/10"
    airr, airr_ex = _by_label(key, "airr"), _by_label(key, "airr_ex")
    selected = SELECTION["selected"]
    omitted = {row["label"]: row for row in SELECTION["omitted"]}
    chosen = [label for segment in "VDJ" for label in selected[segment]]
    assert len(set(chosen)) == len(chosen)
    assert not set(chosen) & omitted.keys()
    assert set(chosen) | omitted.keys() == airr.keys()
    assert SELECTION["counts"] == {"V": 183, "D": 31, "J": 7, "omitted": 28, "inventory": 249}
    assert SELECTION["distinct_genes_selected"] == {"V": 61, "D": 26, "J": 6}
    assert {segment: len(selected[segment]) for segment in "VDJ"} == {"V": 183, "D": 31, "J": 7}
    assert selected["J"] == ["IGHJ1*01", "IGHJ2*01", "IGHJ3*02", "IGHJ4*02", "IGHJ5*02", "IGHJ6*02", "IGHJ6*03"]
    choices = {row["label"]: row for row in SELECTION["record_choices"]}
    for segment in "VDJ":
        assert selected[segment] == sorted(selected[segment], key=_natural)
        for label in selected[segment]:
            record, other = airr[label], airr_ex[label]
            assert record["sequence_type"] == segment
            difference = _differing(record, other)
            if label in CHOSEN_AIRR_RECORDS:
                assert difference == {"coding_sequence"}
            else:
                assert difference == set() and label not in choices
            assert record["functional"] is True
            assert "likely_truncated" not in (record["curational_tags"] or [])
    sequences = [airr[label]["coding_sequence"] for label in chosen]
    assert len(set(sequences)) == len(sequences)
    for label in selected["V"]:
        start = airr[label]["v_gene_delineations"][0]["cdr3_start"] - 4
        assert airr[label]["coding_sequence"][start:start + 3] in CYS
    for label in selected["J"]:
        record = airr[label]
        anchor = record["j_cdr3_end"] - 1
        assert record["coding_sequence"][anchor:anchor + 3] == "TGG"


def test_igh_record_choices_take_the_airr_record_and_keep_the_extension_as_alternative():
    key = "9606.IGH_VDJ/10"
    airr, airr_ex = _by_label(key, "airr"), _by_label(key, "airr_ex")
    choices = {row["label"]: row for row in SELECTION["record_choices"]}
    assert set(choices) == set(CHOSEN_AIRR_RECORDS)
    assert "ends at gene_end" in SELECTION["record_choice_reason"]
    assert "trimming" in SELECTION["record_choice_reason"]
    for label, extension in CHOSEN_AIRR_RECORDS.items():
        row, record = choices[label], airr[label]
        sequence, coding, end = record["sequence"], record["coding_sequence"], record["gene_end"]
        assert row["chosen_export"] == "airr"
        assert row["alternative"] == {"export": "airr_ex", "coding_length": len(coding) + len(extension),
                                      "extension_3prime": extension, "kept_for": "dossier"}
        assert airr_ex[label]["coding_sequence"] == coding + extension
        assert sequence[record["gene_start"] - 1:end] == coding and row["chosen_coding_ends_at_gene_end"] is True
        assert row["chosen_coding_length"] == len(coding)
        assert row["extension_matches_record_context_after_gene_end"] is True
        assert row["record_context_after_gene_end"] == sequence[end:end + len(extension)] == extension
        if label == "IGHV3-7*04":
            assert (end, record["v_rs_start"]) == (455, 457)
            assert row["extension_reaches_rss_minus_1"] is True
            assert sequence[end:record["v_rs_start"] - 1] == "A"
            assert sequence[record["v_rs_start"] - 1:record["v_rs_start"] + 6] == "CACAGTG"
            assert "IGL two-base class" in row["observation"]
        else:
            assert record["v_rs_start"] is None and row["extension_reaches_rss_minus_1"] is None


def test_igh_omissions_list_every_failed_rule():
    omitted = {row["label"]: [f["rule"] for f in row["failed_rules"]] for row in SELECTION["omitted"]}
    expected = {label: ["explicit_record_in_both_variants"] for label in PARALOG_ONLY_IN_AIRR_EX}
    expected.update({label: ["explicit_record_choice", "not_tagged_likely_truncated"]
                     for label in EXTENDED_AND_TRUNCATED})
    expected.update({label: ["provider_functional_true"] for label in FUNCTIONAL_FALSE})
    expected["IGHV3-38*01"] = ["not_tagged_likely_truncated"]
    expected["IGHV1-38-4*01"] = ["v_imgt_104_is_cys"]
    assert omitted == expected
    rules = [rule["rule"] for rule in SELECTION["rules"]]
    assert rules == ["explicit_record_in_both_variants", "explicit_record_choice", "provider_functional_true",
                     "not_tagged_likely_truncated", "v_gapped_record_matches", "v_imgt_104_is_cys",
                     "j_landmark_consistent"]
    cys_rule = next(rule["requirement"] for rule in SELECTION["rules"] if rule["rule"] == "v_imgt_104_is_cys")
    assert "canonical_productive" in cys_rule and "no user view narrowing" in cys_rule
    assert "IGHV3-20*02" in cys_rule and "IGHV3-25*04" in cys_rule
    selected_v = set(SELECTION["selected"]["V"])
    assert {"IGHV3-20*02", "IGHV3-25*04"} <= selected_v   # non-canonical IMGT 23/41 residues are kept
    (tyr,) = [row for row in SELECTION["omitted"] if row["label"] == "IGHV1-38-4*01"]
    assert tyr["failed_rules"][0]["detail"] == "junction anchor codon (IMGT 104) is TAT"


def test_igh_capability_target():
    target = SELECTION["capability_target"]
    assert target["required"] == ["observed_span_simulation"]
    assert target["targeted"] == ["source_annotated_junction", "canonical_productive"]
    assert set(target["required"]) | set(target["targeted"]) <= set(CAPABILITIES)
    basis = target["basis_for_target"]["canonical_productive"]
    assert "below source certification" in basis and "'assertion'" in basis and "'none'" in basis
    assert "canonical_productive_source_certified" in target["not_targeted"]
    assert "original_coding_ends" in target["not_targeted"]


def test_recorded_observations_read_off_the_raw_bytes():
    igh = _by_label("9606.IGH_VDJ/10", "airr")
    record = igh["IGHV1-58*01"]
    sequence = record["sequence"]
    assert (record["gene_start"], record["gene_end"], record["v_rs_start"]) == (143, 438, 324)
    assert sequence[438:445] == "CACAGTG"
    assert sequence[323:330] == "CACAGAA"
    (note,) = SELECTION["selected_record_observations"]
    assert note["label"] == "IGHV1-58*01" and note["observation"] == "v_rss_start_inside_coding_interval"
    assert (note["sequence_after_gene_end"], note["sequence_at_v_rs_start"]) == ("CACAGTG", "CACAGAA")

    igl = _by_label("9606.IGLambda_VJ/4", "airr_ex")["IGLV3-2*03"]
    coding = igl["coding_sequence"]
    assert igl["functional"] is True and igl["curational_tags"] == ["likely_full_length"]
    assert coding[250:271] == "TTACTGTCGCTCGGAACAGCA"
    assert (coding[254:257], coding[258:261], coding[272:275]) == ("TGT", "GCT", "TGC")
    assert 254 % 3 == 272 % 3 == 2 and 258 % 3 == 0

    observations = {row["label"]: row for row in EVIDENCE["record_observations"]}
    assert set(observations) == {"IGLV3-2*03", "IGHV1-58*01", "IGHV3-7*04"}
    assert observations["IGLV3-2*03"]["coding_sequence_250_271"] == coding[250:271]
    assert observations["IGHV1-58*01"]["sequence_after_gene_end"] == sequence[438:445]
    assert observations["IGHV1-58*01"]["sequence_at_v_rs_start"] == sequence[323:330]
    v37 = igh["IGHV3-7*04"]
    assert observations["IGHV3-7*04"]["sequence_after_gene_end"] == v37["sequence"][455:463] == "ACACAGTG"


def test_interpretation_evidence_names_the_binding_requirement_and_zero_coordinates():
    files = {entry["path"]: entry["lines"] for entry in EVIDENCE["producer_code"]["files"]}
    assert "treat 0 as missing" in files["api_v2/api.py"]["560-570"]
    assert {"1107-1108", "2085-2090"} <= set(files["ogrdb/sequence/sequence_routes.py"])
    v_landmark = EVIDENCE["conventions"]["v_landmark"]
    assert v_landmark["scope"].endswith("unsupported until a reviewed profile extension.")
    evidence = " ".join(v_landmark["evidence"])
    assert "sequence_routes.py:1107-1108" in evidence and "sequence_routes.py:2085-2090" in evidence
    origin = next(item for item in v_landmark["evidence"] if item.startswith("Alignment origin"))
    assert "310-312 are not used" in origin and "Cys 104" not in origin
    exceptions = v_landmark["imgt_23_41_exceptions_in_both_gapped_exports"]
    assert all(set(codons) <= {"23", "41"} for per_set in exceptions.values() for codons in per_set.values())
    for key in EXPECT["sets"]:
        short = key.replace("9606.", "").replace("/", ".")
        for export in ("gapped", "gapped_ex"):
            found = {}
            for label, aligned in _fasta(key, export).items():
                if "V" in label[:4]:
                    codons = {"23": aligned[66:69], "41": aligned[120:123]}
                    odd = {pos: c for pos, c in codons.items() if c not in ({"TGT", "TGC"} if pos == "23" else {"TGG"})}
                    if odd:
                        found[label] = odd
            assert found == exceptions[short]
    for key in EXPECT["sets"]:
        for export in ("airr", "airr_ex"):
            for record in _records(key, export):
                for delineation in record["v_gene_delineations"] or []:
                    assert all(value not in (0, None) for name, value in delineation.items()
                               if name.endswith(("_start", "_end")))
