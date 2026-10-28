"""The reviewed human IGH starter packaged with GenAIRR: its recipe, its generator and its installation.

``tests/data/reference_sources/human_igh_ogrdb_10_observed_v1.recipe.json`` is the reviewed recipe:
the candidate selection with a specific reason per allele.  ``tools/build_published_reference.py``
builds it from the pinned exports through the public source API and writes the dossier that is
packaged under ``GenAIRR/data/published_references/``; ``ReferenceBuild.bundled`` loads it by its
exact identifier.

Expected values come from the selection fixture, the pinned raw exports read with the standard
library, or literals pinned from the reviewed build; never from the code under test.
"""
from __future__ import annotations

from collections import Counter
import hashlib
import io
import json
import pathlib
import re
import zipfile

import pytest

from GenAIRR._artifact_io import canonical_json
from GenAIRR._reference_sources.recipe import Recipe
from tests.support.fixtures.ogrdb_inventory import IGH, RECORDS, SELECTION, raw, raw_records
from tests.support.fixtures.source_api import DATA, PROFILE, SOURCES

pytestmark = pytest.mark.behavior

STARTER_ID = "human_igh_ogrdb_10_observed_v1"
RECIPE_PATH = DATA / f"{STARTER_ID}.recipe.json"
# The lock of the pinned IGH release 10 exports and OGRDB statement pages recorded as local files,
# with the data terms they state (from the reviewed build).
LOCK_PIN = "sha256:5f905bfddb0fc15f0d8cb1ff357ba4779dcffcd4810f507ddfd97d103648f88b"
CHOICES = [row["label"] for row in SELECTION["record_choices"]]
REQUIRES = [{"capability": "observed_span_simulation", "evidence": "assertion"},
            {"capability": "source_annotated_junction", "evidence": "source"},
            {"capability": "canonical_productive", "evidence": "assertion"}]
AFFIRMED_UP_TO = re.compile(r"up to and including (?:nucleotide|base) (\d+)")
# Conserved IMGT codons other than the junction anchor, as 0-based gapped columns.
CONSERVED = ((23, 66, "C", "Cys"), (41, 120, "W", "Trp"))
CODONS = {a + b + c: aa for (a, b, c), aa in zip(
    ((a, b, c) for a in "TCAG" for b in "TCAG" for c in "TCAG"),
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")}
RESIDUES = {"C": "Cys", "W": "Trp", "F": "Phe"}
ORF_SCOPE = ("ORF-like record feature (this flag checks only the CAC heptamer core at the declared V RSS start and "
             "the IMGT 23, 41 and 104 codons; it is not a functionality classification): ")
# The IMGT human IGHV set kept in the repository, and the IMGT functionality classes it states
# (FASTA header field 4) that each selected V reason must disclose.
IMGT_IGHV = DATA.parent / "IGHV.fasta"
IMGT_IGHV_SHA256 = "sha256:e0756e5e1cfa6306289b335738a28339f3751cec2ab9f1de8a16a5ae8cb8f934"
IMGT_DISCLOSED = {"ORF": "ORF", "P": "a pseudogene"}
IMGT_CLAUSE = ("This allele is classed as {word} in the IMGT human IGHV set in the GenAIRR source repository "
               "(tests/data/IGHV.fasta, " + IMGT_IGHV_SHA256 + "; its IMGT release is not recorded), whose entry in "
               "that file has the same V-REGION sequence; the provider's functional field is true, it is selected "
               "under the starter's provider-field rule, and no runtime functionality is claimed.")


def sha(sequence: str) -> str:
    return "sha256:" + hashlib.sha256(sequence.encode("ascii")).hexdigest()


def raw_by_label(fmt) -> dict:
    return {record["label"]: record for record in raw_records(IGH, fmt)}


def gapped_by_label() -> dict:
    sequences = {}
    for chunk in raw(IGH, "gapped").decode("ascii").split(">")[1:]:
        header, _, body = chunk.partition("\n")
        sequences[header.strip()] = body.replace("\n", "")
    return sequences


def recipe_document() -> dict:
    return json.loads(RECIPE_PATH.read_bytes().decode("utf-8"))


def reasons() -> dict:
    return {entry["label"]: entry["reason"] for entry in recipe_document()["selection"]}


# ── the reviewed recipe ──────────────────────────────────────────────


def test_the_recipe_is_canonical_lf_json_of_a_valid_recipe():
    data = RECIPE_PATH.read_bytes()
    assert not data.startswith(b"\xef\xbb\xbf") and b"\r" not in data
    document = json.loads(data.decode("utf-8"))
    assert data == (canonical_json(document) + "\n").encode("utf-8")
    assert Recipe.from_json(data.decode("utf-8")).to_dict() == document


def test_the_recipe_names_the_starter_and_requires_its_targeted_capabilities(tmp_path):
    document = recipe_document()
    lock = load_tool().source(STARTER_ID, cache=tmp_path / "cache").lock      # the generator's lock
    assert document["source_lock_hash"] == LOCK_PIN == lock.commitment
    assert document["adapter"] == {"id": "ogrdb_germline_set", "version": 1,
                                   "coordinate_profile": {"id": PROFILE, "version": 1}}
    assert document["projection"] == {"id": "observed_span", "version": 1}
    assert document["reference"] == {"name": STARTER_ID, "version": "1", "build_date": "2026-10-27",
                                     "taxon": {"scientific_name": "Homo sapiens", "taxonomy_id": 9606},
                                     "locus": "IGH", "chain_type": "vdj"}
    assert document["requires"] == REQUIRES
    assert document["defaults_profile"] == "genairr.unfitted_defaults.v1"


def test_the_recipe_selects_exactly_the_reviewed_selection_from_its_own_airr_records():
    records = raw_records(IGH, "airr")
    index = {record["label"]: i for i, record in enumerate(records)}
    expected = [{"segment": segment, "label": label,
                 "record": {"asset": "airr", "pointer": f"{RECORDS}/{index[label]}"},
                 "sequence_sha256": sha(records[index[label]]["coding_sequence"])}
                for segment in ("V", "D", "J") for label in SELECTION["selected"][segment]]
    selection = [{k: v for k, v in entry.items() if k != "reason"} for entry in recipe_document()["selection"]]
    assert selection == expected
    assert Counter(entry["segment"] for entry in selection) == {"V": 190, "D": 33, "J": 7}
    genes = {segment: {entry["label"].split("*")[0] for entry in selection if entry["segment"] == segment}
             for segment in "VDJ"}
    assert {segment: len(names) for segment, names in genes.items()} == SELECTION["distinct_genes_selected"]
    omitted = {row["label"] for row in SELECTION["omitted"]}
    assert len(omitted) == 19 and not omitted & {entry["label"] for entry in selection}


def test_the_recipe_selects_the_airr_record_of_each_record_choice_explicitly():
    airr, entries = raw_by_label("airr"), recipe_document()
    decisions = entries["decisions"]
    assert [decision["sequence_sha256"] for decision in decisions] == [
        sha(airr[label]["coding_sequence"]) for label in CHOICES]
    for label, decision in zip(CHOICES, decisions):
        assert decision["action"] == "select_alternative" and decision["claims"] == [] == decision["evidence"]
        assert len(decision["findings"]) == 1 and decision["reason"] == reasons()[label]


def test_a_record_choice_reason_names_the_chosen_export_its_end_and_the_trim_effect():
    airr, airr_ex, gapped, given = raw_by_label("airr"), raw_by_label("airr_ex"), gapped_by_label(), reasons()
    for label in CHOICES:
        record, text = airr[label], given[label]
        coding = record["coding_sequence"]
        extension = airr_ex[label]["coding_sequence"][len(coding):]
        assert airr_ex[label]["coding_sequence"] == coding + extension
        assert len(gapped[label].rstrip(".")) == 319
        assert text.startswith("Record choice: the 'airr' record is selected;")
        assert f"Its {len(coding)}-nt coding sequence ends at gene_end, at IMGT gapped position 319" in text
        assert f"extends it at the 3' end by {len(extension)} nt ({extension})" in text
        assert f"stops {len(extension)} nt short of the 'airr_ex' end, as a {len(extension)}-nt 3' trim" in text
        assert "not adopted" in text and "does not establish which end is biologically correct" in text
        if record["v_rs_start"] is None:
            assert record["sequence"][record["gene_end"]:] == extension
            assert "of its record sequence (it declares no V RSS start)" in text
        else:
            start = record["v_rs_start"]
            assert label == "IGHV3-7*04" and start == record["gene_end"] + 1 + len(extension)
            assert (f"directly before its declared V RSS start at record position {start} "
                    f"({record['sequence'][start - 1:start + 6]})") in text


def test_identical_dna_reasons_name_the_other_identity_and_keep_each_ones_own_evidence():
    airr, airr_ex, given = raw_by_label("airr"), raw_by_label("airr_ex"), reasons()
    groups = SELECTION["identical_coding_sequence_groups_selected"]
    assert len(groups) == 9
    for group in groups:
        for label, other in (group, group[::-1]):
            assert airr[label]["coding_sequence"] == airr[other]["coding_sequence"]
            text = given[label]
            assert f"Its coding sequence is identical to that of {other}, which is also selected" in text
            assert "separate identities" in text and f"nothing is copied from {other}" in text
            assert ("It has no 'airr_ex' record of its own" in text) == (label not in airr_ex)
    named = [label for label, text in given.items() if "is identical to that of" in text]
    assert len(named) == 18


def test_every_other_reason_states_the_starter_rules_and_each_specific_observation():
    airr, gapped, given = raw_by_label("airr"), gapped_by_label(), reasons()
    for segment in ("V", "D", "J"):
        for label in SELECTION["selected"][segment]:
            text, record = given[label], airr[label]
            assert record["functional"] is True and "likely_truncated" not in (record["curational_tags"] or [])
            if label in CHOICES:
                continue
            assert text.startswith("Selected from its own 'airr' record under the starter rules: provider functional "
                                   "field true")
            if segment == "V":
                start = record["v_gene_delineations"][0]["cdr3_start"] - 4
                codon = record["coding_sequence"][start:start + 3]
                assert CODONS[codon] == "C" and f"{codon} (Cys) at IMGT 104" in text
            elif segment == "J":
                anchor = record["j_cdr3_end"] - 1
                assert record["coding_sequence"][anchor:anchor + 3] == "TGG" and "TGG (Trp) at j_cdr3_end" in text
            else:
                assert "A D allele carries no anchor or coding phase" in text
    # Alleles whose sequence runs past the span their curation note affirms.
    past = {}
    for label in SELECTION["selected"]["V"]:
        match = AFFIRMED_UP_TO.search(airr[label]["curation"] or "")
        end = len(gapped[label].rstrip("."))
        if match and end > int(match.group(1)):
            past[label] = (int(match.group(1)), end)
    assert len(past) == 9
    for label, (column, end) in past.items():
        assert (f"Its curation note affirms the sequence up to IMGT gapped position {column}; its coding sequence runs "
                f"{end - column} nt past it, to {end}, and is used as supplied.") in given[label]
    # ORF-like alleles under the provider's functional=true: a declared heptamer without CAC, a
    # non-conserved IMGT 23 or 41 residue.
    orf = {}
    for label in SELECTION["selected"]["V"]:
        record, sequence = airr[label], gapped[label]
        start = record["v_rs_start"]
        adjacent = start is not None and start == record["gene_end"] + 1
        if adjacent and not record["sequence"][start - 1:].startswith("CAC"):
            orf[label] = f"the record bases at its declared V RSS start read {record['sequence'][start - 1:start + 6]}"
        for position, column, conserved, name in CONSERVED:
            codon = sequence[column:column + 3]
            if CODONS[codon] != conserved:
                orf[label] = (f"its IMGT {position} codon {codon} encodes {RESIDUES[CODONS[codon]]}, not the conserved "
                              f"{name}")
    assert sorted(orf) == ["IGHV3-16*02", "IGHV3-16*03", "IGHV3-16*04", "IGHV3-20*02", "IGHV3-25*04",
                           "IGHV8-51-1*04", "IGHV8-51-1*05"]
    for label, observation in orf.items():
        assert ORF_SCOPE + observation in given[label]
    assert sum("ORF-like" in text for text in given.values()) == 7
    # The one selected V grouped with the record choices by its end, whose note does not call
    # further bases likely.
    assert ("Its coding sequence ends at IMGT gapped position 320, the last position of the span its curation note "
            "affirms 'based on the convincing evidence of the sequence end'; unlike the notes of the 14 record "
            "choices, it does not call further bases likely.") in given["IGHV3-21*07"]
    assert "based on the convincing evidence of the sequence end" in airr["IGHV3-21*07"]["curation"]
    (observation,) = [row for row in SELECTION["selected_record_observations"]
                      if row["observation"] == "v_rss_start_inside_coding_interval"]
    assert (f"Its declared v_rs_start ({observation['v_rs_start']}) lies inside its coding interval, while the record "
            f"bases after gene_end read {observation['sequence_after_gene_end']}") in given[observation["label"]]


def imgt_classes(path: pathlib.Path = IMGT_IGHV) -> dict:
    """``{allele: (IMGT class, ungapped V-REGION)}`` read from the IMGT FASTA headers (field 4)."""
    data = path.read_bytes()
    found = "sha256:" + hashlib.sha256(data).hexdigest()
    if found != IMGT_IGHV_SHA256:
        lf = "sha256:" + hashlib.sha256(data.replace(b"\r\n", b"\n")).hexdigest()
        if b"\r\n" in data and lf == IMGT_IGHV_SHA256:
            pytest.fail(f"{path.name} was checked out with CRLF line endings ({found}); its bytes are pinned, so "
                        ".gitattributes must keep tests/data/IGHV.fasta byte-exact with a -text rule")
        pytest.fail(f"{path.name} is {found}, not the pinned {IMGT_IGHV_SHA256}")
    entries = {}
    for chunk in data.decode("ascii").split(">")[1:]:
        header, _, body = chunk.partition("\n")
        fields = header.split("|")
        assert fields[1] not in entries
        entries[fields[1]] = (fields[3], body.replace("\n", "").replace(".", "").upper())
    return entries


def test_a_crlf_checkout_of_the_imgt_set_is_named_as_the_cause(tmp_path):
    copy = tmp_path / "IGHV.fasta"
    copy.write_bytes(IMGT_IGHV.read_bytes().replace(b"\n", b"\r\n"))
    with pytest.raises(pytest.fail.Exception, match="checked out with CRLF line endings"):
        imgt_classes(copy)


def test_each_v_that_imgt_classes_as_orf_or_pseudogene_says_so():
    classes, airr, given = imgt_classes(), raw_by_label("airr"), reasons()
    selected = SELECTION["selected"]["V"]
    assert Counter(classes[label][0] if label in classes else "absent" for label in selected) == {
        "F": 153, "(F)": 14, "ORF": 12, "P": 5, "absent": 6}
    expected = {label for label in selected if label in classes and classes[label][0] in IMGT_DISCLOSED}
    assert {label for label, text in given.items() if "is classed as" in text} == expected          # both ways
    assert len(expected) == 17
    for label in expected:
        imgt_class, region = classes[label]
        assert region == airr[label]["coding_sequence"] and airr[label]["functional"] is True
        assert IMGT_CLAUSE.format(word=IMGT_DISCLOSED[imgt_class]) in given[label], label
    assert classes["IGHV8-51-1*05"][0] == "P"


def test_leader_observations_are_stated_exactly_where_the_record_fields_show_them():
    airr, given = raw_by_label("airr"), reasons()
    starts, frames, donors = {}, {}, {}
    for label in SELECTION["selected"]["V"]:
        record = airr[label]
        sequence, first, end = record["sequence"], record["leader_1_start"], record["leader_1_end"]
        if first is None or record["leader_2_start"] is None:
            continue
        if sequence[first - 1:first + 2] != "ATG":
            starts[label] = (first, sequence[first - 1:first + 2])
        lengths = (end - first + 1, record["leader_2_end"] - record["leader_2_start"] + 1)
        if sum(lengths) % 3:
            frames[label] = lengths
        if sequence[end:end + 2] not in ("GT", "GC"):
            donors[label] = sequence[end:end + 2]
    assert starts == {"IGHV8-51-1*05": (1, "CTT")}
    assert frames == {"IGHV3-38-3*01": (48, 11)} and donors == {"IGHV3-38-3*01": "GA"}
    assert ("Its declared leader_1 starts at record position 1 with CTT, not ATG; the record begins at that base, so "
            "an ATG upstream of it is not ruled out.") in given["IGHV8-51-1*05"]
    assert ("Its declared leader exons total 59 nt (leader_1 48 nt, leader_2 11 nt), not a whole number of codons, "
            "and the two record bases after leader_1_end read GA, not the GT of a canonical splice donor.") in \
        given["IGHV3-38-3*01"]
    assert sorted(label for label, text in given.items() if "leader" in text) == ["IGHV3-38-3*01", "IGHV8-51-1*05"]


def test_every_reason_is_one_line_of_ascii_text_that_never_overstates_functionality():
    for text in list(reasons().values()) + [d["reason"] for d in recipe_document()["decisions"]]:
        assert text.isascii() and "\n" not in text
        assert "functional allele" not in text.lower()


# ── the generator ────────────────────────────────────────────────────

REPO = DATA.parent.parent.parent
TOOL_PATH = REPO / "tools" / "build_published_reference.py"


def load_tool():
    """The generator script, loaded from the checkout (it is not part of the package)."""
    import importlib.util
    spec = importlib.util.spec_from_file_location("build_published_reference", TOOL_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def refuse_network(patch) -> None:
    """Any connection attempt, name lookup or fetch transport fails the test."""
    import http.client
    import socket
    from GenAIRR._reference_sources import api
    from GenAIRR._reference_sources.acquisition import HttpsTransport

    def refuse(*args, **kwargs):
        pytest.fail("the generator tried to reach the network")

    for owner, name in ((socket, "create_connection"), (socket, "getaddrinfo"), (socket.socket, "connect"),
                        (http.client.HTTPSConnection, "connect"), (http.client.HTTPConnection, "connect"),
                        (HttpsTransport, "get"), (api, "_transport")):
        patch.setattr(owner, name, refuse)


# Every entry of a packaged archive: stored, at a fixed time, POSIX regular file 0644, nothing extra.
ZIP_ENTRY = {"compress_type": zipfile.ZIP_STORED, "date_time": (1980, 1, 1, 0, 0, 0), "create_system": 3,
             "create_version": 20, "extract_version": 20, "external_attr": 0o100644 << 16, "internal_attr": 0,
             "flag_bits": 0, "extra": b"", "comment": b""}


def archive_members(data: bytes) -> dict:
    """The entries of a zip archive, in archive order (name -> bytes), read with the standard library."""
    with zipfile.ZipFile(io.BytesIO(data)) as archive:
        return {entry.filename: archive.read(entry) for entry in archive.infolist()}


def extracted(data: bytes, directory: pathlib.Path) -> pathlib.Path:
    for name, member in archive_members(data).items():
        path = directory.joinpath(*name.split("/"))
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(member)
    return directory


@pytest.fixture(scope="module")
def generated(tmp_path_factory):
    """The starter archive written by the generator's command line into a new file, offline."""
    import contextlib
    target = tmp_path_factory.mktemp("generated") / f"{STARTER_ID}.zip"
    out = io.StringIO()
    with pytest.MonkeyPatch.context() as patch, contextlib.redirect_stdout(out):
        refuse_network(patch)
        code = load_tool().main(["--output", str(target)])
    return code, target, out.getvalue()


def test_the_generator_writes_one_deterministic_stored_zip_of_the_dossier_members_offline(generated, tmp_path):
    from GenAIRR._reference_sources.digests import manifest_hash
    from GenAIRR.reference_sources import ReferenceBuild
    code, target, printed = generated
    assert code == 0
    data = target.read_bytes()
    with zipfile.ZipFile(io.BytesIO(data)) as archive:
        entries = archive.infolist()
        assert archive.comment == b""
    for entry_ in entries:
        assert {field: getattr(entry_, field) for field in ZIP_ENTRY} == ZIP_ENTRY, entry_.filename
        assert entry_.file_size == entry_.compress_size and not entry_.filename.endswith("/")
    members = archive_members(data)
    assert list(members) == sorted(members)                         # sorted by member name, files only
    checked = ReferenceBuild.read_dossier(extracted(data, tmp_path / "dossier"))
    assert list(members) == sorted([member.name for member in checked.manifest.members] + ["manifest.json"])
    reference = checked.reference
    for line in (f"manifest {manifest_hash(checked.manifest)}", f"runtime execution {reference.execution_hash}",
                 f"runtime catalogue {reference.catalogue_hash}", f"runtime models {reference.models_hash}",
                 f"runtime artifact {reference.artifact_hash}", f"builder GenAIRR {checked.manifest.builder.version}",
                 f"archive sha256:{hashlib.sha256(data).hexdigest()} ({len(data)} bytes)",
                 f"lock {checked.lock.commitment}"):
        assert line in printed.splitlines(), line
    hint = printed.splitlines()[-1]                               # every pin a regeneration may move
    for pin in ("src/GenAIRR/_reference_sources/bundled.py", "MANIFEST_PIN", "ARCHIVE_PIN", "RUNTIME_PINS",
                "LOCK_PIN", "tests/components/reference/test_published_reference_resources.py"):
        assert pin in hint, pin
    assert members["recipe.json"] == RECIPE_PATH.read_bytes()
    assert load_tool().archive(dict(reversed(list(members.items())))) == data       # the same bytes again


# OGRDB's data terms as the fixture manifest records them, and the statement pages kept as fixtures.
TERMS = {"status": "recorded", "assets": ["licence_statement", "citation_statement"],
         "licence": SOURCES["data_licence"]["licence"],
         "citation": SOURCES["citation"]["requested_by_provider"] + " doi:" + SOURCES["citation"]["doi"],
         "reason": None}
STATEMENT_ASSETS = {
    "licence_statement": (SOURCES["data_licence"]["statement_url"],
                          DATA / "ogrdb" / "statements" / "licensing_statement.html"),
    "citation_statement": (SOURCES["citation"]["statement_url"], DATA / "ogrdb" / "statements" / "help.html"),
}


def assert_terms_recorded(lock, blob=None) -> None:
    """The lock records OGRDB's terms from the two statement pages (and ``blob(sha256)`` holds their bytes)."""
    assert TERMS["licence"] == "CC0 1.0 Universal" and TERMS["citation"].startswith("The current landscape")
    assert TERMS["citation"].endswith("Nucleic Acids Research, November 2025. doi:10.1093/nar/gkaf1094")
    assert lock.terms.to_dict() == TERMS
    statements = {asset.format: asset for asset in lock.assets if asset.role == "data_terms"}
    assert sorted(statements) == sorted(STATEMENT_ASSETS)
    for fmt, (url, path) in STATEMENT_ASSETS.items():
        data, asset = path.read_bytes(), statements[fmt]
        assert (asset.declared_url, asset.origin_kind, asset.sha256, asset.bytes) == (
            url, "local_file", "sha256:" + hashlib.sha256(data).hexdigest(), len(data))
        if blob is not None:
            assert blob(asset.sha256) == data


def test_the_generator_records_the_providers_data_terms_and_inspection_shows_them(tmp_path):
    source = load_tool().source(STARTER_ID, cache=tmp_path / "cache")
    assert_terms_recorded(source.lock)
    text = str(source.inspect())
    assert f"Data terms: recorded (licence_statement, citation_statement): {TERMS['licence']}" in text.splitlines()
    assert f"Citation: {TERMS['citation']}" in text.splitlines()


def test_a_written_archive_that_does_not_read_back_unchanged_is_removed(tmp_path, monkeypatch):
    from GenAIRR.reference_sources import ReferenceSourceError
    tool = load_tool()
    members = {"manifest.json": b"{}\n", "inputs/v.fasta": b">a\nACGT\n"}
    # The loader refuses it (this manifest is not a manifest): the archive is removed.
    target = tmp_path / "refused.zip"
    with pytest.raises(ReferenceSourceError):
        tool.write_archive(members, target)
    assert not target.exists()

    # The loader extracts other bytes than were written: the archive is removed.
    def altered(resource, directory, where, check):
        directory.mkdir()
        (directory / "manifest.json").write_bytes(b"[]\n")

    monkeypatch.setattr(tool.bundled, "_extract", altered)
    target = tmp_path / "altered.zip"
    with pytest.raises(tool.GeneratorError, match="does not extract"):
        tool.write_archive(members, target)
    assert not target.exists() and sorted(p.name for p in tmp_path.iterdir()) == []


def test_the_generator_writes_only_a_new_zip_file(tmp_path, capsys):
    existing = tmp_path / "kept.zip"
    existing.write_bytes(b"kept\n")
    assert load_tool().main(["--output", str(existing)]) == 2
    assert "already exists" in capsys.readouterr().err and existing.read_bytes() == b"kept\n"
    assert load_tool().main(["--output", str(tmp_path / "missing-parent" / "starter.zip")]) == 2
    assert load_tool().main(["--output", str(tmp_path / "starter")]) == 2
    assert "a .zip file" in capsys.readouterr().err
    assert sorted(p.name for p in tmp_path.iterdir()) == ["kept.zip"]


# ── the packaged starter ─────────────────────────────────────────────

# Pinned from the reviewed build: the typed manifest hash of the packaged dossier (it covers
# the GenAIRR version that generated it) and the runtime hashes of its reference.
MANIFEST_PIN = "sha256:9ab6ba63490c06af94098152af388141f8f946bd6eab3c2ef86f8e5ddb512550"
RUNTIME_PINS = {"execution": "sha256:bacd155a623189d5acd7c1720a9decbde4d4c5b6e6f6c5581e6ce263683df6dc",
                "catalogue": "sha256:fc8d35b4682385b1ef10c1405c7427ca2d7551ee626f4257f7bee4fadafb5266",
                "models": "sha256:e96b1404c64d36a55c0457163ef3a4aa706e1bc5114914147715cd0f964dffd8"}
# The raw sha256 of the packaged archive (it holds the manifest, so it also changes with the version).
ARCHIVE_PIN = "sha256:52df85b32d417c67f882cce18113237eb4889f54962695853dd398bc8a3e33d8"
RECORD_OUTCOMES = {"selected": 230, "excluded": 259, "rejected": 0}
PACKAGE, RESOURCE = "GenAIRR.data", f"published_references/{STARTER_ID}.zip"
DOSSIER_DIRECTORIES = ["cache", "cache/blobs", "cache/blobs/sha256", "inputs"]


@pytest.fixture(scope="module")
def starter():
    from GenAIRR.reference_sources import ReferenceBuild
    return ReferenceBuild.bundled(STARTER_ID)


def runtime_hashes(reference) -> dict:
    return {"execution": reference.execution_hash, "catalogue": reference.catalogue_hash,
            "models": reference.models_hash}


def packaged_archive() -> bytes:
    import importlib.resources
    return importlib.resources.files(PACKAGE).joinpath("published_references").joinpath(f"{STARTER_ID}.zip") \
        .read_bytes()


def packaged_files() -> tuple:
    """Every file of the packaged dossier (relative POSIX name -> bytes), and its directories."""
    from GenAIRR._reference_sources import bundled
    with bundled.dossier_directory(bundled.entry(STARTER_ID)) as directory:
        files = {p.relative_to(directory).as_posix(): p.read_bytes() for p in sorted(directory.rglob("*"))
                 if p.is_file()}
        folders = sorted(p.relative_to(directory).as_posix() for p in directory.rglob("*") if p.is_dir())
    return files, folders


def test_the_registry_names_exactly_the_starter():
    from GenAIRR._reference_sources import bundled
    from GenAIRR.reference_sources import ReferenceBuild, ReferenceSourceError
    assert dict(bundled.REGISTRY) == {STARTER_ID: bundled.BundledDossier(PACKAGE, RESOURCE, MANIFEST_PIN)}
    for near in ("human_igh", STARTER_ID.upper(), STARTER_ID + " ", "human_igh_ogrdb_10_observed"):
        with pytest.raises(ReferenceSourceError) as info:
            ReferenceBuild.bundled(near)
        assert info.value.kind == "missing_asset" and f"available: {STARTER_ID}" in info.value.message


def test_the_bundled_starter_is_the_reviewed_build_with_its_pinned_hashes(starter):
    import GenAIRR
    from GenAIRR._reference_sources.digests import manifest_hash
    from GenAIRR.reference_sources import ReferenceBuild
    assert manifest_hash(starter.manifest) == MANIFEST_PIN
    assert runtime_hashes(starter.reference) == RUNTIME_PINS
    runtime = starter.report.runtime
    assert {"execution": runtime.execution_hash, "catalogue": runtime.catalogue_hash,
            "models": runtime.models_hash} == RUNTIME_PINS
    outcomes = Counter(outcome.outcome for outcome in starter.report.records)
    assert {name: outcomes[name] for name in RECORD_OUTCOMES} == RECORD_OUTCOMES
    assert len(starter.report.records) == sum(RECORD_OUTCOMES.values())
    assert starter.recipe.to_dict() == recipe_document()
    assert starter.checked_with == GenAIRR.__version__
    assert ReferenceBuild.bundled(STARTER_ID) is starter                 # checked once per process
    alleles = starter.reference.alleles
    assert [(a.segment, a.label) for a in alleles] == [
        (segment, label) for segment in ("V", "D", "J") for label in SELECTION["selected"][segment]]
    differs = {finding.observations.get()["label"]: finding.id for finding in starter.report.inventory.findings
               if finding.code == "export_coding_sequence_differs"}
    assert [decision.findings for decision in starter.recipe.decisions] == [(differs[label],) for label in CHOICES]


def test_the_bundled_starter_records_the_providers_data_terms(starter):
    files, _ = packaged_files()
    assert_terms_recorded(starter.lock, lambda sha256: files["cache/blobs/sha256/" + sha256.split(":", 1)[1]])
    lines = str(starter.report).splitlines()
    assert f"  data terms: recorded (licence_statement, citation_statement): {TERMS['licence']}" in lines
    assert f"  citation: {TERMS['citation']}" in lines


def test_the_bundled_starter_states_its_capabilities_truthfully(starter):
    from GenAIRR._reference_sources.projection import PROJECTION_LOSSES
    counts = SELECTION["counts"]
    blockers = {"original_coding_ends": counts["V"] + counts["J"] + 2 * counts["D"],   # V 3', J 5', both D ends
                "complete_reference_span": 2 * (counts["V"] + counts["D"] + counts["J"])}
    for capabilities in (starter.reference.capabilities(), starter.report.capabilities):
        status = capabilities.status
        assert status("observed_span_simulation").available
        assert (status("source_annotated_junction").available, status("source_annotated_junction").evidence) == (
            True, "source")
        assert (status("canonical_productive").available, status("canonical_productive").evidence) == (
            True, "assertion")
        assert capabilities.phase_support == "anchor_derived"
        for name, count in blockers.items():
            assert not status(name).available and len(status(name).unmet) == count, name
    assert blockers == {"original_coding_ends": 263, "complete_reference_span": 460}
    assert [r.to_dict() for r in starter.report.requirements] == [dict(r, satisfied=True) for r in REQUIRES]
    text = " ".join(str(starter.report).split())
    assert [code for code, _ in PROJECTION_LOSSES] == ["v_subregions", "gapped_sequences", "affirmed_spans",
                                                       "genomic_context"]
    for code, loss in PROJECTION_LOSSES:
        assert f"{code}: {' '.join(loss.split())}" in text
    assert "functional alleles" not in text.lower()


def test_the_bundled_starter_simulates_bound_to_its_execution_hash(starter):
    from GenAIRR import Experiment
    result = Experiment.on(starter.reference).recombine().run_records(n=60, seed=11, validate_records=True)
    assert result.reference_content_hash == starter.reference.execution_hash == RUNTIME_PINS["execution"]
    assert len(result) == 60
    for record in result:
        assert set(record["v_call"].split(",")) <= set(SELECTION["selected"]["V"])
        assert set(record["j_call"].split(",")) <= set(SELECTION["selected"]["J"])


def test_rebuilding_from_the_pinned_bytes_and_the_recipe_reproduces_the_packaged_starter(generated, starter,
                                                                                         tmp_path):
    import GenAIRR
    from GenAIRR.reference_sources import ReferenceBuild
    code, target, _ = generated
    assert code == 0
    rebuilt = ReferenceBuild.read_dossier(extracted(target.read_bytes(), tmp_path / "rebuilt"))
    assert runtime_hashes(rebuilt.reference) == runtime_hashes(starter.reference) == RUNTIME_PINS
    ours, packaged = rebuilt.report.to_dict(), starter.report.to_dict()
    for report in (ours, packaged):              # the generator, and the artifact hash that names it
        report.pop("generator")
        report["runtime"].pop("artifact_hash")
    assert ours == packaged
    builder = starter.manifest.builder.version
    if GenAIRR.__version__ != builder:
        pytest.skip(f"the packaged dossier was generated by GenAIRR {builder} and this is GenAIRR "
                    f"{GenAIRR.__version__}: its members are byte-identical only to a dossier generated by the "
                    "same version; regenerate it with tools/build_published_reference.py")
    assert target.read_bytes() == packaged_archive()


def test_the_packaged_dossier_is_one_archive_of_only_its_members(starter):
    import importlib.resources
    folder = importlib.resources.files(PACKAGE).joinpath("published_references")
    assert sorted(item.name for item in folder.iterdir()) == [f"{STARTER_ID}.zip"]  # no marker, no unpacked copy
    data = packaged_archive()
    assert "sha256:" + hashlib.sha256(data).hexdigest() == ARCHIVE_PIN
    files, folders = packaged_files()
    assert list(archive_members(data)) == sorted(files) == sorted(
        [member.name for member in starter.manifest.members] + ["manifest.json"])
    assert folders == DOSSIER_DIRECTORIES                    # no receipts, no __pycache__, no package marker
    assert files["recipe.json"] == RECIPE_PATH.read_bytes()


def test_the_packaged_archive_is_included_in_the_wheel_and_the_sdist():
    import glob
    tomllib = pytest.importorskip("tomllib", reason="reading pyproject.toml needs Python 3.11 or later")
    source = REPO / "src" / "GenAIRR" / "data" / "published_references" / f"{STARTER_ID}.zip"
    if not source.is_file():
        pytest.skip("the source checkout is not available")
    members = {source.resolve()}
    includes = tomllib.loads((REPO / "pyproject.toml").read_bytes().decode("utf-8"))["tool"]["maturin"]["include"]
    for target in ("wheel", "sdist"):
        matched = set()
        for rule in includes:
            formats = [rule["format"]] if isinstance(rule["format"], str) else rule["format"]
            if target in formats:
                found = glob.glob(str(REPO / rule["path"]), recursive=True)
                matched |= {pathlib.Path(path).resolve() for path in found}
        assert members <= matched, f"{target}: {sorted(str(p) for p in members - matched)}"
        if target == "wheel":                                   # no test fixtures ship in the wheel
            assert not [path for path in matched if (REPO / "tests") in path.parents]


# The built-in human_igh alias as it was before the starter was packaged.
HUMAN_IGH_CARTRIDGE = "HUMAN_IGH_OGRDB"
HUMAN_IGH_POOLS = (198, 33, 7)                                         # V, D, J alleles
HUMAN_IGH_CONTENT_HASH = "sha256:a98f4dad00768fcb7cb28164f1646b538dbb42da51a19d618eeff01b48e05b66"


def test_loading_the_starter_leaves_the_human_igh_alias_unchanged(starter):
    import GenAIRR
    from GenAIRR import Experiment
    from GenAIRR._refdata_resolver import _CONFIG_ALIASES
    from GenAIRR.reference_sources import ReferenceBuild
    assert ReferenceBuild.bundled(STARTER_ID) is starter
    assert sorted(alias for alias, name in _CONFIG_ALIASES.items() if name == HUMAN_IGH_CARTRIDGE) == [
        "human_igh", "human_igh_ogrdb"]
    assert not [alias for alias in _CONFIG_ALIASES if "observed" in alias or "published" in alias]
    with pytest.raises(ValueError, match="Unknown config name"):
        Experiment.on(STARTER_ID)
    alias = Experiment.on("human_igh")
    pools = (alias.refdata.v_pool_size(), alias.refdata.d_pool_size(), alias.refdata.j_pool_size())
    assert pools == HUMAN_IGH_POOLS
    first = alias.recombine().run_records(n=20, seed=3)
    second = Experiment.on(getattr(GenAIRR, HUMAN_IGH_CARTRIDGE)).recombine().run_records(n=20, seed=3)
    assert first.reference_content_hash == second.reference_content_hash == HUMAN_IGH_CONTENT_HASH
    assert HUMAN_IGH_CONTENT_HASH != starter.reference.execution_hash
    assert [record["sequence"] for record in first] == [record["sequence"] for record in second]


def declared_version(text: str) -> str:
    """``[project].version`` of a pyproject.toml text, read by the generator's own parser so that the
    release check and the generator follow one rule; a refusal fails the test with its message."""
    tool = load_tool()
    try:
        return tool.declared_version(text)
    except tool.GeneratorError as exc:
        pytest.fail(str(exc))


@pytest.mark.parametrize("text, expected", [
    ('[project]\nname = "GenAIRR"\nversion = "2.3.1"\n', "2.3.1"),
    ('[project]\r\nname = "GenAIRR"\r\nversion = "2.3.1"\r\n[tool.maturin]\r\n', "2.3.1"),          # a CRLF checkout
    ('[project]\nversion = "2.4.0"  # set by the release workflow\n[tool.maturin]\nversion = "9"\n', "2.4.0"),
    ("[project]\nversion = '2.4.1'\n", "2.4.1"),
])
def test_the_declared_version_is_read_from_the_project_table(text, expected):
    assert declared_version(text) == expected


@pytest.mark.parametrize("text, message", [
    ('[project]\nname = "GenAIRR"\ndynamic = ["version"]\n', "dynamic"),
    ('[tool.maturin]\nversion = "2.3.1"\n', "no [project] table"),
    ('[project]\nname = "GenAIRR"\n', "no version"),
])
def test_an_undeclared_version_fails_with_a_clear_message(text, message):
    with pytest.raises(pytest.fail.Exception, match=re.escape(message)):
        declared_version(text)


def test_a_crlf_checkout_of_pyproject_declares_the_same_version(tmp_path):
    pyproject = REPO / "pyproject.toml"
    if not pyproject.is_file():
        pytest.skip("no checkout pyproject.toml")
    text = pyproject.read_bytes().decode("utf-8").replace("\r\n", "\n")
    copy = tmp_path / "pyproject.toml"
    copy.write_bytes(text.replace("\n", "\r\n").encode("utf-8"))
    assert declared_version(copy.read_bytes().decode("utf-8")) == declared_version(text)


def test_the_packaged_dossier_was_generated_by_the_version_this_checkout_declares(starter):
    """A version change without regenerating the packaged archive fails here (checkout only)."""
    pyproject = REPO / "pyproject.toml"
    if not pyproject.is_file():
        pytest.skip("no checkout pyproject.toml: the installed version is not declared here")
    declared = declared_version(pyproject.read_bytes().decode("utf-8"))
    assert starter.manifest.builder.version == declared, (
        f"the packaged dossier was generated by GenAIRR {starter.manifest.builder.version}, but pyproject.toml "
        f"declares {declared}: regenerate it with tools/build_published_reference.py and pin it again")


def test_the_packaged_archive_and_its_sources_are_checked_out_byte_for_byte_on_every_platform():
    """Their digests cover exact bytes, so no checkout may convert line endings; the archive is binary."""
    attributes = REPO / ".gitattributes"
    if not attributes.is_file():
        pytest.skip("the source checkout is not available")
    rules = [line.split() for line in attributes.read_bytes().decode("utf-8").splitlines()
             if line.strip() and not line.startswith("#")]
    assert ["src/GenAIRR/data/published_references/*.zip", "-text", "-diff"] in rules
    assert ["tests/data/reference_sources/**", "-text"] in rules       # exports, statements and the recipe
    assert ["tests/data/IGHV.fasta", "-text"] in rules                 # the IMGT set whose digest is pinned


# ── installation conditions ──────────────────────────────────────────

BLOCKED = ("pandas", "numpy", "scipy", "graphviz", "tqdm", "fastmcp")


def network_guard() -> str:
    """The offline guard of the dossier tests: an audit hook and patched socket functions that
    count and refuse every connection, datagram and name lookup, with positive controls."""
    from tests.support.fixtures.source_dossier import OFFLINE_SCRIPT
    guard = OFFLINE_SCRIPT[:OFFLINE_SCRIPT.index("from GenAIRR import Experiment")]
    assert "sys.addaudithook(_audit)" in guard and "CONTROL_ATTEMPTS" in guard
    return guard


FRESH_PROCESS = r'''
from GenAIRR import Experiment
from GenAIRR._reference_sources.digests import manifest_hash
from GenAIRR.reference_sources import ReferenceBuild

build = ReferenceBuild.bundled(sys.argv[1])
reference = build.reference
result = Experiment.on(reference).recombine().run_records(n=40, seed=5, validate_records=True)
print(json.dumps({
    "manifest_hash": manifest_hash(build.manifest),
    "execution_hash": reference.execution_hash,
    "reference_content_hash": result.reference_content_hash,
    "records": len(result),
    "again": ReferenceBuild.bundled(sys.argv[1]) is build,
    "loaded_optional": sorted(name for name in sys.argv[2].split(",") if sys.modules.get(name) is not None),
    "network_attempts": ATTEMPTS,
    "controls": CONTROLS,
    "control_attempts": CONTROL_ATTEMPTS,
}))
'''


def test_the_starter_loads_and_simulates_in_a_fresh_process_without_optional_packages_or_network(tmp_path):
    import os
    import subprocess
    import sys
    blocker = "import sys\nfor _name in %r:\n    sys.modules[_name] = None\n" % (BLOCKED,)
    script = tmp_path / "fresh_starter.py"
    script.write_bytes((blocker + network_guard() + FRESH_PROCESS).encode("utf-8"))
    home = tmp_path / "home"
    home.mkdir()
    environment = dict(os.environ, HOME=str(home), USERPROFILE=str(home), XDG_CACHE_HOME=str(home / ".cache"),
                       APPDATA=str(home), LOCALAPPDATA=str(home))
    process = subprocess.run([sys.executable, "-I", str(script), STARTER_ID, ",".join(BLOCKED)],
                             capture_output=True, timeout=600, cwd=tmp_path, env=environment)
    assert process.returncode == 0, process.stderr.decode("utf-8", "replace")
    outcome = json.loads(process.stdout.decode("utf-8"))
    guard = "network access is disabled in this offline check"
    assert [name for name, _ in outcome["controls"]] == ["raw socket connect", "datagram send", "name lookup",
                                                         "reverse lookup", "address lookup", "connection"]
    assert all(message == guard for _, message in outcome["controls"]), outcome["controls"]
    assert outcome["control_attempts"] >= 6 and outcome["network_attempts"] == []
    assert b"ResourceWarning" not in process.stderr
    assert outcome["loaded_optional"] == []
    assert outcome["manifest_hash"] == MANIFEST_PIN
    assert outcome["execution_hash"] == outcome["reference_content_hash"] == RUNTIME_PINS["execution"]
    assert outcome["records"] == 40 and outcome["again"] is True
    assert list(home.iterdir()) == []                          # no cache or settings written


def readonly(path, *, on: bool) -> None:
    """Make ``path`` and everything below it read-only (``on``) or writable again."""
    import os
    import stat
    entries = sorted(path.rglob("*"), reverse=True) + [path]
    for entry in (entries if on else reversed(entries)):
        if os.name == "nt":
            os.chmod(entry, stat.S_IREAD if on else stat.S_IREAD | stat.S_IWRITE)
        else:
            os.chmod(entry, (0o555 if entry.is_dir() else 0o444) if on else (0o755 if entry.is_dir() else 0o644))


def test_the_starter_loads_from_a_read_only_package_location(tmp_path, monkeypatch):
    import sys
    import types
    from GenAIRR._reference_sources import build as build_module
    from GenAIRR._reference_sources import bundled
    from GenAIRR._reference_sources.digests import manifest_hash
    from GenAIRR.reference_sources import ReferenceBuild
    from tests.support.fixtures.source_api import tree
    root = tmp_path / "site"
    dossier = root / "genairr_readonly_pkg" / "published"
    dossier.mkdir(parents=True)
    (root / "genairr_readonly_pkg" / "__init__.py").write_bytes(b"")
    (dossier / f"{STARTER_ID}.zip").write_bytes(packaged_archive())
    monkeypatch.syspath_prepend(str(root))
    monkeypatch.delitem(sys.modules, "genairr_readonly_pkg", raising=False)
    entry = bundled.BundledDossier("genairr_readonly_pkg", f"published/{STARTER_ID}.zip", MANIFEST_PIN)
    monkeypatch.setattr(bundled, "REGISTRY", types.MappingProxyType({STARTER_ID: entry}))
    monkeypatch.setattr(build_module, "_BUNDLED", {})
    readonly(root, on=True)
    try:
        probe = dossier / "probe"
        try:
            probe.write_bytes(b"")
        except OSError:
            pass
        else:
            probe.unlink()
            pytest.skip("this platform or user cannot make a directory read-only (for example Windows, or root)")
        before = tree(root)
        loaded = ReferenceBuild.bundled(STARTER_ID)
        after = tree(root)
    finally:
        readonly(root, on=False)
    assert manifest_hash(loaded.manifest) == MANIFEST_PIN and runtime_hashes(loaded.reference) == RUNTIME_PINS
    assert after == before                                      # nothing was written into the package


def test_the_starter_needs_no_source_cache_or_network(tmp_path, monkeypatch):
    from GenAIRR._reference_sources import acquisition, api
    from GenAIRR._reference_sources import build as build_module
    from GenAIRR._reference_sources import cache as cache_module
    from GenAIRR._reference_sources.digests import manifest_hash
    from GenAIRR.reference_sources import ReferenceBuild

    def no_source(*args, **kwargs):
        pytest.fail("a packaged build opened a source cache or acquired a source")

    monkeypatch.setattr(build_module, "_BUNDLED", {})
    refuse_network(monkeypatch)
    monkeypatch.setattr(cache_module.BlobCache, "__init__", no_source)
    for owner in (acquisition, api):
        for name in ("fetch_source", "ingest_files", "read_locked"):
            monkeypatch.setattr(owner, name, no_source)
    home = tmp_path / "home"
    home.mkdir()
    for variable in ("HOME", "USERPROFILE", "XDG_CACHE_HOME", "APPDATA", "LOCALAPPDATA"):
        monkeypatch.setenv(variable, str(home))
    loaded = ReferenceBuild.bundled(STARTER_ID)
    assert manifest_hash(loaded.manifest) == MANIFEST_PIN and runtime_hashes(loaded.reference) == RUNTIME_PINS
    assert list(home.iterdir()) == []


# ── the starter's reasons in the readable report ─────────────────────


def report_section(text: str, title: str) -> str:
    """The lines from ``title`` up to the next unindented line."""
    lines = text.splitlines()
    start = next(i for i, line in enumerate(lines) if line.startswith(title))
    end = next((i for i in range(start + 1, len(lines)) if lines[i] and not lines[i].startswith(" ")), len(lines))
    return "\n".join(lines[start:end])


def test_the_report_text_states_every_selection_and_decision_reason(starter):
    document = recipe_document()
    text = str(starter.report)
    reasons_section = report_section(text, "Selection reasons")
    named = []
    for entry in document["selection"]:
        assert entry["reason"] in reasons_section, entry["label"]
    for line in reasons_section.splitlines()[1:]:
        if line.startswith("      "):
            named += [item for item in line.strip().rstrip(",").split(", ") if item]
    assert sorted(named) == sorted(f"{e['segment']}:{e['label']}" for e in document["selection"])     # each once
    classes = imgt_classes()
    for label in SELECTION["selected"]["V"]:
        if label in classes and classes[label][0] in IMGT_DISCLOSED:
            assert IMGT_CLAUSE.format(word=IMGT_DISCLOSED[classes[label][0]]) in reasons_section, label
    decisions = report_section(text, "Recipe decisions")
    assert len(document["decisions"]) == 14
    for index, decision in enumerate(document["decisions"]):
        (line,) = [line for line in decisions.splitlines() if line.startswith(f"  {index}. ")]
        assert line.startswith(f"  {index}. {decision['action']}; ")
        assert all(finding in line for finding in decision["findings"])
        assert line.endswith(f"reason: {decision['reason']}")


# ── the generator refuses to run under another version than the checkout declares ─


def test_the_generator_refuses_an_installed_version_the_checkout_does_not_declare(tmp_path, monkeypatch, capsys):
    import GenAIRR
    tool = load_tool()
    declared = declared_version((REPO / "pyproject.toml").read_bytes().decode("utf-8"))
    monkeypatch.setattr(GenAIRR, "__version__", "0.0.0+stale")
    target = tmp_path / "starter.zip"
    assert tool.main(["--output", str(target)]) == 2
    message = capsys.readouterr().err
    assert "0.0.0+stale" in message and declared in message and "pyproject.toml" in message
    assert not target.exists()


def test_the_generator_refuses_to_run_without_a_checkout_pyproject(tmp_path, monkeypatch, capsys):
    tool = load_tool()
    monkeypatch.setattr(tool, "PYPROJECT", tmp_path / "absent" / "pyproject.toml")
    assert tool.main(["--output", str(tmp_path / "starter.zip")]) == 2
    assert "pyproject.toml" in capsys.readouterr().err and not (tmp_path / "starter.zip").exists()
    crlf = tmp_path / "pyproject.toml"
    crlf.write_bytes(b'[project]\r\nname = "GenAIRR"\r\nversion = "' + GenAIRR_version().encode() + b'"\r\n')
    assert tool.declared_version(crlf.read_bytes().decode("utf-8")) == GenAIRR_version()


def GenAIRR_version() -> str:
    import GenAIRR
    return GenAIRR.__version__
