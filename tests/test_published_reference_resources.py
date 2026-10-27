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
import json
import re

import pytest

from GenAIRR._artifact_io import canonical_json
from GenAIRR._reference_sources.recipe import Recipe
from tests.test_ogrdb_source_inventory import IGH, RECORDS, SELECTION, raw, raw_records
from tests.test_reference_source_api import DATA, PROFILE, from_files

STARTER_ID = "human_igh_ogrdb_10_observed_v1"
RECIPE_PATH = DATA / f"{STARTER_ID}.recipe.json"
CHOICES = [row["label"] for row in SELECTION["record_choices"]]
REQUIRES = [{"capability": "observed_span_simulation", "evidence": "assertion"},
            {"capability": "source_annotated_junction", "evidence": "source"},
            {"capability": "canonical_productive", "evidence": "assertion"}]
INTERNAL_WORDS = re.compile(r"\b(?:task|plan|ledger|milestone|agent|C[0-5])\b", re.IGNORECASE)
AFFIRMED_UP_TO = re.compile(r"up to and including (?:nucleotide|base) (\d+)")
# Conserved IMGT codons other than the junction anchor, as 0-based gapped columns.
CONSERVED = ((23, 66, "C", "Cys"), (41, 120, "W", "Trp"))
CODONS = {a + b + c: aa for (a, b, c), aa in zip(
    ((a, b, c) for a in "TCAG" for b in "TCAG" for c in "TCAG"),
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")}
RESIDUES = {"C": "Cys", "W": "Trp", "F": "Phe"}


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
    lock = from_files(tmp_path, IGH).lock
    assert document["source_lock_hash"] == lock.commitment
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
    expected = [{"segment": segment, "label": label, "record": {"asset": "airr", "pointer": f"{RECORDS}/{index[label]}"},
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
        if start is not None and start == record["gene_end"] + 1 and not record["sequence"][start - 1:].startswith("CAC"):
            orf[label] = f"the record bases at its declared V RSS start read {record['sequence'][start - 1:start + 6]}"
        for position, column, conserved, name in CONSERVED:
            codon = sequence[column:column + 3]
            if CODONS[codon] != conserved:
                orf[label] = (f"its IMGT {position} codon {codon} encodes {RESIDUES[CODONS[codon]]}, not the conserved "
                              f"{name}")
    assert sorted(orf) == ["IGHV3-16*02", "IGHV3-16*03", "IGHV3-16*04", "IGHV3-20*02", "IGHV3-25*04",
                           "IGHV8-51-1*04", "IGHV8-51-1*05"]
    for label, observation in orf.items():
        assert "ORF-like although the provider's functional field is true: " + observation in given[label]
        assert "No runtime functionality is claimed." in given[label]
    assert sum("ORF-like" in text for text in given.values()) == 7
    (observation,) = [row for row in SELECTION["selected_record_observations"]
                      if row["observation"] == "v_rss_start_inside_coding_interval"]
    assert (f"Its declared v_rs_start ({observation['v_rs_start']}) lies inside its coding interval, while the record "
            f"bases after gene_end read {observation['sequence_after_gene_end']}") in given[observation["label"]]


def test_no_reason_overstates_functionality_or_names_internal_work():
    for text in list(reasons().values()) + [d["reason"] for d in recipe_document()["decisions"]]:
        assert text.isascii() and "\n" not in text
        assert "functional allele" not in text.lower()
        assert not INTERNAL_WORDS.search(text), text


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


@pytest.fixture(scope="module")
def generated(tmp_path_factory):
    """The starter written by the generator's command line into a new directory, offline."""
    import contextlib
    import io
    target = tmp_path_factory.mktemp("generated") / STARTER_ID
    out = io.StringIO()
    with pytest.MonkeyPatch.context() as patch, contextlib.redirect_stdout(out):
        refuse_network(patch)
        code = load_tool().main(["--output", str(target)])
    return code, target, out.getvalue()


def test_the_generator_writes_exactly_the_dossier_members_offline(generated):
    from GenAIRR._reference_sources.digests import manifest_hash
    from GenAIRR.reference_sources import ReferenceBuild
    code, target, printed = generated
    assert code == 0
    checked = ReferenceBuild.read_dossier(target)
    listed = sorted([member.name for member in checked.manifest.members] + ["manifest.json"])
    present = sorted(p.relative_to(target).as_posix() for p in target.rglob("*") if p.is_file())
    assert present == listed and not (target / "receipts").exists()
    assert sorted(p.name for p in target.iterdir() if p.is_dir()) == ["cache", "inputs"]
    reference = checked.reference
    for line in (f"manifest {manifest_hash(checked.manifest)}", f"runtime execution {reference.execution_hash}",
                 f"runtime catalogue {reference.catalogue_hash}", f"runtime models {reference.models_hash}",
                 f"runtime artifact {reference.artifact_hash}", f"builder GenAIRR {checked.manifest.builder.version}"):
        assert line in printed.splitlines(), line
    assert (target / "recipe.json").read_bytes() == RECIPE_PATH.read_bytes()


def test_the_generator_writes_only_into_a_new_or_empty_directory(tmp_path, capsys):
    target = tmp_path / "occupied"
    target.mkdir()
    (target / "notes.txt").write_bytes(b"kept\n")
    assert load_tool().main(["--output", str(target)]) == 2
    assert "is not empty" in capsys.readouterr().err
    assert [p.name for p in target.iterdir()] == ["notes.txt"] and (target / "notes.txt").read_bytes() == b"kept\n"
    assert load_tool().main(["--output", str(tmp_path / "missing-parent" / "dossier")]) == 2
    assert not (tmp_path / "missing-parent").exists()
