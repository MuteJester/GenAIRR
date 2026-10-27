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
import pathlib
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


# ── the packaged starter ─────────────────────────────────────────────

# Pinned from the reviewed build: the typed manifest hash of the packaged dossier (it names
# the GenAIRR version that generated it) and the runtime hashes of its reference.
MANIFEST_PIN = "sha256:adc1961adbe134648d77cf4e8adec09f81465ead027112f6bca5ed7a53bda440"
RUNTIME_PINS = {"execution": "sha256:7c2550263f5a26ebc1fb5369d8d8916c5e2e24e059f5c6661a47d7ec2cace40e",
                "catalogue": "sha256:f404dd27fee06c95c69a41bb79e6cd17b555039750a9d95012ea3c8231ecd4d0",
                "models": "sha256:e96b1404c64d36a55c0457163ef3a4aa706e1bc5114914147715cd0f964dffd8"}
RECORD_OUTCOMES = {"selected": 230, "excluded": 259, "rejected": 0}
PACKAGE, RESOURCE = "GenAIRR.data", f"published_references/{STARTER_ID}"
DOSSIER_DIRECTORIES = ["cache", "cache/blobs", "cache/blobs/sha256", "inputs"]


@pytest.fixture(scope="module")
def starter():
    from GenAIRR.reference_sources import ReferenceBuild
    return ReferenceBuild.bundled(STARTER_ID)


def runtime_hashes(reference) -> dict:
    return {"execution": reference.execution_hash, "catalogue": reference.catalogue_hash,
            "models": reference.models_hash}


def packaged_files() -> dict:
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
    assert len(PROJECTION_LOSSES) == 4
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


def test_rebuilding_from_the_pinned_bytes_and_the_recipe_reproduces_the_packaged_starter(generated, starter):
    import GenAIRR
    from GenAIRR.reference_sources import ReferenceBuild
    code, target, _ = generated
    assert code == 0
    rebuilt = ReferenceBuild.read_dossier(target)
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
    files, _ = packaged_files()
    assert {p.relative_to(target).as_posix(): p.read_bytes() for p in target.rglob("*") if p.is_file()} == files


def test_the_packaged_dossier_holds_only_its_members(starter):
    files, folders = packaged_files()
    assert sorted(files) == sorted([member.name for member in starter.manifest.members] + ["manifest.json"])
    assert folders == DOSSIER_DIRECTORIES                    # no receipts, no __pycache__, no package marker
    assert files["recipe.json"] == RECIPE_PATH.read_bytes()


def test_every_dossier_member_is_packaged_in_the_wheel_and_the_sdist(starter):
    import glob
    tomllib = pytest.importorskip("tomllib", reason="reading pyproject.toml needs Python 3.11 or later")
    source = REPO / "src" / "GenAIRR" / "data" / "published_references" / STARTER_ID
    if not source.is_dir():
        pytest.skip("the source checkout is not available")
    members = {(source / name).resolve() for name in
               [member.name for member in starter.manifest.members] + ["manifest.json"]}
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


def test_loading_the_starter_leaves_the_human_igh_alias_unchanged(starter):
    import GenAIRR
    from GenAIRR import Experiment
    from GenAIRR._refdata_resolver import _CONFIG_ALIASES
    from GenAIRR.reference_sources import ReferenceBuild
    before = dict(_CONFIG_ALIASES)
    assert ReferenceBuild.bundled(STARTER_ID) is starter
    assert dict(_CONFIG_ALIASES) == before and _CONFIG_ALIASES["human_igh"] == "HUMAN_IGH_OGRDB"
    assert not [alias for alias in _CONFIG_ALIASES if "observed" in alias or "published" in alias]
    with pytest.raises(ValueError, match="Unknown config name"):
        Experiment.on(STARTER_ID)
    alias, cartridge = Experiment.on("human_igh"), Experiment.on(GenAIRR.HUMAN_IGH_OGRDB)
    assert (alias.refdata.v_pool_size(), alias.refdata.j_pool_size()) == (
        cartridge.refdata.v_pool_size(), cartridge.refdata.j_pool_size())
    first = alias.recombine().run_records(n=20, seed=3)
    second = Experiment.on(GenAIRR.HUMAN_IGH_OGRDB).recombine().run_records(n=20, seed=3)
    assert first.reference_content_hash == second.reference_content_hash != starter.reference.execution_hash
    assert [record["sequence"] for record in first] == [record["sequence"] for record in second]



def test_the_packaged_dossier_is_checked_out_byte_for_byte_on_every_platform():
    """Its digests cover exact LF bytes, so no checkout may convert its line endings."""
    attributes = REPO / ".gitattributes"
    if not attributes.is_file():
        pytest.skip("the source checkout is not available")
    rules = [line.split() for line in attributes.read_bytes().decode("utf-8").splitlines()
             if line.strip() and not line.startswith("#")]
    assert ["src/GenAIRR/data/published_references/**", "-text"] in rules

