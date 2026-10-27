"""Explicit recipe resolution and deterministic projection into an unchanged runtime reference.

``projection.project`` regenerates the interpreted inventory of a locked source itself,
resolves an explicit recipe against it (one outcome per source record and per claim), writes
deterministic FASTA and annotation-sidecar bytes, builds the runtime reference through the
public ``Reference.from_fasta`` and returns the complete build report.

Expected values come from the independently authored fixture files
(``human_igh_release10_candidate_selection.json``, ``expectations.json``, ``sources.json``),
from literals, or from reading raw payload fields with the standard library; never from the
code under test.  Synthetic sources (length boundary, frameshift, Unicode labels, line
endings) are registered as validated for one test only, by monkeypatching the profile scope.
"""
from __future__ import annotations

import collections
import dataclasses
import hashlib
import json
import re
import tempfile

import pytest

import GenAIRR
from GenAIRR._reference.artifact import Reference
from GenAIRR._reference.errors import ReferenceArtifactError
from GenAIRR._reference_sources import common, evidence, ogrdb, projection
from GenAIRR._reference_sources import model as m
from tests.test_ogrdb_source_coordinates import V_LABEL, in_scope, small
from tests.test_ogrdb_source_inventory import (
    EXPECT, IGH, IGL, RECORDS, SELECTION, error, pinned, raw, raw_records, seal,
)

J1, J7 = "IGLJ1*01", "IGLJ7*02"
MINI = (V_LABEL, J1)
CHOICES = {row["label"]: row for row in SELECTION["record_choices"]}
OMITTED = {row["label"]: [rule["rule"] for rule in row["failed_rules"]] for row in SELECTION["omitted"]}
# The fixture's required capability, plus its targeted capabilities at their target evidence,
# so that the build itself confirms them.
STARTER_REQUIRES = (("observed_span_simulation", "assertion"), ("source_annotated_junction", "source"),
                    ("canonical_productive", "assertion"))
SIDECAR_ALLELE_KEYS = {"segment", "label", "sequence_sha256", "gene", "functionality", "anchor", "coding_phase",
                       "boundary_5", "boundary_3", "notes"}


# ── helpers ──────────────────────────────────────────────────────────


def sha(sequence: str) -> str:
    return "sha256:" + hashlib.sha256(sequence.encode("ascii")).hexdigest()


def raw_sha(data: bytes) -> str:
    return "sha256:" + hashlib.sha256(data).hexdigest()


def records_in(blobs, fmt) -> list:
    """The records of one record export, read with the standard library."""
    return json.loads(blobs[fmt].decode("utf-8"))["GermlineSet"][0]["allele_descriptions"]


def entry(blobs, fmt, label, reason="selected for this test") -> dict:
    records = records_in(blobs, fmt)
    (index,) = [i for i, r in enumerate(records) if r["label"] == label]
    return {"segment": records[index]["sequence_type"], "label": label,
            "record": {"asset": fmt, "pointer": f"{RECORDS}/{index}"},
            "sequence_sha256": sha(records[index]["coding_sequence"]), "reason": reason}


def recipe(lock, selection, *, requires=(("observed_span_simulation", "assertion"),), decisions=(),
           name="test reference", **reference) -> dict:
    identity = lock.identity
    head = {"name": name, "version": "1", "taxon": identity.taxon.to_dict(), "locus": identity.locus,
            "chain_type": identity.chain_type, "build_date": "2026-10-27"}
    head.update(reference)
    return {"artifact": "genairr_reference_recipe", "schema_version": 1, "source_lock_hash": lock.commitment,
            "adapter": lock.adapter.to_dict(), "projection": {"id": "observed_span", "version": 1},
            "reference": head, "selection": list(selection), "decisions": list(decisions),
            "requires": [{"capability": c, "evidence": e} for c, e in requires],
            "defaults_profile": "genairr.unfitted_defaults.v1"}


def build(lock, blobs, document):
    return projection.project(lock, blobs, m.Recipe.from_dict(document))


def mini(monkeypatch, labels=MINI, **edits):
    """A small synthetic IGL source over ``labels``, registered as validated for this test."""
    blobs = small(IGL, set(labels), **edits)
    lock = seal(IGL, blobs)
    in_scope(monkeypatch, lock)
    return blobs, lock


def mini_recipe(blobs, lock, labels=MINI, fmt="airr", **kwargs) -> dict:
    return recipe(lock, [entry(blobs, fmt, label) for label in labels], **kwargs)


def interpreted(blobs, lock):
    return evidence.interpret(lock, blobs, ogrdb.build_inventory(lock, blobs))


def no_native(monkeypatch):
    """Fail the test if the runtime reference importer (and so the native engine) is entered."""
    def refuse(*args, **kwargs):
        raise AssertionError("the runtime reference importer was entered")
    monkeypatch.setattr(Reference, "from_fasta", refuse)


def expected_fasta(items) -> bytes:
    """``>label``, the sequence in 80-column lines, LF line ends and a trailing LF."""
    lines = []
    for label, sequence in items:
        lines.append(">" + label)
        lines.extend(sequence[i:i + 80] for i in range(0, len(sequence), 80))
    return ("\n".join(lines) + "\n").encode("utf-8")


def artifact(reference) -> dict:
    return json.loads(reference.to_json())


def withhold(claims, digest, reason="withheld for this test") -> dict:
    return {"action": "withhold_claim", "findings": [], "claims": list(claims), "sequence_sha256": digest,
            "reason": reason, "evidence": []}


def claim_id(inventory, fmt, label, predicate) -> str:
    (record,) = [r for r in inventory.records if r.locator.asset == fmt and r.label == label]
    (candidate,) = inventory.candidates_of(record.locator)
    (claim,) = [c for c in inventory.claims if c.subject == candidate.locator and c.predicate == predicate]
    return claim.id


# ── the starter: human IGH release 10 ───────────────────────────────


def _airr_index() -> dict:
    return {r["label"]: i for i, r in enumerate(raw_records(IGH, "airr"))}


def starter_recipe(lock, inventory) -> dict:
    """The selection fixture as a recipe: every selected allele from its own 'airr' record, in
    fixture order, and one select_alternative decision per recorded record choice."""
    index, records = _airr_index(), raw_records(IGH, "airr")
    selection = []
    for segment in ("V", "D", "J"):
        for label in SELECTION["selected"][segment]:
            reason = SELECTION["record_choice_reason"] if label in CHOICES else "selected by the starter rules"
            selection.append({"segment": segment, "label": label,
                              "record": {"asset": "airr", "pointer": f"{RECORDS}/{index[label]}"},
                              "sequence_sha256": sha(records[index[label]]["coding_sequence"]), "reason": reason})
    differs = {f.observations.get()["label"]: f.id for f in inventory.findings
               if f.code == "export_coding_sequence_differs"}
    decisions = [{"action": "select_alternative", "findings": [differs[label]], "claims": [],
                  "sequence_sha256": sha(records[index[label]]["coding_sequence"]),
                  "reason": SELECTION["record_choice_reason"], "evidence": []} for label in CHOICES]
    return recipe(lock, selection, requires=STARTER_REQUIRES, decisions=decisions,
                  name="human IGH starter (OGRDB IGH_VDJ release 10)")


@pytest.fixture(scope="module")
def igh():
    blobs = pinned(IGH)
    lock = seal(IGH, blobs)
    return blobs, lock, interpreted(blobs, lock)


@pytest.fixture(scope="module")
def starter(igh):
    blobs, lock, inventory = igh
    document = starter_recipe(lock, inventory)
    return document, build(lock, blobs, document)


def test_starter_fasta_holds_exactly_the_selected_labels_and_supplied_dna(starter):
    _, result = starter
    index, records = _airr_index(), raw_records(IGH, "airr")
    for segment in ("V", "D", "J"):
        labels = SELECTION["selected"][segment]
        assert len(labels) == SELECTION["counts"][segment]
        expected = expected_fasta([(label, records[index[label]]["coding_sequence"]) for label in labels])
        assert result.fasta_of(segment) == expected
    assert sorted(result.inputs) == ["annotations.json", "d.fasta", "j.fasta", "v.fasta"]
    for data in result.inputs.values():
        assert b"\r" not in data and data.endswith(b"\n") and not data.endswith(b"\n\n")
        data.decode("utf-8")


def test_starter_runtime_alleles_are_exactly_the_selected_records(starter):
    _, result = starter
    index, records = _airr_index(), raw_records(IGH, "airr")
    expected = [(segment, label, sha(records[index[label]]["coding_sequence"]))
                for segment in ("V", "D", "J") for label in SELECTION["selected"][segment]]
    assert [(a.segment, a.label, a.sequence_sha256) for a in result.reference.alleles] == expected
    assert {r.status for r in result.reference.records} == {"accepted"}
    usage = artifact(result.reference)["models"]["allele_usage"]
    for segment in ("V", "D", "J"):     # one usage weight per selected identity, none added by relationships
        assert usage[segment] == SELECTION["selected"][segment]
    report = result.report
    assert report.counts.selected_alleles == m.SegmentCounts(V=SELECTION["counts"]["V"], D=SELECTION["counts"]["D"],
                                                             J=SELECTION["counts"]["J"], unassigned=0)
    assert [s.allele for s in report.selection] == [f"{s}:{label}" for s, label, _ in expected]


def _expected_record_outcomes() -> dict:
    """locator id -> (outcome, allele, exclusion, related locator id), from the selection fixture
    and the raw exports read with the standard library."""
    segment_of = {label: segment for segment in ("V", "D", "J") for label in SELECTION["selected"][segment]}
    airr, airr_ex = raw_records(IGH, "airr"), raw_records(IGH, "airr_ex")
    index = _airr_index()

    def omitted(label):
        rules = OMITTED[label]
        if "explicit_record_choice" in rules:
            return "unresolved_conflicting_sequence"
        if "v_imgt_104_is_cys" in rules:          # canonical_productive is required
            return "missing_required_evidence"
        return "outside_explicit_selection"

    expected = {}
    for fmt, records in (("airr", airr), ("airr_ex", airr_ex)):
        for i, record in enumerate(records):
            label, at = record["label"], f"{fmt}#{RECORDS}/{i}"
            if label not in segment_of:
                expected[at] = ("excluded", None, omitted(label), None)
            elif fmt == "airr":
                expected[at] = ("selected", f"{segment_of[label]}:{label}", None, None)
            else:
                same = record["coding_sequence"] == airr[index[label]]["coding_sequence"]
                expected[at] = ("excluded", None, "duplicate_occurrence" if same else "alternative_for_selected_identity",
                                f"airr#{RECORDS}/{index[label]}")
    return expected


def test_every_source_record_has_exactly_one_outcome_and_the_counts_reconcile(starter):
    _, result = starter
    report = result.report
    expected = _expected_record_outcomes()
    assert len(expected) == EXPECT["sets"]["9606.IGH_VDJ/10"]["airr"]["records"] + \
        EXPECT["sets"]["9606.IGH_VDJ/10"]["airr_ex"]["records"]
    assert [o.record.id for o in report.records] == list(expected)
    for outcome in report.records:
        related = None if outcome.related is None else outcome.related.id
        assert (outcome.outcome, outcome.allele, outcome.exclusion, related) == expected[outcome.record.id]
    exclusions = collections.Counter(o.exclusion for o in report.records if o.outcome == "excluded")
    assert exclusions == {"duplicate_occurrence": 207, "alternative_for_selected_identity": 14,
                          "unresolved_conflicting_sequence": 18, "missing_required_evidence": 2,
                          "outside_explicit_selection": 18}
    assert {report.inventory.record(o.record).label for o in report.records
            if o.exclusion == "alternative_for_selected_identity"} == set(CHOICES)
    assert report.counts.record_outcomes == m.OutcomeCounts(selected=230, excluded=259, rejected=0)
    for counts in report.counts.exports:
        stated = EXPECT["sets"]["9606.IGH_VDJ/10"][counts.asset]
        assert counts.record_occurrences == m.SegmentCounts(V=stated["V"], D=stated["D"], J=stated["J"], unassigned=0)
    (tyr,) = [o for o in report.records if o.record.id == f"airr#{RECORDS}/{_airr_index()['IGHV1-38-4*01']}"]
    assert "canonical_productive" in tyr.detail


def test_the_starter_capabilities_are_confirmed_on_the_built_reference(starter):
    _, result = starter
    for capabilities in (result.reference.capabilities(), result.report.capabilities):
        status = capabilities.status
        assert status("observed_span_simulation").available
        assert (status("source_annotated_junction").available, status("source_annotated_junction").evidence) == (
            True, "source")
        assert (status("canonical_productive").available, status("canonical_productive").evidence) == (
            True, "assertion")
        assert capabilities.phase_support == "anchor_derived"
        assert not status("original_coding_ends").available
    target = SELECTION["capability_target"]["basis_for_target"]["canonical_productive"]
    assert "evidence 'assertion' and phase support 'anchor_derived'" in target
    assert [r.to_dict() for r in result.report.requirements] == [
        {"capability": c, "evidence": e, "satisfied": True} for c, e in STARTER_REQUIRES]


def _imgt_delineation(record) -> tuple:
    (j,) = [j for j, d in enumerate(record["v_gene_delineations"]) if d["delineation_scheme"] == "IMGT"]
    return j, record["v_gene_delineations"][j]


def _gapped_index(fmt="gapped") -> dict:
    headers = [chunk.partition("\n")[0].rstrip("\r") for chunk in raw(IGH, fmt).decode("utf-8").split(">")[1:]]
    return {header: n for n, header in enumerate(headers)}


def test_selected_v_and_j_carry_one_accepted_source_anchor_and_an_anchor_derived_phase(starter):
    _, result = starter
    index, records = _airr_index(), raw_records(IGH, "airr")
    claims = artifact(result.reference)["catalogue"]["claims"]
    assert not [c for c in claims if c["derivation"] == "imgt_alignment"]
    assert not [c for c in claims if c["predicate"] in ("functionality", "boundary_5", "boundary_3", "v_subregions")]
    by_subject = collections.defaultdict(list)
    for claim in claims:
        by_subject[(claim["subject"], claim["predicate"])].append(claim)
    for segment in ("V", "J"):
        for label in SELECTION["selected"][segment]:
            record = records[index[label]]
            if segment == "V":
                anchor = _imgt_delineation(record)[1]["cdr3_start"] - 4
                phase, phase_support = anchor % 3, "assertion"
            else:
                anchor, phase, phase_support = record["j_cdr3_end"] - 1, record["j_codon_frame"] - 1, "source"
            allele = f"{segment}:{label}"
            accepted = [c for c in by_subject[(allele, "anchor")] if c["status"] == "accepted"]
            assert [(c["value"], c["support"], c["derivation"]) for c in accepted] == [(anchor, "source", "stated")]
            others = [c for c in by_subject[(allele, "anchor")] if c["status"] != "accepted"]
            assert all((c["status"], c["derivation"]) == ("candidate", "unique_motif") for c in others)
            assert segment == "J" or not others
            (coding_phase,) = by_subject[(allele, "coding_phase")]
            assert (coding_phase["value"], coding_phase["support"], coding_phase["derivation"], coding_phase["status"]) \
                == (phase, phase_support, "anchor_derived", "accepted")
    for label in SELECTION["selected"]["D"]:
        assert not by_subject[(f"D:{label}", "anchor")] and not by_subject[(f"D:{label}", "coding_phase")]


def test_annotations_cite_the_exact_source_fields_and_notes_cite_record_and_build_inputs(starter):
    _, result = starter
    index, records = _airr_index(), raw_records(IGH, "airr")
    airr_sha, gapped_sha = raw_sha(raw(IGH, "airr")), raw_sha(raw(IGH, "gapped"))
    gapped = _gapped_index()
    alleles = {(a.segment, a.label): a for a in result.reference.alleles}
    methods = {(c["subject"], c["predicate"]): c["method"] for c in artifact(result.reference)["catalogue"]["claims"]
               if c["status"] == "accepted"}
    for segment in ("V", "D", "J"):
        for label in SELECTION["selected"][segment]:
            allele, pointer = alleles[(segment, label)], f"{RECORDS}/{index[label]}"
            notes = " ".join(allele.notes)
            assert all(len(note) <= 4096 for note in allele.notes)
            for part in ("airr", airr_sha, pointer, result.build_input_hash, "observed_span"):
                assert part in notes
            if segment == "V":
                j, _ = _imgt_delineation(records[index[label]])
                for part in (airr_sha, f"{pointer}/v_gene_delineations/{j}/cdr3_start",
                             f"gapped {gapped_sha} /{gapped[label]}"):
                    assert part in allele.anchor.source
                assert "cdr3_start - 4" in methods[(f"V:{label}", "anchor")]
                assert "ogrdb_germline_set" in methods[(f"V:{label}", "anchor")]
            elif segment == "J":
                assert "j_cdr3_end - 1" in methods[(f"J:{label}", "anchor")]
                assert "j_codon_frame - 1" in methods[(f"J:{label}", "coding_phase")]
                assert f"{pointer}/j_cdr3_end" in allele.anchor.source
                assert f"{pointer}/j_cdr3_end" in allele.coding_phase.source
                assert f"{pointer}/j_codon_frame" in allele.coding_phase.source
            for name in ("functionality", "boundary_5", "boundary_3"):
                assert getattr(allele, name) is None


def test_identical_dna_under_distinct_labels_stays_two_identities_with_their_own_evidence(starter):
    _, result = starter
    index = _airr_index()
    alleles = {a.label: a for a in result.reference.alleles}
    record_pointer = re.compile(re.escape(RECORDS) + r"/([0-9]+)")
    for first, second in SELECTION["identical_coding_sequence_groups_selected"]:
        a, b = alleles[first], alleles[second]
        assert a.sequence_sha256 == b.sequence_sha256 and a.label != b.label
        for own in (a, b):
            cited = " ".join(own.notes) + " ".join(
                getattr(own, name).source for name in ("anchor", "coding_phase") if getattr(own, name) is not None)
            assert {int(n) for n in record_pointer.findall(cited)} == {index[own.label]}


def test_the_sidecar_is_canonical_lf_json_in_the_v1_shape_without_provider_fields(starter):
    document, result = starter
    data = result.annotations
    sidecar = json.loads(data.decode("utf-8"))
    canonical = json.dumps(sidecar, ensure_ascii=False, sort_keys=True, separators=(",", ":")) + "\n"
    assert data == canonical.encode("utf-8")
    assert set(sidecar) == {"schema", "schema_version", "reference", "conventions", "alleles", "relationships"}
    assert (sidecar["schema"], sidecar["schema_version"]) == ("genairr_reference_annotations", 1)
    assert sidecar["conventions"] == {"coordinates": "zero_based_half_open", "v_alignment": "none"}
    assert sidecar["reference"] == document["reference"]
    for allele in sidecar["alleles"]:
        assert set(allele) == SIDECAR_ALLELE_KEYS and allele["gene"] is None


def test_the_runtime_reader_accepts_the_generated_files_unchanged(starter, tmp_path):
    _, result = starter
    for name, data in result.inputs.items():
        (tmp_path / name).write_bytes(data)
    rebuilt = Reference.from_fasta(v=tmp_path / "v.fasta", d=tmp_path / "d.fasta", j=tmp_path / "j.fasta",
                                   annotations=tmp_path / "annotations.json")
    assert rebuilt.artifact_hash == result.reference.artifact_hash
    assert Reference.from_json(result.reference.to_json()).artifact_hash == result.reference.artifact_hash


def test_the_report_records_hashes_generator_and_claim_outcomes(starter, igh):
    document, result = starter
    blobs, lock, inventory = igh
    report, reference = result.report, result.reference
    recipe_ = m.Recipe.from_dict(document)
    assert report.inventory == inventory
    assert report.hashes.build_inputs == result.build_input_hash == m.build_input_hash(lock, inventory, recipe_)
    assert report.runtime == m.RuntimeHashes(reference.catalogue_hash, reference.models_hash,
                                             reference.execution_hash, reference.artifact_hash)
    assert report.generator == m.Generator("GenAIRR", GenAIRR.__version__)
    assert report.comparison.complete
    selected = set()
    for entry_ in recipe_.selection:
        selected |= {entry_.record.id, entry_.record.id + "/coding_sequence"}
    outcomes = collections.Counter()
    for claim, outcome in zip(inventory.claims, report.claims):
        assert outcome.claim == claim.id and outcome.decision is None
        if claim.subject.id not in selected:
            assert outcome.outcome == "not_selected"
        elif claim.predicate in ("phase_independence", "coding_interval", "extension"):
            assert outcome.outcome == "not_representable" and outcome.reason
        elif claim.status == "unresolved":
            assert outcome.outcome == "withheld" and outcome.reason
        else:
            assert claim.predicate in ("anchor", "coding_phase") and outcome.outcome == "projected"
            if claim.predicate == "coding_phase" and inventory.subject_record(claim).segment == "J":
                assert "anchor_derived" in outcome.reason
        outcomes[outcome.outcome] += 1 if claim.subject.id in selected else 0
    v, d, j = (SELECTION["counts"][s] for s in ("V", "D", "J"))
    assert outcomes["projected"] == 2 * (v + j)                 # anchor and phase of every V and J
    assert outcomes["withheld"] == 3 * (v + d + j)              # functionality and both boundaries, unresolved
    assert outcomes["not_representable"] == (v + j) + (v + d + j)   # phase independence; coding interval
    assert m.BuildReport.from_dict(report.to_dict()) == report


def test_projection_losses_name_what_the_runtime_reference_does_not_hold():
    assert [code for code, _ in projection.PROJECTION_LOSSES] == [
        "v_subregions", "gapped_sequences", "affirmed_spans", "genomic_context"]
    assert all(text.strip() for _, text in projection.PROJECTION_LOSSES)


# ── what a build refuses before any projection ───────────────────────


def test_project_rejects_arguments_of_the_wrong_type(monkeypatch):
    blobs, lock = mini(monkeypatch)
    document = mini_recipe(blobs, lock)
    for args in ((lock.to_dict(), blobs, m.Recipe.from_dict(document)), (lock, blobs, document)):
        with pytest.raises(TypeError):
            projection.project(*args)


def test_a_build_needs_both_record_exports(monkeypatch):
    blobs, lock = mini(monkeypatch, formats=("airr", "gapped"))
    no_native(monkeypatch)
    found = error("missing_asset", build, lock, blobs, mini_recipe(blobs, lock))
    assert found.path == "/source_lock/assets" and "airr_ex" in found.message


def test_a_lock_outside_the_validated_profile_scope_cannot_build():
    blobs = small(IGL, set(MINI))
    lock = seal(IGL, blobs)
    assert error("unsupported_profile", build, lock, blobs, mini_recipe(blobs, lock)).path == "/source_lock"


def test_a_recipe_builds_only_from_its_own_lock(monkeypatch):
    blobs, lock = mini(monkeypatch)
    document = mini_recipe(blobs, lock)
    document["source_lock_hash"] = "sha256:" + "0" * 64
    assert error("digest_mismatch", build, lock, blobs, document).path == "/source_lock_hash"


@pytest.mark.parametrize("drop, missing", [("J", "J"), ("V", "V")])
def test_an_empty_segment_is_rejected_before_projection(monkeypatch, drop, missing):
    blobs, lock = mini(monkeypatch)
    no_native(monkeypatch)
    labels = [label for label in MINI if not label.startswith("IGL" + drop)]
    found = error("unresolved_selection", build, lock, blobs, mini_recipe(blobs, lock, labels))
    assert found.path == "/selection" and f"no {missing} allele" in found.message


def test_an_igh_recipe_without_a_d_allele_is_rejected(monkeypatch):
    blobs = pinned(IGH)
    lock = seal(IGH, blobs)
    no_native(monkeypatch)
    selection = [entry(blobs, "airr", label) for label in ("IGHV1-2*02", "IGHJ4*02")]
    found = error("unresolved_selection", build, lock, blobs, recipe(lock, selection))
    assert found.path == "/selection" and "no D allele" in found.message


@pytest.mark.parametrize("edit, kind, path", [
    (lambda d: d["reference"]["taxon"].update(scientific_name="Mus musculus"), "identity_mismatch", "/reference/taxon"),
    (lambda d: d["reference"]["taxon"].update(taxonomy_id=10090), "identity_mismatch", "/reference/taxon"),
    (lambda d: d["reference"].update(locus="IGK"), "identity_mismatch", "/reference/locus"),
    (lambda d: d["reference"].update(locus="TRA"), "identity_mismatch", "/reference/locus"),
    (lambda d: d["selection"][0].update(label="IGLV1-36*99"), "identity_mismatch", "/selection/0/label"),
    (lambda d: d["selection"][1].update(segment="V"), "identity_mismatch", "/selection/1/segment"),
])
def test_a_recipe_cannot_relabel_the_source_taxon_locus_segment_or_label(monkeypatch, edit, kind, path):
    blobs, lock = mini(monkeypatch)
    no_native(monkeypatch)
    document = mini_recipe(blobs, lock)
    edit(document)
    assert error(kind, build, lock, blobs, document).path == path


def test_a_chain_type_relabel_is_not_even_a_valid_recipe(monkeypatch):
    blobs, lock = mini(monkeypatch)
    document = mini_recipe(blobs, lock)
    document["reference"]["chain_type"] = "vdj"
    assert error("malformed", m.Recipe.from_dict, document).path == "/reference/chain_type"


def test_one_executable_label_cannot_select_two_records(monkeypatch):
    blobs, lock = mini(monkeypatch, ("IGLV10-54*04", J1))
    document = recipe(lock, [entry(blobs, "airr", "IGLV10-54*04"), entry(blobs, "airr_ex", "IGLV10-54*04"),
                             entry(blobs, "airr", J1)])
    assert error("unresolved_selection", m.Recipe.from_dict, document).path == "/selection"


def test_a_selection_binds_an_exact_record_and_digest_never_the_first_label_match(monkeypatch):
    blobs, lock = mini(monkeypatch, ("IGLV10-54*04", J1))
    no_native(monkeypatch)
    wrong = entry(blobs, "airr_ex", "IGLV10-54*04")
    wrong["sequence_sha256"] = entry(blobs, "airr", "IGLV10-54*04")["sequence_sha256"]
    found = error("unresolved_selection", build, lock, blobs, recipe(lock, [wrong, entry(blobs, "airr", J1)]))
    assert found.path == "/selection/0/sequence_sha256"


def test_stale_decisions_are_rejected(monkeypatch):
    blobs, lock = mini(monkeypatch)
    no_native(monkeypatch)
    inventory = interpreted(blobs, lock)
    anchor = claim_id(inventory, "airr", V_LABEL, "anchor")
    v_digest, j_digest = (entry(blobs, "airr", label)["sequence_sha256"] for label in MINI)
    # another sequence digest than the subject's selected one
    found = error("stale_decision", build, lock, blobs, mini_recipe(blobs, lock, decisions=[withhold([anchor], j_digest)]))
    assert found.path == "/decisions/0/sequence_sha256"
    # a claim about a record the recipe does not select
    other = claim_id(inventory, "airr_ex", V_LABEL, "anchor")
    found = error("stale_decision", build, lock, blobs, mini_recipe(blobs, lock, decisions=[withhold([other], v_digest)]))
    assert found.path == "/decisions/0/claims/0"
    # a claim of the same content under another lock (different source bytes)
    changed = dict(blobs)
    changed["gapped"] = blobs["gapped"] + b"\n"
    other_lock = seal(IGL, changed)
    in_scope(monkeypatch, other_lock)
    document = mini_recipe(changed, other_lock, decisions=[withhold([anchor], v_digest)])
    assert error("stale_decision", build, other_lock, changed, document).path == "/decisions/0/claims/0"


# ── a frameshifted source: a stated phase that is not its landmark modulo 3 ──


def _frameshift(fmt, record):
    """IGLJ7*02 states j_codon_frame 3 (phase 2) at j_cdr3_end 8 (anchor 7, and 7 % 3 = 1): a valid
    pair of source observations that the runtime's single phase scalar cannot hold together."""
    if record["label"] == J7:
        record["j_codon_frame"] = 3


def _frameshifted(monkeypatch):
    blobs, lock = mini(monkeypatch, (V_LABEL, J7), edit=_frameshift)
    inventory = interpreted(blobs, lock)
    (finding,) = [f for f in inventory.findings if f.code == "phase_landmark_projection_unsupported"]
    anchor, phase = (claim_id(inventory, "airr", J7, p) for p in ("anchor", "coding_phase"))
    return blobs, lock, finding, anchor, phase, entry(blobs, "airr", J7)["sequence_sha256"]


def test_a_frameshifted_selected_allele_without_a_withholding_decision_is_refused(monkeypatch):
    blobs, lock, finding, anchor, phase, digest = _frameshifted(monkeypatch)
    no_native(monkeypatch)
    found = error("projection", build, lock, blobs, mini_recipe(blobs, lock, (V_LABEL, J7)))
    assert found.path == "/selection/1"
    for part in (f"J:{J7}", finding.id, "phase_landmark_projection_unsupported", anchor, phase):
        assert part in found.message
    assert "exclude the allele" in found.message.lower() and "withhold_claim" in found.message
    # withholding only one of the two observations is not enough
    for claims in ([anchor], [phase]):
        document = mini_recipe(blobs, lock, (V_LABEL, J7), decisions=[withhold(claims, digest)])
        assert error("projection", build, lock, blobs, document).path == "/selection/1"


def test_a_frameshifted_allele_whose_anchor_and_phase_are_withheld_builds_with_null_annotations(monkeypatch):
    blobs, lock, finding, anchor, phase, digest = _frameshifted(monkeypatch)
    document = mini_recipe(blobs, lock, (V_LABEL, J7), decisions=[withhold([anchor, phase], digest)])
    result = build(lock, blobs, document)
    j = result.reference.alleles[-1]
    assert (j.label, j.anchor, j.coding_phase) == (J7, None, None)
    capabilities = result.reference.capabilities()
    assert capabilities.status("observed_span_simulation").available
    junction = capabilities.status("source_annotated_junction")
    assert not junction.available and f"J:{J7}" in {u.allele for u in junction.unmet}
    outcomes = {o.claim: o for o in result.report.claims}
    for claim in (anchor, phase):
        assert (outcomes[claim].outcome, outcomes[claim].decision) == ("withheld", 0)
    # neither source observation is shifted or replaced: both stay in the report as stated
    (raw_j,) = [r for r in records_in(blobs, "airr") if r["label"] == J7]
    kept = result.report.inventory
    assert kept.claim(anchor).value.get() == raw_j["j_cdr3_end"] - 1
    assert kept.claim(phase).value.get() == raw_j["j_codon_frame"] - 1
    assert finding in kept.findings


# ── relationships ────────────────────────────────────────────────────


def _expected_relationships(records_by_label, selected, pointer_of) -> list:
    """(kind, allele, runtime target, field pointer) for each alias or paralog a selected record
    states, in allele order; a name that is a selected allele of the same segment is that allele
    (a paralog only: an alias of an executable label is not representable), any other name a
    non-executable entity."""
    expected = []
    for segment, label in selected:
        record = records_by_label[label]
        for kind, field in (("alias", "aliases"), ("paralog", "paralogs")):
            for k, name in enumerate(record.get(field) or []):
                if name == label or (kind == "alias" and (segment, name) in selected):
                    continue
                target = f"{segment}:{name}" if (segment, name) in selected else f"entity:{name}"
                expected.append((kind, f"{segment}:{label}", target, f"{pointer_of(label)}/{field}/{k}"))
    return expected


def test_starter_relationships_are_the_selected_records_own_aliases_and_paralogs(starter):
    _, result = starter
    index, records = _airr_index(), raw_records(IGH, "airr")
    selected = [(segment, label) for segment in ("V", "D", "J") for label in SELECTION["selected"][segment]]
    expected = _expected_relationships({r["label"]: r for r in records}, selected,
                                       lambda label: f"{RECORDS}/{index[label]}")
    assert expected, "the starter's selected records state aliases and paralogs"
    catalogue = artifact(result.reference)["catalogue"]
    assert [(r["kind"], r["allele"], r["entity"]) for r in catalogue["relationships"]] == [e[:3] for e in expected]
    for relationship, (*_, pointer) in zip(catalogue["relationships"], expected):
        assert f"airr {raw_sha(raw(IGH, 'airr'))} {pointer}" == relationship["source"]
    entities = {e["name"]: e for e in catalogue["entities"]}
    assert set(entities) == {target.split(":", 1)[1] for _, _, target, _ in expected if target.startswith("entity:")}
    assert not any(e["executable"] for e in entities.values())
    assert len(result.reference.alleles) == len(selected)       # names add no executable allele


def _aliases(names):
    def edit(fmt, record):
        if record["label"] == "IGLV10-54*04":
            record["aliases"] = list(names)
    return edit


def test_a_relationship_the_runtime_cannot_hold_is_reported_on_the_record_not_projected(monkeypatch):
    names = (V_LABEL, "IGLV10-54*04", "IGLV10-54*i01")     # an executable label, its own label, an alias name
    blobs, lock = mini(monkeypatch, (V_LABEL, "IGLV10-54*04", J1), edit=_aliases(names))
    result = build(lock, blobs, mini_recipe(blobs, lock, (V_LABEL, "IGLV10-54*04", J1)))
    relationships = artifact(result.reference)["catalogue"]["relationships"]
    assert [(r["kind"], r["allele"], r["entity"]) for r in relationships] == [
        ("alias", "V:IGLV10-54*04", "entity:IGLV10-54*i01")]
    (outcome,) = [o for o in result.report.records if o.allele == "V:IGLV10-54*04"]
    for name in names[:2]:
        assert f"alias {name}" in outcome.detail
    assert "executable" in outcome.detail and "own label" in outcome.detail
    assert [a.label for a in result.reference.alleles] == [V_LABEL, "IGLV10-54*04", J1]


def test_a_paralog_name_without_a_selected_record_stays_a_named_entity_without_weight(monkeypatch):
    blobs, lock = mini(monkeypatch, (V_LABEL, "IGLJ2*01"))
    result = build(lock, blobs, mini_recipe(blobs, lock, (V_LABEL, "IGLJ2*01")))
    catalogue = artifact(result.reference)["catalogue"]
    assert [(r["kind"], r["allele"], r["entity"]) for r in catalogue["relationships"]] == [
        ("paralog", "J:IGLJ2*01", "entity:IGLJ3*01")]
    assert catalogue["entities"] == [{"id": "entity:IGLJ3*01", "name": "IGLJ3*01", "executable": False}]
    assert artifact(result.reference)["models"]["allele_usage"]["J"] == ["IGLJ2*01"]


# ── what the runtime cannot hold is refused in Python, before the native layer ──


def _coding(label, change):
    def edit(fmt, record):
        if record["label"] == label:
            record["coding_sequence"] = change(record["coding_sequence"])
    return edit


def test_a_selected_coding_sequence_of_65535_bases_builds(monkeypatch):
    blobs, lock = mini(monkeypatch, edit=_coding(J1, lambda s: s + "A" * (65535 - len(s))))
    (long,) = [r["coding_sequence"] for r in records_in(blobs, "airr") if r["label"] == J1]
    assert len(long) == 65535
    result = build(lock, blobs, mini_recipe(blobs, lock))
    assert result.fasta_of("J") == expected_fasta([(J1, long)])
    j = result.reference.alleles[-1]
    assert (j.label, j.length, j.sequence_sha256) == (J1, 65535, sha(long))


@pytest.mark.parametrize("change, problem", [
    (lambda s: s + "A" * (65536 - len(s)), "65536 bases"),
    (lambda s: s[:5] + "N" + s[6:], "'N' at position 5"),
    (lambda s: s.lower(), "at position 0"),
])
def test_a_selected_coding_sequence_the_runtime_cannot_hold_is_refused_before_native_entry(monkeypatch, change,
                                                                                          problem):
    blobs, lock = mini(monkeypatch, edit=_coding(J1, change))
    no_native(monkeypatch)
    document = recipe(lock, [entry(blobs, "airr", V_LABEL), entry(blobs, "airr", J1)])
    found = error("projection", build, lock, blobs, document)
    assert found.path == "/selection/1" and f"J:{J1}" in found.message and problem in found.message


def test_an_unselected_record_the_runtime_cannot_hold_is_a_structural_exclusion(monkeypatch):
    blobs, lock = mini(monkeypatch, (V_LABEL, J1, J7), edit=_coding(J7, lambda s: s[:5] + "N" + s[6:]))
    result = build(lock, blobs, mini_recipe(blobs, lock))
    outcomes = [o for o in result.report.records if result.report.inventory.record(o.record).label == J7]
    assert [o.exclusion for o in outcomes] == ["structural_incompatibility"] * 2
    assert all("'N' at position 5" in o.detail for o in outcomes)


PIPED = "x|IGLJ7*02|x|x|J-REGION"      # an IMGT-style header: read back as the second field


def _relabel(fmt, record):
    if record["label"] == J7:
        record["label"] = PIPED


def test_a_label_that_a_fasta_header_would_not_give_back_is_refused_or_excluded(monkeypatch):
    blobs, lock = mini(monkeypatch, (V_LABEL, J1, J7), edit=_relabel)
    no_native(monkeypatch)
    document = recipe(lock, [entry(blobs, "airr", V_LABEL), entry(blobs, "airr", PIPED)])
    found = error("projection", build, lock, blobs, document)
    assert found.path == "/selection/1" and PIPED in found.message and "'IGLJ7*02'" in found.message
    monkeypatch.undo()
    blobs, lock = mini(monkeypatch, (V_LABEL, J1, J7), edit=_relabel)
    result = build(lock, blobs, mini_recipe(blobs, lock))
    outcomes = [o for o in result.report.records if result.report.inventory.record(o.record).label == PIPED]
    assert [o.exclusion for o in outcomes] == ["structural_incompatibility"] * 2
    assert all("'IGLJ7*02'" in o.detail for o in outcomes)


# ── required capabilities, defaults and decisions ────────────────────


@pytest.mark.parametrize("capability, level, parts", [
    ("original_coding_ends", "assertion", (f"V:{V_LABEL}", f"J:{J1}", 'add "boundary_3"', 'add "boundary_5"')),
    ("canonical_productive", "source", ("available only at evidence assertion", "phase support: anchor_derived")),
    ("genotype_constrained", "assertion", ("not yet supported",)),
])
def test_an_unmet_required_capability_fails_the_build(monkeypatch, capability, level, parts):
    blobs, lock = mini(monkeypatch)
    document = mini_recipe(blobs, lock, requires=(("observed_span_simulation", "assertion"), (capability, level)))
    found = error("capability", build, lock, blobs, document)
    assert found.path == "/requires/1"
    for part in (capability, f"evidence {level}") + parts:
        assert part in found.message


def test_a_withheld_accepted_anchor_withholds_its_phase_and_drops_the_junction(monkeypatch):
    blobs, lock = mini(monkeypatch)
    inventory = interpreted(blobs, lock)
    anchor, phase = (claim_id(inventory, "airr", V_LABEL, p) for p in ("anchor", "coding_phase"))
    v_digest = entry(blobs, "airr", V_LABEL)["sequence_sha256"]
    decisions = [withhold([anchor], v_digest)]
    result = build(lock, blobs, mini_recipe(blobs, lock, decisions=decisions))
    v = result.reference.alleles[0]
    assert (v.label, v.anchor, v.coding_phase) == (V_LABEL, None, None)
    outcomes = {o.claim: o for o in result.report.claims}
    assert (outcomes[anchor].outcome, outcomes[anchor].decision) == ("withheld", 0)
    assert (outcomes[phase].outcome, outcomes[phase].decision) == ("withheld", None)
    assert "anchor" in outcomes[phase].reason
    status = result.reference.capabilities().status("source_annotated_junction")
    assert not status.available and f"V:{V_LABEL}" in {u.allele for u in status.unmet}
    document = mini_recipe(blobs, lock, decisions=decisions,
                           requires=(("source_annotated_junction", "assertion"),))
    found = error("capability", build, lock, blobs, document)
    assert found.path == "/requires/0" and f"V:{V_LABEL}" in found.message and "anchor" in found.message


def test_a_record_whose_own_claims_cannot_meet_a_requirement_is_excluded_for_missing_evidence(monkeypatch):
    def drop_gapped(fmt, header, sequence):           # IGLV1-40*01 then has no bound anchor
        return None if header == "IGLV1-40*01" else (header, sequence)
    labels = (V_LABEL, "IGLV3-32*01", "IGLV1-40*01", J1)
    blobs, lock = mini(monkeypatch, labels, gapped_edit=drop_gapped)
    requires = (("source_annotated_junction", "source"), ("canonical_productive", "assertion"))
    result = build(lock, blobs, mini_recipe(blobs, lock, requires=requires))
    details = {}
    for outcome in result.report.records:
        label = result.report.inventory.record(outcome.record).label
        if outcome.outcome == "excluded" and label != V_LABEL and label != J1:
            assert outcome.exclusion == "missing_required_evidence"
            details.setdefault(label, set()).add(outcome.detail)
    assert set(details) == {"IGLV3-32*01", "IGLV1-40*01"}
    (tyr,) = details["IGLV3-32*01"]                    # TAT at its anchor: not the conserved Cys
    assert "canonical_productive" in tyr and "conserved residue" in tyr and "source_annotated_junction" not in tyr
    (unbound,) = details["IGLV1-40*01"]
    assert "source_annotated_junction at evidence source: no accepted anchor" in unbound
    assert all(r.satisfied for r in result.report.requirements)


def test_the_built_defaults_profile_must_be_the_recipe_profile(monkeypatch):
    blobs, lock = mini(monkeypatch)
    # a second registered defaults profile, so that a recipe can name one the build does not use
    monkeypatch.setattr(common, "DEFAULTS_PROFILES", common.DEFAULTS_PROFILES + ("genairr.other_defaults.v1",))
    document = mini_recipe(blobs, lock)
    document["defaults_profile"] = "genairr.other_defaults.v1"
    found = error("projection", build, lock, blobs, document)
    assert found.path == "/defaults_profile" and "genairr.unfitted_defaults.v1" in found.message


def _set_j(label, field_name, value):
    def edit(fmt, record):
        if record["label"] == label:
            record[field_name] = value
    return edit


def test_an_accepted_claim_named_by_accept_claim_keeps_its_support_and_records_the_decision(monkeypatch):
    blobs, lock = mini(monkeypatch)
    inventory = interpreted(blobs, lock)
    anchor = claim_id(inventory, "airr", J1, "anchor")
    (index,) = [i for i, r in enumerate(records_in(blobs, "airr")) if r["label"] == J1]
    accept = {"action": "accept_claim", "findings": [], "claims": [anchor],
              "sequence_sha256": entry(blobs, "airr", J1)["sequence_sha256"], "reason": "reviewed",
              "evidence": [{"asset": "airr", "pointer": f"{RECORDS}/{index}/j_cdr3_end"}]}
    result = build(lock, blobs, mini_recipe(blobs, lock, decisions=[accept]))
    (outcome,) = [o for o in result.report.claims if o.claim == anchor]
    assert (outcome.outcome, outcome.decision) == ("projected", 0)
    assert result.reference.alleles[-1].anchor.support == "source"


def test_accept_claim_cannot_force_a_phase_without_a_projected_anchor(monkeypatch):
    blobs, lock = mini(monkeypatch, edit=_set_j(J1, "j_cdr3_end", None))
    no_native(monkeypatch)
    inventory = interpreted(blobs, lock)
    phase = claim_id(inventory, "airr", J1, "coding_phase")
    assert inventory.claim(phase).status == "candidate"
    (index,) = [i for i, r in enumerate(records_in(blobs, "airr")) if r["label"] == J1]
    accept = {"action": "accept_claim", "findings": [], "claims": [phase],
              "sequence_sha256": entry(blobs, "airr", J1)["sequence_sha256"], "reason": "reviewed",
              "evidence": [{"asset": "airr", "pointer": f"{RECORDS}/{index}/j_codon_frame"}]}
    found = error("projection", build, lock, blobs, mini_recipe(blobs, lock, decisions=[accept]))
    assert found.path == "/decisions/0" and f"J:{J1}" in found.message and "anchor" in found.message


def test_a_candidate_named_by_accept_claim_is_projected_as_an_assertion_citing_the_decision(monkeypatch):
    blobs, lock = mini(monkeypatch)
    real = evidence.interpret

    def demoted(lock_, assets, inventory):
        """The real interpretation with the V anchor demoted to a candidate (a value whose meaning a
        profile has not established), so a decision can accept it."""
        result = real(lock_, assets, inventory)
        claims = tuple(dataclasses.replace(c, status="candidate", support="none")
                       if c.predicate == "anchor" and result.subject_record(c).label == V_LABEL else c
                       for c in result.claims)
        fields = {name: getattr(result, name) for name in ("source_lock_hash", "adapter", "records", "candidates",
                                                           "mentioned_identities", "observations", "support",
                                                           "mappings", "findings")}
        return m.SourceInventory.assemble(claims=claims, **fields)

    monkeypatch.setattr(evidence, "interpret", demoted)
    inventory = demoted(lock, blobs, ogrdb.build_inventory(lock, blobs))
    anchor = claim_id(inventory, "airr", V_LABEL, "anchor")
    (index,) = [i for i, r in enumerate(records_in(blobs, "airr")) if r["label"] == V_LABEL]
    cdr3 = {"asset": "airr", "pointer": f"{RECORDS}/{index}/v_gene_delineations/0/cdr3_start"}
    accept = {"action": "accept_claim", "findings": [], "claims": [anchor],
              "sequence_sha256": entry(blobs, "airr", V_LABEL)["sequence_sha256"], "reason": "reviewed",
              "evidence": [cdr3]}
    without = build(lock, blobs, mini_recipe(blobs, lock))
    assert without.reference.alleles[0].anchor is None
    result = build(lock, blobs, mini_recipe(blobs, lock, decisions=[accept]))
    v = result.reference.alleles[0]
    assert (v.anchor.value, v.anchor.support) == (inventory.claim(anchor).value.get(), "assertion")
    assert v.anchor.source.startswith("recipe decision 0: airr ") and cdr3["pointer"] in v.anchor.source
    assert v.coding_phase.value == v.anchor.value % 3
    (outcome,) = [o for o in result.report.claims if o.claim == anchor]
    assert (outcome.outcome, outcome.decision) == ("projected", 0)
    status = result.reference.capabilities().status("source_annotated_junction")
    assert not status.available and f"V:{V_LABEL}" in {u.allele for u in status.unmet}


# ── determinism, line endings and Unicode labels ─────────────────────


def test_two_temporary_roots_give_identical_files_runtime_hashes_and_report(monkeypatch, tmp_path):
    blobs, lock = mini(monkeypatch)
    document = mini_recipe(blobs, lock)
    seen, real = [], tempfile.TemporaryDirectory

    class Recording(real):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            seen.append(self.name)

    monkeypatch.setattr(tempfile, "TemporaryDirectory", Recording)
    results = []
    for name in ("first", "a second root with a longer name"):
        root = tmp_path / name
        root.mkdir()
        monkeypatch.setattr(tempfile, "tempdir", str(root))
        results.append(build(lock, blobs, document))
        assert seen[-1].startswith(str(root)) and list(root.iterdir()) == []   # used, then removed
    first, second = results
    assert first.inputs == second.inputs
    assert first.report.runtime == second.report.runtime
    assert first.reference.to_json() == second.reference.to_json()
    assert m.report_hash(first.report) == m.report_hash(second.report)
    assert "a second root" not in second.reference.to_json() and str(tmp_path) not in second.reference.to_json()


def test_a_changed_recipe_changes_the_runtime_identity_although_the_dna_is_unchanged(monkeypatch):
    blobs, lock = mini(monkeypatch)
    first = build(lock, blobs, mini_recipe(blobs, lock))
    document = mini_recipe(blobs, lock)
    document["selection"][0]["reason"] = "the same record, selected for another stated reason"
    second = build(lock, blobs, document)
    assert first.build_input_hash != second.build_input_hash
    assert first.fasta == second.fasta
    assert [(a.segment, a.label, a.sequence_sha256) for a in first.reference.alleles] == \
        [(a.segment, a.label, a.sequence_sha256) for a in second.reference.alleles]
    assert first.annotations != second.annotations        # the notes cite the build inputs
    assert first.reference.catalogue_hash != second.reference.catalogue_hash
    assert first.reference.execution_hash != second.reference.execution_hash


def _line_endings(blobs, newline) -> dict:
    """The same exports with LF or CRLF line ends (record exports pretty-printed)."""
    changed = {}
    for fmt, data in blobs.items():
        text = data.decode("utf-8")
        if fmt in ("airr", "airr_ex"):
            text = json.dumps(json.loads(text), ensure_ascii=False, indent=1) + "\n"
        changed[fmt] = text.replace("\n", newline).encode("utf-8")
    return changed


def _annotation_values(data: bytes) -> list:
    values = []
    for allele in json.loads(data.decode("utf-8"))["alleles"]:
        annotations = {name: None if allele[name] is None else
                       (allele[name]["value"], allele[name]["support"], allele[name].get("derivation"))
                       for name in ("functionality", "anchor", "coding_phase", "boundary_5", "boundary_3")}
        values.append((allele["segment"], allele["label"], allele["sequence_sha256"], annotations))
    return values


def test_crlf_source_bytes_change_lock_and_citations_but_not_dna_labels_or_annotation_values(monkeypatch):
    base = small(IGL, set(MINI))
    built = {}
    for newline in ("\n", "\r\n"):
        blobs = _line_endings(base, newline)
        lock = seal(IGL, blobs)
        in_scope(monkeypatch, lock)
        built[newline] = (blobs, lock, build(lock, blobs, mini_recipe(blobs, lock)))
    (lf, lf_lock, a), (crlf, crlf_lock, b) = built["\n"], built["\r\n"]
    assert lf_lock.commitment != crlf_lock.commitment
    for asset in lf_lock.assets:
        assert asset.sha256 == raw_sha(lf[asset.asset_id]) != crlf_lock.asset(asset.asset_id).sha256
    assert a.fasta == b.fasta
    assert _annotation_values(a.annotations) == _annotation_values(b.annotations)
    assert a.annotations != b.annotations
    for blobs, result in ((lf, a), (crlf, b)):        # every citation names the digest of its own raw bytes
        for allele in json.loads(result.annotations.decode("utf-8"))["alleles"]:
            assert allele["anchor"]["source"].startswith(f"airr {raw_sha(blobs['airr'])} ")
            if allele["segment"] == "V":
                assert f"gapped {raw_sha(blobs['gapped'])} " in allele["anchor"]["source"]
            assert f"airr {raw_sha(blobs['airr'])} " in allele["notes"][0]
    for result in (a, b):
        for data in result.inputs.values():
            assert b"\r" not in data


UNICODE = "IGLV1-36*01_\u00fcn\u00efc\u00f6d\u00e9-\u03bb"     # letters beyond ASCII, valid labels


def test_unicode_labels_are_kept_literally_in_utf8_files_and_the_runtime(monkeypatch):
    def relabel(fmt, record):
        if record["label"] == V_LABEL:
            record["label"] = UNICODE

    def regap(fmt, header, sequence):
        return (UNICODE if header == V_LABEL else header), sequence

    blobs, lock = mini(monkeypatch, edit=relabel, gapped_edit=regap)
    document = recipe(lock, [entry(blobs, "airr", UNICODE), entry(blobs, "airr", J1)],
                      requires=(("source_annotated_junction", "source"),))
    result = build(lock, blobs, document)
    assert result.fasta_of("V").startswith(b">" + UNICODE.encode("utf-8") + b"\n")
    assert UNICODE.encode("utf-8") in result.annotations and b"\\u" not in result.annotations
    v = result.reference.alleles[0]
    assert v.label == UNICODE and v.anchor is not None and v.anchor.support == "source"
    assert Reference.from_json(result.reference.to_json()).alleles[0].label == UNICODE


# ── further refusals and withheld annotations ────────────────────────


@pytest.mark.parametrize("field_name, value", [
    ("species", {"id": "NCBITAXON:10090", "label": "Mus musculus"}), ("locus", "IGK"),
])
def test_a_record_declaring_another_taxon_or_locus_cannot_be_selected(monkeypatch, field_name, value):
    def edit(fmt, record):
        if record["label"] == V_LABEL:
            record[field_name] = value
    blobs, lock = mini(monkeypatch, edit=edit)
    no_native(monkeypatch)
    found = error("identity_mismatch", build, lock, blobs, mini_recipe(blobs, lock))
    assert found.path == "/selection/0/record"


def test_conflicting_claims_stay_null_at_runtime_and_keep_their_values_in_the_report(monkeypatch):
    blobs, lock = mini(monkeypatch)
    real = evidence.interpret

    def conflicting(lock_, assets, inventory):
        """The real interpretation with the V anchor replaced by two disagreeing claims."""
        result = real(lock_, assets, inventory)
        claims = []
        for claim in result.claims:
            if claim.predicate == "anchor" and result.subject_record(claim).label == V_LABEL:
                claims.append(dataclasses.replace(claim, status="conflicting"))
                claims.append(dataclasses.replace(claim, status="conflicting", value=claim.value.get() - 3))
            else:
                claims.append(claim)
        fields = {name: getattr(result, name) for name in ("source_lock_hash", "adapter", "records", "candidates",
                                                           "mentioned_identities", "observations", "support",
                                                           "mappings", "findings")}
        return m.SourceInventory.assemble(claims=tuple(claims), **fields)

    monkeypatch.setattr(evidence, "interpret", conflicting)
    result = build(lock, blobs, mini_recipe(blobs, lock))
    v = result.reference.alleles[0]
    assert (v.label, v.anchor, v.coding_phase) == (V_LABEL, None, None)
    inventory, outcomes = result.report.inventory, {o.claim: o for o in result.report.claims}
    disputed = [c for c in inventory.claims if c.status == "conflicting" and c.subject.asset == "airr"]
    assert len({c.value.get() for c in disputed}) == 2
    for claim in disputed:
        assert outcomes[claim.id].outcome == "withheld" and "conflicting" in outcomes[claim.id].reason


def test_the_starter_reports_honest_blockers_for_coding_ends_it_does_not_establish(starter):
    _, result = starter
    ends = result.report.capabilities.status("original_coding_ends")
    v, d, j = (SELECTION["counts"][s] for s in ("V", "D", "J"))
    assert not ends.available and len(ends.unmet) == v + j + 2 * d     # V 3', J 5', both D ends
    assert all('"complete"' in u.requirement and u.hint for u in ends.unmet)


def test_the_built_runtime_alleles_must_be_exactly_the_selection(monkeypatch):
    blobs, lock = mini(monkeypatch, (V_LABEL, J1, J7))
    other = build(lock, blobs, mini_recipe(blobs, lock, (V_LABEL, J1))).reference
    monkeypatch.setattr(Reference, "from_fasta", lambda **files: other)     # an importer returning another selection
    found = error("projection", build, lock, blobs, mini_recipe(blobs, lock, (V_LABEL, J7)))
    assert found.path == "/selection/1" and f"J:{J7}" in found.message

    def rejects(**files):
        raise ReferenceArtifactError("conflict", "alleles[0].anchor", "rejected for this test")
    monkeypatch.setattr(Reference, "from_fasta", rejects)
    found = error("projection", build, lock, blobs, mini_recipe(blobs, lock, (V_LABEL, J7)))
    assert "conflict" in found.message and "alleles[0].anchor" in found.message
