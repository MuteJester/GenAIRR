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
import hashlib
import json
import re

import pytest

import GenAIRR
from GenAIRR._reference.artifact import Reference
from GenAIRR._reference_sources import evidence, ogrdb, projection
from GenAIRR._reference_sources import model as m
from tests.test_ogrdb_source_coordinates import V_LABEL, in_scope, small
from tests.test_ogrdb_source_inventory import (
    EXPECT, IGH, IGL, RECORDS, SELECTION, error, pinned, raw, raw_records, seal,
)

J1 = "IGLJ1*01"
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

    def withhold(claims, digest):
        return {"action": "withhold_claim", "findings": [], "claims": claims, "sequence_sha256": digest,
                "reason": "test decision", "evidence": []}

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
