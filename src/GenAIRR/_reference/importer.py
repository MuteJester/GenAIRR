"""Local FASTA import: turn V/D/J FASTA files plus an annotation sidecar into a catalogue.

Every FASTA record is accounted for (accepted, duplicate or rejected with a reason).
Labels are never renamed.  Claims come only from the sidecar, from a declared IMGT
alignment, or from a unique motif reported as an unaccepted candidate.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Optional

from GenAIRR._artifact_io import ArtifactIOError, read_text_bounded
from GenAIRR._reference.errors import ReferenceArtifactError
from GenAIRR._reference.model import (
    Allele, Catalogue, Claim, Entity, Finding, RecordEntry, ReferenceHeader, Relationship,
    SourceFile, Taxon, Transform, check_label, sequence_digest,
)
from GenAIRR._reference.sidecar import ANNOTATIONS, MAX_BYTES as SIDECAR_MAX_BYTES, load_sidecar
from GenAIRR.utilities.imgt_regions import imgt_region_boundaries

FASTA_MAX_BYTES = 256 * 1024 * 1024
MAX_SEQUENCE = 65535
_BASES = frozenset("ACGT")
_ALIAS_EXECUTABLE = ("alias names an executable record; keep one record and declare the other "
                     "name as an alias entity")


def _text_digest(text: str) -> str:
    return "sha256:" + hashlib.sha256(text.encode("utf-8")).hexdigest()


def _read(path, max_bytes):
    try:
        return read_text_bounded(path, max_bytes=max_bytes)
    except ArtifactIOError as exc:
        kind = "io" if exc.kind in ("io", "resource") else "malformed"
        raise ReferenceArtifactError(kind, Path(path).name, str(exc)) from exc


def _label_of(header: str) -> str:
    fields = header.split("|")
    if len(fields) >= 5 and fields[4].strip().endswith("-REGION"):
        return fields[1].strip()
    tokens = header.split()
    return tokens[0] if tokens else ""


def _parse_fasta(text: str, name: str):
    """Yield ``(header, raw_sequence)`` pairs in file order."""
    records = []
    header = None
    parts: list = []
    for line in text.split("\n"):
        line = line.rstrip("\r")
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(parts)))
            header, parts = line[1:], []
        elif line.strip():
            if header is None:
                raise ReferenceArtifactError("malformed", name, "sequence data before the first '>' header")
            parts.append("".join(line.split()))
    if header is not None:
        records.append((header, "".join(parts)))
    return records


def _normalize(raw: str, segment: str, v_alignment: str):
    """Return ``(sequence, gapped, transforms, rejection_reason)``."""
    if raw == "":
        return None, None, (), "empty_sequence"
    upper = raw.upper()
    transforms = []
    if upper != raw:
        transforms.append(Transform("uppercase", _text_digest(raw), _text_digest(upper)))
    gapped = None
    ungapped = upper
    if "." in upper:
        if segment != "V" or v_alignment != "imgt_gapped":
            return None, None, (), "gap_without_declared_alignment"
        gapped = upper
        ungapped = upper.replace(".", "")
        transforms.append(Transform("remove_imgt_gaps", _text_digest(upper), _text_digest(ungapped)))
    if ungapped == "":
        return None, None, (), "empty_sequence"
    for position, char in enumerate(upper):
        if char != "." and char not in _BASES:
            return None, None, (), f"unsupported_base: {char!r} at position {position}"
    if len(ungapped) > MAX_SEQUENCE:
        return None, None, (), f"too_long: {len(ungapped)} bases exceeds {MAX_SEQUENCE}"
    return ungapped, gapped, tuple(transforms), None


def _scan_segment(path, segment: str, v_alignment: str):
    """Classify every record of one FASTA file.

    Returns ``(source, entries)`` where each entry is a dict with the record fields
    plus ``sequence``, ``gapped`` and ``transforms`` for accepted/duplicate records.
    """
    text = _read(path, FASTA_MAX_BYTES)
    name = Path(path).name
    source = SourceFile(f"src:{segment}", f"{segment.lower()}_fasta", name, _text_digest(text),
                        len(text.encode("utf-8")))
    entries = []
    for index, (header, raw) in enumerate(_parse_fasta(text, name)):
        label = _label_of(header)
        entry = {"id": f"rec:{segment}:{index}", "index": index, "header": header, "label": label,
                 "status": "rejected", "reason": None, "sequence": None, "gapped": None,
                 "transforms": (), "duplicate_of": None}
        problem = check_label(label)
        if problem:
            entry["label"] = None
            entry["reason"] = f"unrepresentable_label: {problem}"
        else:
            seq, gapped, transforms, reason = _normalize(raw, segment, v_alignment)
            if reason:
                entry["reason"] = reason
            else:
                entry.update(sequence=seq, gapped=gapped, transforms=transforms, status="accepted")
        entries.append(entry)

    by_label: dict = {}
    for entry in entries:
        if entry["status"] == "accepted":
            by_label.setdefault(entry["label"], []).append(entry)
    for label, group in by_label.items():
        if len({(e["sequence"], e["gapped"]) for e in group}) > 1:
            for e in group:
                e.update(status="rejected", reason="conflicting_sequence_for_label")
                e["sequence"] = e["gapped"] = None
                e["transforms"] = ()
        else:
            for e in group[1:]:
                e.update(status="duplicate", duplicate_of=group[0]["id"])
    return source, entries


def _resolve(sequence: str, segment: str, locus: str, gapped: Optional[str], imgt_gapped: bool) -> dict:
    from GenAIRR import _engine
    return json.loads(_engine.resolve_reference_anchor(
        sequence.encode("ascii"), segment, locus=locus,
        gapped_sequence=gapped.encode("ascii") if gapped is not None else None,
        imgt_gapped=imgt_gapped))


def _accepted(entries):
    return [e for e in entries if e["status"] == "accepted"]


def _motif_suggestions(segment_entries: dict, locus: str) -> list:
    suggestions = []
    for entry in _accepted(segment_entries.get("J", [])):
        result = _resolve(entry["sequence"], "J", locus, None, False)
        if result["position"] is not None:
            suggestions.append({
                "allele": f"J:{entry['label']}", "predicate": "anchor", "value": result["position"],
                "derivation": "unique_motif",
                "note": "unique F/W-G-X-G motif position; candidate only, verify before stating it",
            })
    return suggestions


def scan_fasta(*, v, j, d, v_alignment: str, locus: str):
    """Scan FASTA files; return accepted records per segment and motif-candidate suggestions."""
    segment_entries = {}
    for segment, path in (("V", v), ("D", d), ("J", j)):
        if path is None:
            continue
        segment_entries[segment] = _scan_segment(path, segment, v_alignment)[1]
    records = {segment: [{"label": e["label"], "sequence_sha256": sequence_digest(e["sequence"])}
                         for e in _accepted(entries)]
               for segment, entries in segment_entries.items()}
    for segment in ("V", "D", "J"):
        records.setdefault(segment, [])
    return records, _motif_suggestions(segment_entries, locus)


def _gene_of(label: str, declared: Optional[str]):
    if declared is not None:
        return declared, "declared"
    if "*" in label:
        prefix = label.split("*", 1)[0]
        if prefix and check_label(prefix) is None:
            return prefix, "label_prefix"
    return label, "label"


def _stated_claims(allele_id: str, side_entry: dict) -> list:
    claims = []
    for name in ANNOTATIONS:
        ann = side_entry[name]
        if ann is None:
            continue
        derivation = "stated"
        if name == "coding_phase" and ann["derivation"] == "anchor_derived":
            derivation = "anchor_derived"
        claims.append(Claim(allele_id, name, ann["value"], ann["support"], derivation, "accepted",
                            ann["source"], ann["method"]))
    return claims


def _alignment_claims(allele: Allele, locus: str) -> list:
    result = _resolve(allele.sequence, "V", locus, allele.gapped_sequence, True)
    claims = []
    if result["position"] is not None:
        claims.append(Claim(allele.id, "anchor", result["position"], "assertion", "imgt_alignment",
                            "accepted", None, result["method"]))
    spans = {label: [start, end] for label, (start, end)
             in imgt_region_boundaries(allele.gapped_sequence).items() if start < end}
    if spans:
        claims.append(Claim(allele.id, "v_subregions", spans, "assertion", "imgt_alignment",
                            "accepted", None, "imgt_gapped_boundaries"))
    return claims


def _mark_anchor_conflicts(claims: list) -> list:
    by_subject: dict = {}
    for claim in claims:
        if claim.predicate == "anchor" and claim.derivation in ("stated", "imgt_alignment"):
            by_subject.setdefault(claim.subject, []).append(claim)
    clash = set()
    for subject, group in by_subject.items():
        if len({c.value for c in group}) > 1:
            clash.update(id(c) for c in group)
    return [Claim(c.subject, c.predicate, c.value, c.support, c.derivation, "conflict", c.source, c.method)
            if id(c) in clash else c for c in claims]


def import_local(*, v, j, d, annotations) -> Catalogue:
    """Import local FASTA files and an annotation sidecar into a validated catalogue."""
    side = load_sidecar(_read(annotations, SIDECAR_MAX_BYTES))
    ref = side["reference"]
    conventions = side["conventions"]
    v_alignment = conventions["v_alignment"]
    locus = ref["locus"]
    chain_type = ref["chain_type"]
    if chain_type == "vj" and d is not None:
        raise ReferenceArtifactError("malformed", "d", "a D FASTA file is not allowed for chain_type vj")
    if chain_type == "vdj" and d is None:
        raise ReferenceArtifactError("malformed", "d", "chain_type vdj requires a D FASTA file")

    sources, records, alleles = [], [], []
    for segment, path in (("V", v), ("D", d), ("J", j)):
        if path is None:
            continue
        source, entries = _scan_segment(path, segment, v_alignment)
        sources.append(source)
        allele_records: dict = {}
        for e in entries:
            allele_id = None
            if e["status"] in ("accepted", "duplicate"):
                allele_id = f"{segment}:{e['label']}"
                allele_records.setdefault(allele_id, []).append(e)
            records.append(RecordEntry(e["id"], source.id, e["index"], e["header"], e["label"], e["status"],
                                       e["reason"], allele_id, e["duplicate_of"]))
        for allele_id, group in allele_records.items():
            first = group[0]
            gene, basis = _gene_of(first["label"], None)
            alleles.append(Allele(allele_id, segment, first["label"], gene, basis, first["sequence"],
                                  sequence_digest(first["sequence"]), first["gapped"],
                                  first["transforms"], tuple(e["id"] for e in group)))

    by_id = {a.id: a for a in alleles}
    claims, replaced = [], {}
    for i, side_entry in enumerate(side["alleles"]):
        allele_id = f"{side_entry['segment']}:{side_entry['label']}"
        path = f"alleles[{i}]"
        allele = by_id.get(allele_id)
        if allele is None:
            raise ReferenceArtifactError("conflict", path, f"sidecar entry {allele_id} has no accepted "
                                         "allele in the FASTA files")
        if allele.sequence_sha256 != side_entry["sequence_sha256"]:
            raise ReferenceArtifactError(
                "conflict", f"{path}.sequence_sha256",
                f"sequence digest mismatch for {allele_id}: sidecar has {side_entry['sequence_sha256']}, "
                f"FASTA gives {allele.sequence_sha256}; digest is over the uppercase ungapped A/C/G/T sequence")
        if side_entry["gene"] is not None:
            gene, basis = _gene_of(allele.label, side_entry["gene"])
            replaced[allele_id] = Allele(allele.id, allele.segment, allele.label, gene, basis,
                                         allele.sequence, allele.sequence_sha256, allele.gapped_sequence,
                                         allele.transforms, allele.records)
        claims.extend(_stated_claims(allele_id, side_entry))
    alleles = [replaced.get(a.id, a) for a in alleles]

    for allele in alleles:
        if allele.segment == "V" and allele.gapped_sequence is not None:
            claims.extend(_alignment_claims(allele, locus))
        elif allele.segment == "J":
            result = _resolve(allele.sequence, "J", locus, None, False)
            if result["position"] is not None:
                claims.append(Claim(allele.id, "anchor", result["position"], "none", "unique_motif",
                                    "candidate", None, result["method"]))
    claims = _mark_anchor_conflicts(claims)

    labels = {(a.segment, a.label) for a in alleles}
    entities, relationships = {}, []
    for i, rel in enumerate(side["relationships"]):
        segment, _, label = rel["allele"].partition(":")
        if (segment, label) not in labels:
            raise ReferenceArtifactError("conflict", f"relationships[{i}].allele",
                                         f"relationship names {rel['allele']}, which has no accepted allele")
        if rel["kind"] == "alias" and rel["entity"] != label and (segment, rel["entity"]) in labels:
            raise ReferenceArtifactError("conflict", f"relationships[{i}].entity", _ALIAS_EXECUTABLE)
        entity_id = f"entity:{rel['entity']}"
        entities.setdefault(entity_id, Entity(entity_id, rel["entity"], False))
        relationships.append(Relationship(rel["kind"], rel["allele"], entity_id, rel["source"]))

    findings = []
    for segment in ("V", "D", "J"):
        groups: dict = {}
        for a in alleles:
            if a.segment == segment:
                groups.setdefault(a.sequence, []).append(a.id)
        for ids in groups.values():
            if len(ids) > 1:
                findings.append(Finding("identical_sequence_distinct_labels", tuple(ids),
                                        f"{len(ids)} alleles share one sequence under different labels"))

    header = ReferenceHeader(ref["name"], ref["version"],
                             Taxon(ref["taxon"]["scientific_name"], ref["taxon"]["taxonomy_id"]),
                             locus, chain_type, ref["build_date"])
    catalogue = Catalogue(header, dict(conventions), tuple(sources), tuple(records), tuple(alleles),
                          tuple(claims), tuple(entities.values()), tuple(relationships), tuple(findings))
    catalogue.validate()
    return catalogue
