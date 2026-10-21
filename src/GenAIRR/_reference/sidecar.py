"""Annotation sidecar v1: a user-written JSON file of per-allele annotations."""
from __future__ import annotations

import copy
import datetime
import re
from typing import Any

from GenAIRR._artifact_io import ArtifactIOError, strict_loads
from GenAIRR._reference.errors import ReferenceArtifactError
from GenAIRR._reference.model import (
    CHAIN_TYPES, DERIVATIONS, PREDICATES, SEGMENTS, V_ALIGNMENTS,
    _bad, _enum, _int, _list, _obj, _str, check_claim_value, check_label,
)

SCHEMA_NAME = "genairr_reference_annotations"
SCHEMA_VERSION = 1
MAX_BYTES = 16 * 1024 * 1024
ANNOTATIONS = ("functionality", "anchor", "coding_phase", "boundary_5", "boundary_3")
_DIGEST_RE = re.compile(r"sha256:[0-9a-f]{64}")
_ANNOTATION_SUPPORTS = ("source", "assertion")
_PHASE_DERIVATIONS = ("independent", "anchor_derived")


def _annotation(value, name, path):
    if value is None:
        return
    keys = ("value", "support", "source", "method")
    if name == "coding_phase":
        keys += ("derivation",)
    _obj(value, keys, path)
    check_claim_value(name, value["value"], None, f"{path}.value")
    support = _enum(value["support"], _ANNOTATION_SUPPORTS, f"{path}.support")
    source = _str(value["source"], f"{path}.source", optional=True)
    _str(value["method"], f"{path}.method", optional=True)
    if support == "source" and (source is None or not source.strip()):
        raise _bad(f"{path}.source", "support \"source\" requires a non-blank source")
    if name == "coding_phase":
        _enum(value["derivation"], _PHASE_DERIVATIONS, f"{path}.derivation")


def _reference(value, path="reference"):
    _obj(value, ("name", "version", "locus", "chain_type", "taxon", "build_date"), path)
    _str(value["name"], f"{path}.name", nonblank=True)
    _str(value["version"], f"{path}.version", optional=True)
    _str(value["locus"], f"{path}.locus", nonblank=True)
    chain_type = _enum(value["chain_type"], CHAIN_TYPES, f"{path}.chain_type")
    taxon = _obj(value["taxon"], ("scientific_name", "taxonomy_id"), f"{path}.taxon")
    _str(taxon["scientific_name"], f"{path}.taxon.scientific_name", nonblank=True)
    _int(taxon["taxonomy_id"], f"{path}.taxon.taxonomy_id", optional=True, minimum=0)
    build_date = _str(value["build_date"], f"{path}.build_date", optional=True)
    if build_date is not None:
        try:
            if not re.fullmatch(r"\d{4}-\d{2}-\d{2}", build_date):
                raise ValueError
            datetime.date.fromisoformat(build_date)
        except ValueError:
            raise _bad(f"{path}.build_date", "expected YYYY-MM-DD or null") from None
    return chain_type


def _conventions(value, path="conventions"):
    _obj(value, ("coordinates", "v_alignment"), path)
    if value["coordinates"] != "zero_based_half_open":
        raise _bad(f"{path}.coordinates", "only zero_based_half_open is supported", "unsupported_convention")
    _enum(value["v_alignment"], V_ALIGNMENTS, f"{path}.v_alignment")


def _allele_id(value, path):
    _str(value, path)
    segment, _, label = value.partition(":")
    if segment not in SEGMENTS or check_label(label) is not None:
        raise _bad(path, "expected an allele id of the form <V|D|J>:<label>")


def load_sidecar(text: str) -> dict:
    """Parse and validate a sidecar; return the document with optional lists defaulted."""
    try:
        doc = strict_loads(text, max_bytes=MAX_BYTES)
    except ArtifactIOError as exc:
        raise ReferenceArtifactError("malformed", "", str(exc)) from exc
    if not isinstance(doc, dict):
        raise _bad("", "expected a JSON object")
    if doc.get("schema") != SCHEMA_NAME:
        raise _bad("schema", f"expected {SCHEMA_NAME!r}")
    version = doc.get("schema_version")
    if isinstance(version, bool) or not isinstance(version, int):
        raise _bad("schema_version", "expected an integer")
    if version != SCHEMA_VERSION:
        raise _bad("schema_version", f"unsupported schema version {version}", "unsupported_schema")
    unknown = [k for k in doc if k not in
               ("schema", "schema_version", "reference", "conventions", "alleles", "relationships", "suggestions")]
    missing = [k for k in ("reference", "conventions", "alleles") if k not in doc]
    if unknown or missing:
        raise _bad("", f"unknown keys {unknown}, missing keys {missing}")
    chain_type = _reference(doc["reference"])
    _conventions(doc["conventions"])

    seen = set()
    for i, entry in enumerate(_list(doc["alleles"], "alleles")):
        path = f"alleles[{i}]"
        _obj(entry, ("segment", "label", "sequence_sha256", "gene") + ANNOTATIONS + ("notes",), path)
        segment = _enum(entry["segment"], SEGMENTS, f"{path}.segment")
        if segment == "D" and chain_type == "vj":
            raise _bad(f"{path}.segment", "D alleles are not allowed for chain_type vj")
        problem = check_label(entry["label"])
        if problem:
            raise _bad(f"{path}.label", problem)
        if not _DIGEST_RE.fullmatch(_str(entry["sequence_sha256"], f"{path}.sequence_sha256")):
            raise _bad(f"{path}.sequence_sha256", "expected sha256:<64 hex>")
        if entry["gene"] is not None:
            problem = check_label(entry["gene"])
            if problem:
                raise _bad(f"{path}.gene", problem)
        for name in ANNOTATIONS:
            _annotation(entry[name], name, f"{path}.{name}")
        for j, note in enumerate(_list(entry["notes"], f"{path}.notes")):
            _str(note, f"{path}.notes[{j}]")
        key = (segment, entry["label"])
        if key in seen:
            raise _bad(path, f"duplicate allele entry {segment}:{entry['label']}")
        seen.add(key)

    for i, rel in enumerate(_list(doc.get("relationships", []), "relationships")):
        path = f"relationships[{i}]"
        _obj(rel, ("kind", "allele", "entity", "source"), path)
        _enum(rel["kind"], ("alias", "paralog"), f"{path}.kind")
        _allele_id(rel["allele"], f"{path}.allele")
        _str(rel["entity"], f"{path}.entity", nonblank=True)
        _str(rel["source"], f"{path}.source", optional=True)

    for i, sug in enumerate(_list(doc.get("suggestions", []), "suggestions")):
        path = f"suggestions[{i}]"
        _obj(sug, ("allele", "predicate", "value", "derivation", "note"), path)
        _allele_id(sug["allele"], f"{path}.allele")
        predicate = _enum(sug["predicate"], PREDICATES, f"{path}.predicate")
        check_claim_value(predicate, sug["value"], None, f"{path}.value")
        _enum(sug["derivation"], DERIVATIONS, f"{path}.derivation")
        _str(sug["note"], f"{path}.note")

    result = dict(doc)
    result.setdefault("relationships", [])
    result.setdefault("suggestions", [])
    return result


def annotation_template(*, records_by_segment: dict, header: dict, conventions: dict, suggestions: list) -> dict:
    """Build a sidecar with every annotation unknown (``null``) and no support claimed."""
    alleles = []
    for segment in SEGMENTS:
        for record in records_by_segment.get(segment) or []:
            entry: dict[str, Any] = {
                "segment": segment, "label": record["label"],
                "sequence_sha256": record["sequence_sha256"], "gene": None,
            }
            for name in ANNOTATIONS:
                entry[name] = None
            entry["notes"] = []
            alleles.append(entry)
    return {
        "schema": SCHEMA_NAME,
        "schema_version": SCHEMA_VERSION,
        "reference": copy.deepcopy(header),
        "conventions": copy.deepcopy(conventions),
        "alleles": alleles,
        "relationships": [],
        "suggestions": copy.deepcopy(suggestions),
    }
