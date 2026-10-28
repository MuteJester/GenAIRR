"""Interpreting OGRDB coordinates: claims, findings and scope.

Helpers for asking what an export says about one allele's coordinates, and what the
interpreter concluded about it.
"""
from __future__ import annotations

import copy
import functools
import json
import re
from GenAIRR._reference_sources import consistency
from GenAIRR._reference_sources import evidence
from GenAIRR._reference_sources import model as m
from tests.support.fixtures.ogrdb_inventory import (
    ALL, EXPECT, IGH, IGK, IGL, ROOT, SETS, _key, build, encode, error, inventory_of, payload, pinned, raw,
    raw_records,
)

EVIDENCE = json.loads((ROOT / "interpretation_evidence.json").read_bytes().decode("utf-8"))

NOTE_LABELS = tuple(row["label"] for row in EXPECT["igl4_airr_ex_coding_end_observations"]
                    if row["finding"] == "coding_end_before_declared_v_rss")

NOTE_COLUMN = re.compile(r"up to and including (?:nucleotide|base) (\d+)")

NOTE_DOTS = re.compile(r"[Tt]railing [\"\u201c]\.+[\"\u201d]|\b(?:one|two|three) dots at the end|a dot at the 3' end")

@functools.lru_cache(maxsize=None)
def _raw_by_label(source, fmt) -> dict:
    found: dict = {}
    for record in raw_records(source, fmt):
        found.setdefault(record["label"], []).append(record)
    return found

def raw_record(source, fmt, label) -> dict:
    """One record of a raw export, read with the standard library JSON reader (a fresh copy)."""
    (record,) = _raw_by_label(source, fmt)[label]
    return copy.deepcopy(record)

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

def span_limiting_notes(source, fmt="airr") -> dict:
    """label -> affirmed columns stated in the raw curation note (empty for a trailing-dot note only),
    for every record of a raw export whose note limits its affirmed span (regex oracle)."""
    found = {}
    for record in raw_records(source, fmt):
        text = record.get("curation") or ""
        columns = sorted({int(n) for n in NOTE_COLUMN.findall(text)})
        if columns or NOTE_DOTS.search(text):
            found[record["label"]] = columns
    return found

GAPPED = {"airr": "gapped", "airr_ex": "gapped_ex"}

_INTERPRETED: dict = {}

def interpreted(source, formats=ALL):
    """(lock, interpreted inventory) of the pinned bytes, built once per module run."""
    key = (source, tuple(formats))
    if key not in _INTERPRETED:
        lock, inventory = inventory_of(source, formats)
        _INTERPRETED[key] = (lock, evidence.interpret(lock, pinned(source, formats), inventory))
    return _INTERPRETED[key]

def interpret(blobs, source=IGL):
    lock, inventory = build(source, blobs)
    return lock, evidence.interpret(lock, blobs, inventory)

def record_of(inventory, fmt, label):
    (record,) = [r for r in inventory.records if r.locator.asset == fmt and r.label == label]
    return record

def candidate_of(inventory, fmt, label):
    (candidate,) = inventory.candidates_of(record_of(inventory, fmt, label).locator)
    return candidate

def claim_of(inventory, subject, predicate):
    locator = subject.locator
    (claim,) = [x for x in inventory.claims if x.subject == locator and x.predicate == predicate]
    return claim

def claims(inventory, predicate, fmt=None, segment=None):
    found = []
    for claim in inventory.claims:
        if claim.predicate != predicate or (fmt is not None and claim.subject.asset != fmt):
            continue
        if segment is not None and inventory.subject_record(claim).segment != segment:
            continue
        found.append(claim)
    return found

def findings(inventory, code) -> list:
    return [f for f in inventory.findings if f.code == code]

def at(record, *tokens) -> m.Locator:
    return m.Locator(record.locator.asset, "/".join((record.locator.pointer,) + tuple(str(t) for t in tokens)))

def small(source, labels, *, formats=ALL, edit=None, gapped_edit=None) -> dict:
    """The pinned exports restricted to ``labels`` (JSON records and gapped records), then
    edited: ``edit(fmt, record)`` for each kept JSON record, ``gapped_edit(fmt, header, sequence)``
    returning ``(header, sequence)`` or ``None`` (drop) for each kept gapped record."""
    blobs = {}
    for fmt in formats:
        if fmt in GAPPED:
            document = payload(source, fmt)
            germline_set = document["GermlineSet"][0]
            germline_set["allele_descriptions"] = [r for r in germline_set["allele_descriptions"]
                                                   if r["label"] in labels]
            for record in germline_set["allele_descriptions"]:
                if edit is not None:
                    edit(fmt, record)
            blobs[fmt] = encode(document)
            continue
        lines = []
        for chunk in raw(source, fmt).decode("utf-8").split(">")[1:]:
            header, _, body = chunk.partition("\n")
            header = header.rstrip("\r")
            if header not in labels:
                continue
            kept = (header, "".join(line.rstrip("\r") for line in body.split("\n")))
            if gapped_edit is not None:
                kept = gapped_edit(fmt, *kept)
            if kept is not None:
                lines.append(f">{kept[0]}\n{kept[1]}\n")
        blobs[fmt] = "".join(lines).encode("utf-8")
    return blobs

def in_scope(monkeypatch, lock) -> None:
    """Register this (synthetic) lock's exact bytes as validated, as a reviewed profile extension would."""
    pinned_assets = tuple((a.format, a.sha256, a.bytes) for a in lock.assets
                          if a.role in ("record_export", "gapped_alignment"))
    key = ("ogrdb_api_v2_human_ig_observed_v1", 1)
    release = (lock.identity.taxon.taxonomy_id, lock.set_id, lock.release)
    scopes = dict(consistency.PROFILE_SCOPES)
    scopes[key] = tuple(e for e in scopes[key] if e[:3] != release) + (release + (pinned_assets,),)
    monkeypatch.setattr(consistency, "PROFILE_SCOPES", scopes)
    assert m.profile_scope_problem(lock) is None

def interpret_in_scope(monkeypatch, blobs, source=IGL):
    lock, inventory = build(source, blobs)
    in_scope(monkeypatch, lock)
    return lock, evidence.interpret(lock, blobs, inventory)

def assert_status(claim, status, support, derivation, value=None):
    assert (claim.status, claim.support, claim.derivation) == (status, support, derivation), claim
    if value is not None:
        assert claim.value.get() == value

V_LABEL = "IGLV1-36*01"
