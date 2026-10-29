"""The pinned OGRDB export fixtures, and the inventory built from them.

Locus identifiers, export formats, the recorded expectations, and the helpers that turn
pinned bytes into a SourceInventory. The base of the reference-source fixture chain: the
coordinate, projection and API fixtures all build on it.
"""
from __future__ import annotations

from tests.support import paths

import hashlib
import json
import pytest
from GenAIRR._reference_sources import model as m
from GenAIRR._reference_sources import ogrdb
from GenAIRR._reference_sources.errors import ReferenceSourceError

ROOT = paths.data("reference_sources")

SOURCES = json.loads((ROOT / "sources.json").read_bytes().decode("utf-8"))

EXPECT = json.loads((ROOT / "expectations.json").read_bytes().decode("utf-8"))

SELECTION = json.loads((ROOT / "human_igh_release10_candidate_selection.json").read_bytes().decode("utf-8"))

IGH, IGK, IGL = ("9606.IGH_VDJ", 10), ("9606.IGKappa_VJ", 5), ("9606.IGLambda_VJ", 4)

SETS = (IGH, IGK, IGL)

ALL = ("airr", "airr_ex", "gapped", "gapped_ex")

SUFFIX = {"airr": "json", "airr_ex": "json", "gapped": "fasta", "gapped_ex": "fasta"}

CHAIN = {"IGH": "vdj", "IGK": "vj", "IGL": "vj"}

RECORDS = "/GermlineSet/0/allele_descriptions"

ADAPTER = m.AdapterRef("ogrdb_germline_set", 1, m.CoordinateProfileRef("ogrdb_api_v2_human_ig_observed_v1", 1))

def _key(source) -> str:
    return f"{source[0]}/{source[1]}"

def _entry(source) -> dict:
    (entry,) = [s for s in SOURCES["sets"] if (s["set_id"], s["release"]) == source]
    return entry

def raw(source, fmt) -> bytes:
    set_id, release = source
    return (ROOT / "ogrdb" / set_id / str(release) / f"{fmt}.{SUFFIX[fmt]}").read_bytes()

def payload(source, fmt) -> dict:
    """The raw export, read independently with the standard library."""
    return json.loads(raw(source, fmt).decode("utf-8"))

def raw_records(source, fmt) -> list:
    return payload(source, fmt)["GermlineSet"][0]["allele_descriptions"]

def encode(document) -> bytes:
    return json.dumps(document, ensure_ascii=False).encode("utf-8")

def _sha(data: bytes) -> str:
    return "sha256:" + hashlib.sha256(data).hexdigest()

def expected_identity(source) -> m.SourceIdentity:
    entry = _entry(source)
    taxonomy_id = int(entry["species"]["id"].split(":")[1])
    return m.SourceIdentity(m.Taxon(entry["species"]["label"], taxonomy_id), entry["locus"], CHAIN[entry["locus"]])

def seal(source, blobs: dict, *, identity=None, set_id=None, release=None) -> m.SourceLock:
    """A local-file lock over exactly ``blobs`` (format -> bytes)."""
    assets = tuple(m.Asset(asset_id=f, role=m.FORMAT_ROLES[f], format=f, sha256=_sha(blobs[f]),
                           bytes=len(blobs[f]), declared_url=None, origin_kind="local_file")
                   for f in m.FORMATS if f in blobs)
    missing = tuple(m.MissingAsset(format=f, role=m.FORMAT_ROLES[f], reason="not_supplied", detail=None)
                    for f in m.EVIDENCE_FORMATS if f not in blobs)
    terms = m.Terms(status="unresolved", assets=(), licence=None, citation=None,
                    reason="no data-terms statement was supplied")
    return m.SourceLock.seal(provider="ogrdb", set_id=set_id or source[0], release=release or source[1],
                             identity=identity or expected_identity(source), adapter=ADAPTER, assets=assets,
                             missing_assets=missing, terms=terms)

def pinned(source, formats=ALL) -> dict:
    return {f: raw(source, f) for f in formats}

_INVENTORIES: dict = {}

def inventory_of(source, formats=ALL):
    """(lock, inventory) over the pinned bytes, built once per module run."""
    key = (source, tuple(formats))
    if key not in _INVENTORIES:
        blobs = pinned(source, formats)
        lock = seal(source, blobs)
        _INVENTORIES[key] = (lock, ogrdb.build_inventory(lock, blobs))
    return _INVENTORIES[key]

def build(source, blobs):
    lock = seal(source, blobs)
    return lock, ogrdb.build_inventory(lock, blobs)

def error(kind, call, *args, **kwargs) -> ReferenceSourceError:
    with pytest.raises(ReferenceSourceError) as info:
        call(*args, **kwargs)
    assert info.value.kind == kind, f"expected {kind!r}, got {info.value.kind!r}: {info.value}"
    return info.value
