"""Building an OgrdbSource over the pinned exports.

Sources recorded from local files, the full and mini recipes, and helpers for reading a
built reference back. MINI itself lives with the projection fixtures, which define it.
"""
from __future__ import annotations

from tests.support import paths

import hashlib
import http.client
import json
import pathlib
import socket
import pytest
from GenAIRR._reference_sources import api
from GenAIRR._reference_sources.acquisition import HttpsTransport
from GenAIRR._reference_sources.errors import KINDS, ReferenceSourceError
from GenAIRR.reference_sources import OgrdbSource, ReferenceBuild
from tests.support.fixtures.ogrdb_coordinates import in_scope, small
from tests.support.fixtures.source_projection import MINI, entry, recipe

DATA = paths.data("reference_sources")

SOURCES = json.loads((DATA / "sources.json").read_bytes().decode("utf-8"))

SUFFIX = {"airr": "json", "airr_ex": "json", "gapped": "fasta", "gapped_ex": "fasta"}

PROFILE = "ogrdb_api_v2_human_ig_observed_v1"

ALL = ("airr", "airr_ex", "gapped", "gapped_ex")

IGH, IGK, IGL = ("9606.IGH_VDJ", 10), ("9606.IGKappa_VJ", 5), ("9606.IGLambda_VJ", 4)

def fixture_path(source, fmt) -> pathlib.Path:
    return DATA / "ogrdb" / source[0] / str(source[1]) / f"{fmt}.{SUFFIX[fmt]}"

def fixture_bytes(source, fmt) -> bytes:
    return fixture_path(source, fmt).read_bytes()

def pinned_asset(source, fmt) -> dict:
    (entry_,) = [s for s in SOURCES["sets"] if (s["set_id"], s["release"]) == source]
    (asset,) = [a for a in entry_["assets"] if a["format"] == fmt]
    return asset

def sha(data: bytes) -> str:
    return "sha256:" + hashlib.sha256(data).hexdigest()

def from_files(tmp_path, source=IGL, formats=ALL, *, cache="cache", files=None) -> OgrdbSource:
    exports = files if files is not None else {fmt: fixture_path(source, fmt) for fmt in formats}
    return OgrdbSource.from_files(set_id=source[0], release=source[1], profile=PROFILE, exports=exports,
                                  cache=tmp_path / cache)

def written(tmp_path, blobs: dict, folder="snapshot") -> dict:
    """``blobs`` (format -> bytes) written as local files; returns ``{format: path}``."""
    directory = tmp_path / folder
    directory.mkdir(parents=True, exist_ok=True)
    paths = {}
    for fmt, data in blobs.items():
        paths[fmt] = directory / f"{fmt}.{SUFFIX[fmt]}"
        paths[fmt].write_bytes(data)
    return paths

def mini_source(tmp_path, monkeypatch, labels=MINI, **edits):
    """A small synthetic IGL source over ``labels``, registered as validated for this test.

    Returns ``(source, blobs)``.  The lock registered in scope is the one ``from_files`` seals
    (local files: the scope compares digests and sizes only)."""
    from tests.support.fixtures.ogrdb_inventory import seal as seal_lock
    blobs = small(IGL, set(labels), **edits)
    in_scope(monkeypatch, seal_lock(IGL, blobs))
    return from_files(tmp_path, files=written(tmp_path, blobs)), blobs

def full_recipe(tmp_path, source: OgrdbSource, labels=MINI, **kwargs) -> dict:
    """A recipe selecting ``labels`` from the pinned IGL release 4 ``airr`` export."""
    blobs = {fmt: fixture_bytes(IGL, fmt) for fmt in ALL}
    return recipe(source.lock, [entry(blobs, "airr", label) for label in labels], **kwargs)

def error(kind, call, *args, **kwargs) -> ReferenceSourceError:
    with pytest.raises(ReferenceSourceError) as info:
        call(*args, **kwargs)
    assert info.value.kind == kind, f"expected {kind}, got {info.value.kind}: {info.value}"
    return info.value

def tree(directory: pathlib.Path) -> dict:
    """Every file below ``directory`` (relative POSIX name -> bytes)."""
    return {p.relative_to(directory).as_posix(): p.read_bytes() for p in sorted(directory.rglob("*"))
            if p.is_file()}

@pytest.fixture
def offline(monkeypatch):
    """Any attempt to open a connection, or to create the fetch transport, fails the test."""
    def refuse(*args, **kwargs):
        pytest.fail("an offline operation tried to reach the network")

    monkeypatch.setattr(socket, "create_connection", refuse)
    monkeypatch.setattr(socket, "getaddrinfo", refuse)
    monkeypatch.setattr(socket.socket, "connect", refuse)
    monkeypatch.setattr(http.client.HTTPSConnection, "connect", refuse)
    monkeypatch.setattr(http.client.HTTPConnection, "connect", refuse)
    monkeypatch.setattr(HttpsTransport, "get", refuse)
    monkeypatch.setattr(api, "_transport", refuse)
