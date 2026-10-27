"""Dossier manifest ``genairr_reference_source_build_manifest`` version 1 (closed shape).

The manifest lists every deterministic dossier member by relative name, digest and size,
the member schema versions, the adapter and projection, the input and report hashes,
the expected runtime hashes and the actual builder.  Member names are a closed set of
fixed names plus content-addressed blob paths ``cache/blobs/sha256/<hex>``; absolute,
drive, backslash, empty, ``.``/``..`` and over-long forms are rejected, and acquisition
receipts (``receipts/...``) are never members.  The manifest hash is computed externally.
"""
from __future__ import annotations

from dataclasses import dataclass
import re

from GenAIRR._reference.artifact import SCHEMA_VERSION as REFERENCE_SCHEMA_VERSION
from GenAIRR._reference.sidecar import SCHEMA_VERSION as ANNOTATIONS_SCHEMA_VERSION
from GenAIRR._reference_sources._fields import (
    bad, check_canonical_order, check_digest, check_instance, check_int, check_member_path, check_tuple,
    construct, parse_array, parse_json, ptr, read_json, require_object, check_header,
)
from GenAIRR._reference_sources.common import AdapterRef, Generator, ProjectionRef, RuntimeHashes
from GenAIRR._reference_sources.digests import report_hash
from GenAIRR._reference_sources.lock import SCHEMA_VERSION as LOCK_SCHEMA_VERSION, SourceLock
from GenAIRR._reference_sources.recipe import SCHEMA_VERSION as RECIPE_SCHEMA_VERSION
from GenAIRR._reference_sources.report import SCHEMA_VERSION as REPORT_SCHEMA_VERSION, BuildReport

ARTIFACT = "genairr_reference_source_build_manifest"
SCHEMA_VERSION = 1
MANIFEST_NAME = "manifest.json"
MAX_MANIFEST_BYTES = 4 * 1024 * 1024
MAX_MEMBERS = 4096
MAX_MEMBER_NAME = 255
MAX_MEMBER_BYTES = 2 ** 40
FIXED_MEMBERS = ("inputs/annotations.json", "inputs/d.fasta", "inputs/j.fasta", "inputs/v.fasta", "recipe.json",
                 "reference.genairr-ref.json", "report.json", "source.lock.json")
OPTIONAL_MEMBERS = ("inputs/d.fasta",)          # present exactly for VDJ references
BLOB_PREFIX = "cache/blobs/sha256/"
RECEIPTS_PREFIX = "receipts/"
MEMBER_SCHEMA_VERSIONS = {"source_lock": LOCK_SCHEMA_VERSION, "recipe": RECIPE_SCHEMA_VERSION,
                          "report": REPORT_SCHEMA_VERSION, "reference": REFERENCE_SCHEMA_VERSION,
                          "annotations": ANNOTATIONS_SCHEMA_VERSION}

_BLOB_RE = re.compile(r"cache/blobs/sha256/([0-9a-f]{64})")
_KEYS = ("artifact", "schema_version", "builder", "versions", "adapter", "projection", "hashes", "runtime",
         "members")


def check_member_name(name, path: str) -> str:
    """A dossier member name: relative, closed, never a receipt or the manifest itself."""
    check_member_path(name, path, MAX_MEMBER_NAME)
    if name.startswith(RECEIPTS_PREFIX):
        raise bad(path, f"{name!r}: acquisition receipts are audit events outside the manifest")
    if name == MANIFEST_NAME:
        raise bad(path, "the manifest does not list itself")
    if name not in FIXED_MEMBERS and not _BLOB_RE.fullmatch(name):
        raise bad(path, f"{name!r} is not a dossier member; expected one of {list(FIXED_MEMBERS)} or "
                        f"{BLOB_PREFIX}<64 lowercase hex>")
    return name


def _closed(cls, d, path, **nested):
    require_object(d, cls.KEYS, path)
    return construct(cls, path, **{k: (nested[k](d[k], ptr(path, k)) if k in nested else d[k]) for k in cls.KEYS})


@dataclass(frozen=True)
class Member:
    name: str
    sha256: str
    bytes: int

    KEYS = ("name", "sha256", "bytes")

    def __post_init__(self):
        check_member_name(self.name, "/name")
        check_digest(self.sha256, "/sha256")
        check_int(self.bytes, "/bytes", minimum=0, maximum=MAX_MEMBER_BYTES, limit_kind="resource_limit")
        blob = _BLOB_RE.fullmatch(self.name)
        if blob and self.sha256 != "sha256:" + blob.group(1):
            raise bad("/sha256", f"blob {self.name} must have digest sha256:{blob.group(1)}, got {self.sha256}",
                      "digest_mismatch")

    @property
    def is_blob(self) -> bool:
        return self.name.startswith(BLOB_PREFIX)

    def to_dict(self) -> dict:
        return {"name": self.name, "sha256": self.sha256, "bytes": self.bytes}

    @classmethod
    def from_dict(cls, d, path: str = "") -> "Member":
        return _closed(cls, d, path)


@dataclass(frozen=True)
class MemberVersions:
    """Schema versions of the lock, recipe, report, reference artifact and annotation sidecar."""

    source_lock: int
    recipe: int
    report: int
    reference: int
    annotations: int

    KEYS = ("source_lock", "recipe", "report", "reference", "annotations")

    def __post_init__(self):
        for key in self.KEYS:
            found = check_int(getattr(self, key), "/" + key)
            if found != MEMBER_SCHEMA_VERSIONS[key]:
                raise bad("/" + key, f"unsupported {key} schema version {found}; this GenAIRR reads version "
                                     f"{MEMBER_SCHEMA_VERSIONS[key]}", "unsupported_schema")

    def to_dict(self) -> dict:
        return {k: getattr(self, k) for k in self.KEYS}

    @classmethod
    def from_dict(cls, d, path: str = "") -> "MemberVersions":
        return _closed(cls, d, path)


@dataclass(frozen=True)
class ManifestHashes:
    """Input hashes (as in the report) plus the report hash."""

    source_lock: str
    inventory: str
    recipe: str
    build_inputs: str
    report: str

    KEYS = ("source_lock", "inventory", "recipe", "build_inputs", "report")

    def __post_init__(self):
        for key in self.KEYS:
            check_digest(getattr(self, key), "/" + key)

    def to_dict(self) -> dict:
        return {k: getattr(self, k) for k in self.KEYS}

    @classmethod
    def from_dict(cls, d, path: str = "") -> "ManifestHashes":
        return _closed(cls, d, path)


@dataclass(frozen=True)
class BuildManifest:
    builder: Generator
    versions: MemberVersions
    adapter: AdapterRef
    projection: ProjectionRef
    hashes: ManifestHashes
    runtime: RuntimeHashes
    members: tuple

    def __post_init__(self):
        for name, cls in (("builder", Generator), ("versions", MemberVersions), ("adapter", AdapterRef),
                          ("projection", ProjectionRef), ("hashes", ManifestHashes), ("runtime", RuntimeHashes)):
            check_instance(getattr(self, name), cls, "/" + name)
        check_tuple(self.members, "/members", MAX_MEMBERS, Member)
        check_canonical_order([x.name for x in self.members], "/members", "member name")
        names = {x.name for x in self.members}
        missing = [n for n in FIXED_MEMBERS if n not in OPTIONAL_MEMBERS and n not in names]
        if missing:
            raise bad("/members", f"missing required members {missing}")

    def member(self, name: str) -> Member:
        for member in self.members:
            if member.name == name:
                return member
        raise KeyError(name)

    def to_dict(self) -> dict:
        return {"artifact": ARTIFACT, "schema_version": SCHEMA_VERSION, "builder": self.builder.to_dict(),
                "versions": self.versions.to_dict(), "adapter": self.adapter.to_dict(),
                "projection": self.projection.to_dict(), "hashes": self.hashes.to_dict(),
                "runtime": self.runtime.to_dict(), "members": [x.to_dict() for x in self.members]}

    @classmethod
    def from_dict(cls, d, path: str = "") -> "BuildManifest":
        check_header(d, ARTIFACT, SCHEMA_VERSION, path)
        require_object(d, _KEYS, path)
        return construct(
            cls, path, builder=Generator.from_dict(d["builder"], ptr(path, "builder")),
            versions=MemberVersions.from_dict(d["versions"], ptr(path, "versions")),
            adapter=AdapterRef.from_dict(d["adapter"], ptr(path, "adapter")),
            projection=ProjectionRef.from_dict(d["projection"], ptr(path, "projection")),
            hashes=ManifestHashes.from_dict(d["hashes"], ptr(path, "hashes")),
            runtime=RuntimeHashes.from_dict(d["runtime"], ptr(path, "runtime")),
            members=parse_array(d["members"], ptr(path, "members"), MAX_MEMBERS, Member.from_dict))

    @classmethod
    def from_json(cls, text: str) -> "BuildManifest":
        return cls.from_dict(parse_json(text, max_bytes=MAX_MANIFEST_BYTES))

    @classmethod
    def read(cls, path) -> "BuildManifest":
        return cls.from_dict(read_json(path, max_bytes=MAX_MANIFEST_BYTES))


def check_manifest(manifest: BuildManifest, lock: SourceLock, report: BuildReport = None) -> None:
    """The manifest commits to this lock (every locked blob, D input exactly for VDJ) and,
    when given, to this report: its hashes, report hash, runtime hashes, builder and projection.

    ``missing_asset``: a locked asset has no blob member.  ``digest_mismatch``: a hash or a
    blob size differs.  ``identity_mismatch``: adapter, builder or projection differ, a blob
    is not locked, or the D input does not match the chain type.
    """
    if not isinstance(manifest, BuildManifest):
        raise TypeError(f"manifest must be a BuildManifest, got {type(manifest).__name__}")
    if not isinstance(lock, SourceLock):
        raise TypeError(f"lock must be a SourceLock, got {type(lock).__name__}")
    if manifest.hashes.source_lock != lock.commitment:
        raise bad("/hashes/source_lock", f"manifest commits to lock {manifest.hashes.source_lock}, not "
                                         f"{lock.commitment}", "digest_mismatch")
    if manifest.adapter != lock.adapter:
        raise bad("/adapter", "manifest adapter differs from the lock adapter", "identity_mismatch")
    blobs = {x.name: (i, x) for i, x in enumerate(manifest.members) if x.is_blob}
    for asset in lock.assets:
        name = BLOB_PREFIX + asset.blob_name
        if name not in blobs:
            raise bad("/members", f"locked asset {asset.asset_id} ({asset.sha256}) has no blob member {name}",
                      "missing_asset")
        index, member = blobs.pop(name)
        if member.bytes != asset.bytes:
            raise bad(ptr(ptr("/members", index), "bytes"), f"blob of {asset.asset_id} has {member.bytes} bytes; "
                                                            f"the lock says {asset.bytes}", "digest_mismatch")
    for name, (index, _) in sorted(blobs.items()):
        raise bad(ptr(ptr("/members", index), "name"), f"blob {name} is not a locked asset", "identity_mismatch")
    has_d = any(x.name == "inputs/d.fasta" for x in manifest.members)
    if has_d != (lock.identity.chain_type == "vdj"):
        raise bad("/members", f"inputs/d.fasta is {'present' if has_d else 'absent'} for a "
                              f"{lock.identity.chain_type} source", "identity_mismatch")
    if report is None:
        return
    if not isinstance(report, BuildReport):
        raise TypeError(f"report must be a BuildReport, got {type(report).__name__}")
    if report.source_lock != lock:
        raise bad("/hashes/source_lock", "the report was built from a different lock", "identity_mismatch")
    for key in ("source_lock", "inventory", "recipe", "build_inputs"):
        if getattr(manifest.hashes, key) != getattr(report.hashes, key):
            raise bad(ptr("/hashes", key), f"manifest {getattr(manifest.hashes, key)} differs from the report "
                                           f"{getattr(report.hashes, key)}", "digest_mismatch")
    actual = report_hash(report)
    if manifest.hashes.report != actual:
        raise bad("/hashes/report", f"manifest commits to report {manifest.hashes.report}, but the report "
                                    f"hashes to {actual}", "digest_mismatch")
    if manifest.runtime != report.runtime:
        raise bad("/runtime", f"manifest runtime hashes {manifest.runtime.to_dict()} differ from the report "
                              f"{report.runtime.to_dict()}", "digest_mismatch")
    if manifest.builder != report.generator:
        raise bad("/builder", f"manifest builder {manifest.builder.to_dict()} differs from the report generator "
                              f"{report.generator.to_dict()}", "identity_mismatch")
    if manifest.projection != report.recipe.projection:
        raise bad("/projection", "manifest projection differs from the recipe projection", "identity_mismatch")
