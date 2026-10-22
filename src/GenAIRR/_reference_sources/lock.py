"""Source lock ``genairr_reference_source_lock`` version 1.

A lock pins the exact bytes of a published germline set: which provider set and
release, its validated scientific identity, the adapter and coordinate profile that
interpret it, every asset by raw SHA-256 digest and size, the evidence assets that are
explicitly absent, and the data terms.  Its ``commitment`` is the typed digest of every
other field.  Retrieval times, absolute paths, local file names, HTTP headers and
acquisition receipts are never part of a lock.
"""
from __future__ import annotations

from dataclasses import dataclass
import re
from typing import Optional

from GenAIRR._reference_sources._fields import (
    MAX_TEXT, MAX_RELEASE, bad, check_canonical_order, check_digest, check_enum, check_instance, check_int,
    check_text, check_token, check_tuple, check_url, construct, digest, parse_array, parse_json, ptr,
    read_json, require_object, check_header,
)
from GenAIRR._reference_sources.common import (
    EVIDENCE_FORMATS, FORMAT_ROLES, FORMATS, PROVIDERS, ROLES, AdapterRef, SourceIdentity, check_asset_id,
)

ARTIFACT = "genairr_reference_source_lock"
SCHEMA_VERSION = 1
DIGEST_DOMAIN = "genairr_reference_source_lock_v1"
MAX_LOCK_BYTES = 1 << 20
MAX_ASSET_BYTES = 64 * 1024 * 1024
MAX_ASSETS = len(FORMATS)

ORIGIN_KINDS = ("provider_download", "local_file")
MISSING_REASONS = ("not_requested", "not_supplied", "not_offered")
TERMS_STATUSES = ("recorded", "unresolved")

_SET_ID_RE = re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]{0,127}")
_LOCK_KEYS = ("artifact", "schema_version", "provider", "set_id", "release", "identity", "adapter",
              "assets", "missing_assets", "terms", "commitment")


@dataclass(frozen=True)
class Asset:
    """One locked asset: raw digest and size of exact bytes, logical format, role and origin.

    ``asset_id`` is generated (it equals ``format``); ``declared_url`` is required for a
    provider download and optional for a local file, whose bytes GenAIRR did not fetch.
    """

    asset_id: str
    role: str
    format: str
    sha256: str
    bytes: int
    declared_url: Optional[str]
    origin_kind: str

    KEYS = ("asset_id", "role", "format", "sha256", "bytes", "declared_url", "origin_kind")

    def __post_init__(self):
        check_asset_id(self.asset_id, "/asset_id")
        check_enum(self.format, FORMATS, "/format")
        if self.asset_id != self.format:
            raise bad("/asset_id", f"asset id must be the generated id {self.format!r}, got {self.asset_id!r}")
        check_enum(self.role, ROLES, "/role")
        if self.role != FORMAT_ROLES[self.format]:
            raise bad("/role", f"format {self.format!r} has role {FORMAT_ROLES[self.format]!r}, got {self.role!r}")
        check_digest(self.sha256, "/sha256")
        check_int(self.bytes, "/bytes", minimum=1, maximum=MAX_ASSET_BYTES, limit_kind="resource_limit")
        check_enum(self.origin_kind, ORIGIN_KINDS, "/origin_kind")
        check_url(self.declared_url, "/declared_url", nullable=self.origin_kind == "local_file")

    @property
    def blob_name(self) -> str:
        """Content-addressed cache name: the hex digest of the raw bytes."""
        return self.sha256.split(":", 1)[1]

    def to_dict(self) -> dict:
        return {"asset_id": self.asset_id, "role": self.role, "format": self.format, "sha256": self.sha256,
                "bytes": self.bytes, "declared_url": self.declared_url, "origin_kind": self.origin_kind}

    @classmethod
    def from_dict(cls, d, path: str = "") -> "Asset":
        require_object(d, cls.KEYS, path)
        return construct(cls, path, **{k: d[k] for k in cls.KEYS})


@dataclass(frozen=True)
class MissingAsset:
    """An evidence format that is explicitly absent, with the reason; never an empty substitute."""

    format: str
    role: str
    reason: str
    detail: Optional[str]

    KEYS = ("format", "role", "reason", "detail")

    def __post_init__(self):
        check_enum(self.format, FORMATS, "/format")
        check_enum(self.role, ROLES, "/role")
        if self.role != FORMAT_ROLES[self.format]:
            raise bad("/role", f"format {self.format!r} has role {FORMAT_ROLES[self.format]!r}, got {self.role!r}")
        check_enum(self.reason, MISSING_REASONS, "/reason")
        check_text(self.detail, "/detail", nullable=True)

    def to_dict(self) -> dict:
        return {"format": self.format, "role": self.role, "reason": self.reason, "detail": self.detail}

    @classmethod
    def from_dict(cls, d, path: str = "") -> "MissingAsset":
        require_object(d, cls.KEYS, path)
        return construct(cls, path, **{k: d[k] for k in cls.KEYS})


@dataclass(frozen=True)
class Terms:
    """Data terms and citation, separate from the provider's API software licence.

    ``recorded``: ``assets`` (the ``data_terms`` assets holding the statements, in lock
    order) is non-empty, ``licence`` is stated and ``reason`` is null.  ``unresolved``:
    no assets, no licence, and a ``reason``.  ``citation`` may be stated either way.
    """

    status: str
    assets: tuple
    licence: Optional[str]
    citation: Optional[str]
    reason: Optional[str]

    KEYS = ("status", "assets", "licence", "citation", "reason")

    def __post_init__(self):
        check_enum(self.status, TERMS_STATUSES, "/status")
        check_tuple(self.assets, "/assets", MAX_ASSETS, str)
        for i, asset_id in enumerate(self.assets):
            check_asset_id(asset_id, ptr("/assets", i))
        check_canonical_order([FORMATS.index(a) for a in self.assets], "/assets", "terms asset")
        check_text(self.licence, "/licence", nullable=True, max_len=MAX_TEXT)
        check_text(self.citation, "/citation", nullable=True, max_len=MAX_TEXT)
        check_text(self.reason, "/reason", nullable=True, max_len=MAX_TEXT)
        if self.status == "recorded":
            if not self.assets:
                raise bad("/assets", "recorded terms need the data-terms assets that state them")
            if self.licence is None:
                raise bad("/licence", "recorded terms need the stated data licence")
            if self.reason is not None:
                raise bad("/reason", "recorded terms have no unresolved reason")
        else:
            if self.assets:
                raise bad("/assets", "unresolved terms reference no assets")
            if self.licence is not None:
                raise bad("/licence", "unresolved terms state no licence")
            if self.reason is None:
                raise bad("/reason", "unresolved terms need a reason")

    def to_dict(self) -> dict:
        return {"status": self.status, "assets": list(self.assets), "licence": self.licence,
                "citation": self.citation, "reason": self.reason}

    @classmethod
    def from_dict(cls, d, path: str = "") -> "Terms":
        require_object(d, cls.KEYS, path)
        assets = parse_array(d["assets"], ptr(path, "assets"), MAX_ASSETS, lambda v, p: v)
        return construct(cls, path, status=d["status"], assets=assets, licence=d["licence"],
                         citation=d["citation"], reason=d["reason"])


_UNSEALED = "sha256:" + "0" * 64


def _body(provider, set_id, release, identity, adapter, assets, missing_assets, terms) -> dict:
    return {"artifact": ARTIFACT, "schema_version": SCHEMA_VERSION, "provider": provider, "set_id": set_id,
            "release": release, "identity": identity.to_dict(), "adapter": adapter.to_dict(),
            "assets": [a.to_dict() for a in assets], "missing_assets": [a.to_dict() for a in missing_assets],
            "terms": terms.to_dict()}


@dataclass(frozen=True)
class SourceLock:
    """An immutable, validated source lock; construct with :meth:`seal` or :meth:`from_dict`."""

    provider: str
    set_id: str
    release: int
    identity: SourceIdentity
    adapter: AdapterRef
    assets: tuple
    missing_assets: tuple
    terms: Terms
    commitment: str

    def __post_init__(self):
        if type(self.provider) is not str:
            raise bad("/provider", "expected a provider id string")
        if self.provider not in PROVIDERS:
            raise bad("/provider", f"unknown provider {self.provider!r}; supported: {list(PROVIDERS)}",
                      "unsupported_profile")
        check_token(self.set_id, "/set_id", _SET_ID_RE, "a provider set id matching [A-Za-z0-9][A-Za-z0-9._-]*")
        check_int(self.release, "/release", minimum=1, maximum=MAX_RELEASE)
        check_instance(self.identity, SourceIdentity, "/identity")
        check_instance(self.adapter, AdapterRef, "/adapter")
        if self.adapter.provider != self.provider:
            raise bad("/adapter", f"adapter {self.adapter.id!r} reads provider {self.adapter.provider!r}, "
                                  f"not {self.provider!r}", "unsupported_profile")
        check_tuple(self.assets, "/assets", MAX_ASSETS, Asset)
        check_tuple(self.missing_assets, "/missing_assets", MAX_ASSETS, MissingAsset)
        check_instance(self.terms, Terms, "/terms")
        check_canonical_order([FORMATS.index(a.format) for a in self.assets], "/assets", "asset format")
        check_canonical_order([FORMATS.index(a.format) for a in self.missing_assets], "/missing_assets",
                              "missing asset format")
        present = {a.format for a in self.assets}
        if not any(a.role == "record_export" for a in self.assets):
            raise bad("/assets", "a lock needs at least one record export (airr or airr_ex)")
        for i, missing in enumerate(self.missing_assets):
            if missing.format in present:
                raise bad("/missing_assets", f"format {missing.format!r} at index {i} is both present and missing")
        accounted = present | {a.format for a in self.missing_assets}
        unaccounted = [f for f in EVIDENCE_FORMATS if f not in accounted]
        if unaccounted:
            raise bad("/missing_assets", f"evidence formats {unaccounted} are neither locked nor recorded "
                                         "as missing")
        terms_assets = [a.asset_id for a in self.assets if a.role == "data_terms"]
        if list(self.terms.assets) != terms_assets:
            raise bad("/terms/assets", f"terms must reference exactly the data-terms assets {terms_assets}, "
                                       f"got {list(self.terms.assets)}")
        check_digest(self.commitment, "/commitment")
        expected = digest(self.body(), DIGEST_DOMAIN)
        if self.commitment != expected:
            raise bad("/commitment", f"commitment {self.commitment} does not match the lock body "
                                     f"(expected {expected})", "digest_mismatch")

    # -- construction -------------------------------------------------

    @classmethod
    def seal(cls, *, provider, set_id, release, identity, adapter, assets, missing_assets, terms) -> "SourceLock":
        """Build a lock, computing its commitment over the other fields."""
        fields = dict(provider=provider, set_id=set_id, release=release, identity=identity, adapter=adapter,
                      assets=assets, missing_assets=missing_assets, terms=terms)
        try:
            commitment = digest(_body(**fields), DIGEST_DOMAIN)
        except (AttributeError, TypeError):
            commitment = _UNSEALED  # a field has the wrong type; construction names it
        return cls(commitment=commitment, **fields)

    # -- views --------------------------------------------------------

    def body(self) -> dict:
        """Every field the commitment covers (all but ``commitment``)."""
        return _body(self.provider, self.set_id, self.release, self.identity, self.adapter, self.assets,
                     self.missing_assets, self.terms)

    def to_dict(self) -> dict:
        d = self.body()
        d["commitment"] = self.commitment
        return d

    def asset(self, asset_id: str) -> Asset:
        for asset in self.assets:
            if asset.asset_id == asset_id:
                return asset
        raise KeyError(asset_id)

    @property
    def asset_ids(self) -> tuple:
        return tuple(a.asset_id for a in self.assets)

    # -- readers ------------------------------------------------------

    @classmethod
    def from_dict(cls, d, path: str = "") -> "SourceLock":
        check_header(d, ARTIFACT, SCHEMA_VERSION, path)
        require_object(d, _LOCK_KEYS, path)
        return construct(
            cls, path, provider=d["provider"], set_id=d["set_id"], release=d["release"],
            identity=SourceIdentity.from_dict(d["identity"], ptr(path, "identity")),
            adapter=AdapterRef.from_dict(d["adapter"], ptr(path, "adapter")),
            assets=parse_array(d["assets"], ptr(path, "assets"), MAX_ASSETS, Asset.from_dict),
            missing_assets=parse_array(d["missing_assets"], ptr(path, "missing_assets"), MAX_ASSETS,
                                       MissingAsset.from_dict),
            terms=Terms.from_dict(d["terms"], ptr(path, "terms")),
            commitment=d["commitment"])

    @classmethod
    def from_json(cls, text: str) -> "SourceLock":
        return cls.from_dict(parse_json(text, max_bytes=MAX_LOCK_BYTES))

    @classmethod
    def read(cls, path) -> "SourceLock":
        return cls.from_dict(read_json(path, max_bytes=MAX_LOCK_BYTES))
