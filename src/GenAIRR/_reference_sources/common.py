"""Vocabularies, registries and small shared shapes of the published-source model.

Asset identifiers are generated, never copied from source file names: each locked
asset's id is its logical format, so a :class:`Locator` (asset id plus JSON pointer)
names one place in one exact set of locked bytes.
"""
from __future__ import annotations

from dataclasses import dataclass
import re

from GenAIRR._reference.defaults import PROFILE as _UNFITTED_DEFAULTS_PROFILE
from GenAIRR._reference_sources._fields import (
    MAX_NAME, MAX_TAXONOMY_ID, MAX_VERSION, bad, check_canonical_order, check_digest, check_enum,
    check_instance, check_int, check_member_path, check_pointer, check_text, check_token, check_tuple,
    check_unique, construct, pointer_sort_key, ptr, require_object, show,
)

# ── Vocabularies ────────────────────────────────────────────────────

PROVIDERS = ("ogrdb",)
SEGMENTS = ("V", "D", "J")
CHAIN_TYPES = ("vdj", "vj")
LOCUS_CHAIN_TYPES = {"IGH": "vdj", "IGK": "vj", "IGL": "vj", "TRA": "vj", "TRB": "vdj", "TRD": "vdj",
                     "TRG": "vj"}
FORMATS = ("airr", "airr_ex", "gapped", "gapped_ex", "ungapped", "ungapped_ex",
           "licence_statement", "citation_statement")
ROLES = ("record_export", "gapped_alignment", "ungapped_sequences", "data_terms")
FORMAT_ROLES = {"airr": "record_export", "airr_ex": "record_export",
                "gapped": "gapped_alignment", "gapped_ex": "gapped_alignment",
                "ungapped": "ungapped_sequences", "ungapped_ex": "ungapped_sequences",
                "licence_statement": "data_terms", "citation_statement": "data_terms"}
# Formats every lock accounts for, as a present asset or an explicit missing asset.
EVIDENCE_FORMATS = ("airr", "airr_ex", "gapped", "gapped_ex")
RECORD_EXPORT_FORMATS = ("airr", "airr_ex")
# A record export is aligned only against the gapped export of the same variant.
GAPPED_VARIANT = {"airr": "gapped", "airr_ex": "gapped_ex"}

# ── Registries of supported semantic profiles ───────────────────────

ADAPTER_ID = "ogrdb_germline_set"
ADAPTER_VERSION = 1
COORDINATE_PROFILE_ID = "ogrdb_api_v2_coordinates"
COORDINATE_PROFILE_VERSION = 1
PROJECTION_ID = "observed_span"
PROJECTION_VERSION = 1
DEFAULTS_PROFILES = (_UNFITTED_DEFAULTS_PROFILE,)

# (adapter id, version) -> (provider, supported (coordinate profile id, version) pairs)
SUPPORTED_ADAPTERS = {
    (ADAPTER_ID, ADAPTER_VERSION): ("ogrdb", ((COORDINATE_PROFILE_ID, COORDINATE_PROFILE_VERSION),)),
}
SUPPORTED_PROJECTIONS = ((PROJECTION_ID, PROJECTION_VERSION),)

_PROFILE_ID_RE = re.compile(r"[a-z][a-z0-9_.]{0,63}")
_ASSET_ID_RE = re.compile(r"[a-z][a-z0-9_]{0,63}")
_GENERATOR_VERSION_RE = re.compile(r"[0-9A-Za-z][0-9A-Za-z.+_-]{0,63}")
MAX_ASSET_ID = 64


def check_asset_id(value, path: str) -> str:
    """A generated asset id: one of :data:`FORMATS`; path forms are rejected explicitly."""
    check_member_path(value, path, MAX_ASSET_ID)
    if "/" in value:
        raise bad(path, f"asset id {show(value)} contains a path separator")
    if not _ASSET_ID_RE.fullmatch(value):
        raise bad(path, f"asset id {show(value)} must match [a-z][a-z0-9_]*")
    if value not in FORMATS:
        raise bad(path, f"unknown asset id {value!r}; asset ids are generated from formats {list(FORMATS)}")
    return value


# ── Shapes ──────────────────────────────────────────────────────────


@dataclass(frozen=True)
class Locator:
    """One place in one locked asset: ``{"asset": <asset id>, "pointer": <JSON pointer>}``.

    For JSON assets the pointer is an RFC 6901 pointer; for FASTA assets it is ``/<n>``,
    the zero-based record index.  ``id`` is ``"<asset>#<pointer>"``.
    """

    asset: str
    pointer: str

    KEYS = ("asset", "pointer")

    def __post_init__(self):
        check_asset_id(self.asset, "/asset")
        check_pointer(self.pointer, "/pointer")

    @property
    def id(self) -> str:
        return f"{self.asset}#{self.pointer}"

    @property
    def role(self) -> str:
        return FORMAT_ROLES[self.asset]

    def sort_key(self) -> tuple:
        return (FORMATS.index(self.asset), pointer_sort_key(self.pointer))

    def is_under(self, other: "Locator") -> bool:
        """True if this locator is strictly inside ``other`` (same asset, longer pointer)."""
        return self.asset == other.asset and self.pointer.startswith(other.pointer + "/")

    def to_dict(self) -> dict:
        return {"asset": self.asset, "pointer": self.pointer}

    @classmethod
    def from_dict(cls, d, path: str = "") -> "Locator":
        require_object(d, cls.KEYS, path)
        return construct(cls, path, asset=d["asset"], pointer=d["pointer"])


def locator_tuple(value, path: str, max_items: int, *, canonical: bool) -> tuple:
    """Validate a tuple of locators; ``canonical`` requires sorted order without duplicates."""
    check_tuple(value, path, max_items, Locator)
    if canonical:
        check_canonical_order([loc.sort_key() for loc in value], path, "locator")
    else:
        check_unique([loc.id for loc in value], path, "locator")
    return value


@dataclass(frozen=True)
class Taxon:
    scientific_name: str
    taxonomy_id: int

    KEYS = ("scientific_name", "taxonomy_id")

    def __post_init__(self):
        check_text(self.scientific_name, "/scientific_name", max_len=MAX_NAME, multiline=False)
        check_int(self.taxonomy_id, "/taxonomy_id", minimum=1, maximum=MAX_TAXONOMY_ID)

    def to_dict(self) -> dict:
        return {"scientific_name": self.scientific_name, "taxonomy_id": self.taxonomy_id}

    @classmethod
    def from_dict(cls, d, path: str = "") -> "Taxon":
        require_object(d, cls.KEYS, path)
        return construct(cls, path, scientific_name=d["scientific_name"], taxonomy_id=d["taxonomy_id"])


@dataclass(frozen=True)
class SourceIdentity:
    """Validated scientific identity of a source set: taxon, locus and chain type."""

    taxon: Taxon
    locus: str
    chain_type: str

    KEYS = ("taxon", "locus", "chain_type")

    def __post_init__(self):
        check_instance(self.taxon, Taxon, "/taxon")
        check_enum(self.locus, tuple(LOCUS_CHAIN_TYPES), "/locus")
        check_enum(self.chain_type, CHAIN_TYPES, "/chain_type")
        if self.chain_type != LOCUS_CHAIN_TYPES[self.locus]:
            raise bad("/chain_type", f"locus {self.locus} has chain type {LOCUS_CHAIN_TYPES[self.locus]!r}, "
                                     f"got {self.chain_type!r}")

    @property
    def segments(self) -> tuple:
        return SEGMENTS if self.chain_type == "vdj" else ("V", "J")

    def to_dict(self) -> dict:
        return {"taxon": self.taxon.to_dict(), "locus": self.locus, "chain_type": self.chain_type}

    @classmethod
    def from_dict(cls, d, path: str = "") -> "SourceIdentity":
        require_object(d, cls.KEYS, path)
        return construct(cls, path, taxon=Taxon.from_dict(d["taxon"], ptr(path, "taxon")),
                         locus=d["locus"], chain_type=d["chain_type"])


def _versioned(self) -> None:
    check_token(self.id, "/id", _PROFILE_ID_RE, "an identifier matching [a-z][a-z0-9_.]*")
    check_int(self.version, "/version", minimum=1, maximum=MAX_VERSION)


@dataclass(frozen=True)
class CoordinateProfileRef:
    """The coordinate profile: which provider field meanings an adapter has established."""

    id: str
    version: int

    KEYS = ("id", "version")

    def __post_init__(self):
        _versioned(self)

    def to_dict(self) -> dict:
        return {"id": self.id, "version": self.version}

    @classmethod
    def from_dict(cls, d, path: str = "") -> "CoordinateProfileRef":
        require_object(d, cls.KEYS, path)
        return construct(cls, path, id=d["id"], version=d["version"])


@dataclass(frozen=True)
class AdapterRef:
    """``{"id", "version", "coordinate_profile": {"id", "version"}}``; unknown versions are rejected."""

    id: str
    version: int
    coordinate_profile: CoordinateProfileRef

    KEYS = ("id", "version", "coordinate_profile")

    def __post_init__(self):
        _versioned(self)
        check_instance(self.coordinate_profile, CoordinateProfileRef, "/coordinate_profile")
        entry = SUPPORTED_ADAPTERS.get((self.id, self.version))
        if entry is None:
            raise bad("", f"unknown adapter {self.id!r} version {self.version}; supported: "
                          f"{sorted(SUPPORTED_ADAPTERS)}", "unsupported_profile")
        profile = (self.coordinate_profile.id, self.coordinate_profile.version)
        if profile not in entry[1]:
            raise bad("", f"adapter {self.id!r} version {self.version} does not support coordinate profile "
                          f"{profile[0]!r} version {profile[1]}; supported: {list(entry[1])}",
                      "unsupported_profile")

    @property
    def provider(self) -> str:
        return SUPPORTED_ADAPTERS[(self.id, self.version)][0]

    def to_dict(self) -> dict:
        return {"id": self.id, "version": self.version, "coordinate_profile": self.coordinate_profile.to_dict()}

    @classmethod
    def from_dict(cls, d, path: str = "") -> "AdapterRef":
        require_object(d, cls.KEYS, path)
        profile = CoordinateProfileRef.from_dict(d["coordinate_profile"], ptr(path, "coordinate_profile"))
        return construct(cls, path, id=d["id"], version=d["version"], coordinate_profile=profile)


@dataclass(frozen=True)
class ProjectionRef:
    """``{"id", "version"}`` of the projection into a runtime reference; unknown versions are rejected."""

    id: str
    version: int

    KEYS = ("id", "version")

    def __post_init__(self):
        _versioned(self)
        if (self.id, self.version) not in SUPPORTED_PROJECTIONS:
            raise bad("", f"unknown projection {self.id!r} version {self.version}; supported: "
                          f"{list(SUPPORTED_PROJECTIONS)}", "unsupported_profile")

    def to_dict(self) -> dict:
        return {"id": self.id, "version": self.version}

    @classmethod
    def from_dict(cls, d, path: str = "") -> "ProjectionRef":
        require_object(d, cls.KEYS, path)
        return construct(cls, path, id=d["id"], version=d["version"])


def check_defaults_profile(value, path: str) -> str:
    if type(value) is not str:
        raise bad(path, f"expected a defaults profile id, got {show(value)}")
    if value not in DEFAULTS_PROFILES:
        raise bad(path, f"unknown defaults profile {show(value)}; supported: {list(DEFAULTS_PROFILES)}",
                  "unsupported_profile")
    return value


@dataclass(frozen=True)
class Generator:
    """The GenAIRR build that actually produced a document: ``{"name": "GenAIRR", "version"}``."""

    name: str
    version: str

    KEYS = ("name", "version")

    def __post_init__(self):
        if self.name != "GenAIRR":
            raise bad("/name", f"expected 'GenAIRR', got {show(self.name)}")
        check_token(self.version, "/version", _GENERATOR_VERSION_RE, "a version string")

    def to_dict(self) -> dict:
        return {"name": self.name, "version": self.version}

    @classmethod
    def from_dict(cls, d, path: str = "") -> "Generator":
        require_object(d, cls.KEYS, path)
        return construct(cls, path, name=d["name"], version=d["version"])


@dataclass(frozen=True)
class RuntimeHashes:
    """The four identities of a runtime ``Reference``, with their existing meanings."""

    catalogue_hash: str
    models_hash: str
    execution_hash: str
    artifact_hash: str

    KEYS = ("catalogue_hash", "models_hash", "execution_hash", "artifact_hash")

    def __post_init__(self):
        for key in self.KEYS:
            check_digest(getattr(self, key), "/" + key)

    def to_dict(self) -> dict:
        return {k: getattr(self, k) for k in self.KEYS}

    @classmethod
    def from_dict(cls, d, path: str = "") -> "RuntimeHashes":
        require_object(d, cls.KEYS, path)
        return construct(cls, path, **{k: d[k] for k in cls.KEYS})
