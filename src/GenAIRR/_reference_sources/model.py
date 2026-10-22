"""Closed data model of published germline-source import (single import point).

Every class is a frozen dataclass whose constructor validates it; ``from_dict`` also
rejects missing and unknown keys and locates every error with a JSON-pointer-like path.
"""
from __future__ import annotations

from GenAIRR._reference_sources._fields import (  # noqa: F401
    MAX_NAME, MAX_POINTER, MAX_REFS, MAX_SEQUENCE, MAX_TEXT, MAX_URL, MAX_VALUE_DEPTH, MAX_VALUE_NODES,
    MAX_VALUE_STRING, MAX_VERSION, JsonValue,
)
from GenAIRR._reference_sources.common import (  # noqa: F401
    ADAPTER_ID, ADAPTER_VERSION, CHAIN_TYPES, COORDINATE_PROFILE_ID, COORDINATE_PROFILE_VERSION,
    DEFAULTS_PROFILES, EVIDENCE_FORMATS, FORMAT_ROLES, FORMATS, GAPPED_VARIANT, LOCUS_CHAIN_TYPES,
    PROJECTION_ID, PROJECTION_VERSION, PROVIDERS, RECORD_EXPORT_FORMATS, ROLES, SEGMENTS,
    SUPPORTED_ADAPTERS, SUPPORTED_PROJECTIONS, AdapterRef, CoordinateProfileRef, Generator, Locator,
    ProjectionRef, RuntimeHashes, SourceIdentity, Taxon,
)
from GenAIRR._reference_sources.lock import (  # noqa: F401
    MAX_ASSET_BYTES, MAX_LOCK_BYTES, MISSING_REASONS, ORIGIN_KINDS, TERMS_STATUSES, Asset, MissingAsset,
    SourceLock, Terms,
)
from GenAIRR._reference_sources.digests import DOMAINS, lock_hash  # noqa: F401
