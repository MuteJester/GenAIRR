"""Closed data model of published germline-source import (single import point).

Every class is a frozen dataclass whose constructor validates it; ``from_dict`` also
rejects missing and unknown keys and locates every error with a JSON-pointer-like path.
"""
from __future__ import annotations

from GenAIRR._reference_sources._fields import (  # noqa: F401
    MAX_ASSET_BYTES, MAX_ASSETS, MAX_JSON_DEPTH, MAX_NAME, MAX_POINTER, MAX_REFS, MAX_SEQUENCE,
    MAX_SOURCE_PAYLOAD_BYTES, MAX_SOURCE_RECORDS, MAX_TEXT, MAX_URL, MAX_VALUE_NODES, MAX_VALUE_STRING,
    MAX_VERSION, JsonValue,
)
from GenAIRR._reference_sources.common import (  # noqa: F401
    ADAPTER_ID, ADAPTER_VERSION, CHAIN_TYPES, COORDINATE_PROFILE_ID, COORDINATE_PROFILE_VERSION,
    DEFAULTS_PROFILES, EVIDENCE_FORMATS, FORMAT_ROLES, FORMATS, GAPPED_VARIANT, LOCUS_CHAIN_TYPES,
    PROJECTION_ID, PROJECTION_VERSION, PROVIDERS, RECORD_EXPORT_FORMATS, ROLES, SEGMENTS,
    PROFILE_SCOPES, SUPPORTED_ADAPTERS, SUPPORTED_PROJECTIONS, CLAIM_ID_RE, FINDING_ID_RE, MAPPING_ID_RE, AdapterRef, CoordinateProfileRef, Generator, Locator,
    ProjectionRef, RuntimeHashes, SourceIdentity, Taxon,
)
from GenAIRR._reference_sources.lock import (  # noqa: F401
    MAX_LOCK_BYTES, MISSING_REASONS, ORIGIN_KINDS, TERMS_STATUSES, Asset, MissingAsset,
    SourceLock, Terms,
)
from GenAIRR._reference_sources.receipt import (  # noqa: F401
    MAX_RECEIPT_BYTES, RECEIPT_HEADERS, RECEIPT_ORIGINS, AcquisitionReceipt,
)
from GenAIRR._reference_sources.recipe import (  # noqa: F401
    DECISION_ACTIONS, EVIDENCE_LEVELS, MAX_RECIPE_BYTES, Decision, Recipe, RecipeReference, Requirement,
    SelectionEntry,
)
from GenAIRR._reference_sources.inventory import (  # noqa: F401
    BOUNDARIES, COORDINATE_PREDICATES, CLAIM_ID_DOMAIN, CLAIM_STATUSES, CLAIM_SUPPORTS, COORDINATE_SYSTEMS, DERIVATIONS,
    EXTENSION_SIDES, FASTA_VALUE_KEYS, IMGT_GAP, FINDING_CATEGORIES, FINDING_ID_DOMAIN, FUNCTIONALITIES, INTERPRETATIONS, MAPPING_ID_DOMAIN,
    MAPPING_KINDS, MAPPING_SPACES, MAPPING_TRANSFORMS, MAPPING_VALIDATIONS, MAX_INVENTORY_BYTES, PHASE_INDEPENDENCE,
    PREDICATE_DERIVATIONS, PREDICATE_SEGMENTS, PREDICATE_SUBJECTS, PREDICATE_SYSTEMS, PREDICATES, PRESENCES,
    RECORD_STATUSES, RELATIONSHIPS, RUNTIME_PREDICATES, SEVERITIES, SUPPORT_KINDS, SUPPORT_STATUSES, Block, Claim,
    Finding, Mapping, Mention, MentionedIdentity, Observation, SequenceCandidate, SourceInventory, SourceRecord,
    SupportReference, alignment_blocks, check_claim_value, text_digest,
)
from GenAIRR._reference_sources.consistency import (  # noqa: F401
    check_build_inputs, check_decisions, check_inventory, check_profile_scope, check_recipe, check_selection,
    profile_scope_problem,
)
from GenAIRR._reference_sources.consistency import selected_subjects  # noqa: F401
from GenAIRR._reference_sources.digests import (  # noqa: F401
    BUILD_INPUTS_DOMAIN, DOMAINS, INVENTORY_DOMAIN, LOCK_DOMAIN, MANIFEST_DOMAIN, RECIPE_DOMAIN, REPORT_DOMAIN,
    build_input_hash, build_inputs, inventory_hash, lock_hash, manifest_hash, recipe_hash, report_hash,
)
from GenAIRR._reference_sources.report import (  # noqa: F401
    CLAIM_OUTCOMES, EXCLUSION_REASONS, MAX_REPORT_BYTES, PHASE_SUPPORTS, RECORD_OUTCOMES, RELATED_EXCLUSIONS,
    BuildHashes, BuildReport, CapabilityStatus, ClaimOutcome, Comparison, Counts, ExportCounts, OutcomeCounts,
    RecordOutcome, RequirementOutcome, RuntimeCapabilities, SegmentCounts, SelectedAllele, UnmetRequirement,
)
from GenAIRR._reference_sources.manifest import (  # noqa: F401
    BLOB_PREFIX, FIXED_MEMBERS, MANIFEST_NAME, MAX_MANIFEST_BYTES, MAX_MEMBER_NAME, MAX_MEMBERS,
    MEMBER_SCHEMA_VERSIONS, OPTIONAL_MEMBERS, RECEIPTS_PREFIX, BuildManifest, ManifestHashes, Member,
    MemberVersions, check_manifest, check_member_name,
)
from GenAIRR._reference.model import sequence_digest  # noqa: F401
