"""Build recipe ``genairr_reference_recipe`` version 1.

A recipe names the locked source it builds from (``source_lock_hash``), the adapter and
projection versions, the derived reference's authored name, version and build date (its
taxon, locus and chain type must equal the lock identity), an ordered selection of exact
source records, a ledger of decisions, required capabilities and the pinned defaults
profile.  Every list keeps its authored order; duplicates are rejected.  The recipe
digest is computed externally over the complete recipe and is never stored inside it.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

from GenAIRR._reference.capabilities import CAPABILITIES
from GenAIRR._reference_sources._fields import (
    MAX_NAME, MAX_REFS, MAX_TEXT, bad, check_date, check_digest, check_enum, check_instance, check_label,
    check_text, check_token, check_tuple, check_unique, construct, parse_array, parse_json, ptr, read_json,
    require_object, check_header,
)
from GenAIRR._reference_sources.common import (
    CHAIN_TYPES, CLAIM_ID_RE, FINDING_ID_RE, LOCUS_CHAIN_TYPES, SEGMENTS, AdapterRef, Locator, ProjectionRef,
    SourceIdentity, Taxon, check_defaults_profile, locator_tuple,
)

ARTIFACT = "genairr_reference_recipe"
SCHEMA_VERSION = 1
MAX_RECIPE_BYTES = 16 * 1024 * 1024
MAX_SELECTION = 65536
MAX_DECISIONS = 65536
MAX_REQUIRES = len(CAPABILITIES)
EVIDENCE_LEVELS = ("assertion", "source")
DECISION_ACTIONS = ("select_alternative", "withhold_claim", "accept_claim")

_KEYS = ("artifact", "schema_version", "source_lock_hash", "adapter", "projection", "reference", "selection",
         "decisions", "requires", "defaults_profile")


@dataclass(frozen=True)
class RecipeReference:
    """The derived reference: authored ``name``, ``version`` and ``build_date``; the taxon,
    locus and chain type restate the lock identity.  Same shape as the annotation
    sidecar's ``reference`` section."""

    name: str
    version: Optional[str]
    taxon: Taxon
    locus: str
    chain_type: str
    build_date: Optional[str]

    KEYS = ("name", "version", "taxon", "locus", "chain_type", "build_date")

    def __post_init__(self):
        check_text(self.name, "/name", max_len=MAX_NAME, multiline=False)
        check_text(self.version, "/version", nullable=True, max_len=MAX_NAME, multiline=False)
        check_instance(self.taxon, Taxon, "/taxon")
        check_enum(self.locus, tuple(LOCUS_CHAIN_TYPES), "/locus")
        check_enum(self.chain_type, CHAIN_TYPES, "/chain_type")
        if self.chain_type != LOCUS_CHAIN_TYPES[self.locus]:
            raise bad("/chain_type", f"locus {self.locus} has chain type {LOCUS_CHAIN_TYPES[self.locus]!r}, "
                                     f"got {self.chain_type!r}")
        check_date(self.build_date, "/build_date", nullable=True)

    @property
    def identity(self) -> SourceIdentity:
        return SourceIdentity(self.taxon, self.locus, self.chain_type)

    def to_dict(self) -> dict:
        return {"name": self.name, "version": self.version, "taxon": self.taxon.to_dict(), "locus": self.locus,
                "chain_type": self.chain_type, "build_date": self.build_date}

    @classmethod
    def from_dict(cls, d, path: str = "") -> "RecipeReference":
        require_object(d, cls.KEYS, path)
        return construct(cls, path, name=d["name"], version=d["version"],
                         taxon=Taxon.from_dict(d["taxon"], ptr(path, "taxon")), locus=d["locus"],
                         chain_type=d["chain_type"], build_date=d["build_date"])


@dataclass(frozen=True)
class SelectionEntry:
    """One selected executable allele, bound to an exact record and normalized coding-sequence digest."""

    segment: str
    label: str
    record: Locator
    sequence_sha256: str
    reason: str

    KEYS = ("segment", "label", "record", "sequence_sha256", "reason")

    def __post_init__(self):
        check_enum(self.segment, SEGMENTS, "/segment")
        check_label(self.label, "/label")
        check_instance(self.record, Locator, "/record")
        check_digest(self.sequence_sha256, "/sequence_sha256")
        check_text(self.reason, "/reason", max_len=MAX_TEXT)

    @property
    def allele_id(self) -> str:
        """The runtime allele id, ``"<segment>:<label>"``."""
        return f"{self.segment}:{self.label}"

    def to_dict(self) -> dict:
        return {"segment": self.segment, "label": self.label, "record": self.record.to_dict(),
                "sequence_sha256": self.sequence_sha256, "reason": self.reason}

    @classmethod
    def from_dict(cls, d, path: str = "") -> "SelectionEntry":
        require_object(d, cls.KEYS, path)
        return construct(cls, path, segment=d["segment"], label=d["label"],
                         record=Locator.from_dict(d["record"], ptr(path, "record")),
                         sequence_sha256=d["sequence_sha256"], reason=d["reason"])


@dataclass(frozen=True)
class Decision:
    """One ledger entry: an action on exact findings/claims bound to a subject sequence digest.

    ``select_alternative`` needs findings and no claims; ``withhold_claim`` needs claims;
    ``accept_claim`` needs claims and evidence locators.  No other action exists.
    """

    action: str
    findings: tuple
    claims: tuple
    sequence_sha256: str
    reason: str
    evidence: tuple

    KEYS = ("action", "findings", "claims", "sequence_sha256", "reason", "evidence")

    def __post_init__(self):
        check_enum(self.action, DECISION_ACTIONS, "/action")
        check_tuple(self.findings, "/findings", MAX_REFS)
        for i, finding in enumerate(self.findings):
            check_token(finding, ptr("/findings", i), FINDING_ID_RE, "a finding id finding:<64 hex>")
        check_unique(list(self.findings), "/findings", "finding id")
        check_tuple(self.claims, "/claims", MAX_REFS)
        for i, claim in enumerate(self.claims):
            check_token(claim, ptr("/claims", i), CLAIM_ID_RE, "a claim id claim:<64 hex>")
        check_unique(list(self.claims), "/claims", "claim id")
        check_digest(self.sequence_sha256, "/sequence_sha256")
        check_text(self.reason, "/reason", max_len=MAX_TEXT)
        locator_tuple(self.evidence, "/evidence", MAX_REFS, canonical=False)
        if self.action == "select_alternative":
            if not self.findings:
                raise bad("/findings", "select_alternative names the findings it resolves")
            if self.claims:
                raise bad("/claims", "select_alternative chooses a supplied sequence, not claims")
        elif not self.claims:
            raise bad("/claims", f"{self.action} names the exact claims it decides")
        if self.action == "accept_claim" and not self.evidence:
            raise bad("/evidence", "accept_claim needs the evidence locators that support the annotation")

    def to_dict(self) -> dict:
        return {"action": self.action, "findings": list(self.findings), "claims": list(self.claims),
                "sequence_sha256": self.sequence_sha256, "reason": self.reason,
                "evidence": [loc.to_dict() for loc in self.evidence]}

    @classmethod
    def from_dict(cls, d, path: str = "") -> "Decision":
        require_object(d, cls.KEYS, path)
        return construct(
            cls, path, action=d["action"],
            findings=parse_array(d["findings"], ptr(path, "findings"), MAX_REFS, lambda v, p: v),
            claims=parse_array(d["claims"], ptr(path, "claims"), MAX_REFS, lambda v, p: v),
            sequence_sha256=d["sequence_sha256"], reason=d["reason"],
            evidence=parse_array(d["evidence"], ptr(path, "evidence"), MAX_REFS, Locator.from_dict))


@dataclass(frozen=True)
class Requirement:
    """A required capability and minimum evidence, in the runtime evaluator's vocabulary."""

    capability: str
    evidence: str

    KEYS = ("capability", "evidence")

    def __post_init__(self):
        check_enum(self.capability, CAPABILITIES, "/capability")
        check_enum(self.evidence, EVIDENCE_LEVELS, "/evidence")

    def to_dict(self) -> dict:
        return {"capability": self.capability, "evidence": self.evidence}

    @classmethod
    def from_dict(cls, d, path: str = "") -> "Requirement":
        require_object(d, cls.KEYS, path)
        return construct(cls, path, capability=d["capability"], evidence=d["evidence"])


@dataclass(frozen=True)
class Recipe:
    source_lock_hash: str
    adapter: AdapterRef
    projection: ProjectionRef
    reference: RecipeReference
    selection: tuple
    decisions: tuple
    requires: tuple
    defaults_profile: str

    def __post_init__(self):
        check_digest(self.source_lock_hash, "/source_lock_hash")
        check_instance(self.adapter, AdapterRef, "/adapter")
        check_instance(self.projection, ProjectionRef, "/projection")
        check_instance(self.reference, RecipeReference, "/reference")
        check_tuple(self.selection, "/selection", MAX_SELECTION, SelectionEntry)
        check_tuple(self.decisions, "/decisions", MAX_DECISIONS, Decision)
        check_tuple(self.requires, "/requires", MAX_REQUIRES, Requirement)
        check_defaults_profile(self.defaults_profile, "/defaults_profile")
        allowed = self.reference.identity.segments
        for i, entry in enumerate(self.selection):
            if entry.segment not in allowed:
                raise bad(ptr(ptr("/selection", i), "segment"),
                          f"{entry.label}: segment {entry.segment} is not part of a {self.reference.chain_type} "
                          f"reference", "identity_mismatch")
        check_unique([(e.segment, e.label) for e in self.selection], "/selection", "selected label",
                     "unresolved_selection")
        check_unique([e.record.id for e in self.selection], "/selection", "selected record",
                     "unresolved_selection")
        check_unique([r.capability for r in self.requires], "/requires", "required capability")
        resolved: dict = {}
        for i, decision in enumerate(self.decisions):
            if decision.action == "select_alternative":
                for finding in decision.findings:
                    if finding in resolved:
                        raise bad("/decisions", f"finding {finding} is resolved by select_alternative decisions "
                                                f"{resolved[finding]} and {i}")
                    resolved[finding] = i
        decided = {}
        for i, decision in enumerate(self.decisions):
            for claim in decision.claims:
                if claim in decided:
                    raise bad("/decisions", f"claim {claim} is decided by decisions {decided[claim]} and {i}")
                decided[claim] = i

    def to_dict(self) -> dict:
        return {"artifact": ARTIFACT, "schema_version": SCHEMA_VERSION, "source_lock_hash": self.source_lock_hash,
                "adapter": self.adapter.to_dict(), "projection": self.projection.to_dict(),
                "reference": self.reference.to_dict(), "selection": [e.to_dict() for e in self.selection],
                "decisions": [x.to_dict() for x in self.decisions],
                "requires": [r.to_dict() for r in self.requires], "defaults_profile": self.defaults_profile}

    @classmethod
    def from_dict(cls, d, path: str = "") -> "Recipe":
        check_header(d, ARTIFACT, SCHEMA_VERSION, path)
        require_object(d, _KEYS, path)
        return construct(
            cls, path, source_lock_hash=d["source_lock_hash"],
            adapter=AdapterRef.from_dict(d["adapter"], ptr(path, "adapter")),
            projection=ProjectionRef.from_dict(d["projection"], ptr(path, "projection")),
            reference=RecipeReference.from_dict(d["reference"], ptr(path, "reference")),
            selection=parse_array(d["selection"], ptr(path, "selection"), MAX_SELECTION, SelectionEntry.from_dict),
            decisions=parse_array(d["decisions"], ptr(path, "decisions"), MAX_DECISIONS, Decision.from_dict),
            requires=parse_array(d["requires"], ptr(path, "requires"), MAX_REQUIRES, Requirement.from_dict),
            defaults_profile=d["defaults_profile"])

    @classmethod
    def from_json(cls, text: str) -> "Recipe":
        return cls.from_dict(parse_json(text, max_bytes=MAX_RECIPE_BYTES))

    @classmethod
    def read(cls, path) -> "Recipe":
        return cls.from_dict(read_json(path, max_bytes=MAX_RECIPE_BYTES))
