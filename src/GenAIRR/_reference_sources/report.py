"""Build report ``genairr_reference_source_build_report`` version 1 (closed shape).

A report embeds the exact lock, recipe and source inventory it was built from and adds
the build's accounting: the four input hashes, the runtime reference hashes, the actual
generator, endpoint-comparison completeness, distinct counts (record occurrences,
sequence candidates, mentioned-only identities, selected executable alleles), the
ordered selection, one outcome per explicit source record (selected, excluded with a
generated reason, or rejected as malformed/unsupported), one outcome per claim
(projected, withheld, not representable at runtime, or about an unselected sequence),
runtime capabilities with per-allele blockers, and each required capability's status.

Construction verifies every embedded document, every cross-document reference, every
recomputable hash and every count, so a rehashed but inconsistent report is rejected.
The report hash is computed externally and never stored inside the report.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

from GenAIRR._reference.capabilities import CAPABILITIES
from GenAIRR._reference_sources._fields import (
    MAX_NAME, MAX_TEXT, bad, check_bool, check_canonical_order, check_digest, check_enum, check_instance,
    check_int, check_text, check_token, check_tuple, construct, parse_array, parse_json, ptr, read_json,
    require_object, check_header,
)
from GenAIRR._reference_sources.common import (
    CLAIM_ID_RE, FINDING_ID_RE, FORMATS, RECORD_EXPORT_FORMATS, SEGMENTS, Generator, Locator, RuntimeHashes,
)
from GenAIRR._reference_sources.consistency import (
    check_decisions, check_inventory, check_recipe, check_selection, selected_subjects,
)
from GenAIRR._reference_sources.digests import build_input_hash, inventory_hash, lock_hash, recipe_hash
from GenAIRR._reference_sources.errors import ReferenceSourceError
from GenAIRR._reference_sources.inventory import (
    MAX_CLAIMS, MAX_FINDINGS, MAX_RECORDS, RUNTIME_PREDICATES, SourceInventory,
)
from GenAIRR._reference_sources.lock import SourceLock
from GenAIRR._reference_sources.recipe import EVIDENCE_LEVELS, MAX_REQUIRES, MAX_SELECTION, Recipe

ARTIFACT = "genairr_reference_source_build_report"
SCHEMA_VERSION = 1
MAX_REPORT_BYTES = 256 * 1024 * 1024
MAX_UNMET = 65536

RECORD_OUTCOMES = ("selected", "excluded", "rejected")
EXCLUSION_REASONS = ("alternative_for_selected_identity", "duplicate_occurrence", "structural_incompatibility",
                     "missing_required_evidence", "unresolved_conflicting_sequence", "outside_explicit_selection")
# Exclusions that point at the selected record of the same identity.
RELATED_EXCLUSIONS = ("alternative_for_selected_identity", "duplicate_occurrence")
CLAIM_OUTCOMES = ("projected", "withheld", "not_representable", "not_selected")
PHASE_SUPPORTS = ("none", "anchor_derived", "independent_source", "independent_assertion")

_EVIDENCE_RANK = {"assertion": 1, "source": 2}
_SEGMENT_KEYS = ("V", "D", "J", "unassigned")
_KEYS = ("artifact", "schema_version", "generator", "source_lock", "recipe", "inventory", "hashes", "runtime",
         "comparison", "counts", "selection", "records", "claims", "capabilities", "requirements")


def _closed(cls, d, path, **nested):
    """Generic ``from_dict``: closed keys, nested parsers by key, the rest passed through."""
    require_object(d, cls.KEYS, path)
    fields = {k: (nested[k](d[k], ptr(path, k)) if k in nested else d[k]) for k in cls.KEYS}
    return construct(cls, path, **fields)


def _under(prefix: str, check, *args) -> None:
    try:
        check(*args)
    except ReferenceSourceError as exc:
        raise exc.under(prefix) from None


# ── Small shapes ────────────────────────────────────────────────────


@dataclass(frozen=True)
class BuildHashes:
    """The input identities: lock, inventory, recipe and build-input hashes."""

    source_lock: str
    inventory: str
    recipe: str
    build_inputs: str

    KEYS = ("source_lock", "inventory", "recipe", "build_inputs")

    def __post_init__(self):
        for key in self.KEYS:
            check_digest(getattr(self, key), "/" + key)

    def to_dict(self) -> dict:
        return {k: getattr(self, k) for k in self.KEYS}

    @classmethod
    def from_dict(cls, d, path: str = "") -> "BuildHashes":
        return _closed(cls, d, path)


@dataclass(frozen=True)
class Comparison:
    """Endpoint comparison: which record exports are present, which are missing, whether
    the cross-variant comparison is complete, and the findings it produced (sorted ids)."""

    complete: bool
    exports: tuple
    missing: tuple
    findings: tuple

    KEYS = ("complete", "exports", "missing", "findings")

    def __post_init__(self):
        check_bool(self.complete, "/complete")
        for name in ("exports", "missing"):
            items = check_tuple(getattr(self, name), "/" + name, len(RECORD_EXPORT_FORMATS))
            for i, fmt in enumerate(items):
                check_enum(fmt, RECORD_EXPORT_FORMATS, ptr("/" + name, i))
            check_canonical_order([FORMATS.index(f) for f in items], "/" + name, "record export")
        if self.complete != (not self.missing):
            raise bad("/complete", "the comparison is complete exactly when no record export is missing")
        check_tuple(self.findings, "/findings", MAX_FINDINGS)
        for i, finding in enumerate(self.findings):
            check_token(finding, ptr("/findings", i), FINDING_ID_RE, "a finding id finding:<64 hex>")
        check_canonical_order(list(self.findings), "/findings", "finding id")

    def check(self, lock: SourceLock, inventory: Optional[SourceInventory]) -> None:
        exports = tuple(a.asset_id for a in lock.assets if a.role == "record_export")
        if self.exports != exports:
            raise bad("/exports", f"the lock holds record exports {list(exports)}, not {list(self.exports)}")
        missing = tuple(f for f in RECORD_EXPORT_FORMATS if f not in exports)
        if self.missing != missing:
            raise bad("/missing", f"the lock lacks record exports {list(missing)}, not {list(self.missing)}")
        if inventory is not None:
            for i, finding in enumerate(self.findings):
                try:
                    inventory.finding(finding)
                except KeyError:
                    raise bad(ptr("/findings", i), f"unknown finding {finding}") from None

    def to_dict(self) -> dict:
        return {"complete": self.complete, "exports": list(self.exports), "missing": list(self.missing),
                "findings": list(self.findings)}

    @classmethod
    def from_dict(cls, d, path: str = "") -> "Comparison":
        def strings(value, where):
            return parse_array(value, where, MAX_FINDINGS, lambda v, p: v)
        return _closed(cls, d, path, exports=strings, missing=strings, findings=strings)


@dataclass(frozen=True)
class SegmentCounts:
    V: int
    D: int
    J: int
    unassigned: int

    KEYS = _SEGMENT_KEYS

    def __post_init__(self):
        for key in self.KEYS:
            check_int(getattr(self, key), "/" + key, minimum=0)

    def to_dict(self) -> dict:
        return {k: getattr(self, k) for k in self.KEYS}

    @classmethod
    def from_dict(cls, d, path: str = "") -> "SegmentCounts":
        return _closed(cls, d, path)


@dataclass(frozen=True)
class ExportCounts:
    """Per record export: record occurrences, sequence candidates and mentioned-only identities."""

    asset: str
    record_occurrences: SegmentCounts
    sequence_candidates: SegmentCounts
    mentioned_identities: SegmentCounts

    KEYS = ("asset", "record_occurrences", "sequence_candidates", "mentioned_identities")

    def __post_init__(self):
        check_enum(self.asset, RECORD_EXPORT_FORMATS, "/asset")
        for key in self.KEYS[1:]:
            check_instance(getattr(self, key), SegmentCounts, "/" + key)

    def to_dict(self) -> dict:
        return {"asset": self.asset, **{k: getattr(self, k).to_dict() for k in self.KEYS[1:]}}

    @classmethod
    def from_dict(cls, d, path: str = "") -> "ExportCounts":
        parse = SegmentCounts.from_dict
        return _closed(cls, d, path, record_occurrences=parse, sequence_candidates=parse,
                       mentioned_identities=parse)


@dataclass(frozen=True)
class OutcomeCounts:
    selected: int
    excluded: int
    rejected: int

    KEYS = RECORD_OUTCOMES

    def __post_init__(self):
        for key in self.KEYS:
            check_int(getattr(self, key), "/" + key, minimum=0)

    def to_dict(self) -> dict:
        return {k: getattr(self, k) for k in self.KEYS}

    @classmethod
    def from_dict(cls, d, path: str = "") -> "OutcomeCounts":
        return _closed(cls, d, path)


@dataclass(frozen=True)
class Counts:
    """Distinct counts: per-export accounting, selected executable alleles, record outcomes."""

    exports: tuple
    selected_alleles: SegmentCounts
    record_outcomes: OutcomeCounts

    KEYS = ("exports", "selected_alleles", "record_outcomes")

    def __post_init__(self):
        check_tuple(self.exports, "/exports", len(RECORD_EXPORT_FORMATS), ExportCounts)
        check_instance(self.selected_alleles, SegmentCounts, "/selected_alleles")
        check_instance(self.record_outcomes, OutcomeCounts, "/record_outcomes")

    def to_dict(self) -> dict:
        return {"exports": [e.to_dict() for e in self.exports], "selected_alleles": self.selected_alleles.to_dict(),
                "record_outcomes": self.record_outcomes.to_dict()}

    @classmethod
    def from_dict(cls, d, path: str = "") -> "Counts":
        return _closed(cls, d, path,
                       exports=lambda v, p: parse_array(v, p, len(RECORD_EXPORT_FORMATS), ExportCounts.from_dict),
                       selected_alleles=SegmentCounts.from_dict, record_outcomes=OutcomeCounts.from_dict)


@dataclass(frozen=True)
class SelectedAllele:
    """One selection entry's runtime allele id and the exact candidate it uses (recipe order)."""

    allele: str
    candidate: Locator

    KEYS = ("allele", "candidate")

    def __post_init__(self):
        check_text(self.allele, "/allele", max_len=MAX_NAME + 2, multiline=False)
        check_instance(self.candidate, Locator, "/candidate")

    def to_dict(self) -> dict:
        return {"allele": self.allele, "candidate": self.candidate.to_dict()}

    @classmethod
    def from_dict(cls, d, path: str = "") -> "SelectedAllele":
        return _closed(cls, d, path, candidate=Locator.from_dict)


@dataclass(frozen=True)
class RecordOutcome:
    """The build outcome of one explicit source record.

    ``selected`` names its runtime ``allele``; ``excluded`` gives a generated
    ``exclusion`` reason (with ``related``, the selected record of the same identity, for
    alternatives and duplicate occurrences); ``rejected`` marks a malformed or
    unsupported source record.  ``detail`` is optional free text.
    """

    record: Locator
    outcome: str
    allele: Optional[str]
    exclusion: Optional[str]
    related: Optional[Locator]
    detail: Optional[str]

    KEYS = ("record", "outcome", "allele", "exclusion", "related", "detail")

    def __post_init__(self):
        check_instance(self.record, Locator, "/record")
        check_enum(self.outcome, RECORD_OUTCOMES, "/outcome")
        selected, excluded = self.outcome == "selected", self.outcome == "excluded"
        check_text(self.allele, "/allele", nullable=not selected, max_len=MAX_NAME + 2, multiline=False)
        if not selected and self.allele is not None:
            raise bad("/allele", f"a {self.outcome} record has no runtime allele")
        check_enum(self.exclusion, EXCLUSION_REASONS, "/exclusion", nullable=not excluded)
        if not excluded and self.exclusion is not None:
            raise bad("/exclusion", f"a {self.outcome} record has no exclusion reason")
        check_instance(self.related, Locator, "/related", nullable=True)
        if (self.exclusion in RELATED_EXCLUSIONS) != (self.related is not None):
            raise bad("/related", f"exactly the exclusions {list(RELATED_EXCLUSIONS)} name the related selected record")
        check_text(self.detail, "/detail", nullable=True, max_len=MAX_TEXT)

    def to_dict(self) -> dict:
        return {"record": self.record.to_dict(), "outcome": self.outcome, "allele": self.allele,
                "exclusion": self.exclusion, "related": None if self.related is None else self.related.to_dict(),
                "detail": self.detail}

    @classmethod
    def from_dict(cls, d, path: str = "") -> "RecordOutcome":
        return _closed(cls, d, path, record=Locator.from_dict,
                       related=lambda v, p: None if v is None else Locator.from_dict(v, p))


@dataclass(frozen=True)
class ClaimOutcome:
    """What the projection did with one inventory claim; ``decision`` is a recipe decision index."""

    claim: str
    outcome: str
    decision: Optional[int]
    reason: Optional[str]

    KEYS = ("claim", "outcome", "decision", "reason")

    def __post_init__(self):
        check_token(self.claim, "/claim", CLAIM_ID_RE, "a claim id claim:<64 hex>")
        check_enum(self.outcome, CLAIM_OUTCOMES, "/outcome")
        check_int(self.decision, "/decision", nullable=True, minimum=0)
        check_text(self.reason, "/reason", nullable=self.outcome != "withheld", max_len=MAX_TEXT)

    def to_dict(self) -> dict:
        return {"claim": self.claim, "outcome": self.outcome, "decision": self.decision, "reason": self.reason}

    @classmethod
    def from_dict(cls, d, path: str = "") -> "ClaimOutcome":
        return _closed(cls, d, path)


@dataclass(frozen=True)
class UnmetRequirement:
    """A per-allele blocker and the evidence that would resolve it (runtime ``Unmet`` shape)."""

    allele: str
    requirement: str
    hint: str

    KEYS = ("allele", "requirement", "hint")

    def __post_init__(self):
        check_text(self.allele, "/allele", max_len=MAX_NAME + 2, multiline=False)
        check_text(self.requirement, "/requirement", max_len=MAX_TEXT)
        check_text(self.hint, "/hint", max_len=MAX_TEXT)

    def to_dict(self) -> dict:
        return {"allele": self.allele, "requirement": self.requirement, "hint": self.hint}

    @classmethod
    def from_dict(cls, d, path: str = "") -> "UnmetRequirement":
        return _closed(cls, d, path)


@dataclass(frozen=True)
class CapabilityStatus:
    """One runtime capability (runtime ``CapabilityStatus`` shape)."""

    name: str
    available: bool
    evidence: Optional[str]
    unmet: tuple
    note: Optional[str]

    KEYS = ("name", "available", "evidence", "unmet", "note")

    def __post_init__(self):
        check_enum(self.name, CAPABILITIES, "/name")
        check_bool(self.available, "/available")
        check_enum(self.evidence, EVIDENCE_LEVELS, "/evidence", nullable=True)
        if not self.available and self.evidence is not None:
            raise bad("/evidence", "an unavailable capability has no evidence level")
        check_tuple(self.unmet, "/unmet", MAX_UNMET, UnmetRequirement)
        if self.available and self.unmet:
            raise bad("/unmet", "an available capability has no unmet requirements")
        check_text(self.note, "/note", nullable=True, max_len=MAX_TEXT)

    def to_dict(self) -> dict:
        return {"name": self.name, "available": self.available, "evidence": self.evidence,
                "unmet": [u.to_dict() for u in self.unmet], "note": self.note}

    @classmethod
    def from_dict(cls, d, path: str = "") -> "CapabilityStatus":
        return _closed(cls, d, path, unmet=lambda v, p: parse_array(v, p, MAX_UNMET, UnmetRequirement.from_dict))


@dataclass(frozen=True)
class RuntimeCapabilities:
    """Every runtime capability in evaluator order, plus phase support (runtime report shape)."""

    statuses: tuple
    phase_support: str

    KEYS = ("statuses", "phase_support")

    def __post_init__(self):
        check_tuple(self.statuses, "/statuses", len(CAPABILITIES), CapabilityStatus)
        if tuple(s.name for s in self.statuses) != CAPABILITIES:
            raise bad("/statuses", f"statuses list every capability in order {list(CAPABILITIES)}")
        check_enum(self.phase_support, PHASE_SUPPORTS, "/phase_support")

    def status(self, name: str) -> CapabilityStatus:
        return self.statuses[CAPABILITIES.index(name)]

    def to_dict(self) -> dict:
        return {"statuses": [s.to_dict() for s in self.statuses], "phase_support": self.phase_support}

    @classmethod
    def from_dict(cls, d, path: str = "") -> "RuntimeCapabilities":
        return _closed(cls, d, path,
                       statuses=lambda v, p: parse_array(v, p, len(CAPABILITIES), CapabilityStatus.from_dict))


@dataclass(frozen=True)
class RequirementOutcome:
    """A recipe requirement and whether the runtime capabilities satisfy it."""

    capability: str
    evidence: str
    satisfied: bool

    KEYS = ("capability", "evidence", "satisfied")

    def __post_init__(self):
        check_enum(self.capability, CAPABILITIES, "/capability")
        check_enum(self.evidence, EVIDENCE_LEVELS, "/evidence")
        check_bool(self.satisfied, "/satisfied")

    def to_dict(self) -> dict:
        return {"capability": self.capability, "evidence": self.evidence, "satisfied": self.satisfied}

    @classmethod
    def from_dict(cls, d, path: str = "") -> "RequirementOutcome":
        return _closed(cls, d, path)


# ── The report ──────────────────────────────────────────────────────


@dataclass(frozen=True)
class BuildReport:
    generator: Generator
    source_lock: SourceLock
    recipe: Recipe
    inventory: SourceInventory
    hashes: BuildHashes
    runtime: RuntimeHashes
    comparison: Comparison
    counts: Counts
    selection: tuple
    records: tuple
    claims: tuple
    capabilities: RuntimeCapabilities
    requirements: tuple

    def __post_init__(self):
        for name, cls in (("generator", Generator), ("source_lock", SourceLock), ("recipe", Recipe),
                          ("inventory", SourceInventory), ("hashes", BuildHashes), ("runtime", RuntimeHashes),
                          ("comparison", Comparison), ("counts", Counts), ("capabilities", RuntimeCapabilities)):
            check_instance(getattr(self, name), cls, "/" + name)
        check_tuple(self.selection, "/selection", MAX_SELECTION, SelectedAllele)
        check_tuple(self.records, "/records", MAX_RECORDS, RecordOutcome)
        check_tuple(self.claims, "/claims", MAX_CLAIMS, ClaimOutcome)
        check_tuple(self.requirements, "/requirements", MAX_REQUIRES, RequirementOutcome)
        lock, recipe, inventory = self.source_lock, self.recipe, self.inventory
        _under("/recipe", check_recipe, recipe, lock)
        _under("/inventory", check_inventory, inventory, lock)
        _under("/recipe", check_selection, recipe, inventory)
        _under("/recipe", check_decisions, recipe, inventory)
        expected = {"source_lock": lock_hash(lock), "inventory": inventory_hash(inventory),
                    "recipe": recipe_hash(recipe), "build_inputs": build_input_hash(lock, inventory, recipe)}
        for key, value in expected.items():
            if getattr(self.hashes, key) != value:
                raise bad(ptr("/hashes", key), f"stored {getattr(self.hashes, key)} differs from the recomputed "
                                               f"{value}", "digest_mismatch")
        _under("/comparison", self.comparison.check, lock, inventory)
        self._check_selection()
        self._check_records()
        self._check_claims()
        self._check_counts()
        self._check_requirements()

    # -- consistency of the build accounting ----------------------------

    def _check_selection(self) -> None:
        if len(self.selection) != len(self.recipe.selection):
            raise bad("/selection", f"{len(self.selection)} selected alleles for {len(self.recipe.selection)} "
                                    "recipe selection entries")
        for i, (chosen, entry) in enumerate(zip(self.selection, self.recipe.selection)):
            if chosen.allele != entry.allele_id:
                raise bad(ptr(ptr("/selection", i), "allele"), f"expected allele {entry.allele_id!r}")
            try:
                candidate = self.inventory.candidate(chosen.candidate)
            except KeyError:
                candidate = None
            if candidate is None or candidate.record != entry.record or candidate.sha256 != entry.sequence_sha256:
                raise bad(ptr(ptr("/selection", i), "candidate"), f"{chosen.candidate.id} is not the candidate of "
                                                                 f"{entry.record.id} with digest {entry.sequence_sha256}")

    def _check_records(self) -> None:
        inventory = self.inventory
        if len(self.records) != len(inventory.records):
            raise bad("/records", f"{len(self.records)} outcomes for {len(inventory.records)} source records")
        entries = {e.record.id: e for e in self.recipe.selection}
        identities = {(e.segment, e.label): e for e in self.recipe.selection}
        for i, (outcome, record) in enumerate(zip(self.records, inventory.records)):
            path = ptr("/records", i)
            if outcome.record != record.locator:
                raise bad(ptr(path, "record"), f"expected the outcome of {record.id} (inventory order)")
            entry = entries.get(record.id)
            expected = "rejected" if record.status != "parsed" else "selected" if entry else "excluded"
            if outcome.outcome != expected:
                raise bad(ptr(path, "outcome"), f"record {record.id} ({record.status}) is {expected}, "
                                                f"not {outcome.outcome}")
            if entry is not None and outcome.allele != entry.allele_id:
                raise bad(ptr(path, "allele"), f"expected allele {entry.allele_id!r}")
            if (expected == "excluded" and outcome.exclusion not in RELATED_EXCLUSIONS
                    and (record.segment, record.label) in identities):
                raise bad(ptr(path, "exclusion"), f"{record.label} is the identity selected at "
                                                  f"{identities[(record.segment, record.label)].record.id}; it is "
                                                  f"excluded as one of {list(RELATED_EXCLUSIONS)}")
            if outcome.related is None:
                continue
            related = entries.get(outcome.related.id)
            if related is None:
                raise bad(ptr(path, "related"), f"{outcome.related.id} is not a selected record")
            if (record.segment, record.label) != (related.segment, related.label):
                raise bad(ptr(path, "related"), f"{record.label} is not the identity {related.label} selected "
                                                f"at {related.record.id}")
            same = any(c.sha256 == related.sequence_sha256 for c in inventory.candidates_of(record.locator))
            if same != (outcome.exclusion == "duplicate_occurrence"):
                raise bad(ptr(path, "exclusion"), "a duplicate occurrence supplies the selected sequence; an "
                                                  "alternative record supplies a different one")

    def _check_claims(self) -> None:
        inventory, decisions = self.inventory, self.recipe.decisions
        if len(self.claims) != len(inventory.claims):
            raise bad("/claims", f"{len(self.claims)} outcomes for {len(inventory.claims)} claims")
        selected = selected_subjects(self.recipe, inventory)
        decided = {}
        for i, (outcome, claim) in enumerate(zip(self.claims, inventory.claims)):
            path = ptr("/claims", i)
            if outcome.claim != claim.id:
                raise bad(ptr(path, "claim"), f"expected the outcome of {claim.id} (inventory order)")
            if claim.subject.id not in selected:
                allowed = ("not_selected",)
            elif claim.predicate not in RUNTIME_PREDICATES:
                allowed = ("not_representable",)
            else:
                allowed = ("projected", "withheld")
            if outcome.outcome not in allowed:
                raise bad(ptr(path, "outcome"), f"claim {claim.id} ({claim.predicate} of {claim.subject.id}) can only "
                                                f"be {' or '.join(allowed)}")
            if outcome.outcome == "projected" and claim.status == "unresolved":
                raise bad(ptr(path, "outcome"), f"claim {claim.id} is unresolved and has no value to project")
            if outcome.decision is not None:
                if outcome.decision >= len(decisions) or claim.id not in decisions[outcome.decision].claims:
                    raise bad(ptr(path, "decision"), f"recipe decision {outcome.decision} does not decide {claim.id}")
                action = decisions[outcome.decision].action
                if (action, outcome.outcome) not in (("withhold_claim", "withheld"), ("accept_claim", "projected")):
                    raise bad(ptr(path, "decision"), f"a {action} decision cannot leave the claim {outcome.outcome}")
                decided[claim.id] = outcome.decision
        for index, decision in enumerate(decisions):
            for claim_id in decision.claims:
                if decided.get(claim_id) != index:
                    raise bad("/claims", f"recipe decision {index} ({decision.action}) on {claim_id} is not applied")

    def _check_counts(self) -> None:
        counts, inventory = self.counts, self.inventory
        exports = tuple(ExportCounts(asset, **{k: SegmentCounts(**v) for k, v in inventory.export_counts(asset).items()})
                        for asset in self.comparison.exports)
        if counts.exports != exports:
            raise bad("/counts/exports", "per-export counts do not match the inventory: expected "
                                         f"{[e.to_dict() for e in exports]}")
        selected = dict.fromkeys(_SEGMENT_KEYS, 0)
        for entry in self.recipe.selection:
            selected[entry.segment] += 1
        if counts.selected_alleles != SegmentCounts(**selected):
            raise bad("/counts/selected_alleles", f"the recipe selects {selected}")
        outcomes = {k: sum(r.outcome == k for r in self.records) for k in RECORD_OUTCOMES}
        if counts.record_outcomes != OutcomeCounts(**outcomes):
            raise bad("/counts/record_outcomes", f"the record outcomes are {outcomes}")

    def _check_requirements(self) -> None:
        requires = self.recipe.requires
        if len(self.requirements) != len(requires):
            raise bad("/requirements", f"{len(self.requirements)} outcomes for {len(requires)} recipe requirements")
        for i, (outcome, requirement) in enumerate(zip(self.requirements, requires)):
            path = ptr("/requirements", i)
            if (outcome.capability, outcome.evidence) != (requirement.capability, requirement.evidence):
                raise bad(path, f"expected requirement {requirement.to_dict()}")
            status = self.capabilities.status(requirement.capability)
            satisfied = status.available and (status.evidence is None or _EVIDENCE_RANK[status.evidence]
                                              >= _EVIDENCE_RANK[requirement.evidence])
            if outcome.satisfied != satisfied:
                raise bad(ptr(path, "satisfied"), f"{requirement.capability} at evidence {requirement.evidence} is "
                                                  f"{'satisfied' if satisfied else 'not satisfied'}")

    # -- serialization -------------------------------------------------

    def to_dict(self) -> dict:
        return {"artifact": ARTIFACT, "schema_version": SCHEMA_VERSION, "generator": self.generator.to_dict(),
                "source_lock": self.source_lock.to_dict(), "recipe": self.recipe.to_dict(),
                "inventory": self.inventory.to_dict(), "hashes": self.hashes.to_dict(),
                "runtime": self.runtime.to_dict(), "comparison": self.comparison.to_dict(),
                "counts": self.counts.to_dict(), "selection": [s.to_dict() for s in self.selection],
                "records": [r.to_dict() for r in self.records], "claims": [c.to_dict() for c in self.claims],
                "capabilities": self.capabilities.to_dict(),
                "requirements": [r.to_dict() for r in self.requirements]}

    @classmethod
    def from_dict(cls, d, path: str = "") -> "BuildReport":
        check_header(d, ARTIFACT, SCHEMA_VERSION, path)
        require_object(d, _KEYS, path)
        return construct(
            cls, path,
            generator=Generator.from_dict(d["generator"], ptr(path, "generator")),
            source_lock=SourceLock.from_dict(d["source_lock"], ptr(path, "source_lock")),
            recipe=Recipe.from_dict(d["recipe"], ptr(path, "recipe")),
            inventory=SourceInventory.from_dict(d["inventory"], ptr(path, "inventory")),
            hashes=BuildHashes.from_dict(d["hashes"], ptr(path, "hashes")),
            runtime=RuntimeHashes.from_dict(d["runtime"], ptr(path, "runtime")),
            comparison=Comparison.from_dict(d["comparison"], ptr(path, "comparison")),
            counts=Counts.from_dict(d["counts"], ptr(path, "counts")),
            selection=parse_array(d["selection"], ptr(path, "selection"), MAX_SELECTION, SelectedAllele.from_dict),
            records=parse_array(d["records"], ptr(path, "records"), MAX_RECORDS, RecordOutcome.from_dict),
            claims=parse_array(d["claims"], ptr(path, "claims"), MAX_CLAIMS, ClaimOutcome.from_dict),
            capabilities=RuntimeCapabilities.from_dict(d["capabilities"], ptr(path, "capabilities")),
            requirements=parse_array(d["requirements"], ptr(path, "requirements"), MAX_REQUIRES,
                                     RequirementOutcome.from_dict))

    @classmethod
    def from_json(cls, text: str) -> "BuildReport":
        return cls.from_dict(parse_json(text, max_bytes=MAX_REPORT_BYTES))

    @classmethod
    def read(cls, path) -> "BuildReport":
        return cls.from_dict(read_json(path, max_bytes=MAX_REPORT_BYTES))
