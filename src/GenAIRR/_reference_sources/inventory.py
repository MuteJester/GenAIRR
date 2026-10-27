"""Normalized source inventory: internal model ``genairr_reference_source_inventory`` version 1.

The inventory is generated from locked bytes by an adapter.  It keeps distinct:

* :class:`SourceRecord` - one supplied record occurrence (including malformed and
  unsupported records), identified by asset plus JSON pointer, never by provider id;
* :class:`SequenceCandidate` - a supplied coding-sequence alternative bound to its record;
* :class:`MentionedIdentity` - a name that appears in an export only through alias or
  paralog relationships, with no sequence, metadata or sampling weight;
* :class:`Observation` - an original value, its field locator and interpretation status;
* :class:`SupportReference` - an embedded, unretrieved provider support reference;
* :class:`Mapping` - an exact coordinate mapping between sequences, variant-specific;
* :class:`Claim` - a build-model claim with support, derivation and status;
* :class:`Finding` - a stable-coded finding with subjects, locators and affected capabilities.

Claims, mappings and findings carry content-derived ids (``claim:``/``mapping:``/
``finding:`` plus 64 hex), verified on reading.  The hashed content includes the
``source_lock_hash`` each item was derived under, so any change to the item or to the
locked source bytes changes its id and a decision naming it goes stale.
Every list is held in canonical order (see :meth:`SourceInventory.assemble`).
"""
from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
import re
from typing import Any, Optional

from GenAIRR._reference.capabilities import CAPABILITIES
from GenAIRR._reference.model import PREDICATES as RUNTIME_PREDICATES, SUBREGION_LABELS, sequence_digest
from GenAIRR._reference_sources._fields import (
    MAX_SOURCE_RECORDS, MAX_NAME, MAX_REFS, MAX_SEQUENCE, MAX_TAXONOMY_ID, MAX_TEXT, JsonValue, bad, check_bool,
    check_canonical_order, check_digest, check_enum, check_instance, check_int, check_label, check_text,
    check_token, check_tuple, construct, digest, parse_array, parse_json, ptr, read_json, require_object,
    check_header, show,
)
from GenAIRR._reference_sources.common import (
    CLAIM_ID_RE, FINDING_ID_RE, FORMAT_ROLES, FORMATS, GAPPED_VARIANT, MAPPING_ID_RE, SEGMENTS, AdapterRef, Locator,
    check_asset_id, locator_tuple,
)

ARTIFACT = "genairr_reference_source_inventory"
SCHEMA_VERSION = 1
MAX_INVENTORY_BYTES = 128 * 1024 * 1024
MAX_RECORDS = MAX_SOURCE_RECORDS
MAX_CANDIDATES = 65536
MAX_MENTIONED = 65536
MAX_OBSERVATIONS = 524288
MAX_SUPPORT = 262144
MAX_CLAIMS = 262144
MAX_MAPPINGS = 65536
MAX_FINDINGS = 65536
MAX_LOCUS = 64

CLAIM_ID_DOMAIN = "genairr_reference_source_claim_id_v1"
FINDING_ID_DOMAIN = "genairr_reference_source_finding_id_v1"
MAPPING_ID_DOMAIN = "genairr_reference_source_mapping_id_v1"

RECORD_STATUSES = ("parsed", "malformed", "unsupported")
PRESENCES = ("present", "absent")
INTERPRETATIONS = ("interpreted", "preserved", "unsupported", "invalid")
RELATIONSHIPS = ("alias", "paralog")
SUPPORT_KINDS = ("unrearranged", "rearranged")
SUPPORT_STATUSES = ("embedded_unretrieved",)

PREDICATES = ("anchor", "coding_phase", "phase_independence", "boundary_5", "boundary_3", "functionality",
              "v_subregions", "coding_interval", "extension")
CLAIM_SUPPORTS = ("source", "assertion", "none")
DERIVATIONS = ("stated", "imgt_alignment", "unique_motif", "anchor_derived", "variant_comparison")
CLAIM_STATUSES = ("accepted", "candidate", "conflicting", "unresolved")
COORDINATE_SYSTEMS = ("coding_sequence", "record_sequence", "imgt_gapped")
PHASE_INDEPENDENCE = ("independent", "anchor_derived")
FUNCTIONALITIES = ("functional", "orf", "pseudogene")
BOUNDARIES = ("complete", "partial")
EXTENSION_SIDES = ("5", "3")

# Coordinate systems a predicate's value may be expressed in; absent -> must be null.
PREDICATE_SYSTEMS = {
    "anchor": ("coding_sequence", "record_sequence", "imgt_gapped"),
    "coding_phase": ("coding_sequence", "record_sequence"),
    "v_subregions": ("coding_sequence", "imgt_gapped"),
    "coding_interval": ("record_sequence",),
}
PREDICATE_DERIVATIONS = {
    "anchor": ("stated", "imgt_alignment", "unique_motif"),
    "coding_phase": ("stated", "imgt_alignment", "anchor_derived"),
    "phase_independence": ("stated", "anchor_derived"),
    "boundary_5": ("stated", "variant_comparison"),
    "boundary_3": ("stated", "variant_comparison"),
    "functionality": ("stated",),
    "v_subregions": ("imgt_alignment",),
    "coding_interval": ("stated",),
    "extension": ("variant_comparison",),
}
# Claims whose meaning rests on the interpretation profile's established coordinates.
COORDINATE_PREDICATES = ("anchor", "coding_phase", "phase_independence", "v_subregions", "coding_interval")
# Subject kind: "record" or "candidate".
PREDICATE_SUBJECTS = {p: ("record" if p == "functionality" else "candidate") for p in PREDICATES}
# Segments a predicate may describe (no universal D coding phase or anchor).
PREDICATE_SEGMENTS = {"anchor": ("V", "J"), "coding_phase": ("V", "J"), "phase_independence": ("V", "J"),
                      "v_subregions": ("V",)}

MAPPING_KINDS = ("gapped_alignment", "coding_interval")
MAPPING_VALIDATIONS = ("exact", "mismatch")
MAPPING_TRANSFORMS = {"gapped_alignment": "remove_imgt_gaps", "coding_interval": "extract_interval"}
MAPPING_SPACES = {"gapped_alignment": ("coding_sequence", "imgt_gapped"),
                  "coding_interval": ("record_sequence", "coding_sequence")}

SEVERITIES = ("info", "warning", "error")
FINDING_CATEGORIES = ("unavailable_evidence", "contradictory_evidence", "unsupported_interpretation",
                      "invalid_source_syntax", "informational")

_SEQ_RE = re.compile(r"[ACGT]+")
_FASTA_POINTER_RE = re.compile(r"/(0|[1-9][0-9]*)")
_CODE_RE = re.compile(r"[a-z][a-z0-9_]{0,63}")
_FASTA_ROLES = ("gapped_alignment", "ungapped_sequences")
_MAX_COORDINATE = 2 ** 31 - 1


def text_digest(text: str) -> str:
    """``sha256:<hex>`` of UTF-8 text (e.g. a gapped sequence or a record's ``sequence`` field)."""
    return "sha256:" + hashlib.sha256(text.encode("utf-8")).hexdigest()


def _content_id(prefix: str, body: dict, domain: str) -> str:
    return prefix + ":" + digest(body, domain).split(":", 1)[1]


def _segment_key(segment: Optional[str]) -> str:
    return segment if segment in SEGMENTS else "unassigned"


def _check_record_locator(locator, path: str) -> None:
    check_instance(locator, Locator, path)
    if locator.role != "record_export":
        raise bad(path, f"{locator.id} is not in a record export")


def _check_id_field(d: dict, item, path: str, pattern) -> None:
    given = d["id"]
    check_token(given, ptr(path, "id"), pattern, f"an id {pattern.pattern}")
    if given != item.id:
        raise bad(ptr(path, "id"), f"id {given} does not match the content (expected {item.id})", "digest_mismatch")


# ── Records, candidates, identities ─────────────────────────────────


@dataclass(frozen=True)
class SourceRecord:
    """One supplied record occurrence.  ``parsed`` records have a valid label and segment;
    ``malformed``/``unsupported`` ones carry a ``problem`` and keep whatever was readable."""

    locator: Locator
    status: str
    problem: Optional[str]
    label: Optional[str]
    provider_id: Optional[str]
    provider_ref: Optional[str]
    segment: Optional[str]
    taxonomy_id: Optional[int]
    locus: Optional[str]

    KEYS = ("locator", "status", "problem", "label", "provider_id", "provider_ref", "segment", "taxonomy_id",
            "locus")

    def __post_init__(self):
        _check_record_locator(self.locator, "/locator")
        check_enum(self.status, RECORD_STATUSES, "/status")
        parsed = self.status == "parsed"
        if parsed and self.problem is not None:
            raise bad("/problem", "a parsed record has no problem")
        check_text(self.problem, "/problem", nullable=parsed)
        if parsed:
            check_label(self.label, "/label")
        else:
            check_text(self.label, "/label", nullable=True, max_len=MAX_TEXT)
        check_text(self.provider_id, "/provider_id", nullable=True, max_len=MAX_NAME, multiline=False)
        check_text(self.provider_ref, "/provider_ref", nullable=True, max_len=MAX_NAME, multiline=False)
        check_enum(self.segment, SEGMENTS, "/segment", nullable=not parsed)
        check_int(self.taxonomy_id, "/taxonomy_id", nullable=True, minimum=1, maximum=MAX_TAXONOMY_ID)
        check_text(self.locus, "/locus", nullable=True, max_len=MAX_LOCUS, multiline=False)

    @property
    def id(self) -> str:
        return self.locator.id

    def to_dict(self) -> dict:
        d = {k: getattr(self, k) for k in self.KEYS}
        d["locator"] = self.locator.to_dict()
        return d

    @classmethod
    def from_dict(cls, d, path: str = "") -> "SourceRecord":
        require_object(d, cls.KEYS, path)
        fields = {k: d[k] for k in cls.KEYS}
        fields["locator"] = Locator.from_dict(d["locator"], ptr(path, "locator"))
        return construct(cls, path, **fields)


@dataclass(frozen=True)
class SequenceCandidate:
    """A supplied coding sequence (normalized, ``[ACGT]``) at ``locator``, inside ``record``."""

    locator: Locator
    record: Locator
    sequence: str
    sha256: str

    KEYS = ("locator", "record", "sequence", "sha256")

    def __post_init__(self):
        _check_record_locator(self.locator, "/locator")
        _check_record_locator(self.record, "/record")
        if not self.locator.is_under(self.record):
            raise bad("/locator", f"{self.locator.id} is not a field of record {self.record.id}")
        if type(self.sequence) is not str:
            raise bad("/sequence", "expected a string")
        if len(self.sequence) > MAX_SEQUENCE:
            raise bad("/sequence", f"sequence of {len(self.sequence)} bases exceeds {MAX_SEQUENCE}", "resource_limit")
        if not _SEQ_RE.fullmatch(self.sequence):
            raise bad("/sequence", "a normalized candidate sequence matches [ACGT]+")
        check_digest(self.sha256, "/sha256")
        if self.sha256 != sequence_digest(self.sequence):
            raise bad("/sha256", f"digest {self.sha256} does not match the sequence "
                                 f"({sequence_digest(self.sequence)})", "digest_mismatch")

    @property
    def id(self) -> str:
        return self.locator.id

    def to_dict(self) -> dict:
        return {"locator": self.locator.to_dict(), "record": self.record.to_dict(), "sequence": self.sequence,
                "sha256": self.sha256}

    @classmethod
    def from_dict(cls, d, path: str = "") -> "SequenceCandidate":
        require_object(d, cls.KEYS, path)
        return construct(cls, path, locator=Locator.from_dict(d["locator"], ptr(path, "locator")),
                         record=Locator.from_dict(d["record"], ptr(path, "record")), sequence=d["sequence"],
                         sha256=d["sha256"])


@dataclass(frozen=True)
class Mention:
    """One occurrence of a name inside a record's alias or paralog field."""

    record: Locator
    locator: Locator
    relationship: str

    KEYS = ("record", "locator", "relationship")

    def __post_init__(self):
        _check_record_locator(self.record, "/record")
        check_instance(self.locator, Locator, "/locator")
        if not self.locator.is_under(self.record):
            raise bad("/locator", f"{self.locator.id} is not a field of record {self.record.id}")
        check_enum(self.relationship, RELATIONSHIPS, "/relationship")

    def to_dict(self) -> dict:
        return {"record": self.record.to_dict(), "locator": self.locator.to_dict(),
                "relationship": self.relationship}

    @classmethod
    def from_dict(cls, d, path: str = "") -> "Mention":
        require_object(d, cls.KEYS, path)
        return construct(cls, path, record=Locator.from_dict(d["record"], ptr(path, "record")),
                         locator=Locator.from_dict(d["locator"], ptr(path, "locator")),
                         relationship=d["relationship"])


@dataclass(frozen=True)
class MentionedIdentity:
    """A label named in ``asset`` only through relationships: no record, sequence or weight."""

    asset: str
    segment: str
    label: str
    mentions: tuple

    KEYS = ("asset", "segment", "label", "mentions")

    def __post_init__(self):
        check_asset_id(self.asset, "/asset")
        if FORMAT_ROLES[self.asset] != "record_export":
            raise bad("/asset", f"{self.asset} is not a record export")
        check_enum(self.segment, SEGMENTS, "/segment")
        check_label(self.label, "/label")
        check_tuple(self.mentions, "/mentions", MAX_REFS, Mention)
        if not self.mentions:
            raise bad("/mentions", f"{self.label} needs at least one mention")
        check_canonical_order([x.locator.sort_key() for x in self.mentions], "/mentions", "mention")
        for i, mention in enumerate(self.mentions):
            if mention.record.asset != self.asset:
                raise bad(ptr(ptr("/mentions", i), "record"), f"mention is in {mention.record.asset}, not {self.asset}")

    def to_dict(self) -> dict:
        return {"asset": self.asset, "segment": self.segment, "label": self.label,
                "mentions": [x.to_dict() for x in self.mentions]}

    @classmethod
    def from_dict(cls, d, path: str = "") -> "MentionedIdentity":
        require_object(d, cls.KEYS, path)
        return construct(cls, path, asset=d["asset"], segment=d["segment"], label=d["label"],
                         mentions=parse_array(d["mentions"], ptr(path, "mentions"), MAX_REFS, Mention.from_dict))


@dataclass(frozen=True)
class Observation:
    """An original value at a field locator.  ``absent`` (no such field) differs from a
    ``present`` JSON null; an absent observation's value is null.  ``record`` is the owning
    record, or null for set-level fields and FASTA records."""

    locator: Locator
    record: Optional[Locator]
    presence: str
    value: Any
    interpretation: str

    KEYS = ("locator", "record", "presence", "value", "interpretation")

    def __post_init__(self):
        check_instance(self.locator, Locator, "/locator")
        if self.record is not None:
            _check_record_locator(self.record, "/record")
            if not self.locator.is_under(self.record):
                raise bad("/locator", f"{self.locator.id} is not a field of record {self.record.id}")
        if self.locator.role in _FASTA_ROLES:
            if not _FASTA_POINTER_RE.fullmatch(self.locator.pointer):
                raise bad("/locator", f"a FASTA record locator is /<index>, got {self.locator.pointer!r}")
            if self.record is not None:
                raise bad("/record", "a FASTA record observation has no owning JSON record")
        check_enum(self.presence, PRESENCES, "/presence")
        if not isinstance(self.value, JsonValue):
            object.__setattr__(self, "value", JsonValue(self.value, "/value"))
        if self.presence == "absent" and self.value.get() is not None:
            raise bad("/value", "an absent field's value must be null")
        check_enum(self.interpretation, INTERPRETATIONS, "/interpretation")

    @property
    def id(self) -> str:
        return self.locator.id

    def to_dict(self) -> dict:
        return {"locator": self.locator.to_dict(), "record": None if self.record is None else self.record.to_dict(),
                "presence": self.presence, "value": self.value.get(), "interpretation": self.interpretation}

    @classmethod
    def from_dict(cls, d, path: str = "") -> "Observation":
        require_object(d, cls.KEYS, path)
        record = None if d["record"] is None else Locator.from_dict(d["record"], ptr(path, "record"))
        return construct(cls, path, locator=Locator.from_dict(d["locator"], ptr(path, "locator")), record=record,
                         presence=d["presence"], value=JsonValue(d["value"], ptr(path, "value")),
                         interpretation=d["interpretation"])


@dataclass(frozen=True)
class SupportReference:
    """A provider support entry as embedded in the export; never retrieved or verified here."""

    locator: Locator
    record: Locator
    kind: str
    repository: Optional[str]
    reference: Optional[str]
    status: str

    KEYS = ("locator", "record", "kind", "repository", "reference", "status")

    def __post_init__(self):
        check_instance(self.locator, Locator, "/locator")
        _check_record_locator(self.record, "/record")
        if not self.locator.is_under(self.record):
            raise bad("/locator", f"{self.locator.id} is not a field of record {self.record.id}")
        check_enum(self.kind, SUPPORT_KINDS, "/kind")
        check_text(self.repository, "/repository", nullable=True, max_len=MAX_NAME, multiline=False)
        check_text(self.reference, "/reference", nullable=True, max_len=MAX_NAME, multiline=False)
        check_enum(self.status, SUPPORT_STATUSES, "/status")

    def to_dict(self) -> dict:
        return {"locator": self.locator.to_dict(), "record": self.record.to_dict(), "kind": self.kind,
                "repository": self.repository, "reference": self.reference, "status": self.status}

    @classmethod
    def from_dict(cls, d, path: str = "") -> "SupportReference":
        require_object(d, cls.KEYS, path)
        return construct(cls, path, locator=Locator.from_dict(d["locator"], ptr(path, "locator")),
                         record=Locator.from_dict(d["record"], ptr(path, "record")), kind=d["kind"],
                         repository=d["repository"], reference=d["reference"], status=d["status"])


# ── Mappings ────────────────────────────────────────────────────────


def _check_interval(value, path: str) -> tuple:
    if not isinstance(value, tuple) or len(value) != 2:
        raise bad(path, "expected a (start, end) pair")
    start = check_int(value[0], path, minimum=0, maximum=_MAX_COORDINATE)
    end = check_int(value[1], path, minimum=0, maximum=_MAX_COORDINATE)
    if start >= end:
        raise bad(path, f"interval [{start}, {end}) is empty or reversed")
    return value


def _interval_from_json(value, path: str) -> tuple:
    if not isinstance(value, list) or len(value) != 2:
        raise bad(path, "expected [start, end]")
    return tuple(value)


@dataclass(frozen=True)
class Block:
    """An exact interval correspondence: ``input[start, end)`` maps onto ``output[start, end)``."""

    input: tuple
    output: tuple

    KEYS = ("input", "output")

    def __post_init__(self):
        a = _check_interval(self.input, "/input")
        b = _check_interval(self.output, "/output")
        if a[1] - a[0] != b[1] - b[0]:
            raise bad("/output", f"output length {b[1] - b[0]} differs from input length {a[1] - a[0]}")

    def to_dict(self) -> dict:
        return {"input": list(self.input), "output": list(self.output)}

    @classmethod
    def from_dict(cls, d, path: str = "") -> "Block":
        require_object(d, cls.KEYS, path)
        return construct(cls, path, input=_interval_from_json(d["input"], ptr(path, "input")),
                         output=_interval_from_json(d["output"], ptr(path, "output")))


@dataclass(frozen=True)
class Mapping:
    """An exact coordinate mapping with its transform profile and validation result.

    ``gapped_alignment``: a candidate (``coding_sequence``) of an ``airr``/``airr_ex``
    record onto one record of the same variant's gapped FASTA (``imgt_gapped``); blocks
    tile the whole candidate.  ``coding_interval``: a record's ``sequence`` field
    (``record_sequence``) onto its candidate, as one block starting at 0.  ``mismatch``
    mappings carry no blocks.
    """

    source_lock_hash: str
    kind: str
    input: Locator
    output: Locator
    input_sha256: str
    output_sha256: str
    input_space: str
    output_space: str
    blocks: tuple
    transform: str
    validation: str
    id: str = field(init=False, compare=False, repr=False)

    KEYS = ("id", "source_lock_hash", "kind", "input", "output", "input_sha256", "output_sha256", "input_space", "output_space",
            "blocks", "transform", "validation")

    def __post_init__(self):
        check_digest(self.source_lock_hash, "/source_lock_hash")
        check_enum(self.kind, MAPPING_KINDS, "/kind")
        _check_record_locator(self.input, "/input")
        check_instance(self.output, Locator, "/output")
        check_digest(self.input_sha256, "/input_sha256")
        check_digest(self.output_sha256, "/output_sha256")
        spaces = MAPPING_SPACES[self.kind]
        check_enum(self.input_space, (spaces[0],), "/input_space")
        check_enum(self.output_space, (spaces[1],), "/output_space")
        check_enum(self.transform, (MAPPING_TRANSFORMS[self.kind],), "/transform")
        check_enum(self.validation, MAPPING_VALIDATIONS, "/validation")
        check_tuple(self.blocks, "/blocks", MAX_REFS, Block)
        if self.kind == "gapped_alignment":
            expected = GAPPED_VARIANT[self.input.asset]
            if self.output.asset != expected:
                raise bad("/output", f"{self.input.asset} records align only to {expected}, never across variants; "
                                     f"got {self.output.asset}")
            if not _FASTA_POINTER_RE.fullmatch(self.output.pointer):
                raise bad("/output", f"a FASTA record locator is /<index>, got {self.output.pointer!r}")
        else:
            if self.output.role != "record_export" or self.output.asset != self.input.asset:
                raise bad("/output", "a coding interval maps within one record export")
            if self.output == self.input:
                raise bad("/output", "a coding interval maps a record field onto a different candidate field")
        if (self.validation == "exact") != bool(self.blocks):
            raise bad("/blocks", "an exact mapping needs blocks; a mismatch mapping must have none")
        for previous, block in zip(self.blocks, self.blocks[1:]):
            if block.input[0] < previous.input[1] or block.output[0] < previous.output[1]:
                raise bad("/blocks", "blocks must increase without overlap on both sides")
        if self.blocks and self.kind == "gapped_alignment":
            if self.blocks[0].input[0] != 0 or any(b.input[0] != a.input[1]
                                                   for a, b in zip(self.blocks, self.blocks[1:])):
                raise bad("/blocks", "gapped-alignment blocks tile the candidate contiguously from 0")
        if self.blocks and self.kind == "coding_interval":
            if len(self.blocks) != 1 or self.blocks[0].output[0] != 0:
                raise bad("/blocks", "a coding interval is one block whose output starts at 0")
        object.__setattr__(self, "id", _content_id("mapping", self._body(), MAPPING_ID_DOMAIN))

    def _body(self) -> dict:
        return {"source_lock_hash": self.source_lock_hash, "kind": self.kind, "input": self.input.to_dict(),
                "output": self.output.to_dict(),
                "input_sha256": self.input_sha256, "output_sha256": self.output_sha256,
                "input_space": self.input_space, "output_space": self.output_space,
                "blocks": [b.to_dict() for b in self.blocks], "transform": self.transform,
                "validation": self.validation}

    def sort_key(self) -> tuple:
        return (self.input.sort_key(), MAPPING_KINDS.index(self.kind), self.output.sort_key(), self.id)

    def to_dict(self) -> dict:
        d = self._body()
        d["id"] = self.id
        return d

    @classmethod
    def from_dict(cls, d, path: str = "") -> "Mapping":
        require_object(d, cls.KEYS, path)
        item = construct(cls, path, source_lock_hash=d["source_lock_hash"], kind=d["kind"],
                         input=Locator.from_dict(d["input"], ptr(path, "input")),
                         output=Locator.from_dict(d["output"], ptr(path, "output")),
                         input_sha256=d["input_sha256"], output_sha256=d["output_sha256"],
                         input_space=d["input_space"], output_space=d["output_space"],
                         blocks=parse_array(d["blocks"], ptr(path, "blocks"), MAX_REFS, Block.from_dict),
                         transform=d["transform"], validation=d["validation"])
        _check_id_field(d, item, path, MAPPING_ID_RE)
        return item


# ── Claims ──────────────────────────────────────────────────────────


def _check_span(value, path: str) -> tuple:
    if not isinstance(value, list) or len(value) != 2:
        raise bad(path, "expected [start, end]")
    start = check_int(value[0], path, minimum=0, maximum=_MAX_COORDINATE)
    end = check_int(value[1], path, minimum=0, maximum=_MAX_COORDINATE)
    if start >= end:
        raise bad(path, f"span [{start}, {end}) is empty or reversed")
    return start, end


def check_claim_value(predicate: str, value: Any, path: str) -> None:
    """Validate a resolved claim value for a build-model predicate."""
    if predicate == "anchor":
        check_int(value, path, minimum=0, maximum=_MAX_COORDINATE)
    elif predicate == "coding_phase":
        check_int(value, path, minimum=0, maximum=2)
    elif predicate == "phase_independence":
        check_enum(value, PHASE_INDEPENDENCE, path)
    elif predicate in ("boundary_5", "boundary_3"):
        check_enum(value, BOUNDARIES, path)
    elif predicate == "functionality":
        check_enum(value, FUNCTIONALITIES, path)
    elif predicate == "v_subregions":
        if not isinstance(value, dict) or not value:
            raise bad(path, "expected a non-empty object of subregion spans")
        spans = []
        for label, span in value.items():
            if label not in SUBREGION_LABELS:
                raise bad(path, f"unknown subregion {show(label)}; expected {list(SUBREGION_LABELS)}")
            start, end = _check_span(span, ptr(path, label))
            spans.append((start, end, label))
        spans.sort()
        for (_, end, first), (start, _, second) in zip(spans, spans[1:]):
            if start < end:
                raise bad(path, f"subregions {first} and {second} overlap")
    elif predicate == "coding_interval":
        _check_span(value, path)
    elif predicate == "extension":
        require_object(value, ("side", "length", "matches_record_context"), path)
        check_enum(value["side"], EXTENSION_SIDES, ptr(path, "side"))
        check_int(value["length"], ptr(path, "length"), minimum=1, maximum=MAX_SEQUENCE)
        check_bool(value["matches_record_context"], ptr(path, "matches_record_context"))
    else:  # pragma: no cover - guarded by the predicate enum
        raise bad(path, f"unknown predicate {predicate!r}")


@dataclass(frozen=True)
class Claim:
    """A build-model claim about a candidate (or, for ``functionality``, a record).

    ``status``: ``accepted`` (support source/assertion), ``candidate`` (support none),
    ``conflicting`` (support source/assertion; disagrees with a peer) or ``unresolved``
    (support none, value null).  Phase value (``coding_phase``) and whether that phase
    is independent of the landmark (``phase_independence``) are separate claims.
    ``evidence`` lists the observed source fields (evidence roots) in canonical order;
    ``mappings`` the coordinate mappings used.
    """

    source_lock_hash: str
    subject: Locator
    predicate: str
    value: Any
    coordinate_system: Optional[str]
    support: str
    derivation: str
    status: str
    evidence: tuple
    mappings: tuple
    id: str = field(init=False, compare=False, repr=False)

    KEYS = ("id", "source_lock_hash", "subject", "predicate", "value", "coordinate_system", "support", "derivation", "status",
            "evidence", "mappings")

    def __post_init__(self):
        check_digest(self.source_lock_hash, "/source_lock_hash")
        _check_record_locator(self.subject, "/subject")
        check_enum(self.predicate, PREDICATES, "/predicate")
        check_enum(self.status, CLAIM_STATUSES, "/status")
        check_enum(self.support, CLAIM_SUPPORTS, "/support")
        check_enum(self.derivation, DERIVATIONS, "/derivation")
        check_enum(self.coordinate_system, COORDINATE_SYSTEMS, "/coordinate_system", nullable=True)
        systems = PREDICATE_SYSTEMS.get(self.predicate)
        if systems is None and self.coordinate_system is not None:
            raise bad("/coordinate_system", f"{self.predicate} is not positional; coordinate system must be null")
        if systems is not None and self.coordinate_system not in systems:
            raise bad("/coordinate_system", f"{self.predicate} is expressed in one of {list(systems)}")
        if self.derivation not in PREDICATE_DERIVATIONS[self.predicate]:
            raise bad("/derivation", f"{self.predicate} derives only by {list(PREDICATE_DERIVATIONS[self.predicate])}")
        unsupported = self.status in ("candidate", "unresolved")
        if unsupported != (self.support == "none"):
            raise bad("/support", f"a {self.status} claim has support "
                                  + ("'none'" if unsupported else "'source' or 'assertion'"))
        if not isinstance(self.value, JsonValue):
            object.__setattr__(self, "value", JsonValue(self.value, "/value"))
        value = self.value.get()
        if self.status == "unresolved":
            if value is not None:
                raise bad("/value", "an unresolved claim's value must be null")
        else:
            if value is None:
                raise bad("/value", f"a {self.status} claim needs a value")
            check_claim_value(self.predicate, value, "/value")
        locator_tuple(self.evidence, "/evidence", MAX_REFS, canonical=True)
        if not self.evidence:
            raise bad("/evidence", "a claim names its evidence roots")
        check_tuple(self.mappings, "/mappings", MAX_REFS)
        for i, mapping in enumerate(self.mappings):
            check_token(mapping, ptr("/mappings", i), MAPPING_ID_RE, "a mapping id mapping:<64 hex>")
        check_canonical_order(list(self.mappings), "/mappings", "mapping id")
        if self.derivation == "imgt_alignment" and self.status != "unresolved" and not self.mappings:
            raise bad("/mappings", "an alignment-derived claim names the gapped-alignment mapping it uses")
        object.__setattr__(self, "id", _content_id("claim", self._body(), CLAIM_ID_DOMAIN))

    def _body(self) -> dict:
        return {"source_lock_hash": self.source_lock_hash, "subject": self.subject.to_dict(),
                "predicate": self.predicate, "value": self.value.get(),
                "coordinate_system": self.coordinate_system, "support": self.support,
                "derivation": self.derivation, "status": self.status,
                "evidence": [e.to_dict() for e in self.evidence], "mappings": list(self.mappings)}

    def sort_key(self) -> tuple:
        return (self.subject.sort_key(), PREDICATES.index(self.predicate), self.id)

    def to_dict(self) -> dict:
        d = self._body()
        d["id"] = self.id
        return d

    @classmethod
    def from_dict(cls, d, path: str = "") -> "Claim":
        require_object(d, cls.KEYS, path)
        item = construct(
            cls, path, source_lock_hash=d["source_lock_hash"],
            subject=Locator.from_dict(d["subject"], ptr(path, "subject")), predicate=d["predicate"],
            value=JsonValue(d["value"], ptr(path, "value")), coordinate_system=d["coordinate_system"],
            support=d["support"], derivation=d["derivation"], status=d["status"],
            evidence=parse_array(d["evidence"], ptr(path, "evidence"), MAX_REFS, Locator.from_dict),
            mappings=parse_array(d["mappings"], ptr(path, "mappings"), MAX_REFS, lambda v, p: v))
        _check_id_field(d, item, path, CLAIM_ID_RE)
        return item


# ── Findings ────────────────────────────────────────────────────────


@dataclass(frozen=True)
class Finding:
    """A stable-coded finding: ``subjects`` are records/candidates, ``locators`` the source
    fields involved, ``observations`` a JSON object of structured values, ``capabilities``
    the runtime capabilities it affects (canonical capability order)."""

    source_lock_hash: str
    code: str
    severity: str
    category: str
    subjects: tuple
    locators: tuple
    observations: Any
    capabilities: tuple
    id: str = field(init=False, compare=False, repr=False)

    KEYS = ("id", "source_lock_hash", "code", "severity", "category", "subjects", "locators", "observations", "capabilities")

    def __post_init__(self):
        check_digest(self.source_lock_hash, "/source_lock_hash")
        check_token(self.code, "/code", _CODE_RE, "a code matching [a-z][a-z0-9_]*")
        check_enum(self.severity, SEVERITIES, "/severity")
        check_enum(self.category, FINDING_CATEGORIES, "/category")
        locator_tuple(self.subjects, "/subjects", MAX_REFS, canonical=True)
        locator_tuple(self.locators, "/locators", MAX_REFS, canonical=True)
        if not self.subjects and not self.locators:
            raise bad("/subjects", "a finding names subjects or source locators")
        if not isinstance(self.observations, JsonValue):
            object.__setattr__(self, "observations", JsonValue(self.observations, "/observations"))
        if not isinstance(self.observations.get(), dict):
            raise bad("/observations", "expected an object of structured observations")
        check_tuple(self.capabilities, "/capabilities", len(CAPABILITIES))
        for i, name in enumerate(self.capabilities):
            check_enum(name, CAPABILITIES, ptr("/capabilities", i))
        check_canonical_order([CAPABILITIES.index(c) for c in self.capabilities], "/capabilities", "capability")
        object.__setattr__(self, "id", _content_id("finding", self._body(), FINDING_ID_DOMAIN))

    def _body(self) -> dict:
        return {"source_lock_hash": self.source_lock_hash, "code": self.code, "severity": self.severity,
                "category": self.category,
                "subjects": [s.to_dict() for s in self.subjects], "locators": [s.to_dict() for s in self.locators],
                "observations": self.observations.get(), "capabilities": list(self.capabilities)}

    def sort_key(self) -> tuple:
        return (self.code, self.id)

    def to_dict(self) -> dict:
        d = self._body()
        d["id"] = self.id
        return d

    @classmethod
    def from_dict(cls, d, path: str = "") -> "Finding":
        require_object(d, cls.KEYS, path)
        item = construct(
            cls, path, source_lock_hash=d["source_lock_hash"], code=d["code"], severity=d["severity"],
            category=d["category"],
            subjects=parse_array(d["subjects"], ptr(path, "subjects"), MAX_REFS, Locator.from_dict),
            locators=parse_array(d["locators"], ptr(path, "locators"), MAX_REFS, Locator.from_dict),
            observations=JsonValue(d["observations"], ptr(path, "observations")),
            capabilities=parse_array(d["capabilities"], ptr(path, "capabilities"), len(CAPABILITIES),
                                     lambda v, p: v))
        _check_id_field(d, item, path, FINDING_ID_RE)
        return item


# ── The inventory ───────────────────────────────────────────────────

_LISTS = (  # field, class, limit, canonical key
    ("records", SourceRecord, MAX_RECORDS, lambda x: x.locator.sort_key()),
    ("candidates", SequenceCandidate, MAX_CANDIDATES, lambda x: x.locator.sort_key()),
    ("mentioned_identities", MentionedIdentity, MAX_MENTIONED, lambda x: (FORMATS.index(x.asset), x.label)),
    ("observations", Observation, MAX_OBSERVATIONS, lambda x: x.locator.sort_key()),
    ("support", SupportReference, MAX_SUPPORT, lambda x: x.locator.sort_key()),
    ("claims", Claim, MAX_CLAIMS, lambda x: x.sort_key()),
    ("mappings", Mapping, MAX_MAPPINGS, lambda x: x.sort_key()),
    ("findings", Finding, MAX_FINDINGS, lambda x: x.sort_key()),
)
_ITEM_NAMES = {"records": "record", "candidates": "candidate", "mentioned_identities": "mentioned identity",
               "observations": "observation", "support": "support reference", "claims": "claim",
               "mappings": "mapping", "findings": "finding"}
_KEYS = ("artifact", "schema_version", "source_lock_hash", "adapter") + tuple(name for name, *_ in _LISTS)


@dataclass(frozen=True)
class SourceInventory:
    """The complete normalized inventory of one lock, every list in canonical order."""

    source_lock_hash: str
    adapter: AdapterRef
    records: tuple
    candidates: tuple
    mentioned_identities: tuple
    observations: tuple
    support: tuple
    claims: tuple
    mappings: tuple
    findings: tuple

    def __post_init__(self):
        check_digest(self.source_lock_hash, "/source_lock_hash")
        check_instance(self.adapter, AdapterRef, "/adapter")
        for name, cls, limit, key in _LISTS:
            items = check_tuple(getattr(self, name), "/" + name, limit, cls)
            check_canonical_order([key(x) for x in items], "/" + name, _ITEM_NAMES[name])
        for name in ("claims", "mappings", "findings"):
            for i, item in enumerate(getattr(self, name)):
                if item.source_lock_hash != self.source_lock_hash:
                    raise bad(ptr(ptr("/" + name, i), "source_lock_hash"),
                              f"{item.id} was derived under lock {item.source_lock_hash}, not this inventory's "
                              f"{self.source_lock_hash}", "digest_mismatch")
        records = {r.id: r for r in self.records}
        candidates = {c.id: c for c in self.candidates}
        observations = {o.id: o for o in self.observations}
        mappings = {x.id: x for x in self.mappings}
        index = {"records": records, "candidates": candidates, "observations": observations, "mappings": mappings,
                 "claims": {c.id: c for c in self.claims}, "findings": {f.id: f for f in self.findings}}
        object.__setattr__(self, "_index", index)
        self._check_items(records, candidates, observations)
        self._check_mappings(candidates, observations)
        self._check_claims(records, candidates, observations, mappings)
        for i, finding in enumerate(self.findings):
            for j, subject in enumerate(finding.subjects):
                if subject.id not in records and subject.id not in candidates:
                    raise bad(ptr(ptr(ptr("/findings", i), "subjects"), j), f"unknown record or candidate {subject.id}")
            for j, locator in enumerate(finding.locators):
                if locator.id not in observations and locator.id not in records and locator.id not in candidates:
                    raise bad(ptr(ptr(ptr("/findings", i), "locators"), j),
                              f"{locator.id} is not an observation, record or candidate")

    # -- checks -------------------------------------------------------

    def _check_items(self, records, candidates, observations) -> None:
        seen = {}
        for i, candidate in enumerate(self.candidates):
            path = ptr("/candidates", i)
            if candidate.record.id not in records:
                raise bad(ptr(path, "record"), f"unknown record {candidate.record.id}")
            key = (candidate.record.id, candidate.sha256)
            if key in seen:
                raise bad("/candidates", f"record {candidate.record.id} has two candidates with digest "
                                         f"{candidate.sha256} (indexes {seen[key]} and {i})")
            seen[key] = i
        for name in ("observations", "support"):
            for i, item in enumerate(getattr(self, name)):
                if item.record is not None and item.record.id not in records:
                    raise bad(ptr(ptr("/" + name, i), "record"), f"unknown record {item.record.id}")
        labels = {(r.locator.asset, r.label) for r in self.records}
        for i, identity in enumerate(self.mentioned_identities):
            path = ptr("/mentioned_identities", i)
            if (identity.asset, identity.label) in labels:
                raise bad(path, f"{identity.label} has its own record in {identity.asset}; it is not a "
                                "mentioned-only identity")
            for j, mention in enumerate(identity.mentions):
                record = records.get(mention.record.id)
                where = ptr(ptr(path, "mentions"), j)
                if record is None:
                    raise bad(ptr(where, "record"), f"unknown record {mention.record.id}")
                if record.segment != identity.segment:
                    raise bad(path, f"{identity.label} is mentioned by {record.segment} record {record.label}, "
                                    f"not as segment {identity.segment}")

    def _check_mappings(self, candidates, observations) -> None:
        inputs, outputs = set(), set()
        for i, mapping in enumerate(self.mappings):
            path = ptr("/mappings", i)
            if mapping.kind == "gapped_alignment":
                candidate, own, other = candidates.get(mapping.input.id), "input", mapping.output
                if candidate is None:
                    raise bad(ptr(path, "input"), f"{mapping.input.id} is not a sequence candidate")
                if other.id not in observations:
                    raise bad(ptr(path, "output"), f"gapped record {other.id} is not an observation")
                sha, end = mapping.input_sha256, (mapping.blocks[-1].input[1] if mapping.blocks else None)
                if other.id in outputs:
                    raise bad(ptr(path, "output"), f"gapped record {other.id} is mapped more than once")
                outputs.add(other.id)
            else:
                candidate, own = candidates.get(mapping.output.id), "output"
                if candidate is None:
                    raise bad(ptr(path, "output"), f"{mapping.output.id} is not a sequence candidate")
                if mapping.input.id not in observations or not mapping.input.is_under(candidate.record):
                    raise bad(ptr(path, "input"), f"{mapping.input.id} is not an observed field of "
                                                  f"record {candidate.record.id}")
                sha, end = mapping.output_sha256, (mapping.blocks[-1].output[1] if mapping.blocks else None)
            if sha != candidate.sha256:
                raise bad(ptr(path, own + "_sha256"), f"digest {sha} differs from candidate {candidate.id} "
                                                      f"({candidate.sha256})", "digest_mismatch")
            if end is not None and end != len(candidate.sequence):
                raise bad(ptr(path, "blocks"), f"blocks cover {end} bases of candidate {candidate.id} "
                                               f"({len(candidate.sequence)} bases)")
            key = (mapping.kind, candidate.id)
            if key in inputs:
                raise bad(path, f"candidate {candidate.id} has more than one {mapping.kind} mapping")
            inputs.add(key)

    def _check_claims(self, records, candidates, observations, mappings) -> None:
        groups: dict = {}
        for i, claim in enumerate(self.claims):
            path = ptr("/claims", i)
            kind = PREDICATE_SUBJECTS[claim.predicate]
            pool = records if kind == "record" else candidates
            subject = pool.get(claim.subject.id)
            if subject is None:
                raise bad(ptr(path, "subject"), f"{claim.predicate} describes a {kind}; {claim.subject.id} is not one")
            record = subject if kind == "record" else records[subject.record.id]
            allowed = PREDICATE_SEGMENTS.get(claim.predicate)
            if allowed is not None and record.segment not in allowed:
                raise bad(ptr(path, "predicate"), f"{claim.predicate} does not apply to a {record.segment} "
                                                  f"segment ({record.label})")
            for j, locator in enumerate(claim.evidence):
                if locator.id not in observations:
                    raise bad(ptr(ptr(path, "evidence"), j), f"evidence root {locator.id} is not an observation")
            for j, mapping in enumerate(claim.mappings):
                if mapping not in mappings:
                    raise bad(ptr(ptr(path, "mappings"), j), f"unknown mapping {mapping}")
            value = claim.value.get()
            if value is not None and kind == "candidate" and claim.coordinate_system == "coding_sequence":
                length = len(subject.sequence)
                if claim.predicate == "anchor" and value > length - 3:
                    raise bad(ptr(path, "value"), f"anchor {value} leaves no codon in {length} bases of {subject.id}")
                if claim.predicate == "v_subregions" and max(span[1] for span in value.values()) > length:
                    raise bad(ptr(path, "value"), f"subregions end beyond {length} bases of {subject.id}")
            groups.setdefault((claim.subject.id, claim.predicate, claim.coordinate_system), []).append((i, claim))
        for (subject, predicate, _), members in groups.items():
            accepted = [(i, c) for i, c in members if c.status == "accepted"]
            conflicting = [(i, c) for i, c in members if c.status == "conflicting"]
            values = {c.value for _, c in accepted}
            if len(values) > 1:
                raise bad(ptr(ptr("/claims", accepted[-1][0]), "value"),
                          f"accepted {predicate} claims about {subject} disagree; disagreeing claims are conflicting")
            if conflicting and accepted:
                raise bad(ptr(ptr("/claims", accepted[0][0]), "status"),
                          f"{predicate} of {subject} is in conflict; no claim of it is accepted")
            if conflicting and len({c.value for _, c in conflicting}) < 2:
                raise bad(ptr(ptr("/claims", conflicting[0][0]), "status"),
                          f"a conflicting {predicate} claim about {subject} needs a disagreeing peer")
        derived = {c.subject.id for c in self.claims
                   if (c.predicate, c.derivation, c.status) == ("coding_phase", "anchor_derived", "accepted")}
        for i, claim in enumerate(self.claims):
            if (claim.predicate == "phase_independence" and claim.status == "accepted"
                    and claim.subject.id in derived and claim.value.get() != "anchor_derived"):
                raise bad(ptr(ptr("/claims", i), "value"),
                          f"the accepted coding phase of {claim.subject.id} is anchor-derived, not independent")

    # -- construction and views ---------------------------------------

    @classmethod
    def assemble(cls, **fields) -> "SourceInventory":
        """Build an inventory, putting every list into canonical order first."""
        for name, _, _, key in _LISTS:
            fields[name] = tuple(sorted(fields[name], key=key))
        return cls(**fields)

    def record(self, locator: Locator) -> SourceRecord:
        return self._index["records"][locator.id]

    def candidate(self, locator: Locator) -> SequenceCandidate:
        return self._index["candidates"][locator.id]

    def observation(self, locator: Locator) -> Observation:
        return self._index["observations"][locator.id]

    def claim(self, claim_id: str) -> Claim:
        return self._index["claims"][claim_id]

    def finding(self, finding_id: str) -> Finding:
        return self._index["findings"][finding_id]

    def mapping(self, mapping_id: str) -> Mapping:
        return self._index["mappings"][mapping_id]

    def candidates_of(self, record: Locator) -> tuple:
        return tuple(c for c in self.candidates if c.record == record)

    def export_counts(self, asset_id: str) -> dict:
        """Record occurrences, sequence candidates and mentioned-only identities of one
        record export, by segment (``unassigned`` for records without a valid segment)."""
        def empty():
            return {"V": 0, "D": 0, "J": 0, "unassigned": 0}
        counts = {"record_occurrences": empty(), "sequence_candidates": empty(), "mentioned_identities": empty()}
        records = self._index["records"]
        for record in self.records:
            if record.locator.asset == asset_id:
                counts["record_occurrences"][_segment_key(record.segment)] += 1
        for candidate in self.candidates:
            if candidate.locator.asset == asset_id:
                counts["sequence_candidates"][_segment_key(records[candidate.record.id].segment)] += 1
        for identity in self.mentioned_identities:
            if identity.asset == asset_id:
                counts["mentioned_identities"][identity.segment] += 1
        return counts

    def asset_references(self):
        """Every ``(path, asset id)`` the inventory refers to, for cross-checks against a lock."""
        for name, *_ in _LISTS:
            for i, item in enumerate(getattr(self, name)):
                base = ptr("/" + name, i)
                if isinstance(item, MentionedIdentity):
                    yield ptr(base, "asset"), item.asset
                for key, value in vars(item).items():
                    if isinstance(value, Locator):
                        yield ptr(ptr(base, key), "asset"), value.asset
                    elif isinstance(value, tuple):
                        for j, entry in enumerate(value):
                            if isinstance(entry, Locator):
                                yield ptr(ptr(ptr(base, key), j), "asset"), entry.asset
                            elif isinstance(entry, Mention):
                                for part in ("record", "locator"):
                                    yield ptr(ptr(ptr(ptr(base, key), j), part), "asset"), getattr(entry, part).asset

    def to_dict(self) -> dict:
        d = {"artifact": ARTIFACT, "schema_version": SCHEMA_VERSION, "source_lock_hash": self.source_lock_hash,
             "adapter": self.adapter.to_dict()}
        for name, *_ in _LISTS:
            d[name] = [x.to_dict() for x in getattr(self, name)]
        return d

    @classmethod
    def from_dict(cls, d, path: str = "") -> "SourceInventory":
        check_header(d, ARTIFACT, SCHEMA_VERSION, path)
        require_object(d, _KEYS, path)
        fields = {name: parse_array(d[name], ptr(path, name), limit, kind.from_dict)
                  for name, kind, limit, _ in _LISTS}
        return construct(cls, path, source_lock_hash=d["source_lock_hash"],
                         adapter=AdapterRef.from_dict(d["adapter"], ptr(path, "adapter")), **fields)

    @classmethod
    def from_json(cls, text: str) -> "SourceInventory":
        return cls.from_dict(parse_json(text, max_bytes=MAX_INVENTORY_BYTES))

    @classmethod
    def read(cls, path) -> "SourceInventory":
        return cls.from_dict(read_json(path, max_bytes=MAX_INVENTORY_BYTES))

