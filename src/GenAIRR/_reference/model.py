"""Closed data model for a portable reference catalogue.

Every class is a frozen dataclass with ``to_dict``/``from_dict``; ``from_dict``
rejects unknown or missing keys.  ``Catalogue.from_dict`` additionally enforces
referential and semantic integrity.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import re
import unicodedata
from typing import Any, Optional

from GenAIRR._artifact_io import typed_digest
from GenAIRR._reference.errors import ReferenceArtifactError

DIGEST_DOMAIN = "genairr_reference_catalogue_v1"

SEGMENTS = ("V", "D", "J")
CHAIN_TYPES = ("vj", "vdj")
V_ALIGNMENTS = ("none", "imgt_gapped")
PREDICATES = ("anchor", "coding_phase", "boundary_5", "boundary_3", "functionality", "v_subregions")
SUPPORTS = ("source", "assertion", "none")
DERIVATIONS = ("stated", "imgt_alignment", "unique_motif", "anchor_derived")
CLAIM_STATUSES = ("accepted", "candidate", "conflict")
FUNCTIONALITIES = ("functional", "orf", "pseudogene", "unknown")
BOUNDARIES = ("complete", "partial")
SUBREGION_LABELS = ("FWR1", "CDR1", "FWR2", "CDR2", "FWR3")
GENE_BASES = ("declared", "label_prefix", "label")
RECORD_STATUSES = ("accepted", "duplicate", "rejected")
RELATIONSHIP_KINDS = ("alias", "paralog")

_SEQ_RE = re.compile(r"[ACGT]{1,65535}")
_DIGEST_RE = re.compile(r"sha256:[0-9a-f]{64}")
_REC_ID_RE = re.compile(r"rec:([VDJ]):(\d+)")


def sequence_digest(sequence: str) -> str:
    return "sha256:" + hashlib.sha256(sequence.encode("ascii")).hexdigest()


def check_label(text: Any) -> Optional[str]:
    """Return ``None`` if ``text`` is a valid label, else a diagnostic."""
    if not isinstance(text, str):
        return "must be a string"
    if not 1 <= len(text) <= 256:
        return f"must be 1-256 characters, got {len(text)}"
    for position, char in enumerate(text):
        if char == ",":
            return f"contains ',' at position {position}"
        if char.isspace():
            return f"contains whitespace U+{ord(char):04X} at position {position}"
        category = unicodedata.category(char)
        if category == "Cc":
            return f"contains control character U+{ord(char):04X} at position {position}"
        if category == "Cf":
            return f"contains format character U+{ord(char):04X} at position {position}"
    return None


def gene_from_label(label: str) -> tuple:
    """``(gene, basis)`` derived from a label: the text before ``*`` if valid, else the label."""
    if "*" in label:
        prefix = label.split("*", 1)[0]
        if prefix and check_label(prefix) is None:
            return prefix, "label_prefix"
    return label, "label"


def anchor_phase_disagreements(claims) -> list:
    """``(index, message)`` for each accepted anchor-derived coding phase that differs from
    its allele's accepted anchor modulo 3 (only when one agreed, unconflicted anchor exists)."""
    anchors: dict = {}
    for c in claims:
        if c.predicate == "anchor":
            anchors.setdefault(c.subject, []).append(c)
    found = []
    for index, c in enumerate(claims):
        if (c.predicate, c.derivation, c.status) != ("coding_phase", "anchor_derived", "accepted"):
            continue
        group = anchors.get(c.subject, [])
        values = {a.value for a in group if a.status == "accepted"}
        if any(a.status == "conflict" for a in group) or len(values) != 1:
            continue
        (anchor,) = values
        if c.value != anchor % 3:
            found.append((index, f"{c.subject}: anchor-derived coding phase {c.value} disagrees with "
                                 f"accepted anchor {anchor} ({anchor} % 3 = {anchor % 3})"))
    return found


def _bad(path: str, message: str, kind: str = "malformed"):
    return ReferenceArtifactError(kind, path, message)


def _obj(value, keys, path):
    if not isinstance(value, dict):
        raise _bad(path, "expected an object")
    missing = [k for k in keys if k not in value]
    unknown = [k for k in value if k not in keys]
    if missing:
        raise _bad(path, f"missing keys {missing}")
    if unknown:
        raise _bad(path, f"unknown keys {unknown}")
    return value


def _str(value, path, *, optional=False, nonblank=False):
    if value is None and optional:
        return None
    if not isinstance(value, str):
        raise _bad(path, "expected a string" + (" or null" if optional else ""))
    if nonblank and not value.strip():
        raise _bad(path, "must not be blank")
    return value


def _int(value, path, *, optional=False, minimum=None):
    if value is None and optional:
        return None
    if type(value) is not int:
        raise _bad(path, "expected an integer" + (" or null" if optional else ""))
    if minimum is not None and value < minimum:
        raise _bad(path, f"must be >= {minimum}")
    return value


def _bool(value, path):
    if not isinstance(value, bool):
        raise _bad(path, "expected a Boolean")
    return value


def _enum(value, allowed, path):
    if not isinstance(value, str) or value not in allowed:
        raise _bad(path, f"must be one of {list(allowed)}")
    return value


def _list(value, path):
    if not isinstance(value, list):
        raise _bad(path, "expected an array")
    return value


def _str_tuple(value, path):
    return tuple(_str(v, f"{path}[{i}]") for i, v in enumerate(_list(value, path)))


def check_claim_value(predicate: str, value: Any, length: Optional[int], path: str) -> None:
    """Validate a claim/suggestion value for ``predicate``.

    ``length`` is the ungapped allele length; when ``None`` bounds are not checked.
    """
    if predicate == "anchor":
        _int(value, path, minimum=0)
        if length is not None and value > length - 3:
            raise _bad(path, f"anchor {value} out of bounds for length {length}", "integrity")
    elif predicate == "coding_phase":
        _int(value, path)
        if value not in (0, 1, 2):
            raise _bad(path, "coding phase must be 0, 1 or 2")
    elif predicate in ("boundary_5", "boundary_3"):
        _enum(value, BOUNDARIES, path)
    elif predicate == "functionality":
        _enum(value, FUNCTIONALITIES, path)
    elif predicate == "v_subregions":
        if not isinstance(value, dict) or not value:
            raise _bad(path, "expected a non-empty object of subregions")
        spans = []
        for label, span in value.items():
            if label not in SUBREGION_LABELS:
                raise _bad(path, f"unknown subregion {label!r}")
            if not isinstance(span, (list, tuple)) or len(span) != 2:
                raise _bad(f"{path}.{label}", "expected [start, end]")
            start = _int(span[0], f"{path}.{label}[0]", minimum=0)
            end = _int(span[1], f"{path}.{label}[1]", minimum=0)
            if start >= end:
                raise _bad(f"{path}.{label}", "start must be less than end", "integrity")
            if length is not None and end > length:
                raise _bad(f"{path}.{label}", f"end {end} out of bounds for length {length}", "integrity")
            spans.append((start, end, label))
        spans.sort()
        for (_, e1, l1), (s2, _, l2) in zip(spans, spans[1:]):
            if s2 < e1:
                raise _bad(path, f"subregions {l1} and {l2} overlap", "integrity")
    else:
        raise _bad(path, f"unknown predicate {predicate!r}")


def _tuple_dicts(items):
    return [i.to_dict() for i in items]


@dataclass(frozen=True)
class Taxon:
    scientific_name: str
    taxonomy_id: Optional[int]

    def to_dict(self):
        return {"scientific_name": self.scientific_name, "taxonomy_id": self.taxonomy_id}

    @classmethod
    def from_dict(cls, d, path="taxon"):
        _obj(d, ("scientific_name", "taxonomy_id"), path)
        return cls(_str(d["scientific_name"], f"{path}.scientific_name", nonblank=True),
                   _int(d["taxonomy_id"], f"{path}.taxonomy_id", optional=True, minimum=0))


@dataclass(frozen=True)
class ReferenceHeader:
    name: str
    version: Optional[str]
    taxon: Taxon
    locus: str
    chain_type: str
    build_date: Optional[str]

    def to_dict(self):
        return {"name": self.name, "version": self.version, "taxon": self.taxon.to_dict(),
                "locus": self.locus, "chain_type": self.chain_type, "build_date": self.build_date}

    @classmethod
    def from_dict(cls, d, path="header"):
        _obj(d, ("name", "version", "taxon", "locus", "chain_type", "build_date"), path)
        build_date = _str(d["build_date"], f"{path}.build_date", optional=True)
        if build_date is not None:
            import datetime
            try:
                if not re.fullmatch(r"\d{4}-\d{2}-\d{2}", build_date):
                    raise ValueError
                datetime.date.fromisoformat(build_date)
            except ValueError:
                raise _bad(f"{path}.build_date", "expected YYYY-MM-DD") from None
        return cls(_str(d["name"], f"{path}.name", nonblank=True),
                   _str(d["version"], f"{path}.version", optional=True),
                   Taxon.from_dict(d["taxon"], f"{path}.taxon"),
                   _str(d["locus"], f"{path}.locus", nonblank=True),
                   _enum(d["chain_type"], CHAIN_TYPES, f"{path}.chain_type"),
                   build_date)


@dataclass(frozen=True)
class SourceFile:
    id: str
    role: str
    file_name: str
    sha256: str
    bytes: int

    def to_dict(self):
        return {"id": self.id, "role": self.role, "file_name": self.file_name,
                "sha256": self.sha256, "bytes": self.bytes}

    @classmethod
    def from_dict(cls, d, path="source"):
        _obj(d, ("id", "role", "file_name", "sha256", "bytes"), path)
        sha = _str(d["sha256"], f"{path}.sha256")
        if not _DIGEST_RE.fullmatch(sha):
            raise _bad(f"{path}.sha256", "expected sha256:<64 hex>")
        return cls(_str(d["id"], f"{path}.id", nonblank=True), _str(d["role"], f"{path}.role"),
                   _str(d["file_name"], f"{path}.file_name"), sha, _int(d["bytes"], f"{path}.bytes", minimum=0))


@dataclass(frozen=True)
class Transform:
    kind: str
    input_sha256: str
    output_sha256: str

    def to_dict(self):
        return {"kind": self.kind, "input_sha256": self.input_sha256, "output_sha256": self.output_sha256}

    @classmethod
    def from_dict(cls, d, path="transform"):
        _obj(d, ("kind", "input_sha256", "output_sha256"), path)
        for key in ("input_sha256", "output_sha256"):
            if not _DIGEST_RE.fullmatch(_str(d[key], f"{path}.{key}")):
                raise _bad(f"{path}.{key}", "expected sha256:<64 hex>")
        return cls(_str(d["kind"], f"{path}.kind", nonblank=True), d["input_sha256"], d["output_sha256"])


@dataclass(frozen=True)
class RecordEntry:
    id: str
    source: str
    index: int
    header: str
    label: Optional[str]
    status: str
    reason: Optional[str]
    allele: Optional[str]
    duplicate_of: Optional[str]

    def to_dict(self):
        return {"id": self.id, "source": self.source, "index": self.index, "header": self.header,
                "label": self.label, "status": self.status, "reason": self.reason,
                "allele": self.allele, "duplicate_of": self.duplicate_of}

    @classmethod
    def from_dict(cls, d, path="record"):
        keys = ("id", "source", "index", "header", "label", "status", "reason", "allele", "duplicate_of")
        _obj(d, keys, path)
        return cls(_str(d["id"], f"{path}.id", nonblank=True), _str(d["source"], f"{path}.source"),
                   _int(d["index"], f"{path}.index", minimum=0), _str(d["header"], f"{path}.header"),
                   _str(d["label"], f"{path}.label", optional=True),
                   _enum(d["status"], RECORD_STATUSES, f"{path}.status"),
                   _str(d["reason"], f"{path}.reason", optional=True),
                   _str(d["allele"], f"{path}.allele", optional=True),
                   _str(d["duplicate_of"], f"{path}.duplicate_of", optional=True))


@dataclass(frozen=True)
class Allele:
    id: str
    segment: str
    label: str
    gene: str
    gene_basis: str
    sequence: str
    sequence_sha256: str
    gapped_sequence: Optional[str]
    transforms: tuple
    records: tuple

    def to_dict(self):
        return {"id": self.id, "segment": self.segment, "label": self.label, "gene": self.gene,
                "gene_basis": self.gene_basis, "sequence": self.sequence,
                "sequence_sha256": self.sequence_sha256, "gapped_sequence": self.gapped_sequence,
                "transforms": _tuple_dicts(self.transforms), "records": list(self.records)}

    @classmethod
    def from_dict(cls, d, path="allele"):
        keys = ("id", "segment", "label", "gene", "gene_basis", "sequence", "sequence_sha256",
                "gapped_sequence", "transforms", "records")
        _obj(d, keys, path)
        transforms = tuple(Transform.from_dict(t, f"{path}.transforms[{i}]")
                           for i, t in enumerate(_list(d["transforms"], f"{path}.transforms")))
        return cls(_str(d["id"], f"{path}.id"), _enum(d["segment"], SEGMENTS, f"{path}.segment"),
                   _str(d["label"], f"{path}.label"), _str(d["gene"], f"{path}.gene"),
                   _enum(d["gene_basis"], GENE_BASES, f"{path}.gene_basis"),
                   _str(d["sequence"], f"{path}.sequence"), _str(d["sequence_sha256"], f"{path}.sequence_sha256"),
                   _str(d["gapped_sequence"], f"{path}.gapped_sequence", optional=True),
                   transforms, _str_tuple(d["records"], f"{path}.records"))


@dataclass(frozen=True)
class Claim:
    subject: str
    predicate: str
    value: Any
    support: str
    derivation: str
    status: str
    source: Optional[str]
    method: Optional[str]

    def to_dict(self):
        return {"subject": self.subject, "predicate": self.predicate, "value": self.value,
                "support": self.support, "derivation": self.derivation, "status": self.status,
                "source": self.source, "method": self.method}

    @classmethod
    def from_dict(cls, d, path="claim"):
        keys = ("subject", "predicate", "value", "support", "derivation", "status", "source", "method")
        _obj(d, keys, path)
        return cls(_str(d["subject"], f"{path}.subject"), _enum(d["predicate"], PREDICATES, f"{path}.predicate"),
                   d["value"], _enum(d["support"], SUPPORTS, f"{path}.support"),
                   _enum(d["derivation"], DERIVATIONS, f"{path}.derivation"),
                   _enum(d["status"], CLAIM_STATUSES, f"{path}.status"),
                   _str(d["source"], f"{path}.source", optional=True),
                   _str(d["method"], f"{path}.method", optional=True))


@dataclass(frozen=True)
class Entity:
    id: str
    name: str
    executable: bool

    def to_dict(self):
        return {"id": self.id, "name": self.name, "executable": self.executable}

    @classmethod
    def from_dict(cls, d, path="entity"):
        _obj(d, ("id", "name", "executable"), path)
        return cls(_str(d["id"], f"{path}.id", nonblank=True), _str(d["name"], f"{path}.name", nonblank=True),
                   _bool(d["executable"], f"{path}.executable"))


@dataclass(frozen=True)
class Relationship:
    kind: str
    allele: str
    entity: str
    source: Optional[str]

    def to_dict(self):
        return {"kind": self.kind, "allele": self.allele, "entity": self.entity, "source": self.source}

    @classmethod
    def from_dict(cls, d, path="relationship"):
        _obj(d, ("kind", "allele", "entity", "source"), path)
        return cls(_enum(d["kind"], RELATIONSHIP_KINDS, f"{path}.kind"), _str(d["allele"], f"{path}.allele"),
                   _str(d["entity"], f"{path}.entity"), _str(d["source"], f"{path}.source", optional=True))


@dataclass(frozen=True)
class Finding:
    kind: str
    subjects: tuple
    detail: str

    def to_dict(self):
        return {"kind": self.kind, "subjects": list(self.subjects), "detail": self.detail}

    @classmethod
    def from_dict(cls, d, path="finding"):
        _obj(d, ("kind", "subjects", "detail"), path)
        return cls(_str(d["kind"], f"{path}.kind", nonblank=True), _str_tuple(d["subjects"], f"{path}.subjects"),
                   _str(d["detail"], f"{path}.detail"))


def _unique(items, path):
    seen = set()
    for item in items:
        if item.id in seen:
            raise _bad(path, f"duplicate id {item.id!r}", "integrity")
        seen.add(item.id)
    return seen


@dataclass(frozen=True)
class Catalogue:
    header: ReferenceHeader
    conventions: dict
    sources: tuple
    records: tuple
    alleles: tuple
    claims: tuple
    entities: tuple
    relationships: tuple
    findings: tuple

    def to_dict(self):
        return {"header": self.header.to_dict(), "conventions": dict(self.conventions),
                "sources": _tuple_dicts(self.sources), "records": _tuple_dicts(self.records),
                "alleles": _tuple_dicts(self.alleles), "claims": _tuple_dicts(self.claims),
                "entities": _tuple_dicts(self.entities), "relationships": _tuple_dicts(self.relationships),
                "findings": _tuple_dicts(self.findings)}

    def digest(self) -> str:
        return typed_digest(self.to_dict(), DIGEST_DOMAIN)

    def allele(self, allele_id: str) -> Allele:
        for allele in self.alleles:
            if allele.id == allele_id:
                return allele
        raise KeyError(allele_id)

    def claims_for(self, allele_id: str, predicate: str) -> tuple:
        return tuple(c for c in self.claims if c.subject == allele_id and c.predicate == predicate)

    @classmethod
    def from_dict(cls, d) -> "Catalogue":
        keys = ("header", "conventions", "sources", "records", "alleles", "claims",
                "entities", "relationships", "findings")
        _obj(d, keys, "catalogue")
        conventions = _obj(d["conventions"], ("coordinates", "v_alignment"), "conventions")
        if conventions["coordinates"] != "zero_based_half_open":
            raise _bad("conventions.coordinates", "unsupported coordinate convention", "unsupported_convention")
        _enum(conventions["v_alignment"], V_ALIGNMENTS, "conventions.v_alignment")

        def many(name, kind):
            return tuple(kind.from_dict(item, f"{name}[{i}]") for i, item in enumerate(_list(d[name], name)))

        catalogue = cls(ReferenceHeader.from_dict(d["header"]), dict(conventions), many("sources", SourceFile),
                        many("records", RecordEntry), many("alleles", Allele), many("claims", Claim),
                        many("entities", Entity), many("relationships", Relationship), many("findings", Finding))
        catalogue.validate()
        return catalogue

    def validate(self) -> None:
        """Check referential and semantic integrity; raise ``ReferenceArtifactError("integrity")``."""
        allowed_segments = ("V", "J") if self.header.chain_type == "vj" else SEGMENTS
        source_ids = _unique(self.sources, "sources")
        _unique(self.records, "records")
        _unique(self.alleles, "alleles")
        _unique(self.entities, "entities")
        records = {r.id: r for r in self.records}
        alleles = {a.id: a for a in self.alleles}
        gapped_ok = self.conventions["v_alignment"] == "imgt_gapped"

        for i, a in enumerate(self.alleles):
            p = f"alleles[{i}]"
            if a.segment not in allowed_segments:
                raise _bad(f"{p}.segment", f"segment {a.segment} not allowed for chain type "
                           f"{self.header.chain_type}", "integrity")
            if a.id != f"{a.segment}:{a.label}":
                raise _bad(f"{p}.id", f"id must be {a.segment}:{a.label}", "integrity")
            for field in ("label", "gene"):
                problem = check_label(getattr(a, field))
                if problem:
                    raise _bad(f"{p}.{field}", problem, "integrity")
            if a.gene_basis != "declared" and (a.gene, a.gene_basis) != gene_from_label(a.label):
                expected, basis = gene_from_label(a.label)
                raise _bad(f"{p}.gene", f"gene {a.gene!r} ({a.gene_basis}) does not follow from label "
                           f"{a.label!r}; expected {expected!r} ({basis})", "integrity")
            if not _SEQ_RE.fullmatch(a.sequence):
                raise _bad(f"{p}.sequence", "sequence must match [ACGT]{1,65535}", "integrity")
            if a.sequence_sha256 != sequence_digest(a.sequence):
                raise _bad(f"{p}.sequence_sha256", "digest does not match sequence", "integrity")
            if a.gapped_sequence is not None:
                if a.segment != "V" or not gapped_ok:
                    raise _bad(f"{p}.gapped_sequence", "gapped sequence only allowed for V under imgt_gapped",
                               "integrity")
                if a.gapped_sequence.replace(".", "") != a.sequence:
                    raise _bad(f"{p}.gapped_sequence", "removing dots does not give the sequence", "integrity")
            if a.transforms and a.transforms[-1].output_sha256 != a.sequence_sha256:
                raise _bad(f"{p}.transforms", "last transform output does not match sequence digest", "integrity")
            for rid in a.records:
                if rid not in records:
                    raise _bad(f"{p}.records", f"unknown record {rid!r}", "integrity")
                if records[rid].allele != a.id:
                    raise _bad(f"{p}.records", f"record {rid!r} does not point back at {a.id!r}", "integrity")

        for i, r in enumerate(self.records):
            p = f"records[{i}]"
            match = _REC_ID_RE.fullmatch(r.id)
            if not match or int(match.group(2)) != r.index:
                raise _bad(f"{p}.id", "record id must be rec:<segment>:<index>", "integrity")
            if match.group(1) not in allowed_segments:
                raise _bad(f"{p}.id", f"segment {match.group(1)} not allowed for chain type "
                           f"{self.header.chain_type}", "integrity")
            if r.source not in source_ids:
                raise _bad(f"{p}.source", f"unknown source {r.source!r}", "integrity")
            if r.allele is not None and r.allele in alleles:
                linked = alleles[r.allele]
                if r.id not in linked.records:
                    raise _bad(f"{p}.allele", f"allele {linked.id!r} does not list record {r.id!r}", "integrity")
                if match.group(1) != linked.segment:
                    raise _bad(f"{p}.id", "record segment differs from its allele's segment", "integrity")
            if r.status == "accepted":
                if r.allele not in alleles:
                    raise _bad(f"{p}.allele", "accepted record must point at an existing allele", "integrity")
            elif r.status == "duplicate":
                target = records.get(r.duplicate_of)
                if target is None or target.status != "accepted":
                    raise _bad(f"{p}.duplicate_of", "must point at an accepted record", "integrity")
            else:
                if not r.reason or not r.reason.strip():
                    raise _bad(f"{p}.reason", "rejected record needs a reason", "integrity")
                if r.allele is not None:
                    raise _bad(f"{p}.allele", "rejected record must not have an allele", "integrity")

        for i, c in enumerate(self.claims):
            p = f"claims[{i}]"
            allele = alleles.get(c.subject)
            if allele is None:
                raise _bad(f"{p}.subject", f"unknown allele {c.subject!r}", "integrity")
            if c.predicate == "v_subregions" and allele.segment != "V":
                raise _bad(f"{p}.predicate", "v_subregions applies to V alleles only", "integrity")
            check_claim_value(c.predicate, c.value, len(allele.sequence), f"{p}.value")
            if c.support == "source" and (c.source is None or not c.source.strip()):
                raise _bad(f"{p}.source", "source-supported claim needs a source", "integrity")
            if c.derivation == "imgt_alignment" and not (
                    self.conventions["v_alignment"] == "imgt_gapped" and allele.segment == "V"
                    and allele.gapped_sequence is not None):
                raise _bad(f"{p}.derivation", "imgt_alignment requires a V allele with a gapped sequence "
                           "under a declared imgt_gapped alignment", "integrity")
            if c.predicate == "v_subregions" and c.derivation != "imgt_alignment":
                raise _bad(f"{p}.derivation", "v_subregions must derive from imgt_alignment", "integrity")
            if c.derivation == "anchor_derived" and c.predicate != "coding_phase":
                raise _bad(f"{p}.derivation", "anchor_derived applies only to coding_phase", "integrity")
            if (c.support == "none") != (c.status == "candidate"):
                raise _bad(f"{p}.status", "support none holds exactly for candidate claims", "integrity")
        disagreements = anchor_phase_disagreements(self.claims)
        if disagreements:
            index, message = disagreements[0]
            raise _bad(f"claims[{index}].value", message, "integrity")

        for i, f in enumerate(self.findings):
            for subject in f.subjects:
                if subject not in alleles:
                    raise _bad(f"findings[{i}].subjects", f"unknown allele {subject!r}", "integrity")

        entities = {e.id: e for e in self.entities}
        for i, e in enumerate(self.entities):
            if e.executable:
                raise _bad(f"entities[{i}].executable", "entities are not executable", "integrity")
            problem = check_label(e.name)
            if problem:
                raise _bad(f"entities[{i}].name", problem, "integrity")
        for i, rel in enumerate(self.relationships):
            p = f"relationships[{i}]"
            allele = alleles.get(rel.allele)
            if allele is None:
                raise _bad(f"{p}.allele", f"unknown allele {rel.allele!r}", "integrity")
            # The target is a named entity, or an executable allele of the same segment.
            if rel.entity in entities:
                name = entities[rel.entity].name
            elif rel.entity in alleles and alleles[rel.entity].segment == allele.segment:
                name = alleles[rel.entity].label
            else:
                raise _bad(f"{p}.entity", f"unknown entity {rel.entity!r}; expected an entity id or an "
                           f"executable {allele.segment} allele id", "integrity")
            if name == allele.label:
                raise _bad(f"{p}.entity", f"{rel.kind} of {allele.id} names its own label", "integrity")
            if rel.kind == "alias" and any(o.segment == allele.segment and o.label == name
                                           for o in self.alleles):
                raise _bad(f"{p}.entity", f"alias {name!r} targets another executable allele label",
                           "integrity")
