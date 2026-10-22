"""The portable ``Reference``: build, inspect, write and reload a local germline reference.

A reference is built from local V/(D)/J FASTA files and an annotation sidecar, and is
written as a closed, canonical ``*.genairr-ref.json`` artifact.  Loading an artifact
regenerates every derived section (derived catalogue content such as alignment and motif
claims, conflict statuses, entities and findings; hashes, view, defaults, admission, native
projection, context and execution hash) from base data and compares each one for exact
equality, so an artifact whose hashes were recomputed over inconsistent content is still
rejected.

Because a stored artifact loads only if this regeneration reproduces it exactly, the
outputs of the derivation are part of the version 1 file format: the motif resolver's
anchor positions and method strings, the IMGT region boundary table, the kinds of native
validation findings, the serialization of the native reference spec (``ReferenceSpecV1``)
and the unfitted defaults profile.  Changing any of them changes which version 1
artifacts load, so it is a schema change and must be treated as one; the golden artifact
read test (``tests/golden/reference``) guards against changing them by accident.
"""
from __future__ import annotations

import copy
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Optional

from GenAIRR._artifact_io import (
    ArtifactIOError, canonical_json, read_text_bounded, strict_loads, typed_digest, write_atomic,
)
from GenAIRR._reference.capabilities import CAPABILITIES, CapabilityReport, _usable, evaluate
from GenAIRR._reference.defaults import models_hash, recombine_defaults, unfitted_defaults_v1
from GenAIRR._reference.defaults import PROFILE as DEFAULTS_PROFILE
from GenAIRR._reference.errors import ReferenceArtifactError, ReferenceCapabilityError
from GenAIRR._reference.importer import import_local, regenerate, scan_fasta
from GenAIRR._reference.lowering import PROFILE_ID, default_view, lower, seal
from GenAIRR._reference.model import Catalogue, Taxon
from GenAIRR._reference.sidecar import annotation_template as _sidecar_template
from GenAIRR._reference.sidecar import load_sidecar

ARTIFACT = "genairr_reference"
SCHEMA_VERSION = 1
DIGEST_DOMAIN = "genairr_reference_artifact_v1"
MAX_BYTES = 64 * 1024 * 1024
VIEW_ID = "default"
ANNOTATIONS = ("anchor", "functionality", "coding_phase", "boundary_5", "boundary_3")

_TOP_KEYS = ("artifact", "schema_version", "generator", "catalogue", "catalogue_hash", "view",
             "models", "models_hash", "native", "commitment")
_SECTION_KEYS = (("generator", ("name", "version")),
                 ("view", ("id", "profile", "alleles", "admission")),
                 ("native", ("reference_spec", "context", "execution_hash")))


# ── Inspection records ───────────────────────────────────────────────


@dataclass(frozen=True)
class ClaimSummary:
    """One annotation exactly as recorded: value, support, derivation and source.

    When several accepted claims agree, the weakest support is shown; evidence is
    never reported stronger than its weakest recorded claim.
    """
    value: Any
    support: str
    derivation: str
    source: Optional[str]

    def to_dict(self):
        return {"value": self.value, "support": self.support, "derivation": self.derivation,
                "source": self.source}


@dataclass(frozen=True)
class ReferenceAllele:
    """An executable allele of the reference with its accepted annotations.

    Each annotation is ``None`` when no accepted claim establishes it (unknown,
    candidate-only or conflicting).  ``notes`` are the author's free-text notes from
    the annotation sidecar, in sidecar order.
    """
    segment: str
    label: str
    gene: str
    length: int
    sequence_sha256: str
    anchor: Optional[ClaimSummary]
    functionality: Optional[ClaimSummary]
    coding_phase: Optional[ClaimSummary]
    boundary_5: Optional[ClaimSummary]
    boundary_3: Optional[ClaimSummary]
    notes: tuple = ()

    def to_dict(self):
        out = {"segment": self.segment, "label": self.label, "gene": self.gene,
               "length": self.length, "sequence_sha256": self.sequence_sha256}
        for name in ANNOTATIONS:
            summary = getattr(self, name)
            out[name] = None if summary is None else summary.to_dict()
        out["notes"] = list(self.notes)
        return out


def _summary(catalogue: Catalogue, allele_id: str, predicate: str) -> Optional[ClaimSummary]:
    claims = _usable(catalogue, allele_id, predicate)
    if not claims:
        return None
    weakest = next((c for c in claims if c.support != "source"), claims[0])
    return ClaimSummary(weakest.value, weakest.support, weakest.derivation, weakest.source)


def _conflicted(catalogue: Catalogue, allele_id: str, predicate: str) -> bool:
    claims = catalogue.claims_for(allele_id, predicate)
    accepted = {canonical_json(c.value) for c in claims if c.status == "accepted"}
    return any(c.status == "conflict" for c in claims) or len(accepted) > 1


def _record_counts(records) -> Dict[str, Dict[str, int]]:
    counts: Dict[str, Dict[str, int]] = {}
    for segment in ("V", "D", "J"):
        statuses = [r.status for r in records if r.id.split(":")[1] == segment]
        if statuses:
            counts[segment] = {s: statuses.count(s) for s in ("accepted", "duplicate", "rejected")}
    return counts


def _support_text(values) -> str:
    if (len(values) > 2 and all(type(v) is int for v in values)
            and list(values) == list(range(values[0], values[-1] + 1))):
        return f"{values[0]}..{values[-1]}"
    return ", ".join(str(v) for v in values)


def _describe_defaults(models: dict) -> tuple:
    usage = models["allele_usage"]
    counts = ", ".join(f"{s} {len(usage[s])}" for s in ("V", "D", "J") if s in usage)
    planes = [("allele_usage", f"{usage['kind']} over {usage['unit']}s ({counts})")]
    for plane in ("np_lengths", "np_bases"):
        planes.append((plane, "; ".join(f"{key} {spec['kind']} over {_support_text(spec['support'])}"
                                        for key, spec in models[plane].items())))
    for plane in ("trims", "p_nucleotides"):
        planes.append((plane, "; ".join(f"{key} {spec['kind']}" for key, spec in models[plane].items())))
    return tuple(planes)


def _cell(summary: Optional[ClaimSummary], conflicted: bool) -> str:
    if conflicted:
        return "conflict"
    if summary is None:
        return "unknown"
    return f"{summary.value} {summary.support}/{summary.derivation}"


@dataclass(frozen=True)
class ReferenceReport:
    """A human- and machine-readable account of a reference.

    ``str(report)`` covers identity, hashes, FASTA record accounting, the allele table
    with recorded support levels, findings, the defaults profile and the capability
    explanation; ``to_dict()`` gives the same content as plain JSON-compatible data.
    """
    name: str
    version: Optional[str]
    taxon: Taxon
    locus: str
    chain_type: str
    build_date: Optional[str]
    catalogue_hash: str
    models_hash: str
    execution_hash: str
    artifact_hash: str
    records: tuple
    alleles: tuple
    conflicts: tuple
    entities: tuple
    relationships: tuple
    findings: tuple
    admitted_findings: tuple
    defaults_profile: str
    defaults: tuple
    capabilities: CapabilityReport

    def to_dict(self) -> dict:
        return {
            "identity": {"name": self.name, "version": self.version, "taxon": self.taxon.to_dict(),
                         "locus": self.locus, "chain_type": self.chain_type, "build_date": self.build_date},
            "hashes": {"catalogue_hash": self.catalogue_hash, "models_hash": self.models_hash,
                       "execution_hash": self.execution_hash, "artifact_hash": self.artifact_hash},
            "records": {"counts": _record_counts(self.records),
                        "not_accepted": [r.to_dict() for r in self.records if r.status != "accepted"]},
            "alleles": [a.to_dict() for a in self.alleles],
            "conflicts": [{"allele": a, "predicate": p} for a, p in self.conflicts],
            "entities": [e.to_dict() for e in self.entities],
            "relationships": [r.to_dict() for r in self.relationships],
            "findings": [f.to_dict() for f in self.findings],
            "admitted_findings": [{"kind": k, "segment": s, "allele": a} for k, s, a in self.admitted_findings],
            "defaults": {"profile": self.defaults_profile, "planes": dict(self.defaults)},
            "capabilities": self.capabilities.to_dict(),
        }

    def _allele_table(self) -> list:
        conflicts = set(self.conflicts)
        rows = [("segment", "label", "gene", "length") + ANNOTATIONS]
        for a in self.alleles:
            allele_id = f"{a.segment}:{a.label}"
            rows.append((a.segment, a.label, a.gene, str(a.length)) + tuple(
                _cell(getattr(a, p), (allele_id, p) in conflicts) for p in ANNOTATIONS))
        widths = [max(len(row[i]) for row in rows) for i in range(len(rows[0]))]
        return ["  " + "  ".join(c.ljust(w) for c, w in zip(row, widths)).rstrip() for row in rows]

    def __str__(self) -> str:
        identity = self.name if self.version is None else f"{self.name}@{self.version}"
        taxon = self.taxon.scientific_name
        if self.taxon.taxonomy_id is not None:
            taxon += f" (taxonomy {self.taxon.taxonomy_id})"
        lines = [f"Reference {identity}",
                 f"  taxon: {taxon}",
                 f"  locus: {self.locus}  chain type: {self.chain_type}  "
                 f"build date: {self.build_date or 'not stated'}",
                 "Hashes",
                 f"  catalogue: {self.catalogue_hash}",
                 f"  models:    {self.models_hash}",
                 f"  execution: {self.execution_hash}",
                 f"  artifact:  {self.artifact_hash}",
                 "Records"]
        for segment, c in _record_counts(self.records).items():
            lines.append(f"  {segment}: {sum(c.values())} records: {c['accepted']} accepted, "
                         f"{c['duplicate']} duplicate, {c['rejected']} rejected")
        for r in self.records:
            name = r.label if r.label is not None else repr(r.header)
            if r.status == "duplicate":
                lines.append(f"  {r.id} {name}: duplicate of {r.duplicate_of}")
            elif r.status == "rejected":
                lines.append(f"  {r.id} {name}: rejected ({r.reason})")
        lines.append("Alleles (view order; each annotation shows value support/derivation as recorded)")
        lines.extend(self._allele_table())
        sources = [f"  {a.segment}:{a.label} {p}: {getattr(a, p).source}" for a in self.alleles
                   for p in ANNOTATIONS if getattr(a, p) is not None and getattr(a, p).source]
        if sources:
            lines.append("Sources")
            lines.extend(sources)
        notes = [f"  {a.segment}:{a.label}: {note}" for a in self.alleles for note in a.notes]
        if notes:
            lines.append("Notes")
            lines.extend(notes)
        lines.append("Entities (named, not executable)")
        for e in self.entities:
            links = [f"{r.kind} of {r.allele}" for r in self.relationships if r.entity == e.id]
            lines.append(f"  {e.name}: {', '.join(links) or 'no relationship'}")
        if not self.entities:
            lines.append("  (none)")
        entity_ids = {e.id for e in self.entities}
        allele_links = [r for r in self.relationships if r.entity not in entity_ids]
        if allele_links:
            lines.append("Relationships between alleles (both executable)")
            lines.extend(f"  {r.allele} {r.kind} of {r.entity}" for r in allele_links)
        lines.append("Findings")
        lines.extend(f"  {f.kind}: {', '.join(f.subjects)} — {f.detail}" for f in self.findings)
        if not self.findings:
            lines.append("  (none)")
        lines.append(f"Admitted engine findings (profile {PROFILE_ID})")
        lines.extend(f"  {kind}: {allele or segment}" for kind, segment, allele in self.admitted_findings)
        if not self.admitted_findings:
            lines.append("  (none)")
        lines.append(f"Defaults (profile {self.defaults_profile}; placeholders, not fitted to data)")
        lines.extend(f"  {plane}: {text}" for plane, text in self.defaults)
        lines.append("Capabilities")
        lines.append(self.capabilities.explain())
        return "\n".join(lines)


# ── Loading helpers ──────────────────────────────────────────────────


def _same(a, b) -> bool:
    """Exact JSON equality: types and list order matter (1 differs from 1.0 and true)."""
    return canonical_json(a) == canonical_json(b)


def _require_keys(value, keys, path):
    if not isinstance(value, dict):
        raise ReferenceArtifactError("malformed", path, "expected an object")
    missing = [k for k in keys if k not in value]
    unknown = sorted(k for k in value if k not in keys)
    if missing:
        raise ReferenceArtifactError("malformed", path, f"missing keys {missing}")
    if unknown:
        raise ReferenceArtifactError("malformed", path, f"unknown keys {unknown}")


def _check(ok: bool, path: str, message: str) -> None:
    if not ok:
        raise ReferenceArtifactError("integrity", path, message)


def _compare_derived(section: str, stored: tuple, regenerated: tuple) -> None:
    """Require a stored derived section to equal its regeneration, item by item and in order."""
    noun = {"claims": "claim", "relationships": "relationship", "entities": "entity",
            "findings": "finding"}[section]
    for index, (have, want) in enumerate(zip(stored, regenerated)):
        if not _same(have.to_dict(), want.to_dict()):
            raise ReferenceArtifactError(
                "integrity", f"catalogue.{section}[{index}]",
                f"{noun} {canonical_json(have.to_dict())} differs from its regeneration from base data: "
                f"{canonical_json(want.to_dict())}")
    if len(stored) > len(regenerated):
        extra = stored[len(regenerated)]
        raise ReferenceArtifactError(
            "integrity", f"catalogue.{section}[{len(regenerated)}]",
            f"{noun} {canonical_json(extra.to_dict())} does not follow from the catalogue's base data")
    if len(stored) < len(regenerated):
        missing = regenerated[len(stored)]
        raise ReferenceArtifactError(
            "integrity", f"catalogue.{section}",
            f"{noun} {canonical_json(missing.to_dict())} follows from the catalogue's base data but is missing")


def _require_selection(catalogue: Catalogue, view: dict) -> None:
    """Every segment the chain type needs must have at least one executable allele."""
    for segment, ids in view.items():
        if ids:
            continue
        records = [r for r in catalogue.records if r.id.split(":")[1] == segment]
        counts = {s: sum(r.status == s for r in records) for s in ("accepted", "duplicate", "rejected")}
        reasons: Dict[str, int] = {}
        for r in records:
            if r.status == "rejected":
                key = r.reason.split(":", 1)[0].split(" ", 1)[0]
                reasons[key] = reasons.get(key, 0) + 1
        detail = ("rejection reasons: " + ", ".join(f"{k} ({n})" for k, n in reasons.items())
                  if reasons else "the FASTA file has no records")
        raise ReferenceArtifactError(
            "empty_selection", segment,
            f"no executable {segment} alleles: {counts['accepted']} accepted, {counts['duplicate']} duplicate, "
            f"{counts['rejected']} rejected; {detail}")


def _context(catalogue_hash: str, models_hash_: str) -> dict:
    return {"schema_version": 1, "catalogue_hash": catalogue_hash, "view": VIEW_ID,
            "profile": PROFILE_ID, "models_hash": models_hash_}


def _generator() -> dict:
    import GenAIRR
    return {"name": "GenAIRR", "version": GenAIRR.__version__}


# ── The Reference ────────────────────────────────────────────────────


class Reference:
    """A portable, immutable germline reference built from local FASTA and annotations.

    Create one with :meth:`from_fasta` (from FASTA files and an annotation sidecar),
    :meth:`read` or :meth:`from_json` (from a ``*.genairr-ref.json`` artifact).  A
    reference never changes after construction: getters return tuples or copies, and the
    native configuration it hands to the engine is sealed.
    """

    __slots__ = ("_catalogue", "_view", "_models", "_models_hash", "_catalogue_hash", "_admission",
                 "_refdata", "_execution_hash", "_alleles", "_capabilities", "_text", "_commitment")

    def __init__(self, *args, **kwargs):
        raise TypeError("create a Reference with Reference.from_fasta, Reference.read or Reference.from_json")

    def __setattr__(self, name, value):
        raise AttributeError("Reference is immutable")

    def __delattr__(self, name):
        raise AttributeError("Reference is immutable")

    def __copy__(self):
        return self

    def __deepcopy__(self, memo):
        return self

    def __repr__(self) -> str:
        pools = ", ".join(f"{segment} {len(ids)}" for segment, ids in self._view.items())
        return f"<Reference {self.name!r} {self.locus} {self.chain_type} ({pools}) {self._execution_hash}>"

    # Construction ------------------------------------------------------

    @classmethod
    def _assemble(cls, catalogue, view, lowered, models, models_hash_, context, refdata, generator):
        execution_hash = refdata.content_hash()
        body = {
            "artifact": ARTIFACT, "schema_version": SCHEMA_VERSION, "generator": generator,
            "catalogue": catalogue.to_dict(), "catalogue_hash": context["catalogue_hash"],
            "view": {"id": VIEW_ID, "profile": PROFILE_ID, "alleles": view, "admission": lowered.admission},
            "models": models, "models_hash": models_hash_,
            "native": {"reference_spec": lowered.spec, "context": context, "execution_hash": execution_hash},
        }
        commitment = typed_digest(body, DIGEST_DOMAIN)
        alleles = tuple(
            ReferenceAllele(a.segment, a.label, a.gene, len(a.sequence), a.sequence_sha256,
                            *(_summary(catalogue, a.id, p) for p in ANNOTATIONS), notes=a.notes)
            for a in (catalogue.allele(aid) for ids in view.values() for aid in ids))
        self = object.__new__(cls)
        for name, value in (
                ("_catalogue", catalogue), ("_view", copy.deepcopy(view)), ("_models", copy.deepcopy(models)),
                ("_models_hash", models_hash_), ("_catalogue_hash", context["catalogue_hash"]),
                ("_admission", copy.deepcopy(lowered.admission)), ("_refdata", refdata),
                ("_execution_hash", execution_hash), ("_alleles", alleles),
                ("_capabilities", evaluate(catalogue, view)),
                ("_text", canonical_json({**body, "commitment": commitment}) + "\n"),
                ("_commitment", commitment)):
            object.__setattr__(self, name, value)
        return self

    @classmethod
    def from_fasta(cls, *, v, j, annotations, d=None) -> "Reference":
        """Build a reference from V, J (and, for ``vdj``, D) FASTA files and an annotation sidecar.

        Every FASTA record is accounted for; the default view holds every accepted allele in
        file order, and the explicit unfitted defaults profile supplies the models.
        """
        catalogue = import_local(v=v, j=j, d=d, annotations=annotations)
        view = default_view(catalogue)
        _require_selection(catalogue, view)
        lowered = lower(catalogue, view)
        models = unfitted_defaults_v1(catalogue.header.chain_type, view)
        models_hash_ = models_hash(models)
        context = _context(catalogue.digest(), models_hash_)
        return cls._assemble(catalogue, view, lowered, models, models_hash_, context,
                             seal(lowered, context), _generator())

    @classmethod
    def read(cls, path) -> "Reference":
        """Read and fully check a ``*.genairr-ref.json`` artifact (at most 64 MiB)."""
        try:
            text = read_text_bounded(path, max_bytes=MAX_BYTES)
        except ArtifactIOError as exc:
            kind = "malformed" if exc.kind == "format" else "io"
            raise ReferenceArtifactError(kind, Path(path).name, str(exc)) from exc
        return cls.from_json(text)

    @classmethod
    def from_json(cls, text: str) -> "Reference":
        """Load an artifact from text, regenerating and comparing every derived section."""
        # 1. Size bound, then strict parsing (duplicate keys, non-finite numbers, trailing data).
        if not isinstance(text, str):
            raise ReferenceArtifactError("malformed", "", "expected JSON text")
        if len(text.encode("utf-8", "surrogatepass")) > MAX_BYTES:
            raise ReferenceArtifactError("io", "", f"artifact exceeds the size limit of {MAX_BYTES} bytes")
        try:
            doc = strict_loads(text, max_bytes=MAX_BYTES)
        except ArtifactIOError as exc:
            raise ReferenceArtifactError("io" if exc.kind == "resource" else "malformed", "", str(exc)) from exc

        # 2. Artifact kind, then schema version, before any key-set check, so that an artifact
        #    written under a newer schema is reported as unsupported rather than malformed.
        if not isinstance(doc, dict):
            raise ReferenceArtifactError("malformed", "", "expected a JSON object")
        if doc.get("artifact") != ARTIFACT:
            raise ReferenceArtifactError("malformed", "artifact",
                                         f"not a GenAIRR reference artifact (expected artifact {ARTIFACT!r})")
        version = doc.get("schema_version")
        if type(version) is not int or version != SCHEMA_VERSION:
            raise ReferenceArtifactError("unsupported_schema", "schema_version",
                                         f"unsupported schema version {version!r}; "
                                         f"this GenAIRR reads version {SCHEMA_VERSION}")

        # 3. Closed key sets, then the generator.
        _require_keys(doc, _TOP_KEYS, "")
        for section, keys in _SECTION_KEYS:
            _require_keys(doc[section], keys, section)
        generator = doc["generator"]
        if generator["name"] != "GenAIRR":
            raise ReferenceArtifactError("malformed", "generator.name", "expected 'GenAIRR'")
        if not isinstance(generator["version"], str) or not generator["version"].strip():
            raise ReferenceArtifactError("malformed", "generator.version", "expected a non-blank string")

        # 4. Commitment over every other key.
        body = {k: v for k, v in doc.items() if k != "commitment"}
        _check(doc["commitment"] == typed_digest(body, DIGEST_DOMAIN), "commitment",
               "commitment does not match the artifact contents")

        # 5. Closed, referential and semantic validation of the catalogue.
        try:
            catalogue = Catalogue.from_dict(doc["catalogue"])
        except ReferenceArtifactError as exc:
            path = exc.path if exc.path.startswith("catalogue") else (
                f"catalogue.{exc.path}" if exc.path else "catalogue")
            raise ReferenceArtifactError(exc.kind, path, exc.message) from exc

        # 5b. Derived catalogue content (claims with their conflict statuses, relationships,
        #     entities, findings) must equal its regeneration from the catalogue's base data.
        derived = regenerate(catalogue)
        for section in ("claims", "relationships", "entities", "findings"):
            _compare_derived(section, getattr(catalogue, section), getattr(derived, section))

        # 6. Catalogue hash.
        catalogue_hash = catalogue.digest()
        _check(doc["catalogue_hash"] == catalogue_hash, "catalogue_hash", "does not match the catalogue")

        # 7. View: identity, profile and the regenerated default view.
        stored_view = doc["view"]
        _check(stored_view["id"] == VIEW_ID, "view.id", f"view id must be {VIEW_ID!r}")
        _check(stored_view["profile"] == PROFILE_ID, "view.profile", f"view profile must be {PROFILE_ID!r}")
        view = default_view(catalogue)
        _check(_same(stored_view["alleles"], view), "view.alleles",
               "view differs from the catalogue's default view (every accepted allele, in record order)")

        # 8. Models: exactly the unfitted defaults profile for this view.
        models = unfitted_defaults_v1(catalogue.header.chain_type, view)
        _check(_same(doc["models"], models), "models", f"models differ from profile {DEFAULTS_PROFILE}")
        models_hash_ = models_hash(models)
        _check(doc["models_hash"] == models_hash_, "models_hash", "does not match the models")

        # 9. Lowering: admission and native projection.
        lowered = lower(catalogue, view)
        _check(_same(stored_view["admission"], lowered.admission), "view.admission",
               "admission differs from catalogue lowering")
        _check(_same(doc["native"]["reference_spec"], lowered.spec), "native.reference_spec",
               "native projection differs from catalogue lowering")

        # 10. Native context.
        context = _context(catalogue_hash, models_hash_)
        _check(_same(doc["native"]["context"], context), "native.context",
               "context differs from the catalogue hash, view, profile and models hash")

        # 11. Seal, then the execution hash.
        refdata = seal(lowered, context)
        _check(doc["native"]["execution_hash"] == refdata.content_hash(), "native.execution_hash",
               "does not match the sealed native reference")

        reference = cls._assemble(catalogue, view, lowered, models, models_hash_, context, refdata,
                                  {"name": generator["name"], "version": generator["version"]})
        _check(reference._text == canonical_json(doc) + "\n", "", "artifact differs from its regenerated form")
        return reference

    @staticmethod
    def annotation_template(*, v, j, d=None, reference: dict, v_alignment: str = "none") -> dict:
        """Return an annotation sidecar for these FASTA files with every annotation unknown.

        ``reference`` is the sidecar's ``reference`` header (name, version, locus,
        chain_type, taxon, build_date).  No support is claimed; unique J motif positions
        are listed under ``suggestions`` as candidates to verify, never as annotations.
        """
        conventions = {"coordinates": "zero_based_half_open", "v_alignment": v_alignment}
        # Check the header and conventions before reading any FASTA file.
        load_sidecar(canonical_json(_sidecar_template(
            records_by_segment={}, header=reference, conventions=conventions, suggestions=[])))
        chain_type = reference["chain_type"]
        if chain_type == "vj" and d is not None:
            raise ReferenceArtifactError("malformed", "d", "a D FASTA file is not allowed for chain_type vj")
        if chain_type == "vdj" and d is None:
            raise ReferenceArtifactError("malformed", "d", "chain_type vdj requires a D FASTA file")
        records, suggestions = scan_fasta(v=v, j=j, d=d, v_alignment=v_alignment, locus=reference["locus"])
        template = _sidecar_template(records_by_segment=records, header=reference,
                                     conventions=conventions, suggestions=suggestions)
        load_sidecar(canonical_json(template))
        return template

    # Serialization -----------------------------------------------------

    def to_json(self) -> str:
        """The artifact text: canonical JSON followed by a newline."""
        return self._text

    def write(self, path, *, replace: bool = False) -> None:
        """Write the artifact atomically; an existing file is kept unless ``replace=True``."""
        try:
            write_atomic(path, self._text.encode("utf-8"), replace=replace)
        except ArtifactIOError as exc:
            raise ReferenceArtifactError("io", Path(path).name, str(exc)) from exc

    # Identity ----------------------------------------------------------

    @property
    def name(self) -> str:
        return self._catalogue.header.name

    @property
    def chain_type(self) -> str:
        return self._catalogue.header.chain_type

    @property
    def locus(self) -> str:
        return self._catalogue.header.locus

    @property
    def taxon(self) -> Taxon:
        return self._catalogue.header.taxon

    @property
    def catalogue_hash(self) -> str:
        """Digest of the catalogue (records, alleles, claims, entities, findings)."""
        return self._catalogue_hash

    @property
    def models_hash(self) -> str:
        """Digest of the defaults record used for simulation."""
        return self._models_hash

    @property
    def execution_hash(self) -> str:
        """Content hash of the sealed native reference the engine runs against."""
        return self._execution_hash

    @property
    def artifact_hash(self) -> str:
        """The artifact's commitment over all of its other contents."""
        return self._commitment

    # Contents ----------------------------------------------------------

    @property
    def alleles(self) -> tuple:
        """Executable alleles (:class:`ReferenceAllele`) in view order: V, then D, then J."""
        return self._alleles

    @property
    def records(self) -> tuple:
        """Every FASTA record with its status (accepted, duplicate or rejected with a reason)."""
        return self._catalogue.records

    @property
    def entities(self) -> tuple:
        """Named, non-executable entities such as aliases and paralog names."""
        return self._catalogue.entities

    @property
    def findings(self) -> tuple:
        """Catalogue findings, such as identical sequences under different labels."""
        return self._catalogue.findings

    @property
    def admission(self) -> dict:
        """A copy of the engine findings admitted by the observed-span profile."""
        return copy.deepcopy(self._admission)

    # Capabilities and reporting ---------------------------------------

    def capabilities(self) -> CapabilityReport:
        """Which workflows this reference supports, and at which evidence level."""
        return self._capabilities

    def require(self, *names: str, evidence: str = "assertion") -> None:
        """Raise :class:`ReferenceCapabilityError` unless every named capability is available.

        ``evidence="source"`` additionally requires that the annotations behind each
        capability are source-supported rather than the author's assertions.
        """
        report = self._capabilities
        available = report.available(evidence=evidence)
        for name in names:
            if name not in CAPABILITIES:
                raise ValueError(f"unknown capability {name!r}; expected one of {list(CAPABILITIES)}")
        for name in names:
            if name not in available:
                raise ReferenceCapabilityError("capability", name, report.explain())

    def report(self) -> ReferenceReport:
        """Identity, hashes, record accounting, alleles with their evidence, defaults and capabilities."""
        catalogue = self._catalogue
        header = catalogue.header
        view_ids = [aid for ids in self._view.values() for aid in ids]
        return ReferenceReport(
            name=header.name, version=header.version, taxon=header.taxon, locus=header.locus,
            chain_type=header.chain_type, build_date=header.build_date,
            catalogue_hash=self._catalogue_hash, models_hash=self._models_hash,
            execution_hash=self._execution_hash, artifact_hash=self._commitment,
            records=catalogue.records, alleles=self._alleles,
            conflicts=tuple((aid, p) for aid in view_ids for p in ANNOTATIONS
                            if _conflicted(catalogue, aid, p)),
            entities=catalogue.entities, relationships=catalogue.relationships,
            findings=catalogue.findings,
            admitted_findings=tuple((f["kind"], f["segment"], f["allele"])
                                    for f in self._admission["accepted_curatable_findings"]),
            defaults_profile=self._models["profile"], defaults=_describe_defaults(self._models),
            capabilities=self._capabilities)

    # Runtime -----------------------------------------------------------

    def native_refdata(self):
        """The sealed native reference configuration (built once; it rejects every mutation)."""
        return self._refdata

    def recombine_defaults(self) -> dict:
        """A fresh copy of the unfitted defaults in the form a recombination step consumes."""
        return recombine_defaults(self._models, chain_type=self.chain_type)
