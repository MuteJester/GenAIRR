"""OGRDB germline-set adapter: set identity and the normalized source inventory.

OGRDB API v2 serves a germline set as AIRR germline-set JSON (``airr``, and ``airr_ex``
with the provider's ``_ex`` content) and as FASTA.  A JSON export is an envelope
``{"Info": {...}, "GermlineSet": [<one germline set>]}`` whose set holds the taxon
(``species``), the ``locus``, set identifiers and the records (``allele_descriptions``).

Provider payloads are third-party data: they are decoded as UTF-8 and parsed only with
the strict JSON loader (duplicate keys, non-finite numbers, trailing data and nesting
deeper than 64 are rejected), and never executed.

:func:`identify` reads the scientific identity of the set (taxon, locus, chain type)
from the payload envelope itself.  It is usable as the acquisition layer's identity
reader: it never trusts the requested set id or URL, checks the payload's own set and
release identifiers against the request wherever the payload states them, and requires
every record export present to agree.

:func:`build_inventory` reads a verified lock's bytes into the normalized inventory.  It
interprets no coordinates.  Every record occurrence becomes a :class:`SourceRecord`
(malformed and unsupported ones included, with a problem and a finding), every supplied
coding sequence a :class:`SequenceCandidate` of its own record, every name that appears
only as a paralog a :class:`MentionedIdentity`, every embedded support entry a
:class:`SupportReference`, and every provider field an :class:`Observation` at its
locator.  Fields the inventory normalizes are ``interpreted`` (or ``invalid``), other
known fields are ``preserved`` uninterpreted, and fields this adapter does not know are
``unsupported``: kept with their values and recorded by a finding, never an error.
Records are identified by asset and JSON pointer, so duplicate or missing provider ids
never merge records, and evidence never moves between identities because their DNA
matches.  The two record exports are compared label by label without preferring either;
a lock with one record export is inspected with the comparison marked incomplete.  Only
a malformed global envelope or an exceeded budget raises.

:func:`account` reconciles an inventory to the exports it was built from, and
:func:`comparison` gives the build report's comparison section.
"""
from __future__ import annotations

import hashlib
import math
import re
from dataclasses import dataclass, field
from typing import Any, Mapping, Optional

from GenAIRR._artifact_io import canonical_json
from GenAIRR._reference_sources._fields import (
    MAX_ASSET_BYTES, MAX_NAME, MAX_REFS, MAX_SEQUENCE, MAX_SOURCE_RECORDS, MAX_TAXONOMY_ID, MAX_TEXT, bad,
    check_label, check_text, parse_json, ptr, require_type, show,
)
from GenAIRR._reference_sources.common import (
    ADAPTER_ID, LOCUS_CHAIN_TYPES, RECORD_EXPORT_FORMATS, SEGMENTS, Locator, SourceIdentity, Taxon,
)
from GenAIRR._reference_sources.digests import lock_hash
from GenAIRR._reference_sources.errors import ReferenceSourceError
from GenAIRR._reference_sources.inventory import (
    MAX_LOCUS, MAX_OBSERVATIONS, Finding, MentionedIdentity, Mention, Observation, SequenceCandidate,
    SourceInventory, SourceRecord, SupportReference, text_digest,
)
from GenAIRR._reference_sources.lock import SourceLock
from GenAIRR._reference_sources.report import Comparison
from GenAIRR._reference.model import sequence_digest

# ── Provider conventions ────────────────────────────────────────────

SET_POINTER = "/GermlineSet/0"
RECORDS_POINTER = SET_POINTER + "/allele_descriptions"
TAXON_PREFIX = "NCBITAXON:"
# The extended export names its set with this suffix (``IGH_VDJ`` / ``IGH_VDJ_extended``).
EXTENDED_SET_SUFFIX = "_extended"
# The provider's suffix on record refs of the extended export.
EXT_SUFFIX = "_ext"

_TAXON_RE = re.compile(r"NCBITAXON:([1-9][0-9]{0,9})")
_SET_ID_RE = re.compile(r"([1-9][0-9]{0,9})\.(.+)")


def _require_assets(assets) -> Mapping:
    if not isinstance(assets, Mapping):
        raise TypeError(f"assets must be a mapping of format to bytes, got {type(assets).__name__}")
    for name, data in assets.items():
        if type(data) is not bytes:
            raise TypeError(f"assets[{name!r}] must be bytes, got {type(data).__name__}")
    return assets


# ── Envelope ────────────────────────────────────────────────────────


@dataclass(frozen=True)
class _Envelope:
    """One parsed record export: its root object, its single germline set and set identity."""

    asset: str
    root: dict
    germline_set: dict
    identity: SourceIdentity
    set_id: Optional[str]          # germline_set_id, where stated
    set_ref: Optional[str]         # germline_set_ref, where stated
    set_name: Optional[str]        # germline_set_name, where stated
    release: Any                   # release_version, where stated (an int or float)

    @property
    def path(self) -> str:
        return ptr("/assets", self.asset)

    def at(self, pointer: str) -> str:
        """The error path of ``pointer`` inside this export."""
        return self.path + pointer


def _load(asset: str, data: bytes) -> Any:
    """Strictly parse the bytes of one JSON export."""
    path = ptr("/assets", asset)
    if len(data) > MAX_ASSET_BYTES:
        raise bad(path, f"{asset}: {len(data)} bytes exceed the limit of {MAX_ASSET_BYTES}", "resource_limit")
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise bad(path, f"{asset}: the export is not UTF-8 text ({exc})") from None
    try:
        return parse_json(text, max_bytes=MAX_ASSET_BYTES)
    except ReferenceSourceError as exc:
        raise bad(path, f"{asset}: {exc.message}", exc.kind) from None


def taxonomy_id(value: Any) -> Optional[int]:
    """The taxonomy id of an OGRDB species id ``NCBITAXON:<n>``, or ``None``."""
    if type(value) is str:
        match = _TAXON_RE.fullmatch(value)
        if match and int(match.group(1)) <= MAX_TAXONOMY_ID:
            return int(match.group(1))
    return None


def _optional_text(germline_set: dict, key: str, path: str) -> Optional[str]:
    value = germline_set.get(key)
    if value is None:
        return None
    try:
        return check_text(value, ptr(path, key), max_len=MAX_NAME, multiline=False)
    except ReferenceSourceError as exc:
        raise bad(exc.path, f"{key}: {exc.message}", exc.kind) from None


def _envelope(asset: str, data: bytes) -> _Envelope:
    """Parse one record export and validate its global envelope and set identity."""
    root = _load(asset, data)
    path = ptr("/assets", asset)
    if not isinstance(root, dict):
        raise bad(path, f"{asset}: expected a JSON object envelope, got {type(root).__name__}")
    sets = root.get("GermlineSet")
    if not isinstance(sets, list) or len(sets) != 1:
        found = "absent" if "GermlineSet" not in root else (
            f"{len(sets)} sets" if isinstance(sets, list) else type(sets).__name__)
        raise bad(ptr(path, "GermlineSet"), f"{asset}: expected a GermlineSet array holding exactly one germline "
                                            f"set, got {found}")
    germline_set = sets[0]
    where = path + SET_POINTER
    if not isinstance(germline_set, dict):
        raise bad(where, f"{asset}: the germline set is not an object ({type(germline_set).__name__})")
    if not isinstance(germline_set.get("allele_descriptions"), list):
        raise bad(ptr(where, "allele_descriptions"), f"{asset}: allele_descriptions must be an array of records")
    species = germline_set.get("species")
    if not isinstance(species, dict):
        raise bad(ptr(where, "species"), f"{asset}: the set names no species object")
    taxon_id = taxonomy_id(species.get("id"))
    if taxon_id is None:
        raise bad(ptr(ptr(where, "species"), "id"), f"{asset}: species id {show(species.get('id'))} is not "
                                                    f"{TAXON_PREFIX}<taxonomy id>")
    try:
        name = check_text(species.get("label"), ptr(ptr(where, "species"), "label"), max_len=MAX_NAME,
                          multiline=False)
    except ReferenceSourceError as exc:
        raise bad(exc.path, f"{asset}: species label: {exc.message}", exc.kind) from None
    locus = germline_set.get("locus")
    if type(locus) is not str or locus not in LOCUS_CHAIN_TYPES:
        raise bad(ptr(where, "locus"), f"{asset}: locus {show(locus)} is not one of {sorted(LOCUS_CHAIN_TYPES)}")
    release = germline_set.get("release_version")
    if release is not None and (type(release) not in (int, float) or not math.isfinite(release)):
        raise bad(ptr(where, "release_version"), f"{asset}: release_version {show(release)} is not a number")
    identity = SourceIdentity(Taxon(name, taxon_id), locus, LOCUS_CHAIN_TYPES[locus])
    return _Envelope(asset=asset, root=root, germline_set=germline_set, identity=identity,
                     set_id=_optional_text(germline_set, "germline_set_id", where),
                     set_ref=_optional_text(germline_set, "germline_set_ref", where),
                     set_name=_optional_text(germline_set, "germline_set_name", where),
                     release=release)


def _mismatch(envelope: _Envelope, name: str, message: str) -> ReferenceSourceError:
    return bad(envelope.at(SET_POINTER + "/" + name), f"{envelope.asset}: {message}", "identity_mismatch")


def _check_request(envelope: _Envelope, set_id: str, release: int) -> None:
    """The payload's own set and release identifiers, where stated, must name the request."""
    match = _SET_ID_RE.fullmatch(set_id)
    taxon = envelope.identity.taxon.taxonomy_id
    if match is None or int(match.group(1)) != taxon:
        raise _mismatch(envelope, "species/id", f"the payload's taxon {TAXON_PREFIX}{taxon} is not the taxonomy "
                                                f"prefix of the requested set id {set_id!r}")
    name = match.group(2)
    if envelope.release is not None and envelope.release != release:
        raise _mismatch(envelope, "release_version", f"the payload is release {show(envelope.release)}, not the "
                                                     f"requested release {release}")
    if envelope.set_ref is not None and envelope.set_ref.rsplit(":", 1)[-1] != f"{name}.{release}":
        raise _mismatch(envelope, "germline_set_ref", f"the payload is set {envelope.set_ref!r}, not the requested "
                                                      f"set {name!r} release {release}")
    if envelope.set_name is not None:
        names = (name, name + EXTENDED_SET_SUFFIX) if envelope.asset == "airr_ex" else (name,)
        if envelope.set_name not in names:
            raise _mismatch(envelope, "germline_set_name", f"the payload names set {envelope.set_name!r}, not the "
                                                           f"requested set {name!r}")


def _check_agreement(first: _Envelope, other: _Envelope) -> None:
    """Every record export of one lock describes the same set."""
    pairs = (("species/id", first.identity.taxon.taxonomy_id, other.identity.taxon.taxonomy_id),
             ("species/label", first.identity.taxon.scientific_name, other.identity.taxon.scientific_name),
             ("locus", first.identity.locus, other.identity.locus),
             ("germline_set_id", first.set_id, other.set_id),
             ("germline_set_ref", first.set_ref, other.set_ref),
             ("release_version", first.release, other.release))
    for name, mine, theirs in pairs:
        if mine is not None and theirs is not None and mine != theirs:
            raise _mismatch(other, name, f"{name} {show(theirs)} disagrees with {show(mine)} in {first.asset}; "
                                          "every record export must describe the same set")


def _identity(envelopes: list, set_id: str, release: int) -> SourceIdentity:
    if not envelopes:
        raise bad("/assets", "no record export (airr or airr_ex) to read the set identity from", "missing_asset")
    for envelope in envelopes:
        _check_request(envelope, set_id, release)
    for other in envelopes[1:]:
        _check_agreement(envelopes[0], other)
    return envelopes[0].identity


def identify(*, set_id: str, release: int, assets: Mapping[str, bytes]) -> SourceIdentity:
    """The set identity (taxon, locus, chain type) read from the record exports' envelopes.

    ``assets`` maps formats to verified bytes; only the record exports (``airr``,
    ``airr_ex``) carry an identity and every one present must agree.  The payload's own
    taxon, set and release identifiers must name ``set_id`` (``<taxonomy id>.<set name>``)
    and ``release`` wherever the payload states them (``identity_mismatch`` otherwise).
    A malformed envelope is ``malformed`` (``resource_limit`` beyond the parser's bounds);
    no record export is ``missing_asset``.  Error paths are ``/assets/<format>/<pointer>``.
    """
    if type(set_id) is not str or type(release) is not int:
        raise TypeError(f"set_id must be a str and release an int, got {type(set_id).__name__} and "
                        f"{type(release).__name__}")
    _require_assets(assets)
    envelopes = [_envelope(fmt, assets[fmt]) for fmt in RECORD_EXPORT_FORMATS if fmt in assets]
    return _identity(envelopes, set_id, release)


# ── Field registries ────────────────────────────────────────────────
#
# The fields of the pinned producer's exports.  A field outside these registries is an
# unknown provider field: it is kept as an ``unsupported`` observation and recorded.

TOP_FIELDS = ("Info", "GermlineSet")
SET_FIELDS = ("germline_set_id", "author", "lab_name", "lab_address", "acknowledgements", "release_version",
              "release_description", "release_date", "germline_set_name", "germline_set_ref", "pub_ids",
              "species", "species_subgroup", "species_subgroup_type", "locus", "allele_descriptions",
              "curation")
RECORD_FIELDS = (
    "allele_description_id", "allele_description_ref", "maintainer", "acknowledgements", "lab_address",
    "release_version", "release_date", "release_description", "label", "sequence", "coding_sequence", "aliases",
    "locus", "chromosome", "sequence_type", "functional", "inference_type", "species", "species_subgroup",
    "species_subgroup_type", "status", "subgroup_designation", "gene_designation", "allele_designation",
    "allele_similarity_cluster_designation", "allele_similarity_cluster_member_id", "j_codon_frame",
    "gene_start", "gene_end", "utr_5_prime_start", "utr_5_prime_end", "leader_1_start", "leader_1_end",
    "leader_2_start", "leader_2_end", "v_rs_start", "v_rs_end", "d_rs_3_prime_start", "d_rs_3_prime_end",
    "d_rs_5_prime_start", "d_rs_5_prime_end", "j_cdr3_end", "j_rs_start", "j_rs_end", "j_donor_splice",
    "v_gene_delineations", "unrearranged_support", "rearranged_support", "paralogs", "curation",
    "curational_tags",
)
DELINEATION_FIELDS = ("sequence_delineation_id", "delineation_scheme", "unaligned_sequence", "aligned_sequence",
                      "fwr1_start", "fwr1_end", "cdr1_start", "cdr1_end", "fwr2_start", "fwr2_end", "cdr2_start",
                      "cdr2_end", "fwr3_start", "fwr3_end", "cdr3_start", "alignment_labels")
SUPPORT_FIELDS = {
    "unrearranged_support": ("sequence_id", "sequence", "curation", "repository_name", "repository_ref",
                             "patch_no", "gff_seqid", "gff_start", "gff_end", "strand"),
    "rearranged_support": ("sequence_id", "sequence", "derivation", "observation_type", "curation",
                           "repository_name", "repository_ref", "deposited_version", "sequence_start",
                           "sequence_end"),
}
SUPPORT_KIND = {"unrearranged_support": "unrearranged", "rearranged_support": "rearranged"}
SPECIES_FIELDS = ("id", "label")
RELATIONSHIP_FIELDS = {"aliases": "alias", "paralogs": "paralog"}
# Fields the inventory normalizes (identity facts); every other known field is kept uninterpreted.
INTERPRETED_SET_FIELDS = ("germline_set_id", "release_version", "germline_set_name", "germline_set_ref",
                          "species", "locus")
INTERPRETED_RECORD_FIELDS = ("allele_description_id", "allele_description_ref", "label", "coding_sequence",
                             "aliases", "locus", "sequence_type", "species", "unrearranged_support",
                             "rearranged_support", "paralogs")
# A record without a valid one of these has no identity or no candidate.
REQUIRED_RECORD_FIELDS = ("label", "sequence_type", "coding_sequence")
# Export-level findings cite the export's set identifier observation.
EXPORT_ANCHOR = SET_POINTER + "/germline_set_ref"

_SEQUENCE_RE = re.compile(r"[ACGT]+")
_JSON_TYPES = {dict: "object", list: "array", str: "string", int: "number", float: "number", bool: "Boolean",
               type(None): "null"}
_FASTA_FORMATS = ("gapped", "gapped_ex", "ungapped", "ungapped_ex")


def _json_type(value: Any) -> str:
    return _JSON_TYPES.get(type(value), type(value).__name__)


def _text_or_none(value: Any, max_len: int) -> Optional[str]:
    """``value`` if it is valid one-line text of at most ``max_len`` characters, else ``None``."""
    try:
        return check_text(value, "", max_len=max_len, multiline=False)
    except ReferenceSourceError:
        return None


# ── Findings ────────────────────────────────────────────────────────
#
# code -> (severity, category, affected capabilities, meaning).

FINDING_CODES = {
    # Records and fields.
    "malformed_record": (
        "error", "invalid_source_syntax", (),
        "A record that is not an object, or has no valid label or sequence_type; it is kept with its problem "
        "and cannot be selected."),
    "unsupported_record": (
        "warning", "unsupported_interpretation", (),
        "A record whose segment is not V, D or J, or not a segment of the set's chain type; it is kept and "
        "cannot be selected."),
    "record_without_sequence_candidate": (
        "warning", "invalid_source_syntax", (),
        "A record whose coding_sequence is absent, null or not [ACGT] text within the length limit; it has no "
        "sequence candidate."),
    "invalid_field_value": (
        "warning", "invalid_source_syntax", (),
        "A field the inventory normalizes holds a value of the wrong form; the value is kept and not used."),
    "record_identity_differs_from_set": (
        "warning", "contradictory_evidence", (),
        "A record's taxon or locus differs from its set's."),
    "unknown_provider_field": (
        "info", "informational", (),
        "A provider field this adapter does not know; its values are kept as unsupported observations (or, "
        "where no locator can name it, in the locked bytes)."),
    # Identities within one export.
    "duplicate_provider_id": (
        "warning", "contradictory_evidence", (),
        "Records of one export share a provider id; they stay separate records."),
    "duplicate_record_occurrence": (
        "info", "informational", (),
        "An identity supplied more than once with the same coding sequence in one export."),
    "same_label_sequence_alternatives": (
        "warning", "contradictory_evidence", (),
        "An identity supplied with different coding sequences in one export; using it needs an explicit "
        "choice of one."),
    "identical_sequence_distinct_identities": (
        "info", "informational", (),
        "Distinct identities supply the same coding sequence in one export; each keeps only its own records "
        "and evidence."),
    "alias_names_explicit_record": (
        "warning", "contradictory_evidence", (),
        "A declared alias names another record's label in the same export, so it is not a plain alias."),
    "mentioned_identity_segment_conflict": (
        "warning", "contradictory_evidence", (),
        "A relationship-only name is named by records of different segments; no mentioned identity is formed."),
    # Endpoint comparison between the two record exports (no export is preferred).
    "endpoint_comparison_incomplete": (
        "warning", "unavailable_evidence", (),
        "Only one record export is locked, so the two exports could not be compared."),
    "mentioned_identity_without_record": (
        "info", "unavailable_evidence", (),
        "A name given only as a paralog has no record in any locked record export."),
    "export_label_absent": (
        "info", "informational", (),
        "An identity has records in one export only; the other may name it as a paralog or alias."),
    "export_coding_sequence_differs": (
        "warning", "contradictory_evidence", ("original_coding_ends",),
        "The two exports supply different coding sequences for one identity."),
    "export_record_sequence_differs": (
        "warning", "contradictory_evidence", (),
        "The two exports' records of one identity hold different sequence fields."),
    "export_support_differs": (
        "info", "informational", (),
        "The two exports' records of one identity embed different support entries."),
    "export_provider_ref_differs": (
        "info", "informational", (),
        "The two exports' records of one identity have provider refs that differ other than by the extended "
        "export's suffix."),
    "export_annotation_differs": (
        "info", "informational", (),
        "The two exports' records of one identity differ in other fields, listed by name."),
    "export_provider_ref_suffix": (
        "info", "informational", (),
        "How many record pairs' provider refs differ only by the extended export's suffix, and how many refs "
        "of each export carry it."),
}
# Findings that make up the endpoint comparison.
COMPARISON_CODES = ("endpoint_comparison_incomplete", "mentioned_identity_without_record", "export_label_absent",
                    "export_coding_sequence_differs", "export_record_sequence_differs", "export_support_differs",
                    "export_provider_ref_differs", "export_annotation_differs", "export_provider_ref_suffix")
# Record fields compared on their own; every other field is compared as an annotation.
_COMPARED_APART = ("label", "sequence_type", "allele_description_ref", "coding_sequence", "sequence",
                   "unrearranged_support", "rearranged_support")


@dataclass
class _Entry:
    """One record occurrence while the inventory is built."""

    record: SourceRecord
    item: Any
    candidate: Optional[SequenceCandidate] = None
    relationships: dict = field(default_factory=lambda: {"alias": [], "paralog": []})


class _Builder:
    """Accumulates the items of one inventory, in any order (``assemble`` sorts them)."""

    def __init__(self, lock: SourceLock, identity: SourceIdentity):
        self.lock = lock
        self.lock_hash = lock_hash(lock)
        self.identity = identity
        self.records: list = []
        self.candidates: list = []
        self.observations: list = []
        self.support: list = []
        self.mentioned: list = []
        self.findings: list = []
        self.entries: dict = {}          # asset -> [_Entry]
        self.unknown: dict = {}          # (asset, field pattern) -> [locator]
        self.observed: set = set()       # locators of every observation

    # -- items ---------------------------------------------------------

    def observe(self, asset: str, pointer: str, record: Optional[Locator], value: Any, interpretation: str,
                presence: str = "present") -> Locator:
        if len(self.observations) >= MAX_OBSERVATIONS:
            raise bad("/assets", f"the inventory would hold more than {MAX_OBSERVATIONS} observations",
                      "resource_limit")
        locator = Locator(asset, pointer)
        try:
            self.observations.append(Observation(locator=locator, record=record, presence=presence,
                                                 value=value if presence == "present" else None,
                                                 interpretation=interpretation))
        except ReferenceSourceError as exc:
            raise bad(ptr("/assets", asset) + pointer, f"{asset}: {exc.message}", exc.kind) from None
        self.observed.add(locator)
        return locator

    def finding(self, code: str, *, subjects=(), locators=(), observations: dict) -> None:
        """Add a finding; one naming more than ``MAX_REFS`` subjects or locators is split into parts."""
        severity, category, capabilities, _ = FINDING_CODES[code]
        subjects = sorted(set(subjects), key=Locator.sort_key)
        locators = sorted(set(locators), key=Locator.sort_key)
        parts = max(1, -(-len(subjects) // MAX_REFS), -(-len(locators) // MAX_REFS))
        for part in range(parts):
            details = dict(observations)
            if parts > 1:
                details.update(part=part + 1, parts=parts)
            window = slice(part * MAX_REFS, (part + 1) * MAX_REFS)
            self.findings.append(Finding(source_lock_hash=self.lock_hash, code=code, severity=severity,
                                         category=category, subjects=tuple(subjects[window]),
                                         locators=tuple(locators[window]), observations=details,
                                         capabilities=tuple(capabilities)))

    def remember_unknown(self, asset: str, pattern: str, key: str, locator: Locator) -> None:
        """Record an unknown provider field ``key`` (at ``pattern``) kept at ``locator``."""
        self.unknown.setdefault((asset, ptr(pattern, key)), []).append(locator)

    def unknown_field(self, asset: str, pattern: str, base: str, key: str, fallback: Locator,
                      record: Optional[Locator], value: Any) -> None:
        """Keep unknown provider field ``key`` of the object at ``base`` as an ``unsupported`` observation.

        ``pattern`` names the field's place for every record (``*`` for an index).  A key no
        pointer can name (too long, or with control characters) is recorded at ``fallback``;
        its value stays in the locked bytes.
        """
        pointer = ptr(base, key)
        try:
            Locator(asset, pointer)
        except ReferenceSourceError:
            locator = fallback
        else:
            locator = self.observe(asset, pointer, record, value, "unsupported")
        self.remember_unknown(asset, pattern, key, locator)

    # -- one export ----------------------------------------------------

    def export(self, envelope: _Envelope) -> None:
        asset = envelope.asset
        anchor = Locator(asset, EXPORT_ANCHOR)
        for key, value in envelope.root.items():
            if key == "GermlineSet":
                continue
            if key in TOP_FIELDS:
                self.observe(asset, ptr("", key), None, value, "preserved")
            else:
                self.unknown_field(asset, "", "", key, anchor, None, value)
        if "Info" not in envelope.root:
            self.observe(asset, "/Info", None, None, "preserved", "absent")
        germline_set = envelope.germline_set
        for key, value in germline_set.items():
            if key == "allele_descriptions":
                continue
            pointer = ptr(SET_POINTER, key)
            if key in SET_FIELDS:
                self.observe(asset, pointer, None, value,
                             "interpreted" if key in INTERPRETED_SET_FIELDS else "preserved")
            else:
                self.unknown_field(asset, SET_POINTER, SET_POINTER, key, anchor, None, value)
        for key in SET_FIELDS:
            if key not in germline_set and key != "allele_descriptions":
                self.observe(asset, ptr(SET_POINTER, key), None, None,
                             "interpreted" if key in INTERPRETED_SET_FIELDS else "preserved", "absent")
        for key, value in germline_set["species"].items():
            if key not in SPECIES_FIELDS:
                self.remember_unknown(asset, SET_POINTER + "/species", key, Locator(asset, SET_POINTER + "/species"))
        self.entries[asset] = [self.record(asset, index, item)
                               for index, item in enumerate(germline_set["allele_descriptions"])]

    def record(self, asset: str, index: int, item: Any) -> _Entry:
        locator = Locator(asset, ptr(RECORDS_POINTER, index))
        if not isinstance(item, dict):
            problem = f"the record is a JSON {_json_type(item)}, not an object"
            record = SourceRecord(locator=locator, status="malformed", problem=problem, label=None,
                                  provider_id=None, provider_ref=None, segment=None, taxonomy_id=None, locus=None)
            self.records.append(record)
            self.finding("malformed_record", subjects=[locator], observations={"problem": problem})
            return _Entry(record, item)
        fields = {}
        for key, value in item.items():
            if key in RECORD_FIELDS:
                fields[key] = self.observe(asset, ptr(locator.pointer, key), locator, value,
                                           self._interpretation(key, value))
            else:
                self.unknown_field(asset, RECORDS_POINTER + "/*", locator.pointer, key, locator, locator, value)
        for key in RECORD_FIELDS:
            if key not in item:
                fields[key] = self.observe(asset, ptr(locator.pointer, key), locator, None,
                                           self._interpretation(key, None, absent=True), "absent")
        record = self._source_record(locator, item, fields)
        self.records.append(record)
        entry = _Entry(record, item)
        self._candidate(entry, fields["coding_sequence"])
        self._nested(entry, fields)
        return entry

    @staticmethod
    def _interpretation(key: str, value: Any, absent: bool = False) -> str:
        if key not in INTERPRETED_RECORD_FIELDS:
            return "preserved"
        if absent:
            return "invalid" if key in REQUIRED_RECORD_FIELDS else "interpreted"
        if key == "label":
            return "invalid" if _label_problem(value) else "interpreted"
        if key == "sequence_type":
            return "interpreted" if type(value) is str and value else "invalid"
        if key == "coding_sequence":
            return "interpreted" if _sequence_problem(value) is None else "invalid"
        if value is None:
            return "interpreted"
        if key in ("allele_description_id", "allele_description_ref"):
            return "interpreted" if _text_or_none(value, MAX_NAME) is not None else "invalid"
        if key == "locus":
            return "interpreted" if _text_or_none(value, MAX_LOCUS) is not None else "invalid"
        if key == "species":
            return "interpreted" if isinstance(value, dict) and taxonomy_id(value.get("id")) else "invalid"
        return "interpreted" if isinstance(value, list) else "invalid"   # relationships and support

    def _invalid(self, record: Locator, locator: Locator, field_name: str, value: Any, expected: str) -> None:
        self.finding("invalid_field_value", subjects=[record], locators=[locator],
                     observations={"field": field_name, "found": _json_type(value), "expected": expected})

    def _status(self, item: dict) -> tuple:
        """``(status, problem, segment)`` of a record object."""
        label, sequence_type = item.get("label"), item.get("sequence_type")
        segment = sequence_type if type(sequence_type) is str and sequence_type in SEGMENTS else None
        if "label" not in item:
            return "malformed", "the record has no label", segment
        problem = _label_problem(label)
        if problem:
            return "malformed", problem, segment
        if "sequence_type" not in item:
            return "malformed", "the record has no sequence_type", segment
        if type(sequence_type) is not str or not sequence_type:
            return "malformed", f"sequence_type {show(sequence_type)} is not a segment name", segment
        if segment is None:
            return "unsupported", f"sequence_type {show(sequence_type)} is not a V, D or J segment", segment
        if segment not in self.identity.segments:
            return "unsupported", (f"a {segment} record cannot belong to a {self.identity.locus} "
                                   f"({self.identity.chain_type}) set"), segment
        return "parsed", None, segment

    def _source_record(self, locator: Locator, item: dict, fields: dict) -> SourceRecord:
        status, problem, segment = self._status(item)
        label = item.get("label")
        if status != "parsed":
            try:   # whatever is readable of the label is kept
                label = check_text(label, "", max_len=MAX_TEXT)
            except ReferenceSourceError:
                label = None
        species = item.get("species")
        taxon = taxonomy_id(species.get("id")) if isinstance(species, dict) else None
        locus = _text_or_none(item.get("locus"), MAX_LOCUS)
        record = SourceRecord(locator=locator, status=status, problem=problem, label=label,
                              provider_id=_text_or_none(item.get("allele_description_id"), MAX_NAME),
                              provider_ref=_text_or_none(item.get("allele_description_ref"), MAX_NAME),
                              segment=segment, taxonomy_id=taxon, locus=locus)
        if status != "parsed":
            self.finding(f"{status}_record", subjects=[locator], observations={"problem": problem})
        for key, expected in (("allele_description_id", "text of at most 256 characters"),
                              ("allele_description_ref", "text of at most 256 characters"),
                              ("locus", f"text of at most {MAX_LOCUS} characters"),
                              ("species", f"an object with id {TAXON_PREFIX}<taxonomy id>")):
            if item.get(key) is not None and self._interpretation(key, item[key]) == "invalid":
                self._invalid(locator, fields[key], key, item[key], expected)
        for key in species if isinstance(species, dict) else ():
            if key not in SPECIES_FIELDS:
                self.remember_unknown(locator.asset, RECORDS_POINTER + "/*/species", key, fields["species"])
        differs, cited = {}, []
        if taxon is not None and taxon != self.identity.taxon.taxonomy_id:
            differs["taxonomy_id"] = {"record": taxon, "set": self.identity.taxon.taxonomy_id}
            cited.append(fields["species"])
        if locus is not None and locus != self.identity.locus:
            differs["locus"] = {"record": locus, "set": self.identity.locus}
            cited.append(fields["locus"])
        if differs:
            self.finding("record_identity_differs_from_set", subjects=[locator], locators=cited,
                         observations=differs)
        return record

    def _candidate(self, entry: _Entry, locator: Locator) -> None:
        record = entry.record
        value = entry.item.get("coding_sequence")
        problem = _sequence_problem(value) if "coding_sequence" in entry.item else "absent"
        if problem is None:
            entry.candidate = SequenceCandidate(locator=locator, record=record.locator, sequence=value,
                                                sha256=sequence_digest(value))
            self.candidates.append(entry.candidate)
        else:
            self.finding("record_without_sequence_candidate", subjects=[record.locator], locators=[locator],
                         observations={"label": record.label, "reason": problem})

    def _nested(self, entry: _Entry, fields: dict) -> None:
        record, item, asset = entry.record.locator, entry.item, entry.record.locator.asset
        for key, relationship in RELATIONSHIP_FIELDS.items():
            value = item.get(key)
            if value is not None and not isinstance(value, list):
                self._invalid(record, fields[key], key, value, "an array of labels or null")
            for j, name in enumerate(value if isinstance(value, list) else ()):
                valid = _label_problem(name) is None
                locator = self.observe(asset, ptr(fields[key].pointer, j), record, name,
                                       "interpreted" if valid else "invalid")
                if valid:
                    entry.relationships[relationship].append((name, locator))
                else:
                    self._invalid(record, locator, key, name, "a label")
        value = item.get("v_gene_delineations")
        for j, delineation in enumerate(value if isinstance(value, list) else ()):
            base = ptr(fields["v_gene_delineations"].pointer, j)
            if not isinstance(delineation, dict):
                self.observe(asset, base, record, delineation, "preserved")
                continue
            for key, part in delineation.items():
                if key in DELINEATION_FIELDS:
                    self.observe(asset, ptr(base, key), record, part, "preserved")
                else:
                    self.unknown_field(asset, RECORDS_POINTER + "/*/v_gene_delineations/*", base, key, record,
                                       record, part)
            for key in DELINEATION_FIELDS:
                if key not in delineation:
                    self.observe(asset, ptr(base, key), record, None, "preserved", "absent")
        for key, kind in SUPPORT_KIND.items():
            value = item.get(key)
            if value is not None and not isinstance(value, list):
                self._invalid(record, fields[key], key, value, "an array of support entries or null")
            for j, support in enumerate(value if isinstance(value, list) else ()):
                valid = isinstance(support, dict)
                locator = self.observe(asset, ptr(fields[key].pointer, j), record, support,
                                       "interpreted" if valid else "invalid")
                if not valid:
                    self._invalid(record, locator, key, support, "a support entry object")
                    continue
                for part in support:
                    if part not in SUPPORT_FIELDS[key]:
                        self.remember_unknown(asset, RECORDS_POINTER + f"/*/{key}/*", part, locator)
                stated = {}
                for part in ("repository_name", "repository_ref"):
                    stated[part] = _text_or_none(support.get(part), MAX_NAME)
                    if support.get(part) is not None and stated[part] is None:
                        self._invalid(record, locator, f"{key}/{part}", support[part],
                                      "text of at most 256 characters")
                self.support.append(SupportReference(locator=locator, record=record, kind=kind,
                                                     repository=stated["repository_name"],
                                                     reference=stated["repository_ref"],
                                                     status="embedded_unretrieved"))

    # -- relationships and identities within one export ----------------

    def relationships(self, asset: str) -> None:
        entries = self.entries[asset]
        labels = {e.record.label for e in entries if e.record.label is not None}
        records_by_label: dict = {}
        for entry in entries:
            if entry.record.status == "parsed":
                records_by_label.setdefault(entry.record.label, []).append(entry.record.locator)
        mentioned: dict = {}
        for entry in entries:
            if entry.record.status != "parsed":
                continue
            own = entry.record.label
            for name, locator in entry.relationships["alias"]:
                if name != own and name in labels:
                    self.finding("alias_names_explicit_record",
                                 subjects=[entry.record.locator] + records_by_label.get(name, []),
                                 locators=[locator], observations={"label": own, "alias": name})
            for name, locator in entry.relationships["paralog"]:
                if name not in labels:
                    mentioned.setdefault(name, []).append((entry.record, locator))
        for name, mentions in sorted(mentioned.items()):
            segments = sorted({record.segment for record, _ in mentions})
            if len(segments) != 1:
                self.finding("mentioned_identity_segment_conflict",
                             subjects=[record.locator for record, _ in mentions],
                             locators=[locator for _, locator in mentions],
                             observations={"label": name, "segments": segments})
                continue
            try:
                self.mentioned.append(MentionedIdentity(
                    asset=asset, segment=segments[0], label=name,
                    mentions=tuple(sorted((Mention(record=record.locator, locator=locator, relationship="paralog")
                                           for record, locator in mentions), key=lambda x: x.locator.sort_key()))))
            except ReferenceSourceError as exc:
                raise bad(ptr("/assets", asset), f"{asset}: paralog name {show(name)}: {exc.message}",
                          exc.kind) from None

    def identities(self, asset: str) -> None:
        entries = [e for e in self.entries[asset] if e.record.status == "parsed"]
        by_provider_id: dict = {}
        for entry in self.entries[asset]:
            if entry.record.provider_id is not None:
                by_provider_id.setdefault(entry.record.provider_id, []).append(entry)
        for provider_id, group in sorted(by_provider_id.items()):
            if len(group) > 1:
                self.finding("duplicate_provider_id", subjects=[e.record.locator for e in group],
                             locators=[Locator(asset, ptr(e.record.locator.pointer, "allele_description_id"))
                                       for e in group],
                             observations={"provider_id": provider_id, "records": len(group)})
        by_identity: dict = {}
        by_sequence: dict = {}
        for entry in entries:
            by_identity.setdefault((entry.record.segment, entry.record.label), []).append(entry)
            if entry.candidate is not None:
                by_sequence.setdefault(entry.candidate.sha256, []).append(entry)
        for (segment, label), group in sorted(by_identity.items()):
            variants: dict = {}
            for entry in group:
                if entry.candidate is not None:
                    variants.setdefault(entry.candidate.sha256, []).append(entry)
            for sha256, same in sorted(variants.items()):
                if len(same) > 1:
                    self.finding("duplicate_record_occurrence", subjects=_subjects(same),
                                 observations={"label": label, "segment": segment, "sha256": sha256,
                                               "occurrences": len(same)})
            if len(variants) > 1:
                self.finding("same_label_sequence_alternatives", subjects=_subjects(group),
                             observations={"label": label, "segment": segment, "alternatives": [
                                 {"sha256": sha256, "length": len(same[0].candidate.sequence),
                                  "records": len(same)} for sha256, same in sorted(variants.items())]})
        for sha256, group in sorted(by_sequence.items()):
            identities = sorted({(e.record.segment, e.record.label) for e in group})
            if len(identities) > 1:
                self.finding("identical_sequence_distinct_identities", subjects=_subjects(group),
                             observations={"sha256": sha256, "length": len(group[0].candidate.sequence),
                                           "identities": [{"segment": s, "label": l} for s, l in identities]})

    # -- endpoint comparison ------------------------------------------

    def compare(self) -> None:
        """Compare the record exports label by label, preferring neither."""
        present = [a for a in RECORD_EXPORT_FORMATS if a in self.entries]
        if len(present) < len(RECORD_EXPORT_FORMATS):
            missing = [{"format": x.format, "reason": x.reason} for x in self.lock.missing_assets
                       if x.format in RECORD_EXPORT_FORMATS]
            self.finding("endpoint_comparison_incomplete", locators=[Locator(a, EXPORT_ANCHOR) for a in present],
                         observations={"exports": present, "missing": missing})
        else:
            self._compare_exports()
        identities = {key for a in present for key in _parsed_identities(self.entries[a])}
        for identity in self.mentioned:
            if (identity.segment, identity.label) not in identities:
                self.finding("mentioned_identity_without_record", subjects=[x.record for x in identity.mentions],
                             locators=[x.locator for x in identity.mentions],
                             observations={"asset": identity.asset, "label": identity.label,
                                           "segment": identity.segment})

    def _compare_exports(self) -> None:
        first, second = RECORD_EXPORT_FORMATS
        sides = {a: _parsed_identities(self.entries[a]) for a in RECORD_EXPORT_FORMATS}
        refs = {"pairs": 0, "suffix_only": 0, "identical": 0, "other": 0}
        for key in sorted(sides[first].keys() | sides[second].keys()):
            if key not in sides[first] or key not in sides[second]:
                self._label_absent(key, sides)
                continue
            mine, theirs = sides[first][key], sides[second][key]
            self._compare_coding(key, mine, theirs)
            if len(mine) == len(theirs) == 1:
                refs["pairs"] += 1
                refs[self._compare_pair(key, mine[0], theirs[0])] += 1
        exports = {a: {"records": len(self.entries[a]),
                       "refs_with_suffix": sum(1 for e in self.entries[a]
                                               if (e.record.provider_ref or "").endswith(EXT_SUFFIX))}
                   for a in RECORD_EXPORT_FORMATS}
        self.finding("export_provider_ref_suffix", locators=[Locator(a, EXPORT_ANCHOR) for a in RECORD_EXPORT_FORMATS],
                     observations={"suffix": EXT_SUFFIX, **refs, **exports})

    def _label_absent(self, key: tuple, sides: dict) -> None:
        segment, label = key
        present = next(a for a in RECORD_EXPORT_FORMATS if key in sides[a])
        absent = next(a for a in RECORD_EXPORT_FORMATS if key not in sides[a])
        named, locators = [], []
        for entry in self.entries[absent]:
            if entry.record.status != "parsed":
                continue
            for relationship, names in entry.relationships.items():
                for name, locator in names:
                    if name == label:
                        named.append({"relationship": relationship, "by": entry.record.label})
                        locators.append(locator)
        self.finding("export_label_absent", subjects=_subjects(sides[present][key]), locators=locators,
                     observations={"label": label, "segment": segment, "present_in": present, "absent_from": absent,
                                   "named_in_absent_export": sorted(named, key=lambda x: (x["by"], x["relationship"]))})

    def _compare_coding(self, key: tuple, mine: list, theirs: list) -> None:
        first, second = RECORD_EXPORT_FORMATS
        if _field_values(mine, "coding_sequence") == _field_values(theirs, "coding_sequence"):
            return
        relation, extra = "not_comparable", None
        if len(mine) == len(theirs) == 1 and mine[0].candidate and theirs[0].candidate:
            relation, extra = _extension(mine[0].candidate.sequence, theirs[0].candidate.sequence, first, second)
        segment, label = key
        self.finding("export_coding_sequence_differs", subjects=_subjects(mine + theirs),
                     locators=[self._field(e, "coding_sequence") for e in mine + theirs],
                     observations={"label": label, "segment": segment, first: _coding_side(mine),
                                   second: _coding_side(theirs), "relation": relation, "extra": extra})

    def _field(self, entry: _Entry, key: str) -> Locator:
        return Locator(entry.record.locator.asset, ptr(entry.record.locator.pointer, key))

    def _compare_pair(self, key: tuple, mine: _Entry, theirs: _Entry) -> str:
        """Compare one record of each export; return the relation of their provider refs."""
        first, second = RECORD_EXPORT_FORMATS
        segment, label = key
        pair, base = [mine, theirs], {"label": label, "segment": segment}
        subjects = _subjects(pair)
        if _field_values([mine], "sequence") != _field_values([theirs], "sequence"):
            self.finding("export_record_sequence_differs", subjects=subjects,
                         locators=[self._field(e, "sequence") for e in pair],
                         observations={**base, first: _text_summary(mine.item, "sequence"),
                                       second: _text_summary(theirs.item, "sequence")})
        support = {k: {first: _count(mine.item, k), second: _count(theirs.item, k),
                       "equal": _field_values([mine], k) == _field_values([theirs], k)} for k in SUPPORT_KIND}
        if not all(v["equal"] for v in support.values()):
            self.finding("export_support_differs", subjects=subjects,
                         locators=[self._field(e, k) for k in SUPPORT_KIND if not support[k]["equal"] for e in pair],
                         observations={**base, **support})
        keys = (set(mine.item) | set(theirs.item)) - set(_COMPARED_APART)
        differing = sorted(k for k in keys if _field_values([mine], k) != _field_values([theirs], k))
        if differing:
            located = [self._field(e, k) for k in differing for e in pair]
            self.finding("export_annotation_differs", subjects=subjects,
                         locators=[x for x in located if x in self.observed],
                         observations={**base, "fields": differing})
        refs = (mine.record.provider_ref, theirs.record.provider_ref)
        if refs[0] == refs[1]:
            return "identical"
        if refs[0] is not None and refs[1] == refs[0] + EXT_SUFFIX:
            return "suffix_only"
        self.finding("export_provider_ref_differs", subjects=subjects,
                     locators=[self._field(e, "allele_description_ref") for e in pair],
                     observations={**base, first: refs[0], second: refs[1]})
        return "other"

    def fasta(self, asset: str, data: bytes) -> None:
        for index, (header, sequence) in enumerate(_fasta_records(asset, data)):
            self.observe(asset, f"/{index}", None, {"header": header, "sequence": sequence}, "preserved")

    def unknown_findings(self) -> None:
        for (asset, pattern), locators in sorted(self.unknown.items()):
            self.finding("unknown_provider_field", locators=locators,
                         observations={"asset": asset, "field": pattern, "occurrences": len(locators)})

    def inventory(self) -> SourceInventory:
        try:
            return SourceInventory.assemble(
                source_lock_hash=self.lock_hash, adapter=self.lock.adapter, records=self.records,
                candidates=self.candidates, mentioned_identities=self.mentioned, observations=self.observations,
                support=self.support, claims=(), mappings=(), findings=self.findings)
        except ReferenceSourceError as exc:   # an inventory bound such as the number of findings
            raise exc.under("/inventory") from None


def _label_problem(value: Any) -> Optional[str]:
    """``None`` for a valid runtime label, else the problem (as the runtime label rule states it)."""
    try:
        check_label(value, "")
    except ReferenceSourceError as exc:
        return exc.message
    return None


def _sequence_problem(value: Any) -> Optional[str]:
    if value is None:
        return "null"
    if type(value) is not str:
        return f"a JSON {_json_type(value)}, not text"
    if len(value) > MAX_SEQUENCE:
        return f"{len(value)} bases exceed the limit of {MAX_SEQUENCE}"
    if not _SEQUENCE_RE.fullmatch(value):
        return "not a sequence of A, C, G and T"
    return None


def _subjects(entries) -> list:
    """The records of ``entries`` and their candidates."""
    found = []
    for entry in entries:
        found.append(entry.record.locator)
        if entry.candidate is not None:
            found.append(entry.candidate.locator)
    return found


def _parsed_identities(entries) -> dict:
    """``(segment, label) -> [entry]`` of the parsed records among ``entries``, in document order."""
    found: dict = {}
    for entry in entries:
        if entry.record.status == "parsed":
            found.setdefault((entry.record.segment, entry.record.label), []).append(entry)
    return found


_ABSENT = object()


def _field_values(entries, key: str) -> list:
    """The distinct canonical values of field ``key`` of ``entries`` (absent apart from null), sorted."""
    return sorted({"" if e.item.get(key, _ABSENT) is _ABSENT else canonical_json(e.item[key]) for e in entries})


def _count(item: dict, key: str) -> Optional[int]:
    value = item.get(key)
    return len(value) if isinstance(value, list) else None


def _text_summary(item: dict, key: str):
    """Length and digest of a text field; its JSON type otherwise; null when absent."""
    if key not in item:
        return None
    value = item[key]
    if type(value) is str:
        return {"length": len(value), "sha256": text_digest(value)}
    return {"type": _json_type(value)}


def _coding_side(entries) -> dict:
    sequences = {e.candidate.sha256: len(e.candidate.sequence) for e in entries if e.candidate is not None}
    return {"records": len(entries),
            "sequences": [{"sha256": sha, "length": length} for sha, length in sorted(sequences.items())]}


def _extension(mine: str, theirs: str, first: str, second: str) -> tuple:
    """How two coding sequences relate: one extends the other at its 3' or 5' end, or neither."""
    for longer, shorter, name in ((theirs, mine, second), (mine, theirs, first)):
        if len(longer) > len(shorter):
            if longer.startswith(shorter):
                return f"{name}_extends_3prime", longer[len(shorter):]
            if longer.endswith(shorter):
                return f"{name}_extends_5prime", longer[:len(longer) - len(shorter)]
    return "different", None


def _fasta_records(asset: str, data: bytes) -> list:
    """``(header, sequence)`` per FASTA record, as supplied.

    The header is the header line after ``>``; the sequence is the record's other lines
    joined without their line terminators (``\\n`` or ``\\r\\n``), with case, gap and every
    other character kept.  One leading byte-order mark is ignored.  Blank lines before the
    first header are allowed; other text there is ``malformed``.
    """
    path = ptr("/assets", asset)
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise bad(path, f"{asset}: the export is not UTF-8 text ({exc})") from None
    if text.startswith("﻿"):
        text = text[1:]
    records, header, lines = [], None, []
    for number, line in enumerate(text.split("\n"), start=1):
        if line.endswith("\r"):
            line = line[:-1]
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(lines)))
            header, lines = line[1:], []
        elif header is not None:
            lines.append(line)
        elif line.strip():
            raise bad(path, f"{asset}: line {number} holds sequence data before the first '>' header")
    if header is not None:
        records.append((header, "".join(lines)))
    return records


# ── Inventory ───────────────────────────────────────────────────────


def _locked_bytes(lock: SourceLock, assets: Mapping) -> dict:
    """The bytes of every locked asset, each exactly the locked digest and size."""
    for name in assets:
        if name not in lock.asset_ids:
            raise bad(ptr("/assets", name), f"{show(name)} is not an asset of this lock (assets "
                                            f"{list(lock.asset_ids)})", "identity_mismatch")
    blobs = {}
    for asset in lock.assets:
        if asset.asset_id not in assets:
            raise bad(ptr("/assets", asset.asset_id), f"no bytes supplied for locked asset {asset.asset_id} "
                                                      f"({asset.sha256}, {asset.bytes} bytes)", "missing_asset")
        data = assets[asset.asset_id]
        digest = "sha256:" + hashlib.sha256(data).hexdigest()
        if (digest, len(data)) != (asset.sha256, asset.bytes):
            raise bad(ptr("/assets", asset.asset_id), f"{asset.asset_id}: the bytes are {digest} ({len(data)} "
                                                      f"bytes), not the locked {asset.sha256} ({asset.bytes} "
                                                      "bytes)", "digest_mismatch")
        blobs[asset.asset_id] = data
    return blobs


def build_inventory(lock: SourceLock, assets: Mapping[str, bytes]) -> SourceInventory:
    """The normalized inventory of a verified lock and its bytes (asset id -> bytes).

    The bytes must be exactly the locked assets (``missing_asset``, ``digest_mismatch``,
    ``identity_mismatch`` for bytes the lock does not hold), and the identity read from
    the record exports must be the lock's (``identity_mismatch`` at ``/identity``).  The
    record exports together may hold at most :data:`MAX_SOURCE_RECORDS` records
    (``resource_limit``).  A malformed envelope or FASTA export is ``malformed`` at
    ``/assets/<asset>...``; a malformed record never raises.  Data-terms statements are
    not read.  Claims and mappings are left empty: coordinates are not interpreted here.
    """
    require_type(lock, SourceLock, "lock")
    _require_assets(assets)
    if lock.adapter.id != ADAPTER_ID:
        raise bad("/adapter", f"adapter {lock.adapter.id!r} is not the OGRDB adapter {ADAPTER_ID!r}",
                  "unsupported_profile")
    blobs = _locked_bytes(lock, assets)
    envelopes = [_envelope(fmt, blobs[fmt]) for fmt in RECORD_EXPORT_FORMATS if fmt in blobs]
    identity = _identity(envelopes, lock.set_id, lock.release)
    if identity != lock.identity:
        raise bad("/identity", f"the record exports identify as {identity.to_dict()}, but the lock pins "
                               f"{lock.identity.to_dict()}", "identity_mismatch")
    total = sum(len(e.germline_set["allele_descriptions"]) for e in envelopes)
    if total > MAX_SOURCE_RECORDS:
        raise bad("/assets", f"the record exports hold {total} records together; the limit is "
                             f"{MAX_SOURCE_RECORDS}", "resource_limit")
    builder = _Builder(lock, identity)
    for envelope in envelopes:
        builder.export(envelope)
        builder.relationships(envelope.asset)
        builder.identities(envelope.asset)
    builder.compare()
    for asset in lock.assets:
        if asset.format in _FASTA_FORMATS:
            builder.fasta(asset.asset_id, blobs[asset.asset_id])
    builder.unknown_findings()
    return builder.inventory()


# ── Accounting ──────────────────────────────────────────────────────

_RELATIONSHIP_ITEM_RE = re.compile(r"/(aliases|paralogs)/(?:0|[1-9][0-9]*)")


@dataclass(frozen=True)
class ExportAccount:
    """Counts of one record export, each reconciling to its raw records.

    Record occurrences, sequence candidates, identities and relationship-only names are
    counted apart; executable alleles are chosen later by an explicit selection and are
    not counted here.  For every export these hold:

    * ``record_occurrences = parsed_records + malformed_records + unsupported_records``
    * ``record_occurrences = sequence_candidates + records_without_candidate``
    * ``parsed_records = identity_sequences + duplicate_occurrences + parsed_without_candidate``
    * ``identity_sequences = record_identities - identities_without_candidate + alternative_sequences``
    * ``named_identities = record_identities + mentioned_only_identities``
    * ``paralog_names = paralog_names_of_records + paralog_names_of_mentioned_only + paralog_names_unresolved``

    An identity is a parsed record's ``(segment, label)``.  ``identity_sequences`` counts
    distinct supplied coding sequences per identity; ``duplicate_occurrences`` the parsed
    records repeating one; ``alternative_sequences`` the sequences beyond the first of an
    identity supplied with several.  ``shared_sequences`` lists, per coding sequence that
    distinct identities share, their labels: they stay distinct identities.
    ``named_identities`` are the identities an export names at all: its own record
    identities plus the names it gives only as paralogs.  Declared aliases name their own
    record and add no identity.  ``segments`` is :meth:`SourceInventory.export_counts`.
    """

    asset: str
    record_occurrences: int
    parsed_records: int
    malformed_records: int
    unsupported_records: int
    sequence_candidates: int
    records_without_candidate: int
    record_identities: int
    identities_without_candidate: int
    identity_sequences: int
    duplicate_occurrences: int
    parsed_without_candidate: int
    alternative_identities: int
    alternative_sequences: int
    shared_sequences: tuple
    mentioned_only_identities: int
    named_identities: int
    alias_names: int
    paralog_names: int
    paralog_names_of_records: int
    paralog_names_of_mentioned_only: int
    paralog_names_unresolved: int
    support_references: int
    provider_refs_with_ext_suffix: int
    segments: Mapping

    def to_dict(self) -> dict:
        d = {k: getattr(self, k) for k in self.__dataclass_fields__}
        d["shared_sequences"] = [list(group) for group in self.shared_sequences]
        d["segments"] = {k: dict(v) for k, v in self.segments.items()}
        return d


@dataclass(frozen=True)
class EndpointComparison:
    """The comparison of the two record exports, preferring neither.

    ``complete`` exactly when both record exports are locked.  ``identities_in_both`` and
    ``only_in`` (asset -> labels with records in that export only) are filled for a
    complete comparison; with ``account.export(a).record_identities`` they satisfy
    ``record_identities = identities_in_both + len(only_in[a])`` for each export.
    ``findings`` are the ids of the comparison findings (:data:`COMPARISON_CODES`), sorted.
    """

    complete: bool
    exports: tuple
    missing: tuple
    identities_in_both: int
    only_in: Mapping
    findings: tuple

    def report(self) -> Comparison:
        """The build report's comparison section."""
        return Comparison(complete=self.complete, exports=self.exports, missing=self.missing,
                          findings=self.findings)

    def to_dict(self) -> dict:
        return {"complete": self.complete, "exports": list(self.exports), "missing": list(self.missing),
                "identities_in_both": self.identities_in_both,
                "only_in": {k: list(v) for k, v in self.only_in.items()}, "findings": list(self.findings)}


@dataclass(frozen=True)
class SourceAccount:
    """The accounting of one inventory: per record export, per FASTA export, and the comparison."""

    exports: tuple                  # ExportAccount per record export of the lock, in lock order
    fasta: tuple                    # (asset, records) per FASTA export of the lock, in lock order
    comparison: EndpointComparison

    def export(self, asset: str) -> ExportAccount:
        for account in self.exports:
            if account.asset == asset:
                return account
        raise KeyError(asset)

    def to_dict(self) -> dict:
        return {"exports": [a.to_dict() for a in self.exports], "fasta": [list(x) for x in self.fasta],
                "comparison": self.comparison.to_dict()}


def _export_account(inventory: SourceInventory, asset: str) -> ExportAccount:
    records = [r for r in inventory.records if r.locator.asset == asset]
    statuses = {s: sum(1 for r in records if r.status == s) for s in ("parsed", "malformed", "unsupported")}
    digests = {r.locator: [c.sha256 for c in inventory.candidates_of(r.locator)] for r in records}
    parsed = [r for r in records if r.status == "parsed"]
    identities: dict = {}
    by_sequence: dict = {}
    for record in parsed:
        identities.setdefault((record.segment, record.label), []).append(record)
        for sha256 in digests[record.locator]:
            by_sequence.setdefault(sha256, set()).add((record.segment, record.label))
    sequences = duplicates = without = alternatives = extra = bare = 0
    for group in identities.values():
        supplied = [d for r in group for d in digests[r.locator]]
        distinct = len(set(supplied))
        sequences += distinct
        duplicates += len(supplied) - distinct
        without += sum(1 for r in group if not digests[r.locator])
        bare += distinct == 0
        alternatives += distinct > 1
        extra += max(0, distinct - 1)
    shared = tuple(sorted(tuple(sorted(label for _, label in names)) for names in by_sequence.values()
                          if len(names) > 1))
    mentioned = [i for i in inventory.mentioned_identities if i.asset == asset]
    by_locator = {r.locator: r for r in records}
    labels = {r.label for r in records if r.label is not None}
    mention_locators = {x.locator for i in mentioned for x in i.mentions}
    aliases = paralogs = paralogs_of_records = unresolved = 0
    for observation in inventory.observations:
        record = observation.record
        if record is None or record.asset != asset:
            continue
        match = _RELATIONSHIP_ITEM_RE.fullmatch(observation.locator.pointer[len(record.pointer):])
        if match is None:
            continue
        if match.group(1) == "aliases":
            aliases += 1
            continue
        paralogs += 1
        if observation.interpretation != "interpreted" or by_locator[record].status != "parsed":
            unresolved += 1
        elif observation.value.get() in labels:
            paralogs_of_records += 1
        elif observation.locator not in mention_locators:
            unresolved += 1
    of_mentioned = sum(len(i.mentions) for i in mentioned)
    return ExportAccount(
        asset=asset, record_occurrences=len(records), parsed_records=statuses["parsed"],
        malformed_records=statuses["malformed"], unsupported_records=statuses["unsupported"],
        sequence_candidates=sum(len(d) for d in digests.values()),
        records_without_candidate=sum(1 for d in digests.values() if not d),
        record_identities=len(identities), identities_without_candidate=bare, identity_sequences=sequences,
        duplicate_occurrences=duplicates, parsed_without_candidate=without, alternative_identities=alternatives,
        alternative_sequences=extra, shared_sequences=shared, mentioned_only_identities=len(mentioned),
        named_identities=len(identities) + len(mentioned), alias_names=aliases, paralog_names=paralogs,
        paralog_names_of_records=paralogs_of_records, paralog_names_of_mentioned_only=of_mentioned,
        paralog_names_unresolved=unresolved,
        support_references=sum(1 for s in inventory.support if s.locator.asset == asset),
        provider_refs_with_ext_suffix=sum(1 for r in records if (r.provider_ref or "").endswith(EXT_SUFFIX)),
        segments=inventory.export_counts(asset))


def _check_pair(lock: SourceLock, inventory: SourceInventory) -> None:
    require_type(lock, SourceLock, "lock")
    require_type(inventory, SourceInventory, "inventory")
    if inventory.source_lock_hash != lock_hash(lock):
        raise bad("/source_lock_hash", f"the inventory was built under lock {inventory.source_lock_hash}, not "
                                       f"this lock {lock_hash(lock)}", "digest_mismatch")


def account(lock: SourceLock, inventory: SourceInventory) -> SourceAccount:
    """The accounting of ``inventory``, built from ``lock`` (``digest_mismatch`` otherwise)."""
    _check_pair(lock, inventory)
    exports = tuple(_export_account(inventory, a.asset_id) for a in lock.assets if a.format in RECORD_EXPORT_FORMATS)
    fasta = tuple((a.asset_id, sum(1 for o in inventory.observations if o.locator.asset == a.asset_id))
                  for a in lock.assets if a.format in _FASTA_FORMATS)
    return SourceAccount(exports=exports, fasta=fasta, comparison=_endpoint_comparison(lock, inventory))


def _endpoint_comparison(lock: SourceLock, inventory: SourceInventory) -> EndpointComparison:
    exports = tuple(a.asset_id for a in lock.assets if a.format in RECORD_EXPORT_FORMATS)
    missing = tuple(f for f in RECORD_EXPORT_FORMATS if f not in exports)
    findings = tuple(sorted(f.id for f in inventory.findings if f.code in COMPARISON_CODES))
    both, only_in = 0, {}
    if not missing:
        identities = {a: {(r.segment, r.label) for r in inventory.records
                          if r.locator.asset == a and r.status == "parsed"} for a in exports}
        first, second = (identities[a] for a in RECORD_EXPORT_FORMATS)
        both = len(first & second)
        only_in = {a: tuple(sorted(label for _, label in identities[a] - identities[b]))
                   for a, b in zip(RECORD_EXPORT_FORMATS, reversed(RECORD_EXPORT_FORMATS))}
    return EndpointComparison(complete=not missing, exports=exports, missing=missing, identities_in_both=both,
                              only_in=only_in, findings=findings)


def comparison(lock: SourceLock, inventory: SourceInventory) -> Comparison:
    """The build report's comparison section for ``inventory`` of ``lock``: complete exactly when
    both record exports are locked, with the ids of the comparison findings."""
    _check_pair(lock, inventory)
    return _endpoint_comparison(lock, inventory).report()
