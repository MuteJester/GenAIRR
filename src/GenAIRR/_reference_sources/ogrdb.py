"""OGRDB germline-set adapter: set identity read from the locked export bytes.

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
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass
from typing import Any, Mapping, Optional

from GenAIRR._reference_sources._fields import (
    MAX_ASSET_BYTES, MAX_NAME, MAX_TAXONOMY_ID, bad, check_text, parse_json, ptr, show,
)
from GenAIRR._reference_sources.common import (
    LOCUS_CHAIN_TYPES, RECORD_EXPORT_FORMATS, SourceIdentity, Taxon,
)
from GenAIRR._reference_sources.errors import ReferenceSourceError

# ── Provider conventions ────────────────────────────────────────────

SET_POINTER = "/GermlineSet/0"
RECORDS_POINTER = SET_POINTER + "/allele_descriptions"
TAXON_PREFIX = "NCBITAXON:"
# The extended export names its set with this suffix (``IGH_VDJ`` / ``IGH_VDJ_extended``).
EXTENDED_SET_SUFFIX = "_extended"

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


def _mismatch(envelope: _Envelope, field: str, message: str) -> ReferenceSourceError:
    return bad(envelope.at(SET_POINTER + "/" + field), f"{envelope.asset}: {message}", "identity_mismatch")


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
    for field, mine, theirs in pairs:
        if mine is not None and theirs is not None and mine != theirs:
            raise _mismatch(other, field, f"{field} {show(theirs)} disagrees with {show(mine)} in {first.asset}; "
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
