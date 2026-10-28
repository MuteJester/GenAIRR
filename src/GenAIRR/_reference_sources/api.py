"""``OgrdbSource``: one explicit OGRDB germline-set release, locked to its exact bytes.

A source comes from an explicit download (:meth:`OgrdbSource.fetch`, the only operation that
uses the network), from local files (:meth:`OgrdbSource.from_files`) or from an existing lock
and a verified cache (:meth:`OgrdbSource.read`).  It never changes and never exposes its bytes
for modification.  :meth:`OgrdbSource.build` projects an explicit recipe over it into a
:class:`~GenAIRR._reference_sources.build.ReferenceBuild`.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
import json
import os
from types import MappingProxyType
from typing import Mapping

from GenAIRR._artifact_io import ArtifactIOError, canonical_json, write_atomic
from GenAIRR._reference_sources import evidence, ogrdb
from GenAIRR._reference_sources._fields import read_json, short_digest
from GenAIRR._reference_sources.acquisition import (
    DATA_TERMS_FORMATS, AcquiredSource, HttpsTransport, Transport, fetch_source, ingest_files, read_locked,
)
from GenAIRR._reference_sources.build import RECIPE_KEYS, ReferenceBuild, document_bytes
from GenAIRR._reference_sources.consistency import profile_scope_problem
from GenAIRR._reference_sources.errors import ReferenceSourceError
from GenAIRR._reference_sources.inventory import (
    CLAIM_STATUSES, COORDINATE_PREDICATES, PREDICATES, SEVERITIES, SourceInventory,
)
from GenAIRR._reference_sources.lock import SourceLock
from GenAIRR._reference_sources.projection import project
from GenAIRR._reference_sources.recipe import MAX_RECIPE_BYTES, Recipe
from GenAIRR._reference_sources.report import asset_lines, segment_counts_text, terms_text

# The receipt clock: ``None`` records the current UTC time.  A test may substitute a callable
# returning a fixed timezone-aware datetime; receipts never enter an identity.
_now = None


def _transport() -> Transport:
    """The transport of an explicit fetch: the standard-library HTTPS client."""
    return HttpsTransport()


def _is_path(value) -> bool:
    return isinstance(value, (str, os.PathLike)) and not isinstance(value, bool)


def _named(exc: ReferenceSourceError, origin) -> ReferenceSourceError:
    """``exc`` with a message naming the file it is about (once)."""
    if str(origin) in exc.message:
        return exc
    return ReferenceSourceError(exc.kind, exc.path, f"recipe {origin}: {exc.message}")


def _first_token(path: str) -> str:
    return path.split("/")[1] if path.startswith("/") else ""


def _parse_recipe(recipe) -> tuple:
    """``(Recipe, file it was read from or None)`` from a recipe, a recipe document or a JSON file path."""
    if isinstance(recipe, Recipe):
        return recipe, None
    if isinstance(recipe, dict):
        return Recipe.from_dict(recipe), None
    if not _is_path(recipe):
        raise TypeError(f"recipe must be a recipe (such as ReferenceBuild.recipe), a recipe document (dict) or a "
                        f"JSON file path, got {type(recipe).__name__}")
    try:
        return Recipe.from_dict(read_json(recipe, max_bytes=MAX_RECIPE_BYTES)), recipe
    except ReferenceSourceError as exc:
        raise _named(exc, recipe) from None


def _relocated(exc: ReferenceSourceError, old: str, new: str) -> ReferenceSourceError:
    if exc.path == old or exc.path.startswith(old + "/"):
        return ReferenceSourceError(exc.kind, new + exc.path[len(old):], exc.message)
    return exc


# ── inspection ───────────────────────────────────────────────────────


@dataclass(frozen=True)
class SourceInspection:
    """What a locked source holds, as interpreted under its profile: identity, assets, scope,
    endpoint comparison, per-export counts, finding and claim counts, and every record.

    ``to_dict()`` returns a fresh JSON-safe copy; ``str()`` renders it as text.  An inspection is
    a view of the source and is not part of any identity.
    """

    _document: str

    @classmethod
    def _of(cls, lock: SourceLock, inventory: SourceInventory) -> "SourceInspection":
        comparison = ogrdb.comparison(lock, inventory)
        problem = profile_scope_problem(lock)
        findings = Counter((f.code, f.severity) for f in inventory.findings)
        claims = Counter((c.predicate, c.status) for c in inventory.claims)
        document = {
            "provider": lock.provider, "set_id": lock.set_id, "release": lock.release,
            "identity": lock.identity.to_dict(), "adapter": lock.adapter.to_dict(), "lock_hash": lock.commitment,
            "scope": {"in_scope": problem is None, "problem": problem},
            "assets": [{"asset_id": a.asset_id, "role": a.role, "sha256": a.sha256, "bytes": a.bytes,
                        "origin_kind": a.origin_kind, "declared_url": a.declared_url} for a in lock.assets],
            "missing_assets": [m.to_dict() for m in lock.missing_assets],
            "terms": lock.terms.to_dict(),
            "comparison": {"complete": comparison.complete, "exports": list(comparison.exports),
                           "missing": list(comparison.missing), "findings": len(comparison.findings)},
            "exports": [{"asset": asset, **inventory.export_counts(asset)} for asset in comparison.exports],
            "findings": [{"code": code, "severity": severity, "count": findings[(code, severity)]}
                         for code, severity in sorted(findings, key=lambda k: (k[0], SEVERITIES.index(k[1])))],
            "claims": [{"predicate": predicate, "status": status, "count": claims[(predicate, status)]}
                       for predicate in PREDICATES for status in CLAIM_STATUSES if claims[(predicate, status)]],
            "accepted_coordinate_claims": sum(1 for c in inventory.claims
                                              if c.status == "accepted" and c.predicate in COORDINATE_PREDICATES),
            "records": [{"record": r.id, "status": r.status, "segment": r.segment, "label": r.label,
                         "candidates": len(inventory.candidates_of(r.locator)), "problem": r.problem}
                        for r in inventory.records],
        }
        return cls(canonical_json(document))

    def to_dict(self) -> dict:
        """A fresh, JSON-safe copy of the inspection."""
        return json.loads(self._document)

    @property
    def set_id(self) -> str:
        return self.to_dict()["set_id"]

    @property
    def release(self) -> int:
        return self.to_dict()["release"]

    @property
    def in_scope(self) -> bool:
        """Whether the interpretation profile has been validated for these exact bytes."""
        return self.to_dict()["scope"]["in_scope"]

    @property
    def comparison_complete(self) -> bool:
        """Whether both record exports are present, so that their differences are compared."""
        return self.to_dict()["comparison"]["complete"]

    def __repr__(self) -> str:
        d = self.to_dict()
        scope = "in scope" if d["scope"]["in_scope"] else "out of scope"
        comparison = "complete" if d["comparison"]["complete"] else "incomplete"
        return (f"<SourceInspection {d['provider']} {d['set_id']} release {d['release']}: {len(d['records'])} records, "
                f"{scope}, comparison {comparison}; lock {short_digest(d['lock_hash'])}>")

    def __str__(self) -> str:
        d = self.to_dict()
        identity, adapter = d["identity"], d["adapter"]
        profile = adapter["coordinate_profile"]
        lines = [f"{d['provider'].upper()} source {d['set_id']} release {d['release']} "
                 f"({identity['taxon']['scientific_name']}, taxonomy {identity['taxon']['taxonomy_id']}, "
                 f"{identity['locus']}, {identity['chain_type']})",
                 f"Lock: {d['lock_hash']}",
                 f"Adapter: {adapter['id']} v{adapter['version']}; interpretation profile {profile['id']} "
                 f"v{profile['version']}"]
        if d["scope"]["in_scope"]:
            lines.append("Profile scope: validated for these exact bytes")
        else:
            lines += ["Profile scope: outside the validated scope",
                      f"  {d['scope']['problem']}",
                      f"  No coordinate claim is accepted ({d['accepted_coordinate_claims']} accepted anchor, phase or "
                      "interval claims): their values stay candidates, and a build from this source is refused"]
        lines.append("Assets:")
        lines += asset_lines(d["assets"], d["missing_assets"], "  ")
        terms = d["terms"]
        lines.append("Data terms: " + terms_text(terms))
        if terms["citation"]:
            lines.append(f"Citation: {terms['citation']}")
        comparison = d["comparison"]
        if comparison["complete"]:
            lines.append(f"Endpoint comparison: complete ({', '.join(comparison['exports'])}); "
                         f"{comparison['findings']} comparison finding{'' if comparison['findings'] == 1 else 's'}")
        else:
            lines.append(f"Endpoint comparison: incomplete: {', '.join(comparison['missing'])} "
                         f"{'is' if len(comparison['missing']) == 1 else 'are'} not locked, so differences between "
                         "the record exports are not compared")
        lines.append("Exports:")
        for counts in d["exports"]:
            lines.append(f"  {counts['asset']}: records {segment_counts_text(counts['record_occurrences'])}; sequence "
                         f"candidates {segment_counts_text(counts['sequence_candidates'])}; mentioned-only identities "
                         f"{segment_counts_text(counts['mentioned_identities'])}")
        lines.append("Findings (code, severity: count):")
        lines += [f"  {f['code']}, {f['severity']}: {f['count']}" for f in d["findings"]] or ["  none"]
        lines.append("Claims (predicate, status: count):")
        lines += [f"  {c['predicate']}, {c['status']}: {c['count']}" for c in d["claims"]] or ["  none"]
        lines.append(f"Records ({len(d['records'])}):")
        for record in d["records"]:
            line = (f"  {record['record']}  {record['status']}  {record['segment'] or '-'}  {record['label'] or '-'}"
                    f"  {record['candidates']} candidate{'' if record['candidates'] == 1 else 's'}")
            lines.append(line + (f"  ({record['problem']})" if record["problem"] else ""))
        return "\n".join(lines)


class OgrdbSource:
    """One explicit OGRDB germline-set release, locked to the exact bytes of its exports.

    ``set_id``, ``release`` and ``profile`` identify it; ``lock`` is its source lock (asset
    digests and sizes, set identity, interpretation profile, data terms); ``receipts`` are the
    acquisition events of this session (outside every identity).
    """

    __slots__ = ("_lock", "_assets", "_receipts")

    def __init__(self, *args, **kwargs):
        raise TypeError("create an OgrdbSource with OgrdbSource.fetch, OgrdbSource.from_files or OgrdbSource.read")

    def __setattr__(self, name, value):
        raise AttributeError("OgrdbSource is immutable")

    def __delattr__(self, name):
        raise AttributeError("OgrdbSource is immutable")

    def __copy__(self):
        return self

    def __deepcopy__(self, memo):
        return self

    def __repr__(self) -> str:
        return (f"<OgrdbSource {self.set_id} release {self.release} ({self.profile}): {len(self._lock.assets)} "
                f"assets; lock {short_digest(self._lock.commitment)}>")

    @classmethod
    def _of(cls, acquired: AcquiredSource) -> "OgrdbSource":
        self = object.__new__(cls)
        assets = MappingProxyType({asset.asset_id: bytes(acquired.assets[asset.asset_id])
                                   for asset in acquired.lock.assets})
        for name, value in (("_lock", acquired.lock), ("_assets", assets), ("_receipts", tuple(acquired.receipts))):
            object.__setattr__(self, name, value)
        return self

    # Acquisition -------------------------------------------------------

    @classmethod
    def fetch(cls, *, set_id: str, release: int, profile: str, cache,
              formats=("airr", "airr_ex", "gapped", "gapped_ex"), expected_lock=None) -> "OgrdbSource":
        """Download the requested exports of one explicit set release and lock their bytes.

        ``release`` is an explicit positive integer (never ``"latest"``); ``profile`` is a known
        interpretation profile id, such as ``"ogrdb_api_v2_human_ig_observed_v1"``; ``formats``
        lists exports among ``airr``, ``airr_ex``, ``gapped`` and ``gapped_ex``, by default all
        four: a build compares both record exports, and V junction anchors are accepted only when
        bound to the gapped exports; ``cache`` is the content-addressed blob directory.  This is the only operation that uses the network: one
        HTTPS request per format to the generated OGRDB endpoint, at most three redirects on the
        same origin, bounded sizes and a deadline, no automatic retries.

        ``expected_lock`` (a lock or a lock file path) reacquires every requested export and
        requires exactly its pinned bytes: bytes changed under the same release are reported as
        drift and nothing is stored or replaced.

        Downloads use the standard-library HTTPS client, which does not use proxy settings
        (``HTTPS_PROXY`` and similar).  Behind a proxy, download the exports with your own tools
        and record them with :meth:`from_files`.
        """
        if expected_lock is not None and not isinstance(expected_lock, SourceLock):
            if not _is_path(expected_lock):
                raise TypeError(f"expected_lock must be a source lock or a lock file path, got "
                                f"{type(expected_lock).__name__}")
            try:
                expected_lock = SourceLock.read(expected_lock)
            except ReferenceSourceError as exc:
                raise exc.under("/expected_lock") from None
        return cls._of(fetch_source(set_id=set_id, release=release, profile=profile, cache=cache,
                                    transport=_transport(), identify=ogrdb.identify, formats=formats,
                                    expected_lock=expected_lock, now=_now))

    @classmethod
    def from_files(cls, *, set_id: str, release: int, profile: str, exports, cache) -> "OgrdbSource":
        """Record local copies of the exports of one explicit set release, never touching the network.

        ``exports`` maps an export format (``airr``, ``airr_ex``, ``gapped``, ``gapped_ex``,
        ``ungapped`` or ``ungapped_ex``) to a file path.  A build needs ``airr`` and ``airr_ex``, and
        V junction anchors need ``gapped`` and ``gapped_ex``.  The bytes are stored in ``cache`` and
        locked as local files: the lock does not claim that GenAIRR fetched them, and it depends on
        the bytes only, never on paths or file names.  The provider's data-terms pages cannot be
        recorded here: the lock records the data terms as unresolved.
        """
        if isinstance(exports, Mapping):
            pages = [fmt for fmt in DATA_TERMS_FORMATS if fmt in exports]
            if pages:
                raise ReferenceSourceError("malformed", "/exports/" + pages[0], (
                    f"{', '.join(pages)}: data-terms pages cannot be recorded through from_files; omit them (a lock "
                    "recorded from local files states its data terms as unresolved)"))
        try:
            acquired = ingest_files(set_id=set_id, release=release, profile=profile, files=exports, cache=cache,
                                    identify=ogrdb.identify, now=_now)
        except ReferenceSourceError as exc:
            raise _relocated(exc, "/files", "/exports") from None
        return cls._of(acquired)

    @classmethod
    def read(cls, lock_path, *, cache) -> "OgrdbSource":
        """Read a lock file and verify every blob it needs in ``cache``, never touching the network.

        A missing or corrupt blob is an error naming the asset and its pinned digest; reading never
        downloads.
        """
        return cls._of(read_locked(lock_path, cache=cache, identify=ogrdb.identify, now=_now))

    # Identity ----------------------------------------------------------

    @property
    def set_id(self) -> str:
        return self._lock.set_id

    @property
    def release(self) -> int:
        return self._lock.release

    @property
    def profile(self) -> str:
        """The interpretation profile id."""
        return self._lock.adapter.coordinate_profile.id

    @property
    def lock(self) -> SourceLock:
        return self._lock

    @property
    def receipts(self) -> tuple:
        """Acquisition receipts of this session, one per locked asset (never part of an identity)."""
        return self._receipts

    # Operations --------------------------------------------------------

    def write_lock(self, path, *, replace: bool = False) -> None:
        """Write the canonical lock (canonical JSON and one LF) atomically.

        An existing file is an error unless ``replace=True``.
        """
        try:
            write_atomic(path, document_bytes(self._lock.to_dict()), replace=replace)
        except ArtifactIOError as exc:
            if not replace and os.path.lexists(path):
                raise ReferenceSourceError("io", "", f"cannot write the source lock to {path}: the lock file exists; "
                                                     "pass replace=True to replace it") from exc
            raise ReferenceSourceError("io", "", f"cannot write the source lock to {path}: {exc}") from exc

    def inspect(self) -> SourceInspection:
        """What the source holds, interpreted under its profile; never touches the network.

        Works for a single export (the endpoint comparison is then incomplete) and for bytes
        outside the profile's validated scope (no coordinate claim is then accepted).
        """
        inventory = evidence._interpret_locked(self._lock, self._assets)      # reads the bytes once
        return SourceInspection._of(self._lock, inventory)

    def build(self, recipe) -> ReferenceBuild:
        """Project ``recipe`` (a recipe document or a JSON file path) over this source.

        Never touches the network.  The recipe is validated as a closed document; errors keep their
        location in it and, for a recipe file, name the file.  A build compares both record exports,
        so it needs ``airr`` and ``airr_ex`` of the same release.
        """
        parsed, origin = _parse_recipe(recipe)
        try:
            projection = project(self._lock, self._assets, parsed)
        except ReferenceSourceError as exc:
            if origin is not None and _first_token(exc.path) in RECIPE_KEYS:     # located in the recipe
                raise _named(exc, origin) from None
            raise
        return ReferenceBuild._create(lock=self._lock, recipe=parsed, assets=self._assets,
                                      inputs=projection.inputs, reference=projection.reference,
                                      report=projection.report, receipts=self._receipts)
