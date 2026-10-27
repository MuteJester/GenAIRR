"""``OgrdbSource``: one explicit OGRDB germline-set release, locked to its exact bytes.

A source comes from an explicit download (:meth:`OgrdbSource.fetch`, the only operation that
uses the network), from local files (:meth:`OgrdbSource.from_files`) or from an existing lock
and a verified cache (:meth:`OgrdbSource.read`).  It never changes and never exposes its bytes
for modification.  :meth:`OgrdbSource.build` projects an explicit recipe over it into a
:class:`~GenAIRR._reference_sources.build.ReferenceBuild`.
"""
from __future__ import annotations

import os
from types import MappingProxyType
from typing import Optional

from GenAIRR._artifact_io import ArtifactIOError, write_atomic
from GenAIRR._reference_sources import ogrdb
from GenAIRR._reference_sources._fields import read_json
from GenAIRR._reference_sources.acquisition import (
    AcquiredSource, HttpsTransport, Transport, fetch_source, ingest_files, read_locked,
)
from GenAIRR._reference_sources.build import ReferenceBuild, document_bytes
from GenAIRR._reference_sources.errors import ReferenceSourceError
from GenAIRR._reference_sources.lock import SourceLock
from GenAIRR._reference_sources.projection import project
from GenAIRR._reference_sources.recipe import MAX_RECIPE_BYTES, Recipe

# The receipt clock: ``None`` records the current UTC time.  A test may substitute a callable
# returning a fixed timezone-aware datetime; receipts never enter an identity.
_now = None

# Top-level recipe keys: an error located below one of them is about the recipe document.
_RECIPE_KEYS = ("artifact", "schema_version", "source_lock_hash", "adapter", "projection", "reference",
                "selection", "decisions", "requires", "defaults_profile")


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
    """``(Recipe, file it was read from or None)`` from a recipe document or a JSON file path."""
    if isinstance(recipe, dict):
        return Recipe.from_dict(recipe), None
    if not _is_path(recipe):
        raise TypeError(f"recipe must be a recipe document (dict) or a JSON file path, got {type(recipe).__name__}")
    try:
        return Recipe.from_dict(read_json(recipe, max_bytes=MAX_RECIPE_BYTES)), recipe
    except ReferenceSourceError as exc:
        raise _named(exc, recipe) from None


def _relocated(exc: ReferenceSourceError, old: str, new: str) -> ReferenceSourceError:
    if exc.path == old or exc.path.startswith(old + "/"):
        return ReferenceSourceError(exc.kind, new + exc.path[len(old):], exc.message)
    return exc


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
        return f"<OgrdbSource {self.set_id} release {self.release} ({self.profile}) {self._lock.commitment}>"

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
    def fetch(cls, *, set_id: str, release: int, profile: str, cache, formats=("airr", "airr_ex"),
              expected_lock=None) -> "OgrdbSource":
        """Download the requested exports of one explicit set release and lock their bytes.

        ``release`` is an explicit positive integer (never ``"latest"``); ``profile`` is a known
        interpretation profile id, such as ``"ogrdb_api_v2_human_ig_observed_v1"``; ``formats``
        lists exports among ``airr``, ``airr_ex``, ``gapped`` and ``gapped_ex``; ``cache`` is the
        content-addressed blob directory.  This is the only operation that uses the network: one
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

        ``exports`` maps a format (``airr``, ``airr_ex``, ``gapped``, ``gapped_ex``, ...) to a file
        path.  The bytes are stored in ``cache`` and locked as local files: the lock does not claim
        that GenAIRR fetched them, and it depends on the bytes only, never on paths or file names.
        """
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
            hint = "; pass replace=True to replace it" if not replace and os.path.lexists(path) else ""
            raise ReferenceSourceError("io", "", f"cannot write the source lock to {path}: {exc}{hint}") from exc

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
            if origin is not None and _first_token(exc.path) in _RECIPE_KEYS:
                raise _named(exc, origin) from None
            raise
        return ReferenceBuild._create(lock=self._lock, recipe=parsed, assets=self._assets,
                                      inputs=projection.inputs, reference=projection.reference,
                                      report=projection.report, receipts=self._receipts)
