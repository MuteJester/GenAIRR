"""Reference builds and their portable build dossiers.

A :class:`ReferenceBuild` holds a runtime :class:`~GenAIRR._reference.artifact.Reference`, its
complete build report, the source lock and recipe it was built from, and the manifest of its
dossier.  The manifest is computed when the build is made; :meth:`ReferenceBuild.write_dossier`
publishes exactly the members it lists.
"""
from __future__ import annotations

import hashlib
from types import MappingProxyType
from typing import Mapping, Optional

from GenAIRR._artifact_io import canonical_json
from GenAIRR._reference.artifact import Reference
from GenAIRR._reference_sources._fields import require_type
from GenAIRR._reference_sources.digests import report_hash
from GenAIRR._reference_sources.lock import SourceLock
from GenAIRR._reference_sources.manifest import (
    BLOB_PREFIX, MEMBER_SCHEMA_VERSIONS, BuildManifest, ManifestHashes, Member, MemberVersions, check_manifest,
)
from GenAIRR._reference_sources.projection import ANNOTATIONS_MEMBER, FASTA_MEMBERS
from GenAIRR._reference_sources.recipe import Recipe
from GenAIRR._reference_sources.report import BuildReport

LOCK_MEMBER = "source.lock.json"
RECIPE_MEMBER = "recipe.json"
REPORT_MEMBER = "report.json"
REFERENCE_MEMBER = "reference.genairr-ref.json"
INPUTS_PREFIX = "inputs/"


def document_bytes(document: dict) -> bytes:
    """A generated JSON member: canonical JSON, one LF, UTF-8."""
    return (canonical_json(document) + "\n").encode("utf-8")


def _sha256(data: bytes) -> str:
    return "sha256:" + hashlib.sha256(data).hexdigest()


def _input_names(lock: SourceLock) -> tuple:
    """The projected input files of a build of ``lock``: FASTA in V, (D,) J order, then the sidecar."""
    return tuple(FASTA_MEMBERS[segment] for segment in lock.identity.segments) + (ANNOTATIONS_MEMBER,)


def _members(lock: SourceLock, recipe: Recipe, assets: Mapping[str, bytes], inputs: Mapping[str, bytes],
             reference: Reference, report: BuildReport) -> dict:
    """``{member name: exact bytes}`` of every deterministic member, in publication order.

    Assets with identical bytes share one content-addressed blob.  Receipts and the manifest
    are not members."""
    files = {}
    for asset in lock.assets:
        files.setdefault(BLOB_PREFIX + asset.blob_name, assets[asset.asset_id])
    files[LOCK_MEMBER] = document_bytes(lock.to_dict())
    files[RECIPE_MEMBER] = document_bytes(recipe.to_dict())
    for name in _input_names(lock):
        files[INPUTS_PREFIX + name] = inputs[name]
    files[REFERENCE_MEMBER] = reference.to_json().encode("utf-8")
    files[REPORT_MEMBER] = document_bytes(report.to_dict())
    return files


def _manifest(lock: SourceLock, recipe: Recipe, report: BuildReport, files: Mapping[str, bytes]) -> BuildManifest:
    """The manifest of exactly ``files``: member digests and sizes, schema versions, adapter and
    projection, the input and report hashes, the runtime hashes and the builder of the report."""
    hashes = report.hashes
    return BuildManifest(
        builder=report.generator, versions=MemberVersions(**MEMBER_SCHEMA_VERSIONS), adapter=lock.adapter,
        projection=recipe.projection,
        hashes=ManifestHashes(source_lock=hashes.source_lock, inventory=hashes.inventory, recipe=hashes.recipe,
                              build_inputs=hashes.build_inputs, report=report_hash(report)),
        runtime=report.runtime,
        members=tuple(Member(name, _sha256(data), len(data)) for name, data in sorted(files.items())))


class ReferenceBuild:
    """A runtime reference built from one locked published source and an explicit recipe.

    ``reference`` is the runtime :class:`~GenAIRR.Reference`; ``report`` is the complete build
    report (``print(build.report)`` for a readable account, ``build.report.to_dict()`` for the
    data); ``lock`` and ``recipe`` are what it was built from; ``manifest`` lists every
    deterministic dossier member with its digest and size.  A build never changes.

    Builds come from ``OgrdbSource.build``, :meth:`read_dossier` and :meth:`bundled`.
    """

    __slots__ = ("_reference", "_report", "_lock", "_recipe", "_manifest", "_files", "_receipts", "_checked_with")

    def __init__(self, *args, **kwargs):
        raise TypeError("a ReferenceBuild comes from OgrdbSource.build, ReferenceBuild.read_dossier or "
                        "ReferenceBuild.bundled")

    def __setattr__(self, name, value):
        raise AttributeError("ReferenceBuild is immutable")

    def __delattr__(self, name):
        raise AttributeError("ReferenceBuild is immutable")

    def __copy__(self):
        return self

    def __deepcopy__(self, memo):
        return self

    def __repr__(self) -> str:
        lock = self._lock
        return (f"<ReferenceBuild {self._reference.name!r} from {lock.provider} {lock.set_id} release "
                f"{lock.release} {self._reference.execution_hash}>")

    @classmethod
    def _create(cls, *, lock: SourceLock, recipe: Recipe, assets: Mapping[str, bytes], inputs: Mapping[str, bytes],
                reference: Reference, report: BuildReport, receipts: tuple = (),
                checked_with: Optional[str] = None) -> "ReferenceBuild":
        require_type(lock, SourceLock, "lock")
        require_type(recipe, Recipe, "recipe")
        require_type(reference, Reference, "reference")
        require_type(report, BuildReport, "report")
        files = _members(lock, recipe, assets, inputs, reference, report)
        manifest = _manifest(lock, recipe, report, files)
        check_manifest(manifest, lock, report)
        self = object.__new__(cls)
        for name, value in (("_reference", reference), ("_report", report), ("_lock", lock), ("_recipe", recipe),
                            ("_manifest", manifest), ("_files", MappingProxyType(files)),
                            ("_receipts", tuple(receipts)), ("_checked_with", checked_with)):
            object.__setattr__(self, name, value)
        return self

    # Results -----------------------------------------------------------

    @property
    def reference(self) -> Reference:
        """The runtime reference: simulate with ``Experiment.on(build.reference)``."""
        return self._reference

    @property
    def report(self) -> BuildReport:
        """The complete build report: the source inventory, the selection, every record and claim
        outcome, the runtime hashes and capabilities.  ``reference.report()`` describes only the
        runtime catalogue."""
        return self._report

    @property
    def lock(self) -> SourceLock:
        """The source lock the build was made from (exact asset digests and sizes)."""
        return self._lock

    @property
    def recipe(self) -> Recipe:
        """The explicit recipe the build was made from."""
        return self._recipe

    @property
    def manifest(self) -> BuildManifest:
        """The dossier manifest: every deterministic member's name, digest and size, the member
        schema versions, the input, report and runtime hashes and the builder."""
        return self._manifest

    @property
    def checked_with(self) -> Optional[str]:
        """The GenAIRR version that checked this build's dossier, or ``None`` for a build made in
        this process.  It is a property of the check, never recorded in a dossier."""
        return self._checked_with
