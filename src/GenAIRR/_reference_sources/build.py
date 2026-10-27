"""Reference builds and their portable build dossiers.

A :class:`ReferenceBuild` holds a runtime :class:`~GenAIRR._reference.artifact.Reference`, its
complete build report, the source lock and recipe it was built from, and the manifest of its
dossier.  The manifest is computed when the build is made; :meth:`ReferenceBuild.write_dossier`
publishes exactly the members it lists.
"""
from __future__ import annotations

from collections import deque
import hashlib
import json
import os
from pathlib import Path
import re
import stat
import threading
from types import MappingProxyType
from typing import Callable, Mapping, Optional

from GenAIRR._artifact_io import ArtifactIOError, canonical_json, write_atomic
from GenAIRR._reference.artifact import MAX_BYTES as REFERENCE_MAX_BYTES, Reference
from GenAIRR._reference.errors import ReferenceArtifactError
from GenAIRR._reference.importer import FASTA_MAX_BYTES
from GenAIRR._reference.sidecar import MAX_BYTES as SIDECAR_MAX_BYTES
from GenAIRR._reference_sources import bundled
from GenAIRR._reference_sources._fields import (
    MAX_ASSET_BYTES, MAX_ASSETS, bad, parse_json, ptr, require_type, show,
)
from GenAIRR._reference_sources.digests import manifest_hash, recipe_hash, report_hash
from GenAIRR._reference_sources.errors import ReferenceSourceError
from GenAIRR._reference_sources.lock import MAX_LOCK_BYTES, SourceLock
from GenAIRR._reference_sources.manifest import (
    BLOB_PREFIX, MANIFEST_NAME, MAX_MANIFEST_BYTES, MAX_MEMBERS, MEMBER_SCHEMA_VERSIONS, RECEIPTS_PREFIX,
    BuildManifest, ManifestHashes, Member, MemberVersions, check_manifest,
)
from GenAIRR._reference_sources.projection import ANNOTATIONS_MEMBER, FASTA_MEMBERS, project, report_for
from GenAIRR._reference_sources.receipt import MAX_RECEIPT_BYTES, AcquisitionReceipt
from GenAIRR._reference_sources.recipe import MAX_RECIPE_BYTES, Recipe
from GenAIRR._reference_sources.report import MAX_REPORT_BYTES, BuildReport

LOCK_MEMBER = "source.lock.json"
RECIPE_MEMBER = "recipe.json"
REPORT_MEMBER = "report.json"
REFERENCE_MEMBER = "reference.genairr-ref.json"
INPUTS_PREFIX = "inputs/"

# Reading budget: every member's declared size is checked against its own limit and the sum
# against the dossier limit before any member is read.  Receipts: at most one per lock asset.
MAX_DOSSIER_BYTES = 1 << 30
MAX_RECEIPTS = MAX_ASSETS
_MEMBER_LIMITS = {LOCK_MEMBER: MAX_LOCK_BYTES, RECIPE_MEMBER: MAX_RECIPE_BYTES, REPORT_MEMBER: MAX_REPORT_BYTES,
                  REFERENCE_MEMBER: REFERENCE_MAX_BYTES, INPUTS_PREFIX + ANNOTATIONS_MEMBER: SIDECAR_MAX_BYTES,
                  **{INPUTS_PREFIX + name: FASTA_MAX_BYTES for name in FASTA_MEMBERS.values()}}
# The only directories a dossier holds; anything else is not part of it.
_DIRECTORIES = ("cache", "cache/blobs", "cache/blobs/sha256", "inputs", "receipts")
_RECEIPT_RE = re.compile(r"receipts/[A-Za-z0-9][A-Za-z0-9._-]{0,63}\.json")
_ARTIFACT_KINDS = {"io": "io", "unsupported_schema": "unsupported_schema", "integrity": "digest_mismatch"}
# Packaged builds already checked in this process, by (identifier, registry entry); one lock per
# key makes concurrent first loads run a single check.
_BUNDLED: dict = {}
_BUNDLED_LOCKS: dict = {}
_BUNDLED_LOCK = threading.Lock()


def document_bytes(document: dict) -> bytes:
    """A generated JSON member: canonical JSON, one LF, UTF-8."""
    return (canonical_json(document) + "\n").encode("utf-8")


def _sha256(data: bytes) -> str:
    return "sha256:" + hashlib.sha256(data).hexdigest()


def _path(directory) -> Path:
    if isinstance(directory, bool) or not isinstance(directory, (str, os.PathLike)):
        raise TypeError(f"directory must be a filesystem path, got {type(directory).__name__}")
    return Path(directory)


def _receipt_name(index: int) -> str:
    return f"{RECEIPTS_PREFIX}{index}.json"


def _receipt_order(name: str) -> tuple:
    stem = name[len(RECEIPTS_PREFIX):-len(".json")]
    return (0, int(stem), "") if stem.isdigit() else (1, 0, stem)


# ── reading a dossier ────────────────────────────────────────────────


def _is_link(entry: os.DirEntry) -> bool:
    """A symbolic link, or on Windows any reparse point (junctions included, on every version)."""
    if entry.is_symlink():
        return True
    attributes = getattr(entry.stat(follow_symlinks=False), "st_file_attributes", 0)   # Windows only
    return bool(attributes & stat.FILE_ATTRIBUTE_REPARSE_POINT)


def _scan(root: Path) -> dict:
    """Every regular file below ``root`` as ``{relative name: size}``, without following anything.

    A symbolic link (or junction) anywhere, a special file or a directory that is not part of a
    dossier is refused; the number of entries is bounded."""
    files: dict = {}
    limit = MAX_MEMBERS + MAX_RECEIPTS + len(_DIRECTORIES) + 1
    seen = 0
    pending = deque([""])
    while pending:
        folder = pending.popleft()
        location = root.joinpath(*folder.split("/")) if folder else root
        entries = []
        try:
            with os.scandir(location) as listing:
                for item in listing:
                    seen += 1
                    if seen > limit:
                        raise bad("", f"{root} holds more than {limit} entries; a dossier holds at most {MAX_MEMBERS} "
                                      f"members and {MAX_RECEIPTS} receipts", "resource_limit")
                    entries.append(item)
        except OSError as exc:
            raise bad("/" + folder if folder else "", f"cannot list {location}: {exc}", "io") from exc
        for item in sorted(entries, key=lambda e: e.name):
            name = f"{folder}/{item.name}" if folder else item.name
            if _is_link(item):
                raise bad("/" + name, f"{name} is a symbolic link or another reparse point; dossier members are "
                                      "regular files and links are never followed", "io")
            if item.is_dir(follow_symlinks=False):
                if name not in _DIRECTORIES:
                    raise bad("/" + name, f"unexpected directory {name}: it is not part of a dossier", "malformed")
                pending.append(name)
            elif item.is_file(follow_symlinks=False):
                files[name] = item.stat(follow_symlinks=False).st_size
            else:
                raise bad("/" + name, f"{name} is not a regular file", "io")
    return files


def _read_file(root: Path, name: str, limit: int) -> bytes:
    """At most ``limit + 1`` bytes of the regular file ``name`` below ``root``, never through a link."""
    path = root.joinpath(*name.split("/"))
    flags = os.O_RDONLY | getattr(os, "O_BINARY", 0) | getattr(os, "O_NOFOLLOW", 0)
    try:
        descriptor = os.open(path, flags)
    except FileNotFoundError as exc:
        raise bad("/" + name, f"{name} is listed in {MANIFEST_NAME} but absent from {root}", "missing_asset") from exc
    except OSError as exc:
        raise bad("/" + name, f"cannot read {name}: {exc}", "io") from exc
    with os.fdopen(descriptor, "rb") as handle:
        if not stat.S_ISREG(os.fstat(handle.fileno()).st_mode):
            raise bad("/" + name, f"{name} is not a regular file", "io")
        return handle.read(limit + 1)


def _text(name: str, data: bytes) -> str:
    try:
        return data.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise bad("/" + name, f"{name} is not UTF-8 text: {exc}") from exc


def _document(name: str, data: bytes, parse: Callable):
    """A generated JSON member, parsed strictly; its bytes must be its canonical form."""
    try:
        document = parse(_text(name, data))
    except ReferenceSourceError as exc:
        raise exc.under("/" + name) from None
    if document_bytes(document.to_dict()) != data:
        raise bad("/" + name, f"{name} is not in canonical form (canonical JSON, UTF-8, one LF); a dossier stores "
                              "each generated document exactly as written")
    return document


def _difference(stored, expected, path: str = "") -> Optional[str]:
    """The JSON pointer of the first difference between two JSON values, or ``None``."""
    if isinstance(stored, dict) and isinstance(expected, dict):
        for key in sorted(set(stored) | set(expected)):
            if key not in stored or key not in expected:
                return ptr(path, key)
            found = _difference(stored[key], expected[key], ptr(path, key))
            if found is not None:
                return found
        return None
    if isinstance(stored, list) and isinstance(expected, list):
        for i, (a, b) in enumerate(zip(stored, expected)):
            found = _difference(a, b, ptr(path, i))
            if found is not None:
                return found
        return None if len(stored) == len(expected) else ptr(path, min(len(stored), len(expected)))
    return None if type(stored) is type(expected) and stored == expected else path


def _value_at(document, pointer: str):
    for token in pointer.split("/")[1:]:
        token = token.replace("~1", "/").replace("~0", "~")
        try:
            document = document[int(token)] if isinstance(document, list) else document[token]
        except (KeyError, IndexError, ValueError, TypeError):
            return "(absent)"
    return document


def _compared(name: str, stored: dict, expected: dict, what: str) -> None:
    """``digest_mismatch`` at the first difference unless ``stored`` equals ``expected``."""
    where = _difference(stored, expected)
    if where is not None:
        raise bad("/" + name + where, f"{name} differs from {what} at {where or 'the top level'}: stored "
                                      f"{show(_value_at(stored, where))}, regenerated "
                                      f"{show(_value_at(expected, where))}", "digest_mismatch")


def _regeneration_failed(exc: ReferenceSourceError, lock: SourceLock) -> ReferenceSourceError:
    """A projection error located in the dossier member it concerns."""
    path = exc.path
    tokens = path.split("/") if path.startswith("/") else ["", ""]
    blobs = {a.asset_id: BLOB_PREFIX + a.blob_name for a in lock.assets}
    if tokens[1] == "source_lock":                        # /source_lock/... -> the lock member
        path = "/" + LOCK_MEMBER + path[len("/source_lock"):]
    elif tokens[1] == "assets" and len(tokens) > 2 and tokens[2] in blobs:     # /assets/<id>... -> its blob
        path = "/" + blobs[tokens[2]] + path[len("/assets/" + tokens[2]):]
    elif tokens[1] in ("artifact", "schema_version", "source_lock_hash", "adapter", "projection", "reference",
                       "selection", "decisions", "requires", "defaults_profile"):
        path = "/" + RECIPE_MEMBER + path
    return ReferenceSourceError(exc.kind, path, "the dossier's lock, blobs and recipe do not regenerate a build "
                                                f"under this GenAIRR: {exc.message}")


def _read_manifest(root: Path) -> BuildManifest:
    """The dossier's manifest, read without following links, bounded and parsed strictly.

    No manifest is ``incomplete_dossier``: it is published last, so its absence means an
    unfinished write (or not a dossier)."""
    if not root.is_dir():
        raise bad("", f"{root} is not a directory; a dossier is a directory written by write_dossier", "io")
    try:
        info = os.lstat(root / MANIFEST_NAME)
    except FileNotFoundError:
        raise bad("", f"{root} has no {MANIFEST_NAME}: the dossier is incomplete (its writing did not finish) or "
                      "is not a dossier, and it is never read", "incomplete_dossier") from None
    except OSError as exc:
        raise bad("/" + MANIFEST_NAME, f"cannot read {MANIFEST_NAME}: {exc}", "io") from exc
    if stat.S_ISLNK(info.st_mode) or not stat.S_ISREG(info.st_mode):
        raise bad("/" + MANIFEST_NAME, f"{MANIFEST_NAME} is a symbolic link or not a regular file; links are "
                                       "never followed", "io")
    data = _read_file(root, MANIFEST_NAME, MAX_MANIFEST_BYTES)
    if len(data) > MAX_MANIFEST_BYTES:
        raise bad("/" + MANIFEST_NAME, f"{MANIFEST_NAME} exceeds {MAX_MANIFEST_BYTES} bytes", "resource_limit")
    return _document(MANIFEST_NAME, data, BuildManifest.from_json)


def _input_names(lock: SourceLock) -> tuple:
    """The projected input files of a build of ``lock``: FASTA in V, (D,) J order, then the sidecar."""
    return tuple(FASTA_MEMBERS[segment] for segment in lock.identity.segments) + (ANNOTATIONS_MEMBER,)


def _member_names(lock: SourceLock) -> tuple:
    """Every deterministic member name of a build of ``lock``, in publication order.

    Assets with identical bytes share one content-addressed blob.  Receipts and the manifest
    are not members."""
    blobs = tuple(dict.fromkeys(BLOB_PREFIX + asset.blob_name for asset in lock.assets))
    return blobs + (LOCK_MEMBER, RECIPE_MEMBER) + tuple(INPUTS_PREFIX + name for name in _input_names(lock)) + (
        REFERENCE_MEMBER, REPORT_MEMBER)


def _members(lock: SourceLock, recipe: Recipe, assets: Mapping[str, bytes], inputs: Mapping[str, bytes],
             reference: Reference, report: BuildReport) -> dict:
    """``{member name: exact bytes}`` of every deterministic member, in publication order."""
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


def _manifest(lock: SourceLock, recipe: Recipe, report: BuildReport, members) -> BuildManifest:
    """The manifest of a build: ``members`` (``(name, digest, size)``), schema versions, adapter and
    projection, the input and report hashes, the runtime hashes and the builder of the report."""
    hashes = report.hashes
    return BuildManifest(
        builder=report.generator, versions=MemberVersions(**MEMBER_SCHEMA_VERSIONS), adapter=lock.adapter,
        projection=recipe.projection,
        hashes=ManifestHashes(source_lock=hashes.source_lock, inventory=hashes.inventory, recipe=hashes.recipe,
                              build_inputs=hashes.build_inputs, report=report_hash(report)),
        runtime=report.runtime,
        members=tuple(Member(name, digest, size) for name, digest, size in sorted(members)))


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
        manifest = _manifest(lock, recipe, report, ((n, _sha256(d), len(d)) for n, d in files.items()))
        check_manifest(manifest, lock, report)
        return cls._assemble(lock=lock, recipe=recipe, reference=reference, report=report, manifest=manifest,
                             files=files, receipts=receipts, checked_with=checked_with)

    @classmethod
    def _assemble(cls, *, lock: SourceLock, recipe: Recipe, reference: Reference, report: BuildReport,
                  manifest: BuildManifest, files: dict, receipts: tuple,
                  checked_with: Optional[str]) -> "ReferenceBuild":
        """A build from parts that are already consistent (``files`` in publication order)."""
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

    # Dossiers ----------------------------------------------------------

    def write_dossier(self, directory) -> None:
        """Write the build as a portable dossier into the new directory ``directory``.

        The directory must not exist (its parent must).  Every member is published with the
        atomic file writer: each locked source asset under ``cache/blobs/sha256/<digest>``, then
        ``source.lock.json``, ``recipe.json``, the exact projected ``inputs/``,
        ``reference.genairr-ref.json``, ``report.json`` and the acquisition receipts
        (``receipts/<n>.json``, outside the manifest), and ``manifest.json`` last.  A directory
        without ``manifest.json`` is incomplete and is never read.  If writing fails, the files
        already written are kept and nothing is removed or replaced.
        """
        root = _path(directory)
        try:
            os.mkdir(root)
        except OSError as exc:
            raise ReferenceSourceError("io", "", f"cannot create the dossier directory {root}: {exc}; a dossier is "
                                                 "written into a new directory only, never into or over an existing "
                                                 "one") from exc
        files = dict(self._files)
        for index, receipt in enumerate(self._receipts):
            files[_receipt_name(index)] = document_bytes(receipt.to_dict())
        files[MANIFEST_NAME] = document_bytes(self._manifest.to_dict())
        created = set()
        for name, data in files.items():
            try:
                parts = name.split("/")
                for depth in range(1, len(parts)):
                    folder = "/".join(parts[:depth])
                    if folder not in created:
                        os.mkdir(root.joinpath(*parts[:depth]))
                        created.add(folder)
                target = root.joinpath(*parts)
                if name == REFERENCE_MEMBER:
                    self._reference.write(target)          # the artifact's own writer: the same bytes
                else:
                    write_atomic(target, data, replace=False)
            except (OSError, ArtifactIOError, ReferenceArtifactError) as exc:
                raise ReferenceSourceError("io", "/" + name, (
                    f"cannot write dossier member {name} into {root}: {exc}. The dossier is incomplete: it has no "
                    f"{MANIFEST_NAME}, so it cannot be read. The files already written are kept; write the build "
                    "into a new directory")) from exc

    @classmethod
    def read_dossier(cls, directory) -> "ReferenceBuild":
        """Check the dossier in ``directory`` fully offline and return the build it stores.

        The dossier must be complete (``manifest.json`` is written last; without it nothing is
        read).  Every file must be a listed member (or an acquisition receipt), a regular file
        and never a link; each member's size and digest must be the manifest's.  The inventory,
        the projected inputs, the runtime reference and the report are then regenerated from the
        locked blobs and the recipe and compared with what is stored: the inputs byte for byte,
        the reference in every runtime field but its generator and commitment, and the report
        exactly, stated for the stored reference and its recorded builder.  A dossier whose
        hashes were recomputed after an edit is therefore still refused.

        Nothing is written.  The result holds the stored reference, report and manifest;
        ``checked_with`` names the GenAIRR version that checked it.
        """
        import GenAIRR
        checking = GenAIRR.__version__
        root = _path(directory)
        manifest = _read_manifest(root)

        # Exactly the listed members (plus receipts) are present, as regular files.
        present = _scan(root)
        listed = {member.name: member for member in manifest.members}
        extra = sorted(name for name in present if name not in listed and name != MANIFEST_NAME)
        receipts = sorted((name for name in extra if _RECEIPT_RE.fullmatch(name)), key=_receipt_order)
        unlisted = [name for name in extra if name not in receipts]
        if unlisted:
            raise bad("/" + unlisted[0], f"{', '.join(unlisted)} {'is' if len(unlisted) == 1 else 'are'} not listed in "
                                         f"{MANIFEST_NAME}; a dossier holds only its listed members and acquisition "
                                         "receipts")
        absent = [name for name in listed if name not in present]
        if absent:
            raise bad("/" + absent[0], f"{', '.join(absent)} {'is' if len(absent) == 1 else 'are'} listed in "
                                       f"{MANIFEST_NAME} but absent from {root}", "missing_asset")
        if len(receipts) > MAX_RECEIPTS:
            raise bad("/receipts", f"{len(receipts)} receipts exceed the limit of {MAX_RECEIPTS}", "resource_limit")

        # Declared and actual sizes, before any member is read.
        total = 0
        for member in manifest.members:
            limit = MAX_ASSET_BYTES if member.is_blob else _MEMBER_LIMITS[member.name]
            if member.bytes > limit:
                raise bad("/" + member.name, f"{member.name} is declared with {member.bytes} bytes; the limit for this "
                                             f"member is {limit}", "resource_limit")
            total += member.bytes
        if total > MAX_DOSSIER_BYTES:
            raise bad("/" + MANIFEST_NAME + "/members", f"the members total {total} bytes; a dossier is read up to "
                                                        f"{MAX_DOSSIER_BYTES}", "resource_limit")
        for member in manifest.members:
            if present[member.name] != member.bytes:
                raise bad("/" + member.name, f"{member.name} holds {present[member.name]} bytes; {MANIFEST_NAME} lists "
                                             f"{member.bytes}", "digest_mismatch")
        for name in receipts:
            if present[name] > MAX_RECEIPT_BYTES:
                raise bad("/" + name, f"{name} holds {present[name]} bytes; a receipt holds at most "
                                      f"{MAX_RECEIPT_BYTES}", "resource_limit")

        # Every member's digest.
        files, observed_digests = {}, {}
        for member in manifest.members:
            data = _read_file(root, member.name, member.bytes)
            observed = observed_digests[member.name] = _sha256(data)
            if len(data) != member.bytes or observed != member.sha256:
                raise bad("/" + member.name, f"{member.name} holds {len(data)} bytes with digest {observed}; "
                                             f"{MANIFEST_NAME} lists {member.bytes} bytes with {member.sha256}",
                          "digest_mismatch")
            files[member.name] = data

        # The lock and the recipe, and the manifest's commitments to them.
        lock = _document(LOCK_MEMBER, files[LOCK_MEMBER], SourceLock.from_json)
        recipe = _document(RECIPE_MEMBER, files[RECIPE_MEMBER], Recipe.from_json)
        try:
            check_manifest(manifest, lock)          # lock commitment, adapter, every locked blob, D input
        except ReferenceSourceError as exc:
            raise exc.under("/" + MANIFEST_NAME) from None
        found = recipe_hash(recipe)
        if found != manifest.hashes.recipe:
            raise bad("/" + RECIPE_MEMBER, f"{RECIPE_MEMBER} has recipe hash {found}, but {MANIFEST_NAME} commits to "
                                           f"{manifest.hashes.recipe}", "digest_mismatch")
        locked = {(a.format, a.sha256, a.bytes) for a in lock.assets}
        stored_receipts = []
        for name in receipts:
            data = _read_file(root, name, MAX_RECEIPT_BYTES)
            receipt = _document(name, data, AcquisitionReceipt.from_json)
            if (receipt.format, receipt.sha256, receipt.bytes) not in locked:
                raise bad("/" + name, f"{name} records {receipt.format} {receipt.sha256} ({receipt.bytes} bytes), "
                                      "which is not a locked asset", "identity_mismatch")
            stored_receipts.append(receipt)

        # Regenerate from the locked blobs and the recipe, and compare.
        assets = {a.asset_id: files[BLOB_PREFIX + a.blob_name] for a in lock.assets}
        try:
            projection = project(lock, assets, recipe)
        except ReferenceSourceError as exc:
            raise _regeneration_failed(exc, lock) from None
        regenerated = projection.inputs
        for name in _input_names(lock):
            stored = files[INPUTS_PREFIX + name]
            if regenerated.get(name) != stored:
                raise bad("/" + INPUTS_PREFIX + name, (
                    f"{INPUTS_PREFIX}{name} ({_sha256(stored)}) is not what the dossier's lock, blobs and recipe "
                    f"regenerate ({_sha256(regenerated.get(name, b''))})"), "digest_mismatch")

        try:
            stored = Reference.from_json(_text(REFERENCE_MEMBER, files[REFERENCE_MEMBER]))
        except ReferenceArtifactError as exc:
            raise bad("/" + REFERENCE_MEMBER, f"the stored reference artifact does not load: {exc}",
                      _ARTIFACT_KINDS.get(exc.kind, "malformed")) from exc
        if stored.to_json().encode("utf-8") != files[REFERENCE_MEMBER]:
            raise bad("/" + REFERENCE_MEMBER, f"{REFERENCE_MEMBER} is not in the artifact's own canonical form")
        stored_document = json.loads(stored.to_json())
        generator = stored_document.pop("generator")
        stored_document.pop("commitment")
        expected_document = json.loads(projection.reference.to_json())
        for key in ("generator", "commitment"):
            expected_document.pop(key)
        _compared(REFERENCE_MEMBER, stored_document, expected_document,
                  "the runtime reference the dossier's lock, blobs and recipe regenerate (every field but the "
                  "generator and the commitment)")

        # The report: the regenerated one, stated for the stored artifact and the recorded builder, must
        # be byte for byte the stored report.  Equal canonical bytes of a validated report prove the
        # stored report valid, so it is parsed only to locate a difference.
        expected = report_for(projection.report, stored, manifest.builder)
        expected_dict = expected.to_dict()
        if document_bytes(expected_dict) != files[REPORT_MEMBER]:
            try:
                stored_report = parse_json(_text(REPORT_MEMBER, files[REPORT_MEMBER]), max_bytes=MAX_REPORT_BYTES)
            except ReferenceSourceError as exc:
                raise exc.under("/" + REPORT_MEMBER) from None
            _compared(REPORT_MEMBER, stored_report, expected_dict,
                      "the report regenerated for the stored reference and the recorded builder")
            raise bad("/" + REPORT_MEMBER, f"{REPORT_MEMBER} is not in canonical form (canonical JSON, UTF-8, one LF); "
                                           "a dossier stores each generated document exactly as written")
        if generator != manifest.builder.to_dict():
            raise bad("/" + REFERENCE_MEMBER + "/generator", (
                f"the stored reference names its generator GenAIRR {show(generator.get('version'))}, but "
                f"{REPORT_MEMBER} and {MANIFEST_NAME} record the builder GenAIRR {manifest.builder.version}"),
                "identity_mismatch")
        try:
            check_manifest(manifest, lock, expected)
        except ReferenceSourceError as exc:
            raise exc.under("/" + MANIFEST_NAME) from None
        # The manifest of exactly the members this build has, with the digests observed above, and of the
        # expected report.
        members = ((name, observed_digests[name], len(files[name])) for name in _member_names(lock))
        _compared(MANIFEST_NAME, manifest.to_dict(), _manifest(lock, recipe, expected, members).to_dict(),
                  "the manifest of its stored members")
        return cls._assemble(lock=lock, recipe=recipe, reference=stored, report=expected, manifest=manifest,
                             files={name: files[name] for name in _member_names(lock)},
                             receipts=tuple(stored_receipts), checked_with=checking)

    @classmethod
    def bundled(cls, exact_id: str) -> "ReferenceBuild":
        """The reference build packaged with GenAIRR under exactly ``exact_id``.

        The packaged dossier's manifest must first have its registered manifest hash, so a
        modified installation fails at once; the dossier is then checked with
        :meth:`read_dossier` (offline, from local package resources).  There are no aliases and
        no downloads; an unknown identifier is an error listing the available ones.  A packaged
        build is checked once per process, concurrent first loads included, and the same build
        is returned afterwards.
        """
        entry = bundled.entry(exact_id)
        key = (exact_id, entry)
        with _BUNDLED_LOCK:
            cached = _BUNDLED.get(key)
            if cached is not None:
                return cached
            key_lock = _BUNDLED_LOCKS.setdefault(key, threading.Lock())
        with key_lock:
            with _BUNDLED_LOCK:
                cached = _BUNDLED.get(key)
            if cached is not None:
                return cached
            with bundled.dossier_directory(entry) as directory:
                found = manifest_hash(_read_manifest(Path(directory)))       # a modified install fails at once
                if found != entry.manifest_hash:
                    raise bad("/" + MANIFEST_NAME, f"the dossier packaged as {exact_id!r} has manifest {found}, not "
                                                   f"the registered {entry.manifest_hash}", "digest_mismatch")
                build = cls.read_dossier(directory)
            if manifest_hash(build.manifest) != entry.manifest_hash:
                raise bad("/" + MANIFEST_NAME, f"the dossier packaged as {exact_id!r} changed while it was checked",
                          "digest_mismatch")
            with _BUNDLED_LOCK:
                _BUNDLED[key] = build
            return build
