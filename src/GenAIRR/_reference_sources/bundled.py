"""The closed registry of reference build dossiers packaged with GenAIRR.

Each entry maps an exact identifier to a dossier packaged as a resource of a Python package
(a directory, or a single zip archive of its members) and to the manifest hash that dossier
must have.  There are no aliases, no pattern matching and no downloads: an identifier is
either listed here or unknown.
"""
from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass
import importlib.resources
import io
from pathlib import Path
import tempfile
import stat
import struct
from types import MappingProxyType
from typing import Callable, Iterator, Mapping, Optional
import unicodedata
import zipfile

from GenAIRR._reference_sources._fields import bad, check_digest, check_member_path
from GenAIRR._reference_sources.errors import ReferenceSourceError
from GenAIRR._reference_sources.manifest import (
    MANIFEST_NAME, MAX_MANIFEST_BYTES, MAX_MEMBER_NAME, MAX_MEMBERS, BuildManifest,
)

MAX_RESOURCE = 255
# A packaged dossier is copied out of a non-file package (for example a zip archive) within the
# same bounds the dossier reader enforces: entries and total bytes.
_MAX_ENTRIES = 2 * MAX_MEMBERS
_MAX_BYTES = 1 << 30
# A packaged archive: a resource named ``*.zip`` holding only the dossier's member files, stored
# without compression.  Its own size may exceed the member bytes by its headers and directory.
ARCHIVE_SUFFIX = ".zip"
_ARCHIVE_OVERHEAD = 64 << 20
_ZIP_ERRORS = (zipfile.BadZipFile, zipfile.LargeZipFile, EOFError, ValueError, NotImplementedError, RuntimeError)
# The end-of-central-directory record, and the zip64 locator that would precede it.  zipfile reads
# the whole central directory (its size taken from this record) before any entry can be checked,
# so the declared entry count and directory size are capped first.  An entry's central record is
# 46 bytes plus its name (at most MAX_MEMBER_NAME characters, 4 UTF-8 bytes each); packaged
# entries carry no extra field and no comment.
_EOCD = b"PK\x05\x06"
_ZIP64_LOCATOR = b"PK\x06\x07"
_MAX_RECORD = 46 + 4 * MAX_MEMBER_NAME
# Names a Windows filesystem cannot hold: device names (with any extension), and these characters.
_RESERVED = frozenset({"CON", "PRN", "AUX", "NUL", *(f"COM{i}" for i in range(1, 10)),
                       *(f"LPT{i}" for i in range(1, 10))})
_WINDOWS_CHARACTERS = frozenset('<>:"|?*')


@dataclass(frozen=True)
class BundledDossier:
    """A packaged dossier: the directory or ``.zip`` archive ``resource`` ('/'-separated) inside
    ``package``, and the ``genairr_reference_source_build_manifest_v1`` hash its manifest must have."""

    package: str
    resource: str
    manifest_hash: str

    def __post_init__(self):
        if type(self.package) is not str or not all(part.isidentifier() for part in self.package.split(".")):
            raise bad("/package", f"expected an importable package name, got {self.package!r}")
        check_member_path(self.resource, "/resource", MAX_RESOURCE)
        check_digest(self.manifest_hash, "/manifest_hash")


# Exact identifier -> packaged dossier.  Each is one zip archive written by
# ``tools/build_published_reference.py``: a single short package path, whatever the member names.
# The manifest hash covers the GenAIRR version that generated the dossier, so an archive
# regenerated after a version change is pinned again here.
REGISTRY: Mapping[str, BundledDossier] = MappingProxyType({
    # Human IGH from OGRDB IGH_VDJ release 10: 190 V, 33 D and 7 J alleles over their observed
    # coding spans, with source-annotated junction anchors.
    "human_igh_ogrdb_10_observed_v1": BundledDossier(
        "GenAIRR.data", "published_references/human_igh_ogrdb_10_observed_v1.zip",
        "sha256:7c306e4c3586c66639a5f37ecc56c9bcb1c75c1e297c3be465d51b22b186fafe"),
})


def entry(exact_id: str) -> BundledDossier:
    """The registry entry of exactly ``exact_id``; an unknown identifier lists the available ones."""
    if type(exact_id) is not str:
        raise TypeError(f"a bundled build identifier is a string, got {type(exact_id).__name__}")
    try:
        return REGISTRY[exact_id]
    except KeyError:
        available = ", ".join(sorted(REGISTRY)) or "none"
        raise bad("", f"no bundled reference build is named exactly {exact_id!r}; available: {available}",
                  "missing_asset") from None


def _entry_name(item, seen: set) -> str:
    """A packaged entry's name, safe as one path component and not repeated in its directory."""
    name = item.name
    if (not isinstance(name, str) or name in ("", ".", "..") or "/" in name or "\\" in name or ":" in name
            or any(ord(char) < 32 for char in name)):
        raise bad("", f"the packaged dossier holds an entry named {name!r}, which is not a plain file or directory "
                      "name", "malformed")
    if name in seen:
        raise bad("", f"the packaged dossier holds an entry named {name!r} twice", "malformed")
    seen.add(name)
    return name


def _copy(source, target: Path, budget: list) -> None:
    """Copy a packaged directory tree (any importlib resource) into ``target``, within bounds: no
    entry is read beyond the remaining byte budget."""
    target.mkdir()
    seen: set = set()
    for item in sorted(source.iterdir(), key=lambda x: str(x.name)):
        name = _entry_name(item, seen)
        budget[0] -= 1
        if budget[0] < 0:
            raise bad("", f"the packaged dossier holds more than {_MAX_ENTRIES} entries", "resource_limit")
        if item.is_dir():
            _copy(item, target / name, budget)
            continue
        with item.open("rb") as handle:
            data = handle.read(budget[1] + 1)
        if len(data) > budget[1]:
            raise bad("", f"the packaged dossier holds more than {_MAX_BYTES} bytes", "resource_limit")
        budget[1] -= len(data)
        (target / name).write_bytes(data)


def _check_directory_record(data: bytes, where: str) -> None:
    """Cap the entry count and central-directory size an archive declares, before zipfile reads it.

    The record is found as zipfile finds it (the last signature within the final 64 KiB); a zip64
    archive is refused, so no other count or size can be in force."""
    at = data.rfind(_EOCD, max(0, len(data) - 22 - 0xFFFF))
    if at < 0 or len(data) - at < 22:
        raise bad("", f"{where} is not a readable zip archive: it has no end-of-central-directory record")
    _, _, on_disk, total, size, offset, _ = struct.unpack("<HHHHIIH", data[at + 4:at + 22])
    if 0xFFFF in (on_disk, total) or 0xFFFFFFFF in (size, offset) or data[max(0, at - 20):at - 16] == _ZIP64_LOCATOR:
        raise bad("", f"{where} is a zip64 archive; a packaged archive is a plain zip")
    if total > _MAX_ENTRIES or on_disk > _MAX_ENTRIES:
        raise bad("", f"{where} declares {max(total, on_disk)} entries; it may hold at most {_MAX_ENTRIES}",
                  "resource_limit")
    if size > _MAX_ENTRIES * _MAX_RECORD:
        raise bad("", f"{where} declares a {size}-byte central directory; at most {_MAX_ENTRIES * _MAX_RECORD} bytes "
                      "fit its entries", "resource_limit")


def _component_problem(part: str) -> Optional[str]:
    """Why a name component cannot exist on every platform, or ``None``."""
    if part.endswith((".", " ")):
        return "ends with a dot or a space"
    if any(char in _WINDOWS_CHARACTERS for char in part):
        return f"contains one of {''.join(sorted(_WINDOWS_CHARACTERS))}"
    if part.split(".", 1)[0].rstrip(" ").upper() in _RESERVED:
        return "is a reserved device name on Windows"
    return None


def _fold(name: str) -> str:
    """The key under which two names would be one file somewhere: Unicode NFC, letter case ignored."""
    return unicodedata.normalize("NFC", unicodedata.normalize("NFC", name).casefold())


def _archive_entries(archive: zipfile.ZipFile, where: str) -> list:
    """The member entries of a packaged archive, refused unless every one is a safely named,
    unique, regular, plain (no extra field or comment), stored and unencrypted file within the
    entry and byte budgets, and no entry's name is another's directory."""
    try:
        entries = archive.infolist()
    except _ZIP_ERRORS as exc:
        raise bad("", f"{where} is not a readable zip archive: {exc}") from exc
    if len(entries) > _MAX_ENTRIES:
        raise bad("", f"{where} holds more than {_MAX_ENTRIES} entries", "resource_limit")
    names, folded, directories, total = set(), {}, {}, 0
    for entry_ in entries:
        name = entry_.orig_filename                     # as stored: zipfile cuts a NUL and, on Windows, maps '\\'
        if name != entry_.filename:
            raise bad("", f"{where} holds an entry stored as {name!r}, which is read back as "
                          f"{entry_.filename!r} (a NUL or a backslash in its name)")
        if name.endswith("/"):
            raise bad("", f"{where} holds a directory entry {name!r}; it holds only the dossier's member files")
        try:
            check_member_path(name, "", MAX_MEMBER_NAME)
        except ReferenceSourceError as exc:
            raise bad("", f"{where} holds an entry whose name is unsafe: {exc.message}", exc.kind) from None
        for part in name.split("/"):
            problem = _component_problem(part)
            if problem:
                raise bad("", f"{where} holds an entry named {name!r}, whose component {part!r} {problem}")
        key = _fold(name)
        if name in names or key in folded:
            raise bad("", f"{where} holds the entries {folded.get(key, name)!r} and {name!r}, one file once letter "
                          "case and Unicode normalisation are ignored")
        names.add(name)
        folded[key] = name
        parts = key.split("/")
        for depth in range(1, len(parts)):
            directories.setdefault("/".join(parts[:depth]), name)
        kind = stat.S_IFMT(entry_.external_attr >> 16)
        if kind not in (0, stat.S_IFREG):
            raise bad("", f"{where} holds the entry {name!r} with file type {kind:o}; it holds only regular files")
        if entry_.extra or entry_.comment:
            raise bad("", f"{where} holds the entry {name!r} with an extra field or a comment; its entries are plain")
        if entry_.flag_bits & 0x1:
            raise bad("", f"{where} holds the encrypted entry {name!r}")
        if entry_.compress_type != zipfile.ZIP_STORED or entry_.compress_size != entry_.file_size:
            raise bad("", f"{where} holds the compressed entry {name!r}; its entries are stored uncompressed")
        total += entry_.file_size
        if total > _MAX_BYTES:
            raise bad("", f"{where} holds more than {_MAX_BYTES} bytes", "resource_limit")
    clash = sorted(key for key in directories if key in folded)
    if clash:
        raise bad("", f"{where} holds {folded[clash[0]]!r} both as a file and as the directory of "
                      f"{directories[clash[0]]!r}")
    return entries


def _extract_entry(archive: zipfile.ZipFile, entry_: zipfile.ZipInfo, target: Path, where: str) -> bytes:
    try:
        with archive.open(entry_) as handle:
            data = handle.read(entry_.file_size + 1)       # reading to the end checks the CRC
    except _ZIP_ERRORS as exc:
        raise bad("", f"{where} entry {entry_.filename!r} cannot be read: {exc}") from exc
    if len(data) != entry_.file_size:
        raise bad("", f"{where} entry {entry_.filename!r} holds {len(data)} bytes, not the {entry_.file_size} its "
                      "header states")
    path = target.joinpath(*entry_.filename.split("/"))
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "xb") as handle:
        handle.write(data)
    return data


def _check_members(manifest: bytes, names: set, where: str) -> None:
    """The archive's entries must be exactly the members its manifest lists, and the manifest."""
    if len(manifest) > MAX_MANIFEST_BYTES:
        raise bad("/" + MANIFEST_NAME, f"{where} holds a {MANIFEST_NAME} of more than {MAX_MANIFEST_BYTES} bytes",
                  "resource_limit")
    try:
        listed = {member.name for member in BuildManifest.from_json(manifest.decode("utf-8")).members}
    except UnicodeDecodeError as exc:
        raise bad("/" + MANIFEST_NAME, f"{where} holds a {MANIFEST_NAME} that is not UTF-8: {exc}") from exc
    except ReferenceSourceError as exc:
        raise exc.under("/" + MANIFEST_NAME) from None
    listed.add(MANIFEST_NAME)
    extra, missing = sorted(names - listed), sorted(listed - names)
    if extra:
        raise bad("", f"{where} holds {', '.join(extra)}, which {MANIFEST_NAME} does not list; a packaged archive "
                      "holds exactly the dossier's members")
    if missing:
        raise bad("/" + missing[0], f"{where} lacks {', '.join(missing)}, which {MANIFEST_NAME} lists",
                  "missing_asset")


def _extract(resource, target: Path, where: str, check: Optional[Callable[[Path], None]]) -> None:
    """Extract a packaged archive into the new directory ``target``: ``manifest.json`` first, then
    ``check(target)`` (when given), then every other member."""
    limit = _MAX_BYTES + _ARCHIVE_OVERHEAD
    with resource.open("rb") as handle:
        data = handle.read(limit + 1)
    if len(data) > limit:
        raise bad("", f"{where} holds more than {limit} bytes", "resource_limit")
    _check_directory_record(data, where)
    try:
        archive = zipfile.ZipFile(io.BytesIO(data))
    except _ZIP_ERRORS as exc:
        raise bad("", f"{where} is not a readable zip archive: {exc}") from exc
    with archive:
        entries = _archive_entries(archive, where)
        target.mkdir()
        first = [e for e in entries if e.filename == MANIFEST_NAME]
        if not first:
            raise bad("", f"{where} holds no {MANIFEST_NAME}: it is not a complete dossier", "incomplete_dossier")
        manifest = _extract_entry(archive, first[0], target, where)
        if check is not None:
            check(target)
        _check_members(manifest, {e.filename for e in entries}, where)
        for entry_ in entries:
            if entry_.filename != MANIFEST_NAME:
                _extract_entry(archive, entry_, target, where)


@contextmanager
def dossier_directory(dossier: BundledDossier, check: Optional[Callable[[Path], None]] = None) -> Iterator[Path]:
    """A filesystem directory holding the packaged dossier, for reading only.

    Where the package's resources are plain files and the dossier is a directory, this is the
    packaged directory itself; otherwise the dossier is copied (a directory) or extracted (a
    ``.zip`` archive) into a private temporary directory, removed afterwards.  ``check``, when
    given, is called with the directory as soon as its ``manifest.json`` is there: for an
    archive, before any other member is extracted.  Nothing is ever written into the package.
    An absent dossier is ``missing_asset``; an unreadable, unsafe or compressed archive is
    ``malformed``.
    """
    try:
        resource = importlib.resources.files(dossier.package)
    except (ImportError, TypeError) as exc:
        raise bad("", f"the package {dossier.package!r} of the bundled dossier cannot be found: {exc}",
                  "missing_asset") from exc
    for part in dossier.resource.split("/"):         # one name per call, as resources join on Python 3.9
        resource = resource.joinpath(part)
    archived = dossier.resource.endswith(ARCHIVE_SUFFIX)
    if not (resource.is_file() if archived else resource.is_dir()):
        raise bad("", f"the bundled dossier {dossier.resource} is not packaged in {dossier.package}", "missing_asset")
    if isinstance(resource, Path) and not archived:
        if check is not None:
            check(resource)
        yield resource
        return
    with tempfile.TemporaryDirectory(prefix="genairr_bundled_") as temporary:
        target = Path(temporary) / "dossier"
        where = f"the bundled dossier archive {dossier.resource} of {dossier.package}"
        try:
            if archived:
                _extract(resource, target, where, check)
            else:
                _copy(resource, target, [_MAX_ENTRIES, _MAX_BYTES])
                if check is not None:
                    check(target)
        except OSError as exc:
            action = "extracted from its archive" if archived else "copied out of its package"
            raise bad("", f"the bundled dossier {dossier.resource} of {dossier.package} cannot be {action} into a "
                          f"private temporary directory: {exc}", "io") from exc
        yield target
