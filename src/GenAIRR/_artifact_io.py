"""Shared strict JSON and atomic file I/O for portable GenAIRR artifacts.

Canonical serialization, typed digests, a bounded duplicate-key-rejecting
parser (native) and no-clobber atomic writes live here so every artifact
format uses one implementation.
"""
from __future__ import annotations

import errno
import io
import json
import os
from pathlib import Path
from secrets import token_hex
from typing import Any

from GenAIRR import _engine


class ArtifactIOError(ValueError):
    """An artifact read/write failure; ``kind`` is "io", "format" or "resource"."""

    def __init__(self, kind: str, message: str):
        self.kind = kind
        super().__init__(message)


def canonical_json(value) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"), allow_nan=False)


def typed_digest(value, domain: str) -> str:
    return _engine._family_replay_digest(canonical_json(value), domain)


def strict_loads(text: str, *, max_bytes: int) -> Any:
    if not isinstance(text, str):
        raise ArtifactIOError("format", "expected JSON text")
    if text.startswith("\ufeff"):
        text = text[1:]
    try:
        normalized = _engine._family_replay_parse_json(text, json.dumps({"max_file_bytes": max_bytes}))
    except Exception as exc:
        message = str(exc)
        kind = "resource" if "resource" in message.lower() else "format"
        raise ArtifactIOError(kind, message) from exc
    return json.loads(normalized)


_READ_CHUNK = 1 << 20


def read_bounded(handle, max_bytes: int) -> bytes:
    """At most ``max_bytes + 1`` bytes of ``handle``, allocated as they arrive.

    The first read asks for the file's own size when it is known (one more byte shows that it
    grew), later reads at most 1 MiB, so a small file never costs a buffer the size of the
    limit.  A result longer than ``max_bytes`` means the source is over the limit.
    """
    try:
        size = os.fstat(handle.fileno()).st_size
    except (AttributeError, OSError, ValueError, io.UnsupportedOperation):
        size = 0
    want = min(size, max_bytes) + 1 if size > 0 else _READ_CHUNK
    parts, total = [], 0
    while total <= max_bytes:
        asked = min(want, max_bytes + 1 - total)
        chunk = handle.read(asked)
        if not chunk:
            break
        parts.append(chunk)
        total += len(chunk)
        if size > 0 and len(parts) == 1 and len(chunk) < asked:
            break                                   # the file's own size, and no more: it is complete
        want = _READ_CHUNK
    return parts[0] if len(parts) == 1 else b"".join(parts)


def read_text_bounded(path, *, max_bytes: int) -> str:
    try:
        with open(path, "rb") as handle:
            data = read_bounded(handle, max_bytes)
    except OSError as exc:
        raise ArtifactIOError("io", str(exc)) from exc
    if len(data) > max_bytes:
        raise ArtifactIOError("resource", "file byte limit exceeded")
    try:
        return data.decode("utf-8")
    except UnicodeError as exc:
        raise ArtifactIOError("format", str(exc)) from exc


_TEMPORARY_ATTEMPTS = 100


def _create_temporary(destination: Path):
    """Create ``.<name>.<random>.tmp`` beside ``destination``; return ``(descriptor, path)``.

    The file is created with mode 0o666 so the user's umask applies, as for any
    ordinary file write (``tempfile.mkstemp`` would force 0o600).  ``O_EXCL`` never
    opens an existing file, so a name already taken is skipped, never reused.
    """
    flags = os.O_CREAT | os.O_EXCL | os.O_WRONLY | getattr(os, "O_BINARY", 0)
    for _ in range(_TEMPORARY_ATTEMPTS):
        candidate = destination.parent / f".{destination.name}.{token_hex(8)}.tmp"
        try:
            return os.open(candidate, flags, 0o666), candidate
        except FileExistsError:
            continue
    raise FileExistsError(errno.EEXIST, "no unused temporary file name", str(destination.parent))


# ``os.link`` errors meaning the filesystem has no hard links (FAT/exFAT, some SMB shares).
_NO_HARD_LINK_ERRNOS = frozenset(
    getattr(errno, name) for name in ("EPERM", "ENOTSUP", "EOPNOTSUPP", "ENOSYS", "EXDEV")
    if hasattr(errno, name))

# Windows ``os.rename`` refuses to overwrite an existing file: an atomic no-replace publish.
_WINDOWS = os.name == "nt"

# Windows ``os.link`` errors meaning no hard links: ERROR_INVALID_FUNCTION (FAT/exFAT) and
# ERROR_NOT_SUPPORTED (many SMB shares).  CPython maps both to EINVAL, which elsewhere
# means something else, so they count only on Windows.
_WINDOWS_NO_HARD_LINK_ERRORS = frozenset({1, 50})


def _hard_links_unsupported(exc: OSError) -> bool:
    if isinstance(exc, FileExistsError):
        return False
    if _WINDOWS and getattr(exc, "winerror", None) in _WINDOWS_NO_HARD_LINK_ERRORS:
        return True
    return isinstance(exc, PermissionError) or exc.errno in _NO_HARD_LINK_ERRNOS


def _publish_without_replacing(temporary: Path, destination: Path) -> None:
    """Publish the complete temporary file at ``destination`` only if nothing is there.

    The destination is never created before this single atomic step: a hard link, or on
    Windows a rename.  Elsewhere, a filesystem without hard links offers no atomic
    no-overwrite step, so the write fails instead.
    """
    try:
        os.link(temporary, destination)
    except OSError as exc:
        if not _hard_links_unsupported(exc):
            raise
        if not _WINDOWS:
            raise ArtifactIOError(
                "io", f"cannot write {destination.name}: the destination filesystem does not support "
                      f"atomic no-overwrite publication (hard links are unavailable: {exc}); write to a "
                      "local disk and copy the file, or pass replace=True if overwriting is intended") from exc
        os.rename(temporary, destination)
    else:
        os.unlink(temporary)


def write_atomic(path, data: bytes, *, replace: bool) -> None:
    """Publish atomically, refusing to replace an existing file unless ``replace``.

    The published file's permissions follow the user's umask.  Without ``replace`` the
    complete file is published by a hard link (by a rename on Windows); where that is
    impossible the write fails with ``ArtifactIOError("io")`` and the destination is left
    untouched.  Only this call's temporary file is ever removed.
    """
    if type(replace) is not bool:
        raise TypeError("replace must be Boolean")
    destination = Path(path)
    temporary = None
    try:
        descriptor, temporary = _create_temporary(destination)
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(data)
            handle.flush()
            os.fsync(handle.fileno())
        if replace:
            os.replace(temporary, destination)
        else:
            _publish_without_replacing(temporary, destination)
        temporary = None
        if hasattr(os, "O_DIRECTORY"):
            directory = os.open(destination.parent, os.O_RDONLY | os.O_DIRECTORY)
            try:
                os.fsync(directory)
            finally:
                os.close(directory)
    except OSError as exc:
        raise ArtifactIOError("io", str(exc)) from exc
    finally:
        if temporary is not None:
            try:
                os.unlink(temporary)
            except FileNotFoundError:
                pass
