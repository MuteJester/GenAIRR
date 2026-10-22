"""Shared strict JSON and atomic file I/O for portable GenAIRR artifacts.

Canonical serialization, typed digests, a bounded duplicate-key-rejecting
parser (native) and no-clobber atomic writes live here so every artifact
format uses one implementation.
"""
from __future__ import annotations

import errno
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
    if text.startswith("﻿"):
        text = text[1:]
    try:
        normalized = _engine._family_replay_parse_json(text, json.dumps({"max_file_bytes": max_bytes}))
    except Exception as exc:
        message = str(exc)
        kind = "resource" if "resource" in message.lower() else "format"
        raise ArtifactIOError(kind, message) from exc
    return json.loads(normalized)


def read_text_bounded(path, *, max_bytes: int) -> str:
    try:
        with open(path, "rb") as handle:
            data = handle.read(max_bytes + 1)
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


def write_atomic(path, data: bytes, *, replace: bool) -> None:
    """Publish atomically, refusing to replace an existing file unless ``replace``.

    The published file's permissions follow the user's umask.
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
            os.link(temporary, destination)
            os.unlink(temporary)
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
