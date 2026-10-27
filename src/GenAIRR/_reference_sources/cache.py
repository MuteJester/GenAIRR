"""Content-addressed cache of exact source bytes: ``<cache>/blobs/sha256/<64 hex>``.

A blob's file name is the SHA-256 digest of its bytes and nothing else, so it never
derives from a provider or local file name.  Every read verifies the exact size and
digest, including blobs that were already present; a symbolic link or any non-regular
file at a blob path is refused.  Writes stream into a temporary file owned by the
writer, in the blob directory, and publish the complete bytes with the shared atomic
no-replace step: an existing identical blob (for example from a concurrent writer) is
accepted only after its full bytes are verified, and a conflicting file is reported and
never overwritten.  A failure or cancellation removes only the writer's own temporary
file.  Reading never creates directories and never acquires anything.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import os
from pathlib import Path
import stat
from typing import Optional

from GenAIRR import _artifact_io
from GenAIRR._reference_sources._fields import MAX_ASSET_BYTES, bad, check_digest, check_int

BLOB_DIRECTORY = ("blobs", "sha256")
_CHUNK = 1 << 16
_TEMPORARY_STEM = "incoming"
_READ_FLAGS = os.O_RDONLY | getattr(os, "O_BINARY", 0) | getattr(os, "O_NOFOLLOW", 0)


@dataclass(frozen=True)
class BlobRef:
    """The identity of one blob: ``sha256`` (``"sha256:<hex>"``) and size in bytes."""

    sha256: str
    bytes: int


def _io(message: str, exc: Optional[BaseException] = None):
    return bad("", message + (f" ({exc})" if exc is not None else ""), "io")


def _describe(asset: Optional[str], sha256: str, size: int) -> str:
    return f"{'asset ' + asset + ' ' if asset else 'blob '}({sha256}, {size} bytes)"


def _hash_open_file(handle, size: int, keep: bool = False):
    """Read at most ``size + 1`` bytes; return ``(digest, count, data or None)``."""
    hasher, count, parts = hashlib.sha256(), 0, []
    while count <= size:
        chunk = handle.read(min(_CHUNK, size + 1 - count))
        if not chunk:
            break
        hasher.update(chunk)
        count += len(chunk)
        if keep:
            parts.append(chunk)
    return "sha256:" + hasher.hexdigest(), count, (b"".join(parts) if keep else None)


class BlobCache:
    """A content-addressed blob store rooted at ``root`` (no I/O until used)."""

    def __init__(self, root):
        if isinstance(root, bool) or not isinstance(root, (str, os.PathLike)):
            raise TypeError(f"cache must be a filesystem path, got {type(root).__name__}")
        root = Path(root)
        if not isinstance(os.fspath(root), str):
            raise TypeError("cache must be a text filesystem path")
        self.root = root
        self.blob_dir = root.joinpath(*BLOB_DIRECTORY)

    def __repr__(self) -> str:
        return f"BlobCache({str(self.root)!r})"

    def blob_path(self, sha256: str) -> Path:
        """The only place a blob with digest ``sha256`` may live."""
        check_digest(sha256, "")
        return self.blob_dir / sha256.split(":", 1)[1]

    # -- verified reads ----------------------------------------------

    def _check(self, sha256: str, size: int, asset: Optional[str], *, keep: bool):
        """Verify the blob; ``None`` when absent, else its bytes (``b""`` unless ``keep``)."""
        path = self.blob_path(sha256)
        check_int(size, "", minimum=1, maximum=MAX_ASSET_BYTES, limit_kind="resource_limit")
        what = _describe(asset, sha256, size)
        try:
            status = os.lstat(path)
        except FileNotFoundError:
            return None
        except OSError as exc:
            raise _io(f"cannot inspect the cached {what} at {path}", exc) from exc
        if stat.S_ISLNK(status.st_mode):
            raise _io(f"refusing a symbolic link at the cache path of {what}: {path}")
        if not stat.S_ISREG(status.st_mode):
            raise _io(f"refusing the cache path of {what}: {path} is not a regular file")
        try:
            descriptor = os.open(path, _READ_FLAGS)
        except FileNotFoundError:
            return None
        except OSError as exc:
            raise _io(f"cannot open the cached {what} at {path}", exc) from exc
        try:
            with os.fdopen(descriptor, "rb") as handle:
                opened = os.fstat(handle.fileno())
                if not stat.S_ISREG(opened.st_mode):
                    raise _io(f"refusing the cache path of {what}: {path} is not a regular file")
                observed_size = opened.st_size
                if observed_size != size:
                    raise bad("", f"the cached {what} at {path} is corrupt: it holds {observed_size} bytes, "
                                  f"expected {size} bytes of {sha256}; it is never overwritten: remove it "
                                  "and acquire the asset again", "digest_mismatch")
                observed, count, data = _hash_open_file(handle, size, keep)
        except OSError as exc:
            raise _io(f"cannot read the cached {what} at {path}", exc) from exc
        if count != size or observed != sha256:
            raise bad("", f"the cached {what} at {path} is corrupt: observed {observed} ({count} bytes), "
                          f"expected {sha256} ({size} bytes); it is never overwritten: remove it and acquire "
                          "the asset again", "digest_mismatch")
        return data if keep else b""

    def _missing(self, sha256: str, size: int, asset: Optional[str]):
        return bad("", f"{_describe(asset, sha256, size)} is not in the cache {self.root} (expected at "
                       f"{self.blob_path(sha256)}); reading never downloads: fetch or import the source again "
                       "with this cache", "missing_asset")

    def read(self, sha256: str, size: int, *, asset: Optional[str] = None) -> bytes:
        """The blob's exact bytes after verifying size and digest.

        Absent: ``missing_asset`` naming ``asset`` and the expected digest.  Corrupt:
        ``digest_mismatch`` with observed and expected values.  A symbolic link or a
        non-regular file: ``io``.
        """
        data = self._check(sha256, size, asset, keep=True)
        if data is None:
            raise self._missing(sha256, size, asset)
        return data

    def verify(self, sha256: str, size: int, *, asset: Optional[str] = None) -> None:
        """Like :meth:`read` without keeping the bytes."""
        if self._check(sha256, size, asset, keep=False) is None:
            raise self._missing(sha256, size, asset)

    def contains(self, sha256: str, size: int) -> bool:
        """True if a blob with exactly these bytes is present (fully verified), False if absent.

        A corrupt, symbolic-link or non-regular entry raises, as for :meth:`read`.
        """
        return self._check(sha256, size, None, keep=False) is not None

    # -- writes -------------------------------------------------------

    def writer(self) -> "BlobWriter":
        """A new owned temporary file in the blob directory (created on demand)."""
        return BlobWriter(self)

    def put(self, data: bytes) -> BlobRef:
        """Stage and publish ``data``; return its :class:`BlobRef`."""
        if not isinstance(data, (bytes, bytearray, memoryview)):
            raise TypeError(f"data must be bytes, got {type(data).__name__}")
        with self.writer() as writer:
            writer.write(bytes(data))
            staged = writer.finish()
        self.publish(staged)
        return BlobRef(staged.sha256, staged.bytes)

    def publish(self, staged: "StagedBlob") -> Path:
        """Publish a staged blob at its content address without ever replacing a file.

        An existing identical blob is accepted after full verification; a conflicting
        file is ``digest_mismatch`` and stays in place.  The staged temporary file is
        always consumed or removed.
        """
        if not isinstance(staged, StagedBlob) or staged.cache_dir != self.blob_dir:
            raise TypeError("publish takes a StagedBlob written by this cache")
        destination = self.blob_path(staged.sha256)
        try:
            if self.contains(staged.sha256, staged.bytes):
                return destination
            try:
                _artifact_io._publish_without_replacing(staged.path, destination)
            except FileExistsError:
                if self.contains(staged.sha256, staged.bytes):
                    return destination
                raise _io(f"the cache entry {destination} changed while publishing; try again") from None
            except _artifact_io.ArtifactIOError as exc:
                raise _io(f"cannot publish blob {staged.sha256} in the cache {self.root}: the cache filesystem "
                          "does not support atomic no-overwrite publication (hard links are unavailable); "
                          "place the cache on a local disk") from exc
            except OSError as exc:
                raise _io(f"cannot publish blob {staged.sha256} at {destination}", exc) from exc
            staged._published()
            _sync_directory(self.blob_dir)
            return destination
        finally:
            staged.discard()


def _sync_directory(directory: Path) -> None:
    if not hasattr(os, "O_DIRECTORY"):
        return
    try:
        descriptor = os.open(directory, os.O_RDONLY | os.O_DIRECTORY)
        try:
            os.fsync(descriptor)
        finally:
            os.close(descriptor)
    except OSError as exc:
        raise _io(f"cannot synchronize the cache directory {directory}", exc) from exc


class StagedBlob:
    """A complete, fsynced temporary file awaiting publication; owns only its own path."""

    def __init__(self, cache_dir: Path, path: Path, sha256: str, size: int):
        self.cache_dir = cache_dir
        self.path = path
        self.sha256 = sha256
        self.bytes = size
        self._owned = True

    @property
    def ref(self) -> BlobRef:
        return BlobRef(self.sha256, self.bytes)

    def _published(self) -> None:
        self._owned = False

    def discard(self) -> None:
        """Remove this blob's own temporary file if it is still there (idempotent)."""
        if self._owned:
            self._owned = False
            try:
                os.unlink(self.path)
            except FileNotFoundError:
                pass

    def __repr__(self) -> str:
        return f"StagedBlob({self.sha256}, {self.bytes} bytes)"


class BlobWriter:
    """Streams bytes into an owned temporary file, hashing as it goes.

    Use as a context manager: leaving the block without :meth:`finish` (an error or a
    cancellation) removes the temporary file.  At most ``MAX_ASSET_BYTES`` are accepted.
    """

    def __init__(self, cache: BlobCache):
        self._cache_dir = cache.blob_dir
        try:
            self._cache_dir.mkdir(parents=True, exist_ok=True)
            descriptor, self._path = _artifact_io._create_temporary(self._cache_dir / _TEMPORARY_STEM)
        except OSError as exc:
            raise _io(f"cannot create a temporary file in the cache directory {self._cache_dir}", exc) from exc
        self._file = os.fdopen(descriptor, "wb")
        self._hasher = hashlib.sha256()
        self.size = 0
        self._finished = False
        self._discarded = False

    def __enter__(self) -> "BlobWriter":
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self.discard()

    def write(self, chunk: bytes) -> None:
        if self._finished or self._discarded:
            raise ValueError("this writer is finished")
        if self.size + len(chunk) > MAX_ASSET_BYTES:
            raise bad("", f"a source asset of more than {MAX_ASSET_BYTES} bytes ({self.size + len(chunk)} "
                          f"bytes so far) exceeds the per-asset limit of {MAX_ASSET_BYTES} bytes",
                      "resource_limit")
        try:
            self._file.write(chunk)
        except OSError as exc:
            raise _io(f"cannot write the temporary file {self._path}", exc) from exc
        self._hasher.update(chunk)
        self.size += len(chunk)

    def finish(self) -> StagedBlob:
        """Flush, fsync and close; return the staged blob (which now owns the file)."""
        if self._finished or self._discarded:
            raise ValueError("this writer is finished")
        if self.size == 0:
            raise bad("", "a source asset must not be empty (0 bytes)")
        try:
            self._file.flush()
            os.fsync(self._file.fileno())
            self._file.close()
        except OSError as exc:
            raise _io(f"cannot write the temporary file {self._path}", exc) from exc
        self._finished = True
        return StagedBlob(self._cache_dir, self._path, "sha256:" + self._hasher.hexdigest(), self.size)

    def discard(self) -> None:
        """Close and remove this writer's own temporary file unless finished (idempotent)."""
        if self._finished or self._discarded:
            return
        self._discarded = True
        try:
            self._file.close()
        except OSError:
            pass
        try:
            os.unlink(self._path)
        except FileNotFoundError:
            pass
