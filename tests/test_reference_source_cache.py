"""Content-addressed blob cache for published germline-source bytes.

Blobs live at ``<cache>/blobs/sha256/<64 hex>``.  Every read verifies the exact size
and SHA-256 digest, symbolic links and non-regular files are refused, writes stream
into an owned temporary file beside the blobs and publish without ever replacing an
existing file, and a failure or cancellation removes only that temporary file.
"""
from __future__ import annotations

import errno
import hashlib
import os
import pathlib
import threading

import pytest

from GenAIRR import _artifact_io
from GenAIRR._reference_sources import cache as cache_module
from GenAIRR._reference_sources.cache import BLOB_DIRECTORY, BlobCache, BlobRef, StagedBlob
from GenAIRR._reference_sources.errors import ReferenceSourceError
from GenAIRR._reference_sources.model import MAX_ASSET_BYTES

ROOT = pathlib.Path(__file__).resolve().parent / "data" / "reference_sources" / "ogrdb"
GAPPED = (ROOT / "9606.IGH_VDJ" / "10" / "gapped.fasta").read_bytes()
GAPPED_SHA = "sha256:681db69932757bd09253cbef990a64959f8c8745e46f364bcb35963ea4f16a65"
GAPPED_SIZE = 72614


def _sha(data: bytes) -> str:
    return "sha256:" + hashlib.sha256(data).hexdigest()


def _error(kind, call, *args, **kwargs) -> ReferenceSourceError:
    with pytest.raises(ReferenceSourceError) as info:
        call(*args, **kwargs)
    assert info.value.kind == kind, f"expected {kind!r}, got {info.value.kind!r}: {info.value}"
    return info.value


def _temporaries(cache: BlobCache) -> list:
    if not cache.blob_dir.exists():
        return []
    return sorted(p.name for p in cache.blob_dir.iterdir() if p.name.endswith(".tmp"))


def _snapshot(directory: pathlib.Path) -> dict:
    if not directory.exists():
        return {}
    return {str(p.relative_to(directory)): (p.read_bytes() if p.is_file() else None)
            for p in sorted(directory.rglob("*"))}


@pytest.fixture
def cache(tmp_path) -> BlobCache:
    return BlobCache(tmp_path / "cache")


def _symlink_or_skip(target: pathlib.Path, link: pathlib.Path, *, directory: bool = False) -> None:
    try:
        os.symlink(target, link, target_is_directory=directory)
    except (OSError, NotImplementedError) as exc:  # e.g. Windows without the symlink privilege
        pytest.skip(f"symbolic links unavailable: {exc}")


# ── layout ────────────────────────────────────────────────────────────


def test_blobs_are_content_addressed_under_blobs_sha256(tmp_path):
    cache = BlobCache(tmp_path / "c")
    assert BLOB_DIRECTORY == ("blobs", "sha256")
    assert cache.root == tmp_path / "c"
    assert cache.blob_dir == tmp_path / "c" / "blobs" / "sha256"
    assert cache.blob_path(GAPPED_SHA) == cache.blob_dir / GAPPED_SHA.split(":")[1]
    assert BlobCache(str(tmp_path / "c")).blob_dir == cache.blob_dir


@pytest.mark.parametrize("value", ["sha256:" + "A" * 64, "sha256:" + "0" * 63, "0" * 64, "md5:" + "0" * 32,
                                   "sha256:../" + "0" * 61, None, 7])
def test_blob_names_come_only_from_validated_digests(cache, value):
    _error("malformed", cache.blob_path, value)


@pytest.mark.parametrize("root", [None, 3, b"bytes-path", True])
def test_cache_root_must_be_a_path(root):
    with pytest.raises(TypeError):
        BlobCache(root)


def test_constructing_or_reading_a_cache_never_creates_it(tmp_path):
    cache = BlobCache(tmp_path / "absent")
    _error("missing_asset", cache.read, GAPPED_SHA, GAPPED_SIZE)
    assert cache.contains(GAPPED_SHA, GAPPED_SIZE) is False
    assert not (tmp_path / "absent").exists()


# ── put, read and verify ──────────────────────────────────────────────


def test_put_publishes_exact_bytes_and_read_returns_them(cache):
    ref = cache.put(GAPPED)
    assert ref == BlobRef(GAPPED_SHA, GAPPED_SIZE)
    assert cache.blob_path(GAPPED_SHA).read_bytes() == GAPPED
    assert cache.read(GAPPED_SHA, GAPPED_SIZE) == GAPPED
    cache.verify(GAPPED_SHA, GAPPED_SIZE)
    assert cache.contains(GAPPED_SHA, GAPPED_SIZE) is True
    assert _temporaries(cache) == []


def test_putting_identical_bytes_twice_keeps_one_verified_blob(cache):
    assert cache.put(GAPPED) == cache.put(GAPPED)
    assert sorted(p.name for p in cache.blob_dir.iterdir()) == [GAPPED_SHA.split(":")[1]]


def test_a_missing_blob_names_the_asset_and_expected_digest(cache):
    error = _error("missing_asset", cache.read, GAPPED_SHA, GAPPED_SIZE, asset="gapped")
    assert "gapped" in error.message and GAPPED_SHA in error.message and str(GAPPED_SIZE) in error.message
    _error("missing_asset", cache.verify, GAPPED_SHA, GAPPED_SIZE, asset="gapped")


def _flip_first_byte(data: bytes) -> bytes:
    return bytes([data[0] ^ 1]) + data[1:]


CORRUPTIONS = {
    "same-size-different-bytes": _flip_first_byte,
    "truncated": lambda data: data[:-1],
    "extended": lambda data: data + b"\n",
    "emptied": lambda data: b"",
}


@pytest.mark.parametrize("corrupt", list(CORRUPTIONS.values()), ids=list(CORRUPTIONS))
def test_a_corrupt_cache_blob_is_rejected_on_every_read(cache, corrupt):
    cache.put(GAPPED)
    path = cache.blob_path(GAPPED_SHA)
    path.write_bytes(corrupt(GAPPED))
    for _ in range(2):  # pre-existing blobs are verified every time, not once
        error = _error("digest_mismatch", cache.read, GAPPED_SHA, GAPPED_SIZE, asset="gapped")
        assert "gapped" in error.message and GAPPED_SHA in error.message
        assert str(GAPPED_SIZE) in error.message and str(path) in error.message
        _error("digest_mismatch", cache.verify, GAPPED_SHA, GAPPED_SIZE)
        _error("digest_mismatch", cache.contains, GAPPED_SHA, GAPPED_SIZE)
    assert path.read_bytes() == corrupt(GAPPED)  # reported, never repaired silently


def test_a_same_size_corruption_reports_the_observed_digest(cache):
    cache.put(GAPPED)
    corrupted = _flip_first_byte(GAPPED)
    cache.blob_path(GAPPED_SHA).write_bytes(corrupted)
    error = _error("digest_mismatch", cache.read, GAPPED_SHA, GAPPED_SIZE)
    assert _sha(corrupted) in error.message


def test_a_symbolic_link_at_a_blob_path_is_refused_even_to_the_right_bytes(cache, tmp_path):
    elsewhere = tmp_path / "elsewhere.fasta"
    elsewhere.write_bytes(GAPPED)
    cache.blob_dir.mkdir(parents=True)
    _symlink_or_skip(elsewhere, cache.blob_path(GAPPED_SHA))
    for call in (cache.read, cache.verify, cache.contains):
        error = _error("io", call, GAPPED_SHA, GAPPED_SIZE)
        assert "symbolic link" in error.message


def test_a_non_regular_file_at_a_blob_path_is_refused(cache):
    cache.blob_path(GAPPED_SHA).mkdir(parents=True)
    for call in (cache.read, cache.verify, cache.contains):
        error = _error("io", call, GAPPED_SHA, GAPPED_SIZE)
        assert "regular file" in error.message


@pytest.mark.parametrize("size", [0, -1, MAX_ASSET_BYTES + 1, True, 1.0])
def test_expected_sizes_are_checked_before_reading(cache, size):
    kind = "resource_limit" if size == MAX_ASSET_BYTES + 1 and type(size) is int else "malformed"
    _error(kind, cache.read, GAPPED_SHA, size)


# ── streaming writes ─────────────────────────────────────────────────


def test_a_writer_streams_into_an_owned_temporary_file_beside_the_blobs(cache):
    with cache.writer() as writer:
        (temporary,) = _temporaries(cache)
        assert temporary.startswith(".incoming.")
        for start in range(0, len(GAPPED), 1000):
            writer.write(GAPPED[start:start + 1000])
        staged = writer.finish()
    assert isinstance(staged, StagedBlob)
    assert (staged.sha256, staged.bytes) == (GAPPED_SHA, GAPPED_SIZE)
    assert staged.path.parent == cache.blob_dir and staged.path.name == temporary
    assert not cache.blob_path(GAPPED_SHA).exists()   # nothing is published before validation
    assert cache.publish(staged) == cache.blob_path(GAPPED_SHA)
    assert cache.read(GAPPED_SHA, GAPPED_SIZE) == GAPPED
    assert _temporaries(cache) == []


def test_discarding_a_staged_blob_removes_only_its_own_temporary_file(cache):
    cache.blob_dir.mkdir(parents=True)
    other = cache.blob_dir / ".incoming.0123456789abcdef.tmp"   # another writer's temporary file
    other.write_bytes(b"not ours")
    with cache.writer() as writer:
        writer.write(GAPPED)
        staged = writer.finish()
    staged.discard()
    staged.discard()   # idempotent
    assert _temporaries(cache) == [other.name] and other.read_bytes() == b"not ours"
    assert not cache.blob_path(GAPPED_SHA).exists()


@pytest.mark.parametrize("interruption", [RuntimeError("failure"), KeyboardInterrupt()],
                         ids=["failure", "cancellation"])
def test_failure_or_cancellation_while_streaming_removes_the_temporary_file(cache, interruption):
    with pytest.raises(type(interruption)):
        with cache.writer() as writer:
            writer.write(GAPPED[:100])
            raise interruption
    assert _temporaries(cache) == []
    assert list(cache.blob_dir.iterdir()) == []


def test_an_empty_blob_is_malformed_and_leaves_nothing(cache):
    with cache.writer() as writer:
        _error("malformed", writer.finish)
    assert _temporaries(cache) == []
    _error("malformed", cache.put, b"")


def test_a_writer_refuses_more_than_one_asset_of_bytes(cache, monkeypatch):
    monkeypatch.setattr(cache_module, "MAX_ASSET_BYTES", 10)
    with cache.writer() as writer:
        writer.write(b"0123456789")
        error = _error("resource_limit", writer.write, b"x")
        assert "11" in error.message and "10" in error.message
    assert _temporaries(cache) == []


def test_the_writer_limit_is_the_shared_asset_budget():
    assert cache_module.MAX_ASSET_BYTES == MAX_ASSET_BYTES == 16 * 1024 * 1024


# ── publication: existing, conflicting and concurrent blobs ───────────


def test_an_existing_identical_blob_is_accepted_only_after_full_verification(cache, monkeypatch):
    cache.put(GAPPED)
    hashed = []
    real = cache_module._hash_open_file

    def spy(handle, size, *args, **kwargs):
        hashed.append(size)
        return real(handle, size, *args, **kwargs)

    monkeypatch.setattr(cache_module, "_hash_open_file", spy)
    with cache.writer() as writer:
        writer.write(GAPPED)
        staged = writer.finish()
    assert cache.publish(staged) == cache.blob_path(GAPPED_SHA)
    assert hashed == [GAPPED_SIZE]          # the existing file was read in full
    assert _temporaries(cache) == []


@pytest.mark.parametrize("corrupt", list(CORRUPTIONS.values()), ids=list(CORRUPTIONS))
def test_a_conflicting_existing_file_is_an_error_and_is_never_overwritten(cache, corrupt):
    cache.blob_dir.mkdir(parents=True)
    path = cache.blob_path(GAPPED_SHA)
    path.write_bytes(corrupt(GAPPED))
    with cache.writer() as writer:
        writer.write(GAPPED)
        staged = writer.finish()
    error = _error("digest_mismatch", cache.publish, staged)
    assert str(path) in error.message and "never overwritten" in error.message
    assert path.read_bytes() == corrupt(GAPPED)
    assert _temporaries(cache) == []


def test_a_symbolic_link_at_the_destination_is_never_replaced(cache, tmp_path):
    elsewhere = tmp_path / "elsewhere"
    elsewhere.write_bytes(b"someone else's file")
    cache.blob_dir.mkdir(parents=True)
    _symlink_or_skip(elsewhere, cache.blob_path(GAPPED_SHA))
    _error("io", cache.put, GAPPED)
    assert os.path.islink(cache.blob_path(GAPPED_SHA))
    assert elsewhere.read_bytes() == b"someone else's file"
    assert _temporaries(cache) == []


def _another_writer_publishes_before_link(monkeypatch, cache, theirs: bytes):
    """Patch ``os.link`` so another writer completes the blob just before this one links."""
    real_link, calls = os.link, []

    def link(source, destination, *args, **kwargs):
        calls.append(pathlib.Path(destination))
        staged = cache.root / "theirs.staged"
        staged.write_bytes(theirs)
        os.replace(staged, destination)
        return real_link(source, destination, *args, **kwargs)

    monkeypatch.setattr(os, "link", link)
    return calls


def test_a_concurrent_identical_writer_winning_the_race_is_accepted(cache, monkeypatch):
    cache.blob_dir.mkdir(parents=True)
    links = _another_writer_publishes_before_link(monkeypatch, cache, GAPPED)
    assert cache.put(GAPPED) == BlobRef(GAPPED_SHA, GAPPED_SIZE)
    assert links == [cache.blob_path(GAPPED_SHA)]
    assert cache.blob_path(GAPPED_SHA).read_bytes() == GAPPED
    assert _temporaries(cache) == []


def test_a_concurrent_conflicting_writer_is_an_error_and_its_file_survives(cache, monkeypatch):
    cache.blob_dir.mkdir(parents=True)
    _another_writer_publishes_before_link(monkeypatch, cache, b"conflicting bytes")
    error = _error("digest_mismatch", cache.put, GAPPED)
    assert "never overwritten" in error.message
    assert cache.blob_path(GAPPED_SHA).read_bytes() == b"conflicting bytes"
    assert _temporaries(cache) == []


def test_two_threads_publishing_the_same_blob_both_succeed_with_one_blob(cache):
    for _ in range(10):
        barrier = threading.Barrier(2)
        results, errors = [], []
        staged = []
        for _ in range(2):
            with cache.writer() as writer:
                writer.write(GAPPED)
                staged.append(writer.finish())

        def publish(blob):
            try:
                barrier.wait(timeout=10)
                results.append(cache.publish(blob))
            except BaseException as exc:  # reported below
                errors.append(exc)

        threads = [threading.Thread(target=publish, args=(blob,)) for blob in staged]
        for thread in threads:
            thread.start()
        for thread in threads:
            thread.join(timeout=30)
        assert errors == []
        assert results == [cache.blob_path(GAPPED_SHA)] * 2
        assert sorted(p.name for p in cache.blob_dir.iterdir()) == [GAPPED_SHA.split(":")[1]]
        assert cache.read(GAPPED_SHA, GAPPED_SIZE) == GAPPED
        cache.blob_path(GAPPED_SHA).unlink()


# ── filesystems without hard links ───────────────────────────────────


def _no_hard_links(monkeypatch):
    calls = []

    def link(source, destination, *args, **kwargs):
        calls.append(pathlib.Path(destination))
        raise PermissionError(errno.EPERM, "Operation not permitted")

    monkeypatch.setattr(os, "link", link)
    return calls


def _recorded_opens(monkeypatch):
    real, calls = os.open, []

    def wrapper(path, *args, **kwargs):
        calls.append(pathlib.Path(path))
        return real(path, *args, **kwargs)

    monkeypatch.setattr(os, "open", wrapper)
    return calls


@pytest.fixture
def posix(monkeypatch):
    monkeypatch.setattr(_artifact_io, "_WINDOWS", False)


@pytest.fixture
def windows(monkeypatch):
    """Simulated Windows: the platform hook and an ``os.rename`` that refuses an existing file."""
    monkeypatch.setattr(_artifact_io, "_WINDOWS", True)
    real_rename, calls = os.rename, []

    def rename(source, destination, *args, **kwargs):
        calls.append(pathlib.Path(destination))
        if os.path.lexists(destination):
            raise FileExistsError(errno.EEXIST, "Cannot create a file when that file already exists",
                                  os.fspath(destination))
        real_rename(source, destination, *args, **kwargs)

    monkeypatch.setattr(os, "rename", rename)
    return calls


def test_without_hard_links_publication_fails_typed_and_never_creates_the_blob(cache, monkeypatch, posix):
    cache.blob_dir.mkdir(parents=True)
    links = _no_hard_links(monkeypatch)
    opened = _recorded_opens(monkeypatch)
    error = _error("io", cache.put, GAPPED)
    assert "atomic no-overwrite publication" in error.message and "local disk" in error.message
    assert "replace=True" not in error.message
    destination = cache.blob_path(GAPPED_SHA)
    assert links == [destination]
    assert not os.path.lexists(destination)
    assert destination not in opened          # never created, not even as an empty placeholder
    assert _temporaries(cache) == []


def test_without_hard_links_an_existing_verified_blob_is_still_accepted(cache, monkeypatch, posix):
    cache.put(GAPPED)
    links = _no_hard_links(monkeypatch)
    assert cache.put(GAPPED) == BlobRef(GAPPED_SHA, GAPPED_SIZE)
    assert links == []
    assert _temporaries(cache) == []


def test_without_hard_links_on_windows_publication_is_a_no_replace_rename(cache, monkeypatch, windows):
    _no_hard_links(monkeypatch)
    assert cache.put(GAPPED) == BlobRef(GAPPED_SHA, GAPPED_SIZE)
    assert windows == [cache.blob_path(GAPPED_SHA)]
    assert cache.read(GAPPED_SHA, GAPPED_SIZE) == GAPPED
    assert _temporaries(cache) == []


def test_without_hard_links_on_windows_a_racing_conflicting_file_is_never_replaced(cache, monkeypatch, windows):
    _no_hard_links(monkeypatch)
    cache.blob_dir.mkdir(parents=True)
    real_contains = BlobCache.contains
    calls = []

    def contains(self, sha256, size):
        calls.append(sha256)
        if len(calls) == 1:   # the destination appears after this writer's check
            self.blob_path(sha256).write_bytes(b"conflicting bytes")
            return False
        return real_contains(self, sha256, size)

    monkeypatch.setattr(BlobCache, "contains", contains)
    _error("digest_mismatch", cache.put, GAPPED)
    assert windows == [cache.blob_path(GAPPED_SHA)]
    assert cache.blob_path(GAPPED_SHA).read_bytes() == b"conflicting bytes"
    assert _temporaries(cache) == []


def test_publication_failures_leave_every_existing_blob_untouched(cache, monkeypatch, posix):
    cache.put(b"an earlier valid blob\n")
    before = _snapshot(cache.root)
    _no_hard_links(monkeypatch)
    _error("io", cache.put, GAPPED)
    assert _snapshot(cache.root) == before


# ── cleanup that cannot complete ─────────────────────────────────────


def _failing_unlink(monkeypatch, *, times=None):
    """Make removing this cache's temporary files fail (``times`` times, or always)."""
    real, failures = os.unlink, []

    def unlink(path, *args, **kwargs):
        if pathlib.Path(path).name.startswith(".incoming.") and (times is None or len(failures) < times):
            failures.append(pathlib.Path(path))
            raise PermissionError(errno.EACCES, "The process cannot access the file because it is being used "
                                                "by another process")
        return real(path, *args, **kwargs)

    monkeypatch.setattr(os, "unlink", unlink)
    return failures


def test_publication_survives_a_failed_removal_of_its_temporary_link(cache, monkeypatch):
    # The blob is linked into place, then removing the temporary name fails once, as a
    # Windows sharing violation does while another process has the new blob open.
    cache.blob_dir.mkdir(parents=True)
    failures = _failing_unlink(monkeypatch, times=1)
    assert cache.put(GAPPED) == BlobRef(GAPPED_SHA, GAPPED_SIZE)
    assert len(failures) == 1
    assert sorted(p.name for p in cache.blob_dir.iterdir()) == [GAPPED_SHA.split(":")[1]]
    assert cache.read(GAPPED_SHA, GAPPED_SIZE) == GAPPED


def test_a_temporary_that_cannot_be_removed_never_fails_a_verified_publication(cache, monkeypatch):
    _failing_unlink(monkeypatch)
    assert cache.put(GAPPED) == BlobRef(GAPPED_SHA, GAPPED_SIZE)
    assert cache.read(GAPPED_SHA, GAPPED_SIZE) == GAPPED
    # the undeletable temporary name stays behind under its own name, never as a blob
    assert [p.name for p in cache.blob_dir.iterdir() if not p.name.startswith(".incoming.")] == [
        GAPPED_SHA.split(":")[1]]


def test_discards_are_best_effort_and_never_hide_the_original_error(cache, monkeypatch):
    failures = _failing_unlink(monkeypatch)
    with pytest.raises(RuntimeError, match="original"):
        with cache.writer() as writer:
            writer.write(GAPPED)
            raise RuntimeError("original")
    with cache.writer() as writer:
        writer.write(GAPPED)
        staged = writer.finish()
    staged.discard()
    assert len(failures) == 2
    assert not cache.blob_path(GAPPED_SHA).exists()


def test_a_failed_link_with_no_blob_in_place_is_still_an_io_error(cache, monkeypatch):
    def link(source, destination, *args, **kwargs):
        raise OSError(errno.EIO, "Input/output error")

    monkeypatch.setattr(os, "link", link)
    _error("io", cache.put, GAPPED)
    assert not os.path.lexists(cache.blob_path(GAPPED_SHA))
    assert _temporaries(cache) == []
