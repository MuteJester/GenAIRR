import contextlib
import errno
import json
import os
import pathlib
import stat

import pytest

from GenAIRR import _artifact_io
from GenAIRR._artifact_io import (
    ArtifactIOError, canonical_json, strict_loads, typed_digest, write_atomic,
)
from GenAIRR._family_replay import FamilyReplayError, FamilyReplayBundle


def test_strict_loads_rejects_duplicate_keys():
    with pytest.raises(ArtifactIOError) as info:
        strict_loads('{"a":1,"a":2}', max_bytes=1000)
    assert info.value.kind == "format"


def test_strict_loads_rejects_oversize():
    with pytest.raises(ArtifactIOError) as info:
        strict_loads('{"a":"' + "x" * 100 + '"}', max_bytes=10)
    assert info.value.kind == "resource"


def test_strict_loads_rejects_trailing_and_nonfinite():
    for text in ('{"a":1} x', '{"a":NaN}'):
        with pytest.raises(ArtifactIOError) as info:
            strict_loads(text, max_bytes=1000)
        assert info.value.kind == "format"


def test_strict_loads_strips_one_bom():
    assert strict_loads("\ufeff" + '{"a":1}', max_bytes=1000) == {"a": 1}
    with pytest.raises(ArtifactIOError):
        strict_loads("\ufeff\ufeff" + '{"a":1}', max_bytes=1000)


def test_typed_digest_ignores_key_order():
    a = typed_digest({"x": 1, "y": [1, 2]}, "test.domain")
    b = typed_digest({"y": [1, 2], "x": 1}, "test.domain")
    assert a == b
    assert a != typed_digest({"x": 1, "y": [1, 2]}, "other.domain")
    assert canonical_json({"b": 1, "a": 2}) == '{"a":2,"b":1}'


def test_write_atomic_no_clobber(tmp_path):
    target = tmp_path / "out.json"
    write_atomic(target, b"first", replace=False)
    with pytest.raises(ArtifactIOError) as info:
        write_atomic(target, b"second", replace=False)
    assert info.value.kind == "io"
    assert target.read_bytes() == b"first"
    assert list(tmp_path.glob("*.tmp")) == [] and list(tmp_path.glob(".*.tmp")) == []


def test_write_atomic_replace(tmp_path):
    target = tmp_path / "out.json"
    write_atomic(target, b"first", replace=False)
    write_atomic(target, b"second", replace=True)
    assert target.read_bytes() == b"second"
    assert list(tmp_path.glob(".*.tmp")) == []


def test_write_atomic_requires_bool_replace(tmp_path):
    with pytest.raises(TypeError):
        write_atomic(tmp_path / "x", b"a", replace=1)


def test_family_replay_write_errors_stay_family_replay_errors(tmp_path):
    bundle = FamilyReplayBundle.__new__(FamilyReplayBundle)
    bundle._text = "{}"
    target = tmp_path / "b.json"
    target.write_text("existing", encoding="utf-8")
    with pytest.raises(FamilyReplayError) as info:
        bundle.write_to(target)
    assert "family_replay[io]" in str(info.value)
    assert target.read_text(encoding="utf-8") == "existing"


# ── File permissions follow the umask, as an ordinary file write would ─────

posix_only = pytest.mark.skipif(os.name == "nt", reason="POSIX permission bits")


@contextlib.contextmanager
def umask(mask):
    previous = os.umask(mask)
    try:
        yield
    finally:
        os.umask(previous)


@posix_only
@pytest.mark.parametrize("replace", [False, True])
@pytest.mark.parametrize("mask,mode", [(0o022, 0o644), (0o077, 0o600)])
def test_write_atomic_mode_follows_umask(tmp_path, mask, mode, replace):
    target = tmp_path / "out.json"
    if replace:
        target.write_bytes(b"old")
        target.chmod(0o600 if mode == 0o644 else 0o644)
    with umask(mask):
        write_atomic(target, b"data", replace=replace)
    assert stat.S_IMODE(target.stat().st_mode) == mode
    assert target.read_bytes() == b"data"
    assert list(tmp_path.glob(".*.tmp")) == []


@posix_only
@pytest.mark.parametrize("mask,mode", [(0o022, 0o644), (0o077, 0o600)])
def test_family_replay_write_to_mode_follows_umask(tmp_path, mask, mode):
    bundle = FamilyReplayBundle.__new__(FamilyReplayBundle)
    bundle._text = "{}"
    target = tmp_path / "b.json"
    with umask(mask):
        bundle.write_to(target)
    assert stat.S_IMODE(target.stat().st_mode) == mode


def test_write_atomic_skips_an_existing_temporary_name(tmp_path, monkeypatch):
    suffixes = iter(["taken", "taken", "fresh"])
    monkeypatch.setattr(_artifact_io, "token_hex", lambda nbytes: next(suffixes))
    other = tmp_path / ".out.json.taken.tmp"      # another writer's temporary file
    other.write_bytes(b"not ours")
    target = tmp_path / "out.json"
    write_atomic(target, b"data", replace=False)
    assert target.read_bytes() == b"data"
    assert other.read_bytes() == b"not ours"       # never opened, truncated or removed
    assert sorted(p.name for p in tmp_path.iterdir()) == [".out.json.taken.tmp", "out.json"]



def test_package_sources_have_no_literal_byte_order_mark():
    # An invisible U+FEFF in source is easy to lose or duplicate; it is spelled as an escape.
    package = pathlib.Path(_artifact_io.__file__).parent
    sources = sorted(package.rglob("*.py"))
    assert pathlib.Path(_artifact_io.__file__) in sources
    assert [str(p.relative_to(package)) for p in sources if b"\xef\xbb\xbf" in p.read_bytes()] == []


# ── Filesystems without hard links (FAT/exFAT, some SMB shares) ────────
#
# Without hard links POSIX has no portable atomic no-overwrite publish, so ``replace=False``
# fails with a typed error and never creates the destination.  Windows ``os.rename`` refuses
# to overwrite an existing file, so it publishes there.


NO_LINK_ERRORS = [PermissionError(errno.EPERM, "Operation not permitted"),
                  OSError(errno.ENOTSUP, "Operation not supported"),
                  OSError(errno.EOPNOTSUPP, "Operation not supported"),
                  OSError(errno.ENOSYS, "Function not implemented"),
                  OSError(errno.EXDEV, "Invalid cross-device link")]
NO_LINK_IDS = ["EPERM", "ENOTSUP", "EOPNOTSUPP", "ENOSYS", "EXDEV"]
NO_ATOMIC_NO_OVERWRITE = "does not support atomic no-overwrite publication"


def _no_hard_links(monkeypatch, error):
    calls = []

    def link(source, destination, *args, **kwargs):
        calls.append(destination)
        raise error

    monkeypatch.setattr(os, "link", link)
    return calls


def _temporaries(directory):
    return [p.name for p in directory.iterdir() if p.name.endswith(".tmp")]


def _recorded(monkeypatch, name):
    """Wrap ``os.<name>`` so every call's positional arguments are recorded."""
    real, calls = getattr(os, name), []

    def wrapper(*args, **kwargs):
        calls.append(args)
        return real(*args, **kwargs)

    monkeypatch.setattr(os, name, wrapper)
    return calls


def _names(args, path):
    return any(isinstance(a, (str, os.PathLike)) and os.fspath(a) == os.fspath(path) for a in args)


@pytest.fixture
def posix(monkeypatch):
    monkeypatch.setattr(_artifact_io, "_WINDOWS", False, raising=False)


@pytest.fixture
def windows(monkeypatch):
    """Simulated Windows: the platform hook, and an ``os.rename`` that refuses an existing file."""
    monkeypatch.setattr(_artifact_io, "_WINDOWS", True, raising=False)
    real_rename, calls = os.rename, []

    def rename(source, destination, *args, **kwargs):
        calls.append((pathlib.Path(source), pathlib.Path(destination)))
        if os.path.lexists(destination):
            raise FileExistsError(errno.EEXIST, "Cannot create a file when that file already exists",
                                  os.fspath(destination))
        real_rename(source, destination, *args, **kwargs)

    monkeypatch.setattr(os, "rename", rename)
    return calls


@pytest.mark.parametrize("error", NO_LINK_ERRORS, ids=NO_LINK_IDS)
def test_write_atomic_without_hard_links_fails_typed_and_never_creates_the_destination(
        tmp_path, monkeypatch, posix, error):
    links = _no_hard_links(monkeypatch, error)
    opened = _recorded(monkeypatch, "open")
    replaced = _recorded(monkeypatch, "replace")
    target = tmp_path / "out.json"
    with pytest.raises(ArtifactIOError) as info:
        write_atomic(target, b"data", replace=False)
    assert info.value.kind == "io"
    message = str(info.value)
    assert NO_ATOMIC_NO_OVERWRITE in message
    assert "local disk" in message and "replace=True" in message
    assert links == [target]
    assert not target.exists()
    assert not any(_names(call, target) for call in opened)    # never created, not even empty
    assert replaced == []
    assert _temporaries(tmp_path) == []


def test_write_atomic_without_hard_links_never_removes_another_writers_artifact(tmp_path, monkeypatch, posix):
    # Another writer completes the destination while this call is running.  This call
    # fails, and the other writer's complete artifact must survive it.
    target = tmp_path / "out.json"
    real_replace, real_unlink = os.replace, os.unlink
    theirs = b"another writer's complete artifact"
    armed, removed = [], []

    def other_writer_finishes():
        if armed:
            armed.clear()
            staged = tmp_path / "theirs.staged"
            staged.write_bytes(theirs)
            real_replace(staged, target)

    def link(source, destination, *args, **kwargs):
        armed.append(True)                  # from here the other writer may finish at any step
        raise PermissionError(errno.EPERM, "Operation not permitted")

    def replace(source, destination, *args, **kwargs):
        other_writer_finishes()
        if pathlib.Path(source).name.startswith(".out.json."):
            raise OSError(errno.EIO, "Input/output error")      # publishing our file fails
        return real_replace(source, destination, *args, **kwargs)

    def unlink(path, *args, **kwargs):
        other_writer_finishes()
        removed.append(pathlib.Path(path))
        return real_unlink(path, *args, **kwargs)

    monkeypatch.setattr(os, "link", link)
    monkeypatch.setattr(os, "replace", replace)
    monkeypatch.setattr(os, "unlink", unlink)
    with pytest.raises(ArtifactIOError) as info:
        write_atomic(target, b"ours", replace=False)
    assert info.value.kind == "io"
    assert target.read_bytes() == theirs
    assert target not in removed
    assert _temporaries(tmp_path) == []


def test_write_atomic_without_hard_links_on_windows_publishes_by_rename(tmp_path, monkeypatch, windows):
    _no_hard_links(monkeypatch, NO_LINK_ERRORS[0])
    replaced = _recorded(monkeypatch, "replace")
    target = tmp_path / "out.json"
    write_atomic(target, b"data", replace=False)
    assert target.read_bytes() == b"data"
    ((source, destination),) = windows
    assert destination == target and source.name.startswith(".out.json.") and source.name.endswith(".tmp")
    assert replaced == []
    assert _temporaries(tmp_path) == []


def test_write_atomic_without_hard_links_on_windows_refuses_an_existing_file(tmp_path, monkeypatch, windows):
    _no_hard_links(monkeypatch, NO_LINK_ERRORS[0])
    target = tmp_path / "out.json"
    target.write_bytes(b"existing")
    with pytest.raises(ArtifactIOError) as info:
        write_atomic(target, b"second", replace=False)
    assert info.value.kind == "io"
    assert [destination for _, destination in windows] == [target]
    assert target.read_bytes() == b"existing"
    assert _temporaries(tmp_path) == []


@pytest.mark.parametrize("platform", ["posix", "windows"])
def test_write_atomic_replace_works_without_hard_links(tmp_path, monkeypatch, request, platform):
    request.getfixturevalue(platform)
    links = _no_hard_links(monkeypatch, NO_LINK_ERRORS[0])
    target = tmp_path / "out.json"
    write_atomic(target, b"first", replace=True)
    assert target.read_bytes() == b"first"
    write_atomic(target, b"second", replace=True)
    assert target.read_bytes() == b"second"
    assert links == []
    assert _temporaries(tmp_path) == []


def test_write_atomic_other_link_errors_are_io_errors(tmp_path, monkeypatch):
    _no_hard_links(monkeypatch, OSError(errno.EIO, "Input/output error"))
    target = tmp_path / "out.json"
    with pytest.raises(ArtifactIOError) as info:
        write_atomic(target, b"data", replace=False)
    assert info.value.kind == "io"
    assert not target.exists()
    assert _temporaries(tmp_path) == []


def test_family_replay_write_without_hard_links_is_a_family_replay_io_error(tmp_path, monkeypatch, posix):
    _no_hard_links(monkeypatch, NO_LINK_ERRORS[0])
    bundle = FamilyReplayBundle.__new__(FamilyReplayBundle)
    bundle._text = "{}"
    target = tmp_path / "b.json"
    with pytest.raises(FamilyReplayError) as info:
        bundle.write_to(target)
    assert info.value.category == "io"
    assert NO_ATOMIC_NO_OVERWRITE in str(info.value)
    assert not target.exists()
    assert _temporaries(tmp_path) == []
    bundle.write_to(target, replace=True)
    assert target.read_text(encoding="utf-8") == "{}"
