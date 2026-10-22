import contextlib
import json
import os
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
    assert strict_loads("﻿" + '{"a":1}', max_bytes=1000) == {"a": 1}
    with pytest.raises(ArtifactIOError):
        strict_loads("﻿﻿" + '{"a":1}', max_bytes=1000)


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
    target.write_text("existing")
    with pytest.raises(FamilyReplayError) as info:
        bundle.write_to(target)
    assert "family_replay[io]" in str(info.value)
    assert target.read_text() == "existing"


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

