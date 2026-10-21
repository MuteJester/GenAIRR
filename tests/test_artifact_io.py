import json
import pytest

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
