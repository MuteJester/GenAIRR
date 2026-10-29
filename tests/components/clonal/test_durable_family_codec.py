"""Independent value/wire boundaries for durable family artifacts."""
import json
import math
import struct

import pytest

from GenAIRR import _engine

pytestmark = pytest.mark.behavior


@pytest.mark.parametrize("value,tag", [
    (None,"null"), (False,"bool"), (0,"int"), (0.0,"f64"), (-0.0,"f64"),
    (-(2**63),"int"), (2**64-1,"int"), (float.fromhex("0x0.0000000000001p-1022"),"f64"),
    ("8000000000000000","str"), ({"t":"f64","v":"8000000000000000"},"object"),
    ([False,0,0.0,-0.0,"é\n\"\u0000"],"array"),
])
def test_lossless_dynamic_values_match_literal_wire_types(value,tag):
    text=_engine._family_replay_pack(value)
    wire=json.loads(text)
    assert wire["t"]==tag
    recovered=_engine._family_replay_unpack(text)
    assert _engine._family_replay_pack(recovered)==text
    assert type(recovered) is type(value)
    if isinstance(value,float):
        assert struct.pack(">d",recovered)==struct.pack(">d",value)
    if value==0.0 and type(value) is float:
        assert wire["v"] == ("8000000000000000" if math.copysign(1,value)<0 else "0000000000000000")


@pytest.mark.parametrize("value", [float("nan"),float("inf"),-float("inf"),2**64,-2**63-1,b"A",{1:"key"},(1,2),{1,2},object()])
def test_non_json_or_unrepresentable_metadata_rejects(value):
    with pytest.raises(ValueError):
        _engine._family_replay_pack(value)


def test_cycles_reject_but_repeated_noncyclic_values_remain_valid():
    cyclic=[]
    cyclic.append(cyclic)
    with pytest.raises(ValueError,match="cyclic"):
        _engine._family_replay_pack(cyclic)
    same={"ok":True}
    assert _engine._family_replay_unpack(_engine._family_replay_pack([same,same]))==[same,same]


@pytest.mark.parametrize("text", [
    '{"x":1,"x":2}', '{"nested":{"x":1,"\\u0078":2}}',
    '[NaN]', '[Infinity]', '1 2', '{"x":18446744073709551616}',
])
def test_strict_parser_rejects_duplicate_keys_nonfinite_and_oversized_integers(text):
    with pytest.raises(ValueError):
        _engine._family_replay_parse_json(text)


@pytest.mark.parametrize("text", [
    '{"t":"bool","v":0}', '{"t":"int","v":"01"}',
    '{"t":"int","v":"-0"}', '{"t":"int","v":"18446744073709551616"}',
    '{"t":"f64","v":"7ff0000000000000"}', '{"t":"f64","v":"000000000000000A"}',
    '{"t":"null","v":null}', '{"t":"str","v":1}',
    '{"t":"object","v":[["x",{"t":"null"}],["x",{"t":"null"}]]}',
])
def test_dynamic_wire_tags_are_closed_and_typed(text):
    with pytest.raises(ValueError):
        _engine._family_replay_unpack(text)


def test_native_canonical_hash_distinguishes_types_domain_and_delimiter_aliases():
    def digest(value,domain="test"):
        return _engine._family_replay_digest(json.dumps(value),domain)
    assert len({digest(v) for v in [False,0,0.0,-0.0,"0000000000000000",None]})==6
    assert digest(["a|b","c"])!=digest(["a","b|c"])
    assert digest({"species":"a\nidentity.locus=b","locus":"c"}) != digest({"species":"a","locus":"b\nidentity.locus=c"})
    assert digest({"b":2,"a":1})==digest({"a":1,"b":2})
    assert digest(1,"first")!=digest(1,"second")


@pytest.mark.parametrize("text,limit", [
    ('[[[0]]]',{"max_depth":2}), ('[1,2,3]',{"max_values":3}),
    ('"abcd"',{"max_string_bytes":3}), ('[1,2]',{"max_file_bytes":4}),
])
def test_parser_enforces_caller_owned_limits(text,limit):
    with pytest.raises(ValueError,match="resource limit"):
        _engine._family_replay_parse_json(text,json.dumps(limit))


def test_shared_budget_is_checked_and_cannot_wrap_counters():
    budget=_engine.FamilyReplayBudget(json.dumps({"max_work":2}))
    budget.charge_work(2)
    with pytest.raises(ValueError,match="resource limit"):
        budget.charge_work(1)
    assert json.loads(budget.snapshot())["work"]==2
    with pytest.raises(ValueError):
        _engine.FamilyReplayBudget('{"unknown_limit":1}')
