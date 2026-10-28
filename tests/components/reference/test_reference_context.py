import json

import pytest

from GenAIRR import data
from GenAIRR import _engine
from GenAIRR._refdata_resolver import dataconfig_to_refdata

CTX = {"schema_version": 1, "catalogue_hash": "sha256:" + "a" * 64, "view": "default",
       "profile": "observed_span_simulation", "models_hash": "sha256:" + "b" * 64}


def _tiny_vj():  # V "TGTAAA" anchor 0, J "TTTGGGAAAGGG" anchor 0, chain "vj"
    r = _engine.RefDataConfig.vj()
    r.add_v_allele("V1", "V", b"TGTAAA", anchor=0)
    r.add_j_allele("J1", "J", b"TTTGGGAAAGGG", anchor=0)
    return r


def test_context_changes_hash_and_seals():
    r = _tiny_vj(); before = r.content_hash()
    r.set_reference_context(json.dumps(CTX))
    assert r.is_sealed() and r.content_hash() != before
    for call in (lambda: r.add_v_allele("V2", "V", b"TGTCCC", anchor=0),
                 lambda: r.set_identity(name="x"), lambda: r.set_allowed_bases(["A"]),
                 lambda: r.curate("raw"), lambda: r.set_reference_context(json.dumps(CTX))):
        with pytest.raises(ValueError, match="sealed"):
            call()


@pytest.mark.parametrize("bad", [{**CTX, "schema_version": 2}, {**CTX, "view": " "},
                                 {**CTX, "catalogue_hash": "sha256:xyz"}, {**CTX, "extra": 1},
                                 {k: v for k, v in CTX.items() if k != "profile"}])
def test_context_rejects_malformed(bad):
    with pytest.raises(ValueError):
        _tiny_vj().set_reference_context(json.dumps(bad))


def test_curated_copy_drops_context():
    r = _tiny_vj(); r.set_reference_context(json.dumps(CTX))
    copy = r.curated("raw")
    assert copy.reference_context_json() is None and not copy.is_sealed()


def test_reference_spec_json_preserves_push_order_and_refuses_sealed():
    r = _tiny_vj(); r.add_v_allele("A0", "A", b"TGTGGG", anchor=0)
    assert [a["name"] for a in json.loads(r.reference_spec_json())["v"]] == ["V1", "A0"]
    r.set_reference_context(json.dumps(CTX))
    with pytest.raises(ValueError, match="reference_spec_v1"):
        r.reference_spec_json()


def test_hash_without_context_is_unchanged_for_bundled_config():
    assert dataconfig_to_refdata(data.HUMAN_IGK_OGRDB).content_hash() == \
        "sha256:6d126831c55da4f7a7f74550123709866a2835c427e16ae1d1b30c1068157b5d"
