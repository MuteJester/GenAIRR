import hashlib

import pytest

from GenAIRR import data
from GenAIRR._dataconfig_extract import extract_recombine_defaults
from GenAIRR._reference import defaults as d

pytestmark = pytest.mark.behavior

VIEW_VDJ = {"V": ["V:IGHV1-2*01", "V:IGHV1-3*01"], "D": ["D:IGHD1-1*01"], "J": ["J:IGHJ4*01", "J:IGHJ5*01"]}
VIEW_VJ = {"V": ["V:IGKV1-5*01"], "J": ["J:IGKJ1*01"]}


def test_vdj_profile_shape():
    rec = d.unfitted_defaults_v1("vdj", VIEW_VDJ)
    assert rec["profile"] == d.PROFILE == "genairr.unfitted_defaults.v1"
    assert rec["allele_usage"] == {
        "unit": "allele", "kind": "uniform",
        "V": ["IGHV1-2*01", "IGHV1-3*01"], "D": ["IGHD1-1*01"], "J": ["IGHJ4*01", "IGHJ5*01"],
    }
    assert rec["np_lengths"] == {k: {"kind": "uniform", "support": [0, 1, 2, 3, 4, 5, 6]} for k in ("NP1", "NP2")}
    assert rec["np_bases"] == {k: {"kind": "uniform", "support": ["A", "C", "G", "T"]} for k in ("NP1", "NP2")}
    assert rec["trims"] == {k: {"kind": "disabled"} for k in ("V_3", "D_5", "D_3", "J_5")}
    assert rec["p_nucleotides"] == rec["trims"]


def test_vj_profile_omits_d_planes():
    rec = d.unfitted_defaults_v1("vj", VIEW_VJ)
    assert "D" not in rec["allele_usage"]
    assert set(rec["np_lengths"]) == {"NP1"}
    assert set(rec["np_bases"]) == {"NP1"}
    assert set(rec["trims"]) == {"V_3", "J_5"}
    assert set(rec["p_nucleotides"]) == {"V_3", "J_5"}


def test_recombine_defaults_values():
    rec = d.unfitted_defaults_v1("vdj", VIEW_VDJ)
    out = d.recombine_defaults(rec, chain_type="vdj")
    assert set(out) == set(extract_recombine_defaults(data.HUMAN_IGH_OGRDB))
    assert out["np1"] == out["np2"] == [(v, 1.0) for v in range(7)]
    assert out["allele_usage_v"] == {"IGHV1-2*01": 1.0, "IGHV1-3*01": 1.0}
    assert out["allele_usage_d"] == {"IGHD1-1*01": 1.0}
    for key, value in out.items():
        if not key.startswith(("np1", "np2", "allele_usage")) or key.endswith(("bases", "transitions")):
            assert value is None, key
    vj = d.recombine_defaults(d.unfitted_defaults_v1("vj", VIEW_VJ), chain_type="vj")
    assert vj["np2"] is None and vj["allele_usage_d"] is None
    assert vj["np1"] == [(v, 1.0) for v in range(7)]


def test_models_hash_stable_and_sensitive():
    a = d.unfitted_defaults_v1("vdj", VIEW_VDJ)
    b = d.unfitted_defaults_v1("vdj", VIEW_VDJ)
    assert d.models_hash(a) == d.models_hash(b)
    assert d.models_hash(a).startswith("sha256:")
    c = d.unfitted_defaults_v1("vdj", {**VIEW_VDJ, "J": ["J:IGHJ4*01"]})
    assert d.models_hash(a) != d.models_hash(c)


def test_extract_recombine_defaults_unchanged_for_bundled():
    out = extract_recombine_defaults(data.HUMAN_IGH_OGRDB)
    assert list(out) == [
        "np1", "np2", "np1_bases", "np2_bases", "np1_markov_transitions", "np2_markov_transitions",
        "p_v_3_lengths", "p_d_5_lengths", "p_d_3_lengths", "p_j_5_lengths",
        "trim_v_3", "trim_d_5", "trim_d_3", "trim_j_5",
        "allele_usage_v", "allele_usage_d", "allele_usage_j",
    ]
    assert out["np1"][:3] == [(0, 0.05265143785167222), (1, 0.04075110575704825), (2, 0.05608972515275564)]
    assert out["trim_v_3"][:3] == [(0, 0.434587338005721), (1, 0.27755841545908777), (2, 0.11297970871321363)]
    assert hashlib.sha256(repr(out).encode()).hexdigest() == (
        "f294d973df0a9fbe622fd5b287582501446fd079203f3662202ebec000c2cd44"
    )
