"""Backward-compat: the genotype machinery is purely additive, so a
run with NO genotype attached must be byte-identical to master.

``MASTER_DIGEST`` was captured on ``master`` (pre-genotype) for
``Experiment.on(HUMAN_IGH_OGRDB).recombine().run_records(n=100, seed=12345)``
and re-verified identical on this branch's no-genotype path.
"""
import hashlib

import GenAIRR as ga
import GenAIRR.data as gdata

MASTER_DIGEST = "3be8e5ea124e1dfff256f93b5ddbb925fdb738d4d6f2eeb5763a81e5b6213460"


def _digest(records):
    h = hashlib.sha256()
    for rec in records:
        h.update(repr(sorted(rec.items())).encode())
    return h.hexdigest()


def test_no_genotype_output_matches_master_baseline():
    res = ga.Experiment.on(gdata.HUMAN_IGH_OGRDB).recombine().run_records(
        n=100, seed=12345
    )
    assert _digest(res) == MASTER_DIGEST


def test_no_genotype_output_is_deterministic():
    a = ga.Experiment.on(gdata.HUMAN_IGH_OGRDB).recombine().run_records(n=100, seed=12345)
    b = ga.Experiment.on(gdata.HUMAN_IGH_OGRDB).recombine().run_records(n=100, seed=12345)
    assert _digest(a) == _digest(b)
