"""A duplex repair fixture with a positionally independent reference.

Builds the refdata, repair script and seed together, so a test can assert on a duplex
observation without restating the whole arrangement.
"""
from __future__ import annotations

import hashlib
import json
from GenAIRR import _engine as ge
from GenAIRR.experimental import CoupledDuplexRepair, DuplexJoin, project_duplex_simulation

def _opened(name):
    return DuplexJoin(name).open("left", offset=0).open("right", offset=0)

def _fixture(*, order=("vd", "dj"), inverted=False, d="CCC", v="AAATGT", verified=False):
    refs = ge.RefDataConfig.vdj()
    refs.add_v_allele("V*01", "V", v.encode(), anchor=3)
    refs.add_d_allele("D*01", "D", d.encode())
    refs.add_j_allele("J*01", "J", b"TTTAAA", anchor=0)
    if verified:
        refs.set_allele_evidence("V", 0, json.dumps({
            "schema_version": 1,
            "sequence_sha256": "sha256:" + hashlib.sha256(v.encode()).hexdigest(),
            "coding_phase": None,
            "anchor_status": "verified",
            "boundary_5": "complete",
            "boundary_3": "complete",
            "source": "independent positional fixture",
            "method": "authored_coordinate",
            "notes": [],
        }))
    scripts = {
        "vd": (_opened("vd").resect("left", "top", count=1)
               .synapse(right_start=6).fill("top").ligate()),
        "dj": _opened("dj").synapse(right_start=3).ligate(),
    }
    repair = CoupledDuplexRepair(v, d, "TTTAAA", d_inverted=inverted)
    for name in order:
        repair = repair.join(scripts[name])
    seed = repair.seed(
        refs, v_call="V*01", d_call="D*01", j_call="J*01",
        reference_ends_are_coding_ends=True,
    )
    return refs, repair, seed

def _compile(plan, refs, *, productive=False, strict=True):
    return plan.compile(
        refdata=refs, respect=ge.productive() if productive else None,
        strict=strict, allow_curatable_refdata=True,
    )
