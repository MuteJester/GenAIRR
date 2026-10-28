"""Minimal hand-built refdata for junction and call provenance.

Small enough that every anchor and coordinate in it is stated here rather
than looked up, so a test asserting on a junction can be read on its own.
"""
from __future__ import annotations

from GenAIRR import _engine as ge

def _vj_basic() -> "ge.RefDataConfig":
    cfg = ge.RefDataConfig.vj()
    # V anchor codon = V[6:9] = TGT (Cys)
    cfg.add_v_allele("v1*01", "v1", b"AAACCCTGT", anchor=6)
    # J anchor codon = J[0:3] = TGG (Trp)
    cfg.add_j_allele("j1*01", "j1", b"TGGAAA", anchor=0)
    return cfg
