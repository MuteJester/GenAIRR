"""Refdata and experiments the DSL tests build on.

Small hand-built configurations whose anchors and coordinates are stated here
rather than looked up, so a test asserting on DSL behavior can be read on its
own.
"""
from __future__ import annotations
from GenAIRR import _engine as ge
from GenAIRR import Experiment
import re as _re

def _vj_refdata(*, canonical=False) -> ge.RefDataConfig:
    cfg = ge.RefDataConfig.vj()
    cfg.add_v_allele("v1*01", "v1", b"AAACCCTGT" if canonical else b"AAACCCGGG", anchor=6)
    cfg.add_j_allele("j1*01", "j1", b"TTTAAA", anchor=0)
    return cfg
def _vdj_refdata(*, canonical=False) -> ge.RefDataConfig:
    cfg = ge.RefDataConfig.vdj()
    cfg.add_v_allele("v1*01", "v1", b"AAACCCTGT" if canonical else b"AAACCCGGG", anchor=6)
    cfg.add_d_allele("d1*01", "d1", b"TTTTTT")
    cfg.add_j_allele("j1*01", "j1", b"TTTAAA", anchor=0)
    return cfg
def _human_igh_exp():
    """Convenience: default human-IGH recombination experiment."""
    return Experiment.on("human_igh").recombine()
def _human_igh_records(n: int = 5, seed: int = 42, **kw):
    """Convenience: run a small human-IGH batch through ``run_records``
    and return the :class:`SimulationResult`."""
    return _human_igh_exp().run_records(n=n, seed=seed, **kw)
def _alignment_invariants_hold(rec) -> list:
    """Return a list of issue strings for a record's alignment triplet
    — empty list if every invariant passes."""
    sa = rec["sequence_alignment"]
    ga_str = rec["germline_alignment"]
    dm = rec["germline_alignment_d_mask"]
    seq = rec["sequence"]
    issues = []
    if not (len(sa) == len(ga_str) == len(dm)):
        issues.append(
            f"length mismatch: sa={len(sa)} ga={len(ga_str)} dm={len(dm)}"
        )
    if len(sa) != len(seq) + sa.count("-"):
        issues.append(
            f"sa len {len(sa)} != seq len {len(seq)} + sa gaps {sa.count('-')}"
        )
    if sa.replace("-", "") != seq.upper():
        issues.append("sa minus gaps != upper(seq)")
    return issues
def _cigar_total_ops(cigar: str) -> int:
    """Sum the run-length counts in a CIGAR string. Used to check
    that the CIGAR's coverage equals the segment's alignment-column
    count."""
    return sum(int(n) for n in _re.findall(r"(\d+)[MID=X]", cigar))
def _cigar_op_counts(cigar: str) -> dict:
    """Return ``{op: total_count}`` so tests can assert e.g. how many
    insertions vs deletions appear."""
    counts: dict = {"M": 0, "I": 0, "D": 0}
    for n, op in _re.findall(r"(\d+)([MID=X])", cigar):
        counts[op] = counts.get(op, 0) + int(n)
    return counts
