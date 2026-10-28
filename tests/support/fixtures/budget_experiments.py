"""VDJ experiments configured for the wall-time budget checks.

Uniform NP-length distributions keep the generated work predictable, so a
budget measures the engine rather than the draw.
"""
from __future__ import annotations

import GenAIRR as ga
from GenAIRR import _engine as ge

def _vdj(*, canonical=False) -> "ge.RefDataConfig":
    # Productive benchmarks need usable original landmarks. Keep the exact
    # lengths, coordinates and flanks; replace only the two declared codons.
    v = b"GAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTACAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTC"
    j = b"TACTACTACGGTATGGACGTCTGGGGCCAAGGGACCACGGTCACCGTCTCCTCAG"
    if canonical:
        v = v[:78] + b"TGT" + v[81:]
        j = j[:8] + b"TGG" + j[11:]
    cfg = ge.RefDataConfig.vdj()
    cfg.add_v_allele(
        "v1*01", "v1",
        v,
        anchor=78,
    )
    cfg.add_d_allele("d1*01", "d1", b"GGGTATAGCAGCAGCTGGTAC")
    cfg.add_j_allele(
        "j1*01", "j1",
        j,
        anchor=8,
    )
    return cfg

def _vdj_recombine(*, canonical=False) -> "ga.Experiment":
    return (
        ga.Experiment.on(_vdj(canonical=canonical)).allow_curatable_refdata()
        .recombine(
            np1_lengths=[(L, 1.0) for L in range(20)],
            np2_lengths=[(L, 1.0) for L in range(20)],
        )
        .trim(
            v_3=[(L, 1.0) for L in range(5)],
            d_5=[(L, 1.0) for L in range(5)],
            d_3=[(L, 1.0) for L in range(5)],
            j_5=[(L, 1.0) for L in range(5)],
        )
    )
