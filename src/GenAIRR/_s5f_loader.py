"""Builtin S5F kernel loaders.

Loads the empirical S5F mutability + substitution tables from the
bundled ``data/mutation_model_parameters/*.pkl`` pickles and
reshapes them into the flat ``(mutability, substitution)`` lists
that ``genairr_engine.PassPlan.push_mutate_s5f`` expects.

The bundled pickles are 3-tuples of ``(mutability, substitution,
combined)`` dicts keyed by 5-mer strings. ``mutability[5mer]`` is a
float; ``substitution[5mer]`` is a per-destination dict
``{"A": p, "C": p, "G": p, "T": p}`` (the source base entry is
absent — S5F is "given a mutation, what does the central base
become"). Keys can include ``"N"`` characters; these get filtered
out since the engine's kernel only handles canonical A/C/G/T
contexts.

The Rust kernel encodes a context as
``b1<<8 | b2<<6 | b3<<4 | b4<<2 | b5`` with A=0, C=1, G=2, T=3.
The substitution Vec is indexed at ``context*4 + dest`` with the
same base ordering.
"""
from __future__ import annotations

import pickle
from functools import lru_cache
from importlib.resources import files
from typing import List, Tuple

# A/C/G/T ↔ 0/1/2/3, matching crate::ir::encode_base in Rust.
_BASE_INDEX = {"A": 0, "C": 1, "G": 2, "T": 3}
_S5F_NUM_CONTEXTS = 1024     # 4⁵
_S5F_SUBSTITUTION_LEN = 4096  # 1024 × 4


def _encode_context(fivemer: str) -> int:
    """Map a 5-mer of canonical A/C/G/T to a kernel context index.

    Returns -1 for any 5-mer containing a non-canonical base (e.g. N).
    """
    idx = 0
    for ch in fivemer:
        b = _BASE_INDEX.get(ch.upper())
        if b is None:
            return -1
        idx = (idx << 2) | b
    return idx


def _build_kernel_lists(
    mutability_dict, substitution_dict
) -> Tuple[List[float], List[float]]:
    """Translate the dict-form S5F tables into flat lists for Rust.

    Missing entries default to zero — both for unmutable contexts
    (mutability) and for the source-base destination (which S5F's
    substitution table doesn't include since the source can't equal
    the destination after a mutation).
    """
    mutability = [0.0] * _S5F_NUM_CONTEXTS
    substitution = [0.0] * _S5F_SUBSTITUTION_LEN

    for fivemer, mu in mutability_dict.items():
        ctx = _encode_context(fivemer)
        if ctx < 0:
            continue
        mutability[ctx] = float(mu)

    for fivemer, dest_probs in substitution_dict.items():
        ctx = _encode_context(fivemer)
        if ctx < 0:
            continue
        for dest_base, prob in dest_probs.items():
            dest_idx = _BASE_INDEX.get(dest_base.upper())
            if dest_idx is None:
                continue
            substitution[ctx * 4 + dest_idx] = float(prob)

    return mutability, substitution


# Map of canonical lowercase model name → bundled pickle filename.
_BUILTIN_S5F_MODELS = {
    "hh_s5f": "HH_S5F_META.pkl",
    "hh_s5f_60": "HH_S5F_60_META.pkl",
    "hh_s5f_opposite": "HH_S5F_Opposite_META.pkl",
    "hkl_s5f": "HKL_S5F_META.pkl",
}


@lru_cache(maxsize=None)
def load_builtin_s5f_kernel(name: str) -> Tuple[List[float], List[float]]:
    """Load a bundled S5F kernel by short name.

    Available names: ``"hh_s5f"`` (default human heavy chain),
    ``"hh_s5f_60"`` (the 60% subset), ``"hh_s5f_opposite"``,
    ``"hkl_s5f"`` (human kappa/lambda). Names are case-insensitive.

    Returns ``(mutability, substitution)`` flat lists of length 1024
    and 4096 respectively, ready for
    ``genairr_engine.PassPlan.push_mutate_s5f``.
    """
    key = name.lower().strip()
    if key not in _BUILTIN_S5F_MODELS:
        avail = ", ".join(sorted(_BUILTIN_S5F_MODELS))
        raise ValueError(
            f"Unknown S5F model {name!r}. Available: {avail}"
        )
    pkg = "GenAIRR.data.mutation_model_parameters"
    fname = _BUILTIN_S5F_MODELS[key]
    raw = files(pkg).joinpath(fname).read_bytes()
    obj = pickle.loads(raw)
    if not isinstance(obj, tuple) or len(obj) < 2:
        raise ValueError(
            f"S5F pickle {fname!r} has unexpected shape {type(obj).__name__}"
        )
    mutability_dict, substitution_dict = obj[0], obj[1]
    return _build_kernel_lists(mutability_dict, substitution_dict)
