"""Strategies for experiment configurations the suite can actually build."""
from __future__ import annotations

from hypothesis import strategies as st

#: Loci with packaged reference data, so a drawn locus always builds.
LOCI = ("human_igh", "human_igk", "human_tcrb")


def loci() -> st.SearchStrategy[str]:
    """One of the loci the suite can build without extra fixtures."""
    return st.sampled_from(LOCI)


def seeds() -> st.SearchStrategy[int]:
    """A seed in the range the engine accepts."""
    return st.integers(min_value=0, max_value=2**31 - 1)


def recombinations(max_sequences: int = 3) -> st.SearchStrategy[tuple]:
    """``(locus, seed, n)`` for a recombination run small enough to stay quick."""
    return st.tuples(loci(), seeds(), st.integers(min_value=1, max_value=max_sequences))
