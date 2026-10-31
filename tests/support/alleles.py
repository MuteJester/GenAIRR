"""Resolve an allele out of a bound reference.

An oracle may not import GenAIRR, so every test that checks a record against
the reference has to look the allele up itself and hand the bytes over. Three
modules had grown their own copy of this walk, and the fourth would have had
to import one of theirs -- which the no-test-to-test-imports rule forbids, for
the reason this module exists.

Kept deliberately dumb: a linear scan over the pool by name. The pools are a
few hundred entries and the tests that use this run at most a few thousand
lookups, so an index would buy nothing and could disagree with the pool.
"""
from __future__ import annotations


def by_name(refdata, segment: str, name: str):
    """The allele called ``name`` in ``segment``'s pool.

    Raises rather than returning None: a call that names no allele in the
    bound reference means the record cannot be checked against it at all,
    which is a finding and not a case to skip.
    """
    size, getter = _pool(refdata, segment)
    for index in range(size):
        allele = getter(index)
        if allele.name == name:
            return allele
    raise AssertionError(
        f"{name} names no {segment} allele in the bound reference, so the "
        "record cannot be checked against it")


def sequence(allele) -> bytes:
    """The allele's bases, whether ``seq`` is a method or an attribute."""
    raw = allele.seq() if callable(allele.seq) else allele.seq
    return bytes(raw)


def v_subregions(refdata, v_call):
    """The V allele's own IMGT subregion table, as the oracle wants it."""
    return by_name(refdata, "V", v_call).subregions


def _pool(refdata, segment: str):
    pools = {
        "V": (refdata.v_pool_size, refdata.v_allele),
        "D": (refdata.d_pool_size, refdata.d_allele),
        "J": (refdata.j_pool_size, refdata.j_allele),
    }
    key = segment.upper()
    if key not in pools:
        raise ValueError(f"{segment!r} is not a segment; expected V, D or J")
    size, getter = pools[key]
    return size(), getter
