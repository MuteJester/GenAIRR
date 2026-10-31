"""Experiments in their commonest shapes.

``Experiment.on("human_igh")`` appears in the suite a few hundred times, usually
followed by ``.recombine()``. Building it here means a change to the default
arrangement is one edit rather than hundreds.
"""
from __future__ import annotations

from typing import Any

import GenAIRR as ga

DEFAULT_LOCUS = "human_igh"


#: Loci whose packaged reference carries alleles with non-canonical anchors
#: (pseudogene / ORF), which the compiler refuses unless curation is allowed.
NEEDS_CURATABLE = frozenset({"human_tcrb"})


def recombined(locus: str = DEFAULT_LOCUS, *, allow_curatable: bool | None = None,
               **recombine: Any):
    """An experiment bound to ``locus`` with recombination configured.

    Keywords go straight to ``recombine``, so a caller needing particular
    allele weights or trim distributions passes them here.

    ``allow_curatable`` defaults to whatever the locus requires: some packaged
    references carry pseudogene or ORF alleles whose anchor codon does not
    translate as expected, and the compiler refuses those unless curation is
    allowed. Defaulting per locus keeps a cross-locus test from having to know
    which ones.
    """
    experiment = ga.Experiment.on(locus)
    if allow_curatable or (allow_curatable is None and locus in NEEDS_CURATABLE):
        experiment = experiment.allow_curatable_refdata()
    return experiment.recombine(**recombine)


def with_p_nucleotides(locus: str = DEFAULT_LOCUS, **recombine: Any):
    """A recombined experiment that reliably produces P nucleotides.

    P addition only happens at a coding end that was *not* trimmed, so the
    coding-end distributions are pinned to negative values -- the convention
    for "extend rather than trim" -- at the V 3' and J 5' ends.

    Shared because a P fixture is needed both to assert the palindrome truth
    and to corrupt it: two fixtures would let the mutant drift away from the
    thing it is supposed to be a mutant of.
    """
    return (recombined(locus=locus, **recombine)
            .process_coding_ends(v_3=[(-3, 1.0)], j_5=[(-2, 1.0)],
                                 reference_ends_are_coding_ends=True))


def with_p_nucleotides_at_every_end(locus: str = DEFAULT_LOCUS, **recombine: Any):
    """As :func:`with_p_nucleotides`, extended to the two D coding ends.

    The V-and-J fixture produces no D P-addition at all -- measured 0 runs
    over 60 records -- so three of the four recorded P lengths had nothing to
    compare and were checked only against the letter case the same pass
    writes. Extending both D ends gives all four boundaries a run on every
    record (60 of 60 each).
    """
    return (recombined(locus=locus, **recombine)
            .process_coding_ends(v_3=[(-3, 1.0)], d_5=[(-2, 1.0)],
                                 d_3=[(-2, 1.0)], j_5=[(-2, 1.0)],
                                 reference_ends_are_coding_ends=True))
