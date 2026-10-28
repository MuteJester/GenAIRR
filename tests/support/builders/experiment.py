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
