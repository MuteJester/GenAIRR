"""Experiments in their commonest shapes.

``Experiment.on("human_igh")`` appears in the suite a few hundred times, usually
followed by ``.recombine()``. Building it here means a change to the default
arrangement is one edit rather than hundreds.
"""
from __future__ import annotations

from typing import Any

import GenAIRR as ga

DEFAULT_LOCUS = "human_igh"


def recombined(locus: str = DEFAULT_LOCUS, **recombine: Any):
    """An experiment bound to ``locus`` with recombination configured.

    Keywords go straight to ``recombine``, so a caller needing particular
    allele weights or trim distributions passes them here.
    """
    return ga.Experiment.on(locus).recombine(**recombine)
