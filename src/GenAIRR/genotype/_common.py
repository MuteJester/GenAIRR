"""Shared leaf primitives for the genotype package.

The canonical V/D/J segment tuple and the reference-allele accessor live
here, imported by both :mod:`.model` and :mod:`._sampling`, so there is a
single definition (no duplication) and no ``model`` <-> ``_sampling``
import cycle.
"""
from __future__ import annotations

from typing import Dict, List

_SEGMENTS = ("V", "D", "J")


def _alleles_by_gene(cfg, segment: str) -> Dict[str, List]:
    return {
        "V": cfg.v_alleles,
        "D": cfg.d_alleles,
        "J": cfg.j_alleles,
    }[segment] or {}
