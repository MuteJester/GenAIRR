"""Population genotype model — a donor-population germline prior carried as a
top-level cartridge plane (``DataConfig.genotype_priors``).

This is NOT an empirical recombination model (those live on
``ReferenceEmpiricalModels``); it is a per-gene carriage/deletion prior plus a
catalogue of population novel/private alleles. ``Genotype.sample`` consumes it to
draw a per-individual diploid genotype. See
``.private/specs/2026-06-17-genotype-cartridge-plane-design.md``.

The plane stays decoupled from a specific catalogue: ``validate()`` does shape +
numeric/DNA sanity only. Catalogue-aware checks (gene/allele existence, novel
functional validation, viability) run at attach time
(``ReferenceCartridgeBuilder.set_genotype_priors``) and at sample time
(``Genotype.sample``).
"""
from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

_SEGMENTS = ("V", "D", "J")
_DNA = set("ACGT")
SCHEMA_TAG = "population_genotype_model/1"


@dataclass
class PopulationNovelAllele:
    """A population novel/private allele carried on a cartridge prior.

    ``name``'s gene token (text before ``"*"``) must equal ``base_allele``'s
    gene (validated catalogue-aware); ``sequence`` is substitution-only (same
    length as the base, validated catalogue-aware). ``frequency`` is the
    population weight the allele competes with inside its base gene's draw.
    """

    name: str
    segment: str
    base_allele: str
    sequence: str
    frequency: float
    allow_nonfunctional: bool = False


def _check_weight(w, where, *, strict_positive):
    if isinstance(w, bool) or not isinstance(w, (int, float)):
        raise ValueError(f"{where}: weight must be a finite number, got {type(w).__name__}")
    wf = float(w)
    if not math.isfinite(wf):
        raise ValueError(f"{where}: weight must be finite, got {w!r}")
    if strict_positive and wf <= 0.0:
        raise ValueError(f"{where}: weight must be > 0, got {wf}")
    if not strict_positive and wf < 0.0:
        raise ValueError(f"{where}: weight must be >= 0, got {wf}")
    return wf


@dataclass
class PopulationGenotypeModel:
    """A donor-population germline prior (see module docstring)."""

    allele_frequencies: Dict[str, Dict[str, Dict[str, float]]] = field(default_factory=dict)
    haplotype_deletion_prob: Dict[str, Dict[str, float]] = field(default_factory=dict)
    chromosome_weights: Tuple[float, float] = (0.5, 0.5)
    novel_alleles: List[PopulationNovelAllele] = field(default_factory=list)
    model_id: str = ""
    source: str = ""
    description: str = ""
    version: str = ""

    def validate(self, chain_type: Optional[str] = None, *, name: str = "genotype_priors") -> None:
        """Catalogue-free shape + numeric/DNA sanity. Raises ``ValueError`` on
        the first violation. ``chain_type="vj"`` rejects any D entry."""
        ct = chain_type.lower() if isinstance(chain_type, str) else None
        if not isinstance(self.model_id, str) or not self.model_id:
            raise ValueError(f"{name}.model_id must be a non-empty string")
        if not isinstance(self.source, str) or not self.source:
            raise ValueError(f"{name}.source must be a non-empty string")
        for label in ("description", "version"):
            if not isinstance(getattr(self, label), str):
                raise ValueError(f"{name}.{label} must be a string")

        # allele_frequencies: seg -> gene -> {allele: weight >= 0}, >=1 positive
        if not isinstance(self.allele_frequencies, dict):
            raise ValueError(f"{name}.allele_frequencies must be a dict")
        for seg, genes in self.allele_frequencies.items():
            self._check_segment(seg, ct, f"{name}.allele_frequencies")
            if not isinstance(genes, dict):
                raise ValueError(f"{name}.allele_frequencies[{seg!r}] must be a dict")
            for gene, alleles in genes.items():
                if not isinstance(alleles, dict) or not alleles:
                    raise ValueError(
                        f"{name}.allele_frequencies[{seg!r}][{gene!r}] must be a "
                        f"non-empty mapping of allele -> weight")
                positive = 0
                for allele, w in alleles.items():
                    if not isinstance(allele, str) or not allele:
                        raise ValueError(
                            f"{name}.allele_frequencies[{seg!r}][{gene!r}]: allele names "
                            f"must be non-empty strings, got {allele!r}")
                    wf = _check_weight(w, f"{name}.allele_frequencies[{seg}][{gene}][{allele}]",
                                       strict_positive=False)
                    if wf > 0:
                        positive += 1
                if positive == 0:
                    raise ValueError(
                        f"{name}.allele_frequencies[{seg!r}][{gene!r}]: at least one "
                        f"allele weight must be > 0")

        # haplotype_deletion_prob: seg -> gene -> prob in [0, 1]
        if not isinstance(self.haplotype_deletion_prob, dict):
            raise ValueError(f"{name}.haplotype_deletion_prob must be a dict")
        for seg, genes in self.haplotype_deletion_prob.items():
            self._check_segment(seg, ct, f"{name}.haplotype_deletion_prob")
            if not isinstance(genes, dict):
                raise ValueError(f"{name}.haplotype_deletion_prob[{seg!r}] must be a dict")
            for gene, p in genes.items():
                if isinstance(p, bool) or not isinstance(p, (int, float)) or not math.isfinite(p):
                    raise ValueError(
                        f"{name}.haplotype_deletion_prob[{seg}][{gene}] must be finite, got {p!r}")
                if not (0.0 <= float(p) <= 1.0):
                    raise ValueError(
                        f"{name}.haplotype_deletion_prob[{seg}][{gene}] must be in [0, 1], got {p}")

        # chromosome_weights: 2-tuple, finite, non-negative, >=1 positive
        cw = self.chromosome_weights
        if not (isinstance(cw, (tuple, list)) and len(cw) == 2):
            raise ValueError(f"{name}.chromosome_weights must be a 2-tuple, got {cw!r}")
        for w in cw:
            if isinstance(w, bool) or not isinstance(w, (int, float)) or not math.isfinite(w):
                raise ValueError(f"{name}.chromosome_weights must be finite numbers, got {cw!r}")
        if cw[0] < 0 or cw[1] < 0 or (cw[0] + cw[1]) <= 0:
            raise ValueError(
                f"{name}.chromosome_weights must be non-negative and sum>0, got {tuple(cw)}")

        # novel_alleles
        seen = set()
        for nv in self.novel_alleles:
            if not isinstance(nv, PopulationNovelAllele):
                raise ValueError(f"{name}.novel_alleles entries must be PopulationNovelAllele")
            for attr in ("name", "base_allele", "sequence"):
                v = getattr(nv, attr)
                if not isinstance(v, str) or not v:
                    raise ValueError(f"{name}.novel_alleles: {attr} must be a non-empty string")
            self._check_segment(nv.segment, ct, f"{name}.novel_alleles[{nv.name!r}]")
            if "*" not in nv.name or not nv.name.split("*")[0]:
                raise ValueError(
                    f"{name}.novel_alleles: name {nv.name!r} must contain a gene token before '*'")
            if any(b not in _DNA for b in nv.sequence.upper()):
                raise ValueError(
                    f"{name}.novel_alleles[{nv.name!r}]: sequence must be DNA (A/C/G/T only)")
            _check_weight(nv.frequency, f"{name}.novel_alleles[{nv.name!r}].frequency",
                          strict_positive=True)
            if nv.name in seen:
                raise ValueError(f"{name}.novel_alleles: duplicate novel name {nv.name!r}")
            seen.add(nv.name)

    @staticmethod
    def _check_segment(seg, ct, where):
        if seg not in _SEGMENTS:
            raise ValueError(f"{where}: segment {seg!r} must be one of {_SEGMENTS}")
        if ct == "vj" and seg == "D":
            raise ValueError(
                f"{where}: D-segment prior on a VJ chain is meaningless (no D pool). "
                f"Drop the D entries or use a VDJ cartridge.")
