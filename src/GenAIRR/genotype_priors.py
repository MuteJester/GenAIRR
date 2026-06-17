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

    def validate(self, chain_type=None, *, name: str = "genotype_priors") -> None:
        """Catalogue-free shape + numeric/DNA sanity. Raises ``ValueError`` on
        the first violation. ``chain_type`` rejects any D entry on a VJ chain;
        it accepts ``"vj"`` / ``"vdj"`` strings or a ``ChainType``-like object
        exposing ``has_d``."""
        ct = self._normalize_chain_type(chain_type)
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
                if not isinstance(gene, str) or not gene:
                    raise ValueError(
                        f"{name}.allele_frequencies[{seg!r}]: gene names must be "
                        f"non-empty strings, got {gene!r}")
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
                if not isinstance(gene, str) or not gene:
                    raise ValueError(
                        f"{name}.haplotype_deletion_prob[{seg!r}]: gene names must be "
                        f"non-empty strings, got {gene!r}")
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
            if not isinstance(nv.allow_nonfunctional, bool):
                raise ValueError(
                    f"{name}.novel_alleles[{nv.name!r}].allow_nonfunctional must be a "
                    f"bool, got {type(nv.allow_nonfunctional).__name__}")
            if nv.name in seen:
                raise ValueError(f"{name}.novel_alleles: duplicate novel name {nv.name!r}")
            seen.add(nv.name)

    @classmethod
    def from_genotypes(cls, genotypes, *, cfg=None, pseudocount=0.0,
                       include_novel=True, min_subjects=None,
                       subject_id_policy="require_unique", segments=None,
                       model_id="estimated", source="from_genotypes",
                       description="", version="") -> "PopulationGenotypeModel":
        """Estimate a population genotype model from observed ``Genotype`` objects.

        Allele frequencies are counted **per carried chromosome** (homozygous=2,
        hemizygous=1, deleted=0); ``pseudocount`` is added to every catalogue
        allele per gene (and to a novel only when ``include_novel`` and the novel
        belongs to the gene). Deletion probability is per haplotype:
        ``deleted_haplotypes / (2 * n_subjects)`` (NO pseudocount). Genotypes
        carrying a duplicated gene (copy-number > 1) are rejected — the plane is
        deletion-only. ``subject_id_policy`` ('require_unique' | 'allow_duplicates')
        guards against double-counting; ``min_subjects`` guards tiny estimates.
        """
        gts = list(genotypes)
        if subject_id_policy not in ("require_unique", "allow_duplicates"):
            raise ValueError(
                f"subject_id_policy must be 'require_unique' or 'allow_duplicates', "
                f"got {subject_id_policy!r}")
        if subject_id_policy == "require_unique":
            ids = [g.subject_id for g in gts]
            none_count = sum(1 for i in ids if i is None)
            if 0 < none_count < len(ids):
                raise ValueError(
                    "subject_id_policy='require_unique': some genotypes have a "
                    "subject_id and others don't")
            if none_count == len(ids) and len(ids) > 1:
                # all-None: uniqueness cannot be verified, so double-counting
                # cannot be ruled out — fail loudly rather than silently allow it.
                raise ValueError(
                    "subject_id_policy='require_unique': none of the genotypes have a "
                    "subject_id, so uniqueness cannot be verified; set with_subject(...) "
                    "or pass subject_id_policy='allow_duplicates'")
            present = [i for i in ids if i is not None]
            if len(present) != len(set(present)):
                raise ValueError(
                    "subject_id_policy='require_unique': duplicate subject_id found")
        n = len(gts)
        if n == 0:
            raise ValueError("from_genotypes: need at least one genotype")
        if min_subjects is not None and n < min_subjects:
            raise ValueError(
                f"from_genotypes: min_subjects={min_subjects} but only {n} given")

        if cfg is None:
            cfg = gts[0]._cfg
        seg_list = segments if segments is not None else cls._segments_for(cfg)

        by_seg = {"V": cfg.v_alleles, "D": cfg.d_alleles, "J": cfg.j_alleles}

        # reject duplicated genes (copy-number > 1 on any slot)
        for g in gts:
            for seg in _SEGMENTS:
                for gene, haps in g._slots[seg].items():
                    for hap in haps:
                        if len(hap) > 1:
                            raise ValueError(
                                f"from_genotypes: genotype {g.subject_id!r} carries a "
                                f"duplicated gene ({seg} {gene}, copy-number > 1); the "
                                f"population plane is deletion-only — remove duplications "
                                f"or collapse them before estimating")

        # Require completeness: a gene absent from _slots is UNSPECIFIED (unknown),
        # NOT deleted. Counting absence as deletion silently inflates p_del, so we
        # require every estimated-segment catalogue gene to be specified (call
        # complete_from_reference() to fill, or delete_gene() to mark absence).
        for g in gts:
            for seg in seg_list:
                for gene in (by_seg[seg] or {}):
                    if gene not in g._slots[seg]:
                        raise ValueError(
                            f"from_genotypes: genotype {g.subject_id!r} does not specify "
                            f"{seg} gene {gene!r}; an unspecified gene is unknown, not "
                            f"deleted. Call complete_from_reference() (or delete_gene to "
                            f"mark it absent) before estimating")

        freqs: Dict[str, Dict[str, Dict[str, float]]] = {}
        dele: Dict[str, Dict[str, float]] = {}
        novel_freq: Dict[str, float] = {}
        novel_spec: Dict[str, "PopulationNovelAllele"] = {}

        for seg in seg_list:
            catalogue = by_seg[seg] or {}
            freqs[seg] = {}
            dele[seg] = {}
            for gene, alleles in catalogue.items():
                counts = {a.name: 0.0 for a in alleles}
                deleted_haps = 0
                for g in gts:
                    haps = g._slots[seg].get(gene, [[], []])
                    for hap in haps:
                        names = {a for (a, _c, _w) in hap}
                        if not names:
                            deleted_haps += 1
                            continue
                        for nm in names:
                            if nm in counts:
                                counts[nm] += 1.0
                            elif include_novel:
                                novel_freq[nm] = novel_freq.get(nm, 0.0) + 1.0
                if pseudocount:
                    for nm in counts:
                        counts[nm] += float(pseudocount)
                kept = {nm: c for nm, c in counts.items() if c > 0}
                if kept:
                    freqs[seg][gene] = kept
                dele[seg][gene] = deleted_haps / (2.0 * n)

        novels: List[PopulationNovelAllele] = []
        if include_novel:
            for g in gts:
                for name, info in getattr(g, "_novel", {}).items():
                    if name in novel_freq and name not in novel_spec:
                        novel_spec[name] = PopulationNovelAllele(
                            name=name, segment=info["segment"],
                            base_allele=info["base"],
                            sequence=info["allele"].ungapped_seq.upper(),
                            frequency=novel_freq[name],
                            allow_nonfunctional=not info.get("functional", True))
            # Novel frequency is carried on the PopulationNovelAllele itself, NOT
            # duplicated into allele_frequencies — putting a non-catalogue name
            # there would (a) fail set_genotype_priors' catalogue-aware check and
            # (b) collide with sample-time novel synthesis. sample() combines the
            # catalogue table with novel frequencies via _augment_freqs_with_novels.
            novels.extend(novel_spec.values())

        return cls(allele_frequencies=freqs, haplotype_deletion_prob=dele,
                   chromosome_weights=(0.5, 0.5), novel_alleles=novels,
                   model_id=model_id, source=source, description=description,
                   version=version)

    @staticmethod
    def _segments_for(cfg):
        segs = ["V"]
        if getattr(cfg, "d_alleles", None):
            segs.append("D")
        segs.append("J")
        return segs

    @staticmethod
    def _normalize_chain_type(chain_type):
        """Return ``"vj"`` / ``"vdj"`` / ``None`` from a string or a
        ``ChainType``-like object (one exposing ``has_d``)."""
        if chain_type is None:
            return None
        if isinstance(chain_type, str):
            return chain_type.lower()
        has_d = getattr(chain_type, "has_d", None)
        if has_d is None:
            return None
        return "vdj" if has_d else "vj"

    @staticmethod
    def _check_segment(seg, ct, where):
        if seg not in _SEGMENTS:
            raise ValueError(f"{where}: segment {seg!r} must be one of {_SEGMENTS}")
        if ct == "vj" and seg == "D":
            raise ValueError(
                f"{where}: D-segment prior on a VJ chain is meaningless (no D pool). "
                f"Drop the D entries or use a VDJ cartridge.")

    def content_checksum(self) -> str:
        """Canonical sha256 of the plane's semantic content — stable across dict
        insertion order, int-vs-float spelling, and DNA case. This (not the
        pickle-based DataConfig checksum) is what manifest/provenance report."""

        def _num(x):
            return repr(float(x))

        def _freqs(d):
            return {seg: {gene: {al: _num(w) for al, w in sorted(alleles.items())}
                          for gene, alleles in sorted(genes.items())}
                    for seg, genes in sorted(d.items())}

        def _del(d):
            return {seg: {gene: _num(p) for gene, p in sorted(genes.items())}
                    for seg, genes in sorted(d.items())}

        novels = sorted(
            ({"name": nv.name, "segment": nv.segment, "base_allele": nv.base_allele,
              "sequence": nv.sequence.upper(), "frequency": _num(nv.frequency),
              "allow_nonfunctional": bool(nv.allow_nonfunctional)}
             for nv in self.novel_alleles),
            key=lambda r: r["name"],
        )
        canon = {
            "schema": SCHEMA_TAG,
            "model_id": self.model_id,
            "source": self.source,
            "description": self.description,
            "version": self.version,
            "allele_frequencies": _freqs(self.allele_frequencies),
            "haplotype_deletion_prob": _del(self.haplotype_deletion_prob),
            "chromosome_weights": [_num(self.chromosome_weights[0]), _num(self.chromosome_weights[1])],
            "novel_alleles": novels,
        }
        blob = json.dumps(canon, sort_keys=True, separators=(",", ":")).encode("utf-8")
        return hashlib.sha256(blob).hexdigest()
