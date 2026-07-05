"""Population-sampling engine for :class:`Genotype`, extracted verbatim as a
mixin (behavior-preserving). ``Genotype`` inherits these class/static methods;
``cls`` resolves to ``Genotype`` at call time, so every call site is unchanged.
"""
from __future__ import annotations

import math
from collections.abc import Mapping
from typing import Dict, List, Optional, Tuple

from ._common import _SEGMENTS, _alleles_by_gene


class _GenotypeSampling:
    @staticmethod
    def _required_segments(cfg) -> List[str]:
        req = ["V", "J"]
        if _alleles_by_gene(cfg, "D"):
            req.insert(1, "D")  # V, D, J
        return req

    @classmethod
    def _resolve_sample_segments(cls, cfg, segments_to_sample) -> List[str]:
        required = cls._required_segments(cfg)
        if segments_to_sample is None:
            return required
        segs = list(segments_to_sample)
        seen = set()
        for s in segs:
            if s not in _SEGMENTS:
                raise ValueError(f"unknown segment {s!r}; expected one of {_SEGMENTS}")
            if s in seen:
                raise ValueError(f"segments_to_sample contains duplicate segment {s!r}")
            seen.add(s)
            if not _alleles_by_gene(cfg, s):
                raise ValueError(f"cartridge has no {s} segment")
        missing = [r for r in required if r not in seen]
        if missing:
            raise ValueError(
                f"segments_to_sample must cover the chain's required segments "
                f"{required}; missing {missing}. Partial sampling is not supported "
                f"by Genotype.sample (it must return a runnable genotype)."
            )
        # canonical _SEGMENTS order so the result is independent of input order
        return [s for s in _SEGMENTS if s in seen]

    @classmethod
    def _gene_segment_index(cls, cfg, segs) -> Dict[str, List[str]]:
        """gene name -> [segments it appears in] (for flat-shape disambiguation)."""
        idx: Dict[str, List[str]] = {}
        for seg in segs:
            for gene in _alleles_by_gene(cfg, seg):
                idx.setdefault(gene, []).append(seg)
        return idx

    @staticmethod
    def _weighted_pick(rng, pairs):
        names = [n for (n, _w) in pairs]
        weights = [w for (_n, w) in pairs]
        return rng.choices(names, weights=weights, k=1)[0]

    @classmethod
    def _normalize_freq_spec(cls, cfg, spec, segs):
        """Return nested ``{seg: {gene: {allele: weight}}}`` from a nested or flat
        ``allele_frequencies`` spec, fully validating segment/gene/allele
        addressing and shapes (unknown names and malformed values raise)."""
        if spec is None:
            return {}
        if not isinstance(spec, Mapping):
            raise ValueError(
                "allele_frequencies must be a mapping (or 'usage_as_prior' / None), "
                f"got {type(spec).__name__}"
            )
        nested: Dict[str, Dict[str, Dict[str, float]]] = {}
        keys = set(spec)
        if keys and keys <= set(_SEGMENTS):  # segment-keyed (nested) shape
            for seg, genes in spec.items():
                if seg not in segs:
                    raise ValueError(
                        f"allele_frequencies: segment {seg!r} is not being sampled "
                        f"(sampling {segs})"
                    )
                if not isinstance(genes, Mapping):
                    raise ValueError(
                        f"allele_frequencies[{seg!r}] must be a mapping of "
                        f"gene -> {{allele: weight}}, got {type(genes).__name__}"
                    )
                catalogue = _alleles_by_gene(cfg, seg)
                for gene, alleles in genes.items():
                    if gene not in catalogue:
                        raise ValueError(
                            f"allele_frequencies: {seg} gene {gene!r} is not in the cartridge"
                        )
                    if not isinstance(alleles, Mapping) or not alleles:
                        raise ValueError(
                            f"allele_frequencies[{seg!r}][{gene!r}] must be a non-empty "
                            f"mapping of allele -> weight"
                        )
                    nested.setdefault(seg, {})[gene] = alleles
        else:  # flat {gene: {...}} shape
            gidx = cls._gene_segment_index(cfg, segs)
            for gene, alleles in spec.items():
                segs_for = gidx.get(gene)
                if not segs_for:
                    raise ValueError(f"allele_frequencies: unknown gene {gene!r}")
                if len(segs_for) > 1:
                    raise ValueError(
                        f"allele_frequencies: gene {gene!r} is ambiguous across segments "
                        f"{segs_for}; use the {{segment: {{gene: ...}}}} shape"
                    )
                if not isinstance(alleles, Mapping) or not alleles:
                    raise ValueError(
                        f"allele_frequencies[{gene!r}] must be a non-empty mapping of "
                        f"allele -> weight"
                    )
                nested.setdefault(segs_for[0], {})[gene] = alleles
        return nested

    @classmethod
    def _usage_frequencies(cls, cfg, segs):
        rm = getattr(cfg, "reference_models", None)
        usage = getattr(rm, "allele_usage", None) if rm else None
        if usage is None:
            raise ValueError(
                "allele_frequencies='usage_as_prior' requires a cartridge with a "
                "typed reference_models.allele_usage; this cartridge has none"
            )
        nested: Dict[str, Dict[str, Dict[str, float]]] = {}
        seg_attr = {"V": "v", "D": "d", "J": "j"}
        for seg in segs:
            table = getattr(usage, seg_attr[seg], None) or {}
            if not table:
                raise ValueError(
                    f"allele_frequencies='usage_as_prior': cartridge allele_usage has no "
                    f"entries for requested segment {seg!r}"
                )
            catalogue = _alleles_by_gene(cfg, seg)
            for allele_name, w in table.items():
                gene = allele_name.split("*")[0]
                if gene not in catalogue:
                    raise ValueError(
                        f"allele_frequencies='usage_as_prior': usage allele {allele_name!r} "
                        f"maps to {seg} gene {gene!r}, which is not in the cartridge catalogue"
                    )
                nested.setdefault(seg, {}).setdefault(gene, {})[allele_name] = w
        return nested

    @classmethod
    def _resolve_allele_frequencies(cls, cfg, spec, segs):
        if spec == "usage_as_prior":
            nested = cls._usage_frequencies(cfg, segs)
        else:
            nested = cls._normalize_freq_spec(cfg, spec, segs)
        out: Dict[str, Dict[str, List[Tuple[str, float]]]] = {}
        for seg in segs:
            out[seg] = {}
            for gene, alleles in _alleles_by_gene(cfg, seg).items():
                names = {a.name for a in alleles}
                supplied = nested.get(seg, {}).get(gene)
                if supplied is None:
                    out[seg][gene] = [(a.name, 1.0) for a in alleles]  # uniform fallback
                    continue
                pairs: List[Tuple[str, float]] = []
                total = 0.0
                # Sort by allele name so the draw is independent of the spec's
                # dict insertion order — two content-equal models (same
                # content_checksum) then produce identical draws at a given seed.
                for nm, w in sorted(supplied.items()):
                    if nm not in names:
                        raise ValueError(f"{seg} gene {gene}: {nm!r} is not a known allele")
                    if isinstance(w, bool) or not isinstance(w, (int, float)) or not math.isfinite(w) or w < 0:
                        raise ValueError(
                            f"{seg} gene {gene}: weight for {nm!r} must be finite and >= 0, got {w!r}"
                        )
                    if w > 0:
                        pairs.append((nm, float(w)))
                    total += w
                if total <= 0 or not pairs:
                    raise ValueError(
                        f"{seg} gene {gene}: at least one allele weight must be > 0"
                    )
                out[seg][gene] = pairs
        return out

    @classmethod
    def _resolve_haplotype_deletion(cls, cfg, spec, segs):
        def _check(p, where):
            if isinstance(p, bool) or not isinstance(p, (int, float)) or not math.isfinite(p):
                raise ValueError(f"{where}: deletion probability must be a finite number, got {p!r}")
            if not (0.0 <= p <= 1.0):
                raise ValueError(f"{where}: deletion probability must be in [0, 1], got {p}")
            return float(p)

        out = {seg: {gene: 0.0 for gene in _alleles_by_gene(cfg, seg)} for seg in segs}
        if isinstance(spec, (int, float)) and not isinstance(spec, bool):
            p = _check(spec, "haplotype_deletion_prob")
            for seg in segs:
                for gene in out[seg]:
                    out[seg][gene] = p
            return out
        if not isinstance(spec, Mapping):
            raise ValueError(
                f"haplotype_deletion_prob must be a float or a mapping, got {type(spec).__name__}"
            )
        keys = set(spec)
        if keys and keys <= set(_SEGMENTS):  # nested {seg: {gene: prob}}
            for seg, genes in spec.items():
                if seg not in segs:
                    raise ValueError(f"haplotype_deletion_prob: segment {seg!r} not being sampled")
                if not isinstance(genes, Mapping):
                    raise ValueError(
                        f"haplotype_deletion_prob[{seg!r}] must be a mapping of "
                        f"gene -> probability, got {type(genes).__name__}"
                    )
                for gene, p in genes.items():
                    if gene not in out[seg]:
                        raise ValueError(f"haplotype_deletion_prob: unknown {seg} gene {gene!r}")
                    out[seg][gene] = _check(p, f"haplotype_deletion_prob[{seg}][{gene}]")
        else:  # flat {gene: prob}
            gidx = cls._gene_segment_index(cfg, segs)
            for gene, p in spec.items():
                segs_for = gidx.get(gene)
                if not segs_for:
                    raise ValueError(f"haplotype_deletion_prob: unknown gene {gene!r}")
                if len(segs_for) > 1:
                    raise ValueError(
                        f"haplotype_deletion_prob: gene {gene!r} is ambiguous across "
                        f"segments {segs_for}; use the {{segment: {{gene: prob}}}} shape"
                    )
                out[segs_for[0]][gene] = _check(p, f"haplotype_deletion_prob[{gene}]")
        return out

    @classmethod
    def sample(
        cls,
        cfg,
        *,
        seed: int = 0,
        allele_frequencies=None,
        haplotype_deletion_prob=None,
        segments_to_sample=None,
        chromosome_weights: Optional[Tuple[float, float]] = None,
        subject_id: Optional[str] = None,
        ensure_viable: bool = True,
        max_resamples: int = 1000,
        use_cartridge_priors: bool = True,
        include_cartridge_novel_alleles="auto",
    ) -> "Genotype":
        """Sample a fully-specified diploid genotype from population priors.

        Independent per-gene, per-chromosome Hardy-Weinberg model: each gene on
        each chromosome is independently deleted or assigned one allele drawn
        from the gene's allele frequencies. Homozygous/heterozygous/hemizygous/
        deleted states emerge at Hardy-Weinberg rates.

        When ``cfg`` carries a ``genotype_priors`` plane and an argument is left
        at its default, the plane supplies it (``use_cartridge_priors=False``
        disables ALL plane consumption — a clean uniform catalogue-only draw).
        Each input is sourced independently and recorded in
        ``g.prior_provenance`` (explicit / cartridge / uniform / default).

        This is NOT a population haplotype model — no linkage disequilibrium, gene
        co-deletion blocks, ancestry, or donor-specific haplotype structure.

        With ``ensure_viable=True`` (default), the draw is repeated (with a
        deterministic sub-seed) up to ``max_resamples`` times until at least one
        **positive-weight** chromosome carries every required segment, raising
        ``ValueError`` if that is impossible under the given deletion settings.

        NOTE: with the default ``ensure_viable=True`` the result is Hardy-Weinberg
        **conditioned on viability** (draws with no complete usable haplotype are
        rejected), not the unconditional HW distribution. Use
        ``ensure_viable=False`` for the raw (possibly infeasible) HW draw.
        """
        import random

        if isinstance(max_resamples, bool) or not isinstance(max_resamples, int) or max_resamples < 1:
            raise ValueError(f"max_resamples must be an int >= 1, got {max_resamples!r}")
        if not isinstance(use_cartridge_priors, bool):
            raise ValueError(
                f"use_cartridge_priors must be a bool, got {use_cartridge_priors!r}")
        # Identity checks (not `in (...)`): Python's `1 == True` / `0 == False`
        # would otherwise let integers slip through and silently act like False.
        _icna = include_cartridge_novel_alleles
        if not (_icna is True or _icna is False or _icna == "auto"):
            raise ValueError(
                "include_cartridge_novel_alleles must be 'auto', True, or False, "
                f"got {include_cartridge_novel_alleles!r}")

        plane = cls._resolve_plane(cfg, use_cartridge_priors)
        if plane is not None:
            # A plane attached via the builder is already validated, but a plane
            # set directly on the DataConfig bypasses that — validate before use so
            # a malformed prior (empty model_id, NaN weights, D-on-VJ) fails loudly
            # rather than being silently sampled and stamped into provenance.
            plane.validate(chain_type=getattr(getattr(cfg, "metadata", None), "chain_type", None))

        # Per-input source resolution.
        if allele_frequencies is not None:
            freq_spec, freq_src = allele_frequencies, "explicit"
        elif plane is not None and plane.allele_frequencies:
            freq_spec, freq_src = plane.allele_frequencies, "cartridge"
        else:
            freq_spec, freq_src = None, "uniform"

        if haplotype_deletion_prob is not None:
            del_spec, del_src = haplotype_deletion_prob, "explicit"
        elif plane is not None and plane.haplotype_deletion_prob:
            del_spec, del_src = plane.haplotype_deletion_prob, "cartridge"
        else:
            del_spec, del_src = 0.0, "default"

        if chromosome_weights is not None:
            cw_in, cw_src = chromosome_weights, "explicit"
        elif plane is not None:
            cw_in, cw_src = plane.chromosome_weights, "cartridge"
        else:
            cw_in, cw_src = (0.5, 0.5), "default"

        cw = cls._check_chromosome_weights(*cw_in)
        segs = cls._resolve_sample_segments(cfg, segments_to_sample)

        # Candidate-novel injection. Plane novels are draw CANDIDATES (injected
        # into the sampling cfg) only when their source is active; the returned
        # genotype still exports only CARRIED novels via effective_dataconfig().
        inject_novels = False
        if plane is not None and plane.novel_alleles:
            if include_cartridge_novel_alleles is True:
                inject_novels = True
            elif include_cartridge_novel_alleles == "auto":
                inject_novels = freq_src in ("cartridge", "uniform")
            # False -> never
        novel_src = "none"
        sampling_cfg = cfg
        novel_helper = None
        if inject_novels:
            novel_helper = cls._register_plane_novels(cfg, plane)
            sampling_cfg = cls._dataconfig_injecting_all_novels(novel_helper)
            freq_spec = cls._augment_freqs_with_novels(sampling_cfg, freq_spec, plane, segs)
            novel_src = "cartridge"

        freqs = cls._resolve_allele_frequencies(sampling_cfg, freq_spec, segs)
        delp = cls._resolve_haplotype_deletion(sampling_cfg, del_spec, segs)
        # Compute the cartridge content hash ONCE (it rebuilds refdata + hashes);
        # reuse it across all draws instead of recomputing per attempt.
        source_hash = cfg.cartridge_manifest()["hashes"]["refdata_content_hash"]
        # Derive each attempt's sub-seed from a base RNG so a failed attempt at
        # `seed` cannot collide with a direct draw at `seed + 1`.
        base_rng = random.Random(seed)

        provenance = {
            "allele_frequencies": freq_src,
            "haplotype_deletion_prob": del_src,
            "chromosome_weights": cw_src,
            "novel_alleles": novel_src,
            "model_id": plane.model_id if plane is not None else None,
            "model_checksum": plane.content_checksum() if plane is not None else None,
        }

        attempts = max_resamples if ensure_viable else 1
        for _attempt in range(attempts):
            sub_seed = base_rng.getrandbits(63)
            g = cls._draw_one(sampling_cfg, sub_seed, segs, freqs, delp, cw, subject_id, source_hash)
            if not ensure_viable or g._is_viable(cfg, cw):
                cls._rebind_to_base(g, cfg, novel_helper)
                g.prior_provenance = dict(provenance)
                return g
        raise ValueError(
            f"could not sample a viable genotype after {max_resamples} attempts; "
            f"haplotype_deletion_prob is too high to leave a complete, positive-weight "
            f"haplotype for required segments {cls._required_segments(cfg)} "
            f"(chromosome_weights={cw})"
        )

    @staticmethod
    def _resolve_plane(cfg, use_cartridge_priors):
        if not use_cartridge_priors:
            return None
        return getattr(cfg, "genotype_priors", None)

    @classmethod
    def _register_plane_novels(cls, cfg, plane):
        """Build a throwaway Genotype carrying every plane novel as a registered
        novel allele. This is where catalogue-aware + functional validation of
        plane novels happens (via add_novel_allele). Returns the helper."""
        helper = cls.from_dataconfig(cfg)
        for nv in plane.novel_alleles:
            helper.add_novel_allele(
                nv.name, base=nv.base_allele, sequence=nv.sequence.upper(),
                segment=nv.segment, allow_nonfunctional=nv.allow_nonfunctional)
        return helper

    @staticmethod
    def _dataconfig_injecting_all_novels(helper):
        """A cfg copy with ALL of the helper's registered novels injected as
        catalogue alleles — the *sampling* reference (candidates), distinct from
        effective_dataconfig() which injects only carried novels."""
        import copy as _copy
        cfg = _copy.deepcopy(helper._cfg)
        by_seg = {"V": cfg.v_alleles, "D": cfg.d_alleles, "J": cfg.j_alleles}
        for name, info in helper._novel.items():
            d = by_seg[info["segment"]]
            existing = list(d.get(info["gene"], []))
            existing.append(_copy.deepcopy(info["allele"]))
            d[info["gene"]] = existing
        return cfg

    @classmethod
    def _augment_freqs_with_novels(cls, sampling_cfg, freq_spec, plane, segs):
        """Return a nested freq spec that includes each plane novel in its gene's
        table per the synthesis rule: authored table -> preserve + add novel;
        no table -> catalogue alleles 1.0 + novel frequency. Collisions raise."""
        if freq_spec == "usage_as_prior":
            nested = cls._usage_frequencies(sampling_cfg, segs)
            nested = {seg: {g: dict(al) for g, al in genes.items()}
                      for seg, genes in nested.items()}
        else:
            nested = {seg: {g: dict(al) for g, al in genes.items()}
                      for seg, genes in cls._normalize_freq_spec(sampling_cfg, freq_spec, segs).items()}
        for nv in plane.novel_alleles:
            seg = nv.segment
            if seg not in segs:
                continue
            gene = nv.name.split("*")[0]
            gene_tbl = nested.setdefault(seg, {}).get(gene)
            if gene_tbl is None:
                # no authored table for this gene: catalogue 1.0 + novel
                catalogue = _alleles_by_gene(sampling_cfg, seg).get(gene, [])
                gene_tbl = {a.name: 1.0 for a in catalogue if a.name != nv.name}
                nested[seg][gene] = gene_tbl
            if nv.name in gene_tbl:
                raise ValueError(
                    f"plane novel {nv.name!r} collides with an existing allele weight")
            gene_tbl[nv.name] = float(nv.frequency)
        return nested

    @staticmethod
    def _rebind_to_base(g, base_cfg, novel_helper):
        """Point a drawn genotype back at the base cfg and register only the
        novels it actually carries, so effective_dataconfig() injects carried
        novels (and nothing else). Task 7 supplies ``novel_helper``."""
        import copy as _copy
        g._cfg = base_cfg
        if novel_helper is None:
            return
        carried = g._carried_allele_names()
        for name, info in novel_helper._novel.items():
            if name in carried:
                g._novel[name] = _copy.deepcopy(info)

    @classmethod
    def _draw_one(cls, cfg, seed, segs, freqs, delp, cw, subject_id, source_hash):
        import random

        rng = random.Random(seed)
        # Build a bare Genotype directly (avoid from_dataconfig, which recomputes
        # the cartridge hash on every draw); reuse the precomputed source_hash.
        g = cls.__new__(cls)
        g._cfg = cfg
        g._permissive = False
        g.subject_id = subject_id
        g._chromosome_weights = cw
        g._slots = {s: {} for s in _SEGMENTS}
        g._novel = {}
        g._source_hash = source_hash
        g.prior_provenance = cls._manual_provenance()  # sample() overwrites with real sources
        for seg in segs:
            for gene in _alleles_by_gene(cfg, seg):
                pdel = delp[seg][gene]
                slots: List[List[Tuple[str, int, float]]] = [[], []]
                for h in (0, 1):
                    if rng.random() < pdel:
                        continue  # deleted on this chromosome
                    allele = cls._weighted_pick(rng, freqs[seg][gene])
                    slots[h] = [(allele, 1, 1.0)]
                g._slots[seg][gene] = slots
        return g
