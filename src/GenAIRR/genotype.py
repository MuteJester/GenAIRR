"""Per-individual diploid genotype builder (PR1: known reference alleles).

A :class:`Genotype` is an editable, narrowed view of a ``DataConfig``'s
reference alleles. Attach it to an experiment with
``Experiment.with_genotype(g)`` to make V(D)J recombination
haplotype-phased: V, D and J of each rearrangement are drawn from a single
chromosome, honouring presence/absence, zygosity, and copy-number/deletion.

Default construction (:meth:`Genotype.from_dataconfig`) is **strict**: a gene
that is used during recombination but never assigned here is an error at
attach time. Use :meth:`complete_from_reference` to fill unspecified genes
with a valid diploid state. :meth:`Genotype.permissive` is a separate,
explicitly non-diploid fallback (see its docstring).
"""
from __future__ import annotations

import math
from collections.abc import Mapping
from typing import Dict, List, Optional, Set, Tuple

_SEGMENTS = ("V", "D", "J")


def _alleles_by_gene(cfg, segment: str) -> Dict[str, List]:
    return {
        "V": cfg.v_alleles,
        "D": cfg.d_alleles,
        "J": cfg.j_alleles,
    }[segment] or {}


class Genotype:
    """A diploid genotype over a ``DataConfig``'s reference alleles."""

    def __init__(self, cfg, *, permissive: bool = False):
        self._cfg = cfg
        self._permissive = bool(permissive)
        self.subject_id: Optional[str] = None  # plain attribute (read anywhere)
        self._chromosome_weights: Tuple[float, float] = (0.5, 0.5)
        # segment -> gene -> [hap0 list[(allele, copies, weight)], hap1 ...]
        self._slots: Dict[str, Dict[str, List[List[Tuple[str, int, float]]]]] = {
            s: {} for s in _SEGMENTS
        }
        # Private/novel alleles defined on this individual:
        # name -> {"allele": <Allele obj>, "gene": str, "segment": str,
        #          "base": str, "mutations": list}
        self._novel: Dict[str, Dict] = {}
        self._source_hash: str = cfg.cartridge_manifest()["hashes"]["refdata_content_hash"]

    # ── constructors ──────────────────────────────────────────────
    @classmethod
    def from_dataconfig(cls, cfg) -> "Genotype":
        """A strict genotype: unspecified-but-used genes error at attach."""
        return cls(cfg, permissive=False)

    @classmethod
    def permissive(cls, cfg) -> "Genotype":
        """NOT a biological diploid genotype. A reference-wide fallback:
        unspecified genes sample over ALL reference alleles WITHOUT
        phasing for those genes. Use only to constrain a few genes; never
        mistake this for 'complete genotype from refdata'."""
        return cls(cfg, permissive=True)

    # ── editing ───────────────────────────────────────────────────
    def with_subject(self, sid: str) -> "Genotype":
        self.subject_id = str(sid)
        return self

    @staticmethod
    def _check_chromosome_weights(w0, w1) -> Tuple[float, float]:
        for w in (w0, w1):
            if isinstance(w, bool) or not isinstance(w, (int, float)) or not math.isfinite(w):
                raise ValueError(f"chromosome_weights must be finite numbers, got {(w0, w1)}")
        if w0 < 0 or w1 < 0 or (w0 + w1) <= 0:
            raise ValueError(
                f"chromosome_weights must be non-negative and sum>0, got {(w0, w1)}"
            )
        return (float(w0), float(w1))

    def chromosome_weights(self, w0: float, w1: float) -> "Genotype":
        self._chromosome_weights = self._check_chromosome_weights(w0, w1)
        return self

    def _check_allele(self, segment: str, gene: str, allele: str) -> None:
        by_gene = _alleles_by_gene(self._cfg, segment)
        if gene not in by_gene:
            raise ValueError(f"{gene!r} is not a known {segment} gene in this cartridge")
        names = {a.name for a in by_gene[gene]}
        # also accept novel alleles defined on this genotype for this gene
        names |= {
            n
            for n, info in self._novel.items()
            if info["gene"] == gene and info["segment"] == segment
        }
        if allele not in names:
            raise ValueError(
                f"{allele!r} is not a known allele of {gene!r} "
                f"(define novel alleles with add_novel_allele first)"
            )

    def homozygous(self, gene: str, allele: str, segment: str = "V") -> "Genotype":
        self._check_allele(segment, gene, allele)
        self._slots[segment][gene] = [[(allele, 1, 1.0)], [(allele, 1, 1.0)]]
        return self

    def heterozygous(
        self, gene: str, allele0: str, allele1: str, segment: str = "V"
    ) -> "Genotype":
        self._check_allele(segment, gene, allele0)
        self._check_allele(segment, gene, allele1)
        self._slots[segment][gene] = [[(allele0, 1, 1.0)], [(allele1, 1, 1.0)]]
        return self

    def delete_gene(self, gene: str, haplotype="both", segment: str = "V") -> "Genotype":
        if haplotype not in ("both", 0, 1):
            raise ValueError(
                f"haplotype must be 'both', 0, or 1, got {haplotype!r}"
            )
        # One-haplotype (hemizygous) deletion requires the gene to be
        # specified first, otherwise the *other* haplotype would also be
        # empty — silently producing a full deletion that
        # complete_from_reference() then skips. Full ("both") deletion of
        # an unspecified gene is fine.
        if haplotype != "both" and gene not in self._slots[segment]:
            raise ValueError(
                f"specify {segment} gene {gene!r} (homozygous/heterozygous) before "
                f"deleting one haplotype; deleting a single haplotype of an "
                f"unspecified gene would delete both"
            )
        cur = self._slots[segment].get(gene, [[], []])
        # copy to avoid aliasing if the slot was shared
        cur = [list(cur[0]), list(cur[1])]
        if haplotype in ("both", 0):
            cur[0] = []
        if haplotype in ("both", 1):
            cur[1] = []
        self._slots[segment][gene] = cur
        return self

    def duplicate_gene(
        self, gene: str, alleles: List[str], haplotype: int, segment: str = "V"
    ) -> "Genotype":
        if haplotype not in (0, 1):
            raise ValueError(f"haplotype must be 0 or 1, got {haplotype!r}")
        for a in alleles:
            self._check_allele(segment, gene, a)
        cur = self._slots[segment].get(gene, [[], []])
        cur = [list(cur[0]), list(cur[1])]
        cur[haplotype] = [(a, 1, 1.0) for a in alleles]
        self._slots[segment][gene] = cur
        return self

    # ── novel / private alleles ───────────────────────────────────
    def _find_ref_allele(self, segment: str, name: str):
        for alleles in _alleles_by_gene(self._cfg, segment).values():
            for a in alleles:
                if a.name == name:
                    return a
        return None

    def add_novel_allele(
        self,
        name: str,
        *,
        base: str,
        mutations: Optional[List[Tuple[int, str]]] = None,
        sequence: Optional[str] = None,
        segment: str = "V",
        allow_nonfunctional: bool = False,
    ) -> "Genotype":
        """Define a private/novel allele not present in the reference.

        Derive it from a reference ``base`` allele by either applying point
        ``mutations`` (a list of ``(0-based position, base)``) or supplying
        an explicit same-length ``sequence`` (substitutions only — no
        indels — so the gene's reading frame and the conserved-anchor
        position are preserved). The novel allele's **gene is taken from
        its name** and must equal the base allele's gene; sub-regions and
        the anchor position are inherited (valid for same-length variants).

        The synthesized coding sequence is **validated**: for V/J the
        conserved anchor codon must still encode the conserved residue
        (Cys for V, Trp/Phe for J) and the coding frame must contain no
        internal stop codon. A variant that breaks either is rejected
        unless ``allow_nonfunctional=True`` (in which case it is kept and
        marked non-functional).

        Registers the novel allele under ``name`` so it can then be placed
        with :meth:`homozygous` / :meth:`heterozygous` / :meth:`duplicate_gene`.
        At compile time it is injected as a real entry in an *effective*
        reference, so it flows through alignment and AIRR output like a
        catalogue allele.
        """
        import copy as _copy

        from .utilities.misc import translate

        base_allele = self._find_ref_allele(segment, base)
        if base_allele is None:
            raise ValueError(f"base allele {base!r} not found in {segment} reference")
        # Gene identity comes from the name; it must match the base's gene
        # (no cross-gene synthesis — that would corrupt gene identity).
        name_gene = name.split("*")[0]
        if name_gene != base_allele.gene:
            raise ValueError(
                f"novel name {name!r} implies gene {name_gene!r} but base {base!r} "
                f"belongs to gene {base_allele.gene!r}; a novel allele must belong "
                f"to its base allele's gene"
            )
        gene = base_allele.gene
        # Name must be unique across the WHOLE catalogue (all segments) and
        # all previously-defined novel alleles.
        if name in self._novel:
            raise ValueError(f"novel allele name {name!r} already defined")
        for seg in _SEGMENTS:
            if self._find_ref_allele(seg, name) is not None:
                raise ValueError(f"novel allele name {name!r} collides with a catalogue allele")
        if (mutations is None) == (sequence is None):
            raise ValueError("provide exactly one of `mutations` or `sequence`")

        base_ungapped = base_allele.ungapped_seq.upper()
        gapped = list(base_allele.gapped_seq or "")
        # ungapped index -> gapped index (positions of non-gap characters).
        # Some custom cartridges carry no (or inconsistent) gapped sequence;
        # in that case we can't project onto gaps, so fall back to an
        # ungapped novel sequence (no gap-derived metadata).
        ung_to_gap = [i for i, ch in enumerate(base_allele.gapped_seq or "") if ch != "."]
        project_gaps = len(ung_to_gap) == len(base_ungapped)
        seq = list(base_ungapped)
        if sequence is not None:
            sequence = sequence.upper()
            if len(sequence) != len(seq):
                raise ValueError(
                    f"explicit sequence length {len(sequence)} != base length {len(seq)}; "
                    "novel alleles are substitution-only (same length as the base)"
                )
            if any(b not in "ACGT" for b in sequence):
                raise ValueError("sequence must contain only A/C/G/T")
            seq = list(sequence)
        else:
            if not mutations:
                raise ValueError("`mutations` must be a non-empty list of (position, base)")
            for pos, b in mutations:
                if not isinstance(pos, int) or isinstance(pos, bool):
                    raise ValueError(f"mutation position must be an int, got {pos!r}")
                if not (isinstance(b, str) and len(b) == 1):
                    raise ValueError(f"mutation base must be a single character, got {b!r}")
                if not (0 <= pos < len(seq)):
                    raise ValueError(f"mutation position {pos} out of range [0,{len(seq)})")
                if b.upper() not in "ACGT":
                    raise ValueError(f"mutation base {b!r} must be A/C/G/T")
                seq[pos] = b.upper()
        new_ungapped = "".join(seq)
        if new_ungapped == base_ungapped:
            raise ValueError("novel allele is identical to its base allele")
        # Project the substitutions onto the gapped sequence too, so
        # gap-dependent metadata stays consistent. If the base has no
        # usable gapped sequence, fall back to the ungapped form.
        if project_gaps:
            for k, b in enumerate(seq):
                gapped[ung_to_gap[k]] = b
            new_gapped = "".join(gapped)
        else:
            new_gapped = new_ungapped

        novel = _copy.deepcopy(base_allele)
        novel.name = name
        novel.gene = gene
        novel.ungapped_seq = new_ungapped
        novel.gapped_seq = new_gapped
        if hasattr(novel, "ungapped_len"):
            novel.ungapped_len = len(new_ungapped)

        # Functional validation (V/J have a conserved coding frame).
        functional, reason = True, None
        anchor = getattr(novel, "anchor", None)
        if segment in ("V", "J") and anchor is not None:
            conserved = {"V": {"C"}, "J": {"W", "F"}}[segment]
            anchor_aa = translate(new_ungapped[anchor : anchor + 3])
            if anchor_aa not in conserved:
                functional = False
                reason = (
                    f"conserved anchor codon now encodes {anchor_aa!r}, "
                    f"expected one of {sorted(conserved)}"
                )
            coding = new_ungapped[:anchor] if segment == "V" else new_ungapped[anchor:]
            if "*" in translate(coding):
                reason = (reason + "; " if reason else "") + "internal stop codon in coding frame"
                functional = False
        if not functional and not allow_nonfunctional:
            raise ValueError(
                f"novel allele {name!r} is non-functional ({reason}); pass "
                f"allow_nonfunctional=True to keep it anyway"
            )
        if not functional:
            try:
                novel.functional_status = "pseudogene"
            except Exception:
                pass

        self._novel[name] = {
            "allele": novel,
            "gene": gene,
            "segment": segment,
            "base": base,
            "mutations": list(mutations) if mutations else None,
            "functional": functional,
        }
        return self

    def has_novel(self) -> bool:
        return bool(self._novel)

    def novel_allele_names(self) -> Set[str]:
        return set(self._novel)

    def _carried_allele_names(self) -> Set[str]:
        """Allele names actually placed on a haplotype (across segments)."""
        names: Set[str] = set()
        for seg in _SEGMENTS:
            for haps in self._slots[seg].values():
                for hap in haps:
                    names.update(a for (a, _c, _w) in hap)
        return names

    def effective_dataconfig(self):
        """Return a copy of the source ``DataConfig`` with this genotype's
        **carried** novel alleles appended to their genes' allele lists —
        the reference the engine actually runs against when novel alleles
        are present. A novel allele that was defined but never placed on a
        haplotype is NOT injected (it would otherwise pollute the aligner
        reference and could surface in ``v_call`` despite being absent from
        the ground truth)."""
        import copy as _copy

        cfg = _copy.deepcopy(self._cfg)
        by_seg = {"V": cfg.v_alleles, "D": cfg.d_alleles, "J": cfg.j_alleles}
        carried = self._carried_allele_names()
        for name, info in self._novel.items():
            if name not in carried:
                continue
            d = by_seg[info["segment"]]
            existing = list(d.get(info["gene"], []))
            existing.append(_copy.deepcopy(info["allele"]))
            d[info["gene"]] = existing
        return cfg

    def complete_from_reference(
        self, policy: str = "homozygous_first_reference"
    ) -> "Genotype":
        """Fill every UNspecified gene with a valid diploid state.

        ``policy``:
        - ``"homozygous_first_reference"`` (default): each unspecified
          gene becomes homozygous for its **first cartridge allele** (NOT
          a population-frequency-common allele — there is no frequency
          prior in PR1; the name says exactly what it does).
        - ``"heterozygous_first_two"``: first two cartridge alleles, one
          per haplotype (homozygous if the gene has a single allele).
        """
        for seg in _SEGMENTS:
            for gene, allele_objs in _alleles_by_gene(self._cfg, seg).items():
                if gene in self._slots[seg] or not allele_objs:
                    continue
                names = [a.name for a in allele_objs]
                if policy == "homozygous_first_reference":
                    self.homozygous(gene, names[0], segment=seg)
                elif policy == "heterozygous_first_two":
                    if len(names) >= 2:
                        self.heterozygous(gene, names[0], names[1], segment=seg)
                    else:
                        self.homozygous(gene, names[0], segment=seg)
                else:
                    raise ValueError(f"unknown policy {policy!r}")
        return self

    # ── snapshot ──────────────────────────────────────────────────
    def _snapshot(self) -> "Genotype":
        """Return an independent copy for attachment to an experiment, so
        that mutating the builder after ``with_genotype()``/``compile()``
        cannot desync ``result.genotypes`` from the compiled engine
        genotype. Shares the (immutable) cartridge reference; deep-copies
        the editable slot state."""
        import copy as _copy

        g = Genotype.__new__(Genotype)
        g._cfg = self._cfg
        g._permissive = self._permissive
        g.subject_id = self.subject_id
        g._chromosome_weights = self._chromosome_weights
        g._slots = _copy.deepcopy(self._slots)
        g._novel = _copy.deepcopy(self._novel)
        g._source_hash = self._source_hash
        return g

    # ── population sampling ───────────────────────────────────────
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
                for nm, w in supplied.items():
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
        haplotype_deletion_prob=0.0,
        segments_to_sample=None,
        chromosome_weights: Tuple[float, float] = (0.5, 0.5),
        subject_id: Optional[str] = None,
        ensure_viable: bool = True,
        max_resamples: int = 1000,
    ) -> "Genotype":
        """Sample a fully-specified diploid genotype from population priors.

        Independent per-gene, per-chromosome Hardy-Weinberg model: each gene on
        each chromosome is independently deleted (prob
        ``haplotype_deletion_prob``) or assigned one allele drawn from the gene's
        allele frequencies. Homozygous/heterozygous/hemizygous/deleted states
        emerge at Hardy-Weinberg rates.

        This is NOT a population haplotype model — no linkage disequilibrium, gene
        co-deletion blocks, ancestry, or donor-specific haplotype structure. It
        samples catalogue alleles only (no novel alleles) and deletion only (no
        copy-number duplication). The default prior is uniform within each gene;
        supply ``allele_frequencies`` for realistic per-gene frequencies.

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
        cw = cls._check_chromosome_weights(*chromosome_weights)
        segs = cls._resolve_sample_segments(cfg, segments_to_sample)
        freqs = cls._resolve_allele_frequencies(cfg, allele_frequencies, segs)
        delp = cls._resolve_haplotype_deletion(cfg, haplotype_deletion_prob, segs)
        # Compute the cartridge content hash ONCE (it rebuilds refdata + hashes);
        # reuse it across all draws instead of recomputing per attempt.
        source_hash = cfg.cartridge_manifest()["hashes"]["refdata_content_hash"]
        # Derive each attempt's sub-seed from a base RNG so a failed attempt at
        # `seed` cannot collide with a direct draw at `seed + 1`.
        base_rng = random.Random(seed)

        attempts = max_resamples if ensure_viable else 1
        for _attempt in range(attempts):
            sub_seed = base_rng.getrandbits(63)
            g = cls._draw_one(cfg, sub_seed, segs, freqs, delp, cw, subject_id, source_hash)
            if not ensure_viable or g._is_viable(cfg, cw):
                return g
        raise ValueError(
            f"could not sample a viable genotype after {max_resamples} attempts; "
            f"haplotype_deletion_prob is too high to leave a complete, positive-weight "
            f"haplotype for required segments {cls._required_segments(cfg)} "
            f"(chromosome_weights={cw})"
        )

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

    def _is_viable(self, cfg, chromosome_weights=None) -> bool:
        """A genotype is viable iff at least one chromosome that can actually be
        expressed (positive chromosome weight) carries every required segment.
        A complete haplotype on a zero-weight chromosome does NOT count — the
        engine would never draw it."""
        cw = chromosome_weights if chromosome_weights is not None else self._chromosome_weights
        for c in (0, 1):
            if cw[c] <= 0:
                continue
            if all(
                any(self._slots[seg].get(gene, [[], []])[c] for gene in self._slots[seg])
                for seg in self._required_segments(cfg)
            ):
                return True
        return False

    # ── queries / export ──────────────────────────────────────────
    @property
    def is_permissive(self) -> bool:
        return self._permissive

    def is_specified(self, segment: str, gene: str) -> bool:
        return gene in self._slots[segment]

    def carried_alleles(self, segment: str, gene: str) -> Set[str]:
        out: Set[str] = set()
        for hap in self._slots[segment].get(gene, [[], []]):
            out.update(a for (a, _c, _w) in hap)
        return out

    @staticmethod
    def _zygosity(h0: List, h1: List) -> str:
        s0 = {a for (a, _, _) in h0}
        s1 = {a for (a, _, _) in h1}
        if not s0 and not s1:
            return "deleted"
        if bool(s0) != bool(s1):  # exactly one haplotype carries the gene
            return "hemizygous"
        if s0 == s1 and len(s0) == 1:
            return "homozygous"
        return "heterozygous"

    def to_table(self) -> List[Dict]:
        """One row per (segment, gene) with full diploid truth: zygosity
        (incl. ``hemizygous`` / ``deleted``), the carried alleles per
        haplotype, and per-haplotype copy/weight detail. Suitable as a
        ground-truth genotype table for inference benchmarks."""
        rows = []
        for seg in _SEGMENTS:
            for gene, haps in self._slots[seg].items():
                h0, h1 = haps[0], haps[1]
                carried = {a for (a, _, _) in h0} | {a for (a, _, _) in h1}
                novel_here = sorted(carried & set(self._novel))
                rows.append(
                    {
                        "subject_id": self.subject_id,
                        "segment": seg,
                        "gene": gene,
                        "zygosity": self._zygosity(h0, h1),
                        "haplotype_0": sorted(a for (a, _, _) in h0),
                        "haplotype_1": sorted(a for (a, _, _) in h1),
                        # per-haplotype (allele, copies, weight) detail
                        "haplotype_0_detail": sorted(h0),
                        "haplotype_1_detail": sorted(h1),
                        "novel": novel_here,  # carried alleles that are private/novel
                        "permissive": self._permissive,
                    }
                )
        return rows

    def to_tsv(self, path: str) -> None:
        import csv

        def _fmt(detail):
            # allele:copies:weight ; ...
            return ";".join(f"{a}:{c}:{w}" for (a, c, w) in detail)

        rows = self.to_table()
        with open(path, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(
                [
                    "subject_id",
                    "segment",
                    "gene",
                    "zygosity",
                    "haplotype_0",
                    "haplotype_1",
                    "novel",
                    "permissive",
                ]
            )
            for r in rows:
                w.writerow(
                    [
                        r["subject_id"],
                        r["segment"],
                        r["gene"],
                        r["zygosity"],
                        _fmt(r["haplotype_0_detail"]),
                        _fmt(r["haplotype_1_detail"]),
                        ";".join(r["novel"]),
                        r["permissive"],
                    ]
                )
