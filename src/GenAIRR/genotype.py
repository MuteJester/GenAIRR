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

    def chromosome_weights(self, w0: float, w1: float) -> "Genotype":
        import math

        if not (math.isfinite(w0) and math.isfinite(w1)):
            raise ValueError(
                f"chromosome_weights must be finite, got {(w0, w1)}"
            )
        if w0 < 0 or w1 < 0 or (w0 + w1) <= 0:
            raise ValueError(
                f"chromosome_weights must be non-negative and sum>0, got {(w0, w1)}"
            )
        self._chromosome_weights = (float(w0), float(w1))
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
        for a in alleles:
            self._check_allele(segment, gene, a)
        cur = self._slots[segment].get(gene, [[], []])
        cur = [list(cur[0]), list(cur[1])]
        cur[int(haplotype)] = [(a, 1, 1.0) for a in alleles]
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
        gapped = list(base_allele.gapped_seq)
        # ungapped index -> gapped index (positions of non-gap characters)
        ung_to_gap = [i for i, ch in enumerate(base_allele.gapped_seq) if ch != "."]
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
        # gap-dependent metadata stays consistent.
        for k, b in enumerate(seq):
            gapped[ung_to_gap[k]] = b

        novel = _copy.deepcopy(base_allele)
        novel.name = name
        novel.gene = gene
        novel.ungapped_seq = new_ungapped
        novel.gapped_seq = "".join(gapped)
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

    def effective_dataconfig(self):
        """Return a copy of the source ``DataConfig`` with this genotype's
        novel alleles appended to their genes' allele lists — the reference
        the engine actually runs against when novel alleles are present."""
        import copy as _copy

        cfg = _copy.deepcopy(self._cfg)
        by_seg = {"V": cfg.v_alleles, "D": cfg.d_alleles, "J": cfg.j_alleles}
        for info in self._novel.values():
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
