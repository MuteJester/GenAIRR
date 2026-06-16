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
        if allele not in names:
            raise ValueError(f"{allele!r} is not a known allele of {gene!r}")

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

    def complete_from_reference(self, policy: str = "homozygous_common") -> "Genotype":
        """Fill every UNspecified gene with a valid diploid state."""
        for seg in _SEGMENTS:
            for gene, allele_objs in _alleles_by_gene(self._cfg, seg).items():
                if gene in self._slots[seg] or not allele_objs:
                    continue
                names = [a.name for a in allele_objs]
                if policy == "homozygous_common":
                    self.homozygous(gene, names[0], segment=seg)
                elif policy == "heterozygous_first_two":
                    if len(names) >= 2:
                        self.heterozygous(gene, names[0], names[1], segment=seg)
                    else:
                        self.homozygous(gene, names[0], segment=seg)
                else:
                    raise ValueError(f"unknown policy {policy!r}")
        return self

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

    def to_table(self) -> List[Dict]:
        rows = []
        for seg in _SEGMENTS:
            for gene, haps in self._slots[seg].items():
                h0 = {a for (a, _, _) in haps[0]}
                h1 = {a for (a, _, _) in haps[1]}
                carried = sorted(h0 | h1)
                if not carried:
                    zyg = "deleted"
                elif h0 == h1 and len(h0) == 1:
                    zyg = "homozygous"
                else:
                    zyg = "heterozygous"
                rows.append(
                    {
                        "segment": seg,
                        "gene": gene,
                        "zygosity": zyg,
                        "haplotype_0": sorted(h0),
                        "haplotype_1": sorted(h1),
                        "permissive": self._permissive,
                    }
                )
        return rows

    def to_tsv(self, path: str) -> None:
        import csv

        rows = self.to_table()
        with open(path, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(
                ["segment", "gene", "zygosity", "haplotype_0", "haplotype_1", "permissive"]
            )
            for r in rows:
                w.writerow(
                    [
                        r["segment"],
                        r["gene"],
                        r["zygosity"],
                        ";".join(r["haplotype_0"]),
                        ";".join(r["haplotype_1"]),
                        r["permissive"],
                    ]
                )
