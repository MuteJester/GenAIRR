//! Per-individual diploid genotype model (PR1: known reference alleles).
use std::collections::HashMap;

use crate::ir::Segment;
use crate::refdata::{AlleleId, GeneId};

/// One carried allele in a haplotype gene slot. `copies` encodes
/// gene-copy multiplicity for the *same* allele; two different alleles
/// in a slot are two `GeneCopy` entries. `weight` is relative
/// within-slot expression.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GeneCopy {
    pub allele: AlleleId,
    pub copies: u8,
    pub weight: f32,
}

/// One chromosome's carried alleles, per V/D/J gene. An absent or empty
/// slot means the gene is deleted on this chromosome.
#[derive(Clone, Debug, Default)]
pub struct Haplotype {
    v: HashMap<GeneId, Vec<GeneCopy>>,
    d: HashMap<GeneId, Vec<GeneCopy>>,
    j: HashMap<GeneId, Vec<GeneCopy>>,
}

impl Haplotype {
    pub fn new() -> Self {
        Self::default()
    }

    fn map(&self, seg: Segment) -> &HashMap<GeneId, Vec<GeneCopy>> {
        match seg {
            Segment::V => &self.v,
            Segment::D => &self.d,
            Segment::J => &self.j,
            _ => panic!("Haplotype: segment must be V/D/J, got {seg:?}"),
        }
    }
    fn map_mut(&mut self, seg: Segment) -> &mut HashMap<GeneId, Vec<GeneCopy>> {
        match seg {
            Segment::V => &mut self.v,
            Segment::D => &mut self.d,
            Segment::J => &mut self.j,
            _ => panic!("Haplotype: segment must be V/D/J, got {seg:?}"),
        }
    }

    /// Set (replace) the copies carried for a gene on this chromosome.
    /// An empty `copies` vec means the gene is deleted here.
    pub fn set(&mut self, seg: Segment, gene: GeneId, copies: Vec<GeneCopy>) {
        self.map_mut(seg).insert(gene, copies);
    }
    pub fn slot(&self, seg: Segment, gene: GeneId) -> &[GeneCopy] {
        self.map(seg).get(&gene).map(Vec::as_slice).unwrap_or(&[])
    }
    pub fn is_deleted(&self, seg: Segment, gene: GeneId) -> bool {
        self.slot(seg, gene).is_empty()
    }
    /// Genes with at least one carried copy on this chromosome, in
    /// ascending GeneId order (deterministic).
    pub fn present_genes(&self, seg: Segment) -> impl Iterator<Item = GeneId> + '_ {
        let mut genes: Vec<GeneId> = self
            .map(seg)
            .iter()
            .filter(|(_, v)| !v.is_empty())
            .map(|(g, _)| *g)
            .collect();
        genes.sort_by_key(|g| g.index());
        genes.into_iter()
    }
    /// (GeneId, usage-weight) for each present gene, weight from `usage`.
    pub fn gene_weights<F: Fn(GeneId) -> f64>(&self, seg: Segment, usage: &F) -> Vec<(GeneId, f64)> {
        self.present_genes(seg).map(|g| (g, usage(g))).collect()
    }
    /// All carried allele ids for a segment across all present genes.
    pub fn carried_alleles(&self, seg: Segment) -> Vec<AlleleId> {
        let mut out = Vec::new();
        for g in self.present_genes(seg) {
            for c in self.slot(seg, g) {
                out.push(c.allele);
            }
        }
        out
    }
}

/// A diploid genotype: two chromosomes + draw weights + provenance.
#[derive(Clone, Debug)]
pub struct Genotype {
    haplotypes: [Haplotype; 2],
    chromosome_weights: [f32; 2],
    subject_id: Option<String>,
    source_refdata_hash: String,
}

impl Genotype {
    pub fn new(
        haplotypes: [Haplotype; 2],
        chromosome_weights: [f32; 2],
        subject_id: Option<String>,
        source_refdata_hash: String,
    ) -> Self {
        Self {
            haplotypes,
            chromosome_weights,
            subject_id,
            source_refdata_hash,
        }
    }
    pub fn haplotype(&self, c: usize) -> &Haplotype {
        &self.haplotypes[c]
    }
    pub fn chromosome_weights(&self) -> [f32; 2] {
        self.chromosome_weights
    }
    pub fn subject_id(&self) -> Option<&str> {
        self.subject_id.as_deref()
    }
    pub fn source_refdata_hash(&self) -> &str {
        &self.source_refdata_hash
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::refdata::AlleleId;

    fn copy(id: u32) -> GeneCopy {
        GeneCopy {
            allele: AlleleId::new(id),
            copies: 1,
            weight: 1.0,
        }
    }

    #[test]
    fn haplotype_reports_carried_alleles_per_gene_with_deletion_as_empty() {
        let mut h = Haplotype::new();
        h.set(Segment::V, GeneId::new(0), vec![copy(10)]); // carried
        h.set(Segment::V, GeneId::new(1), vec![]); // deleted
        assert_eq!(h.slot(Segment::V, GeneId::new(0)).len(), 1);
        assert!(h.is_deleted(Segment::V, GeneId::new(1)));
        assert!(h.is_deleted(Segment::V, GeneId::new(2))); // absent == deleted
        let genes: Vec<GeneId> = h.present_genes(Segment::V).collect();
        assert_eq!(genes, vec![GeneId::new(0)]); // only non-empty slots
    }

    #[test]
    fn genotype_carries_two_haplotypes_and_chromosome_weights() {
        let mut h0 = Haplotype::new();
        let mut h1 = Haplotype::new();
        h0.set(Segment::V, GeneId::new(0), vec![copy(10)]);
        h1.set(Segment::V, GeneId::new(0), vec![copy(11)]); // heterozygous
        let g = Genotype::new([h0, h1], [0.5, 0.5], Some("S1".into()), "sha256:x".into());
        assert_eq!(g.chromosome_weights(), [0.5, 0.5]);
        assert_eq!(g.subject_id(), Some("S1"));
        assert_eq!(
            g.haplotype(0).slot(Segment::V, GeneId::new(0))[0].allele,
            AlleleId::new(10)
        );
        assert_eq!(
            g.haplotype(1).slot(Segment::V, GeneId::new(0))[0].allele,
            AlleleId::new(11)
        );
    }

    #[test]
    fn gene_weights_restrict_to_present_genes_and_apply_usage() {
        // chromosome 0 carries genes 0 and 1; usage favors gene 1.
        let mut h = Haplotype::new();
        h.set(Segment::V, GeneId::new(0), vec![copy(10)]);
        h.set(Segment::V, GeneId::new(1), vec![copy(20)]);
        let usage = |g: GeneId| if g.index() == 1 { 3.0 } else { 1.0 };
        let w = h.gene_weights(Segment::V, &usage);
        assert_eq!(w, vec![(GeneId::new(0), 1.0), (GeneId::new(1), 3.0)]);
    }
}
