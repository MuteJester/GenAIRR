//! Immutable donor allele support. IDs retain their reference-pool indices.

use super::AlleleBitSet;
use crate::genotype::Genotype;
use crate::ir::Segment;
use crate::refdata::AlleleId;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct LiveCallSupport {
    v: Vec<AlleleId>,
    d: Vec<AlleleId>,
    j: Vec<AlleleId>,
}

impl LiveCallSupport {
    /// Inference can use either donor chromosome, irrespective of which
    /// chromosome or expression weight generated this rearrangement.
    pub fn from_genotype(genotype: &Genotype) -> Self {
        let union = |segment| {
            let mut ids = genotype.haplotype(0).carried_alleles(segment);
            ids.extend(genotype.haplotype(1).carried_alleles(segment));
            ids.sort_unstable_by_key(|id| id.index());
            ids.dedup();
            ids
        };
        Self {
            v: union(Segment::V),
            d: union(Segment::D),
            j: union(Segment::J),
        }
    }

    pub fn alleles(&self, segment: Segment) -> Option<&[AlleleId]> {
        match segment {
            Segment::V => Some(&self.v),
            Segment::D => Some(&self.d),
            Segment::J => Some(&self.j),
            _ => None,
        }
    }
}

pub(super) fn supported(id: AlleleId, allowed: Option<&[AlleleId]>) -> bool {
    allowed.map_or(true, |ids| {
        ids.binary_search_by_key(&id.index(), |a| a.index()).is_ok()
    })
}

pub(super) fn max_supported_score(scores: &[u32], allowed: Option<&[AlleleId]>) -> u32 {
    match allowed {
        Some(ids) => ids
            .iter()
            .filter_map(|id| scores.get(id.as_usize()))
            .copied()
            .max()
            .unwrap_or(0),
        None => scores.iter().copied().max().unwrap_or(0),
    }
}

pub(super) fn best_supported_alleles(scores: &[u32], allowed: Option<&[AlleleId]>) -> AlleleBitSet {
    let maximum = max_supported_score(scores, allowed);
    AlleleBitSet::from_ids(
        scores.len(),
        scores.iter().enumerate().filter_map(|(index, score)| {
            let id = AlleleId::new(index as u32);
            (*score == maximum && supported(id, allowed)).then_some(id)
        }),
    )
}
