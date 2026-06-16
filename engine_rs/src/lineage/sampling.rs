//! Sampling observed cells from a grown lineage and collapsing identical
//! genotypes into abundance-bearing observed nodes.

use std::collections::HashMap;

use crate::rng::Rng;

use super::tree::LineageTree;

/// Sample `n_sample` cells uniformly (with replacement) from `sampleable_ids`
/// (the LIVING population at sampling time — the final-generation cells, NOT all
/// tree leaves) and collapse identical genotypes: the first cell *drawn* carrying
/// a given genotype becomes the observed representative and accumulates the
/// abundance; later draws of the same genotype fold into it. Mutates `tree` in
/// place.
///
/// If `sampleable_ids` is empty (an extinct family — no living cells), nothing is
/// marked observed: an unobserved extinct clone yields zero records.
pub fn sample_and_collapse(
    tree: &mut LineageTree,
    sampleable_ids: &[u32],
    n_sample: u32,
    rng: &mut Rng,
) {
    if sampleable_ids.is_empty() || n_sample == 0 {
        return;
    }

    // genotype -> representative node id (first cell drawn with that genotype)
    let mut rep_by_genotype: HashMap<Vec<u8>, u32> = HashMap::new();
    // representative node id -> accumulated abundance
    let mut abundance: HashMap<u32, u32> = HashMap::new();

    for _ in 0..n_sample {
        // Unbiased uniform index (Lemire), matching the engine's RNG convention;
        // avoids the modulo bias of `next_u64() % len`.
        let idx = rng.range_u32(sampleable_ids.len() as u32) as usize;
        let cell_id = sampleable_ids[idx];
        let genotype = tree.get(cell_id).unwrap().genotype.clone();
        let rep = *rep_by_genotype.entry(genotype).or_insert(cell_id);
        *abundance.entry(rep).or_insert(0) += 1;
    }

    for (id, count) in abundance {
        // Write-back relies on the `node.id == index into nodes` invariant
        // (held by construction in the branching loop); guard it in debug builds.
        debug_assert_eq!(
            tree.nodes[id as usize].id, id,
            "id-index invariant violated during sample write-back"
        );
        let node = &mut tree.nodes[id as usize];
        node.abundance = count;
        node.observed = true;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lineage::tree::{LineageNode, LineageTree};
    use crate::rng::Rng;

    // Tree with duplicate genotypes among leaves to exercise collapse.
    fn tree_with_dupes() -> LineageTree {
        LineageTree {
            nodes: vec![
                LineageNode { id: 0, parent_id: None,    generation: 0, genotype: b"AAAA".to_vec(), mutations_from_parent: 0, affinity: 0.0, abundance: 0, observed: false },
                LineageNode { id: 1, parent_id: Some(0), generation: 1, genotype: b"AAAC".to_vec(), mutations_from_parent: 1, affinity: 0.0, abundance: 0, observed: false },
                LineageNode { id: 2, parent_id: Some(0), generation: 1, genotype: b"AAAC".to_vec(), mutations_from_parent: 1, affinity: 0.0, abundance: 0, observed: false },
                LineageNode { id: 3, parent_id: Some(1), generation: 2, genotype: b"AAGC".to_vec(), mutations_from_parent: 1, affinity: 0.0, abundance: 0, observed: false },
            ],
        }
    }

    // Sampleable set with duplicate genotypes among the living cells (nodes 1 & 2
    // share "AAAC") to exercise genotype-collapse.
    fn sampleable() -> Vec<u32> {
        vec![1, 2, 3]
    }

    #[test]
    fn sampling_sets_abundances_summing_to_n_sample() {
        let mut tree = tree_with_dupes();
        let mut rng = Rng::new(5);
        let n_sample = 3;
        sample_and_collapse(&mut tree, &sampleable(), n_sample, &mut rng);

        let total: u32 = tree.nodes.iter().map(|n| n.abundance).sum();
        assert_eq!(total, n_sample);

        assert!(tree.nodes.iter().any(|n| n.observed));
        for n in &tree.nodes {
            assert_eq!(n.observed, n.abundance > 0);
        }
    }

    #[test]
    fn identical_genotypes_collapse_into_one_observed_node() {
        let mut tree = tree_with_dupes();
        let mut rng = Rng::new(1);
        sample_and_collapse(&mut tree, &sampleable(), 3, &mut rng);
        use std::collections::HashSet;
        let mut seen: HashSet<Vec<u8>> = HashSet::new();
        for n in tree.nodes.iter().filter(|n| n.observed) {
            assert!(seen.insert(n.genotype.clone()),
                "genotype observed in more than one node (collapse failed)");
        }
    }

    #[test]
    fn empty_sampleable_set_marks_nothing_observed() {
        // An extinct family (no living cells) yields zero observed cells.
        let mut tree = tree_with_dupes();
        let mut rng = Rng::new(3);
        sample_and_collapse(&mut tree, &[], 5, &mut rng);
        assert!(tree.nodes.iter().all(|n| !n.observed && n.abundance == 0));
    }
}
