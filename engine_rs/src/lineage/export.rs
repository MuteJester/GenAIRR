//! Ground-truth export for a `LineageTree`: node-table TSV, FASTA of all node
//! (ancestral + observed) sequences, and a Newick tree with per-edge mutation
//! counts as branch lengths. Consumed by lineage-inference benchmark tools.

use std::fmt::Write as _;

use super::tree::LineageTree;

/// Tab-separated node table, one row per node plus a header. `parent_id` is
/// `NA` for the root. Columns:
/// `node_id, parent_id, generation, mutations_from_parent, abundance, observed, sequence`.
pub fn to_node_table_tsv(tree: &LineageTree) -> String {
    let mut out = String::new();
    out.push_str(
        "node_id\tparent_id\tgeneration\tmutations_from_parent\tabundance\tobserved\tsequence\n",
    );
    for n in &tree.nodes {
        let parent = match n.parent_id {
            Some(p) => p.to_string(),
            None => "NA".to_string(),
        };
        let seq = String::from_utf8_lossy(&n.genotype);
        let _ = writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            n.id, parent, n.generation, n.mutations_from_parent, n.abundance, n.observed, seq
        );
    }
    out
}

/// FASTA of every node (ancestral + observed). Header carries node id,
/// generation, abundance, and observed flag so downstream tools can map a
/// record back to its ground-truth node.
pub fn to_fasta(tree: &LineageTree) -> String {
    let mut out = String::new();
    for n in &tree.nodes {
        let seq = String::from_utf8_lossy(&n.genotype);
        let _ = writeln!(
            out,
            ">node{}|gen={}|abundance={}|observed={}",
            n.id, n.generation, n.abundance, n.observed
        );
        let _ = writeln!(out, "{}", seq);
    }
    out
}

/// Recursively render the subtree rooted at `id` in Newick form, with the edge
/// to this node labelled by its `mutations_from_parent` count as branch length.
/// Recursion depth is bounded by the number of generations (small), not node
/// count, so this is safe from stack overflow at realistic family sizes.
fn newick_subtree(tree: &LineageTree, id: u32) -> String {
    let node = tree.get(id).expect("newick: node id out of range");
    let children = tree.children_of(id);
    let label = format!("node{id}");
    if children.is_empty() {
        format!("{label}:{}", node.mutations_from_parent)
    } else {
        let inner: Vec<String> = children
            .iter()
            .map(|c| newick_subtree(tree, c.id))
            .collect();
        format!("({}){label}:{}", inner.join(","), node.mutations_from_parent)
    }
}

/// Newick string for the whole tree. Branch lengths are per-edge mutation
/// counts. The founder is the named Newick root (its children are wrapped in
/// parentheses, its label follows) and carries no branch length, since it is
/// the origin. Always terminated with `;`. This is standard rooted Newick
/// (e.g. `((node3:1)node1:1,node2:1)node0;`) — no phantom outer node — so
/// ete3 / dendropy / Bio.Phylo parse `node0` as the actual root.
pub fn to_newick(tree: &LineageTree) -> String {
    let root = tree.root();
    let children = tree.children_of(root.id);
    let root_label = format!("node{}", root.id);
    let body = if children.is_empty() {
        root_label
    } else {
        let inner: Vec<String> = children
            .iter()
            .map(|c| newick_subtree(tree, c.id))
            .collect();
        format!("({}){root_label}", inner.join(","))
    };
    format!("{body};")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dist::{EmpiricalLengthDist, UniformBase};
    use crate::ir::{Nucleotide, NucHandle, Region, Segment, Simulation};
    use crate::lineage::{simulate_family, BranchingParams};
    use crate::lineage::tree::{LineageNode, LineageTree};
    use crate::passes::UniformMutationPass;

    fn grown_founder() -> Simulation {
        let mut sim = Simulation::new();
        for (i, b) in b"AAAAAAAA".iter().enumerate() {
            let (next, _) = sim.with_nucleotide_pushed(
                Nucleotide::germline(*b, i as u16, Segment::V));
            sim = next;
        }
        sim.with_region_added(Region::new(Segment::V, NucHandle::new(0), NucHandle::new(8)))
    }

    #[test]
    fn exports_a_grown_family_consistently() {
        let params = BranchingParams {
            lambda_base: 1.5, lambda_mut: 0.0, max_generations: 6,
            n_max: 300, n_sample: 20, seed: 2024,
        };
        let mutator = UniformMutationPass::new(
            Box::new(EmpiricalLengthDist::from_pairs(vec![(1, 1.0)])),
            Box::new(UniformBase),
        );
        let tree = simulate_family(&grown_founder(), &params, &mutator);
        assert!(tree.validate().is_ok());

        let tsv = to_node_table_tsv(&tree);
        let fasta = to_fasta(&tree);
        let nwk = to_newick(&tree);

        assert_eq!(tsv.lines().count(), tree.len() + 1);
        assert_eq!(fasta.lines().count(), tree.len() * 2);
        assert!(nwk.ends_with(';'));
        assert_eq!(nwk.matches('(').count(), nwk.matches(')').count());
        for n in &tree.nodes {
            assert!(nwk.contains(&format!("node{}", n.id)),
                "newick missing node{}", n.id);
        }
    }

    // root(0,"AAAA") -> 1("AAAC", abundance 2, observed), 2("AAAG"); 1 -> 3("ATAC", abundance 1, observed)
    fn sample_tree() -> LineageTree {
        LineageTree {
            nodes: vec![
                LineageNode { id: 0, parent_id: None,    generation: 0, genotype: b"AAAA".to_vec(), mutations_from_parent: 0, abundance: 0, observed: false },
                LineageNode { id: 1, parent_id: Some(0), generation: 1, genotype: b"AAAC".to_vec(), mutations_from_parent: 1, abundance: 2, observed: true },
                LineageNode { id: 2, parent_id: Some(0), generation: 1, genotype: b"AAAG".to_vec(), mutations_from_parent: 1, abundance: 0, observed: false },
                LineageNode { id: 3, parent_id: Some(1), generation: 2, genotype: b"ATAC".to_vec(), mutations_from_parent: 1, abundance: 1, observed: true },
            ],
        }
    }

    #[test]
    fn fasta_emits_every_node_with_metadata_header() {
        let fasta = to_fasta(&sample_tree());
        let lines: Vec<&str> = fasta.lines().collect();
        assert_eq!(lines.len(), 8); // 4 nodes => header+seq each
        assert_eq!(lines[0], ">node0|gen=0|abundance=0|observed=false");
        assert_eq!(lines[1], "AAAA");
        assert_eq!(lines[2], ">node1|gen=1|abundance=2|observed=true");
        assert_eq!(lines[3], "AAAC");
        assert!(fasta.contains(">node3|gen=2|abundance=1|observed=true"));
    }

    #[test]
    fn node_table_tsv_has_header_and_one_row_per_node() {
        let tsv = to_node_table_tsv(&sample_tree());
        let lines: Vec<&str> = tsv.lines().collect();
        assert_eq!(lines.len(), 5); // 1 header + 4 nodes
        assert_eq!(
            lines[0],
            "node_id\tparent_id\tgeneration\tmutations_from_parent\tabundance\tobserved\tsequence"
        );
        assert_eq!(lines[1], "0\tNA\t0\t0\t0\tfalse\tAAAA");
        assert_eq!(lines[2], "1\t0\t1\t1\t2\ttrue\tAAAC");
        assert_eq!(lines[4], "3\t1\t2\t1\t1\ttrue\tATAC");
    }

    #[test]
    fn newick_encodes_topology_and_branch_lengths() {
        let nwk = to_newick(&sample_tree());
        assert!(nwk.ends_with(';'), "newick must end with ';': {nwk}");
        let opens = nwk.matches('(').count();
        let closes = nwk.matches(')').count();
        assert_eq!(opens, closes, "unbalanced parens: {nwk}");
        assert!(nwk.contains("node1:1"));
        assert!(nwk.contains("node2:1"));
        assert!(nwk.contains("node3:1"));
        assert!(nwk.contains("node0"));
        assert!(!nwk.contains("node0:"), "root must not have a branch length: {nwk}");
        assert_eq!(nwk, "((node3:1)node1:1,node2:1)node0;");
    }

    #[test]
    fn newick_single_node_tree() {
        let tree = LineageTree {
            nodes: vec![LineageNode {
                id: 0, parent_id: None, generation: 0, genotype: b"AAAA".to_vec(),
                mutations_from_parent: 0, abundance: 1, observed: true,
            }],
        };
        assert_eq!(to_newick(&tree), "node0;");
    }
}
