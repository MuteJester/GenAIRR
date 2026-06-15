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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lineage::tree::{LineageNode, LineageTree};

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
}
