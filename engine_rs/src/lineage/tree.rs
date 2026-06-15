//! Clonal lineage tree data structures.
//!
//! A `LineageTree` is a flat arena of `LineageNode`s addressed by `id`. The
//! root is the founder (naive rearrangement); every other node is a somatic
//! descendant produced by one cell division. `genotype` is the node's
//! nucleotide sequence (pool bases) captured at creation, used for
//! genotype-collapse and downstream export.

/// One cell in a clonal lineage.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LineageNode {
    /// Stable index of this node within `LineageTree::nodes` (also its arena position).
    pub id: u32,
    /// Parent node id; `None` only for the root founder.
    pub parent_id: Option<u32>,
    /// Generation depth from the root (root == 0).
    pub generation: u32,
    /// Nucleotide bases of this cell's `Simulation` pool at creation time.
    pub genotype: Vec<u8>,
    /// Mutations introduced on the edge from the parent to this node
    /// (the number of substitutions on the parent→child edge; 0 for the root).
    pub mutations_from_parent: u32,
    /// Observation count after sampling + genotype-collapse. 0 until sampled.
    pub abundance: u32,
    /// Whether this node was observed: the representative node a sampled
    /// genotype collapses onto. (Internal-ancestor observation is a later concern.)
    pub observed: bool,
}

/// A clonal lineage as a flat arena of nodes (node `id` == index into `nodes`).
#[derive(Clone, Debug, Default)]
pub struct LineageTree {
    pub nodes: Vec<LineageNode>,
}

impl LineageTree {
    /// Total node count (founder + all descendants).
    pub fn len(&self) -> usize {
        self.nodes.len()
    }

    /// True when the tree has no nodes.
    pub fn is_empty(&self) -> bool {
        self.nodes.is_empty()
    }

    /// The founder node (the unique node with no parent).
    ///
    /// Panics if the tree has no root — call only on a grown tree.
    pub fn root(&self) -> &LineageNode {
        self.nodes
            .iter()
            .find(|n| n.parent_id.is_none())
            .expect("LineageTree has no root node")
    }

    /// Node by id, or `None` if out of range.
    pub fn get(&self, id: u32) -> Option<&LineageNode> {
        self.nodes.get(id as usize)
    }

    /// Direct children of `id`, in ascending id order.
    pub fn children_of(&self, id: u32) -> Vec<&LineageNode> {
        let mut v: Vec<&LineageNode> = self.nodes
            .iter()
            .filter(|n| n.parent_id == Some(id))
            .collect();
        v.sort_unstable_by_key(|n| n.id);
        v
    }

    /// Nodes with no children (tips), in ascending id order.
    pub fn leaves(&self) -> Vec<&LineageNode> {
        let mut v: Vec<&LineageNode> = self.nodes
            .iter()
            .filter(|n| self.children_of(n.id).is_empty())
            .collect();
        v.sort_unstable_by_key(|n| n.id);
        v
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Build a tiny tree by hand: root(0) -> {1, 2}; 1 -> {3}
    fn hand_tree() -> LineageTree {
        LineageTree {
            nodes: vec![
                LineageNode { id: 0, parent_id: None,    generation: 0, genotype: b"AAAA".to_vec(), mutations_from_parent: 0, abundance: 0, observed: false },
                LineageNode { id: 1, parent_id: Some(0), generation: 1, genotype: b"AAAC".to_vec(), mutations_from_parent: 1, abundance: 0, observed: false },
                LineageNode { id: 2, parent_id: Some(0), generation: 1, genotype: b"AAAG".to_vec(), mutations_from_parent: 1, abundance: 0, observed: false },
                LineageNode { id: 3, parent_id: Some(1), generation: 2, genotype: b"AATC".to_vec(), mutations_from_parent: 1, abundance: 0, observed: false },
            ],
        }
    }

    #[test]
    fn tree_accessors_report_structure() {
        let t = hand_tree();
        assert_eq!(t.root().id, 0);
        assert_eq!(t.len(), 4);

        let root_children: Vec<u32> = t.children_of(0).iter().map(|n| n.id).collect();
        assert_eq!(root_children, vec![1, 2]);

        let leaves: Vec<u32> = t.leaves().iter().map(|n| n.id).collect();
        assert_eq!(leaves, vec![2, 3]); // nodes with no children

        assert_eq!(t.get(3).unwrap().parent_id, Some(1));
        assert!(t.get(99).is_none());
    }
}
