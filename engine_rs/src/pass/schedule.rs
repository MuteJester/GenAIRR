//! `Schedule` — declarative dependency graph over passes.
//!
//! Replaces the previous `Vec<Box<dyn Pass>>` plan-as-list with a
//! proper dep graph in the style of Bevy's `Schedule` or Unity DOTS
//! `SystemGroup`. Edges are auto-derived from `Pass::requirements()`
//! and `Pass::effects()`. The compile step Kahn-topo-sorts; cycles
//! and missing producers are structured errors instead of "the
//! pipeline silently runs in a wrong order."
//!
//! ## What Stage 1 enforces
//!
//! - `AlleleAssignment(seg)` requirement → edge from the producer
//!   pass whose `AssignAllele(seg)` effect satisfies it.
//! - `TrimAllele(seg)` effect → edge to any later
//!   `AssembleSegment(seg)` effect. (Trim before assemble of the
//!   same segment.)
//!
//! Insertion order is the tiebreaker when the graph leaves passes
//! independent — Kahn's algorithm uses a min-heap keyed by insertion
//! index so the canonical V→Np1→D→Np2→J pipeline survives unchanged.
//!
//! ## What later stages add
//!
//! - Stage 2: effect hooks consume `PassEffect` post-execution
//!   (replacing the hand-coded `apply_live_call_updates` match).
//! - Stage 3: explicit `Schedule::add_edge(a, b)` for cross-cutting
//!   constraints not expressible via requirements/effects, plus
//!   named `SystemSet`s for grouping.

use std::collections::{BinaryHeap, HashMap};
use std::cmp::Reverse;

use crate::ir::Segment;

use super::{Pass, PassEffect, PassRequirement};

/// Stable identifier for a pass within a `Schedule`. The numeric
/// value is the insertion index — sorting node IDs ascending yields
/// insertion order, which is the topo-sort tiebreaker.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct NodeId(pub(crate) u32);

impl NodeId {
    #[inline]
    pub fn index(self) -> usize {
        self.0 as usize
    }
}

/// Failure mode for `Schedule::compile`.
#[derive(Debug, Clone)]
pub enum ScheduleError {
    /// Two or more passes form a cycle in the dep graph. The IDs in
    /// the vec form one such cycle (in traversal order).
    Cycle(Vec<NodeId>),
    /// A pass declared `RefData` requirement but no `refdata` was
    /// supplied to the compiler. Reported here so the dep-graph
    /// compile step is the single source of truth for plan-level
    /// errors; the older `analyze.rs` ad-hoc check is gone.
    MissingRefData { pass: NodeId, name: String },
    /// A pass requires an allele assignment for a segment but no
    /// upstream pass produces `AssignAllele(seg)`.
    MissingAlleleAssignment {
        pass: NodeId,
        name: String,
        segment: Segment,
    },
}

/// Edge between two pass nodes. `from` must execute before `to`.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub(crate) struct Edge {
    pub(crate) from: NodeId,
    pub(crate) to: NodeId,
}

/// A dependency-graph schedule of passes. Nodes are owned passes,
/// edges come from three sources:
///
/// 1. **Auto-derived edges** from each pass's `requirements()` /
///    `effects()` metadata. These cover the common producer/consumer
///    cases (e.g. trim before assemble of the same segment).
/// 2. **Explicit edges** added via [`Self::add_edge`] for
///    cross-cutting constraints that aren't naturally expressible
///    through metadata.
/// 3. **System-set edges** added via [`Self::add_to_set`] +
///    [`Self::add_set_edge`]. A set is a named bag of nodes; an edge
///    between sets `A → B` expands to N × M node-level edges at
///    compile time.
pub struct Schedule {
    nodes: Vec<Box<dyn Pass>>,
    /// Cached requirements per node — computed once at `push` time so
    /// `compile()` doesn't re-call the trait method per edge check.
    requirements: Vec<Vec<PassRequirement>>,
    /// Cached effects per node — same rationale.
    effects: Vec<Vec<PassEffect>>,
    /// User-declared explicit edges (Stage 3). Each entry forces
    /// `position(from) < position(to)` in the topo-sorted order.
    explicit_edges: Vec<Edge>,
    /// Named system sets. Each set maps a `&'static str` label to
    /// the list of nodes that belong to it.
    sets: HashMap<&'static str, Vec<NodeId>>,
    /// User-declared set-to-set edges. Expanded into node-level
    /// edges at compile time.
    set_edges: Vec<(&'static str, &'static str)>,
}

impl Schedule {
    /// An empty schedule.
    pub fn new() -> Self {
        Self {
            nodes: Vec::new(),
            requirements: Vec::new(),
            effects: Vec::new(),
            explicit_edges: Vec::new(),
            sets: HashMap::new(),
            set_edges: Vec::new(),
        }
    }

    /// Append a pass and return its stable `NodeId`. Insertion order
    /// drives the topo-sort tiebreaker, so for the canonical
    /// recombination pipeline `push` calls in the existing canonical
    /// order yield the existing canonical execution order.
    pub fn push(&mut self, pass: Box<dyn Pass>) -> NodeId {
        let id = NodeId(self.nodes.len() as u32);
        let reqs = pass.requirements();
        let effs = pass.effects();
        self.nodes.push(pass);
        self.requirements.push(reqs);
        self.effects.push(effs);
        id
    }

    /// Number of passes in the schedule.
    pub fn len(&self) -> usize {
        self.nodes.len()
    }

    /// Whether the schedule contains zero passes.
    pub fn is_empty(&self) -> bool {
        self.nodes.is_empty()
    }

    /// Read-only view of the underlying pass vector. Order is
    /// insertion order — for the topo-sorted execution order, call
    /// [`Self::compile`].
    pub fn passes(&self) -> &[Box<dyn Pass>] {
        &self.nodes
    }

    /// Declare an explicit ordering edge: `from` must execute before
    /// `to`. Use this for cross-cutting constraints that aren't
    /// expressible through `requirements()` / `effects()` — e.g.
    /// "this corruption pass must come after that mutation pass"
    /// when neither produces an effect the other consumes.
    ///
    /// Both `NodeId`s must point at nodes already pushed onto this
    /// schedule. Panics otherwise. Idempotent: redeclaring the same
    /// edge is a no-op (Kahn's tolerates duplicates).
    pub fn add_edge(&mut self, from: NodeId, to: NodeId) -> &mut Self {
        assert!(
            from.index() < self.nodes.len(),
            "Schedule::add_edge: `from` node {:?} not in schedule",
            from
        );
        assert!(
            to.index() < self.nodes.len(),
            "Schedule::add_edge: `to` node {:?} not in schedule",
            to
        );
        assert_ne!(
            from, to,
            "Schedule::add_edge: self-edge ({:?} -> {:?}) would create a cycle",
            from, to
        );
        self.explicit_edges.push(Edge { from, to });
        self
    }

    /// Shorthand for `add_edge(other, self_id)` — make `self_id`
    /// depend on `other`. Mirrors Bevy's `.after()`.
    pub fn after(&mut self, self_id: NodeId, other: NodeId) -> &mut Self {
        self.add_edge(other, self_id)
    }

    /// Shorthand for `add_edge(self_id, other)` — make `self_id`
    /// run before `other`. Mirrors Bevy's `.before()`.
    pub fn before(&mut self, self_id: NodeId, other: NodeId) -> &mut Self {
        self.add_edge(self_id, other)
    }

    /// Assign `node` to a named system set. Sets are bags of
    /// `NodeId`s identified by a `&'static str` label. Sets only
    /// affect the topo-sort when paired with [`Self::add_set_edge`].
    pub fn add_to_set(&mut self, node: NodeId, set: &'static str) -> &mut Self {
        assert!(
            node.index() < self.nodes.len(),
            "Schedule::add_to_set: node {:?} not in schedule",
            node
        );
        self.sets.entry(set).or_default().push(node);
        self
    }

    /// Every node in `from` must execute before every node in `to`.
    /// Expanded into individual node-level edges at compile time.
    /// Either set may be empty (the constraint is trivially
    /// satisfied) or absent from the schedule (treated as empty).
    pub fn add_set_edge(&mut self, from: &'static str, to: &'static str) -> &mut Self {
        self.set_edges.push((from, to));
        self
    }

    /// Cached requirements for a node, in declaration order.
    /// Used by Stage 2's effect-hook dispatcher.
    #[allow(dead_code)]
    pub(crate) fn requirements_of(&self, id: NodeId) -> &[PassRequirement] {
        &self.requirements[id.index()]
    }

    /// Cached effects for a node, in declaration order. Used by
    /// Stage 2's effect-hook dispatcher.
    #[allow(dead_code)]
    pub(crate) fn effects_of(&self, id: NodeId) -> &[PassEffect] {
        &self.effects[id.index()]
    }

    /// Compile the schedule into a topologically-sorted execution
    /// order. Each entry of the returned vec is a `NodeId`; the slice
    /// is a valid linear schedule (every edge `from → to` satisfies
    /// `position(from) < position(to)`).
    ///
    /// `has_refdata` lets the compile step verify `PassRequirement::RefData`
    /// without coupling the schedule to the refdata type itself.
    pub fn compile(&self, has_refdata: bool) -> Result<Vec<NodeId>, ScheduleError> {
        let edges = self.derive_edges();
        self.validate_hard_requirements(has_refdata)?;
        self.kahn_topo_sort(&edges)
    }

    /// Derive directed edges from the cached requirements/effects.
    ///
    /// Two edge sources:
    /// 1. **Producer/consumer**: each `AlleleAssignment(seg)`
    ///    requirement on a node creates an edge from every prior node
    ///    whose effects include `AssignAllele(seg)`.
    /// 2. **Implicit ordering**: any node with effect
    ///    `TrimAllele(seg)` gets an edge to every node with effect
    ///    `AssembleSegment(seg)`. This replaces the ad-hoc
    ///    `invalid_pass_order("trim_after_assembly.{seg}")` check.
    fn derive_edges(&self) -> Vec<Edge> {
        let mut edges = Vec::new();

        for to_idx in 0..self.nodes.len() {
            let to = NodeId(to_idx as u32);
            for req in &self.requirements[to_idx] {
                match req {
                    PassRequirement::AlleleAssignment(seg) => {
                        for from_idx in 0..self.nodes.len() {
                            if from_idx == to_idx {
                                continue;
                            }
                            if self.effects[from_idx]
                                .iter()
                                .any(|e| matches!(e, PassEffect::AssignAllele(s) if s == seg))
                            {
                                edges.push(Edge {
                                    from: NodeId(from_idx as u32),
                                    to,
                                });
                            }
                        }
                    }
                    PassRequirement::RefData => {
                        // RefData isn't produced by a pass; checked
                        // separately in `validate_hard_requirements`.
                    }
                }
            }
        }

        // Implicit: trim(seg) → assemble(seg).
        for from_idx in 0..self.nodes.len() {
            let trims_for: Vec<Segment> = self.effects[from_idx]
                .iter()
                .filter_map(|e| match e {
                    PassEffect::TrimAllele(s) => Some(*s),
                    _ => None,
                })
                .collect();
            if trims_for.is_empty() {
                continue;
            }
            for to_idx in 0..self.nodes.len() {
                if to_idx == from_idx {
                    continue;
                }
                for seg in &trims_for {
                    if self.effects[to_idx]
                        .iter()
                        .any(|e| matches!(e, PassEffect::AssembleSegment(s) if s == seg))
                    {
                        edges.push(Edge {
                            from: NodeId(from_idx as u32),
                            to: NodeId(to_idx as u32),
                        });
                    }
                }
            }
        }

        // Stage 3 — explicit user-declared edges.
        edges.extend(self.explicit_edges.iter().copied());

        // Stage 3 — expand set-level edges into node-level edges.
        // A → B with |A| = m, |B| = n produces m × n node edges.
        // Unknown set names expand to zero edges (treated as empty).
        for (from_set, to_set) in &self.set_edges {
            let empty: Vec<NodeId> = Vec::new();
            let from_members = self.sets.get(from_set).unwrap_or(&empty);
            let to_members = self.sets.get(to_set).unwrap_or(&empty);
            for &from in from_members {
                for &to in to_members {
                    if from != to {
                        edges.push(Edge { from, to });
                    }
                }
            }
        }

        edges
    }

    /// Check that every `RefData` / `AlleleAssignment` requirement
    /// is satisfied (either by the runtime or by an in-plan producer).
    fn validate_hard_requirements(&self, has_refdata: bool) -> Result<(), ScheduleError> {
        for idx in 0..self.nodes.len() {
            let node = NodeId(idx as u32);
            for req in &self.requirements[idx] {
                match req {
                    PassRequirement::RefData => {
                        if !has_refdata {
                            return Err(ScheduleError::MissingRefData {
                                pass: node,
                                name: self.nodes[idx].name().to_string(),
                            });
                        }
                    }
                    PassRequirement::AlleleAssignment(seg) => {
                        let produced = (0..self.nodes.len()).any(|i| {
                            i != idx
                                && self.effects[i]
                                    .iter()
                                    .any(|e| matches!(e, PassEffect::AssignAllele(s) if s == seg))
                        });
                        if !produced {
                            return Err(ScheduleError::MissingAlleleAssignment {
                                pass: node,
                                name: self.nodes[idx].name().to_string(),
                                segment: *seg,
                            });
                        }
                    }
                }
            }
        }
        Ok(())
    }

    /// Kahn's algorithm with a min-heap keyed by insertion index.
    /// The min-heap ensures that when the dep graph allows multiple
    /// passes to run in parallel (no edges between them), we pick the
    /// one inserted first — preserving the user's declared order as
    /// the tiebreaker.
    fn kahn_topo_sort(&self, edges: &[Edge]) -> Result<Vec<NodeId>, ScheduleError> {
        let n = self.nodes.len();
        let mut in_degree: Vec<u32> = vec![0; n];
        let mut adjacency: Vec<Vec<NodeId>> = vec![Vec::new(); n];
        for edge in edges {
            in_degree[edge.to.index()] += 1;
            adjacency[edge.from.index()].push(edge.to);
        }

        // Min-heap of ready nodes ordered by insertion index so the
        // canonical pipeline emerges first when ties exist.
        let mut ready: BinaryHeap<Reverse<NodeId>> = BinaryHeap::new();
        for idx in 0..n {
            if in_degree[idx] == 0 {
                ready.push(Reverse(NodeId(idx as u32)));
            }
        }

        let mut order = Vec::with_capacity(n);
        while let Some(Reverse(node)) = ready.pop() {
            order.push(node);
            for &next in &adjacency[node.index()] {
                let entry = &mut in_degree[next.index()];
                *entry -= 1;
                if *entry == 0 {
                    ready.push(Reverse(next));
                }
            }
        }

        if order.len() != n {
            // Cycle: nodes with nonzero remaining in-degree form one.
            let stuck: Vec<NodeId> = (0..n)
                .filter(|i| in_degree[*i] > 0)
                .map(|i| NodeId(i as u32))
                .collect();
            return Err(ScheduleError::Cycle(stuck));
        }

        Ok(order)
    }
}

impl Default for Schedule {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ir::Simulation;
    use crate::pass::{PassContext, PassEffect, PassRequirement};

    /// Test-only stub pass that records its own name + declared
    /// requirements/effects without doing any simulation work.
    struct StubPass {
        name: String,
        reqs: Vec<PassRequirement>,
        effs: Vec<PassEffect>,
    }

    impl StubPass {
        fn new(name: &str, reqs: Vec<PassRequirement>, effs: Vec<PassEffect>) -> Self {
            Self {
                name: name.to_string(),
                reqs,
                effs,
            }
        }
    }

    impl Pass for StubPass {
        fn name(&self) -> &str {
            &self.name
        }
        fn execute(&self, sim: &Simulation, _ctx: &mut PassContext) -> Simulation {
            sim.clone()
        }
        fn requirements(&self) -> Vec<PassRequirement> {
            self.reqs.clone()
        }
        fn effects(&self) -> Vec<PassEffect> {
            self.effs.clone()
        }
    }

    fn assemble_v_chain() -> Schedule {
        // Canonical fragment: sample_v → trim_v → assemble_v.
        let mut s = Schedule::new();
        s.push(Box::new(StubPass::new(
            "sample.v",
            vec![PassRequirement::RefData],
            vec![PassEffect::AssignAllele(Segment::V)],
        )));
        s.push(Box::new(StubPass::new(
            "trim.v",
            vec![PassRequirement::AlleleAssignment(Segment::V)],
            vec![PassEffect::TrimAllele(Segment::V)],
        )));
        s.push(Box::new(StubPass::new(
            "assemble.v",
            vec![
                PassRequirement::RefData,
                PassRequirement::AlleleAssignment(Segment::V),
            ],
            vec![PassEffect::AssembleSegment(Segment::V)],
        )));
        s
    }

    #[test]
    fn topo_sort_preserves_insertion_order_for_linear_chain() {
        let s = assemble_v_chain();
        let order = s.compile(true).unwrap();
        assert_eq!(
            order,
            vec![NodeId(0), NodeId(1), NodeId(2)],
            "linear dep chain must execute in insertion order"
        );
    }

    #[test]
    fn topo_sort_picks_insertion_order_for_independent_passes() {
        // Two independent allele samples: V then D. No edges between
        // them; insertion order should resolve the tie.
        let mut s = Schedule::new();
        s.push(Box::new(StubPass::new(
            "sample.v",
            vec![PassRequirement::RefData],
            vec![PassEffect::AssignAllele(Segment::V)],
        )));
        s.push(Box::new(StubPass::new(
            "sample.d",
            vec![PassRequirement::RefData],
            vec![PassEffect::AssignAllele(Segment::D)],
        )));
        let order = s.compile(true).unwrap();
        assert_eq!(order, vec![NodeId(0), NodeId(1)]);
    }

    #[test]
    fn missing_refdata_is_a_structured_error() {
        let s = assemble_v_chain();
        let err = s.compile(false).unwrap_err();
        assert!(
            matches!(err, ScheduleError::MissingRefData { .. }),
            "expected MissingRefData, got {:?}",
            err
        );
    }

    #[test]
    fn missing_allele_assignment_is_a_structured_error() {
        let mut s = Schedule::new();
        s.push(Box::new(StubPass::new(
            "trim.d",
            vec![PassRequirement::AlleleAssignment(Segment::D)],
            vec![PassEffect::TrimAllele(Segment::D)],
        )));
        let err = s.compile(true).unwrap_err();
        match err {
            ScheduleError::MissingAlleleAssignment { segment, .. } => {
                assert_eq!(segment, Segment::D);
            }
            other => panic!("expected MissingAlleleAssignment, got {:?}", other),
        }
    }

    #[test]
    fn trim_inserted_after_assemble_is_reordered_before_it() {
        // User pushes trim AFTER assemble — old engine rejected this
        // as `invalid_pass_order`. New engine reorders via the
        // implicit trim→assemble edge.
        let mut s = Schedule::new();
        s.push(Box::new(StubPass::new(
            "sample.v",
            vec![PassRequirement::RefData],
            vec![PassEffect::AssignAllele(Segment::V)],
        )));
        s.push(Box::new(StubPass::new(
            "assemble.v",
            vec![
                PassRequirement::RefData,
                PassRequirement::AlleleAssignment(Segment::V),
            ],
            vec![PassEffect::AssembleSegment(Segment::V)],
        )));
        s.push(Box::new(StubPass::new(
            "trim.v",
            vec![PassRequirement::AlleleAssignment(Segment::V)],
            vec![PassEffect::TrimAllele(Segment::V)],
        )));
        let order = s.compile(true).unwrap();
        // sample(0) → trim(2) → assemble(1): trim emerges before
        // assemble even though it was pushed last.
        assert_eq!(order, vec![NodeId(0), NodeId(2), NodeId(1)]);
    }

    #[test]
    fn explicit_edge_forces_ordering_between_otherwise_independent_passes() {
        // Two independent samples (no producer/consumer link). Push
        // V first, D second — without an explicit edge the topo-sort
        // honors insertion order. Add an explicit `D before V` edge
        // and the order flips.
        let mut s = Schedule::new();
        let v = s.push(Box::new(StubPass::new(
            "sample.v",
            vec![PassRequirement::RefData],
            vec![PassEffect::AssignAllele(Segment::V)],
        )));
        let d = s.push(Box::new(StubPass::new(
            "sample.d",
            vec![PassRequirement::RefData],
            vec![PassEffect::AssignAllele(Segment::D)],
        )));
        s.add_edge(d, v); // D before V
        let order = s.compile(true).unwrap();
        assert_eq!(order, vec![d, v], "explicit edge must override insertion-order tiebreak");
    }

    #[test]
    fn after_and_before_helpers_are_symmetric() {
        let mut s_after = Schedule::new();
        let a1 = s_after.push(Box::new(StubPass::new(
            "a",
            vec![PassRequirement::RefData],
            vec![PassEffect::AssignAllele(Segment::V)],
        )));
        let b1 = s_after.push(Box::new(StubPass::new(
            "b",
            vec![PassRequirement::RefData],
            vec![PassEffect::AssignAllele(Segment::D)],
        )));
        s_after.after(b1, a1); // b runs after a

        let mut s_before = Schedule::new();
        let a2 = s_before.push(Box::new(StubPass::new(
            "a",
            vec![PassRequirement::RefData],
            vec![PassEffect::AssignAllele(Segment::V)],
        )));
        let b2 = s_before.push(Box::new(StubPass::new(
            "b",
            vec![PassRequirement::RefData],
            vec![PassEffect::AssignAllele(Segment::D)],
        )));
        s_before.before(a2, b2); // a runs before b

        assert_eq!(
            s_after.compile(true).unwrap(),
            s_before.compile(true).unwrap(),
            "`.after()` and `.before()` declare the same edge"
        );
    }

    #[test]
    fn cycle_via_explicit_edges_is_detected() {
        let mut s = Schedule::new();
        let a = s.push(Box::new(StubPass::new(
            "a",
            vec![PassRequirement::RefData],
            vec![PassEffect::AssignAllele(Segment::V)],
        )));
        let b = s.push(Box::new(StubPass::new(
            "b",
            vec![PassRequirement::RefData],
            vec![PassEffect::AssignAllele(Segment::D)],
        )));
        s.add_edge(a, b);
        s.add_edge(b, a); // cycle.
        match s.compile(true).unwrap_err() {
            ScheduleError::Cycle(stuck) => {
                assert!(stuck.contains(&a) && stuck.contains(&b));
            }
            other => panic!("expected ScheduleError::Cycle, got {:?}", other),
        }
    }

    #[test]
    fn system_set_edge_orders_groups_of_passes() {
        // Two sets: "first", "second". The `first → second` edge
        // forces every member of "first" before every member of
        // "second", regardless of insertion order.
        let mut s = Schedule::new();
        // Push the "second" group first so insertion order doesn't
        // accidentally satisfy the constraint.
        let s1 = s.push(Box::new(StubPass::new(
            "second.a",
            vec![PassRequirement::RefData],
            vec![PassEffect::AssignAllele(Segment::V)],
        )));
        let s2 = s.push(Box::new(StubPass::new(
            "second.b",
            vec![PassRequirement::RefData],
            vec![PassEffect::AssignAllele(Segment::D)],
        )));
        let f1 = s.push(Box::new(StubPass::new(
            "first.a",
            vec![PassRequirement::RefData],
            vec![PassEffect::AssignAllele(Segment::J)],
        )));
        s.add_to_set(s1, "second")
            .add_to_set(s2, "second")
            .add_to_set(f1, "first")
            .add_set_edge("first", "second");
        let order = s.compile(true).unwrap();
        let pos = |id: NodeId| order.iter().position(|n| *n == id).unwrap();
        assert!(
            pos(f1) < pos(s1) && pos(f1) < pos(s2),
            "first.a must run before every member of `second`, got {:?}",
            order
        );
    }

    #[test]
    fn full_canonical_vdj_pipeline_topo_sorts_in_canonical_order() {
        // Mirrors the canonical Python push order:
        // sample(V,D,J) → trim(V,D,J) → assemble(V) → assemble(D) → assemble(J).
        let mut s = Schedule::new();
        for seg in [Segment::V, Segment::D, Segment::J] {
            s.push(Box::new(StubPass::new(
                &format!("sample.{:?}", seg),
                vec![PassRequirement::RefData],
                vec![PassEffect::AssignAllele(seg)],
            )));
        }
        for seg in [Segment::V, Segment::D, Segment::J] {
            s.push(Box::new(StubPass::new(
                &format!("trim.{:?}", seg),
                vec![PassRequirement::AlleleAssignment(seg)],
                vec![PassEffect::TrimAllele(seg)],
            )));
        }
        for seg in [Segment::V, Segment::D, Segment::J] {
            s.push(Box::new(StubPass::new(
                &format!("assemble.{:?}", seg),
                vec![
                    PassRequirement::RefData,
                    PassRequirement::AlleleAssignment(seg),
                ],
                vec![PassEffect::AssembleSegment(seg)],
            )));
        }
        let order = s.compile(true).unwrap();
        let names: Vec<String> = order
            .iter()
            .map(|id| s.nodes[id.index()].name().to_string())
            .collect();
        assert_eq!(
            names,
            vec![
                "sample.V",
                "sample.D",
                "sample.J",
                "trim.V",
                "trim.D",
                "trim.J",
                "assemble.V",
                "assemble.D",
                "assemble.J",
            ],
            "canonical canonical-order push must topo-sort to canonical order"
        );
    }
}
