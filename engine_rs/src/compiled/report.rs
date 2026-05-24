//! Compile-time analysis report types.
//!
//! Output of `analyze_plan` — pass-level summaries, declared stochastic
//! choices, active contracts, and non-fatal warnings.

use super::ExecutionPolicy;
use crate::pass::{NodeId, PassEffect, PassRequirement};

/// One pass-level compile summary.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PassSummary {
    pub index: usize,
    pub name: String,
    pub declared_choices: Vec<String>,
    pub requirements: Vec<PassRequirement>,
    pub effects: Vec<PassEffect>,
}

/// One declared stochastic choice, annotated with the pass that owns it.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct DeclaredChoice {
    pub pass_index: usize,
    pub pass_name: String,
    pub address: String,
}

/// Compile-time analysis output for a plan.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct CompileReport {
    pub policy: ExecutionPolicy,
    pub pass_summaries: Vec<PassSummary>,
    pub declared_choices: Vec<DeclaredChoice>,
    pub active_contracts: Vec<String>,
    pub warnings: Vec<CompileWarning>,
}

impl CompileReport {
    pub fn pass_names(&self) -> Vec<String> {
        self.pass_summaries
            .iter()
            .map(|summary| summary.name.clone())
            .collect()
    }

    /// Reorder `pass_summaries` from insertion order into the
    /// topo-sorted execution order produced by `Schedule::compile`,
    /// and rewrite each `PassSummary.index` to its post-sort
    /// position so `report.pass_names()` matches what actually runs.
    /// `declared_choices`' `pass_index` is updated to match.
    pub(super) fn reorder_to_execution(&mut self, sorted_order: &[NodeId]) {
        // Build a map: original insertion index → execution position.
        let mut exec_pos = vec![0usize; self.pass_summaries.len()];
        for (exec_idx, node) in sorted_order.iter().enumerate() {
            exec_pos[node.index()] = exec_idx;
        }

        // Reorder pass_summaries into execution order, rewriting
        // each `index` to its execution position.
        let mut reordered: Vec<PassSummary> = Vec::with_capacity(self.pass_summaries.len());
        for node in sorted_order {
            let mut summary = self.pass_summaries[node.index()].clone();
            summary.index = reordered.len();
            reordered.push(summary);
        }
        self.pass_summaries = reordered;

        // Rewrite `pass_index` on declared_choices to match new
        // execution positions.
        for choice in &mut self.declared_choices {
            choice.pass_index = exec_pos[choice.pass_index];
        }
    }
}

/// Non-fatal compile-time diagnostic.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct CompileWarning {
    pub pass_index: Option<usize>,
    pub pass_name: Option<String>,
    pub message: String,
}
