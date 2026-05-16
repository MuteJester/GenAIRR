//! Compile-time analysis report types.
//!
//! Output of `analyze_plan` — pass-level summaries, declared stochastic
//! choices, active contracts, and non-fatal warnings.

use super::ExecutionPolicy;
use crate::pass::{PassEffect, PassRequirement};

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
}

/// Non-fatal compile-time diagnostic.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct CompileWarning {
    pub pass_index: Option<usize>,
    pub pass_name: Option<String>,
    pub message: String,
}
