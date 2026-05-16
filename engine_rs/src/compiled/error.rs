//! Fatal compile-time error types + builder helpers.
//!
//! `analyze_plan` accumulates `CompileError`s into a `CompileErrors`
//! batch and returns them to the caller. The builder fns are
//! `pub(super)` so neighbour submodules (analyze, execute,
//! feasibility_builder) can emit consistently-shaped errors.

use crate::ir::Segment;
use std::fmt;

/// Fatal compile-time error.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct CompileError {
    pub pass_index: Option<usize>,
    pub pass_name: Option<String>,
    pub kind: CompileErrorKind,
}

/// Stable class of compile-time plan error.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum CompileErrorKind {
    MissingRefData,
    MissingAssignment {
        segment: Segment,
    },
    InvalidPassOrder {
        reason: String,
    },
    InvalidParameterSupport {
        address: String,
        reason: String,
    },
    ContractPrecondition {
        contract_name: String,
        reason: String,
    },
}

/// All fatal compile-time errors found in one plan scan.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct CompileErrors {
    pub errors: Vec<CompileError>,
}

impl CompileErrors {
    pub fn is_empty(&self) -> bool {
        self.errors.is_empty()
    }
}

impl fmt::Display for CompileErrors {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.errors.is_empty() {
            return write!(f, "plan compilation failed with no errors");
        }

        write!(f, "plan compilation failed")?;
        for error in &self.errors {
            let location = match (error.pass_index, error.pass_name.as_deref()) {
                (Some(index), Some(name)) => format!("pass {index} ({name})"),
                (Some(index), None) => format!("pass {index}"),
                (None, Some(name)) => name.to_string(),
                (None, None) => "plan".to_string(),
            };
            match &error.kind {
                CompileErrorKind::MissingRefData => {
                    write!(f, "; {location} requires reference data")?
                }
                CompileErrorKind::MissingAssignment { segment } => write!(
                    f,
                    "; {location} requires a prior {:?} allele assignment",
                    segment
                )?,
                CompileErrorKind::InvalidPassOrder { reason } => {
                    write!(f, "; {location} has invalid pass order: {reason}")?
                }
                CompileErrorKind::InvalidParameterSupport { address, reason } => write!(
                    f,
                    "; {location} has invalid parameter support at {address}: {reason}",
                )?,
                CompileErrorKind::ContractPrecondition {
                    contract_name,
                    reason,
                } => write!(
                    f,
                    "; contract {contract_name} precondition failed: {reason}",
                )?,
            }
        }
        Ok(())
    }
}

impl std::error::Error for CompileErrors {}

// ── builder helpers ────────────────────────────────────────────────

pub(super) fn pass_scoped_error(index: usize, name: &str, kind: CompileErrorKind) -> CompileError {
    CompileError {
        pass_index: Some(index),
        pass_name: Some(name.to_string()),
        kind,
    }
}

pub(super) fn plan_scoped_error(kind: CompileErrorKind) -> CompileError {
    CompileError {
        pass_index: None,
        pass_name: None,
        kind,
    }
}

pub(super) fn invalid_parameter(
    address: impl Into<String>,
    reason: impl Into<String>,
) -> CompileErrorKind {
    CompileErrorKind::InvalidParameterSupport {
        address: address.into(),
        reason: reason.into(),
    }
}

pub(super) fn invalid_pass_order(reason: impl Into<String>) -> CompileErrorKind {
    CompileErrorKind::InvalidPassOrder {
        reason: reason.into(),
    }
}

pub(super) fn contract_precondition(
    contract_name: impl Into<String>,
    reason: impl Into<String>,
) -> CompileErrorKind {
    CompileErrorKind::ContractPrecondition {
        contract_name: contract_name.into(),
        reason: reason.into(),
    }
}
