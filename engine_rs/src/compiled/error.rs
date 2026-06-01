//! Fatal compile-time error types + builder helpers.
//!
//! `analyze_plan` accumulates `CompileError`s into a `CompileErrors`
//! batch and returns them to the caller. The builder fns are
//! `pub(super)` so neighbour submodules (analyze, execute,
//! feasibility_builder) can emit consistently-shaped errors.

use crate::ir::Segment;
use crate::pass::ScheduleError;
use crate::refdata::RefDataValidationErrors;
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
    /// Reference data failed structural validation (empty required
    /// pools, malformed bytes, anchor mismatches, etc.). Carries the
    /// full [`RefDataValidationErrors`] so consumers see every
    /// problem in one diagnostic — both Fatal (always-blocking) and
    /// Curatable (pseudogene-shaped) issues — plus the mode the gate
    /// was run under so the remediation message can suggest
    /// `allow_curatable_refdata` when only curatable issues remain.
    /// Emitted at the very start of compile, before schedule or
    /// precondition checks; bad refdata is a root input problem, not
    /// a pass-order problem.
    RefDataValidation {
        errors: RefDataValidationErrors,
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

    /// Wrap a refdata-validation error aggregate as a single
    /// plan-scoped compile error. The aggregated form is intentional:
    /// users get every malformed-refdata issue in one diagnostic —
    /// fatal and curatable alike — instead of having to fix-and-
    /// recompile per issue.
    pub fn from_refdata_validation(errors: RefDataValidationErrors) -> Self {
        Self {
            errors: vec![plan_scoped_error(CompileErrorKind::RefDataValidation { errors })],
        }
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
                CompileErrorKind::RefDataValidation { errors } => {
                    // Delegate to `RefDataValidationErrors::Display`
                    // so the per-issue severity tags and the
                    // pseudogene/ORF remediation hint stay in one
                    // place. Indented under the per-error semicolon
                    // for readability.
                    write!(f, "; {errors}")?;
                }
            }
        }
        Ok(())
    }
}

impl std::error::Error for CompileErrors {}

// ── builder helpers ────────────────────────────────────────────────

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

pub(super) fn contract_precondition(
    contract_name: impl Into<String>,
    reason: impl Into<String>,
) -> CompileErrorKind {
    CompileErrorKind::ContractPrecondition {
        contract_name: contract_name.into(),
        reason: reason.into(),
    }
}

/// Convert a `ScheduleError` from the topo-sort step into the
/// corresponding `CompileError`, preserving pass index/name so the
/// human-readable message points at the offending pass.
pub(super) fn schedule_error_into_compile_error(
    err: ScheduleError,
    name_for: impl Fn(usize) -> String,
) -> CompileError {
    match err {
        ScheduleError::MissingRefData { pass, name } => CompileError {
            pass_index: Some(pass.index()),
            pass_name: Some(name),
            kind: CompileErrorKind::MissingRefData,
        },
        ScheduleError::MissingAlleleAssignment {
            pass,
            name,
            segment,
        } => CompileError {
            pass_index: Some(pass.index()),
            pass_name: Some(name),
            kind: CompileErrorKind::MissingAssignment { segment },
        },
        ScheduleError::Cycle(stuck) => {
            let names: Vec<String> = stuck.iter().map(|n| name_for(n.index())).collect();
            let reason = format!("cycle:{}", names.join("->"));
            CompileError {
                pass_index: stuck.first().map(|n| n.index()),
                pass_name: stuck.first().map(|n| name_for(n.index())),
                kind: CompileErrorKind::InvalidPassOrder { reason },
            }
        }
    }
}
