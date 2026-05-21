use std::fmt;

use crate::contract::ContractViolation;
use crate::dist::FilteredSampleError;
use crate::ir::Segment;

/// Structured failure from a fallible pass execution.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum PassError {
    /// A contract-aware sampling pass could not produce an admissible candidate.
    ConstraintSampling {
        pass_name: String,
        address: String,
        reason: FilteredSampleError,
    },
    /// A pass that requires reference data was run without it.
    MissingRefData { pass_name: String },
    /// A pass requires an allele assignment that is absent.
    MissingAssignment { pass_name: String, segment: Segment },
    /// The IR references an allele id absent from the supplied reference data.
    MissingAllele {
        pass_name: String,
        segment: Segment,
        allele_id: u32,
    },
    /// A distribution emitted a value outside the pass's valid domain.
    InvalidDistributionOutput {
        pass_name: String,
        address: String,
        value: i64,
        reason: String,
    },
    /// The current IR and plan state are inconsistent for this pass.
    InvalidPlanState { pass_name: String, reason: String },
    /// A pass produced a post-state that violates active contracts.
    ContractViolation {
        pass_name: String,
        violations: Vec<ContractViolation>,
    },
}

impl PassError {
    pub fn constraint_sampling(
        pass_name: impl Into<String>,
        address: impl Into<String>,
        reason: FilteredSampleError,
    ) -> Self {
        Self::ConstraintSampling {
            pass_name: pass_name.into(),
            address: address.into(),
            reason,
        }
    }

    pub fn missing_refdata(pass_name: impl Into<String>) -> Self {
        Self::MissingRefData {
            pass_name: pass_name.into(),
        }
    }

    pub fn missing_assignment(pass_name: impl Into<String>, segment: Segment) -> Self {
        Self::MissingAssignment {
            pass_name: pass_name.into(),
            segment,
        }
    }

    pub fn missing_allele(pass_name: impl Into<String>, segment: Segment, allele_id: u32) -> Self {
        Self::MissingAllele {
            pass_name: pass_name.into(),
            segment,
            allele_id,
        }
    }

    pub fn invalid_distribution_output(
        pass_name: impl Into<String>,
        address: impl Into<String>,
        value: i64,
        reason: impl Into<String>,
    ) -> Self {
        Self::InvalidDistributionOutput {
            pass_name: pass_name.into(),
            address: address.into(),
            value,
            reason: reason.into(),
        }
    }

    pub fn invalid_plan_state(pass_name: impl Into<String>, reason: impl Into<String>) -> Self {
        Self::InvalidPlanState {
            pass_name: pass_name.into(),
            reason: reason.into(),
        }
    }

    pub fn contract_violation(
        pass_name: impl Into<String>,
        violations: Vec<ContractViolation>,
    ) -> Self {
        Self::ContractViolation {
            pass_name: pass_name.into(),
            violations,
        }
    }

    pub fn pass_name(&self) -> &str {
        match self {
            Self::ConstraintSampling { pass_name, .. } => pass_name,
            Self::MissingRefData { pass_name } => pass_name,
            Self::MissingAssignment { pass_name, .. } => pass_name,
            Self::MissingAllele { pass_name, .. } => pass_name,
            Self::InvalidDistributionOutput { pass_name, .. } => pass_name,
            Self::InvalidPlanState { pass_name, .. } => pass_name,
            Self::ContractViolation { pass_name, .. } => pass_name,
        }
    }

    pub fn address(&self) -> &str {
        match self {
            Self::ConstraintSampling { address, .. } => address,
            Self::InvalidDistributionOutput { address, .. } => address,
            Self::MissingRefData { .. }
            | Self::MissingAssignment { .. }
            | Self::MissingAllele { .. }
            | Self::InvalidPlanState { .. }
            | Self::ContractViolation { .. } => "",
        }
    }

    pub fn constraint_reason(&self) -> Option<FilteredSampleError> {
        match self {
            Self::ConstraintSampling { reason, .. } => Some(*reason),
            _ => None,
        }
    }
}

impl fmt::Display for PassError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ConstraintSampling {
                pass_name,
                address,
                reason,
            } => {
                let detail = match reason {
                    FilteredSampleError::SupportUnavailable => {
                        "distribution support is not enumerable"
                    }
                    FilteredSampleError::EmptyAdmissibleSupport => {
                        "no admissible candidates remain after contract filtering"
                    }
                    FilteredSampleError::InvalidFilteredSupport => {
                        "filtered support has invalid weights"
                    }
                };
                write!(
                    f,
                    "{} failed strict constrained sampling at {}: {}",
                    pass_name, address, detail
                )
            }
            Self::MissingRefData { pass_name } => {
                write!(
                    f,
                    "{} requires reference data in strict execution",
                    pass_name
                )
            }
            Self::MissingAssignment { pass_name, segment } => write!(
                f,
                "{} requires a {:?} allele assignment before execution",
                pass_name, segment
            ),
            Self::MissingAllele {
                pass_name,
                segment,
                allele_id,
            } => write!(
                f,
                "{} could not resolve {:?} allele id {} in reference data",
                pass_name, segment, allele_id
            ),
            Self::InvalidDistributionOutput {
                pass_name,
                address,
                value,
                reason,
            } => write!(
                f,
                "{} received invalid distribution output at {}: {} ({})",
                pass_name, address, value, reason
            ),
            Self::InvalidPlanState { pass_name, reason } => {
                write!(
                    f,
                    "{} cannot execute in current plan state: {}",
                    pass_name, reason
                )
            }
            Self::ContractViolation {
                pass_name,
                violations,
            } => {
                write!(f, "{} produced a contract-invalid state", pass_name)?;
                for violation in violations {
                    write!(f, "; {}: {}", violation.contract_name, violation.reason)?;
                }
                Ok(())
            }
        }
    }
}

impl std::error::Error for PassError {}
