use crate::assignment::TrimEnd;
use crate::ir::Segment;
use crate::refdata::AlleleId;

/// Compile-time view of an integer distribution's support.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum IntegerSupport {
    Enumerated(Vec<i64>),
    Unavailable,
    Invalid(String),
}

impl IntegerSupport {
    pub(crate) fn from_weighted_pairs(pairs: Option<Vec<(i64, f64)>>) -> Self {
        let Some(pairs) = pairs else {
            return Self::Unavailable;
        };
        if pairs.is_empty() {
            return Self::Invalid("empty_support".to_string());
        }

        let mut values = Vec::with_capacity(pairs.len());
        for (index, (value, weight)) in pairs.into_iter().enumerate() {
            if weight <= 0.0 || !weight.is_finite() {
                return Self::Invalid(format!("invalid_weight_at_index_{index}"));
            }
            values.push(value);
        }
        Self::Enumerated(values)
    }
}

/// Compile-time view of an allele-id distribution's support.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum AlleleIdSupport {
    Enumerated(Vec<AlleleId>),
    Unavailable,
    Invalid(String),
}

impl AlleleIdSupport {
    pub(crate) fn from_weighted_pairs(pairs: Option<Vec<(AlleleId, f64)>>) -> Self {
        let Some(pairs) = pairs else {
            return Self::Unavailable;
        };
        if pairs.is_empty() {
            return Self::Invalid("empty_support".to_string());
        }

        let mut values = Vec::with_capacity(pairs.len());
        for (index, (value, weight)) in pairs.into_iter().enumerate() {
            if weight <= 0.0 || !weight.is_finite() {
                return Self::Invalid(format!("invalid_weight_at_index_{index}"));
            }
            values.push(value);
        }
        Self::Enumerated(values)
    }
}

/// Positive-weight support of a canonical first-order N-base generator.
/// Bits are A/C/G/T; transitions are indexed by the preceding N base.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub struct NpBaseSupport {
    pub first: u8,
    pub transitions: [u8; 4],
}

/// Typed facts exposed by passes for semantic compile-time analysis.
#[derive(Clone, Debug, PartialEq)]
pub enum PassCompileFact {
    AlleleSampleSupport {
        segment: Segment,
        support: AlleleIdSupport,
    },
    TrimSupport {
        segment: Segment,
        end: TrimEnd,
        support: IntegerSupport,
    },
    NpLengthSupport {
        segment: Segment,
        support: IntegerSupport,
    },
    NpBaseSupport {
        segment: Segment,
        support: Option<NpBaseSupport>,
    },
    DOrientationSupport {
        orientations: Vec<bool>,
    },
    /// Per-end P-nucleotide length integer support, declared
    /// by `PAdditionPass`. Symmetric with `NpLengthSupport` for
    /// downstream consumers (schedule analyser, feasibility
    /// pre-flight, plan-signature folding) that want to
    /// inspect the per-end length vocabulary.
    PLengthSupport {
        end: crate::address::PEnd,
        support: IntegerSupport,
    },
    /// Joint, canonical signed outcomes at the original coding ends.
    CodingEndSupport {
        segment: Segment,
        outcomes: Vec<crate::end_processing::EndOutcome>,
    },
    /// Deterministic emission from a committed original-end outcome.
    CodingEndEmission {
        end: crate::address::PEnd,
    },
}
