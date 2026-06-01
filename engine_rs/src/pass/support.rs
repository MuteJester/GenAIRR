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

/// Typed facts exposed by passes for semantic compile-time analysis.
#[derive(Clone, Debug, Eq, PartialEq)]
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
    /// Per-end P-nucleotide length integer support, declared
    /// by `PAdditionPass`. Symmetric with `NpLengthSupport` for
    /// downstream consumers (schedule analyser, feasibility
    /// pre-flight, plan-signature folding) that want to
    /// inspect the per-end length vocabulary.
    PLengthSupport {
        end: crate::address::PEnd,
        support: IntegerSupport,
    },
}
