//! AIRR Rearrangement record builder.
//!
//! Builds a fully-populated AIRR-format record from an `Outcome` and
//! its `RefDataConfig` in one walk. This is the Rust replacement for
//! the Python `_airr_record.outcome_to_airr_record` builder; the
//! Python side is now a thin wrapper around the PyO3 method this
//! module backs.
//!
//! Field semantics match exactly what the Python code produced.
//! Convention is **0-based half-open** for every coordinate field;
//! the `airr_strict=True` export flag in `result.py` does the
//! 1-based-inclusive conversion at TSV/CSV/DataFrame time.

mod builder;
mod junction;
mod locus;
mod projection;
mod record;
mod segment_projection;
mod sequence;
mod trace_fields;
mod validate;
mod walk;

pub use builder::build_airr_record;
pub use record::AirrRecord;
pub use validate::{
    validate_airr_record, AlleleOrderReason, PairedEndField, PairedEndRead,
    ProductiveDecidedBy, RecordValidationIssue,
};

#[cfg(test)]
mod tests;
