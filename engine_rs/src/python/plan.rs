//! `PyPassPlan` — Python-facing builder for [`crate::pass::PassPlan`].
//!
//! Rather than wrapping each concrete pass type as its own
//! `#[pyclass]` (which would explode into a dozen wrappers and a
//! parallel set of `Distribution` wrappers), this file exposes a
//! single builder with typed `push_*` methods. Each method takes
//! the parameters needed to construct the matching Rust pass and
//! appends it to the inner `PassPlan` in one step.
//!
//! Trade-off: Python users can only construct passes the engine
//! already knows how to assemble (recombination + the canonical NP
//! / trim / assemble triples). Custom passes are not addressable
//! from Python at this phase — they require Rust changes. This is
//! the right scope for F.4: the experiment DSL we'll layer on top
//! in F.5 is exactly the set of recombination passes the engine
//! cares about.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use crate::assignment::TrimEnd;
use crate::dist::{AllelePoolDist, EmpiricalLengthDist, UniformBase};
use crate::ir::Segment;
use crate::pass::PassPlan;
use crate::passes::{AssembleSegmentPass, GenerateNPPass, SampleAllelePass, TrimPass};

use super::refdata::PyRefDataConfig;

/// A `PassPlan` constructed from Python.
///
/// Pushed passes are owned by the plan and the plan is borrowed by
/// runners (`PassRuntime::execute*` takes `&PassPlan`). Marked
/// `unsendable` because the `Pass` trait isn't `Send` — the
/// boxed-dyn passes inside the plan can't be moved between
/// threads. Plain Python use is single-threaded so this is fine.
#[pyclass(name = "PassPlan", module = "genairr_engine", unsendable)]
pub struct PyPassPlan {
    inner: PassPlan,
}

impl PyPassPlan {
    /// Borrow the inner `PassPlan`. Used by runners.
    pub(crate) fn inner(&self) -> &PassPlan {
        &self.inner
    }
}

#[pymethods]
impl PyPassPlan {
    /// An empty plan. Use `push_*` methods to populate.
    #[new]
    fn new() -> Self {
        Self {
            inner: PassPlan::new(),
        }
    }

    /// Number of passes currently in the plan.
    fn __len__(&self) -> usize {
        self.inner.len()
    }

    /// Whether the plan contains zero passes.
    fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }

    /// Append a `SampleAllelePass` for `segment`. Requires the active
    /// `refdata` so the per-segment uniform allele distribution can
    /// be built at push time.
    ///
    /// Errors:
    /// - `ValueError` when `segment` isn't `"V"`, `"D"`, or `"J"`.
    /// - `ValueError` when the corresponding pool in `refdata` is empty.
    fn push_sample_allele(
        &mut self,
        segment: &str,
        refdata: &PyRefDataConfig,
    ) -> PyResult<()> {
        let seg = parse_recombinable_segment(segment)?;
        let cfg = refdata.inner();
        let pool = cfg.pool_for(seg).ok_or_else(|| {
            PyValueError::new_err(format!("no pool for segment {:?}", segment))
        })?;
        if pool.is_empty() {
            return Err(PyValueError::new_err(format!(
                "{} pool is empty in supplied refdata",
                segment
            )));
        }
        self.inner.push(Box::new(SampleAllelePass::new(
            seg,
            Box::new(AllelePoolDist::uniform(pool)),
        )));
        Ok(())
    }

    /// Append an `AssembleSegmentPass` for `segment`. The matching
    /// `SampleAllelePass` must already be earlier in the plan
    /// (otherwise the assembler will fail at execute time with a
    /// "no allele assigned" panic).
    ///
    /// Errors: `ValueError` when `segment` isn't `"V"`, `"D"`, or `"J"`.
    fn push_assemble(&mut self, segment: &str) -> PyResult<()> {
        let seg = parse_recombinable_segment(segment)?;
        self.inner.push(Box::new(AssembleSegmentPass::new(seg)));
        Ok(())
    }

    /// Append a `GenerateNPPass` for `np_segment` (`"NP1"` or `"NP2"`).
    /// `length_pairs` is a list of `(length, weight)` tuples that
    /// defines the empirical NP-length distribution; bases are drawn
    /// from `UniformBase`.
    ///
    /// Errors:
    /// - `ValueError` when `np_segment` isn't `"NP1"` or `"NP2"`.
    /// - `ValueError` when `length_pairs` is empty.
    /// - `ValueError` when any pair has a negative length or
    ///   non-finite / non-positive weight (propagated from
    ///   `EmpiricalLengthDist::from_pairs`).
    fn push_generate_np(
        &mut self,
        np_segment: &str,
        length_pairs: Vec<(i64, f64)>,
    ) -> PyResult<()> {
        let seg = parse_np_segment(np_segment)?;
        if length_pairs.is_empty() {
            return Err(PyValueError::new_err(
                "length_pairs must contain at least one (length, weight) entry",
            ));
        }
        // EmpiricalLengthDist::from_pairs panics on bad input — we
        // accept that for now; the Rust constructor's validation is
        // shared with all paths.
        self.inner.push(Box::new(GenerateNPPass::new(
            seg,
            Box::new(EmpiricalLengthDist::from_pairs(length_pairs)),
            Box::new(UniformBase),
        )));
        Ok(())
    }

    /// Append a `TrimPass` for `(segment, end)`. `end` is `"5"` or
    /// `"3"` (the prime end being trimmed). `length_pairs` defines
    /// the trim-amount distribution.
    ///
    /// Errors:
    /// - `ValueError` when `segment` isn't `"V"`, `"D"`, or `"J"`.
    /// - `ValueError` when `end` isn't `"5"` or `"3"`.
    /// - `ValueError` when `length_pairs` is empty.
    fn push_trim(
        &mut self,
        segment: &str,
        end: &str,
        length_pairs: Vec<(i64, f64)>,
    ) -> PyResult<()> {
        let seg = parse_recombinable_segment(segment)?;
        let trim_end = parse_trim_end(end)?;
        if length_pairs.is_empty() {
            return Err(PyValueError::new_err(
                "length_pairs must contain at least one (length, weight) entry",
            ));
        }
        self.inner.push(Box::new(TrimPass::new(
            seg,
            trim_end,
            Box::new(EmpiricalLengthDist::from_pairs(length_pairs)),
        )));
        Ok(())
    }

    fn __repr__(&self) -> String {
        format!("<PassPlan len={}>", self.inner.len())
    }
}

fn parse_recombinable_segment(s: &str) -> PyResult<Segment> {
    match s.to_ascii_uppercase().as_str() {
        "V" => Ok(Segment::V),
        "D" => Ok(Segment::D),
        "J" => Ok(Segment::J),
        other => Err(PyValueError::new_err(format!(
            "segment must be 'V', 'D', or 'J' (got {:?})",
            other
        ))),
    }
}

fn parse_np_segment(s: &str) -> PyResult<Segment> {
    match s.to_ascii_uppercase().as_str() {
        "NP1" => Ok(Segment::Np1),
        "NP2" => Ok(Segment::Np2),
        other => Err(PyValueError::new_err(format!(
            "np_segment must be 'NP1' or 'NP2' (got {:?})",
            other
        ))),
    }
}

fn parse_trim_end(s: &str) -> PyResult<TrimEnd> {
    match s {
        "5" | "5'" | "five" | "FIVE" | "Five" => Ok(TrimEnd::Five),
        "3" | "3'" | "three" | "THREE" | "Three" => Ok(TrimEnd::Three),
        other => Err(PyValueError::new_err(format!(
            "end must be '5' or '3' (got {:?})",
            other
        ))),
    }
}
