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
use crate::passes::{
    AssembleSegmentPass, ContaminantPass, GenerateNPPass, IndelPass, PCRErrorPass,
    QualityErrorPass, S5FMutationPass, SampleAllelePass, TrimPass, UniformMutationPass,
};
use crate::s5f::{S5FKernel, S5F_NUM_CONTEXTS, S5F_SUBSTITUTION_LEN};

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
    /// `refdata` so the per-segment allele distribution can be built
    /// at push time.
    ///
    /// When `allowed_ids` is `None`, the distribution is uniform over
    /// the full pool. When `allowed_ids` is `Some(list)`, the
    /// distribution is uniform over only the listed IDs — the
    /// allele-locking path used by `Experiment.using(...)`.
    ///
    /// Errors:
    /// - `ValueError` when `segment` isn't `"V"`, `"D"`, or `"J"`.
    /// - `ValueError` when the corresponding pool in `refdata` is empty.
    /// - `ValueError` when `allowed_ids` is empty, contains
    ///   out-of-range IDs, or contains duplicates.
    #[pyo3(signature = (segment, refdata, allowed_ids = None))]
    fn push_sample_allele(
        &mut self,
        segment: &str,
        refdata: &PyRefDataConfig,
        allowed_ids: Option<Vec<u32>>,
    ) -> PyResult<()> {
        use crate::refdata::AlleleId;

        let seg = parse_recombinable_segment(segment)?;
        let cfg = refdata.inner();
        let pool = cfg
            .pool_for(seg)
            .ok_or_else(|| PyValueError::new_err(format!("no pool for segment {:?}", segment)))?;
        if pool.is_empty() {
            return Err(PyValueError::new_err(format!(
                "{} pool is empty in supplied refdata",
                segment
            )));
        }

        let dist: Box<AllelePoolDist> = match allowed_ids {
            None => Box::new(AllelePoolDist::uniform(pool)),
            Some(ids) => {
                if ids.is_empty() {
                    return Err(PyValueError::new_err(format!(
                        "{} allowed_ids must contain at least one id",
                        segment
                    )));
                }
                let pool_len = pool.len() as u32;
                let mut seen = std::collections::HashSet::with_capacity(ids.len());
                for raw in &ids {
                    if *raw >= pool_len {
                        return Err(PyValueError::new_err(format!(
                            "{} allowed_id {} out of range (pool size {})",
                            segment, raw, pool_len
                        )));
                    }
                    if !seen.insert(*raw) {
                        return Err(PyValueError::new_err(format!(
                            "{} allowed_id {} appears more than once",
                            segment, raw
                        )));
                    }
                }
                let allele_ids: Vec<AlleleId> = ids.into_iter().map(AlleleId::new).collect();
                Box::new(AllelePoolDist::restricted_uniform(pool, allele_ids))
            }
        };

        self.inner.push(Box::new(SampleAllelePass::new(seg, dist)));
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

    /// Append a `UniformMutationPass`. ``count_pairs`` is a
    /// ``(count, weight)`` distribution over the number of mutations
    /// applied; each mutation samples a position uniformly across
    /// the assembled pool and substitutes the base with a uniformly
    /// drawn A/C/G/T.
    ///
    /// Errors: ``ValueError`` when ``count_pairs`` is empty.
    fn push_mutate_uniform(&mut self, count_pairs: Vec<(i64, f64)>) -> PyResult<()> {
        if count_pairs.is_empty() {
            return Err(PyValueError::new_err(
                "count_pairs must contain at least one (count, weight) entry",
            ));
        }
        self.inner.push(Box::new(UniformMutationPass::new(
            Box::new(EmpiricalLengthDist::from_pairs(count_pairs)),
            Box::new(UniformBase),
        )));
        Ok(())
    }

    /// Append an `S5FMutationPass` with the given S5F kernel data.
    ///
    /// ``mutability`` is a flat ``Vec<f64>`` of length 1024 (one entry
    /// per 4⁵ canonical 5-mer context, indexed via
    /// ``S5FKernel::encode_context``). ``substitution`` is a flat
    /// ``Vec<f64>`` of length 4096 (1024 contexts × 4 destination
    /// bases A/C/G/T). ``count_pairs`` follows the same shape as
    /// ``push_mutate_uniform``.
    ///
    /// Errors: ``ValueError`` when ``count_pairs`` is empty or kernel
    /// dimensions are wrong.
    fn push_mutate_s5f(
        &mut self,
        count_pairs: Vec<(i64, f64)>,
        mutability: Vec<f64>,
        substitution: Vec<f64>,
    ) -> PyResult<()> {
        if count_pairs.is_empty() {
            return Err(PyValueError::new_err(
                "count_pairs must contain at least one (count, weight) entry",
            ));
        }
        if mutability.len() != S5F_NUM_CONTEXTS {
            return Err(PyValueError::new_err(format!(
                "mutability must have length {} (got {})",
                S5F_NUM_CONTEXTS,
                mutability.len()
            )));
        }
        if substitution.len() != S5F_SUBSTITUTION_LEN {
            return Err(PyValueError::new_err(format!(
                "substitution must have length {} (got {})",
                S5F_SUBSTITUTION_LEN,
                substitution.len()
            )));
        }
        let kernel = S5FKernel::new(mutability, substitution);
        self.inner.push(Box::new(S5FMutationPass::new(
            kernel,
            Box::new(EmpiricalLengthDist::from_pairs(count_pairs)),
        )));
        Ok(())
    }

    /// Append a `PCRErrorPass`. ``count_pairs`` is the empirical
    /// distribution over the number of PCR-induced base substitution
    /// errors per simulation; each error samples a uniform position
    /// and replaces the base with a uniform A/C/G/T draw.
    fn push_corrupt_pcr(&mut self, count_pairs: Vec<(i64, f64)>) -> PyResult<()> {
        if count_pairs.is_empty() {
            return Err(PyValueError::new_err(
                "count_pairs must contain at least one (count, weight) entry",
            ));
        }
        self.inner.push(Box::new(PCRErrorPass::new(
            Box::new(EmpiricalLengthDist::from_pairs(count_pairs)),
            Box::new(UniformBase),
        )));
        Ok(())
    }

    /// Append a `QualityErrorPass`. Same shape as PCR but each
    /// substitution writes the destination base as **lowercase** to
    /// preserve the V5 sequencing-error convention (uppercase =
    /// germline, lowercase = corrupted).
    fn push_corrupt_quality(&mut self, count_pairs: Vec<(i64, f64)>) -> PyResult<()> {
        if count_pairs.is_empty() {
            return Err(PyValueError::new_err(
                "count_pairs must contain at least one (count, weight) entry",
            ));
        }
        self.inner.push(Box::new(QualityErrorPass::new(
            Box::new(EmpiricalLengthDist::from_pairs(count_pairs)),
            Box::new(UniformBase),
        )));
        Ok(())
    }

    /// Append an `IndelPass`. ``count_pairs`` is the empirical
    /// distribution over the total number of indel events per
    /// simulation; each event independently chooses insertion vs.
    /// deletion with probability ``insertion_prob`` (defaults to
    /// 0.5). Insertions sample a uniform A/C/G/T base.
    ///
    /// Errors: ``ValueError`` when ``insertion_prob`` is outside
    /// ``[0.0, 1.0]`` or non-finite, or when ``count_pairs`` is empty.
    #[pyo3(signature = (count_pairs, *, insertion_prob=0.5))]
    fn push_corrupt_indel(
        &mut self,
        count_pairs: Vec<(i64, f64)>,
        insertion_prob: f64,
    ) -> PyResult<()> {
        if count_pairs.is_empty() {
            return Err(PyValueError::new_err(
                "count_pairs must contain at least one (count, weight) entry",
            ));
        }
        if !insertion_prob.is_finite() || !(0.0..=1.0).contains(&insertion_prob) {
            return Err(PyValueError::new_err(format!(
                "insertion_prob must be a finite number in [0.0, 1.0], got {}",
                insertion_prob
            )));
        }
        self.inner.push(Box::new(IndelPass::new(
            Box::new(EmpiricalLengthDist::from_pairs(count_pairs)),
            insertion_prob,
            Box::new(UniformBase),
        )));
        Ok(())
    }

    /// Append a `ContaminantPass`. With probability ``apply_prob``
    /// the entire assembled pool is overwritten with uniform
    /// A/C/G/T bases (modelling primer dimers, bacterial DNA, or
    /// any non-receptor sequence in the library).
    ///
    /// Errors: ``ValueError`` when ``apply_prob`` is outside
    /// ``[0.0, 1.0]`` or non-finite.
    fn push_corrupt_contaminant(&mut self, apply_prob: f64) -> PyResult<()> {
        if !apply_prob.is_finite() || !(0.0..=1.0).contains(&apply_prob) {
            return Err(PyValueError::new_err(format!(
                "apply_prob must be a finite number in [0.0, 1.0], got {}",
                apply_prob
            )));
        }
        self.inner.push(Box::new(ContaminantPass::new(
            apply_prob,
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
