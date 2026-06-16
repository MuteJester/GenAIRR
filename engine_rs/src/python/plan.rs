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
use crate::compiled::ExecutionPolicy;
use crate::dist::{AllelePoolDist, EmpiricalLengthDist, UniformBase};
use crate::ir::Segment;
use crate::pass::PassPlan;
use crate::passes::mutate::{SegmentRateWeights, VSubregionRateWeights};
use crate::passes::{
    AssembleSegmentPass, ContaminantPass, EndLossPass, GenerateNPPass, IndelPass, InvertDPass,
    LossEnd, NCorruptionPass, PCRErrorPass, PairedEndLayoutSpec, PairedEndSamplingPass,
    QualityErrorPass, ReceptorRevisionPass, RevCompPass, S5FMutationPass, SampleAllelePass,
    TrimPass, UniformMutationPass,
};
use crate::s5f::{S5FKernel, S5F_NUM_CONTEXTS, S5F_SUBSTITUTION_LEN};

/// Materialise the canonical Python ``[v, d, j, np]`` rate vector
/// into a [`SegmentRateWeights`]. ``None`` → flat defaults (all
/// 1.0); ``Some(vec)`` must be exactly four positive entries in
/// (v, d, j, np) order. The DSL-boundary validator in
/// ``experiment.py`` enforces these invariants before calling
/// across the bridge; this helper re-validates defensively so a
/// custom Python wrapper that bypasses the DSL still surfaces a
/// clear ``ValueError`` rather than panicking inside the pass.
fn build_segment_rate_weights(
    rates: Option<Vec<f64>>,
) -> PyResult<SegmentRateWeights> {
    let Some(vec) = rates else {
        return Ok(SegmentRateWeights::default());
    };
    if vec.len() != 4 {
        return Err(PyValueError::new_err(format!(
            "segment_rates must be a length-4 list [v, d, j, np], got {} entries",
            vec.len()
        )));
    }
    for (i, &r) in vec.iter().enumerate() {
        if !r.is_finite() || r < 0.0 {
            return Err(PyValueError::new_err(format!(
                "segment_rates[{}] must be a finite non-negative number, got {}",
                i, r
            )));
        }
    }
    if vec.iter().copied().sum::<f64>() <= 0.0 {
        return Err(PyValueError::new_err(
            "segment_rates: at least one bucket must have a positive rate",
        ));
    }
    Ok(SegmentRateWeights::from_tuple(vec[0], vec[1], vec[2], vec[3]))
}

/// Materialise the canonical Python ``[fwr1, cdr1, fwr2, cdr2,
/// fwr3]`` rate vector into a [`VSubregionRateWeights`].
/// ``None`` → flat defaults (all 1.0); ``Some(vec)`` must be
/// exactly five finite non-negative entries with at least one
/// strictly positive in
/// (FWR1, CDR1, FWR2, CDR2, FWR3) order. Same defence-in-depth
/// re-validation as [`build_segment_rate_weights`] — the
/// `_validate_v_subregion_rates` boundary in `experiment.py` is
/// the canonical authoritative validator; this helper guards a
/// custom Python wrapper that bypasses the DSL.
fn build_v_subregion_rate_weights(
    rates: Option<Vec<f64>>,
) -> PyResult<VSubregionRateWeights> {
    let Some(vec) = rates else {
        return Ok(VSubregionRateWeights::default());
    };
    if vec.len() != 5 {
        return Err(PyValueError::new_err(format!(
            "v_subregion_rates must be a length-5 list [fwr1, cdr1, fwr2, cdr2, fwr3], got {} entries",
            vec.len()
        )));
    }
    for (i, &r) in vec.iter().enumerate() {
        if !r.is_finite() || r < 0.0 {
            return Err(PyValueError::new_err(format!(
                "v_subregion_rates[{}] must be a finite non-negative number, got {}",
                i, r
            )));
        }
    }
    if vec.iter().copied().sum::<f64>() <= 0.0 {
        return Err(PyValueError::new_err(
            "v_subregion_rates: at least one label must have a positive rate",
        ));
    }
    Ok(VSubregionRateWeights::from_tuple(
        vec[0], vec[1], vec[2], vec[3], vec[4],
    ))
}

use super::compiled::PyCompiledSimulator;
use super::contract::PyContractSet;
use super::refdata::PyRefDataConfig;

/// A `PassPlan` constructed from Python.
///
/// Pushed passes are owned by the plan and the plan is borrowed by
/// runners (`PassRuntime::execute*` takes `&PassPlan`). Marked
/// `unsendable` because the `Pass` trait isn't `Send` — the
/// boxed-dyn passes inside the plan can't be moved between
/// threads. Plain Python use is single-threaded so this is fine.
#[pyclass(name = "PassPlan", module = "GenAIRR._engine", unsendable)]
pub struct PyPassPlan {
    inner: Option<PassPlan>,
}

impl PyPassPlan {
    /// Borrow the inner `PassPlan`. Used by runners.
    pub(crate) fn inner(&self) -> PyResult<&PassPlan> {
        self.inner.as_ref().ok_or_else(consumed_plan_err)
    }

    fn inner_mut(&mut self) -> PyResult<&mut PassPlan> {
        self.inner.as_mut().ok_or_else(consumed_plan_err)
    }

    pub(crate) fn take_inner(&mut self) -> PyResult<PassPlan> {
        self.inner.take().ok_or_else(consumed_plan_err)
    }
}

fn consumed_plan_err() -> PyErr {
    PyValueError::new_err(
        "PassPlan has been consumed by compile(); create a new PassPlan for additional edits or one-shot runs",
    )
}

#[pymethods]
impl PyPassPlan {
    /// An empty plan. Use `push_*` methods to populate.
    #[new]
    fn new() -> Self {
        Self {
            inner: Some(PassPlan::new()),
        }
    }

    /// Number of passes currently in the plan.
    fn __len__(&self) -> PyResult<usize> {
        Ok(self.inner()?.len())
    }

    /// Whether the plan contains zero passes.
    fn is_empty(&self) -> PyResult<bool> {
        Ok(self.inner()?.is_empty())
    }

    /// Stable index of the most recently pushed pass. Use as a
    /// handle for [`Self::add_edge`] / [`Self::after`] /
    /// [`Self::before`]. Equal to `len() - 1` right after any
    /// `push_*` call. Errors if the plan is empty or has already
    /// been consumed by `compile()`.
    fn last_pushed_node(&self) -> PyResult<u32> {
        let plan = self.inner()?;
        let n = plan.len();
        if n == 0 {
            return Err(PyValueError::new_err(
                "plan is empty; push a pass first before calling last_pushed_node()",
            ));
        }
        Ok((n - 1) as u32)
    }

    /// Declare an explicit ordering edge: the pass at `from_idx`
    /// must execute before the pass at `to_idx`. The indices are
    /// node handles obtained from [`Self::last_pushed_node`] (or
    /// equivalently positions in push order, 0-based).
    ///
    /// Used for cross-cutting constraints that the auto-derived
    /// `Pass::requirements()` / `Pass::effects()` edges don't
    /// capture — e.g. "this corruption pass must come after that
    /// mutation pass" when neither produces an effect the other
    /// consumes.
    fn add_edge(&mut self, from_idx: u32, to_idx: u32) -> PyResult<()> {
        use crate::pass::NodeId;
        let plan = self.inner_mut()?;
        let n = plan.len() as u32;
        if from_idx >= n || to_idx >= n {
            return Err(PyValueError::new_err(format!(
                "add_edge: node index out of range (have {} passes, got {}/{})",
                n, from_idx, to_idx
            )));
        }
        if from_idx == to_idx {
            return Err(PyValueError::new_err(
                "add_edge: self-edge would create a cycle",
            ));
        }
        plan.add_edge(NodeId(from_idx), NodeId(to_idx));
        Ok(())
    }

    /// Shorthand for `add_edge(other_idx, self_idx)`: make
    /// `self_idx` execute after `other_idx`.
    fn after(&mut self, self_idx: u32, other_idx: u32) -> PyResult<()> {
        self.add_edge(other_idx, self_idx)
    }

    /// Shorthand for `add_edge(self_idx, other_idx)`: make
    /// `self_idx` execute before `other_idx`.
    fn before(&mut self, self_idx: u32, other_idx: u32) -> PyResult<()> {
        self.add_edge(self_idx, other_idx)
    }

    /// Return the node index of the first pass whose `name()`
    /// equals ``name``, or ``None`` if no such pass is in the plan.
    ///
    /// Used by lowering steps that need to attach an explicit
    /// schedule edge against a pass pushed by an *earlier* step in
    /// the pipeline IR — e.g. `_InvertDStep` finds the previously-
    /// pushed `assemble.d` and adds a `before(invert_d, assemble.d)`
    /// edge so InvertDPass commits its orientation decision before
    /// AssembleSegmentPass reads it.
    ///
    /// Linear scan; pass plans are tens-of-passes deep, the cost is
    /// negligible against pipeline-build time.
    fn find_node_by_name(&self, name: &str) -> PyResult<Option<u32>> {
        let plan = self.inner()?;
        for (i, pass) in plan.passes().iter().enumerate() {
            if pass.name() == name {
                return Ok(Some(i as u32));
            }
        }
        Ok(None)
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
    #[pyo3(signature = (segment, refdata, allowed_ids = None, weights = None))]
    fn push_sample_allele(
        &mut self,
        segment: &str,
        refdata: &PyRefDataConfig,
        allowed_ids: Option<Vec<u32>>,
        weights: Option<Vec<f64>>,
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

        // `allowed_ids` and `weights` are mutually exclusive.
        // The Python DSL never combines them; reject the combination
        // explicitly so a future caller doesn't silently ignore one.
        if allowed_ids.is_some() && weights.is_some() {
            return Err(PyValueError::new_err(format!(
                "{} push_sample_allele: pass either allowed_ids OR weights, not both",
                segment
            )));
        }

        let dist: Box<AllelePoolDist> = if let Some(w) = weights {
            if w.len() != pool.len() {
                return Err(PyValueError::new_err(format!(
                    "{} weights length ({}) must match pool size ({})",
                    segment,
                    w.len(),
                    pool.len()
                )));
            }
            for (i, weight) in w.iter().enumerate() {
                if !weight.is_finite() || *weight <= 0.0 {
                    return Err(PyValueError::new_err(format!(
                        "{} weight at index {} must be finite and > 0, got {}",
                        segment, i, weight
                    )));
                }
            }
            Box::new(AllelePoolDist::from_weights(pool, w))
        } else {
            match allowed_ids {
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
            }
        };

        self.inner_mut()?
            .push(Box::new(SampleAllelePass::new(seg, dist)));
        Ok(())
    }

    /// Append a single `SampleGenotypePass` that replaces the three flat
    /// `SampleAllelePass` passes when a genotype is attached. The pass
    /// draws the chromosome once, then V/(D)/J alleles from that
    /// chromosome's carried set (phased).
    ///
    /// `v`/`d`/`j` are flat rows `(haplotype, allele_id, copies, weight)`
    /// already resolved to this refdata's allele ids; rows are grouped by
    /// the gene the allele belongs to (via the segment's `GeneIndex`).
    #[pyo3(signature = (refdata, chromosome_weights, subject_id, source_hash, v, d, j, d_required))]
    #[allow(clippy::too_many_arguments)]
    fn push_genotype_recombine(
        &mut self,
        refdata: &PyRefDataConfig,
        chromosome_weights: (f32, f32),
        subject_id: Option<String>,
        source_hash: String,
        v: Vec<(u8, u32, u8, f32)>,
        d: Vec<(u8, u32, u8, f32)>,
        j: Vec<(u8, u32, u8, f32)>,
        d_required: bool,
    ) -> PyResult<()> {
        use crate::genotype::{GeneCopy, Genotype, Haplotype};
        use crate::ir::Segment;
        use crate::refdata::{AlleleId, GeneId, GeneIndex};
        use std::collections::HashMap;

        let cfg = refdata.inner();

        let build_index = |seg: Segment| -> PyResult<GeneIndex> {
            let pool = cfg
                .pool_for(seg)
                .ok_or_else(|| PyValueError::new_err(format!("no pool for segment {:?}", seg)))?;
            Ok(GeneIndex::build(pool))
        };
        let v_index = build_index(Segment::V)?;
        let j_index = build_index(Segment::J)?;
        let d_index = if d_required {
            Some(build_index(Segment::D)?)
        } else {
            None
        };

        let mut haps = [Haplotype::new(), Haplotype::new()];
        let fill = |haps: &mut [Haplotype; 2],
                    seg: Segment,
                    idx: &GeneIndex,
                    rows: &[(u8, u32, u8, f32)]|
         -> PyResult<()> {
            let pool_len = cfg.pool_for(seg).map(|p| p.len() as u32).unwrap_or(0);
            let mut grouped: HashMap<(u8, u32), Vec<GeneCopy>> = HashMap::new();
            for (h, aid, copies, weight) in rows {
                if *h > 1 {
                    return Err(PyValueError::new_err(format!(
                        "haplotype index must be 0 or 1, got {}",
                        h
                    )));
                }
                if *aid >= pool_len {
                    return Err(PyValueError::new_err(format!(
                        "{:?} allele id {} out of range (pool size {})",
                        seg, aid, pool_len
                    )));
                }
                let allele = AlleleId::new(*aid);
                let gene = idx.gene_of(allele);
                grouped.entry((*h, gene.index())).or_default().push(GeneCopy {
                    allele,
                    copies: *copies,
                    weight: *weight,
                });
            }
            for ((h, gene_idx), copies) in grouped {
                haps[h as usize].set(seg, GeneId::new(gene_idx), copies);
            }
            Ok(())
        };
        fill(&mut haps, Segment::V, &v_index, &v)?;
        fill(&mut haps, Segment::J, &j_index, &j)?;
        if let Some(di) = &d_index {
            fill(&mut haps, Segment::D, di, &d)?;
        }

        let genotype = std::sync::Arc::new(Genotype::new(
            haps,
            [chromosome_weights.0, chromosome_weights.1],
            subject_id,
            source_hash,
        ));
        self.inner_mut()?
            .push(Box::new(crate::passes::sample_genotype::SampleGenotypePass::new(
                genotype,
                d_required,
                Vec::new(),
                Vec::new(),
                Vec::new(),
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
        self.inner_mut()?
            .push(Box::new(AssembleSegmentPass::new(seg)));
        Ok(())
    }

    /// Append a `PAdditionPass` for `end` (`"V_3"`, `"D_5"`,
    /// `"D_3"`, or `"J_5"`). `length_pairs` is a list of
    /// `(length, weight)` tuples that defines the empirical
    /// per-end P-length distribution. Bases are deterministic
    /// from the source allele's post-trim, post-orientation
    /// coding flank (see
    /// [`docs/p_nucleotide_design.md`](../../../../docs/p_nucleotide_design.md)).
    ///
    /// Errors:
    /// - `ValueError` when `end` isn't one of the four canonical labels.
    /// - `ValueError` when `length_pairs` is empty.
    fn push_p_addition(
        &mut self,
        end: &str,
        length_pairs: Vec<(i64, f64)>,
    ) -> PyResult<()> {
        let p_end = parse_p_end(end)?;
        if length_pairs.is_empty() {
            return Err(PyValueError::new_err(
                "length_pairs must contain at least one (length, weight) entry",
            ));
        }
        let length_dist = Box::new(EmpiricalLengthDist::from_pairs(length_pairs));
        self.inner_mut()?.push(Box::new(
            crate::passes::PAdditionPass::new(p_end, length_dist),
        ));
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
    /// - compile-time validation rejects negative lengths when the
    ///   plan is compiled.
    /// - non-finite / non-positive weights panic in
    ///   `EmpiricalLengthDist::from_pairs` at construction time.
    #[pyo3(signature = (np_segment, length_pairs, *, base_pairs=None, markov_transitions=None))]
    fn push_generate_np(
        &mut self,
        np_segment: &str,
        length_pairs: Vec<(i64, f64)>,
        base_pairs: Option<Vec<(u8, f64)>>,
        markov_transitions: Option<Vec<Vec<(u8, f64)>>>,
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
        //
        // Dispatch (Slice — Markov NP Base Generator):
        //
        //   * ``base_pairs=None`` and ``markov_transitions=None`` →
        //     pre-slice ``UniformBase`` / ``UniformNpGenerator``
        //     (byte-identical to legacy plan signatures).
        //   * ``base_pairs=Some(...)`` and
        //     ``markov_transitions=None`` → cartridge-owned
        //     ``empirical_first_base`` (Slice — Typed NP base model)
        //     via ``CategoricalNpGenerator``.
        //   * ``markov_transitions=Some(rows)`` → full 1-step
        //     ``MarkovBaseGenerator``. ``base_pairs`` is required
        //     in this case (carries the first-base row); the four
        //     rows are in canonical A/C/G/T from-base order; each
        //     row is a ``(to_byte, weight)`` pair list covering
        //     the full A/C/G/T to-alphabet. The Python spec layer
        //     enforces row completeness — defence-in-depth checks
        //     here re-validate the bridge surface.
        let length_dist = Box::new(EmpiricalLengthDist::from_pairs(length_pairs));
        if let Some(rows) = markov_transitions {
            let first_pairs = base_pairs.ok_or_else(|| {
                PyValueError::new_err(
                    "markov_transitions requires base_pairs (the first-base row)",
                )
            })?;
            let first_base = parse_canonical_base_weights(
                &first_pairs,
                "base_pairs",
            )?;
            if rows.len() != 4 {
                return Err(PyValueError::new_err(format!(
                    "markov_transitions must have 4 rows (A/C/G/T from-bases), got {}",
                    rows.len()
                )));
            }
            let mut transition_rows: Vec<Vec<f64>> = Vec::with_capacity(4);
            for (from_idx, row) in rows.iter().enumerate() {
                let label = format!("markov_transitions[{}]", "ACGT".chars().nth(from_idx).unwrap());
                let row_weights = parse_canonical_base_weights(row, &label)?;
                transition_rows.push(row_weights.to_vec());
            }
            let generator = crate::passes::generate_np::MarkovBaseGenerator::from_first_and_rows(
                first_base.to_vec(),
                transition_rows,
            );
            self.inner_mut()?.push(Box::new(
                crate::passes::generate_np::GenerateNPPass::with_generator(
                    seg,
                    length_dist,
                    Box::new(generator),
                ),
            ));
            return Ok(());
        }
        let base_dist: Box<dyn crate::dist::Distribution<Output = u8>> =
            match base_pairs {
                None => Box::new(UniformBase),
                Some(pairs) => {
                    if pairs.is_empty() {
                        return Err(PyValueError::new_err(
                            "base_pairs must contain at least one (base, weight) entry",
                        ));
                    }
                    // Defence-in-depth re-validation. The Python DSL
                    // boundary (`NpBaseModelSpec`) is the authoritative
                    // validator; this catch defends against callers
                    // that bypass the DSL.
                    for (base, weight) in &pairs {
                        if !matches!(*base, b'A' | b'C' | b'G' | b'T') {
                            return Err(PyValueError::new_err(format!(
                                "base_pairs: only canonical A/C/G/T are allowed, got byte {} (0x{:02X})",
                                base, base
                            )));
                        }
                        if !weight.is_finite() || *weight < 0.0 {
                            return Err(PyValueError::new_err(format!(
                                "base_pairs: weights must be finite and non-negative, got {}",
                                weight
                            )));
                        }
                    }
                    let total: f64 = pairs.iter().map(|(_, w)| w).sum();
                    if total <= 0.0 {
                        return Err(PyValueError::new_err(
                            "base_pairs: at least one weight must be strictly positive",
                        ));
                    }
                    Box::new(crate::dist::CategoricalBase::from_pairs(pairs))
                }
            };
        self.inner_mut()?.push(Box::new(GenerateNPPass::new(
            seg,
            length_dist,
            base_dist,
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
    #[pyo3(signature = (count_pairs, *, segment_rates=None, v_subregion_rates=None))]
    fn push_mutate_uniform(
        &mut self,
        count_pairs: Vec<(i64, f64)>,
        segment_rates: Option<Vec<f64>>,
        v_subregion_rates: Option<Vec<f64>>,
    ) -> PyResult<()> {
        if count_pairs.is_empty() {
            return Err(PyValueError::new_err(
                "count_pairs must contain at least one (count, weight) entry",
            ));
        }
        let rates = build_segment_rate_weights(segment_rates)?;
        let v_sub_rates = build_v_subregion_rate_weights(v_subregion_rates)?;
        self.inner_mut()?.push(Box::new(
            UniformMutationPass::new(
                Box::new(EmpiricalLengthDist::from_pairs(count_pairs)),
                Box::new(UniformBase),
            )
            .with_segment_rates(rates)
            .with_v_subregion_rates(v_sub_rates),
        ));
        Ok(())
    }

    /// Append a rate-driven ``UniformMutationPass``. ``rate`` is the
    /// per-base mutation rate (e.g. ``0.03`` for 3 % SHM); the count
    /// per record is drawn from ``Poisson(rate × pool_len)`` so the
    /// number of mutations scales with each record's own sequence
    /// length. Companion to :meth:`push_mutate_uniform`.
    ///
    /// Errors: ``ValueError`` when ``rate`` is outside ``[0.0, 1.0]``
    /// or non-finite.
    #[pyo3(signature = (rate, *, segment_rates=None, v_subregion_rates=None))]
    fn push_mutate_uniform_rate(
        &mut self,
        rate: f64,
        segment_rates: Option<Vec<f64>>,
        v_subregion_rates: Option<Vec<f64>>,
    ) -> PyResult<()> {
        if !rate.is_finite() || !(0.0..=1.0).contains(&rate) {
            return Err(PyValueError::new_err(format!(
                "rate must be a finite value in [0.0, 1.0], got {}",
                rate
            )));
        }
        let rates = build_segment_rate_weights(segment_rates)?;
        let v_sub_rates = build_v_subregion_rate_weights(v_subregion_rates)?;
        self.inner_mut()?.push(Box::new(
            UniformMutationPass::new_rate(rate, Box::new(UniformBase))
                .with_segment_rates(rates)
                .with_v_subregion_rates(v_sub_rates),
        ));
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
    /// ``kernel_name`` is the short identifier of the S5F kernel
    /// being installed (``"hh_s5f"`` / ``"hkl_s5f"`` / …). Threaded
    /// through from the DSL's ``_MutateStep.s5f_model_name`` so the
    /// plan signature can distinguish two passes that share rate +
    /// count but disagree on which kernel's mutability table was
    /// loaded. Optional for backwards compatibility; ``None`` is
    /// recorded in the plan signature as ``kernel=unnamed``.
    ///
    /// Errors: ``ValueError`` when ``count_pairs`` is empty or kernel
    /// dimensions are wrong.
    #[pyo3(signature = (count_pairs, mutability, substitution, *, segment_rates=None, v_subregion_rates=None, kernel_name=None))]
    fn push_mutate_s5f(
        &mut self,
        count_pairs: Vec<(i64, f64)>,
        mutability: Vec<f64>,
        substitution: Vec<f64>,
        segment_rates: Option<Vec<f64>>,
        v_subregion_rates: Option<Vec<f64>>,
        kernel_name: Option<String>,
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
        let rates = build_segment_rate_weights(segment_rates)?;
        let v_sub_rates = build_v_subregion_rate_weights(v_subregion_rates)?;
        let kernel = S5FKernel::new(mutability, substitution);
        let mut pass = S5FMutationPass::new(
            kernel,
            Box::new(EmpiricalLengthDist::from_pairs(count_pairs)),
        )
        .with_segment_rates(rates)
        .with_v_subregion_rates(v_sub_rates);
        if let Some(name) = kernel_name {
            pass = pass.with_kernel_name(name);
        }
        self.inner_mut()?.push(Box::new(pass));
        Ok(())
    }

    /// Append a rate-driven ``S5FMutationPass`` with the given S5F
    /// kernel. ``rate`` is the per-base mutation rate (e.g. ``0.03``
    /// for 3 % SHM); the count per record is drawn from
    /// ``Poisson(rate × pool_len)`` so the number of mutations
    /// scales with each record's own sequence length. Companion to
    /// :meth:`push_mutate_s5f`.
    ///
    /// Errors: ``ValueError`` when ``rate`` is outside ``[0.0, 1.0]``
    /// or kernel dimensions are wrong.
    #[pyo3(signature = (rate, mutability, substitution, *, segment_rates=None, v_subregion_rates=None, kernel_name=None))]
    fn push_mutate_s5f_rate(
        &mut self,
        rate: f64,
        mutability: Vec<f64>,
        substitution: Vec<f64>,
        segment_rates: Option<Vec<f64>>,
        v_subregion_rates: Option<Vec<f64>>,
        kernel_name: Option<String>,
    ) -> PyResult<()> {
        if !rate.is_finite() || !(0.0..=1.0).contains(&rate) {
            return Err(PyValueError::new_err(format!(
                "rate must be a finite value in [0.0, 1.0], got {}",
                rate
            )));
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
        let rates = build_segment_rate_weights(segment_rates)?;
        let v_sub_rates = build_v_subregion_rate_weights(v_subregion_rates)?;
        let kernel = S5FKernel::new(mutability, substitution);
        let mut pass = S5FMutationPass::new_rate(kernel, rate)
            .with_segment_rates(rates)
            .with_v_subregion_rates(v_sub_rates);
        if let Some(name) = kernel_name {
            pass = pass.with_kernel_name(name);
        }
        self.inner_mut()?.push(Box::new(pass));
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
        self.inner_mut()?.push(Box::new(PCRErrorPass::new(
            Box::new(EmpiricalLengthDist::from_pairs(count_pairs)),
            Box::new(UniformBase),
        )));
        Ok(())
    }

    /// Append a rate-driven ``PCRErrorPass``. ``rate`` is the per-base
    /// PCR substitution error probability (e.g. ``1e-4``); the count
    /// per record is drawn from ``Poisson(rate × pool_len)`` so the
    /// number of errors scales with each record's own sequence
    /// length. Companion to :meth:`push_corrupt_pcr`.
    ///
    /// Errors: ``ValueError`` when ``rate`` is outside ``[0.0, 1.0]``
    /// or non-finite.
    fn push_corrupt_pcr_rate(&mut self, rate: f64) -> PyResult<()> {
        if !rate.is_finite() || !(0.0..=1.0).contains(&rate) {
            return Err(PyValueError::new_err(format!(
                "rate must be a finite value in [0.0, 1.0], got {}",
                rate
            )));
        }
        self.inner_mut()?.push(Box::new(PCRErrorPass::new_rate(
            rate,
            Box::new(UniformBase),
        )));
        Ok(())
    }

    /// Append a `QualityErrorPass`. Same shape as PCR but each
    /// substitution writes the destination base as **lowercase** to
    /// preserve the sequencing-error convention (uppercase =
    /// germline, lowercase = corrupted).
    fn push_corrupt_quality(&mut self, count_pairs: Vec<(i64, f64)>) -> PyResult<()> {
        if count_pairs.is_empty() {
            return Err(PyValueError::new_err(
                "count_pairs must contain at least one (count, weight) entry",
            ));
        }
        self.inner_mut()?.push(Box::new(QualityErrorPass::new(
            Box::new(EmpiricalLengthDist::from_pairs(count_pairs)),
            Box::new(UniformBase),
        )));
        Ok(())
    }

    /// Append a rate-driven ``QualityErrorPass``. ``rate`` is the
    /// per-base sequencing error probability (e.g. ``0.001`` for
    /// Q30); the count per record is drawn from
    /// ``Poisson(rate × pool_len)`` so the number of errors scales
    /// with each record's own sequence length. Companion to
    /// :meth:`push_corrupt_quality`.
    ///
    /// Errors: ``ValueError`` when ``rate`` is outside ``[0.0, 1.0]``
    /// or non-finite.
    fn push_corrupt_quality_rate(&mut self, rate: f64) -> PyResult<()> {
        if !rate.is_finite() || !(0.0..=1.0).contains(&rate) {
            return Err(PyValueError::new_err(format!(
                "rate must be a finite value in [0.0, 1.0], got {}",
                rate
            )));
        }
        self.inner_mut()?.push(Box::new(QualityErrorPass::new_rate(
            rate,
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
        self.inner_mut()?.push(Box::new(IndelPass::new(
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
        self.inner_mut()?.push(Box::new(ContaminantPass::new(
            apply_prob,
            Box::new(UniformBase),
        )));
        Ok(())
    }

    /// Append an `NCorruptionPass`. ``count_pairs`` gives the
    /// per-simulation distribution over the number of `N`
    /// substitutions. Each substitution picks a uniform pool
    /// position and overwrites the base with `N`.
    ///
    /// Errors: ``ValueError`` when ``count_pairs`` is empty.
    fn push_corrupt_ns(&mut self, count_pairs: Vec<(i64, f64)>) -> PyResult<()> {
        if count_pairs.is_empty() {
            return Err(PyValueError::new_err(
                "count_pairs must contain at least one (count, weight) entry",
            ));
        }
        self.inner_mut()?
            .push(Box::new(NCorruptionPass::new(Box::new(
                EmpiricalLengthDist::from_pairs(count_pairs),
            ))));
        Ok(())
    }

    /// Append a `EndLossPass` for the 5' end. The
    /// `length_pairs` distribution gives the number of bases to
    /// strip from the start of the assembled sequence; the actual
    /// loss is clamped to the pool length.
    ///
    /// Errors: ``ValueError`` when ``length_pairs`` is empty.
    fn push_corrupt_5prime_loss(&mut self, length_pairs: Vec<(i64, f64)>) -> PyResult<()> {
        if length_pairs.is_empty() {
            return Err(PyValueError::new_err(
                "length_pairs must contain at least one (length, weight) entry",
            ));
        }
        self.inner_mut()?.push(Box::new(EndLossPass::new(
            LossEnd::Five,
            Box::new(EmpiricalLengthDist::from_pairs(length_pairs)),
        )));
        Ok(())
    }

    /// Append a `EndLossPass` for the 3' end.
    ///
    /// Errors: ``ValueError`` when ``length_pairs`` is empty.
    fn push_corrupt_3prime_loss(&mut self, length_pairs: Vec<(i64, f64)>) -> PyResult<()> {
        if length_pairs.is_empty() {
            return Err(PyValueError::new_err(
                "length_pairs must contain at least one (length, weight) entry",
            ));
        }
        self.inner_mut()?.push(Box::new(EndLossPass::new(
            LossEnd::Three,
            Box::new(EmpiricalLengthDist::from_pairs(length_pairs)),
        )));
        Ok(())
    }

    /// Append a `RevCompPass`. With probability
    /// ``apply_prob`` the AIRR record-builder will reverse-complement
    /// the `sequence`, `np1`, `np2`, `junction` fields and flip
    /// the corresponding pool-position coords; alignment / germline
    /// fields stay in forward orientation per the AIRR spec.
    ///
    /// Errors: ``ValueError`` when ``apply_prob`` is outside
    /// ``[0.0, 1.0]`` or non-finite.
    fn push_corrupt_rev_comp(&mut self, apply_prob: f64) -> PyResult<()> {
        if !apply_prob.is_finite() || !(0.0..=1.0).contains(&apply_prob) {
            return Err(PyValueError::new_err(format!(
                "apply_prob must be a finite number in [0.0, 1.0], got {}",
                apply_prob
            )));
        }
        self.inner_mut()?
            .push(Box::new(RevCompPass::new(apply_prob)));
        Ok(())
    }

    /// Append an `InvertDPass`. Records a per-simulation Bool at
    /// ``sample_allele.d.inverted`` and commits
    /// ``SegmentOrientation::ReverseComplement`` on the D allele
    /// when the draw is true. AssembleSegmentPass(D) then emits
    /// the reverse-complemented D slice (Slice B).
    ///
    /// Must be inserted **after** the SampleAllele(D) push and
    /// **before** the AssembleSegment(D) push. The schedule
    /// analyser auto-derives the lower edge through
    /// ``Pass::requirements()``; the upper edge against
    /// AssembleSegment(D) is the caller's responsibility — the
    /// Python lowering wires it via
    /// :meth:`PassPlan.before` against the
    /// :meth:`PassPlan.find_node_by_name` lookup.
    ///
    /// Errors: ``ValueError`` when ``prob`` is outside
    /// ``[0.0, 1.0]`` or non-finite.
    fn push_invert_d(&mut self, prob: f64) -> PyResult<()> {
        if !prob.is_finite() || !(0.0..=1.0).contains(&prob) {
            return Err(PyValueError::new_err(format!(
                "prob must be a finite number in [0.0, 1.0], got {}",
                prob
            )));
        }
        self.inner_mut()?.push(Box::new(InvertDPass::new(prob)));
        Ok(())
    }

    /// Append a `ReceptorRevisionPass`. Records a per-simulation Bool
    /// at `receptor_revision.applied`; on `true`, additionally records
    /// the replacement V allele at `receptor_revision.v_allele` and
    /// the derived 3' trim at `receptor_revision.v_trim_3`, then
    /// commits AssignmentChanged + TrimChanged + SegmentReplaced
    /// against V.
    ///
    /// The replacement V allele distribution is uniform over the V
    /// pool in `refdata`. Slice C's pass filters those candidates to
    /// the same-retained-length subset (`allele.len() >= old_v_len`)
    /// internally; this push site doesn't pre-narrow.
    ///
    /// Must be inserted **after** the AssembleSegment(V) push (and
    /// in v1, after AssembleSegment(J) so receptor revision sees the
    /// full V-NP1-D-NP2-J pool). The Python lowering wires this via
    /// the canonical V-NP1-D-NP2-J recombine sequence and inserts
    /// the call immediately after `push_assemble("J")`.
    ///
    /// Errors:
    /// - `ValueError` when `prob` is non-finite or outside `[0.0, 1.0]`.
    /// - `ValueError` when the V pool in `refdata` is empty.
    /// Append a `PairedEndSamplingPass` configured with the three
    /// integer distributions Slice D's `Experiment.paired_end(…)`
    /// lowers to. The pass records three Int trace records
    /// (`paired_end.r1_length`, `.r2_length`, `.insert_size`) per
    /// simulation; the AIRR builder reads them back at projection
    /// time and applies the
    /// `crate::airr_record::sequence::project_paired_end_layout`
    /// kernel that landed in Slice B.
    ///
    /// `*_pairs` are `(value, weight)` empirical distributions —
    /// the same shape `push_mutate_uniform` / `push_corrupt_*`
    /// already accept. Single-value distributions are passed as
    /// one-element lists.
    ///
    /// Errors:
    /// - `ValueError` when any of the three pair lists is empty.
    ///
    /// The pass itself runs additional per-sample
    /// relationship checks (positive lengths, non-negative insert,
    /// `r1/r2 <= insert`); those errors surface at execute time
    /// as `PassError::InvalidDistributionOutput`. The DSL boundary
    /// already pre-checks the fixed-value cases — see
    /// `Experiment.paired_end`.
    fn push_paired_end(
        &mut self,
        r1_length_pairs: Vec<(i64, f64)>,
        r2_length_pairs: Vec<(i64, f64)>,
        insert_size_pairs: Vec<(i64, f64)>,
    ) -> PyResult<()> {
        if r1_length_pairs.is_empty() {
            return Err(PyValueError::new_err(
                "paired_end r1_length distribution must contain at least one entry",
            ));
        }
        if r2_length_pairs.is_empty() {
            return Err(PyValueError::new_err(
                "paired_end r2_length distribution must contain at least one entry",
            ));
        }
        if insert_size_pairs.is_empty() {
            return Err(PyValueError::new_err(
                "paired_end insert_size distribution must contain at least one entry",
            ));
        }
        let spec = PairedEndLayoutSpec::new(
            Box::new(EmpiricalLengthDist::from_pairs(r1_length_pairs)),
            Box::new(EmpiricalLengthDist::from_pairs(r2_length_pairs)),
            Box::new(EmpiricalLengthDist::from_pairs(insert_size_pairs)),
        );
        self.inner_mut()?
            .push(Box::new(PairedEndSamplingPass::new(spec)));
        Ok(())
    }

    fn push_receptor_revision(
        &mut self,
        prob: f64,
        refdata: &PyRefDataConfig,
    ) -> PyResult<()> {
        use crate::ir::Segment;
        if !prob.is_finite() || !(0.0..=1.0).contains(&prob) {
            return Err(PyValueError::new_err(format!(
                "prob must be a finite number in [0.0, 1.0], got {}",
                prob
            )));
        }
        let cfg = refdata.inner();
        let pool = cfg.pool_for(Segment::V).ok_or_else(|| {
            PyValueError::new_err("receptor_revision requires a V pool in refdata")
        })?;
        if pool.is_empty() {
            return Err(PyValueError::new_err(
                "receptor_revision requires a non-empty V pool",
            ));
        }
        self.inner_mut()?.push(Box::new(ReceptorRevisionPass::new(
            prob,
            Box::new(AllelePoolDist::uniform(pool)),
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
        self.inner_mut()?.push(Box::new(TrimPass::new(
            seg,
            trim_end,
            Box::new(EmpiricalLengthDist::from_pairs(length_pairs)),
        )));
        Ok(())
    }

    /// Compile this builder into an owning `CompiledSimulator`.
    ///
    /// Compilation consumes the plan on success. This is intentional:
    /// the compiled simulator owns the pass IR and can run repeatedly
    /// without re-validating or re-borrowing the builder.
    ///
    /// ``allow_curatable_refdata=True`` runs the refdata gate under
    /// the lenient `AllowCuratable` mode: Fatal issues (empty
    /// required pool, duplicate names, invalid bytes, anchors out of
    /// bounds) still reject, but Curatable issues (V anchor not Cys,
    /// J anchor unexpected AA, missing V/J anchor) pass. Use when
    /// sampling from a real catalogue (bundled mouse_igh,
    /// human_tcrb) that includes pseudogene/ORF alleles. Default
    /// ``False`` — strict validation, every issue gates compile.
    #[pyo3(signature = (*, refdata=None, respect=None, strict=false, allow_curatable_refdata=false))]
    fn compile(
        &mut self,
        refdata: Option<&PyRefDataConfig>,
        respect: Option<&PyContractSet>,
        strict: bool,
        allow_curatable_refdata: bool,
    ) -> PyResult<PyCompiledSimulator> {
        let policy = if strict {
            ExecutionPolicy::Strict
        } else {
            ExecutionPolicy::Permissive
        };

        let refdata = refdata.map(|r| r.inner().clone());
        let contracts = respect.map(|c| c.inner().clone());
        PyCompiledSimulator::compile_from_plan(
            self,
            refdata,
            contracts,
            policy,
            allow_curatable_refdata,
        )
    }

    fn __repr__(&self) -> String {
        match &self.inner {
            Some(plan) => format!("<PassPlan len={}>", plan.len()),
            None => "<PassPlan consumed=true>".to_string(),
        }
    }
}

/// Parse a `(base_byte, weight)` pair list covering the canonical
/// A/C/G/T alphabet into a 4-element `[f64; 4]` indexed in A/C/G/T
/// canonical order. The Python `NpBaseModelSpec` is the
/// authoritative validator (per-base alphabet, weight finiteness,
/// row positivity); this re-validation is defence-in-depth for
/// callers that bypass the DSL boundary.
///
/// Requires every canonical base to be present at least once;
/// duplicates accumulate (matches `CategoricalBase::from_pairs`).
fn parse_canonical_base_weights(
    pairs: &[(u8, f64)],
    label: &str,
) -> PyResult<[f64; 4]> {
    let mut weights = [0.0f64; 4];
    let mut seen = [false; 4];
    for (base, weight) in pairs.iter() {
        let idx = match base {
            b'A' | b'a' => 0,
            b'C' | b'c' => 1,
            b'G' | b'g' => 2,
            b'T' | b't' => 3,
            _ => {
                return Err(PyValueError::new_err(format!(
                    "{label}: only canonical A/C/G/T are allowed, got byte {} (0x{:02X})",
                    base, base
                )));
            }
        };
        if !weight.is_finite() || *weight < 0.0 {
            return Err(PyValueError::new_err(format!(
                "{label}: weights must be finite and non-negative, got {}",
                weight
            )));
        }
        weights[idx] += weight;
        seen[idx] = true;
    }
    let total: f64 = weights.iter().sum();
    if total <= 0.0 {
        return Err(PyValueError::new_err(format!(
            "{label}: at least one weight must be strictly positive",
        )));
    }
    // Markov rows must cover every from→to entry (or the
    // generator silently treats the missing entries as zero
    // weight, which would mask an authoring bug). The
    // Python spec layer already enforces this for the typed
    // path; the defence-in-depth here is for callers that
    // bypass the DSL.
    if !seen.iter().all(|&b| b) {
        let missing: Vec<&str> = ["A", "C", "G", "T"]
            .iter()
            .enumerate()
            .filter_map(|(i, l)| if !seen[i] { Some(*l) } else { None })
            .collect();
        return Err(PyValueError::new_err(format!(
            "{label}: missing canonical base(s) {:?} — the row must cover the full A/C/G/T alphabet",
            missing
        )));
    }
    Ok(weights)
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

/// Parse a P-end label (`"V_3"`, `"D_5"`, `"D_3"`, `"J_5"`)
/// into the typed [`crate::address::PEnd`]. Same on-disk
/// spelling as the trace addresses, so the labels round-trip.
fn parse_p_end(s: &str) -> PyResult<crate::address::PEnd> {
    use crate::address::PEnd;
    match s.to_ascii_uppercase().as_str() {
        "V_3" => Ok(PEnd::V3),
        "D_5" => Ok(PEnd::D5),
        "D_3" => Ok(PEnd::D3),
        "J_5" => Ok(PEnd::J5),
        other => Err(PyValueError::new_err(format!(
            "p_addition end must be 'V_3', 'D_5', 'D_3', or 'J_5' (got {:?})",
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
