//! `GenerateNPPass` — TdT-like N-nucleotide region generation (C.7).

use crate::contract::ChoiceContext;
use crate::dist::{sample_filtered_result, Distribution, FilteredSampleError};
use crate::ir::{flag, NucHandle, Nucleotide, Region, Segment, Simulation};
use crate::pass::{IntegerSupport, Pass, PassCompileFact, PassContext, PassEffect, PassError};
use crate::trace::ChoiceValue;

/// Generate one NP (non-template) region — a stretch of bases
/// interpolated between adjacent V/D/J segments by the TdT enzyme
/// during recombination.
///
/// Single pass produces single IR revision (D2 boundary), but
/// internally records *each* atomic random choice to the trace at
/// its own indexed address. Length, then N bases. The trace stays
/// faithful at the choice level even though the IR revision count
/// stays at the biological-event level.
///
/// Parameterized by:
/// - **np_segment** — `Segment::Np1` or `Segment::Np2`. Other
///   segments rejected at construction.
/// - **length_dist** — `Box<dyn Distribution<Output = i64>>`.
///   Empirical NP length distribution. Negative or oversized
///   lengths are caller bugs and panic at execute time.
/// - **base_dist** — `Box<dyn Distribution<Output = u8>>`.
///   Per-base distribution (typically `UniformBase` until empirical
///   TdT models arrive in Phase E). When contracts are active and the
///   distribution exposes finite support, each base draw is filtered
///   through `ContractSet::admits` before being recorded.
///
/// On execute:
/// 1. Sample length L from `length_dist`. Validate `0 <= L <= u32::MAX`.
/// 2. Record `ChoiceValue::Int(L)` at `"np.{np1|np2}.length"`.
/// 3. For each of L positions: sample a base, record at
///    `"np.{np1|np2}.bases[i]"`, push a synthetic nucleotide
///    (`Nucleotide::synthetic(base, np_segment, N_NUC)`) onto the
///    pool.
/// 4. Construct `Region(np_segment, [start, start+L))` with
///    `frame_phase = 0` (cross-region frame chaining is C.8's
///    responsibility) and recompute its codon rail against the
///    new pool.
/// 5. Append the region to the sequence.
pub struct GenerateNPPass {
    np_segment: Segment,
    length_dist: Box<dyn Distribution<Output = i64>>,
    base_dist: Box<dyn Distribution<Output = u8>>,
}

impl GenerateNPPass {
    /// Construct an NP-generation pass.
    ///
    /// Panics if `np_segment` is not `Segment::Np1` or
    /// `Segment::Np2`.
    pub fn new(
        np_segment: Segment,
        length_dist: Box<dyn Distribution<Output = i64>>,
        base_dist: Box<dyn Distribution<Output = u8>>,
    ) -> Self {
        match np_segment {
            Segment::Np1 | Segment::Np2 => {}
            _ => panic!(
                "GenerateNPPass: np_segment must be Np1 or Np2 — got {:?}",
                np_segment
            ),
        }
        Self {
            np_segment,
            length_dist,
            base_dist,
        }
    }

    pub fn np_segment(&self) -> Segment {
        self.np_segment
    }

    fn length_address(&self) -> &'static str {
        match self.np_segment {
            Segment::Np1 => "np.np1.length",
            Segment::Np2 => "np.np2.length",
            _ => unreachable!("GenerateNPPass with non-NP segment"),
        }
    }

    fn bases_prefix(&self) -> &'static str {
        match self.np_segment {
            Segment::Np1 => "np.np1.bases",
            Segment::Np2 => "np.np2.bases",
            _ => unreachable!("GenerateNPPass with non-NP segment"),
        }
    }

    fn pass_name(&self) -> &'static str {
        match self.np_segment {
            Segment::Np1 => "generate_np.np1",
            Segment::Np2 => "generate_np.np2",
            _ => unreachable!("GenerateNPPass with non-NP segment"),
        }
    }

    /// Constraint-aware length sample (D.6).
    ///
    /// When `ctx.contracts` is `Some`, we use `sample_filtered_result`
    /// to draw only from values the contract set admits at the
    /// length address. In permissive mode, falls back to plain
    /// `sample` when:
    /// - no contracts are active,
    /// - the distribution doesn't expose `support()` (continuous
    ///   or too-large categorical), or
    /// - the filtered support is empty (no admissible candidate).
    ///
    /// In strict mode, the last two cases return a structured
    /// `PassError` instead of silently sampling outside the contract.
    fn sample_length(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
        address: &'static str,
        strict: bool,
    ) -> Result<i64, PassError> {
        // Capture borrows of fields we'll read in the closure so
        // the borrow checker doesn't see overlapping `&mut ctx`
        // and `&ctx.contracts`.
        let refdata = ctx.refdata;
        let contracts = ctx.contracts;

        if let Some(contracts) = contracts {
            match sample_filtered_result(ctx.rng, self.length_dist.as_ref(), |&candidate: &i64| {
                contracts
                    .admits(sim, refdata, address, &ChoiceValue::Int(candidate))
                    .is_ok()
            }) {
                Ok(value) => return Ok(value),
                Err(reason) if strict => {
                    return Err(self.constraint_sampling_error(address, reason));
                }
                Err(_) => {
                    // Permissive legacy path: empty/non-enumerable filters
                    // fall through to unconstrained sampling.
                }
            }
        }
        Ok(self.length_dist.sample(ctx.rng))
    }

    /// Constraint-aware base sample (D.6).
    ///
    /// This is the base-level counterpart to `sample_length`: contracts
    /// can reject a candidate base before it is written to the IR. The
    /// immediate hardening target is `NoStopCodonInJunction`, which
    /// rejects NP bases that would complete a junction stop codon.
    fn sample_base(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
        address: &str,
        index: u32,
        total_len: u32,
        strict: bool,
    ) -> Result<u8, PassError> {
        let refdata = ctx.refdata;
        let contracts = ctx.contracts;

        if let Some(contracts) = contracts {
            let context = ChoiceContext::indexed(index, total_len);
            match sample_filtered_result(ctx.rng, self.base_dist.as_ref(), |candidate: &u8| {
                contracts
                    .admits_with_context(
                        sim,
                        refdata,
                        address,
                        &ChoiceValue::Base(*candidate),
                        context,
                    )
                    .is_ok()
            }) {
                Ok(value) => return Ok(value),
                Err(reason) if strict => {
                    return Err(self.constraint_sampling_error(address, reason));
                }
                Err(_) => {}
            }
        }

        Ok(self.base_dist.sample(ctx.rng))
    }

    fn constraint_sampling_error(
        &self,
        address: impl Into<String>,
        reason: FilteredSampleError,
    ) -> PassError {
        PassError::constraint_sampling(self.pass_name(), address, reason)
    }

    fn execute_with_sampling_mode(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
        strict: bool,
    ) -> Result<Simulation, PassError> {
        // 1. Sample length, possibly constrained by active contracts.
        //
        // **Constraint-aware sampling (D.6):** if `ctx.contracts`
        // is set, we filter the length distribution's support to
        // values that the contract set admits at this address. In
        // strict mode, inability to draw an admissible value is a
        // structured error. In permissive mode, we preserve the
        // original Phase D.6 fallback behavior.
        let address = self.length_address();
        let length = self.sample_length(sim, ctx, address, strict)?;
        if strict && length < 0 {
            return Err(PassError::invalid_distribution_output(
                self.pass_name(),
                address,
                length,
                "negative_length",
            ));
        }
        if strict && length > u32::MAX as i64 {
            return Err(PassError::invalid_distribution_output(
                self.pass_name(),
                address,
                length,
                "length_exceeds_u32",
            ));
        }
        assert!(
            length >= 0,
            "GenerateNPPass({}): length distribution returned negative value {}",
            address,
            length
        );
        assert!(
            length <= u32::MAX as i64,
            "GenerateNPPass({}): length distribution returned {} > u32::MAX",
            address,
            length
        );
        ctx.trace.record(address, ChoiceValue::Int(length));

        let length = length as u32;
        let region_start = NucHandle::new(sim.pool.len() as u32);

        // 2. Sample bases and push them into the pool. Each push
        //    produces a one-step persistent IR update; the loop
        //    accumulates into `current`.
        let mut current = sim.clone();
        for i in 0..length {
            let base_address = format!("{}[{}]", self.bases_prefix(), i);
            let base = self.sample_base(&current, ctx, &base_address, i, length, strict)?;
            ctx.trace.record(base_address, ChoiceValue::Base(base));
            let (next, _h) = current.with_nucleotide_pushed(Nucleotide::synthetic(
                base,
                self.np_segment,
                flag::N_NUC,
            ));
            current = next;
        }

        // 3. Build the region and append it. Codon rail recomputes
        //    against the updated pool. Frame phase chains from
        //    cumulative prior region length, same convention as
        //    AssembleSegmentPass — keeps the codon frame coherent
        //    across V/NP/D/NP/J boundaries.
        let cumulative_len: u32 = sim.sequence.regions.iter().map(|r| r.len()).sum();
        let frame_phase = (cumulative_len % 3) as u8;

        let region_end = NucHandle::new(current.pool.len() as u32);
        let region = Region::new(self.np_segment, region_start, region_end)
            .with_frame_phase(frame_phase)
            .with_codon_rail_recomputed(&current.pool);
        Ok(current.with_region_added(region))
    }
}

impl Pass for GenerateNPPass {
    fn name(&self) -> &str {
        self.pass_name()
    }

    fn execute(&self, sim: &Simulation, ctx: &mut PassContext) -> Simulation {
        self.execute_with_sampling_mode(sim, ctx, false)
            .expect("GenerateNPPass permissive execution must not return PassError")
    }

    fn execute_checked(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
    ) -> Result<Simulation, PassError> {
        self.execute_with_sampling_mode(sim, ctx, true)
    }

    fn declared_choices(&self) -> Vec<String> {
        // Length is fixed-address. Bases are variable-count; per
        // D3 we use the literal `[0..n]` form to indicate runtime
        // expansion.
        vec![
            self.length_address().to_string(),
            format!("{}[0..n]", self.bases_prefix()),
        ]
    }

    fn effects(&self) -> Vec<PassEffect> {
        vec![PassEffect::AppendRegion(self.np_segment)]
    }

    fn compile_facts(&self) -> Vec<PassCompileFact> {
        vec![PassCompileFact::NpLengthSupport {
            segment: self.np_segment,
            support: IntegerSupport::from_weighted_pairs(self.length_dist.support()),
        }]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dist::{EmpiricalLengthDist, UniformBase};
    use crate::pass::{PassError, PassPlan, PassRuntime};
    use crate::passes::EchoPass;

    #[test]
    #[should_panic(expected = "np_segment must be Np1 or Np2")]
    fn generate_np_pass_rejects_v_segment() {
        let _ = GenerateNPPass::new(
            Segment::V,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)])),
            Box::new(UniformBase),
        );
    }

    #[test]
    #[should_panic(expected = "np_segment must be Np1 or Np2")]
    fn generate_np_pass_rejects_d_segment() {
        let _ = GenerateNPPass::new(
            Segment::D,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)])),
            Box::new(UniformBase),
        );
    }

    #[test]
    fn generate_np_pass_zero_length_creates_empty_region() {
        let mut plan = PassPlan::new();
        plan.push(Box::new(GenerateNPPass::new(
            Segment::Np1,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(0, 1.0)])),
            Box::new(UniformBase),
        )));
        let outcome = PassRuntime::execute(&plan, Simulation::new(), 0);

        let sim = outcome.final_simulation();
        // Pool is unchanged.
        assert_eq!(sim.pool.len(), 0);
        // One region was added (NP1, empty).
        assert_eq!(sim.sequence.region_count(), 1);
        let r = &sim.sequence.regions[0];
        assert_eq!(r.segment, Segment::Np1);
        assert!(r.is_empty());
        // Trace recorded length, no bases.
        assert_eq!(outcome.trace.len(), 1);
        assert_eq!(
            outcome.trace.find("np.np1.length").unwrap().value,
            ChoiceValue::Int(0)
        );
    }

    #[test]
    fn generate_np_pass_pushes_n_bases_with_correct_metadata() {
        let mut plan = PassPlan::new();
        plan.push(Box::new(GenerateNPPass::new(
            Segment::Np1,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(5, 1.0)])),
            Box::new(UniformBase),
        )));
        let outcome = PassRuntime::execute(&plan, Simulation::new(), 1234);

        let sim = outcome.final_simulation();
        assert_eq!(sim.pool.len(), 5);
        assert_eq!(sim.sequence.region_count(), 1);

        // Region covers all 5 nucleotides.
        let r = &sim.sequence.regions[0];
        assert_eq!(r.segment, Segment::Np1);
        assert_eq!(r.start, NucHandle::new(0));
        assert_eq!(r.end, NucHandle::new(5));

        // Each pushed nucleotide is a synthetic NP1 N-nuc.
        for i in 0..5 {
            let n = sim.pool.get(NucHandle::new(i)).unwrap();
            assert_eq!(n.segment, Segment::Np1);
            assert_eq!(n.germline_pos, Nucleotide::NO_GERMLINE_POS);
            assert!(n.flags.contains(flag::N_NUC));
            assert!(matches!(n.base, b'A' | b'C' | b'G' | b'T'));
        }
    }

    #[test]
    fn generate_np_pass_records_length_and_indexed_bases() {
        let mut plan = PassPlan::new();
        plan.push(Box::new(GenerateNPPass::new(
            Segment::Np2,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)])),
            Box::new(UniformBase),
        )));
        let outcome = PassRuntime::execute(&plan, Simulation::new(), 7);

        // Length entry + 3 base entries = 4 trace records.
        assert_eq!(outcome.trace.len(), 4);
        assert_eq!(
            outcome.trace.find("np.np2.length").unwrap().value,
            ChoiceValue::Int(3)
        );
        for i in 0..3 {
            let addr = format!("np.np2.bases[{}]", i);
            let rec = outcome
                .trace
                .find(&addr)
                .unwrap_or_else(|| panic!("missing {}", addr));
            assert!(matches!(rec.value, ChoiceValue::Base(_)));
        }
    }

    #[test]
    fn generate_np_pass_trace_value_matches_pushed_base_at_each_position() {
        // Faithfulness: trace[i] equals pool[i].base for each NP base.
        let mut plan = PassPlan::new();
        plan.push(Box::new(GenerateNPPass::new(
            Segment::Np1,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(8, 1.0)])),
            Box::new(UniformBase),
        )));
        let outcome = PassRuntime::execute(&plan, Simulation::new(), 99);
        let sim = outcome.final_simulation();

        for i in 0..8 {
            let addr = format!("np.np1.bases[{}]", i);
            let recorded = match outcome.trace.find(&addr).unwrap().value {
                ChoiceValue::Base(b) => b,
                _ => panic!("wrong variant"),
            };
            let pushed = sim.pool.get(NucHandle::new(i)).unwrap().base;
            assert_eq!(
                recorded, pushed,
                "trace at {} = {} but pool[{}] = {}",
                addr, recorded as char, i, pushed as char
            );
        }
    }

    #[test]
    fn generate_np_pass_is_deterministic_under_same_seed() {
        let plan = || {
            let mut p = PassPlan::new();
            p.push(Box::new(GenerateNPPass::new(
                Segment::Np1,
                Box::new(EmpiricalLengthDist::from_pairs(vec![
                    (1, 1.0),
                    (3, 2.0),
                    (5, 1.0),
                ])),
                Box::new(UniformBase),
            )));
            p
        };
        let oa = PassRuntime::execute(&plan(), Simulation::new(), 0xc0ff_ee);
        let ob = PassRuntime::execute(&plan(), Simulation::new(), 0xc0ff_ee);

        assert_eq!(
            oa.final_simulation().pool.len(),
            ob.final_simulation().pool.len()
        );
        assert_eq!(oa.trace.choices(), ob.trace.choices());
    }

    #[test]
    fn generate_np_pass_declared_choices_lists_length_and_bases_pattern() {
        let pass = GenerateNPPass::new(
            Segment::Np1,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)])),
            Box::new(UniformBase),
        );
        let declared = pass.declared_choices();
        assert_eq!(declared.len(), 2);
        assert!(declared.contains(&"np.np1.length".to_string()));
        assert!(declared.contains(&"np.np1.bases[0..n]".to_string()));
    }

    #[test]
    fn generate_np_pass_appends_region_after_existing_pool_state() {
        // Pre-populate the pool with some nucleotides, then run NP
        // generation. The new region should start at the current
        // pool boundary, not at handle 0.
        let mut plan = PassPlan::new();
        // Three Echo passes push 3 nucleotides as germline.
        plan.push(Box::new(EchoPass::new(b'A', 0, Segment::V)));
        plan.push(Box::new(EchoPass::new(b'C', 1, Segment::V)));
        plan.push(Box::new(EchoPass::new(b'G', 2, Segment::V)));
        // Then NP1 of length 4.
        plan.push(Box::new(GenerateNPPass::new(
            Segment::Np1,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(4, 1.0)])),
            Box::new(UniformBase),
        )));

        let outcome = PassRuntime::execute(&plan, Simulation::new(), 5);
        let sim = outcome.final_simulation();

        // Total pool: 3 V + 4 NP1 = 7.
        assert_eq!(sim.pool.len(), 7);
        // One region (NP1) — Echo doesn't create regions.
        assert_eq!(sim.sequence.region_count(), 1);
        let r = &sim.sequence.regions[0];
        assert_eq!(r.start, NucHandle::new(3));
        assert_eq!(r.end, NucHandle::new(7));
        assert_eq!(r.segment, Segment::Np1);
    }

    #[test]
    #[should_panic(expected = "negative value")]
    fn generate_np_pass_panics_on_negative_length() {
        use crate::dist::UniformInt;
        let mut plan = PassPlan::new();
        plan.push(Box::new(GenerateNPPass::new(
            Segment::Np1,
            Box::new(UniformInt::new(-3, -1)),
            Box::new(UniformBase),
        )));
        let _ = PassRuntime::execute(&plan, Simulation::new(), 0);
    }

    #[test]
    fn generate_np_pass_strict_errors_on_negative_length() {
        use crate::dist::UniformInt;
        let mut plan = PassPlan::new();
        plan.push(Box::new(GenerateNPPass::new(
            Segment::Np1,
            Box::new(UniformInt::new(-3, -1)),
            Box::new(UniformBase),
        )));

        let err = PassRuntime::execute_strict_with_context(&plan, Simulation::new(), 0, None, None)
            .unwrap_err();

        assert_eq!(err.pass_name(), "generate_np.np1");
        assert_eq!(err.address(), "np.np1.length");
        assert!(matches!(
            err,
            PassError::InvalidDistributionOutput { value, .. } if value < 0
        ));
    }
}
