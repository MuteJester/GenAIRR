//! `SampleBasePass` — the sampling reference pass.

use crate::dist::Distribution;
use crate::ir::{NucFlags, Nucleotide, Segment, Simulation};
use crate::pass::{Pass, PassContext, PassEffect};
use crate::trace::ChoiceValue;

/// A sampling pass that draws one base byte from a `Distribution`,
/// records the choice to the trace at a configured address, and
/// appends a synthetic nucleotide carrying that base.
///
/// `SampleBasePass` is the reference implementation for any pass
/// that consumes RNG and must record its choice for the addressed-
/// trace contract (D3 + D9). Future biology-specific passes — N-nuc
/// generation, S5F substitution, PCR error injection — all follow
/// this shape: sample → record → mutate IR.
///
/// **Address discipline:** the address string is stored on the pass
/// instance and passed by reference into `trace.record`. For passes
/// that emit multiple choices per execution (e.g., a future
/// `GenerateNPBases` pass that draws N bases), the address would be
/// constructed per-draw via `format!("...[{}]", i)`. For
/// `SampleBasePass` (one draw per execution) the address is fixed.
pub struct SampleBasePass {
    address: String,
    distribution: Box<dyn Distribution<Output = u8>>,
    segment: Segment,
    flags: NucFlags,
}

impl SampleBasePass {
    /// Construct a sampling pass that draws from `distribution` and
    /// records to `address` in the trace. The pushed nucleotide will
    /// carry the given `segment` and `flags` (typically `flag::N_NUC`
    /// for TdT-like samples, `flag::P_NUC` for P-nucleotides, or
    /// `NucFlags::empty()` for plain test-only synthetic bases).
    pub fn new(
        address: impl Into<String>,
        distribution: Box<dyn Distribution<Output = u8>>,
        segment: Segment,
        flags: NucFlags,
    ) -> Self {
        Self {
            address: address.into(),
            distribution,
            segment,
            flags,
        }
    }

    /// The configured trace address.
    pub fn address(&self) -> &str {
        &self.address
    }

    /// The segment the pushed nucleotide will be tagged with.
    pub fn segment(&self) -> Segment {
        self.segment
    }
}

impl Pass for SampleBasePass {
    fn name(&self) -> &str {
        "sample_base"
    }

    fn execute(&self, sim: &Simulation, ctx: &mut PassContext) -> Simulation {
        let base = self.distribution.sample(ctx.rng);
        ctx.trace.record(&self.address[..], ChoiceValue::Base(base));
        let (next, _h) =
            sim.with_nucleotide_pushed(Nucleotide::synthetic(base, self.segment, self.flags));
        next
    }

    fn declared_choices(&self) -> Vec<String> {
        // SampleBasePass makes exactly one draw per execution at its
        // configured address. Phase D's upstream-bound propagation
        // and build-time validator both consume this list.
        vec![self.address.clone()]
    }

    fn effects(&self) -> Vec<PassEffect> {
        vec![PassEffect::AppendNucleotides]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dist::UniformBase;
    use crate::ir::{flag, NucHandle};
    use crate::pass::{PassPlan, PassRuntime};
    use crate::passes::EchoPass;

    #[test]
    fn sample_base_pass_appends_one_synthetic_nucleotide() {
        let pass = SampleBasePass::new(
            "test.np1.bases[0]",
            Box::new(UniformBase),
            Segment::Np1,
            flag::N_NUC,
        );
        let mut plan = PassPlan::new();
        plan.push(Box::new(pass));

        let outcome = PassRuntime::execute(&plan, Simulation::new(), 7);

        let final_sim = outcome.final_simulation();
        assert_eq!(final_sim.pool.len(), 1);
        let n = final_sim.pool.get(NucHandle::new(0)).unwrap();

        // Synthetic constructor: no germline provenance.
        assert_eq!(n.germline_pos, Nucleotide::NO_GERMLINE_POS);
        assert_eq!(n.segment, Segment::Np1);
        assert!(n.flags.contains(flag::N_NUC));
        // Base is one of A/C/G/T (UniformBase).
        assert!(matches!(n.base, b'A' | b'C' | b'G' | b'T'));
    }

    #[test]
    fn sample_base_pass_records_choice_at_configured_address() {
        let pass = SampleBasePass::new(
            "test.np1.bases[0]",
            Box::new(UniformBase),
            Segment::Np1,
            flag::N_NUC,
        );
        let mut plan = PassPlan::new();
        plan.push(Box::new(pass));

        let outcome = PassRuntime::execute(&plan, Simulation::new(), 7);

        // Trace contains exactly one entry, at the configured address.
        assert_eq!(outcome.trace.len(), 1);
        let rec = outcome.trace.find("test.np1.bases[0]").unwrap();
        match &rec.value {
            ChoiceValue::Base(b) => {
                assert!(matches!(*b, b'A' | b'C' | b'G' | b'T'));
            }
            other => panic!("expected ChoiceValue::Base, got {:?}", other),
        }
    }

    #[test]
    fn sample_base_pass_declares_its_address() {
        // Sampling passes override declared_choices() to report the
        // address(es) they will draw at. The runtime + Phase D
        // validator depend on this introspection.
        let pass = SampleBasePass::new(
            "test.np1.bases[0]",
            Box::new(UniformBase),
            Segment::Np1,
            flag::N_NUC,
        );
        assert_eq!(
            pass.declared_choices(),
            vec!["test.np1.bases[0]".to_string()]
        );
    }

    #[test]
    fn sample_base_pass_recorded_value_matches_pushed_nucleotide() {
        // Faithfulness: the trace must record exactly what was written
        // to the IR. Replay against a divergent record would produce
        // a different sequence; this test pins the invariant.
        let mut plan = PassPlan::new();
        for i in 0..10 {
            plan.push(Box::new(SampleBasePass::new(
                format!("test.bases[{}]", i),
                Box::new(UniformBase),
                Segment::Np1,
                flag::N_NUC,
            )));
        }

        let outcome = PassRuntime::execute(&plan, Simulation::new(), 0xdead_beef);
        let final_sim = outcome.final_simulation();

        for i in 0..10 {
            let addr = format!("test.bases[{}]", i);
            let recorded = match outcome.trace.find(&addr).unwrap().value {
                ChoiceValue::Base(b) => b,
                _ => panic!("wrong variant"),
            };
            let pushed = final_sim.pool.get(NucHandle::new(i)).unwrap().base;
            assert_eq!(
                recorded, pushed,
                "trace entry at {} = {} but IR has {} at handle {}",
                addr, recorded as char, pushed as char, i
            );
        }
    }

    // ── Phase B milestone: replay determinism ──────────────────────

    /// Build a representative mixed plan: alternating EchoPass (no
    /// RNG) and SampleBasePass (RNG-consuming). Used by the replay
    /// determinism test below.
    fn mixed_plan() -> PassPlan {
        let mut plan = PassPlan::new();
        for i in 0..8 {
            plan.push(Box::new(EchoPass::new(b'A', i as u16, Segment::V)));
            plan.push(Box::new(SampleBasePass::new(
                format!("np.np1.bases[{}]", i),
                Box::new(UniformBase),
                Segment::Np1,
                flag::N_NUC,
            )));
        }
        plan
    }

    #[test]
    fn replay_determinism_same_seed_same_trace_same_ir() {
        // The Phase B success criterion: two independent runs of the
        // same plan with the same seed must produce identical traces
        // and identical final IR revisions, byte for byte.
        let oa = PassRuntime::execute(&mixed_plan(), Simulation::new(), 0xc0ff_ee);
        let ob = PassRuntime::execute(&mixed_plan(), Simulation::new(), 0xc0ff_ee);

        // Trace equality: same addresses, same values, same order.
        assert_eq!(oa.trace.len(), ob.trace.len());
        assert_eq!(oa.trace.choices().len(), 8); // 8 SampleBase + 8 Echo (no record)
        for (a, b) in oa.trace.choices().iter().zip(ob.trace.choices().iter()) {
            assert_eq!(a.address, b.address);
            assert_eq!(a.value, b.value);
        }

        // Final IR equality: same pool size, identical bases at every
        // handle.
        let pa = &oa.final_simulation().pool;
        let pb = &ob.final_simulation().pool;
        assert_eq!(pa.len(), pb.len());
        for i in 0..pa.len() {
            let h = NucHandle::new(i as u32);
            let na = pa.get(h).unwrap();
            let nb = pb.get(h).unwrap();
            assert_eq!(na.base, nb.base);
            assert_eq!(na.germline, nb.germline);
            assert_eq!(na.germline_pos, nb.germline_pos);
            assert_eq!(na.segment, nb.segment);
            assert_eq!(na.flags, nb.flags);
        }

        // Pass-name list also identical.
        assert_eq!(oa.pass_names, ob.pass_names);
    }

    #[test]
    fn replay_determinism_different_seed_diverges() {
        let oa = PassRuntime::execute(&mixed_plan(), Simulation::new(), 1);
        let ob = PassRuntime::execute(&mixed_plan(), Simulation::new(), 2);

        // Trace addresses should still match (plan is identical).
        for (a, b) in oa.trace.choices().iter().zip(ob.trace.choices().iter()) {
            assert_eq!(a.address, b.address);
        }

        // At least one sampled base must differ between the two seeds —
        // 8 independent draws have ~vanishing probability of full agreement.
        let any_diff = oa
            .trace
            .choices()
            .iter()
            .zip(ob.trace.choices().iter())
            .any(|(a, b)| a.value != b.value);
        assert!(
            any_diff,
            "different seeds produced byte-identical traces — RNG plumbing broken"
        );
    }
}
