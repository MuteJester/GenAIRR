//! `TrimPass` — recombination-stage trim sampling (C.6).

use crate::assignment::TrimEnd;
use crate::dist::Distribution;
use crate::ir::{Segment, Simulation};
use crate::pass::{Pass, PassContext};
use crate::trace::ChoiceValue;

/// Sample a trim amount from a distribution and apply it to the
/// assigned allele on the given segment / end.
///
/// The pass is parameterized by:
/// - **Segment** — must be V, D, or J (NP and C are rejected at
///   construction).
/// - **End** — `TrimEnd::Five` or `TrimEnd::Three`. All six
///   `(segment, end)` combinations are syntactically valid; whether
///   biology uses them is up to the plan author.
/// - **Distribution** — any `Box<dyn Distribution<Output = i64>>`,
///   typically an `EmpiricalLengthDist` constructed at plan-build
///   time with the empirical trim distribution for that segment.
///
/// On execute the pass:
/// 1. Draws one `i64` trim value from the distribution via
///    `ctx.rng`.
/// 2. Validates `0 <= value <= u16::MAX` — distributions that
///    return negative or oversized trims are caller bugs and
///    panic loudly with the offending value.
/// 3. Records `ChoiceValue::Int(value)` to the trace at address
///    `"trim.{segment}_{end}"` (e.g. `"trim.v_3"`, `"trim.j_5"`).
/// 4. Returns `sim.with_trim(segment, end, value as u16)`.
///
/// **Pre-condition:** an allele must already be assigned to
/// `segment` (i.e., the matching `SampleAllelePass` ran earlier in
/// the plan). `Simulation::with_trim` panics if the slot is empty,
/// and so does this pass by extension.
pub struct TrimPass {
    segment: Segment,
    end: TrimEnd,
    distribution: Box<dyn Distribution<Output = i64>>,
}

impl TrimPass {
    /// Construct a trim pass.
    ///
    /// Panics if `segment` is anything other than V, D, or J.
    pub fn new(
        segment: Segment,
        end: TrimEnd,
        distribution: Box<dyn Distribution<Output = i64>>,
    ) -> Self {
        match segment {
            Segment::V | Segment::D | Segment::J => {}
            _ => panic!("TrimPass: segment must be V, D, or J — got {:?}", segment),
        }
        Self {
            segment,
            end,
            distribution,
        }
    }

    pub fn segment(&self) -> Segment {
        self.segment
    }

    pub fn end(&self) -> TrimEnd {
        self.end
    }

    /// The hierarchical-string address (D3) at which this pass
    /// records its choice. Same string as `name()`.
    fn address(&self) -> &'static str {
        match (self.segment, self.end) {
            (Segment::V, TrimEnd::Five) => "trim.v_5",
            (Segment::V, TrimEnd::Three) => "trim.v_3",
            (Segment::D, TrimEnd::Five) => "trim.d_5",
            (Segment::D, TrimEnd::Three) => "trim.d_3",
            (Segment::J, TrimEnd::Five) => "trim.j_5",
            (Segment::J, TrimEnd::Three) => "trim.j_3",
            _ => unreachable!("TrimPass with non-V/D/J segment"),
        }
    }
}

impl Pass for TrimPass {
    fn name(&self) -> &str {
        self.address()
    }

    fn execute(&self, sim: &Simulation, ctx: &mut PassContext) -> Simulation {
        let value = self.distribution.sample(ctx.rng);
        assert!(
            value >= 0,
            "TrimPass({}): distribution returned negative trim {}",
            self.address(),
            value
        );
        assert!(
            value <= u16::MAX as i64,
            "TrimPass({}): distribution returned trim {} > u16::MAX",
            self.address(),
            value
        );

        ctx.trace.record(self.address(), ChoiceValue::Int(value));
        sim.with_trim(self.segment, self.end, value as u16)
    }

    fn declared_choices(&self) -> Vec<String> {
        vec![self.address().to_string()]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dist::{AllelePoolDist, EmpiricalLengthDist};
    use crate::pass::{PassPlan, PassRuntime};
    use crate::passes::sample_allele::test_support::make_test_pool;
    use crate::passes::SampleAllelePass;

    #[test]
    #[should_panic(expected = "segment must be V, D, or J")]
    fn trim_pass_rejects_np1() {
        let _ = TrimPass::new(
            Segment::Np1,
            TrimEnd::Three,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(0, 1.0)])),
        );
    }

    #[test]
    fn trim_pass_addresses_cover_all_six_combinations() {
        let dist = || Box::new(EmpiricalLengthDist::from_pairs(vec![(0, 1.0)]));

        let cases = [
            (Segment::V, TrimEnd::Five, "trim.v_5"),
            (Segment::V, TrimEnd::Three, "trim.v_3"),
            (Segment::D, TrimEnd::Five, "trim.d_5"),
            (Segment::D, TrimEnd::Three, "trim.d_3"),
            (Segment::J, TrimEnd::Five, "trim.j_5"),
            (Segment::J, TrimEnd::Three, "trim.j_3"),
        ];
        for (seg, end, expected_addr) in cases {
            let pass = TrimPass::new(seg, end, dist());
            assert_eq!(pass.name(), expected_addr);
            assert_eq!(pass.declared_choices(), vec![expected_addr.to_string()]);
        }
    }

    #[test]
    fn trim_pass_applies_trim_and_records_choice() {
        // V allele assigned, then trim_3 sampled.
        let v_pool = make_test_pool(1, Segment::V);
        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(&v_pool)),
        )));
        plan.push(Box::new(TrimPass::new(
            Segment::V,
            TrimEnd::Three,
            // Single-value dist: always returns 4.
            Box::new(EmpiricalLengthDist::from_pairs(vec![(4, 1.0)])),
        )));

        let outcome = PassRuntime::execute(&plan, Simulation::new(), 7);

        // Trim recorded at canonical address.
        let rec = outcome.trace.find("trim.v_3").expect("trim recorded");
        assert_eq!(rec.value, ChoiceValue::Int(4));

        // Trim applied to the V allele instance.
        let v = outcome.final_simulation().assignments.v.unwrap();
        assert_eq!(v.trim_3, 4);
        assert_eq!(v.trim_5, 0);
    }

    #[test]
    #[should_panic(expected = "no instance assigned to segment")]
    fn trim_pass_panics_when_no_allele_assigned() {
        // The trim pass relies on `Simulation::with_trim`, which
        // panics when no instance is assigned to the segment. This
        // test pins the propagation.
        let mut plan = PassPlan::new();
        plan.push(Box::new(TrimPass::new(
            Segment::V,
            TrimEnd::Three,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(2, 1.0)])),
        )));
        let _ = PassRuntime::execute(&plan, Simulation::new(), 0);
    }

    #[test]
    #[should_panic(expected = "negative trim")]
    fn trim_pass_panics_on_negative_distribution_output() {
        // Runtime defensive check. UniformInt::new(-5, -1) always
        // returns a negative trim, which TrimPass must reject loudly.
        use crate::dist::UniformInt;
        let v_pool = make_test_pool(1, Segment::V);
        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(&v_pool)),
        )));
        plan.push(Box::new(TrimPass::new(
            Segment::V,
            TrimEnd::Three,
            Box::new(UniformInt::new(-5, -1)),
        )));
        let _ = PassRuntime::execute(&plan, Simulation::new(), 1);
    }

    #[test]
    fn trim_pass_full_vdj_chain_records_six_choices() {
        // Heavy chain shape: sample V/D/J, then trim each side.
        let v_pool = make_test_pool(1, Segment::V);
        let d_pool = make_test_pool(1, Segment::D);
        let j_pool = make_test_pool(1, Segment::J);
        let dist = || Box::new(EmpiricalLengthDist::from_pairs(vec![(2, 1.0)]));

        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(&v_pool)),
        )));
        plan.push(Box::new(SampleAllelePass::new(
            Segment::D,
            Box::new(AllelePoolDist::uniform(&d_pool)),
        )));
        plan.push(Box::new(SampleAllelePass::new(
            Segment::J,
            Box::new(AllelePoolDist::uniform(&j_pool)),
        )));
        plan.push(Box::new(TrimPass::new(Segment::V, TrimEnd::Three, dist())));
        plan.push(Box::new(TrimPass::new(Segment::D, TrimEnd::Five, dist())));
        plan.push(Box::new(TrimPass::new(Segment::D, TrimEnd::Three, dist())));
        plan.push(Box::new(TrimPass::new(Segment::J, TrimEnd::Five, dist())));

        let outcome = PassRuntime::execute(&plan, Simulation::new(), 42);

        // 3 sampling + 4 trim = 7 trace entries.
        assert_eq!(outcome.trace.len(), 7);

        // All four trim addresses present with value 2.
        for addr in ["trim.v_3", "trim.d_5", "trim.d_3", "trim.j_5"] {
            assert_eq!(outcome.trace.find(addr).unwrap().value, ChoiceValue::Int(2));
        }

        // Trims applied to the assignments.
        let sim = outcome.final_simulation();
        assert_eq!(sim.assignments.v.unwrap().trim_3, 2);
        assert_eq!(sim.assignments.d.unwrap().trim_5, 2);
        assert_eq!(sim.assignments.d.unwrap().trim_3, 2);
        assert_eq!(sim.assignments.j.unwrap().trim_5, 2);
    }

    #[test]
    fn trim_pass_is_deterministic_under_same_seed() {
        let v_pool = make_test_pool(1, Segment::V);
        let plan = || {
            let mut p = PassPlan::new();
            p.push(Box::new(SampleAllelePass::new(
                Segment::V,
                Box::new(AllelePoolDist::uniform(&v_pool)),
            )));
            p.push(Box::new(TrimPass::new(
                Segment::V,
                TrimEnd::Three,
                Box::new(EmpiricalLengthDist::from_pairs(vec![
                    (0, 1.0),
                    (1, 2.0),
                    (2, 3.0),
                    (5, 1.0),
                ])),
            )));
            p
        };

        let oa = PassRuntime::execute(&plan(), Simulation::new(), 0xc0ff_ee);
        let ob = PassRuntime::execute(&plan(), Simulation::new(), 0xc0ff_ee);
        assert_eq!(
            oa.final_simulation().assignments.v.unwrap().trim_3,
            ob.final_simulation().assignments.v.unwrap().trim_3
        );
        assert_eq!(oa.trace.choices(), ob.trace.choices());
    }
}
