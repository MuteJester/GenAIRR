//! `EchoPass` — the deterministic transform reference pass.

use crate::ir::{Nucleotide, Segment, Simulation};
use crate::pass::{Pass, PassContext};

/// A deterministic transform pass that appends one configured
/// nucleotide to the simulation's pool when executed.
///
/// `EchoPass` makes no random choices — it consumes no RNG words
/// and writes nothing to the trace. The IR revision it produces
/// differs from the input only by one extra nucleotide at the end
/// of the pool.
///
/// Useful for:
/// - Building synthetic test plans where the input is fully
///   determined.
/// - Anchoring the "transform pass" pattern as a code reference
///   for future deterministic passes.
/// - Smoke-testing the runtime end to end without involving the
///   RNG or the contract system.
#[derive(Clone, Debug)]
pub struct EchoPass {
    base: u8,
    germline_pos: u16,
    segment: Segment,
}

impl EchoPass {
    /// Construct an `EchoPass` that, when executed, pushes a single
    /// `Nucleotide::germline(base, germline_pos, segment)` onto the
    /// simulation pool.
    pub fn new(base: u8, germline_pos: u16, segment: Segment) -> Self {
        Self {
            base,
            germline_pos,
            segment,
        }
    }

    /// Inspect the configured base.
    pub fn base(&self) -> u8 {
        self.base
    }

    /// Inspect the configured germline position.
    pub fn germline_pos(&self) -> u16 {
        self.germline_pos
    }

    /// Inspect the configured segment.
    pub fn segment(&self) -> Segment {
        self.segment
    }
}

impl Pass for EchoPass {
    fn name(&self) -> &str {
        "echo"
    }

    fn execute(&self, sim: &Simulation, _ctx: &mut PassContext) -> Simulation {
        let (next, _h) = sim.with_nucleotide_pushed(Nucleotide::germline(
            self.base,
            self.germline_pos,
            self.segment,
        ));
        next
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ir::NucHandle;
    use crate::pass::{PassPlan, PassRuntime};

    #[test]
    fn echo_pass_appends_configured_nucleotide() {
        let plan = {
            let mut p = PassPlan::new();
            p.push(Box::new(EchoPass::new(b'C', 42, Segment::V)));
            p
        };
        let outcome = PassRuntime::execute(&plan, Simulation::new(), 0);

        assert_eq!(outcome.revisions.len(), 2);
        let final_sim = outcome.final_simulation();
        assert_eq!(final_sim.pool.len(), 1);

        let n = final_sim.pool.get(NucHandle::new(0)).unwrap();
        assert_eq!(n.base, b'C');
        assert_eq!(n.germline, b'C'); // germline constructor: base == germline
        assert_eq!(n.germline_pos, 42);
        assert_eq!(n.segment, Segment::V);
    }

    #[test]
    fn echo_pass_records_nothing_to_trace() {
        // Determinism contract for transform passes: no RNG draws,
        // no trace writes. Verifies the design distinction between
        // transform and sampling passes.
        let plan = {
            let mut p = PassPlan::new();
            p.push(Box::new(EchoPass::new(b'A', 0, Segment::V)));
            p.push(Box::new(EchoPass::new(b'C', 1, Segment::V)));
            p.push(Box::new(EchoPass::new(b'G', 2, Segment::V)));
            p
        };
        let outcome = PassRuntime::execute(&plan, Simulation::new(), 0);

        assert!(
            outcome.trace.is_empty(),
            "EchoPass produced trace entries: {:?}",
            outcome.trace.choices()
        );
    }

    #[test]
    fn echo_pass_is_deterministic_under_any_seed() {
        // Two runs with different seeds should produce identical
        // outcomes when only EchoPass instances are in the plan,
        // because EchoPass doesn't consume the RNG.
        let plan = || {
            let mut p = PassPlan::new();
            p.push(Box::new(EchoPass::new(b'T', 5, Segment::J)));
            p.push(Box::new(EchoPass::new(b'A', 6, Segment::J)));
            p
        };

        let o1 = PassRuntime::execute(&plan(), Simulation::new(), 1);
        let o2 = PassRuntime::execute(&plan(), Simulation::new(), 0xdead_beef);

        assert_eq!(o1.revisions.len(), o2.revisions.len());
        for i in 0..o1.revisions[2].pool.len() {
            let h = NucHandle::new(i as u32);
            assert_eq!(
                o1.revisions[2].pool.get(h).unwrap().base,
                o2.revisions[2].pool.get(h).unwrap().base
            );
        }
    }

    #[test]
    fn echo_pass_composes_with_pass_plan_in_order() {
        let plan = {
            let mut p = PassPlan::new();
            p.push(Box::new(EchoPass::new(b'A', 0, Segment::V)));
            p.push(Box::new(EchoPass::new(b'C', 1, Segment::V)));
            p.push(Box::new(EchoPass::new(b'G', 2, Segment::V)));
            p.push(Box::new(EchoPass::new(b'T', 3, Segment::V)));
            p
        };
        let outcome = PassRuntime::execute(&plan, Simulation::new(), 0);

        // 5 revisions: initial + 4 echo passes.
        assert_eq!(outcome.revisions.len(), 5);
        assert_eq!(outcome.pass_names, vec!["echo"; 4]);

        // Each successive revision has one more nucleotide.
        for (i, rev) in outcome.revisions.iter().enumerate() {
            assert_eq!(rev.pool.len(), i);
        }

        // Bases pushed in order ACGT.
        let final_pool = &outcome.final_simulation().pool;
        for (i, &expected) in [b'A', b'C', b'G', b'T'].iter().enumerate() {
            assert_eq!(
                final_pool.get(NucHandle::new(i as u32)).unwrap().base,
                expected
            );
        }
    }

    #[test]
    fn echo_pass_accessors_round_trip() {
        let p = EchoPass::new(b'G', 17, Segment::Np1);
        assert_eq!(p.base(), b'G');
        assert_eq!(p.germline_pos(), 17);
        assert_eq!(p.segment(), Segment::Np1);
        assert_eq!(p.name(), "echo");
    }

    #[test]
    fn echo_pass_declares_no_choices() {
        // Transform passes (no RNG draws) inherit the default empty
        // declared_choices(). Phase D's validator + upstream-bound
        // propagation will skip them by virtue of getting an empty
        // address list.
        let p = EchoPass::new(b'A', 0, Segment::V);
        assert!(p.declared_choices().is_empty());
    }

    #[test]
    fn echo_pass_preserves_history_chain_under_persistent_ir() {
        // D1 contract at the integration level: every revision retains
        // its own pool size after subsequent passes have run.
        let plan = {
            let mut p = PassPlan::new();
            p.push(Box::new(EchoPass::new(b'A', 0, Segment::V)));
            p.push(Box::new(EchoPass::new(b'C', 1, Segment::V)));
            p.push(Box::new(EchoPass::new(b'G', 2, Segment::V)));
            p
        };
        let outcome = PassRuntime::execute(&plan, Simulation::new(), 0);

        // Initial revision: empty pool, even after 3 passes ran.
        assert_eq!(outcome.revisions[0].pool.len(), 0);
        // After first echo: 1 nucleotide, never grows.
        assert_eq!(outcome.revisions[1].pool.len(), 1);
        assert_eq!(
            outcome.revisions[1]
                .pool
                .get(NucHandle::new(0))
                .unwrap()
                .base,
            b'A'
        );
        // After second echo: 2 nucleotides; the original 'A' is unchanged.
        assert_eq!(outcome.revisions[2].pool.len(), 2);
        assert_eq!(
            outcome.revisions[2]
                .pool
                .get(NucHandle::new(0))
                .unwrap()
                .base,
            b'A'
        );
        assert_eq!(
            outcome.revisions[2]
                .pool
                .get(NucHandle::new(1))
                .unwrap()
                .base,
            b'C'
        );
    }
}
