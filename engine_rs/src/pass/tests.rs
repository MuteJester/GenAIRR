use super::*;
use crate::ir::{NucFlags, NucHandle, Nucleotide, Segment, Simulation};
use crate::trace::ChoiceValue;

struct TestEchoPass {
    base: u8,
    germline_pos: u16,
}

impl Pass for TestEchoPass {
    fn name(&self) -> &str {
        "test_echo"
    }

    fn execute(&self, sim: &Simulation, ctx: &mut PassContext) -> Simulation {
        let (next, _h) = sim.with_nucleotide_pushed(Nucleotide::germline(
            self.base,
            self.germline_pos,
            Segment::V,
        ));
        ctx.trace
            .record("test_echo.base", ChoiceValue::Base(self.base));
        next
    }
}

struct TestRandomBasePass;

impl Pass for TestRandomBasePass {
    fn name(&self) -> &str {
        "test_random_base"
    }

    fn execute(&self, sim: &Simulation, ctx: &mut PassContext) -> Simulation {
        const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
        let idx = ctx.rng.range_u32(4) as usize;
        let base = BASES[idx];
        let (next, _h) = sim.with_nucleotide_pushed(Nucleotide::synthetic(
            base,
            Segment::Np1,
            NucFlags::empty(),
        ));
        ctx.trace
            .record("test_random_base.base", ChoiceValue::Base(base));
        next
    }
}

#[test]
fn empty_plan_returns_initial_unchanged() {
    let initial = Simulation::new();
    let plan = PassPlan::new();
    let outcome = PassRuntime::execute(&plan, initial.clone(), 42);

    assert_eq!(outcome.revisions.len(), 1);
    assert!(outcome.pass_names.is_empty());
    assert!(outcome.trace.is_empty());
    assert_eq!(outcome.final_simulation().pool.len(), 0);
}

#[test]
fn single_pass_plan_produces_two_revisions() {
    let initial = Simulation::new();
    let mut plan = PassPlan::new();
    plan.push(Box::new(TestEchoPass {
        base: b'A',
        germline_pos: 0,
    }));

    let outcome = PassRuntime::execute(&plan, initial.clone(), 0);

    assert_eq!(outcome.revisions.len(), 2);
    assert_eq!(outcome.pass_names, vec!["test_echo".to_string()]);
    assert_eq!(outcome.revisions[0].pool.len(), 0);
    assert_eq!(outcome.revisions[1].pool.len(), 1);
    assert_eq!(
        outcome.revisions[1]
            .pool
            .get(NucHandle::new(0))
            .unwrap()
            .base,
        b'A'
    );
    assert_eq!(outcome.trace.len(), 1);
    assert_eq!(
        outcome.trace.find("test_echo.base").unwrap().value,
        ChoiceValue::Base(b'A')
    );
}

#[test]
fn multi_pass_plan_applies_in_order_and_history_is_consistent() {
    let initial = Simulation::new();
    let mut plan = PassPlan::new();
    plan.push(Box::new(TestEchoPass {
        base: b'A',
        germline_pos: 0,
    }));
    plan.push(Box::new(TestEchoPass {
        base: b'C',
        germline_pos: 1,
    }));
    plan.push(Box::new(TestEchoPass {
        base: b'G',
        germline_pos: 2,
    }));

    let outcome = PassRuntime::execute(&plan, initial, 0);

    assert_eq!(outcome.revisions.len(), 4);
    assert_eq!(outcome.pass_names.len(), 3);
    assert_eq!(outcome.revisions[0].pool.len(), 0);
    assert_eq!(outcome.revisions[1].pool.len(), 1);
    assert_eq!(outcome.revisions[2].pool.len(), 2);
    assert_eq!(outcome.revisions[3].pool.len(), 3);

    let pool = &outcome.revisions[3].pool;
    assert_eq!(pool.get(NucHandle::new(0)).unwrap().base, b'A');
    assert_eq!(pool.get(NucHandle::new(1)).unwrap().base, b'C');
    assert_eq!(pool.get(NucHandle::new(2)).unwrap().base, b'G');
    assert_eq!(outcome.revisions[1].pool.len(), 1);
    assert!(outcome.revisions[1].pool.get(NucHandle::new(1)).is_none());
}

#[test]
fn runtime_threads_rng_through_passes_deterministically() {
    let mut plan_a = PassPlan::new();
    for _ in 0..10 {
        plan_a.push(Box::new(TestRandomBasePass));
    }
    let mut plan_b = PassPlan::new();
    for _ in 0..10 {
        plan_b.push(Box::new(TestRandomBasePass));
    }

    let oa = PassRuntime::execute(&plan_a, Simulation::new(), 0xc0ff_ee);
    let ob = PassRuntime::execute(&plan_b, Simulation::new(), 0xc0ff_ee);

    assert_eq!(oa.trace.len(), 10);
    assert_eq!(ob.trace.len(), 10);
    for (a, b) in oa.trace.choices().iter().zip(ob.trace.choices().iter()) {
        assert_eq!(a.address, b.address);
        assert_eq!(a.value, b.value);
    }

    let pa = &oa.final_simulation().pool;
    let pb = &ob.final_simulation().pool;
    assert_eq!(pa.len(), pb.len());
    for i in 0..pa.len() {
        let h = NucHandle::new(i as u32);
        assert_eq!(pa.get(h).unwrap().base, pb.get(h).unwrap().base);
    }
}

#[test]
fn runtime_with_different_seeds_produces_different_traces() {
    let mut plan = PassPlan::new();
    for _ in 0..10 {
        plan.push(Box::new(TestRandomBasePass));
    }

    let o42 = PassRuntime::execute(&plan, Simulation::new(), 42);
    let o43 = PassRuntime::execute(&plan, Simulation::new(), 43);

    let any_diff = o42
        .trace
        .choices()
        .iter()
        .zip(o43.trace.choices().iter())
        .any(|(a, b)| a.value != b.value);
    assert!(
        any_diff,
        "different seeds produced identical traces — RNG plumbing broken"
    );
}

#[test]
fn outcome_revision_after_returns_correct_revision() {
    let mut plan = PassPlan::new();
    plan.push(Box::new(TestEchoPass {
        base: b'A',
        germline_pos: 0,
    }));
    plan.push(Box::new(TestEchoPass {
        base: b'C',
        germline_pos: 1,
    }));

    let outcome = PassRuntime::execute(&plan, Simulation::new(), 0);

    let after_first_echo = outcome.revision_after("test_echo").unwrap();
    assert_eq!(after_first_echo.pool.len(), 1);
    assert!(outcome.revision_after("does_not_exist").is_none());
}

#[test]
fn pass_plan_construction_helpers() {
    let mut plan = PassPlan::new();
    assert!(plan.is_empty());
    assert_eq!(plan.len(), 0);

    plan.push(Box::new(TestEchoPass {
        base: b'A',
        germline_pos: 0,
    }));
    assert!(!plan.is_empty());
    assert_eq!(plan.len(), 1);
    assert_eq!(plan.passes()[0].name(), "test_echo");
}
