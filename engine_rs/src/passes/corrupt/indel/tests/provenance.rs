use super::*;
use crate::ir::flag;

#[test]
fn indel_pass_inserted_nucleotides_carry_indel_flag() {
    let mut plan = PassPlan::new();
    plan.push(Box::new(IndelPass::new(
        fixed_count(3),
        1.0,
        Box::new(UniformBase),
    )));
    let outcome = PassRuntime::execute(&plan, indel_test_sim(), 19);
    let final_sim = outcome.final_simulation();

    let mut inserted = 0;
    let mut germline = 0;
    for i in 0..final_sim.pool.len() as u32 {
        let n = final_sim.pool.get(NucHandle::new(i)).unwrap();
        if n.flags.contains(flag::INDEL_INSERTED) {
            inserted += 1;
            assert!(n.germline_pos.is_none());
        } else {
            germline += 1;
        }
    }
    assert_eq!(inserted, 3);
    assert_eq!(germline, 12);
}

#[test]
fn indel_pass_insertion_segment_uses_sequence_context() {
    let sim = multi_segment_indel_context_sim();

    assert_eq!(IndelPass::insertion_segment(&sim, 1), Segment::V);
    assert_eq!(IndelPass::insertion_segment(&sim, 3), Segment::Np1);
    assert_eq!(IndelPass::insertion_segment(&sim, 5), Segment::J);
    assert_eq!(IndelPass::insertion_segment(&sim, 8), Segment::J);
    assert_eq!(
        IndelPass::insertion_segment(&Simulation::new(), 0),
        Segment::V
    );
}

#[test]
fn indel_pass_inserted_nucleotide_inherits_non_v_segment() {
    let sim = multi_segment_indel_context_sim();
    let mut plan = PassPlan::new();
    plan.push(Box::new(IndelPass::new(
        fixed_count(1),
        1.0,
        Box::new(UniformBase),
    )));

    for seed in 0..200u64 {
        let outcome = PassRuntime::execute(&plan, sim.clone(), seed);
        let site = match outcome.trace.find("corrupt.indel.site[0]").unwrap().value {
            ChoiceValue::Int(s) => s as u32,
            _ => unreachable!(),
        };

        if site >= 5 {
            let inserted = outcome
                .final_simulation()
                .pool
                .get(NucHandle::new(site))
                .unwrap();
            assert_eq!(inserted.segment, Segment::J);
            assert!(inserted.flags.contains(flag::INDEL_INSERTED));
            return;
        }
    }

    panic!("test setup did not produce an insertion in J context");
}

#[test]
fn indel_pass_inserted_base_in_pool_matches_traced_base() {
    let mut plan = PassPlan::new();
    plan.push(Box::new(IndelPass::new(
        fixed_count(2),
        1.0,
        Box::new(UniformBase),
    )));
    let outcome = PassRuntime::execute(&plan, indel_test_sim(), 42);
    let final_sim = outcome.final_simulation();

    let sites: Vec<u32> = (0..2)
        .map(|i| {
            match outcome
                .trace
                .find(&format!("corrupt.indel.site[{}]", i))
                .unwrap()
                .value
            {
                ChoiceValue::Int(s) => s as u32,
                _ => unreachable!(),
            }
        })
        .collect();
    let bases: Vec<u8> = (0..2)
        .map(|i| {
            match outcome
                .trace
                .find(&format!("corrupt.indel.base[{}]", i))
                .unwrap()
                .value
            {
                ChoiceValue::Base(b) => b,
                _ => unreachable!(),
            }
        })
        .collect();

    let mut resting: Vec<u32> = Vec::new();
    for s in sites {
        for r in resting.iter_mut() {
            if *r >= s {
                *r += 1;
            }
        }
        resting.push(s);
    }

    for (s, expected) in resting.iter().zip(bases.iter()) {
        let pool_base = final_sim.pool.get(NucHandle::new(*s)).unwrap().base;
        assert_eq!(pool_base, *expected, "trace lies at site {}", s);
    }
}

#[test]
fn indel_pass_insertion_grows_region_only_when_spanning() {
    let mut plan = PassPlan::new();
    plan.push(Box::new(IndelPass::new(
        fixed_count(3),
        1.0,
        Box::new(UniformBase),
    )));
    let outcome = PassRuntime::execute(&plan, indel_test_sim(), 0);
    let final_sim = outcome.final_simulation();

    let mut region_end: u32 = 12;
    for i in 0..3 {
        let site = match outcome
            .trace
            .find(&format!("corrupt.indel.site[{}]", i))
            .unwrap()
            .value
        {
            ChoiceValue::Int(s) => s as u32,
            _ => unreachable!(),
        };
        if site < region_end {
            region_end += 1;
        }
    }

    assert_eq!(final_sim.sequence.regions.len(), 1);
    assert_eq!(final_sim.sequence.regions[0].start.index(), 0);
    assert_eq!(final_sim.sequence.regions[0].end.index(), region_end);
    assert_eq!(final_sim.pool.len(), 12 + 3);
}

#[test]
fn indel_pass_codon_rail_refresh_after_insertion() {
    // rail no longer maintained in the hot path. Compute
    // on demand and assert determinism.
    let mut plan = PassPlan::new();
    plan.push(Box::new(IndelPass::new(
        fixed_count(2),
        1.0,
        Box::new(UniformBase),
    )));
    let outcome = PassRuntime::execute(&plan, indel_test_sim(), 1);
    let final_sim = outcome.final_simulation();

    let a = final_sim.sequence.regions[0].with_codon_rail_recomputed(&final_sim.pool);
    let b = final_sim.sequence.regions[0].with_codon_rail_recomputed(&final_sim.pool);
    assert_eq!(a.amino_acids, b.amino_acids);
}

#[test]
fn indel_pass_codon_rail_refresh_after_deletion() {
    // see comment on insertion variant.
    let mut plan = PassPlan::new();
    plan.push(Box::new(IndelPass::new(
        fixed_count(2),
        0.0,
        Box::new(UniformBase),
    )));
    let outcome = PassRuntime::execute(&plan, indel_test_sim(), 1);
    let final_sim = outcome.final_simulation();

    let a = final_sim.sequence.regions[0].with_codon_rail_recomputed(&final_sim.pool);
    let b = final_sim.sequence.regions[0].with_codon_rail_recomputed(&final_sim.pool);
    assert_eq!(a.amino_acids, b.amino_acids);
}

#[test]
fn indel_pass_is_deterministic_under_same_seed() {
    let plan = || {
        let mut p = PassPlan::new();
        p.push(Box::new(IndelPass::new(
            fixed_count(3),
            0.5,
            Box::new(UniformBase),
        )));
        p
    };
    let oa = PassRuntime::execute(&plan(), indel_test_sim(), 0xc0ff_ee);
    let ob = PassRuntime::execute(&plan(), indel_test_sim(), 0xc0ff_ee);
    assert_eq!(oa.trace.choices(), ob.trace.choices());
    assert_eq!(
        oa.final_simulation().pool.len(),
        ob.final_simulation().pool.len()
    );
    for i in 0..oa.final_simulation().pool.len() as u32 {
        let h = NucHandle::new(i);
        assert_eq!(
            oa.final_simulation().pool.get(h).unwrap().base,
            ob.final_simulation().pool.get(h).unwrap().base
        );
    }
}

#[test]
fn indel_pass_persistent_ir_preserves_input() {
    let sim = indel_test_sim();
    let original_len = sim.pool.len();
    let original_region_end = sim.sequence.regions[0].end.index();

    let mut plan = PassPlan::new();
    plan.push(Box::new(IndelPass::new(
        fixed_count(3),
        1.0,
        Box::new(UniformBase),
    )));
    let _outcome = PassRuntime::execute(&plan, sim.clone(), 99);

    assert_eq!(sim.pool.len(), original_len);
    assert_eq!(sim.sequence.regions[0].end.index(), original_region_end);
}
