use super::common::{assembled_v_sim, SEED_RANGE};
use genairr_engine::dist::{EmpiricalLengthDist, UniformBase};
use genairr_engine::ir::{flag, NucHandle};
use genairr_engine::pass::{PassPlan, PassRuntime};
use genairr_engine::passes::IndelPass;
use genairr_engine::trace::ChoiceValue;

#[test]
fn property_indel_pool_length_matches_trace_arithmetic() {
    let mut plan = PassPlan::new();
    plan.push(Box::new(IndelPass::new(
        Box::new(EmpiricalLengthDist::from_pairs(vec![(5, 1.0)])),
        0.5,
        Box::new(UniformBase),
    )));

    for seed in 0..SEED_RANGE {
        let initial_len = assembled_v_sim().pool.len() as i64;
        let outcome = PassRuntime::execute(&plan, assembled_v_sim(), seed);

        let count = match outcome.trace.find("corrupt.indel.count").unwrap().value {
            ChoiceValue::Int(n) => n,
            _ => unreachable!(),
        };
        assert_eq!(count, 5, "seed {} had unexpected count", seed);

        let mut insertions = 0i64;
        let mut deletions = 0i64;
        for i in 0..5 {
            let kind = matches!(
                outcome
                    .trace
                    .find(&format!("corrupt.indel.kind[{}]", i))
                    .unwrap()
                    .value,
                ChoiceValue::Bool(true)
            );
            let site = match outcome
                .trace
                .find(&format!("corrupt.indel.site[{}]", i))
                .unwrap()
                .value
            {
                ChoiceValue::Int(s) => s,
                _ => unreachable!(),
            };
            if kind {
                insertions += 1;
            } else if site >= 0 {
                deletions += 1;
            }
        }

        let expected_len = initial_len + insertions - deletions;
        assert_eq!(
            outcome.final_simulation().pool.len() as i64,
            expected_len,
            "seed {}: expected pool len {} (initial {} + {} ins - {} del), got {}",
            seed,
            expected_len,
            initial_len,
            insertions,
            deletions,
            outcome.final_simulation().pool.len()
        );
    }
}

#[test]
fn property_indel_inserted_flag_count_at_most_traced_insertions() {
    let mut plan = PassPlan::new();
    plan.push(Box::new(IndelPass::new(
        Box::new(EmpiricalLengthDist::from_pairs(vec![(6, 1.0)])),
        0.5,
        Box::new(UniformBase),
    )));

    for seed in 0..SEED_RANGE {
        let outcome = PassRuntime::execute(&plan, assembled_v_sim(), seed);
        let final_sim = outcome.final_simulation();

        let mut flagged = 0u32;
        for i in 0..final_sim.pool.len() as u32 {
            let n = final_sim.pool.get(NucHandle::new(i)).unwrap();
            if n.flags.contains(flag::INDEL_INSERTED) {
                flagged += 1;
                assert!(
                    n.germline_pos.is_none(),
                    "seed {}: flagged nucleotide at {} has germline pos {:?}",
                    seed,
                    i,
                    n.germline_pos.get()
                );
            }
        }

        let traced_insertions = (0..6)
            .filter(|i| {
                matches!(
                    outcome
                        .trace
                        .find(&format!("corrupt.indel.kind[{}]", i))
                        .unwrap()
                        .value,
                    ChoiceValue::Bool(true)
                )
            })
            .count() as u32;

        assert!(
            flagged <= traced_insertions,
            "seed {}: flagged count {} exceeds traced insertions {}",
            seed,
            flagged,
            traced_insertions
        );
    }
}

#[test]
fn property_indel_inserted_flag_count_equals_insertions_when_no_deletions() {
    let mut plan = PassPlan::new();
    plan.push(Box::new(IndelPass::new(
        Box::new(EmpiricalLengthDist::from_pairs(vec![(5, 1.0)])),
        1.0,
        Box::new(UniformBase),
    )));

    for seed in 0..SEED_RANGE {
        let outcome = PassRuntime::execute(&plan, assembled_v_sim(), seed);
        let final_sim = outcome.final_simulation();

        let flagged = (0..final_sim.pool.len() as u32)
            .filter(|i| {
                final_sim
                    .pool
                    .get(NucHandle::new(*i))
                    .unwrap()
                    .flags
                    .contains(flag::INDEL_INSERTED)
            })
            .count();

        assert_eq!(
            flagged, 5,
            "seed {}: expected 5 flagged nucleotides, got {}",
            seed, flagged
        );
    }
}
