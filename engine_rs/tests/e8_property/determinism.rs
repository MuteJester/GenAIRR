use super::common::{assembled_v_sim, uniform_s5f, SEED_RANGE};
use genairr_engine::dist::{EmpiricalLengthDist, UniformBase};
use genairr_engine::ir::NucHandle;
use genairr_engine::pass::PassPlan;
use genairr_engine::pass::testing::PassRuntime;
use genairr_engine::passes::{
    ContaminantPass, IndelPass, PCRErrorPass, QualityErrorPass, S5FMutationPass,
};

#[test]
fn property_determinism_full_corruption_stack() {
    let plan = || {
        let mut p = PassPlan::new();
        p.push(Box::new(S5FMutationPass::new(
            uniform_s5f(),
            Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)])),
        )));
        p.push(Box::new(PCRErrorPass::new(
            Box::new(EmpiricalLengthDist::from_pairs(vec![(2, 1.0)])),
            Box::new(UniformBase),
        )));
        p.push(Box::new(IndelPass::new(
            Box::new(EmpiricalLengthDist::from_pairs(vec![(2, 1.0)])),
            0.5,
            Box::new(UniformBase),
        )));
        p.push(Box::new(QualityErrorPass::new(
            Box::new(EmpiricalLengthDist::from_pairs(vec![(1, 1.0)])),
            Box::new(UniformBase),
        )));
        p.push(Box::new(ContaminantPass::new(0.2, Box::new(UniformBase))));
        p
    };

    for seed in 0..SEED_RANGE {
        let oa = PassRuntime::execute(&plan(), assembled_v_sim(), seed);
        let ob = PassRuntime::execute(&plan(), assembled_v_sim(), seed);

        assert_eq!(
            oa.trace.choices(),
            ob.trace.choices(),
            "seed {} non-deterministic trace",
            seed
        );

        let pa = &oa.final_simulation().pool;
        let pb = &ob.final_simulation().pool;
        assert_eq!(pa.len(), pb.len(), "seed {} non-deterministic len", seed);
        for i in 0..pa.len() as u32 {
            let na = pa.get(NucHandle::new(i)).unwrap();
            let nb = pb.get(NucHandle::new(i)).unwrap();
            assert_eq!(
                na.base, nb.base,
                "seed {} non-deterministic base at {}",
                seed, i
            );
            assert_eq!(
                na.flags.bits(),
                nb.flags.bits(),
                "seed {} non-deterministic flags at {}",
                seed,
                i
            );
        }
    }
}
