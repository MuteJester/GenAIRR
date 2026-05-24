use super::*;
use crate::assignment::AlleleInstance;
use crate::dist::EmpiricalLengthDist;
use crate::ir::{NucHandle, Nucleotide, Region, Segment};
use crate::pass::PassPlan;
    use crate::pass::testing::PassRuntime;
use crate::refdata::{Allele, AlleleId, ChainType, RefDataConfig};
use crate::s5f::{S5F_NUM_CONTEXTS, S5F_SUBSTITUTION_LEN};
use crate::trace::ChoiceValue;

mod constraints;
mod core;
mod provenance;

fn s5f_uniform_kernel() -> S5FKernel {
    S5FKernel::new(
        vec![1.0; S5F_NUM_CONTEXTS],
        vec![0.25; S5F_SUBSTITUTION_LEN],
    )
}

fn s5f_zero_kernel() -> S5FKernel {
    S5FKernel::new(
        vec![0.0; S5F_NUM_CONTEXTS],
        vec![0.25; S5F_SUBSTITUTION_LEN],
    )
}

fn s5f_stop_filter_kernel(include_safe_base: bool) -> S5FKernel {
    let context =
        S5FKernel::encode_context(b'T', b'A', b'C', b'A', b'A').expect("canonical context");
    let mut mutability = vec![0.0; S5F_NUM_CONTEXTS];
    mutability[context as usize] = 1.0;

    let mut substitution = vec![0.0; S5F_SUBSTITUTION_LEN];
    let offset = context as usize * 4;
    substitution[offset] = 1.0;
    if include_safe_base {
        substitution[offset + 1] = 1.0;
    }

    S5FKernel::new(mutability, substitution)
}

fn s5f_single_mutation_plan(kernel: S5FKernel) -> PassPlan {
    let mut plan = PassPlan::new();
    plan.push(Box::new(S5FMutationPass::new(
        kernel,
        Box::new(EmpiricalLengthDist::from_pairs(vec![(1, 1.0)])),
    )));
    plan
}

fn s5f_productive_vj_fixture() -> (RefDataConfig, Simulation) {
    let mut cfg = RefDataConfig::empty(ChainType::Vj);
    let _ = cfg.v_pool.push(Allele {
        name: "v_s5f*01".into(),
        gene: "v_s5f".into(),
        seq: b"TACAAA".to_vec(),
        segment: Segment::V,
        anchor: Some(0),
    });
    let _ = cfg.j_pool.push(Allele {
        name: "j_s5f*01".into(),
        gene: "j_s5f".into(),
        seq: b"TGG".to_vec(),
        segment: Segment::J,
        anchor: Some(0),
    });

    let mut sim = Simulation::new();
    for (i, &b) in b"TACAAA".iter().enumerate() {
        let (next, _) = sim.with_nucleotide_pushed(Nucleotide::germline(b, i as u16, Segment::V));
        sim = next;
    }
    let v_region = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(6));
    sim = sim.with_region_added(v_region);

    for (i, &b) in b"TGG".iter().enumerate() {
        let (next, _) = sim.with_nucleotide_pushed(Nucleotide::germline(b, i as u16, Segment::J));
        sim = next;
    }
    let j_region = Region::new(Segment::J, NucHandle::new(6), NucHandle::new(9));
    sim = sim.with_region_added(j_region);

    sim = sim
        .with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(0)))
        .with_allele_assigned(Segment::J, AlleleInstance::new(AlleleId::new(0)));

    (cfg, sim)
}

fn find_seed_for_s5f_unconstrained_base(
    kernel: &S5FKernel,
    sim: &Simulation,
    target_base: u8,
) -> u64 {
    for seed in 0..512u64 {
        let outcome =
            PassRuntime::execute(&s5f_single_mutation_plan(kernel.clone()), sim.clone(), seed);
        if let Some(rec) = outcome.trace.find("mutate.s5f.base[0]") {
            if rec.value == ChoiceValue::Base(target_base) {
                return seed;
            }
        }
    }
    panic!(
        "no seed in search range produced S5F base {}",
        target_base as char
    );
}

fn s5f_test_sim() -> Simulation {
    let mut sim = Simulation::new();
    for (i, b) in b"AAACCCGGGTTTACGTACGT".iter().enumerate() {
        let (next, _) = sim.with_nucleotide_pushed(Nucleotide::germline(*b, i as u16, Segment::V));
        sim = next;
    }
    let region = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(20));
    sim.with_region_added(region)
}
