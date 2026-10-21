use super::*;
use crate::contract;
use crate::dist::AllelePoolDist;
use crate::ir::{NucHandle, Region};
use crate::refdata::ChainType;
use crate::replay::TraceCursor;
use crate::rng::Rng;
use crate::trace::Trace;

fn fixture() -> (RefDataConfig, Simulation) {
    let mut cfg = RefDataConfig::empty(ChainType::Vdj);
    for (name, sequence) in [
        ("original", b"AAACCCTGT".as_slice()),
        ("stop", b"TAACCCTGT"),
        ("valid_one", b"CCCAAATGT"),
        ("valid_three", b"GGGAAATGT"),
    ] {
        let _ = cfg.v_pool.push(Allele {
            name: name.into(),
            gene: name.into(),
            seq: sequence.to_vec(),
            segment: Segment::V,
            anchor: Some(6),
            functional_status: None,
            subregions: vec![],
        });
    }
    for (segment, sequence, anchor) in [
        (Segment::D, b"GGG".as_slice(), None),
        (Segment::J, b"TTTAAA".as_slice(), Some(0)),
    ] {
        let allele = Allele {
            name: format!("{segment:?}"),
            gene: format!("{segment:?}"),
            seq: sequence.to_vec(),
            segment,
            anchor,
            functional_status: None,
            subregions: vec![],
        };
        match segment {
            Segment::D => cfg.d_pool.push(allele),
            Segment::J => cfg.j_pool.push(allele),
            _ => unreachable!(),
        };
    }
    let mut builder = SimulationBuilder::from_simulation(Simulation::new());
    for segment in [Segment::V, Segment::D, Segment::J] {
        let allele = cfg.get(segment, AlleleId::new(0)).unwrap();
        let start = builder.peek().pool.len();
        for (i, &base) in allele.seq.iter().enumerate() {
            builder.push_nucleotide(Nucleotide::germline(base, i as u16, segment));
        }
        let end = builder.peek().pool.len();
        builder.assign_allele(
            segment,
            AlleleInstance::new(AlleleId::new(0)).with_haplotype(0),
        );
        builder.add_region(Region::new(
            segment,
            NucHandle::new(start as u32),
            NucHandle::new(end as u32),
        ));
    }
    let sim = builder.seal();
    assert_eq!(
        crate::productivity::assess(&sim, &cfg).productive,
        Some(true)
    );
    (cfg, sim)
}

fn pass(cfg: &RefDataConfig, genotype: bool, weights: [f64; 4]) -> ReceptorRevisionPass {
    let distribution = if weights.contains(&0.0) {
        AllelePoolDist::restricted_uniform(
            &cfg.v_pool,
            weights
                .iter()
                .enumerate()
                .filter(|(_, w)| **w > 0.0)
                .map(|(i, _)| AlleleId::new(i as u32))
                .collect(),
        )
    } else {
        AllelePoolDist::from_weights(&cfg.v_pool, weights.to_vec())
    };
    let pass = ReceptorRevisionPass::new(1.0, Box::new(distribution));
    if genotype {
        pass.with_genotype_constraint(GenotypeVConstraint {
            per_hap: [
                weights
                    .into_iter()
                    .enumerate()
                    .map(|(i, weight)| (AlleleId::new(i as u32), weight))
                    .collect(),
                vec![],
            ],
            same_haplotype: true,
        })
    } else {
        pass
    }
}

#[test]
fn replacement_filter_preserves_supported_weights_and_emits_only_committed_events() {
    let (cfg, sim) = fixture();
    let contracts = contract::productive();
    for genotype in [false, true] {
        let pass = pass(&cfg, genotype, [1., 100., 1., 3.]);
        for strict in [false, true] {
            for seed in 0..64 {
                let mut expected_rng = Rng::new(seed);
                expected_rng.next_f64(); // application coin
                let expected = if expected_rng.next_f64() * 4. <= 1. {
                    2
                } else {
                    3
                };
                let mut trace = Trace::new();
                let mut rng = Rng::new(seed);
                let mut events = Vec::new();
                let mut ctx = PassContext {
                    trace: &mut trace,
                    rng: &mut rng,
                    pass_index: 0,
                    refdata: Some(&cfg),
                    contracts: Some(&contracts),
                    feasibility: None,
                    reference_index: None,
                    replay_cursor: None,
                    event_log_sink: Some(&mut events),
                };
                let result = pass.execute_with_policy(&sim, &mut ctx, strict).unwrap();
                assert_eq!(
                    result
                        .assignments
                        .get(Segment::V)
                        .unwrap()
                        .allele_id
                        .index(),
                    expected
                );
                assert_eq!(
                    crate::productivity::assess(&result, &cfg).productive,
                    Some(true)
                );
                assert_eq!(trace.choices().len(), 3);
                assert_eq!(events.len(), 3);
                assert_eq!(
                    sim.assignments.get(Segment::V).unwrap().allele_id.index(),
                    0
                );
            }
        }
    }
}

#[test]
fn biologically_rejected_replacement_trace_is_rejected_in_both_modes() {
    let (cfg, sim) = fixture();
    let contracts = contract::productive();
    for genotype in [false, true] {
        let pass = pass(&cfg, genotype, [1., 100., 1., 3.]);
        for strict in [false, true] {
            let mut input = Trace::new();
            input.record_choice(
                address::ChoiceAddress::ReceptorRevisionApplied,
                ChoiceValue::Bool(true),
            );
            input.record_choice(
                address::ChoiceAddress::ReceptorRevisionVAllele,
                ChoiceValue::AlleleId(1),
            );
            input.record_choice(
                address::ChoiceAddress::ReceptorRevisionVTrim3,
                ChoiceValue::Int(0),
            );
            let mut cursor = TraceCursor::from_owned(input.choices().to_vec());
            let mut trace = Trace::new();
            let mut rng = Rng::new(9);
            let mut events = Vec::new();
            let mut ctx = PassContext {
                trace: &mut trace,
                rng: &mut rng,
                pass_index: 0,
                refdata: Some(&cfg),
                contracts: Some(&contracts),
                feasibility: None,
                reference_index: None,
                replay_cursor: Some(&mut cursor),
                event_log_sink: Some(&mut events),
            };
            assert!(matches!(
                pass.execute_with_policy(&sim, &mut ctx, strict),
                Err(PassError::ContractViolation { .. })
            ));
            assert!(events.is_empty());
        }
    }
}

#[test]
fn no_productive_replacement_is_explicit_noop_or_strict_empty_support() {
    let (cfg, sim) = fixture();
    let contracts = contract::productive();
    for genotype in [false, true] {
        let pass = pass(&cfg, genotype, [0., 1., 0., 0.]);
        for strict in [false, true] {
            let mut trace = Trace::new();
            let mut rng = Rng::new(9);
            let mut events = Vec::new();
            let mut ctx = PassContext {
                trace: &mut trace,
                rng: &mut rng,
                pass_index: 0,
                refdata: Some(&cfg),
                contracts: Some(&contracts),
                feasibility: None,
                reference_index: None,
                replay_cursor: None,
                event_log_sink: Some(&mut events),
            };
            let result = pass.execute_with_policy(&sim, &mut ctx, strict);
            if strict {
                assert!(matches!(result, Err(PassError::ConstraintSampling { .. })));
            } else {
                assert_eq!(
                    result
                        .unwrap()
                        .assignments
                        .get(Segment::V)
                        .unwrap()
                        .allele_id
                        .index(),
                    0
                );
                assert_eq!(trace.choices().len(), 1);
                assert_eq!(trace.choices()[0].value, ChoiceValue::Bool(false));
            }
            assert!(events.is_empty());
        }
    }
}
