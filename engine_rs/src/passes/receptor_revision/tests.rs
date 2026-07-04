    use super::*;
    use crate::assignment::AlleleInstance;
    use crate::contract::ContractSet;
    use crate::dist::AllelePoolDist;
    use crate::ir::{NucHandle, Region, Segment, SimulationEvent};
    use crate::pass::testing::PassRuntime;
    use crate::pass::{PassError, PassPlan};
    use crate::refdata::{Allele, AllelePool, ChainType, RefDataConfig};
    use crate::replay::TraceCursor;
    use crate::rng::Rng;
    use crate::trace::Trace;

    fn allele(name: &str, seq: &[u8]) -> Allele {
        Allele {
            name: name.to_string(),
            gene: name.split('*').next().unwrap_or(name).to_string(),
            seq: seq.to_vec(),
            segment: Segment::V,
            anchor: None,
            functional_status: None,
            subregions: Vec::new(),
        }
    }

    /// Reference data with two V alleles whose retained 6-byte
    /// prefix (under 0 trim) differs. The original assembled V is
    /// `AAAAAA` matching allele V0 exactly; replacement candidate
    /// V1 carries `GGGGGGCC` — length 8, retains 6 bytes at
    /// `trim_3 = 2`.
    fn two_v_refdata() -> (RefDataConfig, AlleleId, AlleleId) {
        let mut cfg = RefDataConfig::empty(ChainType::Vdj);
        let v0 = cfg.v_pool.push(allele("V1*01", b"AAAAAA"));
        let v1 = cfg.v_pool.push(allele("V1*02", b"GGGGGGCC"));
        (cfg, v0, v1)
    }

    /// Build a sim with the V slot already assigned to V0 and a
    /// V region of length 6 (`AAAAAA`). Doubles for tests as the
    /// post-recombination starting state.
    fn sim_v_assembled(v0: AlleleId) -> Simulation {
        let mut sim = Simulation::new();
        for (i, &b) in b"AAAAAA".iter().enumerate() {
            let (next, _) = sim.with_nucleotide_pushed(Nucleotide::germline(
                b,
                i as u16,
                Segment::V,
            ));
            sim = next;
        }
        sim.with_allele_assigned(Segment::V, AlleleInstance::new(v0))
            .with_region_added(Region::new(
                Segment::V,
                NucHandle::new(0),
                NucHandle::new(6),
            ))
    }

    fn v_only_pool(cfg: &RefDataConfig) -> AllelePool {
        cfg.v_pool.clone()
    }

    fn run_pass(
        prob: f64,
        seed: u64,
        cfg: RefDataConfig,
        sim: Simulation,
    ) -> (Trace, Simulation) {
        let mut plan = PassPlan::new();
        let pool = v_only_pool(&cfg);
        plan.push(Box::new(ReceptorRevisionPass::new(
            prob,
            Box::new(AllelePoolDist::uniform(&pool)),
        )));
        let outcome = PassRuntime::execute_with_refdata(&plan, sim, seed, &cfg);
        let final_sim = outcome.final_simulation().clone();
        (outcome.trace, final_sim)
    }

    fn run_with_ctx(
        pass: &ReceptorRevisionPass,
        cfg: &RefDataConfig,
        contracts: Option<&ContractSet>,
        initial: Simulation,
        cursor: Option<&mut TraceCursor>,
        event_sink: Option<&mut Vec<SimulationEvent>>,
    ) -> Result<(Trace, Simulation), PassError> {
        let mut trace = Trace::new();
        let mut rng = Rng::new(0xc0ff_ee);
        let mut ctx = PassContext {
            trace: &mut trace,
            rng: &mut rng,
            pass_index: 0,
            refdata: Some(cfg),
            contracts,
            feasibility: None,
            reference_index: None,
            replay_cursor: cursor,
            event_log_sink: event_sink,
        };
        let next = pass.execute_checked(&initial, &mut ctx)?;
        Ok((trace, next))
    }

    // ── Construction guards ──────────────────────────────────────

    #[test]
    #[should_panic(expected = "prob must be in [0.0, 1.0]")]
    fn receptor_revision_rejects_out_of_range_prob() {
        let (cfg, _, _) = two_v_refdata();
        let _ = ReceptorRevisionPass::new(
            1.5,
            Box::new(AllelePoolDist::uniform(&cfg.v_pool)),
        );
    }

    #[test]
    #[should_panic(expected = "prob must be in [0.0, 1.0]")]
    fn receptor_revision_rejects_negative_prob() {
        let (cfg, _, _) = two_v_refdata();
        let _ = ReceptorRevisionPass::new(
            -0.5,
            Box::new(AllelePoolDist::uniform(&cfg.v_pool)),
        );
    }

    // ── genotype-aware candidate selection ──────────────────────

    fn constraint(per_hap: [Vec<(AlleleId, f64)>; 2], same: bool) -> GenotypeVConstraint {
        GenotypeVConstraint { per_hap, same_haplotype: same }
    }

    fn geno_sim(v0: AlleleId) -> Simulation {
        sim_v_assembled(v0)
            .with_allele_assigned(Segment::V, AlleleInstance::new(v0).with_haplotype(0))
    }

    fn geno_pass(cfg: &RefDataConfig, per_hap: [Vec<(AlleleId, f64)>; 2], same: bool) -> ReceptorRevisionPass {
        ReceptorRevisionPass::new(1.0, Box::new(AllelePoolDist::uniform(&cfg.v_pool)))
            .with_genotype_constraint(constraint(per_hap, same))
    }

    #[test]
    fn genotype_same_haplotype_excludes_current_and_restricts_to_chromosome() {
        let (cfg, v0, v1) = two_v_refdata();
        // hap0 carries {V0, V1}; hap1 carries {V0}. Current is V0 on hap0.
        let pass = geno_pass(&cfg, [vec![(v0, 1.0), (v1, 1.0)], vec![(v0, 1.0)]], true);
        let (trace, after) = run_with_ctx(&pass, &cfg, None, geno_sim(v0), None, None).unwrap();
        assert_eq!(trace.find("receptor_revision.applied").unwrap().value, ChoiceValue::Bool(true));
        // exclude-current: the only eligible alternate on hap0 is V1
        assert_eq!(after.assignments.get(Segment::V).unwrap().allele_id, v1);
        assert_eq!(after.assignments.get(Segment::V).unwrap().receptor_revision_original_id, Some(v0));
        assert_eq!(after.assignments.get(Segment::V).unwrap().haplotype, Some(0));
    }

    fn run_permissive(pass: &ReceptorRevisionPass, cfg: &RefDataConfig, initial: Simulation) -> (Trace, Simulation) {
        let mut trace = Trace::new();
        let mut rng = Rng::new(0xc0ff_ee);
        let mut ctx = PassContext {
            trace: &mut trace,
            rng: &mut rng,
            pass_index: 0,
            refdata: Some(cfg),
            contracts: None,
            feasibility: None,
            reference_index: None,
            replay_cursor: None,
            event_log_sink: None,
        };
        let next = pass.execute(&initial, &mut ctx);
        (trace, next)
    }

    #[test]
    fn genotype_no_eligible_alternate_permissive_applied_false() {
        let (cfg, v0, _v1) = two_v_refdata();
        let pass = geno_pass(&cfg, [vec![(v0, 1.0)], vec![(v0, 1.0)]], true);
        let (trace, after) = run_permissive(&pass, &cfg, geno_sim(v0));
        assert_eq!(trace.find("receptor_revision.applied").unwrap().value, ChoiceValue::Bool(false));
        assert!(trace.find("receptor_revision.v_allele").is_none());
        assert_eq!(after.assignments.get(Segment::V).unwrap().allele_id, v0);
    }

    #[test]
    fn genotype_no_eligible_alternate_strict_errors() {
        let (cfg, v0, _v1) = two_v_refdata();
        let pass = geno_pass(&cfg, [vec![(v0, 1.0)], vec![(v0, 1.0)]], true);
        let err = run_with_ctx(&pass, &cfg, None, geno_sim(v0), None, None).unwrap_err();
        assert!(matches!(err, PassError::ConstraintSampling { .. }), "got {err:?}");
    }

    #[test]
    fn genotype_both_haplotypes_admits_other_chromosome_allele() {
        let (cfg, v0, v1) = two_v_refdata();
        // V1 carried only on hap1; same_haplotype=false aggregates both.
        let pass = geno_pass(&cfg, [vec![(v0, 1.0)], vec![(v1, 1.0)]], false);
        let (_t, after) = run_with_ctx(&pass, &cfg, None, geno_sim(v0), None, None).unwrap();
        assert_eq!(after.assignments.get(Segment::V).unwrap().allele_id, v1);
        // haplotype provenance preserved as the original rearrangement chromosome (0)
        assert_eq!(after.assignments.get(Segment::V).unwrap().haplotype, Some(0));
    }

    #[test]
    fn genotype_missing_haplotype_stamp_errors() {
        let (cfg, v0, v1) = two_v_refdata();
        let pass = geno_pass(&cfg, [vec![(v0, 1.0), (v1, 1.0)], vec![(v0, 1.0)]], true);
        // sim_v_assembled assigns V0 WITHOUT a haplotype stamp.
        let err = run_with_ctx(&pass, &cfg, None, sim_v_assembled(v0), None, None).unwrap_err();
        assert!(matches!(err, PassError::InvalidPlanState { .. }), "got {err:?}");
    }

    #[test]
    fn genotype_strict_post_event_contract_rejection_errors() {
        use crate::contract::{Contract, ContractViolation};

        // A contract whose verify() always rejects the post-event state.
        struct RejectAll;
        impl Contract for RejectAll {
            fn name(&self) -> &str {
                "reject_all_test"
            }
            fn verify(
                &self,
                _sim: &Simulation,
                _refdata: Option<&RefDataConfig>,
            ) -> Result<(), ContractViolation> {
                Err(ContractViolation::new(self.name(), "rejected by test"))
            }
        }

        let (cfg, v0, v1) = two_v_refdata();
        let pass = geno_pass(&cfg, [vec![(v0, 1.0), (v1, 1.0)], vec![(v0, 1.0)]], true);
        let contracts = ContractSet::new().with(Box::new(RejectAll));
        // Strict mode (execute_checked) with a rejecting contract must surface a
        // contract violation — parity with the non-genotype path.
        let err =
            run_with_ctx(&pass, &cfg, Some(&contracts), geno_sim(v0), None, None).unwrap_err();
        assert!(matches!(err, PassError::ContractViolation { .. }), "got {err:?}");
    }

    #[test]
    fn genotype_replay_strict_post_event_contract_rejection_errors() {
        use crate::contract::{Contract, ContractViolation};

        struct RejectAll;
        impl Contract for RejectAll {
            fn name(&self) -> &str {
                "reject_all_test"
            }
            fn verify(
                &self,
                _sim: &Simulation,
                _refdata: Option<&RefDataConfig>,
            ) -> Result<(), ContractViolation> {
                Err(ContractViolation::new(self.name(), "rejected by test"))
            }
        }

        let (cfg, v0, v1) = two_v_refdata();
        let pass = geno_pass(&cfg, [vec![(v0, 1.0), (v1, 1.0)], vec![(v0, 1.0)]], true);
        let contracts = ContractSet::new().with(Box::new(RejectAll));
        // A valid applied=true replay (v1, trim 2) must STILL be rejected by the
        // post-event contract in strict mode — replay parity with the fresh path.
        let mut cursor = TraceCursor::from_owned(replay_records(true, Some(v1.index()), Some(2)));
        let err = run_with_ctx(&pass, &cfg, Some(&contracts), geno_sim(v0), Some(&mut cursor), None)
            .unwrap_err();
        assert!(matches!(err, PassError::ContractViolation { .. }), "got {err:?}");
    }

    // ── genotype-aware replay validation ────────────────────────

    fn replay_records(applied: bool, allele: Option<u32>, trim: Option<i64>) -> Vec<crate::trace::ChoiceRecord> {
        let mut t = Trace::new();
        t.record_choice(address::ChoiceAddress::ReceptorRevisionApplied, ChoiceValue::Bool(applied));
        if let Some(a) = allele {
            t.record_choice(address::ChoiceAddress::ReceptorRevisionVAllele, ChoiceValue::AlleleId(a));
        }
        if let Some(tr) = trim {
            t.record_choice(address::ChoiceAddress::ReceptorRevisionVTrim3, ChoiceValue::Int(tr));
        }
        t.choices().to_vec()
    }

    #[test]
    fn replay_genotype_valid_replacement_reproduces() {
        let (cfg, v0, v1) = two_v_refdata();
        let pass = geno_pass(&cfg, [vec![(v0, 1.0), (v1, 1.0)], vec![(v0, 1.0)]], true);
        let mut cursor = TraceCursor::from_owned(replay_records(true, Some(v1.index()), Some(2)));
        let (_t, after) = run_with_ctx(&pass, &cfg, None, geno_sim(v0), Some(&mut cursor), None).unwrap();
        assert_eq!(after.assignments.get(Segment::V).unwrap().allele_id, v1);
        assert!(cursor.is_drained());
    }

    #[test]
    fn replay_genotype_allele_not_carried_on_haplotype_errors() {
        let (cfg, v0, v1) = two_v_refdata();
        // hap0 carries only V0; recorded V1 is not carried on the drawn chromosome.
        let pass = geno_pass(&cfg, [vec![(v0, 1.0)], vec![(v0, 1.0), (v1, 1.0)]], true);
        let mut cursor = TraceCursor::from_owned(replay_records(true, Some(v1.index()), Some(2)));
        let err = run_with_ctx(&pass, &cfg, None, geno_sim(v0), Some(&mut cursor), None).unwrap_err();
        assert!(matches!(err, PassError::InvalidDistributionOutput { .. }), "got {err:?}");
    }

    #[test]
    fn replay_genotype_unresolvable_allele_id_reports_missing_allele() {
        let (cfg, v0, _v1) = two_v_refdata();
        let pass = geno_pass(&cfg, [vec![(v0, 1.0)], vec![(v0, 1.0)]], true);
        // out-of-range recorded v_allele -> missing_allele (not "not_carried"),
        // matching the non-genotype replay diagnostic contract.
        let mut cursor = TraceCursor::from_owned(replay_records(true, Some(999), Some(0)));
        let err = run_with_ctx(&pass, &cfg, None, geno_sim(v0), Some(&mut cursor), None).unwrap_err();
        match err {
            PassError::MissingAllele { allele_id, .. } => assert_eq!(allele_id, 999),
            other => panic!("expected MissingAllele, got {other:?}"),
        }
    }

    #[test]
    fn replay_genotype_equals_current_allele_errors() {
        let (cfg, v0, v1) = two_v_refdata();
        let pass = geno_pass(&cfg, [vec![(v0, 1.0), (v1, 1.0)], vec![(v0, 1.0)]], true);
        let mut cursor = TraceCursor::from_owned(replay_records(true, Some(v0.index()), Some(0)));
        let err = run_with_ctx(&pass, &cfg, None, geno_sim(v0), Some(&mut cursor), None).unwrap_err();
        assert!(matches!(err, PassError::InvalidDistributionOutput { .. }), "got {err:?}");
    }

    #[test]
    fn replay_genotype_trim_length_mismatch_errors() {
        let (cfg, v0, v1) = two_v_refdata();
        let pass = geno_pass(&cfg, [vec![(v0, 1.0), (v1, 1.0)], vec![(v0, 1.0)]], true);
        // V1 len 8, old 6 -> trim must be 2; record 3 (retains 5) => mismatch.
        let mut cursor = TraceCursor::from_owned(replay_records(true, Some(v1.index()), Some(3)));
        let err = run_with_ctx(&pass, &cfg, None, geno_sim(v0), Some(&mut cursor), None).unwrap_err();
        assert!(matches!(err, PassError::InvalidPlanState { .. }), "got {err:?}");
    }

    #[test]
    fn genotype_signature_differs_by_candidate_set_and_same_haplotype() {
        let (cfg, v0, v1) = two_v_refdata();
        let no_geno = ReceptorRevisionPass::new(0.5, Box::new(AllelePoolDist::uniform(&cfg.v_pool)))
            .parameter_signature();
        let g_true = geno_pass_prob(&cfg, [vec![(v0, 1.0), (v1, 1.0)], vec![(v0, 1.0)]], true).parameter_signature();
        let g_false = geno_pass_prob(&cfg, [vec![(v0, 1.0), (v1, 1.0)], vec![(v0, 1.0)]], false).parameter_signature();
        let g_other = geno_pass_prob(&cfg, [vec![(v0, 2.0)], vec![(v1, 1.0)]], true).parameter_signature();
        assert_ne!(g_true, no_geno);
        assert_ne!(g_true, g_false);
        assert_ne!(g_true, g_other);
    }

    fn geno_pass_prob(cfg: &RefDataConfig, per_hap: [Vec<(AlleleId, f64)>; 2], same: bool) -> ReceptorRevisionPass {
        ReceptorRevisionPass::new(0.5, Box::new(AllelePoolDist::uniform(&cfg.v_pool)))
            .with_genotype_constraint(constraint(per_hap, same))
    }

    // ── prob=0: no replacement ──────────────────────────────────

    #[test]
    fn prob_zero_records_applied_false_and_no_mutation_events() {
        let (cfg, v0, _) = two_v_refdata();
        let pass = ReceptorRevisionPass::new(0.0, Box::new(AllelePoolDist::uniform(&cfg.v_pool)));
        let mut events: Vec<SimulationEvent> = Vec::new();
        let (trace, after) = run_with_ctx(
            &pass,
            &cfg,
            None,
            sim_v_assembled(v0),
            None,
            Some(&mut events),
        )
        .unwrap();

        // Applied=false recorded.
        let rec = trace
            .find("receptor_revision.applied")
            .expect("applied Bool must be recorded even when prob = 0");
        assert_eq!(rec.value, ChoiceValue::Bool(false));
        // No allele/trim records.
        assert!(trace.find("receptor_revision.v_allele").is_none());
        assert!(trace.find("receptor_revision.v_trim_3").is_none());

        // No state-changing events.
        assert!(events.iter().all(|e| !matches!(
            e,
            SimulationEvent::AssignmentChanged { .. }
                | SimulationEvent::TrimChanged { .. }
                | SimulationEvent::SegmentReplaced { .. }
        )));

        // V assignment + pool bytes unchanged.
        assert_eq!(after.assignments.get(Segment::V).unwrap().allele_id, v0);
        let bases: Vec<u8> = after.pool.as_slice().iter().map(|n| n.base).collect();
        assert_eq!(&bases, b"AAAAAA");
    }

    // ── prob=1: replacement ─────────────────────────────────────

    #[test]
    fn prob_one_records_three_choices_and_emits_three_events() {
        let (cfg, v0, v1) = two_v_refdata();
        let pass = ReceptorRevisionPass::new(1.0, Box::new(AllelePoolDist::uniform(&cfg.v_pool)));
        let mut events: Vec<SimulationEvent> = Vec::new();
        let (trace, after) = run_with_ctx(
            &pass,
            &cfg,
            None,
            sim_v_assembled(v0),
            None,
            Some(&mut events),
        )
        .unwrap();

        // applied = true
        assert_eq!(
            trace.find("receptor_revision.applied").unwrap().value,
            ChoiceValue::Bool(true),
        );
        // The only length-eligible candidate at old_v_len=6 is V1
        // (length 8). V0 is length 6 too — same length is eligible
        // (trim_3=0), so RNG decides. Either way, the recorded
        // pair must satisfy `len - trim_3 == 6`.
        let recorded_id = match trace.find("receptor_revision.v_allele").unwrap().value {
            ChoiceValue::AlleleId(id) => id,
            _ => panic!("expected AlleleId record"),
        };
        let recorded_trim = match trace.find("receptor_revision.v_trim_3").unwrap().value {
            ChoiceValue::Int(t) => t,
            _ => panic!("expected Int record"),
        };
        let recorded_allele = cfg
            .get(Segment::V, AlleleId::new(recorded_id))
            .expect("recorded allele must resolve");
        assert_eq!(
            recorded_allele.len() as i64 - recorded_trim,
            6,
            "retained length must equal old V region length"
        );
        let _unused = v1;

        // Exactly one of each of the three replacement events.
        let assignments = events
            .iter()
            .filter(|e| matches!(e, SimulationEvent::AssignmentChanged { segment: Segment::V, .. }))
            .count();
        let trims = events
            .iter()
            .filter(|e| matches!(e, SimulationEvent::TrimChanged { segment: Segment::V, end: TrimEnd::Three, .. }))
            .count();
        let replaces = events
            .iter()
            .filter(|e| matches!(e, SimulationEvent::SegmentReplaced { segment: Segment::V, .. }))
            .count();
        assert_eq!(assignments, 1);
        assert_eq!(trims, 1);
        assert_eq!(replaces, 1);

        // The committed V assignment matches the recorded id.
        assert_eq!(
            after.assignments.get(Segment::V).unwrap().allele_id.index(),
            recorded_id,
        );
        // Pool length unchanged (same-length constraint).
        assert_eq!(after.pool.len(), 6);
        // The replacement bytes equal the recorded allele's
        // 6-byte prefix.
        let bases: Vec<u8> = after.pool.as_slice().iter().map(|n| n.base).collect();
        assert_eq!(bases, recorded_allele.seq[..6].to_vec());
    }

    // ── Replay ──────────────────────────────────────────────────

    #[test]
    fn replay_applied_true_reproduces_replacement_without_rng() {
        let (cfg, v0, v1) = two_v_refdata();
        // prob=0.0 would normally fire applied=false; the trace
        // overrides via cursor. Pins "trace is the source of truth".
        let pass = ReceptorRevisionPass::new(0.0, Box::new(AllelePoolDist::uniform(&cfg.v_pool)));

        let mut input = Trace::new();
        input.record_choice(
            address::ChoiceAddress::ReceptorRevisionApplied,
            ChoiceValue::Bool(true),
        );
        input.record_choice(
            address::ChoiceAddress::ReceptorRevisionVAllele,
            ChoiceValue::AlleleId(v1.index()),
        );
        // V1 length 8, old V length 6 → trim_3 = 2.
        input.record_choice(
            address::ChoiceAddress::ReceptorRevisionVTrim3,
            ChoiceValue::Int(2),
        );
        let mut cursor = TraceCursor::from_owned(input.choices().to_vec());

        let (trace, after) =
            run_with_ctx(&pass, &cfg, None, sim_v_assembled(v0), Some(&mut cursor), None).unwrap();

        assert_eq!(
            trace.find("receptor_revision.applied").unwrap().value,
            ChoiceValue::Bool(true),
        );
        assert_eq!(
            after.assignments.get(Segment::V).unwrap().allele_id,
            v1,
        );
        assert_eq!(
            after.assignments.get(Segment::V).unwrap().trim_3,
            2,
        );
        // V1's retained 6-byte prefix is "GGGGGG".
        let bases: Vec<u8> = after.pool.as_slice().iter().map(|n| n.base).collect();
        assert_eq!(&bases, b"GGGGGG");
        assert!(cursor.is_drained());
    }

    #[test]
    fn replay_applied_false_consumes_only_bool_record() {
        let (cfg, v0, _) = two_v_refdata();
        // prob=1.0 would normally fire applied=true; trace says
        // false → no allele/trim consumption.
        let pass = ReceptorRevisionPass::new(1.0, Box::new(AllelePoolDist::uniform(&cfg.v_pool)));

        let mut input = Trace::new();
        input.record_choice(
            address::ChoiceAddress::ReceptorRevisionApplied,
            ChoiceValue::Bool(false),
        );
        let mut cursor = TraceCursor::from_owned(input.choices().to_vec());

        let (trace, after) =
            run_with_ctx(&pass, &cfg, None, sim_v_assembled(v0), Some(&mut cursor), None).unwrap();

        assert_eq!(
            trace.find("receptor_revision.applied").unwrap().value,
            ChoiceValue::Bool(false),
        );
        assert!(trace.find("receptor_revision.v_allele").is_none());
        // Unchanged.
        assert_eq!(after.assignments.get(Segment::V).unwrap().allele_id, v0);
        assert!(cursor.is_drained());
    }

    #[test]
    fn replay_missing_allele_record_after_applied_true_errors() {
        let (cfg, v0, _) = two_v_refdata();
        let pass = ReceptorRevisionPass::new(0.0, Box::new(AllelePoolDist::uniform(&cfg.v_pool)));

        let mut input = Trace::new();
        input.record_choice(
            address::ChoiceAddress::ReceptorRevisionApplied,
            ChoiceValue::Bool(true),
        );
        // Missing allele + trim records.
        let mut cursor = TraceCursor::from_owned(input.choices().to_vec());

        let err = run_with_ctx(&pass, &cfg, None, sim_v_assembled(v0), Some(&mut cursor), None)
            .unwrap_err();
        match err {
            PassError::Replay { pass_name, reason } => {
                assert_eq!(pass_name, "receptor_revision");
                let msg = format!("{reason}");
                assert!(
                    msg.contains("receptor_revision.v_allele") || msg.contains("exhausted"),
                    "expected replay error to mention v_allele or exhaustion, got: {msg}",
                );
            }
            other => panic!("expected PassError::Replay, got {other:?}"),
        }
    }

    #[test]
    fn replay_unknown_allele_id_errors() {
        let (cfg, v0, _) = two_v_refdata();
        let pass = ReceptorRevisionPass::new(0.0, Box::new(AllelePoolDist::uniform(&cfg.v_pool)));

        let mut input = Trace::new();
        input.record_choice(
            address::ChoiceAddress::ReceptorRevisionApplied,
            ChoiceValue::Bool(true),
        );
        input.record_choice(
            address::ChoiceAddress::ReceptorRevisionVAllele,
            ChoiceValue::AlleleId(999),
        );
        input.record_choice(
            address::ChoiceAddress::ReceptorRevisionVTrim3,
            ChoiceValue::Int(0),
        );
        let mut cursor = TraceCursor::from_owned(input.choices().to_vec());

        let err = run_with_ctx(&pass, &cfg, None, sim_v_assembled(v0), Some(&mut cursor), None)
            .unwrap_err();
        match err {
            PassError::MissingAllele { pass_name, segment, allele_id } => {
                assert_eq!(pass_name, "receptor_revision");
                assert_eq!(segment, Segment::V);
                assert_eq!(allele_id, 999);
            }
            other => panic!("expected PassError::MissingAllele, got {other:?}"),
        }
    }

    #[test]
    fn replay_trim_causing_length_mismatch_errors() {
        let (cfg, v0, v1) = two_v_refdata();
        let pass = ReceptorRevisionPass::new(0.0, Box::new(AllelePoolDist::uniform(&cfg.v_pool)));

        let mut input = Trace::new();
        input.record_choice(
            address::ChoiceAddress::ReceptorRevisionApplied,
            ChoiceValue::Bool(true),
        );
        input.record_choice(
            address::ChoiceAddress::ReceptorRevisionVAllele,
            ChoiceValue::AlleleId(v1.index()),
        );
        // V1 has length 8, old V region is 6. A trim_3 of 3 would
        // retain 5 bytes — mismatch.
        input.record_choice(
            address::ChoiceAddress::ReceptorRevisionVTrim3,
            ChoiceValue::Int(3),
        );
        let mut cursor = TraceCursor::from_owned(input.choices().to_vec());

        let err = run_with_ctx(&pass, &cfg, None, sim_v_assembled(v0), Some(&mut cursor), None)
            .unwrap_err();
        match err {
            PassError::InvalidPlanState { pass_name, reason } => {
                assert_eq!(pass_name, "receptor_revision");
                assert!(reason.contains("length mismatch"), "got: {reason}");
            }
            other => panic!("expected PassError::InvalidPlanState, got {other:?}"),
        }
    }

    // ── Plan-state guards ───────────────────────────────────────

    #[test]
    fn missing_v_assignment_errors_in_checked_path() {
        let (cfg, _, _) = two_v_refdata();
        let pass = ReceptorRevisionPass::new(0.5, Box::new(AllelePoolDist::uniform(&cfg.v_pool)));
        let err = run_with_ctx(&pass, &cfg, None, Simulation::new(), None, None).unwrap_err();
        match err {
            PassError::MissingAssignment { pass_name, segment } => {
                assert_eq!(pass_name, "receptor_revision");
                assert_eq!(segment, Segment::V);
            }
            other => panic!("expected MissingAssignment, got {other:?}"),
        }
    }

    #[test]
    fn missing_refdata_errors_in_checked_path() {
        // Build a sim with V assigned but execute without refdata.
        let (cfg, v0, _) = two_v_refdata();
        let pass = ReceptorRevisionPass::new(0.5, Box::new(AllelePoolDist::uniform(&cfg.v_pool)));

        let mut trace = Trace::new();
        let mut rng = Rng::new(0);
        let mut ctx = PassContext {
            trace: &mut trace,
            rng: &mut rng,
            pass_index: 0,
            refdata: None,
            contracts: None,
            feasibility: None,
            reference_index: None,
            replay_cursor: None,
            event_log_sink: None,
        };
        let err = pass
            .execute_checked(&sim_v_assembled(v0), &mut ctx)
            .unwrap_err();
        match err {
            PassError::MissingRefData { pass_name } => {
                assert_eq!(pass_name, "receptor_revision");
            }
            other => panic!("expected MissingRefData, got {other:?}"),
        }
    }

    #[test]
    fn no_v_region_errors_in_checked_path() {
        let (cfg, v0, _) = two_v_refdata();
        let pass = ReceptorRevisionPass::new(0.5, Box::new(AllelePoolDist::uniform(&cfg.v_pool)));
        // Assignment without a region — assembly never ran.
        let sim = Simulation::new().with_allele_assigned(Segment::V, AlleleInstance::new(v0));
        let err = run_with_ctx(&pass, &cfg, None, sim, None, None).unwrap_err();
        match err {
            PassError::InvalidPlanState { pass_name, reason } => {
                assert_eq!(pass_name, "receptor_revision");
                assert!(reason.contains("no V region"), "got: {reason}");
            }
            other => panic!("expected InvalidPlanState, got {other:?}"),
        }
    }

    // ── Pass metadata ───────────────────────────────────────────

    #[test]
    fn declares_three_choice_patterns_in_order() {
        let (cfg, _, _) = two_v_refdata();
        let pass = ReceptorRevisionPass::new(0.5, Box::new(AllelePoolDist::uniform(&cfg.v_pool)));
        assert_eq!(
            pass.declared_choice_patterns(),
            vec![
                address::ChoiceAddressPattern::ReceptorRevisionApplied,
                address::ChoiceAddressPattern::ReceptorRevisionVAllele,
                address::ChoiceAddressPattern::ReceptorRevisionVTrim3,
            ],
        );
    }

    #[test]
    fn declares_refdata_and_v_assignment_requirements() {
        let (cfg, _, _) = two_v_refdata();
        let pass = ReceptorRevisionPass::new(0.5, Box::new(AllelePoolDist::uniform(&cfg.v_pool)));
        assert_eq!(
            pass.requirements(),
            vec![
                PassRequirement::RefData,
                PassRequirement::AlleleAssignment(Segment::V),
            ],
        );
    }

    #[test]
    fn declares_no_compile_effects() {
        // Pinning the empty declaration — receptor revision is
        // event-driven, the schedule analyser must not treat it
        // as initial recombination.
        let (cfg, _, _) = two_v_refdata();
        let pass = ReceptorRevisionPass::new(0.5, Box::new(AllelePoolDist::uniform(&cfg.v_pool)));
        assert!(pass.effects().is_empty());
    }

    #[test]
    fn pass_name_is_stable() {
        // The frozen pass name flows through `pass_plan_signature`
        // — a rename here breaks every existing trace.
        let (cfg, _, _) = two_v_refdata();
        let pass = ReceptorRevisionPass::new(0.5, Box::new(AllelePoolDist::uniform(&cfg.v_pool)));
        assert_eq!(pass.name(), "receptor_revision");
    }

    // ── Full plan round-trip (V assembled → replacement) ────────

    #[test]
    fn full_runtime_round_trip_with_prob_one() {
        let (cfg, v0, _) = two_v_refdata();
        let (trace, after) = run_pass(1.0, 0, cfg.clone(), sim_v_assembled(v0));

        // Three records present.
        assert!(trace.find("receptor_revision.applied").is_some());
        assert!(trace.find("receptor_revision.v_allele").is_some());
        assert!(trace.find("receptor_revision.v_trim_3").is_some());

        // V region length unchanged (same-length constraint).
        let v_region = after
            .sequence
            .regions
            .iter()
            .find(|r| r.segment == Segment::V)
            .expect("V region must remain after replacement");
        assert_eq!(v_region.len(), 6);
        assert_eq!(after.pool.len(), 6);
    }
