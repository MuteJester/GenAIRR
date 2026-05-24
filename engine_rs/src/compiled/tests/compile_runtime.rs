use super::*;

#[test]
fn compiled_simulator_collects_report_metadata() {
    let cfg = vj_refdata();
    let plan = vj_plan(&cfg);
    let contracts = crate::contract::productive();
    let compiled =
        CompiledSimulator::compile(&plan, Some(&cfg), Some(&contracts), ExecutionPolicy::Strict)
            .expect("plan should compile");

    assert_eq!(compiled.policy(), ExecutionPolicy::Strict);
    assert_eq!(
        compiled.report().pass_names(),
        vec![
            "sample_allele.v",
            "sample_allele.j",
            "assemble.v",
            "generate_np.np1",
            "assemble.j",
        ]
    );
    assert!(compiled
        .report()
        .declared_choices
        .iter()
        .any(|choice| choice.address == "np.np1.length"));
    assert_eq!(
        compiled.report().active_contracts,
        vec![
            "productive_junction_frame",
            "no_stop_codon_in_junction",
            "anchor_preserved.v",
            "anchor_preserved.j",
        ]
    );
}

#[test]
fn compiled_simulators_own_reference_match_index_when_refdata_exists() {
    let cfg = vj_refdata();
    let plan = vj_plan(&cfg);
    let compiled = CompiledSimulator::compile(&plan, Some(&cfg), None, ExecutionPolicy::Permissive)
        .expect("plan should compile");

    let index = compiled
        .reference_index()
        .expect("refdata-backed compile should build a reference index");
    assert_eq!(index.v.allele_count(), cfg.v_pool.len());
    assert_eq!(index.j.allele_count(), cfg.j_pool.len());
    assert_eq!(index.d.allele_count(), 0);
    assert_eq!(
        index.v.kmer_hits(b"AAACCCG").unwrap(),
        &[crate::live_call::KmerHit {
            allele_id: AlleleId::new(0),
            ref_pos: 0,
        }]
    );

    let owned = OwnedCompiledSimulator::compile(plan, Some(cfg), None, ExecutionPolicy::Permissive)
        .expect("owned compile should also build the index");
    assert!(owned.reference_index().is_some());

    let empty_plan = PassPlan::new();
    let no_refdata =
        CompiledSimulator::compile(&empty_plan, None, None, ExecutionPolicy::Permissive)
            .expect("empty plan without refdata should compile");
    assert!(no_refdata.reference_index().is_none());
}

#[test]
fn compiled_assembly_initializes_exact_live_call_from_current_region() {
    let mut cfg = RefDataConfig::empty(ChainType::Vdj);
    let v0 = cfg.v_pool.push(Allele {
        name: "v1*01".into(),
        gene: "v1".into(),
        seq: b"GGAAACCC".to_vec(),
        segment: Segment::V,
        anchor: None,
    });
    let v1 = cfg.v_pool.push(Allele {
        name: "v2*01".into(),
        gene: "v2".into(),
        seq: b"TTAAACCC".to_vec(),
        segment: Segment::V,
        anchor: None,
    });

    let mut plan = PassPlan::new();
    plan.push(Box::new(SampleAllelePass::new(
        Segment::V,
        Box::new(AllelePoolDist::restricted_uniform(&cfg.v_pool, vec![v0])),
    )));
    plan.push(Box::new(TrimPass::new(
        Segment::V,
        TrimEnd::Five,
        Box::new(EmpiricalLengthDist::from_pairs(vec![(2, 1.0)])),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));

    let compiled = CompiledSimulator::compile(&plan, Some(&cfg), None, ExecutionPolicy::Permissive)
        .expect("plan should compile");
    let outcome = compiled.run_one(0).expect("run should succeed");
    let final_sim = outcome.final_simulation();
    let v_call = final_sim
        .segment_calls
        .get(Segment::V)
        .cloned()
        .expect("V live call should exist");

    assert_eq!(final_sim.assignments.get(Segment::V).unwrap().allele_id, v0);
    assert_eq!(v_call.allele_call.to_ids(), vec![v0, v1]);
    assert_eq!(
        v_call.confidence,
        crate::live_call::LiveCallConfidence::ExactAmbiguous
    );
    assert_eq!(
        v_call.boundary_summary.ref_start,
        crate::live_call::BoundaryValue::Single(2)
    );
    assert_eq!(
        v_call.boundary_summary.ref_end,
        crate::live_call::BoundaryValue::Single(8)
    );
    assert_eq!(
        v_call.boundary_summary.seq_start,
        crate::live_call::BoundaryValue::Single(0)
    );
    assert_eq!(
        v_call.boundary_summary.seq_end,
        crate::live_call::BoundaryValue::Single(6)
    );
}

#[test]
fn compiled_v_three_prime_trim_widens_live_and_airr_call() {
    let mut cfg = RefDataConfig::empty(ChainType::Vdj);
    let shared_prefix = b"AAACCCGGGTTT";
    let suffixes: [&[u8]; 5] = [
        b"AAAAAAAAAA",
        b"CCCCCCCCCC",
        b"GGGGGGGGGG",
        b"TTTTTTTTTT",
        b"ACGTACGTAC",
    ];
    let mut ids = Vec::new();
    for (index, suffix) in suffixes.iter().enumerate() {
        let mut seq = shared_prefix.to_vec();
        seq.extend_from_slice(suffix);
        ids.push(cfg.v_pool.push(Allele {
            name: format!("IGHVtrim*0{}", index + 1),
            gene: "IGHVtrim".into(),
            seq,
            segment: Segment::V,
            anchor: None,
        }));
    }

    let mut plan = PassPlan::new();
    plan.push(Box::new(SampleAllelePass::new(
        Segment::V,
        Box::new(AllelePoolDist::restricted_uniform(
            &cfg.v_pool,
            vec![ids[0]],
        )),
    )));
    plan.push(Box::new(TrimPass::new(
        Segment::V,
        TrimEnd::Three,
        Box::new(EmpiricalLengthDist::from_pairs(vec![(10, 1.0)])),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));

    let compiled = CompiledSimulator::compile(&plan, Some(&cfg), None, ExecutionPolicy::Permissive)
        .expect("trim-before-assembly plan should compile");
    let outcome = compiled.run_one(0).expect("run should succeed");
    let final_sim = outcome.final_simulation();
    let v_call = final_sim
        .segment_calls
        .get(Segment::V)
        .cloned()
        .expect("V live call exists");

    assert_eq!(
        final_sim.assignments.get(Segment::V).unwrap().allele_id,
        ids[0]
    );
    assert_eq!(v_call.allele_call.to_ids(), ids);
    assert_eq!(
        v_call.confidence,
        crate::live_call::LiveCallConfidence::ExactAmbiguous
    );
    assert_eq!(
        v_call.boundary_summary.ref_end,
        crate::live_call::BoundaryValue::Single(shared_prefix.len() as u32)
    );

    let rec = build_airr_record(&outcome, &cfg, "trim-ambiguity");
    assert_eq!(
        rec.v_call,
        "IGHVtrim*01,IGHVtrim*02,IGHVtrim*03,IGHVtrim*04,IGHVtrim*05"
    );
}

#[test]
fn compiled_simulator_auto_reorders_trim_pushed_after_assemble() {
    // Before the dep-graph scheduler, pushing trim AFTER assemble of
    // the same segment was a hard compile error. The scheduler now
    // derives the implicit `TrimAllele → AssembleSegment` edge from
    // the pass metadata and topo-sorts trim back ahead of assemble.
    let cfg = vj_refdata();
    let mut plan = PassPlan::new();
    plan.push(Box::new(SampleAllelePass::new(
        Segment::V,
        Box::new(AllelePoolDist::uniform(&cfg.v_pool)),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
    plan.push(Box::new(TrimPass::new(
        Segment::V,
        TrimEnd::Three,
        Box::new(EmpiricalLengthDist::from_pairs(vec![(1, 1.0)])),
    )));

    let compiled =
        CompiledSimulator::compile(&plan, Some(&cfg), None, ExecutionPolicy::Permissive)
            .expect("scheduler must auto-fix misordered trim/assemble");

    // Report's pass_names() reads execution order from the topo-sort;
    // trim should appear before assemble even though it was pushed
    // last.
    let names = compiled.report().pass_names();
    let pos = |needle: &str| {
        names
            .iter()
            .position(|n| n == needle)
            .unwrap_or_else(|| panic!("pass {:?} not in {:?}", needle, names))
    };
    assert!(
        pos("trim.v_3") < pos("assemble.v"),
        "trim.v_3 must run before assemble.v, got order {:?}",
        names
    );
}

#[test]
fn compiled_simulator_rejects_assemble_without_refdata() {
    let cfg = vj_refdata();
    let plan = vj_plan(&cfg);
    let err = match CompiledSimulator::compile(&plan, None, None, ExecutionPolicy::Permissive) {
        Ok(_) => panic!("assembly without refdata should fail at compile time"),
        Err(err) => err,
    };

    assert!(err
        .errors
        .iter()
        .any(|e| matches!(e.kind, CompileErrorKind::MissingRefData)));
}

#[test]
fn compiled_simulator_rejects_trim_before_sample_allele() {
    let mut plan = PassPlan::new();
    plan.push(Box::new(TrimPass::new(
        Segment::V,
        TrimEnd::Three,
        Box::new(EmpiricalLengthDist::from_pairs(vec![(1, 1.0)])),
    )));

    let err = match CompiledSimulator::compile(&plan, None, None, ExecutionPolicy::Permissive) {
        Ok(_) => panic!("trim before sample allele should fail"),
        Err(err) => err,
    };
    assert_eq!(
        err.errors[0].kind,
        CompileErrorKind::MissingAssignment {
            segment: Segment::V
        }
    );
}

#[test]
fn compiled_simulator_rejects_negative_np_length_support() {
    let cfg = vj_refdata();
    let plan = vj_plan_with_np_lengths(&cfg, vec![(-1, 1.0), (0, 1.0)]);

    let err = match CompiledSimulator::compile(&plan, Some(&cfg), None, ExecutionPolicy::Permissive)
    {
        Ok(_) => panic!("negative NP length support should fail at compile time"),
        Err(err) => err,
    };

    assert!(err.errors.iter().any(|e| matches!(
        &e.kind,
        CompileErrorKind::InvalidParameterSupport { address, reason }
            if address == "np.np1.length" && reason == "negative_length"
    )));
}

#[test]
fn compiled_simulator_rejects_productive_when_v_support_has_no_anchor() {
    let mut cfg = RefDataConfig::empty(ChainType::Vj);
    let _ = cfg.v_pool.push(Allele {
        name: "v_anchorless*01".into(),
        gene: "v_anchorless".into(),
        seq: b"AAACCCGGG".to_vec(),
        segment: Segment::V,
        anchor: None,
    });
    let _ = cfg.j_pool.push(Allele {
        name: "j*01".into(),
        gene: "j".into(),
        seq: b"TTTAAA".to_vec(),
        segment: Segment::J,
        anchor: Some(0),
    });
    let plan = vj_plan_with_np_lengths(&cfg, vec![(0, 1.0), (3, 1.0)]);
    let contracts = crate::contract::productive();

    let err = match CompiledSimulator::compile(
        &plan,
        Some(&cfg),
        Some(&contracts),
        ExecutionPolicy::Permissive,
    ) {
        Ok(_) => panic!("productive compile should reject anchorless V support"),
        Err(err) => err,
    };

    assert!(err.errors.iter().any(|e| matches!(
        &e.kind,
        CompileErrorKind::ContractPrecondition {
            contract_name,
            reason,
        } if contract_name == "productive_junction_frame"
            && reason.contains("V sample support has no alleles with valid anchors")
    )));
}

#[test]
fn compiled_simulator_rejects_productive_without_in_frame_np_mass() {
    let cfg = vj_refdata();
    let plan = vj_plan_with_np_lengths(&cfg, vec![(1, 1.0), (2, 1.0)]);
    let contracts = crate::contract::productive();

    let err = match CompiledSimulator::compile(
        &plan,
        Some(&cfg),
        Some(&contracts),
        ExecutionPolicy::Permissive,
    ) {
        Ok(_) => panic!("productive compile should reject NP support with no in-frame mass"),
        Err(err) => err,
    };

    assert!(err.errors.iter().any(|e| matches!(
        &e.kind,
        CompileErrorKind::ContractPrecondition {
            contract_name,
            reason,
        } if contract_name == "productive_junction_frame"
            && reason.contains("NP1 length support has no in-frame mass")
    )));
}

#[test]
fn productive_runtime_filters_anchorless_v_allele_candidates() {
    let mut cfg = RefDataConfig::empty(ChainType::Vj);
    let _ = cfg.v_pool.push(Allele {
        name: "v_anchorless*01".into(),
        gene: "v_anchorless".into(),
        seq: b"AAACCCGGG".to_vec(),
        segment: Segment::V,
        anchor: None,
    });
    let anchored_v = cfg.v_pool.push(Allele {
        name: "v_anchored*01".into(),
        gene: "v_anchored".into(),
        seq: b"AAACCCGGG".to_vec(),
        segment: Segment::V,
        anchor: Some(6),
    });
    let _ = cfg.j_pool.push(Allele {
        name: "j*01".into(),
        gene: "j".into(),
        seq: b"TTTAAA".to_vec(),
        segment: Segment::J,
        anchor: Some(0),
    });
    let plan = vj_plan(&cfg);
    let contracts = crate::contract::productive();
    let compiled =
        CompiledSimulator::compile(&plan, Some(&cfg), Some(&contracts), ExecutionPolicy::Strict)
            .expect("mixed anchored/anchorless support should compile");

    let outcome = compiled.run_one(0).expect("run should succeed");

    assert_eq!(
        outcome.final_simulation().assignments.get(Segment::V).copied().unwrap().allele_id,
        anchored_v
    );
}

#[test]
fn vj_productive_feasibility_filters_j_before_np_length_sampling() {
    let (cfg, j_good, j_bad) = vj_refdata_with_runtime_j_residue();
    let plan = vj_plan_with_np_lengths(&cfg, vec![(0, 1.0)]);
    let contracts = crate::contract::productive();
    let compiled =
        CompiledSimulator::compile(&plan, Some(&cfg), Some(&contracts), ExecutionPolicy::Strict)
            .expect("mixed J support has one feasible completion");

    for seed in 0..32 {
        let outcome = compiled.run_one(seed).expect("run should stay feasible");
        let final_sim = outcome.final_simulation();
        assert_eq!(final_sim.assignments.get(Segment::J).copied().unwrap().allele_id, j_good);
        assert_ne!(final_sim.assignments.get(Segment::J).copied().unwrap().allele_id, j_bad);
        assert!(
            contracts.verify(final_sim, Some(&cfg)).is_ok(),
            "seed {seed} should satisfy productive()"
        );
    }
}

#[test]
fn vj_productive_feasibility_respects_future_trim_support() {
    let mut cfg = RefDataConfig::empty(ChainType::Vj);
    let _ = cfg.v_pool.push(Allele {
        name: "v*01".into(),
        gene: "v".into(),
        seq: b"AAACCCGGGAAA".to_vec(),
        segment: Segment::V,
        anchor: Some(6),
    });
    let j_good = cfg.j_pool.push(Allele {
        name: "j_good*01".into(),
        gene: "j_good".into(),
        seq: b"TTTAAA".to_vec(),
        segment: Segment::J,
        anchor: Some(0),
    });
    let j_needs_future_v_trim = cfg.j_pool.push(Allele {
        name: "j_future_trim*01".into(),
        gene: "j_future_trim".into(),
        seq: b"TTTAAA".to_vec(),
        segment: Segment::J,
        anchor: Some(1),
    });

    let mut plan = PassPlan::new();
    plan.push(Box::new(SampleAllelePass::new(
        Segment::V,
        Box::new(AllelePoolDist::uniform(&cfg.v_pool)),
    )));
    plan.push(Box::new(SampleAllelePass::new(
        Segment::J,
        Box::new(AllelePoolDist::uniform(&cfg.j_pool)),
    )));
    plan.push(Box::new(TrimPass::new(
        Segment::V,
        TrimEnd::Three,
        Box::new(EmpiricalLengthDist::from_pairs(vec![(1, 1.0)])),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
    plan.push(Box::new(GenerateNPPass::new(
        Segment::Np1,
        Box::new(EmpiricalLengthDist::from_pairs(vec![(0, 1.0)])),
        Box::new(UniformBase),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::J)));

    let contracts = crate::contract::productive();
    let compiled =
        CompiledSimulator::compile(&plan, Some(&cfg), Some(&contracts), ExecutionPolicy::Strict)
            .expect("future V trim support makes one J allele feasible");

    let outcome = compiled.run_one(0).expect("run should succeed");
    let final_sim = outcome.final_simulation();

    assert_eq!(
        final_sim.assignments.get(Segment::J).copied().unwrap().allele_id,
        j_needs_future_v_trim
    );
    assert_ne!(final_sim.assignments.get(Segment::J).copied().unwrap().allele_id, j_good);
    assert_eq!(final_sim.assignments.get(Segment::V).copied().unwrap().trim_3, 1);
    contracts
        .verify(final_sim, Some(&cfg))
        .expect("final simulation should satisfy productive()");
}

#[test]
fn vj_productive_feasibility_filters_upstream_known_stop() {
    let mut cfg = RefDataConfig::empty(ChainType::Vj);
    let _v_stop = cfg.v_pool.push(Allele {
        name: "v_stop*01".into(),
        gene: "v_stop".into(),
        seq: b"AAATAAGGG".to_vec(),
        segment: Segment::V,
        anchor: Some(3),
    });
    let v_productive = cfg.v_pool.push(Allele {
        name: "v_productive*01".into(),
        gene: "v_productive".into(),
        seq: b"AAATGTGGG".to_vec(),
        segment: Segment::V,
        anchor: Some(3),
    });
    let _ = cfg.j_pool.push(Allele {
        name: "j*01".into(),
        gene: "j".into(),
        seq: b"TGGAAA".to_vec(),
        segment: Segment::J,
        anchor: Some(0),
    });
    let plan = vj_plan_with_np_lengths(&cfg, vec![(0, 1.0)]);
    let contracts = crate::contract::productive();
    let compiled =
        CompiledSimulator::compile(&plan, Some(&cfg), Some(&contracts), ExecutionPolicy::Strict)
            .expect("mixed V support has one stop-free productive completion");

    for seed in 0..32 {
        let outcome = compiled.run_one(seed).expect("run should stay feasible");
        let final_sim = outcome.final_simulation();
        assert_eq!(final_sim.assignments.get(Segment::V).copied().unwrap().allele_id, v_productive);
        contracts
            .verify(final_sim, Some(&cfg))
            .unwrap_or_else(|_| panic!("seed {seed} should satisfy productive()"));
    }
}

#[test]
fn compiled_run_one_matches_direct_runtime_for_valid_plan() {
    let cfg = vj_refdata();
    let plan = vj_plan(&cfg);
    let compiled = CompiledSimulator::compile(&plan, Some(&cfg), None, ExecutionPolicy::Permissive)
        .expect("plan should compile");

    let direct = PassRuntime::execute_with_refdata(&plan, Simulation::new(), 42, &cfg);
    let via_compiled = compiled.run_one(42).expect("run should succeed");

    assert_eq!(direct.trace.choices(), via_compiled.trace.choices());
    assert!(direct.events.is_empty());
    assert_eq!(via_compiled.events.len(), plan.len());
    assert_eq!(
        direct.final_simulation().pool.as_slice(),
        via_compiled.final_simulation().pool.as_slice()
    );
}

#[test]
fn owned_compiled_simulator_owns_inputs_and_runs_repeatedly() {
    let cfg = vj_refdata();
    let plan = vj_plan(&cfg);
    let contracts = ContractSet::new().with(Box::new(MaxPoolLen::new(64)));
    let borrowed_report = {
        let borrowed = CompiledSimulator::compile(
            &plan,
            Some(&cfg),
            Some(&contracts),
            ExecutionPolicy::Strict,
        )
        .expect("borrowed simulator should compile");
        borrowed.report().clone()
    };

    let compiled = OwnedCompiledSimulator::compile(
        plan,
        Some(cfg.clone()),
        Some(contracts.clone()),
        ExecutionPolicy::Strict,
    )
    .expect("owned simulator should compile");

    let first = compiled.run_one(7).expect("first run should succeed");
    let second = compiled.run_one(7).expect("second run should succeed");

    assert_eq!(compiled.report(), &borrowed_report);
    assert_eq!(compiled.policy(), ExecutionPolicy::Strict);
    assert_eq!(compiled.plan().len(), 5);
    assert_eq!(compiled.refdata().unwrap().v_pool.len(), 1);
    assert_eq!(compiled.contracts().unwrap().len(), 1);
    assert_eq!(first.trace.choices(), second.trace.choices());
    assert_eq!(first.events.len(), compiled.plan().len());
    assert_eq!(second.events.len(), compiled.plan().len());
}

#[test]
fn compiled_execution_uses_per_pass_trace_deltas() {
    for policy in [ExecutionPolicy::Permissive, ExecutionPolicy::Strict] {
        let mut plan = PassPlan::new();
        plan.push(Box::new(TraceProbePass::new(
            "trace_probe_one",
            "trace.probe.one",
            b'A',
        )));
        plan.push(Box::new(TraceProbePass::new(
            "trace_probe_two",
            "trace.probe.two",
            b'C',
        )));
        let compiled =
            CompiledSimulator::compile(&plan, None, None, policy).expect("plan should compile");

        let outcome = compiled.run_one(0).expect("run should succeed");

        assert_eq!(outcome.trace.len(), 2);
        assert_eq!(outcome.trace.choices()[0].address, "trace.probe.one");
        assert_eq!(outcome.trace.choices()[1].address, "trace.probe.two");
        assert_eq!(outcome.events.len(), 2);
        assert_eq!(outcome.events[0].pass_index, 0);
        assert_eq!(outcome.events[0].pass_name, "trace_probe_one");
        assert_eq!(outcome.events[0].trace_span, TraceSpan::new(0, 1));
        assert_eq!(outcome.events[0].pre.pool_len, 0);
        assert_eq!(outcome.events[0].post.pool_len, 1);
        assert_eq!(outcome.events[1].pass_index, 1);
        assert_eq!(outcome.events[1].pass_name, "trace_probe_two");
        assert_eq!(outcome.events[1].trace_span, TraceSpan::new(1, 2));
        assert_eq!(outcome.events[1].pre.pool_len, 1);
        assert_eq!(outcome.events[1].post.pool_len, 2);
        assert_eq!(
            outcome.pass_names,
            vec!["trace_probe_one", "trace_probe_two"]
        );
        assert_eq!(outcome.final_simulation().pool.len(), 2);
    }
}

#[test]
fn strict_compiled_run_verifies_contracts_after_each_pass() {
    let mut plan = PassPlan::new();
    plan.push(Box::new(AppendBasePass::new("append_one", b'A')));
    plan.push(Box::new(AppendBasePass::new("append_two", b'C')));
    let contracts = ContractSet::new().with(Box::new(MaxPoolLen::new(1)));
    let compiled =
        CompiledSimulator::compile(&plan, None, Some(&contracts), ExecutionPolicy::Strict)
            .expect("plan should compile");

    let err = compiled
        .run_one(0)
        .expect_err("second pass should violate the strict contract fence");

    match err {
        PassError::ContractViolation {
            pass_name,
            violations,
        } => {
            assert_eq!(pass_name, "append_two");
            assert_eq!(violations.len(), 1);
            assert_eq!(violations[0].contract_name, "test.max_pool_len");
        }
        other => panic!("expected contract violation, got {other:?}"),
    }
}

#[test]
fn owned_compiled_simulator_allows_runtime_policy_override() {
    let mut plan = PassPlan::new();
    plan.push(Box::new(AppendBasePass::new("append_one", b'A')));
    plan.push(Box::new(AppendBasePass::new("append_two", b'C')));
    let contracts = ContractSet::new().with(Box::new(MaxPoolLen::new(1)));
    let compiled =
        OwnedCompiledSimulator::compile(plan, None, Some(contracts), ExecutionPolicy::Permissive)
            .expect("plan should compile");

    assert!(compiled.run_one(0).is_ok());
    let err = compiled
        .run_one_with_policy(0, ExecutionPolicy::Strict)
        .expect_err("strict override should enforce contract fence");

    assert!(matches!(err, PassError::ContractViolation { .. }));
}

#[test]
fn strict_abort_does_not_commit_rejected_pass_trace_or_state() {
    let mut plan = PassPlan::new();
    plan.push(Box::new(TraceProbePass::new(
        "append_one",
        "trace.append.one",
        b'A',
    )));
    plan.push(Box::new(TraceProbePass::new(
        "append_two",
        "trace.append.two",
        b'C',
    )));
    let contracts = ContractSet::new().with(Box::new(MaxPoolLen::new(1)));
    let compiled =
        CompiledSimulator::compile(&plan, None, Some(&contracts), ExecutionPolicy::Strict)
            .expect("plan should compile");

    let abort = compiled
        .execute_transactional(Simulation::new(), 0)
        .expect_err("second pass should abort before commit");

    match abort.error {
        PassError::ContractViolation {
            pass_name,
            violations,
        } => {
            assert_eq!(pass_name, "append_two");
            assert_eq!(violations[0].contract_name, "test.max_pool_len");
        }
        other => panic!("expected contract violation, got {other:?}"),
    }

    assert_eq!(abort.committed.pass_names, vec!["append_one"]);
    assert_eq!(abort.committed.revisions.len(), 2);
    assert_eq!(abort.committed.final_simulation().pool.len(), 1);
    assert_eq!(abort.committed.trace.len(), 1);
    assert_eq!(abort.committed.events.len(), 1);
    assert_eq!(abort.committed.events[0].pass_name, "append_one");
    assert_eq!(abort.committed.events[0].trace_span, TraceSpan::new(0, 1));
    assert!(abort.committed.trace.find("trace.append.one").is_some());
    assert!(abort.committed.trace.find("trace.append.two").is_none());
}
