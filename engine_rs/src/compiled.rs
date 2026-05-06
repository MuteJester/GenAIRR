//! Compiled simulator boundary.
//!
//! `PassPlan` is the engine's intermediate representation: an ordered
//! sequence of concrete simulation steps. It should not be the final
//! execution artifact. A `CompiledSimulator` binds that IR to the
//! reference data, active contracts, execution policy, and compile-time
//! analysis report needed to run safely and introspectably.

use crate::contract::ContractSet;
use crate::event::{EventRecord, StateSummary, TraceSpan};
use crate::ir::{Segment, Simulation};
use crate::pass::{Outcome, PassContext, PassEffect, PassError, PassPlan, PassRequirement};
use crate::refdata::RefDataConfig;
use crate::rng::Rng;
use crate::trace::Trace;
use std::collections::HashSet;
use std::fmt;

/// Runtime failure policy selected at compile time.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ExecutionPolicy {
    /// Preserve permissive runtime behavior: constrained draws may
    /// fall back to unconstrained sampling when no support is available.
    Permissive,
    /// Fail loudly with structured errors when a pass cannot execute
    /// safely or cannot sample an admissible candidate.
    Strict,
}

/// One pass-level compile summary.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PassSummary {
    pub index: usize,
    pub name: String,
    pub declared_choices: Vec<String>,
    pub requirements: Vec<PassRequirement>,
    pub effects: Vec<PassEffect>,
}

/// One declared stochastic choice, annotated with the pass that owns it.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct DeclaredChoice {
    pub pass_index: usize,
    pub pass_name: String,
    pub address: String,
}

/// Compile-time analysis output for a plan.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct CompileReport {
    pub policy: ExecutionPolicy,
    pub pass_summaries: Vec<PassSummary>,
    pub declared_choices: Vec<DeclaredChoice>,
    pub active_contracts: Vec<String>,
    pub warnings: Vec<CompileWarning>,
}

impl CompileReport {
    pub fn pass_names(&self) -> Vec<String> {
        self.pass_summaries
            .iter()
            .map(|summary| summary.name.clone())
            .collect()
    }
}

/// Non-fatal compile-time diagnostic.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct CompileWarning {
    pub pass_index: Option<usize>,
    pub pass_name: Option<String>,
    pub message: String,
}

/// Fatal compile-time error.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct CompileError {
    pub pass_index: usize,
    pub pass_name: String,
    pub kind: CompileErrorKind,
}

/// Stable class of compile-time plan error.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum CompileErrorKind {
    MissingRefData,
    MissingAssignment { segment: Segment },
}

/// All fatal compile-time errors found in one plan scan.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct CompileErrors {
    pub errors: Vec<CompileError>,
}

impl CompileErrors {
    pub fn is_empty(&self) -> bool {
        self.errors.is_empty()
    }
}

impl fmt::Display for CompileErrors {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.errors.is_empty() {
            return write!(f, "plan compilation failed with no errors");
        }

        write!(f, "plan compilation failed")?;
        for error in &self.errors {
            match &error.kind {
                CompileErrorKind::MissingRefData => write!(
                    f,
                    "; pass {} ({}) requires reference data",
                    error.pass_index, error.pass_name
                )?,
                CompileErrorKind::MissingAssignment { segment } => write!(
                    f,
                    "; pass {} ({}) requires a prior {:?} allele assignment",
                    error.pass_index, error.pass_name, segment
                )?,
            }
        }
        Ok(())
    }
}

impl std::error::Error for CompileErrors {}

/// A compiled simulation artifact ready for deterministic execution.
///
/// This first implementation borrows the plan/refdata/contracts from
/// the owning layer. That keeps the current Python API stable while
/// moving all execution through a real compile boundary. A future
/// owning `PyCompiledSimulator` can reuse this same analysis surface.
pub struct CompiledSimulator<'a> {
    plan: &'a PassPlan,
    refdata: Option<&'a RefDataConfig>,
    contracts: Option<&'a ContractSet>,
    policy: ExecutionPolicy,
    report: CompileReport,
}

#[derive(Debug)]
struct ExecutionAbort {
    error: PassError,
    /// Already-committed prefix at the point of failure. Kept private:
    /// public runners expose only the structured pass error, while
    /// internal tests can assert transaction commit boundaries.
    #[cfg_attr(not(test), allow(dead_code))]
    committed: Outcome,
}

impl<'a> CompiledSimulator<'a> {
    /// Compile a pass plan into an executable simulator.
    ///
    /// Validation is intentionally semantic and typed: it consumes
    /// `Pass::requirements()` / `Pass::effects()` rather than parsing
    /// pass names. This is the seam where deeper contract precondition
    /// checks and resource validation should be added.
    pub fn compile(
        plan: &'a PassPlan,
        refdata: Option<&'a RefDataConfig>,
        contracts: Option<&'a ContractSet>,
        policy: ExecutionPolicy,
    ) -> Result<Self, CompileErrors> {
        let (report, errors) = analyze_plan(plan, refdata, contracts, policy);
        if !errors.is_empty() {
            return Err(CompileErrors { errors });
        }

        Ok(Self {
            plan,
            refdata,
            contracts,
            policy,
            report,
        })
    }

    pub fn report(&self) -> &CompileReport {
        &self.report
    }

    pub fn policy(&self) -> ExecutionPolicy {
        self.policy
    }

    pub fn plan(&self) -> &PassPlan {
        self.plan
    }

    pub fn refdata(&self) -> Option<&RefDataConfig> {
        self.refdata
    }

    pub fn contracts(&self) -> Option<&ContractSet> {
        self.contracts
    }

    /// Run one simulation from an empty initial IR.
    pub fn run_one(&self, seed: u64) -> Result<Outcome, PassError> {
        self.run_one_from(Simulation::new(), seed)
    }

    /// Run one simulation from a caller-supplied initial IR.
    pub fn run_one_from(&self, initial: Simulation, seed: u64) -> Result<Outcome, PassError> {
        self.execute_transactional(initial, seed)
            .map_err(|abort| abort.error)
    }

    /// Run a deterministic seed-stitched batch. The i-th simulation
    /// uses `seed + i`, matching the current Python `CompiledExperiment`
    /// contract. Parallel seed splitting can replace this behind the
    /// same method later.
    pub fn run_batch(&self, n: usize, seed: u64) -> Vec<Result<Outcome, PassError>> {
        (0..n)
            .map(|i| self.run_one(seed.wrapping_add(i as u64)))
            .collect()
    }

    /// Compiled execution owns the step lifecycle so the compiled
    /// artifact, not `PassPlan`, is the authority that decides when
    /// trace deltas and state revisions are safe to commit.
    fn execute_transactional(
        &self,
        initial: Simulation,
        seed: u64,
    ) -> Result<Outcome, ExecutionAbort> {
        let mut trace = Trace::new();
        let mut rng = Rng::new(seed);

        let mut revisions: Vec<Simulation> = Vec::with_capacity(self.plan.len() + 1);
        let mut pass_names: Vec<String> = Vec::with_capacity(self.plan.len());
        let mut events: Vec<EventRecord> = Vec::with_capacity(self.plan.len());

        revisions.push(initial);

        for pass in self.plan.passes() {
            let prev = revisions.last().expect("non-empty by construction");
            let pass_index = pass_names.len();
            let pass_name = pass.name().to_string();
            let effects = pass.effects();
            let trace_start = trace.len();
            let pre = StateSummary::from_simulation(prev);
            let mut trace_delta = Trace::new();
            let mut ctx = PassContext {
                trace: &mut trace_delta,
                rng: &mut rng,
                refdata: self.refdata,
                contracts: self.contracts,
            };

            let next = match self.policy {
                ExecutionPolicy::Permissive => pass.execute(prev, &mut ctx),
                ExecutionPolicy::Strict => match pass.execute_checked(prev, &mut ctx) {
                    Ok(next) => next,
                    Err(error) => {
                        return Err(Self::abort(error, revisions, pass_names, trace, events));
                    }
                },
            };
            drop(ctx);

            if self.policy == ExecutionPolicy::Strict {
                if let Some(contracts) = self.contracts {
                    if let Err(violations) = contracts.verify(&next, self.refdata) {
                        return Err(Self::abort(
                            PassError::contract_violation(pass_name, violations),
                            revisions,
                            pass_names,
                            trace,
                            events,
                        ));
                    }
                }
            }

            let trace_end = trace_start + trace_delta.len();
            let post = StateSummary::from_simulation(&next);
            let event = EventRecord::pass_committed(
                pass_index,
                pass_name.clone(),
                effects,
                TraceSpan::new(trace_start, trace_end),
                pre,
                post,
            );

            trace.append_delta(trace_delta);
            pass_names.push(pass_name);
            revisions.push(next);
            events.push(event);
        }

        Ok(Outcome {
            revisions,
            pass_names,
            trace,
            events,
        })
    }

    fn abort(
        error: PassError,
        revisions: Vec<Simulation>,
        pass_names: Vec<String>,
        trace: Trace,
        events: Vec<EventRecord>,
    ) -> ExecutionAbort {
        ExecutionAbort {
            error,
            committed: Outcome {
                revisions,
                pass_names,
                trace,
                events,
            },
        }
    }
}

fn analyze_plan(
    plan: &PassPlan,
    refdata: Option<&RefDataConfig>,
    contracts: Option<&ContractSet>,
    policy: ExecutionPolicy,
) -> (CompileReport, Vec<CompileError>) {
    let mut assigned: HashSet<Segment> = HashSet::new();
    let mut pass_summaries = Vec::with_capacity(plan.len());
    let mut declared_choices = Vec::new();
    let mut errors = Vec::new();

    for (index, pass) in plan.passes().iter().enumerate() {
        let name = pass.name().to_string();
        let choices = pass.declared_choices();
        let requirements = pass.requirements();
        let effects = pass.effects();

        for req in &requirements {
            match req {
                PassRequirement::RefData if refdata.is_none() => {
                    errors.push(CompileError {
                        pass_index: index,
                        pass_name: name.clone(),
                        kind: CompileErrorKind::MissingRefData,
                    });
                }
                PassRequirement::AlleleAssignment(segment) if !assigned.contains(segment) => {
                    errors.push(CompileError {
                        pass_index: index,
                        pass_name: name.clone(),
                        kind: CompileErrorKind::MissingAssignment { segment: *segment },
                    });
                }
                _ => {}
            }
        }

        for address in &choices {
            declared_choices.push(DeclaredChoice {
                pass_index: index,
                pass_name: name.clone(),
                address: address.clone(),
            });
        }

        for effect in &effects {
            if let PassEffect::AssignAllele(segment) = effect {
                assigned.insert(*segment);
            }
        }

        pass_summaries.push(PassSummary {
            index,
            name,
            declared_choices: choices,
            requirements,
            effects,
        });
    }

    let active_contracts = contracts
        .map(|set| {
            set.iter()
                .map(|contract| contract.name().to_string())
                .collect()
        })
        .unwrap_or_default();

    (
        CompileReport {
            policy,
            pass_summaries,
            declared_choices,
            active_contracts,
            warnings: Vec::new(),
        },
        errors,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::assignment::TrimEnd;
    use crate::contract::{Contract, ContractViolation};
    use crate::dist::{AllelePoolDist, EmpiricalLengthDist, UniformBase};
    use crate::ir::{NucFlags, Nucleotide, Segment};
    use crate::pass::{Pass, PassRuntime};
    use crate::passes::{AssembleSegmentPass, GenerateNPPass, SampleAllelePass, TrimPass};
    use crate::refdata::{Allele, ChainType};
    use crate::trace::ChoiceValue;

    fn vj_refdata() -> RefDataConfig {
        let mut cfg = RefDataConfig::empty(ChainType::Vj);
        let _ = cfg.v_pool.push(Allele {
            name: "v*01".into(),
            gene: "v".into(),
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
        cfg
    }

    fn vj_plan(cfg: &RefDataConfig) -> PassPlan {
        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(&cfg.v_pool)),
        )));
        plan.push(Box::new(SampleAllelePass::new(
            Segment::J,
            Box::new(AllelePoolDist::uniform(&cfg.j_pool)),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
        plan.push(Box::new(GenerateNPPass::new(
            Segment::Np1,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(0, 1.0)])),
            Box::new(UniformBase),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::J)));
        plan
    }

    struct AppendBasePass {
        name: &'static str,
        base: u8,
    }

    impl AppendBasePass {
        fn new(name: &'static str, base: u8) -> Self {
            Self { name, base }
        }
    }

    impl Pass for AppendBasePass {
        fn name(&self) -> &str {
            self.name
        }

        fn execute(&self, sim: &Simulation, _ctx: &mut PassContext) -> Simulation {
            let (next, _handle) = sim.with_nucleotide_pushed(Nucleotide::synthetic(
                self.base,
                Segment::Np1,
                NucFlags::empty(),
            ));
            next
        }

        fn effects(&self) -> Vec<PassEffect> {
            vec![PassEffect::AppendNucleotides]
        }
    }

    struct TraceProbePass {
        name: &'static str,
        address: &'static str,
        base: u8,
    }

    impl TraceProbePass {
        fn new(name: &'static str, address: &'static str, base: u8) -> Self {
            Self {
                name,
                address,
                base,
            }
        }
    }

    impl Pass for TraceProbePass {
        fn name(&self) -> &str {
            self.name
        }

        fn execute(&self, sim: &Simulation, ctx: &mut PassContext) -> Simulation {
            assert!(
                ctx.trace.is_empty(),
                "{} saw previously committed trace records",
                self.name
            );
            ctx.trace.record(self.address, ChoiceValue::Base(self.base));
            let (next, _handle) = sim.with_nucleotide_pushed(Nucleotide::synthetic(
                self.base,
                Segment::Np1,
                NucFlags::empty(),
            ));
            next
        }

        fn declared_choices(&self) -> Vec<String> {
            vec![self.address.to_string()]
        }

        fn effects(&self) -> Vec<PassEffect> {
            vec![PassEffect::AppendNucleotides]
        }
    }

    struct MaxPoolLen {
        max: usize,
    }

    impl MaxPoolLen {
        fn new(max: usize) -> Self {
            Self { max }
        }
    }

    impl Contract for MaxPoolLen {
        fn name(&self) -> &str {
            "test.max_pool_len"
        }

        fn verify(
            &self,
            sim: &Simulation,
            _refdata: Option<&RefDataConfig>,
        ) -> Result<(), ContractViolation> {
            if sim.pool.len() <= self.max {
                Ok(())
            } else {
                Err(ContractViolation::new(
                    self.name(),
                    format!("pool length {} exceeds {}", sim.pool.len(), self.max),
                ))
            }
        }
    }

    #[test]
    fn compiled_simulator_collects_report_metadata() {
        let cfg = vj_refdata();
        let plan = vj_plan(&cfg);
        let contracts = crate::contract::productive();
        let compiled = CompiledSimulator::compile(
            &plan,
            Some(&cfg),
            Some(&contracts),
            ExecutionPolicy::Strict,
        )
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
    fn compiled_run_one_matches_direct_runtime_for_valid_plan() {
        let cfg = vj_refdata();
        let plan = vj_plan(&cfg);
        let compiled =
            CompiledSimulator::compile(&plan, Some(&cfg), None, ExecutionPolicy::Permissive)
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
}
