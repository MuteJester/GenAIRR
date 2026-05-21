//! Compile-time plan analysis: validates declared supports, builds the
//! `CompileReport`, and derives the optional `FeasibilityContext`.
//!
//! The `CompileFactIndex` is the shared data structure that holds
//! supports gathered from each pass. `analyze_plan` is the single
//! entry point used by both `CompiledSimulator::compile` and
//! `OwnedCompiledSimulator::compile`.

use super::error::{
    invalid_parameter, invalid_pass_order, pass_scoped_error, plan_scoped_error, CompileError,
    CompileErrorKind,
};
use super::feasibility_builder::{
    anchor_head_residues, anchor_tail_residues, build_feasibility_context, d_length_residues,
    optional_np_length_residues, required_np_length_residues, residues_admit_frame,
};
use super::report::{CompileReport, DeclaredChoice, PassSummary};
use super::{np_length_address, sample_allele_address, trim_address, ExecutionPolicy};
use crate::assignment::TrimEnd;
use crate::contract::{ContractKind, ContractSet};
use crate::feasibility::FeasibilityContext;
use crate::ir::Segment;
use crate::pass::{
    AlleleIdSupport, IntegerSupport, PassCompileFact, PassEffect, PassPlan, PassRequirement,
};
use crate::refdata::RefDataConfig;
use std::collections::{HashMap, HashSet};

#[derive(Clone, Debug)]
pub(super) struct Located<T> {
    pub(super) value: T,
    pub(super) pass_index: usize,
    pub(super) pass_name: String,
}

impl<T> Located<T> {
    fn new(value: T, pass_index: usize, pass_name: &str) -> Self {
        Self {
            value,
            pass_index,
            pass_name: pass_name.to_string(),
        }
    }

    fn error(&self, kind: CompileErrorKind) -> CompileError {
        CompileError {
            pass_index: Some(self.pass_index),
            pass_name: Some(self.pass_name.clone()),
            kind,
        }
    }
}

#[derive(Default)]
pub(super) struct CompileFactIndex {
    pub(super) allele_supports: HashMap<Segment, Located<AlleleIdSupport>>,
    pub(super) trim_supports: HashMap<(Segment, TrimEnd), Located<IntegerSupport>>,
    pub(super) np_length_supports: HashMap<Segment, Located<IntegerSupport>>,
    pub(super) assigned_segments: HashSet<Segment>,
}

impl CompileFactIndex {
    fn record_fact(
        &mut self,
        fact: PassCompileFact,
        pass_index: usize,
        pass_name: &str,
        refdata: Option<&RefDataConfig>,
        errors: &mut Vec<CompileError>,
    ) {
        match fact {
            PassCompileFact::AlleleSampleSupport { segment, support } => {
                let located = Located::new(support, pass_index, pass_name);
                validate_allele_support(segment, &located, refdata, errors);
                self.allele_supports.insert(segment, located);
            }
            PassCompileFact::TrimSupport {
                segment,
                end,
                support,
            } => {
                let located = Located::new(support, pass_index, pass_name);
                validate_trim_support(segment, end, &located, errors);
                self.trim_supports.insert((segment, end), located);
            }
            PassCompileFact::NpLengthSupport { segment, support } => {
                let located = Located::new(support, pass_index, pass_name);
                validate_np_length_support(segment, &located, errors);
                self.np_length_supports.insert(segment, located);
            }
        }
    }
}

fn validate_allele_support(
    segment: Segment,
    support: &Located<AlleleIdSupport>,
    refdata: Option<&RefDataConfig>,
    errors: &mut Vec<CompileError>,
) {
    match &support.value {
        AlleleIdSupport::Invalid(reason) => {
            errors.push(support.error(invalid_parameter(
                sample_allele_address(segment),
                reason.clone(),
            )));
        }
        AlleleIdSupport::Enumerated(ids) => {
            if ids.is_empty() {
                errors.push(support.error(invalid_parameter(
                    sample_allele_address(segment),
                    "empty_support",
                )));
                return;
            }
            let Some(refdata) = refdata else {
                return;
            };
            let Some(pool) = refdata.pool_for(segment) else {
                return;
            };
            for id in ids {
                if id.as_usize() >= pool.len() {
                    errors.push(support.error(invalid_parameter(
                        sample_allele_address(segment),
                        format!("allele_id_{}_out_of_range", id.index()),
                    )));
                }
            }
        }
        AlleleIdSupport::Unavailable => {}
    }
}

fn validate_trim_support(
    segment: Segment,
    end: TrimEnd,
    support: &Located<IntegerSupport>,
    errors: &mut Vec<CompileError>,
) {
    let address = trim_address(segment, end);
    match &support.value {
        IntegerSupport::Invalid(reason) => {
            errors.push(support.error(invalid_parameter(address, reason.clone())));
        }
        IntegerSupport::Enumerated(values) => {
            if values.is_empty() {
                errors.push(support.error(invalid_parameter(address, "empty_support")));
                return;
            }
            for value in values {
                if *value < 0 {
                    errors.push(support.error(invalid_parameter(address, "negative_trim")));
                } else if *value > u16::MAX as i64 {
                    errors.push(support.error(invalid_parameter(address, "trim_exceeds_u16")));
                }
            }
        }
        IntegerSupport::Unavailable => {}
    }
}

fn validate_np_length_support(
    segment: Segment,
    support: &Located<IntegerSupport>,
    errors: &mut Vec<CompileError>,
) {
    let address = np_length_address(segment);
    match &support.value {
        IntegerSupport::Invalid(reason) => {
            errors.push(support.error(invalid_parameter(address, reason.clone())));
        }
        IntegerSupport::Enumerated(values) => {
            if values.is_empty() {
                errors.push(support.error(invalid_parameter(address, "empty_support")));
                return;
            }
            for value in values {
                if *value < 0 {
                    errors.push(support.error(invalid_parameter(address, "negative_length")));
                } else if *value > u32::MAX as i64 {
                    errors.push(support.error(invalid_parameter(address, "length_exceeds_u32")));
                }
            }
        }
        IntegerSupport::Unavailable => {}
    }
}

fn validate_contract_preconditions(
    contracts: Option<&ContractSet>,
    refdata: Option<&RefDataConfig>,
    facts: &CompileFactIndex,
    errors: &mut Vec<CompileError>,
) {
    let Some(contracts) = contracts else {
        return;
    };

    if contracts.contains_kind(ContractKind::ProductiveJunctionFrame) {
        validate_productive_frame_preconditions(refdata, facts, errors);
    }
}

fn validate_productive_frame_preconditions(
    refdata: Option<&RefDataConfig>,
    facts: &CompileFactIndex,
    errors: &mut Vec<CompileError>,
) {
    let contract_name = "productive_junction_frame";
    let Some(refdata) = refdata else {
        errors.push(plan_scoped_error(super::error::contract_precondition(
            contract_name,
            "productive frame validation requires reference data",
        )));
        return;
    };

    let Some(v_residues) = anchor_tail_residues(Segment::V, refdata, facts, errors) else {
        return;
    };
    let Some(j_residues) = anchor_head_residues(Segment::J, refdata, facts, errors) else {
        return;
    };

    let uses_d = facts.assigned_segments.contains(&Segment::D);
    if uses_d {
        let Some(np1_residues) = optional_np_length_residues(Segment::Np1, facts, errors) else {
            return;
        };
        let Some(d_residues) = d_length_residues(refdata, facts, errors) else {
            return;
        };
        let Some(np2_residues) = required_np_length_residues(Segment::Np2, facts, errors) else {
            return;
        };

        if residues_admit_frame(&[
            &v_residues,
            &np1_residues,
            &d_residues,
            &np2_residues,
            &j_residues,
        ]) {
            return;
        }

        errors.push(plan_scoped_error(super::error::contract_precondition(
            contract_name,
            "NP2 length support has no in-frame mass for the active V/D/J anchor and trim supports",
        )));
    } else {
        let Some(np1_residues) = required_np_length_residues(Segment::Np1, facts, errors) else {
            return;
        };

        if residues_admit_frame(&[&v_residues, &np1_residues, &j_residues]) {
            return;
        }

        errors.push(plan_scoped_error(super::error::contract_precondition(
            contract_name,
            "NP1 length support has no in-frame mass for the active V/J anchor and trim supports",
        )));
    }
}

pub(super) fn analyze_plan(
    plan: &PassPlan,
    refdata: Option<&RefDataConfig>,
    contracts: Option<&ContractSet>,
    policy: ExecutionPolicy,
) -> (CompileReport, Vec<CompileError>, Option<FeasibilityContext>) {
    let mut assigned: HashSet<Segment> = HashSet::new();
    let mut pass_summaries = Vec::with_capacity(plan.len());
    let mut declared_choices = Vec::new();
    let mut errors = Vec::new();
    let mut fact_index = CompileFactIndex::default();
    let mut assembled: HashSet<Segment> = HashSet::new();

    for (index, pass) in plan.passes().iter().enumerate() {
        let name = pass.name().to_string();
        let choices = pass.declared_choices();
        let requirements = pass.requirements();
        let effects = pass.effects();
        let compile_facts = pass.compile_facts();

        for req in &requirements {
            match req {
                PassRequirement::RefData if refdata.is_none() => {
                    errors.push(pass_scoped_error(
                        index,
                        &name,
                        CompileErrorKind::MissingRefData,
                    ));
                }
                PassRequirement::AlleleAssignment(segment) if !assigned.contains(segment) => {
                    errors.push(pass_scoped_error(
                        index,
                        &name,
                        CompileErrorKind::MissingAssignment { segment: *segment },
                    ));
                }
                _ => {}
            }
        }

        for fact in compile_facts {
            fact_index.record_fact(fact, index, &name, refdata, &mut errors);
        }

        for address in &choices {
            declared_choices.push(DeclaredChoice {
                pass_index: index,
                pass_name: name.clone(),
                address: address.clone(),
            });
        }

        for effect in &effects {
            match effect {
                PassEffect::AssignAllele(segment) => {
                    assigned.insert(*segment);
                    fact_index.assigned_segments.insert(*segment);
                }
                PassEffect::TrimAllele(segment) if assembled.contains(segment) => {
                    errors.push(pass_scoped_error(
                        index,
                        &name,
                        invalid_pass_order(format!("trim_after_assembly.{:?}", segment)),
                    ));
                }
                PassEffect::AssembleSegment(segment) => {
                    assembled.insert(*segment);
                }
                _ => {}
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

    validate_contract_preconditions(contracts, refdata, &fact_index, &mut errors);
    let feasibility = if errors.is_empty() {
        build_feasibility_context(contracts, refdata, &fact_index)
    } else {
        None
    };

    (
        CompileReport {
            policy,
            pass_summaries,
            declared_choices,
            active_contracts,
            warnings: Vec::new(),
        },
        errors,
        feasibility,
    )
}
