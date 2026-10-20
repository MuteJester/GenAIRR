//! Compile-time plan analysis: validates declared supports, builds the
//! `CompileReport`, and derives the optional `FeasibilityContext`.
//!
//! The `CompileFactIndex` is the shared data structure that holds
//! supports gathered from each pass. `analyze_plan` is the single
//! entry point used by both `CompiledSimulator::compile` and
//! `OwnedCompiledSimulator::compile`.

use super::error::{invalid_parameter, plan_scoped_error, CompileError, CompileErrorKind};
use super::feasibility_builder::{
    anchor_head_residues, anchor_tail_residues, build_feasibility_context, d_length_residues,
    optional_np_length_residues, required_np_length_residues, residues_admit_frame,
    validate_coding_supports,
};
use super::report::{CompileReport, DeclaredChoice, PassSummary};
use super::{np_length_address, sample_allele_address, trim_address, ExecutionPolicy};
use crate::address::PEnd;
use crate::assignment::TrimEnd;
use crate::contract::{ContractKind, ContractSet};
use crate::feasibility::FeasibilityContext;
use crate::ir::Segment;
use crate::pass::{
    AlleleIdSupport, IntegerSupport, NpBaseSupport, PassCompileFact, PassEffect, PassPlan,
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
    pub(super) fn new(value: T, pass_index: usize, pass_name: &str) -> Self {
        Self {
            value,
            pass_index,
            pass_name: pass_name.to_string(),
        }
    }

    pub(super) fn error(&self, kind: CompileErrorKind) -> CompileError {
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
    pub(super) p_length_supports: Vec<(PEnd, Located<IntegerSupport>)>,
    pub(super) coding_end_supports:
        HashMap<Segment, Located<Vec<crate::end_processing::EndOutcome>>>,
    pub(super) coding_end_emissions: Vec<(PEnd, Located<()>)>,
    pub(super) np_base_supports: HashMap<Segment, Located<Option<NpBaseSupport>>>,
    pub(super) d_orientation_supports: Vec<Located<Vec<bool>>>,
    pub(super) append_passes: HashMap<Segment, Vec<usize>>,
    pub(super) assigned_segments: HashSet<Segment>,
    pub(super) editing_passes: Vec<usize>,
}

impl CompileFactIndex {
    /// True iff the plan contains the recombination passes that
    /// produce the V/J anchor + NP-length facts the
    /// `productive_junction_frame` precondition check requires.
    /// Post-fork clonal plans inherit a pre-recombined IR and won't
    /// have any of these facts; we use this to skip preconditions
    /// that can't be evaluated rather than fail compile.
    pub(super) fn has_recombination_support(&self) -> bool {
        !self.allele_supports.is_empty() || !self.np_length_supports.is_empty()
    }

    pub(super) fn has_nonzero_p_support(&self) -> bool {
        !self.coding_end_supports.is_empty()
            || self.p_length_supports.iter().any(|(_, support)| {
                !matches!(&support.value, IntegerSupport::Enumerated(values)
                if values.iter().all(|value| *value == 0))
            })
    }

    /// Full productive generation needs the same complete finite continuation
    /// solver with or without P additions. Frame residues alone cannot reject
    /// future D/junction stops. Partial inherited plans keep their runtime
    /// contracts and have no new recombination choices to solve here.
    pub(super) fn needs_complete_productive_support(&self, contracts: &ContractSet) -> bool {
        let complete = contracts.requires_productive_emission()
            && [Segment::V, Segment::J].into_iter().all(|segment| {
                self.allele_supports.contains_key(&segment)
                    && self
                        .append_passes
                        .get(&segment)
                        .is_some_and(|p| p.len() == 1)
            });
        complete
            || (self.has_nonzero_p_support()
                && (contracts.contains_kind(ContractKind::ProductiveJunctionFrame)
                    || contracts.contains_kind(ContractKind::NoStopCodonInJunction)
                    || contracts.contains_kind(ContractKind::CodingSequenceIntegrity)))
    }
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
            PassCompileFact::PLengthSupport { end, support } => {
                // Retain every pass, including repeated additions at one
                // end, so no nonzero support can hide behind a later zero.
                self.p_length_supports
                    .push((end, Located::new(support, pass_index, pass_name)));
            }
            PassCompileFact::NpBaseSupport { segment, support } => {
                self.np_base_supports
                    .insert(segment, Located::new(support, pass_index, pass_name));
            }
            PassCompileFact::DOrientationSupport { orientations } => {
                self.d_orientation_supports
                    .push(Located::new(orientations, pass_index, pass_name));
            }
            PassCompileFact::CodingEndSupport { segment, outcomes } => {
                let entry = Located::new(outcomes, pass_index, pass_name);
                if self.coding_end_supports.contains_key(&segment) {
                    errors.push(entry.error(invalid_parameter(
                        pass_name,
                        "coding-end processing must be declared once per segment",
                    )));
                }
                self.coding_end_supports.insert(segment, entry);
            }
            PassCompileFact::CodingEndEmission { end } => {
                self.coding_end_emissions
                    .push((end, Located::new((), pass_index, pass_name)));
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

/// Reject provably impossible trim geometry independently of productivity.
/// Mixed supports can still select an incompatible allele at runtime; assembly
/// checks those realized values through the fallible execution path.
fn validate_trim_geometry(
    refdata: Option<&RefDataConfig>,
    facts: &CompileFactIndex,
    errors: &mut Vec<CompileError>,
) {
    let Some(refdata) = refdata else {
        return;
    };
    for segment in [Segment::V, Segment::D, Segment::J] {
        if !facts.append_passes.contains_key(&segment) {
            continue;
        }
        let minimum = |end| match facts.trim_supports.get(&(segment, end)) {
            None => Some(0),
            Some(entry) => match &entry.value {
                IntegerSupport::Enumerated(values) => values
                    .iter()
                    .copied()
                    .filter(|v| (0..=u16::MAX as i64).contains(v))
                    .min(),
                _ => None,
            },
        };
        let (Some(five), Some(three)) = (minimum(TrimEnd::Five), minimum(TrimEnd::Three)) else {
            continue;
        };
        let maximum_length =
            super::feasibility_builder::allele_ids_for_segment(segment, refdata, facts)
                .into_iter()
                .filter_map(|id| refdata.get(segment, id).map(|a| i64::from(a.len())))
                .max();
        let Some(maximum_length) = maximum_length else {
            continue;
        };
        if five + three <= maximum_length {
            continue;
        }
        let Some(((segment, end), support)) = facts
            .trim_supports
            .iter()
            .filter(|((s, _), _)| *s == segment)
            .max_by_key(|(_, support)| support.pass_index)
        else {
            continue;
        };
        errors.push(support.error(invalid_parameter(
            trim_address(*segment, *end),
            format!(
                "trim_exceeds_allele_length: minimum combined trim {} exceeds every sampled allele (maximum length {})",
                five + three, maximum_length,
            ),
        )));
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

    if contracts.contains_kind(ContractKind::CodingSequenceIntegrity) {
        if facts.has_recombination_support() {
            validate_coding_supports(refdata, facts, errors);
        }
    }

    if facts.needs_complete_productive_support(contracts) {
        super::feasibility_builder::validate_p_productive_supports(
            contracts, refdata, facts, errors,
        );
        return;
    }

    if contracts.contains_kind(ContractKind::ProductiveJunctionFrame) {
        // Only validate productive-frame preconditions when the plan
        // actually contains the recombination passes that produce the
        // anchor / NP-length facts the check consumes. Post-fork plans
        // (the per-descendant half of a clonal fork) inherit the
        // assembled IR from the parent and only run mutation /
        // corruption passes, so those facts will be absent. The
        // contract bundle itself stays active at runtime — every
        // substitution and indel still gets contract-filtered against
        // the productive shape — but the build-time precondition
        // can't be evaluated without recombination support.
        if facts.has_recombination_support() {
            validate_productive_frame_preconditions(refdata, facts, errors);
        }
    }
}

/// Keep the opt-in outcome model unambiguous even for direct native plans.
fn validate_coding_end_geometry(
    refdata: Option<&RefDataConfig>,
    facts: &CompileFactIndex,
    errors: &mut Vec<CompileError>,
) {
    use crate::assignment::SegmentOrientation;
    if !facts.coding_end_supports.is_empty() {
        if !facts.p_length_supports.is_empty() {
            errors.push(plan_scoped_error(invalid_parameter(
                "end_processing.v1",
                "original-end processing cannot be mixed with legacy P-addition passes",
            )));
        }
        let mut ordered = Vec::new();
        for (&segment, passes) in &facts.append_passes {
            let rank = match segment {
                Segment::V => 0,
                Segment::Np1 => 2,
                Segment::D => 5,
                Segment::Np2 => 7,
                Segment::J => 9,
            };
            if passes.len() != 1 {
                errors.push(plan_scoped_error(invalid_parameter("end_processing.v1",
                    "original-end processing requires at most one assembly per segment and one N generation per gap")));
            }
            ordered.extend(passes.iter().map(|&index| (index, rank)));
        }
        for (end, entry) in &facts.coding_end_emissions {
            ordered.push((
                entry.pass_index,
                match end {
                    PEnd::V3 => 1,
                    PEnd::D5 => 4,
                    PEnd::D3 => 6,
                    PEnd::J5 => 8,
                },
            ));
        }
        for entry in &facts.d_orientation_supports {
            ordered.push((entry.pass_index, 3));
        }
        ordered.sort_unstable();
        if facts.d_orientation_supports.len() > 1
            || ordered.windows(2).any(|pair| pair[0].1 > pair[1].1)
        {
            errors.push(plan_scoped_error(invalid_parameter("end_processing.v1",
                "original-end processing requires coding-order assembly and P emission, with D orientation before D5-P and D assembly")));
        }
        if let Some(&(completion, _)) = ordered.last() {
            if facts.editing_passes.iter().any(|&index| index < completion) {
                errors.push(plan_scoped_error(invalid_parameter("end_processing.v1",
                    "mutation and observation edits must follow completion of original-end recombination")));
            }
        }
    }
    for (&segment, entry) in &facts.coding_end_supports {
        let invalid = |reason| entry.error(invalid_parameter(&entry.pass_name, reason));
        if facts.trim_supports.keys().any(|(s, _)| *s == segment)
            || facts
                .p_length_supports
                .iter()
                .any(|(end, _)| end.source_segment() == segment)
        {
            errors.push(invalid(
                "coding-end processing replaces legacy trim and P choices for its segment",
            ));
        }
        let expected: &[PEnd] = match segment {
            Segment::V => &[PEnd::V3],
            Segment::D => &[PEnd::D5, PEnd::D3],
            Segment::J => &[PEnd::J5],
            _ => {
                errors.push(invalid(
                    "coding-end processing requires a V, D or J segment",
                ));
                continue;
            }
        };
        let assembly = facts
            .append_passes
            .get(&segment)
            .filter(|passes| passes.len() == 1)
            .map(|passes| passes[0]);
        if assembly.is_none_or(|index| index <= entry.pass_index) {
            errors.push(invalid(
                "coding-end processing requires one later segment assembly",
            ));
        }
        for &end in expected {
            let emissions: Vec<_> = facts
                .coding_end_emissions
                .iter()
                .filter(|(e, _)| *e == end)
                .collect();
            if emissions.len() != 1 || emissions[0].1.pass_index <= entry.pass_index {
                errors.push(invalid("coding-end processing requires exactly one later P emission for each coding end, including zero-length outcomes"));
            } else if let Some(assembly) = assembly {
                let emission = emissions[0].1.pass_index;
                let is_three = matches!(end, PEnd::V3 | PEnd::D3);
                if (is_three && emission <= assembly) || (!is_three && emission >= assembly) {
                    errors.push(invalid(
                        "coding-end P emission must bracket segment assembly in coding order",
                    ));
                }
            }
        }
        if entry.value.is_empty()
            || entry.value.iter().any(|row| {
                !row.weight.is_finite()
                    || row.weight <= 0.0
                    || (segment == Segment::V && row.five != 0)
                    || (segment == Segment::J && row.three != 0)
            })
        {
            errors.push(invalid(
                "coding-end outcomes require positive finite weights and valid segment ends",
            ));
            continue;
        }
        let Some(refdata) = refdata else {
            errors.push(invalid("coding-end processing requires reference data"));
            continue;
        };
        let possible = super::feasibility_builder::allele_ids_for_segment(segment, refdata, facts)
            .into_iter()
            .any(|id| {
                refdata.coding_ends_are_not_known_partial(segment, id)
                    && refdata.get(segment, id).is_some_and(|allele| {
                        entry.value.iter().any(|row| {
                            crate::end_processing::materialize(
                                &allele.seq,
                                row.five,
                                row.three,
                                SegmentOrientation::Forward,
                            )
                            .is_ok()
                        })
                    })
            });
        if !possible {
            errors.push(invalid(
                "coding-end outcomes have no legal geometry for any sampled allele (known partial coding boundaries are excluded)",
            ));
        }
    }
    for (end, emission) in &facts.coding_end_emissions {
        if !facts
            .coding_end_supports
            .contains_key(&end.source_segment())
        {
            errors.push(emission.error(invalid_parameter(
                &emission.pass_name,
                "coding-end P emission requires an original-end processing choice",
            )));
        }
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
    let mut pass_summaries = Vec::with_capacity(plan.len());
    let mut declared_choices = Vec::new();
    let mut errors = Vec::new();
    let mut fact_index = CompileFactIndex::default();

    for (index, pass) in plan.passes().iter().enumerate() {
        let name = pass.name().to_string();
        let choices: Vec<String> = pass
            .declared_choice_patterns()
            .into_iter()
            .map(String::from)
            .collect();
        let requirements = pass.requirements();
        let effects = pass.effects();
        let compile_facts = pass.compile_facts();
        // Revision rewrites a segment through runtime events and deliberately
        // has no static pool effect; it still follows the original join.
        if name == crate::address::RECEPTOR_REVISION {
            fact_index.editing_passes.push(index);
        }

        // Pre-`Schedule` we checked `RefData` and `AlleleAssignment`
        // requirements here, and rejected `TrimAllele(seg)` if it
        // came after `AssembleSegment(seg)` in declared order. The
        // `Schedule::compile` topo-sort now owns all three: missing
        // requirements come back as `ScheduleError`, and a misordered
        // trim is silently reordered before assemble via the
        // implicit `TrimAllele → AssembleSegment` edge. `analyze_plan`
        // keeps only the facts-and-feasibility responsibilities.

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
            if matches!(effect, PassEffect::EditBases | PassEffect::StructuralIndel) {
                fact_index.editing_passes.push(index);
            }
            if let PassEffect::AssignAllele(segment) = effect {
                fact_index.assigned_segments.insert(*segment);
            }
            if let PassEffect::AppendRegion(segment) | PassEffect::AssembleSegment(segment) = effect
            {
                fact_index
                    .append_passes
                    .entry(*segment)
                    .or_default()
                    .push(index);
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

    validate_trim_geometry(refdata, &fact_index, &mut errors);
    validate_coding_end_geometry(refdata, &fact_index, &mut errors);
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
