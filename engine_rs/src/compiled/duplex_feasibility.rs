//! Compile a finite authored repair vocabulary before any biological masking.
use super::analyze::CompileFactIndex;
use super::error::{invalid_parameter, plan_scoped_error, CompileError, CompileErrorKind};
use crate::contract::{ContractKind, ContractSet};
use crate::duplex_repair::model::DuplexContextKey;
use crate::feasibility::{DuplexCompletionFeasibility, FeasibilityContext};
use crate::ir::Segment;
use crate::pass::{AlleleIdSupport, PassPlan};
use crate::passes::sample_allele::validate_finite_weights;
use crate::refdata::RefDataConfig;
use std::any::TypeId;
use std::collections::{BTreeMap, BTreeSet};

pub(super) fn build(
    plan: &PassPlan,
    refdata: Option<&RefDataConfig>,
    contracts: Option<&ContractSet>,
    facts: &CompileFactIndex,
    errors: &mut Vec<CompileError>,
) -> Option<FeasibilityContext> {
    if facts.finite_duplex_sources.len() != 1 {
        errors.push(plan_scoped_error(invalid_parameter(
            "duplex.recombine",
            "exactly one atomic duplex source is supported",
        )));
        return None;
    }
    let entry = &facts.finite_duplex_sources[0];
    let invalid = |reason: String| entry.error(invalid_parameter("duplex.recombine", reason));
    let Some(refdata) = refdata else {
        errors.push(entry.error(CompileErrorKind::MissingRefData));
        return None;
    };
    let model = &entry.value;
    let resource = |name: &str, required, limit| {
        entry.error(CompileErrorKind::ResourceLimit {
            resource: name.into(),
            required,
            limit,
        })
    };
    if let Err(reason) = model.validate_schema() {
        errors.push(invalid(reason));
        return None;
    }
    if model.reference_content_hash != crate::trace_file::refdata_content_hash(refdata) {
        errors.push(invalid(
            "finite duplex model reference content hash mismatch".into(),
        ));
        return None;
    }
    if contracts.is_some_and(|c| {
        c.iter().any(|c| {
            matches!(
                c.kind(),
                ContractKind::Custom | ContractKind::AnchorPreserved { .. }
            )
        })
    }) {
        errors.push(entry.error(CompileErrorKind::UnsupportedContractCombo {
            contract_name: "finite_duplex".into(),
            reason: "finite duplex generation supports canonical built-in contracts; custom and preservation-only contracts are unsupported".into(),
        }));
    }
    let segments = [Segment::V, Segment::D, Segment::J];
    let mut source_passes = [usize::MAX; 3];
    let mut marginals = [Vec::new(), Vec::new(), Vec::new()];
    for (index, segment) in segments.into_iter().enumerate() {
        crate::family_replay::primitive_work(1);
        let Some(passes) = facts
            .allele_sample_passes
            .get(&segment)
            .filter(|p| p.len() == 1)
        else {
            errors.push(invalid(format!(
                "finite duplex generation requires exactly one {segment:?} sampling source"
            )));
            continue;
        };
        source_passes[index] = passes[0];
        if passes[0] >= entry.pass_index {
            errors.push(invalid(
                "all allele choices must precede atomic duplex repair".into(),
            ));
        }
        match facts.allele_supports.get(&segment).map(|s| &s.value) {
            Some(AlleleIdSupport::Enumerated(ids)) if !ids.is_empty() => {
                marginals[index] = ids.clone();
                // A custom distribution may expose several weighted atoms
                // with the same outcome. Preserve those atoms in its CDF, but
                // enumerate each biological source context only once.
                marginals[index].sort_by_key(|id| id.index());
                marginals[index].dedup();
            }
            _ => errors.push(invalid(format!(
                "finite duplex generation requires finite positive {segment:?} allele support"
            ))),
        }
    }
    if !facts.trim_supports.is_empty()
        || !facts.np_length_supports.is_empty()
        || !facts.np_base_supports.is_empty()
        || !facts.p_length_supports.is_empty()
        || !facts.coding_end_supports.is_empty()
        || !facts.coding_end_emissions.is_empty()
        || facts
            .append_passes
            .values()
            .flatten()
            .any(|&i| i != entry.pass_index)
    {
        errors.push(invalid(
            "atomic duplex repair replaces legacy trimming, N/P generation and segment assembly"
                .into(),
        ));
    }
    let mut orientation_pass = None;
    let orientations = match facts.d_orientation_supports.as_slice() {
        [] => vec![false],
        [support]
            if support.pass_index < entry.pass_index
                && support.pass_index > source_passes[1]
                && !support.value.is_empty() =>
        {
            orientation_pass = Some(support.pass_index);
            support.value.clone()
        }
        _ => {
            errors.push(invalid("finite duplex generation permits one D orientation choice after D selection and before repair".into()));
            Vec::new()
        }
    };
    let mut flat_weights = Vec::new();
    let mut numeric_checks = 0u64;
    for (index, pass) in plan.passes().iter().enumerate() {
        crate::family_replay::primitive_work(1);
        let id = pass.as_ref().type_id();
        if let Some(source) =
            (pass.as_ref() as &dyn std::any::Any).downcast_ref::<crate::passes::SampleAllelePass>()
        {
            if let Some(weights) = source.finite_law_support() {
                numeric_checks = numeric_checks.saturating_add(weights.len() as u64);
                if numeric_checks > model.max_preflight_operations {
                    errors.push(resource(
                        "finite_duplex_preflight_operations",
                        numeric_checks,
                        model.max_preflight_operations,
                    ));
                    return None;
                }
                if let Err(reason) = validate_finite_weights(
                    &weights
                        .iter()
                        .map(|(_, weight)| *weight)
                        .collect::<Vec<_>>(),
                ) {
                    errors.push(invalid(reason));
                }
                flat_weights.push((index, source.segment(), weights));
            }
        }
        if let Some(genotype) = (pass.as_ref() as &dyn std::any::Any)
            .downcast_ref::<crate::passes::sample_genotype::SampleGenotypePass>(
        ) {
            match genotype.validate_finite_law() {
                Err(reason) => errors.push(invalid(reason)),
                Ok(checks) => {
                    numeric_checks = numeric_checks.saturating_add(checks);
                    if numeric_checks > model.max_preflight_operations {
                        errors.push(resource(
                            "finite_duplex_preflight_operations",
                            numeric_checks,
                            model.max_preflight_operations,
                        ));
                        return None;
                    }
                }
            }
        }
        let supported = if index == entry.pass_index {
            id == TypeId::of::<crate::passes::DuplexRecombinePass>()
        } else if index < entry.pass_index {
            (id == TypeId::of::<crate::passes::SampleAllelePass>()
                && source_passes.contains(&index))
                || (id == TypeId::of::<crate::passes::sample_genotype::SampleGenotypePass>()
                    && source_passes.iter().all(|&i| i == index))
                || (id == TypeId::of::<crate::passes::InvertDPass>()
                    && orientation_pass == Some(index))
        } else {
            super::execute::duplex_pass_supported(pass.as_ref())
        };
        if !supported {
            errors.push(invalid(format!(
                "pass {} is unsupported in a finite duplex generation plan",
                pass.name()
            )));
        }
    }
    if plan.compile(true).is_ok_and(|order| {
        order
            .iter()
            .take(entry.pass_index + 1)
            .enumerate()
            .any(|(i, node)| node.index() != i)
    }) {
        errors.push(invalid(
            "finite duplex source choices must execute in declared order".into(),
        ));
    }
    let groups = match facts.joint_allele_supports.as_slice() {
        [] => vec![marginals],
        [joint] if source_passes.iter().all(|&index| index == joint.pass_index) => {
            joint.value.clone()
        }
        _ => {
            errors.push(invalid("finite duplex source linkage is ambiguous".into()));
            Vec::new()
        }
    };
    if !errors.is_empty() {
        return None;
    }
    let mut contexts = BTreeSet::new();
    let mut operations = numeric_checks;
    for group in &groups {
        for &v in &group[0] {
            for &d in &group[1] {
                for &j in &group[2] {
                    for &d_inverted in &orientations {
                        crate::family_replay::primitive_work(1);
                        if contexts.insert(DuplexContextKey {
                            v,
                            d,
                            j,
                            d_inverted,
                        }) {
                            operations = operations.saturating_add(1);
                            if operations > model.max_preflight_operations {
                                errors.push(resource(
                                    "finite_duplex_preflight_operations",
                                    operations,
                                    model.max_preflight_operations,
                                ));
                                return None;
                            }
                        }
                    }
                }
            }
        }
    }
    if contexts.is_empty() {
        errors.push(invalid(
            "finite duplex source has no positive natural joint context".into(),
        ));
        return None;
    }
    let mut covered = BTreeSet::new();
    let mut context_rows: BTreeMap<DuplexContextKey, Vec<usize>> = BTreeMap::new();
    let mut active_rows = Vec::new();
    let mut reserved_bases = 0u64;
    for (index, row) in model
        .rows
        .iter()
        .enumerate()
        .filter(|(_, row)| row.weight > 0.0)
    {
        let key = match model.context_key(index, refdata) {
            Ok(key) => key,
            Err(reason) => {
                errors.push(invalid(reason));
                continue;
            }
        };
        if !contexts.contains(&key) {
            continue;
        }
        covered.insert(key);
        context_rows.entry(key).or_default().push(index);
        active_rows.push(index);
        operations = row.script.joins.iter().fold(operations, |total, join| {
            total.saturating_add(join.operations.len() as u64)
        });
        reserved_bases = reserved_bases.saturating_add(row.script.max_bases as u64);
    }
    if operations > model.max_preflight_operations {
        errors.push(resource(
            "finite_duplex_preflight_operations",
            operations,
            model.max_preflight_operations,
        ));
    }
    if reserved_bases > model.max_preflight_bases {
        errors.push(resource(
            "finite_duplex_preflight_base_reservation",
            reserved_bases,
            model.max_preflight_bases,
        ));
    }
    if let Some(missing) = contexts.difference(&covered).next() {
        errors.push(invalid(format!(
            "finite duplex model lacks positive history support for natural context V={},D={},J={},d_inverted={}",
            missing.v.index(), missing.d.index(), missing.j.index(), missing.d_inverted,
        )));
    }
    if !errors.is_empty() {
        return None;
    }
    for indices in context_rows.values() {
        crate::family_replay::primitive_work(indices.len() as u64);
        if let Err(reason) = model.validate_weight_support(indices) {
            errors.push(invalid(reason));
        }
    }
    if !errors.is_empty() {
        return None;
    }
    // Complete physical/capability preflight precedes every productivity mask.
    let mut prepared = Vec::new();
    let model_digest = model.digest();
    for index in active_rows {
        crate::family_replay::primitive_work(1);
        match model.prepare_row_with_digest(index, refdata, &model_digest) {
            Ok(row) => prepared.push(row),
            Err(reason) => errors.push(invalid(reason)),
        }
    }
    if !errors.is_empty() {
        return None;
    }
    let mut proven_false = 0usize;
    let mut unresolved = 0usize;
    let mut contract_rejected = 0usize;
    prepared.retain(|row| {
        contracts.is_none_or(|contracts| {
            let contract_ok = contracts.verify(&row.seed, Some(refdata)).is_ok();
            if contracts.requires_productive_emission() {
                match crate::productivity::assess(&row.seed, refdata).productive {
                    Some(false) => {
                        proven_false += 1;
                        return false;
                    }
                    None => {
                        unresolved += 1;
                        return false;
                    }
                    Some(true) => {}
                }
            }
            if !contract_ok {
                contract_rejected += 1;
            }
            contract_ok
        })
    });
    let mut accepted_context_rows: BTreeMap<DuplexContextKey, Vec<usize>> = BTreeMap::new();
    for row in &prepared {
        crate::family_replay::primitive_work(1);
        accepted_context_rows
            .entry(row.context)
            .or_default()
            .push(row.row_index);
    }
    for indices in accepted_context_rows.values() {
        crate::family_replay::primitive_work(indices.len() as u64);
        if let Err(reason) = model.validate_weight_support(indices) {
            errors.push(invalid(reason));
        }
    }
    if !errors.is_empty() {
        return None;
    }
    // Each accepted row witnesses a reachable prefix. Histories with the same
    // prefix induce precisely the IDs that the local allele sampler retains.
    // Full-support CDF advancement alone does not prove every masked support:
    // removing a preceding half-ULP atom can change round-to-even tie parity.
    for (owner, segment, weights) in flat_weights {
        crate::family_replay::primitive_work(weights.len() as u64 + 1);
        let slot = segments
            .iter()
            .position(|s| *s == segment)
            .expect("V/D/J source");
        let mut by_prefix: BTreeMap<Vec<u32>, BTreeSet<u32>> = BTreeMap::new();
        for row in &prepared {
            let ids = [
                row.context.v.index(),
                row.context.d.index(),
                row.context.j.index(),
            ];
            let mut prefix: Vec<_> = source_passes
                .iter()
                .enumerate()
                .filter(|(_, source)| **source < owner)
                .map(|(index, _)| ids[index])
                .collect();
            if orientation_pass.is_some_and(|source| source < owner) {
                prefix.push(u32::from(row.context.d_inverted));
            }
            by_prefix.entry(prefix).or_default().insert(ids[slot]);
        }
        let masks: BTreeSet<Vec<u32>> = by_prefix
            .into_values()
            .map(|ids| ids.into_iter().collect())
            .collect();
        for allowed in masks {
            crate::family_replay::primitive_work(allowed.len() as u64 + 1);
            // Charge every advertised entry examined, including excluded
            // entries, once for each distinct candidate mask at this source.
            operations = operations.saturating_add(weights.len() as u64);
            if operations > model.max_preflight_operations {
                errors.push(resource(
                    "finite_duplex_preflight_operations",
                    operations,
                    model.max_preflight_operations,
                ));
                return None;
            }
            let selected: Vec<_> = weights
                .iter()
                .filter(|(id, _)| allowed.binary_search(&id.index()).is_ok())
                .map(|(_, weight)| *weight)
                .collect();
            if let Err(reason) = validate_finite_weights(&selected) {
                errors.push(invalid(format!(
                    "{reason} in reachable {segment:?} allele support after completion masking"
                )));
            }
        }
    }
    if !errors.is_empty() {
        return None;
    }
    if prepared.is_empty() {
        errors.push(entry.error(CompileErrorKind::ContractPrecondition {
            contract_name: "duplex_observed_vdj_v1".into(),
            reason: format!(
                "finite authored histories have no positive admissible completion under the active contracts \
                 (proven_false={proven_false}, unresolved={unresolved}, contract_rejected={contract_rejected})"
            ),
        }));
        return None;
    }
    Some(
        FeasibilityContext::empty().with_duplex(DuplexCompletionFeasibility {
            source_passes,
            orientation_pass,
            repair_pass: entry.pass_index,
            rows: prepared,
        }),
    )
}
