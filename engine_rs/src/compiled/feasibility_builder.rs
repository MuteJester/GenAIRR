//! Build-time helpers that derive `FeasibilityContext` and productive-frame
//! residue sets from a plan's declared supports.
//!
//! These free functions operate over the `CompileFactIndex` collected
//! by `analyze.rs`. They are kept `pub(super)` so the analyze layer can
//! call them during precondition validation and feasibility context
//! construction, but they are not part of the compiled-simulator public
//! API.

use super::analyze::{CompileFactIndex, Located};
use super::error::{contract_precondition, plan_scoped_error, CompileError};
use super::{np_length_address, trim_address};
use crate::assignment::TrimEnd;
use crate::contract::{ContractKind, ContractSet};
use crate::feasibility::{ChoiceDomain, FeasibilityContext, VjProductiveFeasibility};
use crate::ir::Segment;
use crate::pass::{AlleleIdSupport, IntegerSupport};
use crate::refdata::{AlleleId, RefDataConfig};
use std::collections::HashSet;

pub(super) fn allele_ids_for_segment(
    segment: Segment,
    refdata: &RefDataConfig,
    facts: &CompileFactIndex,
) -> Vec<AlleleId> {
    let Some(pool) = refdata.pool_for(segment) else {
        return Vec::new();
    };

    match facts
        .allele_supports
        .get(&segment)
        .map(|entry| &entry.value)
    {
        Some(AlleleIdSupport::Enumerated(ids)) => ids
            .iter()
            .copied()
            .filter(|id| id.as_usize() < pool.len())
            .collect(),
        Some(AlleleIdSupport::Invalid(_)) => Vec::new(),
        Some(AlleleIdSupport::Unavailable) | None => pool.iter().map(|(id, _)| id).collect(),
    }
}

pub(super) fn trim_values(
    segment: Segment,
    end: TrimEnd,
    facts: &CompileFactIndex,
    errors: &mut Vec<CompileError>,
) -> Option<Vec<u32>> {
    let Some(entry) = facts.trim_supports.get(&(segment, end)) else {
        return Some(vec![0]);
    };

    match &entry.value {
        IntegerSupport::Enumerated(values) => {
            let values: Vec<u32> = values
                .iter()
                .copied()
                .filter(|value| (0..=u16::MAX as i64).contains(value))
                .map(|value| value as u32)
                .collect();
            if values.is_empty() {
                None
            } else {
                Some(values)
            }
        }
        IntegerSupport::Unavailable => {
            errors.push(plan_scoped_error(contract_precondition(
                "productive_junction_frame",
                format!(
                    "{} support is not enumerable, so productive frame preconditions cannot be validated",
                    trim_address(segment, end)
                ),
            )));
            None
        }
        IntegerSupport::Invalid(_) => None,
    }
}

pub(super) fn required_np_length_residues(
    segment: Segment,
    facts: &CompileFactIndex,
    errors: &mut Vec<CompileError>,
) -> Option<HashSet<u8>> {
    let Some(entry) = facts.np_length_supports.get(&segment) else {
        errors.push(plan_scoped_error(contract_precondition(
            "productive_junction_frame",
            format!(
                "{} support is required for productive frame validation",
                np_length_address(segment)
            ),
        )));
        return None;
    };
    np_length_residues_from_entry(segment, entry, true, errors)
}

pub(super) fn optional_np_length_residues(
    segment: Segment,
    facts: &CompileFactIndex,
    errors: &mut Vec<CompileError>,
) -> Option<HashSet<u8>> {
    let Some(entry) = facts.np_length_supports.get(&segment) else {
        return Some(HashSet::from([0]));
    };
    np_length_residues_from_entry(segment, entry, false, errors)
}

fn np_length_residues_from_entry(
    segment: Segment,
    entry: &Located<IntegerSupport>,
    required: bool,
    errors: &mut Vec<CompileError>,
) -> Option<HashSet<u8>> {
    match &entry.value {
        IntegerSupport::Enumerated(values) => {
            let residues: HashSet<u8> = values
                .iter()
                .copied()
                .filter(|value| (0..=u32::MAX as i64).contains(value))
                .map(|value| (value % 3) as u8)
                .collect();
            if residues.is_empty() {
                if required {
                    errors.push(plan_scoped_error(contract_precondition(
                        "productive_junction_frame",
                        format!(
                            "{} has no valid non-negative support",
                            np_length_address(segment)
                        ),
                    )));
                }
                None
            } else {
                Some(residues)
            }
        }
        IntegerSupport::Unavailable => {
            errors.push(plan_scoped_error(contract_precondition(
                "productive_junction_frame",
                format!(
                    "{} support is not enumerable, so productive frame preconditions cannot be validated",
                    np_length_address(segment)
                ),
            )));
            None
        }
        IntegerSupport::Invalid(_) => None,
    }
}

pub(super) fn anchor_tail_residues(
    segment: Segment,
    refdata: &RefDataConfig,
    facts: &CompileFactIndex,
    errors: &mut Vec<CompileError>,
) -> Option<HashSet<u8>> {
    let ids = allele_ids_for_segment(segment, refdata, facts);
    let anchored_count = ids
        .iter()
        .filter(|id| {
            refdata
                .get(segment, **id)
                .is_some_and(|allele| valid_anchor(allele.anchor, allele.len()))
        })
        .count();
    if anchored_count == 0 {
        errors.push(plan_scoped_error(contract_precondition(
            "productive_junction_frame",
            format!(
                "{:?} sample support has no alleles with valid anchors",
                segment
            ),
        )));
        return None;
    }

    let trim_5 = trim_values(segment, TrimEnd::Five, facts, errors)?;
    let trim_3 = trim_values(segment, TrimEnd::Three, facts, errors)?;
    let mut residues = HashSet::new();
    for id in ids {
        let Some(allele) = refdata.get(segment, id) else {
            continue;
        };
        let Some(anchor) = allele.anchor.map(|anchor| anchor as u32) else {
            continue;
        };
        if !valid_anchor(Some(anchor as u16), allele.len()) {
            continue;
        }
        for five in &trim_5 {
            for three in &trim_3 {
                let retained_end = allele.len().saturating_sub(*three);
                if *five <= anchor && anchor + 3 <= retained_end {
                    residues.insert(((retained_end - anchor) % 3) as u8);
                }
            }
        }
    }

    if residues.is_empty() {
        errors.push(plan_scoped_error(contract_precondition(
            "productive_junction_frame",
            format!("{:?} trim support removes every valid anchor", segment),
        )));
        None
    } else {
        Some(residues)
    }
}

pub(super) fn anchor_head_residues(
    segment: Segment,
    refdata: &RefDataConfig,
    facts: &CompileFactIndex,
    errors: &mut Vec<CompileError>,
) -> Option<HashSet<u8>> {
    let ids = allele_ids_for_segment(segment, refdata, facts);
    let anchored_count = ids
        .iter()
        .filter(|id| {
            refdata
                .get(segment, **id)
                .is_some_and(|allele| valid_anchor(allele.anchor, allele.len()))
        })
        .count();
    if anchored_count == 0 {
        errors.push(plan_scoped_error(contract_precondition(
            "productive_junction_frame",
            format!(
                "{:?} sample support has no alleles with valid anchors",
                segment
            ),
        )));
        return None;
    }

    let trim_5 = trim_values(segment, TrimEnd::Five, facts, errors)?;
    let trim_3 = trim_values(segment, TrimEnd::Three, facts, errors)?;
    let mut residues = HashSet::new();
    for id in ids {
        let Some(allele) = refdata.get(segment, id) else {
            continue;
        };
        let Some(anchor) = allele.anchor.map(|anchor| anchor as u32) else {
            continue;
        };
        if !valid_anchor(Some(anchor as u16), allele.len()) {
            continue;
        }
        for five in &trim_5 {
            for three in &trim_3 {
                let retained_end = allele.len().saturating_sub(*three);
                if *five <= anchor && anchor + 3 <= retained_end {
                    residues.insert(((anchor + 3 - *five) % 3) as u8);
                }
            }
        }
    }

    if residues.is_empty() {
        errors.push(plan_scoped_error(contract_precondition(
            "productive_junction_frame",
            format!("{:?} trim support removes every valid anchor", segment),
        )));
        None
    } else {
        Some(residues)
    }
}

pub(super) fn d_length_residues(
    refdata: &RefDataConfig,
    facts: &CompileFactIndex,
    errors: &mut Vec<CompileError>,
) -> Option<HashSet<u8>> {
    let ids = allele_ids_for_segment(Segment::D, refdata, facts);
    if ids.is_empty() {
        errors.push(plan_scoped_error(contract_precondition(
            "productive_junction_frame",
            "D sample support is empty for a D-containing productive plan",
        )));
        return None;
    }

    let trim_5 = trim_values(Segment::D, TrimEnd::Five, facts, errors)?;
    let trim_3 = trim_values(Segment::D, TrimEnd::Three, facts, errors)?;
    let mut residues = HashSet::new();
    for id in ids {
        let Some(allele) = refdata.get(Segment::D, id) else {
            continue;
        };
        for five in &trim_5 {
            for three in &trim_3 {
                if *five + *three <= allele.len() {
                    residues.insert(((allele.len() - *five - *three) % 3) as u8);
                }
            }
        }
    }

    if residues.is_empty() {
        errors.push(plan_scoped_error(contract_precondition(
            "productive_junction_frame",
            "D trim support leaves no valid retained D lengths",
        )));
        None
    } else {
        Some(residues)
    }
}

pub(super) fn valid_anchor(anchor: Option<u16>, allele_len: u32) -> bool {
    anchor.is_some_and(|anchor| (anchor as u32) + 3 <= allele_len)
}

pub(super) fn residues_admit_frame(sets: &[&HashSet<u8>]) -> bool {
    fn walk(sets: &[&HashSet<u8>], index: usize, residue: u8) -> bool {
        if index == sets.len() {
            return residue % 3 == 0;
        }
        sets[index]
            .iter()
            .any(|next| walk(sets, index + 1, (residue + *next) % 3))
    }
    walk(sets, 0, 0)
}

pub(super) fn build_feasibility_context(
    contracts: Option<&ContractSet>,
    refdata: Option<&RefDataConfig>,
    facts: &CompileFactIndex,
) -> Option<FeasibilityContext> {
    let contracts = contracts?;
    if !contracts.contains_kind(ContractKind::ProductiveJunctionFrame) {
        return None;
    }
    if facts.assigned_segments.contains(&Segment::D) {
        return None;
    }

    let refdata = refdata?;
    let vj_productive = build_vj_productive_feasibility(refdata, facts)?;
    let context = FeasibilityContext::empty().with_vj_productive(vj_productive);
    if context.is_empty() {
        None
    } else {
        Some(context)
    }
}

fn build_vj_productive_feasibility(
    refdata: &RefDataConfig,
    facts: &CompileFactIndex,
) -> Option<VjProductiveFeasibility> {
    Some(VjProductiveFeasibility::new(
        allele_domain_for_feasibility(Segment::V, refdata, facts)?,
        allele_domain_for_feasibility(Segment::J, refdata, facts)?,
        trim_domain_for_feasibility(Segment::V, TrimEnd::Five, facts)?,
        trim_domain_for_feasibility(Segment::V, TrimEnd::Three, facts)?,
        trim_domain_for_feasibility(Segment::J, TrimEnd::Five, facts)?,
        trim_domain_for_feasibility(Segment::J, TrimEnd::Three, facts)?,
        np_length_domain_for_feasibility(Segment::Np1, facts)?,
    ))
}

fn allele_domain_for_feasibility(
    segment: Segment,
    refdata: &RefDataConfig,
    facts: &CompileFactIndex,
) -> Option<ChoiceDomain<AlleleId>> {
    let support = facts.allele_supports.get(&segment)?;
    let values = unique_allele_ids(allele_ids_for_segment(segment, refdata, facts));
    if values.is_empty() {
        None
    } else {
        Some(ChoiceDomain::new(Some(support.pass_index), values))
    }
}

fn trim_domain_for_feasibility(
    segment: Segment,
    end: TrimEnd,
    facts: &CompileFactIndex,
) -> Option<ChoiceDomain<u32>> {
    let Some(entry) = facts.trim_supports.get(&(segment, end)) else {
        return Some(ChoiceDomain::new(None, vec![0]));
    };

    let IntegerSupport::Enumerated(values) = &entry.value else {
        return None;
    };

    let values = unique_u32_values(
        values
            .iter()
            .copied()
            .filter(|value| (0..=u16::MAX as i64).contains(value))
            .map(|value| value as u32),
    );
    if values.is_empty() {
        None
    } else {
        Some(ChoiceDomain::new(Some(entry.pass_index), values))
    }
}

fn np_length_domain_for_feasibility(
    segment: Segment,
    facts: &CompileFactIndex,
) -> Option<Vec<u32>> {
    let entry = facts.np_length_supports.get(&segment)?;
    let IntegerSupport::Enumerated(values) = &entry.value else {
        return None;
    };

    let lengths = unique_u32_values(
        values
            .iter()
            .copied()
            .filter(|value| (0..=u32::MAX as i64).contains(value))
            .map(|value| value as u32),
    );

    if lengths.is_empty() {
        None
    } else {
        Some(lengths)
    }
}

fn unique_allele_ids(values: Vec<AlleleId>) -> Vec<AlleleId> {
    let mut seen = HashSet::new();
    let mut unique = Vec::new();
    for value in values {
        if seen.insert(value) {
            unique.push(value);
        }
    }
    unique
}

fn unique_u32_values(values: impl IntoIterator<Item = u32>) -> Vec<u32> {
    let mut seen = HashSet::new();
    let mut unique = Vec::new();
    for value in values {
        if seen.insert(value) {
            unique.push(value);
        }
    }
    unique
}
