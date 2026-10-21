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
use crate::address::PEnd;
use crate::assignment::TrimEnd;
use crate::contract::{ContractKind, ContractSet};
use crate::feasibility::{
    ChoiceDomain, CodingEndGeometry, CodingSegmentFeasibility, FeasibilityContext,
    NpCompletionDomain, PProductiveFeasibility, SegmentCompletionDomain, VjProductiveFeasibility,
};
use crate::ir::Segment;
use crate::pass::{AlleleIdSupport, IntegerSupport, NpBaseSupport};
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
    let refdata = refdata?;
    let mut context = FeasibilityContext::empty();
    for (&segment, entry) in &facts.coding_end_supports {
        context = context.with_coding_end_geometry(CodingEndGeometry {
            segment,
            allele_pass: facts
                .allele_supports
                .get(&segment)
                .map(|entry| entry.pass_index),
            outcomes: ChoiceDomain::new(Some(entry.pass_index), entry.value.clone()),
        });
    }
    let Some(contracts) = contracts else {
        return (!context.is_empty()).then_some(context);
    };
    if contracts.contains_kind(ContractKind::CodingSequenceIntegrity) {
        for segment in [Segment::V, Segment::J] {
            if facts.coding_end_supports.contains_key(&segment) {
                continue;
            }
            if let (Some(alleles), Some(trim_5), Some(trim_3)) = (
                allele_domain_for_feasibility(segment, refdata, facts),
                trim_domain_for_feasibility(segment, TrimEnd::Five, facts),
                trim_domain_for_feasibility(segment, TrimEnd::Three, facts),
            ) {
                context = context.with_coding_segment(CodingSegmentFeasibility::new(
                    segment, alleles, trim_5, trim_3,
                ));
            }
        }
    }
    let p_aware = facts.needs_complete_productive_support(contracts);
    if p_aware {
        if let Some(feasibility) =
            build_p_productive_feasibility(contracts, refdata, facts, &mut Vec::new())
        {
            context = context.with_p_productive(feasibility);
        }
    } else if contracts.contains_kind(ContractKind::ProductiveJunctionFrame)
        && !facts.assigned_segments.contains(&Segment::D)
    {
        if let Some(mut vj_productive) = build_vj_productive_feasibility(refdata, facts) {
            if contracts.contains_kind(ContractKind::CodingSequenceIntegrity) {
                vj_productive = vj_productive.with_coding_integrity();
            }
            context = context.with_vj_productive(vj_productive);
        }
    }
    if context.is_empty() {
        None
    } else {
        Some(context)
    }
}

pub(super) fn validate_p_productive_supports(
    contracts: &ContractSet,
    refdata: Option<&RefDataConfig>,
    facts: &CompileFactIndex,
    errors: &mut Vec<CompileError>,
) {
    let Some(refdata) = refdata else {
        errors.push(plan_scoped_error(contract_precondition(
            "productive_junction_frame",
            "P-aware productive feasibility requires reference data",
        )));
        return;
    };
    if let Some(feasibility) = build_p_productive_feasibility(contracts, refdata, facts, errors) {
        if !feasibility.has_completion(&crate::ir::Simulation::new(), refdata) {
            errors.push(plan_scoped_error(contract_precondition(
                "productive_junction_frame",
                "declared allele, trim, N and P supports have no productive completion",
            )));
        }
    }
}

fn build_p_productive_feasibility(
    contracts: &ContractSet,
    refdata: &RefDataConfig,
    facts: &CompileFactIndex,
    errors: &mut Vec<CompileError>,
) -> Option<PProductiveFeasibility> {
    let initial_errors = errors.len();
    let unsupported = |reason: String| super::error::CompileErrorKind::UnsupportedContractCombo {
        contract_name: "productive_junction_frame".into(),
        reason,
    };
    let mut p_lengths = Vec::new();
    let processing = facts
        .coding_end_supports
        .iter()
        .map(|(&segment, entry)| {
            (
                segment,
                ChoiceDomain::new(Some(entry.pass_index), entry.value.clone()),
            )
        })
        .collect();
    let mut ends = HashSet::new();
    let mut ordered = Vec::new();
    for (end, entry) in &facts.p_length_supports {
        if !ends.insert(*end) {
            errors.push(entry.error(unsupported(format!("P-aware productive feasibility requires one P-length pass per coding end; repeated {end:?}"))));
            continue;
        }
        let Some(values) = completion_integer_values(&entry.value) else {
            errors.push(entry.error(unsupported(format!("P-aware productive feasibility requires finite nonnegative enumerable P-length support at {end:?}"))));
            continue;
        };
        p_lengths.push((*end, ChoiceDomain::new(Some(entry.pass_index), values)));
        ordered.push((
            entry.pass_index,
            match end {
                PEnd::V3 => 1,
                PEnd::D5 => 4,
                PEnd::D3 => 6,
                PEnd::J5 => 8,
            },
        ));
        if !facts.allele_supports.contains_key(&end.source_segment()) {
            errors.push(entry.error(unsupported(format!(
                "P-aware productive feasibility requires an allele-sampling pass for {:?}",
                end.source_segment()
            ))));
        }
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
    let uses_d = facts.assigned_segments.contains(&Segment::D);
    let mut segments = Vec::new();
    for segment in [Segment::V, Segment::D, Segment::J] {
        if segment == Segment::D && !uses_d {
            continue;
        }
        let (Some(alleles), Some(trim_5), Some(trim_3)) = (
            allele_domain_for_feasibility(segment, refdata, facts),
            trim_domain_for_feasibility(segment, TrimEnd::Five, facts),
            trim_domain_for_feasibility(segment, TrimEnd::Three, facts),
        ) else {
            errors.push(plan_scoped_error(unsupported(format!("P-aware productive feasibility requires finite allele and trim support for {segment:?}"))));
            continue;
        };
        let assembly = facts
            .append_passes
            .get(&segment)
            .filter(|passes| passes.len() == 1)
            .map(|passes| passes[0]);
        if let Some(index) = assembly {
            ordered.push((
                index,
                match segment {
                    Segment::V => 0,
                    Segment::D => 5,
                    _ => 9,
                },
            ));
            for (end, trim) in [(TrimEnd::Five, &trim_5), (TrimEnd::Three, &trim_3)] {
                let earliest_p = facts
                    .p_length_supports
                    .iter()
                    .filter(|(e, _)| e.source_segment() == segment)
                    .map(|(_, e)| e.pass_index)
                    .min()
                    .unwrap_or(index)
                    .min(index);
                if trim.pass_index().is_some_and(|t| t >= earliest_p) {
                    errors.push(plan_scoped_error(unsupported(format!(
                        "{segment:?} {end:?} trimming must precede P templating and assembly"
                    ))));
                }
            }
        } else {
            errors.push(plan_scoped_error(unsupported(format!(
                "P-aware productive feasibility requires one assembly pass for {segment:?}"
            ))));
        }
        segments.push(SegmentCompletionDomain {
            segment,
            alleles,
            trim_5,
            trim_3,
        });
    }
    let mut np = Vec::new();
    for segment in [Segment::Np1, Segment::Np2] {
        let Some(entry) = facts.np_length_supports.get(&segment) else {
            continue;
        };
        if !facts
            .append_passes
            .get(&segment)
            .is_some_and(|passes| passes.as_slice() == [entry.pass_index])
        {
            errors.push(entry.error(unsupported(format!(
                "P-aware productive feasibility requires one N-generation pass for {segment:?}"
            ))));
        }
        if segment == Segment::Np2 && !uses_d {
            errors.push(entry.error(unsupported(
                "NP2 is not part of a VJ productive junction".into(),
            )));
        }
        let Some(values) = completion_integer_values(&entry.value) else {
            errors.push(entry.error(unsupported(format!("P-aware productive feasibility requires finite enumerable {segment:?} length support"))));
            continue;
        };
        let bases = facts
            .np_base_supports
            .get(&segment)
            .and_then(|entry| entry.value)
            .or_else(|| {
                values.iter().all(|n| *n == 0).then_some(NpBaseSupport {
                    first: 15,
                    transitions: [15; 4],
                })
            });
        let Some(bases) = bases.filter(|b| {
            b.first != 0
                && b.first & !15 == 0
                && b.transitions.iter().all(|m| *m != 0 && *m & !15 == 0)
        }) else {
            errors.push(entry.error(unsupported(format!("P-aware productive feasibility requires declared canonical first-order N-base support for {segment:?}"))));
            continue;
        };
        ordered.push((
            entry.pass_index,
            if segment == Segment::Np1 { 2 } else { 7 },
        ));
        np.push(NpCompletionDomain {
            segment,
            lengths: ChoiceDomain::new(Some(entry.pass_index), values),
            bases,
        });
    }
    let d_orientation = match facts.d_orientation_supports.as_slice() {
        [] => ChoiceDomain::new(None, vec![false]),
        [entry] if !entry.value.is_empty() => {
            ordered.push((entry.pass_index, 3));
            ChoiceDomain::new(Some(entry.pass_index), entry.value.clone())
        }
        _ => {
            errors.push(plan_scoped_error(unsupported(
                "P-aware productive feasibility requires at most one D-orientation decision".into(),
            )));
            ChoiceDomain::new(None, vec![false])
        }
    };
    ordered.sort_unstable();
    if ordered.windows(2).any(|pair| pair[0].1 > pair[1].1) {
        errors.push(plan_scoped_error(unsupported("P-aware productive feasibility requires coding-order assembly: V, V3-P, N1, D orientation/D5-P/D/D3-P, N2, J5-P, J".into())));
    }
    if errors.len() != initial_errors {
        return None;
    }
    Some(PProductiveFeasibility {
        segments,
        processing,
        p_lengths,
        np,
        d_orientation,
        check_frame: contracts.contains_kind(ContractKind::ProductiveJunctionFrame)
            || contracts.contains_kind(ContractKind::CodingSequenceIntegrity),
        check_stops: contracts.contains_kind(ContractKind::NoStopCodonInJunction)
            || contracts.contains_kind(ContractKind::CodingSequenceIntegrity),
        check_coding: contracts.contains_kind(ContractKind::CodingSequenceIntegrity),
        check_v_anchor: contracts.contains_kind(ContractKind::CanonicalAnchor {
            segment: Segment::V,
        }),
        check_j_anchor: contracts.contains_kind(ContractKind::CanonicalAnchor {
            segment: Segment::J,
        }),
    })
}

fn completion_integer_values(support: &IntegerSupport) -> Option<Vec<u32>> {
    let IntegerSupport::Enumerated(values) = support else {
        return None;
    };
    if values.is_empty() || values.iter().any(|n| !(0..=u32::MAX as i64).contains(n)) {
        return None;
    }
    Some(unique_u32_values(values.iter().map(|n| *n as u32)))
}

pub(super) fn validate_coding_supports(
    refdata: Option<&RefDataConfig>,
    facts: &CompileFactIndex,
    errors: &mut Vec<CompileError>,
) {
    let Some(refdata) = refdata else { return };
    for segment in [Segment::V, Segment::J] {
        // The joint end-outcome solver validates these without reducing
        // correlated outcomes to independent trim marginals.
        if facts.coding_end_supports.contains_key(&segment) {
            continue;
        }
        if !facts.allele_supports.contains_key(&segment) {
            continue;
        }
        let (Some(five), Some(three)) = (
            trim_values(segment, TrimEnd::Five, facts, errors),
            trim_values(segment, TrimEnd::Three, facts, errors),
        ) else {
            continue;
        };
        let feasible = allele_ids_for_segment(segment, refdata, facts)
            .into_iter()
            .any(|id| {
                let Some(allele) = refdata.get(segment, id) else {
                    return false;
                };
                let Some(anchor) = allele.anchor.map(u32::from) else {
                    return false;
                };
                five.iter().any(|&t5| {
                    three.iter().any(|&t3| {
                        t5 <= anchor
                            && allele
                                .seq
                                .len()
                                .checked_sub(t3 as usize)
                                .is_some_and(|end| anchor as usize + 3 <= end)
                            && crate::productivity::retained_stop_free(refdata, segment, id, t5, t3)
                    })
                })
            });
        if !feasible {
            errors.push(plan_scoped_error(contract_precondition(
                "coding_sequence_integrity",
                format!(
                    "{segment:?} allele/trim support has no stop-free retained coding sequence"
                ),
            )));
        }
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
