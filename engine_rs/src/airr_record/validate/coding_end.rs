//! Independent checks of the public retained-outcome extension against source
//! alleles, persistent choices, observed pool flags, and the P ancestry ledger.

use std::collections::HashSet;

use serde_json::{json, Value};

use crate::address::PEnd;
use crate::ir::{complement_base, flag, NucHandle, Segment, Simulation};
use crate::refdata::RefDataConfig;

use super::{AirrRecord, RecordValidationIssue};

fn mismatch(issues: &mut Vec<RecordValidationIssue>, field: &str, reason: &str) {
    issues.push(RecordValidationIssue::CodingEndProcessingMismatch {
        field: field.to_owned(),
        reason: reason.to_owned(),
    });
}

fn check_value(
    actual: Option<&Value>,
    expected: Value,
    field: &str,
    issues: &mut Vec<RecordValidationIssue>,
) {
    if actual != Some(&expected) {
        mismatch(
            issues,
            field,
            "reported value disagrees with persistent ancestry",
        );
    }
}

pub(super) fn check(
    record: &AirrRecord,
    sim: &Simulation,
    refdata: &RefDataConfig,
    issues: &mut Vec<RecordValidationIssue>,
) {
    let Some(provenance) = sim.coding_end_provenance.as_ref() else {
        if record.coding_end_processing.is_some() {
            mismatch(issues, "coding_end_processing", "model was not used");
        }
        return;
    };
    let Some(metadata) = &record.coding_end_processing else {
        mismatch(
            issues,
            "coding_end_processing",
            "model provenance is missing",
        );
        return;
    };
    for (field, expected) in [
        ("model", "original_coding_end_v1"),
        ("boundary_contract", "provided_allele_coding_ends"),
        ("coordinate_system", "zero_based_coding_oriented"),
    ] {
        check_value(metadata.get(field), json!(expected), field, issues);
    }
    let complete = [Segment::V, Segment::D, Segment::J]
        .into_iter()
        .filter(|segment| *segment != Segment::D || sim.assignments.has(Segment::D))
        .all(|segment| {
            sim.sequence
                .regions
                .iter()
                .any(|region| region.segment == segment)
        });
    check_value(
        metadata.get("recombination_complete"),
        json!(complete),
        "recombination_complete",
        issues,
    );
    let choices = metadata.get("choices").and_then(Value::as_array);
    if choices.map(Vec::len) != Some(provenance.choices.len()) {
        mismatch(issues, "choices", "choice history length differs");
    }
    for (index, choice) in provenance.choices.iter().enumerate() {
        let Some(allele) = refdata.get(choice.segment, choice.allele_id) else {
            mismatch(issues, "choices", "original source allele is missing");
            continue;
        };
        let value = choices.and_then(|values| values.get(index));
        let expected = json!({
            "choice_id": index,
            "segment": format!("{:?}", choice.segment),
            "original_allele": allele.name,
            "outcome_index": choice.outcome_index,
            "five_prime_outcome": choice.five,
            "three_prime_outcome": choice.three,
            "net_deletion_5": if choice.five > 0 { choice.five } else { 0 },
            "net_deletion_3": if choice.three > 0 { choice.three } else { 0 },
            "retained_p_5": if choice.five < 0 { choice.five.unsigned_abs() } else { 0 },
            "retained_p_3": if choice.three < 0 { choice.three.unsigned_abs() } else { 0 },
        });
        check_value(value, expected, &format!("choices[{index}]"), issues);
        if choice.choice_id != index
            || choice.five.unsigned_abs() > allele.len() as u64
            || choice.three.unsigned_abs() > allele.len() as u64
            || i128::from(choice.five.max(0)) + i128::from(choice.three.max(0))
                > i128::from(allele.len())
            || (choice.segment == Segment::V && choice.five != 0)
            || (choice.segment == Segment::J && choice.three != 0)
        {
            mismatch(issues, "choices", "invalid original coding-end geometry");
        }
    }

    let origins = metadata.get("p_origins").and_then(Value::as_array);
    if origins.map(Vec::len) != Some(provenance.p_origins.len()) {
        mismatch(issues, "p_origins", "P ancestry history length differs");
    }
    let mut observed_positions = HashSet::new();
    let mut original_templates = HashSet::new();
    let mut retained = [0usize; 4];
    let mut observed = [0usize; 4];
    for (index, origin) in provenance.p_origins.iter().enumerate() {
        let (end_index, assembled_name, physical_five) = match origin.end {
            PEnd::V3 => (0, "V_3", false),
            PEnd::D5 => (1, "D_5", true),
            PEnd::D3 => (2, "D_3", false),
            PEnd::J5 => (3, "J_5", true),
        };
        retained[end_index] += 1;
        let canonical_five = physical_five != origin.orientation.is_reverse();
        if !original_templates.insert((origin.choice_id, canonical_five, origin.template_position))
        {
            mismatch(
                issues,
                "p_origins",
                "duplicate original coding-tip template ancestry",
            );
        }
        let source = refdata.get(origin.source_segment, origin.allele_id);
        let current = origin
            .query_position
            .and_then(|position| sim.pool.get(NucHandle::new(position)));
        let expected = json!({
            "origin_id": index,
            "choice_id": origin.choice_id,
            "source_segment": format!("{:?}", origin.source_segment),
            "original_allele": source.map(|a| a.name.as_str()).unwrap_or(""),
            "canonical_end": if canonical_five { "5" } else { "3" },
            "assembled_end": assembled_name,
            "source_orientation": if origin.orientation.is_reverse() { "reverse_complement" } else { "forward" },
            "template_position": origin.template_position,
            "query_position": origin.query_position,
            "original_base": (origin.original_base as char).to_string(),
            "current_base": current.map(|n| (n.base as char).to_string()),
        });
        check_value(
            origins.and_then(|values| values.get(index)),
            expected,
            &format!("p_origins[{index}]"),
            issues,
        );
        let choice = provenance.choices.get(origin.choice_id);
        let ancestry_valid = choice.is_some_and(|choice| {
            choice.allele_id == origin.allele_id
                && choice.segment == origin.source_segment
                && origin.end.source_segment() == origin.source_segment
                && source.is_some_and(|allele| {
                    let outcome = if canonical_five {
                        choice.five
                    } else {
                        choice.three
                    };
                    let template = origin.template_position as usize;
                    let in_tip = outcome < 0
                        && template < allele.seq.len()
                        && if canonical_five {
                            (template as u64) < outcome.unsigned_abs()
                        } else {
                            (allele.seq.len() - template) as u64 <= outcome.unsigned_abs()
                        };
                    in_tip
                        && matches!(
                            allele.seq[template].to_ascii_uppercase(),
                            b'A' | b'C' | b'G' | b'T'
                        )
                        && origin.original_base
                            == if origin.orientation.is_reverse() {
                                allele.seq[template]
                            } else {
                                complement_base(allele.seq[template])
                            }
                })
        });
        if !ancestry_valid {
            mismatch(
                issues,
                "p_origins",
                "P ancestry is inconsistent with the original coding tip",
            );
        }
        if let Some(position) = origin.query_position {
            observed[end_index] += 1;
            if !observed_positions.insert(position)
                || current.is_none_or(|n| {
                    !n.flags.contains(flag::P_NUC)
                        || n.flags.contains(flag::N_NUC)
                        || n.germline != origin.original_base
                })
            {
                mismatch(
                    issues,
                    "p_origins",
                    "current P position or nucleotide ancestry is invalid",
                );
            }
        }
    }
    for (index, end) in ["V_3", "D_5", "D_3", "J_5"].iter().enumerate() {
        for (field, count) in [
            ("retained_p_at_join", retained[index]),
            ("observed_p", observed[index]),
        ] {
            check_value(
                metadata.get(field).and_then(|counts| counts.get(*end)),
                json!(count),
                &format!("{field}.{end}"),
                issues,
            );
        }
    }
    let flagged_positions: HashSet<u32> = sim
        .pool
        .as_slice()
        .iter()
        .enumerate()
        .filter_map(|(position, n)| n.flags.contains(flag::P_NUC).then_some(position as u32))
        .collect();
    if flagged_positions != observed_positions {
        mismatch(
            issues,
            "p_origins",
            "observed P flags and ancestry positions differ",
        );
    }

    // A complete recombination has emitted every P side before the final J
    // region closes. Partial native outcomes may legitimately have a choice
    // whose emit pass has not run yet. Regions survive full observation loss.
    if complete {
        let mut expected_templates = HashSet::new();
        for choice in &provenance.choices {
            let Some(allele) = refdata.get(choice.segment, choice.allele_id) else {
                continue;
            };
            for (five, outcome) in [(true, choice.five), (false, choice.three)] {
                if outcome >= 0 || outcome.unsigned_abs() > allele.len() as u64 {
                    continue;
                }
                let p_len = outcome.unsigned_abs() as u32;
                let start = if five { 0 } else { allele.len() - p_len };
                for template in start..start + p_len {
                    expected_templates.insert((choice.choice_id, five, template));
                }
            }
        }
        if original_templates != expected_templates {
            mismatch(
                issues,
                "p_origins",
                "completed emission does not account for every selected original-tip P base",
            );
        }
        for (index, reported) in [
            record.p_v_3_length,
            record.p_d_5_length,
            record.p_d_3_length,
            record.p_j_5_length,
        ]
        .iter()
        .enumerate()
        {
            if *reported != retained[index] as i64 {
                mismatch(
                    issues,
                    "retained_p_at_join",
                    "historical P counters disagree with emitted ancestry",
                );
            }
        }
    }
}
