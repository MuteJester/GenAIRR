//! Projection of retained outcomes and P ancestry; no microscopic processing
//! events are inferred from these final outcomes.

use serde_json::{json, Value};

use crate::address::PEnd;
use crate::ir::{CodingEndPOrigin, NucHandle, Segment, Simulation};
use crate::refdata::RefDataConfig;

pub(crate) const MODEL: &str = "original_coding_end_v1";
pub(crate) const BOUNDARY_CONTRACT: &str = "provided_allele_coding_ends";
pub(crate) const COORDINATE_SYSTEM: &str = "zero_based_coding_oriented";

pub(crate) fn end_name(end: PEnd) -> &'static str {
    match end {
        PEnd::V3 => "V_3",
        PEnd::D5 => "D_5",
        PEnd::D3 => "D_3",
        PEnd::J5 => "J_5",
    }
}

pub(crate) fn segment_name(segment: Segment) -> &'static str {
    match segment {
        Segment::V => "V",
        Segment::D => "D",
        Segment::J => "J",
        Segment::Np1 => "NP1",
        Segment::Np2 => "NP2",
    }
}

pub(crate) fn canonical_end(origin: &CodingEndPOrigin) -> &'static str {
    let assembled_five = matches!(origin.end, PEnd::D5 | PEnd::J5);
    if assembled_five ^ origin.orientation.is_reverse() {
        "5"
    } else {
        "3"
    }
}

pub(crate) fn project(sim: &Simulation, refdata: &RefDataConfig) -> Option<Value> {
    let provenance = sim.coding_end_provenance.as_ref()?;
    let recombination_complete = [Segment::V, Segment::J]
        .into_iter()
        .chain(sim.assignments.has(Segment::D).then_some(Segment::D))
        .all(|segment| {
            sim.sequence
                .regions
                .iter()
                .any(|region| region.segment == segment)
        });
    let choices: Vec<Value> = provenance
        .choices
        .iter()
        .map(|choice| {
            let original_allele = refdata
                .get(choice.segment, choice.allele_id)
                .map(|allele| allele.name.as_str())
                .unwrap_or("");
            json!({
                "choice_id": choice.choice_id,
                "segment": segment_name(choice.segment),
                "original_allele": original_allele,
                "outcome_index": choice.outcome_index,
                "five_prime_outcome": choice.five,
                "three_prime_outcome": choice.three,
                "net_deletion_5": choice.five.max(0),
                "net_deletion_3": choice.three.max(0),
                "retained_p_5": choice.five.min(0).unsigned_abs(),
                "retained_p_3": choice.three.min(0).unsigned_abs(),
            })
        })
        .collect();
    let origins: Vec<Value> = provenance
        .p_origins
        .iter()
        .enumerate()
        .map(|(origin_id, origin)| {
            let allele = refdata.get(origin.source_segment, origin.allele_id);
            let current_base = origin
                .query_position
                .and_then(|position| sim.pool.get(NucHandle::new(position)))
                .map(|n| (n.base as char).to_string());
            json!({
                "origin_id": origin_id,
                "choice_id": origin.choice_id,
                "source_segment": segment_name(origin.source_segment),
                "original_allele": allele.map(|a| a.name.as_str()).unwrap_or(""),
                "canonical_end": canonical_end(origin),
                "assembled_end": end_name(origin.end),
                "source_orientation": if origin.orientation.is_reverse() { "reverse_complement" } else { "forward" },
                "template_position": origin.template_position,
                "query_position": origin.query_position,
                "original_base": (origin.original_base as char).to_string(),
                "current_base": current_base,
            })
        })
        .collect();
    let mut retained = serde_json::Map::new();
    let mut observed = serde_json::Map::new();
    for end in PEnd::all() {
        retained.insert(
            end_name(end).into(),
            json!(provenance.p_origins.iter().filter(|o| o.end == end).count()),
        );
        observed.insert(
            end_name(end).into(),
            json!(provenance
                .p_origins
                .iter()
                .filter(|o| o.end == end && o.query_position.is_some())
                .count()),
        );
    }
    Some(json!({
        "model": MODEL,
        "boundary_contract": BOUNDARY_CONTRACT,
        "coordinate_system": COORDINATE_SYSTEM,
        "recombination_complete": recombination_complete,
        "choices": choices,
        "retained_p_at_join": retained,
        "observed_p": observed,
        "p_origins": origins,
    }))
}
