//! Projection of a current observation without rewriting its founder duplex.
pub(crate) mod layout;

use super::projection::{self, BodyPoint, Path, ProjectionAnalysis};
use super::Substrate;
use crate::ir::{
    complement_base, DuplexEditAction, DuplexQueryLink, DuplexRemovalCause, Segment, Simulation,
};
use crate::refdata::RefDataConfig;
use serde_json::{json, Map, Value};
use sha2::{Digest, Sha256};
use std::collections::BTreeSet;

fn label(segment: Segment) -> &'static str {
    match segment {
        Segment::V => "v",
        Segment::D => "d",
        Segment::J => "j",
        _ => "np",
    }
}

fn append(runs: &mut Vec<(usize, char)>, op: char) {
    if let Some((n, last)) = runs.last_mut() {
        if *last == op {
            *n += 1;
            return;
        }
    }
    runs.push((1, op));
}

fn current_path(sim: &Simulation, segment: Segment, reference: &[u8]) -> Result<Path, String> {
    let provenance = sim
        .duplex_provenance
        .as_ref()
        .ok_or("simulation has no duplex ancestry")?;
    let inst = sim
        .assignments
        .get(segment)
        .ok_or("duplex assignment is missing")?;
    let region = sim
        .sequence
        .regions
        .iter()
        .find(|r| r.segment == segment)
        .ok_or("duplex structural region is missing")?;
    if sim
        .sequence
        .regions
        .iter()
        .filter(|r| r.segment == segment)
        .count()
        != 1
    {
        return Err("duplex projection requires one region per source".into());
    }
    let (start, end) = (region.start.index() as usize, region.end.index() as usize);
    if start > end || end > sim.pool.len() {
        return Err("duplex structural region is out of bounds".into());
    }
    let mut path = Path {
        available: false,
        reason: None,
        body: Vec::new(),
        query_start: None,
        query_end: None,
        reference_start: None,
        reference_end: None,
        oriented_reference_start: None,
        oriented_reference_end: None,
        cigar: None,
        sequence_alignment: None,
        germline_alignment: None,
        identity: None,
        empty_query_boundary: None,
    };
    let gstart = inst.observed_ref_start.unwrap_or(u32::from(inst.trim_5)) as usize;
    let gend = inst
        .observed_ref_end
        .map(|x| x as usize)
        .unwrap_or(reference.len() - usize::from(inst.trim_3));
    if gstart >= gend {
        path.reason = Some("source_reference_not_observed");
        return Ok(path);
    }
    if gend > reference.len() {
        return Err("duplex observed reference interval is out of bounds".into());
    }
    let inverted = inst.orientation.is_reverse();
    let rstart = if inverted {
        reference.len() - gend
    } else {
        gstart
    };
    let rend = if inverted {
        reference.len() - gstart
    } else {
        gend
    };
    let base_at = |r: usize| {
        let g = if inverted { reference.len() - 1 - r } else { r };
        if inverted {
            complement_base(reference[g])
        } else {
            reference[g].to_ascii_uppercase()
        }
    };
    let mut r = rstart;
    let mut qa = Vec::new();
    let mut ra = Vec::new();
    let mut runs = Vec::new();
    let mut matches = 0;
    for q in start..end {
        let n = &sim.pool.as_slice()[q];
        if let Some(g) = provenance.body_reference_position(segment, q) {
            let expected = if inverted { reference.len() - 1 - g } else { g };
            if expected < r || expected >= rend {
                return Err("duplex observed body path is nonmonotone or outside coverage".into());
            }
            while r < expected {
                qa.push(b'-');
                ra.push(base_at(r));
                append(&mut runs, 'D');
                r += 1;
            }
            let germline = base_at(r);
            qa.push(n.base);
            ra.push(germline);
            append(&mut runs, 'M');
            matches += usize::from(n.base.eq_ignore_ascii_case(&germline));
            r += 1;
            let DuplexQueryLink::FounderColumn(column) = provenance.links[q] else {
                return Err("inserted observation acquired a founder reference coordinate".into());
            };
            let claim = provenance.founder.body_claims[column]
                .iter()
                .find(|c| {
                    c.reference_position == g
                        && match (segment, c.segment) {
                            (Segment::V, Substrate::V)
                            | (Segment::D, Substrate::D)
                            | (Segment::J, Substrate::J) => true,
                            _ => false,
                        }
                })
                .ok_or("body support is missing")?;
            path.body.push(BodyPoint {
                query_position: q,
                reference_position: g,
                support_ids: claim.support_ids.clone(),
            });
        } else {
            qa.push(n.base);
            ra.push(b'-');
            append(&mut runs, 'I');
        }
    }
    while r < rend {
        qa.push(b'-');
        ra.push(base_at(r));
        append(&mut runs, 'D');
        r += 1;
    }
    path.available = true;
    if start < end {
        path.query_start = Some(start);
        path.query_end = Some(end);
    } else {
        path.empty_query_boundary = Some(start);
    }
    path.reference_start = Some(gstart);
    path.reference_end = Some(gend);
    path.oriented_reference_start = Some(rstart);
    path.oriented_reference_end = Some(rend);
    path.cigar = Some(runs.iter().map(|(n, op)| format!("{n}{op}")).collect());
    path.identity = (!qa.is_empty()).then_some(matches as f64 / qa.len() as f64);
    path.sequence_alignment = Some(String::from_utf8(qa).map_err(|_| "non-ASCII query")?);
    path.germline_alignment = Some(String::from_utf8(ra).map_err(|_| "non-ASCII reference")?);
    Ok(path)
}

/// Project a persistent duplex-backed simulation after audited continuations.
pub fn project_duplex_simulation(
    sim: &Simulation,
    refdata: &RefDataConfig,
) -> Result<Value, String> {
    let p = sim
        .duplex_provenance
        .as_ref()
        .ok_or("simulation has no duplex repair provenance")?;
    p.founder
        .bindings
        .validate(&p.founder.annotated.product, refdata)?;
    p.validate_pool(sim.pool.as_slice())?;
    layout::validate_layout(sim)?;
    let founder_analysis =
        projection::analyze_coupled(&p.founder.annotated, &p.founder.bindings, refdata)?;
    if !founder_analysis.layout_reasons.is_empty() {
        return Err("duplex founder is not a supported execution layout".into());
    }
    let query: Vec<_> = sim.pool.as_slice().iter().map(|n| n.base).collect();
    if query.iter().any(|b| {
        !matches!(
            b.to_ascii_uppercase(),
            b'A' | b'C'
                | b'G'
                | b'T'
                | b'N'
                | b'R'
                | b'Y'
                | b'S'
                | b'W'
                | b'K'
                | b'M'
                | b'B'
                | b'D'
                | b'H'
                | b'V'
        )
    }) {
        return Err("duplex query contains a non-DNA nucleotide".into());
    }
    let substrates = [Substrate::V, Substrate::D, Substrate::J];
    let segments = [Segment::V, Segment::D, Segment::J];
    let mut built = Vec::new();
    for (i, s) in segments.into_iter().enumerate() {
        let binding = p.founder.bindings.for_substrate(substrates[i]);
        let inst = sim
            .assignments
            .get(s)
            .ok_or("duplex generating assignment is missing")?;
        let source = &founder_analysis.paths[i];
        if inst.allele_id != p.founder.bindings.allele_id(substrates[i])
            || inst.receptor_revision_original_id.is_some()
            || inst.orientation.is_reverse() != (i == 1 && p.founder.annotated.product.d_inverted)
            || Some(usize::from(inst.trim_5)) != source.reference_start
            || source.reference_end != binding.sequence.len().checked_sub(usize::from(inst.trim_3))
        {
            return Err("duplex assignment no longer matches its generating source".into());
        }
        built.push(current_path(sim, s, binding.sequence.as_bytes())?);
    }
    let paths: [Path; 3] = built.try_into().map_err(|_| "missing duplex source path")?;
    let mut reasons = BTreeSet::new();
    for (s, path) in segments.into_iter().zip(&paths) {
        if !path.available {
            reasons.insert(format!("{}:{}", label(s), path.reason.unwrap()));
        }
    }
    let mut previous = 0;
    for path in &paths {
        if let Some((start, end)) = path
            .query_start
            .or(path.empty_query_boundary)
            .zip(path.query_end.or(path.empty_query_boundary))
        {
            if start < previous {
                return Err("duplex current source regions overlap or are reordered".into());
            }
            previous = end;
        }
    }
    let current_claims:Vec<_>=p.links.iter().enumerate().map(|(q,link)| {
        let claims=match link { DuplexQueryLink::FounderColumn(column)=>p.founder.body_claims[*column].iter().map(|c|json!({
            "segment":c.segment,"reference_position":c.reference_position,"support_ids":c.support_ids,
            "source_instance":p.founder.bindings.for_substrate(c.segment).source_instance})).collect::<Vec<_>>(),_=>Vec::new() };
        json!({"query_position":q,"claims":claims})
    }).collect();
    let np_bounds: [Option<(usize, usize)>; 2] = std::array::from_fn(|i| {
        paths[i].query_end.or(paths[i].empty_query_boundary).zip(
            paths[i + 1]
                .query_start
                .or(paths[i + 1].empty_query_boundary),
        )
    });
    let analysis = ProjectionAnalysis {
        paths,
        column_claims: current_claims,
        overlaps: founder_analysis.overlaps,
        layout_reasons: reasons,
    };
    let mut value = projection::render_projection(
        &p.founder.annotated,
        &p.founder.bindings,
        refdata,
        &query,
        analysis,
    )?;
    let assessment = crate::productivity::assess(sim, refdata);
    let record = value["record"]
        .as_object_mut()
        .ok_or("missing duplex record")?;
    let mut availability = record["duplex_repair"]["availability"]
        .as_object()
        .cloned()
        .unwrap_or_default();
    let mut put = |field: &str, v: Value, reason: &str| {
        availability.remove(field);
        if v.is_null() {
            availability.insert(field.into(), json!(reason));
        }
        record.insert(field.into(), v);
    };
    for (field, v) in [
        ("productive", json!(assessment.productive)),
        ("vj_in_frame", json!(assessment.vj_in_frame)),
        ("stop_codon", json!(assessment.stop_codon)),
        ("junction_start", json!(assessment.junction.map(|x| x.0))),
        ("junction_end", json!(assessment.junction.map(|x| x.1))),
        (
            "junction_length",
            json!(assessment.junction.map(|(a, b)| b - a)),
        ),
        (
            "junction",
            json!(assessment
                .junction
                .map(|(a, b)| String::from_utf8_lossy(&query[a..b]).into_owned())),
        ),
        (
            "junction_aa",
            json!(assessment.junction.map(|(a, b)| query[a..b]
                .chunks_exact(3)
                .map(|c| crate::codon::translate_codon(c[0], c[1], c[2]) as char)
                .collect::<String>())),
        ),
        (
            "sequence_aa",
            json!(assessment.frame.map(|f| query
                .get(f.start..)
                .unwrap_or_default()
                .chunks_exact(3)
                .map(|c| crate::codon::translate_codon(c[0], c[1], c[2]) as char)
                .collect::<String>())),
        ),
    ] {
        put(field, v, "productivity_evidence_unresolved");
    }
    put("productivity_reasons", json!(assessment.reasons), "");
    for (i, s) in segments.into_iter().enumerate() {
        let source = &founder_analysis.paths[i];
        put(
            &format!("{}_trim_5", label(s)),
            json!(source.reference_start),
            "founder_trim_unresolved",
        );
        put(
            &format!("{}_trim_3", label(s)),
            json!(source.reference_end.map(|g| p
                .founder
                .bindings
                .for_substrate(substrates[i])
                .sequence
                .len()
                - g)),
            "founder_trim_unresolved",
        );
    }
    for (i, np) in ["np1", "np2"].into_iter().enumerate() {
        let bounds = np_bounds[i];
        put(
            &format!("{np}_aa"),
            json!(bounds
                .zip(assessment.frame)
                .map(|((a, b), f)| projection::translate_interval(&query, a, b, f.phase))),
            "coding_phase_unresolved",
        );
    }
    for segment in segments {
        let recalled = crate::airr_record::recalled_call_name(
            refdata,
            sim,
            segment,
            sim.assignments.get(segment).map(|a| a.allele_id),
        );
        put(
            &format!("recalled_{}_call", label(segment)),
            if recalled.is_empty() {
                Value::Null
            } else {
                json!(recalled)
            },
            "recall_not_computed",
        );
    }
    drop(put);
    record.insert(
        "productivity_policy".into(),
        json!(crate::productivity::policy_name(sim)),
    );
    record.insert("rev_comp".into(), json!(p.rev_comp));
    if p.rev_comp {
        record.insert(
            "sequence".into(),
            json!(
                String::from_utf8(query.iter().rev().map(|b| complement_base(*b)).collect())
                    .map_err(|_| "non-ASCII reverse read")?
            ),
        );
    }
    let founder_id = record["sequence_id"].clone();
    let identity = serde_json::to_vec(&(
        &founder_id,
        &p.links,
        &p.edits,
        p.rev_comp,
        &p.attempt_stages,
        &query,
    ))
    .map_err(|e| e.to_string())?;
    let identity = if let Some(selection) = &p.founder.selection {
        serde_json::to_vec(&(identity, selection)).map_err(|e| e.to_string())?
    } else {
        identity
    };
    record.insert(
        "sequence_id".into(),
        json!(format!("duplex-observation:{:x}", Sha256::digest(identity))),
    );
    project_counters(sim, record, &mut availability);
    let observation = json!({"model":"duplex_observation_v1","founder_id":founder_id,
        "coordinate_system":"zero_based_coding_oriented","links":p.links,"removed":p.removed,
        "edits":p.edits,"attempt_stages":p.attempt_stages,"rev_comp":p.rev_comp,"insertion_identity_scope":"observation_branch_local"});
    record["duplex_repair"]["observation"] = observation;
    if p.founder.selection.is_some() {
        record["duplex_repair"]["selection"] = json!(p.founder.selection);
    }
    record["duplex_repair"]["founder_column_body_claims"] = json!(founder_analysis.column_claims);
    record["duplex_repair"]["body_claims_timepoint"] = json!("current_observation");
    record["duplex_repair"]["availability"] = json!(availability);
    value["availability"] = json!(availability);
    value["capabilities"]["persistent_simulation"] = json!(true);
    value["capabilities"]["experiment_integration"] = json!(p.founder.selection.is_some());
    Ok(value)
}

fn project_counters(
    sim: &Simulation,
    record: &mut Map<String, Value>,
    availability: &mut Map<String, Value>,
) {
    let p = sim.duplex_provenance.as_ref().unwrap();
    let mut counts = Map::<String, Value>::new();
    let mut unknown = BTreeSet::new();
    let mut shm_written: Vec<(u8, u16)> = Vec::new();
    let mut increment = |field: String| {
        let n = counts.get(&field).and_then(Value::as_u64).unwrap_or(0);
        counts.insert(field, json!(n + 1));
    };
    for edit in &p.edits {
        let pass = edit.source_pass.as_deref();
        match &edit.action {
            DuplexEditAction::Substitution {
                segment,
                old_base,
                new_base,
                germline_pos,
                ..
            } => match pass {
                // A write whose drawn base equalled the base already there is
                // in the ledger -- it is a physical record of every write --
                // but it changed nothing and is not a mutation, so it must
                // not reach these counters. Leaving it in put the per-segment
                // sum one above `n_mutations` in the same projected record:
                // 7 against 6 on a grown duplex family. Case-insensitive, as
                // everywhere: a templated base is stored lowercase and an SHM
                // write emits uppercase.
                Some("mutate.uniform" | "mutate.s5f")
                    if old_base.eq_ignore_ascii_case(new_base) => {}
                Some("mutate.uniform" | "mutate.s5f") => {
                    let s = match segment.as_str() {
                        "V" => "v",
                        "D" => "d",
                        "J" => "j",
                        _ => "np",
                    };
                    increment(format!("n_{s}_mutations"));
                    // Remember the germline position for the recoverable
                    // count below. The NP segments carry no germline
                    // coordinate, so the three templated segments are the
                    // whole key space; `germline_pos` is `None` on the rest
                    // and nothing is collected for them.
                    if let Some(pos) = germline_pos {
                        let code = match segment.as_str() {
                            "V" => crate::ir::Segment::V as u8,
                            "D" => crate::ir::Segment::D as u8,
                            "J" => crate::ir::Segment::J as u8,
                            _ => continue,
                        };
                        shm_written.push((code, *pos));
                    }
                }
                Some("corrupt.pcr") => increment("applied_pcr_substitutions".into()),
                Some("corrupt.quality") => increment("applied_quality_substitutions".into()),
                Some("corrupt.ns") => (),
                _ => {
                    unknown.insert("unattributed_substitution");
                }
            },
            DuplexEditAction::Insertion { segment, .. }
            | DuplexEditAction::Deletion {
                segment,
                cause: DuplexRemovalCause::StructuralDeletion,
                ..
            } => {
                if pass == Some("corrupt.indel") {
                    increment("n_indels".into());
                    if matches!(segment.as_str(), "V" | "D" | "J") {
                        increment(format!("n_{}_indels", segment.to_ascii_lowercase()));
                    }
                } else {
                    unknown.insert("unattributed_structural_edit");
                }
            }
            DuplexEditAction::Deletion {
                cause: DuplexRemovalCause::ObservationLoss,
                ..
            } => match pass {
                Some("corrupt.end_loss.5") => increment("end_loss_5_length".into()),
                Some("corrupt.end_loss.3") => increment("end_loss_3_length".into()),
                _ => {
                    unknown.insert("unattributed_observation_loss");
                }
            },
            DuplexEditAction::ReverseComplement { .. } => (),
            _ => {
                unknown.insert("unsupported_replacement");
            }
        }
    }
    drop(increment);
    for field in [
        "n_pcr_errors",
        "n_quality_errors",
        "n_indels",
        "n_v_indels",
        "n_d_indels",
        "n_j_indels",
        "n_v_mutations",
        "n_d_mutations",
        "n_j_mutations",
        "n_np_mutations",
        "end_loss_5_length",
        "end_loss_3_length",
    ] {
        record.insert(field.into(), counts.remove(field).unwrap_or(json!(0)));
    }
    for (field, stage) in [
        ("n_pcr_errors", crate::ir::DuplexAttemptStage::Pcr),
        ("n_quality_errors", crate::ir::DuplexAttemptStage::Quality),
    ] {
        record.insert(field.into(), json!(p.sampled_attempts(stage)));
    }
    record["duplex_repair"]["applied_error_substitutions"] = json!({
        "pcr":counts.remove("applied_pcr_substitutions").unwrap_or(json!(0)),
        "quality":counts.remove("applied_quality_substitutions").unwrap_or(json!(0))});
    record.insert("n_mutations".into(), json!(sim.mutation_count));
    // Resolved against the observed molecule exactly as the native builder
    // does it, from the SHM sites this ledger attributes to a mutation pass.
    // The field is in the default column order, so a duplex record that
    // omitted it wrote a blank column into every export.
    record.insert(
        "n_observable_mutations".into(),
        json!(crate::airr_record::builder::observable_mutations(
            sim,
            &mut shm_written
        )),
    );
    record.insert("mutation_rate".into(), json!(sim.mutation_burden));
    if !unknown.is_empty() {
        for field in [
            "n_mutations",
            "mutation_rate",
            "n_pcr_errors",
            "n_quality_errors",
            "n_indels",
            "n_v_indels",
            "n_d_indels",
            "n_j_indels",
            "n_v_mutations",
            "n_d_mutations",
            "n_j_mutations",
            "n_np_mutations",
            "end_loss_5_length",
            "end_loss_3_length",
        ] {
            record.insert(field.into(), Value::Null);
            availability.insert(field.into(), json!("unattributed_observation_history"));
        }
        record["duplex_repair"]["applied_error_substitutions"] =
            json!({"pcr":Value::Null,"quality":Value::Null});
        record["duplex_repair"]["observation_counter_reasons"] = json!(unknown);
    }
}

#[cfg(test)]
mod tests;
