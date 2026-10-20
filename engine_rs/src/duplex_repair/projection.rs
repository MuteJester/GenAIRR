//! Standalone truth projection. Physical ancestry and alignment placement are
//! separate; no aligner, motif search, or post-mutation re-call chooses a path.
use super::bindings::DuplexReferenceBindings;
use super::*;
use crate::codon::translate_codon;
use crate::ir::Segment;
use crate::refdata::RefDataConfig;
use serde::Serialize;
use serde_json::{json, Map, Value};
use sha2::{Digest, Sha256};
use std::collections::BTreeSet;

pub const PROJECTION_MODEL: &str = "duplex_airr_projection_v1";
pub const PRODUCTIVITY_POLICY: &str = crate::productivity::DUPLEX_POLICY;

#[derive(Clone, Debug, Serialize)]
pub(super) struct BodyPoint {
    pub(super) query_position: usize,
    pub(super) reference_position: usize,
    pub(super) support_ids: Vec<u64>,
}

#[derive(Clone, Debug, Serialize)]
pub(super) struct Path {
    pub(super) available: bool,
    pub(super) reason: Option<&'static str>,
    pub(super) body: Vec<BodyPoint>,
    pub(super) query_start: Option<usize>,
    pub(super) query_end: Option<usize>,
    pub(super) reference_start: Option<usize>,
    pub(super) reference_end: Option<usize>,
    pub(super) oriented_reference_start: Option<usize>,
    pub(super) oriented_reference_end: Option<usize>,
    pub(super) cigar: Option<String>,
    pub(super) sequence_alignment: Option<String>,
    pub(super) germline_alignment: Option<String>,
    pub(super) identity: Option<f64>,
    #[serde(rename = "query_boundary", skip_serializing_if = "Option::is_none")]
    pub(super) empty_query_boundary: Option<usize>,
}

fn label(segment: Substrate) -> &'static str {
    match segment {
        Substrate::V => "v",
        Substrate::D => "d",
        Substrate::J => "j",
    }
}
fn segment(segment: Substrate) -> Segment {
    match segment {
        Substrate::V => Segment::V,
        Substrate::D => Segment::D,
        Substrate::J => Segment::J,
    }
}
fn oriented(g: usize, length: usize, inverted: bool) -> usize {
    if inverted {
        length - 1 - g
    } else {
        g
    }
}
fn translated(sequence: &[u8], phase: usize) -> String {
    sequence
        .get(phase..)
        .unwrap_or_default()
        .chunks_exact(3)
        .map(|c| translate_codon(c[0], c[1], c[2]) as char)
        .collect()
}
/// Translate only complete coding-frame codons wholly inside an interval.
pub(super) fn translate_interval(query: &[u8], start: usize, end: usize, phase: usize) -> String {
    let first = start + (phase + 3 - start % 3) % 3;
    if first >= end {
        return String::new();
    }
    translated(&query[first..end], 0)
}

fn push_run(runs: &mut Vec<(usize, char)>, op: char) {
    if let Some((n, previous)) = runs.last_mut() {
        if *previous == op {
            *n += 1;
            return;
        }
    }
    runs.push((1, op));
}

fn build_path(body: Vec<BodyPoint>, reference: &[u8], query: &[u8], inverted: bool) -> Path {
    let mut path = Path {
        available: false,
        reason: None,
        body,
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
    if path.body.is_empty() {
        path.reason = Some("body_absent_after_repair");
        return path;
    }
    let pairs: Vec<_> = path
        .body
        .iter()
        .map(|b| {
            (
                b.query_position,
                oriented(b.reference_position, reference.len(), inverted),
            )
        })
        .collect();
    if pairs
        .windows(2)
        .any(|p| p[0].0 >= p[1].0 || p[0].1 >= p[1].1)
    {
        path.reason = Some("nonunique_or_nonmonotone_body_path");
        return path;
    }
    let (qstart, rstart) = pairs[0];
    let (qend, rend) = pairs[pairs.len() - 1];
    let mut qa = Vec::new();
    let mut ra = Vec::new();
    let mut runs = Vec::new();
    let mut r = rstart;
    let mut next = 0;
    let mut matches = 0usize;
    for (q, &base) in query.iter().enumerate().take(qend + 1).skip(qstart) {
        if next < pairs.len() && pairs[next].0 == q {
            while r < pairs[next].1 {
                let canonical = oriented(r, reference.len(), inverted);
                qa.push(b'-');
                ra.push(if inverted {
                    complement(reference[canonical] as char) as u8
                } else {
                    reference[canonical]
                });
                push_run(&mut runs, 'D');
                r += 1;
            }
            let canonical = oriented(r, reference.len(), inverted);
            let germline = if inverted {
                complement(reference[canonical] as char) as u8
            } else {
                reference[canonical]
            };
            qa.push(base);
            ra.push(germline);
            push_run(&mut runs, 'M');
            matches += usize::from(base.eq_ignore_ascii_case(&germline));
            r += 1;
            next += 1;
        } else {
            qa.push(base);
            ra.push(b'-');
            push_run(&mut runs, 'I');
        }
    }
    path.available = true;
    path.query_start = Some(qstart);
    path.query_end = Some(qend + 1);
    path.reference_start = path.body.iter().map(|b| b.reference_position).min();
    path.reference_end = path.body.iter().map(|b| b.reference_position + 1).max();
    path.oriented_reference_start = Some(rstart);
    path.oriented_reference_end = Some(rend + 1);
    path.cigar = Some(runs.iter().map(|(n, op)| format!("{n}{op}")).collect());
    path.identity = Some(matches as f64 / qa.len() as f64);
    path.sequence_alignment = Some(String::from_utf8(qa).expect("ASCII alignment"));
    path.germline_alignment = Some(String::from_utf8(ra).expect("ASCII alignment"));
    path
}

/// Verify the complete physical/evidence snapshot against deterministic replay.
/// Bounds use actual material and work; the original larger limits do not
/// authorize any operation missing from the recorded history.
fn verify_snapshot(annotated: &AnnotatedCoupledProduct) -> Result<(), String> {
    let p = &annotated.product;
    let operations = p
        .joins
        .iter()
        .try_fold(0usize, |n, j| n.checked_add(j.operations.len()))
        .ok_or("operation count overflow")?;
    let replay = repair_coupled_annotated(
        &p.v_reference,
        &p.d_reference,
        &p.j_reference,
        p.d_inverted,
        &p.joins,
        operations,
        p.ledger.len(),
    )?;
    if &replay != annotated {
        return Err("duplex snapshot or placement evidence disagrees with exact replay".into());
    }
    Ok(())
}

fn set(
    record: &mut Map<String, Value>,
    availability: &mut Map<String, Value>,
    field: &str,
    value: Value,
    reason: &str,
) {
    if value.is_null() {
        availability.insert(field.into(), json!(reason));
    }
    record.insert(field.into(), value);
}

#[derive(Debug)]
struct Landmark {
    position: Option<usize>,
    accepted: Option<bool>,
    reason: &'static str,
}

fn landmark(
    path: &Path,
    substrate: Substrate,
    bindings: &DuplexReferenceBindings,
    query: &[u8],
) -> Landmark {
    let source = bindings.for_substrate(substrate);
    let unknown = |reason| Landmark {
        position: None,
        accepted: None,
        reason,
    };
    if !source.landmark_ready {
        return unknown("reference_anchor_unresolved");
    }
    let Some(anchor) = source.anchor else {
        return unknown("reference_anchor_unresolved");
    };
    let mut positions = Vec::new();
    for g in anchor..anchor + 3 {
        let matches: Vec<_> = path
            .body
            .iter()
            .filter(|b| b.reference_position == g)
            .collect();
        if matches.is_empty() {
            return Landmark {
                position: None,
                accepted: Some(false),
                reason: "anchor_removed",
            };
        }
        if matches.len() != 1 {
            return unknown("anchor_placement_ambiguous");
        }
        let q = matches[0].query_position;
        if path
            .body
            .iter()
            .any(|b| b.query_position == q && b.reference_position != g)
        {
            return unknown("anchor_placement_ambiguous");
        }
        positions.push(q);
    }
    if positions[1] != positions[0] + 1 || positions[2] != positions[0] + 2 {
        return Landmark {
            position: Some(positions[0]),
            accepted: Some(false),
            reason: "anchor_provenance_disrupted",
        };
    }
    let q = positions[0];
    let accepted = crate::productivity::anchor_class(
        segment(substrate),
        [query[q], query[q + 1], query[q + 2]],
    );
    Landmark {
        position: Some(q),
        accepted,
        reason: match accepted {
            Some(true) => "anchor_accepted",
            Some(false) => "anchor_residue_rejected",
            None => "anchor_ambiguous",
        },
    }
}

fn assess(
    paths: &[Path; 3],
    bindings: &DuplexReferenceBindings,
    query: &[u8],
    record: &mut Map<String, Value>,
    availability: &mut Map<String, Value>,
) -> Option<usize> {
    let v = landmark(&paths[0], Substrate::V, bindings, query);
    let j = landmark(&paths[2], Substrate::J, bindings, query);
    let junction = v
        .position
        .zip(j.position)
        .filter(|&(v, j)| j >= v && j + 3 <= query.len())
        .map(|(v, j)| (v, j + 3));
    let frame = bindings
        .v
        .coding_phase
        .zip(paths[0].body.first())
        .map(|(native, b)| (b.query_position + 3 - ((b.reference_position + 3 - native) % 3)) % 3);
    let coding_span = frame
        .zip(paths[0].query_start)
        .zip(paths[2].query_end)
        .map(|((phase, start), end)| (start + (phase + 3 - start % 3) % 3, end));
    let stop = coding_span.and_then(|(start, end)| {
        if start + 3 > end {
            return None;
        }
        let mut unknown = false;
        for c in query[start..end].chunks_exact(3) {
            match crate::productivity::codon_stop([c[0], c[1], c[2]]) {
                Some(true) => return Some(true),
                None => unknown = true,
                Some(false) => (),
            }
        }
        (!unknown).then_some(false)
    });
    let in_frame = junction
        .zip(frame)
        .map(|((start, end), phase)| (end - start) % 3 == 0 && start % 3 == phase);
    let tail = frame
        .zip(bindings.j.coding_phase)
        .zip(paths[2].body.last())
        .zip(paths[2].query_end)
        .zip(paths[2].reference_end)
        .map(|((((phase, native), b), end), reference_end)| {
            (reference_end + 3 - native) % 3 == (end + 3 - phase) % 3
                && (b.reference_position + 3 - native) % 3 == (b.query_position + 3 - phase) % 3
                && (end - b.query_position - 1) % 3 == 0
        });
    let predicates = [v.accepted, j.accepted, in_frame, stop.map(|x| !x), tail];
    let productive = if predicates.contains(&Some(false)) {
        Some(false)
    } else if predicates.iter().all(|p| *p == Some(true)) {
        Some(true)
    } else {
        None
    };
    let mut reasons = Vec::new();
    if v.accepted != Some(true) {
        reasons.push(v.reason);
    }
    if j.accepted != Some(true) {
        reasons.push(j.reason);
    }
    if frame.is_none() {
        reasons.push("coding_phase_unresolved");
    }
    match in_frame {
        Some(false) => reasons.push("vj_out_of_frame"),
        None => reasons.push("vj_frame_unresolved"),
        _ => (),
    }
    match stop {
        Some(true) => reasons.push("coding_stop"),
        None => reasons.push("coding_stop_unresolved"),
        _ => (),
    }
    match tail {
        Some(false) => reasons.push("j_tail_out_of_frame"),
        None => reasons.push("j_tail_unresolved"),
        _ => (),
    }
    if productive == Some(true) {
        reasons.push("accepted");
    }
    for (field, value) in [
        ("productive", json!(productive)),
        ("vj_in_frame", json!(in_frame)),
        ("stop_codon", json!(stop)),
        ("junction_start", json!(junction.map(|j| j.0))),
        ("junction_end", json!(junction.map(|j| j.1))),
        ("junction_length", json!(junction.map(|(a, b)| b - a))),
        (
            "junction",
            json!(junction.map(|(a, b)| String::from_utf8_lossy(&query[a..b]).into_owned())),
        ),
        (
            "junction_aa",
            json!(junction.map(|(a, b)| translated(&query[a..b], 0))),
        ),
    ] {
        set(
            record,
            availability,
            field,
            value,
            "productivity_evidence_unresolved",
        );
    }
    record.insert("productivity_policy".into(), json!(PRODUCTIVITY_POLICY));
    record.insert(
        "productivity_scope".into(),
        json!(crate::productivity::SCOPE),
    );
    record.insert("productivity_reasons".into(), json!(reasons));
    set(
        record,
        availability,
        "sequence_aa",
        json!(frame.map(|p| translated(query, p))),
        "coding_phase_unresolved",
    );
    frame
}

/// Project a complete coupled repair snapshot using exact bound generating
/// alleles. Valid overlapping products produce partial records, not errors.
pub(super) struct ProjectionAnalysis {
    pub(super) paths: [Path; 3],
    pub(super) column_claims: Vec<Value>,
    pub(super) overlaps: Vec<Value>,
    pub(super) layout_reasons: BTreeSet<String>,
}

/// Limits of the current execution IR, separate from whether a complete
/// standalone alignment exists. This is shared with the atomic seed importer.
pub(super) fn seed_capability_reasons(
    paths: &[Path; 3],
    layout_reasons: &BTreeSet<String>,
    bindings: &DuplexReferenceBindings,
    query_len: usize,
) -> BTreeSet<String> {
    let mut reasons = layout_reasons.clone();
    if query_len > u32::MAX as usize {
        reasons.insert("query_length_outside_native_handle_range".into());
    }
    for (source, path) in [Substrate::V, Substrate::D, Substrate::J]
        .into_iter()
        .zip(paths)
    {
        let prefix = label(source);
        if [path.query_start, path.query_end, path.empty_query_boundary]
            .into_iter()
            .flatten()
            .any(|q| q > u32::MAX as usize)
            || path
                .body
                .iter()
                .any(|p| p.query_position > u32::MAX as usize)
        {
            reasons.insert(format!(
                "{prefix}:query_coordinate_outside_native_handle_range"
            ));
        }
        // u16::MAX is reserved by GermlinePos::NONE, not a usable coordinate.
        if path
            .body
            .iter()
            .any(|p| p.reference_position >= u16::MAX as usize)
        {
            reasons.insert(format!(
                "{prefix}:body_reference_coordinate_outside_native_ir"
            ));
        }
        if path.reference_start.is_some_and(|g| g > u16::MAX as usize) {
            reasons.insert(format!("{prefix}:trim_5_outside_native_ir"));
        }
        if let Some(end) = path.reference_end {
            match bindings
                .for_substrate(source)
                .sequence
                .len()
                .checked_sub(end)
            {
                Some(trim) if trim <= u16::MAX as usize => {}
                _ => {
                    reasons.insert(format!("{prefix}:trim_3_outside_native_ir"));
                }
            }
        }
    }
    reasons
}

pub(super) fn analyze_coupled(
    annotated: &AnnotatedCoupledProduct,
    bindings: &DuplexReferenceBindings,
    refdata: &RefDataConfig,
) -> Result<ProjectionAnalysis, String> {
    let product = &annotated.product;
    bindings.validate(product, refdata)?;
    verify_snapshot(annotated)?;
    let query = product.sequence.as_bytes();
    let mut bodies: [Vec<BodyPoint>; 3] = std::array::from_fn(|_| Vec::new());
    let substrates = [Substrate::V, Substrate::D, Substrate::J];
    let mut column_claims = Vec::new();
    let mut overlaps = Vec::new();
    for (q, column) in product.columns.iter().enumerate() {
        let mut claims: BTreeMap<(Substrate, usize), Vec<u64>> = BTreeMap::new();
        for base in [&column.top, &column.bottom] {
            let placement = annotated
                .evidence
                .placements
                .get(&base.id)
                .ok_or("missing physical placement")?;
            if let Some(body) = &placement.body {
                let binding = bindings.for_substrate(body.segment);
                if body.reference_position >= binding.sequence.len() {
                    return Err("body reference coordinate out of bounds".into());
                }
                let offset = oriented(
                    body.reference_position,
                    binding.sequence.len(),
                    body.segment == Substrate::D && product.d_inverted,
                );
                let expected = product.origin_positions[&body.segment]
                    .checked_add(offset as i64)
                    .ok_or("body placement overflow")?;
                if expected != column.position {
                    return Err("body certificate is not in its native column".into());
                }
                claims
                    .entry((body.segment, body.reference_position))
                    .or_default()
                    .push(base.id);
            }
        }
        let mut serialized = Vec::new();
        for ((s, g), mut ids) in claims {
            ids.sort_unstable();
            let index = substrates.iter().position(|&x| x == s).unwrap();
            bodies[index].push(BodyPoint {
                query_position: q,
                reference_position: g,
                support_ids: ids.clone(),
            });
            serialized.push(
                json!({"segment":s,"source_instance":bindings.for_substrate(s).source_instance,
                "reference_position":g,"support_ids":ids}),
            );
        }
        if serialized.len() > 1 {
            overlaps.push(json!({"query_position":q,"claims":serialized}));
        }
        column_claims.push(json!({"query_position":q,"claims":serialized}));
    }
    let paths: [Path; 3] = std::array::from_fn(|i| {
        build_path(
            std::mem::take(&mut bodies[i]),
            bindings
                .for_substrate(substrates[i])
                .sequence
                .to_ascii_uppercase()
                .as_bytes(),
            query,
            i == 1 && product.d_inverted,
        )
    });
    let mut layout_reasons = BTreeSet::new();
    if !overlaps.is_empty() {
        layout_reasons.insert("overlapping_body_claims".to_string());
    }
    for (s, path) in substrates.iter().zip(&paths) {
        if !path.available {
            layout_reasons.insert(format!("{}:{}", label(*s), path.reason.unwrap()));
        }
    }
    for pair in paths.windows(2) {
        if pair[0]
            .query_end
            .zip(pair[1].query_start)
            .is_some_and(|(a, b)| a > b)
        {
            layout_reasons.insert("overlapping_or_reordered_segment_paths".to_string());
        }
    }
    if paths[0].query_start.is_some_and(|p| p != 0)
        || paths[2].query_end.is_some_and(|p| p != query.len())
    {
        layout_reasons.insert("outer_nonbody_query_columns".to_string());
    }
    Ok(ProjectionAnalysis {
        paths,
        column_claims,
        overlaps,
        layout_reasons,
    })
}

pub fn project_coupled(
    annotated: &AnnotatedCoupledProduct,
    bindings: &DuplexReferenceBindings,
    refdata: &RefDataConfig,
) -> Result<Value, String> {
    let analysis = analyze_coupled(annotated, bindings, refdata)?;
    render_projection(
        annotated,
        bindings,
        refdata,
        annotated.product.sequence.as_bytes(),
        analysis,
    )
}

pub(super) fn render_projection(
    annotated: &AnnotatedCoupledProduct,
    bindings: &DuplexReferenceBindings,
    refdata: &RefDataConfig,
    query: &[u8],
    analysis: ProjectionAnalysis,
) -> Result<Value, String> {
    let ProjectionAnalysis {
        paths,
        column_claims,
        overlaps,
        layout_reasons,
    } = analysis;
    let product = &annotated.product;
    let substrates = [Substrate::V, Substrate::D, Substrate::J];
    let representable = layout_reasons.is_empty();
    let seed_reasons = seed_capability_reasons(&paths, &layout_reasons, bindings, query.len());
    let mut record = Map::new();
    let mut availability = Map::new();
    // Replayed history has the same stable identity; distinct generating
    // bindings or histories do not collapse just because DNA matches.
    let identity = serde_json::to_vec(&(PROJECTION_MODEL, bindings, &product.joins))
        .map_err(|e| e.to_string())?;
    record.insert(
        "sequence_id".into(),
        json!(format!("duplex:{:x}", Sha256::digest(identity))),
    );
    record.insert("sequence".into(), json!(String::from_utf8_lossy(query)));
    record.insert("sequence_length".into(), json!(query.len()));
    record.insert("rev_comp".into(), json!(false));
    record.insert("d_inverted".into(), json!(product.d_inverted));
    set(
        &mut record,
        &mut availability,
        "locus",
        json!(refdata.identity.locus),
        "reference_locus_unspecified",
    );
    for (s, path) in substrates.iter().zip(&paths) {
        let l = label(*s);
        let binding = bindings.for_substrate(*s);
        record.insert(format!("{l}_call"), json!(binding.allele_name));
        set(
            &mut record,
            &mut availability,
            &format!("recalled_{l}_call"),
            Value::Null,
            "recall_not_computed",
        );
        for (suffix, value) in [
            ("cigar", json!(path.cigar)),
            ("identity", json!(path.identity)),
            ("sequence_start", json!(path.query_start)),
            ("sequence_end", json!(path.query_end)),
            ("germline_start", json!(path.reference_start)),
            ("germline_end", json!(path.reference_end)),
            ("trim_5", json!(path.reference_start)),
            (
                "trim_3",
                json!(path.reference_end.map(|g| binding.sequence.len() - g)),
            ),
        ] {
            set(
                &mut record,
                &mut availability,
                &format!("{l}_{suffix}"),
                value,
                path.reason.unwrap_or("segment_alignment_unresolved"),
            );
        }
    }
    let mut aligned_query = String::new();
    let mut aligned_reference = String::new();
    let mut d_mask = String::new();
    let mut cursor = 0;
    for (i, path) in paths.iter().enumerate() {
        let l = label(substrates[i]);
        if representable {
            let start = path.query_start.or(path.empty_query_boundary).unwrap();
            let gap = &query[cursor..start];
            aligned_query.push_str(&String::from_utf8_lossy(gap));
            aligned_reference.push_str(&"N".repeat(gap.len()));
            d_mask.push_str(&"N".repeat(gap.len()));
            record.insert(format!("{l}_alignment_start"), json!(aligned_query.len()));
            aligned_query.push_str(path.sequence_alignment.as_ref().unwrap());
            let ra = path.germline_alignment.as_ref().unwrap();
            aligned_reference.push_str(ra);
            if i == 1 {
                d_mask.extend(ra.bytes().map(|b| if b == b'-' { '-' } else { 'N' }));
            } else {
                d_mask.push_str(ra);
            }
            record.insert(format!("{l}_alignment_end"), json!(aligned_query.len()));
            cursor = path.query_end.or(path.empty_query_boundary).unwrap();
        } else {
            for suffix in ["alignment_start", "alignment_end"] {
                set(
                    &mut record,
                    &mut availability,
                    &format!("{l}_{suffix}"),
                    Value::Null,
                    "combined_layout_unrepresentable",
                );
            }
        }
    }
    for (field, value) in [
        ("sequence_alignment", &aligned_query),
        ("germline_alignment", &aligned_reference),
        ("germline_alignment_d_mask", &d_mask),
    ] {
        set(
            &mut record,
            &mut availability,
            field,
            if representable {
                json!(value)
            } else {
                Value::Null
            },
            "combined_layout_unrepresentable",
        );
    }
    let frame = assess(&paths, bindings, query, &mut record, &mut availability);
    for (i, np) in ["np1", "np2"].iter().enumerate() {
        let bounds = paths[i]
            .query_end
            .or(paths[i].empty_query_boundary)
            .zip(
                paths[i + 1]
                    .query_start
                    .or(paths[i + 1].empty_query_boundary),
            )
            .map(|(a, b)| (a, b.max(a)));
        let tract =
            bounds.map(|(a, b)| std::str::from_utf8(&query[a..b]).expect("ASCII coding query"));
        for (field, value) in [
            (np.to_string(), json!(tract)),
            (format!("{np}_length"), json!(tract.map(str::len))),
            (
                format!("{np}_aa"),
                json!(bounds
                    .zip(frame)
                    .map(|((a, b), phase)| translate_interval(query, a, b, phase))),
            ),
        ] {
            set(
                &mut record,
                &mut availability,
                &field,
                value,
                if field.ends_with("_aa") && bounds.is_some() {
                    "coding_phase_unresolved"
                } else {
                    "junction_boundary_unresolved"
                },
            );
        }
    }
    for end in ["v_3", "d_5", "d_3", "j_5"] {
        set(
            &mut record,
            &mut availability,
            &format!("p_{end}_length"),
            Value::Null,
            "legacy_p_emission_not_defined_for_duplex",
        );
    }
    // Repair is not SHM, PCR damage or a later polymerase indel.
    for field in [
        "n_mutations",
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
        record.insert(field.into(), json!(0));
    }
    record.insert("mutation_rate".into(), json!(0.0));
    let mut roles = Map::new();
    for column in &product.columns {
        let role = serde_json::to_value(annotated.evidence.placements[&column.top.id].role)
            .map_err(|e| e.to_string())?;
        let name = role.as_str().ok_or("invalid placement role")?;
        let n = roles.get(name).and_then(Value::as_u64).unwrap_or(0);
        roles.insert(name.into(), json!(n + 1));
    }
    let counts = json!({"physical":{
        "original_nucleotides":2*(product.v_reference.len()+product.d_reference.len()+product.j_reference.len()),
        "resected_nucleotides":product.removed.len(),
        "untemplated_additions":product.ledger.iter().filter(|b| b.kind=="n_addition").count(),
        "templated_synthesis":product.ledger.iter().filter(|b| b.kind=="templated_synthesis").count(),
        "hairpin_transfer_events":annotated.evidence.events.iter().filter(|e| matches!(e,RepairEvidenceEvent::Transferred{..})).count(),
        "retained_nucleotides":2*product.sequence.len()},
        "founder_query":{"timepoint":"after_both_joins","strand":"top","columns":product.sequence.len(),"placement_roles":roles}});
    record.insert("duplex_repair".into(), json!({"model":"duplex_provenance_v1",
        "projection_policy":PROJECTION_MODEL,"coordinate_system":"zero_based_coding_oriented",
        "bindings":bindings,"product":product,"evidence":annotated.evidence,
        "column_body_claims":column_claims,"overlaps":overlaps,"counts":counts,"availability":availability}));
    let mut path_json = Map::new();
    for (s, path) in substrates.iter().zip(&paths) {
        path_json.insert(
            label(*s).into(),
            serde_json::to_value(path).map_err(|e| e.to_string())?,
        );
    }
    Ok(
        json!({"model":PROJECTION_MODEL,"record":record,"paths":path_json,
        "availability":availability,"capabilities":{"simulation_seed":seed_reasons.is_empty(),
            "experiment_integration":false,"reasons":seed_reasons}}),
    )
}

#[cfg(test)]
mod tests;
