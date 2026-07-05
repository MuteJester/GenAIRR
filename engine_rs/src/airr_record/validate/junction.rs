//! C3: Junction truth.

use super::*;

use super::super::junction::{anchor_amino_acid_preserved, anchor_pool_position, junction_has_stop};
use super::super::projection::lookup_allele;
use crate::ir::Segment;

pub(super) fn check_junction(
    record: &AirrRecord,
    sim: &crate::ir::Simulation,
    refdata: &RefDataConfig,
    issues: &mut Vec<RecordValidationIssue>,
) {
    // Junction recomputation only meaningful for non-rev-comp records;
    // rev-comp flips coordinates and re-translates AA after the
    // junction is sliced. Skip when record is reverse-complemented —
    // the §5 audit covers that path with dedicated tests.
    if record.rev_comp {
        return;
    }

    // Re-derive the junction window the same way builder.rs does:
    // germline_pos scan over V/J regions, not offset arithmetic.
    // (The contract-side `crate::junction::compute_junction` uses
    // offsets and diverges under indels/end-loss; the AIRR builder
    // uses the scan, so the validator must match the builder.)
    let v_region = sim
        .sequence
        .regions
        .iter()
        .find(|r| r.segment == Segment::V);
    let j_region = sim
        .sequence
        .regions
        .iter()
        .find(|r| r.segment == Segment::J);
    let v_id = sim.assignments.get(Segment::V).map(|i| i.allele_id);
    let j_id = sim.assignments.get(Segment::J).map(|i| i.allele_id);
    let v_anchor = lookup_allele(refdata, Segment::V, v_id).and_then(|a| a.anchor);
    let j_anchor = lookup_allele(refdata, Segment::J, j_id).and_then(|a| a.anchor);

    let (Some(vr), Some(jr), Some(va), Some(ja)) = (v_region, j_region, v_anchor, j_anchor) else {
        // Junction undefined; builder leaves fields as defaults.
        return;
    };
    let (Some(vap), Some(jap)) = (
        anchor_pool_position(sim, vr, va as u32),
        anchor_pool_position(sim, jr, ja as u32),
    ) else {
        return;
    };
    let v_anchor_in_pool = vap as i64;
    let j_anchor_in_pool = jap as i64;
    if j_anchor_in_pool + 3 <= v_anchor_in_pool {
        return;
    }
    let jstart = v_anchor_in_pool;
    let jend = j_anchor_in_pool + 3;
    let seq_len = record.sequence_length;
    let safe_start = jstart.clamp(0, seq_len) as usize;
    let safe_end = jend.clamp(0, seq_len) as usize;
    let recomputed_content: String = if safe_end > safe_start {
        record.sequence[safe_start..safe_end].to_string()
    } else {
        String::new()
    };
    let recomputed_len = recomputed_content.len() as u32;

    let reported_len = record.junction_length;
    if reported_len != Some(recomputed_len as i64) {
        issues.push(RecordValidationIssue::JunctionLengthMismatch {
            reported: reported_len,
            recomputed: recomputed_len,
        });
    }

    if !record.junction.eq_ignore_ascii_case(&recomputed_content) {
        issues.push(RecordValidationIssue::JunctionContentMismatch {
            reported: record.junction.clone(),
            recomputed: recomputed_content.clone(),
        });
    }

    // Frame.
    let recomputed_in_frame = recomputed_len % 3 == 0;
    if record.vj_in_frame != Some(recomputed_in_frame) {
        issues.push(RecordValidationIssue::VjInFrameMismatch {
            reported: record.vj_in_frame,
            recomputed: recomputed_in_frame,
        });
    }

    // Stop codon (only meaningful when in-frame).
    let recomputed_stop = recomputed_in_frame && junction_has_stop(&record.junction);
    if record.stop_codon != Some(recomputed_stop) {
        issues.push(RecordValidationIssue::StopCodonMismatch {
            reported: record.stop_codon,
            recomputed: recomputed_stop,
        });
    }

    // Productive triad: in-frame ∧ no stop ∧ V/J anchor amino acids preserved.
    let v_region = sim
        .sequence
        .regions
        .iter()
        .find(|r| r.segment == Segment::V);
    let j_region = sim
        .sequence
        .regions
        .iter()
        .find(|r| r.segment == Segment::J);
    let v_anchor_ok = v_region
        .map(|r| {
            anchor_amino_acid_preserved(
                sim,
                refdata,
                Segment::V,
                r,
                sim.assignments.get(Segment::V).map(|i| i.allele_id),
                record.v_trim_5,
            )
        })
        .unwrap_or(true);
    let j_anchor_ok = j_region
        .map(|r| {
            anchor_amino_acid_preserved(
                sim,
                refdata,
                Segment::J,
                r,
                sim.assignments.get(Segment::J).map(|i| i.allele_id),
                record.j_trim_5,
            )
        })
        .unwrap_or(true);

    let (recomputed_productive, reason) = if !recomputed_in_frame {
        (false, ProductiveDecidedBy::OutOfFrame)
    } else if recomputed_stop {
        (false, ProductiveDecidedBy::JunctionStopCodon)
    } else if !v_anchor_ok {
        (false, ProductiveDecidedBy::VAnchorAaChanged)
    } else if !j_anchor_ok {
        (false, ProductiveDecidedBy::JAnchorAaChanged)
    } else {
        (true, ProductiveDecidedBy::InFrameAndAnchorsPreserved)
    };
    if record.productive != Some(recomputed_productive) {
        issues.push(RecordValidationIssue::ProductiveMismatch {
            reported: record.productive,
            recomputed: recomputed_productive,
            reason,
        });
    }
}
