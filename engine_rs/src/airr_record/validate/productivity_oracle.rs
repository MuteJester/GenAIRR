//! Deliberately independent, slow assessment for record validation. It does not
//! call the projection or production productivity/frame/anchor helpers.
use crate::ir::{Segment, Simulation};
use crate::refdata::{AnchorEvidenceStatus, RefDataConfig};

pub(super) struct Oracle {
    pub productive: Option<bool>,
    pub frame: Option<bool>,
    pub stop: Option<bool>,
    pub junction: Option<(usize, usize)>,
    pub reasons: Vec<&'static str>,
}

fn resolutions(b: u8) -> &'static [u8] {
    match b.to_ascii_uppercase() {
        b'A' => b"A",
        b'C' => b"C",
        b'G' => b"G",
        b'T' | b'U' => b"T",
        b'R' => b"AG",
        b'Y' => b"CT",
        b'S' => b"CG",
        b'W' => b"AT",
        b'K' => b"GT",
        b'M' => b"AC",
        b'B' => b"CGT",
        b'D' => b"AGT",
        b'H' => b"ACT",
        b'V' => b"ACG",
        b'N' => b"ACGT",
        _ => b"",
    }
}

fn predicate(codon: &[u8], accepts: impl Fn(u8) -> bool) -> Option<bool> {
    if codon.len() != 3 {
        return None;
    }
    let mut yes = false;
    let mut no = false;
    for &a in resolutions(codon[0]) {
        for &b in resolutions(codon[1]) {
            for &c in resolutions(codon[2]) {
                if accepts(crate::codon::translate_codon(a, b, c)) {
                    yes = true;
                } else {
                    no = true;
                }
            }
        }
    }
    match (yes, no) {
        (true, false) => Some(true),
        (false, true) => Some(false),
        _ => None,
    }
}

struct Landmark {
    ready: bool,
    phase: Option<usize>,
    accepted: Option<bool>,
    position: Option<usize>,
    /// The first surviving codon base, and one past the last, found on their own from each
    /// base's placement: this codon's share of the junction window.
    codon_start: Option<usize>,
    codon_end: Option<usize>,
    reason: &'static str,
}

/// Independent reading of certified placements. Do not use the production
/// SequenceView or landmark helpers: a synthetic cache entry alone is not
/// evidence that the founder's reference body was restored.
fn placed_reference(sim: &Simulation, segment: Segment, position: usize) -> Option<usize> {
    if let Some(p) = &sim.duplex_provenance {
        let crate::ir::DuplexQueryLink::FounderColumn(column) = p.links.get(position)? else {
            return None;
        };
        let source = match segment {
            Segment::V => crate::duplex_repair::Substrate::V,
            Segment::D => crate::duplex_repair::Substrate::D,
            Segment::J => crate::duplex_repair::Substrate::J,
            _ => return None,
        };
        let mut claims = p
            .founder
            .body_claims
            .get(*column)?
            .iter()
            .filter(|c| c.segment == source);
        let coordinate = claims.next()?.reference_position;
        return claims
            .all(|c| c.reference_position == coordinate)
            .then_some(coordinate);
    }
    let n = sim.pool.as_slice().get(position)?;
    (n.segment == segment)
        .then(|| n.germline_pos.get())
        .flatten()
        .map(usize::from)
}

fn landmark(sim: &Simulation, refdata: &RefDataConfig, segment: Segment) -> Landmark {
    let mut out = Landmark {
        ready: false,
        phase: None,
        accepted: None,
        position: None,
        codon_start: None,
        codon_end: None,
        reason: "assignment_missing",
    };
    let Some(inst) = sim.assignments.get(segment) else {
        return out;
    };
    out.reason = "reference_missing";
    let Some(allele) = refdata.get(segment, inst.allele_id) else {
        return out;
    };
    let evidence = refdata.allele_evidence(segment, inst.allele_id);
    let evidence_ok = evidence.is_none_or(|e| {
        matches!(
            e.anchor_status,
            AnchorEvidenceStatus::Verified | AnchorEvidenceStatus::Inferred
        )
    });
    let accepts = |aa| match segment {
        Segment::V => aa == b'C',
        Segment::J => aa == b'F' || aa == b'W',
        _ => false,
    };
    let source_codon = allele
        .anchor
        .and_then(|a| allele.seq.get(a as usize..a as usize + 3));
    let canonical_source =
        source_codon.is_some_and(|codon| predicate(codon, accepts) == Some(true));
    // An independently verified coordinate can identify a rejected residue.
    // Legacy or inferred positions require the canonical reference residue.
    let verified_position =
        evidence.is_some_and(|e| e.anchor_status == AnchorEvidenceStatus::Verified);
    out.ready = source_codon.is_some() && (verified_position || (evidence_ok && canonical_source));
    out.phase = evidence
        .and_then(|e| e.coding_phase)
        .filter(|&p| p < 3)
        .map(usize::from)
        .or_else(|| (evidence_ok && canonical_source).then(|| allele.anchor.unwrap() as usize % 3));
    out.reason = "reference_anchor_unresolved";
    if !out.ready {
        return out;
    }
    let a = allele.anchor.unwrap() as usize;
    out.reason = "segment_unassembled";
    let Some(region) = sim.sequence.regions.iter().find(|r| r.segment == segment) else {
        return out;
    };
    let begin = region.start.index() as usize;
    let end = region.end.index() as usize;
    let retained_end = allele.seq.len().saturating_sub(inst.trim_3 as usize);
    // Whatever survives of the codon at its own reference coordinates bounds the junction, as long
    // as the codon keeps its outer base (V: the first, J: the last; only trimming from the
    // segment's outer end removes it) and the read covers every codon base that recombination
    // kept: a trimmed or deleted base is lost, an uncovered one is unknown.
    let kept = inst.trim_5 as usize..retained_end;
    let observed = inst.observed_ref_start.unwrap_or(u32::from(inst.trim_5)) as usize
        ..inst.observed_ref_end.unwrap_or(retained_end as u32) as usize;
    let outer = if segment == Segment::J { a + 2 } else { a };
    if kept.contains(&outer) && (a..a + 3).all(|g| !kept.contains(&g) || observed.contains(&g)) {
        let placed: Vec<usize> = (a..a + 3)
            .filter_map(|g| (begin..end).find(|&q| placed_reference(sim, segment, q) == Some(g)))
            .collect();
        out.codon_start = placed.first().copied();
        out.codon_end = placed.last().map(|q| q + 1);
    }
    if inst.trim_5 as usize > a || a + 3 > retained_end {
        out.accepted = Some(false);
        out.reason = "anchor_removed";
        return out;
    }
    if a < inst.observed_ref_start.unwrap_or(u32::from(inst.trim_5)) as usize
        || a + 3 > inst.observed_ref_end.unwrap_or(retained_end as u32) as usize
        || begin == end
    {
        out.reason = "anchor_not_observed";
        return out;
    }
    let pool = sim.pool.as_slice();
    out.position = (begin..end).find(|&p| placed_reference(sim, segment, p) == Some(a));
    let Some(p) = out.position else {
        out.accepted = Some(false);
        out.reason = "anchor_removed";
        return out;
    };
    let mut codon = Vec::new();
    for offset in 0..3 {
        let Some(n) = pool.get(p + offset) else {
            out.position = None;
            out.reason = "anchor_not_observed";
            return out;
        };
        if p + offset >= end || placed_reference(sim, segment, p + offset) != Some(a + offset) {
            out.accepted = Some(false);
            out.reason = "anchor_provenance_disrupted";
            return out;
        }
        codon.push(n.base);
    }
    out.accepted = predicate(&codon, accepts);
    out.reason = match out.accepted {
        Some(true) => "anchor_accepted",
        Some(false) => "anchor_residue_rejected",
        None => "anchor_ambiguous",
    };
    out
}

pub(super) fn derive(sim: &Simulation, refdata: &RefDataConfig) -> Oracle {
    let v = landmark(sim, refdata, Segment::V);
    let j = landmark(sim, refdata, Segment::J);
    let pool = sim.pool.as_slice();
    let vr = sim
        .sequence
        .regions
        .iter()
        .find(|r| r.segment == Segment::V);
    let jr = sim
        .sequence
        .regions
        .iter()
        .find(|r| r.segment == Segment::J);
    let phase = vr
        .filter(|r| !r.is_empty())
        .zip(v.phase)
        .and_then(|(r, phase)| {
            let inst = sim.assignments.get(Segment::V)?;
            Some(
                (r.start.index() as usize + phase + 3
                    - (inst.trim_5 as usize + inst.observation_loss_5 as usize) % 3)
                    % 3,
            )
        });
    let mut stop = None;
    if let (Some(phase), Some(vr)) = (phase, vr) {
        let begin = vr.start.index() as usize;
        let start = begin + (phase + 3 - begin % 3) % 3;
        let end = jr.unwrap_or(vr).end.index() as usize;
        let mut unknown = false;
        let mut seen = false;
        for p in (start..end.min(pool.len()).saturating_sub(2)).step_by(3) {
            seen = true;
            let codon = [pool[p].base, pool[p + 1].base, pool[p + 2].base];
            match predicate(&codon, |aa| aa == b'*') {
                Some(true) => {
                    stop = Some(true);
                    break;
                }
                None => unknown = true,
                Some(false) => {}
            }
        }
        if stop != Some(true) && seen && !unknown {
            stop = Some(false);
        }
    }
    let junction = v
        .codon_start
        .zip(j.codon_end)
        .filter(|&(start, end)| end > start && end <= pool.len());
    let frame = junction
        .zip(phase)
        .map(|((start, end), phase)| (end - start) % 3 == 0 && start % 3 == phase);
    let tail = (|| {
        let phase = phase?;
        let jp = j.phase?;
        let jr = jr.filter(|r| !r.is_empty())?;
        let inst = sim.assignments.get(Segment::J)?;
        let allele = refdata.get(Segment::J, inst.allele_id)?;
        let reference_end = if sim.duplex_provenance.is_some() {
            let retained = allele.seq.len().checked_sub(inst.trim_3 as usize)?;
            let observed = inst
                .observed_ref_end
                .map(|x| x as usize)
                .unwrap_or(retained);
            (observed <= retained).then_some(observed)?
        } else {
            allele
                .seq
                .len()
                .checked_sub(inst.trim_3 as usize + inst.observation_loss_3 as usize)?
        };
        let end = jr.end.index() as usize;
        if (reference_end + 3 - jp) % 3 != (end + 3 - phase) % 3 {
            return Some(false);
        }
        let (p, g) = (jr.start.index() as usize..end)
            .rev()
            .find_map(|p| placed_reference(sim, Segment::J, p).map(|g| (p, g)))?;
        Some((g + 3 - jp) % 3 == (p + 3 - phase) % 3 && (end - p - 1) % 3 == 0)
    })();
    let predicates = [v.accepted, j.accepted, frame, stop.map(|s| !s), tail];
    let productive = if predicates.contains(&Some(false)) {
        Some(false)
    } else if predicates.contains(&None) {
        None
    } else {
        Some(true)
    };
    let mut reasons = Vec::new();
    if v.accepted != Some(true) {
        reasons.push(v.reason);
    }
    if j.accepted != Some(true) {
        reasons.push(j.reason);
    }
    if phase.is_none() {
        reasons.push("coding_phase_unresolved");
    }
    match frame {
        Some(false) => reasons.push("vj_out_of_frame"),
        None => reasons.push("vj_frame_unresolved"),
        _ => {}
    }
    match stop {
        Some(true) => reasons.push("coding_stop"),
        None => reasons.push("coding_stop_unresolved"),
        _ => {}
    }
    match tail {
        Some(false) => reasons.push("j_tail_out_of_frame"),
        None => reasons.push("j_tail_unresolved"),
        _ => {}
    }
    if productive == Some(true) {
        reasons.push("accepted");
    }
    Oracle {
        productive,
        frame,
        stop,
        junction,
        reasons,
    }
}
