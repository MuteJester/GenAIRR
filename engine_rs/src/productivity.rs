//! Versioned sequence-level assessment of the observed coding span.
//!
//! Catalogue admission is deliberately separate: accepting a catalogue cannot
//! make an invalid landmark usable. This policy predicts sequence acceptance,
//! not expression, splicing, folding, or antigen binding.
use crate::assignment::AlleleInstance;
use crate::codon::translate_codon;
use crate::ir::{NucHandle, Nucleotide, Segment, Simulation};
use crate::refdata::{Allele, AlleleId, AnchorEvidenceStatus, RefDataConfig};

pub const POLICY: &str = "canonical_observed_vdj_v1";
pub const DUPLEX_POLICY: &str = "duplex_observed_vdj_v1";
pub const SCOPE: &str = "observed_coding_span";

#[derive(Clone, Copy, Debug, Default)]
pub struct ProductivityPolicy;

impl ProductivityPolicy {
    pub fn accepts_residue(self, segment: Segment, residue: u8) -> bool {
        match segment {
            Segment::V => residue == b'C',
            Segment::J => matches!(residue, b'F' | b'W'),
            _ => false,
        }
    }
}

/// Read-only adapter for ordinary simulations and temporary structural edits.
pub trait SequenceView {
    fn len(&self) -> usize;
    fn nucleotide(&self, position: usize) -> Option<Nucleotide>;
    fn region(&self, segment: Segment) -> Option<(usize, usize)>;
    fn assignment(&self, segment: Segment) -> Option<AlleleInstance>;
    /// The landmark-placement policy carried by this observed molecule.
    /// Temporary adapters must preserve this value or decline the fast path.
    fn productivity_policy(&self) -> &'static str {
        POLICY
    }
    /// Certified reference placement, independent of the physical identity
    /// of the query nucleotide. Legacy molecules retain their existing
    /// nucleotide provenance rule; duplex molecules consult their sidecar.
    fn body_reference_position(&self, segment: Segment, position: usize) -> Option<u16> {
        let n = self.nucleotide(position)?;
        (n.segment == segment)
            .then(|| n.germline_pos.get())
            .flatten()
    }
    fn germline_position(&self, segment: Segment, germline: u16) -> Option<usize> {
        let (start, end) = self.region(segment)?;
        (start..end).find(|&p| self.body_reference_position(segment, p) == Some(germline))
    }
    /// The observed position of the first surviving base of the anchor codon at germline
    /// `anchor..anchor + 3`, read from each base's recorded germline position.
    fn anchor_codon_start(&self, segment: Segment, anchor: u16) -> Option<usize> {
        (0..3u16).find_map(|offset| self.germline_position(segment, anchor.checked_add(offset)?))
    }
    /// One past the observed position of the last surviving base of the anchor codon at germline
    /// `anchor..anchor + 3`, read from each base's recorded germline position: a base inserted
    /// inside the codon lengthens the span and a deleted one shortens it.
    fn anchor_codon_end(&self, segment: Segment, anchor: u16) -> Option<usize> {
        (0..3u16)
            .rev()
            .find_map(|offset| self.germline_position(segment, anchor.checked_add(offset)?))
            .map(|p| p + 1)
    }
    fn last_germline(&self, segment: Segment) -> Option<(usize, u16)> {
        let (start, end) = self.region(segment)?;
        (start..end)
            .rev()
            .find_map(|p| self.body_reference_position(segment, p).map(|g| (p, g)))
    }
    /// Reference endpoint used to assess the observed J tail. Legacy
    /// molecules retain their historical query-loss accounting.
    fn j_reference_end_for_phase(&self, reference_len: usize) -> Option<usize> {
        let inst = self.assignment(Segment::J)?;
        reference_len.checked_sub(usize::from(inst.trim_3) + inst.observation_loss_3 as usize)
    }
    /// Whether any complete codon is definitely a stop. Unknown means a
    /// possible stop or no complete codon, not evidence of a stop-free span.
    fn coding_stop(&self, start: usize, end: usize) -> Option<bool> {
        if end.min(self.len()) < start.saturating_add(3) {
            return None;
        }
        let mut unknown = false;
        for p in (start..end.min(self.len()).saturating_sub(2)).step_by(3) {
            let codon =
                [p, p + 1, p + 2].map(|i| self.nucleotide(i).map(|n| n.base).unwrap_or(b'N'));
            match codon_stop(codon) {
                Some(true) => return Some(true),
                None => unknown = true,
                Some(false) => {}
            }
        }
        (!unknown).then_some(false)
    }
}

impl SequenceView for Simulation {
    fn len(&self) -> usize {
        self.pool.len() as usize
    }
    fn nucleotide(&self, position: usize) -> Option<Nucleotide> {
        self.pool.get(NucHandle::new(position as u32)).copied()
    }
    fn region(&self, segment: Segment) -> Option<(usize, usize)> {
        self.sequence
            .regions
            .iter()
            .find(|r| r.segment == segment)
            .map(|r| (r.start.index() as usize, r.end.index() as usize))
    }
    fn assignment(&self, segment: Segment) -> Option<AlleleInstance> {
        self.assignments.get(segment).copied()
    }
    fn productivity_policy(&self) -> &'static str {
        if self.duplex_provenance.is_some() {
            DUPLEX_POLICY
        } else {
            POLICY
        }
    }
    fn body_reference_position(&self, segment: Segment, position: usize) -> Option<u16> {
        if let Some(duplex) = &self.duplex_provenance {
            return duplex
                .body_reference_position(segment, position)
                .and_then(|g| u16::try_from(g).ok());
        }
        let n = self.nucleotide(position)?;
        (n.segment == segment)
            .then(|| n.germline_pos.get())
            .flatten()
    }
    fn j_reference_end_for_phase(&self, reference_len: usize) -> Option<usize> {
        let inst = self.assignment(Segment::J)?;
        if self.duplex_provenance.is_some() {
            // Duplex observation coverage is validated against the removal
            // ledger. Clipping an inserted suffix consumes no reference;
            // clipping past an earlier deletion can consume several reference
            // positions. Neither operation is described by a query-base count.
            let retained_end = reference_len.checked_sub(usize::from(inst.trim_3))?;
            let observed_end = inst
                .observed_ref_end
                .map(|p| p as usize)
                .unwrap_or(retained_end);
            return (observed_end <= retained_end).then_some(observed_end);
        }
        reference_len.checked_sub(usize::from(inst.trim_3) + inst.observation_loss_3 as usize)
    }
}

pub fn policy_name(view: &impl SequenceView) -> &'static str {
    view.productivity_policy()
}

/// IUPAC code -> the set of A/C/G/T it admits, as a 4-bit mask, for every
/// byte value in both cases.
///
/// Built once at compile time instead of matching per call. `base_mask` is the
/// hottest function on the `productive_only()` path -- 12.1% of samples in a
/// native profile of the indel-plus-productive workload, 7.6% of them in the
/// function itself -- because it runs per base inside the productivity checks
/// AND inside the feasibility solver's inner loop, through `codon_stop`. Most
/// of that was the `to_ascii_uppercase()` call: its `is_ascii_lowercase` test
/// alone accounted for 4.7% of all samples. Folding the case into the table
/// removes both the branch chain and the conversion.
const BASE_MASK: [u8; 256] = {
    let mut table = [0u8; 256];
    let codes: [(u8, u8); 16] = [
        (b'A', 1),
        (b'C', 2),
        (b'G', 4),
        (b'T', 8),
        (b'U', 8),
        (b'R', 5),
        (b'Y', 10),
        (b'S', 6),
        (b'W', 9),
        (b'K', 12),
        (b'M', 3),
        (b'B', 14),
        (b'D', 13),
        (b'H', 11),
        (b'V', 7),
        (b'N', 15),
    ];
    let mut i = 0;
    while i < codes.len() {
        let (code, mask) = codes[i];
        table[code as usize] = mask;
        // Case is a separate channel in this engine: a templated base is
        // stored lowercase and an SHM write emits uppercase, so both spellings
        // must map to the same mask.
        table[code.to_ascii_lowercase() as usize] = mask;
        i += 1;
    }
    table
};

fn base_mask(base: u8) -> u8 {
    BASE_MASK[base as usize]
}

/// Amino-acid possibility bitset (ASCII bit positions), without interpreting
/// an unresolved translation X as evidence that a stop is absent.
pub fn codon_amino_acids(codon: [u8; 3]) -> Option<u128> {
    let masks = codon.map(base_mask);
    if masks.contains(&0) {
        return None;
    }
    let bases = *b"ACGT";
    let mut result = 0;
    for a in 0..4 {
        if masks[0] & (1 << a) == 0 {
            continue;
        }
        for b in 0..4 {
            if masks[1] & (1 << b) == 0 {
                continue;
            }
            for c in 0..4 {
                if masks[2] & (1 << c) != 0 {
                    result |= 1u128 << translate_codon(bases[a], bases[b], bases[c]);
                }
            }
        }
    }
    Some(result)
}

pub fn codon_stop(codon: [u8; 3]) -> Option<bool> {
    let [first, second, third] = codon.map(base_mask);
    if first == 0 || second == 0 || third == 0 {
        return None;
    }
    // Standard-code stops are exactly TAA, TAG and TGA. Count compatible
    // stop realizations directly; this is also the hot codon transition in
    // finite junction-completion search.
    let stops = if first & 8 == 0 {
        0
    } else {
        u32::from(second & 1 != 0) * (third & 5).count_ones()
            + u32::from(second & 4 != 0) * u32::from(third & 1 != 0)
    };
    if stops == 0 {
        Some(false)
    } else if stops == first.count_ones() * second.count_ones() * third.count_ones() {
        Some(true)
    } else {
        None
    }
}

pub fn anchor_class(segment: Segment, codon: [u8; 3]) -> Option<bool> {
    let aas = codon_amino_acids(codon)?;
    let allowed = match segment {
        Segment::V => 1u128 << b'C',
        Segment::J => (1u128 << b'F') | (1u128 << b'W'),
        _ => return None,
    };
    if aas & !allowed == 0 {
        Some(true)
    } else if aas & allowed == 0 {
        Some(false)
    } else {
        None
    }
}

/// Compatibility readiness for legacy reference annotations. A canonical
/// residue is necessary, not proof that the declared position was curated.
pub fn reference_anchor_ready(allele: &Allele) -> bool {
    let Some(a) = allele.anchor.map(usize::from) else {
        return false;
    };
    let Some(codon) = allele.seq.get(a..a + 3) else {
        return false;
    };
    anchor_class(allele.segment, [codon[0], codon[1], codon[2]]) == Some(true)
}

pub fn reference_phase(allele: &Allele) -> Option<usize> {
    reference_anchor_ready(allele).then(|| usize::from(allele.anchor.unwrap()) % 3)
}

pub fn reference_ready(refdata: &RefDataConfig, segment: Segment, id: AlleleId) -> bool {
    let Some(allele) = refdata.get(segment, id) else {
        return false;
    };
    if let Some(evidence) = refdata.allele_evidence(segment, id) {
        // Verified positional evidence establishes the landmark even when
        // its original residue is biologically rejected. Keep that rejection
        // distinct from an unresolved annotation. It does not establish a
        // coding phase unless one is supplied independently or the source
        // codon itself is canonical.
        if evidence.anchor_status == AnchorEvidenceStatus::Verified {
            return allele
                .anchor
                .is_some_and(|a| usize::from(a) + 3 <= allele.seq.len());
        }
        if !matches!(evidence.anchor_status, AnchorEvidenceStatus::Inferred) {
            return false;
        }
    }
    reference_anchor_ready(allele)
}

/// Evidence is sequence-bound and validated at installation/compilation.
/// Never substitute an invalid landmark for independent coding-phase evidence.
pub fn reference_phase_for(
    refdata: &RefDataConfig,
    segment: Segment,
    id: AlleleId,
) -> Option<usize> {
    let allele = refdata.get(segment, id)?;
    if let Some(evidence) = refdata.allele_evidence(segment, id) {
        if let Some(phase) = evidence.coding_phase {
            return (phase < 3).then_some(phase as usize);
        }
        if !reference_ready(refdata, segment, id) {
            return None;
        }
    }
    reference_phase(allele)
}

pub fn retained_stop_free(
    refdata: &RefDataConfig,
    segment: Segment,
    id: AlleleId,
    trim_5: u32,
    trim_3: u32,
) -> bool {
    if !reference_ready(refdata, segment, id) {
        return false;
    }
    let Some(allele) = refdata.get(segment, id) else {
        return false;
    };
    // Productive generation additionally requires an accepted source
    // landmark; metadata readiness alone can include a verified
    // noncanonical source landmark.
    if !reference_anchor_ready(allele) {
        return false;
    }
    let Some(phase) = reference_phase_for(refdata, segment, id) else {
        return false;
    };
    let Some(end) = allele.seq.len().checked_sub(trim_3 as usize) else {
        return false;
    };
    let begin = trim_5 as usize;
    if begin > end {
        return false;
    }
    let start = begin + (phase + 3 - begin % 3) % 3;
    (start..end.saturating_sub(2))
        .step_by(3)
        .all(|p| codon_stop([allele.seq[p], allele.seq[p + 1], allele.seq[p + 2]]) == Some(false))
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CodingSpan {
    pub phase: usize,
    pub start: usize,
    pub end: usize,
}

pub fn coding_span(view: &impl SequenceView, refdata: &RefDataConfig) -> Option<CodingSpan> {
    let inst = view.assignment(Segment::V)?;
    let native_phase = reference_phase_for(refdata, Segment::V, inst.allele_id)?;
    let (v_start, v_end) = view.region(Segment::V)?;
    if v_start >= v_end {
        return None;
    }
    let removed = usize::from(inst.trim_5) + inst.observation_loss_5 as usize;
    let phase = (v_start + native_phase + 3 - removed % 3) % 3;
    let start = v_start + (phase + 3 - v_start % 3) % 3;
    let end = view.region(Segment::J).map(|(_, end)| end).unwrap_or(v_end);
    Some(CodingSpan { phase, start, end })
}

#[derive(Clone, Debug)]
pub struct AnchorAssessment {
    pub accepted: Option<bool>,
    pub position: Option<usize>,
    pub reason: &'static str,
}

pub fn anchor(
    view: &impl SequenceView,
    refdata: &RefDataConfig,
    segment: Segment,
) -> AnchorAssessment {
    let unknown = |reason| AnchorAssessment {
        accepted: None,
        position: None,
        reason,
    };
    let Some(inst) = view.assignment(segment) else {
        return unknown("assignment_missing");
    };
    let Some(allele) = refdata.get(segment, inst.allele_id) else {
        return unknown("reference_missing");
    };
    if !reference_ready(refdata, segment, inst.allele_id) {
        return unknown("reference_anchor_unresolved");
    }
    let a = usize::from(allele.anchor.unwrap());
    let Some((start, end)) = view.region(segment) else {
        return unknown("segment_unassembled");
    };
    let retained_end = allele.seq.len().saturating_sub(usize::from(inst.trim_3));
    if usize::from(inst.trim_5) > a || a + 3 > retained_end {
        return AnchorAssessment {
            accepted: Some(false),
            position: None,
            reason: "anchor_removed",
        };
    }
    let observed_start = inst.observed_ref_start.unwrap_or(u32::from(inst.trim_5)) as usize;
    let observed_end = inst.observed_ref_end.unwrap_or(retained_end as u32) as usize;
    if a < observed_start || a + 3 > observed_end || start == end {
        return unknown("anchor_not_observed");
    }
    let position = view.germline_position(segment, a as u16);
    let Some(p) = position else {
        return AnchorAssessment {
            accepted: Some(false),
            position: None,
            reason: "anchor_removed",
        };
    };
    let mut codon = [b'N'; 3];
    for (offset, byte) in codon.iter_mut().enumerate() {
        let Some(n) = view.nucleotide(p + offset) else {
            return unknown("anchor_not_observed");
        };
        if p + offset >= end
            || view
                .body_reference_position(segment, p + offset)
                .map(usize::from)
                != Some(a + offset)
        {
            return AnchorAssessment {
                accepted: Some(false),
                position: Some(p),
                reason: "anchor_provenance_disrupted",
            };
        }
        *byte = n.base;
    }
    let accepted = anchor_class(segment, codon);
    AnchorAssessment {
        accepted,
        position: Some(p),
        reason: match accepted {
            Some(true) => "anchor_accepted",
            Some(false) => "anchor_residue_rejected",
            None => "anchor_ambiguous",
        },
    }
}

/// The junction window, from each base's recorded germline position: from the first surviving
/// base of the V Cys codon to one past the last surviving base of the J W/F codon.  Exonuclease
/// trimming eats a codon only from its inner side, so the outer limits (the V codon's first base,
/// the J codon's last) stay put and N/P bases filling trimmed codon positions fall inside; an indel
/// that deletes an outer base moves the limit to the next surviving codon base, and inserted bases
/// inside a codon widen the window.  It is undefined when either codon is lost entirely, when a
/// codon's outer base is trimmed from the segment's outer end, or when a codon base that
/// recombination kept lies outside the observed read.
pub(crate) fn junction_bounds(
    view: &impl SequenceView,
    refdata: &RefDataConfig,
) -> Option<(usize, usize)> {
    let start = codon_span(view, refdata, Segment::V)?.0;
    let end = codon_span(view, refdata, Segment::J)?.1;
    (end > start && end <= view.len()).then_some((start, end))
}

/// The observed span of whatever survives of the segment's anchor codon: its first surviving base
/// and one past its last.
fn codon_span(
    view: &impl SequenceView,
    refdata: &RefDataConfig,
    segment: Segment,
) -> Option<(usize, usize)> {
    let inst = view.assignment(segment)?;
    let allele = refdata.get(segment, inst.allele_id)?;
    let anchor = allele.anchor?;
    let a = usize::from(anchor);
    let kept = usize::from(inst.trim_5)..allele.seq.len().saturating_sub(usize::from(inst.trim_3));
    // The codon's outer base marks the junction's limit; trimming from the segment's outer end
    // past it leaves no limit in the molecule.
    let outer = if segment == Segment::J { a + 2 } else { a };
    if !kept.contains(&outer) {
        return None;
    }
    // A kept codon base outside the read is unknown, not lost.
    let observed = inst.observed_ref_start.unwrap_or(kept.start as u32) as usize
        ..inst.observed_ref_end.unwrap_or(kept.end as u32) as usize;
    if (a..a + 3).any(|g| kept.contains(&g) && !observed.contains(&g)) {
        return None;
    }
    Some((
        view.anchor_codon_start(segment, anchor)?,
        view.anchor_codon_end(segment, anchor)?,
    ))
}

pub fn j_tail(
    view: &impl SequenceView,
    refdata: &RefDataConfig,
    frame: CodingSpan,
) -> Option<bool> {
    let inst = view.assignment(Segment::J)?;
    let allele = refdata.get(Segment::J, inst.allele_id)?;
    let native_phase = reference_phase_for(refdata, Segment::J, inst.allele_id)?;
    let (start, end) = view.region(Segment::J)?;
    if start == end {
        return None;
    }
    let reference_end = view.j_reference_end_for_phase(allele.seq.len())?;
    if (reference_end + 3 - native_phase) % 3 != (end + 3 - frame.phase) % 3 {
        return Some(false);
    }
    let (p, g) = view.last_germline(Segment::J)?;
    let g = usize::from(g);
    Some((g + 3 - native_phase) % 3 == (p + 3 - frame.phase) % 3 && (end - p - 1) % 3 == 0)
}

#[derive(Clone, Debug)]
pub struct ProductivityAssessment {
    pub productive: Option<bool>,
    pub vj_in_frame: Option<bool>,
    pub stop_codon: Option<bool>,
    pub junction: Option<(usize, usize)>,
    pub frame: Option<CodingSpan>,
    /// The fifth predicate of [`Self::productive`]: whether the J tail joins
    /// the constant region in frame. Carried out of the assessment because a
    /// consumer cannot derive it -- it is germline-phase arithmetic over the
    /// J's reference end and its last germline position, and nothing in the
    /// projected record comes close. Without it, a record that satisfies the
    /// four documented conditions and is reported unproductive looks like a
    /// defect rather than a fifth condition.
    pub j_tail_in_frame: Option<bool>,
    /// The two anchor predicates, carried out for the same reason as the tail:
    /// their tri-state depends on provenance a reader does not have.
    pub v_anchor_accepted: Option<bool>,
    pub j_anchor_accepted: Option<bool>,
    pub reasons: Vec<&'static str>,
}

pub fn assess(view: &impl SequenceView, refdata: &RefDataConfig) -> ProductivityAssessment {
    let v = anchor(view, refdata, Segment::V);
    let j = anchor(view, refdata, Segment::J);
    let frame = coding_span(view, refdata);
    let stop_codon = frame.and_then(|f| view.coding_stop(f.start, f.end));
    // Landmarks only bound the junction when the reference anchors are trusted.
    let trusted = |segment| {
        view.assignment(segment)
            .is_some_and(|inst| reference_ready(refdata, segment, inst.allele_id))
    };
    let junction = (trusted(Segment::V) && trusted(Segment::J))
        .then(|| junction_bounds(view, refdata))
        .flatten();
    let vj_in_frame = junction
        .zip(frame)
        .map(|((start, end), f)| (end - start) % 3 == 0 && start % 3 == f.phase);
    let tail = frame.and_then(|f| j_tail(view, refdata, f));
    let predicates = [
        v.accepted,
        j.accepted,
        vj_in_frame,
        stop_codon.map(|s| !s),
        tail,
    ];
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
    match vj_in_frame {
        Some(false) => reasons.push("vj_out_of_frame"),
        None => reasons.push("vj_frame_unresolved"),
        _ => {}
    }
    match stop_codon {
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
    ProductivityAssessment {
        productive,
        vj_in_frame,
        stop_codon,
        junction,
        frame,
        j_tail_in_frame: tail,
        v_anchor_accepted: v.accepted,
        j_anchor_accepted: j.accepted,
        reasons,
    }
}

#[cfg(test)]
mod duplex_tests;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::airr_record::{build_airr_record, validate_airr_record, RecordValidationIssue};
    use crate::contract::{productive, AnchorPreserved, Contract};
    use crate::ir::Region;
    use crate::pass::Outcome;
    use crate::refdata::ChainType;
    use crate::trace::Trace;

    fn fixture(v: &[u8], j: &[u8]) -> (RefDataConfig, Simulation) {
        let mut cfg = RefDataConfig::empty(ChainType::Vj);
        let mut sim = Simulation::new();
        for (segment, seq, anchor) in [(Segment::V, v, 6), (Segment::J, j, 0)] {
            let allele = Allele {
                name: format!("{segment:?}*01"),
                gene: format!("{segment:?}"),
                seq: seq.to_vec(),
                segment,
                anchor: Some(anchor),
                functional_status: None,
                subregions: vec![],
            };
            let id = if segment == Segment::V {
                cfg.v_pool.push(allele)
            } else {
                cfg.j_pool.push(allele)
            };
            sim = sim.with_allele_assigned(segment, AlleleInstance::new(id));
            let begin = sim.pool.len() as u32;
            for (g, &base) in seq.iter().enumerate() {
                sim = sim
                    .with_nucleotide_pushed(Nucleotide::germline(base, g as u16, segment))
                    .0;
            }
            sim = sim.with_region_added(Region::new(
                segment,
                NucHandle::new(begin),
                NucHandle::new(sim.pool.len() as u32),
            ));
        }
        (cfg, sim)
    }
    fn outcome(sim: Simulation) -> Outcome {
        Outcome {
            revisions: vec![sim],
            pass_names: vec![],
            trace: Trace::new(),
            events: vec![],
        }
    }
    fn assert_record(cfg: &RefDataConfig, sim: Simulation, expected: (Option<bool>, Option<bool>)) {
        let out = outcome(sim);
        let rec = build_airr_record(&out, cfg, "assessment");
        assert_eq!((rec.productive, rec.stop_codon), expected);
        assert!(validate_airr_record(&rec, &out, cfg).is_empty());
    }

    #[test]
    fn stop_certainty_matches_explicit_realizations_for_every_iupac_codon() {
        // Independent expansion table and literal standard-code stops keep
        // the optimized mask arithmetic from validating itself.
        let symbols: &[(u8, &[u8])] = &[
            (b'A', b"A"),
            (b'C', b"C"),
            (b'G', b"G"),
            (b'T', b"T"),
            (b'R', b"AG"),
            (b'Y', b"CT"),
            (b'S', b"CG"),
            (b'W', b"AT"),
            (b'K', b"GT"),
            (b'M', b"AC"),
            (b'B', b"CGT"),
            (b'D', b"AGT"),
            (b'H', b"ACT"),
            (b'V', b"ACG"),
            (b'N', b"ACGT"),
        ];
        for &(a, a_bases) in symbols {
            for &(b, b_bases) in symbols {
                for &(c, c_bases) in symbols {
                    let mut stops = 0;
                    let mut total = 0;
                    for &x in a_bases {
                        for &y in b_bases {
                            for &z in c_bases {
                                total += 1;
                                stops +=
                                    usize::from(matches!(&[x, y, z], b"TAA" | b"TAG" | b"TGA"));
                            }
                        }
                    }
                    let expected = if stops == 0 {
                        Some(false)
                    } else if stops == total {
                        Some(true)
                    } else {
                        None
                    };
                    assert_eq!(codon_stop([a, b, c]), expected, "{a} {b} {c}");
                    assert_eq!(
                        codon_stop([a, b, c].map(|base| base.to_ascii_lowercase())),
                        expected
                    );
                }
            }
        }
        assert_eq!(codon_stop(*b"UAR"), Some(true));
        assert_eq!(codon_stop(*b"TA?"), None);
    }

    #[test]
    fn invalid_reference_is_unknown_not_a_biological_false_or_invented_frame() {
        let (cfg, sim) = fixture(b"AAACCCGGG", b"TTTGGG");
        let a = assess(&sim, &cfg);
        assert_eq!(
            (a.productive, a.stop_codon, a.vj_in_frame, a.junction),
            (None, None, None, None)
        );
        assert!(productive().verify(&sim, Some(&cfg)).is_err());
        assert_record(&cfg, sim, (None, None));
    }

    #[test]
    fn ambiguity_is_not_mistaken_for_absence_of_stops() {
        for (v, expected) in [
            (b"AAANNNTGT".as_slice(), (None, None)),
            (b"AAAGCNTGT".as_slice(), (Some(true), Some(false))),
            (b"AAATARTGT".as_slice(), (Some(false), Some(true))),
            (b"TAANNNTGT".as_slice(), (Some(false), Some(true))),
        ] {
            let (cfg, sim) = fixture(v, b"TTTGGG");
            assert_record(&cfg, sim, expected);
        }
    }

    #[test]
    fn canonical_j_class_allows_f_to_w_while_preservation_remains_optional() {
        let (cfg, sim) = fixture(b"AAACCCTGT", b"TTTGGG");
        // TTT -> TGG at the same three original J nucleotides.
        let sim = sim
            .with_base_changed(NucHandle::new(10), b'G')
            .with_base_changed(NucHandle::new(11), b'G');
        assert_eq!(assess(&sim, &cfg).productive, Some(true));
        assert!(productive().verify(&sim, Some(&cfg)).is_ok());
        assert!(AnchorPreserved::new(Segment::J)
            .verify(&sim, Some(&cfg))
            .is_err());
        assert_record(&cfg, sim, (Some(true), Some(false)));
    }

    #[test]
    fn observed_loss_is_unknown_but_original_anchor_deletion_is_false() {
        let (cfg, sim) = fixture(b"AAACCCTGT", b"TTTGGG");
        let deleted = sim.with_indel_deleted(6);
        assert_eq!(assess(&deleted, &cfg).productive, Some(false));
        let mut lost = sim;
        for _ in 0..9 {
            lost = lost.with_indel_deleted(0);
        }
        let mut inst = lost.assignment(Segment::V).unwrap();
        inst.observation_loss_5 = 9;
        inst.observed_ref_start = Some(9);
        lost = lost.with_allele_assigned(Segment::V, inst);
        assert_record(&cfg, lost, (None, None));
    }

    #[test]
    fn clipping_before_anchor_retains_an_observed_span_prediction() {
        let (cfg, sim) = fixture(b"AAACCCTGT", b"TTTGGG");
        let mut lost = sim.with_indel_deleted(0);
        let mut inst = lost.assignment(Segment::V).unwrap();
        inst.observation_loss_5 = 1;
        inst.observed_ref_start = Some(1);
        lost = lost.with_allele_assigned(Segment::V, inst);
        assert_record(&cfg, lost, (Some(true), Some(false)));
    }

    #[test]
    fn independent_phase_can_prove_a_stop_without_repairing_a_landmark() {
        use crate::refdata::{AlleleEvidence, CodingBoundaryStatus};
        let (mut cfg, sim) = fixture(b"AAATAAGGG", b"TTTGGG");
        let allele = cfg.get(Segment::V, AlleleId::new(0)).unwrap();
        let evidence = AlleleEvidence {
            schema_version: 1,
            sequence_sha256: crate::refdata::evidence::sequence_sha256(&allele.seq),
            coding_phase: Some(0),
            anchor_status: AnchorEvidenceStatus::Invalid,
            boundary_5: CodingBoundaryStatus::Unknown,
            boundary_3: CodingBoundaryStatus::Unknown,
            source: Some("fixture annotation".into()),
            method: "independent_phase".into(),
            notes: vec![],
        };
        cfg.set_allele_evidence(Segment::V, AlleleId::new(0), evidence)
            .unwrap();
        let a = assess(&sim, &cfg);
        assert_eq!(a.productive, Some(false));
        assert_eq!(a.stop_codon, Some(true));
        assert_eq!(a.vj_in_frame, None);
        assert_eq!(a.junction, None);
        assert_record(&cfg, sim, (Some(false), Some(true)));
    }

    #[test]
    fn final_emission_guard_rejects_failed_and_unknown_initial_states_in_both_modes() {
        use crate::compiled::{CompileOptions, CompiledSimulator, ExecutionPolicy};
        use crate::pass::PassPlan;
        for v in [b"AAACCCTGT".as_slice(), b"AAANNNTGT".as_slice()] {
            let (cfg, mut sim) = fixture(v, b"TTTGGG");
            if v == b"AAACCCTGT" {
                sim = sim.with_indel_deleted(6);
            }
            let plan = PassPlan::new();
            let contracts = productive();
            for policy in [ExecutionPolicy::Permissive, ExecutionPolicy::Strict] {
                let compiled = CompiledSimulator::compile_with_options(
                    &plan,
                    Some(&cfg),
                    Some(&contracts),
                    policy,
                    CompileOptions::skip_refdata_validation(),
                )
                .unwrap();
                let error = compiled.run_one_from(sim.clone(), 0).unwrap_err();
                assert_eq!(error.pass_name(), "productive.final");
            }
        }
    }

    #[test]
    fn validator_checks_nullable_claims_and_absent_junction_fields() {
        let (cfg, sim) = fixture(b"AAACCCGGG", b"TTTGGG");
        let out = outcome(sim);
        let mut rec = build_airr_record(&out, &cfg, "tampered");
        rec.productive = Some(true);
        rec.stop_codon = Some(false);
        rec.vj_in_frame = Some(true);
        rec.junction = "GGGTTT".into();
        rec.junction_aa = "GF".into();
        rec.junction_start = Some(6);
        rec.junction_end = Some(12);
        rec.junction_length = Some(6);
        let issues = validate_airr_record(&rec, &out, &cfg);
        for field in [
            "productive",
            "stop_codon",
            "vj_in_frame",
            "junction_start",
            "junction_end",
            "junction_length",
        ] {
            assert!(issues.iter().any(|i| matches!(i,RecordValidationIssue::ProductivityAssessmentMismatch{field:f,..} if f==field)), "{field}: {issues:?}");
        }
        assert!(issues
            .iter()
            .any(|i| matches!(i, RecordValidationIssue::JunctionContentMismatch { .. })));
        assert!(issues
            .iter()
            .any(|i| matches!(i, RecordValidationIssue::JunctionAaMismatch { .. })));
    }

    /// Observed V and J bases: each entry is a base and the germline position it copies, or `None`
    /// for a base an indel inserted.  The V germline is `AAACCCTGT` (Cys codon at germline 6), the J
    /// germline `TGGGGCAAA` (Trp codon at germline 0); `v_inst` and `j_inst` record trims and read
    /// coverage on the assignments.
    fn edited_fixture(
        v: &[(u8, Option<u16>)],
        j: &[(u8, Option<u16>)],
        v_inst: impl FnOnce(AlleleInstance) -> AlleleInstance,
        j_inst: impl FnOnce(AlleleInstance) -> AlleleInstance,
    ) -> (RefDataConfig, Simulation) {
        let (cfg, sim) = fixture(b"AAACCCTGT", b"TGGGGCAAA");
        let mut edited = Simulation::new()
            .with_allele_assigned(
                Segment::V,
                v_inst(sim.assignments.get(Segment::V).copied().unwrap()),
            )
            .with_allele_assigned(
                Segment::J,
                j_inst(sim.assignments.get(Segment::J).copied().unwrap()),
            );
        for (segment, layout) in [(Segment::V, v), (Segment::J, j)] {
            let begin = edited.pool.len() as u32;
            for &(base, germline) in layout {
                let nucleotide = match germline {
                    Some(g) => Nucleotide::germline(base, g, segment),
                    None => Nucleotide::synthetic(base, segment, crate::ir::flag::INDEL_INSERTED),
                };
                edited = edited.with_nucleotide_pushed(nucleotide).0;
            }
            let end = edited.pool.len() as u32;
            edited = edited.with_region_added(Region::new(
                segment,
                NucHandle::new(begin),
                NucHandle::new(end),
            ));
        }
        (cfg, edited)
    }

    /// The intact V (`AAACCCTGT`, germline 0..9) and the J bases `j`.
    fn edited_j_fixture(j: &[(u8, Option<u16>)]) -> (RefDataConfig, Simulation) {
        let v: Vec<_> = b"AAACCCTGT"
            .iter()
            .enumerate()
            .map(|(g, &b)| (b, Some(g as u16)))
            .collect();
        edited_fixture(&v, j, |i| i, |i| i)
    }

    /// The V bases `AAACCCTGT` at the given germline positions.
    fn v_bases(germline: impl IntoIterator<Item = u16>) -> Vec<(u8, Option<u16>)> {
        germline
            .into_iter()
            .map(|g| (b"AAACCCTGT"[g as usize], Some(g)))
            .collect()
    }

    /// The J bases `TGGGGCAAA` at the given germline positions.
    fn j_bases(germline: impl IntoIterator<Item = u16>) -> Vec<(u8, Option<u16>)> {
        germline
            .into_iter()
            .map(|g| (b"TGGGGCAAA"[g as usize], Some(g)))
            .collect()
    }

    /// Whatever survives of each conserved codon bounds the junction: it starts at the first
    /// surviving base of the V Cys codon and ends one past the last surviving base of the J Trp
    /// codon, whether the other bases were lost to exonuclease trimming or to an indel.  Only a
    /// codon lost entirely, or one the read does not cover, leaves the junction undefined.
    /// Expected bounds are written out by hand from the layouts.
    #[test]
    fn junction_spans_whatever_survives_of_each_conserved_codon() {
        type Inst = fn(AlleleInstance) -> AlleleInstance;
        let keep: Inst = |i| i;
        let cases: Vec<(
            &str,
            Vec<(u8, Option<u16>)>,
            Inst,
            Vec<(u8, Option<u16>)>,
            Inst,
            Option<(usize, usize)>,
        )> = vec![
            (
                "V codon's first base deleted",
                v_bases([0, 1, 2, 3, 4, 5, 7, 8]),
                keep,
                j_bases(0..9),
                keep,
                Some((6, 11)),
            ),
            (
                "V 3' trim removes the codon's last base",
                v_bases(0..8),
                |i| i.with_trim_3(1),
                j_bases(0..9),
                keep,
                Some((6, 11)),
            ),
            (
                "V 3' trim removes the codon's last two bases",
                v_bases(0..7),
                |i| i.with_trim_3(2),
                j_bases(0..9),
                keep,
                Some((6, 10)),
            ),
            (
                "V 3' trim removes the whole codon",
                v_bases(0..6),
                |i| i.with_trim_3(3),
                j_bases(0..9),
                keep,
                None,
            ),
            (
                "J 5' trim removes the codon's first base",
                v_bases(0..9),
                keep,
                j_bases(1..9),
                |i| i.with_trim_5(1),
                Some((6, 11)),
            ),
            (
                "J 5' trim removes the codon's first two bases",
                v_bases(0..9),
                keep,
                j_bases(2..9),
                |i| i.with_trim_5(2),
                Some((6, 10)),
            ),
            (
                "J 5' trim removes the whole codon",
                v_bases(0..9),
                keep,
                j_bases(3..9),
                |i| i.with_trim_5(3),
                None,
            ),
            (
                "only the J codon's middle base survives",
                v_bases(0..9),
                keep,
                j_bases([1, 3, 4, 5, 6, 7, 8]),
                keep,
                Some((6, 10)),
            ),
            (
                "V 5' trim removes the codon's first base",
                v_bases([7, 8]),
                |i| i.with_trim_5(7),
                j_bases(0..9),
                keep,
                None,
            ),
            (
                "J 3' trim removes the codon's last base",
                v_bases(0..9),
                keep,
                j_bases([0, 1]),
                |i| i.with_trim_3(7),
                None,
            ),
            (
                "read ends exactly at the J codon's end",
                v_bases(0..9),
                keep,
                j_bases(0..3),
                |i| AlleleInstance {
                    observed_ref_end: Some(3),
                    observation_loss_3: 6,
                    ..i
                },
                Some((6, 12)),
            ),
            (
                "read starts inside the V codon",
                v_bases([7, 8]),
                |i| AlleleInstance {
                    observed_ref_start: Some(7),
                    observation_loss_5: 7,
                    ..i
                },
                j_bases(0..9),
                keep,
                None,
            ),
            (
                "read stops inside the J codon",
                v_bases(0..9),
                keep,
                j_bases(0..2),
                |i| AlleleInstance {
                    observed_ref_end: Some(2),
                    observation_loss_3: 7,
                    ..i
                },
                None,
            ),
        ];
        for (name, v, v_inst, j, j_inst, expected) in cases {
            let (cfg, sim) = edited_fixture(&v, &j, v_inst, j_inst);
            assert_eq!(
                assess(&sim, &cfg).junction,
                expected,
                "{name}: assessed junction"
            );
            let contract = crate::junction::compute_junction(&sim, &cfg)
                .map(|j| (j.start.index() as usize, j.end.index() as usize));
            assert_eq!(contract, expected, "{name}: contract junction");
            let out = outcome(sim);
            let rec = build_airr_record(&out, &cfg, name);
            assert_eq!(
                (rec.junction_start, rec.junction_end),
                (expected.map(|e| e.0 as i64), expected.map(|e| e.1 as i64)),
                "{name}: record junction"
            );
            if let Some((start, end)) = expected {
                assert_eq!(
                    rec.junction,
                    rec.sequence[start..end],
                    "{name}: junction bases"
                );
            }
            let issues = validate_airr_record(&rec, &out, &cfg);
            assert!(issues.is_empty(), "{name}: {issues:?}");
        }
    }

    /// The junction ends one past the observed position of the last surviving base of the J Trp
    /// codon, read from each base's recorded germline position; an indel inside the codon moves
    /// that end by its own length.  Expected ends are written out by hand from the layouts.
    #[test]
    fn junction_ends_after_the_last_surviving_base_of_the_j_codon() {
        let tail = [
            (b'G', Some(3)),
            (b'G', Some(4)),
            (b'C', Some(5)),
            (b'A', Some(6)),
        ];
        let cases: Vec<(&str, Vec<(u8, Option<u16>)>, usize)> = vec![
            (
                "intact codon",
                vec![(b'T', Some(0)), (b'G', Some(1)), (b'G', Some(2))],
                12,
            ),
            (
                "base inserted inside the codon",
                vec![
                    (b'T', Some(0)),
                    (b'A', None),
                    (b'G', Some(1)),
                    (b'G', Some(2)),
                ],
                13,
            ),
            (
                "middle codon base deleted",
                vec![(b'T', Some(0)), (b'G', Some(2))],
                11,
            ),
            (
                "last codon base deleted",
                vec![(b'T', Some(0)), (b'G', Some(1))],
                11,
            ),
            (
                "base inserted between the codon's second and third bases",
                vec![
                    (b'T', Some(0)),
                    (b'G', Some(1)),
                    (b'A', None),
                    (b'G', Some(2)),
                ],
                13,
            ),
            (
                "two bases inserted inside the codon",
                vec![
                    (b'T', Some(0)),
                    (b'A', None),
                    (b'G', Some(1)),
                    (b'C', None),
                    (b'G', Some(2)),
                ],
                14,
            ),
            (
                "second and third codon bases deleted",
                vec![(b'T', Some(0))],
                10,
            ),
            (
                "base inserted just before the codon",
                vec![
                    (b'A', None),
                    (b'T', Some(0)),
                    (b'G', Some(1)),
                    (b'G', Some(2)),
                ],
                13,
            ),
            (
                "base inserted just after the codon",
                vec![
                    (b'T', Some(0)),
                    (b'G', Some(1)),
                    (b'G', Some(2)),
                    (b'A', None),
                ],
                12,
            ),
            (
                "substitution inside the codon",
                vec![(b'T', Some(0)), (b'C', Some(1)), (b'G', Some(2))],
                12,
            ),
        ];
        for (name, codon, end) in cases {
            let mut j = codon;
            j.extend_from_slice(&tail);
            let (cfg, sim) = edited_j_fixture(&j);
            assert_eq!(
                assess(&sim, &cfg).junction,
                Some((6, end)),
                "{name}: assessed junction"
            );
            let contract = crate::junction::compute_junction(&sim, &cfg).expect(name);
            assert_eq!(
                (
                    contract.start.index() as usize,
                    contract.end.index() as usize
                ),
                (6, end),
                "{name}: contract junction"
            );
            let out = outcome(sim);
            let rec = build_airr_record(&out, &cfg, name);
            assert_eq!(
                (rec.junction_start, rec.junction_end),
                (Some(6), Some(end as i64)),
                "{name}"
            );
            assert_eq!(rec.junction, rec.sequence[6..end], "{name}: junction bases");
            assert!(
                validate_airr_record(&rec, &out, &cfg).is_empty(),
                "{name}: record validates"
            );
        }
    }

    /// With every base of a conserved codon deleted, nothing marks the junction's end, so no
    /// junction is defined, in the assessment, in the contract junction and in the validated
    /// record alike.
    #[test]
    fn junction_is_undefined_only_when_a_whole_conserved_codon_is_deleted() {
        let (cfg, sim) = edited_j_fixture(&j_bases(3..9));
        assert_eq!(assess(&sim, &cfg).junction, None);
        assert!(crate::junction::compute_junction(&sim, &cfg).is_none());
        let out = outcome(sim);
        let rec = build_airr_record(&out, &cfg, "whole J codon deleted");
        assert_eq!((rec.junction_start, rec.junction_end), (None, None));
        assert!(rec.junction.is_empty());
        assert!(validate_airr_record(&rec, &out, &cfg).is_empty());
    }

    /// The frame contract classifies an edit by whether it shifts the J anchor's first base, so an
    /// indel already inside the codon leaves that base where its recorded germline position puts
    /// it: the next insertion just after it does not move the anchor.
    #[test]
    fn frame_contract_reads_the_j_anchor_from_its_own_base_after_an_edit_inside_the_codon() {
        use crate::contract::{IndelEventClass, IndelKindHint, ProductiveJunctionFrame};
        let mut j = vec![
            (b'T', Some(0)),
            (b'A', None),
            (b'G', Some(1)),
            (b'G', Some(2)),
        ];
        j.extend_from_slice(&[
            (b'G', Some(3)),
            (b'G', Some(4)),
            (b'C', Some(5)),
            (b'A', Some(6)),
        ]);
        let (cfg, sim) = edited_j_fixture(&j);
        let frame = ProductiveJunctionFrame::new();
        // The J anchor's first base is at pool 9; the junction now ends at 13.
        assert_eq!(
            frame.admissible_indel_class_at(&sim, Some(&cfg), 9, IndelKindHint::Insertion),
            IndelEventClass::FrameDelta(1)
        );
        assert_eq!(
            frame.admissible_indel_class_at(&sim, Some(&cfg), 10, IndelKindHint::Insertion),
            IndelEventClass::FrameNeutral
        );
        assert_eq!(
            frame.admissible_indel_class_at(&sim, Some(&cfg), 8, IndelKindHint::Deletion),
            IndelEventClass::FrameDelta(-1)
        );
        assert_eq!(
            frame.admissible_indel_class_at(&sim, Some(&cfg), 9, IndelKindHint::Deletion),
            IndelEventClass::FrameNeutral
        );
    }

    /// The validator derives the junction end on its own: a record whose end was computed by
    /// adding 3 to the first codon base, when an indel sits inside the codon, is reported.
    #[test]
    fn validator_reports_a_junction_end_that_ignores_an_indel_inside_the_j_codon() {
        let mut j = vec![
            (b'T', Some(0)),
            (b'A', None),
            (b'G', Some(1)),
            (b'G', Some(2)),
        ];
        j.extend_from_slice(&[
            (b'G', Some(3)),
            (b'G', Some(4)),
            (b'C', Some(5)),
            (b'A', Some(6)),
        ]);
        let (cfg, sim) = edited_j_fixture(&j);
        let out = outcome(sim);
        let mut rec = build_airr_record(&out, &cfg, "first-base-plus-three");
        rec.junction_end = Some(12);
        rec.junction_length = Some(6);
        rec.junction = rec.sequence[6..12].to_string();
        rec.junction_aa = crate::codon::translate_seq(&rec.junction);
        let issues = validate_airr_record(&rec, &out, &cfg);
        assert!(issues.iter().any(|i| matches!(i,
            RecordValidationIssue::ProductivityAssessmentMismatch { field, .. } if *field == "junction_end")),
            "{issues:?}");
    }
}

#[cfg(test)]
mod base_mask_tests {
    use super::base_mask;

    /// Characterization over the whole byte range, written so a lookup-table
    /// rewrite of `base_mask` has to reproduce it exactly rather than
    /// approximately. It passes before and after that rewrite on purpose:
    /// there is no new behaviour here, only a faster way to compute the same
    /// mapping, and the only thing worth asserting is that it IS the same.
    ///
    /// The IUPAC rows are spelled out rather than generated from the function
    /// under test, which would make the test a tautology.
    #[test]
    fn base_mask_is_the_iupac_table_in_both_cases_and_zero_elsewhere() {
        let expected: [(u8, u8); 16] = [
            (b'A', 1),
            (b'C', 2),
            (b'G', 4),
            (b'T', 8),
            (b'U', 8),
            (b'R', 5),
            (b'Y', 10),
            (b'S', 6),
            (b'W', 9),
            (b'K', 12),
            (b'M', 3),
            (b'B', 14),
            (b'D', 13),
            (b'H', 11),
            (b'V', 7),
            (b'N', 15),
        ];
        for byte in 0u8..=255 {
            let want = expected
                .iter()
                .find(|(code, _)| *code == byte.to_ascii_uppercase())
                .map(|(_, mask)| *mask)
                .unwrap_or(0);
            assert_eq!(
                base_mask(byte),
                want,
                "base_mask({byte}) = {:?} as a char",
                byte as char
            );
        }
        // The case channel specifically: a templated base is stored lowercase
        // and an SHM write emits uppercase, so every code must agree across
        // case or the productivity checks would read mutated bases differently
        // from germline ones.
        for (code, mask) in expected {
            assert_eq!(base_mask(code.to_ascii_lowercase()), mask);
            assert_eq!(base_mask(code), mask);
        }
    }
}
