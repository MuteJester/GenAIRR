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
    fn germline_position(&self, segment: Segment, germline: u16) -> Option<usize> {
        let (start, end) = self.region(segment)?;
        (start..end).find(|&p| {
            self.nucleotide(p)
                .is_some_and(|n| n.segment == segment && n.germline_pos.get() == Some(germline))
        })
    }
    fn last_germline(&self, segment: Segment) -> Option<(usize, u16)> {
        let (start, end) = self.region(segment)?;
        (start..end).rev().find_map(|p| {
            let n = self.nucleotide(p)?;
            (n.segment == segment)
                .then(|| n.germline_pos.get().map(|g| (p, g)))
                .flatten()
        })
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
}

fn base_mask(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'A' => 1,
        b'C' => 2,
        b'G' => 4,
        b'T' | b'U' => 8,
        b'R' => 5,
        b'Y' => 10,
        b'S' => 6,
        b'W' => 9,
        b'K' => 12,
        b'M' => 3,
        b'B' => 14,
        b'D' => 13,
        b'H' => 11,
        b'V' => 7,
        b'N' => 15,
        _ => 0,
    }
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
            || n.segment != segment
            || n.germline_pos.get().map(usize::from) != Some(a + offset)
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
    let reference_end = allele
        .seq
        .len()
        .checked_sub(usize::from(inst.trim_3) + inst.observation_loss_3 as usize)?;
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
    pub reasons: Vec<&'static str>,
}

pub fn assess(view: &impl SequenceView, refdata: &RefDataConfig) -> ProductivityAssessment {
    let v = anchor(view, refdata, Segment::V);
    let j = anchor(view, refdata, Segment::J);
    let frame = coding_span(view, refdata);
    let stop_codon = frame.and_then(|f| view.coding_stop(f.start, f.end));
    let junction = v
        .position
        .zip(j.position)
        .filter(|&(v, j)| j + 3 > v && j + 3 <= view.len())
        .map(|(v, j)| (v, j + 3));
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
        reasons,
    }
}

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
}
