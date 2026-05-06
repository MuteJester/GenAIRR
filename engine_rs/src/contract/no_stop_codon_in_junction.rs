//! `NoStopCodonInJunction` — junction codons must translate to amino acids.

use crate::ir::{translate_codon, NucHandle, Segment, Simulation, AMINO_STOP};
use crate::junction::compute_junction;
use crate::refdata::RefDataConfig;
use crate::trace::ChoiceValue;

use super::{ChoiceContext, ChoiceKind, Contract, ContractViolation};

/// Verifies that no codon inside the junction translates to a stop
/// (TAA, TAG, TGA). Walks codons from `junction.start` in steps of
/// 3 and checks each translation.
///
/// Vacuously satisfied (returns `Ok`) when:
/// - no `RefDataConfig` is provided,
/// - the junction is undefined (`compute_junction` returns `None`),
/// - the junction is out of frame — that's a different contract's
///   concern (`ProductiveJunctionFrame`); we pass on rather than
///   double-reporting.
///
/// Returns a `ContractViolation` with the first stop's position
/// and codon string when a stop is detected.
///
/// **Note:** the contract reads from `sim.pool`, not from the
/// allele reference. Mutations / NP base substitutions are
/// reflected — this is the *current* pool's stop-codon status, not
/// the germline's.
pub struct NoStopCodonInJunction;

#[derive(Copy, Clone, Debug, Eq, PartialEq)]
enum NpFilterCandidate {
    Length {
        segment: Segment,
        length: u32,
    },
    Base {
        segment: Segment,
        index: u32,
        total_len: Option<u32>,
        base: u8,
    },
}

impl NoStopCodonInJunction {
    pub fn new() -> Self {
        Self
    }

    fn parse_np_base_address(address: &str) -> Option<(Segment, u32)> {
        let (segment, prefix) = if address.starts_with("np.np1.bases[") {
            (Segment::Np1, "np.np1.bases[")
        } else if address.starts_with("np.np2.bases[") {
            (Segment::Np2, "np.np2.bases[")
        } else {
            return None;
        };
        let rest = address.strip_prefix(prefix)?;
        let idx = rest.strip_suffix(']')?.parse::<u32>().ok()?;
        Some((segment, idx))
    }

    fn parse_np_length_address(address: &str) -> Option<Segment> {
        match address {
            "np.np1.length" => Some(Segment::Np1),
            "np.np2.length" => Some(Segment::Np2),
            _ => None,
        }
    }

    fn parse_targeted_substitution_base_address(address: &str) -> Option<u32> {
        for prefix in [
            "mutate.uniform.base[",
            "mutate.s5f.base[",
            "corrupt.pcr.error_base[",
            "corrupt.contaminant.bases[",
        ] {
            if let Some(rest) = address.strip_prefix(prefix) {
                return rest.strip_suffix(']')?.parse::<u32>().ok();
            }
        }
        None
    }

    fn is_targeted_substitution_candidate(address: &str, context: ChoiceContext) -> bool {
        context.kind == ChoiceKind::TargetedBaseSubstitution
            || Self::parse_targeted_substitution_base_address(address).is_some()
    }

    fn parse_filter_candidate(
        address: &str,
        candidate: &ChoiceValue,
        context: ChoiceContext,
    ) -> Option<NpFilterCandidate> {
        if let Some(segment) = Self::parse_np_length_address(address) {
            let length = match candidate {
                ChoiceValue::Int(n) if *n >= 0 => *n as u32,
                _ => return None,
            };
            return Some(NpFilterCandidate::Length { segment, length });
        }

        let (segment, index_from_address) = Self::parse_np_base_address(address)?;
        let base = match candidate {
            ChoiceValue::Base(b) => *b,
            _ => return None,
        };
        let index = context.draw_index.unwrap_or(index_from_address);
        Some(NpFilterCandidate::Base {
            segment,
            index,
            total_len: context.draw_count,
            base,
        })
    }

    fn append_pool_range(out: &mut Vec<Option<u8>>, sim: &Simulation, start: u32, end: u32) {
        for pos in start..end {
            out.push(sim.pool.get(NucHandle::new(pos)).map(|n| n.base));
        }
    }

    fn retained_bounds(allele_len: u32, trim_5: u16, trim_3: u16) -> Option<(u32, u32)> {
        let trim_5 = trim_5 as u32;
        let trim_3 = trim_3 as u32;
        if trim_5 + trim_3 > allele_len {
            return None;
        }
        Some((trim_5, allele_len - trim_3))
    }

    fn append_ref_slice(out: &mut Vec<Option<u8>>, seq: &[u8], start: u32, end: u32) {
        let end = end.min(seq.len() as u32);
        for pos in start..end {
            out.push(seq.get(pos as usize).copied());
        }
    }

    fn append_fixed_segment(
        out: &mut Vec<Option<u8>>,
        sim: &Simulation,
        refdata: &RefDataConfig,
        segment: Segment,
        allele_end_exclusive: Option<u32>,
    ) -> Option<()> {
        let inst = sim.assignments.get(segment).copied()?;
        let allele = refdata.get(segment, inst.allele_id)?;
        let (retained_start, retained_end) =
            Self::retained_bounds(allele.len(), inst.trim_5, inst.trim_3)?;
        let allele_end = allele_end_exclusive
            .unwrap_or(retained_end)
            .min(retained_end);
        if allele_end <= retained_start {
            return Some(());
        }

        if let Some(region) = sim.sequence.regions.iter().find(|r| r.segment == segment) {
            let pool_len = allele_end.saturating_sub(retained_start);
            let pool_end = region
                .start
                .index()
                .saturating_add(pool_len)
                .min(region.end.index());
            Self::append_pool_range(out, sim, region.start.index(), pool_end);
        } else {
            Self::append_ref_slice(out, &allele.seq, retained_start, allele_end);
        }

        Some(())
    }

    fn append_np_segment(
        out: &mut Vec<Option<u8>>,
        sim: &Simulation,
        segment: Segment,
        np_candidate: NpFilterCandidate,
    ) -> Option<bool> {
        if let Some(region) = sim.sequence.regions.iter().find(|r| r.segment == segment) {
            Self::append_pool_range(out, sim, region.start.index(), region.end.index());
            return Some(true);
        }

        match np_candidate {
            NpFilterCandidate::Length {
                segment: candidate_segment,
                length,
            } if candidate_segment == segment => {
                out.extend((0..length).map(|_| None));
                Some(true)
            }
            NpFilterCandidate::Base {
                segment: candidate_segment,
                index,
                total_len,
                base,
            } if candidate_segment == segment => {
                let count_to_append = total_len.unwrap_or(index.saturating_add(1));
                if count_to_append <= index {
                    return Some(false);
                }
                let current_np_start = (sim.pool.len() as u32).checked_sub(index)?;
                for i in 0..count_to_append {
                    if i < index {
                        let pos = current_np_start + i;
                        out.push(sim.pool.get(NucHandle::new(pos)).map(|n| n.base));
                    } else if i == index {
                        out.push(Some(base));
                    } else {
                        out.push(None);
                    }
                }
                Some(total_len.is_some())
            }
            _ => Some(false),
        }
    }

    fn hypothetical_junction_bases(
        sim: &Simulation,
        refdata: &RefDataConfig,
        np_candidate: NpFilterCandidate,
    ) -> Option<Vec<Option<u8>>> {
        let v_inst = sim.assignments.v?;
        let j_inst = sim.assignments.j?;
        let v_allele = refdata.get(Segment::V, v_inst.allele_id)?;
        let j_allele = refdata.get(Segment::J, j_inst.allele_id)?;
        let v_anchor = v_allele.anchor? as u32;
        let j_anchor = j_allele.anchor? as u32;

        let (v_retained_start, v_retained_end) =
            Self::retained_bounds(v_allele.len(), v_inst.trim_5, v_inst.trim_3)?;
        let (j_retained_start, j_retained_end) =
            Self::retained_bounds(j_allele.len(), j_inst.trim_5, j_inst.trim_3)?;
        if v_anchor < v_retained_start || v_anchor >= v_retained_end {
            return None;
        }
        if j_anchor < j_retained_start {
            return None;
        }

        let v_region = sim
            .sequence
            .regions
            .iter()
            .find(|r| r.segment == Segment::V)?;
        let v_anchor_pool = v_region.start.index() + (v_anchor - v_retained_start);

        let mut bases = Vec::new();
        Self::append_pool_range(&mut bases, sim, v_anchor_pool, v_region.end.index());

        if sim.assignments.d.is_some() {
            if !Self::append_np_segment(&mut bases, sim, Segment::Np1, np_candidate)? {
                return Some(bases);
            }
            Self::append_fixed_segment(&mut bases, sim, refdata, Segment::D, None)?;
            if !Self::append_np_segment(&mut bases, sim, Segment::Np2, np_candidate)? {
                return Some(bases);
            }
        } else if !Self::append_np_segment(&mut bases, sim, Segment::Np1, np_candidate)? {
            return Some(bases);
        }

        let j_junction_end = j_anchor.saturating_add(3).min(j_retained_end);
        Self::append_fixed_segment(&mut bases, sim, refdata, Segment::J, Some(j_junction_end))?;

        Some(bases)
    }

    fn reject_known_stop(
        &self,
        address: &str,
        bases: &[Option<u8>],
    ) -> Result<(), ContractViolation> {
        let mut offset = 0usize;
        while offset + 3 <= bases.len() {
            if let (Some(b1), Some(b2), Some(b3)) =
                (bases[offset], bases[offset + 1], bases[offset + 2])
            {
                if translate_codon(b1, b2, b3) == AMINO_STOP {
                    return Err(ContractViolation::new(
                        self.name(),
                        format!(
                            "candidate at {} would force stop codon '{}{}{}' at junction offset {}",
                            address, b1 as char, b2 as char, b3 as char, offset
                        ),
                    ));
                }
            }
            offset += 3;
        }
        Ok(())
    }

    fn admits_np_candidate(
        &self,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
        address: &str,
        candidate: &ChoiceValue,
        context: ChoiceContext,
    ) -> Result<(), ContractViolation> {
        let np_candidate = match Self::parse_filter_candidate(address, candidate, context) {
            None => return Ok(()),
            Some(c) => c,
        };
        let refdata = match refdata {
            None => return Ok(()),
            Some(r) => r,
        };
        let bases = match Self::hypothetical_junction_bases(sim, refdata, np_candidate) {
            None => return Ok(()),
            Some(b) => b,
        };
        self.reject_known_stop(address, &bases)
    }

    fn pool_base_with_candidate(
        sim: &Simulation,
        handle: NucHandle,
        target: NucHandle,
        candidate_base: u8,
    ) -> Option<u8> {
        if handle == target {
            Some(candidate_base)
        } else {
            sim.pool.get(handle).map(|n| n.base)
        }
    }

    fn admits_targeted_substitution_candidate(
        &self,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
        address: &str,
        candidate: &ChoiceValue,
        context: ChoiceContext,
    ) -> Result<(), ContractViolation> {
        if !Self::is_targeted_substitution_candidate(address, context) {
            return Ok(());
        }

        let candidate_base = match candidate {
            ChoiceValue::Base(b) => *b,
            _ => return Ok(()),
        };
        let target = match context.target {
            Some(h) => h,
            None => return Ok(()),
        };
        let refdata = match refdata {
            None => return Ok(()),
            Some(r) => r,
        };
        let junction = match compute_junction(sim, refdata) {
            None => return Ok(()),
            Some(j) => j,
        };

        // Frame violations are owned by ProductiveJunctionFrame. This
        // contract only reasons over codon-shaped triples.
        if !junction.is_in_frame() {
            return Ok(());
        }

        let start = junction.start.index();
        let end = junction.end.index();
        let target_idx = target.index();
        if target_idx < start || target_idx >= end {
            return Ok(());
        }

        let mut pos = start;
        while pos + 3 <= end {
            let b1 =
                Self::pool_base_with_candidate(sim, NucHandle::new(pos), target, candidate_base);
            let b2 = Self::pool_base_with_candidate(
                sim,
                NucHandle::new(pos + 1),
                target,
                candidate_base,
            );
            let b3 = Self::pool_base_with_candidate(
                sim,
                NucHandle::new(pos + 2),
                target,
                candidate_base,
            );

            if let (Some(b1), Some(b2), Some(b3)) = (b1, b2, b3) {
                if translate_codon(b1, b2, b3) == AMINO_STOP {
                    return Err(ContractViolation::new(
                        self.name(),
                        format!(
                            "candidate at {} for target {} would leave stop codon '{}{}{}' at junction position {}",
                            address,
                            target_idx,
                            b1 as char,
                            b2 as char,
                            b3 as char,
                            pos
                        ),
                    ));
                }
            }
            pos += 3;
        }

        Ok(())
    }
}

impl Default for NoStopCodonInJunction {
    fn default() -> Self {
        Self::new()
    }
}

impl Contract for NoStopCodonInJunction {
    fn name(&self) -> &str {
        "no_stop_codon_in_junction"
    }

    fn verify(
        &self,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
    ) -> Result<(), ContractViolation> {
        let refdata = match refdata {
            None => return Ok(()),
            Some(r) => r,
        };

        let junction = match compute_junction(sim, refdata) {
            None => return Ok(()),
            Some(j) => j,
        };

        // Out of frame → not our concern. ProductiveJunctionFrame
        // is the contract that owns the frame check; here we'd just
        // be walking codon-shaped triples that aren't real codons.
        if !junction.is_in_frame() {
            return Ok(());
        }

        let start = junction.start.index();
        let end = junction.end.index();

        let mut pos = start;
        while pos + 3 <= end {
            // Defensive lookups — pool pointers from the junction
            // window should always be valid for a well-formed
            // simulation, but we don't want to panic on a malformed
            // input. Treat missing handles as ambiguous.
            let b1 = sim
                .pool
                .get(NucHandle::new(pos))
                .map(|n| n.base)
                .unwrap_or(b'N');
            let b2 = sim
                .pool
                .get(NucHandle::new(pos + 1))
                .map(|n| n.base)
                .unwrap_or(b'N');
            let b3 = sim
                .pool
                .get(NucHandle::new(pos + 2))
                .map(|n| n.base)
                .unwrap_or(b'N');

            let aa = translate_codon(b1, b2, b3);
            if aa == AMINO_STOP {
                return Err(ContractViolation::new(
                    self.name(),
                    format!(
                        "stop codon at junction position {}: '{}{}{}'",
                        pos, b1 as char, b2 as char, b3 as char
                    ),
                ));
            }
            pos += 3;
        }

        Ok(())
    }

    fn admits(
        &self,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
        address: &str,
        candidate: &ChoiceValue,
    ) -> Result<(), ContractViolation> {
        self.admits_np_candidate(sim, refdata, address, candidate, ChoiceContext::none())
    }

    fn admits_with_context(
        &self,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
        address: &str,
        candidate: &ChoiceValue,
        context: ChoiceContext,
    ) -> Result<(), ContractViolation> {
        self.admits_np_candidate(sim, refdata, address, candidate, context)?;
        self.admits_targeted_substitution_candidate(sim, refdata, address, candidate, context)
    }
}

#[cfg(test)]
mod tests {
    use super::super::test_support::{make_assembled_sim_from_refdata, make_vj_with_anchor_codons};
    use super::*;
    use crate::assignment::AlleleInstance;
    use crate::ir::{Nucleotide, Region};
    use crate::refdata::{Allele, AlleleId, AllelePool, ChainType, RefDataConfig};

    #[test]
    fn no_stop_codon_name_is_canonical() {
        assert_eq!(
            NoStopCodonInJunction::new().name(),
            "no_stop_codon_in_junction"
        );
    }

    #[test]
    fn no_stop_codon_no_refdata_passes_vacuously() {
        let c = NoStopCodonInJunction::new();
        assert!(c.verify(&Simulation::new(), None).is_ok());
    }

    #[test]
    fn no_stop_codon_undefined_junction_passes_vacuously() {
        let cfg = make_vj_with_anchor_codons(b"GGG", b"TTT");
        let c = NoStopCodonInJunction::new();
        // No assignments at all → junction undefined.
        assert!(c.verify(&Simulation::new(), Some(&cfg)).is_ok());
    }

    #[test]
    fn no_stop_codon_in_frame_no_stops_passes() {
        // V anchor codon GGG (Gly), J anchor codon TTT (Phe).
        // Junction codons: GGG, TTT → no stops.
        let cfg = make_vj_with_anchor_codons(b"GGG", b"TTT");
        let c = NoStopCodonInJunction::new();
        let sim = make_assembled_sim_from_refdata(&cfg);

        assert!(c.verify(&sim, Some(&cfg)).is_ok());
    }

    #[test]
    fn no_stop_codon_admits_rejects_uniform_mutation_base_that_creates_stop() {
        // Junction codons start as TAC, GGG. Mutating the third base
        // of TAC to A would create TAA.
        let cfg = make_vj_with_anchor_codons(b"TAC", b"GGG");
        let c = NoStopCodonInJunction::new();
        let sim = make_assembled_sim_from_refdata(&cfg);

        let err = c
            .admits_with_context(
                &sim,
                Some(&cfg),
                "mutate.uniform.base[0]",
                &ChoiceValue::Base(b'A'),
                ChoiceContext::indexed_target(0, 1, NucHandle::new(8)),
            )
            .unwrap_err();

        assert_eq!(err.contract_name, "no_stop_codon_in_junction");
        assert!(err.reason.contains("TAA"), "{}", err.reason);
        assert!(c
            .admits_with_context(
                &sim,
                Some(&cfg),
                "mutate.uniform.base[0]",
                &ChoiceValue::Base(b'C'),
                ChoiceContext::indexed_target(0, 1, NucHandle::new(8)),
            )
            .is_ok());
    }

    #[test]
    fn no_stop_codon_admits_rejects_s5f_mutation_base_that_creates_stop() {
        // Same target-site contract path as uniform mutation, but
        // dispatched through the biologically important S5F address.
        let cfg = make_vj_with_anchor_codons(b"TAC", b"GGG");
        let c = NoStopCodonInJunction::new();
        let sim = make_assembled_sim_from_refdata(&cfg);

        let err = c
            .admits_with_context(
                &sim,
                Some(&cfg),
                "mutate.s5f.base[0]",
                &ChoiceValue::Base(b'A'),
                ChoiceContext::indexed_target(0, 1, NucHandle::new(8)),
            )
            .unwrap_err();

        assert_eq!(err.contract_name, "no_stop_codon_in_junction");
        assert!(err.reason.contains("TAA"), "{}", err.reason);
    }

    #[test]
    fn no_stop_codon_admits_rejects_pcr_error_base_that_creates_stop() {
        // Observation-stage substitutions use the same target-site
        // contract path as biological mutations.
        let cfg = make_vj_with_anchor_codons(b"TAC", b"GGG");
        let c = NoStopCodonInJunction::new();
        let sim = make_assembled_sim_from_refdata(&cfg);

        let err = c
            .admits_with_context(
                &sim,
                Some(&cfg),
                "corrupt.pcr.error_base[0]",
                &ChoiceValue::Base(b'A'),
                ChoiceContext::indexed_target(0, 1, NucHandle::new(8)),
            )
            .unwrap_err();

        assert_eq!(err.contract_name, "no_stop_codon_in_junction");
        assert!(err.reason.contains("TAA"), "{}", err.reason);
    }

    #[test]
    fn no_stop_codon_admits_rejects_contaminant_base_that_creates_stop() {
        // Contaminant replacement is also a target-site substitution.
        let cfg = make_vj_with_anchor_codons(b"AAA", b"GGG");
        let c = NoStopCodonInJunction::new();
        let sim = make_assembled_sim_from_refdata(&cfg);

        let err = c
            .admits_with_context(
                &sim,
                Some(&cfg),
                "corrupt.contaminant.bases[6]",
                &ChoiceValue::Base(b'T'),
                ChoiceContext::indexed_target(6, 12, NucHandle::new(6)),
            )
            .unwrap_err();

        assert_eq!(err.contract_name, "no_stop_codon_in_junction");
        assert!(err.reason.contains("TAA"), "{}", err.reason);
    }

    #[test]
    fn no_stop_codon_v_anchor_taa_violates() {
        // V anchor codon TAA (stop).
        let cfg = make_vj_with_anchor_codons(b"TAA", b"GGG");
        let c = NoStopCodonInJunction::new();
        let sim = make_assembled_sim_from_refdata(&cfg);

        let err = c.verify(&sim, Some(&cfg)).unwrap_err();
        assert_eq!(err.contract_name, "no_stop_codon_in_junction");
        assert!(err.reason.contains("stop"), "{}", err.reason);
        assert!(err.reason.contains("TAA"), "{}", err.reason);
    }

    #[test]
    fn no_stop_codon_v_anchor_tag_violates() {
        let cfg = make_vj_with_anchor_codons(b"TAG", b"GGG");
        let c = NoStopCodonInJunction::new();
        let sim = make_assembled_sim_from_refdata(&cfg);

        let err = c.verify(&sim, Some(&cfg)).unwrap_err();
        assert!(err.reason.contains("TAG"), "{}", err.reason);
    }

    #[test]
    fn no_stop_codon_v_anchor_tga_violates() {
        let cfg = make_vj_with_anchor_codons(b"TGA", b"GGG");
        let c = NoStopCodonInJunction::new();
        let sim = make_assembled_sim_from_refdata(&cfg);

        let err = c.verify(&sim, Some(&cfg)).unwrap_err();
        assert!(err.reason.contains("TGA"), "{}", err.reason);
    }

    #[test]
    fn no_stop_codon_j_anchor_stop_violates() {
        // V anchor GGG (fine), J anchor TAA (stop).
        let cfg = make_vj_with_anchor_codons(b"GGG", b"TAA");
        let c = NoStopCodonInJunction::new();
        let sim = make_assembled_sim_from_refdata(&cfg);

        let err = c.verify(&sim, Some(&cfg)).unwrap_err();
        assert!(err.reason.contains("TAA"), "{}", err.reason);
    }

    #[test]
    fn no_stop_codon_first_stop_reported() {
        // V anchor TAA (first), J anchor TAG (second). Should report V's TAA.
        let cfg = make_vj_with_anchor_codons(b"TAA", b"TAG");
        let c = NoStopCodonInJunction::new();
        let sim = make_assembled_sim_from_refdata(&cfg);

        let err = c.verify(&sim, Some(&cfg)).unwrap_err();
        assert!(err.reason.contains("TAA"), "{}", err.reason);
        // Position 6 is the V anchor pool position.
        assert!(err.reason.contains("position 6"), "{}", err.reason);
    }

    #[test]
    fn no_stop_codon_out_of_frame_junction_passes_vacuously() {
        // J anchor at 2 → J anchor in pool = 9 + 2 = 11.
        // V anchor at 6 → V anchor pool = 6. Junction = [6, 14). Length 8.
        // 8 % 3 != 0 → out of frame → contract returns Ok (deferred to
        // ProductiveJunctionFrame).
        let mut cfg = make_vj_with_anchor_codons(b"GGG", b"TTT");
        // Move J anchor to 2 to break frame.
        let j_allele = cfg.j_pool.get(AlleleId::new(0)).unwrap().clone();
        cfg.j_pool = AllelePool::from_vec(vec![Allele {
            anchor: Some(2),
            ..j_allele
        }]);

        let c = NoStopCodonInJunction::new();
        let sim = make_assembled_sim_from_refdata(&cfg);
        assert!(c.verify(&sim, Some(&cfg)).is_ok());
    }

    #[test]
    fn no_stop_codon_works_through_box_dyn() {
        let cfg = make_vj_with_anchor_codons(b"GGG", b"TTT");
        let sim = make_assembled_sim_from_refdata(&cfg);
        let c: Box<dyn Contract> = Box::new(NoStopCodonInJunction::new());
        assert!(c.verify(&sim, Some(&cfg)).is_ok());
        assert_eq!(c.name(), "no_stop_codon_in_junction");
    }

    fn make_partial_np_stop_filter_case() -> (RefDataConfig, Simulation) {
        let mut cfg = RefDataConfig::empty(ChainType::Vj);
        let _ = cfg.v_pool.push(Allele {
            name: "v_test*01".into(),
            gene: "v_test".into(),
            seq: b"GGGTA".to_vec(),
            segment: Segment::V,
            anchor: Some(0),
        });
        let _ = cfg.j_pool.push(Allele {
            name: "j_test*01".into(),
            gene: "j_test".into(),
            seq: b"TTTAAA".to_vec(),
            segment: Segment::J,
            anchor: Some(0),
        });

        let mut sim = Simulation::new();
        for (i, &b) in b"GGGTA".iter().enumerate() {
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(b, i as u16, Segment::V));
            sim = next;
        }
        sim = sim.with_region_added(Region::new(
            Segment::V,
            NucHandle::new(0),
            NucHandle::new(5),
        ));
        sim = sim
            .with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(0)))
            .with_allele_assigned(Segment::J, AlleleInstance::new(AlleleId::new(0)));

        (cfg, sim)
    }

    fn make_partial_np1_future_d_stop_case(v_seq: &[u8]) -> (RefDataConfig, Simulation) {
        let mut cfg = RefDataConfig::empty(ChainType::Vdj);
        let _ = cfg.v_pool.push(Allele {
            name: "v_test*01".into(),
            gene: "v_test".into(),
            seq: v_seq.to_vec(),
            segment: Segment::V,
            anchor: Some(0),
        });
        let _ = cfg.d_pool.push(Allele {
            name: "d_test*01".into(),
            gene: "d_test".into(),
            seq: b"AAC".to_vec(),
            segment: Segment::D,
            anchor: None,
        });
        let _ = cfg.j_pool.push(Allele {
            name: "j_test*01".into(),
            gene: "j_test".into(),
            seq: b"TTTAAA".to_vec(),
            segment: Segment::J,
            anchor: Some(0),
        });

        let mut sim = Simulation::new();
        for (i, &b) in v_seq.iter().enumerate() {
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(b, i as u16, Segment::V));
            sim = next;
        }
        sim = sim.with_region_added(Region::new(
            Segment::V,
            NucHandle::new(0),
            NucHandle::new(v_seq.len() as u32),
        ));
        sim = sim
            .with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(0)))
            .with_allele_assigned(Segment::D, AlleleInstance::new(AlleleId::new(0)))
            .with_allele_assigned(Segment::J, AlleleInstance::new(AlleleId::new(0)));

        (cfg, sim)
    }

    #[test]
    fn no_stop_codon_admits_rejects_np_base_that_completes_stop() {
        let (cfg, sim) = make_partial_np_stop_filter_case();
        let c = NoStopCodonInJunction::new();

        let err = c
            .admits(
                &sim,
                Some(&cfg),
                "np.np1.bases[0]",
                &ChoiceValue::Base(b'A'),
            )
            .unwrap_err();
        assert_eq!(err.contract_name, "no_stop_codon_in_junction");
        assert!(err.reason.contains("TAA"), "{}", err.reason);

        assert!(c
            .admits(
                &sim,
                Some(&cfg),
                "np.np1.bases[0]",
                &ChoiceValue::Base(b'C')
            )
            .is_ok());
    }

    #[test]
    fn no_stop_codon_admits_ignores_np_base_until_codon_complete() {
        let (cfg, sim) = make_partial_np_stop_filter_case();
        let c = NoStopCodonInJunction::new();
        let (sim, _) = sim.with_nucleotide_pushed(Nucleotide::synthetic(
            b'T',
            Segment::Np1,
            crate::ir::flag::N_NUC,
        ));

        // Candidate position is offset 6 from the junction start, so it
        // starts a new codon rather than completing one. The future third
        // base will be filtered when that draw is reached.
        assert!(c
            .admits(
                &sim,
                Some(&cfg),
                "np.np1.bases[1]",
                &ChoiceValue::Base(b'A')
            )
            .is_ok());
    }

    #[test]
    fn no_stop_codon_admits_rejects_np1_base_that_forces_future_d_stop() {
        let (cfg, sim) = make_partial_np1_future_d_stop_case(b"GGG");
        let c = NoStopCodonInJunction::new();

        let err = c
            .admits_with_context(
                &sim,
                Some(&cfg),
                "np.np1.bases[0]",
                &ChoiceValue::Base(b'T'),
                ChoiceContext::indexed(0, 1),
            )
            .unwrap_err();
        assert_eq!(err.contract_name, "no_stop_codon_in_junction");
        assert!(err.reason.contains("TAA"), "{}", err.reason);

        assert!(c
            .admits_with_context(
                &sim,
                Some(&cfg),
                "np.np1.bases[0]",
                &ChoiceValue::Base(b'C'),
                ChoiceContext::indexed(0, 1),
            )
            .is_ok());
    }

    #[test]
    fn no_stop_codon_admits_rejects_zero_np1_length_that_forces_future_d_stop() {
        let (cfg, sim) = make_partial_np1_future_d_stop_case(b"GGGT");
        let c = NoStopCodonInJunction::new();

        let err = c
            .admits(&sim, Some(&cfg), "np.np1.length", &ChoiceValue::Int(0))
            .unwrap_err();
        assert_eq!(err.contract_name, "no_stop_codon_in_junction");
        assert!(err.reason.contains("TAA"), "{}", err.reason);

        assert!(c
            .admits(&sim, Some(&cfg), "np.np1.length", &ChoiceValue::Int(1),)
            .is_ok());
    }
}
