use crate::address;
use crate::ir::{translate_codon, NucHandle, Segment, Simulation, AMINO_STOP};
use crate::refdata::RefDataConfig;
use crate::trace::ChoiceValue;

use super::NoStopCodonInJunction;
use crate::contract::{ChoiceContext, Contract, ContractViolation};

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
    pub(super) fn parse_np_base_address(address: &str) -> Option<(Segment, u32)> {
        address::parse_np_base(address)
    }

    fn parse_np_length_address(address: &str) -> Option<Segment> {
        address::parse_np_length(address)
    }

    fn parse_filter_candidate(
        address: &str,
        candidate: &ChoiceValue,
        context: ChoiceContext<'_>,
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
        let v_inst = sim.assignments.get(Segment::V).copied()?;
        let j_inst = sim.assignments.get(Segment::J).copied()?;
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

        if sim.assignments.has(Segment::D) {
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

    pub(super) fn reject_known_stop(
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

    pub(super) fn admits_np_candidate(
        &self,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
        address: &str,
        candidate: &ChoiceValue,
        context: ChoiceContext<'_>,
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
}
