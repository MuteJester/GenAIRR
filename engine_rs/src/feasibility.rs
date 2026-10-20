//! Compiled downstream feasibility checks.
//!
//! Contracts answer "does this concrete candidate violate the current
//! state?" Feasibility answers the complementary compiled-runtime
//! question: "if this candidate is committed now, does any admissible
//! completion of the remaining stochastic choices still exist?"
//!
//! This module intentionally starts with the bounded case that matters
//! for VJ productive recombination. It is not a retry loop and it is not
//! a second pass scheduler: it is a typed, compile-derived support
//! artifact consulted by sampling passes before committing choices.

use crate::address::{self, ChoiceAddress, VdjSegment};
use crate::assignment::TrimEnd;
use crate::ir::{Segment, Simulation};
use crate::refdata::{AlleleId, RefDataConfig};
use crate::trace::ChoiceValue;

mod p_productive;
pub(crate) use p_productive::{
    NpCompletionDomain, PProductiveFeasibility, SegmentCompletionDomain,
};

/// Runtime view of compile-derived feasibility domains.
#[derive(Clone, Debug, Default)]
pub struct FeasibilityContext {
    vj_productive: Option<VjProductiveFeasibility>,
    coding_segments: Vec<CodingSegmentFeasibility>,
    p_productive: Option<PProductiveFeasibility>,
    coding_end_geometry: Vec<CodingEndGeometry>,
}

impl FeasibilityContext {
    pub(crate) fn empty() -> Self {
        Self::default()
    }

    pub(crate) fn with_vj_productive(mut self, feasibility: VjProductiveFeasibility) -> Self {
        self.vj_productive = Some(feasibility);
        self
    }

    pub(crate) fn is_empty(&self) -> bool {
        self.vj_productive.is_none()
            && self.coding_segments.is_empty()
            && self.p_productive.is_none()
            && self.coding_end_geometry.is_empty()
    }

    pub(crate) fn with_p_productive(mut self, feasibility: PProductiveFeasibility) -> Self {
        self.p_productive = Some(feasibility);
        self
    }

    pub(crate) fn has_p_productive(&self) -> bool {
        self.p_productive.is_some()
    }

    pub(crate) fn has_coding_end_geometry(&self) -> bool {
        !self.coding_end_geometry.is_empty()
    }

    pub(crate) fn with_coding_end_geometry(mut self, geometry: CodingEndGeometry) -> Self {
        self.coding_end_geometry.push(geometry);
        self
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn admits_with_allele_restrictions(
        &self,
        pass_index: usize,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
        address: &str,
        candidate: &ChoiceValue,
        allowed_alleles: &[(Segment, Vec<AlleleId>)],
    ) -> bool {
        self.admits(pass_index, sim, refdata, address, candidate)
            && self.p_productive.as_ref().map_or(true, |p| {
                p.admits_with_allele_restrictions(
                    pass_index,
                    sim,
                    refdata,
                    address,
                    candidate,
                    allowed_alleles,
                )
            })
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn admits_np_base(
        &self,
        pass_index: usize,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
        segment: Segment,
        index: u32,
        total_len: u32,
        candidate: u8,
    ) -> bool {
        self.p_productive.as_ref().map_or(true, |p| {
            refdata.is_some_and(|r| {
                p.admits_np_base(pass_index, sim, r, segment, index, total_len, candidate)
            })
        })
    }

    pub(crate) fn with_coding_segment(mut self, feasibility: CodingSegmentFeasibility) -> Self {
        self.coding_segments.push(feasibility);
        self
    }

    /// Return whether `candidate` at `address` keeps at least one
    /// compiled downstream path feasible.
    pub(crate) fn admits(
        &self,
        pass_index: usize,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
        address: &str,
        candidate: &ChoiceValue,
    ) -> bool {
        if self
            .coding_end_geometry
            .iter()
            .any(|geometry| !geometry.admits(pass_index, sim, refdata, address, candidate))
        {
            return false;
        }
        if let Some(p) = &self.p_productive {
            if !p.admits(pass_index, sim, refdata, address, candidate) {
                return false;
            }
        }
        for coding in &self.coding_segments {
            if !coding.admits(pass_index, sim, refdata, address, candidate) {
                return false;
            }
        }
        if let Some(vj_productive) = &self.vj_productive {
            if !vj_productive.admits(pass_index, sim, refdata, address, candidate) {
                return false;
            }
        }
        true
    }
}

/// Physical end geometry applies even when no productivity contracts are active.
/// It filters original allele choices before a too-short allele can be committed.
#[derive(Clone, Debug)]
pub(crate) struct CodingEndGeometry {
    pub(crate) segment: Segment,
    pub(crate) allele_pass: Option<usize>,
    pub(crate) outcomes: ChoiceDomain<crate::end_processing::EndOutcome>,
}

impl CodingEndGeometry {
    fn admits(
        &self,
        pass_index: usize,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
        address: &str,
        candidate: &ChoiceValue,
    ) -> bool {
        let parsed = ChoiceAddress::parse(address);
        let allele_choice = matches!(parsed, Some(ChoiceAddress::SampleAllele(s))
            if Segment::from(s) == self.segment)
            && self.allele_pass == Some(pass_index);
        let end_choice = matches!(parsed, Some(ChoiceAddress::CodingEndOutcome(s))
            if Segment::from(s) == self.segment)
            && self.outcomes.pass_index == Some(pass_index);
        if !allele_choice && !end_choice {
            return true;
        }
        let Some(refdata) = refdata else { return false };
        let id = if allele_choice {
            match candidate {
                ChoiceValue::AlleleId(id) => AlleleId::new(*id),
                _ => return false,
            }
        } else {
            let Some(instance) = sim.assignments.get(self.segment) else {
                return false;
            };
            instance.allele_id
        };
        if !refdata.coding_ends_are_not_known_partial(self.segment, id) {
            return false;
        }
        let Some(allele) = refdata.get(self.segment, id) else {
            return false;
        };
        let valid = |outcome: &crate::end_processing::EndOutcome| {
            crate::end_processing::materialize(
                &allele.seq,
                outcome.five,
                outcome.three,
                crate::assignment::SegmentOrientation::Forward,
            )
            .is_ok()
        };
        if end_choice {
            let ChoiceValue::Int(index) = candidate else {
                return false;
            };
            usize::try_from(*index)
                .ok()
                .and_then(|index| self.outcomes.values.get(index))
                .is_some_and(valid)
        } else {
            self.outcomes.values.iter().any(valid)
        }
    }
}

/// Stop-free V/J reference support conditioned on the actual future trim
/// domains. Filtering alleles before those domains are known would incorrectly
/// reject an allele whose stop is removed by an explicitly requested trim.
#[derive(Clone, Debug)]
pub(crate) struct CodingSegmentFeasibility {
    segment: Segment,
    alleles: ChoiceDomain<AlleleId>,
    trim_5: ChoiceDomain<u32>,
    trim_3: ChoiceDomain<u32>,
}

impl CodingSegmentFeasibility {
    pub(crate) fn new(
        segment: Segment,
        alleles: ChoiceDomain<AlleleId>,
        trim_5: ChoiceDomain<u32>,
        trim_3: ChoiceDomain<u32>,
    ) -> Self {
        Self {
            segment,
            alleles,
            trim_5,
            trim_3,
        }
    }

    fn admits(
        &self,
        pass_index: usize,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
        address: &str,
        candidate: &ChoiceValue,
    ) -> bool {
        let sample_address = sample_allele_address(self.segment);
        let five_address = trim_address(self.segment, TrimEnd::Five);
        let three_address = trim_address(self.segment, TrimEnd::Three);
        if address != sample_address && address != five_address && address != three_address {
            return true;
        }
        let Some(refdata) = refdata else { return true };
        let allele_ids = if address == sample_address {
            let ChoiceValue::AlleleId(id) = candidate else {
                return false;
            };
            vec![AlleleId::new(*id)]
        } else if let Some(inst) = sim.assignments.get(self.segment) {
            vec![inst.allele_id]
        } else {
            self.alleles.values.clone()
        };
        let trim_values = |domain: &ChoiceDomain<u32>, end: TrimEnd, target_address: &str| {
            if address == target_address {
                match candidate {
                    ChoiceValue::Int(n) if (0..=u16::MAX as i64).contains(n) => vec![*n as u32],
                    _ => Vec::new(),
                }
            } else if domain.is_committed_before(pass_index) {
                sim.assignments
                    .get(self.segment)
                    .map(|inst| {
                        vec![u32::from(match end {
                            TrimEnd::Five => inst.trim_5,
                            TrimEnd::Three => inst.trim_3,
                        })]
                    })
                    .unwrap_or_default()
            } else {
                domain.values.clone()
            }
        };
        let five = trim_values(&self.trim_5, TrimEnd::Five, five_address);
        let three = trim_values(&self.trim_3, TrimEnd::Three, three_address);
        allele_ids.into_iter().any(|id| {
            let Some(allele) = refdata.get(self.segment, id) else {
                return false;
            };
            let Some(anchor) = allele.anchor.map(u32::from) else {
                return false;
            };
            five.iter().any(|&t5| {
                three.iter().any(|&t3| {
                    t5 <= anchor
                        && allele
                            .seq
                            .len()
                            .checked_sub(t3 as usize)
                            .is_some_and(|end| anchor as usize + 3 <= end)
                        && crate::productivity::retained_stop_free(
                            refdata,
                            self.segment,
                            id,
                            t5,
                            t3,
                        )
                })
            })
        })
    }
}

/// Finite support for one stochastic choice plus the plan index that
/// owns it. `None` means the choice is absent from the plan and the
/// support is a fixed structural default, e.g. zero trim.
#[derive(Clone, Debug)]
pub(crate) struct ChoiceDomain<T> {
    pass_index: Option<usize>,
    values: Vec<T>,
}

impl<T> ChoiceDomain<T> {
    pub(crate) fn new(pass_index: Option<usize>, values: Vec<T>) -> Self {
        Self { pass_index, values }
    }

    pub(crate) fn pass_index(&self) -> Option<usize> {
        self.pass_index
    }

    fn is_committed_before(&self, pass_index: usize) -> bool {
        self.pass_index.is_some_and(|index| index < pass_index)
    }
}

/// Exact VJ productive-frame feasibility over V/J allele, V/J trim,
/// and NP1-length supports.
#[derive(Clone, Debug)]
pub(crate) struct VjProductiveFeasibility {
    v_alleles: ChoiceDomain<AlleleId>,
    j_alleles: ChoiceDomain<AlleleId>,
    v_trim_5: ChoiceDomain<u32>,
    v_trim_3: ChoiceDomain<u32>,
    j_trim_5: ChoiceDomain<u32>,
    j_trim_3: ChoiceDomain<u32>,
    np1_lengths: Vec<u32>,
    require_coding_integrity: bool,
}

impl VjProductiveFeasibility {
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn new(
        v_alleles: ChoiceDomain<AlleleId>,
        j_alleles: ChoiceDomain<AlleleId>,
        v_trim_5: ChoiceDomain<u32>,
        v_trim_3: ChoiceDomain<u32>,
        j_trim_5: ChoiceDomain<u32>,
        j_trim_3: ChoiceDomain<u32>,
        np1_lengths: Vec<u32>,
    ) -> Self {
        Self {
            v_alleles,
            j_alleles,
            v_trim_5,
            v_trim_3,
            j_trim_5,
            j_trim_3,
            np1_lengths,
            require_coding_integrity: false,
        }
    }

    pub(crate) fn with_coding_integrity(mut self) -> Self {
        self.require_coding_integrity = true;
        self
    }

    fn admits(
        &self,
        pass_index: usize,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
        address: &str,
        candidate: &ChoiceValue,
    ) -> bool {
        if !is_vj_productive_choice(address) {
            return true;
        }

        let Some(refdata) = refdata else {
            return true;
        };

        let Some(v_alleles) = self.allele_values(
            Segment::V,
            &self.v_alleles,
            pass_index,
            sim,
            address,
            candidate,
        ) else {
            return false;
        };
        let Some(j_alleles) = self.allele_values(
            Segment::J,
            &self.j_alleles,
            pass_index,
            sim,
            address,
            candidate,
        ) else {
            return false;
        };

        let Some(v_trim_5) = self.trim_values(
            Segment::V,
            TrimEnd::Five,
            &self.v_trim_5,
            pass_index,
            sim,
            address,
            candidate,
        ) else {
            return false;
        };
        let Some(v_trim_3) = self.trim_values(
            Segment::V,
            TrimEnd::Three,
            &self.v_trim_3,
            pass_index,
            sim,
            address,
            candidate,
        ) else {
            return false;
        };
        let Some(j_trim_5) = self.trim_values(
            Segment::J,
            TrimEnd::Five,
            &self.j_trim_5,
            pass_index,
            sim,
            address,
            candidate,
        ) else {
            return false;
        };
        let Some(j_trim_3) = self.trim_values(
            Segment::J,
            TrimEnd::Three,
            &self.j_trim_3,
            pass_index,
            sim,
            address,
            candidate,
        ) else {
            return false;
        };

        self.has_feasible_completion(
            refdata, &v_alleles, &j_alleles, &v_trim_5, &v_trim_3, &j_trim_5, &j_trim_3,
        )
    }

    fn allele_values(
        &self,
        segment: Segment,
        domain: &ChoiceDomain<AlleleId>,
        pass_index: usize,
        sim: &Simulation,
        address: &str,
        candidate: &ChoiceValue,
    ) -> Option<Vec<AlleleId>> {
        if address == sample_allele_address(segment) {
            return match candidate {
                ChoiceValue::AlleleId(index) => Some(vec![AlleleId::new(*index)]),
                _ => None,
            };
        }

        if domain.is_committed_before(pass_index) {
            return sim
                .assignments
                .get(segment)
                .map(|instance| vec![instance.allele_id]);
        }

        // Even at the *same* pass index, if the segment is already
        // assigned in the partial simulation, treat it as committed. The
        // consolidated `SampleGenotypePass` assigns V (and D) before
        // sampling J within one pass index; without this, J feasibility
        // would ignore the V already chosen and over-accept. The flat
        // one-segment-per-pass path never assigns another segment at the
        // same index, so this is a no-op there.
        if let Some(instance) = sim.assignments.get(segment) {
            return Some(vec![instance.allele_id]);
        }

        Some(domain.values.clone())
    }

    #[allow(clippy::too_many_arguments)]
    fn trim_values(
        &self,
        segment: Segment,
        end: TrimEnd,
        domain: &ChoiceDomain<u32>,
        pass_index: usize,
        sim: &Simulation,
        address: &str,
        candidate: &ChoiceValue,
    ) -> Option<Vec<u32>> {
        if address == trim_address(segment, end) {
            return match candidate {
                ChoiceValue::Int(value) if (0..=u16::MAX as i64).contains(value) => {
                    Some(vec![*value as u32])
                }
                _ => None,
            };
        }

        if domain.is_committed_before(pass_index) {
            let instance = sim.assignments.get(segment)?;
            let value = match end {
                TrimEnd::Five => instance.trim_5,
                TrimEnd::Three => instance.trim_3,
            };
            return Some(vec![value as u32]);
        }

        Some(domain.values.clone())
    }

    #[allow(clippy::too_many_arguments)]
    fn has_feasible_completion(
        &self,
        refdata: &RefDataConfig,
        v_alleles: &[AlleleId],
        j_alleles: &[AlleleId],
        v_trim_5: &[u32],
        v_trim_3: &[u32],
        j_trim_5: &[u32],
        j_trim_3: &[u32],
    ) -> bool {
        for v_id in v_alleles {
            let Some(v_allele) = refdata.get(Segment::V, *v_id) else {
                continue;
            };
            let Some(v_anchor) = v_allele.anchor.map(|anchor| anchor as u32) else {
                continue;
            };

            for v_five in v_trim_5 {
                for v_three in v_trim_3 {
                    let Some(v_tail) =
                        anchor_tail(v_allele.seq.as_slice(), v_anchor, *v_five, *v_three)
                    else {
                        continue;
                    };
                    if self.require_coding_integrity
                        && !crate::productivity::retained_stop_free(
                            refdata,
                            Segment::V,
                            *v_id,
                            *v_five,
                            *v_three,
                        )
                    {
                        continue;
                    }

                    for j_id in j_alleles {
                        let Some(j_allele) = refdata.get(Segment::J, *j_id) else {
                            continue;
                        };
                        let Some(j_anchor) = j_allele.anchor.map(|anchor| anchor as u32) else {
                            continue;
                        };

                        for j_five in j_trim_5 {
                            for j_three in j_trim_3 {
                                let Some(j_head) = anchor_head(
                                    j_allele.seq.as_slice(),
                                    j_anchor,
                                    *j_five,
                                    *j_three,
                                ) else {
                                    continue;
                                };

                                if self.require_coding_integrity
                                    && !crate::productivity::retained_stop_free(
                                        refdata,
                                        Segment::J,
                                        *j_id,
                                        *j_five,
                                        *j_three,
                                    )
                                {
                                    continue;
                                }

                                if self.np1_lengths.iter().any(|np_len| {
                                    productive_length_is_feasible(v_tail, *np_len, j_head)
                                }) {
                                    return true;
                                }
                            }
                        }
                    }
                }
            }
        }

        false
    }
}

fn is_vj_productive_choice(address: &str) -> bool {
    matches!(
        ChoiceAddress::parse(address),
        Some(ChoiceAddress::SampleAllele(VdjSegment::V | VdjSegment::J))
            | Some(ChoiceAddress::Trim {
                segment: VdjSegment::V | VdjSegment::J,
                ..
            })
    )
}

fn sample_allele_address(segment: Segment) -> &'static str {
    match segment {
        Segment::V | Segment::J => address::sample_allele_vdj(segment),
        _ => address::SAMPLE_ALLELE_UNSUPPORTED,
    }
}

fn trim_address(segment: Segment, end: TrimEnd) -> &'static str {
    match (segment, end) {
        (Segment::V | Segment::J, _) => address::trim_vdj(segment, end),
        _ => address::TRIM_UNSUPPORTED,
    }
}

fn anchor_tail(seq: &[u8], anchor: u32, trim_5: u32, trim_3: u32) -> Option<&[u8]> {
    let retained_end = (seq.len() as u32).checked_sub(trim_3)?;
    if trim_5 <= anchor && anchor + 3 <= retained_end {
        Some(&seq[anchor as usize..retained_end as usize])
    } else {
        None
    }
}

fn anchor_head(seq: &[u8], anchor: u32, trim_5: u32, trim_3: u32) -> Option<&[u8]> {
    let retained_end = (seq.len() as u32).checked_sub(trim_3)?;
    if trim_5 <= anchor && anchor + 3 <= retained_end {
        Some(&seq[trim_5 as usize..(anchor + 3) as usize])
    } else {
        None
    }
}

fn productive_length_is_feasible(v_tail: &[u8], np_len: u32, j_head: &[u8]) -> bool {
    let junction_len = v_tail.len() as u32 + np_len + j_head.len() as u32;
    if junction_len % 3 != 0 {
        return false;
    }

    no_known_stop_for_length(v_tail, np_len, j_head)
}

fn no_known_stop_for_length(v_tail: &[u8], np_len: u32, j_head: &[u8]) -> bool {
    let mut bases = Vec::with_capacity(v_tail.len() + np_len as usize + j_head.len());
    bases.extend(v_tail.iter().copied().map(Some));
    bases.extend((0..np_len).map(|_| None));
    bases.extend(j_head.iter().copied().map(Some));

    let mut offset = 0usize;
    while offset + 3 <= bases.len() {
        if let (Some(a), Some(b), Some(c)) = (bases[offset], bases[offset + 1], bases[offset + 2]) {
            if crate::productivity::codon_stop([a, b, c]) != Some(false) {
                return false;
            }
        }
        offset += 3;
    }

    true
}

#[cfg(test)]
mod tests;
