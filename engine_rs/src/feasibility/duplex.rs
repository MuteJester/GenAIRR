//! Boolean completion support for the finite authored duplex law.
//! Weights remain in their original local samplers; no history mass is folded
//! into upstream allele, gene, chromosome, or orientation probabilities.
use crate::address::ChoiceAddress;
use crate::duplex_repair::model::PreparedDuplexHistory;
use crate::ir::{Segment, Simulation};
use crate::refdata::AlleleId;
use crate::trace::ChoiceValue;

#[derive(Clone, Debug)]
pub(crate) struct DuplexCompletionFeasibility {
    pub(crate) source_passes: [usize; 3],
    pub(crate) orientation_pass: Option<usize>,
    pub(crate) repair_pass: usize,
    pub(crate) rows: Vec<PreparedDuplexHistory>,
}

impl DuplexCompletionFeasibility {
    pub(crate) fn owns_source_pass(&self, index: usize) -> bool {
        self.source_passes.contains(&index) || self.orientation_pass == Some(index)
    }

    pub(crate) fn admits(
        &self,
        pass_index: usize,
        sim: &Simulation,
        address: &str,
        candidate: &ChoiceValue,
        restricted: &[(Segment, Vec<AlleleId>)],
    ) -> bool {
        let parsed = ChoiceAddress::parse(address);
        let candidate_segment = match parsed {
            Some(ChoiceAddress::SampleAllele(s)) => Some(Segment::from(s)),
            _ => None,
        };
        let segments = [Segment::V, Segment::D, Segment::J];
        let owner = candidate_segment
            .and_then(|segment| {
                segments
                    .iter()
                    .position(|s| *s == segment)
                    .map(|i| self.source_passes[i])
            })
            .or_else(|| {
                matches!(parsed, Some(ChoiceAddress::SampleAlleleDInverted))
                    .then_some(self.orientation_pass)
                    .flatten()
            });
        if owner != Some(pass_index) {
            return true;
        }
        self.rows.iter().any(|row| {
            let ids = [row.context.v, row.context.d, row.context.j];
            for (index, segment) in segments.into_iter().enumerate() {
                if let Some((_, allowed)) = restricted.iter().find(|(s, _)| *s == segment) {
                    if !allowed.contains(&ids[index]) {
                        return false;
                    }
                }
                if candidate_segment == Some(segment) {
                    if !matches!(candidate, ChoiceValue::AlleleId(id) if *id == ids[index].index()) {
                        return false;
                    }
                } else if self.source_passes[index] <= pass_index {
                    // Genotype commits V/D/J within one pass. Its already
                    // assigned slots bind the prefix, while absent slots are
                    // still future choices on the same chromosome.
                    if let Some(assigned) = sim.assignments.get(segment) {
                        if assigned.allele_id != ids[index] {
                            return false;
                        }
                    } else if self.source_passes[index] < pass_index {
                        return false;
                    }
                }
            }
            if matches!(parsed, Some(ChoiceAddress::SampleAlleleDInverted)) {
                return matches!(candidate, ChoiceValue::Bool(value) if *value == row.context.d_inverted);
            }
            if self.orientation_pass.is_none() {
                !row.context.d_inverted
            } else if self.orientation_pass.is_some_and(|owner| owner < pass_index) {
                sim.assignments.get(Segment::D)
                    .is_some_and(|d| d.orientation.is_reverse() == row.context.d_inverted)
            } else {
                true
            }
        })
    }
}
