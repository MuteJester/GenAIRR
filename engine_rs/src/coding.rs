//! Coding coordinates follow the retained V coding origin, not a displaced
//! junction anchor. This distinguishes terminal coverage loss from an internal
//! insertion/deletion that changes the translation of the observed molecule.

use crate::ir::{NucHandle, Region, Segment, Simulation};
use crate::refdata::RefDataConfig;

#[derive(Copy, Clone, Debug)]
pub(crate) struct CodingFrame {
    /// Pool-coordinate phase of complete codon starts.
    pub(crate) phase: usize,
    /// First complete observed V codon; excludes a partial terminal codon.
    pub(crate) start: usize,
    /// Exclusive end of the observed V(D)J coding span.
    pub(crate) end: usize,
}

pub(crate) fn anchor_position(sim: &Simulation, region: &Region, anchor: u32) -> Option<u32> {
    if sim.duplex_provenance.is_some() {
        use crate::productivity::SequenceView;
        return (region.start.index()..region.end.index()).find(|&p| {
            sim.body_reference_position(region.segment, p as usize)
                .map(u32::from)
                == Some(anchor)
        });
    }
    (region.start.index()..region.end.index()).find(|&p| {
        sim.pool
            .get(NucHandle::new(p))
            .is_some_and(|n| n.germline_pos.get().map(u32::from) == Some(anchor))
    })
}

pub(crate) fn coding_frame(sim: &Simulation, refdata: &RefDataConfig) -> Option<CodingFrame> {
    crate::productivity::coding_span(sim, refdata).map(|frame| CodingFrame {
        phase: frame.phase,
        start: frame.start,
        end: frame.end,
    })
}

/// Pending assembly is allowed by incremental contracts; a fully assembled J
/// must have a definitely valid endpoint phase. Final emission is assessed by
/// the nullable observed-span policy.
pub(crate) fn j_tail_in_frame(
    sim: &Simulation,
    refdata: &RefDataConfig,
    frame: CodingFrame,
) -> bool {
    if !sim.sequence.regions.iter().any(|r| r.segment == Segment::J) {
        return true;
    }
    crate::productivity::j_tail(
        sim,
        refdata,
        crate::productivity::CodingSpan {
            phase: frame.phase,
            start: frame.start,
            end: frame.end,
        },
    ) == Some(true)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::airr_record::build_airr_record;
    use crate::assignment::AlleleInstance;
    use crate::ir::{flag, Nucleotide};
    use crate::pass::Outcome;
    use crate::refdata::{Allele, ChainType};
    use crate::trace::Trace;

    fn fixture(v: &[u8], np: &[u8], j: &[u8]) -> (RefDataConfig, Simulation) {
        let mut cfg = RefDataConfig::empty(ChainType::Vj);
        let mut sim = Simulation::new();
        for (segment, seq, anchor) in [(Segment::V, v, 6), (Segment::J, j, 3)] {
            let allele = Allele {
                name: format!("{segment:?}*01"),
                gene: format!("{segment:?}"),
                seq: seq.to_vec(),
                segment,
                anchor: Some(anchor),
                functional_status: None,
                subregions: Vec::new(),
            };
            let id = match segment {
                Segment::V => cfg.v_pool.push(allele),
                _ => cfg.j_pool.push(allele),
            };
            sim = sim.with_allele_assigned(segment, AlleleInstance::new(id));
            if segment == Segment::J && !np.is_empty() {
                let start = sim.pool.len() as u32;
                for &base in np {
                    sim = sim
                        .with_nucleotide_pushed(Nucleotide::synthetic(
                            base,
                            Segment::Np1,
                            flag::N_NUC,
                        ))
                        .0;
                }
                sim = sim.with_region_added(Region::new(
                    Segment::Np1,
                    NucHandle::new(start),
                    NucHandle::new(sim.pool.len() as u32),
                ));
            }
            let start = sim.pool.len() as u32;
            for (g, &base) in seq.iter().enumerate() {
                sim = sim
                    .with_nucleotide_pushed(Nucleotide::germline(base, g as u16, segment))
                    .0;
            }
            sim = sim.with_region_added(Region::new(
                segment,
                NucHandle::new(start),
                NucHandle::new(sim.pool.len() as u32),
            ));
        }
        (cfg, sim)
    }

    fn record(cfg: &RefDataConfig, sim: Simulation) -> crate::airr_record::AirrRecord {
        build_airr_record(
            &Outcome {
                revisions: vec![sim],
                pass_names: Vec::new(),
                trace: Trace::new(),
                events: Vec::new(),
            },
            cfg,
            "coding-test",
        )
    }

    #[test]
    fn stops_before_v_anchor_and_after_j_anchor_are_reported() {
        for (v, j) in [
            (b"ATGTAATGTAAA".as_slice(), b"AAATGGGGG".as_slice()),
            (b"ATGAAATGTAAA".as_slice(), b"AAATGGTAA".as_slice()),
        ] {
            let (cfg, sim) = fixture(v, b"", j);
            let rec = record(&cfg, sim);
            assert_eq!(rec.stop_codon, Some(true));
            assert_eq!(rec.productive, Some(false));
            assert_eq!(rec.vj_in_frame, Some(true));
            assert!(!rec.junction_aa.contains('*'));
        }
    }

    #[test]
    fn out_of_frame_junction_does_not_suppress_stop_status() {
        let (cfg, sim) = fixture(b"ATGAAATGTAAA", b"T", b"AAATGGGGG");
        let rec = record(&cfg, sim);
        assert_eq!(rec.vj_in_frame, Some(false));
        assert_eq!(rec.stop_codon, Some(true));
        assert_eq!(rec.productive, Some(false));
    }

    #[test]
    fn internal_and_first_base_insertions_keep_original_translation_origin() {
        let (cfg, sim) = fixture(b"ATGAAATGTAAA", b"", b"AAATGGGGG");
        for (position, expected) in [(2, "IEM*KMG"), (0, "HEM*KMG")] {
            let changed = sim.with_indel_inserted(
                position,
                Nucleotide::synthetic(b'C', Segment::V, flag::INDEL_INSERTED),
            );
            let rec = record(&cfg, changed);
            assert_eq!(rec.sequence_aa, expected);
            assert_eq!(rec.vj_in_frame, Some(false));
            assert_eq!(rec.productive, Some(false));
            assert_eq!(rec.stop_codon, Some(true));
        }
    }

    #[test]
    fn first_base_deletion_differs_from_explicit_terminal_coverage_loss() {
        let (cfg, sim) = fixture(b"ATGAAATGTAAA", b"", b"AAATGGGGG");
        let deleted = sim.with_indel_deleted(0);
        let rec = record(&cfg, deleted.clone());
        assert_eq!(rec.sequence_aa, "*NVKNG");
        assert_eq!(rec.productive, Some(false));
        let mut assignment = *deleted.assignments.get(Segment::V).unwrap();
        assignment.observation_loss_5 = 1;
        let rec = record(&cfg, deleted.with_allele_assigned(Segment::V, assignment));
        assert_eq!(rec.sequence_aa, "KCKKWG");
        assert_eq!(rec.productive, Some(true));
        assert_eq!(rec.stop_codon, Some(false));
    }

    #[test]
    fn j_tail_frameshift_is_unproductive_without_a_new_stop() {
        let (cfg, sim) = fixture(b"ATGAAATGTAAA", b"", b"AAATGGGGG");
        let changed = sim.with_indel_inserted(
            19,
            Nucleotide::synthetic(b'A', Segment::J, flag::INDEL_INSERTED),
        );
        let rec = record(&cfg, changed);
        assert_eq!(rec.vj_in_frame, Some(true));
        assert_eq!(rec.stop_codon, Some(false));
        assert_eq!(rec.productive, Some(false));
    }

    #[test]
    fn last_j_base_deletion_differs_from_explicit_terminal_coverage_loss() {
        use crate::contract::{CodingSequenceIntegrity, Contract};
        let (cfg, sim) = fixture(b"ATGAAATGTAAA", b"", b"AAATGGGGG");
        let contract = CodingSequenceIntegrity;
        for count in [1, 2] {
            let mut deleted = sim.clone();
            for _ in 0..count {
                deleted = deleted.with_indel_deleted(deleted.pool.len() as u32 - 1);
            }
            let rec = record(&cfg, deleted.clone());
            assert_eq!(rec.vj_in_frame, Some(true));
            assert_eq!(rec.stop_codon, Some(false));
            assert_eq!(rec.productive, Some(false));
            assert!(contract.verify(&deleted, Some(&cfg)).is_err());

            let mut assignment = *deleted.assignments.get(Segment::J).unwrap();
            assignment.observation_loss_3 = count;
            let lost = deleted.with_allele_assigned(Segment::J, assignment);
            let rec = record(&cfg, lost.clone());
            assert_eq!(rec.productive, Some(true));
            assert_eq!(rec.stop_codon, Some(false));
            assert!(contract.verify(&lost, Some(&cfg)).is_ok());
        }
    }

    #[test]
    fn productive_mutation_support_blocks_stops_on_both_sides_of_junction() {
        use crate::contract::{CodingSequenceIntegrity, Contract};
        let (cfg, sim) = fixture(b"ATGAAATGTAAA", b"", b"AAATGGAAA");
        let contract = CodingSequenceIntegrity;
        for position in [3, 18] {
            assert!(!contract.admits_fixed_base_at(
                &sim,
                Some(&cfg),
                NucHandle::new(position),
                b'T',
            ));
            assert!(contract.admits_fixed_base_at(
                &sim,
                Some(&cfg),
                NucHandle::new(position + 2),
                b'G',
            ));
        }
    }
}
