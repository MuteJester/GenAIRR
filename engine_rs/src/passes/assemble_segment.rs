//! `AssembleSegmentPass` — copy a germline allele slice into the pool (C.8).

use crate::ir::{NucHandle, Nucleotide, Region, Segment, Simulation};
use crate::pass::{Pass, PassContext};

/// Assemble one germline segment (V, D, or J) from its assigned
/// `AlleleInstance` into the simulation pool.
///
/// Reads the allele's reference bases from
/// `ctx.refdata.expect(...)`, slices off `trim_5` from the start
/// and `trim_3` from the end, and pushes the retained bases as
/// germline-derived nucleotides into the pool. Then constructs a
/// `Region` for the segment with `frame_phase` chained from the
/// cumulative length of all prior regions, and recomputes its
/// codon rail.
///
/// **Determinism:** the pass makes no random choices.
/// `declared_choices()` returns the empty vector — there is
/// nothing to record to the trace.
///
/// **Pre-conditions:**
/// - An allele must be assigned to `segment` (panics otherwise —
///   propagated from `Simulation::assignments.get`).
/// - `ctx.refdata` must be `Some(...)` (panics otherwise —
///   `PassContext::refdata` is `Option` so this pass can be
///   constructed in any plan, but it requires the data at
///   execute time).
/// - `trim_5 + trim_3 <= allele.len()` (panics otherwise — caller
///   bug, asserted defensively).
///
/// **Frame phase:** computed as `(cumulative prior region length)
/// % 3`. This anchors the codon frame at the start of the *first*
/// region. Real biology anchors at the V Cys; that anchor-aware
/// offset is a future refinement (Phase D / E). For C.8 the raw
/// cumulative-length-mod-3 is correct as long as the user isn't
/// reading frame-relative properties yet.
pub struct AssembleSegmentPass {
    segment: Segment,
}

impl AssembleSegmentPass {
    /// Construct an assembly pass for the given segment.
    /// Panics if `segment` is not V, D, or J.
    pub fn new(segment: Segment) -> Self {
        match segment {
            Segment::V | Segment::D | Segment::J => {}
            _ => panic!(
                "AssembleSegmentPass: segment must be V, D, or J — got {:?}",
                segment
            ),
        }
        Self { segment }
    }

    pub fn segment(&self) -> Segment {
        self.segment
    }
}

impl Pass for AssembleSegmentPass {
    fn name(&self) -> &str {
        match self.segment {
            Segment::V => "assemble.v",
            Segment::D => "assemble.d",
            Segment::J => "assemble.j",
            _ => unreachable!("AssembleSegmentPass with non-V/D/J segment"),
        }
    }

    fn execute(&self, sim: &Simulation, ctx: &mut PassContext) -> Simulation {
        let refdata = ctx.refdata.unwrap_or_else(|| {
            panic!(
                "AssembleSegmentPass({:?}): PassContext.refdata is None — \
                 use PassRuntime::execute_with_refdata for plans containing \
                 assembly passes",
                self.segment
            )
        });

        let inst = sim
            .assignments
            .get(self.segment)
            .copied()
            .unwrap_or_else(|| {
                panic!(
                    "AssembleSegmentPass({:?}): no allele assigned — \
                 SampleAllelePass for this segment must run before assembly",
                    self.segment
                )
            });

        let allele = refdata
            .get(self.segment, inst.allele_id)
            .unwrap_or_else(|| {
                panic!(
                    "AssembleSegmentPass({:?}): allele_id {:?} out of bounds \
                 for refdata pool",
                    self.segment, inst.allele_id
                )
            });

        let trim_5 = inst.trim_5 as u32;
        let trim_3 = inst.trim_3 as u32;
        let allele_len = allele.len();

        assert!(
            trim_5 + trim_3 <= allele_len,
            "AssembleSegmentPass({:?}): trim_5 ({}) + trim_3 ({}) exceeds \
             allele length ({}) for {}",
            self.segment,
            trim_5,
            trim_3,
            allele_len,
            allele.name
        );

        let slice_start = trim_5;
        let slice_end = allele_len - trim_3;
        let slice_len = slice_end - slice_start;

        // Frame phase: cumulative length of prior regions, mod 3.
        // The frame is implicitly anchored at the start of the first
        // assembled region. Anchor-aware framing arrives in Phase D/E.
        let cumulative_len: u32 = sim.sequence.regions.iter().map(|r| r.len()).sum();
        let frame_phase = (cumulative_len % 3) as u8;

        let region_start = NucHandle::new(sim.pool.len() as u32);

        // Push the post-trim slice into the pool as germline nucleotides.
        // Each base carries its position-in-original-allele as `germline_pos`.
        let mut current = sim.clone();
        for i in 0..slice_len {
            let allele_pos = slice_start + i;
            let base = allele.seq[allele_pos as usize];
            let (next, _h) = current.with_nucleotide_pushed(Nucleotide::germline(
                base,
                allele_pos as u16,
                self.segment,
            ));
            current = next;
        }

        let region_end = NucHandle::new(current.pool.len() as u32);
        let region = Region::new(self.segment, region_start, region_end)
            .with_frame_phase(frame_phase)
            .with_codon_rail_recomputed(&current.pool);

        current.with_region_added(region)
    }

    fn declared_choices(&self) -> Vec<String> {
        // Deterministic — no choices.
        Vec::new()
    }
}

#[cfg(test)]
pub(super) mod test_support {
    use crate::refdata::{Allele, ChainType, RefDataConfig};

    use crate::ir::Segment;

    /// Build a tiny reference config with one V, D, and J allele
    /// each of known length / sequence for assembly tests.
    pub fn make_test_refdata() -> RefDataConfig {
        let mut cfg = RefDataConfig::empty(ChainType::Vdj);
        let _ = cfg.v_pool.push(Allele {
            name: "v_test*01".to_string(),
            gene: "v_test".to_string(),
            seq: b"AAACCCGGG".to_vec(), // 9 bases
            segment: Segment::V,
            anchor: Some(6), // C of CCC at position 6
        });
        let _ = cfg.d_pool.push(Allele {
            name: "d_test*01".to_string(),
            gene: "d_test".to_string(),
            seq: b"TTTTTT".to_vec(), // 6 bases
            segment: Segment::D,
            anchor: None,
        });
        let _ = cfg.j_pool.push(Allele {
            name: "j_test*01".to_string(),
            gene: "j_test".to_string(),
            seq: b"GGGCCC".to_vec(), // 6 bases
            segment: Segment::J,
            anchor: Some(0),
        });
        cfg
    }
}

#[cfg(test)]
mod tests {
    use super::test_support::make_test_refdata;
    use super::*;
    use crate::assignment::TrimEnd;
    use crate::dist::{AllelePoolDist, EmpiricalLengthDist};
    use crate::pass::{PassPlan, PassRuntime};
    use crate::passes::sample_allele::test_support::make_test_pool;
    use crate::passes::{SampleAllelePass, TrimPass};

    #[test]
    #[should_panic(expected = "segment must be V, D, or J")]
    fn assemble_segment_pass_rejects_np1() {
        let _ = AssembleSegmentPass::new(Segment::Np1);
    }

    #[test]
    #[should_panic(expected = "PassContext.refdata is None")]
    fn assemble_segment_pass_panics_without_refdata() {
        // Use plain `execute` (refdata = None) — should panic.
        let mut plan = PassPlan::new();
        let v_pool = make_test_pool(1, Segment::V);
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(&v_pool)),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
        let _ = PassRuntime::execute(&plan, Simulation::new(), 0);
    }

    #[test]
    #[should_panic(expected = "no allele assigned")]
    fn assemble_segment_pass_panics_without_allele_assignment() {
        let refdata = make_test_refdata();
        let mut plan = PassPlan::new();
        plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
        let _ = PassRuntime::execute_with_refdata(&plan, Simulation::new(), 0, &refdata);
    }

    #[test]
    #[should_panic(expected = "exceeds allele length")]
    fn assemble_segment_pass_panics_when_trim_exceeds_length() {
        let refdata = make_test_refdata();
        // V is 9 bases. Trim 5 + 3 = 8, leaves 1 base — fine.
        // Trim 5' = 5, 3' = 5 → sum 10 > 9 → panic.
        let v_pool = make_test_pool(1, Segment::V);
        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(&v_pool)),
        )));
        plan.push(Box::new(TrimPass::new(
            Segment::V,
            TrimEnd::Five,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(5, 1.0)])),
        )));
        plan.push(Box::new(TrimPass::new(
            Segment::V,
            TrimEnd::Three,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(5, 1.0)])),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));

        // Note: SampleAllelePass uses make_test_pool, but AssembleSegmentPass
        // reads from refdata. They MUST be the same allele for the
        // test to succeed; we rely on test_refdata having allele 0
        // be the V whose length is 9.
        let _ = PassRuntime::execute_with_refdata(&plan, Simulation::new(), 0, &refdata);
    }

    #[test]
    fn assemble_segment_pass_v_no_trim_pushes_full_allele() {
        let refdata = make_test_refdata();
        // Use the V pool from refdata directly so SampleAllele sees
        // the same alleles as Assemble. We sample id 0 (the only V).
        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(&refdata.v_pool)),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));

        let outcome = PassRuntime::execute_with_refdata(&plan, Simulation::new(), 0, &refdata);
        let sim = outcome.final_simulation();

        // V allele was 9 bases; no trim → 9 nucleotides in pool.
        assert_eq!(sim.pool.len(), 9);
        // One region (V) of length 9.
        assert_eq!(sim.sequence.region_count(), 1);
        let r = &sim.sequence.regions[0];
        assert_eq!(r.segment, Segment::V);
        assert_eq!(r.len(), 9);
        assert_eq!(r.frame_phase, 0); // first region

        // Bases match the allele sequence exactly, with germline_pos
        // matching position-in-allele.
        let expected = b"AAACCCGGG";
        for i in 0..9 {
            let n = sim.pool.get(NucHandle::new(i)).unwrap();
            assert_eq!(n.base, expected[i as usize]);
            assert_eq!(n.germline, expected[i as usize]);
            assert_eq!(n.germline_pos, i as u16);
            assert_eq!(n.segment, Segment::V);
        }
    }

    #[test]
    fn assemble_segment_pass_v_with_trim_pushes_post_trim_slice() {
        let refdata = make_test_refdata();
        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(&refdata.v_pool)),
        )));
        plan.push(Box::new(TrimPass::new(
            Segment::V,
            TrimEnd::Five,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(2, 1.0)])),
        )));
        plan.push(Box::new(TrimPass::new(
            Segment::V,
            TrimEnd::Three,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)])),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));

        let outcome = PassRuntime::execute_with_refdata(&plan, Simulation::new(), 0, &refdata);
        let sim = outcome.final_simulation();

        // Allele = "AAACCCGGG" (9 bases); trim_5=2, trim_3=3.
        // Slice = "ACCCG" (positions 2..6, length 4 — wait let me recompute).
        // Actually 9 - 2 - 3 = 4 → slice is positions [2, 6) = "ACCC".
        assert_eq!(sim.pool.len(), 4);
        let expected = b"ACCC";
        for i in 0..4 {
            let n = sim.pool.get(NucHandle::new(i)).unwrap();
            assert_eq!(n.base, expected[i as usize]);
            assert_eq!(n.germline_pos, (i + 2) as u16); // shifted by trim_5
        }
    }

    #[test]
    fn assemble_segment_pass_chains_frame_phase_across_regions() {
        // V (9 bases, frame_phase=0) + assemble J (6 bases,
        // frame_phase = 9 % 3 = 0). So both regions have phase 0.
        // Test that the chain is correct after assembly.
        let refdata = make_test_refdata();
        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(&refdata.v_pool)),
        )));
        plan.push(Box::new(SampleAllelePass::new(
            Segment::J,
            Box::new(AllelePoolDist::uniform(&refdata.j_pool)),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::J)));

        let outcome = PassRuntime::execute_with_refdata(&plan, Simulation::new(), 0, &refdata);
        let sim = outcome.final_simulation();

        assert_eq!(sim.pool.len(), 15); // 9 V + 6 J
        assert_eq!(sim.sequence.region_count(), 2);
        assert_eq!(sim.sequence.regions[0].segment, Segment::V);
        assert_eq!(sim.sequence.regions[0].frame_phase, 0);
        assert_eq!(sim.sequence.regions[1].segment, Segment::J);
        // J starts after 9 V bases; 9 % 3 = 0 → frame_phase 0.
        assert_eq!(sim.sequence.regions[1].frame_phase, 0);
    }

    #[test]
    fn assemble_segment_pass_chains_frame_phase_with_offset() {
        // Build a plan where V has 2 bases (after heavy trim), so
        // J's frame_phase = 2 % 3 = 2. Tests the non-zero phase chain.
        let refdata = make_test_refdata();
        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(&refdata.v_pool)),
        )));
        plan.push(Box::new(SampleAllelePass::new(
            Segment::J,
            Box::new(AllelePoolDist::uniform(&refdata.j_pool)),
        )));
        // V is 9 bases. Trim 5'=4, 3'=3 → slice has 9-4-3 = 2 bases.
        plan.push(Box::new(TrimPass::new(
            Segment::V,
            TrimEnd::Five,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(4, 1.0)])),
        )));
        plan.push(Box::new(TrimPass::new(
            Segment::V,
            TrimEnd::Three,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)])),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::J)));

        let outcome = PassRuntime::execute_with_refdata(&plan, Simulation::new(), 0, &refdata);
        let sim = outcome.final_simulation();

        assert_eq!(sim.sequence.regions[0].len(), 2);
        assert_eq!(sim.sequence.regions[0].frame_phase, 0);
        // J's first base sits at position 2 of the cumulative chain;
        // 2 % 3 = 2 → frame_phase 2.
        assert_eq!(sim.sequence.regions[1].frame_phase, 2);
    }

    #[test]
    fn assemble_segment_pass_does_not_record_to_trace() {
        // Determinism contract for transform passes: no trace
        // entries.
        let refdata = make_test_refdata();
        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(&refdata.v_pool)),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));

        let outcome = PassRuntime::execute_with_refdata(&plan, Simulation::new(), 7, &refdata);

        // Trace contains only the SampleAllele entry, none from Assemble.
        assert_eq!(outcome.trace.len(), 1);
        assert_eq!(outcome.trace.choices()[0].address, "sample_allele.v");
    }

    #[test]
    fn assemble_segment_pass_declared_choices_is_empty() {
        let pass = AssembleSegmentPass::new(Segment::V);
        assert!(pass.declared_choices().is_empty());
    }

    #[test]
    fn assemble_segment_pass_codon_rail_computed_on_assembled_region() {
        // V allele "AAACCCGGG" (9 bases), no trim, frame_phase 0
        // → codons AAA, CCC, GGG → K, P, G.
        let refdata = make_test_refdata();
        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(&refdata.v_pool)),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));

        let outcome = PassRuntime::execute_with_refdata(&plan, Simulation::new(), 0, &refdata);
        let sim = outcome.final_simulation();
        let r = &sim.sequence.regions[0];

        assert_eq!(r.amino_acids, b"KPG");
        assert!(r.stop_codon_positions.is_empty());
    }
}
