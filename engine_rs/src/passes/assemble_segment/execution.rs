use crate::ir::{NucHandle, Nucleotide, Region, Simulation, SimulationBuilder};
use crate::pass::{Pass, PassContext, PassError};

use super::AssembleSegmentPass;

impl AssembleSegmentPass {
    pub(super) fn execute_with_validation(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
        strict: bool,
    ) -> Result<Simulation, PassError> {
        let refdata = match ctx.refdata {
            Some(refdata) => refdata,
            None if strict => return Err(PassError::missing_refdata(self.name())),
            None => {
                panic!(
                    "AssembleSegmentPass({:?}): PassContext.refdata is None — \
                     use PassRuntime::execute_with_refdata for plans containing \
                     assembly passes",
                    self.segment
                )
            }
        };

        let inst = match sim.assignments.get(self.segment).copied() {
            Some(inst) => inst,
            None if strict => return Err(PassError::missing_assignment(self.name(), self.segment)),
            None => {
                panic!(
                    "AssembleSegmentPass({:?}): no allele assigned — \
                     SampleAllelePass for this segment must run before assembly",
                    self.segment
                )
            }
        };

        let allele = match refdata.get(self.segment, inst.allele_id) {
            Some(allele) => allele,
            None if strict => {
                return Err(PassError::missing_allele(
                    self.name(),
                    self.segment,
                    inst.allele_id.index(),
                ));
            }
            None => {
                panic!(
                    "AssembleSegmentPass({:?}): allele_id {:?} out of bounds \
                     for refdata pool",
                    self.segment, inst.allele_id
                )
            }
        };

        let trim_5 = inst.trim_5 as u32;
        let trim_3 = inst.trim_3 as u32;
        let allele_len = allele.len();

        if strict && trim_5 + trim_3 > allele_len {
            return Err(PassError::invalid_plan_state(
                self.name(),
                "trim_exceeds_allele_length",
            ));
        }
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
        let cumulative_len: u32 = sim.sequence.regions.iter().map(|r| r.len()).sum();
        let frame_phase = (cumulative_len % 3) as u8;

        let seg = self.segment;
        let slice = &allele.seq[slice_start as usize..(slice_start + slice_len) as usize];

        // No codon-rail data is stored on `Region`; the pool is the
        // authoritative source. `NoStopCodonInJunction` reads straight
        // from `sim.pool`, AIRR's `sequence_aa` / `junction_aa`
        // re-translate from raw bytes (the per-region frame doesn't
        // match the junction frame in general). Callers that need
        // the rail call `crate::ir::compute_codon_rail(region, pool)`
        // against the current pool.
        //
        // The walker observer still runs when a `ReferenceMatchIndex`
        // is available — that one drives the live-call call set,
        // which IS read by AIRR's v_call / d_call / j_call projection.
        let seq_start = sim.pool.len() as u32;

        let mut builder = SimulationBuilder::from_simulation(sim.clone());

        let walker_attached = if let Some(reference_index) = ctx.reference_index {
            if let Some(segment_index) = reference_index.get(seg) {
                builder.attach_walker_observer(segment_index, seq_start);
                true
            } else {
                false
            }
        } else {
            false
        };

        for (i, &base) in slice.iter().enumerate() {
            let allele_pos = slice_start + i as u32;
            builder.push_nucleotide(Nucleotide::germline(base, allele_pos as u16, seg));
        }

        let seq_end = builder.peek().pool.len() as u32;

        let segment_call = if walker_attached {
            let sealed = builder.seal_walker_observer(seq_end);
            let base_version = builder.peek().segment_calls.version;
            let evidence_version = base_version.saturating_add(1);
            let reference_index = ctx
                .reference_index
                .expect("walker_attached implies ctx.reference_index is Some");
            let segment_index = reference_index
                .get(seg)
                .expect("walker_attached implies ref_index.get(seg) is Some");
            Some(sealed.finalize_with_extensions(
                builder.peek(),
                segment_index,
                evidence_version,
                seq_start,
                seq_end,
            ))
        } else {
            None
        };

        let current = builder.seal();

        let region_start = NucHandle::new(seq_start);
        let region_end = NucHandle::new(seq_end);
        // Region carries no codon-rail data. On-demand consumers call
        // `crate::ir::compute_codon_rail(&region, &pool)`.
        let region = Region::new(seg, region_start, region_end).with_frame_phase(frame_phase);
        let current = current.with_region_added(region);

        // If a walker observer was attached, stash the
        // observer-produced call on the segment_calls sidecar without
        // bumping the version. The post-pass `LiveCallRefreshHook`
        // calls `with_assembled_segment_live_call`, whose fast path
        // notices the pre-staged call (evidence_version ==
        // base.version + 1) and absorbs it by bumping the version to
        // match — the final version trajectory matches the original
        // two-walk path.
        if let Some(segment_call) = segment_call {
            let mut next_calls = (*current.segment_calls).clone();
            next_calls.stage_segment_call(segment_call);
            Ok(current.with_segment_calls(next_calls))
        } else {
            Ok(current)
        }
    }
}
