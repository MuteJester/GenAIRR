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

        // Phase 19: single canonical assembly path. Build the new
        // region's bytes via `SimulationBuilder` so observers see each
        // base as it lands. The codon-rail observer ALWAYS runs (it
        // populates the new region's `amino_acids` /
        // `stop_codon_positions` byte-for-byte equivalent to a
        // post-push `with_codon_rail_recomputed`). The walker observer
        // runs ONLY when a `ReferenceMatchIndex` is available
        // (compiled execution path) — test harnesses like
        // `PassRuntime::execute_with_refdata` don't compile a
        // reference index, so they pay nothing extra for assembly
        // beyond the codon rail.
        //
        // When the walker observer is attached, the post-pass
        // `with_assembled_segment_live_call` hook detects the
        // pre-staged call via its `evidence_version == version + 1`
        // fast path and reuses it instead of running
        // `call_from_region` from scratch.
        let seq_start = sim.pool.len() as u32;

        let mut builder = SimulationBuilder::from_simulation(sim.clone());
        builder.attach_codon_rail_observer(seg, seq_start, frame_phase);

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
        let sealed_rail = builder.seal_codon_rail_observer();

        let segment_call = if walker_attached {
            let sealed = builder.seal_walker_observer(seq_end);
            // The `evidence_version` we stamp into the observer-produced
            // call must equal `base_state.version + 1` so the post-pass
            // `with_assembled_segment_live_call` fast-path check fires.
            let base_version = builder
                .peek()
                .live_calls
                .as_ref()
                .map(|s| s.version)
                .unwrap_or(0);
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
        // Region's codon rail comes from the sealed observer,
        // byte-for-byte equivalent to running
        // `with_codon_rail_recomputed(&current.pool)` on this range.
        let region =
            Region::from_sealed_codon_rail(seg, region_start, region_end, frame_phase, sealed_rail);
        let current = current.with_region_added(region);

        // If a walker observer was attached, stash the
        // observer-produced call on the live-calls sidecar without
        // bumping the LiveCallState version. The post-pass
        // `apply_live_call_updates` hook (in `compiled/execute.rs`)
        // calls `with_assembled_segment_live_call`, whose fast path
        // notices the pre-staged call (evidence_version == base.version + 1)
        // and absorbs it by bumping the version to match — the final
        // version trajectory matches the original two-walk path.
        if let Some(segment_call) = segment_call {
            let base_state = current
                .live_calls
                .as_ref()
                .map(|s| (**s).clone())
                .unwrap_or_default();
            let next_state = base_state.with_segment_call_observed(segment_call);
            Ok(current.with_live_calls(next_state))
        } else {
            Ok(current)
        }
    }
}
