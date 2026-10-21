//! Explicit version-1 semantic field inventory; derived caches are rebuilt.
use super::{canonical_digest, F64Bits};
use crate::assignment::AlleleInstance;
use crate::ir::{DuplexProvenance, Region, Simulation, SimulationEvent};
use crate::pass::Outcome;
use serde_json::{json, Value};

fn region(r: &Region) -> Value {
    json!({"segment":format!("{:?}",r.segment),"start":r.start.index(),"end":r.end.index(),"frame_phase":r.frame_phase})
}
fn assignment(a: &AlleleInstance) -> Value {
    json!({
        "allele_id":a.allele_id.index(),"trim_5":a.trim_5,"trim_3":a.trim_3,
        "p_length_5":a.p_length_5,"p_length_3":a.p_length_3,
        "observation_loss_5":a.observation_loss_5,"observation_loss_3":a.observation_loss_3,
        "observed_ref_start":a.observed_ref_start,"observed_ref_end":a.observed_ref_end,
        "orientation":format!("{:?}",a.orientation),
        "receptor_revision_original_id":a.receptor_revision_original_id.map(|id|id.index()),"haplotype":a.haplotype
    })
}
fn duplex(p: &DuplexProvenance) -> Value {
    json!({
        "founder":{"annotated":p.founder.annotated.as_ref(),"bindings":p.founder.bindings.as_ref(),
            "body_claims":p.founder.body_claims,"p_metrics":p.founder.p_metrics,"selection":p.founder.selection},
        "links":p.links,"removed":p.removed,"next_edit_ordinal":p.next_edit_ordinal,
        "edits":p.edits,"attempt_stages":p.attempt_stages,"rev_comp":p.rev_comp
    })
}
pub fn semantic_state_digest(sim: &Simulation) -> Result<String, String> {
    if sim.coding_end_provenance.is_some() {
        return Err("coding-end snapshots are outside family replay v1".into());
    }
    let mut history_work = sim.pool.len() as u64 + 1;
    if let Some(p) = &sim.duplex_provenance {
        for n in [
            p.links.len(),
            p.removed.len(),
            p.edits.len(),
            p.attempt_stages.len(),
            p.founder.body_claims.len(),
            p.founder.annotated.product.ledger.len(),
            p.founder.annotated.evidence.events.len(),
            p.founder.annotated.evidence.placements.len(),
        ] {
            history_work = history_work
                .checked_add(n as u64)
                .ok_or("history commitment work overflow")?;
        }
    }
    super::primitive_work(history_work);
    let pool: Vec<_> = sim
        .pool
        .as_slice()
        .iter()
        .map(|n| {
            json!([
                n.base,
                n.germline,
                n.germline_pos.get(),
                format!("{:?}", n.segment),
                n.flags.bits()
            ])
        })
        .collect();
    let assignments: Vec<_> = sim
        .assignments
        .iter()
        .map(|(seg, a)| json!([format!("{seg:?}"), assignment(a)]))
        .collect();
    let support = sim.live_call_support.as_ref().map(|support| {
        [
            crate::ir::Segment::V,
            crate::ir::Segment::D,
            crate::ir::Segment::J,
        ]
        .iter()
        .map(|s| {
            json!([
                format!("{s:?}"),
                support
                    .alleles(*s)
                    .map(|ids| ids.iter().map(|id| id.index()).collect::<Vec<_>>())
            ])
        })
        .collect::<Vec<_>>()
    });
    canonical_digest(
        "semantic_state_v1",
        &json!({"pool":pool,
        "regions":sim.sequence.regions.iter().map(region).collect::<Vec<_>>(),"assignments":assignments,
        "live_call_support":support,"coding_end_provenance":Value::Null,
        "duplex_provenance":sim.duplex_provenance.as_ref().map(|p|duplex(p)),
        "mutation_count":sim.mutation_count,"mutation_burden":F64Bits::new(sim.mutation_burden)?}),
    )
}
fn event(e: &SimulationEvent) -> Value {
    use SimulationEvent::*;
    match e {
        DuplexAttemptStageRecorded { summary } => json!(["duplex_attempt_stage", summary]),
        DuplexProvenanceInstalled { provenance } => json!(["duplex_installed", duplex(provenance)]),
        DuplexEditContextRecorded { source_pass } => json!(["duplex_context", source_pass]),
        // These variants are not admitted by the duplex source whitelist. Still
        // commit their complete explicit scalar payload rather than dropping them.
        CodingEndChoiceRecorded { choice: c } => json!([
            "coding_choice",
            c.choice_id,
            format!("{:?}", c.segment),
            c.allele_id.index(),
            c.five,
            c.three,
            c.outcome_index
        ]),
        CodingEndPOriginRecorded { origin: o } => json!([
            "coding_p",
            o.choice_id,
            format!("{:?}", o.source_segment),
            o.allele_id.index(),
            format!("{:?}", o.end),
            format!("{:?}", o.orientation),
            o.template_position,
            o.original_base,
            o.query_position
        ]),
        BasePushed {
            handle,
            base,
            segment,
            germline_pos,
            flags,
        } => json!([
            "push",
            handle.index(),
            base,
            format!("{segment:?}"),
            germline_pos,
            flags.bits()
        ]),
        BaseChanged {
            handle,
            old_base,
            new_base,
            segment,
            germline_pos,
        } => json!([
            "change",
            handle.index(),
            old_base,
            new_base,
            format!("{segment:?}"),
            germline_pos
        ]),
        BaseDeleted { at, removed_base } => json!(["observation_delete", at, removed_base]),
        IndelInserted {
            at,
            base,
            segment,
            flags,
        } => json!(["insert", at, base, format!("{segment:?}"), flags.bits()]),
        IndelDeleted {
            at,
            removed_base,
            segment,
        } => json!(["delete", at, removed_base, format!("{segment:?}")]),
        AssignmentChanged { segment, old, new } => json!([
            "assignment",
            format!("{segment:?}"),
            old.as_ref().map(assignment),
            assignment(new)
        ]),
        TrimChanged {
            segment,
            end,
            old,
            new,
        } => json!(["trim", format!("{segment:?}"), format!("{end:?}"), old, new]),
        RegionAdded { region: r } => json!(["region_add", region(r)]),
        PRegionAdded { end, region: r } => json!(["p_region_add", format!("{end:?}"), region(r)]),
        RegionReplaced { old, new } => json!(["region_replace", region(old), region(new)]),
        ReverseComplementFlagRecorded { applied } => json!(["rev_comp", applied]),
        MutationCountChanged { old, new, delta } => json!(["mutation_count", old, new, delta]),
        OrientationChanged { segment, old, new } => json!([
            "orientation",
            format!("{segment:?}"),
            format!("{old:?}"),
            format!("{new:?}")
        ]),
        SegmentReplaced {
            segment,
            old_region,
            new_region,
            bytes_delta,
        } => json!([
            "segment_replace",
            format!("{segment:?}"),
            region(old_region),
            region(new_region),
            bytes_delta
        ]),
    }
}
pub fn outcome_events_digest(outcome: &Outcome) -> Result<String, String> {
    let mut records = Vec::new();
    for r in &outcome.events {
        super::primitive_work(r.simulation_events.len() as u64 + 1);
        records.push(
            json!({"pass_index":r.pass_index,"pass_name":r.pass_name,"kind":format!("{:?}",r.kind),
            "compile_effects":r.compile_effects.iter().map(|e|format!("{e:?}")).collect::<Vec<_>>(),
            "trace_span":[r.trace_span.start,r.trace_span.end],
            "pre":[r.pre.pool_len,r.pre.region_count,r.pre.assigned_allele_count],
            "post":[r.post.pool_len,r.post.region_count,r.post.assigned_allele_count],
            "simulation_events":r.simulation_events.iter().map(event).collect::<Vec<_>>() }),
        );
    }
    canonical_digest(
        "committed_events_v1",
        &json!({"pass_names":outcome.pass_names,"events":records}),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ir::{NucHandle, Nucleotide, Segment};
    #[test]
    fn state_inventory_includes_molecule_assignment_frame_and_burden_not_cache_versions() {
        let sim = Simulation::new()
            .with_nucleotide_pushed(Nucleotide::germline(b'A', 0, Segment::V))
            .0
            .with_region_added(Region::new(
                Segment::V,
                NucHandle::new(0),
                NucHandle::new(1),
            ));
        let before = semantic_state_digest(&sim).unwrap();
        assert_ne!(
            before,
            semantic_state_digest(&sim.with_base_changed(NucHandle::new(0), b'C')).unwrap()
        );
        assert_ne!(
            before,
            semantic_state_digest(&sim.with_mutation_count(1)).unwrap()
        );
        assert_ne!(
            before,
            semantic_state_digest(&sim.with_mutation_burden(-0.0)).unwrap()
        );
        let mut calls = sim.segment_calls.as_ref().clone();
        calls.version += 1;
        assert_eq!(
            before,
            semantic_state_digest(&sim.with_segment_calls(calls)).unwrap()
        );
    }
}
