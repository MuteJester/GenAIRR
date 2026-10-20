//! Original-coding-end processing: one joint choice, deterministic P emission.
//!
//! The authored row index is the stochastic choice. Both strict and permissive
//! execution reject impossible support; no zero-length replacement is invented.
//! The pure kernel is also used by compile-time productive completion.

use crate::address::{ChoiceAddress, ChoiceAddressPattern, PEnd, VdjSegment};
use crate::assignment::TrimEnd;
use crate::contract::ChoiceContext;
use crate::dist::{sample_filtered_result, Distribution, FilteredSampleError};
use crate::end_processing::{materialize, EndOutcome};
use crate::ir::{NucHandle, Region, Segment, Simulation, SimulationBuilder};
use crate::pass::{Pass, PassCompileFact, PassContext, PassEffect, PassError, PassRequirement};
use crate::rng::Rng;
use crate::trace::ChoiceValue;

/// Choose one joint canonical-reference signed net-deletion outcome.
pub struct ProcessCodingEndsPass {
    segment: Segment,
    outcomes: Vec<EndOutcome>,
}

/// Keep the exact authored weights and row indices, including duplicate
/// geometries. An empirical CDF reconstructed by differencing cumulative
/// weights could lose very small authored weights.
struct OutcomeIndices<'a>(&'a [EndOutcome]);

impl Distribution for OutcomeIndices<'_> {
    type Output = i64;

    fn sample(&self, rng: &mut Rng) -> i64 {
        let target = rng.next_f64() * self.0.iter().map(|o| o.weight).sum::<f64>();
        let mut cumulative = 0.0;
        for (index, row) in self.0.iter().enumerate() {
            cumulative += row.weight;
            if target < cumulative {
                return index as i64;
            }
        }
        (self.0.len() - 1) as i64
    }

    fn support(&self) -> Option<Vec<(i64, f64)>> {
        Some(
            self.0
                .iter()
                .enumerate()
                .map(|(i, row)| (i as i64, row.weight))
                .collect(),
        )
    }
}

impl ProcessCodingEndsPass {
    pub fn new(segment: Segment, outcomes: Vec<EndOutcome>) -> Result<Self, String> {
        if !matches!(segment, Segment::V | Segment::D | Segment::J) {
            return Err("coding_end_processing_requires_v_d_or_j".into());
        }
        if outcomes.is_empty() || i64::try_from(outcomes.len()).is_err() {
            return Err("coding_end_processing_requires_nonempty_finite_support".into());
        }
        let mut total = 0.0;
        for (index, row) in outcomes.iter().enumerate() {
            if !row.weight.is_finite() || row.weight <= 0.0 {
                return Err(format!("coding_end_invalid_weight_at_row_{index}"));
            }
            total += row.weight;
            if (segment == Segment::V && row.five != 0) || (segment == Segment::J && row.three != 0)
            {
                return Err(format!("coding_end_nonjunction_end_at_row_{index}"));
            }
            if row.five > i64::from(u16::MAX) || row.three > i64::from(u16::MAX) {
                return Err(format!(
                    "coding_end_trim_exceeds_ir_capacity_at_row_{index}"
                ));
            }
        }
        if !total.is_finite() || total <= 0.0 {
            return Err("coding_end_invalid_total_weight".into());
        }
        Ok(Self { segment, outcomes })
    }

    pub fn segment(&self) -> Segment {
        self.segment
    }

    pub fn outcomes(&self) -> &[EndOutcome] {
        &self.outcomes
    }

    fn typed_segment(&self) -> VdjSegment {
        self.segment
            .try_into()
            .expect("constructor validates V/D/J")
    }

    fn apply(&self, sim: &Simulation, ctx: &mut PassContext) -> Result<Simulation, PassError> {
        let refdata = ctx
            .refdata
            .ok_or_else(|| PassError::missing_refdata(self.name()))?;
        let inst = sim
            .assignments
            .get(self.segment)
            .copied()
            .ok_or_else(|| PassError::missing_assignment(self.name(), self.segment))?;
        if self.segment != Segment::D && inst.orientation.is_reverse() {
            return Err(PassError::invalid_plan_state(
                self.name(),
                "coding_end_v1_supports_only_d_inversion",
            ));
        }
        if !refdata.coding_ends_are_not_known_partial(self.segment, inst.allele_id) {
            return Err(PassError::invalid_plan_state(
                self.name(),
                "coding_end_reference_has_known_partial_boundary",
            ));
        }
        let allele = refdata.get(self.segment, inst.allele_id).ok_or_else(|| {
            PassError::missing_allele(self.name(), self.segment, inst.allele_id.index())
        })?;
        if sim.coding_end_choice(self.segment).is_some()
            || inst.trim_5 != 0
            || inst.trim_3 != 0
            || inst.p_length_5 != 0
            || inst.p_length_3 != 0
            || sim
                .sequence
                .regions
                .iter()
                .any(|r| r.segment == self.segment)
        {
            return Err(PassError::invalid_plan_state(
                self.name(),
                "coding_end_processing_requires_unprocessed_unassembled_segment",
            ));
        }
        let choice_address = ChoiceAddress::CodingEndOutcome(self.typed_segment());
        let address = choice_address.to_string();
        let feasibility = ctx.feasibility;
        let contracts = ctx.contracts;
        let pass_index = ctx.pass_index;
        let admits = |index: &i64| {
            let Some(row) = usize::try_from(*index)
                .ok()
                .and_then(|i| self.outcomes.get(i))
            else {
                return false;
            };
            let choice = ChoiceValue::Int(*index);
            materialize(&allele.seq, row.five, row.three, inst.orientation).is_ok()
                && feasibility
                    .is_none_or(|f| f.admits(pass_index, sim, Some(refdata), &address, &choice))
                && contracts.is_none_or(|k| {
                    k.admits_typed(
                        sim,
                        Some(refdata),
                        ChoiceContext::none().with_address(choice_address),
                        &choice,
                    )
                    .is_ok()
                })
        };
        let index = if let Some(cursor) = ctx.replay_cursor.as_deref_mut() {
            cursor
                .expect_int(choice_address)
                .map_err(|reason| PassError::replay(self.name(), reason))?
        } else {
            sample_filtered_result(ctx.rng, &OutcomeIndices(&self.outcomes), admits)
                .map_err(|reason| PassError::constraint_sampling(self.name(), &address, reason))?
        };
        if !admits(&index) {
            return Err(PassError::constraint_sampling(
                self.name(),
                address,
                FilteredSampleError::EmptyAdmissibleSupport,
            ));
        }
        let outcome_index = index as usize;
        let row = self.outcomes[outcome_index];
        let mut builder = SimulationBuilder::from_simulation(sim.clone());
        if ctx.event_log_sink.is_some() {
            builder.attach_event_log_observer();
        }
        builder.update_trim(self.segment, TrimEnd::Five, row.five.max(0) as u16);
        builder.update_trim(self.segment, TrimEnd::Three, row.three.max(0) as u16);
        builder.record_coding_end_choice(
            self.segment,
            inst.allele_id,
            row.five,
            row.three,
            outcome_index,
        );
        ctx.trace
            .record_choice(choice_address, ChoiceValue::Int(index));
        if let Some(sink) = ctx.event_log_sink.as_deref_mut() {
            sink.extend(builder.seal_event_log_observer());
        }
        Ok(builder.seal())
    }
}

impl Pass for ProcessCodingEndsPass {
    fn name(&self) -> &str {
        match self.segment {
            Segment::V => "end_processing.v1.v",
            Segment::D => "end_processing.v1.d",
            Segment::J => "end_processing.v1.j",
            _ => unreachable!(),
        }
    }

    fn parameter_signature(&self) -> String {
        let rows = self
            .outcomes
            .iter()
            .map(|r| format!("{}:{}:{:016x}", r.five, r.three, r.weight.to_bits()))
            .collect::<Vec<_>>()
            .join(",");
        format!("model=original_coding_end_v1;outcomes=[{rows}]")
    }

    fn execute(&self, sim: &Simulation, ctx: &mut PassContext) -> Simulation {
        self.apply(sim, ctx)
            .expect("original coding-end processing requires valid support")
    }

    fn execute_checked(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
    ) -> Result<Simulation, PassError> {
        self.apply(sim, ctx)
    }

    fn requirements(&self) -> Vec<PassRequirement> {
        vec![
            PassRequirement::RefData,
            PassRequirement::AlleleAssignment(self.segment),
        ]
    }

    fn effects(&self) -> Vec<PassEffect> {
        vec![PassEffect::TrimAllele(self.segment)]
    }

    fn declared_choice_patterns(&self) -> Vec<ChoiceAddressPattern> {
        vec![ChoiceAddressPattern::CodingEndOutcome(self.typed_segment())]
    }

    fn compile_facts(&self) -> Vec<PassCompileFact> {
        vec![PassCompileFact::CodingEndSupport {
            segment: self.segment,
            outcomes: self.outcomes.clone(),
        }]
    }
}

/// Emit one assembled P flank from the already committed original-end choice.
///
/// Zero-length ends still have an emission pass, but create neither a base nor
/// a stochastic P-length address. Compile facts expose every end to validation.
pub struct EmitCodingEndPPass {
    end: PEnd,
}

impl EmitCodingEndPPass {
    pub fn new(end: PEnd) -> Self {
        Self { end }
    }

    pub fn end(&self) -> PEnd {
        self.end
    }

    fn apply(&self, sim: &Simulation, ctx: &mut PassContext) -> Result<Simulation, PassError> {
        let segment = self.end.source_segment();
        let refdata = ctx
            .refdata
            .ok_or_else(|| PassError::missing_refdata(self.name()))?;
        let inst = sim
            .assignments
            .get(segment)
            .copied()
            .ok_or_else(|| PassError::missing_assignment(self.name(), segment))?;
        if segment != Segment::D && inst.orientation.is_reverse() {
            return Err(PassError::invalid_plan_state(
                self.name(),
                "coding_end_v1_supports_only_d_inversion",
            ));
        }
        let choice = sim.coding_end_choice(segment).ok_or_else(|| {
            PassError::invalid_plan_state(self.name(), "coding_end_emission_missing_choice")
        })?;
        if choice.allele_id != inst.allele_id
            || choice.five.max(0) != i64::from(inst.trim_5)
            || choice.three.max(0) != i64::from(inst.trim_3)
        {
            return Err(PassError::invalid_plan_state(
                self.name(),
                "coding_end_emission_assignment_changed",
            ));
        }
        if !refdata.coding_ends_are_not_known_partial(segment, choice.allele_id) {
            return Err(PassError::invalid_plan_state(
                self.name(),
                "coding_end_reference_has_known_partial_boundary",
            ));
        }
        let allele = refdata.get(segment, choice.allele_id).ok_or_else(|| {
            PassError::missing_allele(self.name(), segment, choice.allele_id.index())
        })?;
        let processed = materialize(&allele.seq, choice.five, choice.three, inst.orientation)
            .map_err(|reason| PassError::invalid_plan_state(self.name(), reason))?;
        let is_five = matches!(self.end, PEnd::D5 | PEnd::J5);
        let (bytes, positions, previous) = if is_five {
            (
                &processed.p5,
                &processed.p5_template_positions,
                inst.p_length_5,
            )
        } else {
            (
                &processed.p3,
                &processed.p3_template_positions,
                inst.p_length_3,
            )
        };
        if previous != 0 {
            return Err(PassError::invalid_plan_state(
                self.name(),
                "coding_end_emission_already_materialized",
            ));
        }
        if bytes.is_empty() {
            return Ok(sim.clone());
        }
        let start = sim.pool.len() as u32;
        let mut builder = SimulationBuilder::from_simulation(sim.clone());
        if ctx.event_log_sink.is_some() {
            builder.attach_event_log_observer();
        }
        for (&base, &position) in bytes.iter().zip(positions) {
            builder.push_coding_end_p(base, segment, self.end, position);
        }
        let updated = if is_five {
            inst.with_p_length_5(bytes.len() as u32)
        } else {
            inst.with_p_length_3(bytes.len() as u32)
        };
        builder.assign_allele(segment, updated);
        let end = builder.peek().pool.len() as u32;
        builder.record_p_region(
            self.end,
            Region::new(segment, NucHandle::new(start), NucHandle::new(end))
                .with_frame_phase((start % 3) as u8),
        );
        if let Some(sink) = ctx.event_log_sink.as_deref_mut() {
            sink.extend(builder.seal_event_log_observer());
        }
        Ok(builder.seal())
    }
}

impl Pass for EmitCodingEndPPass {
    fn name(&self) -> &str {
        match self.end {
            PEnd::V3 => "emit_coding_end_p.v1.v_3",
            PEnd::D5 => "emit_coding_end_p.v1.d_5",
            PEnd::D3 => "emit_coding_end_p.v1.d_3",
            PEnd::J5 => "emit_coding_end_p.v1.j_5",
        }
    }

    fn parameter_signature(&self) -> String {
        "model=original_coding_end_v1".into()
    }

    fn execute(&self, sim: &Simulation, ctx: &mut PassContext) -> Simulation {
        self.apply(sim, ctx)
            .expect("coding-end emission requires a valid committed outcome")
    }

    fn execute_checked(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
    ) -> Result<Simulation, PassError> {
        self.apply(sim, ctx)
    }

    fn requirements(&self) -> Vec<PassRequirement> {
        vec![
            PassRequirement::RefData,
            PassRequirement::AlleleAssignment(self.end.source_segment()),
        ]
    }

    fn effects(&self) -> Vec<PassEffect> {
        vec![PassEffect::AppendNucleotides]
    }

    fn compile_facts(&self) -> Vec<PassCompileFact> {
        vec![PassCompileFact::CodingEndEmission { end: self.end }]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::assignment::{AlleleInstance, SegmentOrientation};
    use crate::ir::SimulationEvent;
    use crate::passes::test_support::{run_pass_capturing_events, run_pass_with_replay_records};
    use crate::refdata::{Allele, AlleleId, ChainType, RefDataConfig};
    use crate::trace::ChoiceRecord;

    fn fixture(segment: Segment, orientation: SegmentOrientation) -> (RefDataConfig, Simulation) {
        let mut cfg = RefDataConfig::empty(ChainType::Vdj);
        let allele = Allele {
            name: "source*01".into(),
            gene: "source".into(),
            seq: b"AACGTC".to_vec(),
            segment,
            anchor: None,
            functional_status: None,
            subregions: Vec::new(),
        };
        match segment {
            Segment::V => cfg.v_pool.push(allele),
            Segment::D => cfg.d_pool.push(allele),
            Segment::J => cfg.j_pool.push(allele),
            _ => unreachable!(),
        };
        let sim = Simulation::new().with_allele_assigned(
            segment,
            AlleleInstance::new(AlleleId::new(0)).with_orientation(orientation),
        );
        (cfg, sim)
    }

    fn row(five: i64, three: i64, weight: f64) -> EndOutcome {
        EndOutcome {
            five,
            three,
            weight,
        }
    }

    #[test]
    fn constructor_rejects_invalid_tables_and_nonjunction_ends() {
        for outcomes in [
            vec![],
            vec![row(0, 0, 0.0)],
            vec![row(0, 0, -1.0)],
            vec![row(0, 0, f64::NAN)],
            vec![row(0, 0, f64::INFINITY)],
            vec![row(0, 0, f64::MAX), row(0, 0, f64::MAX)],
        ] {
            assert!(ProcessCodingEndsPass::new(Segment::D, outcomes).is_err());
        }
        assert!(ProcessCodingEndsPass::new(Segment::V, vec![row(-1, 0, 1.0)]).is_err());
        assert!(ProcessCodingEndsPass::new(Segment::J, vec![row(0, 1, 1.0)]).is_err());
        assert!(ProcessCodingEndsPass::new(Segment::Np1, vec![row(0, 0, 1.0)]).is_err());
    }

    #[test]
    fn physical_filter_preserves_authored_row_index_and_records_trim_events() {
        let (cfg, sim) = fixture(Segment::D, SegmentOrientation::Forward);
        let pass = ProcessCodingEndsPass::new(
            Segment::D,
            vec![row(4, 3, 10.0), row(2, -2, 1.0), row(-7, 0, 50.0)],
        )
        .unwrap();
        let run = run_pass_capturing_events(&pass, &sim, Some(&cfg), None).unwrap();
        assert_eq!(run.trace.choices()[0].value, ChoiceValue::Int(1));
        assert_eq!(run.sim.assignments.get(Segment::D).unwrap().trim_5, 2);
        assert_eq!(run.sim.assignments.get(Segment::D).unwrap().trim_3, 0);
        assert_eq!(
            run.sim.coding_end_choice(Segment::D).unwrap().outcome_index,
            1
        );
        assert!(run
            .events
            .iter()
            .any(|e| matches!(e, SimulationEvent::TrimChanged { .. })));
        assert!(run
            .events
            .iter()
            .any(|e| matches!(e, SimulationEvent::CodingEndChoiceRecorded { .. })));
    }

    #[test]
    fn replay_rejects_unknown_or_impossible_indices_before_any_commit() {
        let (cfg, sim) = fixture(Segment::D, SegmentOrientation::Forward);
        let pass =
            ProcessCodingEndsPass::new(Segment::D, vec![row(1, 1, 1.0), row(4, 3, 1.0)]).unwrap();
        for index in [-1, 1, 2, i64::MAX] {
            let (result, trace, events) = run_pass_with_replay_records(
                &pass,
                &sim,
                vec![ChoiceRecord::new(
                    ChoiceAddress::CodingEndOutcome(VdjSegment::D),
                    ChoiceValue::Int(index),
                )],
                Some(&cfg),
                None,
            );
            assert!(result.is_err(), "index={index}");
            assert!(trace.choices().is_empty());
            assert!(events.is_empty());
        }
    }

    #[test]
    fn deterministic_reverse_d_emission_tracks_original_templates_and_has_no_draw() {
        let (cfg, sim) = fixture(Segment::D, SegmentOrientation::ReverseComplement);
        let process = ProcessCodingEndsPass::new(Segment::D, vec![row(-2, -3, 1.0)]).unwrap();
        let selected = run_pass_capturing_events(&process, &sim, Some(&cfg), None).unwrap();
        let emit = EmitCodingEndPPass::new(PEnd::D5);
        let emitted = run_pass_capturing_events(&emit, &selected.sim, Some(&cfg), None).unwrap();
        assert_eq!(
            emitted
                .sim
                .pool
                .as_slice()
                .iter()
                .map(|n| n.base)
                .collect::<Vec<_>>(),
            b"GTC"
        );
        assert_eq!(
            emitted.sim.assignments.get(Segment::D).unwrap().p_length_5,
            3
        );
        assert!(emitted.trace.choices().is_empty());
        assert!(emit.declared_choice_patterns().is_empty());
        assert!(emitted
            .events
            .iter()
            .any(|e| matches!(e, SimulationEvent::PRegionAdded { end: PEnd::D5, .. })));
        let origins = &emitted
            .sim
            .coding_end_provenance
            .as_ref()
            .unwrap()
            .p_origins;
        assert_eq!(
            origins
                .iter()
                .map(|p| p.template_position)
                .collect::<Vec<_>>(),
            [3, 4, 5]
        );
        assert!(origins
            .iter()
            .all(|p| p.allele_id == AlleleId::new(0) && p.end == PEnd::D5));
        assert!(run_pass_capturing_events(&emit, &emitted.sim, Some(&cfg), None).is_err());
    }

    #[test]
    fn zero_extension_is_noop_but_missing_choice_is_error() {
        let (cfg, sim) = fixture(Segment::V, SegmentOrientation::Forward);
        let emit = EmitCodingEndPPass::new(PEnd::V3);
        assert!(run_pass_capturing_events(&emit, &sim, Some(&cfg), None).is_err());
        let process = ProcessCodingEndsPass::new(Segment::V, vec![row(0, 2, 1.0)]).unwrap();
        let selected = run_pass_capturing_events(&process, &sim, Some(&cfg), None).unwrap();
        let run = run_pass_capturing_events(&emit, &selected.sim, Some(&cfg), None).unwrap();
        assert!(run.sim.pool.is_empty());
        assert!(run.events.is_empty());
        assert!(run.trace.choices().is_empty());
    }

    #[test]
    fn row_order_and_weights_are_part_of_replay_signature() {
        let a =
            ProcessCodingEndsPass::new(Segment::D, vec![row(1, 0, 1.0), row(0, -1, 2.0)]).unwrap();
        let reordered =
            ProcessCodingEndsPass::new(Segment::D, vec![row(0, -1, 2.0), row(1, 0, 1.0)]).unwrap();
        let reweighted =
            ProcessCodingEndsPass::new(Segment::D, vec![row(1, 0, 1.0), row(0, -1, 3.0)]).unwrap();
        assert_ne!(a.parameter_signature(), reordered.parameter_signature());
        assert_ne!(a.parameter_signature(), reweighted.parameter_signature());
        let address = ChoiceAddress::CodingEndOutcome(VdjSegment::D);
        assert_eq!(ChoiceAddress::parse(&address.to_string()), Some(address));
        assert_eq!(
            ChoiceAddressPattern::parse(&address.to_string()),
            Some(ChoiceAddressPattern::CodingEndOutcome(VdjSegment::D)),
        );
    }

    #[test]
    fn physical_filter_preserves_one_to_three_weights_without_invalid_mass() {
        let (cfg, sim) = fixture(Segment::D, SegmentOrientation::Forward);
        let pass = ProcessCodingEndsPass::new(
            Segment::D,
            vec![row(1, 0, 1.0), row(4, 3, 1_000.0), row(0, -1, 3.0)],
        )
        .unwrap();
        let trials = 2048;
        let mut high_weight = 0;
        for seed in 0..trials {
            let mut rng = Rng::new(seed as u64);
            let mut trace = crate::trace::Trace::new();
            let mut ctx = PassContext {
                trace: &mut trace,
                rng: &mut rng,
                pass_index: 0,
                refdata: Some(&cfg),
                contracts: None,
                feasibility: None,
                reference_index: None,
                replay_cursor: None,
                event_log_sink: None,
            };
            let next = pass.execute_checked(&sim, &mut ctx).unwrap();
            let index = next.coding_end_choice(Segment::D).unwrap().outcome_index;
            assert!(index == 0 || index == 2, "invalid geometry was sampled");
            high_weight += usize::from(index == 2);
        }
        let expectation = trials as f64 * 0.75;
        let sigma = (trials as f64 * 0.75 * 0.25).sqrt();
        assert!(
            (high_weight as f64 - expectation).abs() <= 5.0 * sigma,
            "valid row 2 count {high_weight}; expected {expectation} ± {}",
            5.0 * sigma,
        );
    }

    #[test]
    fn a_valid_replayed_index_preserves_choice_and_derived_emission() {
        let (cfg, sim) = fixture(Segment::D, SegmentOrientation::Forward);
        let process =
            ProcessCodingEndsPass::new(Segment::D, vec![row(2, 0, 1000.0), row(0, -2, 1.0)])
                .unwrap();
        let (selected, trace, _) = run_pass_with_replay_records(
            &process,
            &sim,
            vec![ChoiceRecord::new(
                ChoiceAddress::CodingEndOutcome(VdjSegment::D),
                ChoiceValue::Int(1),
            )],
            Some(&cfg),
            None,
        );
        let selected = selected.unwrap();
        assert_eq!(trace.choices()[0].value, ChoiceValue::Int(1));
        assert_eq!(selected.coding_end_choice(Segment::D).unwrap().three, -2);
        let emitted = run_pass_capturing_events(
            &EmitCodingEndPPass::new(PEnd::D3),
            &selected,
            Some(&cfg),
            None,
        )
        .unwrap();
        assert_eq!(
            emitted
                .sim
                .pool
                .as_slice()
                .iter()
                .map(|n| n.base)
                .collect::<Vec<_>>(),
            b"GA",
        );
        assert!(emitted.trace.choices().is_empty());
    }

    #[test]
    fn unsupported_v_j_inversion_cannot_silently_drop_a_forced_extension() {
        for (segment, outcome, end) in [
            (Segment::V, row(0, -2, 1.0), PEnd::V3),
            (Segment::J, row(-2, 0, 1.0), PEnd::J5),
        ] {
            let (cfg, reversed) = fixture(segment, SegmentOrientation::ReverseComplement);
            let process = ProcessCodingEndsPass::new(segment, vec![outcome]).unwrap();
            assert!(run_pass_capturing_events(&process, &reversed, Some(&cfg), None).is_err());
            let forward = reversed.with_allele_orientation(segment, SegmentOrientation::Forward);
            let selected = run_pass_capturing_events(&process, &forward, Some(&cfg), None).unwrap();
            let changed = selected
                .sim
                .with_allele_orientation(segment, SegmentOrientation::ReverseComplement);
            assert!(run_pass_capturing_events(
                &EmitCodingEndPPass::new(end),
                &changed,
                Some(&cfg),
                None
            )
            .is_err());
        }
    }
}
