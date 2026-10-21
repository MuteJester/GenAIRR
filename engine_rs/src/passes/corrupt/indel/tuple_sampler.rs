//! Natural sequential indel histories conditioned on the complete final
//! contract. Complete independent proposals precede a bounded exact search,
//! optimized for counts <= 2. Exhaustion never means that biological support
//! is empty.
use crate::address;
use crate::contract::{ChoiceContext, ContractKind, ContractSet};
use crate::ir::{flag, NucHandle, Nucleotide, Segment, Simulation};
use crate::pass::{Pass, PassContext, PassError};
use crate::productivity::SequenceView;
use crate::refdata::RefDataConfig;
use crate::rng::Rng;

use super::edit_view::{EditIndex, EditView};
use super::event::IndelEvent;
use super::IndelPass;

pub(super) const SAMPLER_MODEL: &str = "natural_sequential_v1";

#[derive(Clone, Copy, Debug)]
pub struct IndelSamplingOptions {
    /// Complete independent proposals before the exact fallback.
    pub rejection_attempts: u32,
    /// Maximum fallback work nodes, including exhaustive nodes and restart-probe edits.
    pub exact_node_budget: u64,
    /// Explicit compatibility for historical all-NoOp permissive traces.
    pub legacy_replay_sentinels: bool,
}

impl Default for IndelSamplingOptions {
    fn default() -> Self {
        Self {
            rejection_attempts: 256,
            exact_node_budget: 20_000_000,
            legacy_replay_sentinels: false,
        }
    }
}

#[derive(Clone, Debug)]
pub(super) struct SamplingReport {
    pub strategy: &'static str,
    pub attempts: u64,
    pub explored_states: u64,
}

pub(super) struct TupleSample {
    /// None is used only after a complete proof of empty support.
    pub events: Option<Vec<IndelEvent>>,
    pub report: SamplingReport,
}

#[derive(Clone, Copy)]
struct WeightedEvent {
    event: IndelEvent,
    log_weight: f64,
}

/// Stable weighted reservoir. Reachability is represented by `selected`, not
/// by a floating-point sum that could underflow to zero.
struct Reservoir {
    log_mass: f64,
    selected: Option<Vec<IndelEvent>>,
}
impl Reservoir {
    fn new() -> Self {
        Self {
            log_mass: f64::NEG_INFINITY,
            selected: None,
        }
    }
    fn add(&mut self, events: &[IndelEvent], log_weight: f64, rng: &mut Rng) {
        debug_assert!(log_weight.is_finite());
        if self.selected.is_none() {
            self.log_mass = log_weight;
            self.selected = Some(events.to_vec());
            return;
        }
        let high = self.log_mass.max(log_weight);
        let low = self.log_mass.min(log_weight);
        let total = high + (low - high).exp().ln_1p();
        if rng.next_f64() < (log_weight - total).exp() {
            self.selected = Some(events.to_vec());
        }
        self.log_mass = total;
    }
}

struct Work {
    count: u32,
    attempts: u64,
    states: u64,
    limit: u64,
    /// Count only exhaustive nodes for the deterministic restart schedule.
    search_nodes: u64,
    restart_probes: bool,
    accepted_probe: Option<Vec<IndelEvent>>,
}
impl Work {
    fn visit(&mut self) -> Result<(), PassError> {
        if self.states >= self.limit {
            return Err(PassError::sampling_budget_exceeded(
                address::CORRUPT_INDEL,
                address::ChoiceAddress::CorruptIndelCount.to_string(),
                self.count,
                self.attempts,
                self.states,
                self.limit,
                "natural_sequential_exact_search",
            ));
        }
        self.states += 1;
        Ok(())
    }
}

pub(super) fn tuple_admitted(
    initial: &Simulation,
    final_state: &Simulation,
    contracts: &ContractSet,
    refdata: Option<&RefDataConfig>,
    count: u32,
) -> bool {
    contracts
        .admits_post_event(
            initial,
            final_state,
            refdata,
            ChoiceContext::indel_insertion(count.saturating_sub(1), count, NucHandle::new(0)),
        )
        .is_ok()
        && (!contracts.requires_productive_emission()
            || refdata.is_some_and(|cfg| {
                crate::productivity::assess(final_state, cfg).productive == Some(true)
            }))
}

fn standard_productive_bundle(contracts: &ContractSet) -> bool {
    let kinds = [
        ContractKind::ProductiveJunctionFrame,
        ContractKind::NoStopCodonInJunction,
        ContractKind::CodingSequenceIntegrity,
        ContractKind::CanonicalAnchor {
            segment: Segment::V,
        },
        ContractKind::CanonicalAnchor {
            segment: Segment::J,
        },
    ];
    contracts.len() == kinds.len() && kinds.into_iter().all(|kind| contracts.contains_kind(kind))
}

impl IndelPass {
    fn proposal_with_view<const CAPACITY: usize>(
        &self,
        index: &EditIndex<'_>,
        count: u32,
        bases: &[(u8, f64)],
        rng: &mut Rng,
        cfg: &RefDataConfig,
    ) -> Option<Vec<IndelEvent>> {
        let mut state = index.view_with_capacity::<CAPACITY>();
        let mut events = Vec::with_capacity(count as usize);
        for _ in 0..count {
            let event = self.natural_proposal(state.len(), bases, rng);
            state = state.edited(event);
            events.push(event);
        }
        (crate::productivity::assess(&state, cfg).productive == Some(true)).then_some(events)
    }

    /// Each proposal starts at the unchanged root and consumes exactly the
    /// same natural draws, independent of the representation used to assess it.
    #[allow(clippy::too_many_arguments)]
    fn independent_proposal(
        &self,
        initial: &Simulation,
        index: Option<&EditIndex<'_>>,
        contracts: &ContractSet,
        refdata: Option<&RefDataConfig>,
        count: u32,
        bases: &[(u8, f64)],
        rng: &mut Rng,
    ) -> Option<Vec<IndelEvent>> {
        if let Some(index) = index {
            if count <= 2 {
                return self.proposal_with_view::<5>(index, count, bases, rng, refdata.unwrap());
            }
            if count <= 64 {
                return self.proposal_with_view::<129>(index, count, bases, rng, refdata.unwrap());
            }
        }
        let mut state = initial.clone();
        let mut events = Vec::with_capacity(count as usize);
        for _ in 0..count {
            let event = self.natural_proposal(state.pool.len(), bases, rng);
            state = apply_event_hypothetical(&state, event);
            events.push(event);
        }
        tuple_admitted(initial, &state, contracts, refdata, count).then_some(events)
    }

    /// Draw the same normalized natural law that exact enumeration uses.
    /// In particular, do not delegate to a categorical sampler whose finite
    /// raw weights may overflow when summed.
    fn natural_proposal(&self, length: usize, bases: &[(u8, f64)], rng: &mut Rng) -> IndelEvent {
        if rng.next_f64() < self.insertion_prob {
            let site = if length == 0 {
                0
            } else {
                rng.range_u32(length as u32 + 1)
            };
            let draw = rng.next_f64();
            let mut cumulative = 0.0;
            let mut base = bases.last().expect("positive insertion support checked").0;
            for &(candidate, log_weight) in bases {
                cumulative += log_weight.exp();
                if draw < cumulative {
                    base = candidate;
                    break;
                }
            }
            IndelEvent::Insertion { site, base }
        } else {
            IndelEvent::Deletion {
                site: if length == 0 {
                    None
                } else {
                    Some(rng.range_u32(length as u32))
                },
            }
        }
    }

    fn natural_events(&self, length: usize, bases: &[(u8, f64)]) -> Vec<WeightedEvent> {
        let mut events = Vec::new();
        if self.insertion_prob > 0.0 {
            let log_site = self.insertion_prob.ln() - ((length + 1) as f64).ln();
            for site in 0..=length {
                for &(base, log_base) in bases {
                    events.push(WeightedEvent {
                        event: IndelEvent::Insertion {
                            site: site as u32,
                            base,
                        },
                        log_weight: log_site + log_base,
                    });
                }
            }
        }
        if self.insertion_prob < 1.0 {
            let log_kind = (-self.insertion_prob).ln_1p();
            if length == 0 {
                events.push(WeightedEvent {
                    event: IndelEvent::Deletion { site: None },
                    log_weight: log_kind,
                });
            } else {
                let log_site = log_kind - (length as f64).ln();
                for site in 0..length {
                    events.push(WeightedEvent {
                        event: IndelEvent::Deletion {
                            site: Some(site as u32),
                        },
                        log_weight: log_site,
                    });
                }
            }
        }
        events
    }

    fn view_admitted(
        &self,
        initial: &Simulation,
        view: &EditView<'_>,
        events: &[IndelEvent],
        contracts: &ContractSet,
        refdata: Option<&RefDataConfig>,
        fast: bool,
    ) -> bool {
        if fast {
            return crate::productivity::assess(view, refdata.unwrap()).productive == Some(true);
        }
        let mut state = initial.clone();
        for &event in events {
            state = apply_event_hypothetical(&state, event);
        }
        tuple_admitted(initial, &state, contracts, refdata, events.len() as u32)
    }

    fn exact_short(
        &self,
        initial: &Simulation,
        contracts: &ContractSet,
        refdata: Option<&RefDataConfig>,
        count: u32,
        bases: &[(u8, f64)],
        rng: &mut Rng,
        work: &mut Work,
    ) -> Result<Reservoir, PassError> {
        let index = EditIndex::new(initial);
        let view = index.view();
        // EditView currently indexes the legacy nucleotide cache. A duplex
        // molecule instead needs its persistent certified-placement map,
        // updated by each materialized hypothetical event.
        let fast = initial.duplex_provenance.is_none()
            && refdata.is_some()
            && standard_productive_bundle(contracts);
        let mut reservoir = Reservoir::new();
        if count == 0 {
            work.visit()?;
            if tuple_admitted(initial, initial, contracts, refdata, 0) {
                reservoir.add(&[], 0.0, rng);
            }
            return Ok(reservoir);
        }
        for first in self.natural_events(view.len(), bases) {
            work.visit()?;
            let after_first = view.edited(first.event);
            if count == 1 {
                if self.view_admitted(
                    initial,
                    &after_first,
                    &[first.event],
                    contracts,
                    refdata,
                    fast,
                ) {
                    reservoir.add(&[first.event], first.log_weight, rng);
                }
                continue;
            }
            for second in self.natural_events(after_first.len(), bases) {
                work.visit()?;
                let events = [first.event, second.event];
                let after_second = after_first.edited(second.event);
                if self.view_admitted(initial, &after_second, &events, contracts, refdata, fast) {
                    reservoir.add(&events, first.log_weight + second.log_weight, rng);
                }
            }
        }
        Ok(reservoir)
    }

    #[allow(clippy::too_many_arguments)]
    fn exact_deeper(
        &self,
        initial: &Simulation,
        state: &Simulation,
        contracts: &ContractSet,
        refdata: Option<&RefDataConfig>,
        remaining: u32,
        bases: &[(u8, f64)],
        events: &mut Vec<IndelEvent>,
        log_weight: f64,
        rng: &mut Rng,
        work: &mut Work,
        reservoir: &mut Reservoir,
        index: Option<&EditIndex<'_>>,
    ) -> Result<(), PassError> {
        work.visit()?;
        work.search_nodes += 1;
        // Exhaustive traversal has enormous branching for realistic receptors.
        // Independent restart probes can finish promptly without sampling the
        // incomplete reservoir. Their schedule depends only on visited search
        // nodes, never on which histories entered that reservoir.
        if work.restart_probes
            && work.search_nodes % 16 == 0
            && work.limit - work.states >= u64::from(work.count)
        {
            for _ in 0..work.count {
                work.visit()?;
            }
            work.attempts += 1;
            work.accepted_probe = self
                .independent_proposal(initial, index, contracts, refdata, work.count, bases, rng);
            if work.accepted_probe.is_some() {
                return Ok(());
            }
        }
        if remaining == 0 {
            if tuple_admitted(initial, state, contracts, refdata, work.count) {
                reservoir.add(events, log_weight, rng);
            }
            return Ok(());
        }
        for candidate in self.natural_events(state.pool.len(), bases) {
            let next = apply_event_hypothetical(state, candidate.event);
            events.push(candidate.event);
            self.exact_deeper(
                initial,
                &next,
                contracts,
                refdata,
                remaining - 1,
                bases,
                events,
                log_weight + candidate.log_weight,
                rng,
                work,
                reservoir,
                index,
            )?;
            events.pop();
            if work.accepted_probe.is_some() {
                return Ok(());
            }
        }
        Ok(())
    }

    pub(super) fn sample_admissible_tuple(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
        count: u32,
    ) -> Result<TupleSample, PassError> {
        let contracts = ctx
            .contracts
            .expect("constrained sampler requires contracts");
        if count > 4096 {
            return Err(PassError::sampling_budget_exceeded(
                self.name(),
                address::ChoiceAddress::CorruptIndelCount.to_string(),
                count,
                0,
                0,
                self.sampling_options.exact_node_budget,
                "tuple_length_limit_4096",
            ));
        }
        let bases = if self.insertion_prob > 0.0 {
            self.base_log_support().map_err(|reason| {
                PassError::constraint_sampling(
                    self.name(),
                    address::ChoiceAddress::CorruptIndelCount.to_string(),
                    reason,
                )
            })?
        } else {
            Vec::new()
        };
        let mut work = Work {
            count,
            attempts: 0,
            states: 0,
            limit: self.sampling_options.exact_node_budget,
            search_nodes: 0,
            restart_probes: self.sampling_options.rejection_attempts > 0,
            accepted_probe: None,
        };
        let index = (sim.duplex_provenance.is_none()
            && ctx.refdata.is_some()
            && standard_productive_bundle(contracts))
        .then(|| EditIndex::new(sim));
        if count > 0 {
            // Each attempt is complete and starts at the same input. A failed
            // candidate does not alter the root or redraw the requested count.
            for _ in 0..self.sampling_options.rejection_attempts {
                work.attempts += 1;
                if let Some(events) = self.independent_proposal(
                    sim,
                    index.as_ref(),
                    contracts,
                    ctx.refdata,
                    count,
                    &bases,
                    ctx.rng,
                ) {
                    return Ok(TupleSample {
                        events: Some(events),
                        report: SamplingReport {
                            strategy: "independent_tuple_rejection",
                            attempts: work.attempts,
                            explored_states: 0,
                        },
                    });
                }
            }
        }
        let reservoir = if count <= 2 {
            self.exact_short(
                sim,
                contracts,
                ctx.refdata,
                count,
                &bases,
                ctx.rng,
                &mut work,
            )?
        } else {
            // The stack depth is itself a deterministic resource limit. No
            // partial sample is returned if exhaustive completion is too deep.
            if count > 64 {
                return Err(PassError::sampling_budget_exceeded(
                    self.name(),
                    address::ChoiceAddress::CorruptIndelCount.to_string(),
                    count,
                    work.attempts,
                    0,
                    work.limit,
                    "exact_search_depth_limit",
                ));
            }
            let mut reservoir = Reservoir::new();
            self.exact_deeper(
                sim,
                sim,
                contracts,
                ctx.refdata,
                count,
                &bases,
                &mut Vec::new(),
                0.0,
                ctx.rng,
                &mut work,
                &mut reservoir,
                index.as_ref(),
            )?;
            if let Some(events) = work.accepted_probe.take() {
                return Ok(TupleSample {
                    events: Some(events),
                    report: SamplingReport {
                        strategy: "exact_search_independent_probe",
                        attempts: work.attempts,
                        explored_states: work.states,
                    },
                });
            }
            reservoir
        };
        Ok(TupleSample {
            events: reservoir.selected,
            report: SamplingReport {
                strategy: if count <= 2 {
                    "exact_short_tuple"
                } else {
                    "exact_fallback"
                },
                attempts: work.attempts,
                explored_states: work.states,
            },
        })
    }

    /// Replay validation for the historical all-no-op encoding in new-model
    /// traces: permit it only after a complete proof of empty support.
    pub(super) fn prove_empty_support(
        &self,
        sim: &Simulation,
        contracts: &ContractSet,
        refdata: Option<&RefDataConfig>,
        count: u32,
    ) -> Result<(bool, u64), PassError> {
        let bases = if self.insertion_prob > 0.0 {
            self.base_log_support().map_err(|reason| {
                PassError::constraint_sampling(
                    self.name(),
                    address::ChoiceAddress::CorruptIndelCount.to_string(),
                    reason,
                )
            })?
        } else {
            Vec::new()
        };
        let mut work = Work {
            count,
            attempts: 0,
            states: 0,
            limit: self.sampling_options.exact_node_budget,
            search_nodes: 0,
            restart_probes: false,
            accepted_probe: None,
        };
        let mut rng = Rng::new(0);
        let reservoir = if count <= 2 {
            self.exact_short(sim, contracts, refdata, count, &bases, &mut rng, &mut work)?
        } else {
            if count > 64 {
                return Err(PassError::sampling_budget_exceeded(
                    self.name(),
                    address::ChoiceAddress::CorruptIndelCount.to_string(),
                    count,
                    0,
                    0,
                    work.limit,
                    "exact_search_depth_limit",
                ));
            }
            let mut reservoir = Reservoir::new();
            self.exact_deeper(
                sim,
                sim,
                contracts,
                refdata,
                count,
                &bases,
                &mut Vec::new(),
                0.0,
                &mut rng,
                &mut work,
                &mut reservoir,
                None,
            )?;
            reservoir
        };
        Ok((reservoir.selected.is_none(), work.states))
    }
}

pub(super) fn apply_event_hypothetical(sim: &Simulation, event: IndelEvent) -> Simulation {
    match event {
        IndelEvent::Insertion { site, base } => {
            let segment = IndelPass::insertion_segment(sim, site);
            sim.with_indel_inserted(
                site,
                Nucleotide::synthetic(base, segment, flag::INDEL_INSERTED),
            )
        }
        IndelEvent::Deletion { site: Some(site) } => sim.with_indel_deleted(site),
        IndelEvent::Deletion { site: None } | IndelEvent::NoOp => sim.clone(),
    }
}
