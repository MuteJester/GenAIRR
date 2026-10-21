//! Closed per-attempt receipts. Choices reconstruct states; no state snapshot
//! is accepted from the wire. Seed labels are provenance, never replay input.
use super::branching::BranchingParams;
use super::tree::LineageTree;
use crate::family_replay::{
    outcome_events_digest, semantic_state_digest, F64Bits, FamilyReplayBudget,
};
use crate::ir::Simulation;
use crate::pass::{Outcome, PassError};
use crate::rng::Rng;
use crate::trace::Trace;
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AttemptConfig {
    pub rate: F64Bits,
    pub lambda_base: F64Bits,
    pub lambda_mut: F64Bits,
    pub max_generations: u32,
    pub n_max: u32,
    pub n_sample: u32,
    pub seed_label: String,
    pub selection_strength: F64Bits,
    pub beta: F64Bits,
    pub target_aa: Option<String>,
    pub mature_substitutions: u32,
    pub mutability: Vec<F64Bits>,
    pub substitution: Vec<F64Bits>,
    pub contracts: Option<Vec<String>>,
    pub strict: bool,
}

pub fn error(reason: impl Into<String>) -> PassError {
    PassError::invalid_plan_state("family.replay", reason)
}
pub fn label(value: &str) -> Result<u64, String> {
    let n: u64 = value.parse().map_err(|_| "invalid u64 seed label")?;
    if n.to_string() != value {
        return Err("noncanonical u64 seed label".into());
    }
    Ok(n)
}

impl AttemptConfig {
    pub fn validate(&self) -> Result<(), String> {
        if [
            self.rate,
            self.lambda_base,
            self.lambda_mut,
            self.selection_strength,
            self.beta,
        ]
        .iter()
        .any(|x| !x.get().is_finite())
            || self
                .mutability
                .iter()
                .chain(&self.substitution)
                .any(|x| !x.get().is_finite())
        {
            return Err("nonfinite family capture parameter".into());
        }
        let rate = self.rate.get();
        let lambda = self.lambda_base.get();
        let strength = self.selection_strength.get();
        if !(0.0..=1.0).contains(&rate)
            || lambda < 0.0
            || self.lambda_mut.get() < 0.0
            || strength < 0.0
            || self.beta.get() < 0.0
            || !(lambda * (1.0 + strength)).is_finite()
            || lambda * (1.0 + strength) > 64.0
        {
            return Err("unsupported family capture numerical domain; require lambda_base*(1+selection_strength)<=64".into());
        }
        if self.max_generations == 0
            || self.max_generations > 1000
            || self.n_max == 0
            || self.n_sample == 0
        {
            return Err(
                "family capture requires positive generation/capacity/sample controls".into(),
            );
        }
        label(&self.seed_label)?;
        if self.mutability.len() != 1024
            || self.substitution.len() != 4096
            || self
                .mutability
                .iter()
                .chain(&self.substitution)
                .any(|x| x.get() < 0.0)
        {
            return Err("invalid captured S5F kernel".into());
        }
        if let Some(target) = &self.target_aa {
            if target.is_empty()
                || !target
                    .bytes()
                    .all(|aa| super::affinity::AA_ORDER.contains(&aa))
            {
                return Err(
                    "explicit captured target must contain standard uppercase amino acids".into(),
                );
            }
        }
        self.contract_set()?;
        Ok(())
    }
    pub fn params(&self) -> Result<BranchingParams, String> {
        Ok(BranchingParams {
            lambda_base: self.lambda_base.get(),
            lambda_mut: self.lambda_mut.get(),
            max_generations: self.max_generations,
            n_max: self.n_max,
            n_sample: self.n_sample,
            seed: label(&self.seed_label)?,
        })
    }
    pub fn contract_set(&self) -> Result<Option<crate::contract::ContractSet>, String> {
        use crate::contract::*;
        use crate::ir::Segment;
        self.contracts
            .as_ref()
            .map(|names| {
                let mut set = ContractSet::new();
                for name in names {
                    let c: Box<dyn Contract> = match name.as_str() {
                        "coding_sequence_integrity" => Box::new(CodingSequenceIntegrity),
                        "productive_junction_frame" => Box::new(ProductiveJunctionFrame::new()),
                        "no_stop_codon_in_junction" => Box::new(NoStopCodonInJunction::new()),
                        "canonical_anchor.v1.v" => Box::new(AnchorPreserved::canonical(Segment::V)),
                        "canonical_anchor.v1.j" => Box::new(AnchorPreserved::canonical(Segment::J)),
                        _ => return Err(format!("unsupported captured contract {name}")),
                    };
                    set.add(c);
                }
                Ok(set)
            })
            .transpose()
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TargetChoice {
    pub position: u32,
    pub residue: String,
}
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TargetReceipt {
    pub choices: Vec<TargetChoice>,
    pub resolved_target: Option<String>,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ParentProposal {
    pub generation: u32,
    pub parent_id: u32,
    pub rate: F64Bits,
    pub proposed_offspring: u32,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DivisionReceipt {
    pub node_id: u32,
    pub parent_id: u32,
    pub generation: u32,
    pub seed_label: String,
    pub trace: Trace,
    pub parent_state_digest: String,
    pub child_state_digest: String,
    pub events_digest: String,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SampleChoice {
    pub live_index: u32,
    pub node_id: u32,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct NodeSummary {
    pub node_id: u32,
    pub parent_id: Option<u32>,
    pub generation: u32,
    pub mutations_from_parent: u32,
    pub affinity: F64Bits,
    pub abundance: u32,
    pub observed: bool,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FamilyAttemptReceipt {
    pub schema_version: u32,
    pub model: String,
    pub config: AttemptConfig,
    pub founder_state_digest: String,
    pub founder_events_digest: String,
    pub target: TargetReceipt,
    pub growth: Vec<ParentProposal>,
    pub edges: Vec<DivisionReceipt>,
    pub final_live: Vec<u32>,
    pub samples: Vec<SampleChoice>,
    pub nodes: Vec<NodeSummary>,
    pub completed_generations: u32,
    pub peak_live: usize,
    pub termination: String,
}

pub struct FamilyReplayDriver<'a> {
    pub receipt: FamilyAttemptReceipt,
    expected: Option<&'a FamilyAttemptReceipt>,
    pub budget: &'a mut FamilyReplayBudget,
    growth_cursor: usize,
    edge_cursor: usize,
    sample_cursor: usize,
    target_cursor: usize,
}
impl<'a> FamilyReplayDriver<'a> {
    pub fn new(
        config: AttemptConfig,
        founder: &Outcome,
        expected: Option<&'a FamilyAttemptReceipt>,
        budget: &'a mut FamilyReplayBudget,
    ) -> Result<Self, PassError> {
        config.validate().map_err(error)?;
        budget
            .charge_history(founder.final_simulation())
            .map_err(error)?;
        let receipt = FamilyAttemptReceipt {
            schema_version: 1,
            model: "duplex_lineage_attempt_v1".into(),
            config,
            founder_state_digest: semantic_state_digest(founder.final_simulation())
                .map_err(error)?,
            founder_events_digest: outcome_events_digest(founder).map_err(error)?,
            target: TargetReceipt::default(),
            growth: vec![],
            edges: vec![],
            final_live: vec![],
            samples: vec![],
            nodes: vec![],
            completed_generations: 0,
            peak_live: 1,
            termination: String::new(),
        };
        if let Some(e) = expected {
            if e.schema_version != 1
                || e.model != receipt.model
                || serde_json::to_value(&e.config).map_err(|e| error(e.to_string()))?
                    != serde_json::to_value(&receipt.config).map_err(|e| error(e.to_string()))?
                || e.founder_state_digest != receipt.founder_state_digest
                || e.founder_events_digest != receipt.founder_events_digest
            {
                return Err(error("attempt founder/configuration/version mismatch"));
            }
        }
        Ok(Self {
            receipt,
            expected,
            budget,
            growth_cursor: 0,
            edge_cursor: 0,
            sample_cursor: 0,
            target_cursor: 0,
        })
    }
    pub fn is_replay(&self) -> bool {
        self.expected.is_some()
    }
    pub fn work(&mut self, n: u64) -> Result<(), PassError> {
        self.budget.charge_work(n).map_err(error)
    }
    pub fn target_choice(&mut self, target: &[u8], rng: &mut Rng) -> Result<(usize, u8), String> {
        self.budget.charge_work(1)?;
        self.budget.charge_choices(1)?;
        let (position, residue) = if let Some(e) = self.expected {
            let c = e
                .target
                .choices
                .get(self.target_cursor)
                .ok_or("missing target choice")?;
            let bytes = c.residue.as_bytes();
            if bytes.len() != 1 {
                return Err("target residue must be one ASCII amino acid".into());
            }
            (c.position as usize, bytes[0])
        } else {
            let before = rng.words_consumed();
            let p = rng.range_u32(target.len() as u32) as usize;
            let cur = target[p];
            let aa = match super::affinity::AA_ORDER.iter().position(|x| *x == cur) {
                Some(index) => {
                    let pick = rng.range_u32(19) as usize;
                    super::affinity::AA_ORDER[if pick >= index { pick + 1 } else { pick }]
                }
                None => super::affinity::AA_ORDER[rng.range_u32(20) as usize],
            };
            self.budget.charge_work(rng.words_consumed() - before)?;
            (p, aa)
        };
        self.receipt.target.choices.push(TargetChoice {
            position: u32::try_from(position).map_err(|_| "target position overflow")?,
            residue: String::from_utf8(vec![residue]).map_err(|_| "non-ASCII target")?,
        });
        self.target_cursor += 1;
        Ok((position, residue))
    }
    pub fn resolved_target(&mut self, target: Option<&[u8]>) -> Result<(), PassError> {
        self.receipt.target.resolved_target =
            target.map(|x| String::from_utf8_lossy(x).into_owned());
        if let Some(e) = self.expected {
            if self.target_cursor != e.target.choices.len()
                || self.receipt.target.resolved_target != e.target.resolved_target
            {
                return Err(error("target realization mismatch or surplus choices"));
            }
        }
        Ok(())
    }
    pub fn offspring(
        &mut self,
        generation: u32,
        parent_id: u32,
        rate: f64,
        rng: &mut Rng,
    ) -> Result<u32, PassError> {
        self.work(1)?;
        self.budget.charge_choices(1).map_err(error)?;
        if !rate.is_finite() || !(0.0..=64.0).contains(&rate) {
            return Err(error("invalid effective captured offspring rate"));
        }
        let proposed_offspring = if let Some(e) = self.expected {
            let p = e
                .growth
                .get(self.growth_cursor)
                .ok_or_else(|| error("missing visited-parent proposal"))?;
            if p.generation != generation
                || p.parent_id != parent_id
                || p.rate.get().to_bits() != rate.to_bits()
                || (rate == 0.0 && p.proposed_offspring != 0)
            {
                return Err(error(format!("offspring scope, numerical rate or support mismatch at generation {generation} parent {parent_id}: recorded ({},{},{:016x}), actual ({generation},{parent_id},{:016x})",p.generation,p.parent_id,p.rate.get().to_bits(),rate.to_bits())));
            }
            p.proposed_offspring
        } else {
            super::poisson::poisson_sample_with_work(rng, rate, || self.budget.charge_work(1))
                .map_err(error)?
        };
        self.receipt.growth.push(ParentProposal {
            generation,
            parent_id,
            rate: F64Bits::new(rate).map_err(error)?,
            proposed_offspring,
        });
        self.growth_cursor += 1;
        Ok(proposed_offspring)
    }
    pub fn begin_edge(
        &mut self,
        node_id: u32,
        parent_id: u32,
        generation: u32,
        parent: &Simulation,
        rng: &mut Rng,
    ) -> Result<(u64, Option<Trace>), PassError> {
        self.budget.charge_nodes(1).map_err(error)?;
        self.budget
            .charge_bases(parent.pool.len() as u64)
            .map_err(error)?;
        self.budget.charge_history(parent).map_err(error)?;
        self.work(1)?;
        self.budget.charge_choices(1).map_err(error)?;
        if let Some(e) = self.expected {
            let edge = e
                .edges
                .get(self.edge_cursor)
                .ok_or_else(|| error("missing division receipt"))?;
            if (edge.node_id, edge.parent_id, edge.generation) != (node_id, parent_id, generation)
                || semantic_state_digest(parent).map_err(error)? != edge.parent_state_digest
            {
                return Err(error("division parent identity mismatch"));
            }
            self.budget
                .charge_choices(edge.trace.len() as u64)
                .map_err(error)?;
            Ok((
                label(&edge.seed_label).map_err(error)?,
                Some(edge.trace.clone()),
            ))
        } else {
            self.work(1)?;
            Ok((rng.next_u64(), None))
        }
    }
    pub fn commit_edge(
        &mut self,
        node_id: u32,
        parent_id: u32,
        generation: u32,
        seed: u64,
        edge: &Outcome,
    ) -> Result<(), PassError> {
        let child = edge.final_simulation();
        self.budget
            .charge_bases(child.pool.len() as u64)
            .map_err(error)?;
        self.budget.charge_history(child).map_err(error)?;
        if self.expected.is_none() {
            self.budget
                .charge_choices(edge.trace.len() as u64)
                .map_err(error)?;
        }
        let record = DivisionReceipt {
            node_id,
            parent_id,
            generation,
            seed_label: seed.to_string(),
            trace: edge.trace.clone(),
            parent_state_digest: semantic_state_digest(&edge.revisions[0]).map_err(error)?,
            child_state_digest: semantic_state_digest(child).map_err(error)?,
            events_digest: outcome_events_digest(edge).map_err(error)?,
        };
        if let Some(e) = self.expected {
            let expected = &e.edges[self.edge_cursor];
            if record.child_state_digest != expected.child_state_digest
                || record.events_digest != expected.events_digest
            {
                return Err(error("division state/event commitment mismatch"));
            }
        }
        self.receipt.edges.push(record);
        self.edge_cursor += 1;
        Ok(())
    }
    pub fn finish_growth(
        &mut self,
        live: &[u32],
        completed: u32,
        peak: usize,
    ) -> Result<(), PassError> {
        self.work(live.len() as u64 + 1)?;
        self.receipt.final_live = live.to_vec();
        self.receipt.completed_generations = completed;
        self.receipt.peak_live = peak;
        self.receipt.termination = if live.is_empty() {
            "extinct"
        } else {
            "max_generations"
        }
        .into();
        if let Some(e) = self.expected {
            if self.growth_cursor != e.growth.len()
                || self.edge_cursor != e.edges.len()
                || self.receipt.final_live != e.final_live
                || completed != e.completed_generations
                || peak != e.peak_live
                || self.receipt.termination != e.termination
            {
                return Err(error(
                    "growth termination, topology or receipt consumption mismatch",
                ));
            }
        }
        Ok(())
    }
    pub fn sample(&mut self, live: &[u32], rng: &mut Rng) -> Result<usize, String> {
        self.budget.charge_work(1)?;
        self.budget.charge_choices(1)?;
        let index = if let Some(e) = self.expected {
            let c = e
                .samples
                .get(self.sample_cursor)
                .ok_or("missing sample choice")?;
            let index = c.live_index as usize;
            if live.get(index) != Some(&c.node_id) {
                return Err("sample scope or eligibility mismatch".into());
            }
            index
        } else {
            let before = rng.words_consumed();
            let index = rng.range_u32(live.len() as u32) as usize;
            self.budget.charge_work(rng.words_consumed() - before)?;
            index
        };
        self.receipt.samples.push(SampleChoice {
            live_index: index as u32,
            node_id: live[index],
        });
        self.sample_cursor += 1;
        Ok(index)
    }
    pub fn finish(&mut self, tree: &LineageTree) -> Result<(), PassError> {
        self.work(tree.len() as u64)?;
        self.receipt.nodes = tree
            .nodes
            .iter()
            .map(|node| {
                Ok(NodeSummary {
                    node_id: node.id,
                    parent_id: node.parent_id,
                    generation: node.generation,
                    mutations_from_parent: node.mutations_from_parent,
                    affinity: F64Bits::new(node.affinity)?,
                    abundance: node.abundance,
                    observed: node.observed,
                })
            })
            .collect::<Result<_, String>>()
            .map_err(error)?;
        if let Some(e) = self.expected {
            if self.sample_cursor != e.samples.len()
                || serde_json::to_value(&self.receipt).map_err(|e| error(e.to_string()))?
                    != serde_json::to_value(e).map_err(|e| error(e.to_string()))?
            {
                return Err(error(
                    "attempt summary, samples or complete receipt mismatch",
                ));
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests;
