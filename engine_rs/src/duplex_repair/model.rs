//! Finite authored repair histories. Binding/preflight is a compile-time task;
//! rows remain distinct even when their completed nucleotide strings agree.
use super::{
    bindings::DuplexReferenceBindings, repair_coupled_annotated, seed::seed_from_coupled,
    JoinScript, COUPLED_MODEL,
};
use crate::ir::{DuplexHistorySelection, Segment, Simulation};
use crate::refdata::{AlleleId, RefDataConfig};
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::collections::HashSet;
use std::sync::Arc;

pub const FINITE_DUPLEX_MODEL: &str = "finite_authored_duplex_history_v1";
pub const FINITE_DUPLEX_LAW: &str = "finite_authored_history_local_conditioning_v1";

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct CoupledDuplexScript {
    pub model: String,
    pub v: String,
    pub d: String,
    pub j: String,
    pub d_inverted: bool,
    pub joins: Vec<JoinScript>,
    pub max_operations: usize,
    pub max_bases: usize,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct DuplexHistoryRow {
    pub id: String,
    pub weight: f64,
    pub v_call: String,
    pub d_call: String,
    pub j_call: String,
    pub script: CoupledDuplexScript,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct FiniteDuplexModel {
    pub model: String,
    pub reference_content_hash: String,
    pub reference_ends_are_coding_ends: bool,
    pub max_preflight_operations: u64,
    /// Aggregate reservation of active positive rows' declared max_bases.
    pub max_preflight_bases: u64,
    pub rows: Vec<DuplexHistoryRow>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct DuplexContextKey {
    pub v: AlleleId,
    pub d: AlleleId,
    pub j: AlleleId,
    pub d_inverted: bool,
}

impl DuplexContextKey {
    pub fn from_simulation(sim: &Simulation) -> Option<Self> {
        Some(Self {
            v: sim.assignments.get(Segment::V)?.allele_id,
            d: sim.assignments.get(Segment::D)?.allele_id,
            j: sim.assignments.get(Segment::J)?.allele_id,
            d_inverted: sim.assignments.get(Segment::D)?.orientation.is_reverse(),
        })
    }
    fn ordered(self) -> (u32, u32, u32, bool) {
        (
            self.v.index(),
            self.d.index(),
            self.j.index(),
            self.d_inverted,
        )
    }
}
impl PartialOrd for DuplexContextKey {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for DuplexContextKey {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.ordered().cmp(&other.ordered())
    }
}

#[derive(Clone, Debug)]
pub struct PreparedDuplexHistory {
    pub row_index: usize,
    pub context: DuplexContextKey,
    pub seed: Simulation,
}

impl FiniteDuplexModel {
    /// Refuse a floating-point support that this version cannot represent.
    /// In particular, a positive authored history must never become zero mass
    /// as a side effect of overflow-safe normalization.
    pub fn validate_weight_support(&self, indices: &[usize]) -> Result<(), String> {
        if indices.is_empty() {
            return Err("finite duplex weight support is empty".into());
        }
        let weights = indices
            .iter()
            .map(|index| {
                self.rows
                    .get(*index)
                    .map(|row| row.weight)
                    .filter(|weight| weight.is_finite() && *weight > 0.0)
                    .ok_or("finite duplex weight support contains an invalid or nonpositive row")
            })
            .collect::<Result<Vec<_>, _>>()?;
        let max = weights.iter().copied().fold(0.0f64, f64::max);
        let scaled: Vec<_> = weights.iter().map(|weight| weight / max).collect();
        if scaled
            .iter()
            .any(|weight| !weight.is_finite() || *weight <= 0.0)
            || !scaled.iter().sum::<f64>().is_finite()
        {
            return Err(
                "finite duplex positive weights exceed representable numerical support".into(),
            );
        }
        let mut cumulative = 0.0;
        for weight in scaled {
            let next = cumulative + weight;
            if next <= cumulative {
                return Err(
                    "finite duplex positive weight is lost at the authored-order CDF numerical resolution".into(),
                );
            }
            cumulative = next;
        }
        Ok(())
    }

    pub fn validate_schema(&self) -> Result<(), String> {
        if self.model != FINITE_DUPLEX_MODEL {
            return Err("unsupported finite duplex model version".into());
        }
        if self.rows.is_empty() || self.rows.len() > 100_000 {
            return Err("finite duplex model requires 1..100000 authored rows".into());
        }
        let hash = self
            .reference_content_hash
            .strip_prefix("sha256:")
            .ok_or("finite duplex model requires a bound sha256 reference content hash")?;
        if hash.len() != 64 || !hash.bytes().all(|b| b.is_ascii_hexdigit()) {
            return Err("invalid finite duplex reference content hash".into());
        }
        if self.max_preflight_operations == 0 || self.max_preflight_bases == 0 {
            return Err("finite duplex preflight budgets must be positive".into());
        }
        let mut ids = HashSet::new();
        for (index, row) in self.rows.iter().enumerate() {
            if row.id.is_empty() || !ids.insert(row.id.as_str()) {
                return Err(format!(
                    "finite duplex row {index} has an empty or duplicate history ID"
                ));
            }
            if !row.weight.is_finite() || row.weight < 0.0 {
                return Err(format!(
                    "finite duplex row {} has an invalid natural weight",
                    row.id
                ));
            }
            if row.script.model != COUPLED_MODEL {
                return Err(format!(
                    "finite duplex row {} has an unsupported script version",
                    row.id
                ));
            }
            if [&row.v_call, &row.d_call, &row.j_call]
                .iter()
                .any(|s| s.is_empty())
            {
                return Err(format!(
                    "finite duplex row {} lacks exact source calls",
                    row.id
                ));
            }
        }
        Ok(())
    }

    /// Includes authored order, IDs, weights, operation order, all budgets,
    /// bindings and the named probability/placement/capability contracts.
    pub fn digest(&self) -> String {
        let bytes = serde_json::to_vec(&(
            self,
            FINITE_DUPLEX_LAW,
            super::projection::PROJECTION_MODEL,
            crate::productivity::DUPLEX_POLICY,
            "native_single_owner_seed_v1",
        ))
        .expect("finite model is serializable");
        format!("sha256:{:x}", Sha256::digest(bytes))
    }

    pub fn context_key(
        &self,
        row_index: usize,
        refdata: &RefDataConfig,
    ) -> Result<DuplexContextKey, String> {
        let row = self
            .rows
            .get(row_index)
            .ok_or("finite duplex row index outside model")?;
        let get = |segment, name: &str| {
            refdata
                .pool_for(segment)
                .and_then(|p| p.find_by_name(name))
                .map(|(id, _)| id)
                .ok_or_else(|| format!("finite duplex row {} has unknown source {name}", row.id))
        };
        Ok(DuplexContextKey {
            v: get(Segment::V, &row.v_call)?,
            d: get(Segment::D, &row.d_call)?,
            j: get(Segment::J, &row.j_call)?,
            d_inverted: row.script.d_inverted,
        })
    }

    /// Materialize exactly this row through the shared kernel and importer.
    /// The compiler must call this for all positive active rows BEFORE masks.
    /// Replay also calls it, never drawing an alternate history on failure.
    pub fn prepare_row(
        &self,
        row_index: usize,
        refdata: &RefDataConfig,
    ) -> Result<PreparedDuplexHistory, String> {
        self.prepare_row_with_digest(row_index, refdata, &self.digest())
    }

    /// Compile batches and this model's source pass already hold its immutable
    /// digest. Reuse that digest without serializing the entire table for every
    /// row. Reference binding and physical replay validation remain per-row.
    pub(crate) fn prepare_row_with_digest(
        &self,
        row_index: usize,
        refdata: &RefDataConfig,
        model_digest: &str,
    ) -> Result<PreparedDuplexHistory, String> {
        let row = self
            .rows
            .get(row_index)
            .ok_or("finite duplex row index outside model")?;
        if self.reference_content_hash != crate::trace_file::refdata_content_hash(refdata) {
            return Err("finite duplex model reference content hash mismatch".into());
        }
        let s = &row.script;
        let bindings = DuplexReferenceBindings::new(
            refdata,
            [&row.v_call, &row.d_call, &row.j_call],
            [&s.v, &s.d, &s.j],
            s.d_inverted,
            self.reference_ends_are_coding_ends,
        )
        .map_err(|e| format!("history {}: {e}", row.id))?;
        let annotated = repair_coupled_annotated(
            &s.v,
            &s.d,
            &s.j,
            s.d_inverted,
            &s.joins,
            s.max_operations,
            s.max_bases,
        )
        .map_err(|e| format!("history {}: {e}", row.id))?;
        let mut seed = seed_from_coupled(annotated, bindings, refdata)
            .map_err(|e| format!("history {}: {e}", row.id))?;
        let p = Arc::make_mut(seed.duplex_provenance.as_mut().expect("seed has founder"));
        Arc::make_mut(&mut p.founder).selection = Some(DuplexHistorySelection {
            model_digest: model_digest.to_owned(),
            history_id: row.id.clone(),
            row_index,
            natural_weight_bits: row.weight.to_bits(),
            conditioning_law: FINITE_DUPLEX_LAW.into(),
        });
        Ok(PreparedDuplexHistory {
            row_index,
            context: self.context_key(row_index, refdata)?,
            seed,
        })
    }
}

#[cfg(test)]
pub(crate) mod tests;
