//! Closed, reconstructible native execution descriptors. Signatures are not codecs.
use super::{
    reference::{ReferenceSpecV1, SegmentSpec},
    F64Bits,
};
use crate::passes::corrupt::indel::IndelSamplingOptions;
use crate::{
    dist::*,
    pass::{Pass, PassPlan},
    passes::*,
    refdata::{AlleleId, RefDataConfig},
    s5f::S5FKernel,
};
use serde::{Deserialize, Serialize};
use std::sync::Arc;

fn finite(v: F64Bits) -> Result<f64, String> {
    let n = v.get();
    if n.is_finite() {
        Ok(n)
    } else {
        Err("nonfinite binary64 parameter".into())
    }
}
fn prob(v: F64Bits) -> Result<f64, String> {
    let n = finite(v)?;
    if (0.0..=1.0).contains(&n) {
        Ok(n)
    } else {
        Err("probability/rate outside [0,1]".into())
    }
}
pub fn bits(v: f64) -> Option<F64Bits> {
    F64Bits::new(v).ok()
}
pub fn vec_bits(v: &[f64]) -> Option<Vec<F64Bits>> {
    v.iter().map(|v| bits(*v)).collect()
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum DistributionSpec {
    UniformBase,
    CategoricalBase {
        weights: [F64Bits; 4],
    },
    UniformInt {
        min: i64,
        max: i64,
    },
    EmpiricalInt {
        values: Vec<i64>,
        cumulative: Vec<F64Bits>,
        total: F64Bits,
    },
    AllelePool {
        ids: Option<Vec<u32>>,
        cumulative: Vec<F64Bits>,
        total: F64Bits,
    },
}
pub fn cdf(values: &[F64Bits], total: F64Bits) -> Result<Vec<f64>, String> {
    crate::family_replay::primitive_work(values.len() as u64);
    let values = values
        .iter()
        .map(|v| finite(*v))
        .collect::<Result<Vec<_>, _>>()?;
    let t = finite(total)?;
    if values.is_empty()
        || values[0] < 0.0
        || values.windows(2).any(|p| p[1] < p[0])
        || t <= 0.0
        || values.last().copied() != Some(t)
    {
        return Err("invalid exact categorical CDF".into());
    }
    Ok(values)
}
impl DistributionSpec {
    pub fn supports_integer(&self, value: i64) -> Result<bool, String> {
        match self {
            Self::UniformInt { min, max } => Ok(*min <= value && value < *max),
            Self::EmpiricalInt {
                values,
                cumulative,
                total,
            } => {
                let c = cdf(cumulative, *total)?;
                if values.len() != c.len() {
                    return Err("integer CDF shape mismatch".into());
                }
                let mut previous = 0.0;
                Ok(values.iter().zip(c).any(|(v, end)| {
                    let positive = end > previous;
                    previous = end;
                    *v == value && positive
                }))
            }
            _ => Err("expected integer distribution".into()),
        }
    }
    /// Positive preimage of the existing end-loss clamp/cast transformation.
    pub fn end_loss_preimage(&self, n: u32, pool_len: u32, greater: bool) -> Result<bool, String> {
        if n > pool_len || (greater && n == pool_len) {
            return Ok(false);
        }
        let qualifies = |raw: i64| {
            let target = (raw.max(0) as u32).min(pool_len);
            if greater {
                target > n
            } else {
                target == n
            }
        };
        match self {
            Self::EmpiricalInt {
                values,
                cumulative,
                total,
            } => {
                let c = cdf(cumulative, *total)?;
                if values.len() != c.len() {
                    return Err("integer CDF shape mismatch".into());
                }
                let mut previous = 0.0;
                Ok(values.iter().zip(c).any(|(raw, end)| {
                    let positive = end > previous;
                    previous = end;
                    positive && qualifies(*raw)
                }))
            }
            Self::UniformInt { min, max } => {
                if max <= min {
                    return Err("empty integer domain".into());
                }
                if !greater && n == 0 && *min <= 0 {
                    return Ok(true);
                }
                let lo = (*min as i128).max(1);
                let hi = *max as i128;
                if lo >= hi {
                    return Ok(false);
                }
                let threshold = if greater { n as i128 + 1 } else { n as i128 };
                let modulo = 1i128 << 32;
                let residue = lo % modulo;
                let first = if greater || n == pool_len {
                    if residue >= threshold {
                        lo
                    } else {
                        lo + threshold - residue
                    }
                } else {
                    lo + (threshold - residue).rem_euclid(modulo)
                };
                Ok(first < hi)
            }
            _ => Err("expected integer distribution".into()),
        }
    }
    /// Counts cannot rely on the observation end-loss negative/cast clamp.
    pub fn count_integer(&self) -> Result<Box<dyn Distribution<Output = i64>>, String> {
        match self {
            Self::UniformInt { min, max }
                if *min >= 0 && *max > *min && *max - 1 <= u32::MAX as i64 => {}
            Self::EmpiricalInt {
                values,
                cumulative,
                total,
            } => {
                let cumulative = cdf(cumulative, *total)?;
                if values.len() != cumulative.len() {
                    return Err("count CDF shape mismatch".into());
                }
                let mut previous = 0.0;
                for (value, end) in values.iter().zip(cumulative) {
                    if end > previous && !(0..=u32::MAX as i64).contains(value) {
                        return Err("positive count support is outside [0,u32::MAX]".into());
                    }
                    previous = end;
                }
            }
            _ => return Err("invalid native count distribution domain".into()),
        }
        self.integer()
    }
    pub fn integer(&self) -> Result<Box<dyn Distribution<Output = i64>>, String> {
        match self {
            Self::UniformInt { min, max } => {
                let span = max
                    .checked_sub(*min)
                    .filter(|n| *n > 0 && *n <= u32::MAX as i64)
                    .ok_or("invalid integer range")?;
                let _ = span;
                Ok(Box::new(UniformInt::new(*min, *max)))
            }
            Self::EmpiricalInt {
                values,
                cumulative,
                total,
            } => Ok(Box::new(EmpiricalLengthDist::from_replay_parts(
                values.clone(),
                cdf(cumulative, *total)?,
                total.get(),
            )?)),
            _ => Err("expected native integer distribution".into()),
        }
    }
    pub fn base(&self) -> Result<Box<dyn Distribution<Output = u8>>, String> {
        match self {
            Self::UniformBase => Ok(Box::new(UniformBase)),
            Self::CategoricalBase { weights } => {
                let mut w = [0.0; 4];
                for (i, b) in weights.iter().enumerate() {
                    w[i] = finite(*b)?;
                    if w[i] < 0.0 {
                        return Err("negative base weight".into());
                    }
                }
                if !w.iter().sum::<f64>().is_finite() || w.iter().sum::<f64>() <= 0.0 {
                    return Err("invalid base total".into());
                }
                Ok(Box::new(CategoricalBase::from_weights(w)))
            }
            _ => Err("expected native base distribution".into()),
        }
    }
    pub fn allele(
        &self,
        pool: &crate::refdata::AllelePool,
    ) -> Result<Box<dyn Distribution<Output = AlleleId>>, String> {
        match self {
            Self::AllelePool {
                ids,
                cumulative,
                total,
            } => Ok(Box::new(AllelePoolDist::from_replay_parts(
                pool,
                ids.as_ref()
                    .map(|ids| ids.iter().copied().map(AlleleId::new).collect()),
                cdf(cumulative, *total)?,
                total.get(),
            )?)),
            _ => Err("expected native allele distribution".into()),
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum CountSpec {
    Distribution { distribution: DistributionSpec },
    Rate { rate: F64Bits },
}
impl CountSpec {
    pub fn capture(c: &crate::passes::count_source::CountSource) -> Option<Self> {
        match c {
            crate::passes::count_source::CountSource::Distribution(d) => Some(Self::Distribution {
                distribution: d.family_replay_spec()?,
            }),
            crate::passes::count_source::CountSource::Rate(v) => {
                Some(Self::Rate { rate: bits(*v)? })
            }
        }
    }
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct KernelSpecV1 {
    pub schema_version: u32,
    pub mutability: Vec<F64Bits>,
    pub substitution: Vec<F64Bits>,
}
impl KernelSpecV1 {
    pub fn capture(k: &S5FKernel) -> Option<Self> {
        crate::family_replay::primitive_work(
            (k.mutability_table().len() + k.substitution_table().len()) as u64,
        );
        Some(Self {
            schema_version: 1,
            mutability: vec_bits(k.mutability_table())?,
            substitution: vec_bits(k.substitution_table())?,
        })
    }
    pub fn reconstruct(&self) -> Result<S5FKernel, String> {
        if self.schema_version != 1
            || self.mutability.len() != 1024
            || self.substitution.len() != 4096
        {
            return Err("invalid numeric S5F kernel shape/version".into());
        }
        let decode = |v: &[F64Bits]| {
            v.iter()
                .map(|v| {
                    let n = finite(*v)?;
                    if n < 0.0 {
                        Err("negative S5F entry".into())
                    } else {
                        Ok(n)
                    }
                })
                .collect::<Result<Vec<_>, String>>()
        };
        Ok(S5FKernel::new(
            decode(&self.mutability)?,
            decode(&self.substitution)?,
        ))
    }
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DuplexRowSpec {
    pub id: String,
    pub weight: F64Bits,
    pub v_call: String,
    pub d_call: String,
    pub j_call: String,
    pub script: crate::duplex_repair::model::CoupledDuplexScript,
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DuplexModelSpecV1 {
    pub model: String,
    pub reference_content_hash: String,
    pub reference_ends_are_coding_ends: bool,
    pub max_preflight_operations: u64,
    pub max_preflight_bases: u64,
    pub rows: Vec<DuplexRowSpec>,
}
impl DuplexModelSpecV1 {
    pub fn capture(m: &crate::duplex_repair::model::FiniteDuplexModel) -> Option<Self> {
        Some(Self {
            model: m.model.clone(),
            reference_content_hash: m.reference_content_hash.clone(),
            reference_ends_are_coding_ends: m.reference_ends_are_coding_ends,
            max_preflight_operations: m.max_preflight_operations,
            max_preflight_bases: m.max_preflight_bases,
            rows: m
                .rows
                .iter()
                .map(|r| {
                    Some(DuplexRowSpec {
                        id: r.id.clone(),
                        weight: bits(r.weight)?,
                        v_call: r.v_call.clone(),
                        d_call: r.d_call.clone(),
                        j_call: r.j_call.clone(),
                        script: r.script.clone(),
                    })
                })
                .collect::<Option<Vec<_>>>()?,
        })
    }
    pub fn reconstruct(&self) -> Result<crate::duplex_repair::model::FiniteDuplexModel, String> {
        use crate::duplex_repair::model::*;
        let m = FiniteDuplexModel {
            model: self.model.clone(),
            reference_content_hash: self.reference_content_hash.clone(),
            reference_ends_are_coding_ends: self.reference_ends_are_coding_ends,
            max_preflight_operations: self.max_preflight_operations,
            max_preflight_bases: self.max_preflight_bases,
            rows: self
                .rows
                .iter()
                .map(|r| {
                    Ok(DuplexHistoryRow {
                        id: r.id.clone(),
                        weight: finite(r.weight)?,
                        v_call: r.v_call.clone(),
                        d_call: r.d_call.clone(),
                        j_call: r.j_call.clone(),
                        script: r.script.clone(),
                    })
                })
                .collect::<Result<Vec<_>, String>>()?,
        };
        m.validate_schema()?;
        Ok(m)
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum ReplayPassSpecV1 {
    SampleAllele {
        segment: SegmentSpec,
        distribution: DistributionSpec,
    },
    InvertD {
        probability: F64Bits,
    },
    DuplexSource {
        model: DuplexModelSpecV1,
    },
    UniformMutation {
        count: CountSpec,
        base: DistributionSpec,
        segment_rates: [F64Bits; 4],
        v_subregion_rates: [F64Bits; 5],
    },
    S5fMutation {
        count: CountSpec,
        kernel: KernelSpecV1,
        kernel_name: Option<String>,
        segment_rates: [F64Bits; 4],
        v_subregion_rates: [F64Bits; 5],
    },
    Pcr {
        count: CountSpec,
        base: DistributionSpec,
    },
    Quality {
        count: CountSpec,
        base: DistributionSpec,
    },
    Ns {
        count: CountSpec,
    },
    Indel {
        count: CountSpec,
        base: DistributionSpec,
        insertion_probability: F64Bits,
        rejection_attempts: u32,
        exact_node_budget: u64,
        legacy_replay_sentinels: bool,
    },
    EndLoss {
        five_prime: bool,
        distribution: DistributionSpec,
    },
    RevComp {
        probability: F64Bits,
    },
}
pub fn rate_bits(r: crate::passes::mutate::SegmentRateWeights) -> Option<[F64Bits; 4]> {
    Some([bits(r.v)?, bits(r.d)?, bits(r.j)?, bits(r.np)?])
}
pub fn subregion_bits(r: crate::passes::mutate::VSubregionRateWeights) -> Option<[F64Bits; 5]> {
    Some([
        bits(r.fwr1)?,
        bits(r.cdr1)?,
        bits(r.fwr2)?,
        bits(r.cdr2)?,
        bits(r.fwr3)?,
    ])
}
fn rates<const N: usize>(b: &[F64Bits; N]) -> Result<[f64; N], String> {
    let mut a = [0.0; N];
    for (i, v) in b.iter().enumerate() {
        a[i] = finite(*v)?;
        if a[i] < 0.0 {
            return Err("negative mutation rate scalar".into());
        }
    }
    if !a.iter().any(|n| *n > 0.0) {
        return Err("empty mutation rate support".into());
    }
    Ok(a)
}
impl ReplayPassSpecV1 {
    pub fn reconstruct(&self, r: &RefDataConfig) -> Result<Box<dyn Pass>, String> {
        macro_rules! count_pass {
            ($c:expr,$ctor:expr,$rate:expr) => {
                match $c {
                    CountSpec::Distribution { distribution } => {
                        ($ctor)(distribution.count_integer()?)
                    }
                    CountSpec::Rate { rate } => ($rate)(prob(*rate)?),
                }
            };
        }
        Ok(match self {
            Self::SampleAllele {
                segment,
                distribution,
            } => {
                let s = segment.native();
                if !matches!(
                    s,
                    crate::ir::Segment::V | crate::ir::Segment::D | crate::ir::Segment::J
                ) {
                    return Err("invalid allele source segment".into());
                }
                Box::new(SampleAllelePass::new(
                    s,
                    distribution.allele(r.pool_for(s).ok_or("missing source pool")?)?,
                ))
            }
            Self::InvertD { probability } => Box::new(InvertDPass::new(prob(*probability)?)),
            Self::DuplexSource { model } => {
                Box::new(DuplexRecombinePass::new(Arc::new(model.reconstruct()?))?)
            }
            Self::UniformMutation {
                count,
                base,
                segment_rates,
                v_subregion_rates,
            } => {
                let mut p = count_pass!(
                    count,
                    |d| Ok::<_, String>(UniformMutationPass::new(d, base.base()?)),
                    |v| Ok::<_, String>(UniformMutationPass::new_rate(v, base.base()?))
                )?;
                let a = rates(segment_rates)?;
                let b = rates(v_subregion_rates)?;
                p = p
                    .with_segment_rates(crate::passes::mutate::SegmentRateWeights::from_tuple(
                        a[0], a[1], a[2], a[3],
                    ))
                    .with_v_subregion_rates(
                        crate::passes::mutate::VSubregionRateWeights::from_tuple(
                            b[0], b[1], b[2], b[3], b[4],
                        ),
                    );
                Box::new(p)
            }
            Self::S5fMutation {
                count,
                kernel,
                kernel_name,
                segment_rates,
                v_subregion_rates,
            } => {
                let mut p = count_pass!(
                    count,
                    |d| Ok::<_, String>(S5FMutationPass::new(kernel.reconstruct()?, d)),
                    |v| Ok::<_, String>(S5FMutationPass::new_rate(kernel.reconstruct()?, v))
                )?;
                let a = rates(segment_rates)?;
                let b = rates(v_subregion_rates)?;
                p = p
                    .with_segment_rates(crate::passes::mutate::SegmentRateWeights::from_tuple(
                        a[0], a[1], a[2], a[3],
                    ))
                    .with_v_subregion_rates(
                        crate::passes::mutate::VSubregionRateWeights::from_tuple(
                            b[0], b[1], b[2], b[3], b[4],
                        ),
                    );
                if let Some(n) = kernel_name {
                    p = p.with_kernel_name(n.clone());
                }
                Box::new(p)
            }
            Self::Pcr { count, base } => Box::new(count_pass!(
                count,
                |d| Ok::<_, String>(PCRErrorPass::new(d, base.base()?)),
                |v| Ok::<_, String>(PCRErrorPass::new_rate(v, base.base()?))
            )?),
            Self::Quality { count, base } => Box::new(count_pass!(
                count,
                |d| Ok::<_, String>(QualityErrorPass::new(d, base.base()?)),
                |v| Ok::<_, String>(QualityErrorPass::new_rate(v, base.base()?))
            )?),
            Self::Ns { count } => {
                let CountSpec::Distribution { distribution } = count else {
                    return Err("N masking has no rate constructor".into());
                };
                Box::new(NCorruptionPass::new(distribution.count_integer()?))
            }
            Self::Indel {
                count,
                base,
                insertion_probability,
                rejection_attempts,
                exact_node_budget,
                legacy_replay_sentinels,
            } => {
                if *legacy_replay_sentinels {
                    return Err(
                        "legacy indel sentinel policy is unsupported in family replay v1".into(),
                    );
                }
                let CountSpec::Distribution { distribution } = count else {
                    return Err("indel has no rate constructor".into());
                };
                Box::new(
                    IndelPass::new(
                        distribution.count_integer()?,
                        prob(*insertion_probability)?,
                        base.base()?,
                    )
                    .with_sampling_options(IndelSamplingOptions {
                        rejection_attempts: *rejection_attempts,
                        exact_node_budget: *exact_node_budget,
                        legacy_replay_sentinels: *legacy_replay_sentinels,
                    }),
                )
            }
            Self::EndLoss {
                five_prime,
                distribution,
            } => Box::new(EndLossPass::new(
                if *five_prime {
                    LossEnd::Five
                } else {
                    LossEnd::Three
                },
                distribution.integer()?,
            )),
            Self::RevComp { probability } => Box::new(RevCompPass::new(prob(*probability)?)),
        })
    }
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum ContractSpec {
    CodingSequenceIntegrity,
    ProductiveJunctionFrame,
    NoStopCodonInJunction,
    CanonicalAnchor { segment: SegmentSpec },
}
impl ContractSpec {
    pub fn capture(c: &dyn crate::contract::Contract) -> Result<Self, String> {
        use crate::contract::ContractKind::*;
        Ok(match c.kind() {
            CodingSequenceIntegrity => Self::CodingSequenceIntegrity,
            ProductiveJunctionFrame => Self::ProductiveJunctionFrame,
            NoStopCodonInJunction => Self::NoStopCodonInJunction,
            CanonicalAnchor { segment } => Self::CanonicalAnchor {
                segment: segment.into(),
            },
            _ => return Err("unsupported family contract".into()),
        })
    }
    pub fn reconstruct(&self) -> Result<Box<dyn crate::contract::Contract>, String> {
        use crate::contract::*;
        Ok(match self {
            Self::CodingSequenceIntegrity => Box::new(CodingSequenceIntegrity),
            Self::ProductiveJunctionFrame => Box::new(ProductiveJunctionFrame::new()),
            Self::NoStopCodonInJunction => Box::new(NoStopCodonInJunction::new()),
            Self::CanonicalAnchor { segment } => {
                let s = segment.native();
                if !matches!(s, crate::ir::Segment::V | crate::ir::Segment::J) {
                    return Err("unsupported canonical anchor segment".into());
                }
                Box::new(AnchorPreserved::canonical(s))
            }
        })
    }
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ReplayPlanSpecV1 {
    pub schema_version: u32,
    pub reference: ReferenceSpecV1,
    pub passes: Vec<ReplayPassSpecV1>,
    pub contracts: Option<Vec<ContractSpec>>,
    pub strict: bool,
    pub allow_curatable_refdata: bool,
    pub productivity_policy: String,
}
impl ReplayPlanSpecV1 {
    pub fn capture(
        plan: &PassPlan,
        r: &RefDataConfig,
        contracts: Option<&crate::contract::ContractSet>,
        strict: bool,
        allow_curatable_refdata: bool,
    ) -> Result<Self, String> {
        Ok(Self {
            schema_version: 1,
            reference: ReferenceSpecV1::capture(r)?,
            passes: plan
                .passes()
                .iter()
                .map(|p| {
                    p.family_replay_spec()
                        .ok_or_else(|| format!("unsupported native replay codec for {}", p.name()))
                })
                .collect::<Result<Vec<_>, _>>()?,
            contracts: contracts
                .map(|c| {
                    c.iter()
                        .map(ContractSpec::capture)
                        .collect::<Result<Vec<_>, _>>()
                })
                .transpose()?,
            strict,
            allow_curatable_refdata,
            productivity_policy: crate::productivity::DUPLEX_POLICY.into(),
        })
    }
    pub fn reconstruct(
        &self,
    ) -> Result<
        (
            PassPlan,
            RefDataConfig,
            Option<crate::contract::ContractSet>,
        ),
        String,
    > {
        if self.schema_version != 1
            || self.productivity_policy != crate::productivity::DUPLEX_POLICY
        {
            return Err("unsupported plan codec/productivity policy".into());
        }
        let r = self.reference.reconstruct()?;
        let mut p = PassPlan::new();
        for s in &self.passes {
            p.push(s.reconstruct(&r)?);
        }
        let c = self
            .contracts
            .as_ref()
            .map(|specs| {
                let mut c = crate::contract::ContractSet::new();
                for s in specs {
                    c.add(s.reconstruct()?);
                }
                Ok::<_, String>(c)
            })
            .transpose()?;
        Ok((p, r, c))
    }
    pub fn digest(&self) -> Result<String, String> {
        super::canonical_digest(
            "replay_plan_spec_v1",
            &serde_json::to_value(self).map_err(|e| e.to_string())?,
        )
    }
}

#[cfg(test)]
mod tests;
