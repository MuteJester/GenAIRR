use super::*;
use crate::airr_record::build_airr_record;
use crate::assignment::TrimEnd;
use crate::contract::{Contract, ContractViolation};
use crate::dist::{AllelePoolDist, EmpiricalLengthDist, UniformBase};
use crate::ir::{NucFlags, Nucleotide, Segment};
use crate::pass::Pass;
use crate::pass::testing::PassRuntime;
use crate::passes::{AssembleSegmentPass, GenerateNPPass, SampleAllelePass, TrimPass};
use crate::refdata::{Allele, ChainType};
use crate::trace::ChoiceValue;

fn vj_refdata() -> RefDataConfig {
    let mut cfg = RefDataConfig::empty(ChainType::Vj);
    let _ = cfg.v_pool.push(Allele {
        name: "v*01".into(),
        gene: "v".into(),
        seq: b"AAACCCGGG".to_vec(),
        segment: Segment::V,
        anchor: Some(6),
    });
    let _ = cfg.j_pool.push(Allele {
        name: "j*01".into(),
        gene: "j".into(),
        seq: b"TTTAAA".to_vec(),
        segment: Segment::J,
        anchor: Some(0),
    });
    cfg
}

fn vj_plan(cfg: &RefDataConfig) -> PassPlan {
    vj_plan_with_np_lengths(cfg, vec![(0, 1.0)])
}

fn vj_plan_with_np_lengths(cfg: &RefDataConfig, np_lengths: Vec<(i64, f64)>) -> PassPlan {
    let mut plan = PassPlan::new();
    plan.push(Box::new(SampleAllelePass::new(
        Segment::V,
        Box::new(AllelePoolDist::uniform(&cfg.v_pool)),
    )));
    plan.push(Box::new(SampleAllelePass::new(
        Segment::J,
        Box::new(AllelePoolDist::uniform(&cfg.j_pool)),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
    plan.push(Box::new(GenerateNPPass::new(
        Segment::Np1,
        Box::new(EmpiricalLengthDist::from_pairs(np_lengths)),
        Box::new(UniformBase),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::J)));
    plan
}

fn vj_refdata_with_runtime_j_residue() -> (RefDataConfig, AlleleId, AlleleId) {
    let mut cfg = RefDataConfig::empty(ChainType::Vj);
    let _ = cfg.v_pool.push(Allele {
        name: "v*01".into(),
        gene: "v".into(),
        seq: b"AAACCCGGG".to_vec(),
        segment: Segment::V,
        anchor: Some(6),
    });
    let j_good = cfg.j_pool.push(Allele {
        name: "j_good*01".into(),
        gene: "j_good".into(),
        seq: b"TTTAAA".to_vec(),
        segment: Segment::J,
        anchor: Some(0),
    });
    let j_bad = cfg.j_pool.push(Allele {
        name: "j_bad*01".into(),
        gene: "j_bad".into(),
        seq: b"TTTAAA".to_vec(),
        segment: Segment::J,
        anchor: Some(1),
    });
    (cfg, j_good, j_bad)
}

struct AppendBasePass {
    name: &'static str,
    base: u8,
}

impl AppendBasePass {
    fn new(name: &'static str, base: u8) -> Self {
        Self { name, base }
    }
}

impl Pass for AppendBasePass {
    fn name(&self) -> &str {
        self.name
    }

    fn execute(&self, sim: &Simulation, _ctx: &mut PassContext) -> Simulation {
        let (next, _handle) = sim.with_nucleotide_pushed(Nucleotide::synthetic(
            self.base,
            Segment::Np1,
            NucFlags::empty(),
        ));
        next
    }

    fn effects(&self) -> Vec<PassEffect> {
        vec![PassEffect::AppendNucleotides]
    }
}

struct TraceProbePass {
    name: &'static str,
    address: &'static str,
    base: u8,
}

impl TraceProbePass {
    fn new(name: &'static str, address: &'static str, base: u8) -> Self {
        Self {
            name,
            address,
            base,
        }
    }
}

impl Pass for TraceProbePass {
    fn name(&self) -> &str {
        self.name
    }

    fn execute(&self, sim: &Simulation, ctx: &mut PassContext) -> Simulation {
        assert!(
            ctx.trace.is_empty(),
            "{} saw previously committed trace records",
            self.name
        );
        ctx.trace.record(self.address, ChoiceValue::Base(self.base));
        let (next, _handle) = sim.with_nucleotide_pushed(Nucleotide::synthetic(
            self.base,
            Segment::Np1,
            NucFlags::empty(),
        ));
        next
    }

    fn declared_choices(&self) -> Vec<String> {
        vec![self.address.to_string()]
    }

    fn effects(&self) -> Vec<PassEffect> {
        vec![PassEffect::AppendNucleotides]
    }
}

struct MaxPoolLen {
    max: usize,
}

impl MaxPoolLen {
    fn new(max: usize) -> Self {
        Self { max }
    }
}

impl Contract for MaxPoolLen {
    fn name(&self) -> &str {
        "test.max_pool_len"
    }

    fn verify(
        &self,
        sim: &Simulation,
        _refdata: Option<&RefDataConfig>,
    ) -> Result<(), ContractViolation> {
        if sim.pool.len() <= self.max {
            Ok(())
        } else {
            Err(ContractViolation::new(
                self.name(),
                format!("pool length {} exceeds {}", sim.pool.len(), self.max),
            ))
        }
    }
}

/// Test-only deterministic base distribution. Always samples the
/// configured byte; reports it as the only point in `support()`.
/// Used by mutation- and NP-pass fixtures to produce known base
/// sequences without seed-pinning.
#[derive(Clone, Debug)]
struct ConstBaseDist(u8);

impl crate::dist::Distribution for ConstBaseDist {
    type Output = u8;
    fn sample(&self, _rng: &mut crate::rng::Rng) -> u8 {
        self.0
    }
    fn support(&self) -> Option<Vec<(u8, f64)>> {
        Some(vec![(self.0, 1.0)])
    }
}

mod compile_runtime;
mod curated;
mod d_extensions;
pub(super) mod indels;
mod j_extensions;
pub(super) mod live_call_edits;
mod overlaps;
mod v_extensions;
