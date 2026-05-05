//! Module-level Python entry points that build and execute plans.
//!
//! This file holds the minimal smoke runners that prove the
//! Cargo → maturin → Python pipeline works end-to-end. The
//! production-style entry points (compiled from the user-facing
//! Experiment DSL) arrive in later F.x steps.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use crate::dist::{AllelePoolDist, EmpiricalLengthDist, UniformBase};
use crate::ir::{flag, Segment, Simulation};
use crate::pass::{PassPlan, PassRuntime};
use crate::passes::{
    AssembleSegmentPass, EchoPass, GenerateNPPass, SampleAllelePass, SampleBasePass,
};
use crate::refdata::{Allele, ChainType, RefDataConfig};

use super::outcome::PyOutcome;
use super::plan::PyPassPlan;
use super::refdata::PyRefDataConfig;

/// Run a small built-in mixed plan (8 × `EchoPass` interleaved with
/// 8 × `SampleBasePass`) and return the resulting [`PyOutcome`].
///
/// Used by the F.1 Python tests to exercise the read-only wrapper
/// types end-to-end. Has no biological meaning — it's the smallest
/// plan that produces a non-empty pool *and* a non-empty trace, so
/// every accessor on every wrapper has something to read back.
#[pyfunction]
fn run_smoke_plan(seed: u64) -> PyOutcome {
    let mut plan = PassPlan::new();
    for i in 0..8 {
        plan.push(Box::new(EchoPass::new(b'A', i as u16, Segment::V)));
        plan.push(Box::new(SampleBasePass::new(
            format!("np.np1.bases[{}]", i),
            Box::new(UniformBase),
            Segment::Np1,
            flag::N_NUC,
        )));
    }

    let outcome = PassRuntime::execute(&plan, Simulation::new(), seed);
    PyOutcome::new(outcome)
}

/// Build the synthetic VJ refdata used by [`run_smoke_vj_recombination`]:
/// V `AAACCCGGG` (9bp, anchor 6) + J `TTTAAA` (6bp, anchor 0). Same
/// fixture shape as the D.7 / E.8 integration tests so the resulting
/// outcomes have predictable structure (junction = `[6, 12)`,
/// in-frame at NP1 lengths divisible by 3).
fn smoke_vj_refdata() -> RefDataConfig {
    let mut cfg = RefDataConfig::empty(ChainType::Vj);
    let _ = cfg.v_pool.push(Allele {
        name: "v_smoke*01".into(),
        gene: "v_smoke".into(),
        seq: b"AAACCCGGG".to_vec(),
        segment: Segment::V,
        anchor: Some(6),
    });
    let _ = cfg.j_pool.push(Allele {
        name: "j_smoke*01".into(),
        gene: "j_smoke".into(),
        seq: b"TTTAAA".to_vec(),
        segment: Segment::J,
        anchor: Some(0),
    });
    cfg
}

/// Run a full VJ recombination plan against a built-in synthetic
/// refdata and return the resulting [`PyOutcome`].
///
/// Pipeline:
/// 1. Sample the V allele (single-allele pool → always id 0).
/// 2. Sample the J allele (single-allele pool → always id 0).
/// 3. Assemble V into the pool (9 bases of `AAACCCGGG`).
/// 4. Generate NP1 (length distribution `0..=6`, base distribution
///    `UniformBase`).
/// 5. Assemble J into the pool (6 bases of `TTTAAA`).
///
/// Used by F.2 Python tests to exercise [`PyRegion`], V/J allele-id
/// accessors, and the `sample_allele.*` / `np.np1.*` trace addresses.
/// Replaced in F.3 by a runner that takes a real `RefDataConfig`
/// from Python.
#[pyfunction]
fn run_smoke_vj_recombination(seed: u64) -> PyOutcome {
    let cfg = smoke_vj_refdata();

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
        Box::new(EmpiricalLengthDist::from_pairs(
            (0..7).map(|i| (i, 1.0)).collect::<Vec<_>>(),
        )),
        Box::new(UniformBase),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::J)));

    let outcome = PassRuntime::execute_with_refdata(&plan, Simulation::new(), seed, &cfg);
    PyOutcome::new(outcome)
}

/// Run a VJ recombination plan against a Python-supplied
/// [`PyRefDataConfig`].
///
/// Plan shape (identical to [`run_smoke_vj_recombination`]):
/// 1. Sample V allele (uniform over `refdata.v_pool`).
/// 2. Sample J allele (uniform over `refdata.j_pool`).
/// 3. Assemble V.
/// 4. Generate NP1 with length distribution `0..np1_max_length`
///    (uniform over the integers in that range), bases drawn from
///    `UniformBase`.
/// 5. Assemble J.
///
/// `np1_max_length` defaults to 7 (i.e. lengths 0..=6 with equal
/// weight), matching the existing baked-in fixture so test
/// expectations carry over between runners.
///
/// Returns an error when:
/// - the refdata's chain type is `vdj` (use a VDJ runner instead),
/// - the V or J pool is empty,
/// - `np1_max_length` is `< 1`.
#[pyfunction]
#[pyo3(signature = (refdata, seed, *, np1_max_length=7))]
fn run_vj_recombination(
    refdata: &PyRefDataConfig,
    seed: u64,
    np1_max_length: i64,
) -> PyResult<PyOutcome> {
    let cfg = refdata.inner();
    if !matches!(cfg.chain_type, ChainType::Vj) {
        return Err(PyValueError::new_err(
            "run_vj_recombination requires a VJ-chain RefDataConfig",
        ));
    }
    if cfg.v_pool.len() == 0 {
        return Err(PyValueError::new_err("V pool is empty"));
    }
    if cfg.j_pool.len() == 0 {
        return Err(PyValueError::new_err("J pool is empty"));
    }
    if np1_max_length < 1 {
        return Err(PyValueError::new_err(
            "np1_max_length must be at least 1 (must include length 0 plus a positive ceiling)",
        ));
    }

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
        Box::new(EmpiricalLengthDist::from_pairs(
            (0..np1_max_length).map(|i| (i, 1.0)).collect::<Vec<_>>(),
        )),
        Box::new(UniformBase),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::J)));

    let outcome = PassRuntime::execute_with_refdata(&plan, Simulation::new(), seed, cfg);
    Ok(PyOutcome::new(outcome))
}

/// Execute a Python-built [`PyPassPlan`] and return the resulting
/// [`PyOutcome`].
///
/// `refdata` is optional: pass it for plans that include
/// `AssembleSegmentPass` or any other pass that requires reference
/// data. Plans that only use NP / sample-base / echo passes can run
/// without one.
#[pyfunction]
#[pyo3(signature = (plan, seed, *, refdata=None))]
fn run(
    plan: &PyPassPlan,
    seed: u64,
    refdata: Option<&PyRefDataConfig>,
) -> PyOutcome {
    let outcome = match refdata {
        Some(r) => PassRuntime::execute_with_refdata(
            plan.inner(),
            Simulation::new(),
            seed,
            r.inner(),
        ),
        None => PassRuntime::execute(plan.inner(), Simulation::new(), seed),
    };
    PyOutcome::new(outcome)
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(run_smoke_plan, m)?)?;
    m.add_function(wrap_pyfunction!(run_smoke_vj_recombination, m)?)?;
    m.add_function(wrap_pyfunction!(run_vj_recombination, m)?)?;
    m.add_function(wrap_pyfunction!(run, m)?)?;
    Ok(())
}
