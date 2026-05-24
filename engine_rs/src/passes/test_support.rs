//! Shared test fixtures used across multiple pass test modules.
//!
//! Compiled only under `#[cfg(test)]`. Lives in its own file so that
//! per-pass test modules can reuse stop-codon-filtering fixtures
//! without duplicating ~50 lines of setup each.

#![cfg(test)]
#![allow(dead_code)]

use crate::assignment::AlleleInstance;
use crate::dist::Distribution;
use crate::ir::{NucHandle, Nucleotide, Region, Segment, Simulation};
use crate::refdata::{Allele, AlleleId, ChainType, RefDataConfig};

/// A base distribution that always samples `A` but reports a
/// 2-element support `{A, C}`. Use to demonstrate that an active
/// productive-contract filter can divert from the default `A` to
/// the safe alternative `C` when `A` would create a stop codon.
#[derive(Clone, Debug)]
pub(crate) struct StopThenSafeMutationBaseDist;

impl Distribution for StopThenSafeMutationBaseDist {
    type Output = u8;

    fn sample(&self, _rng: &mut crate::rng::Rng) -> u8 {
        b'A'
    }

    fn support(&self) -> Option<Vec<(u8, f64)>> {
        Some(vec![(b'A', 1.0), (b'C', 1.0)])
    }
}

/// A degenerate distribution whose support is `{A}` only. When a
/// contract rejects `A` at the filter stage, strict mode must
/// surface `EmptyAdmissibleSupport`. Use to exercise the strict-mode
/// failure path.
#[derive(Clone, Debug)]
pub(crate) struct StopOnlyMutationBaseDist;

impl Distribution for StopOnlyMutationBaseDist {
    type Output = u8;

    fn sample(&self, _rng: &mut crate::rng::Rng) -> u8 {
        b'A'
    }

    fn support(&self) -> Option<Vec<(u8, f64)>> {
        Some(vec![(b'A', 1.0)])
    }
}

/// VJ fixture used by the substitution-style pass tests
/// (uniform-mutation, PCR error, contaminant). The pool already
/// holds `TAC` (V) + `TGG` (J) so position 2's neighbours are
/// `T A · T G G`. With anchors at V[0] and J[0], the junction is
/// `[0, 6)` and codon 0 (`TAC`) becomes `TAA` if position 2 is
/// substituted to `A` — the perfect stop-codon filter probe.
pub(crate) fn make_substitution_productive_vj_fixture() -> (RefDataConfig, Simulation) {
    let mut cfg = RefDataConfig::empty(ChainType::Vj);
    let _ = cfg.v_pool.push(Allele {
        name: "v_mut*01".into(),
        gene: "v_mut".into(),
        seq: b"TAC".to_vec(),
        segment: Segment::V,
        anchor: Some(0),
    });
    let _ = cfg.j_pool.push(Allele {
        name: "j_mut*01".into(),
        gene: "j_mut".into(),
        seq: b"TGG".to_vec(),
        segment: Segment::J,
        anchor: Some(0),
    });

    let mut sim = Simulation::new();
    for (i, &b) in b"TAC".iter().enumerate() {
        let (next, _) = sim.with_nucleotide_pushed(Nucleotide::germline(b, i as u16, Segment::V));
        sim = next;
    }
    let v_region = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(3));
    sim = sim.with_region_added(v_region);

    for (i, &b) in b"TGG".iter().enumerate() {
        let (next, _) = sim.with_nucleotide_pushed(Nucleotide::germline(b, i as u16, Segment::J));
        sim = next;
    }
    let j_region = Region::new(Segment::J, NucHandle::new(3), NucHandle::new(6));
    sim = sim.with_region_added(j_region);

    sim = sim
        .with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(0)))
        .with_allele_assigned(Segment::J, AlleleInstance::new(AlleleId::new(0)));

    (cfg, sim)
}
