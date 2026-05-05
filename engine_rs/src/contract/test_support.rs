//! Shared test fixtures used across contract test modules.
//!
//! Compiled only under `#[cfg(test)]`. Lives in its own file so
//! per-contract test modules can reuse the V/J refdata builders
//! without duplicating setup boilerplate.

#![cfg(test)]
#![allow(dead_code)]

use crate::assignment::AlleleInstance;
use crate::ir::{flag, NucHandle, Nucleotide, Region, Segment, Simulation};
use crate::refdata::{Allele, AlleleId, ChainType, RefDataConfig};

/// Build a `RefDataConfig` whose V pool has a single allele of length
/// `len` with the given anchor. Used by `AnchorPreserved` tests and
/// the `ContractSet` verify tests.
pub fn make_v_anchor_at(len: u32, anchor: Option<u16>) -> RefDataConfig {
    let mut cfg = RefDataConfig::empty(ChainType::Vdj);
    let _ = cfg.v_pool.push(Allele {
        name: "v_test*01".to_string(),
        gene: "v_test".to_string(),
        seq: vec![b'A'; len as usize],
        segment: Segment::V,
        anchor,
    });
    cfg
}

/// Build a complete V+J refdata with controllable anchors.
/// V: "AAACCCGGG" (9bp, anchor at v_anchor)
/// J: "TTTAAA" (6bp, anchor at j_anchor)
pub fn make_vj_for_frame_test(v_anchor: Option<u16>, j_anchor: Option<u16>) -> RefDataConfig {
    let mut cfg = RefDataConfig::empty(ChainType::Vdj);
    let _ = cfg.v_pool.push(Allele {
        name: "v_test*01".into(),
        gene: "v_test".into(),
        seq: b"AAACCCGGG".to_vec(),
        segment: Segment::V,
        anchor: v_anchor,
    });
    let _ = cfg.j_pool.push(Allele {
        name: "j_test*01".into(),
        gene: "j_test".into(),
        seq: b"TTTAAA".to_vec(),
        segment: Segment::J,
        anchor: j_anchor,
    });
    cfg
}

/// Build an assembled sim with V at `[v_pool_start, v_pool_start +
/// v_len)`, optional NP1 padding, and J at `[j_pool_start,
/// j_pool_start + j_len)`. Used by `ProductiveJunctionFrame` tests
/// and `ContractSet` verify tests.
pub fn make_assembled_sim(
    v_pool_start: u32,
    v_len: u32,
    j_pool_start: u32,
    j_len: u32,
    v_inst: AlleleInstance,
    j_inst: AlleleInstance,
) -> Simulation {
    let mut sim = Simulation::new();
    for i in 0..v_pool_start {
        let (next, _) =
            sim.with_nucleotide_pushed(Nucleotide::germline(b'X', i as u16, Segment::V));
        sim = next;
    }
    for i in 0..v_len {
        let (next, _) =
            sim.with_nucleotide_pushed(Nucleotide::germline(b'A', i as u16, Segment::V));
        sim = next;
    }
    for _ in v_pool_start + v_len..j_pool_start {
        let (next, _) =
            sim.with_nucleotide_pushed(Nucleotide::synthetic(b'a', Segment::Np1, flag::N_NUC));
        sim = next;
    }
    for i in 0..j_len {
        let (next, _) =
            sim.with_nucleotide_pushed(Nucleotide::germline(b'T', i as u16, Segment::J));
        sim = next;
    }
    sim = sim.with_region_added(Region::new(
        Segment::V,
        NucHandle::new(v_pool_start),
        NucHandle::new(v_pool_start + v_len),
    ));
    sim = sim.with_region_added(Region::new(
        Segment::J,
        NucHandle::new(j_pool_start),
        NucHandle::new(j_pool_start + j_len),
    ));
    sim.with_allele_assigned(Segment::V, v_inst)
        .with_allele_assigned(Segment::J, j_inst)
}

/// Build a refdata with V allele "AAACCCxxx" (anchor at 6, so
/// junction starts at AAACCCxxx[6..9] = "xxx") and J allele
/// "yyyAAA" (anchor at 0, so junction ends at yyyAAA[0..3] +
/// 3 = first 3 chars). The exact bases at the V anchor codon
/// and J anchor codon control whether the junction has stops.
pub fn make_vj_with_anchor_codons(
    v_anchor_codon: &[u8; 3],
    j_anchor_codon: &[u8; 3],
) -> RefDataConfig {
    let mut v_seq = b"AAACCC".to_vec();
    v_seq.extend_from_slice(v_anchor_codon);
    // V is now 9 bases, anchor at 6.

    let mut j_seq = j_anchor_codon.to_vec();
    j_seq.extend_from_slice(b"AAA");
    // J is now 6 bases, anchor at 0.

    let mut cfg = RefDataConfig::empty(ChainType::Vdj);
    let _ = cfg.v_pool.push(Allele {
        name: "v_test*01".into(),
        gene: "v_test".into(),
        seq: v_seq,
        segment: Segment::V,
        anchor: Some(6),
    });
    let _ = cfg.j_pool.push(Allele {
        name: "j_test*01".into(),
        gene: "j_test".into(),
        seq: j_seq,
        segment: Segment::J,
        anchor: Some(0),
    });
    cfg
}

/// Build a simulation that copies the V/J sequences from `cfg`
/// into the pool with the given pool placement, no NP padding.
/// V anchor codon ends at V[8], J anchor codon at J[0..3], so
/// junction is `[V_anchor_pool, J_anchor_pool + 3) = [6, 12)`,
/// which spans the V anchor codon (3 bases) + 3 bases of J.
/// Length 6 → in frame.
pub fn make_assembled_sim_from_refdata(cfg: &RefDataConfig) -> Simulation {
    let v_allele = cfg.v_pool.get(AlleleId::new(0)).unwrap();
    let j_allele = cfg.j_pool.get(AlleleId::new(0)).unwrap();

    let mut sim = Simulation::new();
    for (i, &b) in v_allele.seq.iter().enumerate() {
        let (next, _) = sim.with_nucleotide_pushed(Nucleotide::germline(b, i as u16, Segment::V));
        sim = next;
    }
    for (i, &b) in j_allele.seq.iter().enumerate() {
        let (next, _) = sim.with_nucleotide_pushed(Nucleotide::germline(b, i as u16, Segment::J));
        sim = next;
    }
    sim = sim.with_region_added(Region::new(
        Segment::V,
        NucHandle::new(0),
        NucHandle::new(v_allele.seq.len() as u32),
    ));
    sim = sim.with_region_added(Region::new(
        Segment::J,
        NucHandle::new(v_allele.seq.len() as u32),
        NucHandle::new((v_allele.seq.len() + j_allele.seq.len()) as u32),
    ));
    sim.with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(0)))
        .with_allele_assigned(Segment::J, AlleleInstance::new(AlleleId::new(0)))
}
