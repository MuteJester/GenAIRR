//! Shared fixtures for the metadata-truth integration suite.
//!
//! All helpers here build VJ / VDJ test simulations from a small
//! curated `RefDataConfig` whose alleles are hand-designed to exercise
//! ambiguity scenarios:
//!
//! - V alleles share suffixes so that 5' trims of the same length
//!   leave overlapping bases that become ambiguous between V01/V02,
//!   and the same trims leave V03 distinguishable by a single byte
//!   so callers can test the "ambiguity narrows back" path when
//!   NP1 recreates the missing suffix.
//! - D alleles share both prefixes and suffixes so left and right
//!   trims independently widen / narrow the D call.
//! - J alleles differ in their 5' prefix so left trims widen J calls
//!   and NP1/NP2 extension can narrow them back.
//!
//! Use these fixtures as a starting point. If your test needs a
//! different distinguishing-byte layout, build a one-off
//! `RefDataConfig` locally — the helpers here are conveniences, not a
//! contract.

#![allow(dead_code)] // not every helper is consumed by every category

use genairr_engine::ir::{Region, Segment, Simulation};
use genairr_engine::refdata::{Allele, AlleleId, ChainType, RefDataConfig};

/// Default seed used by the small smoke tests in each category file.
/// Individual tests can override.
pub const DEFAULT_SEED: u64 = 42;

/// Build a small curated VJ refdata with hand-designed ambiguity.
///
/// V pool (length 12 each):
/// - `v01*01` AAACCCGGGTTT  (V01)
/// - `v01*02` AAACCCGGGTTT  (alias of V01 — same nt, different name; both
///   should appear in the v_call set whenever V01 is the truth)
/// - `v02*01` AAACCCGGGAAA  (last 3 differ from v01)
/// - `v03*01` AAACCCGGGCCC  (last 3 differ from v01 and v02)
///
/// All three real alleles share the first 9 bases. A 3-base 3' trim
/// makes them indistinguishable; NP1 bases can later recreate the
/// trimmed suffix and narrow the call set back to the truth.
///
/// J pool (length 9 each):
/// - `j01*01` TTTAAACCC  (J01)
/// - `j02*01` GGGAAACCC  (last 6 shared with J01; first 3 differ)
///
/// A 3-base 5' trim of J makes both J alleles indistinguishable.
pub fn vj_ambiguous_refdata() -> RefDataConfig {
    let mut cfg = RefDataConfig::empty(ChainType::Vj);
    let _ = cfg.v_pool.push(Allele {
        name: "v01*01".into(),
        gene: "v01".into(),
        seq: b"AAACCCGGGTTT".to_vec(),
        segment: Segment::V,
        anchor: Some(9),
    });
    let _ = cfg.v_pool.push(Allele {
        name: "v01*02".into(),
        gene: "v01".into(),
        seq: b"AAACCCGGGTTT".to_vec(),
        segment: Segment::V,
        anchor: Some(9),
    });
    let _ = cfg.v_pool.push(Allele {
        name: "v02*01".into(),
        gene: "v02".into(),
        seq: b"AAACCCGGGAAA".to_vec(),
        segment: Segment::V,
        anchor: Some(9),
    });
    let _ = cfg.v_pool.push(Allele {
        name: "v03*01".into(),
        gene: "v03".into(),
        seq: b"AAACCCGGGCCC".to_vec(),
        segment: Segment::V,
        anchor: Some(9),
    });
    let _ = cfg.j_pool.push(Allele {
        name: "j01*01".into(),
        gene: "j01".into(),
        seq: b"TTTAAACCC".to_vec(),
        segment: Segment::J,
        anchor: Some(0),
    });
    let _ = cfg.j_pool.push(Allele {
        name: "j02*01".into(),
        gene: "j02".into(),
        seq: b"GGGAAACCC".to_vec(),
        segment: Segment::J,
        anchor: Some(0),
    });
    cfg
}

/// Build a small curated VDJ refdata with ambiguity on every segment.
///
/// V pool: 4 alleles sharing first 9 bytes (same layout as
/// [`vj_ambiguous_refdata`]).
///
/// D pool (length 9 each):
/// - `d01*01` AAATTTGGG  (D01)
/// - `d02*01` AAATTTCCC  (shares 6-byte prefix with D01)
/// - `d03*01` CCCTTTGGG  (shares 6-byte suffix with D01)
///
/// 3-base 5' trim of D makes D01/D02 indistinguishable (both share
/// the surviving suffix TTTGGG / TTTCCC respectively — actually only
/// D01/D03 share the surviving TTTGGG; D02 has CCC). 3-base 3' trim
/// makes D01/D02 indistinguishable.
///
/// J pool: same as [`vj_ambiguous_refdata`].
pub fn vdj_ambiguous_refdata() -> RefDataConfig {
    let mut cfg = RefDataConfig::empty(ChainType::Vdj);
    for v in [
        ("v01*01", "v01", "AAACCCGGGTTT"),
        ("v01*02", "v01", "AAACCCGGGTTT"),
        ("v02*01", "v02", "AAACCCGGGAAA"),
        ("v03*01", "v03", "AAACCCGGGCCC"),
    ] {
        let _ = cfg.v_pool.push(Allele {
            name: v.0.into(),
            gene: v.1.into(),
            seq: v.2.as_bytes().to_vec(),
            segment: Segment::V,
            anchor: Some(9),
        });
    }
    for d in [
        ("d01*01", "d01", "AAATTTGGG"),
        ("d02*01", "d02", "AAATTTCCC"),
        ("d03*01", "d03", "CCCTTTGGG"),
    ] {
        let _ = cfg.d_pool.push(Allele {
            name: d.0.into(),
            gene: d.1.into(),
            seq: d.2.as_bytes().to_vec(),
            segment: Segment::D,
            anchor: None,
        });
    }
    for j in [
        ("j01*01", "j01", "TTTAAACCC"),
        ("j02*01", "j02", "GGGAAACCC"),
    ] {
        let _ = cfg.j_pool.push(Allele {
            name: j.0.into(),
            gene: j.1.into(),
            seq: j.2.as_bytes().to_vec(),
            segment: Segment::J,
            anchor: Some(0),
        });
    }
    cfg
}

/// Look up an allele id by name in the given segment pool.
///
/// Panics if not found — fixtures with a typo should fail loudly.
pub fn allele_id_by_name(refdata: &RefDataConfig, segment: Segment, name: &str) -> AlleleId {
    let pool = match segment {
        Segment::V => &refdata.v_pool,
        Segment::D => &refdata.d_pool,
        Segment::J => &refdata.j_pool,
        Segment::Np1 | Segment::Np2 => panic!("NP segments have no allele pool"),
    };
    for i in 0..pool.len() {
        let a = pool.get(AlleleId::new(i as u32)).unwrap();
        if a.name == name {
            return AlleleId::new(i as u32);
        }
    }
    panic!(
        "allele_id_by_name: {} not found in {:?} pool",
        name, segment
    );
}

/// Find the v_call name set on the assembled simulation. Returns the
/// set as a sorted Vec for easy `==` comparison in tests.
pub fn v_call_names(sim: &Simulation, refdata: &RefDataConfig) -> Vec<String> {
    call_names_for_segment(sim, refdata, Segment::V)
}

pub fn d_call_names(sim: &Simulation, refdata: &RefDataConfig) -> Vec<String> {
    call_names_for_segment(sim, refdata, Segment::D)
}

pub fn j_call_names(sim: &Simulation, refdata: &RefDataConfig) -> Vec<String> {
    call_names_for_segment(sim, refdata, Segment::J)
}

fn call_names_for_segment(
    sim: &Simulation,
    refdata: &RefDataConfig,
    segment: Segment,
) -> Vec<String> {
    let pool = match segment {
        Segment::V => &refdata.v_pool,
        Segment::D => &refdata.d_pool,
        Segment::J => &refdata.j_pool,
        Segment::Np1 | Segment::Np2 => return Vec::new(),
    };
    let Some(seg_call) = sim.segment_calls.get(segment) else {
        return Vec::new();
    };
    let mut names: Vec<String> = Vec::new();
    seg_call.allele_call.for_each_id(|id| {
        if let Some(a) = pool.get(id) {
            names.push(a.name.clone());
        }
    });
    names.sort();
    names
}

/// Helper to look up a region by segment on a sealed simulation.
/// Returns `None` if not assembled.
pub fn region_for_segment(sim: &Simulation, segment: Segment) -> Option<&Region> {
    sim.sequence.regions.iter().find(|r| r.segment == segment)
}

/// Build a VJ refdata where three V alleles share the exact same
/// nucleotide sequence under different names ("triplet alias"), plus
/// an additional distinguishable allele.
///
/// V pool (length 12):
/// - `v_triplet*01` AAACCCGGGTTT
/// - `v_triplet*02` AAACCCGGGTTT  (alias)
/// - `v_triplet*03` AAACCCGGGTTT  (alias)
/// - `v_unique*01`  AAACCCGGGAAA  (distinguishable in tail)
///
/// J pool: single `j_single*01` allele = TTTAAACCC.
///
/// Whenever any of the three triplet alleles is the sampled truth,
/// all three should appear in `v_call`; `v_unique*01` must NOT appear
/// (its distinguishing tail differs in the surviving bytes).
pub fn vj_triplet_alias_refdata() -> RefDataConfig {
    let mut cfg = RefDataConfig::empty(ChainType::Vj);
    for name in ["v_triplet*01", "v_triplet*02", "v_triplet*03"] {
        let _ = cfg.v_pool.push(Allele {
            name: name.into(),
            gene: "v_triplet".into(),
            seq: b"AAACCCGGGTTT".to_vec(),
            segment: Segment::V,
            anchor: Some(9),
        });
    }
    let _ = cfg.v_pool.push(Allele {
        name: "v_unique*01".into(),
        gene: "v_unique".into(),
        seq: b"AAACCCGGGAAA".to_vec(),
        segment: Segment::V,
        anchor: Some(9),
    });
    let _ = cfg.j_pool.push(Allele {
        name: "j_single*01".into(),
        gene: "j_single".into(),
        seq: b"TTTAAACCC".to_vec(),
        segment: Segment::J,
        anchor: Some(0),
    });
    cfg
}
