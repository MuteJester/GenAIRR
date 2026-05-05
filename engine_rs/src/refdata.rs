//! Reference data — immutable shared inputs to the simulation.
//!
//! ## What lives here
//!
//! Allele sequences, gene metadata, chain configuration. Loaded once
//! at simulator construction (Phase F: from `*.v6dat` files per D10),
//! shared across every simulation that uses it. *Never mutated* by the
//! simulator; the per-simulation state lives on `AlleleInstance`
//! (Phase C.4).
//!
//! Empirical distributions over alleles (frequencies) and trim/NP
//! distributions live in `dist.rs` rather than here — they are sampled
//! against the reference data, not part of it.
//!
//! ## Phase C.1 scope
//!
//! Just the typed data shapes plus construction / lookup helpers.
//! No biology, no PyO3, no serde. Phase F will add `serde` derives so
//! the types round-trip through `bincode` (D10).

use crate::ir::Segment;

// ──────────────────────────────────────────────────────────────────
// AlleleId — typed u32 newtype, distinct from NucHandle / RegionHandle
// ──────────────────────────────────────────────────────────────────

/// Index into an `AllelePool`. Stable for the lifetime of the
/// `RefDataConfig` it was issued for.
///
/// Distinct type from `NucHandle` and `RegionHandle` so the compiler
/// catches handle-confusion at call sites — this is the same
/// discipline as the IR handles in §3 of the design doc.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
pub struct AlleleId(u32);

impl AlleleId {
    pub const fn new(idx: u32) -> Self {
        Self(idx)
    }
    pub const fn index(self) -> u32 {
        self.0
    }
    pub const fn as_usize(self) -> usize {
        self.0 as usize
    }
}

// ──────────────────────────────────────────────────────────────────
// ChainType — VJ (light) vs VDJ (heavy)
// ──────────────────────────────────────────────────────────────────

/// The biological chain configuration.
///
/// `Vj` chains (light: kappa, lambda, TCR alpha/gamma) recombine V
/// and J only — there is no D segment, NP1 spans the V→J interval,
/// and NP2 / D pools are unused.
///
/// `Vdj` chains (heavy: IGH, TCR beta/delta) recombine V, D, and J
/// with NP1 between V and D and NP2 between D and J.
///
/// Other axes (allele frequency models, isotype constants, etc.)
/// are independent of `ChainType` and live elsewhere in the config.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
pub enum ChainType {
    Vj,
    Vdj,
}

impl ChainType {
    /// Whether this chain has a D segment (and therefore an NP2 region).
    pub const fn has_d(self) -> bool {
        matches!(self, ChainType::Vdj)
    }
}

// ──────────────────────────────────────────────────────────────────
// Allele — one germline allele entry
// ──────────────────────────────────────────────────────────────────

/// One germline allele in a reference data set.
///
/// Immutable once constructed. Per-simulation state (which allele was
/// sampled, what trim was applied, what ambiguity set the post-trim
/// retained bases project to) lives on `AlleleInstance` (Phase C.4).
///
/// **Field discipline:**
/// - `seq` is uppercase (`b'A'`, `b'C'`, `b'G'`, `b'T'`). Mixed case
///   is not allowed — case carries semantic meaning later in the
///   simulation (lowercase = SHM-mutated, etc.) and should not leak
///   into reference data.
/// - `anchor` is `Some(pos)` when the allele has a conserved codon
///   (Cys for V, W/F for J) at position `pos` (0-indexed within
///   `seq`). `None` for anchorless / partial alleles. Phase C.9's
///   `AnchorPreserved` contract reads this field.
/// - `name` is the canonical allele identifier (e.g.,
///   `"IGHV1-2*01"`). `gene` is the truncation to the gene level
///   (e.g., `"IGHV1-2"`).
#[derive(Clone, Debug)]
pub struct Allele {
    pub name: String,
    pub gene: String,
    pub seq: Vec<u8>,
    pub segment: Segment,
    pub anchor: Option<u16>,
}

impl Allele {
    /// Length of the allele sequence in nucleotides.
    pub fn len(&self) -> u32 {
        self.seq.len() as u32
    }

    /// Whether the allele has a known anchor position.
    pub fn has_anchor(&self) -> bool {
        self.anchor.is_some()
    }
}

// ──────────────────────────────────────────────────────────────────
// AllelePool — flat indexed collection of alleles for one segment
// ──────────────────────────────────────────────────────────────────

/// All alleles available for one segment role (V, D, J, or C).
///
/// Indexed by `AlleleId`. The pool is the reference data; *which*
/// allele a particular simulation sampled is recorded by an
/// `AlleleInstance` referring to a stable `AlleleId` (Phase C.4).
#[derive(Clone, Debug, Default)]
pub struct AllelePool {
    alleles: Vec<Allele>,
}

impl AllelePool {
    pub fn new() -> Self {
        Self::default()
    }

    /// Construct from an existing `Vec<Allele>`. Caller is responsible
    /// for ensuring all alleles have the right segment role for the
    /// pool being built.
    pub fn from_vec(alleles: Vec<Allele>) -> Self {
        Self { alleles }
    }

    /// Number of alleles in the pool.
    pub fn len(&self) -> usize {
        self.alleles.len()
    }

    /// Whether the pool contains zero alleles.
    pub fn is_empty(&self) -> bool {
        self.alleles.is_empty()
    }

    /// Append an allele. Returns the issued `AlleleId` so callers can
    /// reference the just-added allele in subsequent setup. This is
    /// the construction-time API; reference data is sealed before any
    /// simulation runs against it.
    #[must_use = "AllelePool::push returns the issued AlleleId; use \
                  `_ = pool.push(...)` if it isn't needed"]
    pub fn push(&mut self, allele: Allele) -> AlleleId {
        let id = AlleleId::new(self.alleles.len() as u32);
        self.alleles.push(allele);
        id
    }

    /// Look up an allele by id. Returns `None` if `id` is out of
    /// bounds — defensive to allow `RefDataConfig` to expose
    /// fallible lookup at the boundary even though correct callers
    /// always have valid handles.
    pub fn get(&self, id: AlleleId) -> Option<&Allele> {
        self.alleles.get(id.as_usize())
    }

    /// Read-only slice of all alleles in pool order. Phase C's
    /// `AllelePoolDist` (C.3) iterates this to build cumulative
    /// frequency tables.
    pub fn as_slice(&self) -> &[Allele] {
        &self.alleles
    }

    /// Iterator over `(AlleleId, &Allele)` pairs in pool order.
    pub fn iter(&self) -> impl Iterator<Item = (AlleleId, &Allele)> {
        self.alleles
            .iter()
            .enumerate()
            .map(|(i, a)| (AlleleId::new(i as u32), a))
    }

    /// Find the first allele whose name matches exactly. O(N).
    /// Used in tests and tooling; production sampling goes through
    /// `AllelePoolDist` (C.3) and never name-resolves here.
    pub fn find_by_name(&self, name: &str) -> Option<(AlleleId, &Allele)> {
        self.iter().find(|(_, a)| a.name == name)
    }
}

// ──────────────────────────────────────────────────────────────────
// RefDataConfig — top-level container for one species/chain
// ──────────────────────────────────────────────────────────────────

/// Top-level immutable reference configuration for one simulation
/// target (e.g., human IGH).
///
/// In production, this is loaded from a `.v6dat` file (D10) at
/// `Experiment.compile()` time. In tests, build directly via
/// `RefDataConfig::builder()` or by constructing the fields
/// manually.
///
/// **Invariant:** for `chain_type == ChainType::Vj`, the `d_pool` is
/// expected to be empty. Construction does not enforce this — callers
/// of the builder API in tests are expected to honor it; the
/// assembly pass in C.8 will explicitly skip D for VJ chains.
#[derive(Clone, Debug)]
pub struct RefDataConfig {
    pub chain_type: ChainType,
    pub v_pool: AllelePool,
    pub d_pool: AllelePool,
    pub j_pool: AllelePool,
    pub c_pool: AllelePool,
}

impl RefDataConfig {
    /// Empty config for the given chain type. Use the builder /
    /// direct field assignment to populate the pools.
    pub fn empty(chain_type: ChainType) -> Self {
        Self {
            chain_type,
            v_pool: AllelePool::new(),
            d_pool: AllelePool::new(),
            j_pool: AllelePool::new(),
            c_pool: AllelePool::new(),
        }
    }

    /// Resolve an allele by its segment role and id. Returns `None`
    /// if the segment isn't in this config's pools (e.g., asking for
    /// a D allele on a VJ chain) or if the id is out of bounds.
    pub fn get(&self, segment: Segment, id: AlleleId) -> Option<&Allele> {
        match segment {
            Segment::V => self.v_pool.get(id),
            Segment::D => self.d_pool.get(id),
            Segment::J => self.j_pool.get(id),
            // NP regions don't have alleles; return None defensively
            // rather than panicking — a caller asking for an NP
            // allele is misusing the API.
            Segment::Np1 | Segment::Np2 => None,
        }
    }

    /// Pool for the given segment, or `None` for NP segments.
    pub fn pool_for(&self, segment: Segment) -> Option<&AllelePool> {
        match segment {
            Segment::V => Some(&self.v_pool),
            Segment::D => Some(&self.d_pool),
            Segment::J => Some(&self.j_pool),
            Segment::Np1 | Segment::Np2 => None,
        }
    }
}

// ──────────────────────────────────────────────────────────────────
// Tests
// ──────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use std::mem::size_of;

    /// Helper: a tiny synthetic V allele for tests.
    fn make_v(name: &str, gene: &str, seq: &[u8], anchor: Option<u16>) -> Allele {
        Allele {
            name: name.to_string(),
            gene: gene.to_string(),
            seq: seq.to_vec(),
            segment: Segment::V,
            anchor,
        }
    }

    #[test]
    fn allele_id_is_zero_cost_newtype() {
        assert_eq!(size_of::<AlleleId>(), size_of::<u32>());
    }

    #[test]
    fn allele_id_round_trip() {
        let id = AlleleId::new(7);
        assert_eq!(id.index(), 7);
        assert_eq!(id.as_usize(), 7);
    }

    #[test]
    fn chain_type_has_d_distinction() {
        assert!(!ChainType::Vj.has_d());
        assert!(ChainType::Vdj.has_d());
    }

    #[test]
    fn allele_basic_accessors() {
        let a = make_v("IGHV1-2*01", "IGHV1-2", b"ACGTACGT", Some(3));
        assert_eq!(a.len(), 8);
        assert!(a.has_anchor());
        assert_eq!(a.anchor, Some(3));
        assert_eq!(a.segment, Segment::V);
    }

    #[test]
    fn allele_anchorless_round_trip() {
        let a = make_v("IGHV-pseudo*01", "IGHV-pseudo", b"ACGT", None);
        assert!(!a.has_anchor());
        assert_eq!(a.anchor, None);
    }

    #[test]
    fn allele_pool_starts_empty() {
        let p = AllelePool::new();
        assert_eq!(p.len(), 0);
        assert!(p.is_empty());
        assert!(p.get(AlleleId::new(0)).is_none());
        assert!(p.find_by_name("nonexistent").is_none());
    }

    #[test]
    fn allele_pool_push_returns_sequential_ids() {
        let mut p = AllelePool::new();
        let id0 = p.push(make_v("a*01", "a", b"AA", Some(0)));
        let id1 = p.push(make_v("b*01", "b", b"CC", Some(0)));
        let id2 = p.push(make_v("c*01", "c", b"GG", None));

        assert_eq!(id0.index(), 0);
        assert_eq!(id1.index(), 1);
        assert_eq!(id2.index(), 2);
        assert_eq!(p.len(), 3);
    }

    #[test]
    fn allele_pool_get_returns_stored_allele() {
        let mut p = AllelePool::new();
        let id = p.push(make_v("x*01", "x", b"ATGC", Some(2)));
        let got = p.get(id).expect("just pushed");
        assert_eq!(got.name, "x*01");
        assert_eq!(got.seq, b"ATGC");
        assert_eq!(got.anchor, Some(2));
    }

    #[test]
    fn allele_pool_get_out_of_bounds_returns_none() {
        let p = AllelePool::new();
        assert!(p.get(AlleleId::new(99)).is_none());
    }

    #[test]
    fn allele_pool_iter_yields_id_allele_pairs_in_order() {
        let mut p = AllelePool::new();
        let _ = p.push(make_v("a*01", "a", b"AA", None));
        let _ = p.push(make_v("b*01", "b", b"CC", None));
        let _ = p.push(make_v("c*01", "c", b"GG", None));

        let collected: Vec<(u32, String)> = p
            .iter()
            .map(|(id, a)| (id.index(), a.name.clone()))
            .collect();
        assert_eq!(
            collected,
            vec![
                (0, "a*01".to_string()),
                (1, "b*01".to_string()),
                (2, "c*01".to_string()),
            ]
        );
    }

    #[test]
    fn allele_pool_find_by_name_locates_exact_match() {
        let mut p = AllelePool::new();
        let _ = p.push(make_v("IGHV1-2*01", "IGHV1-2", b"AA", None));
        let target_id = p.push(make_v("IGHV1-2*02", "IGHV1-2", b"AC", None));
        let _ = p.push(make_v("IGHV3-23*01", "IGHV3-23", b"GG", None));

        let (id, allele) = p.find_by_name("IGHV1-2*02").expect("name should exist");
        assert_eq!(id, target_id);
        assert_eq!(allele.seq, b"AC");

        // Partial / similar names should not match.
        assert!(p.find_by_name("IGHV1-2").is_none());
        assert!(p.find_by_name("IGHV1-2*03").is_none());
    }

    #[test]
    fn ref_data_config_empty_for_chain_type() {
        let cfg = RefDataConfig::empty(ChainType::Vdj);
        assert_eq!(cfg.chain_type, ChainType::Vdj);
        assert!(cfg.v_pool.is_empty());
        assert!(cfg.d_pool.is_empty());
        assert!(cfg.j_pool.is_empty());
        assert!(cfg.c_pool.is_empty());
    }

    #[test]
    fn ref_data_config_pool_for_segment_routes_correctly() {
        let mut cfg = RefDataConfig::empty(ChainType::Vdj);
        let _ = cfg.v_pool.push(make_v("v*01", "v", b"AA", None));
        let _ = cfg.d_pool.push(Allele {
            name: "d*01".into(),
            gene: "d".into(),
            seq: b"GG".to_vec(),
            segment: Segment::D,
            anchor: None,
        });
        let _ = cfg.j_pool.push(Allele {
            name: "j*01".into(),
            gene: "j".into(),
            seq: b"TT".to_vec(),
            segment: Segment::J,
            anchor: Some(0),
        });

        assert_eq!(cfg.pool_for(Segment::V).unwrap().len(), 1);
        assert_eq!(cfg.pool_for(Segment::D).unwrap().len(), 1);
        assert_eq!(cfg.pool_for(Segment::J).unwrap().len(), 1);
        assert!(cfg.pool_for(Segment::Np1).is_none());
        assert!(cfg.pool_for(Segment::Np2).is_none());
    }

    #[test]
    fn ref_data_config_get_resolves_segment_id_pair() {
        let mut cfg = RefDataConfig::empty(ChainType::Vdj);
        let v_id = cfg.v_pool.push(make_v("v*01", "v", b"AAAT", Some(1)));

        let v = cfg.get(Segment::V, v_id).expect("v*01 should resolve");
        assert_eq!(v.name, "v*01");
        assert_eq!(v.anchor, Some(1));

        // Wrong segment -> None.
        assert!(cfg.get(Segment::J, v_id).is_none());
        // NP segment -> None defensively.
        assert!(cfg.get(Segment::Np1, v_id).is_none());
    }

    #[test]
    fn ref_data_config_supports_vj_chain_with_empty_d_pool() {
        // VJ chains: d_pool is conventionally empty. Construction
        // does not enforce this; assembly (C.8) handles VJ vs VDJ
        // explicitly via chain_type.has_d().
        let mut cfg = RefDataConfig::empty(ChainType::Vj);
        let _ = cfg.v_pool.push(make_v("v*01", "v", b"AA", Some(0)));
        let _ = cfg.j_pool.push(Allele {
            name: "j*01".into(),
            gene: "j".into(),
            seq: b"TT".to_vec(),
            segment: Segment::J,
            anchor: Some(0),
        });

        assert!(!cfg.chain_type.has_d());
        assert!(cfg.d_pool.is_empty());
        assert_eq!(cfg.v_pool.len(), 1);
        assert_eq!(cfg.j_pool.len(), 1);
    }
}
