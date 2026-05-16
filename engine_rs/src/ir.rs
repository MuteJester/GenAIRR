//! Intermediate-representation (IR) type definitions for the
//! simulation engine.
//!
//! This module defines the typed data shapes from §3 of the design
//! document. It deliberately contains **no biological logic** — only
//! the structural definitions of the entities. Logic lives in
//! sibling modules (passes, contracts, queries).
//!
//! ## Architectural commitments reflected here
//!
//! - **Arena allocation (D-§9):** `NucleotidePool` is a flat `Vec<Nucleotide>`.
//!   Cross-entity references are typed `u32` newtype handles, not pointers.
//! - **Persistent IR (D1):** all top-level entities derive `Clone` and have
//!   no interior mutability.
//! - **Entity-attached metadata (D5):** derived state lives on the entity
//!   that owns it (codon rail / amino acids).
//!
//! ## Performance / cost model (be honest about this)
//!
//! D1 commits to a *persistent IR contract* — every `with_*` method
//! takes `&self`, returns a new value, leaves the receiver intact.
//! The contract is honored exactly. The *implementation* is currently
//! **naive deep clone**: every `with_*` call clones the affected
//! `Vec` (the nucleotide pool, the regions list) wholesale.
//!
//! This means:
//!
//! - Cost per write is **O(N)** in pool size, not O(log N) or O(1).
//! - For a typical simulation (~400 nucleotides × ~50–60 revisions)
//!   the clone overhead is ~192 KB live history per active sim —
//!   acceptable for development.
//! - Branching (contracts that retry) multiplies this linearly
//!   with branch depth. We should measure before committing to the
//!   naive path long-term.
//! - The structural-sharing optimization (`Arc<Vec<Nucleotide>>` +
//!   copy-on-write at the whole-vector granularity, or an HAMT-based
//!   persistent vector via the `im` crate) is a forward-compatible
//!   change — it preserves the existing `with_*` API surface and
//!   tightens the cost model to ~O(1) for unchanged regions.
//!
//! When the contract test suite starts driving real
//! simulation patterns we'll measure the cost and decide whether to
//! land structure sharing. Until then, the API is honest about the
//! contract; this comment is honest about the cost.

// ──────────────────────────────────────────────────────────────────
// Handle types — typed u32 newtypes
//
// Distinct types so that a `NucHandle` and a `RegionHandle` cannot be
// confused at a call site. Both wrap `u32`. Index-based references are
// our answer to the "linked structures fight the borrow checker"
// problem — see §3 of the design doc.
// ──────────────────────────────────────────────────────────────────

/// Index into a `NucleotidePool`. Stable for the lifetime of the
/// pool revision it was issued for.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
pub struct NucHandle(u32);

impl NucHandle {
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

/// Index into a `Sequence`'s regions list. Stable for the lifetime of
/// the sequence revision it was issued for.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
pub struct RegionHandle(u32);

impl RegionHandle {
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
// Segment — the biological role a nucleotide participates in
// ──────────────────────────────────────────────────────────────────

/// The biological segment a nucleotide or region belongs to.
///
/// Order matters for assembly: V → NP1 → D → NP2 → J for VDJ chains;
/// V → NP1 → J for VJ chains (NP1 spans the V-J interval and there
/// are no D / NP2 entries).
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
pub enum Segment {
    V,
    Np1,
    D,
    Np2,
    J,
}

// ──────────────────────────────────────────────────────────────────
// NucFlags — per-nucleotide bitflags
//
// Hand-rolled bitset over u8. Constants exposed as `pub const`. We
// don't pull in the `bitflags` crate at this stage; if the flag set
// grows past 8 distinct kinds we revisit.
// ──────────────────────────────────────────────────────────────────

/// Per-nucleotide flags. Bitset over u8 — see `flag` constants.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug, Default)]
pub struct NucFlags(u8);

impl NucFlags {
    /// Empty flag set.
    pub const fn empty() -> Self {
        Self(0)
    }

    /// Construct from a raw u8 bitfield.
    pub const fn from_bits(bits: u8) -> Self {
        Self(bits)
    }

    /// Raw bits (for serialization / introspection).
    pub const fn bits(self) -> u8 {
        self.0
    }

    /// Test whether all bits in `flags` are set.
    pub const fn contains(self, flags: NucFlags) -> bool {
        (self.0 & flags.0) == flags.0
    }

    /// Set the bits in `flags`. Returns the new flags value.
    #[must_use]
    pub const fn with(self, flags: NucFlags) -> Self {
        Self(self.0 | flags.0)
    }

    /// Clear the bits in `flags`. Returns the new flags value.
    #[must_use]
    pub const fn without(self, flags: NucFlags) -> Self {
        Self(self.0 & !flags.0)
    }
}

/// Flag constants. Each is a single-bit `NucFlags` value.
pub mod flag {
    use super::NucFlags;

    /// Nucleotide is a P-nucleotide (templated palindromic complement
    /// of the segment edge during V(D)J recombination).
    pub const P_NUC: NucFlags = NucFlags::from_bits(1 << 0);
    /// Nucleotide is an N-nucleotide (TdT-added at the junction).
    pub const N_NUC: NucFlags = NucFlags::from_bits(1 << 1);
    /// Nucleotide sits inside the junction region (V Cys → J W/F + 3).
    pub const JUNCTION: NucFlags = NucFlags::from_bits(1 << 2);
    /// Nucleotide came from a D segment that was inverted (RSS-symmetric
    /// V(D)J inversion event, biologically real for D segments only).
    pub const INVERTED: NucFlags = NucFlags::from_bits(1 << 3);
    /// Nucleotide was inserted by a downstream indel pass — has no
    /// germline-allele provenance.
    pub const INDEL_INSERTED: NucFlags = NucFlags::from_bits(1 << 4);
}

// ──────────────────────────────────────────────────────────────────
// GermlinePos — typed source-allele position with `NONE` for synthetic
// ──────────────────────────────────────────────────────────────────

/// Position in the source allele, or `NONE` for synthetic bases
/// (NP, P-nuc, contaminant, indel-inserted) with no germline
/// provenance. Replaces the previous `germline_pos: u16` field
/// where `u16::MAX` was a magic sentinel.
///
/// `#[repr(transparent)]` keeps the in-memory layout identical to a
/// raw `u16`, so this type adds zero memory cost on `Nucleotide`.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
#[repr(transparent)]
pub struct GermlinePos(u16);

impl GermlinePos {
    /// Sentinel value meaning "no germline provenance".
    pub const NONE: Self = Self(u16::MAX);

    /// Construct a real allele position. Panics if `v == u16::MAX`
    /// — that value is reserved for `NONE`.
    pub const fn pos(v: u16) -> Self {
        assert!(
            v != u16::MAX,
            "GermlinePos::pos: u16::MAX is reserved for NONE"
        );
        Self(v)
    }

    /// Project to `Option<u16>`. `Some(v)` for a real position,
    /// `None` for `NONE`.
    pub const fn get(self) -> Option<u16> {
        if self.0 == u16::MAX {
            None
        } else {
            Some(self.0)
        }
    }

    /// Returns `true` when this is the `NONE` sentinel.
    pub const fn is_none(self) -> bool {
        self.0 == u16::MAX
    }

    /// Returns `true` when this carries a real allele position.
    pub const fn is_some(self) -> bool {
        !self.is_none()
    }
}

// ──────────────────────────────────────────────────────────────────
// Nucleotide — single base with provenance metadata
// ──────────────────────────────────────────────────────────────────

/// A single nucleotide in the simulation IR.
///
/// Stored in a `NucleotidePool` arena and referenced via `NucHandle`.
/// Layout is `Copy`-able for cheap snapshots — total size is currently
/// 8 bytes, target is to keep this <= 16 bytes.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
pub struct Nucleotide {
    /// Current base — usually one of `b'A'`, `b'C'`, `b'G'`, `b'T'`
    /// (uppercase for germline-derived, lowercase for SHM-mutated /
    /// observation-time edits) plus `b'N'`/`b'n'` for ambiguous /
    /// quality-masked bases.
    pub base: u8,

    /// Wild-type base at this position — the germline reference. Equal
    /// to `base` until a mutation pass changes `base`.
    pub germline: u8,

    /// Position in the source allele (0-indexed), or `GermlinePos::NONE`
    /// for NP / P-nuc / contaminant / indel-inserted bases with no
    /// germline provenance. See `GermlinePos` for the safe accessor API.
    pub germline_pos: GermlinePos,

    /// Biological segment role.
    pub segment: Segment,

    /// Per-nucleotide flags (see `flag` constants).
    pub flags: NucFlags,
}

impl Nucleotide {
    /// Construct a germline-derived nucleotide (base == germline).
    /// `germline_pos` must be a real allele position (`< u16::MAX`);
    /// passing `u16::MAX` panics — use `Nucleotide::synthetic` for
    /// bases with no germline provenance.
    pub const fn germline(base: u8, germline_pos: u16, segment: Segment) -> Self {
        Self {
            base,
            germline: base,
            germline_pos: GermlinePos::pos(germline_pos),
            segment,
            flags: NucFlags::empty(),
        }
    }

    /// Construct a synthetic (NP / P-nuc / inserted) nucleotide with no
    /// germline provenance. The `flags` argument carries the kind tag
    /// (e.g., `flag::N_NUC`).
    pub const fn synthetic(base: u8, segment: Segment, flags: NucFlags) -> Self {
        Self {
            base,
            germline: base,
            germline_pos: GermlinePos::NONE,
            segment,
            flags,
        }
    }
}

// ──────────────────────────────────────────────────────────────────
// NucleotidePool — flat arena of nucleotides
//
// Owned by the Simulation. NucHandles index into `nucleotides`.
// ──────────────────────────────────────────────────────────────────

/// Arena of all nucleotides in a `Simulation`. `NucHandle` indexes
/// into this pool.
#[derive(Clone, Debug, Default)]
pub struct NucleotidePool {
    nucleotides: Vec<Nucleotide>,
}

impl NucleotidePool {
    /// New empty pool.
    pub fn new() -> Self {
        Self::default()
    }

    /// New empty pool pre-allocated for `cap` nucleotides. Use this at
    /// `Simulation` construction to avoid reallocations during a single
    /// simulation's lifetime.
    pub fn with_capacity(cap: usize) -> Self {
        Self {
            nucleotides: Vec::with_capacity(cap),
        }
    }

    /// Number of nucleotides currently in the pool.
    pub fn len(&self) -> usize {
        self.nucleotides.len()
    }

    /// Whether the pool contains zero nucleotides.
    pub fn is_empty(&self) -> bool {
        self.nucleotides.is_empty()
    }

    /// Read-only access to a nucleotide by handle.
    pub fn get(&self, h: NucHandle) -> Option<&Nucleotide> {
        self.nucleotides.get(h.as_usize())
    }

    /// Read-only access to the underlying slice (for bulk iteration
    /// inside metadata computation).
    pub fn as_slice(&self) -> &[Nucleotide] {
        &self.nucleotides
    }

    /// Append a nucleotide to the pool, returning its handle.
    /// Internal — the user-facing path goes through the persistent
    /// `Simulation::with_*` API. Kept `pub` for use during initial
    /// construction (passes build a fresh sim and push nucleotides
    /// directly into a fresh pool before sealing it).
    pub fn push(&mut self, n: Nucleotide) -> NucHandle {
        let h = NucHandle::new(self.nucleotides.len() as u32);
        self.nucleotides.push(n);
        h
    }

    // ── Persistent update API ──────────────────────────
    //
    // Each `with_*` method takes `&self`, returns a new `NucleotidePool`,
    // and never mutates the receiver. This is the persistent-IR
    // contract from D1.
    //
    // **Implementation note:** the current implementation uses naive
    // deep clone (Vec::clone()) for simplicity. Structure sharing
    // (Arc<Vec<…>> with copy-on-write, or `im::Vector`-style HAMT) is
    // a future optimization — the API surface is forwards-compatible
    // with that change.

    /// Return a new pool with `n` appended; receiver unchanged.
    ///
    /// **API note:** this method returns a `(Self, NucHandle)` tuple,
    /// not a plain `Self`. The asymmetry vs `with_base_changed` /
    /// `with_nucleotide_changed` (which return `Self` only) is
    /// deliberate: an *appending* operation creates a *new* handle,
    /// and the caller almost always needs it to refer to the just-
    /// added nucleotide in subsequent updates. Returning the handle
    /// here avoids forcing every caller to recompute it from
    /// `pool.len() - 1`. `#[must_use]` on the return type catches
    /// silent drops at compile time.
    ///
    /// Idiomatic call shape:
    /// ```ignore
    /// let (pool2, h) = pool1.with_pushed(nucleotide);
    /// // use both pool2 and h
    /// ```
    #[must_use = "with_pushed returns (NucleotidePool, NucHandle); both must \
                  be used or destructured. Drop the handle explicitly with \
                  `let (next, _) = ...` if it isn't needed."]
    pub fn with_pushed(&self, n: Nucleotide) -> (Self, NucHandle) {
        let mut next = self.clone();
        let h = next.push(n);
        (next, h)
    }

    /// Return a new pool with the nucleotide at `handle` replaced by
    /// `new_n`. Panics if `handle` is out of bounds — callers are
    /// responsible for handle validity.
    pub fn with_nucleotide_changed(&self, handle: NucHandle, new_n: Nucleotide) -> Self {
        let mut next = self.clone();
        next.nucleotides[handle.as_usize()] = new_n;
        next
    }

    /// Return a new pool with only the `base` field of the nucleotide
    /// at `handle` changed to `new_base`. All other fields (germline,
    /// segment, flags, germline_pos) are preserved. Panics if
    /// `handle` is out of bounds.
    pub fn with_base_changed(&self, handle: NucHandle, new_base: u8) -> Self {
        let mut next = self.clone();
        next.nucleotides[handle.as_usize()].base = new_base;
        next
    }

    // ── Indel API ──────────────────────────────────────
    //
    // Insertions and deletions change pool length and *shift* the
    // positional handles of all nucleotides at-or-after the indel
    // point. Old handles into the pool are NOT stable across an
    // indel — callers that hold handles before an indel must
    // remap them through `shift_handle_for_indel` if they want to
    // keep referring to the same logical nucleotide.

    /// Return a new pool with a fresh nucleotide inserted at
    /// position `at` — the new nucleotide takes that index and
    /// every nucleotide formerly at index ≥ `at` shifts up by one.
    ///
    /// `at == self.len()` is allowed and equivalent to appending
    /// at the end. Panics if `at > self.len()`.
    pub fn with_inserted(&self, at: u32, n: Nucleotide) -> Self {
        assert!(
            (at as usize) <= self.nucleotides.len(),
            "NucleotidePool::with_inserted: at ({}) > pool length ({})",
            at,
            self.nucleotides.len()
        );
        let mut next = self.clone();
        next.nucleotides.insert(at as usize, n);
        next
    }

    /// Return a new pool with the nucleotide at `at` removed —
    /// every nucleotide formerly at index > `at` shifts down by
    /// one. Panics if `at >= self.len()`.
    pub fn with_deleted(&self, at: u32) -> Self {
        assert!(
            (at as usize) < self.nucleotides.len(),
            "NucleotidePool::with_deleted: at ({}) >= pool length ({})",
            at,
            self.nucleotides.len()
        );
        let mut next = self.clone();
        next.nucleotides.remove(at as usize);
        next
    }
}

// ──────────────────────────────────────────────────────────────────
// Genetic code — codon → amino acid translation
//
// Standard NCBI translation table 1 (the universal code). Encoded as
// a 64-entry const array indexed by `(b1*16 + b2*4 + b3)` where each
// base is encoded as A=0, C=1, G=2, T=3.
//
// Stops are represented as `b'*'` in the output. Ambiguous codons
// (containing N or any unrecognized base) translate to `b'X'`.
// ──────────────────────────────────────────────────────────────────

/// Map a DNA base byte to its 0..=3 numeric encoding (A=0, C=1, G=2,
/// T=3). Returns `None` for any other byte (N, gap, ambiguity codes,
/// etc.). Case-insensitive. U is treated as T (RNA-tolerant).
///
/// Public so the S5F kernel can reuse the same encoding
/// to build 5-mer context indices.
pub const fn encode_base(b: u8) -> Option<u8> {
    match b.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' | b'U' => Some(3),
        _ => None,
    }
}

// Codon arithmetic (translate_codon, AMINO_STOP, AMINO_AMBIGUOUS,
// GENETIC_CODE) lives in `crate::codon`. Re-exported here for backwards
// compatibility — every existing `use crate::ir::translate_codon`
// callsite (feasibility, contracts, the codon-rail recompute below)
// keeps compiling unchanged.
pub use crate::codon::{translate_codon, AMINO_AMBIGUOUS, AMINO_STOP};

// ──────────────────────────────────────────────────────────────────
// Region — a contiguous run of nucleotides bound to a segment role
// ──────────────────────────────────────────────────────────────────

/// One region of the assembled sequence. Carries its own segment role,
/// its range of `NucHandle`s in the pool, and its derived codon-rail
/// metadata (D5 — entity-attached metadata, no separate cache).
///
/// The `amino_acids` and `stop_codon_positions` fields are computed
/// from the pool by `with_codon_rail_recomputed` and recomputed on
/// every IR revision that touches the region's bases. The persistent
/// IR contract (D1) means a stale field is structurally impossible —
/// any `Region` value either has up-to-date metadata or was constructed
/// before any metadata was needed.
#[derive(Clone, Debug)]
pub struct Region {
    /// Biological role of this region.
    pub segment: Segment,

    /// Half-open range of nucleotide handles `[start, end)` in the
    /// owning `NucleotidePool`. Empty range means the region exists
    /// but contributes no nucleotides (e.g., NP1 with zero N/P bases).
    pub start: NucHandle,
    pub end: NucHandle,

    /// Position within the codon frame at the first nucleotide of this
    /// region. `0` = this nucleotide is the first base of a codon, `1`
    /// = second base, `2` = third base. Determined by upstream regions'
    /// cumulative lengths during assembly.
    pub frame_phase: u8,

    /// Codon-rail metadata — translated amino acids for codons fully
    /// contained within this region, in 5'→3' order. Cross-region
    /// codons are not included here; sequence-level walks handle those.
    /// Each byte is either a single-letter amino-acid code,
    /// `AMINO_STOP` (`b'*'`) for a stop, or `AMINO_AMBIGUOUS`
    /// (`b'X'`) for a codon containing a non-{A,C,G,T,U} base.
    /// Empty until `with_codon_rail_recomputed` has been called.
    pub amino_acids: Vec<u8>,

    /// Handles of the *first* base of every stop codon in this region.
    /// Subset of `start..end`. Empty under the same conditions as
    /// `amino_acids`.
    pub stop_codon_positions: Vec<NucHandle>,
}

impl Region {
    /// Construct a region with the given segment and nucleotide range.
    /// `frame_phase` defaults to 0; assembly will set it from the
    /// preceding regions. Codon-rail metadata starts empty —
    /// `with_codon_rail_recomputed` populates it once the pool is
    /// available.
    pub fn new(segment: Segment, start: NucHandle, end: NucHandle) -> Self {
        Self {
            segment,
            start,
            end,
            frame_phase: 0,
            amino_acids: Vec::new(),
            stop_codon_positions: Vec::new(),
        }
    }

    /// Number of nucleotides in this region.
    pub fn len(&self) -> u32 {
        self.end.index().saturating_sub(self.start.index())
    }

    /// Whether the region has zero nucleotides.
    pub fn is_empty(&self) -> bool {
        self.start.index() == self.end.index()
    }

    // ── Persistent update API ──────────────────────────

    /// Return a new region with `end` advanced by one nucleotide.
    /// Used when a pass appends to the underlying pool.
    pub fn with_end_extended(&self, new_end: NucHandle) -> Self {
        Self {
            end: new_end,
            ..self.clone()
        }
    }

    /// Return a new region with the given frame phase. Used by
    /// assembly to set the phase based on cumulative upstream length.
    pub fn with_frame_phase(&self, phase: u8) -> Self {
        Self {
            frame_phase: phase,
            ..self.clone()
        }
    }

    // ── Codon-rail metadata ─────────────────────────────

    /// Return a new region with `amino_acids` and
    /// `stop_codon_positions` recomputed from `pool`. Receiver
    /// unchanged.
    ///
    /// Translation walks codons fully contained within the region's
    /// `[start, end)` range, taking `frame_phase` into account.
    /// Cross-region codons (where the third base of a codon spills
    /// into the next region) are NOT included here — they are
    /// handled at the sequence level when both regions are available.
    ///
    /// Encoded amino acids are single-letter codes from
    /// `GENETIC_CODE`; stops emit `AMINO_STOP` (`b'*'`); codons
    /// containing any non-{A,C,G,T,U} base emit `AMINO_AMBIGUOUS`
    /// (`b'X'`).
    pub fn with_codon_rail_recomputed(&self, pool: &NucleotidePool) -> Self {
        // Skip leading bases that complete a codon started in the
        // previous region:
        //   frame_phase=0 → first base of this region is codon[0]   → skip 0
        //   frame_phase=1 → first base of this region is codon[1]   → skip 2
        //                    (the next 2 bases finish the spanning codon)
        //   frame_phase=2 → first base of this region is codon[2]   → skip 1
        let skip = (3 - (self.frame_phase as u64)) % 3;

        // Loop bookkeeping in u64 to make `i + 3` overflow-impossible
        // even at the top of the u32 handle range. NucHandle is u32
        // so any in-bounds index fits in u32 trivially after the
        // addition; we only widen for the loop arithmetic.
        let start_idx_u64 = (self.start.index() as u64).saturating_add(skip);
        let end_idx_u64 = self.end.index() as u64;

        // Defensive: a malformed region with end < start (which
        // `len()` already saturates at zero) yields an empty rail
        // here too. No codons emitted, no panic.
        if start_idx_u64 >= end_idx_u64 {
            return Self {
                amino_acids: Vec::new(),
                stop_codon_positions: Vec::new(),
                ..self.clone()
            };
        }

        // Pre-allocate based on a tight upper bound so we don't reallocate.
        let max_codons = ((end_idx_u64 - start_idx_u64) / 3) as usize;
        let mut amino_acids = Vec::with_capacity(max_codons);
        let mut stops = Vec::new();

        let mut i = start_idx_u64;
        while i + 3 <= end_idx_u64 {
            // i + 2 < end_idx_u64 ≤ u32::MAX (NucHandle is u32), so
            // every cast back to u32 is lossless.
            let h0 = NucHandle::new(i as u32);
            // Out-of-bounds handles here would indicate a bug in the
            // caller's pool/region invariants. Use unwrap intentionally:
            // a panic is the right signal during development.
            let b1 = pool.get(NucHandle::new(i as u32)).unwrap().base;
            let b2 = pool.get(NucHandle::new((i + 1) as u32)).unwrap().base;
            let b3 = pool.get(NucHandle::new((i + 2) as u32)).unwrap().base;
            let aa = translate_codon(b1, b2, b3);
            amino_acids.push(aa);
            if aa == AMINO_STOP {
                stops.push(h0);
            }
            i += 3;
        }

        Self {
            amino_acids,
            stop_codon_positions: stops,
            ..self.clone()
        }
    }

    /// Convenience: number of stop codons in this region's
    /// already-computed metadata.
    pub fn stop_codon_count(&self) -> usize {
        self.stop_codon_positions.len()
    }
}

// ──────────────────────────────────────────────────────────────────
// Sequence — the assembled product, owns its regions
// ──────────────────────────────────────────────────────────────────

/// The assembled sequence — the root structural entity at the
/// "biological product" level.
#[derive(Clone, Debug, Default)]
pub struct Sequence {
    /// Regions in biological assembly order.
    pub regions: Vec<Region>,
}

impl Sequence {
    pub fn new() -> Self {
        Self::default()
    }

    /// Number of regions currently in the sequence.
    pub fn region_count(&self) -> usize {
        self.regions.len()
    }

    // ── Persistent update API ──────────────────────────

    /// Return a new sequence with `region` appended; receiver unchanged.
    pub fn with_region_added(&self, region: Region) -> Self {
        let mut next = self.clone();
        next.regions.push(region);
        next
    }

    /// Return a new sequence with the region at `idx` replaced by
    /// `region`. Panics if `idx` is out of bounds.
    pub fn with_region_replaced(&self, idx: usize, region: Region) -> Self {
        let mut next = self.clone();
        next.regions[idx] = region;
        next
    }

    /// Return a new sequence whose regions have `frame_phase`
    /// recomputed from biological order and cumulative upstream
    /// lengths.
    ///
    /// This is the canonical repair primitive after any operation
    /// that changes one region's length. A downstream region's pool
    /// range may be shifted correctly but still have stale codon-frame
    /// metadata unless its phase is recalculated from the adjusted
    /// upstream region lengths.
    pub fn with_frame_phases_recomputed(&self) -> Self {
        let mut cumulative_len = 0u64;
        let regions: Vec<Region> = self
            .regions
            .iter()
            .map(|r| {
                let next = r.with_frame_phase((cumulative_len % 3) as u8);
                cumulative_len = cumulative_len.saturating_add(next.len() as u64);
                next
            })
            .collect();
        Self { regions }
    }

    // ── Indel adjustment ───────────────────────────────

    /// Return a new sequence with every region's range adjusted for
    /// an indel at pool position `pos` with the given `delta`
    /// (`+1` for insertion, `-1` for deletion).
    ///
    /// Adjustment rules per region:
    /// - **Entirely before** the indel (`region.end <= pos`):
    ///   unchanged.
    /// - **Entirely after** (`region.start > pos`): both `start`
    ///   and `end` shift by `delta`.
    /// - **Spanning the indel** (`region.start <= pos < region.end`):
    ///   `start` unchanged, `end` shifts by `delta` (region grows
    ///   on insertion, shrinks on deletion).
    ///
    /// After range adjustment, `frame_phase` is recomputed for all
    /// regions from their adjusted lengths. Codon rails are *not*
    /// recomputed here — that's `Simulation::with_indel_inserted` /
    /// `with_indel_deleted`'s job, since rails need the new pool not
    /// the old.
    ///
    pub fn with_indel_adjusted(&self, pos: u32, delta: i32) -> Self {
        let new_regions: Vec<Region> = self
            .regions
            .iter()
            .map(|r| {
                let r_start = r.start.index();
                let r_end = r.end.index();

                if r_end <= pos {
                    // Entirely before — unchanged.
                    r.clone()
                } else if r_start > pos {
                    // Entirely after — shift both.
                    Region {
                        start: NucHandle::new(shift_pos(r_start, delta)),
                        end: NucHandle::new(shift_pos(r_end, delta)),
                        ..r.clone()
                    }
                } else {
                    // Spans the indel — extend / shrink.
                    Region {
                        end: NucHandle::new(shift_pos(r_end, delta)),
                        ..r.clone()
                    }
                }
            })
            .collect();
        Sequence {
            regions: new_regions,
        }
        .with_frame_phases_recomputed()
    }
}

/// Shift a position by a signed delta. Saturates at zero on
/// negative deltas — callers should ensure deltas are valid
/// (e.g., a deletion never shifts a position below zero).
fn shift_pos(pos: u32, delta: i32) -> u32 {
    if delta >= 0 {
        pos.saturating_add(delta as u32)
    } else {
        pos.saturating_sub((-delta) as u32)
    }
}

// ──────────────────────────────────────────────────────────────────
// refresh_regions_covering — internal helper for codon-rail consistency
// ──────────────────────────────────────────────────────────────────

/// For every region whose half-open range covers `handle`, return
/// a copy with its codon rail recomputed against `pool`. Other
/// regions are kept unchanged (cloned by-value, but their metadata
/// is still valid because they don't include the mutated handle).
///
/// Internal helper for `Simulation::with_base_changed` and
/// `with_nucleotide_changed`. Without this, those methods would
/// produce IR revisions whose `Region.amino_acids` lies about
/// the pool's actual content — a subtle staleness bug that would
/// surface the moment SHM mutations enter the picture.
fn refresh_regions_covering(seq: &Sequence, pool: &NucleotidePool, handle: NucHandle) -> Sequence {
    let h = handle.index();
    let new_regions: Vec<Region> = seq
        .regions
        .iter()
        .map(|r| {
            if h >= r.start.index() && h < r.end.index() {
                r.with_codon_rail_recomputed(pool)
            } else {
                r.clone()
            }
        })
        .collect();
    Sequence {
        regions: new_regions,
    }
}

// ──────────────────────────────────────────────────────────────────
// Simulation — the root entity, owns the pool and the sequence
// ──────────────────────────────────────────────────────────────────

/// The root of one simulation. Owns the nucleotide pool, the
/// assembled sequence, and the per-simulation allele assignments.
/// `live_calls` is the dynamic allele-call evidence sidecar. Early
/// phases populate exact assembled-region calls from the compiled
/// simulator; AIRR call projection reads V/D/J call sets from it
/// when present.
/// Future phases add: pass history (already collected by the
/// runtime in `Outcome`), contract set, RNG state.
#[derive(Clone, Debug, Default)]
pub struct Simulation {
    /// Arena of all nucleotides.
    pub pool: NucleotidePool,

    /// The assembled sequence (initially empty).
    pub sequence: Sequence,

    /// V/D/J/C allele instances sampled for this simulation.
    /// Pre-recombination this is all-`None`. Each sampling pass
    /// (C.5) populates one slot; trim passes (C.6) update the
    /// instance via `with_trim`. See `assignment.rs` for the type.
    pub assignments: crate::assignment::AlleleAssignments,

    /// Dynamic V/D/J call evidence over the current observed sequence.
    /// Optional because direct/non-compiled execution can still build
    /// structural simulations without the live-call index.
    pub live_calls: Option<crate::live_call::LiveCallState>,
}

impl Simulation {
    /// Empty simulation.
    pub fn new() -> Self {
        Self::default()
    }

    /// Empty simulation with `nuc_capacity` reserved up front.
    pub fn with_capacity(nuc_capacity: usize) -> Self {
        Self {
            pool: NucleotidePool::with_capacity(nuc_capacity),
            sequence: Sequence::new(),
            assignments: crate::assignment::AlleleAssignments::new(),
            live_calls: None,
        }
    }

    /// Return a new simulation with dynamic live-call evidence set.
    ///
    /// The compiled runtime uses this to persist live-call updates
    /// after pass effects.
    pub fn with_live_calls(&self, live_calls: crate::live_call::LiveCallState) -> Self {
        Self {
            pool: self.pool.clone(),
            sequence: self.sequence.clone(),
            assignments: self.assignments,
            live_calls: Some(live_calls),
        }
    }

    // ── Persistent update API ──────────────────────────
    //
    // Every method here takes `&self`, returns a new `Simulation`, and
    // never mutates the receiver. This is the public surface that
    // passes use to evolve the simulation without ever mutating an
    // existing IR revision.

    /// Return a new simulation with `n` pushed onto the pool. The
    /// returned tuple gives the issued `NucHandle` so the caller can
    /// reference the new nucleotide in subsequent updates.
    ///
    /// See `NucleotidePool::with_pushed` for the rationale behind
    /// the tuple return — it's the only `with_*` method on
    /// `Simulation` that does not return `Self` alone, and it does
    /// so for the same reason: an appending operation creates a
    /// new handle that the caller almost always needs.
    #[must_use = "with_nucleotide_pushed returns (Simulation, NucHandle); \
                  both must be used or destructured. Drop the handle \
                  explicitly with `let (next, _) = ...` if it isn't needed."]
    pub fn with_nucleotide_pushed(&self, n: Nucleotide) -> (Self, NucHandle) {
        let (pool, h) = self.pool.with_pushed(n);
        (
            Self {
                pool,
                sequence: self.sequence.clone(),
                assignments: self.assignments,
                live_calls: self.live_calls.clone(),
            },
            h,
        )
    }

    /// Return a new simulation with the nucleotide at `handle`
    /// replaced. Panics if `handle` is out of bounds.
    ///
    /// **Codon-rail consistency (D5 invariant):** if the changed
    /// nucleotide lies inside any assembled region, that region's
    /// codon rail is automatically recomputed against the new pool.
    /// This preserves the entity-attached metadata invariant
    /// promised by D5: a `Region` value either has up-to-date
    /// metadata or has not been recomputed yet — staleness against
    /// the pool the simulation actually carries is impossible.
    pub fn with_nucleotide_changed(&self, handle: NucHandle, new_n: Nucleotide) -> Self {
        let new_pool = self.pool.with_nucleotide_changed(handle, new_n);
        let new_sequence = refresh_regions_covering(&self.sequence, &new_pool, handle);
        Self {
            pool: new_pool,
            sequence: new_sequence,
            assignments: self.assignments,
            live_calls: self.live_calls.clone(),
        }
    }

    /// Return a new simulation with only the base of the nucleotide
    /// at `handle` changed to `new_base`. Panics if `handle` is out
    /// of bounds. This is the canonical mutation primitive that
    /// SHM passes compose against.
    ///
    /// **Codon-rail consistency (D5 invariant):** the affected
    /// region (if any contains `handle`) has its codon rail
    /// recomputed automatically. See `with_nucleotide_changed` for
    /// the same discipline applied to whole-nucleotide replacement.
    pub fn with_base_changed(&self, handle: NucHandle, new_base: u8) -> Self {
        let new_pool = self.pool.with_base_changed(handle, new_base);
        let new_sequence = refresh_regions_covering(&self.sequence, &new_pool, handle);
        Self {
            pool: new_pool,
            sequence: new_sequence,
            assignments: self.assignments,
            live_calls: self.live_calls.clone(),
        }
    }

    /// Return a new simulation with `region` added to the sequence.
    pub fn with_region_added(&self, region: Region) -> Self {
        Self {
            pool: self.pool.clone(),
            sequence: self.sequence.with_region_added(region),
            assignments: self.assignments,
            live_calls: self.live_calls.clone(),
        }
    }

    // ── Allele-assignment persistent API (C.4) ─────────────────────

    /// Return a new simulation with `instance` set on the allele
    /// slot for `segment` (V, D, J, or C). Panics if `segment` is
    /// an NP segment.
    pub fn with_allele_assigned(
        &self,
        segment: Segment,
        instance: crate::assignment::AlleleInstance,
    ) -> Self {
        Self {
            pool: self.pool.clone(),
            sequence: self.sequence.clone(),
            assignments: self.assignments.with_assigned(segment, instance),
            live_calls: self.live_calls.clone(),
        }
    }

    /// Return a new simulation with the trim of the assigned allele
    /// at `segment` updated to `value` on the given `end`. Panics
    /// if no allele is currently assigned to that segment, or if
    /// `segment` is an NP segment.
    pub fn with_trim(&self, segment: Segment, end: crate::assignment::TrimEnd, value: u16) -> Self {
        Self {
            pool: self.pool.clone(),
            sequence: self.sequence.clone(),
            assignments: self.assignments.with_trim(segment, end, value),
            live_calls: self.live_calls.clone(),
        }
    }

    // ── Indel API ──────────────────────────────────────

    /// Return a new simulation with `n` inserted at pool position
    /// `at`. The pool grows by 1; every existing nucleotide at
    /// position ≥ `at` shifts up by 1. Region ranges are adjusted
    /// per `Sequence::with_indel_adjusted`. Codon rails on regions
    /// that span or follow the insertion are recomputed against
    /// the new pool, since their range changed (and the bases
    /// inside may have shifted).
    pub fn with_indel_inserted(&self, at: u32, n: Nucleotide) -> Self {
        let new_pool = self.pool.with_inserted(at, n);
        let adjusted = self.sequence.with_indel_adjusted(at, 1);
        // Recompute rails on regions whose range changed (i.e., they
        // span or follow the insertion point).
        let new_regions: Vec<Region> = adjusted
            .regions
            .iter()
            .enumerate()
            .map(|(i, r)| {
                let original = &self.sequence.regions[i];
                if r.start != original.start
                    || r.end != original.end
                    || r.frame_phase != original.frame_phase
                {
                    r.with_codon_rail_recomputed(&new_pool)
                } else {
                    r.clone()
                }
            })
            .collect();
        Self {
            pool: new_pool,
            sequence: Sequence {
                regions: new_regions,
            },
            assignments: self.assignments,
            live_calls: self.live_calls.clone(),
        }
    }

    /// Return a new simulation with the nucleotide at pool position
    /// `at` removed. The pool shrinks by 1; every existing
    /// nucleotide at position > `at` shifts down by 1. Region
    /// ranges are adjusted per `Sequence::with_indel_adjusted`.
    /// Codon rails on regions that span or follow the deletion are
    /// recomputed against the new pool.
    pub fn with_indel_deleted(&self, at: u32) -> Self {
        let new_pool = self.pool.with_deleted(at);
        let adjusted = self.sequence.with_indel_adjusted(at, -1);
        let new_regions: Vec<Region> = adjusted
            .regions
            .iter()
            .enumerate()
            .map(|(i, r)| {
                let original = &self.sequence.regions[i];
                if r.start != original.start
                    || r.end != original.end
                    || r.frame_phase != original.frame_phase
                {
                    r.with_codon_rail_recomputed(&new_pool)
                } else {
                    r.clone()
                }
            })
            .collect();
        Self {
            pool: new_pool,
            sequence: Sequence {
                regions: new_regions,
            },
            assignments: self.assignments,
            live_calls: self.live_calls.clone(),
        }
    }
}

// ──────────────────────────────────────────────────────────────────
// Tests
// ──────────────────────────────────────────────────────────────────


#[cfg(test)]
#[path = "ir_tests.rs"]
mod tests;
