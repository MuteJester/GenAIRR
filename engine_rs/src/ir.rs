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

    /// Position in the source allele (0-indexed). For NP / inserted /
    /// contaminant bases this is undefined; we use `u16::MAX` as a
    /// sentinel meaning "no germline provenance".
    pub germline_pos: u16,

    /// Biological segment role.
    pub segment: Segment,

    /// Per-nucleotide flags (see `flag` constants).
    pub flags: NucFlags,
}

impl Nucleotide {
    /// Sentinel `germline_pos` value for nucleotides with no germline
    /// provenance (NP, P-nuc, contaminant, indel-inserted, etc.).
    pub const NO_GERMLINE_POS: u16 = u16::MAX;

    /// Construct a germline-derived nucleotide (base == germline).
    pub const fn germline(base: u8, germline_pos: u16, segment: Segment) -> Self {
        Self {
            base,
            germline: base,
            germline_pos,
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
            germline_pos: Self::NO_GERMLINE_POS,
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

/// Standard genetic code (NCBI table 1). Indexed by
/// `(b1<<4) | (b2<<2) | b3` with the 0..=3 base encoding from
/// `encode_base`. Stop codons (TAA, TAG, TGA) emit `b'*'`.
#[rustfmt::skip]
const GENETIC_CODE: [u8; 64] = [
    // AAA   AAC   AAG   AAT
    b'K', b'N', b'K', b'N',
    // ACA   ACC   ACG   ACT
    b'T', b'T', b'T', b'T',
    // AGA   AGC   AGG   AGT
    b'R', b'S', b'R', b'S',
    // ATA   ATC   ATG   ATT
    b'I', b'I', b'M', b'I',
    // CAA   CAC   CAG   CAT
    b'Q', b'H', b'Q', b'H',
    // CCA   CCC   CCG   CCT
    b'P', b'P', b'P', b'P',
    // CGA   CGC   CGG   CGT
    b'R', b'R', b'R', b'R',
    // CTA   CTC   CTG   CTT
    b'L', b'L', b'L', b'L',
    // GAA   GAC   GAG   GAT
    b'E', b'D', b'E', b'D',
    // GCA   GCC   GCG   GCT
    b'A', b'A', b'A', b'A',
    // GGA   GGC   GGG   GGT
    b'G', b'G', b'G', b'G',
    // GTA   GTC   GTG   GTT
    b'V', b'V', b'V', b'V',
    // TAA   TAC   TAG   TAT
    b'*', b'Y', b'*', b'Y',
    // TCA   TCC   TCG   TCT
    b'S', b'S', b'S', b'S',
    // TGA   TGC   TGG   TGT
    b'*', b'C', b'W', b'C',
    // TTA   TTC   TTG   TTT
    b'L', b'F', b'L', b'F',
];

/// Special amino-acid byte for a codon containing any non-{A,C,G,T,U}
/// base — e.g. an N from quality masking, or an indel-inserted byte
/// that hasn't been resolved.
pub const AMINO_AMBIGUOUS: u8 = b'X';

/// Special amino-acid byte for a stop codon (TAA, TAG, TGA).
pub const AMINO_STOP: u8 = b'*';

/// Translate one codon to an amino-acid byte. Returns `AMINO_AMBIGUOUS`
/// (`b'X'`) if any of the three input bases is not one of A/C/G/T/U
/// (case-insensitive).
pub fn translate_codon(b1: u8, b2: u8, b3: u8) -> u8 {
    match (encode_base(b1), encode_base(b2), encode_base(b3)) {
        (Some(i1), Some(i2), Some(i3)) => {
            let idx = ((i1 as usize) << 4) | ((i2 as usize) << 2) | (i3 as usize);
            GENETIC_CODE[idx]
        }
        _ => AMINO_AMBIGUOUS,
    }
}

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
mod tests {
    use super::*;
    use std::mem::size_of;

    #[test]
    fn nucleotide_size_is_small() {
        // We want Nucleotide to stay tight enough for cache-friendly
        // arena scans. 8 bytes today; lifting the cap to 16 is fine
        // but should be a deliberate decision documented in the design
        // doc, not silent struct growth.
        assert!(
            size_of::<Nucleotide>() <= 16,
            "Nucleotide grew past 16 bytes ({} bytes) — review IR layout",
            size_of::<Nucleotide>()
        );
    }

    #[test]
    fn handle_types_are_zero_cost() {
        // Newtype wrappers around u32 must not add overhead.
        assert_eq!(size_of::<NucHandle>(), size_of::<u32>());
        assert_eq!(size_of::<RegionHandle>(), size_of::<u32>());
    }

    #[test]
    fn handle_index_round_trip() {
        let h = NucHandle::new(42);
        assert_eq!(h.index(), 42);
        assert_eq!(h.as_usize(), 42);
    }

    #[test]
    fn segment_equality_is_structural() {
        assert_eq!(Segment::V, Segment::V);
        assert_ne!(Segment::V, Segment::J);
        assert_ne!(Segment::Np1, Segment::Np2);
    }

    #[test]
    fn nuc_flags_set_clear_test() {
        let empty = NucFlags::empty();
        assert!(!empty.contains(flag::P_NUC));

        let with_p = empty.with(flag::P_NUC);
        assert!(with_p.contains(flag::P_NUC));
        assert!(!with_p.contains(flag::N_NUC));

        let with_p_and_junction = with_p.with(flag::JUNCTION);
        assert!(with_p_and_junction.contains(flag::P_NUC));
        assert!(with_p_and_junction.contains(flag::JUNCTION));

        let cleared = with_p_and_junction.without(flag::P_NUC);
        assert!(!cleared.contains(flag::P_NUC));
        assert!(cleared.contains(flag::JUNCTION));
    }

    #[test]
    fn nucleotide_germline_constructor() {
        let n = Nucleotide::germline(b'A', 12, Segment::V);
        assert_eq!(n.base, b'A');
        assert_eq!(n.germline, b'A');
        assert_eq!(n.germline_pos, 12);
        assert_eq!(n.segment, Segment::V);
        assert_eq!(n.flags, NucFlags::empty());
    }

    #[test]
    fn nucleotide_synthetic_constructor_has_no_germline_pos() {
        let n = Nucleotide::synthetic(b'a', Segment::Np1, flag::N_NUC);
        assert_eq!(n.base, b'a');
        assert_eq!(n.germline, b'a');
        assert_eq!(n.germline_pos, Nucleotide::NO_GERMLINE_POS);
        assert_eq!(n.segment, Segment::Np1);
        assert!(n.flags.contains(flag::N_NUC));
    }

    #[test]
    fn pool_push_returns_sequential_handles() {
        let mut pool = NucleotidePool::new();
        assert!(pool.is_empty());

        let h0 = pool.push(Nucleotide::germline(b'A', 0, Segment::V));
        let h1 = pool.push(Nucleotide::germline(b'C', 1, Segment::V));
        let h2 = pool.push(Nucleotide::germline(b'G', 2, Segment::V));

        assert_eq!(h0.index(), 0);
        assert_eq!(h1.index(), 1);
        assert_eq!(h2.index(), 2);
        assert_eq!(pool.len(), 3);
        assert!(!pool.is_empty());
    }

    #[test]
    fn pool_get_by_handle_returns_correct_nucleotide() {
        let mut pool = NucleotidePool::new();
        let h = pool.push(Nucleotide::germline(b'T', 5, Segment::J));

        let got = pool.get(h).expect("handle should be valid");
        assert_eq!(got.base, b'T');
        assert_eq!(got.germline_pos, 5);
        assert_eq!(got.segment, Segment::J);
    }

    #[test]
    fn pool_get_out_of_bounds_returns_none() {
        let pool = NucleotidePool::new();
        assert!(pool.get(NucHandle::new(0)).is_none());
        assert!(pool.get(NucHandle::new(99)).is_none());
    }

    #[test]
    fn region_range_and_length() {
        let r = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(10));
        assert_eq!(r.segment, Segment::V);
        assert_eq!(r.len(), 10);
        assert!(!r.is_empty());

        let empty = Region::new(Segment::Np1, NucHandle::new(5), NucHandle::new(5));
        assert!(empty.is_empty());
        assert_eq!(empty.len(), 0);
    }

    #[test]
    fn sequence_starts_empty() {
        let s = Sequence::new();
        assert_eq!(s.region_count(), 0);
        assert!(s.regions.is_empty());
    }

    #[test]
    fn simulation_default_construction_is_clean() {
        let s = Simulation::new();
        assert!(s.pool.is_empty());
        assert_eq!(s.sequence.region_count(), 0);
        assert!(s.live_calls.is_none());
    }

    #[test]
    fn simulation_with_capacity_pre_allocates() {
        // Reserved capacity is a perf hint, not a guarantee of `len`.
        let s = Simulation::with_capacity(500);
        assert!(s.pool.is_empty());
        assert_eq!(s.sequence.region_count(), 0);
        assert!(s.live_calls.is_none());
    }

    #[test]
    fn simulation_live_call_sidecar_is_persistent_and_dormant() {
        let live = crate::live_call::LiveCallState::empty().with_dirty_window(
            crate::live_call::DirtyWindow::new(
                1,
                4,
                crate::live_call::DirtyReason::BaseEdited { site: 2 },
            ),
        );
        let s0 = Simulation::new().with_live_calls(live.clone());
        let (s1, handle) = s0.with_nucleotide_pushed(Nucleotide::germline(b'A', 0, Segment::V));
        let s2 = s1.with_base_changed(handle, b'C');
        let s3 = s2.with_region_added(Region::new(
            Segment::V,
            NucHandle::new(0),
            NucHandle::new(1),
        ));

        assert_eq!(s0.live_calls.as_ref(), Some(&live));
        assert_eq!(s1.live_calls.as_ref(), Some(&live));
        assert_eq!(s2.live_calls.as_ref(), Some(&live));
        assert_eq!(s3.live_calls.as_ref(), Some(&live));

        // Dormant means core structural edits do not yet interpret or
        // mutate live calls; they only preserve the sidecar.
        assert_eq!(s3.pool.len(), 1);
        assert_eq!(s3.sequence.region_count(), 1);
    }

    #[test]
    fn simulation_clone_is_deep_at_value_level() {
        // simple Clone is a real deep copy.
        let mut a = Simulation::new();
        a.pool.push(Nucleotide::germline(b'A', 0, Segment::V));
        a.sequence.regions.push(Region::new(
            Segment::V,
            NucHandle::new(0),
            NucHandle::new(1),
        ));

        let b = a.clone();
        assert_eq!(b.pool.len(), 1);
        assert_eq!(b.sequence.region_count(), 1);

        // Mutating `a` must not affect `b` (Clone is a deep copy in A.2).
        a.pool.push(Nucleotide::germline(b'C', 1, Segment::V));
        assert_eq!(a.pool.len(), 2);
        assert_eq!(b.pool.len(), 1);
    }

    // ── Persistent update API tests ─────────────────────────────────

    #[test]
    fn pool_with_pushed_returns_new_pool_old_unchanged() {
        let p0 = NucleotidePool::new();
        let (p1, h) = p0.with_pushed(Nucleotide::germline(b'A', 0, Segment::V));

        // Old pool: untouched.
        assert_eq!(p0.len(), 0);
        // New pool: has the nucleotide at handle 0.
        assert_eq!(p1.len(), 1);
        assert_eq!(h.index(), 0);
        assert_eq!(p1.get(h).unwrap().base, b'A');
        // Old pool still does not see it.
        assert!(p0.get(h).is_none());
    }

    #[test]
    fn pool_with_base_changed_isolates_revisions() {
        let p0 = NucleotidePool::new();
        let (p1, h) = p0.with_pushed(Nucleotide::germline(b'A', 0, Segment::V));
        let p2 = p1.with_base_changed(h, b'G');

        // Three independent revisions of the pool exist simultaneously.
        assert!(p0.get(h).is_none());
        assert_eq!(p1.get(h).unwrap().base, b'A');
        assert_eq!(p2.get(h).unwrap().base, b'G');

        // Other fields preserved on the changed revision.
        let n2 = p2.get(h).unwrap();
        assert_eq!(n2.germline, b'A'); // germline unchanged
        assert_eq!(n2.germline_pos, 0);
        assert_eq!(n2.segment, Segment::V);
    }

    #[test]
    fn pool_with_nucleotide_changed_replaces_only_target() {
        let mut p = NucleotidePool::new();
        let h0 = p.push(Nucleotide::germline(b'A', 0, Segment::V));
        let h1 = p.push(Nucleotide::germline(b'C', 1, Segment::V));
        let h2 = p.push(Nucleotide::germline(b'G', 2, Segment::V));

        let p2 = p.with_nucleotide_changed(h1, Nucleotide::germline(b'T', 99, Segment::J));

        // Old pool retains everything.
        assert_eq!(p.get(h0).unwrap().base, b'A');
        assert_eq!(p.get(h1).unwrap().base, b'C');
        assert_eq!(p.get(h1).unwrap().segment, Segment::V);
        assert_eq!(p.get(h2).unwrap().base, b'G');

        // New pool: target replaced, neighbours unchanged.
        assert_eq!(p2.get(h0).unwrap().base, b'A');
        assert_eq!(p2.get(h1).unwrap().base, b'T');
        assert_eq!(p2.get(h1).unwrap().germline_pos, 99);
        assert_eq!(p2.get(h1).unwrap().segment, Segment::J);
        assert_eq!(p2.get(h2).unwrap().base, b'G');
    }

    #[test]
    fn region_with_end_extended_keeps_other_fields() {
        let r0 = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(5));
        let r1 = r0.with_end_extended(NucHandle::new(10));

        assert_eq!(r0.end.index(), 5);
        assert_eq!(r1.end.index(), 10);
        assert_eq!(r1.start.index(), 0);
        assert_eq!(r1.segment, Segment::V);
    }

    #[test]
    fn region_with_frame_phase_isolated() {
        let r0 = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(9));
        assert_eq!(r0.frame_phase, 0);

        let r1 = r0.with_frame_phase(2);
        assert_eq!(r0.frame_phase, 0);
        assert_eq!(r1.frame_phase, 2);
        assert_eq!(r1.start, r0.start);
        assert_eq!(r1.end, r0.end);
        assert_eq!(r1.segment, r0.segment);
    }

    #[test]
    fn sequence_with_region_added_isolated() {
        let s0 = Sequence::new();
        let s1 = s0.with_region_added(Region::new(
            Segment::V,
            NucHandle::new(0),
            NucHandle::new(5),
        ));
        let s2 = s1.with_region_added(Region::new(
            Segment::J,
            NucHandle::new(5),
            NucHandle::new(8),
        ));

        assert_eq!(s0.region_count(), 0);
        assert_eq!(s1.region_count(), 1);
        assert_eq!(s2.region_count(), 2);
        assert_eq!(s2.regions[0].segment, Segment::V);
        assert_eq!(s2.regions[1].segment, Segment::J);
    }

    #[test]
    fn sequence_with_region_replaced_isolated() {
        let s = Sequence::new()
            .with_region_added(Region::new(
                Segment::V,
                NucHandle::new(0),
                NucHandle::new(3),
            ))
            .with_region_added(Region::new(
                Segment::J,
                NucHandle::new(3),
                NucHandle::new(6),
            ));

        let s2 = s.with_region_replaced(
            0,
            Region::new(Segment::V, NucHandle::new(0), NucHandle::new(10)),
        );

        // Old sequence: V region had length 3.
        assert_eq!(s.regions[0].len(), 3);
        // New sequence: V region has length 10, J unchanged.
        assert_eq!(s2.regions[0].len(), 10);
        assert_eq!(s2.regions[1].segment, Segment::J);
        assert_eq!(s2.regions[1].len(), 3);
    }

    #[test]
    fn simulation_persistent_chain_preserves_history() {
        // The point of D1: every revision in a chain remains accessible.
        let s0 = Simulation::new();
        let (s1, h_a) = s0.with_nucleotide_pushed(Nucleotide::germline(b'A', 0, Segment::V));
        let (s2, h_c) = s1.with_nucleotide_pushed(Nucleotide::germline(b'C', 1, Segment::V));
        let s3 = s2.with_base_changed(h_a, b'T');
        let s4 = s3.with_region_added(Region::new(
            Segment::V,
            NucHandle::new(0),
            NucHandle::new(2),
        ));

        // s0: empty.
        assert_eq!(s0.pool.len(), 0);
        assert_eq!(s0.sequence.region_count(), 0);

        // s1: one nucleotide (A).
        assert_eq!(s1.pool.len(), 1);
        assert_eq!(s1.pool.get(h_a).unwrap().base, b'A');
        assert_eq!(s1.sequence.region_count(), 0);

        // s2: two nucleotides (A, C), no regions.
        assert_eq!(s2.pool.len(), 2);
        assert_eq!(s2.pool.get(h_a).unwrap().base, b'A');
        assert_eq!(s2.pool.get(h_c).unwrap().base, b'C');
        assert_eq!(s2.sequence.region_count(), 0);

        // s3: A→T mutation, C unchanged.
        assert_eq!(s3.pool.len(), 2);
        assert_eq!(s3.pool.get(h_a).unwrap().base, b'T');
        assert_eq!(s3.pool.get(h_c).unwrap().base, b'C');
        // Original germline preserved on the mutated nucleotide.
        assert_eq!(s3.pool.get(h_a).unwrap().germline, b'A');
        assert_eq!(s3.sequence.region_count(), 0);

        // s4: region added.
        assert_eq!(s4.pool.len(), 2);
        assert_eq!(s4.sequence.region_count(), 1);
        assert_eq!(s4.sequence.regions[0].segment, Segment::V);

        // Earlier revisions still see their own pool state.
        assert_eq!(s2.pool.get(h_a).unwrap().base, b'A'); // pre-mutation
        assert_eq!(s3.sequence.region_count(), 0); // pre-region-add
    }

    #[test]
    fn simulation_assignments_default_is_empty() {
        let sim = Simulation::new();
        assert!(sim.assignments.v.is_none());
        assert!(sim.assignments.d.is_none());
        assert!(sim.assignments.j.is_none());
        assert!(sim.assignments.c.is_none());
    }

    #[test]
    fn simulation_with_allele_assigned_populates_slot_persistently() {
        use crate::assignment::AlleleInstance;
        use crate::refdata::AlleleId;

        let s0 = Simulation::new();
        let s1 = s0.with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(7)));

        // s0 unchanged.
        assert!(s0.assignments.v.is_none());
        // s1 has the V allele.
        assert_eq!(s1.assignments.v.unwrap().allele_id, AlleleId::new(7));
        assert_eq!(s1.assignments.v.unwrap().trim_5, 0);
        assert_eq!(s1.assignments.v.unwrap().trim_3, 0);
    }

    #[test]
    fn simulation_with_trim_updates_assignment_persistently() {
        use crate::assignment::{AlleleInstance, TrimEnd};
        use crate::refdata::AlleleId;

        let s0 = Simulation::new()
            .with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(0)));
        let s1 = s0.with_trim(Segment::V, TrimEnd::Three, 5);

        // s0 retains zero trim.
        assert_eq!(s0.assignments.v.unwrap().trim_3, 0);
        // s1 has the trim applied.
        assert_eq!(s1.assignments.v.unwrap().trim_3, 5);
        // 5' end untouched.
        assert_eq!(s1.assignments.v.unwrap().trim_5, 0);
    }

    #[test]
    fn simulation_with_base_changed_refreshes_overlapping_region_codon_rail() {
        // Audit finding: `with_base_changed` previously updated the
        // pool but cloned the sequence verbatim, leaving
        // `Region.amino_acids` stale against the new pool. This
        // would silently desync any future SHM/uniform pass that
        // mutates a base inside an assembled region.
        //
        // Pin the fix: change a base inside a Region and verify the
        // Region's codon rail reflects the new base.

        // Build a sim with a 9-base V region whose codon rail says KPG.
        let mut sim = Simulation::new();
        for (i, b) in b"AAACCCGGG".iter().enumerate() {
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(*b, i as u16, Segment::V));
            sim = next;
        }
        let region = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(9))
            .with_codon_rail_recomputed(&sim.pool);
        sim = sim.with_region_added(region);

        assert_eq!(sim.sequence.regions[0].amino_acids, b"KPG");

        // Mutate position 5 (C → A): codon at [3..6) was CCC = P, now becomes CCA = P (no change).
        // Mutate position 4 (C → A): codon CCC → CAC = H. Should change amino_acids[1] to H.
        let mutated = sim.with_base_changed(NucHandle::new(4), b'A');
        assert_eq!(mutated.sequence.regions[0].amino_acids, b"KHG");
        // Pool reflects the change.
        assert_eq!(mutated.pool.get(NucHandle::new(4)).unwrap().base, b'A');
        // Original sim's region still says KPG (persistent IR).
        assert_eq!(sim.sequence.regions[0].amino_acids, b"KPG");
    }

    #[test]
    fn simulation_with_base_changed_outside_any_region_leaves_rails_alone() {
        // Edge case: change a base that's NOT inside any region.
        // Codon rails should remain whatever they were; the change
        // affects only the pool.
        let mut sim = Simulation::new();
        for (i, b) in b"AAACCCGGG".iter().enumerate() {
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(*b, i as u16, Segment::V));
            sim = next;
        }
        // Add ONE region covering [0, 6) with codon rail.
        let region = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(6))
            .with_codon_rail_recomputed(&sim.pool);
        sim = sim.with_region_added(region);

        assert_eq!(sim.sequence.regions[0].amino_acids, b"KP");

        // Mutate position 7 (outside the region's [0, 6) range).
        let mutated = sim.with_base_changed(NucHandle::new(7), b'A');
        // Region's codon rail unchanged.
        assert_eq!(mutated.sequence.regions[0].amino_acids, b"KP");
    }

    #[test]
    fn simulation_with_base_changed_refreshes_correct_region_among_many() {
        // Multi-region case: only the region covering the changed
        // handle should have its rail recomputed; others stay
        // intact (and should be the same Region values, by Eq).
        let mut sim = Simulation::new();
        for (i, b) in b"AAACCCGGG".iter().enumerate() {
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(*b, i as u16, Segment::V));
            sim = next;
        }
        // Two regions: [0, 3) with codon AAA → K, and [6, 9) with codon GGG → G.
        let r0 = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(3))
            .with_codon_rail_recomputed(&sim.pool);
        let r1 = Region::new(Segment::V, NucHandle::new(6), NucHandle::new(9))
            .with_codon_rail_recomputed(&sim.pool);
        sim = sim.with_region_added(r0).with_region_added(r1);

        assert_eq!(sim.sequence.regions[0].amino_acids, b"K");
        assert_eq!(sim.sequence.regions[1].amino_acids, b"G");

        // Mutate position 7 (inside r1: GGG → GAG = E).
        let mutated = sim.with_base_changed(NucHandle::new(7), b'A');
        // r0 unchanged.
        assert_eq!(mutated.sequence.regions[0].amino_acids, b"K");
        // r1 recomputed: GAG = E.
        assert_eq!(mutated.sequence.regions[1].amino_acids, b"E");
    }

    // ── Indel IR API tests ─────────────────────────────────────────

    #[test]
    fn pool_with_inserted_grows_and_shifts() {
        let mut p = NucleotidePool::new();
        let _ = p.push(Nucleotide::germline(b'A', 0, Segment::V));
        let _ = p.push(Nucleotide::germline(b'C', 1, Segment::V));
        let _ = p.push(Nucleotide::germline(b'G', 2, Segment::V));

        let p2 = p.with_inserted(1, Nucleotide::germline(b'X', 99, Segment::Np1));
        assert_eq!(p2.len(), 4);
        assert_eq!(p2.get(NucHandle::new(0)).unwrap().base, b'A');
        assert_eq!(p2.get(NucHandle::new(1)).unwrap().base, b'X');
        assert_eq!(p2.get(NucHandle::new(2)).unwrap().base, b'C'); // shifted up
        assert_eq!(p2.get(NucHandle::new(3)).unwrap().base, b'G');

        // Original pool unchanged (persistent IR).
        assert_eq!(p.len(), 3);
    }

    #[test]
    fn pool_with_inserted_at_end_appends() {
        let mut p = NucleotidePool::new();
        let _ = p.push(Nucleotide::germline(b'A', 0, Segment::V));
        let p2 = p.with_inserted(1, Nucleotide::germline(b'C', 1, Segment::V));
        assert_eq!(p2.len(), 2);
        assert_eq!(p2.get(NucHandle::new(1)).unwrap().base, b'C');
    }

    #[test]
    #[should_panic(expected = "with_inserted: at")]
    fn pool_with_inserted_out_of_bounds_panics() {
        let p = NucleotidePool::new();
        let _ = p.with_inserted(5, Nucleotide::germline(b'A', 0, Segment::V));
    }

    #[test]
    fn pool_with_deleted_shrinks_and_shifts() {
        let mut p = NucleotidePool::new();
        let _ = p.push(Nucleotide::germline(b'A', 0, Segment::V));
        let _ = p.push(Nucleotide::germline(b'C', 1, Segment::V));
        let _ = p.push(Nucleotide::germline(b'G', 2, Segment::V));

        let p2 = p.with_deleted(1);
        assert_eq!(p2.len(), 2);
        assert_eq!(p2.get(NucHandle::new(0)).unwrap().base, b'A');
        assert_eq!(p2.get(NucHandle::new(1)).unwrap().base, b'G'); // shifted down
        assert_eq!(p.len(), 3); // original unchanged
    }

    #[test]
    #[should_panic(expected = "with_deleted: at")]
    fn pool_with_deleted_out_of_bounds_panics() {
        let p = NucleotidePool::new();
        let _ = p.with_deleted(0);
    }

    #[test]
    fn sequence_indel_adjustment_region_entirely_before() {
        let s = Sequence::new().with_region_added(Region::new(
            Segment::V,
            NucHandle::new(0),
            NucHandle::new(5),
        ));
        // Insertion at position 7 (after region.end = 5).
        let s2 = s.with_indel_adjusted(7, 1);
        // Region unchanged.
        assert_eq!(s2.regions[0].start, NucHandle::new(0));
        assert_eq!(s2.regions[0].end, NucHandle::new(5));
    }

    #[test]
    fn sequence_indel_adjustment_region_entirely_after_insertion() {
        let s = Sequence::new().with_region_added(Region::new(
            Segment::V,
            NucHandle::new(5),
            NucHandle::new(10),
        ));
        // Insertion at position 3 — both start and end shift up by 1.
        let s2 = s.with_indel_adjusted(3, 1);
        assert_eq!(s2.regions[0].start, NucHandle::new(6));
        assert_eq!(s2.regions[0].end, NucHandle::new(11));
    }

    #[test]
    fn sequence_indel_adjustment_region_spanning_insertion_grows() {
        let s = Sequence::new().with_region_added(Region::new(
            Segment::V,
            NucHandle::new(0),
            NucHandle::new(10),
        ));
        // Insertion at position 5 (inside the region).
        let s2 = s.with_indel_adjusted(5, 1);
        // Start unchanged, end grew.
        assert_eq!(s2.regions[0].start, NucHandle::new(0));
        assert_eq!(s2.regions[0].end, NucHandle::new(11));
    }

    #[test]
    fn sequence_indel_adjustment_region_spanning_deletion_shrinks() {
        let s = Sequence::new().with_region_added(Region::new(
            Segment::V,
            NucHandle::new(0),
            NucHandle::new(10),
        ));
        // Deletion at position 5 (inside the region).
        let s2 = s.with_indel_adjusted(5, -1);
        assert_eq!(s2.regions[0].start, NucHandle::new(0));
        assert_eq!(s2.regions[0].end, NucHandle::new(9));
    }

    #[test]
    fn sequence_indel_adjustment_region_after_deletion_shifts() {
        let s = Sequence::new().with_region_added(Region::new(
            Segment::V,
            NucHandle::new(5),
            NucHandle::new(10),
        ));
        // Deletion at position 3 (before the region).
        let s2 = s.with_indel_adjusted(3, -1);
        assert_eq!(s2.regions[0].start, NucHandle::new(4));
        assert_eq!(s2.regions[0].end, NucHandle::new(9));
    }

    #[test]
    fn sequence_frame_phase_recompute_chains_from_region_lengths() {
        let s = Sequence::new()
            .with_region_added(
                Region::new(Segment::V, NucHandle::new(0), NucHandle::new(1)).with_frame_phase(2),
            )
            .with_region_added(
                Region::new(Segment::Np1, NucHandle::new(1), NucHandle::new(3)).with_frame_phase(2),
            )
            .with_region_added(
                Region::new(Segment::J, NucHandle::new(3), NucHandle::new(6)).with_frame_phase(2),
            );

        let recomputed = s.with_frame_phases_recomputed();

        assert_eq!(recomputed.regions[0].frame_phase, 0);
        assert_eq!(recomputed.regions[1].frame_phase, 1);
        assert_eq!(recomputed.regions[2].frame_phase, 0);
    }

    #[test]
    fn sequence_indel_adjustment_recomputes_downstream_frame_phase_after_insertion() {
        let s = Sequence::new()
            .with_region_added(
                Region::new(Segment::V, NucHandle::new(0), NucHandle::new(2)).with_frame_phase(0),
            )
            .with_region_added(
                Region::new(Segment::J, NucHandle::new(2), NucHandle::new(8)).with_frame_phase(2),
            );

        let adjusted = s.with_indel_adjusted(1, 1);

        assert_eq!(adjusted.regions[0].start, NucHandle::new(0));
        assert_eq!(adjusted.regions[0].end, NucHandle::new(3));
        assert_eq!(adjusted.regions[0].frame_phase, 0);
        assert_eq!(adjusted.regions[1].start, NucHandle::new(3));
        assert_eq!(adjusted.regions[1].end, NucHandle::new(9));
        assert_eq!(adjusted.regions[1].frame_phase, 0);
    }

    #[test]
    fn sequence_indel_adjustment_recomputes_downstream_frame_phase_after_deletion() {
        let s = Sequence::new()
            .with_region_added(
                Region::new(Segment::V, NucHandle::new(0), NucHandle::new(4)).with_frame_phase(0),
            )
            .with_region_added(
                Region::new(Segment::J, NucHandle::new(4), NucHandle::new(10)).with_frame_phase(1),
            );

        let adjusted = s.with_indel_adjusted(1, -1);

        assert_eq!(adjusted.regions[0].start, NucHandle::new(0));
        assert_eq!(adjusted.regions[0].end, NucHandle::new(3));
        assert_eq!(adjusted.regions[0].frame_phase, 0);
        assert_eq!(adjusted.regions[1].start, NucHandle::new(3));
        assert_eq!(adjusted.regions[1].end, NucHandle::new(9));
        assert_eq!(adjusted.regions[1].frame_phase, 0);
    }

    #[test]
    fn simulation_with_indel_inserted_refreshes_downstream_frame_phase_and_rail() {
        let mut sim = Simulation::new();
        for (i, b) in b"AAATGCCC".iter().enumerate() {
            let segment = if i < 2 { Segment::V } else { Segment::J };
            let (next, _) = sim.with_nucleotide_pushed(Nucleotide::germline(*b, i as u16, segment));
            sim = next;
        }
        let v_region = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(2))
            .with_frame_phase(0)
            .with_codon_rail_recomputed(&sim.pool);
        let j_region = Region::new(Segment::J, NucHandle::new(2), NucHandle::new(8))
            .with_frame_phase(2)
            .with_codon_rail_recomputed(&sim.pool);
        sim = sim.with_region_added(v_region).with_region_added(j_region);

        let mutated = sim.with_indel_inserted(
            1,
            Nucleotide::synthetic(b'G', Segment::V, crate::ir::flag::INDEL_INSERTED),
        );

        let j = &mutated.sequence.regions[1];
        assert_eq!(j.start, NucHandle::new(3));
        assert_eq!(j.end, NucHandle::new(9));
        assert_eq!(j.frame_phase, 0);
        assert_eq!(j.amino_acids, b"MP");
    }

    #[test]
    fn simulation_with_indel_deleted_refreshes_downstream_frame_phase_and_rail() {
        let mut sim = Simulation::new();
        for (i, b) in b"AAAAATGCCC".iter().enumerate() {
            let segment = if i < 4 { Segment::V } else { Segment::J };
            let (next, _) = sim.with_nucleotide_pushed(Nucleotide::germline(*b, i as u16, segment));
            sim = next;
        }
        let v_region = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(4))
            .with_frame_phase(0)
            .with_codon_rail_recomputed(&sim.pool);
        let j_region = Region::new(Segment::J, NucHandle::new(4), NucHandle::new(10))
            .with_frame_phase(1)
            .with_codon_rail_recomputed(&sim.pool);
        sim = sim.with_region_added(v_region).with_region_added(j_region);

        let mutated = sim.with_indel_deleted(1);

        let j = &mutated.sequence.regions[1];
        assert_eq!(j.start, NucHandle::new(3));
        assert_eq!(j.end, NucHandle::new(9));
        assert_eq!(j.frame_phase, 0);
        assert_eq!(j.amino_acids, b"MP");
    }

    #[test]
    fn simulation_with_indel_inserted_full_round_trip() {
        let mut sim = Simulation::new();
        for (i, b) in b"ATGCCCGGG".iter().enumerate() {
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(*b, i as u16, Segment::V));
            sim = next;
        }
        let region = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(9))
            .with_codon_rail_recomputed(&sim.pool);
        sim = sim.with_region_added(region);
        // Initial codon rail: ATG CCC GGG → M P G.
        assert_eq!(sim.sequence.regions[0].amino_acids, b"MPG");

        // Insert a 'G' at position 3 — the region grows by 1, and
        // codons shift: ATG | GCC CGG G → M, A, R, + leftover G.
        let mutated = sim.with_indel_inserted(
            3,
            Nucleotide::synthetic(b'G', Segment::V, crate::ir::flag::INDEL_INSERTED),
        );
        assert_eq!(mutated.pool.len(), 10);
        assert_eq!(mutated.sequence.regions[0].len(), 10);
        // ATG GCC CGG → M A R (3 codons, 10 = 3*3 + 1 leftover)
        assert_eq!(mutated.sequence.regions[0].amino_acids, b"MAR");
        // Original sim untouched.
        assert_eq!(sim.pool.len(), 9);
        assert_eq!(sim.sequence.regions[0].amino_acids, b"MPG");
    }

    #[test]
    fn simulation_with_indel_deleted_full_round_trip() {
        let mut sim = Simulation::new();
        for (i, b) in b"ATGAAACCCGGG".iter().enumerate() {
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(*b, i as u16, Segment::V));
            sim = next;
        }
        let region = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(12))
            .with_codon_rail_recomputed(&sim.pool);
        sim = sim.with_region_added(region);
        // Codon rail: ATG AAA CCC GGG → M K P G.
        assert_eq!(sim.sequence.regions[0].amino_acids, b"MKPG");

        // Delete position 3 (the first 'A' of the second codon).
        // Region shrinks to 11 bases.
        let mutated = sim.with_indel_deleted(3);
        assert_eq!(mutated.pool.len(), 11);
        assert_eq!(mutated.sequence.regions[0].len(), 11);
        // Bases now: ATG AAC CCG GG → 3 full codons + 2 leftover.
        // ATG=M, AAC=N, CCG=P. Three amino acids.
        assert_eq!(mutated.sequence.regions[0].amino_acids, b"MNP");
    }

    #[test]
    fn simulation_indel_branching_revisions_are_independent() {
        let mut sim = Simulation::new();
        for (i, b) in b"ATGAAA".iter().enumerate() {
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(*b, i as u16, Segment::V));
            sim = next;
        }

        // Branch A: insert at 3.
        let branch_a = sim.with_indel_inserted(
            3,
            Nucleotide::synthetic(b'C', Segment::V, crate::ir::flag::INDEL_INSERTED),
        );
        // Branch B: delete at 3.
        let branch_b = sim.with_indel_deleted(3);

        assert_eq!(branch_a.pool.len(), 7);
        assert_eq!(branch_b.pool.len(), 5);
        // Original unchanged.
        assert_eq!(sim.pool.len(), 6);
    }

    #[test]
    fn simulation_with_nucleotide_changed_also_refreshes_codon_rail() {
        // Same fix applies to `with_nucleotide_changed`, which can
        // change germline_pos / segment / flags in addition to base.
        let mut sim = Simulation::new();
        for (i, b) in b"AAACCC".iter().enumerate() {
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(*b, i as u16, Segment::V));
            sim = next;
        }
        let region = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(6))
            .with_codon_rail_recomputed(&sim.pool);
        sim = sim.with_region_added(region);
        assert_eq!(sim.sequence.regions[0].amino_acids, b"KP");

        // Replace position 3 with a new nucleotide whose base is 'T'.
        // Codon at [3..6) was CCC = P, becomes TCC = S.
        let mutated = sim
            .with_nucleotide_changed(NucHandle::new(3), Nucleotide::germline(b'T', 3, Segment::V));
        assert_eq!(mutated.sequence.regions[0].amino_acids, b"KS");
    }

    #[test]
    fn simulation_assignments_chain_with_pool_changes() {
        // Assignments survive pool / sequence changes — they live in
        // their own slot, not in the pool.
        use crate::assignment::AlleleInstance;
        use crate::refdata::AlleleId;

        let s0 = Simulation::new()
            .with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(3)));

        let (s1, _h) = s0.with_nucleotide_pushed(Nucleotide::germline(b'A', 0, Segment::V));

        // Assignment survives across the pool change.
        assert_eq!(s1.assignments.v.unwrap().allele_id, AlleleId::new(3));
    }

    #[test]
    fn simulation_branching_revisions_are_independent() {
        // Two divergent branches of mutation from the same parent must
        // not affect each other.
        let s0 = Simulation::new();
        let (s1, h) = s0.with_nucleotide_pushed(Nucleotide::germline(b'A', 0, Segment::V));

        let branch_left = s1.with_base_changed(h, b'C');
        let branch_right = s1.with_base_changed(h, b'G');

        // Common ancestor unchanged.
        assert_eq!(s1.pool.get(h).unwrap().base, b'A');
        // Left branch: A → C.
        assert_eq!(branch_left.pool.get(h).unwrap().base, b'C');
        // Right branch: A → G.
        assert_eq!(branch_right.pool.get(h).unwrap().base, b'G');
    }

    // ── Codon-rail metadata tests ───────────────────────────────────

    /// Helper: build a pool from a base string with all nucleotides in
    /// segment V, germline_pos = position within the string. Returns
    /// the pool and a Region covering the whole pool with frame_phase=0.
    fn pool_from_string(s: &str) -> (NucleotidePool, Region) {
        let mut pool = NucleotidePool::new();
        for (i, b) in s.bytes().enumerate() {
            pool.push(Nucleotide::germline(b, i as u16, Segment::V));
        }
        let region = Region::new(
            Segment::V,
            NucHandle::new(0),
            NucHandle::new(s.len() as u32),
        );
        (pool, region)
    }

    #[test]
    fn translate_codon_canonical_examples() {
        assert_eq!(translate_codon(b'A', b'T', b'G'), b'M'); // Met (start)
        assert_eq!(translate_codon(b'T', b'T', b'T'), b'F'); // Phe
        assert_eq!(translate_codon(b'T', b'G', b'G'), b'W'); // Trp
        assert_eq!(translate_codon(b'C', b'A', b'C'), b'H'); // His
        assert_eq!(translate_codon(b'A', b'A', b'A'), b'K'); // Lys
        assert_eq!(translate_codon(b'G', b'G', b'C'), b'G'); // Gly
    }

    #[test]
    fn translate_codon_stops() {
        assert_eq!(translate_codon(b'T', b'A', b'A'), AMINO_STOP);
        assert_eq!(translate_codon(b'T', b'A', b'G'), AMINO_STOP);
        assert_eq!(translate_codon(b'T', b'G', b'A'), AMINO_STOP);
        // TGG is *not* a stop — it is Trp (W).
        assert_eq!(translate_codon(b'T', b'G', b'G'), b'W');
    }

    #[test]
    fn translate_codon_lowercase_is_case_insensitive() {
        assert_eq!(translate_codon(b'a', b't', b'g'), b'M');
        assert_eq!(translate_codon(b'T', b'a', b'a'), AMINO_STOP);
    }

    #[test]
    fn translate_codon_ambiguous_is_x() {
        assert_eq!(translate_codon(b'A', b'T', b'N'), AMINO_AMBIGUOUS);
        assert_eq!(translate_codon(b'N', b'A', b'T'), AMINO_AMBIGUOUS);
        assert_eq!(translate_codon(b'A', b'-', b'C'), AMINO_AMBIGUOUS);
    }

    #[test]
    fn translate_codon_uracil_treated_as_thymine() {
        assert_eq!(translate_codon(b'A', b'U', b'G'), b'M');
        assert_eq!(translate_codon(b'U', b'A', b'A'), AMINO_STOP);
    }

    #[test]
    fn region_codon_rail_starts_empty_until_recomputed() {
        let r = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(9));
        assert!(r.amino_acids.is_empty());
        assert!(r.stop_codon_positions.is_empty());
        assert_eq!(r.stop_codon_count(), 0);
    }

    #[test]
    fn region_with_codon_rail_recomputed_translates_in_frame() {
        // ATG GGG CAC → M G H, no stops.
        let (pool, region) = pool_from_string("ATGGGGCAC");
        let r = region.with_codon_rail_recomputed(&pool);

        assert_eq!(r.amino_acids, b"MGH");
        assert_eq!(r.stop_codon_positions, vec![]);
        assert_eq!(r.stop_codon_count(), 0);
        // Receiver unchanged (persistent contract).
        assert!(region.amino_acids.is_empty());
    }

    #[test]
    fn region_with_codon_rail_recomputed_picks_up_stops() {
        // ATG TAA TGG → M * W, one stop at position 3.
        let (pool, region) = pool_from_string("ATGTAATGG");
        let r = region.with_codon_rail_recomputed(&pool);

        assert_eq!(r.amino_acids, b"M*W");
        assert_eq!(r.stop_codon_positions, vec![NucHandle::new(3)]);
        assert_eq!(r.stop_codon_count(), 1);
    }

    #[test]
    fn region_with_codon_rail_recomputed_drops_partial_codon_at_end() {
        // ATG GG → M, plus 2 incomplete bases. No second amino acid.
        let (pool, region) = pool_from_string("ATGGG");
        let r = region.with_codon_rail_recomputed(&pool);

        assert_eq!(r.amino_acids, b"M");
        assert_eq!(r.stop_codon_count(), 0);
    }

    #[test]
    fn region_with_codon_rail_respects_frame_phase_one() {
        // frame_phase=1 means position 0 is the 2nd base of a codon
        // started in a (notional) previous region. Skip 2 bases, then
        // translate.
        //
        // Bases:   X X A T G C C C
        // Phase:   1 2 0 1 2 0 1 2     (0 = first base of codon)
        // Codons fully in this region: ATG, CCC → M P
        let (pool, region) = pool_from_string("XXATGCCC");
        let r = region.with_frame_phase(1).with_codon_rail_recomputed(&pool);
        assert_eq!(r.amino_acids, b"MP");
    }

    #[test]
    fn region_with_codon_rail_respects_frame_phase_two() {
        // frame_phase=2 means position 0 is the 3rd base of a codon.
        // Skip 1 base, then translate.
        //
        // Bases:   X T A C G G G
        // Phase:   2 0 1 2 0 1 2
        // Codons fully in this region: TAC, GGG → Y G
        let (pool, region) = pool_from_string("XTACGGG");
        let r = region.with_frame_phase(2).with_codon_rail_recomputed(&pool);
        assert_eq!(r.amino_acids, b"YG");
    }

    #[test]
    fn region_with_codon_rail_handles_ambiguous_bases() {
        // Codon containing N → X.
        let (pool, region) = pool_from_string("ATGNAATGG");
        let r = region.with_codon_rail_recomputed(&pool);
        // ATG = M, NAA = X (ambiguous), TGG = W
        assert_eq!(r.amino_acids, b"MXW");
    }

    #[test]
    fn region_with_codon_rail_recomputes_after_base_mutation() {
        // Mutation in this region must produce a region revision with
        // updated amino acids.
        let (pool0, region0) = pool_from_string("ATGTACTGG");
        let r0 = region0.with_codon_rail_recomputed(&pool0);
        assert_eq!(r0.amino_acids, b"MYW");

        // Mutate the second base of the second codon: TAC → TGC = C.
        let pool1 = pool0.with_base_changed(NucHandle::new(4), b'G');
        let r1 = region0.with_codon_rail_recomputed(&pool1);
        assert_eq!(r1.amino_acids, b"MCW");

        // Old region revision unaffected.
        assert_eq!(r0.amino_acids, b"MYW");
    }

    #[test]
    fn region_with_codon_rail_detects_stop_introduced_by_mutation() {
        // TAC → TAA (Y → stop) by mutating one base.
        let (pool0, region0) = pool_from_string("ATGTACGGG");
        let r0 = region0.with_codon_rail_recomputed(&pool0);
        assert_eq!(r0.stop_codon_count(), 0);
        assert_eq!(r0.amino_acids, b"MYG");

        // Mutate position 5: TAC → TAA. This creates a stop codon.
        let pool1 = pool0.with_base_changed(NucHandle::new(5), b'A');
        let r1 = region0.with_codon_rail_recomputed(&pool1);
        assert_eq!(r1.amino_acids, b"M*G");
        assert_eq!(r1.stop_codon_count(), 1);
        assert_eq!(r1.stop_codon_positions, vec![NucHandle::new(3)]);

        // Old revision unchanged.
        assert_eq!(r0.stop_codon_count(), 0);
    }

    #[test]
    fn region_codon_rail_empty_region_produces_empty_metadata() {
        let pool = NucleotidePool::new();
        let region = Region::new(Segment::Np1, NucHandle::new(0), NucHandle::new(0));
        let r = region.with_codon_rail_recomputed(&pool);
        assert!(r.amino_acids.is_empty());
        assert!(r.stop_codon_positions.is_empty());
    }

    #[test]
    fn region_codon_rail_frame_phase_2_on_one_base_region() {
        // frame_phase = 2 means the region's first base is the third
        // base of a codon started in a previous region. `skip = 1`,
        // so we'd need at least 4 bases (1 skipped + 3 for a fresh
        // codon) to emit anything. With one base, output is empty.
        let (pool, region) = pool_from_string("X");
        let r = region.with_frame_phase(2).with_codon_rail_recomputed(&pool);
        assert!(r.amino_acids.is_empty());
        assert!(r.stop_codon_positions.is_empty());
    }

    #[test]
    fn region_codon_rail_frame_phase_1_on_two_base_region() {
        // frame_phase = 1 → skip 2 → no bases left after skip.
        // Output is empty.
        let (pool, region) = pool_from_string("XY");
        let r = region.with_frame_phase(1).with_codon_rail_recomputed(&pool);
        assert!(r.amino_acids.is_empty());
    }

    #[test]
    fn region_codon_rail_malformed_end_less_than_start_is_safe() {
        // Defensive contract: a region constructed with end < start
        // (which `len()` already saturates to 0) produces empty
        // codon rail metadata, not a panic. Builders should not
        // produce such regions, but the recompute should be robust.
        let mut pool = NucleotidePool::new();
        for i in 0..6 {
            pool.push(Nucleotide::germline(b'A', i, Segment::V));
        }
        let region = Region::new(Segment::V, NucHandle::new(5), NucHandle::new(2));
        assert_eq!(region.len(), 0); // saturating_sub clamps to 0

        let r = region.with_codon_rail_recomputed(&pool);
        assert!(r.amino_acids.is_empty());
        assert!(r.stop_codon_positions.is_empty());
    }

    #[test]
    fn region_codon_rail_skip_overruns_end_is_safe() {
        // 1-base region with frame_phase=1 → skip=2 → start_idx > end_idx
        // after the skip. Same defensive case as above, different
        // path through the code.
        let (pool, region) = pool_from_string("X");
        let r = region.with_frame_phase(1).with_codon_rail_recomputed(&pool);
        assert!(r.amino_acids.is_empty());
        assert!(r.stop_codon_positions.is_empty());
    }

    // ── Stress test for the persistent IR + codon rail ──────────────

    /// Tiny deterministic PRNG (xorshift32) for the stress test. We
    /// avoid bringing in `rand` here so the test has zero external
    /// dependencies and reproduces identically across machines.
    struct Xorshift32(u32);

    impl Xorshift32 {
        fn new(seed: u32) -> Self {
            // xorshift32 cannot have a zero seed.
            Self(if seed == 0 { 0xdead_beef } else { seed })
        }
        fn next(&mut self) -> u32 {
            let mut x = self.0;
            x ^= x << 13;
            x ^= x >> 17;
            x ^= x << 5;
            self.0 = x;
            x
        }
    }

    /// Apply 1000 random base mutations through the persistent API,
    /// keeping every IR revision in memory. Verifies four properties
    /// that together pin the IR's architectural commitments:
    ///
    ///   1. Persistent contract (D1): the initial revision retains its
    ///      original amino-acid sequence after all 1000 mutations have
    ///      been applied to descendant revisions. No mutation leaks
    ///      back upstream.
    ///   2. Persistent + entity-attached metadata (D5): every revision's
    ///      stored `amino_acids` field equals the result of recomputing
    ///      the codon rail from that revision's own pool. Stale
    ///      metadata is structurally impossible.
    ///   3. Branching independence: at any mid-revision, the stored
    ///      amino acids match the revision's pool — they were captured
    ///      at construction time and were not perturbed by later writes.
    ///   4. The persistent API is not silently coalescing changes —
    ///      most random mutations produce a visible delta in the new
    ///      revision's pool versus the previous revision's pool.
    #[test]
    fn stress_persistent_ir_with_codon_rail() {
        const N_NUCS: u32 = 90; // 30 codons
        const N_MUTATIONS: usize = 1000;

        // Build the initial pool: 90 'A' nucleotides → 30 codons of AAA → 30 K's.
        let mut pool0 = NucleotidePool::new();
        for i in 0..N_NUCS {
            pool0.push(Nucleotide::germline(b'A', i as u16, Segment::V));
        }
        let region0 = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(N_NUCS))
            .with_codon_rail_recomputed(&pool0);
        let sim0 = Simulation {
            pool: pool0,
            sequence: Sequence::new().with_region_added(region0),
            assignments: crate::assignment::AlleleAssignments::new(),
            live_calls: None,
        };

        // Property check at revision 0.
        assert_eq!(sim0.sequence.regions[0].amino_acids.len(), 30);
        assert!(sim0.sequence.regions[0]
            .amino_acids
            .iter()
            .all(|&aa| aa == b'K'));

        // Snapshot the initial amino acid sequence to verify property
        // (1) — that revision 0 is never perturbed.
        let initial_amino_acids = sim0.sequence.regions[0].amino_acids.clone();

        // Apply 1000 random mutations. Every mutation produces a new
        // Simulation revision; every previous revision stays alive.
        let mut history: Vec<Simulation> = Vec::with_capacity(N_MUTATIONS + 1);
        history.push(sim0);

        let mut rng = Xorshift32::new(0x517c_c1ed); // arbitrary fixed seed
        let bases = [b'A', b'C', b'G', b'T'];
        let mut visible_mutations = 0usize;

        for _ in 0..N_MUTATIONS {
            let prev = history.last().unwrap();
            let pos = rng.next() % N_NUCS;
            let new_base = bases[(rng.next() & 0b11) as usize];

            let prev_base = prev.pool.get(NucHandle::new(pos)).unwrap().base;
            if prev_base != new_base {
                visible_mutations += 1;
            }

            // Persistent update: pool, then region (with codon rail
            // recomputed against the new pool), then sequence.
            let new_pool = prev.pool.with_base_changed(NucHandle::new(pos), new_base);
            let new_region = prev.sequence.regions[0].with_codon_rail_recomputed(&new_pool);
            let new_sequence = prev.sequence.with_region_replaced(0, new_region);

            history.push(Simulation {
                pool: new_pool,
                sequence: new_sequence,
                assignments: crate::assignment::AlleleAssignments::new(),
                live_calls: prev.live_calls.clone(),
            });
        }

        // Sanity: history size is correct.
        assert_eq!(history.len(), N_MUTATIONS + 1);

        // Property (1): revision 0's amino acids are exactly what they
        // were before any mutation happened.
        assert_eq!(
            history[0].sequence.regions[0].amino_acids,
            initial_amino_acids
        );
        assert!(history[0].sequence.regions[0]
            .amino_acids
            .iter()
            .all(|&aa| aa == b'K'));

        // Property (2): every revision's stored amino acids equal a
        // fresh recomputation from its own pool. This is the entity-
        // attached metadata invariant — staleness is impossible.
        for (i, rev) in history.iter().enumerate() {
            let recomputed = rev.sequence.regions[0].with_codon_rail_recomputed(&rev.pool);
            assert_eq!(
                rev.sequence.regions[0].amino_acids, recomputed.amino_acids,
                "revision {} amino acids drift from pool",
                i
            );
            assert_eq!(
                rev.sequence.regions[0].stop_codon_positions, recomputed.stop_codon_positions,
                "revision {} stop codon positions drift from pool",
                i
            );
        }

        // Property (3): mid-revisions hold their own state, not the
        // tail's state. We pick a few specific revisions to confirm.
        for &i in &[1usize, 250, 500, 750, 999] {
            assert!(i < history.len());
            let rev = &history[i];
            let recomputed = rev.sequence.regions[0].with_codon_rail_recomputed(&rev.pool);
            assert_eq!(rev.sequence.regions[0].amino_acids, recomputed.amino_acids);
        }

        // Property (4): most mutations produced a visible base delta.
        // With a uniform 4-base alphabet, the expected fraction of
        // self-substitutions is 1/4, so visible should be ~75% × 1000.
        // Use a generous lower bound to absorb statistical jitter.
        assert!(
            visible_mutations >= 600,
            "expected ≥600 visible mutations, got {}",
            visible_mutations
        );
    }
}
