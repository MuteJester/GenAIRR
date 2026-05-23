//! Streaming codon-rail observer attached to `SimulationBuilder`.
//!
//! Today's `Region::with_codon_rail_recomputed` runs as a post-pass
//! after `AssembleSegmentPass`: it re-walks the just-pushed
//! structural region position-by-position to translate each
//! frame-aligned codon and accumulate `amino_acids` /
//! `stop_codon_positions`.
//!
//! `CodonRailObserverState` performs the **same translation state
//! machine** incrementally while `AssembleSegmentPass` is still
//! pushing bases into the pool. Each `on_base_pushed` call slides
//! one byte into a 3-byte codon buffer; the third byte triggers a
//! single `translate_codon` call and updates the rail. At seal
//! time the observer hands back `(amino_acids, stop_codon_positions)`
//! that the assembly pass writes directly into the new `Region`,
//! skipping the post-pass region rebuild entirely.
//!
//! ## Phase 2 scope
//!
//! This observer is the second concrete implementer of
//! [`IrEventObserver`], introduced after the walker observer in
//! Phase 1.5. Its update rule is byte-for-byte equivalent to
//! `Region::with_codon_rail_recomputed`: the observer keeps the
//! original function around (a) as the property-test oracle and
//! (b) for code paths that still rebuild from scratch (post-indel
//! refresh, post-mutation refresh — these will be ported when
//! Phase 4 introduces the corresponding edit events).
//!
//! ## Bit-identical invariant
//!
//! The 500-record `human_tcrb` regression sha256 catches any drift.
//! Three correctness anchors:
//!
//! - The frame-phase skip `(3 - frame_phase) % 3` is computed once
//!   at attach time and frozen for the observer's lifetime; bases at
//!   positions `[start, start+skip)` contribute no codon.
//! - Codons are emitted only when the 3-byte buffer is fully
//!   populated; trailing partial codons are dropped (matches the
//!   `while i + 3 <= end_idx_u64` guard in the from-scratch path).
//! - The `NucHandle` recorded in `stop_codon_positions` is the
//!   handle of the codon's *first* base, not its third. Same as the
//!   from-scratch loop's `h0`.

use super::region::Region;
use super::{translate_codon, NucHandle, Nucleotide, Segment, AMINO_STOP};
use crate::ir::builder::IrEventObserver;

/// State a streaming codon-rail observer accumulates as
/// `SimulationBuilder` pushes one base after another into the pool
/// for an assembled structural region.
///
/// One observer instance covers one region's assembly. The
/// `SimulationBuilder` holds it as `Option<CodonRailObserverState>`.
pub(crate) struct CodonRailObserverState {
    /// Segment this observer's region belongs to. Used to filter
    /// edit events: `on_base_changed` only fires on bases in our
    /// region range. Pure assembly events (`on_base_pushed`) come
    /// in already-filtered by the builder loop, so this is mainly
    /// load-bearing for the Phase 5 edit path.
    segment: Segment,
    /// First pool index of the assembled region — equal to
    /// `region.start.index()`. Captured at attach time; immutable
    /// for the observer's lifetime.
    seq_start: u32,
    /// Number of bytes to discard at region start before the first
    /// frame-aligned codon: `(3 - frame_phase) % 3`. Captured at
    /// attach time. For frame_phase=0 this is 0; for 1 it's 2; for
    /// 2 it's 1.
    skip: u32,
    /// Accumulated amino-acid translations, in codon order.
    /// Equivalent to `region.amino_acids`.
    amino_acids: Vec<u8>,
    /// Handles of the first base of every stop codon observed so
    /// far. Equivalent to `region.stop_codon_positions`.
    stop_codon_positions: Vec<NucHandle>,
    /// In-progress codon buffer. Bytes [0..codon_pos] are populated;
    /// when `codon_pos` reaches 3 we translate and reset.
    codon_bytes: [u8; 3],
    /// Number of bytes filled into `codon_bytes`: 0, 1, or 2.
    codon_pos: u8,
    /// Handle of byte 0 of the in-progress codon. Recorded into
    /// `stop_codon_positions` when the completed codon translates
    /// to `AMINO_STOP`.
    codon_start_handle: NucHandle,
    /// All bytes pushed so far, indexed by `pool_pos - seq_start`.
    /// Maintained to support Phase 5 incremental edits: when
    /// `on_base_changed` fires for a previously-pushed handle, we
    /// need the *other two bytes* of the affected codon to
    /// retranslate. The pool itself isn't readable inside the event
    /// (the observer fires before the mutation lands), so we keep
    /// our own copy.
    ///
    /// Memory cost: one `u8` per pushed base. For a 320-nt TCRB
    /// junction with up to 3 segments this is well under 1 KB per
    /// builder lifetime — negligible.
    pushed_bytes: Vec<u8>,
    /// Phase 16: set when an `on_indel_*` event arrives that the
    /// observer cannot patch in place (an internal indel that
    /// invalidates the codon frame for the affected suffix).
    /// External indels (before the region) shift `seq_start` and
    /// leave this flag clear. The builder consults this flag at
    /// seal time via [`Self::rebuild_if_stale`] and replaces the
    /// observer with a fresh `from_existing_region` rebuild against
    /// the post-mutation simulation when set.
    needs_rebuild: bool,
}

/// Sealed output of a `CodonRailObserverState`: the rail that the
/// assembly pass can write directly into the new `Region` (bypassing
/// the post-pass `with_codon_rail_recomputed` rebuild).
pub(crate) struct SealedCodonRail {
    pub amino_acids: Vec<u8>,
    pub stop_codon_positions: Vec<NucHandle>,
}

impl CodonRailObserverState {
    /// Initialize a fresh codon-rail observer for one region's
    /// assembly. `seq_start` must be `sim.pool.len()` at the moment
    /// of attach (= the soon-to-be-built region's `start.index()`).
    /// `frame_phase` is the region's frame phase (chained from
    /// cumulative-length-mod-3 across prior regions).
    pub(crate) fn new(segment: Segment, seq_start: u32, frame_phase: u8) -> Self {
        let skip = (3u32 - (frame_phase as u32)) % 3;
        Self {
            segment,
            seq_start,
            skip,
            // Capacity hint: structural region of ~300 bases yields
            // ~100 codons. Allocator will round up; saves the first
            // re-grow on typical V/J spans.
            amino_acids: Vec::with_capacity(128),
            stop_codon_positions: Vec::new(),
            codon_bytes: [0; 3],
            codon_pos: 0,
            codon_start_handle: NucHandle::new(0),
            pushed_bytes: Vec::with_capacity(384),
            needs_rebuild: false,
        }
    }

    /// Phase 6: rebuild an observer's internal state from an
    /// already-assembled `Region` plus the current pool slice.
    ///
    /// Used by post-assembly passes (S5F, PCR, quality, contaminant,
    /// indel) that want to attach a codon-rail observer to a
    /// `SimulationBuilder` *after* assembly has finished and torn
    /// down its own observer. The reconstructor:
    ///
    /// 1. Copies the existing `amino_acids` + `stop_codon_positions`
    ///    verbatim (they're already correct).
    /// 2. Copies the region's bytes into the shadow `pushed_bytes`
    ///    buffer so `on_base_changed` can retranslate codons.
    /// 3. Recomputes the in-progress codon state (`codon_bytes`,
    ///    `codon_pos`, `codon_start_handle`) by looking at the trailing
    ///    bytes that haven't yet formed a complete codon — necessary
    ///    if a mutation passes ends up extending the region (not
    ///    typical for S5F but defensive).
    ///
    /// `region.frame_phase` provides the `skip` value; `region.start`
    /// provides `seq_start`. `pool_slice` is the half-open range
    /// `[region.start, region.end)` of bytes the observer needs to
    /// reconstruct internal state. Must have length `region.len()`.
    pub(crate) fn from_existing_region(region: &Region, pool_slice: &[u8]) -> Self {
        assert_eq!(
            pool_slice.len() as u32,
            region.len(),
            "from_existing_region: pool_slice length {} != region.len() {}",
            pool_slice.len(),
            region.len()
        );

        let skip = (3u32 - (region.frame_phase as u32)) % 3;
        let seq_start = region.start.index();
        let amino_acids = region.amino_acids.clone();
        let stop_codon_positions = region.stop_codon_positions.clone();
        let pushed_bytes = pool_slice.to_vec();

        // Reconstruct the in-progress codon buffer: how many bytes
        // past the last complete codon have we seen?
        let coding_bytes = pool_slice.len() as u32;
        let codon_pos = if coding_bytes < skip {
            // Region is shorter than the leading skip — no codons,
            // no in-progress buffer either.
            0u8
        } else {
            ((coding_bytes - skip) % 3) as u8
        };

        let mut codon_bytes = [0u8; 3];
        let mut codon_start_handle = NucHandle::new(seq_start);
        if codon_pos > 0 {
            // Trailing partial codon: copy its bytes into the buffer.
            let codon_byte_0 = coding_bytes - codon_pos as u32;
            codon_start_handle = NucHandle::new(seq_start + codon_byte_0);
            for k in 0..codon_pos {
                codon_bytes[k as usize] =
                    pool_slice[(codon_byte_0 + k as u32) as usize];
            }
        }

        Self {
            segment: region.segment,
            seq_start,
            skip,
            amino_acids,
            stop_codon_positions,
            codon_bytes,
            codon_pos,
            codon_start_handle,
            pushed_bytes,
            needs_rebuild: false,
        }
    }

    /// Drain the observer and return its sealed rail.
    ///
    /// Any trailing partial codon (fewer than 3 bytes after the last
    /// full codon) is dropped, matching the
    /// `while i + 3 <= end_idx_u64` guard in
    /// `Region::with_codon_rail_recomputed`.
    pub(crate) fn seal(self) -> SealedCodonRail {
        SealedCodonRail {
            amino_acids: self.amino_acids,
            stop_codon_positions: self.stop_codon_positions,
        }
    }

    /// Segment the observer is attached to. Used by the builder to
    /// pair a sealed rail with the matching `Region` for write-back.
    pub(crate) fn observed_segment(&self) -> Segment {
        self.segment
    }

    /// Handle at which this observer's region begins. Reflects any
    /// `on_indel_*` shift updates the observer has absorbed.
    pub(crate) fn observed_start_handle(&self) -> NucHandle {
        NucHandle::new(self.seq_start)
    }

    /// Phase 16: was the observer's state invalidated by an internal
    /// indel that the in-place patch path can't safely handle?
    /// `true` means [`Self::rebuild_if_stale`] will replace the
    /// observer with a fresh `from_existing_region` rebuild at
    /// seal time.
    pub(crate) fn needs_rebuild(&self) -> bool {
        self.needs_rebuild
    }

    /// Phase 16: rebuild the observer from the post-mutation
    /// simulation when [`Self::needs_rebuild`] is set, otherwise
    /// return self unchanged. The rebuild uses
    /// [`Self::from_existing_region`] against the region currently
    /// assigned to the observer's segment in `sim` (region.start
    /// has typically shifted from attach time, so we match by
    /// segment rather than by stored start handle).
    pub(crate) fn rebuild_if_stale(
        self,
        sim: &crate::ir::Simulation,
    ) -> Self {
        if !self.needs_rebuild {
            return self;
        }
        let Some(region) = sim
            .sequence
            .regions
            .iter()
            .find(|r| r.segment == self.segment)
            .cloned()
        else {
            // Region disappeared (shouldn't happen for V/D/J in
            // production paths). Leave the observer as-is so the
            // seal-time region-write-back lookup safely fails.
            return self;
        };
        let start = region.start.index() as usize;
        let end = region.end.index() as usize;
        let pool_slice: Vec<u8> = sim
            .pool
            .as_slice()
            .get(start..end)
            .expect("rebuild_if_stale: region range outside pool")
            .iter()
            .map(|n| n.base)
            .collect();
        Self::from_existing_region(&region, &pool_slice)
    }
}

impl IrEventObserver for CodonRailObserverState {
    fn on_base_pushed(&mut self, handle: NucHandle, n: &Nucleotide) {
        // Phase 7: when multiple codon-rail observers are attached
        // (one per existing region during mutation passes), only the
        // observer whose segment matches the incoming nucleotide
        // should react. During single-observer assembly this is a
        // no-op because the caller attaches exactly one observer
        // for one segment and never pushes cross-segment bytes.
        if n.segment != self.segment {
            return;
        }
        // Record the byte in our own buffer so Phase 5
        // `on_base_changed` can retranslate the affected codon
        // without reaching into the pool (which hasn't been
        // mutated yet by the time the event fires).
        self.pushed_bytes.push(n.base);

        let pos = handle.index();
        // Bytes before the first frame-aligned codon contribute
        // nothing. Equivalent to skipping the first `skip` indices
        // of the structural region in the from-scratch loop.
        if pos < self.seq_start + self.skip {
            return;
        }

        // First byte of a new codon: record its handle so we can
        // attribute a stop to the codon's start position later.
        if self.codon_pos == 0 {
            self.codon_start_handle = handle;
        }

        self.codon_bytes[self.codon_pos as usize] = n.base;
        self.codon_pos += 1;

        if self.codon_pos == 3 {
            let aa = translate_codon(self.codon_bytes[0], self.codon_bytes[1], self.codon_bytes[2]);
            self.amino_acids.push(aa);
            if aa == AMINO_STOP {
                self.stop_codon_positions.push(self.codon_start_handle);
            }
            self.codon_pos = 0;
        }
    }

    /// Phase 5: incremental codon-rail patch on a base change.
    ///
    /// When the byte at `handle` flips from `old_n.base` to
    /// `new_base`, only the codon *containing* that byte can
    /// change in the amino-acid rail. We:
    ///
    /// 1. Filter out edits outside our region's segment / range.
    /// 2. Locate the codon index `(rel - skip) / 3` if and only if
    ///    the position is inside a fully-translated codon (i.e.
    ///    not in the leading `skip` bytes, not in a trailing
    ///    partial codon).
    /// 3. Patch our own `pushed_bytes` shadow buffer with the new
    ///    base, retranslate the three bytes of that codon, and
    ///    update `amino_acids[codon_idx]` plus the
    ///    `stop_codon_positions` set accordingly.
    ///
    /// This is the architectural payoff of Phase 5 for the rail:
    /// a single byte change costs one `translate_codon` call and
    /// at most one membership update of `stop_codon_positions`,
    /// instead of `Region::with_codon_rail_recomputed` re-walking
    /// the entire region.
    fn on_base_changed(&mut self, handle: NucHandle, old_n: &Nucleotide, new_base: u8) {
        if old_n.segment != self.segment {
            return;
        }
        let pos = handle.index();
        if pos < self.seq_start {
            return;
        }
        let rel = pos - self.seq_start;
        if (rel as usize) >= self.pushed_bytes.len() {
            // Position is outside the bytes we've observed —
            // either past our seal window or never pushed through us.
            return;
        }
        if old_n.base == new_base {
            return;
        }

        // Update our shadow buffer to reflect the new byte.
        self.pushed_bytes[rel as usize] = new_base;

        // Was this position inside a fully-translated codon? The
        // leading `skip` bytes never form a codon; bytes past the
        // last `(amino_acids.len() * 3)`-th codon byte are still
        // in `codon_pos > 0` partial state and not yet in
        // `amino_acids`.
        if rel < self.skip {
            // The codon_pos partial buffer also tracks this byte
            // if it's currently in-flight, but the skip rule says
            // no — those leading bytes never enter codon_bytes.
            return;
        }
        let codon_relative_byte = rel - self.skip;
        let codon_idx = (codon_relative_byte / 3) as usize;
        if codon_idx >= self.amino_acids.len() {
            // The codon containing this byte hasn't been fully
            // translated yet — it's either in the in-progress
            // `codon_bytes` buffer or a future codon. Patch the
            // in-progress buffer if applicable; otherwise no rail
            // update is needed (the next on_base_pushed will see
            // the new value via our shadow buffer when we get to
            // it).
            let pos_in_codon = (codon_relative_byte % 3) as usize;
            if pos_in_codon < self.codon_pos as usize {
                // The in-progress codon already contains a stale
                // byte; patch it.
                self.codon_bytes[pos_in_codon] = new_base;
            }
            return;
        }

        // Recompute the 3 bytes of the affected codon from the
        // shadow buffer and update amino_acids + stop set.
        let codon_byte_0 = self.skip + (codon_idx as u32) * 3;
        let b0 = self.pushed_bytes[codon_byte_0 as usize];
        let b1 = self.pushed_bytes[(codon_byte_0 + 1) as usize];
        let b2 = self.pushed_bytes[(codon_byte_0 + 2) as usize];
        let new_aa = translate_codon(b0, b1, b2);
        let old_aa = self.amino_acids[codon_idx];
        if new_aa == old_aa {
            return;
        }
        self.amino_acids[codon_idx] = new_aa;

        let codon_start_handle = NucHandle::new(self.seq_start + codon_byte_0);
        match (old_aa == AMINO_STOP, new_aa == AMINO_STOP) {
            (false, false) | (true, true) => {}
            (true, false) => {
                // Codon was a stop, no longer is — remove from set.
                if let Some(idx) = self
                    .stop_codon_positions
                    .iter()
                    .position(|h| *h == codon_start_handle)
                {
                    self.stop_codon_positions.remove(idx);
                }
            }
            (false, true) => {
                // Codon was not a stop, now is — insert maintaining
                // ascending handle order to match the from-scratch
                // walker's emission order.
                let insert_at = self
                    .stop_codon_positions
                    .iter()
                    .position(|h| h.index() > codon_start_handle.index())
                    .unwrap_or(self.stop_codon_positions.len());
                self.stop_codon_positions
                    .insert(insert_at, codon_start_handle);
            }
        }
    }

    /// Phase 16: handle an insertion event.
    ///
    /// Any indel before-or-inside our region invalidates state and
    /// triggers a `from_existing_region` rebuild at seal time —
    /// an external indel shifts upstream cumulative length and so
    /// changes our `frame_phase` / `skip`, and an internal indel
    /// changes our bytes. Only indels strictly after our right
    /// boundary leave our rail unaffected.
    fn on_indel_inserted(&mut self, at: u32, _n: &Nucleotide) {
        if self.needs_rebuild {
            return;
        }
        let region_end = self.seq_start + self.pushed_bytes.len() as u32;
        if at < region_end {
            self.needs_rebuild = true;
        }
    }

    /// Phase 16: handle a deletion event. Same rebuild rule as
    /// insertion — any deletion before-or-inside our region
    /// requires re-deriving the rail.
    fn on_indel_deleted(&mut self, at: u32, _removed: &Nucleotide) {
        if self.needs_rebuild {
            return;
        }
        let region_end = self.seq_start + self.pushed_bytes.len() as u32;
        if at < region_end {
            self.needs_rebuild = true;
        }
    }
}

impl Region {
    /// Construct a `Region` whose codon rail comes from a streaming
    /// observer's sealed state instead of a post-pass rebuild.
    ///
    /// Bit-identical to the chain
    /// `Region::new(...).with_frame_phase(p).with_codon_rail_recomputed(pool)`
    /// on the same range when the observer received every nucleotide
    /// in that range via [`CodonRailObserverState::on_base_pushed`].
    pub(crate) fn from_sealed_codon_rail(
        segment: Segment,
        start: NucHandle,
        end: NucHandle,
        frame_phase: u8,
        rail: SealedCodonRail,
    ) -> Self {
        Self {
            segment,
            start,
            end,
            frame_phase,
            amino_acids: rail.amino_acids,
            stop_codon_positions: rail.stop_codon_positions,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ir::flag;

    fn n(base: u8) -> Nucleotide {
        Nucleotide::synthetic(base, Segment::V, flag::N_NUC)
    }

    /// Drive the observer over a byte sequence with the given
    /// frame_phase starting at `seq_start`, return the sealed rail.
    fn drive(seq_start: u32, frame_phase: u8, bases: &[u8]) -> SealedCodonRail {
        let mut obs = CodonRailObserverState::new(Segment::V, seq_start, frame_phase);
        for (i, &b) in bases.iter().enumerate() {
            let h = NucHandle::new(seq_start + i as u32);
            IrEventObserver::on_base_pushed(&mut obs, h, &n(b));
        }
        obs.seal()
    }

    #[test]
    fn empty_region_yields_empty_rail() {
        let rail = drive(10, 0, b"");
        assert!(rail.amino_acids.is_empty());
        assert!(rail.stop_codon_positions.is_empty());
    }

    #[test]
    fn frame_phase_zero_no_skip() {
        let rail = drive(0, 0, b"ATGAAA");
        assert_eq!(rail.amino_acids, b"MK");
        assert!(rail.stop_codon_positions.is_empty());
    }

    #[test]
    fn frame_phase_one_skips_two_leading_bytes() {
        // skip = (3 - 1) % 3 = 2
        let rail = drive(0, 1, b"XXATGAAA");
        assert_eq!(rail.amino_acids, b"MK");
    }

    #[test]
    fn frame_phase_two_skips_one_leading_byte() {
        // skip = (3 - 2) % 3 = 1
        let rail = drive(0, 2, b"XATGAAA");
        assert_eq!(rail.amino_acids, b"MK");
    }

    #[test]
    fn trailing_partial_codon_is_dropped() {
        let rail = drive(0, 0, b"ATGAA"); // 5 bases: ATG + AA (incomplete)
        assert_eq!(rail.amino_acids, b"M");
    }

    #[test]
    fn stop_codon_position_is_codon_start_handle() {
        // TAA at offset 3 → stop, codon_start_handle = NucHandle::new(13)
        let rail = drive(10, 0, b"ATGTAA");
        assert_eq!(rail.amino_acids, b"M*");
        assert_eq!(rail.stop_codon_positions, vec![NucHandle::new(13)]);
    }

    /// Drive the observer over `bases`, apply `(rel_pos, new_base)`
    /// edits, return the sealed rail.
    fn drive_with_edits(
        seq_start: u32,
        frame_phase: u8,
        bases: &[u8],
        edits: &[(u32, u8)],
    ) -> SealedCodonRail {
        let mut obs = CodonRailObserverState::new(Segment::V, seq_start, frame_phase);
        for (i, &b) in bases.iter().enumerate() {
            let h = NucHandle::new(seq_start + i as u32);
            IrEventObserver::on_base_pushed(&mut obs, h, &n(b));
        }
        for &(rel_pos, new_base) in edits {
            let h = NucHandle::new(seq_start + rel_pos);
            let old = n(bases[rel_pos as usize]);
            IrEventObserver::on_base_changed(&mut obs, h, &old, new_base);
        }
        obs.seal()
    }

    /// Apply the same edits to `bases` and drive a fresh observer
    /// over the mutated sequence. The two should agree.
    fn apply_edits(bases: &[u8], edits: &[(u32, u8)]) -> Vec<u8> {
        let mut mutated = bases.to_vec();
        for &(rel_pos, new_base) in edits {
            mutated[rel_pos as usize] = new_base;
        }
        mutated
    }

    #[test]
    fn change_base_no_op_when_old_equals_new() {
        let bases = b"ATGAAACGTTAA";
        let edits = &[(0, b'A')]; // pos 0 already 'A'
        let rail = drive_with_edits(0, 0, bases, edits);
        let fresh = drive(0, 0, bases);
        assert_eq!(rail.amino_acids, fresh.amino_acids);
        assert_eq!(rail.stop_codon_positions, fresh.stop_codon_positions);
    }

    #[test]
    fn change_base_in_middle_codon_retranslates_one_codon() {
        // ATG AAA CGT TAA → change pos 4 (A in AAA, middle byte)
        // to G → "AGA" → R. Stop codon unchanged.
        let bases = b"ATGAAACGTTAA";
        let edits: &[(u32, u8)] = &[(4, b'G')];
        let rail = drive_with_edits(0, 0, bases, edits);
        let mutated = apply_edits(bases, edits);
        let fresh = drive(0, 0, &mutated);
        assert_eq!(rail.amino_acids, fresh.amino_acids);
        assert_eq!(rail.stop_codon_positions, fresh.stop_codon_positions);
    }

    #[test]
    fn change_base_creates_new_stop() {
        // ATG AAA → change AAA to TAA via two edits → introduces a stop
        // at codon 1.
        let bases = b"ATGAAA";
        let edits: &[(u32, u8)] = &[(3, b'T')]; // AAA → TAA = stop
        let rail = drive_with_edits(0, 0, bases, edits);
        let mutated = apply_edits(bases, edits);
        let fresh = drive(0, 0, &mutated);
        assert_eq!(rail.amino_acids, fresh.amino_acids);
        assert_eq!(
            rail.stop_codon_positions, fresh.stop_codon_positions,
            "stop insertion order diverged"
        );
        assert_eq!(rail.stop_codon_positions.len(), 1);
    }

    #[test]
    fn change_base_removes_stop() {
        // ATG TAA → change pos 3 from T → A so TAA becomes AAA.
        // Stop at codon 1 should disappear.
        let bases = b"ATGTAA";
        let edits: &[(u32, u8)] = &[(3, b'A')];
        let rail = drive_with_edits(0, 0, bases, edits);
        let mutated = apply_edits(bases, edits);
        let fresh = drive(0, 0, &mutated);
        assert_eq!(rail.amino_acids, fresh.amino_acids);
        assert_eq!(rail.stop_codon_positions, fresh.stop_codon_positions);
        assert!(rail.stop_codon_positions.is_empty());
    }

    #[test]
    fn change_base_multiple_stops_inserted_in_handle_order() {
        // ATGAAACGTAAA → change pos 0 (A→T) and pos 6 (C→T):
        //   codon 0: TTG → L
        //   codon 1: AAA → K (unchanged)
        //   codon 2: TGT → C (was R: CGT)
        //   codon 3: AAA → K (unchanged)
        // Then apply third edit: pos 3 (A→T) AAA→TAA = stop at codon 1.
        // Final stops: only codon 1.
        let bases = b"ATGAAACGTAAA";
        let edits: &[(u32, u8)] = &[(0, b'T'), (6, b'T'), (3, b'T')];
        let rail = drive_with_edits(0, 0, bases, edits);
        let mutated = apply_edits(bases, edits);
        let fresh = drive(0, 0, &mutated);
        assert_eq!(rail.amino_acids, fresh.amino_acids);
        assert_eq!(rail.stop_codon_positions, fresh.stop_codon_positions);
    }

    #[test]
    fn change_base_outside_observer_range_ignored() {
        // Observer attached at seq_start=10. Edit at handle 5 should be ignored.
        let bases = b"ATGAAA";
        let mut obs = CodonRailObserverState::new(Segment::V, 10, 0);
        for (i, &b) in bases.iter().enumerate() {
            IrEventObserver::on_base_pushed(&mut obs, NucHandle::new(10 + i as u32), &n(b));
        }
        let irrelevant_old = n(b'A');
        IrEventObserver::on_base_changed(
            &mut obs,
            NucHandle::new(5), // before our seq_start
            &irrelevant_old,
            b'T',
        );
        let rail = obs.seal();
        let fresh = drive(10, 0, bases);
        assert_eq!(rail.amino_acids, fresh.amino_acids);
        assert_eq!(rail.stop_codon_positions, fresh.stop_codon_positions);
    }

    #[test]
    fn change_base_during_in_progress_codon_patches_buffer() {
        // Push 5 bases of ATGAAA (5 of 6 — one codon plus first 2 bytes
        // of second). Edit pos 4 (the partially-pushed byte) before pushing
        // pos 5. Result should match driving fresh over ATGAA(edited)A.
        let mut obs = CodonRailObserverState::new(Segment::V, 0, 0);
        for (i, &b) in b"ATGAA".iter().enumerate() {
            IrEventObserver::on_base_pushed(&mut obs, NucHandle::new(i as u32), &n(b));
        }
        // Edit pos 4 (second 'A' in second codon, in-flight) from A → C.
        IrEventObserver::on_base_changed(&mut obs, NucHandle::new(4), &n(b'A'), b'C');
        // Push pos 5 (final byte of second codon).
        IrEventObserver::on_base_pushed(&mut obs, NucHandle::new(5), &n(b'A'));
        let rail = obs.seal();
        let fresh = drive(0, 0, b"ATGACA");
        assert_eq!(rail.amino_acids, fresh.amino_acids);
        assert_eq!(rail.stop_codon_positions, fresh.stop_codon_positions);
    }

    /// Build a `Region` populated with codon-rail data by running
    /// `with_codon_rail_recomputed` over a freshly-built pool.
    fn region_from(seq_start: u32, frame_phase: u8, bases: &[u8]) -> (Region, Vec<u8>) {
        use super::super::pool::NucleotidePool;
        let mut pool = NucleotidePool::with_capacity(seq_start as usize + bases.len());
        for _ in 0..seq_start {
            pool.push(n(b'A')); // pad with dummy bases so handles align
        }
        for &b in bases {
            pool.push(n(b));
        }
        let region = Region::new(
            Segment::V,
            NucHandle::new(seq_start),
            NucHandle::new(seq_start + bases.len() as u32),
        )
        .with_frame_phase(frame_phase)
        .with_codon_rail_recomputed(&pool);
        (region, bases.to_vec())
    }

    #[test]
    fn from_existing_region_matches_fresh_drive_at_seal() {
        // After rebuilding from a Region that was fully assembled, the
        // observer's seal output should equal the Region's own rail.
        for frame_phase in 0..=2u8 {
            for bases in [
                &b"ATGAAACGTTAA"[..],
                &b"ATGGGG"[..],
                &b"TAGTAATAA"[..],
                &b""[..],
            ] {
                let (region, slice) = region_from(0, frame_phase, bases);
                let obs = CodonRailObserverState::from_existing_region(&region, &slice);
                let rail = obs.seal();
                assert_eq!(
                    rail.amino_acids, region.amino_acids,
                    "frame_phase={} bases={:?}",
                    frame_phase, bases
                );
                assert_eq!(
                    rail.stop_codon_positions, region.stop_codon_positions,
                    "frame_phase={} bases={:?}",
                    frame_phase, bases
                );
            }
        }
    }

    #[test]
    fn from_existing_region_supports_subsequent_on_base_changed() {
        // Rebuild + apply edit → result should match fresh observer
        // driven over the equivalent mutated byte stream.
        let bases = b"ATGAAACGTAAA";
        let edits: &[(u32, u8)] = &[(3, b'T')]; // codon 1 AAA → TAA = stop

        // Rebuild path
        let (region, slice) = region_from(0, 0, bases);
        let mut obs = CodonRailObserverState::from_existing_region(&region, &slice);
        for &(rel_pos, new_base) in edits {
            let h = NucHandle::new(rel_pos);
            let old = n(bases[rel_pos as usize]);
            IrEventObserver::on_base_changed(&mut obs, h, &old, new_base);
        }
        let rail_via_rebuild = obs.seal();

        // Fresh path
        let mut mutated = bases.to_vec();
        for &(rel_pos, new_base) in edits {
            mutated[rel_pos as usize] = new_base;
        }
        let rail_fresh = drive(0, 0, &mutated);

        assert_eq!(rail_via_rebuild.amino_acids, rail_fresh.amino_acids);
        assert_eq!(
            rail_via_rebuild.stop_codon_positions,
            rail_fresh.stop_codon_positions
        );
    }

    #[test]
    fn from_existing_region_preserves_in_progress_codon_state() {
        // Region of 8 bases with frame_phase=0 → 2 complete codons +
        // 2 trailing bytes in the in-progress buffer. After rebuild,
        // pushing one more byte should complete a third codon
        // matching what fresh observation would produce.
        let bases = b"ATGAAACG"; // 2 codons (ATG, AAA) + partial "CG"
        let (region, slice) = region_from(0, 0, bases);
        let mut obs = CodonRailObserverState::from_existing_region(&region, &slice);
        // Push one more byte to complete the third codon as CGT → R
        IrEventObserver::on_base_pushed(&mut obs, NucHandle::new(8), &n(b'T'));
        let rail_via_rebuild = obs.seal();

        let rail_fresh = drive(0, 0, b"ATGAAACGT");
        assert_eq!(rail_via_rebuild.amino_acids, rail_fresh.amino_acids);
        assert_eq!(
            rail_via_rebuild.stop_codon_positions,
            rail_fresh.stop_codon_positions
        );
    }

    #[test]
    fn from_existing_region_at_nonzero_seq_start() {
        // Region not starting at pool index 0 — common after upstream
        // V/D regions have already been assembled.
        let bases = b"ATGAAACGT";
        let (region, slice) = region_from(42, 0, bases);
        assert_eq!(region.start.index(), 42);
        let mut obs = CodonRailObserverState::from_existing_region(&region, &slice);

        // Apply edit at absolute handle 45 (rel_pos 3 inside region)
        IrEventObserver::on_base_changed(
            &mut obs,
            NucHandle::new(45),
            &n(b'A'),
            b'T', // codon 1 AAA → TAA = stop
        );
        let rail = obs.seal();
        let fresh = drive(42, 0, b"ATGTAACGT");
        assert_eq!(rail.amino_acids, fresh.amino_acids);
        assert_eq!(rail.stop_codon_positions, fresh.stop_codon_positions);
    }

    #[test]
    fn observer_matches_with_codon_rail_recomputed() {
        use super::super::pool::NucleotidePool;

        // Build a pool with a recognisable sequence; assemble a
        // region over the whole thing and compare both code paths.
        for frame_phase in 0..=2u8 {
            for bases in [
                &b"ATGAAACGTTAA"[..],
                &b"ATGGGGCCC"[..],
                &b"TAGTAATAA"[..],
                &b"ATG"[..],
                &b""[..],
            ] {
                let mut pool = NucleotidePool::with_capacity(64);
                for &b in bases {
                    pool.push(n(b));
                }
                let region = Region::new(
                    Segment::V,
                    NucHandle::new(0),
                    NucHandle::new(bases.len() as u32),
                )
                .with_frame_phase(frame_phase)
                .with_codon_rail_recomputed(&pool);

                let rail = drive(0, frame_phase, bases);
                assert_eq!(
                    rail.amino_acids, region.amino_acids,
                    "frame_phase={} bases={:?}",
                    frame_phase, bases
                );
                assert_eq!(
                    rail.stop_codon_positions, region.stop_codon_positions,
                    "frame_phase={} bases={:?}",
                    frame_phase, bases
                );
            }
        }
    }
}
