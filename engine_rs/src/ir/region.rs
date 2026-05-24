use super::{translate_codon, NucHandle, NucleotidePool, Segment, AMINO_STOP};

/// One region of the assembled sequence. Carries its segment role,
/// nucleotide range, and codon frame phase.
///
/// **Codon-rail metadata is no longer stored here.** Earlier versions
/// of the engine kept `amino_acids: Vec<u8>` and `stop_codon_positions:
/// Vec<NucHandle>` on the region as a per-pass cache, but no
/// hot-path consumer ever read them — every contract reads pool bytes
/// directly, AIRR projection re-translates from the raw sequence.
/// Carrying the cache made `Region` ~50 bytes heavier per instance,
/// forced two separate `Simulation::with_indel_*_no_rail_recompute`
/// API variants, and was a steady source of "is this stale?"
/// confusion. The cache is now computed on demand via
/// [`compute_codon_rail`] only when an external consumer (the Python
/// boundary, debug tests) actually needs it.
#[derive(Clone, Debug)]
pub struct Region {
    /// Biological role of this region.
    pub segment: Segment,

    /// Half-open range of nucleotide handles `[start, end)`.
    pub start: NucHandle,
    pub end: NucHandle,

    /// Position within the codon frame at the first nucleotide.
    pub frame_phase: u8,
}

impl Region {
    /// Construct a region with the given segment and nucleotide range.
    pub fn new(segment: Segment, start: NucHandle, end: NucHandle) -> Self {
        Self {
            segment,
            start,
            end,
            frame_phase: 0,
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

    /// Return a new region with `end` advanced to `new_end`.
    pub fn with_end_extended(&self, new_end: NucHandle) -> Self {
        Self {
            end: new_end,
            ..self.clone()
        }
    }

    /// Return a new region with the given frame phase.
    pub fn with_frame_phase(&self, phase: u8) -> Self {
        Self {
            frame_phase: phase,
            ..self.clone()
        }
    }
}

/// Computed codon-rail metadata for one region against a specific
/// nucleotide pool snapshot. Produced on demand by
/// [`compute_codon_rail`]; never stored on `Region` itself.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct CodonRail {
    /// Translated amino acids for codons fully contained within the
    /// region. `b'*'` marks a stop codon; `b'X'` marks an ambiguous
    /// codon (non-canonical base).
    pub amino_acids: Vec<u8>,
    /// First-base handles of every stop codon in this region.
    pub stop_codon_positions: Vec<NucHandle>,
}

impl CodonRail {
    /// An empty codon rail (no codons, no stops).
    pub fn empty() -> Self {
        Self {
            amino_acids: Vec::new(),
            stop_codon_positions: Vec::new(),
        }
    }

    /// Number of stop codons in the rail.
    pub fn stop_codon_count(&self) -> usize {
        self.stop_codon_positions.len()
    }
}

/// Compute the codon rail for `region` against `pool`. Codons start
/// at `region.start + ((3 - frame_phase) % 3)` and proceed in steps
/// of 3 until the codon would overflow `region.end`. Bases outside
/// the canonical `ACGT` alphabet translate to `b'X'`.
pub fn compute_codon_rail(region: &Region, pool: &NucleotidePool) -> CodonRail {
    let skip = (3 - (region.frame_phase as u64)) % 3;
    let start_idx_u64 = (region.start.index() as u64).saturating_add(skip);
    let end_idx_u64 = region.end.index() as u64;

    if start_idx_u64 >= end_idx_u64 {
        return CodonRail::empty();
    }

    let max_codons = ((end_idx_u64 - start_idx_u64) / 3) as usize;
    let mut amino_acids = Vec::with_capacity(max_codons);
    let mut stops = Vec::new();

    let mut i = start_idx_u64;
    while i + 3 <= end_idx_u64 {
        let h0 = NucHandle::new(i as u32);
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

    CodonRail {
        amino_acids,
        stop_codon_positions: stops,
    }
}
