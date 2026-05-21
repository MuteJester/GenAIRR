use super::{translate_codon, NucHandle, NucleotidePool, Segment, AMINO_STOP};

/// One region of the assembled sequence. Carries its segment role,
/// nucleotide range, and derived codon-rail metadata.
#[derive(Clone, Debug)]
pub struct Region {
    /// Biological role of this region.
    pub segment: Segment,

    /// Half-open range of nucleotide handles `[start, end)`.
    pub start: NucHandle,
    pub end: NucHandle,

    /// Position within the codon frame at the first nucleotide.
    pub frame_phase: u8,

    /// Translated amino acids for codons fully contained within this region.
    pub amino_acids: Vec<u8>,

    /// Handles of the first base of every stop codon in this region.
    pub stop_codon_positions: Vec<NucHandle>,
}

impl Region {
    /// Construct a region with the given segment and nucleotide range.
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

    /// Return a new region with codon-rail metadata recomputed from `pool`.
    pub fn with_codon_rail_recomputed(&self, pool: &NucleotidePool) -> Self {
        let skip = (3 - (self.frame_phase as u64)) % 3;
        let start_idx_u64 = (self.start.index() as u64).saturating_add(skip);
        let end_idx_u64 = self.end.index() as u64;

        if start_idx_u64 >= end_idx_u64 {
            return Self {
                amino_acids: Vec::new(),
                stop_codon_positions: Vec::new(),
                ..self.clone()
            };
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

        Self {
            amino_acids,
            stop_codon_positions: stops,
            ..self.clone()
        }
    }

    /// Number of stop codons in this region's already-computed metadata.
    pub fn stop_codon_count(&self) -> usize {
        self.stop_codon_positions.len()
    }
}
