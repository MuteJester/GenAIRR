//! Allocation-free, fixed-capacity views of ordered edits. Original nucleotide
//! identities, rather than just bytes, survive the mapping.
use crate::assignment::AlleleInstance;
use crate::ir::{flag, Nucleotide, PoolRange, Segment, Simulation};
use crate::productivity::{codon_stop, SequenceView};

use super::event::IndelEvent;

#[derive(Clone, Copy, Debug)]
enum Piece {
    Original { start: usize, end: usize },
    Inserted(Nucleotide),
}

impl Piece {
    fn len(self) -> usize {
        match self {
            Self::Original { start, end } => end - start,
            Self::Inserted(_) => 1,
        }
    }
}

/// Indexed once for a tuple root, shared by every hypothetical pair.
pub(super) struct EditIndex<'a> {
    original: &'a Simulation,
    by_reference: [Vec<(u16, usize)>; Segment::COUNT],
    by_position: [Vec<(usize, u16)>; Segment::COUNT],
    stops: [Vec<u32>; 3],
    unknown: [Vec<u32>; 3],
}

impl<'a> EditIndex<'a> {
    pub(super) fn new(original: &'a Simulation) -> Self {
        let mut index = Self {
            original,
            by_reference: std::array::from_fn(|_| Vec::new()),
            by_position: std::array::from_fn(|_| Vec::new()),
            stops: std::array::from_fn(|_| vec![0]),
            unknown: std::array::from_fn(|_| vec![0]),
        };
        for (p, n) in original.pool.as_slice().iter().enumerate() {
            if let Some(g) = n.germline_pos.get() {
                index.by_reference[n.segment.index()].push((g, p));
                index.by_position[n.segment.index()].push((p, g));
            }
        }
        for positions in &mut index.by_reference {
            positions.sort_unstable();
        }
        let pool = original.pool.as_slice();
        for phase in 0..3 {
            for p in (phase..pool.len().saturating_sub(2)).step_by(3) {
                let result = codon_stop([pool[p].base, pool[p + 1].base, pool[p + 2].base]);
                index.stops[phase]
                    .push(index.stops[phase].last().unwrap() + u32::from(result == Some(true)));
                index.unknown[phase]
                    .push(index.unknown[phase].last().unwrap() + u32::from(result.is_none()));
            }
        }
        index
    }

    pub(super) fn view(&self) -> EditView<'_> {
        self.view_with_capacity::<5>()
    }

    pub(super) fn view_with_capacity<const CAPACITY: usize>(&self) -> EditView<'_, CAPACITY> {
        assert!(CAPACITY > 0);
        EditView {
            index: self,
            pieces: [Piece::Original {
                start: 0,
                end: self.original.pool.len(),
            }; CAPACITY],
            piece_count: usize::from(!self.original.pool.is_empty()),
            length: self.original.pool.len(),
            regions: std::array::from_fn(|i| {
                self.original
                    .sequence
                    .regions
                    .iter()
                    .find(|r| r.segment.index() == i)
                    .map(|r| PoolRange::new(r.start.index(), r.end.index()))
            }),
        }
    }
}

#[derive(Clone, Copy)]
pub(super) struct EditView<'a, const CAPACITY: usize = 5> {
    index: &'a EditIndex<'a>,
    pieces: [Piece; CAPACITY],
    piece_count: usize,
    length: usize,
    regions: [Option<PoolRange>; Segment::COUNT],
}

impl<const CAPACITY: usize> EditView<'_, CAPACITY> {
    fn push(&mut self, piece: Piece) {
        if piece.len() == 0 {
            return;
        }
        if let Piece::Original { start, end } = piece {
            if self.piece_count > 0 {
                if let Piece::Original {
                    end: previous_end, ..
                } = &mut self.pieces[self.piece_count - 1]
                {
                    if *previous_end == start {
                        *previous_end = end;
                        return;
                    }
                }
            }
        }
        assert!(
            self.piece_count < self.pieces.len(),
            "EditView capacity must hold at most two pieces per edit plus the original span"
        );
        self.pieces[self.piece_count] = piece;
        self.piece_count += 1;
    }

    fn append_range(&mut self, source: &Self, start: usize, end: usize) {
        let mut offset = 0;
        for &piece in &source.pieces[..source.piece_count] {
            let right = offset + piece.len();
            let lo = start.max(offset);
            let hi = end.min(right);
            if lo < hi {
                self.push(match piece {
                    Piece::Original {
                        start: original, ..
                    } => Piece::Original {
                        start: original + lo - offset,
                        end: original + hi - offset,
                    },
                    Piece::Inserted(n) => Piece::Inserted(n),
                });
            }
            offset = right;
        }
    }

    fn insertion_segment(&self, site: usize) -> Segment {
        for segment in [
            Segment::V,
            Segment::Np1,
            Segment::D,
            Segment::Np2,
            Segment::J,
        ] {
            if self.regions[segment.index()].is_some_and(|r| r.contains(site as u32)) {
                return segment;
            }
        }
        self.nucleotide(site)
            .or_else(|| site.checked_sub(1).and_then(|p| self.nucleotide(p)))
            .map(|n| n.segment)
            .unwrap_or(Segment::V)
    }

    pub(super) fn edited(&self, event: IndelEvent) -> Self {
        let (at, removed, inserted) = match event {
            IndelEvent::Insertion { site, base } => (
                site as usize,
                0,
                Some(Nucleotide::synthetic(
                    base,
                    self.insertion_segment(site as usize),
                    flag::INDEL_INSERTED,
                )),
            ),
            IndelEvent::Deletion { site: Some(site) } => (site as usize, 1, None),
            IndelEvent::Deletion { site: None } | IndelEvent::NoOp => return *self,
        };
        let mut result = *self;
        result.piece_count = 0;
        result.append_range(self, 0, at);
        if let Some(n) = inserted {
            result.push(Piece::Inserted(n));
        }
        result.append_range(self, at + removed, self.length);
        result.length = self.length + usize::from(inserted.is_some()) - removed;
        for region in result.regions.iter_mut().flatten() {
            *region = if inserted.is_some() {
                region.after_insertion(at as u32)
            } else {
                region.after_deletion(at as u32)
            };
        }
        result
    }

    fn mapped_original(&self, original: usize) -> Option<usize> {
        let mut offset = 0;
        for &piece in &self.pieces[..self.piece_count] {
            if let Piece::Original { start, end } = piece {
                if start <= original && original < end {
                    return Some(offset + original - start);
                }
            }
            offset += piece.len();
        }
        None
    }
}

impl<const CAPACITY: usize> SequenceView for EditView<'_, CAPACITY> {
    fn len(&self) -> usize {
        self.length
    }
    fn nucleotide(&self, position: usize) -> Option<Nucleotide> {
        let mut offset = 0;
        for &piece in &self.pieces[..self.piece_count] {
            if position < offset + piece.len() {
                return Some(match piece {
                    Piece::Original { start, .. } => {
                        self.index.original.pool.as_slice()[start + position - offset]
                    }
                    Piece::Inserted(n) => n,
                });
            }
            offset += piece.len();
        }
        None
    }
    fn region(&self, segment: Segment) -> Option<(usize, usize)> {
        self.regions[segment.index()].map(|r| (r.start as usize, r.end as usize))
    }
    fn assignment(&self, segment: Segment) -> Option<AlleleInstance> {
        self.index.original.assignments.get(segment).copied()
    }
    fn germline_position(&self, segment: Segment, germline: u16) -> Option<usize> {
        let range = self.regions[segment.index()]?;
        let positions = &self.index.by_reference[segment.index()];
        let mut i = positions.partition_point(|&(g, _)| g < germline);
        while let Some(&(g, p)) = positions.get(i) {
            if g != germline {
                break;
            }
            if let Some(mapped) = self
                .mapped_original(p)
                .filter(|&p| range.contains(p as u32))
            {
                return Some(mapped);
            }
            i += 1;
        }
        None
    }
    fn last_germline(&self, segment: Segment) -> Option<(usize, u16)> {
        let range = self.regions[segment.index()]?;
        let positions = &self.index.by_position[segment.index()];
        let mut offset = self.length;
        for &piece in self.pieces[..self.piece_count].iter().rev() {
            offset -= piece.len();
            let Piece::Original { start, end } = piece else {
                continue;
            };
            let lo = (range.start as usize).max(offset);
            let hi = (range.end as usize).min(offset + end - start);
            if lo >= hi {
                continue;
            }
            let upper = start + hi - offset;
            let lower = start + lo - offset;
            let i = positions.partition_point(|&(p, _)| p < upper);
            if i > 0 {
                let (p, g) = positions[i - 1];
                if p >= lower {
                    return Some((offset + p - start, g));
                }
            }
        }
        None
    }
    fn coding_stop(&self, start: usize, end: usize) -> Option<bool> {
        let end = end.min(self.length);
        if end < start.saturating_add(3) {
            return None;
        }
        let mut p = start;
        let mut unknown = false;
        while p + 2 < end {
            let mut offset = 0;
            let mut unchanged = None;
            for &piece in &self.pieces[..self.piece_count] {
                let right = offset + piece.len();
                if offset <= p && p + 2 < right {
                    if let Piece::Original { start, .. } = piece {
                        unchanged = Some((start + p - offset, right.min(end)));
                    }
                    break;
                }
                offset = right;
            }
            if let Some((original, right)) = unchanged {
                let count = (right - p) / 3;
                let phase = original % 3;
                let lo = original / 3;
                let hi = lo + count;
                if self.index.stops[phase][hi] > self.index.stops[phase][lo] {
                    return Some(true);
                }
                unknown |= self.index.unknown[phase][hi] > self.index.unknown[phase][lo];
                p += count * 3;
            } else {
                match codon_stop([p, p + 1, p + 2].map(|i| self.nucleotide(i).unwrap().base)) {
                    Some(true) => return Some(true),
                    None => unknown = true,
                    Some(false) => {}
                }
                p += 3;
            }
        }
        (!unknown).then_some(false)
    }
}
