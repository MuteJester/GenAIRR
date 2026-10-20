//! Experimental local two-end repair. Both strands use one coordinate axis:
//! top is 5'->3', bottom is 3'->5'. No rates or enzyme identities are inferred.
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

pub const MODEL: &str = "experimental_duplex_repair_v1";
pub mod bindings;
mod coupled;
mod evidence;
pub mod model;
pub mod observation;
pub mod projection;
pub mod seed;
pub use coupled::{
    repair_coupled, repair_coupled_annotated, CoupledOrigin, CoupledProduct, JoinName, JoinScript,
    LigationSeam, Substrate, COUPLED_MODEL,
};
pub use evidence::{
    AnnotatedCoupledProduct, CodingEnd, CoupledEvidence, OperationId, PhysicalLocation,
    PhysicalTip, Placement, PlacementRole, RepairEvidenceEvent,
};
use evidence::{EvidenceObserver, OperationObserver};

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum Side {
    Left,
    Right,
}
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum Strand {
    Top,
    Bottom,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
#[serde(tag = "operation", rename_all = "snake_case", deny_unknown_fields)]
pub enum RepairOperation {
    Open {
        side: Side,
        offset: usize,
    },
    Resect {
        side: Side,
        strand: Strand,
        count: usize,
    },
    AddN {
        side: Side,
        sequence: String,
    },
    Synapse {
        right_start: i64,
    },
    Fill {
        strand: Strand,
    },
    Ligate {},
}

#[derive(Clone, Debug, Serialize, PartialEq, Eq)]
pub struct Origin {
    pub side: Side,
    pub strand: Strand,
    pub reference_position: usize,
}

#[derive(Clone, Debug, Serialize, PartialEq, Eq)]
pub struct Base<O = Origin> {
    pub id: u64,
    pub base: char,
    pub kind: String,
    pub original: Option<O>,
    pub template_id: Option<u64>,
}

#[derive(Clone, Debug)]
struct End<O = Origin> {
    top: BTreeMap<i64, Base<O>>,
    bottom: BTreeMap<i64, Base<O>>,
    length: usize,
    original_start: i64,
    opened: bool,
    /// The previous ligation's two strand boundaries. Operations at the other
    /// end must retain both bases adjacent to those bonds in this local model.
    barrier: Option<(i64, i64)>,
}

#[derive(Clone, Debug, Serialize, PartialEq, Eq)]
pub struct Column<O = Origin> {
    pub position: i64,
    pub top: Base<O>,
    pub bottom: Base<O>,
}

#[derive(Clone, Debug, Serialize, PartialEq, Eq)]
pub struct RepairProduct {
    pub model: String,
    pub top_direction: String,
    pub bottom_direction: String,
    pub sequence: String,
    /// Bottom strand represented in its physical 5'->3' direction.
    pub complementary_sequence: String,
    pub columns: Vec<Column>,
    pub removed: Vec<Base>,
    pub history: Vec<RepairOperation>,
    pub left_reference: String,
    pub right_reference: String,
}

fn complement(base: char) -> char {
    match base {
        'A' => 'T',
        'C' => 'G',
        'G' => 'C',
        'T' => 'A',
        _ => unreachable!(),
    }
}
fn dna(sequence: &str) -> Result<String, String> {
    if !sequence.is_ascii()
        || sequence.is_empty()
        || !sequence
            .bytes()
            .all(|b| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T'))
    {
        return Err("duplex substrates and additions require nonempty resolved A/C/G/T".into());
    }
    Ok(sequence.to_ascii_uppercase())
}
fn bounds<O>(strand: &BTreeMap<i64, Base<O>>) -> Result<(i64, i64), String> {
    let first = *strand
        .first_key_value()
        .ok_or("a required strand is empty")?
        .0;
    let end = strand
        .last_key_value()
        .unwrap()
        .0
        .checked_add(1)
        .ok_or("coordinate overflow")?;
    if end.checked_sub(first) != Some(strand.len() as i64) {
        return Err("strand has an unresolved internal gap".into());
    }
    Ok((first, end))
}

impl End {
    fn new(sequence: &str, side: Side, next: &mut u64) -> Self {
        let mut top = BTreeMap::new();
        let mut bottom = BTreeMap::new();
        for (p, base) in sequence.chars().enumerate() {
            for (strand, value, dest) in [
                (Strand::Top, base, &mut top),
                (Strand::Bottom, complement(base), &mut bottom),
            ] {
                let id = *next;
                *next += 1;
                dest.insert(
                    p as i64,
                    Base {
                        id,
                        base: value,
                        kind: "original_strand".into(),
                        original: Some(Origin {
                            side: side.clone(),
                            strand,
                            reference_position: p,
                        }),
                        template_id: None,
                    },
                );
            }
        }
        Self {
            top,
            bottom,
            length: sequence.len(),
            original_start: 0,
            opened: false,
            barrier: None,
        }
    }
}

impl<O: Clone> End<O> {
    fn check_barrier(&self, side: &Side, strand: &Strand, p: i64) -> Result<(), String> {
        if let Some((top, bottom)) = self.barrier {
            let boundary = match strand {
                Strand::Top => top,
                Strand::Bottom => bottom,
            };
            if match side {
                Side::Left => p <= boundary,
                Side::Right => p >= boundary - 1,
            } {
                return Err("processing across a previous ligation seam is unsupported".into());
            }
        }
        Ok(())
    }

    fn open(
        &mut self,
        side: &Side,
        offset: usize,
        observer: &mut Option<&mut dyn OperationObserver>,
    ) -> Result<(), String> {
        if self.opened {
            return Err("a hairpin may be opened only once".into());
        }
        if self.barrier.is_none() && offset >= self.length {
            return Err(
                "opening must retain a paired stem and remain inside the original coding end"
                    .into(),
            );
        }
        let original_end = self
            .original_start
            .checked_add(self.length as i64)
            .ok_or("coordinate overflow")?;
        let (top_start, top_end) = bounds(&self.top)?;
        let (bottom_start, bottom_end) = bounds(&self.bottom)?;
        let source_len = match side {
            Side::Left => self.bottom.len(),
            Side::Right => self.top.len(),
        };
        if offset >= source_len {
            return Err("opening must retain an actual paired component stem".into());
        }
        if match side {
            Side::Left => top_end != original_end || bottom_end != original_end,
            Side::Right => top_start != self.original_start || bottom_start != self.original_start,
        } {
            return Err("the original hairpin turn is no longer present at this terminus".into());
        }
        if let Some((top_seam, bottom_seam)) = self.barrier {
            let cuts_seam = match side {
                Side::Left => original_end - offset as i64 <= bottom_seam,
                Side::Right => self.original_start + offset as i64 >= top_seam,
            };
            if cuts_seam {
                return Err("opening at or across a previous ligation seam is unsupported".into());
            }
        }
        // Opening moves the actual extant material, which can include synthesis
        // from a previous join. It never recreates a deleted reference base.
        let mut transfers = Vec::with_capacity(offset);
        match side {
            Side::Left => {
                for i in 0..offset {
                    let p = original_end - 1 - i as i64;
                    let q = original_end
                        .checked_add(i as i64)
                        .ok_or("coordinate overflow")?;
                    self.check_barrier(side, &Strand::Bottom, p)?;
                    if !self.bottom.contains_key(&p) || self.top.contains_key(&q) {
                        return Err("hairpin opening lacks its required current substrate".into());
                    }
                    transfers.push((p, q));
                }
            }
            Side::Right => {
                for i in 0..offset {
                    let p = self
                        .original_start
                        .checked_add(i as i64)
                        .ok_or("coordinate overflow")?;
                    let q = self
                        .original_start
                        .checked_sub(1 + i as i64)
                        .ok_or("coordinate overflow")?;
                    self.check_barrier(side, &Strand::Top, p)?;
                    if !self.top.contains_key(&p) || self.bottom.contains_key(&q) {
                        return Err("hairpin opening lacks its required current substrate".into());
                    }
                    transfers.push((p, q));
                }
            }
        }
        let paired_start = match side {
            Side::Left => top_start.max(bottom_start),
            Side::Right => (self.original_start + offset as i64).max(bottom_start),
        };
        let paired_end = match side {
            Side::Left => (original_end - offset as i64).min(top_end),
            Side::Right => top_end.min(bottom_end),
        };
        if paired_start >= paired_end {
            return Err("opening must retain an actual paired component stem".into());
        }
        for (p, q) in transfers {
            match side {
                Side::Left => {
                    let base = self.bottom.remove(&p).ok_or("missing hairpin substrate")?;
                    if let Some(o) = observer.as_mut() {
                        o.transfer(side, base.id, p, q)?;
                    }
                    self.top.insert(q, base);
                }
                Side::Right => {
                    let base = self.top.remove(&p).ok_or("missing hairpin substrate")?;
                    if let Some(o) = observer.as_mut() {
                        o.transfer(side, base.id, p, q)?;
                    }
                    self.bottom.insert(q, base);
                }
            }
        }
        self.opened = true;
        Ok(())
    }
    fn resect(
        &mut self,
        side: &Side,
        strand: &Strand,
        count: usize,
        removed: &mut Vec<Base<O>>,
        observer: &mut Option<&mut dyn OperationObserver>,
    ) -> Result<(), String> {
        if !self.opened {
            return Err("resection requires an opened coding end".into());
        }
        let target = match strand {
            Strand::Top => &self.top,
            Strand::Bottom => &self.bottom,
        };
        if count >= target.len() {
            return Err("resection must retain a strand primer/stem in this model".into());
        }
        let positions: Vec<_> = match side {
            Side::Left => target.keys().rev().take(count).copied().collect(),
            Side::Right => target.keys().take(count).copied().collect(),
        };
        for p in positions {
            self.check_barrier(side, strand, p)?;
        }
        let target = match strand {
            Strand::Top => &mut self.top,
            Strand::Bottom => &mut self.bottom,
        };
        for _ in 0..count {
            let item = match side {
                Side::Left => target.pop_last(),
                Side::Right => target.pop_first(),
            };
            let (at, base) = item.unwrap();
            if let Some(o) = observer.as_mut() {
                o.remove(side, strand, base.id, at)?;
            }
            removed.push(base);
        }
        Ok(())
    }
    fn add_n(
        &mut self,
        side: &Side,
        sequence: &str,
        next: &mut u64,
        ledger: &mut Vec<Base<O>>,
        observer: &mut Option<&mut dyn OperationObserver>,
    ) -> Result<(), String> {
        if !self.opened {
            return Err("TdT extension requires an opened 3-prime terminus".into());
        }
        let sequence = dna(sequence)?;
        // Input additions follow polymerization order, always 5'->3'.
        let (target, step) = match side {
            Side::Left => (&mut self.top, 1),
            Side::Right => (&mut self.bottom, -1),
        };
        let (start, end) = bounds(target)?;
        let mut position = if step == 1 { end } else { start - 1 };
        for base in sequence.chars() {
            let id = *next;
            *next += 1;
            let base = Base {
                id,
                base,
                kind: "n_addition".into(),
                original: None,
                template_id: None,
            };
            if let Some(o) = observer.as_mut() {
                let strand = match side {
                    Side::Left => Strand::Top,
                    Side::Right => Strand::Bottom,
                };
                o.create(side, &strand, id, position, None)?;
            }
            ledger.push(base.clone());
            target.insert(position, base);
            position = position.checked_add(step).ok_or("coordinate overflow")?;
        }
        Ok(())
    }

    fn translate(&mut self, shift: i64) -> Result<(), String> {
        for target in [&mut self.top, &mut self.bottom] {
            *target = std::mem::take(target)
                .into_iter()
                .map(|(p, b)| {
                    p.checked_add(shift)
                        .map(|p| (p, b))
                        .ok_or_else(|| "coordinate overflow".to_string())
                })
                .collect::<Result<_, _>>()?;
        }
        self.original_start = self
            .original_start
            .checked_add(shift)
            .ok_or("coordinate overflow")?;
        if let Some((top, bottom)) = self.barrier {
            self.barrier = Some((
                top.checked_add(shift).ok_or("coordinate overflow")?,
                bottom.checked_add(shift).ok_or("coordinate overflow")?,
            ));
        }
        Ok(())
    }
}

/// Resolve a finite authored repair history. Incomplete or incompatible products
/// are errors; only a fully paired, doubly ligated product is exportable.
pub fn repair(
    left: &str,
    right: &str,
    history: &[RepairOperation],
    max_operations: usize,
    max_bases: usize,
) -> Result<RepairProduct, String> {
    let left = dna(left)?;
    let right = dna(right)?;
    if history.len() > max_operations {
        return Err("repair operation budget exceeded".into());
    }
    let initial = left
        .len()
        .checked_add(right.len())
        .and_then(|n| n.checked_mul(2))
        .ok_or("base budget overflow")?;
    if initial > max_bases || max_bases > 1_000_000 {
        return Err("repair base budget exceeded (maximum 1000000)".into());
    }
    let mut next = 0;
    let mut l = End::new(&left, Side::Left, &mut next);
    let mut r = End::new(&right, Side::Right, &mut next);
    let mut removed = Vec::new();
    let mut ledger = Vec::new();
    execute_join(
        &mut l,
        &mut r,
        history,
        &mut next,
        &mut removed,
        &mut ledger,
        max_bases,
        None,
    )?;
    l.top.extend(r.top);
    l.bottom.extend(r.bottom);
    let columns = paired_columns(l.top, l.bottom)?;
    let sequence = columns.iter().map(|c| c.top.base).collect();
    let complementary_sequence = columns.iter().rev().map(|c| c.bottom.base).collect();
    Ok(RepairProduct {
        model: MODEL.into(),
        top_direction: "5_to_3".into(),
        bottom_direction: "3_to_5".into(),
        sequence,
        complementary_sequence,
        columns,
        removed,
        history: history.to_vec(),
        left_reference: left,
        right_reference: right,
    })
}

struct JoinedGeometry {
    right_translation: i64,
    top_seam: i64,
    bottom_seam: i64,
}

/// Execute one local join on existing physical components. Original coordinates
/// and nucleotide identities survive; only the right component is translated.
fn execute_join<O: Clone>(
    l: &mut End<O>,
    r: &mut End<O>,
    history: &[RepairOperation],
    next: &mut u64,
    removed: &mut Vec<Base<O>>,
    ledger: &mut Vec<Base<O>>,
    max_bases: usize,
    mut observer: Option<&mut dyn OperationObserver>,
) -> Result<JoinedGeometry, String> {
    let mut shift = None;
    let mut sealed = false;
    for (ordinal, operation) in history.iter().enumerate() {
        if sealed {
            return Err("no operation may follow complete ligation".into());
        }
        if let Some(o) = observer.as_mut() {
            o.start(ordinal, operation)?;
        }
        match operation {
            RepairOperation::Open { side, offset } => {
                if shift.is_some() {
                    return Err("opening must precede synapsis".into());
                }
                match side {
                    Side::Left => &mut *l,
                    Side::Right => &mut *r,
                }
                .open(side, *offset, &mut observer)?;
            }
            RepairOperation::Resect {
                side,
                strand,
                count,
            } => {
                if shift.is_some() {
                    return Err("this model requires resection before synapsis".into());
                }
                match side {
                    Side::Left => &mut *l,
                    Side::Right => &mut *r,
                }
                .resect(side, strand, *count, removed, &mut observer)?;
            }
            RepairOperation::AddN { side, sequence } => {
                if shift.is_some() {
                    return Err("this model requires TdT extension before synapsis".into());
                }
                if (*next as usize)
                    .checked_add(sequence.len())
                    .is_none_or(|n| n > max_bases)
                {
                    return Err("repair base budget exceeded".into());
                }
                match side {
                    Side::Left => &mut *l,
                    Side::Right => &mut *r,
                }
                .add_n(side, sequence, next, ledger, &mut observer)?;
            }
            RepairOperation::Synapse { right_start } => {
                if shift.is_some() || !l.opened || !r.opened {
                    return Err("synapsis requires two opened, unjoined ends".into());
                }
                if right_start.unsigned_abs() > max_bases as u64 {
                    return Err("synapsis coordinate exceeds local repair budget".into());
                }
                let translation = l
                    .original_start
                    .checked_add(*right_start)
                    .and_then(|p| p.checked_sub(r.original_start))
                    .ok_or("coordinate overflow")?;
                r.translate(translation)?;
                if let Some(o) = observer.as_mut() {
                    o.translate(
                        translation,
                        r.top
                            .values()
                            .chain(r.bottom.values())
                            .map(|b| b.id)
                            .collect(),
                    )?;
                }
                let (_, lt) = bounds(&l.top)?;
                let (rt, _) = bounds(&r.top)?;
                let (_, lb) = bounds(&l.bottom)?;
                let (rb, _) = bounds(&r.bottom)?;
                if lt > rt || lb > rb {
                    return Err(
                        "overlapping strands require flap processing outside this model's support"
                            .into(),
                    );
                }
                // Annealing is permitted only where paired bases agree. Blunt
                // juxtaposition can have zero overlap; incompatible pairs may
                // not silently become a consensus nucleotide.
                for (p, a) in l.top.iter().chain(r.top.iter()) {
                    if let Some(b) = l.bottom.get(p).or_else(|| r.bottom.get(p)) {
                        if complement(a.base) != b.base {
                            return Err("synapsed ends contain incompatible base pairs".into());
                        }
                    }
                }
                shift = Some(translation);
            }
            RepairOperation::Fill { strand } => {
                if shift.is_none() {
                    return Err("template-dependent synthesis requires paired-end synapsis".into());
                }
                let (begin, end) = match strand {
                    Strand::Top => (bounds(&l.top)?.1, bounds(&r.top)?.0),
                    Strand::Bottom => (bounds(&l.bottom)?.1, bounds(&r.bottom)?.0),
                };
                if begin == end {
                    return Err("there is no gap for this polymerase extension".into());
                }
                let (primer, opposing) = match strand {
                    Strand::Top => {
                        let p = begin - 1;
                        (l.top.get(&p), l.bottom.get(&p).or_else(|| r.bottom.get(&p)))
                    }
                    Strand::Bottom => {
                        let p = end;
                        (r.bottom.get(&p), l.top.get(&p).or_else(|| r.top.get(&p)))
                    }
                };
                if !primer
                    .zip(opposing)
                    .is_some_and(|(p, t)| complement(p.base) == t.base)
                {
                    return Err(
                        "template-dependent synthesis requires a paired 3-prime primer terminus"
                            .into(),
                    );
                }
                if (*next as usize)
                    .checked_add((end - begin) as usize)
                    .is_none_or(|n| n > max_bases)
                {
                    return Err("repair base budget exceeded".into());
                }
                let mut additions = Vec::new();
                for p in begin..end {
                    let template = match strand {
                        Strand::Top => l.bottom.get(&p).or_else(|| r.bottom.get(&p)),
                        Strand::Bottom => l.top.get(&p).or_else(|| r.top.get(&p)),
                    }.ok_or("a gap lacks a continuous opposing template; no fictitious fill-in is allowed")?;
                    let base = Base {
                        id: *next,
                        base: complement(template.base),
                        kind: "templated_synthesis".into(),
                        original: None,
                        template_id: Some(template.id),
                    };
                    if let Some(o) = observer.as_mut() {
                        let side = match strand {
                            Strand::Top => Side::Left,
                            Strand::Bottom => Side::Right,
                        };
                        o.create(&side, strand, base.id, p, Some(template.id))?;
                    }
                    ledger.push(base.clone());
                    additions.push((p, base));
                    *next += 1;
                }
                // Top grows rightward from the left 3' primer; bottom grows
                // leftward from the right 3' primer. A recessed 5' end never
                // serves as a polymerase primer.
                match strand {
                    Strand::Top => l.top.extend(additions),
                    Strand::Bottom => r.bottom.extend(additions),
                };
            }
            RepairOperation::Ligate {} => {
                if shift.is_none() {
                    return Err("ligation requires synapsed ends".into());
                }
                if bounds(&l.top)?.1 != bounds(&r.top)?.0
                    || bounds(&l.bottom)?.1 != bounds(&r.bottom)?.0
                {
                    return Err(
                        "ligation requires two resolved 3-OH/5-phosphate nicks without gaps".into(),
                    );
                }
                sealed = true;
                if let Some(o) = observer.as_mut() {
                    o.ligate(bounds(&l.top)?.1, bounds(&l.bottom)?.1)?;
                }
            }
        }
    }
    if !sealed {
        return Err("repair history has no complete double-strand ligation".into());
    }
    Ok(JoinedGeometry {
        right_translation: shift.ok_or("missing synapsis")?,
        top_seam: bounds(&l.top)?.1,
        bottom_seam: bounds(&l.bottom)?.1,
    })
}

fn paired_columns<O>(
    top: BTreeMap<i64, Base<O>>,
    mut bottom: BTreeMap<i64, Base<O>>,
) -> Result<Vec<Column<O>>, String> {
    if bounds(&top)? != bounds(&bottom)? {
        return Err("completed duplex has unequal strand boundaries".into());
    }
    let mut columns = Vec::with_capacity(top.len());
    for (position, top) in top {
        let bottom = bottom
            .remove(&position)
            .ok_or("completed duplex has an unpaired base")?;
        if complement(top.base) != bottom.base {
            return Err("completed duplex contains an unresolved mismatch".into());
        }
        columns.push(Column {
            position,
            top,
            bottom,
        });
    }
    Ok(columns)
}

#[cfg(test)]
mod tests {
    use super::*;
    fn run(l: &str, r: &str, history: &[RepairOperation]) -> Result<RepairProduct, String> {
        repair(l, r, history, 128, 4096)
    }
    fn opens(offset: usize) -> Vec<RepairOperation> {
        vec![
            RepairOperation::Open {
                side: Side::Left,
                offset,
            },
            RepairOperation::Open {
                side: Side::Right,
                offset,
            },
        ]
    }
    #[test]
    fn blunt_ligation_preserves_both_original_strands() {
        let mut h = opens(0);
        h.extend([
            RepairOperation::Synapse { right_start: 3 },
            RepairOperation::Ligate {},
        ]);
        let p = run("ACG", "TTA", &h).unwrap();
        assert_eq!(p.sequence, "ACGTTA");
        assert_eq!(p.complementary_sequence, "TAACGT");
        assert!(p
            .columns
            .iter()
            .all(|c| c.top.original.is_some() && c.bottom.original.is_some()));
    }
    #[test]
    fn compatible_hairpin_overhangs_join_with_two_contributing_origins() {
        let mut h = opens(1);
        h.extend([
            RepairOperation::Synapse { right_start: 3 },
            RepairOperation::Ligate {},
        ]);
        let p = run("AAG", "CAA", &h).unwrap();
        assert_eq!(p.sequence, "AAGCAA");
        assert_eq!(p.columns[2].top.original.as_ref().unwrap().side, Side::Left);
        assert_eq!(
            p.columns[2].bottom.original.as_ref().unwrap().side,
            Side::Right
        );
        assert_eq!(
            p.columns[3].top.original.as_ref().unwrap().strand,
            Strand::Bottom
        );
        // All twelve input bases remain, nicking itself creates no base.
        let ids: std::collections::BTreeSet<_> = p
            .columns
            .iter()
            .flat_map(|c| [c.top.id, c.bottom.id])
            .collect();
        assert_eq!(ids.len(), 12);
    }
    #[test]
    fn physical_resection_and_template_restoration_have_distinct_history() {
        let mut h = opens(0);
        h.extend([
            RepairOperation::Resect {
                side: Side::Left,
                strand: Strand::Top,
                count: 1,
            },
            RepairOperation::Synapse { right_start: 3 },
            RepairOperation::Fill {
                strand: Strand::Top,
            },
            RepairOperation::Ligate {},
        ]);
        let p = run("ACG", "TTA", &h).unwrap();
        assert_eq!(p.sequence, "ACGTTA");
        assert_eq!(p.removed.len(), 1);
        assert_eq!(p.columns[2].top.kind, "templated_synthesis");
        assert_eq!(p.columns[2].top.template_id, Some(p.columns[2].bottom.id));
    }
    #[test]
    fn opposed_n_additions_can_share_one_final_column() {
        let mut h = opens(0);
        h.extend([
            RepairOperation::AddN {
                side: Side::Left,
                sequence: "G".into(),
            },
            RepairOperation::AddN {
                side: Side::Right,
                sequence: "C".into(),
            },
            RepairOperation::Synapse { right_start: 4 },
            RepairOperation::Ligate {},
        ]);
        let p = run("AAA", "CCC", &h).unwrap();
        assert_eq!(p.sequence, "AAAGCCC");
        assert_eq!(p.columns[3].top.kind, "n_addition");
        assert_eq!(p.columns[3].bottom.kind, "n_addition");
    }
    #[test]
    fn chemistry_gaps_and_unresolved_histories_are_rejected() {
        assert!(run("AAA", "CCC", &opens(0))
            .unwrap_err()
            .contains("ligation"));
        let mut h = opens(1);
        h.extend([
            RepairOperation::Synapse { right_start: 3 },
            RepairOperation::Ligate {},
        ]);
        assert!(run("AAA", "CCC", &h).unwrap_err().contains("incompatible"));
        let mut h = opens(0);
        h.extend([
            RepairOperation::Synapse { right_start: 4 },
            RepairOperation::Fill {
                strand: Strand::Top,
            },
            RepairOperation::Ligate {},
        ]);
        assert!(run("AAA", "CCC", &h).unwrap_err().contains("template"));
        assert!(repair("AAA", "CCC", &opens(0), 1, 100).is_err());
        assert!(repair("AAA", "CCC", &opens(0), 100, 5).is_err());
    }
    #[test]
    fn full_operation_json_replay_is_exact() {
        let mut h = opens(0);
        h.extend([
            RepairOperation::Synapse { right_start: 3 },
            RepairOperation::Ligate {},
        ]);
        let decoded: Vec<RepairOperation> =
            serde_json::from_str(&serde_json::to_string(&h).unwrap()).unwrap();
        assert_eq!(run("ACG", "TTA", &h), run("ACG", "TTA", &decoded));
        assert!(serde_json::from_str::<RepairOperation>(
            r#"{"operation":"ligate","invented":true}"#
        )
        .is_err());
    }
}
