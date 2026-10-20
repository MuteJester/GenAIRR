//! Two authored local joins acting on one shared, physically persistent D.
use super::*;

pub const COUPLED_MODEL: &str = "experimental_coupled_duplex_repair_v1";

#[derive(Clone, Copy, Debug, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "lowercase")]
pub enum JoinName {
    VD,
    DJ,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct JoinScript {
    pub join: JoinName,
    pub operations: Vec<RepairOperation>,
}

#[derive(Clone, Copy, Debug, Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord)]
#[serde(rename_all = "lowercase")]
pub enum Substrate {
    V,
    D,
    J,
}

#[derive(Clone, Debug, Serialize, PartialEq, Eq)]
pub struct CoupledOrigin {
    pub segment: Substrate,
    /// Original canonical strand, independent of current physical placement.
    pub strand: Strand,
    pub reference_position: usize,
}

#[derive(Clone, Debug, Serialize, PartialEq, Eq)]
pub struct LigationSeam {
    pub join: JoinName,
    /// Coordinates of the right nucleotide adjacent to each ligated bond.
    pub top: i64,
    pub bottom: i64,
}

#[derive(Clone, Debug, Serialize, PartialEq, Eq)]
pub struct CoupledProduct {
    pub model: String,
    pub top_direction: String,
    pub bottom_direction: String,
    pub sequence: String,
    pub complementary_sequence: String,
    pub columns: Vec<Column<CoupledOrigin>>,
    pub removed: Vec<Base<CoupledOrigin>>,
    /// Every original and synthesized physical nucleotide, in allocation order.
    /// Template links remain resolvable even if their template is later removed.
    pub ledger: Vec<Base<CoupledOrigin>>,
    pub joins: Vec<JoinScript>,
    pub v_reference: String,
    pub d_reference: String,
    pub j_reference: String,
    pub d_inverted: bool,
    /// Original oriented coordinate-zero anchors in the final V-rooted axis.
    pub origin_positions: BTreeMap<Substrate, i64>,
    pub ligation_seams: Vec<LigationSeam>,
}

struct Component {
    end: End<CoupledOrigin>,
    origins: BTreeMap<Substrate, i64>,
    seams: Vec<LigationSeam>,
}

impl Component {
    fn original(
        sequence: &str,
        segment: Substrate,
        inverted: bool,
        next: &mut u64,
        ledger: &mut Vec<Base<CoupledOrigin>>,
    ) -> Self {
        let mut top = BTreeMap::new();
        let mut bottom = BTreeMap::new();
        for (p, value) in sequence.chars().enumerate() {
            for (strand, base) in [(Strand::Top, value), (Strand::Bottom, complement(value))] {
                let nucleotide = Base {
                    id: *next,
                    base,
                    kind: "original_strand".into(),
                    original: Some(CoupledOrigin {
                        segment,
                        strand: strand.clone(),
                        reference_position: p,
                    }),
                    template_id: None,
                };
                *next += 1;
                ledger.push(nucleotide.clone());
                let physical_top = (strand == Strand::Top) != inverted;
                let position = if inverted { sequence.len() - 1 - p } else { p } as i64;
                if physical_top {
                    top.insert(position, nucleotide);
                } else {
                    bottom.insert(position, nucleotide);
                }
            }
        }
        Self {
            end: End {
                top,
                bottom,
                length: sequence.len(),
                original_start: 0,
                opened: false,
                barrier: None,
            },
            origins: BTreeMap::from([(segment, 0)]),
            seams: Vec::new(),
        }
    }

    fn expose(&mut self, segment: Substrate, original_length: usize) -> Result<(), String> {
        self.end.original_start = *self
            .origins
            .get(&segment)
            .ok_or("the requested substrate is absent from this component")?;
        self.end.length = original_length;
        self.end.opened = false;
        self.end.barrier = self.seams.first().map(|s| (s.top, s.bottom));
        Ok(())
    }

    fn translate_ancestry(&mut self, shift: i64) -> Result<(), String> {
        // execute_join already translated this component's physical strands.
        for p in self.origins.values_mut() {
            *p = p.checked_add(shift).ok_or("coordinate overflow")?;
        }
        for seam in &mut self.seams {
            seam.top = seam.top.checked_add(shift).ok_or("coordinate overflow")?;
            seam.bottom = seam
                .bottom
                .checked_add(shift)
                .ok_or("coordinate overflow")?;
        }
        Ok(())
    }
}

/// Resolve VD and DJ in either authored order using exactly one physical D.
///
/// Each local `right_start` is relative to its named left substrate's original
/// oriented coordinate zero: D relative to V for VD, J relative to D for DJ.
/// The previous ligation remains intact. Cutting/resecting at or across it is
/// outside this bounded model; surviving stems may contain previous synthesis.
pub fn repair_coupled(
    v: &str,
    d: &str,
    j: &str,
    d_inverted: bool,
    joins: &[JoinScript],
    max_operations: usize,
    max_bases: usize,
) -> Result<CoupledProduct, String> {
    let v = dna(v)?;
    let d = dna(d)?;
    let j = dna(j)?;
    if joins.len() != 2 || joins[0].join == joins[1].join {
        return Err(
            "coupled repair requires exactly one VD and one DJ join in either order".into(),
        );
    }
    let operation_count = joins.iter().try_fold(0usize, |n, script| {
        n.checked_add(script.operations.len())
            .ok_or("repair operation budget overflow")
    })?;
    if operation_count > max_operations {
        return Err("repair operation budget exceeded across successive joins".into());
    }
    let initial = v
        .len()
        .checked_add(d.len())
        .and_then(|n| n.checked_add(j.len()))
        .and_then(|n| n.checked_mul(2))
        .ok_or("base budget overflow")?;
    if initial > max_bases || max_bases > 1_000_000 {
        return Err("repair base budget exceeded (maximum 1000000)".into());
    }
    let mut ledger = Vec::with_capacity(initial);
    let mut next = 0;
    // IDs are fixed by canonical substrate/strand coordinates, not join order
    // or orientation. A second join never reinitializes a physical nucleotide.
    let mut components = vec![
        Some(Component::original(
            &v,
            Substrate::V,
            false,
            &mut next,
            &mut ledger,
        )),
        Some(Component::original(
            &d,
            Substrate::D,
            d_inverted,
            &mut next,
            &mut ledger,
        )),
        Some(Component::original(
            &j,
            Substrate::J,
            false,
            &mut next,
            &mut ledger,
        )),
    ];
    let lengths = BTreeMap::from([
        (Substrate::V, v.len()),
        (Substrate::D, d.len()),
        (Substrate::J, j.len()),
    ]);
    let mut removed = Vec::new();
    for script in joins {
        let (left, right) = match script.join {
            JoinName::VD => (Substrate::V, Substrate::D),
            JoinName::DJ => (Substrate::D, Substrate::J),
        };
        let find = |segment| {
            components
                .iter()
                .position(|c| c.as_ref().is_some_and(|c| c.origins.contains_key(&segment)))
                .ok_or("join substrate is absent")
        };
        let li = find(left)?;
        let ri = find(right)?;
        if li == ri {
            return Err("a local join requires two distinct unjoined components".into());
        }
        let mut l = components[li].take().ok_or("missing left component")?;
        let mut r = components[ri].take().ok_or("missing right component")?;
        l.expose(left, lengths[&left])?;
        r.expose(right, lengths[&right])?;
        let geometry = execute_join(
            &mut l.end,
            &mut r.end,
            &script.operations,
            &mut next,
            &mut removed,
            &mut ledger,
            max_bases,
        )?;
        r.translate_ancestry(geometry.right_translation)?;
        l.end.top.extend(r.end.top);
        l.end.bottom.extend(r.end.bottom);
        // Each completed local join must already be a real paired component.
        if bounds(&l.end.top)? != bounds(&l.end.bottom)? {
            return Err("completed duplex has unequal strand boundaries".into());
        }
        for (p, top) in &l.end.top {
            if !l
                .end
                .bottom
                .get(p)
                .is_some_and(|b| complement(top.base) == b.base)
            {
                return Err("completed duplex contains an unresolved mismatch or gap".into());
            }
        }
        l.origins.extend(r.origins);
        l.seams.extend(r.seams);
        l.seams.push(LigationSeam {
            join: script.join,
            top: geometry.top_seam,
            bottom: geometry.bottom_seam,
        });
        components[li] = Some(l);
    }
    let mut remaining = components.into_iter().flatten();
    let component = remaining
        .next()
        .ok_or("completed repair has no component")?;
    if remaining.next().is_some() || component.origins.len() != 3 {
        return Err("completed coupled repair must contain one connected VDJ component".into());
    }
    let columns = paired_columns(component.end.top, component.end.bottom)?;
    Ok(CoupledProduct {
        model: COUPLED_MODEL.into(),
        top_direction: "5_to_3".into(),
        bottom_direction: "3_to_5".into(),
        sequence: columns.iter().map(|c| c.top.base).collect(),
        complementary_sequence: columns.iter().rev().map(|c| c.bottom.base).collect(),
        columns,
        removed,
        ledger,
        joins: joins.to_vec(),
        v_reference: v,
        d_reference: d,
        j_reference: j,
        d_inverted,
        origin_positions: component.origins,
        ligation_seams: component.seams,
    })
}

#[cfg(test)]
mod tests;
