//! Optional, deterministic certificates emitted by the physical repair executor.
//! These describe placement, not an inference from final sequence similarity.
use super::*;

pub const EVIDENCE_MODEL: &str = "duplex_operation_placement_v1";

#[derive(Clone, Copy, Debug, Serialize, PartialEq, Eq)]
pub struct OperationId {
    pub join_ordinal: usize,
    pub operation_ordinal: usize,
}

#[derive(Clone, Copy, Debug, Serialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum CodingEnd {
    Five,
    Three,
}

#[derive(Clone, Copy, Debug, Serialize, PartialEq, Eq)]
pub struct PhysicalTip {
    pub segment: Substrate,
    /// End on the oriented substrate, before any whole-read transformation.
    pub end: CodingEnd,
    /// Corresponding canonical source end; inversion swaps D's two ends.
    pub canonical_end: CodingEnd,
}

#[derive(Clone, Copy, Debug, Serialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum PlacementRole {
    OriginalBody,
    RestoredBody,
    HairpinTransfer,
    NAddition,
    CopiedNonbody,
}

#[derive(Clone, Debug, Serialize, PartialEq, Eq)]
pub struct PhysicalLocation {
    /// Original oriented zero of the component's coordinate frame.
    pub frame: Substrate,
    pub strand: Strand,
    pub position: i64,
}

#[derive(Clone, Debug, Serialize, PartialEq, Eq)]
pub struct Placement {
    /// Live body claim. Its strand is the canonical supporting origin strand,
    /// not necessarily the physical strand occupied by a restored nucleotide.
    /// Removal and hairpin transfer revoke it; historical events retain it.
    pub body: Option<CoupledOrigin>,
    pub role: PlacementRole,
    pub birth_operation: Option<OperationId>,
    pub last_transfer_operation: Option<OperationId>,
    pub removal_operation: Option<OperationId>,
    /// Absent after physical resection. Surviving final locations use V's axis.
    pub location: Option<PhysicalLocation>,
}

#[derive(Clone, Debug, Serialize, PartialEq, Eq)]
#[serde(tag = "event", rename_all = "snake_case")]
pub enum RepairEvidenceEvent {
    Initialized {
        base_id: u64,
        body: CoupledOrigin,
        location: PhysicalLocation,
    },
    OperationStarted {
        operation: OperationId,
        join: JoinName,
        action: RepairOperation,
    },
    Created {
        operation: OperationId,
        base_id: u64,
        location: PhysicalLocation,
        tip: PhysicalTip,
        role: PlacementRole,
        template_id: Option<u64>,
    },
    TemplateUsed {
        operation: OperationId,
        base_id: u64,
        template_id: u64,
        template_location: PhysicalLocation,
        body_at_use: Option<CoupledOrigin>,
    },
    Transferred {
        operation: OperationId,
        base_id: u64,
        from: PhysicalLocation,
        to: PhysicalLocation,
        tip: PhysicalTip,
        previous_body: Option<CoupledOrigin>,
    },
    Removed {
        operation: OperationId,
        base_id: u64,
        location: PhysicalLocation,
        tip: PhysicalTip,
        previous_body: Option<CoupledOrigin>,
    },
    Translated {
        operation: OperationId,
        base_ids: Vec<u64>,
        from_frame: Substrate,
        to_frame: Substrate,
        shift: i64,
    },
    Ligated {
        operation: OperationId,
        join: JoinName,
        frame: Substrate,
        top: i64,
        bottom: i64,
    },
}

#[derive(Clone, Debug, Serialize, PartialEq, Eq)]
pub struct CoupledEvidence {
    pub model: String,
    pub events: Vec<RepairEvidenceEvent>,
    pub placements: BTreeMap<u64, Placement>,
}

#[derive(Clone, Debug, Serialize, PartialEq, Eq)]
pub struct AnnotatedCoupledProduct {
    pub product: CoupledProduct,
    pub evidence: CoupledEvidence,
}

/// Shared-executor callbacks carry physical changes, never reconstructed
/// hypothetical histories. The unannotated executor supplies no observer.
pub(super) trait OperationObserver {
    fn start(&mut self, ordinal: usize, action: &RepairOperation) -> Result<(), String>;
    fn transfer(&mut self, side: &Side, id: u64, from: i64, to: i64) -> Result<(), String>;
    fn remove(&mut self, side: &Side, strand: &Strand, id: u64, at: i64) -> Result<(), String>;
    fn create(
        &mut self,
        side: &Side,
        strand: &Strand,
        id: u64,
        at: i64,
        template: Option<u64>,
    ) -> Result<(), String>;
    fn translate(&mut self, shift: i64, ids: Vec<u64>) -> Result<(), String>;
    fn ligate(&mut self, top: i64, bottom: i64) -> Result<(), String>;
}

pub(super) struct EvidenceObserver {
    evidence: CoupledEvidence,
    join: JoinName,
    operation: OperationId,
    left_frame: Substrate,
    right_frame: Substrate,
    d_inverted: bool,
}

impl EvidenceObserver {
    pub fn new(d_inverted: bool) -> Self {
        Self {
            evidence: CoupledEvidence {
                model: EVIDENCE_MODEL.into(),
                events: Vec::new(),
                placements: BTreeMap::new(),
            },
            join: JoinName::VD,
            operation: OperationId {
                join_ordinal: 0,
                operation_ordinal: 0,
            },
            left_frame: Substrate::V,
            right_frame: Substrate::D,
            d_inverted,
        }
    }

    pub fn finish(self) -> CoupledEvidence {
        self.evidence
    }

    pub fn initialize(
        &mut self,
        base: &Base<CoupledOrigin>,
        frame: Substrate,
        strand: Strand,
        position: i64,
    ) -> Result<(), String> {
        let body = base
            .original
            .clone()
            .ok_or("original nucleotide has no origin")?;
        let location = PhysicalLocation {
            frame,
            strand,
            position,
        };
        if self
            .evidence
            .placements
            .insert(
                base.id,
                Placement {
                    body: Some(body.clone()),
                    role: PlacementRole::OriginalBody,
                    birth_operation: None,
                    last_transfer_operation: None,
                    removal_operation: None,
                    location: Some(location.clone()),
                },
            )
            .is_some()
        {
            return Err("duplicate initialized nucleotide identity".into());
        }
        self.evidence.events.push(RepairEvidenceEvent::Initialized {
            base_id: base.id,
            body,
            location,
        });
        Ok(())
    }

    pub fn begin_join(
        &mut self,
        ordinal: usize,
        join: JoinName,
        left_frame: Substrate,
        right_frame: Substrate,
    ) {
        self.join = join;
        self.operation.join_ordinal = ordinal;
        self.left_frame = left_frame;
        self.right_frame = right_frame;
    }

    fn location(&self, side: &Side, strand: &Strand, position: i64) -> PhysicalLocation {
        PhysicalLocation {
            frame: match side {
                Side::Left => self.left_frame,
                Side::Right => self.right_frame,
            },
            strand: strand.clone(),
            position,
        }
    }

    fn tip(&self, side: &Side) -> PhysicalTip {
        let segment = match (self.join, side) {
            (JoinName::VD, Side::Left) => Substrate::V,
            (JoinName::VD, Side::Right) | (JoinName::DJ, Side::Left) => Substrate::D,
            (JoinName::DJ, Side::Right) => Substrate::J,
        };
        let end = match side {
            Side::Left => CodingEnd::Three,
            Side::Right => CodingEnd::Five,
        };
        let canonical_end = if segment == Substrate::D && self.d_inverted {
            match end {
                CodingEnd::Five => CodingEnd::Three,
                CodingEnd::Three => CodingEnd::Five,
            }
        } else {
            end
        };
        PhysicalTip {
            segment,
            end,
            canonical_end,
        }
    }

    fn require_location(&self, id: u64, expected: &PhysicalLocation) -> Result<(), String> {
        if self
            .evidence
            .placements
            .get(&id)
            .and_then(|p| p.location.as_ref())
            != Some(expected)
        {
            return Err("operation evidence disagrees with the current physical placement".into());
        }
        Ok(())
    }
}

impl OperationObserver for EvidenceObserver {
    fn start(&mut self, ordinal: usize, action: &RepairOperation) -> Result<(), String> {
        self.operation.operation_ordinal = ordinal;
        self.evidence
            .events
            .push(RepairEvidenceEvent::OperationStarted {
                operation: self.operation,
                join: self.join,
                action: action.clone(),
            });
        Ok(())
    }

    fn transfer(&mut self, side: &Side, id: u64, from: i64, to: i64) -> Result<(), String> {
        let (source, destination) = match side {
            Side::Left => (Strand::Bottom, Strand::Top),
            Side::Right => (Strand::Top, Strand::Bottom),
        };
        let from = self.location(side, &source, from);
        let to = self.location(side, &destination, to);
        self.require_location(id, &from)?;
        let tip = self.tip(side);
        let placement = self
            .evidence
            .placements
            .get_mut(&id)
            .ok_or("missing transfer origin")?;
        let previous_body = placement.body.take();
        placement.role = PlacementRole::HairpinTransfer;
        placement.last_transfer_operation = Some(self.operation);
        placement.location = Some(to.clone());
        self.evidence.events.push(RepairEvidenceEvent::Transferred {
            operation: self.operation,
            base_id: id,
            from,
            to,
            tip,
            previous_body,
        });
        Ok(())
    }

    fn remove(&mut self, side: &Side, strand: &Strand, id: u64, at: i64) -> Result<(), String> {
        let location = self.location(side, strand, at);
        self.require_location(id, &location)?;
        let tip = self.tip(side);
        let placement = self
            .evidence
            .placements
            .get_mut(&id)
            .ok_or("missing removed origin")?;
        let previous_body = placement.body.take();
        placement.location = None;
        placement.removal_operation = Some(self.operation);
        self.evidence.events.push(RepairEvidenceEvent::Removed {
            operation: self.operation,
            base_id: id,
            location,
            tip,
            previous_body,
        });
        Ok(())
    }

    fn create(
        &mut self,
        side: &Side,
        strand: &Strand,
        id: u64,
        at: i64,
        template: Option<u64>,
    ) -> Result<(), String> {
        let location = self.location(side, strand, at);
        let (body, role) = if let Some(template_id) = template {
            let template = self
                .evidence
                .placements
                .get(&template_id)
                .ok_or("missing template evidence")?;
            let template_location = template
                .location
                .as_ref()
                .ok_or("template was already removed")?;
            if template_location.position != at
                || template_location.frame != location.frame
                || template_location.strand == *strand
            {
                return Err(
                    "synthesis lacks a current opposite-strand template at its column".into(),
                );
            }
            let body = template.body.clone();
            self.evidence
                .events
                .push(RepairEvidenceEvent::TemplateUsed {
                    operation: self.operation,
                    base_id: id,
                    template_id,
                    template_location: template_location.clone(),
                    body_at_use: body.clone(),
                });
            let role = if body.is_some() {
                PlacementRole::RestoredBody
            } else {
                PlacementRole::CopiedNonbody
            };
            (body, role)
        } else {
            (None, PlacementRole::NAddition)
        };
        if self
            .evidence
            .placements
            .insert(
                id,
                Placement {
                    body,
                    role,
                    birth_operation: Some(self.operation),
                    last_transfer_operation: None,
                    removal_operation: None,
                    location: Some(location.clone()),
                },
            )
            .is_some()
        {
            return Err("duplicate synthesized nucleotide identity".into());
        }
        self.evidence.events.push(RepairEvidenceEvent::Created {
            operation: self.operation,
            base_id: id,
            location,
            tip: self.tip(side),
            role,
            template_id: template,
        });
        Ok(())
    }

    fn translate(&mut self, shift: i64, ids: Vec<u64>) -> Result<(), String> {
        for id in &ids {
            let location = self
                .evidence
                .placements
                .get_mut(id)
                .and_then(|p| p.location.as_mut())
                .ok_or("missing translated placement")?;
            if location.frame != self.right_frame {
                return Err("translated component has inconsistent coordinate frames".into());
            }
            location.position = location
                .position
                .checked_add(shift)
                .ok_or("coordinate overflow")?;
            location.frame = self.left_frame;
        }
        self.evidence.events.push(RepairEvidenceEvent::Translated {
            operation: self.operation,
            base_ids: ids,
            from_frame: self.right_frame,
            to_frame: self.left_frame,
            shift,
        });
        self.right_frame = self.left_frame;
        Ok(())
    }

    fn ligate(&mut self, top: i64, bottom: i64) -> Result<(), String> {
        self.evidence.events.push(RepairEvidenceEvent::Ligated {
            operation: self.operation,
            join: self.join,
            frame: self.left_frame,
            top,
            bottom,
        });
        Ok(())
    }
}

#[cfg(test)]
mod tests;
