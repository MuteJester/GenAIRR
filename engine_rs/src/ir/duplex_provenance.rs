//! Founder repair ancestry and independent descendant observation links.
//! A repair nucleotide ID names founder ancestry after replication, not an
//! assertion that descendants contain the same physical atoms.
use super::{Nucleotide, Segment};
use crate::address::PEnd;
use crate::duplex_repair::bindings::DuplexReferenceBindings;
use crate::duplex_repair::{AnnotatedCoupledProduct, Substrate};
use serde::Serialize;
use std::collections::BTreeMap;
use std::sync::Arc;

#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
#[serde(tag = "status", content = "value", rename_all = "snake_case")]
pub enum DuplexMetric<T> {
    Known(T),
    Unavailable { reason: String },
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
pub struct DuplexPMetrics {
    pub v3: DuplexMetric<u32>,
    pub d5: DuplexMetric<u32>,
    pub d3: DuplexMetric<u32>,
    pub j5: DuplexMetric<u32>,
}

impl Default for DuplexPMetrics {
    fn default() -> Self {
        let missing = || DuplexMetric::Unavailable {
            reason: "legacy_at_join_p_contract_not_determined".into(),
        };
        Self {
            v3: missing(),
            d5: missing(),
            d3: missing(),
            j5: missing(),
        }
    }
}

impl DuplexPMetrics {
    pub fn get(&self, end: PEnd) -> &DuplexMetric<u32> {
        match end {
            PEnd::V3 => &self.v3,
            PEnd::D5 => &self.d5,
            PEnd::D3 => &self.d3,
            PEnd::J5 => &self.j5,
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
pub struct DuplexBodyClaim {
    pub segment: Substrate,
    pub reference_position: usize,
    pub support_ids: Vec<u64>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DuplexFounder {
    pub annotated: Arc<AnnotatedCoupledProduct>,
    pub bindings: Arc<DuplexReferenceBindings>,
    /// Deduplicated column claims, independent of which physical strand
    /// provides the body evidence. The importer separately gates layout.
    pub body_claims: Vec<Vec<DuplexBodyClaim>>,
    pub p_metrics: DuplexPMetrics,
}

impl DuplexFounder {
    /// Construction preserves all claims; it does not substitute for the
    /// reference binding, replay integrity, or layout checks of the importer.
    pub fn new(
        annotated: Arc<AnnotatedCoupledProduct>,
        bindings: Arc<DuplexReferenceBindings>,
    ) -> Result<Self, String> {
        let mut body_claims = Vec::with_capacity(annotated.product.columns.len());
        for column in &annotated.product.columns {
            let mut claims = BTreeMap::<(Substrate, usize), Vec<u64>>::new();
            for base in [&column.top, &column.bottom] {
                let placement = annotated
                    .evidence
                    .placements
                    .get(&base.id)
                    .ok_or("founder nucleotide lacks placement evidence")?;
                if let Some(body) = &placement.body {
                    claims
                        .entry((body.segment, body.reference_position))
                        .or_default()
                        .push(base.id);
                }
            }
            body_claims.push(
                claims
                    .into_iter()
                    .map(|((segment, reference_position), mut support_ids)| {
                        support_ids.sort_unstable();
                        DuplexBodyClaim {
                            segment,
                            reference_position,
                            support_ids,
                        }
                    })
                    .collect(),
            );
        }
        Ok(Self {
            annotated,
            bindings,
            body_claims,
            p_metrics: DuplexPMetrics::default(),
        })
    }

    pub fn with_p_metrics(mut self, metrics: DuplexPMetrics) -> Self {
        self.p_metrics = metrics;
        self
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
#[serde(tag = "kind", content = "origin", rename_all = "snake_case")]
pub enum DuplexQueryLink {
    FounderColumn(usize),
    /// This is a branch-local edit reference, not a globally unique physical
    /// nucleotide ID. Pair it with its observation/lineage record to compare
    /// independently edited branches.
    Insertion {
        event_ordinal: u64,
        offset: u32,
        original_base: u8,
    },
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum DuplexRemovalCause {
    StructuralDeletion,
    ObservationLoss,
    Replacement,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
pub struct DuplexRemovedLink {
    pub link: DuplexQueryLink,
    pub query_position_at_removal: u32,
    pub cause: DuplexRemovalCause,
    pub event_ordinal: u64,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum DuplexEditAction {
    Substitution {
        at: u32,
        old_base: u8,
        new_base: u8,
        segment: String,
        germline_pos: Option<u16>,
        link: DuplexQueryLink,
    },
    Insertion {
        at: u32,
        base: u8,
        segment: String,
        link: DuplexQueryLink,
    },
    Deletion {
        at: u32,
        base: u8,
        segment: String,
        link: DuplexQueryLink,
        cause: DuplexRemovalCause,
    },
    Replacement {
        start: u32,
        end: u32,
        replacement_len: usize,
    },
    ReverseComplement {
        applied: bool,
    },
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
pub struct DuplexEdit {
    pub event_ordinal: u64,
    /// None denotes a low-level edit whose biological mechanism is unknown.
    pub source_pass: Option<String>,
    pub action: DuplexEditAction,
}

/// Sampled attempt counts are not applied substitutions. Every invocation,
/// including an empty query or zero count, contributes an ordered stage.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum DuplexAttemptStage {
    Pcr,
    Quality,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
pub struct DuplexAttemptSummary {
    pub stage_ordinal: u64,
    pub stage: DuplexAttemptStage,
    pub sampled_attempts: u32,
    pub input_pool_len: u32,
    /// Boundary before the first possible applied edit in this stage.
    pub edit_ordinal: u64,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DuplexProvenance {
    pub founder: Arc<DuplexFounder>,
    pub links: Vec<DuplexQueryLink>,
    pub removed: Vec<DuplexRemovedLink>,
    pub next_edit_ordinal: u64,
    pub edits: Vec<DuplexEdit>,
    pub attempt_stages: Vec<DuplexAttemptSummary>,
    /// Latest explicit output-orientation decision, matching the normal
    /// reverse-complement pass (the coding-oriented pool is never reversed).
    pub rev_comp: bool,
}

impl DuplexProvenance {
    pub fn new(founder: Arc<DuplexFounder>) -> Self {
        let links = (0..founder.annotated.product.columns.len())
            .map(DuplexQueryLink::FounderColumn)
            .collect();
        Self {
            founder,
            links,
            removed: Vec::new(),
            next_edit_ordinal: 0,
            edits: Vec::new(),
            attempt_stages: Vec::new(),
            rev_comp: false,
        }
    }

    /// Cumulative sampled attempts across the inherited history and all
    /// subsequent executions. Cloning preserves, rather than re-adds, history.
    pub fn sampled_attempts(&self, stage: DuplexAttemptStage) -> u64 {
        self.attempt_stages
            .iter()
            .filter(|s| s.stage == stage)
            .map(|s| u64::from(s.sampled_attempts))
            .sum()
    }

    pub(crate) fn with_attempt_stage(
        &self,
        stage: DuplexAttemptStage,
        sampled_attempts: u32,
    ) -> Self {
        let mut next = self.clone();
        next.attempt_stages.push(DuplexAttemptSummary {
            stage_ordinal: self.attempt_stages.len() as u64,
            stage,
            sampled_attempts,
            input_pool_len: self.links.len() as u32,
            edit_ordinal: self.next_edit_ordinal,
        });
        next
    }

    pub fn body_reference_position(
        &self,
        segment: Segment,
        query_position: usize,
    ) -> Option<usize> {
        let DuplexQueryLink::FounderColumn(column) = self.links.get(query_position)? else {
            return None;
        };
        let substrate = match segment {
            Segment::V => Substrate::V,
            Segment::D => Substrate::D,
            Segment::J => Substrate::J,
            _ => return None,
        };
        let mut candidates = self
            .founder
            .body_claims
            .get(*column)?
            .iter()
            .filter(|c| c.segment == substrate);
        let result = candidates.next()?.reference_position;
        candidates.next().is_none().then_some(result)
    }

    pub fn validate_mapping(&self, pool_len: usize) -> Result<(), String> {
        self.replay_observation(pool_len).map(|_| ())
    }

    /// Replay ancestry and bytes together: a changed current base is trustworthy
    /// only when the complete recorded history accounts for that change.
    fn replay_observation(&self, pool_len: usize) -> Result<Vec<u8>, String> {
        if self.links.len() != pool_len {
            return Err("duplex observation map length differs from the nucleotide pool".into());
        }
        let rebuilt = DuplexFounder::new(
            self.founder.annotated.clone(),
            self.founder.bindings.clone(),
        )?;
        if rebuilt.body_claims != self.founder.body_claims {
            return Err("duplex founder body claims differ from placement certificates".into());
        }
        if self.next_edit_ordinal != self.edits.len() as u64
            || self
                .edits
                .iter()
                .enumerate()
                .any(|(i, e)| e.event_ordinal != i as u64)
        {
            return Err("duplex edit ordinals are not a complete ordered history".into());
        }
        let mut founder_seen = std::collections::BTreeSet::new();
        let mut inserted_seen = std::collections::BTreeSet::new();
        for link in self
            .links
            .iter()
            .chain(self.removed.iter().map(|r| &r.link))
        {
            match link {
                DuplexQueryLink::FounderColumn(column) => {
                    if *column >= self.founder.body_claims.len() || !founder_seen.insert(*column) {
                        return Err(
                            "duplex observation ledger has an invalid or duplicated founder column"
                                .into(),
                        );
                    }
                }
                DuplexQueryLink::Insertion {
                    event_ordinal,
                    offset,
                    original_base,
                } => {
                    if !inserted_seen.insert((*event_ordinal, *offset)) {
                        return Err(
                            "duplex observation ledger has a duplicated insertion origin".into(),
                        );
                    }
                    let action = &self
                        .edits
                        .get(*event_ordinal as usize)
                        .ok_or("duplex insertion refers to a missing edit")?
                        .action;
                    match action {
                        DuplexEditAction::Insertion {
                            base, link: origin, ..
                        } if *offset == 0 && base == original_base && origin == link => {}
                        DuplexEditAction::Replacement {
                            replacement_len, ..
                        } if (*offset as usize) < *replacement_len => {}
                        _ => {
                            return Err(
                                "duplex insertion origin does not match its birth event".into()
                            )
                        }
                    }
                }
            }
        }
        if founder_seen.len() != self.founder.body_claims.len() {
            return Err("duplex ledger lost a founder column without recording removal".into());
        }
        let mut links: Vec<_> = (0..self.founder.body_claims.len())
            .map(DuplexQueryLink::FounderColumn)
            .collect();
        let mut query = self.founder.annotated.product.sequence.as_bytes().to_vec();
        if query.len() != links.len() {
            return Err("duplex founder query length differs from its column evidence".into());
        }
        let mut removed = Vec::new();
        let mut rev_comp = false;
        let mut query_lengths = vec![links.len()];
        for edit in &self.edits {
            let ordinal = edit.event_ordinal;
            match &edit.action {
                DuplexEditAction::Substitution {
                    at,
                    link,
                    old_base,
                    new_base,
                    ..
                } => {
                    if links.get(*at as usize) != Some(link) {
                        return Err(
                            "duplex substitution does not address its recorded ancestry".into()
                        );
                    }
                    if query.get(*at as usize) != Some(old_base) {
                        return Err(
                            "duplex substitution source byte differs from its prior observation"
                                .into(),
                        );
                    }
                    query[*at as usize] = *new_base;
                }
                DuplexEditAction::Insertion { at, link, base, .. } => {
                    if *at as usize > links.len() {
                        return Err("duplex insertion outside query".into());
                    }
                    links.insert(*at as usize, link.clone());
                    query.insert(*at as usize, *base);
                }
                DuplexEditAction::Deletion {
                    at,
                    link,
                    cause,
                    base,
                    ..
                } => {
                    if links.get(*at as usize) != Some(link) {
                        return Err("duplex deletion ancestry mismatch".into());
                    }
                    if query.get(*at as usize) != Some(base) {
                        return Err("duplex deleted byte differs from its prior observation".into());
                    }
                    links.remove(*at as usize);
                    query.remove(*at as usize);
                    removed.push(DuplexRemovedLink {
                        link: link.clone(),
                        query_position_at_removal: *at,
                        cause: *cause,
                        event_ordinal: ordinal,
                    });
                }
                DuplexEditAction::Replacement {
                    start,
                    end,
                    replacement_len,
                } => {
                    if start > end || *end as usize > links.len() {
                        return Err("duplex replacement outside query".into());
                    }
                    for (offset, link) in links[*start as usize..*end as usize]
                        .iter()
                        .cloned()
                        .enumerate()
                    {
                        removed.push(DuplexRemovedLink {
                            link,
                            query_position_at_removal: *start + offset as u32,
                            cause: DuplexRemovalCause::Replacement,
                            event_ordinal: ordinal,
                        });
                    }
                    let mut replacement = Vec::new();
                    for offset in 0..*replacement_len {
                        let link = self
                            .links
                            .iter()
                            .chain(self.removed.iter().map(|r| &r.link))
                            .find(|l| {
                                matches!(l, DuplexQueryLink::Insertion {event_ordinal,offset:o,..}
                                if *event_ordinal == ordinal && *o as usize == offset)
                            })
                            .ok_or("duplex replacement lost an insertion origin")?;
                        replacement.push(link.clone());
                    }
                    let replacement_bytes: Vec<_> = replacement
                        .iter()
                        .map(|link| match link {
                            DuplexQueryLink::Insertion { original_base, .. } => *original_base,
                            DuplexQueryLink::FounderColumn(_) => {
                                unreachable!("replacement origins were checked above")
                            }
                        })
                        .collect();
                    links.splice(*start as usize..*end as usize, replacement);
                    query.splice(*start as usize..*end as usize, replacement_bytes);
                }
                DuplexEditAction::ReverseComplement { applied } => {
                    rev_comp = *applied;
                }
            }
            query_lengths.push(links.len());
        }
        let mut previous_boundary = 0;
        for (ordinal, stage) in self.attempt_stages.iter().enumerate() {
            if stage.stage_ordinal != ordinal as u64
                || stage.edit_ordinal < previous_boundary
                || query_lengths.get(stage.edit_ordinal as usize).copied()
                    != Some(stage.input_pool_len as usize)
            {
                return Err(
                    "duplex attempt summary does not match its ordered input observation".into(),
                );
            }
            previous_boundary = stage.edit_ordinal;
        }
        if links != self.links || removed != self.removed || rev_comp != self.rev_comp {
            return Err("duplex current observation differs from its edit history".into());
        }
        Ok(query)
    }

    /// Validate the derived single-owner nucleotide cache without mistaking
    /// chemical ancestry for reference placement. Current bytes must be exactly
    /// explained by recorded edits, including case changes and ambiguity bases.
    pub fn validate_pool(&self, pool: &[Nucleotide]) -> Result<(), String> {
        let replayed = self.replay_observation(pool.len())?;
        if replayed.into_iter().ne(pool.iter().map(|n| n.base)) {
            return Err(
                "duplex nucleotide bytes differ from the recorded observation history".into(),
            );
        }
        for (n, link) in pool.iter().zip(&self.links) {
            match link {
                DuplexQueryLink::FounderColumn(column) => {
                    let claims = &self.founder.body_claims[*column];
                    if claims.len() == 1 {
                        let claim = &claims[0];
                        let segment = match claim.segment {
                            Substrate::V => Segment::V,
                            Substrate::D => Segment::D,
                            Substrate::J => Segment::J,
                        };
                        let binding = self.founder.bindings.for_substrate(claim.segment);
                        let mut germline = *binding
                            .sequence
                            .as_bytes()
                            .get(claim.reference_position)
                            .ok_or("duplex body placement is outside its reference")?;
                        if segment == Segment::D && self.founder.bindings.d_inverted {
                            germline = super::complement_base(germline);
                        }
                        if n.segment != segment
                            || n.germline_pos.get().map(usize::from)
                                != Some(claim.reference_position)
                            || n.germline.to_ascii_uppercase() != germline.to_ascii_uppercase()
                        {
                            return Err(
                                "duplex nucleotide reference cache differs from body placement"
                                    .into(),
                            );
                        }
                    } else if !claims.is_empty() {
                        return Err(
                            "duplex execution cache cannot represent overlapping body placements"
                                .into(),
                        );
                    } else if n.germline_pos.get().is_some() {
                        return Err(
                            "duplex nonbody column has an invented germline coordinate".into()
                        );
                    }
                }
                DuplexQueryLink::Insertion { original_base, .. } => {
                    if n.germline_pos.get().is_some() || n.germline != *original_base {
                        return Err(
                            "duplex insertion nucleotide has an invented reference origin".into(),
                        );
                    }
                }
            }
        }
        Ok(())
    }

    pub fn validate_initial_pool(&self, pool: &[Nucleotide]) -> Result<(), String> {
        self.validate_pool(pool)?;
        if !self.attempt_stages.is_empty()
            || !self.removed.is_empty()
            || !self.edits.is_empty()
            || self.next_edit_ordinal != 0
            || self.rev_comp
            || self
                .links
                .iter()
                .enumerate()
                .any(|(p, link)| *link != DuplexQueryLink::FounderColumn(p))
            || pool.iter().map(|n| n.base.to_ascii_uppercase()).ne(self
                .founder
                .annotated
                .product
                .sequence
                .bytes())
        {
            return Err("duplex installation requires the complete unedited founder query".into());
        }
        Ok(())
    }

    pub fn with_inserted(&self, at: u32, nucleotide: Nucleotide) -> Self {
        let mut next =
            self.with_replacement(at, at, &[nucleotide], DuplexRemovalCause::Replacement);
        next.edits.last_mut().unwrap().action = DuplexEditAction::Insertion {
            at,
            base: nucleotide.base,
            segment: format!("{:?}", nucleotide.segment),
            link: next.links[at as usize].clone(),
        };
        next
    }

    pub fn with_removed(&self, at: u32, nucleotide: Nucleotide, cause: DuplexRemovalCause) -> Self {
        let mut next = self.with_replacement(at, at + 1, &[], cause);
        next.edits.last_mut().unwrap().action = DuplexEditAction::Deletion {
            at,
            base: nucleotide.base,
            segment: format!("{:?}", nucleotide.segment),
            link: self.links[at as usize].clone(),
            cause,
        };
        next
    }

    pub fn with_substitution(&self, at: u32, old: Nucleotide, new_base: u8) -> Self {
        let mut next = self.clone();
        let ordinal = next.next_edit_ordinal;
        next.next_edit_ordinal = ordinal
            .checked_add(1)
            .expect("duplex edit ordinal overflow");
        next.edits.push(DuplexEdit {
            event_ordinal: ordinal,
            source_pass: None,
            action: DuplexEditAction::Substitution {
                at,
                old_base: old.base,
                new_base,
                segment: format!("{:?}", old.segment),
                germline_pos: old.germline_pos.get(),
                link: self.links[at as usize].clone(),
            },
        });
        next
    }

    pub fn with_reverse_complement(&self, applied: bool) -> Self {
        let mut next = self.clone();
        let ordinal = next.next_edit_ordinal;
        next.next_edit_ordinal = ordinal
            .checked_add(1)
            .expect("duplex edit ordinal overflow");
        next.rev_comp = applied;
        next.edits.push(DuplexEdit {
            event_ordinal: ordinal,
            source_pass: None,
            action: DuplexEditAction::ReverseComplement { applied },
        });
        next
    }

    pub(crate) fn set_last_edit_source(&mut self, source: Option<&str>) {
        if let Some(edit) = self.edits.last_mut() {
            edit.source_pass = source.map(str::to_owned);
        }
    }

    pub fn with_replacement(
        &self,
        start: u32,
        end: u32,
        replacement: &[Nucleotide],
        cause: DuplexRemovalCause,
    ) -> Self {
        assert!(
            start <= end && end as usize <= self.links.len(),
            "duplex edit outside observation"
        );
        if start == end && replacement.is_empty() {
            return self.clone();
        }
        let mut next = self.clone();
        let ordinal = next.next_edit_ordinal;
        next.next_edit_ordinal = ordinal
            .checked_add(1)
            .expect("duplex edit ordinal overflow");
        for (offset, link) in next.links[start as usize..end as usize]
            .iter()
            .cloned()
            .enumerate()
        {
            next.removed.push(DuplexRemovedLink {
                link,
                query_position_at_removal: start + offset as u32,
                cause,
                event_ordinal: ordinal,
            });
        }
        let insertions =
            replacement
                .iter()
                .enumerate()
                .map(|(offset, n)| DuplexQueryLink::Insertion {
                    event_ordinal: ordinal,
                    offset: offset as u32,
                    original_base: n.base,
                });
        next.links.splice(start as usize..end as usize, insertions);
        next.edits.push(DuplexEdit {
            event_ordinal: ordinal,
            source_pass: None,
            action: DuplexEditAction::Replacement {
                start,
                end,
                replacement_len: replacement.len(),
            },
        });
        next
    }
}

#[cfg(test)]
mod tests;
