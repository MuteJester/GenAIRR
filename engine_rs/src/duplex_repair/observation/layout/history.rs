//! Independently recover observation coverage from recorded terminal clipping.
//! An assignment cache is not evidence that reference material was observed.
use crate::duplex_repair::Substrate;
use crate::ir::{DuplexEditAction, DuplexQueryLink, DuplexRemovalCause, Segment, Simulation};

#[derive(Clone, Copy, Default)]
struct Coverage {
    start: Option<u32>,
    end: Option<u32>,
    lost_five: u32,
    lost_three: u32,
}

#[derive(Clone, Copy)]
struct Span {
    segment: Segment,
    start: usize,
    end: usize,
}

fn source_index(segment: Segment) -> Option<usize> {
    match segment {
        Segment::V => Some(0),
        Segment::D => Some(1),
        Segment::J => Some(2),
        _ => None,
    }
}

fn parse_segment(label: &str) -> Result<Segment, String> {
    match label {
        "V" => Ok(Segment::V),
        "D" => Ok(Segment::D),
        "J" => Ok(Segment::J),
        "Np1" => Ok(Segment::Np1),
        "Np2" => Ok(Segment::Np2),
        _ => Err("duplex edit records an unsupported segment".into()),
    }
}

pub(super) fn validate_history(sim: &Simulation) -> Result<(), String> {
    let p = sim
        .duplex_provenance
        .as_ref()
        .ok_or("simulation has no duplex provenance")?;
    let length = p.founder.body_claims.len();
    let mut spans = Vec::new();
    for (segment, substrate) in [
        (Segment::V, Substrate::V),
        (Segment::D, Substrate::D),
        (Segment::J, Substrate::J),
    ] {
        let columns: Vec<_> = p
            .founder
            .body_claims
            .iter()
            .enumerate()
            .filter(|(_, claims)| claims.iter().any(|claim| claim.segment == substrate))
            .map(|(q, _)| q)
            .collect();
        spans.push(Span {
            segment,
            start: *columns
                .first()
                .ok_or("duplex founder source has no body span")?,
            end: columns.last().unwrap() + 1,
        });
    }
    for (index, np) in [Segment::Np1, Segment::Np2].into_iter().enumerate() {
        let (start, end) = (spans[index].end, spans[index + 1].start);
        if start > end {
            return Err("duplex founder body spans overlap".into());
        }
        if start < end {
            spans.push(Span {
                segment: np,
                start,
                end,
            });
        }
    }
    spans.sort_by_key(|span| span.start);
    let mut roles = vec![None; length];
    for span in &spans {
        for q in span.start..span.end {
            roles[q] = Some(span.segment);
        }
    }
    let mut roles: Vec<_> = roles
        .into_iter()
        .collect::<Option<Vec<_>>>()
        .ok_or("duplex founder has an unclassified column")?;
    let mut links: Vec<_> = (0..length).map(DuplexQueryLink::FounderColumn).collect();
    let mut coverage = [Coverage::default(); 3];
    for edit in &p.edits {
        match &edit.action {
            DuplexEditAction::Substitution {
                at,
                segment,
                germline_pos,
                ..
            } => {
                let at = *at as usize;
                let segment = parse_segment(segment)?;
                if roles.get(at) != Some(&segment) {
                    return Err(
                        "duplex substitution segment differs from its recorded placement".into(),
                    );
                }
                let position = reference_position(p, &links[at], segment);
                if position != germline_pos.map(usize::from) {
                    return Err(
                        "duplex substitution reference cache differs from its placement".into(),
                    );
                }
            }
            DuplexEditAction::Insertion {
                at, segment, link, ..
            } => {
                let at = *at as usize;
                let segment = parse_segment(segment)?;
                if at > roles.len() {
                    return Err("duplex insertion lies outside the current query".into());
                }
                let owner = spans
                    .iter()
                    .position(|span| span.segment == segment && span.start <= at && at <= span.end)
                    .ok_or("duplex insertion has no declared structural owner")?;
                for (index, span) in spans.iter_mut().enumerate() {
                    if index == owner {
                        span.end += 1;
                    } else if at < span.start || (at == span.start && index > owner) {
                        span.start += 1;
                        span.end += 1;
                    }
                }
                roles.insert(at, segment);
                links.insert(at, link.clone());
            }
            DuplexEditAction::Deletion {
                at, segment, cause, ..
            } => {
                let at = *at as usize;
                let segment = parse_segment(segment)?;
                if roles.get(at) != Some(&segment) {
                    return Err(
                        "duplex deletion segment differs from its recorded placement".into(),
                    );
                }
                match cause {
                    DuplexRemovalCause::StructuralDeletion => {}
                    DuplexRemovalCause::Replacement => {
                        return Err("duplex replacement is unsupported".into())
                    }
                    DuplexRemovalCause::ObservationLoss => {
                        let five = match edit.source_pass.as_deref() {
                            Some("corrupt.end_loss.5") if at == 0 => true,
                            Some("corrupt.end_loss.3") if at + 1 == roles.len() => false,
                            _ => return Err(
                                "duplex observation loss lacks an audited terminal clipping event"
                                    .into(),
                            ),
                        };
                        if let Some(index) = source_index(segment) {
                            let c = &mut coverage[index];
                            if segment == Segment::V && five {
                                c.lost_five = c.lost_five.saturating_add(1);
                            } else if segment == Segment::J && !five {
                                c.lost_three = c.lost_three.saturating_add(1);
                            }
                            if let Some(g) = reference_position(p, &links[at], segment) {
                                let reverse =
                                    segment == Segment::D && p.founder.bindings.d_inverted;
                                if five != reverse {
                                    c.start = Some(c.start.unwrap_or(0).max(g as u32 + 1));
                                } else {
                                    c.end = Some(c.end.unwrap_or(u32::MAX).min(g as u32));
                                }
                            }
                            if spans.iter().any(|span| {
                                span.segment == segment && span.start == at && span.end == at + 1
                            }) {
                                c.start = Some(0);
                                c.end = Some(0);
                            }
                        }
                        let crossed = if five { at } else { at + 1 };
                        for span in &spans {
                            if span.start == crossed && span.end == crossed {
                                if let Some(index) = source_index(span.segment) {
                                    coverage[index].start = Some(0);
                                    coverage[index].end = Some(0);
                                }
                            }
                        }
                    }
                }
                for span in &mut spans {
                    if at < span.start {
                        span.start -= 1;
                        span.end -= 1;
                    } else if at < span.end {
                        span.end -= 1;
                    }
                }
                roles.remove(at);
                links.remove(at);
            }
            DuplexEditAction::Replacement { .. } => {
                return Err("duplex replacement is unsupported".into())
            }
            DuplexEditAction::ReverseComplement { .. } => {}
        }
    }
    if spans.len() != sim.sequence.regions.len()
        || spans
            .iter()
            .zip(&sim.sequence.regions)
            .any(|(expected, actual)| {
                expected.segment != actual.segment
                    || expected.start != actual.start.index() as usize
                    || expected.end != actual.end.index() as usize
            })
        || roles
            .iter()
            .copied()
            .ne(sim.pool.as_slice().iter().map(|n| n.segment))
    {
        return Err("duplex structural layout differs from its observation history".into());
    }
    for (index, segment) in [Segment::V, Segment::D, Segment::J].into_iter().enumerate() {
        let actual = sim
            .assignments
            .get(segment)
            .ok_or("duplex generating assignment is missing")?;
        let expected = coverage[index];
        if actual.observed_ref_start != expected.start
            || actual.observed_ref_end != expected.end
            || actual.observation_loss_5 != expected.lost_five
            || actual.observation_loss_3 != expected.lost_three
        {
            return Err(
                "duplex reference coverage cache differs from its terminal clipping history".into(),
            );
        }
    }
    Ok(())
}

fn reference_position(
    p: &crate::ir::DuplexProvenance,
    link: &DuplexQueryLink,
    segment: Segment,
) -> Option<usize> {
    let DuplexQueryLink::FounderColumn(column) = link else {
        return None;
    };
    p.founder.body_claims[*column]
        .iter()
        .find(|claim| {
            matches!(
                (segment, claim.segment),
                (Segment::V, Substrate::V)
                    | (Segment::D, Substrate::D)
                    | (Segment::J, Substrate::J)
            )
        })
        .map(|claim| claim.reference_position)
}
