//! Exact existence of a productive completion with templated P additions.
//!
//! The codon automaton retains at most two bases. N regions traverse their
//! declared positive-weight first-order support; segment/P alternatives are
//! visited lazily. No complete N strings or cross-product of all junction
//! choices is materialized, and no random draws are made by the solver.

use super::ChoiceDomain;
use crate::address::{ChoiceAddress, PEnd};
use crate::assignment::TrimEnd;
use crate::ir::{Segment, Simulation};
use crate::pass::NpBaseSupport;
use crate::passes::PAdditionPass;
use crate::refdata::{AlleleId, RefDataConfig};
use crate::trace::ChoiceValue;

#[derive(Clone, Debug)]
pub(crate) struct SegmentCompletionDomain {
    pub segment: Segment,
    pub alleles: ChoiceDomain<AlleleId>,
    pub trim_5: ChoiceDomain<u32>,
    pub trim_3: ChoiceDomain<u32>,
}

#[derive(Clone, Debug)]
pub(crate) struct NpCompletionDomain {
    pub segment: Segment,
    pub lengths: ChoiceDomain<u32>,
    pub bases: NpBaseSupport,
}

#[derive(Clone, Debug)]
pub(crate) struct PProductiveFeasibility {
    pub segments: Vec<SegmentCompletionDomain>,
    pub p_lengths: Vec<(PEnd, ChoiceDomain<u32>)>,
    pub processing: Vec<(Segment, ChoiceDomain<crate::end_processing::EndOutcome>)>,
    pub np: Vec<NpCompletionDomain>,
    pub d_orientation: ChoiceDomain<bool>,
    pub check_frame: bool,
    pub check_stops: bool,
    pub check_coding: bool,
    pub check_v_anchor: bool,
    pub check_j_anchor: bool,
}

#[derive(Copy, Clone)]
struct Query<'a> {
    pass_index: usize,
    sim: &'a Simulation,
    refdata: &'a RefDataConfig,
    address: Option<ChoiceAddress>,
    candidate: Option<&'a ChoiceValue>,
    np_slot: Option<(Segment, u32, u32)>,
    allowed_alleles: &'a [(Segment, Vec<AlleleId>)],
}

impl PProductiveFeasibility {
    pub(crate) fn has_completion(&self, sim: &Simulation, refdata: &RefDataConfig) -> bool {
        self.solve(Query {
            pass_index: 0,
            sim,
            refdata,
            address: None,
            candidate: None,
            np_slot: None,
            allowed_alleles: &[],
        })
    }

    pub(crate) fn admits(
        &self,
        pass_index: usize,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
        address: &str,
        candidate: &ChoiceValue,
    ) -> bool {
        self.admits_with_allele_restrictions(pass_index, sim, refdata, address, candidate, &[])
    }

    pub(crate) fn admits_with_allele_restrictions(
        &self,
        pass_index: usize,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
        address: &str,
        candidate: &ChoiceValue,
        allowed_alleles: &[(Segment, Vec<AlleleId>)],
    ) -> bool {
        let address = ChoiceAddress::parse(address);
        // Later receptor revision can reuse original sample/trim addresses.
        // Its existing P bytes are historical, so only this compilation's
        // recombination owners may query the template-based completion model.
        let owner = match address {
            Some(ChoiceAddress::SampleAllele(s)) => self
                .segments
                .iter()
                .find(|d| d.segment == Segment::from(s))
                .and_then(|d| d.alleles.pass_index),
            Some(ChoiceAddress::Trim { segment, end }) => self
                .segments
                .iter()
                .find(|d| d.segment == Segment::from(segment))
                .and_then(|d| match end {
                    TrimEnd::Five => d.trim_5.pass_index,
                    TrimEnd::Three => d.trim_3.pass_index,
                }),
            Some(ChoiceAddress::PLength { end }) => self
                .p_lengths
                .iter()
                .find(|(e, _)| *e == end)
                .and_then(|(_, d)| d.pass_index),
            Some(ChoiceAddress::CodingEndOutcome(s)) => self
                .processing
                .iter()
                .find(|(segment, _)| *segment == Segment::from(s))
                .and_then(|(_, domain)| domain.pass_index),
            Some(ChoiceAddress::NpLength(s)) => self
                .np
                .iter()
                .find(|d| d.segment == Segment::from(s))
                .and_then(|d| d.lengths.pass_index),
            Some(ChoiceAddress::SampleAlleleDInverted) => self.d_orientation.pass_index,
            _ => None,
        };
        if owner != Some(pass_index) {
            return true;
        }
        refdata.is_some_and(|refdata| {
            self.solve(Query {
                pass_index,
                sim,
                refdata,
                address,
                candidate: Some(candidate),
                np_slot: None,
                allowed_alleles,
            })
        })
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn admits_np_base(
        &self,
        pass_index: usize,
        sim: &Simulation,
        refdata: &RefDataConfig,
        segment: Segment,
        index: u32,
        total_len: u32,
        base: u8,
    ) -> bool {
        if !self
            .np
            .iter()
            .any(|np| np.segment == segment && np.lengths.pass_index == Some(pass_index))
        {
            return true;
        }
        if index >= total_len {
            return false;
        }
        let candidate = ChoiceValue::Base(base);
        self.solve(Query {
            pass_index,
            sim,
            refdata,
            address: None,
            candidate: Some(&candidate),
            np_slot: Some((segment, index, total_len)),
            allowed_alleles: &[],
        })
    }

    fn solve(&self, query: Query<'_>) -> bool {
        Solver {
            spec: self,
            query,
            memo: [[None; 241]; 6],
        }
        .exists(0, 0)
    }
}

struct Solver<'a> {
    spec: &'a PProductiveFeasibility,
    query: Query<'a>,
    memo: [[Option<bool>; 241]; 6],
}

impl Solver<'_> {
    fn exists(&mut self, stage: usize, carry: u8) -> bool {
        if stage == 5 {
            return !self.spec.check_frame || carry == 0;
        }
        if let Some(answer) = self.memo[stage][carry as usize] {
            return answer;
        }
        let spec = self.spec;
        let query = self.query;
        let answer = match stage {
            0 => visit_segment(spec, query, Segment::V, |body, anchor| {
                let Some(anchor) = anchor else {
                    return false;
                };
                let Some(after_v) = advance(carry, &body[anchor..], spec.check_stops) else {
                    return false;
                };
                p_lengths(spec, query, PEnd::V3).iter().any(|length| {
                    if length as usize > body.len() {
                        return false;
                    }
                    let p = PAdditionPass::p_bytes_iter(PEnd::V3, body, length as usize);
                    advance_iter(after_v, p, spec.check_stops)
                        .is_some_and(|next| self.exists(1, next))
                })
            }),
            1 => self.np_exists(Segment::Np1, carry, 2),
            2 if spec.segments.iter().any(|d| d.segment == Segment::D) => {
                visit_segment(spec, query, Segment::D, |body, _| {
                    p_lengths(spec, query, PEnd::D5).iter().any(|five| {
                        if five as usize > body.len() {
                            return false;
                        }
                        let p5 = PAdditionPass::p_bytes_iter(PEnd::D5, body, five as usize);
                        let Some(after_d) = advance_iter(carry, p5, spec.check_stops)
                            .and_then(|s| advance(s, body, spec.check_stops))
                        else {
                            return false;
                        };
                        p_lengths(spec, query, PEnd::D3).iter().any(|three| {
                            if three as usize > body.len() {
                                return false;
                            }
                            let p3 = PAdditionPass::p_bytes_iter(PEnd::D3, body, three as usize);
                            advance_iter(after_d, p3, spec.check_stops)
                                .is_some_and(|s| self.exists(3, s))
                        })
                    })
                })
            }
            2 => self.exists(3, carry),
            3 => self.np_exists(Segment::Np2, carry, 4),
            4 => visit_segment(spec, query, Segment::J, |body, anchor| {
                let Some(anchor) = anchor else {
                    return false;
                };
                p_lengths(spec, query, PEnd::J5).iter().any(|length| {
                    if length as usize > body.len() {
                        return false;
                    }
                    let p = PAdditionPass::p_bytes_iter(PEnd::J5, body, length as usize);
                    advance_iter(carry, p, spec.check_stops)
                        .and_then(|s| advance(s, &body[..anchor + 3], spec.check_stops))
                        .is_some_and(|s| self.exists(5, s))
                })
            }),
            _ => unreachable!(),
        };
        self.memo[stage][carry as usize] = Some(answer);
        answer
    }

    fn np_exists(&mut self, segment: Segment, carry: u8, next_stage: usize) -> bool {
        let Some(domain) = self.spec.np.iter().find(|np| np.segment == segment) else {
            return self.exists(next_stage, carry);
        };
        let query = self.query;
        if domain.lengths.is_committed_before(query.pass_index) {
            let Some(region) = query
                .sim
                .sequence
                .regions
                .iter()
                .find(|r| r.segment == segment)
            else {
                return false;
            };
            let mut state = Some(carry);
            for position in region.start.index()..region.end.index() {
                state = state.and_then(|s| {
                    push(
                        s,
                        query.sim.pool.as_slice()[position as usize].base,
                        self.spec.check_stops,
                    )
                });
            }
            return state.is_some_and(|s| self.exists(next_stage, s));
        }
        let lengths = if let Some((active, _, total)) = query.np_slot {
            if active == segment {
                vec![total]
            } else {
                domain.lengths.values.clone()
            }
        } else if matches!(query.address, Some(ChoiceAddress::NpLength(s)) if Segment::from(s) == segment)
        {
            integer_candidate(query).into_iter().collect()
        } else {
            domain.lengths.values.clone()
        };
        let bases = domain.bases;
        let check_stops = self.spec.check_stops;
        // The bases already committed in this slot. Invariant across the
        // candidate lengths -- nothing in it reads `length` -- so it is
        // computed once rather than per length.
        let mut initial_carry = carry;
        let mut previous = 4usize;
        let mut consumed = 0;
        if let Some((active, index, _)) = query.np_slot {
            if active == segment {
                let Some(start) = query.sim.pool.len().checked_sub(index as usize) else {
                    return false;
                };
                for offset in 0..=index {
                    let base = if offset == index {
                        match query.candidate {
                            Some(ChoiceValue::Base(b)) => *b,
                            _ => return false,
                        }
                    } else {
                        query.sim.pool.as_slice()[start + offset as usize].base
                    };
                    let Some(state) = push(initial_carry, base, check_stops) else {
                        return false;
                    };
                    initial_carry = state;
                    previous = base_index(base);
                }
                consumed = index + 1;
            }
        }
        // ONE walk over the positions, not one per candidate length.
        //
        // The state set after position `p` does not depend on the total
        // length: the mask at each position is a function of the position and
        // of the previous base's index, both of which the state carries. So
        // the automaton for length L is a strict prefix of the one for L+1,
        // and walking each length from its start recomputed that prefix. With
        // the 51-value support (0..=50) every packaged locus declares, that is
        // 1275 position-steps per call where 50 suffice.
        //
        // A length is answered where the walk reaches it instead: at each
        // position the support admits, the carries reachable there are tested.
        // The function returns `bool`, so answering the lengths in position
        // order rather than domain order cannot change what it returns.
        // Sorted in place and walked with a cursor rather than indexed through
        // a membership array: `lengths` is already owned here, and the array
        // was one more allocation per call on top of the two this function
        // already makes -- measurable on the loci where a single candidate
        // length is the common case.
        let mut lengths = lengths;
        lengths.sort_unstable();
        let Some(&max_length) = lengths.last() else {
            return false;
        };
        if consumed > max_length {
            return false;
        }
        let mut cursor = lengths.partition_point(|&length| length < consumed);
        // Canonical N additions reach only a small fraction of the
        // IUPAC automaton. Keep sparse active states instead of scanning
        // all 1,205 carry/previous-base combinations per nucleotide.
        let mut states = vec![initial_carry as usize * 5 + previous];
        let mut next = Vec::with_capacity(16);
        let mut seen = [false; 1205];
        // Carries already passed to `exists` in this call. The answer is a
        // disjunction over (length, carry) pairs, so a carry that failed at one
        // admitted length fails at every other one -- testing it again only
        // repeats a memo lookup. Walking ascending rather than in domain order
        // meant answering every shorter length on the way, which is what this
        // takes back: each carry is tested at most once per call.
        let mut tested = [0u64; 4];
        // The zero-walk case, which the old code reached through an empty
        // `for position in consumed..length` range.
        if lengths[cursor] == consumed {
            cursor += 1;
            if self.any_new_carry_reaches(&states, next_stage, &mut tested) {
                return true;
            }
        }
        for position in consumed..max_length {
            next.clear();
            for &index in &states {
                let prior = index % 5;
                let mask = if position == 0 || prior == 4 {
                    bases.first
                } else {
                    bases.transitions[prior]
                };
                for bit in 0..4 {
                    if mask & (1 << bit) == 0 {
                        continue;
                    }
                    if let Some(state) = push((index / 5) as u8, b"ACGT"[bit], check_stops) {
                        let target = state as usize * 5 + bit;
                        if !seen[target] {
                            seen[target] = true;
                            next.push(target);
                        }
                    }
                }
            }
            for &index in &next {
                seen[index] = false;
            }
            std::mem::swap(&mut states, &mut next);
            // Terminal: states can only arise from this prefix, so none here
            // means none for every longer length either.
            if states.is_empty() {
                break;
            }
            let walked = position + 1;
            while cursor < lengths.len() && lengths[cursor] < walked {
                cursor += 1;
            }
            if cursor < lengths.len()
                && lengths[cursor] == walked
                && self.any_new_carry_reaches(&states, next_stage, &mut tested)
            {
                return true;
            }
        }
        false
    }

    /// Whether any carry reachable in `states` and not yet tried completes
    /// from `next_stage`. Records what it tried in `tested`.
    ///
    /// `states` are packed `carry * 5 + previous_base`; the previous base is
    /// spent once the N region ends, so the distinct carries are what matters
    /// and several states collapse onto one.
    ///
    /// A four-word bitset over the 241 carries rather than a 241-byte array
    /// scanned end to end: the active states are sparse -- a handful out of
    /// 1205 -- so a scan visiting every carry slot cost more than the carries
    /// it found. It also makes "which of these are new" one mask.
    fn any_new_carry_reaches(
        &mut self,
        states: &[usize],
        next_stage: usize,
        tested: &mut [u64; 4],
    ) -> bool {
        let mut carries = [0u64; 4];
        for &state in states {
            let carry = state / 5;
            carries[carry / 64] |= 1 << (carry % 64);
        }
        for (word, bits) in carries.into_iter().enumerate() {
            let mut fresh = bits & !tested[word];
            tested[word] |= bits;
            while fresh != 0 {
                let bit = fresh.trailing_zeros() as usize;
                fresh &= fresh - 1;
                if self.exists(next_stage, (word * 64 + bit) as u8) {
                    return true;
                }
            }
        }
        false
    }
}

fn integer_candidate(query: Query<'_>) -> Option<u32> {
    match query.candidate {
        Some(ChoiceValue::Int(value)) if (0..=u32::MAX as i64).contains(value) => {
            Some(*value as u32)
        }
        _ => None,
    }
}

/// The candidate P lengths for `end`: no candidates, exactly one, or the
/// declared domain.
///
/// Returned as a borrow rather than a `Vec` because the last case is the
/// common one and `domain.values.clone()` allocated and copied the whole
/// domain on every call, inside the solver's nested candidate loops. The
/// other three cases carry at most one value, which needs no allocation
/// either.
enum PLengths<'a> {
    None,
    One(u32),
    Domain(&'a [u32]),
}

impl PLengths<'_> {
    fn iter(&self) -> impl Iterator<Item = u32> + '_ {
        let (single, slice) = match self {
            PLengths::None => (None, [].as_slice()),
            PLengths::One(value) => (Some(*value), [].as_slice()),
            PLengths::Domain(values) => (None, *values),
        };
        single.into_iter().chain(slice.iter().copied())
    }
}

fn p_lengths<'a>(
    spec: &'a PProductiveFeasibility,
    query: Query<'a>,
    end: PEnd,
) -> PLengths<'a> {
    let Some((_, domain)) = spec.p_lengths.iter().find(|(e, _)| *e == end) else {
        return PLengths::One(0);
    };
    if query.address == Some(ChoiceAddress::PLength { end }) {
        return match integer_candidate(query) {
            Some(value) => PLengths::One(value),
            None => PLengths::None,
        };
    }
    if domain.is_committed_before(query.pass_index) {
        return match query.sim.assignments.get(end.source_segment()) {
            Some(inst) => PLengths::One(match end {
                PEnd::V3 | PEnd::D3 => inst.p_length_3,
                _ => inst.p_length_5,
            }),
            None => PLengths::None,
        };
    }
    PLengths::Domain(&domain.values)
}

fn trim_values(
    query: Query<'_>,
    segment: Segment,
    end: TrimEnd,
    domain: &ChoiceDomain<u32>,
) -> Vec<u32> {
    if matches!(query.address, Some(ChoiceAddress::Trim { segment: s, end: e }) if Segment::from(s) == segment && e == end)
    {
        return integer_candidate(query).into_iter().collect();
    }
    if domain.is_committed_before(query.pass_index) {
        return query
            .sim
            .assignments
            .get(segment)
            .map(|inst| {
                vec![u32::from(match end {
                    TrimEnd::Five => inst.trim_5,
                    TrimEnd::Three => inst.trim_3,
                })]
            })
            .unwrap_or_default();
    }
    domain.values.clone()
}

fn visit_segment(
    spec: &PProductiveFeasibility,
    query: Query<'_>,
    segment: Segment,
    mut visit: impl FnMut(&[u8], Option<usize>) -> bool,
) -> bool {
    let Some(domain) = spec.segments.iter().find(|d| d.segment == segment) else {
        return false;
    };
    let mut alleles = if matches!(query.address, Some(ChoiceAddress::SampleAllele(s)) if Segment::from(s) == segment)
    {
        match query.candidate {
            Some(ChoiceValue::AlleleId(id)) => vec![AlleleId::new(*id)],
            _ => return false,
        }
    } else if let Some(inst) = query.sim.assignments.get(segment) {
        vec![inst.allele_id]
    } else {
        domain.alleles.values.clone()
    };
    if let Some((_, allowed)) = query.allowed_alleles.iter().find(|(s, _)| *s == segment) {
        alleles.retain(|id| allowed.contains(id));
    }
    let five = trim_values(query, segment, TrimEnd::Five, &domain.trim_5);
    let three = trim_values(query, segment, TrimEnd::Three, &domain.trim_3);
    let orientations = if segment != Segment::D {
        vec![false]
    } else if query.address == Some(ChoiceAddress::SampleAlleleDInverted) {
        match query.candidate {
            Some(ChoiceValue::Bool(b)) => vec![*b],
            _ => return false,
        }
    } else if spec.d_orientation.is_committed_before(query.pass_index) {
        query
            .sim
            .assignments
            .get(segment)
            .map(|inst| vec![inst.orientation.is_reverse()])
            .unwrap_or_default()
    } else {
        spec.d_orientation.values.clone()
    };

    for id in alleles {
        let Some(allele) = query.refdata.get(segment, id) else {
            continue;
        };
        if ((segment == Segment::V && spec.check_v_anchor)
            || (segment == Segment::J && spec.check_j_anchor))
            && (!crate::productivity::reference_ready(query.refdata, segment, id)
                || !crate::productivity::reference_anchor_ready(allele))
        {
            continue;
        }
        if let Some((_, processing)) = spec.processing.iter().find(|(s, _)| *s == segment) {
            if !query.refdata.coding_ends_are_not_known_partial(segment, id) {
                continue;
            }
            let indices: Vec<usize> = if matches!(query.address,
                Some(ChoiceAddress::CodingEndOutcome(s)) if Segment::from(s) == segment)
            {
                integer_candidate(query)
                    .map(|index| index as usize)
                    .into_iter()
                    .collect()
            } else if processing.is_committed_before(query.pass_index) {
                query
                    .sim
                    .coding_end_choice(segment)
                    .map(|choice| choice.outcome_index)
                    .into_iter()
                    .collect()
            } else {
                (0..processing.values.len()).collect()
            };
            for index in indices {
                let Some(outcome) = processing.values.get(index) else {
                    continue;
                };
                for &reverse in &orientations {
                    let orientation = if reverse {
                        crate::assignment::SegmentOrientation::ReverseComplement
                    } else {
                        crate::assignment::SegmentOrientation::Forward
                    };
                    let Ok(processed) = crate::end_processing::materialize(
                        &allele.seq,
                        outcome.five,
                        outcome.three,
                        orientation,
                    ) else {
                        continue;
                    };
                    let t5 = processed.retained_start;
                    let end = processed.retained_end;
                    let anchor = if segment == Segment::D {
                        None
                    } else {
                        let Some(anchor) = allele.anchor.map(usize::from) else {
                            continue;
                        };
                        if t5 > anchor || anchor + 3 > end {
                            continue;
                        }
                        if spec.check_coding
                            && !crate::productivity::retained_stop_free(
                                query.refdata,
                                segment,
                                id,
                                t5 as u32,
                                (allele.seq.len() - end) as u32,
                            )
                        {
                            continue;
                        }
                        Some(processed.p5.len() + anchor - t5)
                    };
                    let mut body = processed.p5;
                    body.extend_from_slice(&processed.body);
                    body.extend_from_slice(&processed.p3);
                    if visit(&body, anchor) {
                        return true;
                    }
                }
            }
            continue;
        }
        for &t5 in &five {
            for &t3 in &three {
                let Some(end) = allele.seq.len().checked_sub(t3 as usize) else {
                    continue;
                };
                if t5 as usize > end {
                    continue;
                }
                let anchor = if segment == Segment::D {
                    None
                } else {
                    let Some(anchor) = allele.anchor.map(usize::from) else {
                        continue;
                    };
                    if (t5 as usize) > anchor || anchor + 3 > end {
                        continue;
                    }
                    if spec.check_coding
                        && !crate::productivity::retained_stop_free(
                            query.refdata,
                            segment,
                            id,
                            t5,
                            t3,
                        )
                    {
                        continue;
                    }
                    Some(anchor - t5 as usize)
                };
                for &reverse in &orientations {
                    let body: Vec<u8> = if reverse {
                        allele.seq[t5 as usize..end]
                            .iter()
                            .rev()
                            .map(|&b| crate::ir::complement_base(b))
                            .collect()
                    } else {
                        allele.seq[t5 as usize..end].to_vec()
                    };
                    if visit(&body, anchor) {
                        return true;
                    }
                }
            }
        }
    }
    false
}

fn base_index(base: u8) -> usize {
    match base.to_ascii_uppercase() {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => 4,
    }
}

/// 0 = empty, 1..15 = one IUPAC base, 16..240 = two bases.
fn push(carry: u8, base: u8, check_stops: bool) -> Option<u8> {
    let alphabet = b"ACGTRYSWKMBDHVN";
    let base = match base.to_ascii_uppercase() {
        b'U' => b'T',
        b => b,
    };
    let b = alphabet.iter().position(|&value| value == base)? as u8;
    if carry == 0 {
        Some(1 + b)
    } else if carry <= 15 {
        Some(16 + (carry - 1) * 15 + b)
    } else {
        let pair = carry - 16;
        if check_stops
            && crate::productivity::codon_stop([
                alphabet[(pair / 15) as usize],
                alphabet[(pair % 15) as usize],
                alphabet[b as usize],
            ]) != Some(false)
        {
            None
        } else {
            Some(0)
        }
    }
}

fn advance(carry: u8, bytes: &[u8], check_stops: bool) -> Option<u8> {
    advance_iter(carry, bytes.iter().copied(), check_stops)
}

/// Fold `push` over a stream of bases. The P-addition candidates come through
/// here as an iterator so they never have to be collected into a `Vec` first.
fn advance_iter(
    mut carry: u8,
    bases: impl Iterator<Item = u8>,
    check_stops: bool,
) -> Option<u8> {
    for base in bases {
        carry = push(carry, base, check_stops)?;
    }
    Some(carry)
}

#[cfg(test)]
mod tests;

#[cfg(test)]
mod processing_tests;
