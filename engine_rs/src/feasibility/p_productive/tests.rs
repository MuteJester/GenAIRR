//! Brute-force oracle: enumerate complete N strings, physically concatenate
//! every retained segment and palindrome, and inspect literal stop codons.
//! The oracle does not use the solver's codon automaton, P derivation,
//! retained-reference checks, or runtime productivity contracts.

use super::*;
use crate::assignment::{AlleleInstance, SegmentOrientation};
use crate::ir::{flag, NucHandle, Nucleotide, Region};
use crate::refdata::{Allele, ChainType};
use std::collections::HashSet;

// Choice order follows lowering, with deterministic assembly steps omitted.
const ADDRESSES: [&str; 14] = [
    "sample_allele.v",
    "sample_allele.d",
    "sample_allele.j",
    "trim.v_3",
    "trim.d_5",
    "trim.d_3",
    "trim.j_5",
    "p.v_3.length",
    "np.np1.length",
    "sample_allele.d.inverted",
    "p.d_5.length",
    "p.d_3.length",
    "np.np2.length",
    "p.j_5.length",
];

#[derive(Clone)]
struct Fixture {
    v: Vec<(&'static [u8], usize)>,
    d: Option<&'static [u8]>,
    j: (&'static [u8], usize),
    choices: Vec<Vec<u32>>,
    bases: [NpBaseSupport; 2],
}

impl Fixture {
    fn vj(v: &'static [u8], anchor: usize) -> Self {
        let mut choices = vec![vec![0]; 14];
        choices[7] = vec![0, 1, 2];
        choices[8] = vec![0, 1, 2];
        choices[13] = vec![0, 1];
        Self {
            v: vec![(v, anchor)],
            d: None,
            j: (b"TTTGGG", 0),
            choices,
            bases: [NpBaseSupport {
                first: 15,
                transitions: [15; 4],
            }; 2],
        }
    }

    fn active(&self, index: usize) -> bool {
        self.d.is_some() || ![1, 4, 5, 9, 10, 11, 12].contains(&index)
    }

    fn reference(&self) -> RefDataConfig {
        let mut cfg = RefDataConfig::empty(if self.d.is_some() {
            ChainType::Vdj
        } else {
            ChainType::Vj
        });
        let allele = |segment, seq: &[u8], anchor| Allele {
            name: format!("{segment:?}*01"),
            gene: format!("{segment:?}"),
            seq: seq.to_vec(),
            segment,
            anchor,
            functional_status: None,
            subregions: Vec::new(),
        };
        for &(sequence, anchor) in &self.v {
            let _ = cfg
                .v_pool
                .push(allele(Segment::V, sequence, Some(anchor as u16)));
        }
        if let Some(d) = self.d {
            let _ = cfg.d_pool.push(allele(Segment::D, d, None));
        }
        let _ = cfg
            .j_pool
            .push(allele(Segment::J, self.j.0, Some(self.j.1 as u16)));
        cfg
    }

    fn spec(&self) -> PProductiveFeasibility {
        let domain = |index| ChoiceDomain::new(Some(index), self.choices[index].clone());
        let segment =
            |s, sample, five: Option<usize>, three: Option<usize>| SegmentCompletionDomain {
                segment: s,
                alleles: ChoiceDomain::new(
                    Some(sample),
                    self.choices[sample]
                        .iter()
                        .map(|&n| AlleleId::new(n))
                        .collect(),
                ),
                trim_5: five
                    .map(&domain)
                    .unwrap_or_else(|| ChoiceDomain::new(None, vec![0])),
                trim_3: three
                    .map(&domain)
                    .unwrap_or_else(|| ChoiceDomain::new(None, vec![0])),
            };
        let mut segments = vec![segment(Segment::V, 0, None, Some(3))];
        let mut p_lengths = vec![(PEnd::V3, domain(7)), (PEnd::J5, domain(13))];
        let mut np = vec![NpCompletionDomain {
            segment: Segment::Np1,
            lengths: domain(8),
            bases: self.bases[0],
        }];
        if self.d.is_some() {
            segments.push(segment(Segment::D, 1, Some(4), Some(5)));
            p_lengths.extend([(PEnd::D5, domain(10)), (PEnd::D3, domain(11))]);
            np.push(NpCompletionDomain {
                segment: Segment::Np2,
                lengths: domain(12),
                bases: self.bases[1],
            });
        }
        segments.push(segment(Segment::J, 2, Some(6), None));
        PProductiveFeasibility {
            segments,
            processing: vec![],
            p_lengths,
            np,
            d_orientation: ChoiceDomain::new(
                Some(9),
                self.choices[9].iter().map(|&n| n != 0).collect(),
            ),
            check_frame: true,
            check_stops: true,
            check_coding: true,
            check_v_anchor: true,
            check_j_anchor: true,
        }
    }
}

#[derive(Clone, Debug)]
struct Path {
    choices: Vec<u32>,
    n: [Vec<u8>; 2],
}

fn reverse_complement(bytes: &[u8]) -> Vec<u8> {
    bytes
        .iter()
        .rev()
        .map(|b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => panic!("oracle fixtures contain only canonical bases"),
        })
        .collect()
}

fn stop_free(sequence: &[u8]) -> bool {
    sequence
        .chunks_exact(3)
        .all(|codon| !matches!(codon, b"TAA" | b"TAG" | b"TGA"))
}

/// Materialize all positive-support Markov strings. Tiny fixtures only.
fn strings(length: u32, support: NpBaseSupport) -> Vec<Vec<u8>> {
    let mut strings = vec![vec![]];
    for _ in 0..length {
        strings = strings
            .into_iter()
            .flat_map(|prefix| {
                let mask = prefix.last().map_or(support.first, |base| {
                    support.transitions[b"ACGT".iter().position(|b| b == base).unwrap()]
                });
                (0..4)
                    .filter(move |&bit| mask & (1 << bit) != 0)
                    .map(move |bit| {
                        let mut next = prefix.clone();
                        next.push(b"ACGT"[bit]);
                        next
                    })
            })
            .collect();
    }
    strings
}

fn paths(fixture: &Fixture) -> Vec<Path> {
    let mut combinations = vec![vec![]];
    for domain in &fixture.choices {
        combinations = combinations
            .into_iter()
            .flat_map(|prefix| {
                domain.iter().map(move |&choice| {
                    let mut next = prefix.clone();
                    next.push(choice);
                    next
                })
            })
            .collect();
    }
    combinations
        .into_iter()
        .flat_map(|choices| {
            let first = strings(choices[8], fixture.bases[0]);
            let second = strings(choices[12], fixture.bases[1]);
            first.into_iter().flat_map(move |n1| {
                let choices = choices.clone();
                second.clone().into_iter().map(move |n2| Path {
                    choices: choices.clone(),
                    n: [n1.clone(), n2],
                })
            })
        })
        .collect()
}

fn productive(fixture: &Fixture, path: &Path) -> bool {
    let c = &path.choices;
    let (v_full, v_anchor) = fixture.v[c[0] as usize];
    let Some(v_end) = v_full.len().checked_sub(c[3] as usize) else {
        return false;
    };
    let v = &v_full[..v_end];
    let j_start = c[6] as usize;
    if v_anchor + 3 > v.len() || j_start > fixture.j.1 {
        return false;
    }
    let j = &fixture.j.0[j_start..];
    if !matches!(&v[v_anchor..v_anchor + 3], b"TGT" | b"TGC")
        || !matches!(
            &fixture.j.0[fixture.j.1..fixture.j.1 + 3],
            b"TTT" | b"TTC" | b"TGG"
        )
    {
        return false;
    }
    if c[7] as usize > v.len() || c[13] as usize > j.len() {
        return false;
    }
    let mut sequence = v.to_vec();
    sequence.extend(reverse_complement(&v[v.len() - c[7] as usize..]));
    sequence.extend(&path.n[0]);
    if let Some(d) = fixture.d {
        let Some(end) = d.len().checked_sub(c[5] as usize) else {
            return false;
        };
        if c[4] as usize > end {
            return false;
        }
        let mut retained = d[c[4] as usize..end].to_vec();
        if c[9] != 0 {
            retained = reverse_complement(&retained);
        }
        if c[10] as usize > retained.len() || c[11] as usize > retained.len() {
            return false;
        }
        sequence.extend(reverse_complement(&retained[..c[10] as usize]));
        sequence.extend(&retained);
        sequence.extend(reverse_complement(
            &retained[retained.len() - c[11] as usize..],
        ));
        sequence.extend(&path.n[1]);
    }
    sequence.extend(reverse_complement(&j[..c[13] as usize]));
    let junction_end = sequence.len() + fixture.j.1 - j_start + 3;
    sequence.extend(j);
    (junction_end - v_anchor) % 3 == 0 && stop_free(&sequence[v_anchor % 3..])
}

fn candidate(path: &Path, index: usize) -> ChoiceValue {
    match index {
        0..=2 => ChoiceValue::AlleleId(path.choices[index]),
        9 => ChoiceValue::Bool(path.choices[index] != 0),
        _ => ChoiceValue::Int(i64::from(path.choices[index])),
    }
}

fn prefix_key(fixture: &Fixture, path: &Path, index: usize) -> String {
    let mut key = format!("{index}:");
    for previous in 0..index {
        if !fixture.active(previous) {
            continue;
        }
        if previous == 8 || previous == 12 {
            key.push_str(&format!("{:?};", path.n[usize::from(previous == 12)]));
        } else {
            key.push_str(&format!("{};", path.choices[previous]));
        }
    }
    key
}

fn candidate_key(fixture: &Fixture, path: &Path, index: usize) -> String {
    format!(
        "{}choice={}",
        prefix_key(fixture, path, index),
        path.choices[index]
    )
}

fn base_key(fixture: &Fixture, path: &Path, event: usize, position: usize) -> String {
    let np = usize::from(event == 12);
    format!(
        "{}len={};prefix={:?}",
        prefix_key(fixture, path, event),
        path.n[np].len(),
        &path.n[np][..=position]
    )
}

fn append_np(mut sim: Simulation, segment: Segment, bytes: &[u8], completed: bool) -> Simulation {
    let start = sim.pool.len() as u32;
    for &base in bytes {
        sim = sim
            .with_nucleotide_pushed(Nucleotide::synthetic(base, segment, flag::N_NUC))
            .0;
    }
    if completed {
        let end = sim.pool.len() as u32;
        sim = sim.with_region_added(Region::new(
            segment,
            NucHandle::new(start),
            NucHandle::new(end),
        ));
    }
    sim
}

/// Only assignments and N bytes are read by the completion query. This
/// reconstructs the exact committed state at each stochastic pass, including
/// the zero-length region that a completed empty N pass creates.
fn before(fixture: &Fixture, path: &Path, event: usize) -> Simulation {
    let c = &path.choices;
    let mut sim = Simulation::new();
    for (segment, sample, five, three, p5, p3) in [
        (Segment::V, 0, None, Some(3), None, Some(7)),
        (Segment::D, 1, Some(4), Some(5), Some(10), Some(11)),
        (Segment::J, 2, Some(6), None, Some(13), None),
    ] {
        if sample >= event || (segment == Segment::D && fixture.d.is_none()) {
            continue;
        }
        let committed = |index: Option<usize>| index.filter(|&i| i < event).map_or(0, |i| c[i]);
        let mut instance = AlleleInstance::new(AlleleId::new(c[sample]))
            .with_trim_5(committed(five) as u16)
            .with_trim_3(committed(three) as u16)
            .with_p_length_5(committed(p5))
            .with_p_length_3(committed(p3));
        if segment == Segment::D && 9 < event && c[9] != 0 {
            instance = instance.with_orientation(SegmentOrientation::ReverseComplement);
        }
        sim = sim.with_allele_assigned(segment, instance);
    }
    if 8 < event {
        sim = append_np(sim, Segment::Np1, &path.n[0], true);
    }
    if fixture.d.is_some() && 12 < event {
        sim = append_np(sim, Segment::Np2, &path.n[1], true);
    }
    sim
}

fn check_exhaustive(fixture: &Fixture) -> (usize, usize) {
    let reference = fixture.reference();
    let spec = fixture.spec();
    let all = paths(fixture);
    let accepted: Vec<_> = all
        .iter()
        .filter(|path| productive(fixture, path))
        .collect();
    assert_eq!(
        spec.has_completion(&Simulation::new(), &reference),
        !accepted.is_empty()
    );
    let mut allowed_choices = HashSet::new();
    let mut allowed_bases = HashSet::new();
    for path in &accepted {
        for event in 0..14 {
            if fixture.active(event) {
                allowed_choices.insert(candidate_key(fixture, path, event));
            }
        }
        for (np, event) in [(0, 8), (1, 12)] {
            if !fixture.active(event) {
                continue;
            }
            for position in 0..path.n[np].len() {
                allowed_bases.insert(base_key(fixture, path, event, position));
            }
        }
    }
    let mut checked_choices = HashSet::new();
    let mut checked_bases = HashSet::new();
    for path in &all {
        for (event, address) in ADDRESSES.iter().enumerate() {
            if !fixture.active(event) {
                continue;
            }
            let key = candidate_key(fixture, path, event);
            if checked_choices.insert(key.clone()) {
                assert_eq!(
                    spec.admits(
                        event,
                        &before(fixture, path, event),
                        Some(&reference),
                        address,
                        &candidate(path, event)
                    ),
                    allowed_choices.contains(&key),
                    "candidate mismatch: {key}",
                );
            }
        }
        for (np, event, segment) in [(0, 8, Segment::Np1), (1, 12, Segment::Np2)] {
            if !fixture.active(event) {
                continue;
            }
            for position in 0..path.n[np].len() {
                let key = base_key(fixture, path, event, position);
                if !checked_bases.insert(key.clone()) {
                    continue;
                }
                let sim = append_np(
                    before(fixture, path, event),
                    segment,
                    &path.n[np][..position],
                    false,
                );
                assert_eq!(
                    spec.admits_np_base(
                        event,
                        &sim,
                        &reference,
                        segment,
                        position as u32,
                        path.n[np].len() as u32,
                        path.n[np][position]
                    ),
                    allowed_bases.contains(&key),
                    "N-base mismatch: {key}",
                );
            }
        }
    }
    (all.len(), accepted.len())
}

#[test]
fn vj_all_choices_and_n_prefixes_match_complete_sequence_oracle() {
    let mut fixture = Fixture::vj(b"TGTA", 0);
    fixture.v.push((b"TGTT", 0));
    fixture.choices[0] = vec![0, 1];
    fixture.choices[3] = vec![0, 1];
    fixture.j = (b"CTTTGGG", 1);
    fixture.choices[6] = vec![0, 1];
    let (all, productive) = check_exhaustive(&fixture);
    assert!(productive > 0 && productive < all);
}

#[test]
fn vdj_trim_inversion_p_and_both_n_prefixes_match_complete_sequence_oracle() {
    let mut fixture = Fixture::vj(b"TGT", 0);
    fixture.d = Some(b"TAC");
    fixture.choices[4] = vec![0, 1];
    fixture.choices[5] = vec![0, 1];
    fixture.choices[7] = vec![0, 1];
    fixture.choices[8] = vec![0, 1];
    fixture.choices[9] = vec![0, 1];
    fixture.choices[10] = vec![0, 1];
    fixture.choices[11] = vec![0, 1];
    fixture.choices[12] = vec![0, 1];
    // Categorical zero weights must not become viable base alternatives.
    fixture.bases = [NpBaseSupport {
        first: 3,
        transitions: [3; 4],
    }; 2];
    let (all, productive) = check_exhaustive(&fixture);
    assert!(productive > 0 && productive < all);
}

#[test]
fn markov_lookahead_rejects_prefix_whose_forced_suffix_creates_stop() {
    let mut fixture = Fixture::vj(b"TGTT", 0);
    fixture.choices[7] = vec![0];
    fixture.choices[8] = vec![2];
    fixture.choices[13] = vec![3];
    // First A or C; both must next emit A. V ends T, so AA yields TAA,
    // while CA yields TCA. A one-base/local-stop check admits both firsts.
    fixture.bases[0] = NpBaseSupport {
        first: 3,
        transitions: [1; 4],
    };
    assert_eq!(check_exhaustive(&fixture), (2, 1));
}

#[test]
fn future_p_template_rejects_nonstop_n_prefix_before_completion() {
    let mut fixture = Fixture::vj(b"TGTT", 0);
    fixture.choices[7] = vec![0];
    fixture.choices[8] = vec![1];
    fixture.choices[13] = vec![1];
    // The future J-side P is A; N=A/G completes TAA/TGA.
    assert_eq!(check_exhaustive(&fixture), (4, 2));
}

#[test]
fn invalid_positive_p_length_and_unavoidable_stop_have_no_completion() {
    let mut too_long = Fixture::vj(b"TGT", 0);
    too_long.choices[7] = vec![4];
    too_long.choices[8] = vec![0];
    too_long.choices[13] = vec![0];
    assert_eq!(check_exhaustive(&too_long), (1, 0));
    let mut stop = Fixture::vj(b"TGTT", 0);
    stop.choices[7] = vec![2];
    stop.choices[8] = vec![0];
    stop.choices[13] = vec![0];
    assert_eq!(check_exhaustive(&stop), (1, 0));
}

#[test]
fn complete_coding_flanks_outside_junction_are_part_of_oracle() {
    let mut fixture = Fixture::vj(b"TAATGT", 3);
    fixture.v.push((b"CCCTGT", 3));
    fixture.choices[0] = vec![0, 1];
    let (all, productive) = check_exhaustive(&fixture);
    assert!(productive > 0 && productive < all / 2);
    fixture.j = (b"TTTTAA", 0);
    assert_eq!(check_exhaustive(&fixture).1, 0);
}
