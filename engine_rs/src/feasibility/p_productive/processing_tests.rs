//! Independent literal-sequence enumeration of the original-end solver.
use super::*;
use crate::assignment::AlleleInstance;
use crate::end_processing::EndOutcome;
use crate::refdata::{Allele, ChainType};

fn rc(bytes: &[u8]) -> Vec<u8> {
    bytes
        .iter()
        .rev()
        .map(|base| match base {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => unreachable!(),
        })
        .collect()
}

// Deliberately does not call the production materializer or codon automaton.
fn literal(seq: &[u8], outcome: EndOutcome, reverse: bool) -> Option<(Vec<u8>, usize)> {
    let five = outcome.five.unsigned_abs() as usize;
    let three = outcome.three.unsigned_abs() as usize;
    if five > seq.len() || three > seq.len() {
        return None;
    }
    let start = if outcome.five > 0 { five } else { 0 };
    let end = seq.len() - if outcome.three > 0 { three } else { 0 };
    if start > end {
        return None;
    }
    let mut result = if outcome.five < 0 {
        rc(&seq[..five])
    } else {
        vec![]
    };
    let offset = result.len();
    result.extend_from_slice(&seq[start..end]);
    if outcome.three < 0 {
        result.extend(rc(&seq[seq.len() - three..]));
    }
    Some((if reverse { rc(&result) } else { result }, offset))
}

fn stop_free(bytes: &[u8]) -> bool {
    bytes
        .chunks_exact(3)
        .all(|codon| !matches!(codon, b"TAA" | b"TAG" | b"TGA"))
}

fn words(length: usize, markov: bool) -> Vec<Vec<u8>> {
    let mut result = vec![vec![]];
    for _ in 0..length {
        result = result
            .into_iter()
            .flat_map(|prefix| {
                let alphabet: &[u8] = if !markov {
                    b"ACGT"
                } else if prefix.last() == Some(&b'A') {
                    b"C"
                } else {
                    b"A"
                };
                alphabet.iter().map(move |&base| {
                    let mut next = prefix.clone();
                    next.push(base);
                    next
                })
            })
            .collect();
    }
    result
}

fn run_oracle(with_d: bool, markov: bool) {
    let mut reference = RefDataConfig::empty(if with_d {
        ChainType::Vdj
    } else {
        ChainType::Vj
    });
    let v = b"TGTACG";
    let d = b"AGTC";
    let j = b"AGTTTGGG";
    let add = |segment, seq: &[u8], anchor| Allele {
        name: format!("{segment:?}*01"),
        gene: format!("{segment:?}"),
        seq: seq.to_vec(),
        segment,
        anchor,
        functional_status: None,
        subregions: vec![],
    };
    let _ = reference.v_pool.push(add(Segment::V, v, Some(0)));
    let _ = reference.j_pool.push(add(Segment::J, j, Some(2)));
    if with_d {
        let _ = reference.d_pool.push(add(Segment::D, d, None));
    }
    let row = |five, three| EndOutcome {
        five,
        three,
        weight: 1.0,
    };
    let v_rows = vec![row(0, -2), row(0, -1), row(0, 0), row(0, 1), row(0, 4)];
    let d_rows = if with_d {
        [-2, 0, 2, 4]
            .into_iter()
            .flat_map(|five| [-1, 0, 2, 4].into_iter().map(move |three| row(five, three)))
            .collect()
    } else {
        vec![row(0, 0)]
    };
    let j_rows = vec![row(-2, 0), row(-1, 0), row(0, 0), row(1, 0), row(3, 0)];
    let rows = [v_rows, d_rows, j_rows];
    let segments: Vec<_> = [Segment::V, Segment::D, Segment::J]
        .into_iter()
        .filter(|s| with_d || *s != Segment::D)
        .map(|segment| SegmentCompletionDomain {
            segment,
            alleles: ChoiceDomain::new(Some(0), vec![AlleleId::new(0)]),
            trim_5: ChoiceDomain::new(None, vec![0]),
            trim_3: ChoiceDomain::new(None, vec![0]),
        })
        .collect();
    let bases = if markov {
        NpBaseSupport {
            first: 1,
            transitions: [2, 1, 1, 1],
        }
    } else {
        NpBaseSupport {
            first: 15,
            transitions: [15; 4],
        }
    };
    let np: Vec<_> = [Segment::Np1, Segment::Np2]
        .into_iter()
        .filter(|s| with_d || *s != Segment::Np2)
        .map(|segment| NpCompletionDomain {
            segment,
            lengths: ChoiceDomain::new(
                Some(if segment == Segment::Np1 { 7 } else { 11 }),
                vec![0, 1, 2],
            ),
            bases,
        })
        .collect();
    let spec = PProductiveFeasibility {
        segments,
        p_lengths: vec![],
        processing: [Segment::V, Segment::D, Segment::J]
            .into_iter()
            .enumerate()
            .filter(|(_, s)| with_d || *s != Segment::D)
            .map(|(i, s)| (s, ChoiceDomain::new(Some(3 + i), rows[i].clone())))
            .collect(),
        np,
        d_orientation: ChoiceDomain::new(
            Some(8),
            if with_d {
                vec![false, true]
            } else {
                vec![false]
            },
        ),
        check_frame: true,
        check_stops: true,
        check_coding: true,
        check_v_anchor: true,
        check_j_anchor: true,
    };
    let ns: Vec<_> = (0..=2).flat_map(|n| words(n, markov)).collect();
    let mut valid = Vec::new();
    for vi in 0..rows[0].len() {
        for di in 0..rows[1].len() {
            for ji in 0..rows[2].len() {
                if rows[0][vi].three > 3 || rows[2][ji].five > 2 {
                    continue;
                }
                let Some((vb, _)) = literal(v, rows[0][vi], false) else {
                    continue;
                };
                let Some((jb, jp)) = literal(j, rows[2][ji], false) else {
                    continue;
                };
                for reverse in [false, true] {
                    if reverse && !with_d {
                        continue;
                    }
                    let Some((db, _)) = literal(d, rows[1][di], reverse) else {
                        continue;
                    };
                    for n1 in &ns {
                        for n2 in &ns {
                            if !with_d && !n2.is_empty() {
                                continue;
                            }
                            let mut sequence = vb.clone();
                            sequence.extend(n1);
                            if with_d {
                                sequence.extend(&db);
                                sequence.extend(n2);
                            }
                            let junction_end =
                                sequence.len() + jp + 2 - rows[2][ji].five.max(0) as usize + 3;
                            sequence.extend(&jb);
                            if junction_end % 3 == 0 && stop_free(&sequence) {
                                valid.push([vi, di, ji]);
                            }
                        }
                    }
                }
            }
        }
    }
    assert!(!valid.is_empty());
    assert!(spec.has_completion(&Simulation::new(), &reference));
    let mut sim = Simulation::new();
    for segment in [Segment::V, Segment::D, Segment::J] {
        if with_d || segment != Segment::D {
            sim = sim.with_allele_assigned(segment, AlleleInstance::new(AlleleId::new(0)));
        }
    }
    for vi in 0..rows[0].len() {
        let expected = valid.iter().any(|path| path[0] == vi);
        assert_eq!(
            spec.admits(
                3,
                &sim,
                Some(&reference),
                "end_processing.v1.v.outcome",
                &ChoiceValue::Int(vi as i64)
            ),
            expected
        );
        let v_sim = sim.with_coding_end_choice(
            Segment::V,
            AlleleId::new(0),
            rows[0][vi].five,
            rows[0][vi].three,
            vi,
        );
        for di in 0..rows[1].len() {
            if with_d {
                let expected = valid.iter().any(|path| path[0] == vi && path[1] == di);
                assert_eq!(
                    spec.admits(
                        4,
                        &v_sim,
                        Some(&reference),
                        "end_processing.v1.d.outcome",
                        &ChoiceValue::Int(di as i64)
                    ),
                    expected
                );
            }
            let d_sim = if with_d {
                v_sim.with_coding_end_choice(
                    Segment::D,
                    AlleleId::new(0),
                    rows[1][di].five,
                    rows[1][di].three,
                    di,
                )
            } else {
                v_sim.clone()
            };
            for ji in 0..rows[2].len() {
                let expected = valid.iter().any(|path| *path == [vi, di, ji]);
                assert_eq!(
                    spec.admits(
                        5,
                        &d_sim,
                        Some(&reference),
                        "end_processing.v1.j.outcome",
                        &ChoiceValue::Int(ji as i64)
                    ),
                    expected,
                    "with_d={with_d}, markov={markov}, outcomes={vi}/{di}/{ji}"
                );
            }
        }
    }
    assert!(!spec.admits(
        3,
        &sim,
        Some(&reference),
        "end_processing.v1.v.outcome",
        &ChoiceValue::Int(99)
    ));
}

#[test]
fn original_end_vj_support_matches_literal_enumeration() {
    run_oracle(false, false);
}
#[test]
fn original_end_vdj_support_matches_literal_enumeration() {
    run_oracle(true, false);
}
#[test]
fn original_end_vj_markov_support_matches_literal_enumeration() {
    run_oracle(false, true);
}
#[test]
fn original_end_vdj_markov_support_matches_literal_enumeration() {
    run_oracle(true, true);
}
