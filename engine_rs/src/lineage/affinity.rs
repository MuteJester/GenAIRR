//! Affinity model for clonal selection: BLOSUM62-weighted, region-weighted
//! amino-acid distance to a target antigen, and target generation.

use crate::ir::{compute_codon_rail, Simulation};
use crate::rng::Rng;

/// 20 standard amino acids in the BLOSUM62 matrix's row/column order.
const AA_ORDER: &[u8; 20] = b"ARNDCQEGHILKMFPSTWYV";

/// BLOSUM62 substitution score for two amino-acid bytes (uppercase one-letter
/// codes). Unknown / non-standard residues (including `*` and `X`) score as a
/// strong mismatch so they contribute maximally to distance.
pub fn blosum62(a: u8, b: u8) -> i8 {
    fn idx(aa: u8) -> Option<usize> {
        AA_ORDER.iter().position(|&x| x == aa)
    }
    // Canonical BLOSUM62, rows/cols in AA_ORDER.
    const M: [[i8; 20]; 20] = [
        // A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V
        [ 4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0],
        [-1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3],
        [-2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3],
        [-2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3],
        [ 0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1],
        [-1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2],
        [-1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2],
        [ 0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3],
        [-2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3],
        [-1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3],
        [-1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1],
        [-1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2],
        [-1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1],
        [-2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1],
        [-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2],
        [ 1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2],
        [ 0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0],
        [-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3],
        [-2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1],
        [ 0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4],
    ];
    match (idx(a), idx(b)) {
        (Some(i), Some(j)) => M[i][j],
        _ => -4,
    }
}

/// Per-position substitution cost from BLOSUM62: `S(a,a) - S(a,b)`, which is 0
/// when `b == a` and grows as the substitution becomes less favorable.
fn substitution_cost(a: u8, b: u8) -> f64 {
    (blosum62(a, a) as f64 - blosum62(a, b) as f64).max(0.0)
}

/// Region-weighted BLOSUM62 amino-acid distance between two aa sequences.
/// `weights` is per amino-acid position; positions beyond the shorter sequence
/// are ignored, and weights shorter than the sequences default to 1.0.
pub fn weighted_aa_distance(a: &[u8], b: &[u8], weights: &[f64]) -> f64 {
    let n = a.len().min(b.len());
    let mut d = 0.0;
    for i in 0..n {
        let w = weights.get(i).copied().unwrap_or(1.0);
        d += w * substitution_cost(a[i], b[i]);
    }
    d
}

/// Generate a "mature" target: the naive aa sequence with up to `m` random
/// single-residue substitutions to a different standard amino acid.
/// Deterministic for the given `rng` state.
pub fn make_mature_target(naive_aa: &[u8], m: u32, rng: &mut Rng) -> Vec<u8> {
    let mut target = naive_aa.to_vec();
    if target.is_empty() {
        return target;
    }
    for _ in 0..m {
        // Unbiased uniform draws (Lemire `range_u32`), matching the engine's
        // RNG convention rather than `% len` modulo bias.
        let pos = rng.range_u32(target.len() as u32) as usize;
        let cur = target[pos];
        // Pick a standard amino acid uniformly among the 19 != cur (no skew).
        let pick = match AA_ORDER.iter().position(|&x| x == cur) {
            Some(ci) => {
                let r = rng.range_u32(19) as usize;
                AA_ORDER[if r >= ci { r + 1 } else { r }]
            }
            None => AA_ORDER[rng.range_u32(20) as usize],
        };
        target[pos] = pick;
    }
    target
}

/// Translate a `Simulation` to its amino-acid sequence by concatenating the
/// in-frame translation of each region (respecting each region's frame phase),
/// in biological order. Used to score a cell's affinity to a target antigen.
pub fn sim_to_aa(sim: &Simulation) -> Vec<u8> {
    let mut aa = Vec::new();
    for region in sim.sequence.regions.iter() {
        let rail = compute_codon_rail(region, &sim.pool);
        aa.extend_from_slice(&rail.amino_acids);
    }
    aa
}

/// Affinity-selection model: maps a cell's amino-acid sequence to an affinity in
/// (0, 1] (1 = identical to target) and to a fitness multiplier for its offspring
/// rate. `selection_strength == 0` makes fitness identically 1 (neutral).
#[derive(Clone, Debug)]
pub struct AffinityModel {
    target_aa: Vec<u8>,
    aa_weights: Vec<f64>,
    beta: f64,
    selection_strength: f64,
    /// Affinity of the founder; fitness is measured relative to this baseline so
    /// the founder has fitness ~1 and cells that improve on it exceed 1.
    founder_affinity: f64,
}

impl AffinityModel {
    /// Build a model. `founder_aa` is the naive founder's aa sequence; its
    /// affinity becomes the fitness baseline.
    pub fn new(
        target_aa: Vec<u8>,
        aa_weights: Vec<f64>,
        beta: f64,
        selection_strength: f64,
        founder_aa: &[u8],
    ) -> Self {
        let founder_affinity =
            (-beta * weighted_aa_distance(founder_aa, &target_aa, &aa_weights)).exp();
        Self {
            target_aa,
            aa_weights,
            beta,
            selection_strength,
            founder_affinity,
        }
    }

    /// Affinity in (0, 1]: `exp(-beta * region-weighted BLOSUM distance to target)`.
    pub fn affinity_value(&self, aa: &[u8]) -> f64 {
        (-self.beta * weighted_aa_distance(aa, &self.target_aa, &self.aa_weights)).exp()
    }

    /// Fitness from an already-computed affinity value (avoids re-translation).
    pub fn fitness_from_affinity(&self, affinity_value: f64) -> f64 {
        (1.0 + self.selection_strength * (affinity_value - self.founder_affinity)).max(0.0)
    }

    /// Fitness multiplier for offspring rate: `1 + strength * (affinity - founder_affinity)`,
    /// clamped at 0. `selection_strength == 0` ⇒ always 1 (neutral).
    pub fn fitness(&self, aa: &[u8]) -> f64 {
        self.fitness_from_affinity(self.affinity_value(aa))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ir::{Nucleotide, NucHandle, Region, Segment, Simulation};
    use crate::rng::Rng;

    fn aa_founder(seq: &[u8]) -> Simulation {
        let mut sim = Simulation::new();
        for (i, b) in seq.iter().enumerate() {
            let (next, _) = sim.with_nucleotide_pushed(
                Nucleotide::germline(*b, i as u16, Segment::V));
            sim = next;
        }
        sim.with_region_added(Region::new(Segment::V, NucHandle::new(0), NucHandle::new(seq.len() as u32)))
    }

    #[test]
    fn sim_to_aa_translates_in_frame() {
        // 8 'A' bases -> codons AAA, AAA (last 2 ignored) -> "KK" (Lys, Lys)
        let sim = aa_founder(b"AAAAAAAA");
        assert_eq!(sim_to_aa(&sim), b"KK".to_vec());
    }

    #[test]
    fn affinity_value_is_one_at_target() {
        let m = AffinityModel::new(b"KK".to_vec(), vec![1.0; 2], 1.0, 1.0, b"KK");
        assert!((m.affinity_value(b"KK") - 1.0).abs() < 1e-9);
    }

    #[test]
    fn affinity_value_decreases_with_distance() {
        let m = AffinityModel::new(b"KK".to_vec(), vec![1.0; 2], 1.0, 1.0, b"KK");
        let near = m.affinity_value(b"KK");
        let far = m.affinity_value(b"WW");
        assert!(far < near, "far {far} should be < near {near}");
        assert!(far > 0.0);
    }

    #[test]
    fn fitness_is_neutral_when_strength_zero() {
        // selection_strength = 0 => fitness == 1.0 for ANY sequence
        let m = AffinityModel::new(b"KK".to_vec(), vec![1.0; 2], 1.0, 0.0, b"NN");
        assert!((m.fitness(b"KK") - 1.0).abs() < 1e-9);
        assert!((m.fitness(b"WW") - 1.0).abs() < 1e-9);
        assert!((m.fitness(b"NN") - 1.0).abs() < 1e-9);
    }

    #[test]
    fn fitness_rewards_closer_than_founder_penalizes_farther() {
        // founder = "NN" (far from target "KK"); strength 1.
        let m = AffinityModel::new(b"KK".to_vec(), vec![1.0; 2], 1.0, 1.0, b"NN");
        // a cell AT the target is fitter than the founder baseline (fitness > 1)
        assert!(m.fitness(b"KK") > 1.0);
        // a cell equal to the founder has fitness exactly 1 (baseline)
        assert!((m.fitness(b"NN") - 1.0).abs() < 1e-9);
        // fitness never goes negative
        assert!(m.fitness(b"WWWWWWWW") >= 0.0);
    }

    #[test]
    fn blosum62_known_entries() {
        assert_eq!(blosum62(b'A', b'A'), 4);
        assert_eq!(blosum62(b'W', b'W'), 11);
        assert_eq!(blosum62(b'C', b'C'), 9);
        assert_eq!(blosum62(b'A', b'R'), -1);
        assert_eq!(blosum62(b'W', b'C'), -2);
        assert_eq!(blosum62(b'D', b'E'), 2);
        assert_eq!(blosum62(b'R', b'A'), blosum62(b'A', b'R'));
        let _ = blosum62(b'*', b'A');
        let _ = blosum62(b'X', b'X');
    }

    #[test]
    fn weighted_distance_zero_for_identical() {
        let a = b"ACDEW";
        let w = vec![1.0; 5];
        assert!(weighted_aa_distance(a, a, &w).abs() < 1e-9);
    }

    #[test]
    fn weighted_distance_increases_with_substitutions_and_respects_weights() {
        let a = b"AAAAA";
        let b = b"AAAAW";
        let flat = vec![1.0; 5];
        let d_flat = weighted_aa_distance(a, b, &flat);
        assert!(d_flat > 0.0);
        let mut heavy = vec![1.0; 5];
        heavy[4] = 10.0;
        assert!(weighted_aa_distance(a, b, &heavy) > d_flat);
        let mut zero = vec![1.0; 5];
        zero[4] = 0.0;
        assert!(weighted_aa_distance(a, b, &zero).abs() < 1e-9);
    }

    #[test]
    fn mature_target_applies_m_substitutions_deterministically() {
        let naive = b"ACDEFGHIKLMNPQRSTVWY".to_vec();
        let mut rng = Rng::new(42);
        let t = make_mature_target(&naive, 3, &mut rng);
        assert_eq!(t.len(), naive.len());
        let diffs = naive.iter().zip(&t).filter(|(a, b)| a != b).count();
        assert!(diffs >= 1 && diffs <= 3, "expected up to 3 aa substitutions, got {diffs}");
        let mut rng2 = Rng::new(42);
        assert_eq!(make_mature_target(&naive, 3, &mut rng2), t);
    }
}
