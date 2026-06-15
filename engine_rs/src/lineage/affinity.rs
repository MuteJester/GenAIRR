//! Affinity model for clonal selection: BLOSUM62-weighted, region-weighted
//! amino-acid distance to a target antigen, and target generation.

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
        let pos = (rng.next_u64() % target.len() as u64) as usize;
        let cur = target[pos];
        // pick a standard amino acid different from the current residue
        let start = (rng.next_u64() % 20) as usize;
        let mut pick = AA_ORDER[start];
        if pick == cur {
            pick = AA_ORDER[(start + 1) % 20];
        }
        target[pos] = pick;
    }
    target
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rng::Rng;

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
