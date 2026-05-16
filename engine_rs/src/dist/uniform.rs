//! `UniformBase` — uniform distribution over the four canonical DNA bases.

use super::Distribution;
use crate::rng::Rng;

/// Uniform distribution over the four canonical DNA bases
/// `{b'A', b'C', b'G', b'T'}` (uppercase). Each base has probability
/// 0.25.
///
/// This is the default fallback when no empirical TdT / NP-base
/// model is configured. Empirical transition matrices override at
/// runtime through the same trait surface.
#[derive(Clone, Debug, Default)]
pub struct UniformBase;

impl Distribution for UniformBase {
    type Output = u8;

    fn sample(&self, rng: &mut Rng) -> u8 {
        const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
        BASES[rng.range_u32(4) as usize]
    }

    fn support(&self) -> Option<Vec<(u8, f64)>> {
        Some(vec![(b'A', 1.0), (b'C', 1.0), (b'G', 1.0), (b'T', 1.0)])
    }
}
