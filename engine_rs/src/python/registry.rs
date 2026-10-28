//! The engine's own list of the passes it defines.
//!
//! A completeness check over the biology the engine simulates needs to know
//! what the engine consists of. [`crate::pass::Outcome::pass_names`] reports
//! only the passes a given compiled plan scheduled, so a pass nobody schedules
//! is invisible to it -- precisely the drift such a check exists to catch.
//!
//! The list below names the pass types the crate exports from
//! [`crate::passes`]. Test doubles live behind `cfg(test)` and cannot appear
//! here.

use pyo3::prelude::*;

/// Every pass the engine defines, sorted, with no duplicates.
///
/// Kept beside the re-export block in [`crate::passes`]: a pass added there
/// and not here fails the catalog completeness test on the Python side.
const PASS_NAMES: &[&str] = &[
    "AssembleSegmentPass",
    "ContaminantPass",
    "DuplexRecombinePass",
    "EchoPass",
    "EmitCodingEndPPass",
    "EndLossPass",
    "GenerateNPPass",
    "IndelPass",
    "InvertDPass",
    "NCorruptionPass",
    "PAdditionPass",
    "PCRErrorPass",
    "PairedEndSamplingPass",
    "ProcessCodingEndsPass",
    "QualityErrorPass",
    "ReceptorRevisionPass",
    "RevCompPass",
    "S5FMutationPass",
    "SampleAllelePass",
    "SampleBasePass",
    "SampleGenotypePass",
    "TrimPass",
    "UniformMutationPass",
];

/// The engine's pass registry, as a sorted tuple of names.
#[pyfunction]
pub fn pass_registry() -> Vec<String> {
    PASS_NAMES.iter().map(|s| (*s).to_string()).collect()
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(pyo3::wrap_pyfunction!(pass_registry, m)?)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::PASS_NAMES;

    #[test]
    fn registry_is_sorted_and_unique() {
        let mut sorted = PASS_NAMES.to_vec();
        sorted.sort_unstable();
        sorted.dedup();
        assert_eq!(sorted.as_slice(), PASS_NAMES, "PASS_NAMES must be sorted and unique");
    }

    #[test]
    fn registry_names_types_the_crate_exports() {
        // A compile-time check that each name is a real exported pass type:
        // referencing the type fails to build if it is renamed or removed.
        fn assert_exists<T>() {}
        assert_exists::<crate::passes::AssembleSegmentPass>();
        assert_exists::<crate::passes::ContaminantPass>();
        assert_exists::<crate::passes::DuplexRecombinePass>();
        assert_exists::<crate::passes::EchoPass>();
        assert_exists::<crate::passes::EmitCodingEndPPass>();
        assert_exists::<crate::passes::EndLossPass>();
        assert_exists::<crate::passes::GenerateNPPass>();
        assert_exists::<crate::passes::IndelPass>();
        assert_exists::<crate::passes::InvertDPass>();
        assert_exists::<crate::passes::NCorruptionPass>();
        assert_exists::<crate::passes::PAdditionPass>();
        assert_exists::<crate::passes::PCRErrorPass>();
        assert_exists::<crate::passes::PairedEndSamplingPass>();
        assert_exists::<crate::passes::ProcessCodingEndsPass>();
        assert_exists::<crate::passes::QualityErrorPass>();
        assert_exists::<crate::passes::ReceptorRevisionPass>();
        assert_exists::<crate::passes::RevCompPass>();
        assert_exists::<crate::passes::S5FMutationPass>();
        assert_exists::<crate::passes::SampleAllelePass>();
        assert_exists::<crate::passes::SampleBasePass>();
        // Not re-exported from `passes`, so named by its module path.
        assert_exists::<crate::passes::sample_genotype::SampleGenotypePass>();
        assert_exists::<crate::passes::TrimPass>();
        assert_exists::<crate::passes::UniformMutationPass>();
    }
}
