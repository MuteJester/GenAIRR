//! Corruption passes — observation-stage perturbations applied
//! after the biology has been simulated.
//!
//! - [`PCRErrorPass`] — random base substitutions modeling PCR
//!   amplification errors (Phase E.4).
//! - [`QualityErrorPass`] — sequencing errors written in lowercase
//!   to mark the position as corrupted (Phase E.5).
//! - [`ContaminantPass`] — wholesale read-level contamination
//!   driven by a Bernoulli switch (Phase E.6).
//! - [`IndelPass`] — insertions and deletions that change pool
//!   length (Phase E.7).

pub mod contaminant;
pub mod end_loss;
pub mod indel;
pub mod pcr;
pub mod quality;
pub mod rev_comp;

pub use contaminant::ContaminantPass;
pub use end_loss::{EndLossPass, LossEnd};
pub use indel::IndelPass;
pub use pcr::PCRErrorPass;
pub use quality::QualityErrorPass;
pub use rev_comp::RevCompPass;
