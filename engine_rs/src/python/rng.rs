//! The keyed seed derivation, as Python calls it.
//!
//! The Python layer splits a clonal run into sub-simulations — a
//! clone's founder, a family's retry attempt, a read, an observed node
//! — and each wants a seed of its own. It cannot reach
//! [`StreamRoot`](crate::rng::substream::StreamRoot) directly, and a derivation
//! written a second time in Python would be a second definition of the
//! engine's one keyed primitive, free to drift from this one. So the
//! primitive is exposed instead of copied, and
//! [`crate::rng::discipline`] keeps watching the only place a stream is
//! born.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use crate::rng::substream::{StreamRoot, SubstreamLabel};

/// The seed of the stream named `(role, index)` under `seed`.
///
/// `index` is a `u32` because that is what a
/// [`SubstreamLabel`] occurrence is; an index that does not fit is an
/// error rather than a silent truncation, which would alias two
/// sub-simulations onto one stream.
#[pyfunction]
#[pyo3(signature = (seed, role, index))]
pub fn substream_seed(seed: u64, role: &str, index: u32) -> PyResult<u64> {
    if role.is_empty() {
        return Err(PyValueError::new_err("role must not be empty"));
    }
    Ok(StreamRoot::new(seed).substream_seed(SubstreamLabel::named(role, index)))
}

/// Register this module's functions on `m`.
///
/// There is no `#[cfg(test)]` block here on purpose: a test binary that
/// touches `PyErr` has to link libpython, which the `extension-module`
/// build does not. What this wrapper computes is pinned from both
/// sides instead — `crate::rng::substream`'s tests assert the word for
/// `("clone", 0)`, and `tests/components/engine/test_substream_seed.py`
/// asserts the same word through the binding, so a wrapper that stops
/// calling the primitive fails there.
pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(pyo3::wrap_pyfunction!(substream_seed, m)?)?;
    Ok(())
}
