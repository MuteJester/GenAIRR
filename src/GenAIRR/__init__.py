"""GenAIRR — synthetic immune-receptor-sequence simulator.

The user-facing API is the fluent :class:`Experiment` builder:

    >>> import GenAIRR as ga
    >>> outcomes = ga.Experiment.on("human_igh").recombine().run(n=100, seed=42)

Internally, ``Experiment`` compiles to a Rust ``PassPlan`` (via
:mod:`genairr_engine`) and runs it through the V6 simulation kernel.
The legacy V5 C engine has been retired in favour of the Rust kernel.
"""
from importlib.metadata import (
    PackageNotFoundError as _PackageNotFoundError,
    version as _pkg_version,
)

# Single source of truth for the version is pyproject.
try:
    __version__ = _pkg_version("GenAIRR")
except _PackageNotFoundError:
    # Fallback for editable installs where metadata is missing.
    __version__ = "0.0.0"


# ──────────────────────────────────────────────────────────────────
# Public API
# ──────────────────────────────────────────────────────────────────

# The simulation entry point.
from .experiment import CompiledExperiment, Experiment, dataconfig_to_refdata

# Canonical contract bundle. Re-exported here so users write
# ``GenAIRR.productive()`` without reaching into ``genairr_engine``.
from genairr_engine import StrictSamplingError, productive

# Reference data — DataConfig + species pickles.
from .data import list_configs
from .dataconfig import ChainType, ConfigInfo, DataConfig, DataConfigError, Species
from .dataconfig.enums import Productivity

# Reproducibility helpers.
from .seed import get_seed, reset_seed, set_seed


# ──────────────────────────────────────────────────────────────────
# Lazy DataConfig attribute access
# ──────────────────────────────────────────────────────────────────
#
# Common builtin configs are accessible as ``GenAIRR.HUMAN_IGH_OGRDB``
# etc. ``__getattr__`` defers the pickle load to first access so plain
# ``import GenAIRR`` stays cheap. The full list of 100+ builtins is
# also available via :func:`list_configs` and ``GenAIRR.data.<NAME>``.

_LAZY_DATA_CONFIGS = frozenset(
    {
        "HUMAN_IGH_OGRDB",
        "HUMAN_IGH_EXTENDED",
        "HUMAN_IGK_OGRDB",
        "HUMAN_IGL_OGRDB",
        "HUMAN_TCRB_IMGT",
    }
)


def __getattr__(name):
    if name in _LAZY_DATA_CONFIGS:
        from . import data as _data

        return getattr(_data, name)
    raise AttributeError(f"module 'GenAIRR' has no attribute {name!r}")


__all__ = [
    # Version
    "__version__",
    # Simulation API
    "Experiment",
    "CompiledExperiment",
    "dataconfig_to_refdata",
    "productive",
    "StrictSamplingError",
    # Reference data
    "DataConfig",
    "DataConfigError",
    "ConfigInfo",
    "ChainType",
    "Species",
    "Productivity",
    "list_configs",
    "HUMAN_IGH_OGRDB",
    "HUMAN_IGH_EXTENDED",
    "HUMAN_IGK_OGRDB",
    "HUMAN_IGL_OGRDB",
    "HUMAN_TCRB_IMGT",
    # Reproducibility
    "set_seed",
    "get_seed",
    "reset_seed",
]
