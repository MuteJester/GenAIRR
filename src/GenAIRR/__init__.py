from importlib.metadata import version as _pkg_version, PackageNotFoundError as _PackageNotFoundError

# Expose package version without hardcoding (single source of truth in pyproject)
try:
    __version__ = _pkg_version("GenAIRR")
except _PackageNotFoundError:
    # Fallback for editable installs during development when metadata is missing
    __version__ = "0.0.0"

# =============================================================================
# Top-level exports for convenient access
# =============================================================================

# Experiment DSL — the single entry point for simulation
from .experiment import Experiment

# CompiledSimulator — reusable handle returned by Experiment.compile()
# SimulationStream — streaming interface returned by CompiledSimulator.stream()
from .protocol import CompiledSimulator, SimulationStream

# Step descriptors are an internal IR layer produced by the clause →
# step compiler in protocol.py. They are NOT the user API — pass
# clauses from `GenAIRR.ops` (e.g. `rate(...)`, `model('s5f')`,
# `with_indels()`) to the phase methods (`.mutate(...)`, etc.).
# Advanced/internal users can still import from `GenAIRR.steps.*`.

# Seed management for reproducibility
from .seed import set_seed, get_seed, reset_seed

# Data configs (lazy-loaded) — most commonly used configs.
# All 100+ configs are accessible via GenAIRR.data or Experiment.on("config_name").
# Use GenAIRR.data.list_configs() to see all available configs.
from .data import (
    HUMAN_IGH_OGRDB,
    HUMAN_IGH_EXTENDED,
    HUMAN_IGK_OGRDB,
    HUMAN_IGL_OGRDB,
    HUMAN_TCRB_IMGT,
)
from .data import list_configs

# DataConfig and metadata
from .dataconfig import DataConfig, DataConfigError, ConfigInfo, ChainType, Species
from .dataconfig.enums import Productivity

# Batch result wrapper + narration
from .result import SimulationResult, narrate, narrate_from_record

# Visualization
from .utilities.visualize import visualize_sequence


# Define what gets exported with `from GenAIRR import *`
__all__ = [
    # Version
    "__version__",
    # Core API
    "Experiment",
    "CompiledSimulator",
    "SimulationStream",
    "SimulationResult",
    "narrate",
    "narrate_from_record",
    # Step classes are intentionally NOT exported at the top level —
    # they're an internal IR. Use `GenAIRR.ops` clauses instead, or
    # `from GenAIRR.steps import ...` for advanced/internal use.
    # Seed management
    "set_seed",
    "get_seed",
    "reset_seed",
    # Data configs (common) + discovery
    "HUMAN_IGH_OGRDB",
    "HUMAN_IGH_EXTENDED",
    "HUMAN_IGK_OGRDB",
    "HUMAN_IGL_OGRDB",
    "HUMAN_TCRB_IMGT",
    "list_configs",
    # Visualization
    "visualize_sequence",
    # DataConfig and metadata
    "DataConfig",
    "DataConfigError",
    "ConfigInfo",
    "ChainType",
    "Species",
    "Productivity",
]
