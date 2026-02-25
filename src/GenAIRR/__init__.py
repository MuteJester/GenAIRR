from importlib.metadata import version as _pkg_version, PackageNotFoundError as _PackageNotFoundError
from typing import Union, List

# Expose package version without hardcoding (single source of truth in pyproject)
try:
    __version__ = _pkg_version("GenAIRR")
except _PackageNotFoundError:
    # Fallback for editable installs during development when metadata is missing
    __version__ = "0.0.0"

# =============================================================================
# Top-level exports for convenient access
# =============================================================================

# Pipeline - the main entry point for sequence simulation
from .pipeline import AugmentationPipeline
# Alias for cleaner API
Pipeline = AugmentationPipeline

# Steps namespace - allows `from GenAIRR import steps` then `steps.SimulateSequence`
from . import steps

# Mutation models
from .mutation import S5F, Uniform

# Seed management for reproducibility
from .seed import set_seed, get_seed, reset_seed

# Data configs (lazy-loaded)
from .data import (
    HUMAN_IGH_OGRDB,
    HUMAN_IGH_EXTENDED,
    HUMAN_IGK_OGRDB,
    HUMAN_IGL_OGRDB,
    HUMAN_TCRB_IMGT,
)

# DataConfig and metadata
from .dataconfig import DataConfig, DataConfigError, ConfigInfo, ChainType, Species

# Container for results
from .container import SimulationContainer


def simulate(config, mutation_model, productive: bool = True, n: int = 1):
    """
    One-line convenience function for generating simulated AIRR sequences.

    This creates a minimal pipeline with SimulateSequence and position correction
    steps, suitable for most common use cases.

    Args:
        config: DataConfig to use (e.g., HUMAN_IGH_OGRDB, HUMAN_IGK_OGRDB)
        mutation_model: Mutation model instance (e.g., S5F(0.003, 0.25))
        productive: If True, only generate productive sequences (default: True)
        n: Number of sequences to generate (default: 1)

    Returns:
        SimulationContainer if n=1, or List[SimulationContainer] if n>1

    Example:
        >>> from GenAIRR import simulate, HUMAN_IGH_OGRDB, S5F
        >>> result = simulate(HUMAN_IGH_OGRDB, S5F(0.003, 0.25))
        >>> print(result.sequence)

        >>> # Generate multiple sequences
        >>> results = simulate(HUMAN_IGH_OGRDB, S5F(0.003, 0.25), n=100)
    """
    from .steps import SimulateSequence, FixVPositionAfterTrimmingIndexAmbiguity

    pipeline = Pipeline(
        config=config,
        steps=[
            SimulateSequence(mutation_model, productive),
            FixVPositionAfterTrimmingIndexAmbiguity(),
        ]
    )

    if n == 1:
        return pipeline.execute()
    return [pipeline.execute() for _ in range(n)]


# Define what gets exported with `from GenAIRR import *`
__all__ = [
    # Version
    "__version__",
    # Core classes
    "Pipeline",
    "AugmentationPipeline",
    "SimulationContainer",
    # Steps namespace
    "steps",
    # Mutation models
    "S5F",
    "Uniform",
    # Seed management
    "set_seed",
    "get_seed",
    "reset_seed",
    # Data configs
    "HUMAN_IGH_OGRDB",
    "HUMAN_IGH_EXTENDED",
    "HUMAN_IGK_OGRDB",
    "HUMAN_IGL_OGRDB",
    "HUMAN_TCRB_IMGT",
    # DataConfig and metadata
    "DataConfig",
    "DataConfigError",
    "ConfigInfo",
    "ChainType",
    "Species",
    # Convenience function
    "simulate",
]
