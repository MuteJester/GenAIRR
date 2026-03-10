"""
Seed Management for GenAIRR

Provides centralized control over random number generation for reproducible
sequence simulation. Setting a seed ensures identical results across runs.

Usage:
    >>> from GenAIRR import set_seed
    >>> set_seed(42)
    >>> # All subsequent pipeline executions will be reproducible

    >>> from GenAIRR import Experiment, set_seed
    >>> set_seed(42)
    >>> r1 = Experiment.on("human_igh").somatic_hypermutation().run(n=1, seed=42)
    >>> r2 = Experiment.on("human_igh").somatic_hypermutation().run(n=1, seed=42)
    >>> assert r1[0]["sequence"] == r2[0]["sequence"]  # Always True
"""

import random
from typing import Optional


_current_seed: Optional[int] = None


def set_seed(seed: int) -> None:
    """
    Set the global random seed for reproducible sequence simulation.

    This function seeds both Python's built-in `random` module and NumPy's
    random number generator, ensuring all random operations in GenAIRR
    produce identical results when given the same seed.

    Args:
        seed: Integer seed value. Use the same seed to reproduce exact results.

    Example:
        >>> from GenAIRR import set_seed, Experiment
        >>> set_seed(42)
        >>> result = Experiment.on("human_igh").somatic_hypermutation().run(n=10, seed=42)
        >>> # Running the above with seed 42 will always produce the same sequences

    Note:
        Call set_seed() before any simulation to ensure reproducibility.
        Each call resets the random state, so you can re-seed between runs.
    """
    global _current_seed
    _current_seed = seed
    random.seed(seed)
    try:
        import numpy as np
        np.random.seed(seed)
    except ImportError:
        pass


def get_seed() -> Optional[int]:
    """
    Get the current random seed, if one has been set.

    Returns:
        The seed value passed to set_seed(), or None if set_seed()
        has not been called in this session.

    Example:
        >>> from GenAIRR import set_seed, get_seed
        >>> get_seed()  # Returns None
        >>> set_seed(42)
        >>> get_seed()  # Returns 42
    """
    return _current_seed


def reset_seed() -> None:
    """
    Clear the stored seed and re-randomize the generators.

    This restores non-deterministic behavior by seeding from system entropy.
    Useful when you want to return to random behavior after reproducible runs.

    Example:
        >>> from GenAIRR import set_seed, reset_seed
        >>> set_seed(42)
        >>> # ... reproducible operations ...
        >>> reset_seed()
        >>> # ... back to random behavior ...
    """
    global _current_seed
    _current_seed = None
    # Seed from system entropy (None triggers this behavior)
    random.seed(None)
    try:
        import numpy as np
        np.random.seed(None)
    except ImportError:
        pass
