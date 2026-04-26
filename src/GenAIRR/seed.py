"""
Seed Management for GenAIRR.

Provides a process-wide default seed for reproducible sequence simulation.
The seed flows into the C simulator via Experiment.run/compile/run_to_file:
when the caller does not pass an explicit ``seed=`` argument, the global
value set by ``set_seed()`` is used instead. If neither is set, the C
engine auto-seeds from wall time + simulator address (different streams
per concurrent simulator).

The C simulator owns its own PCG32 RNG state per instance (T0-5), so two
simulators with the same explicit seed produce byte-identical output.

Note:
    set_seed() also seeds Python's ``random`` and ``numpy.random`` for
    any auxiliary code paths that draw outside the C engine.
    An explicit ``seed=0`` does NOT trigger the global fallback — it is
    passed through to the C engine, which interprets 0 as auto-seed.

Usage:
    >>> from GenAIRR import set_seed, Experiment
    >>> set_seed(42)
    >>> r1 = Experiment.on("human_igh").run(n=10)   # uses global seed 42
    >>> r2 = Experiment.on("human_igh").run(n=10)   # also uses 42 → identical
    >>> assert r1[0]["sequence"] == r2[0]["sequence"]
"""

import random
from typing import Optional


_current_seed: Optional[int] = None


def set_seed(seed: int) -> None:
    """
    Set the process-wide default seed for reproducible sequence simulation.

    Subsequent ``Experiment.run/compile/run_to_file`` calls that do NOT
    pass an explicit ``seed=`` argument will use this value when seeding
    the C simulator's PCG32 RNG. Also seeds Python's ``random`` and
    ``numpy.random`` for any auxiliary draws outside the C engine.

    Args:
        seed: Integer seed value. Same seed → identical output.

    Example:
        >>> from GenAIRR import set_seed, Experiment
        >>> set_seed(42)
        >>> r1 = Experiment.on("human_igh").run(n=10)
        >>> r2 = Experiment.on("human_igh").run(n=10)
        >>> assert r1[0]["sequence"] == r2[0]["sequence"]

    Note:
        An explicit ``seed=`` on ``Experiment.run/compile`` always
        overrides the global value. Pass ``seed=0`` to the C engine
        to request auto-seeding (time + simulator address) regardless
        of any global ``set_seed`` value.
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
