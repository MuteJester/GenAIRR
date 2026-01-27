"""
Mutation Result

Immutable dataclass representing the result of applying mutations to a sequence.
This provides a clear, typed return contract for all mutation models.
"""

from dataclasses import dataclass, field
from typing import Dict


@dataclass(frozen=True)
class MutationResult:
    """
    Immutable result from applying mutations to a sequence.

    This dataclass provides a clear contract for what mutation models return,
    replacing the previous untyped tuple(str, dict, float) return type.

    Attributes:
        mutated_sequence: The sequence after mutations have been applied.
        mutations: Dictionary mapping position -> "original>mutated" (e.g., {42: "A>T"}).
        mutation_rate: The actual mutation rate that was applied (sampled from min/max range).
        mutation_count: The number of mutations that were made.

    Example:
        >>> result = MutationResult(
        ...     mutated_sequence="ATCGATCG",
        ...     mutations={2: "A>T", 5: "G>C"},
        ...     mutation_rate=0.05,
        ...     mutation_count=2
        ... )
        >>> result.mutated_sequence
        'ATCGATCG'
        >>> result.mutations
        {2: 'A>T', 5: 'G>C'}
    """

    mutated_sequence: str
    mutations: Dict[int, str] = field(default_factory=dict)
    mutation_rate: float = 0.0
    mutation_count: int = 0

    def __post_init__(self):
        """Validate mutation_count matches mutations dict if not explicitly set."""
        # Since frozen, we can't modify after creation, but we validate consistency
        if self.mutation_count == 0 and self.mutations:
            # Use object.__setattr__ to bypass frozen restriction during init
            object.__setattr__(self, 'mutation_count', len(self.mutations))
