"""
Mutation Model Abstract Base Class

Defines the interface that all mutation models must implement.
"""

from abc import ABC, abstractmethod
from typing import Tuple, Dict, Union

from .context import MutationContext
from .result import MutationResult


class MutationModel(ABC):
    """
    Abstract base class for all mutation models.

    This class defines the interface that mutation models must implement.
    It provides both a new typed interface (`apply`) and backwards compatibility
    with the legacy `apply_mutation` method.

    Subclasses must implement either:
    - `apply(context: MutationContext) -> MutationResult` (preferred)
    - `apply_mutation(sequence_object) -> tuple` (legacy, deprecated)

    If only `apply_mutation` is implemented (for backwards compatibility),
    the `apply` method will wrap it automatically.

    Attributes:
        min_mutation_rate: Minimum mutation rate (0.0 to 1.0).
        max_mutation_rate: Maximum mutation rate (0.0 to 1.0).
        productive: Whether to ensure mutations don't create stop codons
                   or destroy conserved anchor residues.

    Example (new style):
        >>> class MyModel(MutationModel):
        ...     @property
        ...     def name(self) -> str:
        ...         return "My Custom Model"
        ...
        ...     def apply(self, context: MutationContext) -> MutationResult:
        ...         # Apply mutations...
        ...         return MutationResult(mutated_sequence=..., mutations=...)

    Example (legacy style, still supported):
        >>> class LegacyModel(MutationModel):
        ...     def apply_mutation(self, sequence_object):
        ...         # Apply mutations...
        ...         return mutated_seq, mutations_dict, mutation_rate
    """

    def __init__(
        self,
        min_mutation_rate: float = 0.0,
        max_mutation_rate: float = 0.0,
        productive: bool = False,
    ):
        """
        Initialize the mutation model.

        Args:
            min_mutation_rate: Minimum mutation rate (default: 0.0).
            max_mutation_rate: Maximum mutation rate (default: 0.0).
            productive: If True, avoid stop codons and preserve anchor residues.
        """
        self.min_mutation_rate = min_mutation_rate
        self.max_mutation_rate = max_mutation_rate
        self.productive = productive

    # =========================================================================
    # Abstract properties - subclasses must implement
    # =========================================================================

    @property
    @abstractmethod
    def name(self) -> str:
        """Human-readable name of this mutation model."""
        pass

    # =========================================================================
    # Core mutation method - subclasses should implement one of these
    # =========================================================================

    def apply(self, context: MutationContext) -> MutationResult:
        """
        Apply mutations to a sequence.

        This is the preferred method for applying mutations. It takes an
        explicit MutationContext and returns a typed MutationResult.

        Args:
            context: MutationContext containing the sequence and metadata.

        Returns:
            MutationResult with the mutated sequence and mutation details.

        Note:
            Subclasses can override this method directly, or implement
            `_apply_impl` for the actual mutation logic.
        """
        # Default implementation calls legacy method if available
        # Subclasses should override this with proper implementation
        raise NotImplementedError(
            f"{self.__class__.__name__} must implement either 'apply' or 'apply_mutation'"
        )

    def apply_mutation(self, sequence_object) -> Tuple[str, Dict[int, str], float]:
        """
        Apply mutations to a sequence object (legacy interface).

        This method is deprecated but maintained for backwards compatibility.
        New code should use the `apply` method instead.

        Args:
            sequence_object: A sequence object with required attributes
                           (ungapped_seq, v_seq_start, v_seq_end, etc.)

        Returns:
            Tuple of (mutated_sequence, mutations_dict, mutation_rate)
        """
        # Create context from sequence object
        context = MutationContext.from_sequence_object(sequence_object)

        # Call the new interface
        result = self.apply(context)

        # Return in legacy format
        return result.mutated_sequence, result.mutations, result.mutation_rate

    # =========================================================================
    # Metadata properties - subclasses can override these
    # =========================================================================

    @property
    def supports_productive_mode(self) -> bool:
        """
        Whether this model supports productive mode (avoiding stop codons).

        Returns:
            True if the model can ensure sequences remain productive.
        """
        return False

    @property
    def is_context_dependent(self) -> bool:
        """
        Whether mutation probabilities depend on surrounding sequence context.

        Returns:
            True for models like S5F that use k-mer context.
        """
        return False

    @property
    def requires_external_data(self) -> bool:
        """
        Whether this model requires external data files (e.g., mutation matrices).

        Returns:
            True if external data must be loaded before use.
        """
        return False

    # =========================================================================
    # Utility methods
    # =========================================================================

    def __repr__(self) -> str:
        """String representation of the model."""
        return (
            f"{self.__class__.__name__}("
            f"min_rate={self.min_mutation_rate}, "
            f"max_rate={self.max_mutation_rate}, "
            f"productive={self.productive})"
        )
