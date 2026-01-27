"""
Uniform Mutation Model

Applies mutations uniformly across the sequence at a random rate between
min_mutation_rate and max_mutation_rate.
"""

import random
from typing import List, Set

from .mutation_model import MutationModel
from .context import MutationContext
from .result import MutationResult
from .components import StopCodonChecker, AnchorChecker


class Uniform(MutationModel):
    """
    Uniform mutation model that applies mutations at a constant rate.

    This model selects positions randomly from the V, D, and J regions
    (excluding NP regions) and mutates them to a different nucleotide base.

    Args:
        min_mutation_rate: Minimum mutation rate (0.0 to 1.0). Default: 0.0
        max_mutation_rate: Maximum mutation rate (0.0 to 1.0). Default: 0.0
        productive: If True, avoid stop codons and preserve anchor residues.
                   Default: False

    Example:
        >>> from GenAIRR import Uniform, set_seed
        >>> set_seed(42)
        >>> model = Uniform(min_mutation_rate=0.01, max_mutation_rate=0.05)
        >>> # Use with pipeline or apply directly to a context

    Attributes:
        NUCLEOTIDES: Set of valid nucleotide bases.
    """

    NUCLEOTIDES: Set[str] = {"A", "T", "C", "G"}

    def __init__(
        self,
        min_mutation_rate: float = 0.0,
        max_mutation_rate: float = 0.0,
        productive: bool = False,
    ):
        """Initialize the Uniform mutation model."""
        super().__init__(
            min_mutation_rate=min_mutation_rate,
            max_mutation_rate=max_mutation_rate,
            productive=productive,
        )
        self.bases = self.NUCLEOTIDES  # Legacy attribute for compatibility
        # Shared components
        self._stop_checker = StopCodonChecker()
        self._anchor_checker = AnchorChecker()

    @property
    def name(self) -> str:
        """Human-readable name of this mutation model."""
        return "Uniform"

    @property
    def supports_productive_mode(self) -> bool:
        """This model supports productive mode."""
        return True

    @property
    def is_context_dependent(self) -> bool:
        """Uniform model does not use sequence context."""
        return False

    def apply(self, context: MutationContext) -> MutationResult:
        """
        Apply uniform mutations to a sequence.

        Args:
            context: MutationContext containing sequence and region information.

        Returns:
            MutationResult with the mutated sequence and mutation details.
        """
        sequence = context.sequence
        mutation_rate = random.uniform(self.min_mutation_rate, self.max_mutation_rate)
        num_mutations = int(mutation_rate * len(sequence))

        # Get mutable positions (V, D, J regions only)
        mutable_positions = context.mutable_positions

        if num_mutations == 0 or not mutable_positions:
            return MutationResult(
                mutated_sequence=sequence,
                mutations={},
                mutation_rate=mutation_rate,
                mutation_count=0,
            )

        # Select positions to mutate
        positions_to_mutate = random.sample(
            mutable_positions, k=min(num_mutations, len(mutable_positions))
        )

        # Apply mutations
        mutations = {}
        mutated_seq = list(sequence)

        if self.productive:
            mutations, mutated_seq = self._apply_productive_mutations(
                context, positions_to_mutate, mutated_seq, mutable_positions
            )
        else:
            for position in positions_to_mutate:
                new_base = self._mutate_base(mutated_seq[position])
                mutations[position] = f"{mutated_seq[position]}>{new_base}"
                mutated_seq[position] = new_base

        return MutationResult(
            mutated_sequence="".join(mutated_seq),
            mutations=mutations,
            mutation_rate=mutation_rate,
            mutation_count=len(mutations),
        )

    def _apply_productive_mutations(
        self,
        context: MutationContext,
        positions: List[int],
        mutated_seq: List[str],
        mutable_positions: List[int],
    ) -> tuple:
        """
        Apply mutations while ensuring the sequence remains productive.

        This avoids:
        - Creating stop codons
        - Destroying conserved anchor residues (Cysteine at V, F/W at J)
        """
        mutations = {}
        restricted_positions = self._anchor_checker.get_restricted_positions(
            context.junction_start, context.junction_end
        )

        for position in positions:
            position, new_base = self._productive_mutation(
                position,
                mutated_seq,
                restricted_positions,
                mutable_positions,
                max_attempts=100,
            )
            if new_base:
                mutations[position] = f"{mutated_seq[position]}>{new_base}"
                mutated_seq[position] = new_base

        return mutations, mutated_seq

    def _productive_mutation(
        self,
        position: int,
        mutated_seq: List[str],
        restricted_positions: dict,
        mutable_positions: List[int],
        max_attempts: int = 100,
    ) -> tuple:
        """
        Find a valid mutation that doesn't create stop codons or destroy anchors.

        Args:
            position: Initial position to try mutating.
            mutated_seq: Current sequence as list of characters.
            restricted_positions: Dict of anchor positions.
            mutable_positions: List of positions that can be mutated.
            max_attempts: Maximum number of attempts before giving up.

        Returns:
            Tuple of (final_position, new_base) or (position, None) if failed.
        """
        for _ in range(max_attempts):
            new_base = self._mutate_base(mutated_seq[position])

            # Check for stop codon
            if self._stop_checker.would_create_stop(
                mutated_seq, position, new_base
            ):
                position = random.choice(mutable_positions)
                continue

            # Check for anchor destruction
            if position in restricted_positions:
                anchor_type = restricted_positions[position]
                codon = self._get_mutated_codon(
                    mutated_seq, position, new_base, restricted_positions, anchor_type
                )
                if not self._anchor_checker.is_valid_anchor_mutation(codon, anchor_type):
                    position = random.choice(mutable_positions)
                    continue

            return position, new_base

        # Failed to find valid mutation
        return position, None

    def _get_mutated_codon(
        self,
        sequence: List[str],
        position: int,
        new_base: str,
        restricted_positions: dict,
        anchor_type: str,
    ) -> str:
        """Build the codon that would result from this mutation."""
        codon_chars = []
        for pos, tag in restricted_positions.items():
            if tag == anchor_type:
                if pos == position:
                    codon_chars.append((pos, new_base))
                else:
                    codon_chars.append((pos, sequence[pos]))
        # Sort by position and extract characters
        codon_chars.sort(key=lambda x: x[0])
        return "".join(char for _, char in codon_chars)

    def _mutate_base(self, base: str) -> str:
        """
        Mutate a nucleotide to a different base.

        Args:
            base: The original nucleotide.

        Returns:
            A different nucleotide base.
        """
        return random.choice(list(self.NUCLEOTIDES - {base}))

    # =========================================================================
    # Legacy interface (for backwards compatibility)
    # =========================================================================

    def mutable_positions(self, sequence_object) -> List[int]:
        """
        Get mutable positions from a sequence object (legacy interface).

        Deprecated: Use MutationContext.mutable_positions instead.
        """
        positions = []
        positions.extend(range(sequence_object.v_seq_start, sequence_object.v_seq_end))
        if hasattr(sequence_object, "d_seq_start"):
            positions.extend(
                range(sequence_object.d_seq_start, sequence_object.d_seq_end)
            )
        positions.extend(range(sequence_object.j_seq_start, sequence_object.j_seq_end))
        return positions

    def apply_mutation(self, sequence_object):
        """
        Apply mutations to a sequence object (legacy interface).

        This method is maintained for backwards compatibility.
        New code should use the `apply` method with MutationContext instead.

        Args:
            sequence_object: A sequence object with required attributes.

        Returns:
            Tuple of (mutated_sequence, mutations_dict, mutation_rate).
        """
        # Create context from sequence object
        context = MutationContext.from_sequence_object(sequence_object)

        # Call the new interface
        result = self.apply(context)

        # Return in legacy format
        return result.mutated_sequence, result.mutations, result.mutation_rate
