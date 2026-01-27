"""
Mutation Context

Explicit input contract for mutation models. This dataclass contains all the
information a mutation model needs to apply mutations to a sequence.
"""

from dataclasses import dataclass
from typing import Optional, Tuple, List
from enum import Enum, auto


class ChainType(Enum):
    """Enumeration of supported immunoglobulin chain types."""
    BCR_HEAVY = auto()
    BCR_LIGHT_KAPPA = auto()
    BCR_LIGHT_LAMBDA = auto()
    TCR_ALPHA = auto()
    TCR_BETA = auto()
    TCR_GAMMA = auto()
    TCR_DELTA = auto()


@dataclass
class MutationContext:
    """
    Everything a mutation model needs to know about the sequence.

    This provides an explicit input contract, replacing the implicit requirement
    that sequence objects have specific attributes. Models receive this context
    rather than raw sequence objects, making the interface clear and testable.

    Attributes:
        sequence: The nucleotide sequence to mutate.
        v_region: Tuple of (start, end) positions for the V region.
        j_region: Tuple of (start, end) positions for the J region.
        d_region: Optional tuple of (start, end) for D region (None for light chains).
        junction_start: Start position of the CDR3/junction region.
        junction_end: End position of the CDR3/junction region.
        chain_type: The type of immunoglobulin chain.

    Example:
        >>> context = MutationContext(
        ...     sequence="ATCGATCGATCG...",
        ...     v_region=(0, 100),
        ...     j_region=(150, 200),
        ...     d_region=(110, 140),  # None for light chains
        ...     junction_start=95,
        ...     junction_end=155,
        ...     chain_type=ChainType.BCR_HEAVY
        ... )

    Note:
        For light chains (kappa/lambda), d_region should be None since they
        don't have a D segment.
    """

    sequence: str
    v_region: Tuple[int, int]
    j_region: Tuple[int, int]
    junction_start: int
    junction_end: int
    chain_type: ChainType
    d_region: Optional[Tuple[int, int]] = None

    @property
    def has_d_region(self) -> bool:
        """Check if this context has a D region (heavy chains only)."""
        return self.d_region is not None

    @property
    def is_heavy_chain(self) -> bool:
        """Check if this is a heavy chain context."""
        return self.chain_type in (ChainType.BCR_HEAVY, ChainType.TCR_BETA, ChainType.TCR_DELTA)

    @property
    def mutable_positions(self) -> list:
        """
        Get positions eligible for mutation (V, D, J regions only, not NP regions).

        Returns:
            List of positions that can be mutated.
        """
        positions = []
        # V region positions
        positions.extend(range(self.v_region[0], self.v_region[1]))
        # D region positions (if present)
        if self.d_region:
            positions.extend(range(self.d_region[0], self.d_region[1]))
        # J region positions
        positions.extend(range(self.j_region[0], self.j_region[1]))
        return positions

    @property
    def v_anchor_positions(self) -> Tuple[int, int, int]:
        """Get the 3 positions of the conserved V anchor (Cysteine codon)."""
        return (self.junction_start, self.junction_start + 1, self.junction_start + 2)

    @property
    def j_anchor_positions(self) -> Tuple[int, int, int]:
        """Get the 3 positions of the conserved J anchor (F/W codon)."""
        return (self.junction_end - 3, self.junction_end - 2, self.junction_end - 1)

    @classmethod
    def from_sequence_object(cls, sequence_obj) -> "MutationContext":
        """
        Create a MutationContext from a legacy sequence object.

        This factory method provides backwards compatibility by extracting
        the required information from existing sequence objects.

        Args:
            sequence_obj: A HeavyChainSequence or LightChainSequence object.

        Returns:
            MutationContext populated from the sequence object.
        """
        from ..sequence import HeavyChainSequence
        from ..sequence.light_chain import LightChainSequence

        # Determine chain type
        if isinstance(sequence_obj, HeavyChainSequence):
            chain_type = ChainType.BCR_HEAVY
            d_region = (sequence_obj.d_seq_start, sequence_obj.d_seq_end)
        elif isinstance(sequence_obj, LightChainSequence):
            # Check if kappa or lambda based on allele name
            if hasattr(sequence_obj, 'v_allele') and 'IGK' in str(sequence_obj.v_allele):
                chain_type = ChainType.BCR_LIGHT_KAPPA
            else:
                chain_type = ChainType.BCR_LIGHT_LAMBDA
            d_region = None
        else:
            # Default to heavy for unknown types
            chain_type = ChainType.BCR_HEAVY
            d_region = (
                getattr(sequence_obj, 'd_seq_start', 0),
                getattr(sequence_obj, 'd_seq_end', 0)
            )

        return cls(
            sequence=sequence_obj.ungapped_seq,
            v_region=(sequence_obj.v_seq_start, sequence_obj.v_seq_end),
            j_region=(sequence_obj.j_seq_start, sequence_obj.j_seq_end),
            d_region=d_region,
            junction_start=sequence_obj.junction_start,
            junction_end=sequence_obj.junction_end,
            chain_type=chain_type,
        )
