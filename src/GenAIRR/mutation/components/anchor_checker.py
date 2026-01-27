"""
Anchor Checker

Reusable component for checking and preserving conserved CDR3 anchor residues
(V anchor Cysteine, J anchor Phenylalanine/Tryptophan) during mutation.
"""

from typing import Tuple, Set


class AnchorChecker:
    """
    Checks and preserves conserved CDR3 anchor residues during mutation.

    The CDR3 region is bounded by:
    - V anchor: Conserved Cysteine (C) encoded by TGC or TGT
    - J anchor: Conserved Phenylalanine (F) encoded by TTC/TTT, or
                Tryptophan (W) encoded by TGG

    This class consolidates the anchor preservation logic that was duplicated
    in both Uniform and S5F mutation models.

    Attributes:
        V_ANCHOR_CODONS: Codons encoding Cysteine (C).
        J_ANCHOR_CODONS: Codons encoding Phenylalanine (F) or Tryptophan (W).
    """

    # Cysteine codons (V anchor)
    V_ANCHOR_CODONS: frozenset = frozenset(["TGC", "TGT"])

    # Phenylalanine and Tryptophan codons (J anchor)
    J_ANCHOR_CODONS: frozenset = frozenset(["TTT", "TTC", "TGG"])

    def is_v_anchor_position(
        self, position: int, junction_start: int
    ) -> bool:
        """
        Check if a position is part of the V anchor codon.

        Args:
            position: The position being checked (0-indexed).
            junction_start: The start of the junction/CDR3 region.

        Returns:
            True if position is in the V anchor codon (3 nucleotides).
        """
        return junction_start <= position < junction_start + 3

    def is_j_anchor_position(
        self, position: int, junction_end: int
    ) -> bool:
        """
        Check if a position is part of the J anchor codon.

        Args:
            position: The position being checked (0-indexed).
            junction_end: The end of the junction/CDR3 region.

        Returns:
            True if position is in the J anchor codon (last 3 nucleotides).
        """
        return junction_end - 3 <= position < junction_end

    def is_anchor_position(
        self, position: int, junction_start: int, junction_end: int
    ) -> Tuple[bool, str]:
        """
        Check if a position is part of either anchor codon.

        Args:
            position: The position being checked (0-indexed).
            junction_start: The start of the junction/CDR3 region.
            junction_end: The end of the junction/CDR3 region.

        Returns:
            Tuple of (is_anchor, anchor_type) where anchor_type is 'v', 'j', or ''.
        """
        if self.is_v_anchor_position(position, junction_start):
            return True, "v"
        if self.is_j_anchor_position(position, junction_end):
            return True, "j"
        return False, ""

    def would_destroy_anchor(
        self,
        sequence: str,
        position: int,
        new_base: str,
        junction_start: int,
        junction_end: int,
    ) -> bool:
        """
        Check if a mutation would destroy a conserved anchor residue.

        Args:
            sequence: The current sequence (string or list).
            position: The position being mutated.
            new_base: The new nucleotide to place.
            junction_start: Start of junction/CDR3 region.
            junction_end: End of junction/CDR3 region.

        Returns:
            True if mutation would destroy the anchor, False otherwise.
        """
        is_anchor, anchor_type = self.is_anchor_position(
            position, junction_start, junction_end
        )

        if not is_anchor:
            return False

        # Build the mutated codon
        if anchor_type == "v":
            codon_start = junction_start
            valid_codons = self.V_ANCHOR_CODONS
        else:  # j
            codon_start = junction_end - 3
            valid_codons = self.J_ANCHOR_CODONS

        # Create the codon with the mutation applied
        seq = list(sequence)
        seq[position] = new_base
        codon = "".join(seq[codon_start : codon_start + 3])

        # Check if the resulting codon is still valid
        return codon not in valid_codons

    def get_restricted_positions(
        self, junction_start: int, junction_end: int
    ) -> dict:
        """
        Get a mapping of restricted anchor positions to their anchor type.

        This provides the same format used by the legacy Uniform and S5F models.

        Args:
            junction_start: Start of junction/CDR3 region.
            junction_end: End of junction/CDR3 region.

        Returns:
            Dict mapping position -> anchor_type ('v' or 'j').
        """
        return {
            junction_start: "v",
            junction_start + 1: "v",
            junction_start + 2: "v",
            junction_end - 3: "j",
            junction_end - 2: "j",
            junction_end - 1: "j",
        }

    def is_valid_anchor_mutation(
        self,
        codon: str,
        anchor_type: str,
    ) -> bool:
        """
        Check if a mutated codon is still a valid anchor.

        Args:
            codon: The 3-nucleotide codon string.
            anchor_type: Either 'v' or 'j'.

        Returns:
            True if the codon still encodes the correct anchor residue.
        """
        if anchor_type == "v":
            return codon in self.V_ANCHOR_CODONS
        elif anchor_type == "j":
            return codon in self.J_ANCHOR_CODONS
        return True
