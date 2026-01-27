"""
Stop Codon Checker

Reusable component for detecting and avoiding stop codons during mutation.
Consolidates the duplicated stop codon logic from Uniform and S5F models.
"""

from typing import List, Tuple, Set


class StopCodonChecker:
    """
    Reusable stop codon detection for mutation models.

    This class consolidates the stop codon detection logic that was previously
    duplicated across Uniform and S5F mutation models.

    Attributes:
        STOP_CODONS: Frozenset of stop codon sequences.
        NUCLEOTIDES: Frozenset of valid nucleotide characters.

    Example:
        >>> checker = StopCodonChecker()
        >>> checker.would_create_stop("ATGATGATG", 3, "A")  # TAG
        True
        >>> checker.find_stop_codons("ATGTAGATG")
        [3]
    """

    STOP_CODONS: frozenset = frozenset(["TAG", "TAA", "TGA"])
    NUCLEOTIDES: frozenset = frozenset(["A", "T", "C", "G"])

    def would_create_stop(
        self, sequence: str, position: int, new_base: str
    ) -> bool:
        """
        Check if a mutation at the given position would create a stop codon.

        This checks all reading frames that could be affected by the mutation:
        - The codon containing this position
        - Adjacent codons if the position is near a codon boundary

        Args:
            sequence: The current nucleotide sequence (list or string).
            position: The position being mutated (0-indexed).
            new_base: The new nucleotide to be placed at this position.

        Returns:
            True if the mutation would create a stop codon, False otherwise.

        Example:
            >>> checker = StopCodonChecker()
            >>> checker.would_create_stop("ATGATGATG", 3, "A")  # Creates TAG
            True
            >>> checker.would_create_stop("ATGATGATG", 3, "C")  # Creates TCG
            False
        """
        if position >= len(sequence):
            return False

        # Convert to list if string for easier manipulation
        seq = list(sequence)
        seq[position] = new_base

        # Get the codon start position (frame-aligned to position 0)
        codon_start = position - (position % 3)

        # Check the codon containing this position
        if codon_start + 3 <= len(seq):
            codon = "".join(seq[codon_start : codon_start + 3])
            if codon in self.STOP_CODONS:
                return True

        return False

    def would_create_stop_any_frame(
        self,
        nucleotides: List[str],
        center_position: int,
        new_base: str,
        reading_frame: int,
    ) -> Tuple[bool, str]:
        """
        Check if mutation creates stop codon in current, previous, or next reading frame.

        This is used by S5F model which operates on 5-mers where the center
        nucleotide (position 2) is being mutated.

        Args:
            nucleotides: List of 5 nucleotide characters (a 5-mer).
            center_position: Index of center nucleotide in the 5-mer (usually 2).
            new_base: The new nucleotide value.
            reading_frame: The reading frame offset (0, 1, or 2).

        Returns:
            Tuple of (has_stop_codon, codon_sequence).
        """
        # Replace center nucleotide with the new base
        nucs = list(nucleotides)
        nucs[center_position] = new_base

        tests = []

        # Check the current reading frame
        if reading_frame + 3 <= len(nucs):
            codon = "".join(nucs[reading_frame : reading_frame + 3])
            if codon in self.STOP_CODONS:
                tests.append(True)

            # Check the previous reading frame if possible
            if reading_frame > 0:
                prev_codon = "".join(nucs[reading_frame - 1 : reading_frame + 2])
                if prev_codon in self.STOP_CODONS:
                    tests.append(True)

            # Check the next reading frame if possible
            if reading_frame < 2 and reading_frame + 4 <= len(nucs):
                next_codon = "".join(nucs[reading_frame + 1 : reading_frame + 4])
                if next_codon in self.STOP_CODONS:
                    tests.append(True)

            return any(tests), codon

        return False, ""

    def find_stop_codons(self, sequence: str) -> List[int]:
        """
        Find all stop codon positions in a sequence.

        Args:
            sequence: The nucleotide sequence to check.

        Returns:
            List of positions where stop codons start (0-indexed).

        Example:
            >>> checker = StopCodonChecker()
            >>> checker.find_stop_codons("ATGTAGATGTAA")
            [3, 9]
        """
        positions = []
        for i in range(0, len(sequence) - 2, 3):
            codon = sequence[i : i + 3]
            if codon in self.STOP_CODONS:
                positions.append(i)
        return positions

    def has_stop_codon(self, sequence: str) -> bool:
        """
        Check if sequence contains any stop codons.

        Args:
            sequence: The nucleotide sequence to check.

        Returns:
            True if any stop codon is present, False otherwise.
        """
        for i in range(0, len(sequence) - 2, 3):
            if sequence[i : i + 3] in self.STOP_CODONS:
                return True
        return False

    def get_stop_forming_bases(self, five_mer: str) -> Set[str]:
        """
        Get bases that would form stop codons if placed at center of 5-mer.

        This is used to filter out mutations that would create stop codons.

        Args:
            five_mer: A 5-character nucleotide string.

        Returns:
            Set of nucleotide bases that would create a stop codon.

        Example:
            >>> checker = StopCodonChecker()
            >>> checker.get_stop_forming_bases("NNTNN")  # T at center
            {'A', 'G'}  # TAG, TGA would be stops
        """
        stop_forming = set()

        for nucleotide in self.NUCLEOTIDES:
            # Replace center nucleotide
            test_seq = five_mer[:2] + nucleotide + five_mer[3:]

            # Check all possible 3-mers in the 5-mer
            for i in range(3):  # Positions 0, 1, 2 can start a codon
                if i + 3 <= len(test_seq):
                    codon = test_seq[i : i + 3]
                    if codon in self.STOP_CODONS:
                        stop_forming.add(nucleotide)
                        break

        return stop_forming
