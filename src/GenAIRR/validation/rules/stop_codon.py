"""
Stop Codon Rule

Checks for the presence of stop codons (TAG, TAA, TGA) in the reading frame.
"""

from typing import Tuple
from ..base import ValidationRule, ValidationContext, ValidationRuleResult


class StopCodonRule(ValidationRule):
    """
    Check for stop codons (TAG, TAA, TGA) in the reading frame.

    Scans the sequence from position 0, checking every codon (3 nucleotides).
    A stop codon anywhere in the reading frame makes the sequence non-functional.

    Metadata:
        is_fatal: True - Cannot translate if stop codon present
        contributes_to_vj_in_frame: True - Stop codon means not in frame
    """

    STOP_CODONS = frozenset(["TAG", "TAA", "TGA"])

    @property
    def name(self) -> str:
        return "Stop Codon Check"

    @property
    def is_fatal(self) -> bool:
        """Stop codons prevent translation, so this is fatal."""
        return True

    @property
    def contributes_to_vj_in_frame(self) -> bool:
        """A stop codon means the sequence is not in frame."""
        return True

    def check(self, context: ValidationContext) -> ValidationRuleResult:
        has_stop, position = self._find_stop_codon(context.sequence)

        if has_stop:
            return ValidationRuleResult(
                passed=False,
                note="Stop codon present.",
                data={"stop_codon_position": position}
            )

        return ValidationRuleResult(
            passed=True,
            data={"stop_codon_position": -1}
        )

    def update_context(self, context: ValidationContext, result: ValidationRuleResult) -> None:
        """Store stop codon position in context."""
        context.stop_codon_position = result.data.get("stop_codon_position", -1)

    def _find_stop_codon(self, sequence: str) -> Tuple[bool, int]:
        """Find the first stop codon in the sequence."""
        for pos in range(0, len(sequence), 3):
            codon = sequence[pos:pos + 3]
            if codon in self.STOP_CODONS:
                return True, pos
        return False, -1

    def check_sequence(self, sequence: str) -> Tuple[bool, int]:
        """
        Public method to check a sequence for stop codons.

        This method is exposed for use by other components that need
        to check for stop codons without running the full validation.

        Args:
            sequence: Nucleotide sequence to check

        Returns:
            Tuple of (has_stop_codon, position)
        """
        return self._find_stop_codon(sequence)
