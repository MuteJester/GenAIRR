"""
Result dataclass for functionality assessment.
"""

from dataclasses import dataclass


@dataclass
class FunctionalityResult:
    """
    Complete result of functionality assessment.

    Attributes:
        functional: True if the sequence meets all functionality criteria
        stop_codon: True if a stop codon was found in the reading frame
        stop_codon_position: Position of the first stop codon (-1 if none)
        vj_in_frame: True if junction boundaries and length are properly aligned
        junction_aa: Amino acid translation of the junction (empty if not translatable)
        note: Description of why the sequence is not functional (empty if functional)
    """
    functional: bool
    stop_codon: bool
    stop_codon_position: int
    vj_in_frame: bool
    junction_aa: str
    note: str
