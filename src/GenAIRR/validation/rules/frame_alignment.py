"""
Frame Alignment Rule

Checks that the junction is properly frame-aligned for translation.
"""

from ..base import ValidationRule, ValidationContext, ValidationRuleResult


class FrameAlignmentRule(ValidationRule):
    """
    Check that the junction is properly frame-aligned.

    For a junction to be in-frame:
    - Junction start position must be divisible by 3
    - Junction end position must be divisible by 3
    - Junction length must be divisible by 3

    This ensures the junction can be translated into a valid amino acid sequence
    without frameshift issues.

    Metadata:
        contributes_to_vj_in_frame: True - Frame alignment directly affects vj_in_frame
    """

    @property
    def name(self) -> str:
        return "Frame Alignment Check"

    @property
    def contributes_to_vj_in_frame(self) -> bool:
        """Frame alignment is a core component of vj_in_frame."""
        return True

    def check(self, context: ValidationContext) -> ValidationRuleResult:
        start_aligned = context.junction_start % 3 == 0
        end_aligned = context.junction_end % 3 == 0
        length_aligned = context.junction_length % 3 == 0

        all_aligned = start_aligned and end_aligned and length_aligned

        if not all_aligned:
            reasons = []
            if not start_aligned:
                reasons.append("junction start not aligned")
            if not end_aligned:
                reasons.append("junction end not aligned")
            if not length_aligned:
                reasons.append("junction length not divisible by 3")

            return ValidationRuleResult(
                passed=False,
                note=f"Frame misalignment: {', '.join(reasons)}."
            )

        return ValidationRuleResult(passed=True)
