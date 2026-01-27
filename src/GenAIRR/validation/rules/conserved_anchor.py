"""
Conserved Anchor Rule

Checks for the conserved Phenylalanine (F) or Tryptophan (W) at junction end.
"""

from ..base import ValidationRule, ValidationContext, ValidationRuleResult


class ConservedAnchorRule(ValidationRule):
    """
    Check for conserved Phenylalanine (F) or Tryptophan (W) at junction end.

    The J region anchor residue (F or W) should be at the end of the
    junction. This conserved motif is required for functional antibodies
    and marks the boundary of the CDR3 region.

    Note: This rule handles the case where junction_aa is empty by returning
    failure - no need to skip via requires_translation.
    """

    ALLOWED_RESIDUES = frozenset(["F", "W"])

    @property
    def name(self) -> str:
        return "Conserved J Anchor (F/W) Check"

    def check(self, context: ValidationContext) -> ValidationRuleResult:
        # should_skip handles the case where junction_aa is empty
        if not context.junction_aa:
            return ValidationRuleResult(
                passed=False,
                note="Cannot check: junction not translated."
            )

        if context.junction_aa[-1] not in self.ALLOWED_RESIDUES:
            return ValidationRuleResult(
                passed=False,
                note="J anchor (W/F) not present."
            )

        return ValidationRuleResult(passed=True)
