"""
Conserved Cysteine Rule

Checks for the conserved Cysteine (C) at the junction start.
"""

from ..base import ValidationRule, ValidationContext, ValidationRuleResult


class ConservedCysteineRule(ValidationRule):
    """
    Check for conserved Cysteine (C) at the junction start.

    The second conserved cysteine of the V region should be at the
    start of the junction (CDR3). This is a hallmark of functional
    immunoglobulin sequences and is essential for proper antibody folding.

    Note: This rule handles the case where junction_aa is empty by returning
    failure - no need to skip via requires_translation.
    """

    REQUIRED_RESIDUE = "C"

    @property
    def name(self) -> str:
        return "Conserved Cysteine (C) Check"

    def check(self, context: ValidationContext) -> ValidationRuleResult:
        # should_skip handles the case where junction_aa is empty
        if not context.junction_aa:
            return ValidationRuleResult(
                passed=False,
                note="Cannot check: junction not translated."
            )

        if not context.junction_aa.startswith(self.REQUIRED_RESIDUE):
            return ValidationRuleResult(
                passed=False,
                note="V second C not present."
            )

        return ValidationRuleResult(passed=True)
