"""
Junction Translatable Rule

Checks that the junction can be translated to amino acids and performs the translation.
"""

from ..base import ValidationRule, ValidationContext, ValidationRuleResult
from ...utilities import translate


class JunctionTranslatableRule(ValidationRule):
    """
    Check that the junction can be translated to amino acids.

    This rule:
    1. Verifies junction length is divisible by 3
    2. Translates the junction sequence
    3. Stores the translation in the context for subsequent rules

    This rule should run before conserved residue checks, as they
    depend on the translated junction being available in the context.

    Metadata:
        is_fatal: True - If translation fails, subsequent rules can't check residues
    """

    @property
    def name(self) -> str:
        return "Junction Translation Check"

    @property
    def is_fatal(self) -> bool:
        """If we can't translate, conserved residue checks can't run."""
        return True

    def check(self, context: ValidationContext) -> ValidationRuleResult:
        if context.junction_length % 3 != 0:
            return ValidationRuleResult(
                passed=False,
                note="Junction length not divisible by 3."
            )

        # Translate the junction
        junction_aa = translate(context.junction_seq)

        return ValidationRuleResult(
            passed=True,
            data={"junction_aa": junction_aa}
        )

    def update_context(self, context: ValidationContext, result: ValidationRuleResult) -> None:
        """Store translated junction in context for subsequent rules."""
        if result.passed:
            context.junction_aa = result.data.get("junction_aa", "")
