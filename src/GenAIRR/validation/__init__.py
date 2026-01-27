"""
Validation Package

This package provides a modular, composition-based framework for validating
immunoglobulin sequence functionality.

Structure:
    - base.py: Base classes (ValidationRule, ValidationContext, ValidationRuleResult)
    - result.py: FunctionalityResult dataclass
    - validator.py: Main FunctionalityValidator class
    - rules/: Individual validation rules
        - stop_codon.py: StopCodonRule
        - frame_alignment.py: FrameAlignmentRule
        - junction_translatable.py: JunctionTranslatableRule
        - conserved_cysteine.py: ConservedCysteineRule
        - conserved_anchor.py: ConservedAnchorRule

Usage:
    from GenAIRR.validation import FunctionalityValidator

    validator = FunctionalityValidator()
    result = validator.assess(sequence, junction_start, junction_end)

    if result.functional:
        print("Sequence is productive")
    else:
        print(f"Not productive: {result.note}")

Custom Rules:
    from GenAIRR.validation import ValidationRule, ValidationRuleResult

    class MyCustomRule(ValidationRule):
        @property
        def name(self) -> str:
            return "My Custom Check"

        def check(self, context):
            # Your validation logic here
            return ValidationRuleResult(passed=True)

    validator = FunctionalityValidator()
    validator.add_rule(MyCustomRule())
"""

# Base classes
from .base import (
    ValidationRule,
    ValidationContext,
    ValidationRuleResult,
)

# Result dataclass
from .result import FunctionalityResult

# Main validator
from .validator import FunctionalityValidator

# Individual rules (for advanced usage)
from .rules import (
    StopCodonRule,
    FrameAlignmentRule,
    JunctionTranslatableRule,
    ConservedCysteineRule,
    ConservedAnchorRule,
)

__all__ = [
    # Base classes
    "ValidationRule",
    "ValidationContext",
    "ValidationRuleResult",
    # Result
    "FunctionalityResult",
    # Validator
    "FunctionalityValidator",
    # Rules
    "StopCodonRule",
    "FrameAlignmentRule",
    "JunctionTranslatableRule",
    "ConservedCysteineRule",
    "ConservedAnchorRule",
]
