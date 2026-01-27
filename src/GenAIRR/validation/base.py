"""
Base classes for the validation framework.

This module provides the foundational classes for building validation rules:
- ValidationContext: Data container passed to each rule
- ValidationRuleResult: Result from a single rule check
- ValidationRule: Abstract base class for all rules with metadata support
"""

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Dict, Any, Set


@dataclass
class ValidationContext:
    """
    Context object containing all data needed for validation rules.

    This object is passed to each validation rule, providing a consistent
    interface for accessing sequence data. Rules can also store computed
    values here for use by subsequent rules.

    Attributes:
        sequence: The full nucleotide sequence
        junction_start: Start position of the junction
        junction_end: End position of the junction
        junction_length: Length of the junction (calculated)
        junction_seq: The junction nucleotide sequence (extracted)
        junction_aa: The junction amino acid sequence (set by translation rule)
        stop_codon_position: Position of stop codon if found (set by stop codon rule)
    """
    sequence: str
    junction_start: int
    junction_end: int
    junction_length: int = field(init=False)
    junction_seq: str = field(init=False)
    junction_aa: str = ""
    stop_codon_position: int = -1

    def __post_init__(self):
        self.junction_length = self.junction_end - self.junction_start
        self.junction_seq = self.sequence[self.junction_start:self.junction_end].upper()


@dataclass
class ValidationRuleResult:
    """
    Result from a single validation rule check.

    Attributes:
        passed: True if the rule passed
        note: Description of why it failed (empty if passed)
        data: Optional additional data from the check (e.g., stop codon position)
    """
    passed: bool
    note: str = ""
    data: Dict[str, Any] = field(default_factory=dict)


class ValidationRule(ABC):
    """
    Abstract base class for validation rules with metadata support.

    Each rule implements a single, focused validation check. Rules declare
    their behavior through metadata properties, allowing the validator to
    handle them generically without type checking.

    Metadata Properties:
        name: Human-readable name of the rule
        is_fatal: If True, failing this rule stops further validation
        requires_translation: If True, only run if junction_aa is available
        contributes_to_vj_in_frame: If True, this rule's pass/fail affects vj_in_frame

    To create a custom rule:
        1. Subclass ValidationRule
        2. Implement the `name` property
        3. Implement the `check` method
        4. Override metadata properties as needed

    Example:
        class MyCustomRule(ValidationRule):
            @property
            def name(self) -> str:
                return "My Custom Check"

            @property
            def is_fatal(self) -> bool:
                return True  # Stop validation if this fails

            def check(self, context: ValidationContext) -> ValidationRuleResult:
                if some_condition:
                    return ValidationRuleResult(passed=True)
                return ValidationRuleResult(passed=False, note="Failed because...")
    """

    # =========================================================================
    # Abstract methods - must be implemented
    # =========================================================================

    @property
    @abstractmethod
    def name(self) -> str:
        """Human-readable name of this validation rule."""
        pass

    @abstractmethod
    def check(self, context: ValidationContext) -> ValidationRuleResult:
        """
        Execute this validation rule.

        Args:
            context: ValidationContext containing sequence data

        Returns:
            ValidationRuleResult indicating pass/fail and any notes
        """
        pass

    # =========================================================================
    # Metadata properties - override as needed
    # =========================================================================

    @property
    def is_fatal(self) -> bool:
        """
        If True, failing this rule stops further validation.

        Use for rules where subsequent rules cannot meaningfully run
        if this one fails (e.g., stop codon check before translation).

        Default: False
        """
        return False

    @property
    def requires_translation(self) -> bool:
        """
        If True, this rule only runs if junction_aa is populated.

        Use for rules that need the translated junction sequence
        (e.g., conserved residue checks).

        Default: False
        """
        return False

    @property
    def contributes_to_vj_in_frame(self) -> bool:
        """
        If True, this rule's result affects the vj_in_frame calculation.

        vj_in_frame is True only if ALL rules with this flag pass.

        Default: False
        """
        return False

    # =========================================================================
    # Optional hooks - override for special behavior
    # =========================================================================

    def should_skip(self, context: ValidationContext) -> bool:
        """
        Determine if this rule should be skipped based on context.

        Override for custom skip logic beyond requires_translation.

        Args:
            context: Current validation context

        Returns:
            True if this rule should be skipped
        """
        if self.requires_translation and not context.junction_aa:
            return True
        return False

    def update_context(self, context: ValidationContext, result: ValidationRuleResult) -> None:
        """
        Update the context with data computed by this rule.

        Override to store computed values for use by subsequent rules.

        Args:
            context: Validation context to update
            result: The result from this rule's check
        """
        pass
