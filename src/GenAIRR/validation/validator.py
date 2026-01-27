"""
Functionality Validator

The main validator class that orchestrates validation rules using
metadata-driven execution - no type checking required.
"""

from typing import List, Tuple, Optional

from .base import ValidationRule, ValidationContext
from .result import FunctionalityResult
from .rules import (
    StopCodonRule,
    FrameAlignmentRule,
    JunctionTranslatableRule,
    ConservedCysteineRule,
    ConservedAnchorRule,
)


class FunctionalityValidator:
    """
    Validates whether an immunoglobulin sequence is functional (productive).

    This validator uses a metadata-driven approach where each validation
    rule declares its own behavior through properties. The validator
    executes rules generically without any type checking.

    Rule Metadata:
        - is_fatal: If True, failing stops further validation
        - requires_translation: If True, skipped if junction_aa is empty
        - contributes_to_vj_in_frame: If True, affects vj_in_frame calculation

    Default Rules (in order):
        1. StopCodonRule - No stop codons (fatal, contributes to vj_in_frame)
        2. FrameAlignmentRule - Proper frame alignment (contributes to vj_in_frame)
        3. JunctionTranslatableRule - Can translate junction (fatal)
        4. ConservedCysteineRule - Starts with C (requires translation)
        5. ConservedAnchorRule - Ends with F/W (requires translation)

    Customization:
        - Add rules: validator.add_rule(MyCustomRule())
        - Remove rules: validator.remove_rule("Rule Name")
        - Replace all rules: validator.rules = [Rule1(), Rule2()]

    Example:
        validator = FunctionalityValidator()
        result = validator.assess(sequence, junction_start, junction_end)

    Example (custom rule):
        class MyRule(ValidationRule):
            @property
            def name(self): return "My Check"

            @property
            def is_fatal(self): return True  # Stop if this fails

            def check(self, context): ...

        validator.add_rule(MyRule())
    """

    def __init__(self):
        """Initialize with default validation rules."""
        self.rules: List[ValidationRule] = [
            StopCodonRule(),
            FrameAlignmentRule(),
            JunctionTranslatableRule(),
            ConservedCysteineRule(),
            ConservedAnchorRule(),
        ]

    def add_rule(self, rule: ValidationRule, position: Optional[int] = None) -> None:
        """
        Add a validation rule.

        Args:
            rule: The ValidationRule to add
            position: Optional index to insert at (appends if None)
        """
        if position is None:
            self.rules.append(rule)
        else:
            self.rules.insert(position, rule)

    def remove_rule(self, rule_name: str) -> bool:
        """
        Remove a validation rule by name.

        Args:
            rule_name: The name of the rule to remove

        Returns:
            True if a rule was removed, False if not found
        """
        for i, rule in enumerate(self.rules):
            if rule.name == rule_name:
                self.rules.pop(i)
                return True
        return False

    def get_rule_names(self) -> List[str]:
        """Get list of all current rule names in order."""
        return [rule.name for rule in self.rules]

    def check_stop_codons(self, sequence: str) -> Tuple[bool, int]:
        """
        Check for stop codons in a sequence.

        Convenience method for checking stop codons without full validation.
        Finds the StopCodonRule and uses it directly.

        Args:
            sequence: Nucleotide sequence to check

        Returns:
            Tuple of (has_stop_codon, position)
        """
        for rule in self.rules:
            if hasattr(rule, 'check_sequence'):
                return rule.check_sequence(sequence)
        # Fallback if no stop codon rule found
        return False, -1

    def assess(self, sequence: str, junction_start: int, junction_end: int) -> FunctionalityResult:
        """
        Perform complete functionality assessment of a sequence.

        Runs all validation rules in order, respecting their metadata:
        - Skips rules where should_skip() returns True
        - Stops on fatal rules that fail
        - Calculates vj_in_frame from rules with contributes_to_vj_in_frame=True

        Args:
            sequence: The full nucleotide sequence
            junction_start: Start position of the junction
            junction_end: End position of the junction

        Returns:
            FunctionalityResult containing all assessment details
        """
        # Create context for rules
        context = ValidationContext(
            sequence=sequence,
            junction_start=junction_start,
            junction_end=junction_end
        )

        # Track state
        all_passed = True
        first_failure_note = ""
        vj_in_frame_passed = True
        stop_codon_found = False

        # Run each rule
        for rule in self.rules:
            # Check if rule should be skipped
            if rule.should_skip(context):
                continue

            # Execute the rule
            result = rule.check(context)

            # Let rule update context (e.g., store computed values)
            rule.update_context(context, result)

            # Track failures
            if not result.passed:
                all_passed = False

                # Store first failure note
                if not first_failure_note:
                    first_failure_note = result.note

                # Track vj_in_frame contributors
                if rule.contributes_to_vj_in_frame:
                    vj_in_frame_passed = False

                # Check for stop codon (for result field)
                if "stop_codon_position" in result.data:
                    stop_codon_found = True

                # Stop if this rule is fatal
                if rule.is_fatal:
                    break

        # Build result
        return FunctionalityResult(
            functional=all_passed,
            stop_codon=stop_codon_found,
            stop_codon_position=context.stop_codon_position,
            vj_in_frame=vj_in_frame_passed,
            junction_aa=context.junction_aa,
            note=first_failure_note
        )
