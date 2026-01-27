"""
Tests for FunctionalityValidator

This test module verifies that the FunctionalityValidator produces
identical results to the current implementation in HeavyChainSequence
and LightChainSequence classes.
"""

import unittest
from GenAIRR.validation import (
    FunctionalityValidator,
    FunctionalityResult,
    ValidationContext,
    StopCodonRule,
    FrameAlignmentRule,
    ConservedCysteineRule,
    ConservedAnchorRule,
)
from GenAIRR.sequence import HeavyChainSequence, LightChainSequence
from GenAIRR.data import HUMAN_IGH_OGRDB, HUMAN_IGK_OGRDB, HUMAN_IGL_OGRDB
from GenAIRR.utilities import translate


class TestFunctionalityValidatorUnit(unittest.TestCase):
    """Unit tests for FunctionalityValidator and its rules."""

    def setUp(self):
        self.validator = FunctionalityValidator()
        # Individual rules for direct testing
        self.stop_codon_rule = StopCodonRule()
        self.frame_rule = FrameAlignmentRule()
        self.cysteine_rule = ConservedCysteineRule()
        self.anchor_rule = ConservedAnchorRule()

    # =========================================================================
    # check_stop_codons tests
    # =========================================================================

    def test_check_stop_codons_no_stops(self):
        """Test sequence with no stop codons."""
        seq = "ATGAAACCCGGG"  # MKP G - no stops
        has_stop, pos = self.validator.check_stop_codons(seq)
        self.assertFalse(has_stop)
        self.assertEqual(pos, -1)

    def test_check_stop_codons_tag(self):
        """Test sequence with TAG stop codon."""
        seq = "ATGTAGCCC"  # M * P
        has_stop, pos = self.validator.check_stop_codons(seq)
        self.assertTrue(has_stop)
        self.assertEqual(pos, 3)

    def test_check_stop_codons_taa(self):
        """Test sequence with TAA stop codon."""
        seq = "ATGTAACCC"  # M * P
        has_stop, pos = self.validator.check_stop_codons(seq)
        self.assertTrue(has_stop)
        self.assertEqual(pos, 3)

    def test_check_stop_codons_tga(self):
        """Test sequence with TGA stop codon."""
        seq = "ATGTGACCC"  # M * P
        has_stop, pos = self.validator.check_stop_codons(seq)
        self.assertTrue(has_stop)
        self.assertEqual(pos, 3)

    def test_check_stop_codons_at_start(self):
        """Test stop codon at position 0."""
        seq = "TAGATGCCC"
        has_stop, pos = self.validator.check_stop_codons(seq)
        self.assertTrue(has_stop)
        self.assertEqual(pos, 0)

    def test_check_stop_codons_at_end(self):
        """Test stop codon at the end."""
        seq = "ATGCCCTAG"
        has_stop, pos = self.validator.check_stop_codons(seq)
        self.assertTrue(has_stop)
        self.assertEqual(pos, 6)

    def test_check_stop_codons_out_of_frame_ignored(self):
        """Test that stop codons out of frame are correctly ignored."""
        # TAG appears but not at codon boundary (positions 1-3)
        seq = "ATAGAACCC"  # Codons: ATA, GAA, CCC - no stops
        has_stop, pos = self.validator.check_stop_codons(seq)
        self.assertFalse(has_stop)
        self.assertEqual(pos, -1)

    # =========================================================================
    # FrameAlignmentRule tests
    # =========================================================================

    def test_check_frame_alignment_valid(self):
        """Test properly aligned junction."""
        # All divisible by 3
        ctx1 = ValidationContext(sequence="A" * 33, junction_start=0, junction_end=33)
        ctx2 = ValidationContext(sequence="A" * 351, junction_start=288, junction_end=351)
        ctx3 = ValidationContext(sequence="A" * 48, junction_start=9, junction_end=48)

        self.assertTrue(self.frame_rule.check(ctx1).passed)
        self.assertTrue(self.frame_rule.check(ctx2).passed)
        self.assertTrue(self.frame_rule.check(ctx3).passed)

    def test_check_frame_alignment_start_misaligned(self):
        """Test junction with misaligned start."""
        ctx1 = ValidationContext(sequence="A" * 33, junction_start=1, junction_end=33)
        ctx2 = ValidationContext(sequence="A" * 33, junction_start=2, junction_end=33)

        self.assertFalse(self.frame_rule.check(ctx1).passed)
        self.assertFalse(self.frame_rule.check(ctx2).passed)

    def test_check_frame_alignment_end_misaligned(self):
        """Test junction with misaligned end."""
        ctx1 = ValidationContext(sequence="A" * 32, junction_start=0, junction_end=32)
        ctx2 = ValidationContext(sequence="A" * 34, junction_start=0, junction_end=34)

        self.assertFalse(self.frame_rule.check(ctx1).passed)
        self.assertFalse(self.frame_rule.check(ctx2).passed)

    def test_check_frame_alignment_length_misaligned(self):
        """Test junction with length not divisible by 3."""
        ctx1 = ValidationContext(sequence="A" * 31, junction_start=0, junction_end=31)
        ctx2 = ValidationContext(sequence="A" * 32, junction_start=0, junction_end=32)

        self.assertFalse(self.frame_rule.check(ctx1).passed)
        self.assertFalse(self.frame_rule.check(ctx2).passed)

    # =========================================================================
    # ConservedCysteineRule and ConservedAnchorRule tests
    # =========================================================================

    def _make_context_with_junction_aa(self, junction_aa):
        """Helper to create a context with pre-set junction_aa."""
        ctx = ValidationContext(sequence="AAA", junction_start=0, junction_end=3)
        ctx.junction_aa = junction_aa
        return ctx

    def test_check_conserved_residues_valid_f(self):
        """Test valid junction ending with F."""
        ctx = self._make_context_with_junction_aa("CARWDYSSSGYYF")
        self.assertTrue(self.cysteine_rule.check(ctx).passed)
        self.assertTrue(self.anchor_rule.check(ctx).passed)

    def test_check_conserved_residues_valid_w(self):
        """Test valid junction ending with W."""
        ctx = self._make_context_with_junction_aa("CARWDYSSSGYW")
        self.assertTrue(self.cysteine_rule.check(ctx).passed)
        self.assertTrue(self.anchor_rule.check(ctx).passed)

    def test_check_conserved_residues_missing_c(self):
        """Test junction without starting C."""
        ctx = self._make_context_with_junction_aa("AARWDYSSSGYYF")
        result = self.cysteine_rule.check(ctx)
        self.assertFalse(result.passed)
        self.assertEqual(result.note, "V second C not present.")

    def test_check_conserved_residues_missing_fw(self):
        """Test junction without ending F or W."""
        ctx = self._make_context_with_junction_aa("CARWDYSSSGYYG")
        result = self.anchor_rule.check(ctx)
        self.assertFalse(result.passed)
        self.assertEqual(result.note, "J anchor (W/F) not present.")

    def test_check_conserved_residues_empty(self):
        """Test empty junction."""
        ctx = self._make_context_with_junction_aa("")
        result = self.cysteine_rule.check(ctx)
        self.assertFalse(result.passed)
        self.assertEqual(result.note, "Cannot check: junction not translated.")

    # =========================================================================
    # assess (full assessment) tests
    # =========================================================================

    def test_assess_functional_sequence(self):
        """Test a fully functional sequence."""
        # Create a synthetic functional sequence
        # Junction at position 0-33 (length 33, divisible by 3)
        # Junction: TGT...TTC (C...F)
        junction = "TGTGCCAGGTGGGACTACTACTACTACTTTTTC"  # C A R W D Y Y Y Y F F
        sequence = junction + "ATGATGATG"  # Add some extra sequence

        result = self.validator.assess(sequence, 0, 33)

        self.assertTrue(result.functional)
        self.assertFalse(result.stop_codon)
        self.assertEqual(result.stop_codon_position, -1)
        self.assertTrue(result.vj_in_frame)
        self.assertTrue(result.junction_aa.startswith("C"))
        self.assertTrue(result.junction_aa.endswith("F"))
        self.assertEqual(result.note, "")

    def test_assess_stop_codon(self):
        """Test sequence with stop codon."""
        sequence = "ATGTAGCCCGGGAAATTT"  # Has TAG at position 3
        result = self.validator.assess(sequence, 0, 18)

        self.assertFalse(result.functional)
        self.assertTrue(result.stop_codon)
        self.assertEqual(result.stop_codon_position, 3)
        self.assertFalse(result.vj_in_frame)
        self.assertEqual(result.note, "Stop codon present.")

    def test_assess_frame_misalignment(self):
        """Test sequence with frame misalignment."""
        # Junction length not divisible by 3
        sequence = "TGTGCCAGGTGGGACTACTAC"  # 21 bases - C A R W D Y Y (7 AAs)
        result = self.validator.assess(sequence, 0, 20)  # Length 20, not divisible by 3

        self.assertFalse(result.functional)
        self.assertFalse(result.stop_codon)
        self.assertFalse(result.vj_in_frame)  # Length not divisible by 3


class TestFunctionalityValidatorVsHeavyChain(unittest.TestCase):
    """
    Compare FunctionalityValidator results against HeavyChainSequence.

    These tests ensure the validator produces identical results to the
    current implementation embedded in HeavyChainSequence.
    """

    @classmethod
    def setUpClass(cls):
        cls.config = HUMAN_IGH_OGRDB
        cls.validator = FunctionalityValidator()

    def test_heavy_chain_comparison_multiple_sequences(self):
        """
        Generate multiple heavy chain sequences and verify validator
        produces identical results to the embedded check_functionality.
        """
        mismatches = []

        for i in range(100):
            # Generate a random heavy chain sequence
            seq = HeavyChainSequence.create_random(self.config)

            # Get results from current implementation
            current_functional = seq.functional
            current_stop_codon = seq.stop_codon
            current_vj_in_frame = seq.vj_in_frame
            current_note = seq.note
            current_junction_aa = getattr(seq, 'junction_aa', '')

            # Get results from validator
            result = self.validator.assess(
                seq.ungapped_seq,
                seq.junction_start,
                seq.junction_end
            )

            # Compare all fields
            if (result.functional != current_functional or
                result.stop_codon != current_stop_codon or
                result.vj_in_frame != current_vj_in_frame or
                result.junction_aa != current_junction_aa):

                mismatches.append({
                    'index': i,
                    'validator': {
                        'functional': result.functional,
                        'stop_codon': result.stop_codon,
                        'vj_in_frame': result.vj_in_frame,
                        'junction_aa': result.junction_aa,
                        'note': result.note
                    },
                    'current': {
                        'functional': current_functional,
                        'stop_codon': current_stop_codon,
                        'vj_in_frame': current_vj_in_frame,
                        'junction_aa': current_junction_aa,
                        'note': current_note
                    }
                })

        self.assertEqual(len(mismatches), 0,
                        f"Found {len(mismatches)} mismatches: {mismatches[:3]}")

    def test_heavy_chain_stop_codon_position(self):
        """Verify stop codon position matches between implementations."""
        for _ in range(50):
            seq = HeavyChainSequence.create_random(self.config)

            # Current implementation stores position during check
            current_has_stop, current_pos = seq.check_stops(seq.ungapped_seq, return_pos=True)

            # Validator result
            result = self.validator.assess(
                seq.ungapped_seq,
                seq.junction_start,
                seq.junction_end
            )

            self.assertEqual(result.stop_codon, current_has_stop)
            self.assertEqual(result.stop_codon_position, current_pos)


class TestFunctionalityValidatorVsLightChain(unittest.TestCase):
    """
    Compare FunctionalityValidator results against LightChainSequence.
    """

    @classmethod
    def setUpClass(cls):
        cls.kappa_config = HUMAN_IGK_OGRDB
        cls.lambda_config = HUMAN_IGL_OGRDB
        cls.validator = FunctionalityValidator()

    def test_light_chain_kappa_comparison(self):
        """Compare validator against kappa light chain sequences."""
        mismatches = []

        for i in range(50):
            seq = LightChainSequence.create_random(self.kappa_config)

            current_functional = seq.functional
            current_stop_codon = seq.stop_codon
            current_vj_in_frame = seq.vj_in_frame
            current_junction_aa = getattr(seq, 'junction_aa', '')

            result = self.validator.assess(
                seq.ungapped_seq,
                seq.junction_start,
                seq.junction_end
            )

            if (result.functional != current_functional or
                result.stop_codon != current_stop_codon or
                result.vj_in_frame != current_vj_in_frame or
                result.junction_aa != current_junction_aa):

                mismatches.append({
                    'index': i,
                    'validator': result,
                    'current': {
                        'functional': current_functional,
                        'stop_codon': current_stop_codon,
                        'vj_in_frame': current_vj_in_frame,
                        'junction_aa': current_junction_aa
                    }
                })

        self.assertEqual(len(mismatches), 0,
                        f"Kappa mismatches: {mismatches[:3]}")

    def test_light_chain_lambda_comparison(self):
        """Compare validator against lambda light chain sequences."""
        mismatches = []

        for i in range(50):
            seq = LightChainSequence.create_random(self.lambda_config)

            current_functional = seq.functional
            current_stop_codon = seq.stop_codon
            current_vj_in_frame = seq.vj_in_frame
            current_junction_aa = getattr(seq, 'junction_aa', '')

            result = self.validator.assess(
                seq.ungapped_seq,
                seq.junction_start,
                seq.junction_end
            )

            if (result.functional != current_functional or
                result.stop_codon != current_stop_codon or
                result.vj_in_frame != current_vj_in_frame or
                result.junction_aa != current_junction_aa):

                mismatches.append({
                    'index': i,
                    'validator': result,
                    'current': {
                        'functional': current_functional,
                        'stop_codon': current_stop_codon,
                        'vj_in_frame': current_vj_in_frame,
                        'junction_aa': current_junction_aa
                    }
                })

        self.assertEqual(len(mismatches), 0,
                        f"Lambda mismatches: {mismatches[:3]}")


class TestFunctionalityValidatorEdgeCases(unittest.TestCase):
    """Test edge cases and boundary conditions."""

    def setUp(self):
        self.validator = FunctionalityValidator()

    def test_empty_sequence(self):
        """Test with empty sequence."""
        result = self.validator.assess("", 0, 0)
        self.assertFalse(result.functional)

    def test_very_short_sequence(self):
        """Test with sequence shorter than a codon."""
        result = self.validator.assess("AT", 0, 2)
        self.assertFalse(result.functional)

    def test_single_codon_junction(self):
        """Test junction that is exactly one codon."""
        # TGT = Cysteine, but can't end with F/W in single codon
        result = self.validator.assess("TGT", 0, 3)
        self.assertFalse(result.functional)
        self.assertEqual(result.junction_aa, "C")

    def test_lowercase_sequence(self):
        """Test that lowercase sequences are handled."""
        junction = "tgtgccaggtgggactactactactactttttc"
        sequence = junction + "atgatgatg"
        result = self.validator.assess(sequence, 0, 33)
        # Should still work - junction is uppercased internally
        self.assertTrue(result.junction_aa.startswith("C"))


class TestFunctionalityValidatorModularity(unittest.TestCase):
    """Test the modular rule-based design."""

    def test_get_rule_names(self):
        """Test that we can list all rule names."""
        validator = FunctionalityValidator()
        names = validator.get_rule_names()

        self.assertEqual(len(names), 5)
        self.assertIn("Stop Codon Check", names)
        self.assertIn("Frame Alignment Check", names)
        self.assertIn("Junction Translation Check", names)
        self.assertIn("Conserved Cysteine (C) Check", names)
        self.assertIn("Conserved J Anchor (F/W) Check", names)

    def test_remove_rule(self):
        """Test that we can remove a rule."""
        validator = FunctionalityValidator()
        initial_count = len(validator.rules)

        # Remove the anchor check
        removed = validator.remove_rule("Conserved J Anchor (F/W) Check")
        self.assertTrue(removed)
        self.assertEqual(len(validator.rules), initial_count - 1)
        self.assertNotIn("Conserved J Anchor (F/W) Check", validator.get_rule_names())

    def test_remove_nonexistent_rule(self):
        """Test that removing a non-existent rule returns False."""
        validator = FunctionalityValidator()
        removed = validator.remove_rule("Nonexistent Rule")
        self.assertFalse(removed)

    def test_add_custom_rule(self):
        """Test that we can add a custom rule."""
        from GenAIRR.validation import ValidationRule, ValidationRuleResult, ValidationContext

        class AlwaysFailRule(ValidationRule):
            @property
            def name(self):
                return "Always Fail"

            def check(self, context):
                return ValidationRuleResult(passed=False, note="This always fails.")

        validator = FunctionalityValidator()
        validator.add_rule(AlwaysFailRule())

        self.assertIn("Always Fail", validator.get_rule_names())

        # Now any sequence should fail
        result = validator.assess("ATGATGATG", 0, 9)
        self.assertFalse(result.functional)

    def test_validator_without_conserved_checks(self):
        """Test validator with only stop codon and frame checks."""
        validator = FunctionalityValidator()
        validator.remove_rule("Conserved Cysteine (C) Check")
        validator.remove_rule("Conserved J Anchor (F/W) Check")

        # A sequence that has no C at start and no F/W at end
        # but has no stops and is in frame should now pass
        # AAA = Lysine (K), GGG = Glycine (G)
        sequence = "AAAGGGAAA"  # KGK - no C, no F/W
        result = validator.assess(sequence, 0, 9)

        # Should pass because we removed the conserved residue checks
        self.assertTrue(result.functional)


if __name__ == '__main__':
    unittest.main(verbosity=2)
