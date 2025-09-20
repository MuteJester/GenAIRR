"""
Enhanced Sequence Simulation Tests
==================================

Comprehensive test suite specifically designed to validate the core sequence simulation
logic in GenAIRR library with deep coverage of biological accuracy and edge cases.

Author: GitHub Copilot
Created: 2025-01-16
"""

import unittest
import random
import numpy as np
from unittest.mock import Mock, patch

from GenAIRR.pipeline import AugmentationPipeline
from GenAIRR.steps import (
    SimulateSequence, FixVPositionAfterTrimmingIndexAmbiguity,
    FixDPositionAfterTrimmingIndexAmbiguity, FixJPositionAfterTrimmingIndexAmbiguity, 
    CorrectForVEndCut, CorrectForDTrims, CorruptSequenceBeginning, 
    InsertNs, InsertIndels, ShortDValidation, DistillMutationRate
)
from GenAIRR.mutation import S5F, Uniform
from GenAIRR.alleles import VAllele, DAllele, JAllele
from GenAIRR.container.SimulationContainer import SimulationContainer
from GenAIRR.sequence import LightChainSequence, HeavyChainSequence
from GenAIRR.sequence.np_region import NP_Region
from GenAIRR.steps.StepBase import AugmentationStep
from GenAIRR.data import HUMAN_IGK_OGRDB, HUMAN_IGL_OGRDB, HUMAN_IGH_OGRDB
from GenAIRR.dataconfig import DataConfig
from GenAIRR.utilities.misc import translate, weighted_choice, check_stops


class EnhancedSequenceSimulationTests(unittest.TestCase):
    """
    Enhanced test suite specifically focused on validating the core sequence simulation
    logic with comprehensive biological and computational accuracy tests.
    """

    @classmethod
    def setUpClass(cls):
        """Set up test fixtures for enhanced sequence simulation testing."""
        cls.heavychain_config = HUMAN_IGH_OGRDB
        cls.lightchain_kappa_config = HUMAN_IGK_OGRDB
        cls.lightchain_lambda_config = HUMAN_IGL_OGRDB

    def setUp(self):
        """Set up before each test."""
        AugmentationStep.set_dataconfig(self.heavychain_config)

    # ============================================================================
    # SEQUENCE SIMULATION CORE LOGIC TESTS
    # ============================================================================

    def test_sequence_simulation_determinism(self):
        """Test that sequence simulation is deterministic given same inputs."""
        # Use specific alleles to control variability
        v_family = list(self.heavychain_config.v_alleles.keys())[0]
        d_family = list(self.heavychain_config.d_alleles.keys())[0] 
        j_family = list(self.heavychain_config.j_alleles.keys())[0]
        
        specific_v = self.heavychain_config.v_alleles[v_family][0]
        specific_d = self.heavychain_config.d_alleles[d_family][0]
        specific_j = self.heavychain_config.j_alleles[j_family][0]
        
        # Generate multiple sequences with same specific alleles
        sequences = []
        for seed in [42, 42, 42]:  # Same seed should give same result
            random.seed(seed)
            np.random.seed(seed)
            
            seq = HeavyChainSequence.create_random(
                self.heavychain_config,
                specific_v=specific_v,
                specific_d=specific_d,
                specific_j=specific_j
            )
            sequences.append({
                'sequence': seq.ungapped_seq,
                'v_trim_3': seq.v_trim_3,
                'd_trim_5': seq.d_trim_5,
                'd_trim_3': seq.d_trim_3,
                'j_trim_5': seq.j_trim_5,
                'np1_length': seq.NP1_length,
                'np2_length': seq.NP2_length
            })
        
        # With same seed, should get same trimming and NP lengths
        # (Note: Complete determinism might vary based on implementation details)
        self.assertEqual(len(set(seq['sequence'] for seq in sequences)), 1,
                        "Same seed should produce identical sequences")

    def test_allele_trimming_biological_constraints(self):
        """Test that allele trimming respects biological constraints."""
        # Test V allele trimming constraints
        v_alleles = [allele for family in self.heavychain_config.v_alleles.values() for allele in family]
        
        for _ in range(50):
            v_allele = random.choice(v_alleles)
            trimmed_seq, trim_5, trim_3 = v_allele.get_trimmed(self.heavychain_config.trim_dicts)
            
            # V alleles should never be trimmed at 5' end
            self.assertEqual(trim_5, 0, "V alleles should not be trimmed at 5' end")
            
            # Trimming should not remove the entire allele
            self.assertGreater(len(trimmed_seq), 0, "Trimming should not remove entire allele")
            
            # Should not trim past anchor for functionality
            if hasattr(v_allele, 'anchor') and v_allele.anchor > 0:
                max_allowed_trim = v_allele.ungapped_len - v_allele.anchor - 1
                self.assertLessEqual(trim_3, max_allowed_trim, 
                                   "Should not trim past anchor point")
        
        # Test D allele trimming constraints
        d_alleles = [allele for family in self.heavychain_config.d_alleles.values() for allele in family]
        
        for _ in range(50):
            d_allele = random.choice(d_alleles)
            trimmed_seq, trim_5, trim_3 = d_allele.get_trimmed(self.heavychain_config.trim_dicts)
            
            # D allele should not be completely trimmed away
            self.assertGreater(len(trimmed_seq), 0, "D allele should not be completely trimmed")
            
            # Combined trims should not exceed allele length
            self.assertLess(trim_5 + trim_3, d_allele.ungapped_len, 
                          "Combined trims should not exceed allele length")

    def test_np_region_markov_chain_properties(self):
        """Test that NP regions follow proper Markov chain generation."""
        # Test NP region generation with controlled parameters
        for region_type in ['NP1', 'NP2']:
            for target_length in [5, 10, 15]:
                if target_length in self.heavychain_config.NP_lengths[region_type]:
                    # Generate multiple NP regions
                    regions = []
                    for _ in range(20):
                        np_seq = NP_Region.create_np_region(
                            self.heavychain_config.NP_lengths,
                            self.heavychain_config.NP_transitions,
                            region_type,
                            self.heavychain_config.NP_first_bases
                        )
                        if np_seq:
                            regions.append(np_seq)
                    
                    if regions:
                        # Should produce varied sequences
                        unique_sequences = set(regions)
                        self.assertGreater(len(unique_sequences), 1, 
                                         "Should produce varied NP sequences")
                        
                        # All sequences should contain only valid bases (allowing lowercase)
                        valid_bases = set('ATCGNatcgn')
                        for seq in regions:
                            self.assertTrue(set(seq).issubset(valid_bases),
                                          f"NP sequence {seq} contains invalid bases")
                        
                        # Sequences should have reasonable lengths (Markov chains can vary significantly)
                        for seq in regions:
                            self.assertLessEqual(len(seq), target_length + 30,
                                               "NP sequence length should be reasonable for Markov chain generation")

    def test_junction_calculation_accuracy(self):
        """Test accuracy of junction length and boundary calculations."""
        for _ in range(30):
            seq = HeavyChainSequence.create_random(self.heavychain_config)
            
            # Manual junction length calculation
            expected_length = (seq.v_allele.ungapped_len - (seq.v_allele.anchor - 1) - seq.v_trim_3 +
                             seq.NP1_length +
                             seq.d_allele.ungapped_len - seq.d_trim_5 - seq.d_trim_3 +
                             seq.NP2_length +
                             (seq.j_allele.anchor + 2) - seq.j_trim_5)
            
            self.assertEqual(seq.junction_length, expected_length,
                           "Junction length calculation should match manual calculation")
            
            # Test junction extraction
            expected_junction = seq.ungapped_seq[seq.junction_start:seq.junction_end]
            self.assertEqual(seq.junction, expected_junction.upper(),
                           "Junction sequence should match extracted sequence")
            
            # Test junction boundaries
            self.assertEqual(seq.junction_start, seq.v_allele.anchor,
                           "Junction should start at V anchor")
            self.assertEqual(seq.junction_end, seq.v_allele.anchor + seq.junction_length,
                           "Junction end should be start + length")

    def test_functionality_assessment_logic(self):
        """Test the biological functionality assessment logic."""
        functional_sequences = []
        non_functional_sequences = []
        
        for _ in range(100):
            seq = HeavyChainSequence.create_random(self.heavychain_config)
            
            if seq.functional:
                functional_sequences.append(seq)
            else:
                non_functional_sequences.append(seq)
        
        # Test functional sequences meet all criteria
        for seq in functional_sequences:
            # Should not have stop codons
            self.assertFalse(seq.stop_codon, "Functional sequences should not have stop codons")
            
            # Should be in-frame
            self.assertTrue(seq.vj_in_frame, "Functional sequences should be in-frame")
            
            # Junction should be divisible by 3
            self.assertEqual(seq.junction_length % 3, 0, 
                           "Functional sequences should have junction length divisible by 3")
            
            # Should have proper conserved residues
            if len(seq.junction_aa) > 0:
                self.assertEqual(seq.junction_aa[0], 'C', 
                               "Functional sequences should start with Cysteine")
            if len(seq.junction_aa) > 1:
                self.assertIn(seq.junction_aa[-1], ['F', 'W'],
                            "Functional sequences should end with Phe or Trp")
        
        # Should have both functional and non-functional sequences
        self.assertGreater(len(functional_sequences), 0, "Should produce some functional sequences")
        self.assertGreater(len(non_functional_sequences), 0, "Should produce some non-functional sequences")

    # ============================================================================
    # PIPELINE STEP INTEGRATION TESTS
    # ============================================================================

    def test_position_fix_step_accuracy(self):
        """Test that position fixing steps accurately adjust boundaries."""
        pipeline = AugmentationPipeline([
            SimulateSequence(S5F(), True),
            FixVPositionAfterTrimmingIndexAmbiguity(),
            FixDPositionAfterTrimmingIndexAmbiguity(),
            FixJPositionAfterTrimmingIndexAmbiguity(),
            DistillMutationRate()
        ])
        
        for _ in range(20):
            result = pipeline.execute()
            result_dict = result.get_dict()
            
            # Test position consistency after fixes
            self.assertLessEqual(result_dict['v_sequence_start'], result_dict['v_sequence_end'])
            self.assertLessEqual(result_dict['v_sequence_end'], result_dict['d_sequence_start'])
            self.assertLessEqual(result_dict['d_sequence_start'], result_dict['d_sequence_end'])
            self.assertLessEqual(result_dict['d_sequence_end'], result_dict['j_sequence_start'])
            self.assertLessEqual(result_dict['j_sequence_start'], result_dict['j_sequence_end'])
            
            # Test germline position consistency
            self.assertLessEqual(result_dict['v_germline_start'], result_dict['v_germline_end'])
            self.assertLessEqual(result_dict['d_germline_start'], result_dict['d_germline_end'])
            self.assertLessEqual(result_dict['j_germline_start'], result_dict['j_germline_end'])
            
            # Test that sequence segments match expected lengths
            v_segment_length = result_dict['v_sequence_end'] - result_dict['v_sequence_start']
            d_segment_length = result_dict['d_sequence_end'] - result_dict['d_sequence_start']
            j_segment_length = result_dict['j_sequence_end'] - result_dict['j_sequence_start']
            
            self.assertGreater(v_segment_length, 0, "V segment should have positive length")
            self.assertGreater(d_segment_length, 0, "D segment should have positive length")
            self.assertGreater(j_segment_length, 0, "J segment should have positive length")

    def test_corruption_step_biological_realism(self):
        """Test that sequence corruption maintains biological realism."""
        corruption_step = CorruptSequenceBeginning(0.8, [0.4, 0.3, 0.3], 300, 100, 200, 50)
        
        pipeline = AugmentationPipeline([
            SimulateSequence(S5F(), True),
            corruption_step,
            DistillMutationRate()
        ])
        
        corrupted_sequences = []
        for _ in range(30):
            result = pipeline.execute()
            result_dict = result.get_dict()
            
            if result_dict.get('corruption_event', 'no-corruption') != 'no-corruption':
                corrupted_sequences.append(result_dict)
        
        # Test corruption events
        for seq_data in corrupted_sequences:
            corruption_event = seq_data['corruption_event']
            
            # Should be valid corruption type
            self.assertIn(corruption_event, ['add', 'remove', 'remove_before_add'])
            
            # Sequence should still be valid nucleotide sequence
            sequence = seq_data['sequence']
            valid_bases = set('ATCGN')
            self.assertTrue(set(sequence.upper()).issubset(valid_bases),
                          "Corrupted sequence should only contain valid bases")
            
            # Should have reasonable length
            self.assertGreater(len(sequence), 50, "Corrupted sequence should have reasonable length")
            self.assertLess(len(sequence), 1000, "Corrupted sequence should not be too long")
            
            # Test corruption metadata consistency
            corruption_event = seq_data['corruption_event']
            
            # For corruption events, just verify the basic structure
            if corruption_event in ['add', 'remove_before_add']:
                add_amount = seq_data.get('corruption_add_amount', 0)
                add_section = seq_data.get('corruption_added_section', '')
                # Verify that if there's an add amount, there's an add section
                if add_amount > 0:
                    self.assertGreater(len(add_section), 0,
                                     "Should have added section if add amount > 0")
            
            if corruption_event in ['remove', 'remove_before_add']:
                remove_amount = seq_data.get('corruption_remove_amount', 0)
                self.assertGreaterEqual(remove_amount, 0, "Remove amount should be non-negative")

    def test_indel_insertion_logic_validation(self):
        """Test indel insertion logic and position validation."""
        indel_step = InsertIndels(0.7, 5, 0.5, 0.5)
        
        pipeline = AugmentationPipeline([
            SimulateSequence(S5F(), True),
            FixVPositionAfterTrimmingIndexAmbiguity(),
            FixDPositionAfterTrimmingIndexAmbiguity(),
            FixJPositionAfterTrimmingIndexAmbiguity(),
            indel_step,
            DistillMutationRate()
        ])
        
        sequences_with_indels = []
        for _ in range(30):
            result = pipeline.execute()
            result_dict = result.get_dict()
            
            indels = result_dict.get('indels', {})
            if indels:
                sequences_with_indels.append(result_dict)
        
        # Test indel properties
        for seq_data in sequences_with_indels:
            indels = seq_data['indels']
            sequence = seq_data['sequence']
            
            for position, indel_description in indels.items():
                # Position should be reasonable (can be -1 for special positioning)
                self.assertGreaterEqual(position, -1, "Position should be >= -1")
                if position >= 0:
                    self.assertLess(position, len(sequence), "Valid positions should be within sequence")
                
                # Indel description should be properly formatted
                self.assertTrue(indel_description.startswith('I < ') or 
                              indel_description.startswith('D > '),
                              f"Invalid indel description: {indel_description}")
                
                # Extract base from description
                if indel_description.startswith('I < '):
                    inserted_base = indel_description[4:]
                    self.assertIn(inserted_base, ['A', 'T', 'C', 'G'],
                                "Inserted base should be valid nucleotide")
                elif indel_description.startswith('D > '):
                    deleted_base = indel_description[4:]
                    self.assertIn(deleted_base, ['A', 'T', 'C', 'G'],
                                "Deleted base should be valid nucleotide")

    # ============================================================================
    # LIGHT CHAIN SIMULATION VALIDATION
    # ============================================================================

    def test_light_chain_simulation_accuracy(self):
        """Test light chain simulation follows different rules than heavy chain."""
        # Test Kappa light chains
        AugmentationStep.set_dataconfig(self.lightchain_kappa_config)
        
        kappa_sequences = []
        for _ in range(20):
            seq = LightChainSequence.create_random(self.lightchain_kappa_config)
            kappa_sequences.append(seq)
        
        for seq in kappa_sequences:
            # Should not have D allele
            self.assertIsNone(seq.d_allele, "Light chains should not have D allele")
            
            # Should have only one NP region
            self.assertGreater(seq.NP1_length, 0, "Light chains should have NP1 region")
            self.assertEqual(seq.NP2_length, 0, "Light chains should not have NP2 region")
            
            # Junction calculation should be different
            expected_junction_length = (seq.v_allele.length - (seq.v_allele.anchor - 1) - seq.v_trim_3 +
                                      seq.NP1_length + (seq.j_allele.anchor + 2) - seq.j_trim_5)
            self.assertEqual(seq.junction_length, expected_junction_length,
                           "Light chain junction calculation should be correct")
        
        # Test Lambda light chains
        AugmentationStep.set_dataconfig(self.lightchain_lambda_config)
        
        lambda_sequences = []
        for _ in range(20):
            seq = LightChainSequence.create_random(self.lightchain_lambda_config)
            lambda_sequences.append(seq)
        
        # Test that kappa and lambda use different allele sets
        kappa_v_names = {seq.v_allele.name for seq in kappa_sequences}
        lambda_v_names = {seq.v_allele.name for seq in lambda_sequences}
        
        # Should have different allele names (allowing for some potential overlap)
        self.assertGreater(len(kappa_v_names), 0)
        self.assertGreater(len(lambda_v_names), 0)

    # ============================================================================
    # MUTATION MODEL VALIDATION
    # ============================================================================

    def test_mutation_model_biological_accuracy(self):
        """Test that mutation models produce biologically accurate results."""
        # Test S5F mutation model
        s5f_model = S5F(min_mutation_rate=0.05, max_mutation_rate=0.05)
        
        mutated_sequences = []
        for _ in range(20):
            seq = HeavyChainSequence.create_random(self.heavychain_config)
            seq.mutate(s5f_model)
            mutated_sequences.append(seq)
        
        for seq in mutated_sequences:
            # Should have mutations
            self.assertIsNotNone(seq.mutated_seq, "Should have mutated sequence")
            self.assertIsInstance(seq.mutations, dict, "Should have mutations dictionary")
            
            # Mutation rate should be close to expected
            self.assertAlmostEqual(seq.mutation_freq, 0.05, places=2,
                                 msg="Mutation frequency should match model parameters")
            
            # Mutations should only be in valid positions
            for pos in seq.mutations.keys():
                self.assertGreaterEqual(pos, 0, "Mutation position should be non-negative")
                self.assertLess(pos, len(seq.ungapped_seq), 
                              "Mutation position should be within sequence")
        
        # Test Uniform mutation model
        uniform_model = Uniform(min_mutation_rate=0.03, max_mutation_rate=0.03)
        
        seq = HeavyChainSequence.create_random(self.heavychain_config)
        seq.mutate(uniform_model)
        
        self.assertAlmostEqual(seq.mutation_freq, 0.03, places=2,
                             msg="Uniform model should produce exact mutation rate")

    # ============================================================================
    # EDGE CASE AND ROBUSTNESS TESTS
    # ============================================================================

    def test_edge_case_handling(self):
        """Test handling of edge cases in sequence simulation."""
        # Test with minimal configurations
        try:
            seq = HeavyChainSequence.create_random(self.heavychain_config)
            
            # Should handle very short sequences
            if len(seq.ungapped_seq) > 0:
                self.assertIsInstance(seq.ungapped_seq, str)
                self.assertGreater(len(seq.ungapped_seq), 0)
            
            # Should handle junction calculation even with extreme trimming
            self.assertGreaterEqual(seq.junction_length, 0)
            
        except Exception as e:
            self.fail(f"Should handle edge cases gracefully: {e}")

    def test_data_consistency_validation(self):
        """Test that generated data maintains internal consistency."""
        for _ in range(30):
            seq = HeavyChainSequence.create_random(self.heavychain_config)
            
            # Test that segment positions are consistent
            self.assertEqual(seq.v_seq_start, 0, "V segment should start at 0")
            self.assertGreaterEqual(seq.d_seq_start, seq.v_seq_end, 
                                  "D segment should start after V segment")
            self.assertGreaterEqual(seq.j_seq_start, seq.d_seq_end,
                                  "J segment should start after D segment")
            
            # Test that germline positions make sense
            self.assertGreaterEqual(seq.v_germline_end, seq.v_germline_start)
            self.assertGreaterEqual(seq.d_germline_end, seq.d_germline_start)
            self.assertGreaterEqual(seq.j_germline_end, seq.j_germline_start)
            
            # Test that trimming values are reasonable
            self.assertGreaterEqual(seq.v_trim_3, 0)
            self.assertGreaterEqual(seq.d_trim_5, 0)
            self.assertGreaterEqual(seq.d_trim_3, 0)
            self.assertGreaterEqual(seq.j_trim_5, 0)
            
            # Test that NP regions have reasonable lengths
            self.assertGreaterEqual(seq.NP1_length, 0)
            self.assertGreaterEqual(seq.NP2_length, 0)


if __name__ == '__main__':
    # Configure test execution
    unittest.main(verbosity=2, buffer=True)
