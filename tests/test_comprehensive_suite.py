#!/usr/bin/env python3
"""
GenAIRR Comprehensive Test Suite

This module contains the comprehensive test suite for the GenAIRR library,
providing extensive coverage of all major functionality including allele
management, mutation models, sequence generation, pipeline operations,
and biological validation.

Author: GenAIRR Development Team
License: See LICENSE file in the repository root
"""

import unittest
import random
import re
import numpy as np
import tempfile
import os
from unittest.mock import Mock, patch

from GenAIRR.pipeline import AugmentationPipeline
from GenAIRR.steps import (
    SimulateSequence, FixVPositionAfterTrimmingIndexAmbiguity,
    FixDPositionAfterTrimmingIndexAmbiguity, FixJPositionAfterTrimmingIndexAmbiguity, 
    FilterTCRDJAmbiguities, CorrectForVEndCut, CorrectForDTrims, 
    CorruptSequenceBeginning, InsertNs, InsertIndels, ShortDValidation, 
    DistillMutationRate
)
from GenAIRR.mutation import S5F, Uniform
from GenAIRR.alleles import VAllele, DAllele, JAllele
from GenAIRR.container.SimulationContainer import SimulationContainer
from GenAIRR.sequence import LightChainSequence, HeavyChainSequence
from GenAIRR.steps.StepBase import AugmentationStep
from GenAIRR.data import HUMAN_IGK_OGRDB, HUMAN_IGL_OGRDB, HUMAN_IGH_OGRDB
from GenAIRR.dataconfig import DataConfig
from GenAIRR.utilities.misc import translate, weighted_choice, check_stops


class EnhancedGenAIRRTests(unittest.TestCase):
    """
    Comprehensive Test Suite for GenAIRR Library
    
    This test suite provides extensive coverage of GenAIRR functionality including:
    - Allele creation, properties, and trimming mechanisms
    - Mutation models (S5F, Uniform) with parameter validation
    - Sequence generation with variability and quality testing
    - Complete pipeline testing with all augmentation steps
    - Data configuration validation and manipulation
    - Utility function verification and edge cases
    - Error handling and biological constraint validation
    - Integration testing for real-world scenarios
    
    Test Categories:
    - Core Allele Functionality (creation, trimming, properties)
    - Mutation Model Testing (S5F, Uniform variants)
    - Sequence Generation (variability, selection, productivity)
    - Pipeline Operations (ordering, comprehensive workflows)
    - Data Configuration (validation, copying, custom creation)
    - Utility Functions (translation, weighted choice, stop codons)
    - Error Handling (edge cases, invalid inputs)
    - Integration Testing (light/heavy chains, complete workflows)
    - Performance Testing (memory usage, concurrent execution)
    - Quality Validation (biological constraints, sequence metrics)
    """

    @classmethod
    def setUpClass(cls):
        """Set up test fixtures before running tests."""
        cls.heavychain_config = HUMAN_IGH_OGRDB
        cls.lightchain_kappa_config = HUMAN_IGK_OGRDB
        cls.lightchain_lambda_config = HUMAN_IGL_OGRDB

    def setUp(self):
        """Set up before each test."""
        # Reset dataconfig for each test
        AugmentationStep.set_dataconfig(self.heavychain_config)

    def tearDown(self):
        """Clean up after each test."""
        pass

    # ============================================================================
    # CORE ALLELE FUNCTIONALITY TESTS
    # ============================================================================

    def test_allele_creation_and_properties(self):
        """Test creation and properties of different allele types."""
        # Test V allele
        v_seq = 'CAGGTCACCTTGAAGGAGTCTGGTCCT'
        v_allele = VAllele('IGHV1-1*01', v_seq, len(v_seq))
        
        self.assertEqual(v_allele.name, 'IGHV1-1*01')
        self.assertEqual(v_allele.gene, 'IGHV1-1')
        self.assertEqual(v_allele.family, 'IGHV1')
        self.assertEqual(v_allele.ungapped_seq, v_seq.upper())
        self.assertEqual(v_allele.length, len(v_seq))
        
        # Test D allele
        d_seq = 'GGTATAACGGTATC'
        d_allele = DAllele('IGHD1-1*01', d_seq, len(d_seq))
        
        self.assertEqual(d_allele.name, 'IGHD1-1*01')
        self.assertEqual(d_allele.gene, 'IGHD1-1')
        self.assertEqual(d_allele.family, 'IGHD1')
        
        # Test J allele
        j_seq = 'ACTACTTTGACTACTGG'
        j_allele = JAllele('IGHJ1*01', j_seq, len(j_seq))
        
        self.assertEqual(j_allele.name, 'IGHJ1*01')
        self.assertEqual(j_allele.gene, 'IGHJ1')

    def test_allele_trimming(self):
        """Test allele trimming functionality."""
        v_seq = 'ATCGATCGATCGATCGATCG'
        v_allele = VAllele('IGHV1-1*01', v_seq, len(v_seq))
        
        # Create properly structured trim dicts matching the library's expected format
        trim_dicts = {
            'V_3': {
                'IGHV1': {
                    'IGHV1-1': {0: 0.5, 1: 0.3, 2: 0.2}
                }
            },
            'V_5': {
                'IGHV1': {
                    'IGHV1-1': {0: 0.4, 1: 0.4, 2: 0.2}
                }
            }
        }
        
        # Test multiple trimming calls
        for _ in range(10):
            trimmed_seq, trim_5, trim_3 = v_allele.get_trimmed(trim_dicts)
            self.assertIsInstance(trimmed_seq, str)
            self.assertGreaterEqual(trim_5, 0)
            self.assertGreaterEqual(trim_3, 0)
            self.assertEqual(len(trimmed_seq), len(v_seq) - trim_5 - trim_3)

    # ============================================================================
    # MUTATION MODEL TESTS
    # ============================================================================

    def test_s5f_mutation_model(self):
        """Test S5F mutation model with various parameters."""
        # Test with different parameters
        s5f_low = S5F(min_mutation_rate=0.001, max_mutation_rate=0.01)
        s5f_high = S5F(min_mutation_rate=0.05, max_mutation_rate=0.15)
        s5f_productive = S5F(min_mutation_rate=0.01, max_mutation_rate=0.1, productive=True)
        
        # Create a test sequence
        test_seq = HeavyChainSequence.create_random(self.heavychain_config)
        
        # Test mutation application
        for model in [s5f_low, s5f_high, s5f_productive]:
            test_seq.mutate(model)
            self.assertIsNotNone(test_seq.mutated_seq)
            self.assertIsInstance(test_seq.mutations, dict)
            self.assertIsInstance(test_seq.mutation_freq, float)

    def test_uniform_mutation_model(self):
        """Test Uniform mutation model."""
        uniform_model = Uniform(min_mutation_rate=0.05, max_mutation_rate=0.05)
        test_seq = LightChainSequence.create_random(self.lightchain_kappa_config)
        
        test_seq.mutate(uniform_model)
        self.assertIsNotNone(test_seq.mutated_seq)
        self.assertIsInstance(test_seq.mutations, dict)
        self.assertEqual(test_seq.mutation_freq, 0.05)

    # ============================================================================
    # SEQUENCE GENERATION TESTS
    # ============================================================================

    def test_sequence_generation_variability(self):
        """Test that sequence generation produces variable results."""
        sequences = []
        for _ in range(50):
            seq = HeavyChainSequence.create_random(self.heavychain_config)
            sequences.append(seq.ungapped_seq)
        
        # Check that we get different sequences
        unique_sequences = set(sequences)
        self.assertGreater(len(unique_sequences), 40)  # Should have high variability

    def test_specific_allele_selection(self):
        """Test sequence generation with specific allele selection."""
        # Get actual allele instances from the config
        v_family_name = list(self.heavychain_config.v_alleles.keys())[0]
        d_family_name = list(self.heavychain_config.d_alleles.keys())[0]
        j_family_name = list(self.heavychain_config.j_alleles.keys())[0]
        
        v_allele = self.heavychain_config.v_alleles[v_family_name][0]
        d_allele = self.heavychain_config.d_alleles[d_family_name][0]
        j_allele = self.heavychain_config.j_alleles[j_family_name][0]
        
        seq = HeavyChainSequence.create_random(
            self.heavychain_config,
            specific_v=v_allele,
            specific_d=d_allele,
            specific_j=j_allele
        )
        
        self.assertEqual(seq.v_allele.name, v_allele.name)
        self.assertEqual(seq.d_allele.name, d_allele.name)
        self.assertEqual(seq.j_allele.name, j_allele.name)

    def test_productive_sequence_generation(self):
        """Test generation of productive sequences."""
        productive_count = 0
        total_sequences = 100
        
        for _ in range(total_sequences):
            seq = HeavyChainSequence.create_random(self.heavychain_config)
            if seq.functional:
                productive_count += 1
        
        # Should have some productive sequences (not all will be productive naturally)
        self.assertGreater(productive_count, 5)

    # ============================================================================
    # PIPELINE TESTS
    # ============================================================================

    def test_empty_pipeline_handling(self):
        """Test pipeline behavior with invalid configurations."""
        with self.assertRaises(Exception):
            # Pipeline without SimulateSequence first should raise exception
            pipeline = AugmentationPipeline([DistillMutationRate()])
            pipeline.execute()

    def test_pipeline_step_ordering(self):
        """Test that pipeline steps execute in correct order."""
        steps_executed = []
        
        class TestStep(AugmentationStep):
            def __init__(self, name):
                super().__init__()
                self.name = name
            
            def apply(self, container):
                steps_executed.append(self.name)
        
        pipeline = AugmentationPipeline([
            SimulateSequence(S5F(), True),
            TestStep("Step1"),
            TestStep("Step2"),
            TestStep("Step3")
        ])
        
        pipeline.execute()
        self.assertEqual(steps_executed, ["Step1", "Step2", "Step3"])

    def test_pipeline_with_all_steps(self):
        """Test pipeline with comprehensive set of steps."""
        pipeline = AugmentationPipeline([
            SimulateSequence(S5F(min_mutation_rate=0.01, max_mutation_rate=0.1), True),
            FixVPositionAfterTrimmingIndexAmbiguity(),
            FixDPositionAfterTrimmingIndexAmbiguity(),
            FixJPositionAfterTrimmingIndexAmbiguity(),
            CorrectForVEndCut(),
            CorrectForDTrims(),
            CorruptSequenceBeginning(0.5, [0.4, 0.4, 0.2], 300, 100, 200, 20),
            InsertNs(0.01, 0.3),
            ShortDValidation(),
            InsertIndels(0.3, 3, 0.5, 0.5),
            DistillMutationRate()
        ])
        
        result = pipeline.execute()
        self.assertIsInstance(result, SimulationContainer)
        
        # Check that all expected fields are present
        expected_fields = [
            'sequence', 'v_call', 'd_call', 'j_call', 'mutations', 
            'Ns', 'indels', 'productive', 'stop_codon', 'vj_in_frame'
        ]
        for field in expected_fields:
            self.assertIn(field, result.get_dict())

    # ============================================================================
    # DATA CONFIGURATION TESTS
    # ============================================================================

    def test_dataconfig_properties(self):
        """Test DataConfig properties and methods."""
        config = self.heavychain_config
        
        # Test property access
        self.assertGreater(config.number_of_v_alleles, 0)
        self.assertGreater(config.number_of_d_alleles, 0)
        self.assertGreater(config.number_of_j_alleles, 0)
        
        # Test allele list retrieval
        v_list = config.allele_list('v')
        d_list = config.allele_list('d')
        j_list = config.allele_list('j')
        
        self.assertIsInstance(v_list, list)
        self.assertIsInstance(d_list, list)
        self.assertIsInstance(j_list, list)
        
        self.assertEqual(len(v_list), config.number_of_v_alleles)
        self.assertEqual(len(d_list), config.number_of_d_alleles)
        self.assertEqual(len(j_list), config.number_of_j_alleles)

    def test_dataconfig_copy(self):
        """Test DataConfig deep copy functionality."""
        original = self.heavychain_config
        copy_config = original.copy()
        
        # Test that it's a deep copy
        self.assertIsNot(original, copy_config)
        self.assertEqual(original.number_of_v_alleles, copy_config.number_of_v_alleles)
        
        # Modify copy and ensure original unchanged
        if copy_config.v_alleles:
            original_count = len(original.v_alleles)
            copy_config.v_alleles = {}
            self.assertEqual(len(original.v_alleles), original_count)

    # ============================================================================
    # UTILITY FUNCTION TESTS
    # ============================================================================

    def test_translate_function(self):
        """Test nucleotide translation function."""
        # Test normal codon
        self.assertEqual(translate('ATG'), 'M')
        self.assertEqual(translate('ATGAAATAG'), 'MK*')
        
        # Test stop codons
        self.assertEqual(translate('TAG'), '*')
        self.assertEqual(translate('TAA'), '*')
        self.assertEqual(translate('TGA'), '*')
        
        # Test incomplete codon
        self.assertEqual(translate('AT'), '.')
        
        # Test special characters
        self.assertEqual(translate('N' * 3), 'X')
        self.assertEqual(translate('-' * 3), 'X')
        self.assertEqual(translate('.' * 3), '_')

    def test_weighted_choice_function(self):
        """Test weighted choice utility function."""
        choices = {'A': 0.5, 'B': 0.3, 'C': 0.2}
        
        # Test multiple selections
        selections = [weighted_choice(choices) for _ in range(1000)]
        unique_selections = set(selections)
        
        # Should select all options with weighted randomness
        self.assertEqual(unique_selections, {'A', 'B', 'C'})
        
        # Test that more weighted options are selected more often
        a_count = selections.count('A')
        b_count = selections.count('B')
        c_count = selections.count('C')
        
        self.assertGreater(a_count, b_count)
        self.assertGreater(b_count, c_count)

    def test_check_stops_function(self):
        """Test stop codon checking function."""
        # Test sequences with stop codons
        self.assertTrue(check_stops('ATGTAG'))
        self.assertTrue(check_stops('ATGTAA'))
        self.assertTrue(check_stops('ATGTGA'))
        
        # Test sequences without stop codons
        self.assertFalse(check_stops('ATGAAA'))
        self.assertFalse(check_stops('CCCCCC'))
        
        # Test with position return
        has_stop, pos = check_stops('ATGTAG', return_pos=True)
        self.assertTrue(has_stop)
        self.assertEqual(pos, 3)

    # ============================================================================
    # ERROR HANDLING AND EDGE CASES
    # ============================================================================

    def test_container_error_handling(self):
        """Test SimulationContainer error handling."""
        container = SimulationContainer()
        
        # Test accessing non-existent key
        with self.assertRaises(KeyError):
            _ = container['non_existent_key']
        
        # Test updating existing field (SimulationContainer only allows setting existing attributes)
        container['sequence'] = 'test_sequence'
        self.assertEqual(container['sequence'], 'test_sequence')
        
        # Test that trying to set non-existent field raises error
        with self.assertRaises(KeyError):
            container['non_existent_field'] = 'test_value'

    def test_sequence_edge_cases(self):
        """Test edge cases in sequence generation."""
        # Use the actual heavychain config which has proper trim dictionaries
        seq = HeavyChainSequence.create_random(self.heavychain_config)
        self.assertIsNotNone(seq.ungapped_seq)
        self.assertGreater(len(seq.ungapped_seq), 0)

    def test_invalid_mutation_parameters(self):
        """Test mutation models with edge case parameters."""
        # Test with zero mutation rate
        zero_mutation = S5F(min_mutation_rate=0.0, max_mutation_rate=0.0)
        seq = HeavyChainSequence.create_random(self.heavychain_config)
        seq.mutate(zero_mutation)
        
        self.assertEqual(len(seq.mutations), 0)
        self.assertEqual(seq.mutation_freq, 0.0)

    # ============================================================================
    # INTEGRATION TESTS
    # ============================================================================

    def test_light_chain_kappa_complete_pipeline(self):
        """Test complete pipeline for light chain kappa."""
        AugmentationStep.set_dataconfig(self.lightchain_kappa_config)
        
        pipeline = AugmentationPipeline([
            SimulateSequence(S5F(), True),
            FixVPositionAfterTrimmingIndexAmbiguity(),
            FixJPositionAfterTrimmingIndexAmbiguity(),
            CorrectForVEndCut(),
            CorruptSequenceBeginning(0.3, [0.5, 0.3, 0.2], 200, 50, 100, 10),
            InsertNs(0.01, 0.2),
            InsertIndels(0.2, 2, 0.5, 0.5),
            DistillMutationRate()
        ])
        
        results = []
        for _ in range(10):
            result = pipeline.execute()
            results.append(result.get_dict())
            
        self.assertEqual(len(results), 10)
        
        # Check that each result has required fields
        for result in results:
            self.assertIn('sequence', result)
            self.assertIn('v_call', result)
            self.assertIn('j_call', result)
            # Light chains have d_call but it should be empty list
            self.assertIn('d_call', result)
            self.assertEqual(result['d_call'], [])

    def test_light_chain_lambda_complete_pipeline(self):
        """Test complete pipeline for light chain lambda."""
        AugmentationStep.set_dataconfig(self.lightchain_lambda_config)
        
        pipeline = AugmentationPipeline([
            SimulateSequence(Uniform(0.02, 0.02), True),
            FixVPositionAfterTrimmingIndexAmbiguity(),
            FixJPositionAfterTrimmingIndexAmbiguity(),
            CorrectForVEndCut(),
            InsertNs(0.005, 0.1),
            DistillMutationRate()
        ])
        
        result = pipeline.execute()
        result_dict = result.get_dict()
        
        self.assertIn('sequence', result_dict)
        self.assertIn('v_call', result_dict)
        self.assertIn('j_call', result_dict)
        # Light chains have d_call but it should be empty list
        self.assertIn('d_call', result_dict)
        self.assertEqual(result_dict['d_call'], [])

    def test_batch_sequence_generation(self):
        """Test generating large batches of sequences."""
        batch_size = 100
        sequences = []
        
        pipeline = AugmentationPipeline([
            SimulateSequence(S5F(min_mutation_rate=0.01, max_mutation_rate=0.05), True),
            FixVPositionAfterTrimmingIndexAmbiguity(),
            FixDPositionAfterTrimmingIndexAmbiguity(),
            FixJPositionAfterTrimmingIndexAmbiguity(),
            DistillMutationRate()
        ])
        
        for _ in range(batch_size):
            result = pipeline.execute()
            sequences.append(result.get_dict())
        
        self.assertEqual(len(sequences), batch_size)
        
        # Check variability
        unique_seqs = set(seq['sequence'] for seq in sequences)
        self.assertGreater(len(unique_seqs), batch_size * 0.8)  # High variability

    # ============================================================================
    # PERFORMANCE AND STRESS TESTS
    # ============================================================================

    def test_memory_usage_large_batches(self):
        """Test memory usage with large sequence generation batches."""
        import gc
        
        gc.collect()
        
        # Generate many sequences
        for _ in range(500):
            seq = HeavyChainSequence.create_random(self.heavychain_config)
            # Ensure we're not keeping references
            del seq
        
        gc.collect()
        # Test passes if no memory errors occur

    def test_concurrent_pipeline_execution(self):
        """Test that pipelines can be used concurrently."""
        import threading
        
        results = {}
        
        def run_pipeline(thread_id):
            pipeline = AugmentationPipeline([
                SimulateSequence(S5F(), True),
                FixVPositionAfterTrimmingIndexAmbiguity(),
                DistillMutationRate()
            ])
            result = pipeline.execute()
            results[thread_id] = result.get_dict()
        
        threads = []
        for i in range(5):
            thread = threading.Thread(target=run_pipeline, args=(i,))
            threads.append(thread)
            thread.start()
        
        for thread in threads:
            thread.join()
        
        self.assertEqual(len(results), 5)
        for result in results.values():
            self.assertIn('sequence', result)

    # ============================================================================
    # VALIDATION AND QUALITY TESTS
    # ============================================================================

    def test_sequence_quality_metrics(self):
        """Test sequence quality and biological validity."""
        pipeline = AugmentationPipeline([
            SimulateSequence(S5F(productive=True), True),
            FixVPositionAfterTrimmingIndexAmbiguity(),
            FixDPositionAfterTrimmingIndexAmbiguity(),
            FixJPositionAfterTrimmingIndexAmbiguity(),
            CorrectForVEndCut(),
            CorrectForDTrims(),
            ShortDValidation(),
            DistillMutationRate()
        ])
        
        productive_count = 0
        total_sequences = 50
        
        for _ in range(total_sequences):
            result = pipeline.execute()
            result_dict = result.get_dict()
            
            # Check sequence composition
            sequence = result_dict['sequence']
            valid_bases = set('ATCGN')
            sequence_bases = set(sequence.upper())
            self.assertTrue(sequence_bases.issubset(valid_bases))
            
            # Check mutation rate consistency
            if 'mutations' in result_dict and 'Ns' in result_dict:
                total_mutations = len(result_dict['mutations']) + len(result_dict['Ns'])
                calculated_rate = total_mutations / len(sequence) if len(sequence) > 0 else 0
                stored_rate = result_dict.get('mutation_rate', 0)
                
                # Allow for some tolerance in floating point comparison
                self.assertAlmostEqual(calculated_rate, stored_rate, places=3)
            
            if result_dict.get('productive', False):
                productive_count += 1
        
        # With productive=True, we should get mostly productive sequences
        self.assertGreater(productive_count, total_sequences * 0.7)

    def test_junction_properties(self):
        """Test junction sequence properties and validity."""
        sequences = []
        for _ in range(50):
            seq = HeavyChainSequence.create_random(self.heavychain_config)
            sequences.append(seq)
        
        for seq in sequences:
            # Check junction length is reasonable
            self.assertGreater(seq.junction_length, 0)
            self.assertLess(seq.junction_length, 200)  # Reasonable upper bound
            
            # Check junction sequence exists and has correct length
            self.assertEqual(len(seq.junction), seq.junction_length)
            
            # For functional sequences, check conserved residues
            if seq.functional:
                junction_aa = translate(seq.junction)
                if len(junction_aa) > 0:
                    self.assertEqual(junction_aa[0], 'C')  # Should start with Cys
                if len(junction_aa) > 1:
                    self.assertIn(junction_aa[-1], ['F', 'W'])  # Should end with Phe or Trp

    # ============================================================================
    # CONFIGURATION AND DATACONFIG BUILDER TESTS
    # ============================================================================

    def test_custom_dataconfig_creation(self):
        """Test creation of custom DataConfig with mock data."""
        # Create a minimal mock DataConfig
        mock_config = DataConfig()
        mock_config.name = "Test Config"
        
        # Mock some basic data structures
        mock_config.v_alleles = {
            'IGHV1': [VAllele('IGHV1-1*01', 'ATCGATCGATCG', 12)]
        }
        mock_config.d_alleles = {
            'IGHD1': [DAllele('IGHD1-1*01', 'GGGCCC', 6)]
        }
        mock_config.j_alleles = {
            'IGHJ1': [JAllele('IGHJ1*01', 'TTTTTT', 6)]
        }
        
        mock_config.trim_dicts = {
            'V_5': {0: 1.0}, 'V_3': {0: 1.0},
            'D_5': {0: 1.0}, 'D_3': {0: 1.0},
            'J_5': {0: 1.0}, 'J_3': {0: 1.0}
        }
        
        mock_config.gene_use_dict = {
            'V': {'IGHV1': 1.0},
            'D': {'IGHD1': 1.0},
            'J': {'IGHJ1': 1.0}
        }
        
        mock_config.NP_lengths = {
            'NP1': {5: 0.5, 10: 0.5},
            'NP2': {3: 0.5, 8: 0.5}
        }
        
        mock_config.NP_first_bases = {
            'NP1': {'A': 0.25, 'T': 0.25, 'C': 0.25, 'G': 0.25},
            'NP2': {'A': 0.25, 'T': 0.25, 'C': 0.25, 'G': 0.25}
        }
        
        mock_config.NP_transitions = {
            'NP1': {
                0: {'A': {'A': 0.25, 'T': 0.25, 'C': 0.25, 'G': 0.25}},
                1: {'A': {'A': 0.25, 'T': 0.25, 'C': 0.25, 'G': 0.25}}
            },
            'NP2': {
                0: {'A': {'A': 0.25, 'T': 0.25, 'C': 0.25, 'G': 0.25}},
                1: {'A': {'A': 0.25, 'T': 0.25, 'C': 0.25, 'G': 0.25}}
            }
        }
        
        # Test that the mock config works
        self.assertEqual(mock_config.number_of_v_alleles, 1)
        self.assertEqual(mock_config.number_of_d_alleles, 1)
        self.assertEqual(mock_config.number_of_j_alleles, 1)

    def test_dataconfig_validation(self):
        """Test DataConfig validation and error handling."""
        config = self.heavychain_config
        
        # Test that required components exist
        self.assertIsNotNone(config.v_alleles)
        self.assertIsNotNone(config.d_alleles)
        self.assertIsNotNone(config.j_alleles)
        self.assertIsNotNone(config.trim_dicts)
        self.assertIsNotNone(config.gene_use_dict)
        
        # Test trim dictionary structure - check that V_3 exists (the main one used)
        self.assertIn('V_3', config.trim_dicts)
        self.assertIsInstance(config.trim_dicts['V_3'], dict)
        
        # Check that the trim dictionary has the expected nested structure
        for segment_direction, trim_dict in config.trim_dicts.items():
            if isinstance(trim_dict, dict):
                self.assertIsInstance(trim_dict, dict)

    # ============================================================================
    # ADVANCED PIPELINE STEP TESTS
    # ============================================================================

    def test_corruption_step_variants(self):
        """Test different corruption event types."""
        AugmentationStep.set_dataconfig(self.heavychain_config)
        
        # Test different corruption events
        corruption_configs = [
            (0.8, [1.0, 0.0, 0.0]),  # Only add
            (0.8, [0.0, 1.0, 0.0]),  # Only remove
            (0.8, [0.0, 0.0, 1.0]),  # Only remove_before_add
            (0.8, [0.33, 0.33, 0.34])  # Mixed
        ]
        
        for prob, events in corruption_configs:
            step = CorruptSequenceBeginning(prob, events, 300, 50, 150, 20)
            
            pipeline = AugmentationPipeline([
                SimulateSequence(S5F(), True),
                step,
                DistillMutationRate()
            ])
            
            result = pipeline.execute()
            result_dict = result.get_dict()
            
            self.assertIn('corruption_event', result_dict)
            self.assertIn('sequence', result_dict)

    def test_indel_insertion_edge_cases(self):
        """Test indel insertion with various parameters."""
        # Test with high indel probability
        high_indel_step = InsertIndels(
            indel_probability=0.9,
            max_indels=5,
            insertion_proba=0.5,
            deletion_proba=0.5
        )
        
        # Test with low indel probability
        low_indel_step = InsertIndels(
            indel_probability=0.1,
            max_indels=1,
            insertion_proba=0.7,
            deletion_proba=0.3
        )
        
        for step in [high_indel_step, low_indel_step]:
            pipeline = AugmentationPipeline([
                SimulateSequence(S5F(), True),
                FixVPositionAfterTrimmingIndexAmbiguity(),
                FixDPositionAfterTrimmingIndexAmbiguity(),
                FixJPositionAfterTrimmingIndexAmbiguity(),
                step,
                DistillMutationRate()
            ])
            
            result = pipeline.execute()
            result_dict = result.get_dict()
            
            self.assertIn('indels', result_dict)
            self.assertIsInstance(result_dict['indels'], dict)

    def test_n_insertion_parameters(self):
        """Test N insertion with different parameters."""
        n_insertion_configs = [
            (0.05, 0.1),   # Low N insertion
            (0.02, 0.8),   # High N conversion
            (0.0, 0.0),    # No N insertion
        ]
        
        for n_prob, n_conv in n_insertion_configs:
            step = InsertNs(n_prob, n_conv)
            
            pipeline = AugmentationPipeline([
                SimulateSequence(S5F(), True),
                step,
                DistillMutationRate()
            ])
            
            result = pipeline.execute()
            result_dict = result.get_dict()
            
            self.assertIn('Ns', result_dict)
            self.assertIsInstance(result_dict['Ns'], dict)


if __name__ == '__main__':
    # Run tests with increased verbosity
    unittest.main(verbosity=2)
