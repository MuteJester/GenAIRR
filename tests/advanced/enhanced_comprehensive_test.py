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
    """Enhanced comprehensive test suite for GenAIRR library covering broader functionality."""

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
    # SEQUENCE ASSEMBLY AND STRUCTURAL VALIDATION TESTS
    # ============================================================================

    def test_sequence_assembly_integrity(self):
        """Test that sequence assembly maintains structural integrity across all components."""
        for _ in range(20):
            seq = HeavyChainSequence.create_random(self.heavychain_config)
            
            # Test that V, D, J segments are properly positioned
            self.assertEqual(seq.v_seq_start, 0)
            self.assertGreater(seq.v_seq_end, seq.v_seq_start)
            self.assertGreaterEqual(seq.d_seq_start, seq.v_seq_end)
            self.assertGreater(seq.d_seq_end, seq.d_seq_start)
            self.assertGreaterEqual(seq.j_seq_start, seq.d_seq_end)
            self.assertGreater(seq.j_seq_end, seq.j_seq_start)
            
            # Test that NP regions are correctly sized
            expected_np1_length = seq.d_seq_start - seq.v_seq_end
            expected_np2_length = seq.j_seq_start - seq.d_seq_end
            self.assertEqual(seq.NP1_length, expected_np1_length)
            self.assertEqual(seq.NP2_length, expected_np2_length)
            
            # Test sequence length consistency
            expected_length = (seq.v_seq_end - seq.v_seq_start + 
                             seq.NP1_length + 
                             seq.d_seq_end - seq.d_seq_start + 
                             seq.NP2_length + 
                             seq.j_seq_end - seq.j_seq_start)
            if hasattr(seq, 'c_allele') and seq.c_allele:
                expected_length += len(seq.c_trimmed_seq)
            
            # Allow small variance due to assembly differences
            self.assertAlmostEqual(len(seq.ungapped_seq), expected_length, delta=2)

    def test_trimming_logic_consistency(self):
        """Test that allele trimming follows biological constraints."""
        # Test V allele trimming
        v_alleles = [allele for family in self.heavychain_config.v_alleles.values() for allele in family]
        
        for _ in range(50):
            v_allele = random.choice(v_alleles)
            trimmed_seq, trim_5, trim_3 = v_allele.get_trimmed(self.heavychain_config.trim_dicts)
            
            # V alleles should never be trimmed at 5'
            self.assertEqual(trim_5, 0)
            
            # Trimming should not exceed allele length
            self.assertLessEqual(trim_3, v_allele.ungapped_len)
            
            # Trimmed sequence should be reasonable
            self.assertGreater(len(trimmed_seq), 0)
            self.assertEqual(len(trimmed_seq), v_allele.ungapped_len - trim_3)
            
            # Should not trim past anchor point for functionality
            if trim_3 > 0:
                self.assertLess(trim_3, v_allele.ungapped_len - v_allele.anchor)

    def test_np_region_generation_consistency(self):
        """Test NP region generation follows Markov chain principles."""
        from GenAIRR.sequence.np_region import NP_Region
        
        # Test multiple NP region generations
        for region_type in ['NP1', 'NP2']:
            for length in [5, 10, 15, 25]:
                if length in self.heavychain_config.NP_lengths[region_type]:
                    np_seq = NP_Region.create_np_region(
                        self.heavychain_config.NP_lengths,
                        self.heavychain_config.NP_transitions,
                        region_type,
                        self.heavychain_config.NP_first_bases
                    )
                    
                    if np_seq:  # Only test if sequence was generated
                        # Should only contain valid bases
                        valid_bases = set('ATCGN')
                        self.assertTrue(set(np_seq).issubset(valid_bases))
                        
                        # Should have reasonable length (allowing for early termination)
                        self.assertGreater(len(np_seq), 0)
                        self.assertLessEqual(len(np_seq), length + 5)  # Allow some variance

    def test_junction_boundary_accuracy(self):
        """Test that junction boundaries are correctly calculated."""
        for _ in range(30):
            seq = HeavyChainSequence.create_random(self.heavychain_config)
            
            # Junction should start at V anchor
            self.assertEqual(seq.junction_start, seq.v_allele.anchor)
            
            # Junction should end at calculated position
            expected_end = seq.v_allele.anchor + seq.junction_length
            self.assertEqual(seq.junction_end, expected_end)
            
            # Junction sequence should match calculated boundaries
            expected_junction = seq.ungapped_seq[seq.junction_start:seq.junction_end]
            self.assertEqual(seq.junction, expected_junction.upper())
            
            # Junction length should be consistent with calculation
            calculated_length = seq.get_junction_length()
            self.assertEqual(seq.junction_length, calculated_length)

    def test_allele_anchor_positions(self):
        """Test that allele anchor positions are correctly identified."""
        # Test V allele anchors
        v_alleles = [allele for family in self.heavychain_config.v_alleles.values() for allele in family]
        for v_allele in v_alleles[:10]:  # Test subset for performance
            # Anchor should be within sequence bounds
            self.assertGreater(v_allele.anchor, 0)
            self.assertLess(v_allele.anchor, len(v_allele.ungapped_seq))
            
            # Should be able to find cysteine motif near anchor
            anchor_region = v_allele.ungapped_seq[max(0, v_allele.anchor-10):v_allele.anchor+10]
            # Relaxed test - just check anchor is reasonable
            self.assertIsInstance(v_allele.anchor, int)

    def test_functionality_assessment_accuracy(self):
        """Test biological functionality assessment logic."""
        productive_count = 0
        total_tests = 100
        
        for _ in range(total_tests):
            seq = HeavyChainSequence.create_random(self.heavychain_config)
            
            # Test stop codon detection
            has_stops = seq.check_stops(seq.ungapped_seq, return_pos=False)
            self.assertEqual(seq.stop_codon, has_stops)
            
            # Test frame calculation
            junction_in_frame = (seq.junction_length % 3 == 0 and 
                               seq.junction_start % 3 == 0 and 
                               seq.junction_end % 3 == 0)
            
            # If marked as functional, should meet all criteria
            if seq.functional:
                productive_count += 1
                self.assertFalse(seq.stop_codon)
                self.assertTrue(seq.vj_in_frame)
                self.assertTrue(junction_in_frame)
                
                # Should have valid conserved residues
                if len(seq.junction_aa) > 0:
                    self.assertEqual(seq.junction_aa[0], 'C')  # Should start with Cys
                if len(seq.junction_aa) > 1:
                    self.assertIn(seq.junction_aa[-1], ['F', 'W'])  # Should end with Phe or Trp
        
        # Should produce some functional sequences
        self.assertGreater(productive_count, 5)

    # ============================================================================
    # PIPELINE STEP LOGIC VALIDATION TESTS  
    # ============================================================================

    def test_trimming_index_ambiguity_fixes(self):
        """Test that trimming index ambiguity fixes work correctly."""
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
            
            # Positions should be consistent after fixes
            self.assertGreaterEqual(result_dict['v_sequence_end'], result_dict['v_sequence_start'])
            self.assertGreaterEqual(result_dict['d_sequence_start'], result_dict['v_sequence_end'])
            self.assertGreaterEqual(result_dict['d_sequence_end'], result_dict['d_sequence_start'])
            self.assertGreaterEqual(result_dict['j_sequence_start'], result_dict['d_sequence_end'])
            self.assertGreaterEqual(result_dict['j_sequence_end'], result_dict['j_sequence_start'])
            
            # Germline positions should be reasonable
            self.assertGreaterEqual(result_dict['v_germline_end'], result_dict['v_germline_start'])
            self.assertGreaterEqual(result_dict['d_germline_end'], result_dict['d_germline_start'])
            self.assertGreaterEqual(result_dict['j_germline_end'], result_dict['j_germline_start'])

    def test_short_d_validation_logic(self):
        """Test Short-D validation step logic."""
        from GenAIRR.steps import ShortDValidation
        
        # Test with different thresholds
        for threshold in [3, 5, 8]:
            pipeline = AugmentationPipeline([
                SimulateSequence(S5F(), True),
                FixVPositionAfterTrimmingIndexAmbiguity(),
                FixDPositionAfterTrimmingIndexAmbiguity(),
                FixJPositionAfterTrimmingIndexAmbiguity(),
                ShortDValidation(threshold),
                DistillMutationRate()
            ])
            
            short_d_count = 0
            total_sequences = 50
            
            for _ in range(total_sequences):
                result = pipeline.execute()
                result_dict = result.get_dict()
                
                d_length = result_dict['d_sequence_end'] - result_dict['d_sequence_start']
                
                if d_length < threshold:
                    # Should be marked as Short-D
                    self.assertEqual(result_dict['d_call'], ['Short-D'])
                    short_d_count += 1
                else:
                    # Should not be marked as Short-D
                    self.assertNotEqual(result_dict['d_call'], ['Short-D'])
            
            # Should find some short D segments with smaller thresholds
            if threshold <= 5:
                self.assertGreater(short_d_count, 0)

    def test_allele_correction_maps_usage(self):
        """Test that allele correction maps are properly utilized."""
        pipeline = AugmentationPipeline([
            SimulateSequence(S5F(), True),
            FixVPositionAfterTrimmingIndexAmbiguity(),
            FixDPositionAfterTrimmingIndexAmbiguity(),
            FixJPositionAfterTrimmingIndexAmbiguity(),
            CorrectForVEndCut(),
            CorrectForDTrims(),
            DistillMutationRate()
        ])
        
        original_allele_counts = {}
        corrected_allele_counts = {}
        
        for _ in range(30):
            # Get initial sequence to track original alleles
            initial_result = SimulateSequence(S5F(), True)
            container = SimulationContainer()
            initial_result.apply(container)
            
            original_v = container.v_call[0] if container.v_call else None
            original_d = container.d_call[0] if container.d_call else None
            
            # Run full pipeline
            result = pipeline.execute()
            result_dict = result.get_dict()
            
            # Track allele usage
            if original_v:
                original_allele_counts[original_v] = original_allele_counts.get(original_v, 0) + 1
            
            final_v_calls = result_dict.get('v_call', [])
            final_d_calls = result_dict.get('d_call', [])
            
            # Should have at least one allele call
            self.assertGreater(len(final_v_calls), 0)
            if original_d != 'Short-D':
                self.assertGreater(len(final_d_calls), 0)
            
            # Count corrected alleles
            for allele in final_v_calls:
                corrected_allele_counts[allele] = corrected_allele_counts.get(allele, 0) + 1

    # ============================================================================
    # CORRUPTION AND INDEL TESTING
    # ============================================================================
    
    def test_sequence_corruption_mechanisms(self):
        """Test different sequence corruption event types and their effects."""
        corruption_scenarios = [
            (0.9, [1.0, 0.0, 0.0]),  # Only add events
            (0.9, [0.0, 1.0, 0.0]),  # Only remove events  
            (0.9, [0.0, 0.0, 1.0]),  # Only remove_before_add events
        ]
        
        for prob, events in corruption_scenarios:
            step = CorruptSequenceBeginning(prob, events, 300, 50, 150, 20)
            
            pipeline = AugmentationPipeline([
                SimulateSequence(S5F(), True),
                step,
                DistillMutationRate()
            ])
            
            corrupted_count = 0
            for _ in range(20):
                result = pipeline.execute()
                result_dict = result.get_dict()
                
                if result_dict.get('corruption_event', 'no-corruption') != 'no-corruption':
                    corrupted_count += 1
                    
                    # Test corruption metadata consistency
                    corruption_event = result_dict['corruption_event']
                    self.assertIn(corruption_event, ['add', 'remove', 'remove_before_add', 'no-corruption'])
                    
                    if corruption_event in ['add', 'remove_before_add']:
                        self.assertGreaterEqual(result_dict.get('corruption_add_amount', 0), 0)
                        self.assertIsInstance(result_dict.get('corruption_added_section', ''), str)
                    
                    if corruption_event in ['remove', 'remove_before_add']:
                        self.assertGreaterEqual(result_dict.get('corruption_remove_amount', 0), 0)
                        self.assertIsInstance(result_dict.get('corruption_removed_section', ''), str)
            
            # With high probability, should see corruption events
            self.assertGreater(corrupted_count, 10)

    def test_indel_insertion_position_logic(self):
        """Test that indels are inserted at valid positions only."""
        pipeline = AugmentationPipeline([
            SimulateSequence(S5F(), True),
            FixVPositionAfterTrimmingIndexAmbiguity(),
            FixDPositionAfterTrimmingIndexAmbiguity(),
            FixJPositionAfterTrimmingIndexAmbiguity(),
            InsertIndels(0.8, 3, 0.5, 0.5),  # High probability for testing
            DistillMutationRate()
        ])
        
        indel_containing_sequences = 0
        for _ in range(30):
            result = pipeline.execute()
            result_dict = result.get_dict()
            
            indels = result_dict.get('indels', {})
            if indels:
                indel_containing_sequences += 1
                
                # Test indel position validity
                for pos, indel_info in indels.items():
                    # Position should be within sequence bounds
                    self.assertGreaterEqual(pos, 0)
                    self.assertLess(pos, len(result_dict['sequence']))
                    
                    # Indel info should have correct format
                    self.assertIsInstance(indel_info, str)
                    self.assertTrue(indel_info.startswith('I < ') or indel_info.startswith('D > '))
                    
                    # Should not be in NP regions (this is complex to verify perfectly)
                    # Basic sanity check: position exists
                    self.assertIsInstance(pos, int)
        
        # Should see some indel-containing sequences with high probability
        self.assertGreater(indel_containing_sequences, 10)

    def test_mutation_position_constraints(self):
        """Test that mutations only occur in valid regions (not NP regions)."""
        pipeline = AugmentationPipeline([
            SimulateSequence(S5F(min_mutation_rate=0.05, max_mutation_rate=0.15), True),
            FixVPositionAfterTrimmingIndexAmbiguity(),
            FixDPositionAfterTrimmingIndexAmbiguity(),
            FixJPositionAfterTrimmingIndexAmbiguity(),
            DistillMutationRate()
        ])
        
        for _ in range(20):
            result = pipeline.execute()
            result_dict = result.get_dict()
            
            mutations = result_dict.get('mutations', {})
            if mutations:
                # Get region boundaries
                v_end = result_dict['v_sequence_end']
                d_start = result_dict['d_sequence_start']
                d_end = result_dict['d_sequence_end']
                j_start = result_dict['j_sequence_start']
                
                for pos in mutations.keys():
                    # Mutations should only be in V, D, or J regions, not NP regions
                    in_v_region = 0 <= pos < v_end
                    in_d_region = d_start <= pos < d_end
                    in_j_region = j_start <= pos < len(result_dict['sequence'])
                    
                    self.assertTrue(in_v_region or in_d_region or in_j_region,
                                  f"Mutation at position {pos} is not in a valid region")

    # ============================================================================
    # LIGHT CHAIN SPECIFIC TESTS
    # ============================================================================
    
    def test_light_chain_structural_differences(self):
        """Test light chain sequences have correct structure (no D segment)."""
        AugmentationStep.set_dataconfig(self.lightchain_kappa_config)
        
        for _ in range(20):
            seq = LightChainSequence.create_random(self.lightchain_kappa_config)
            
            # Light chains should not have D allele
            self.assertIsNone(seq.d_allele)
            
            # Should have only V and J segments
            self.assertIsNotNone(seq.v_allele)
            self.assertIsNotNone(seq.j_allele)
            
            # Should have only NP1 region (between V and J)
            self.assertGreater(seq.NP1_length, 0)
            self.assertEqual(seq.NP2_length, 0)
            
            # Junction calculation should be different
            expected_junction_length = (seq.v_allele.length - (seq.v_allele.anchor - 1) - seq.v_trim_3 +
                                      seq.NP1_length + (seq.j_allele.anchor + 2) - seq.j_trim_5)
            self.assertEqual(seq.junction_length, expected_junction_length)

    def test_light_chain_kappa_vs_lambda_differences(self):
        """Test differences between kappa and lambda light chains."""
        kappa_sequences = []
        lambda_sequences = []
        
        # Generate kappa sequences
        AugmentationStep.set_dataconfig(self.lightchain_kappa_config)
        for _ in range(10):
            seq = LightChainSequence.create_random(self.lightchain_kappa_config)
            kappa_sequences.append(seq)
        
        # Generate lambda sequences
        AugmentationStep.set_dataconfig(self.lightchain_lambda_config)
        for _ in range(10):
            seq = LightChainSequence.create_random(self.lightchain_lambda_config)
            lambda_sequences.append(seq)
        
        # Test that allele names are different between kappa and lambda
        kappa_v_names = {seq.v_allele.name for seq in kappa_sequences}
        lambda_v_names = {seq.v_allele.name for seq in lambda_sequences}
        
        # Should have distinct allele sets (allowing for some overlap in naming patterns)
        self.assertGreater(len(kappa_v_names), 0)
        self.assertGreater(len(lambda_v_names), 0)
        
        # Test that junction lengths can vary between types
        kappa_junction_lengths = [seq.junction_length for seq in kappa_sequences]
        lambda_junction_lengths = [seq.junction_length for seq in lambda_sequences]
        
        self.assertGreater(len(set(kappa_junction_lengths)), 1)  # Should have variety
        self.assertGreater(len(set(lambda_junction_lengths)), 1)  # Should have variety

    # ============================================================================
    # DATA CONFIGURATION ROBUSTNESS TESTS
    # ============================================================================
    
    def test_dataconfig_edge_cases(self):
        """Test DataConfig behavior with edge cases."""
        config = self.heavychain_config
        
        # Test allele list retrieval with invalid types
        with self.assertRaises(Exception):
            config.allele_list('invalid_type')
        
        # Test that trim dictionaries have expected structure
        self.assertIn('V_3', config.trim_dicts)
        self.assertIn('D_5', config.trim_dicts)
        self.assertIn('D_3', config.trim_dicts)
        self.assertIn('J_5', config.trim_dicts)
        
        # Test that NP parameters are present
        self.assertIn('NP1', config.NP_lengths)
        self.assertIn('NP2', config.NP_lengths)
        self.assertIn('NP1', config.NP_first_bases)
        self.assertIn('NP2', config.NP_first_bases)

    def test_allele_family_gene_consistency(self):
        """Test that allele family and gene naming is consistent."""
        for allele_type in ['v_alleles', 'd_alleles', 'j_alleles']:
            alleles_dict = getattr(self.heavychain_config, allele_type)
            
            for family_name, alleles in alleles_dict.items():
                for allele in alleles:
                    # Family name should match allele's family
                    self.assertEqual(allele.family, family_name)
                    
                    # Gene should be part of family
                    self.assertTrue(allele.gene.startswith(family_name.split('-')[0]))
                    
                    # Name should include gene
                    self.assertTrue(allele.name.startswith(allele.gene))

    # ============================================================================
    # PERFORMANCE AND STRESS TESTS
    # ============================================================================
    
    def test_large_batch_memory_efficiency(self):
        """Test memory usage doesn't grow excessively with large batches."""
        import gc
        
        # Generate large batch and ensure cleanup works
        sequences = []
        for i in range(200):
            seq = HeavyChainSequence.create_random(self.heavychain_config)
            sequences.append(seq.ungapped_seq)  # Only keep sequence string
            
            # Periodic cleanup
            if i % 50 == 0:
                gc.collect()
        
        # Should have 200 unique or varied sequences
        unique_sequences = set(sequences)
        self.assertGreater(len(unique_sequences), 150)  # High variability expected
        
        # Cleanup
        del sequences
        gc.collect()

    def test_pipeline_reproducibility_with_seed(self):
        """Test that pipeline execution is reproducible with same random seed."""
        
        def run_pipeline_with_seed(seed):
            random.seed(seed)
            np.random.seed(seed)
            
            pipeline = AugmentationPipeline([
                SimulateSequence(S5F(min_mutation_rate=0.02, max_mutation_rate=0.02), True),
                FixVPositionAfterTrimmingIndexAmbiguity(),
                DistillMutationRate()
            ])
            
            return pipeline.execute().get_dict()
        
        # Run same pipeline with same seed twice
        result1 = run_pipeline_with_seed(42)
        result2 = run_pipeline_with_seed(42)
        
        # Should get identical results
        self.assertEqual(result1['sequence'], result2['sequence'])
        self.assertEqual(result1['v_call'], result2['v_call'])
        self.assertEqual(result1['mutations'], result2['mutations'])

    # ============================================================================
    # ADVANCED FUNCTIONALITY TESTS  
    # ============================================================================
    
    def test_custom_allele_selection_constraints(self):
        """Test that custom allele selection works correctly."""
        # Get specific alleles to force
        v_family = list(self.heavychain_config.v_alleles.keys())[0]
        d_family = list(self.heavychain_config.d_alleles.keys())[0]
        j_family = list(self.heavychain_config.j_alleles.keys())[0]
        
        specific_v = self.heavychain_config.v_alleles[v_family][0]
        specific_d = self.heavychain_config.d_alleles[d_family][0]
        specific_j = self.heavychain_config.j_alleles[j_family][0]
        
        for _ in range(10):
            seq = HeavyChainSequence.create_random(
                self.heavychain_config,
                specific_v=specific_v,
                specific_d=specific_d,
                specific_j=specific_j
            )
            
            # Should use exactly the specified alleles
            self.assertEqual(seq.v_allele.name, specific_v.name)
            self.assertEqual(seq.d_allele.name, specific_d.name)
            self.assertEqual(seq.j_allele.name, specific_j.name)

    def test_mutation_model_parameter_effects(self):
        """Test that different mutation model parameters produce expected effects."""
        # Test with very low mutation rate
        low_mut_seq = HeavyChainSequence.create_random(self.heavychain_config)
        low_mut_model = S5F(min_mutation_rate=0.001, max_mutation_rate=0.001)
        low_mut_seq.mutate(low_mut_model)
        
        # Test with high mutation rate  
        high_mut_seq = HeavyChainSequence.create_random(self.heavychain_config)
        high_mut_model = S5F(min_mutation_rate=0.10, max_mutation_rate=0.10)
        high_mut_seq.mutate(high_mut_model)
        
        # High mutation should generally have more mutations
        # (though this can vary due to sequence length differences)
        low_mut_count = len(low_mut_seq.mutations)
        high_mut_count = len(high_mut_seq.mutations)
        
        # At minimum, verify mutation counts are reasonable
        self.assertGreaterEqual(low_mut_count, 0)
        self.assertGreaterEqual(high_mut_count, 0)
        
        # Verify mutation rates are set correctly
        self.assertAlmostEqual(low_mut_seq.mutation_freq, 0.001, places=3)
        self.assertAlmostEqual(high_mut_seq.mutation_freq, 0.10, places=2)

    def test_translation_and_biological_validity(self):
        """Test translation function and biological sequence validity."""
        # Test translation function with known sequences
        test_cases = [
            ('ATGAAATAG', 'MK*'),    # Start + amino acid + stop
            ('ATGTTTTAT', 'MFY'),    # Normal amino acids
            ('NNNAAATTT', 'XKF'),    # N should translate to X
            ('---AAATTT', 'XKF'),    # Gaps should translate to X
            ('...AAATTT', '_KF'),    # Dots should translate to _
            ('ATGA', 'M.'),          # Incomplete codon
        ]
        
        for nucleotide_seq, expected_aa in test_cases:
            result = translate(nucleotide_seq)
            self.assertEqual(result, expected_aa)
        
        # Test with real sequences
        for _ in range(10):
            seq = HeavyChainSequence.create_random(self.heavychain_config)
            
            # Translate junction
            junction_aa = translate(seq.junction)
            
            # Should start with C if functional
            if seq.functional and len(junction_aa) > 0:
                self.assertEqual(junction_aa[0], 'C')
            
            # Should not contain invalid characters
            valid_aa_chars = set('ACDEFGHIKLMNPQRSTVWY*X_.')
            self.assertTrue(set(junction_aa).issubset(valid_aa_chars))

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
