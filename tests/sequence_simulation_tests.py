import unittest

from GenAIRR.alleles import VAllele
from GenAIRR.sequence import LightChainSequence, HeavyChainSequence
import pickle
from importlib import resources
import GenAIRR

with resources.path('GenAIRR.data', 'LightChain_KAPPA_DataConfigV2.pkl') as data_path:
    with open(data_path, 'rb') as h:
        lightchain_kappa_config = pickle.load(h)

with resources.path('GenAIRR.data', 'LightChain_LAMBDA_DataConfigV2.pkl') as data_path:
    with open(data_path, 'rb') as h:
        lightchain_lambda_config = pickle.load(h)

with resources.path('GenAIRR.data', 'HeavyChain_DataConfig_OGRDB_V2.pkl') as data_path:
    with open(data_path, 'rb') as h:
        heavychain_config = pickle.load(h)
        print(heavychain_config)


# Define a test case class inheriting from unittest.TestCase
class TestSequenceSimulation(unittest.TestCase):

    # Define a setup method if you need to prepare anything before each test method (optional)
    def setUp(self):
        # Setup code here, e.g., initializing objects
        pass

    # Define a teardown method if you need to clean up after each test method (optional)
    def tearDown(self):
        # Teardown code here, e.g., closing files or connections
        pass

    def test_v_allele_creation(self):
        seq = 'caggtcaccttgaaggagtctggtcct...gtgctggtgaaacccacagagaccctcacgctgacctgcaccgtctctgggttctcactcagc......aatgctagaatgggtgtgagctggatccgtcagcccccagggaaggccctggagtggcttgcacacattttttcgaat.........gacgaaaaatcctacagcacatctctgaag...agcaggctcaccatctccaaggacacctccaaaagccaggtggtccttaccatgaccaatatggaccctgtggacacagccacatattactgtgcatggatac'
        name = 'IGHVF1-G1*01'
        length = 301

        allele = VAllele(name,seq,length)

        self.assertEqual(allele.gene,'IGHVF1-G1')
        self.assertEqual(allele.family,'IGHVF1')
        self.assertEqual(allele.anchor,288)


    # Define your test methods, each starting with 'test_'
    def test_single_sequence_generation(self):
        # Call the function you want to test with the inputs you want to test
        lc1_seq = LightChainSequence.create_random(lightchain_lambda_config)
        lc2_seq = LightChainSequence.create_random(lightchain_kappa_config)
        # print(seq)
        h_seq = HeavyChainSequence.create_random(heavychain_config)

        # Use self.assertEqual(), self.assertTrue(), etc., to assert the expected outcomes
        self.assertIsNotNone(lc1_seq)
        self.assertIsNotNone(lc2_seq)
        self.assertIsNotNone(h_seq)

    def test_lighchain_sequence_simulator(self):
        from GenAIRR.simulation import LightChainKappaLambdaSequenceAugmentor, SequenceAugmentorArguments
        args = SequenceAugmentorArguments()

        aug = LightChainKappaLambdaSequenceAugmentor(lambda_args=args,
                                                     kappa_args=args,
                                                     lambda_dataconfig=lightchain_lambda_config,
                                                     kappa_dataconfig=lightchain_kappa_config)
        generated_seqs = []
        for _ in range(100):
            generated_seqs.append(aug.simulate_augmented_sequence())
        self.assertEqual(len(generated_seqs), 100)

    def test_heavychain_sequence_simulator(self):
        from GenAIRR.simulation import HeavyChainSequenceAugmentor, SequenceAugmentorArguments
        args = SequenceAugmentorArguments(simulate_indels=True)

        aug = HeavyChainSequenceAugmentor(heavychain_config, args)
        generated_seqs = []
        for _ in range(100):
            generated_seqs.append(aug.simulate_augmented_sequence())
        self.assertEqual(len(generated_seqs), 100)


# This conditional ensures that the tests will run only if this script is executed, not when it's imported
if __name__ == '__main__':
    unittest.main()
