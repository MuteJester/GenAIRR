import unittest
import random

import numpy as np

from GenAIRR.alleles import VAllele
from GenAIRR.generateDataConfig import RandomDataConfigGenerator, CustomDataConfigGenerator
from GenAIRR.sequence import LightChainSequence, HeavyChainSequence
from GenAIRR.data import builtin_heavy_chain_data_config,builtin_kappa_chain_data_config,builtin_lambda_chain_data_config
lightchain_kappa_config = builtin_kappa_chain_data_config()
lightchain_lambda_config = builtin_lambda_chain_data_config()
heavychain_config = builtin_heavy_chain_data_config()

class TestSequenceSimulation(unittest.TestCase):

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

    def test_fix_v_position_after_trimming_index_ambiguity(self):
        # CREATE a sequence that has a junction that recreates a part of the reference alleles
        # test wheter the index was correctly updated for the v
        from GenAIRR.simulation import HeavyChainSequenceAugmentor, SequenceAugmentorArguments
        args = SequenceAugmentorArguments(simulate_indels=0.2)
        augmentor = HeavyChainSequenceAugmentor(heavychain_config, args)
        simulated = {
            'sequence': 'TGTNCTTCCAGTACACCTGAAGGGTGGGTACTCGTTACCCTTCTCGCTCTNTTAGGATATTGGGTCAGCAAGTGGATGACAAAGAGGGATGNACTCAGATTGGCGTGTAGTGGAGGAGGTGCAGCTGGTGGAGTCTGGGGGANGCTTGGTACAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTCAGTAGCTATAGCATGAACTNGGTCCGCCAGGCTCNAGGGAAGGGTCTGGAGTGGGTTTCATACATCAGTAGTAGTAGTAATAGCATATACTACGCAGACTCTGTNAAGGGCCGATTCACCATCTCCAGAGACAANGCCAAGAACTCACTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCTATNTATTACTGTGCGAGAGATGTNCTTCCAGTACACCTGAAGGGTGGGTACTCGTTACCCTTCTCGCTCTNTTAGGATATTGGGTCAGCAAGTGGATGACAAAGAGGGATGNACTCAGATTGGCGTGTAGTGGAGGAGGTGCAGCTGGTGGAGTCTGGGGGANGCTTGGTACAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTCAGTAGCTATAGCATGAACTNGGTCCGCCAGGCTCNAGGGAAGGGTCTGGAGTGGGTTTCATACATCAGTAGTAGTAGTAATAGCATATACTACG',
            'v_sequence_start': 115,
            'v_sequence_end': 410,
            'd_sequence_start': 421,
            'd_sequence_end': 431,
            'j_sequence_start': 445,
            'j_sequence_end': 494,
            'v_germline_start': 0,
            'v_germline_end': 295,
            'd_germline_start': 4,
            'd_germline_end': 14,
            'j_germline_start': 2,
            'j_germline_end': 51,
            'junction_sequence_start': 400,
            'junction_sequence_end': 463,
            'v_call': ['IGHVF10-G52*05','IGHVF10-G52*04'],
            'v_trim_5': 0,
            'v_trim_3': 1,
            'd_trim_5': 4,
            'd_trim_3': 5,
            'j_trim_5': 2,
            'j_trim_3': 0,
            'c_trim_3': 25,}
        augmentor.fix_v_position_after_trimming_index_ambiguity(simulated)

        self.assertEqual(simulated['v_trim_3'],0)
        self.assertEqual(simulated['v_sequence_end'],411)
        self.assertEqual(simulated['v_germline_end'],296)
    def test_fix_d_position_after_trimming_index_ambiguity(self):
        # CREATE a sequence that has a junction that recreates a part of the reference alleles
        # test wheter the index was correctly updated for the v
        from GenAIRR.simulation import HeavyChainSequenceAugmentor, SequenceAugmentorArguments
        args = SequenceAugmentorArguments(simulate_indels=0.2)
        augmentor = HeavyChainSequenceAugmentor(heavychain_config, args)
        simulated = {
            'sequence': 'TGTNCTTCCAGTACACCTGAAGGGTGGGTACTCGTTACCCTTCTCGCTCTNTTAGGATATTGGGTCAGCAAGTGGATGACAAAGAGGGATGNACTCAGATTGGCGTGTAGTGGAGGAGGTGCAGCTGGTGGAGTCTGGGGGANGCTTGGTACAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTCAGTAGCTATAGCATGAACTNGGTCCGCCAGGCTCNAGGGAAGGGTCTGGAGTGGGTTTCATACATCAGTAGTAGTAGTAATAGCATATACTACGCAGACTCTGTNAAGGGCCGATTCACCATCTCCAGAGACAANGCCAAGAACTCACTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCTATNTATTACTGTGCGAGAGTTCCAGGGCACTACGGTGGTAACAAAACTGGTTNGACCCCTGGGGCCAGTGAACCCTGGTCACCGTCTCCTCAGGCTTCCACCAAGGGC',
            'v_sequence_start': 115,
            'v_sequence_end': 410,
            'd_sequence_start': 421,
            'd_sequence_end': 431,
            'j_sequence_start': 445,
            'j_sequence_end': 494,
            'v_germline_start': 0,
            'v_germline_end': 295,
            'd_germline_start': 4,
            'd_germline_end': 14,
            'j_germline_start': 2,
            'j_germline_end': 51,
            'junction_sequence_start': 400,
            'junction_sequence_end': 463,
            'v_call': 'IGHVF10-G52*05,IGHVF10-G52*04',
            'd_call': ['IGHD4-23*01','IGHD4-23*01'],
            'j_call': 'IGHJ5*02',
            'c_call': 'IGHG3*09',
            'mutation_rate': 0.03143418467583497,
            'v_trim_5': 0,
            'v_trim_3': 1,
            'd_trim_5': 4,
            'd_trim_3': 5,
            'j_trim_5': 2,
            'j_trim_3': 0,
            'c_trim_3': 25,
            'corruption_event': 'add',
            'corruption_add_amount': 115,
            'corruption_remove_amount': 0,
            'mutations': {246: 'G>T',
                          267: 'T>C',
                          281: 'G>A',
                          284: 'C>G',
                          391: 'G>A',
                          469: 'G>T'},
            'Ns': {3: 'A>N',
                   50: 'A>N',
                   91: 'A>N',
                   142: 'G>N',
                   221: 'G>N',
                   236: 'C>N',
                   306: 'G>N',
                   336: 'T>N',
                   393: 'G>N',
                   453: 'C>N'},
            'indels': {},
            'productive': False,
            'stop_codon': True,
            'vj_in_frame': False,
            'note': 'Stop codon present.'}
        augmentor.fix_d_position_after_trimming_index_ambiguity(simulated)

        self.assertEqual(simulated['d_trim_3'],3)
        self.assertEqual(simulated['d_trim_5'],2)
        self.assertEqual(simulated['d_sequence_start'],419)
        self.assertEqual(simulated['d_sequence_end'],433)
        self.assertEqual(simulated['d_germline_start'], 2)
        self.assertEqual(simulated['d_germline_end'], 16)
    def test_fix_j_position_after_trimming_index_ambiguity(self):
        # CREATE a sequence that has a junction that recreates a part of the reference alleles
        # test wheter the index was correctly updated for the v
        from GenAIRR.simulation import HeavyChainSequenceAugmentor, SequenceAugmentorArguments
        args = SequenceAugmentorArguments(simulate_indels=0.2)
        augmentor = HeavyChainSequenceAugmentor(heavychain_config, args)
        simulated = {
            'sequence': 'TGTNCTTCCAGTACACCTGAAGGGTGGGTACTCGTTACCCTTCTCGCTCTNTTAGGATATTGGGTCAGCAAGTGGATGACAAAGAGGGATGNACTCAGATTGGCGTGTAGTGGAGGAGGTGCAGCTGGTGGAGTCTGGGGGANGCTTGGTACAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTCAGTAGCTATAGCATGAACTNGGTCCGCCAGGCTCNAGGGAAGGGTCTGGAGTGGGTTTCATACATCAGTAGTAGTAGTAATAGCATATACTACGCAGACTCTGTNAAGGGCCGATTCACCATCTCCAGAGACAANGCCAAGAACTCACTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCTATNTATTACTGTGCGAGAGTTCCAGGGCATTACGGTGGTAGGGTCCGCCCGGACTACGGTGGTAGGGTCCGCCCGGAAAACTGGTTNGACCCCTGGGGCCAGTGAACCCTGGTCACCGTCTCCTCAGGCTTCCACCAAGGGC',
            'v_sequence_start': 115,
            'v_sequence_end': 410,
            'd_sequence_start': 421,
            'd_sequence_end': 431,
            'j_sequence_start': 445,
            'j_sequence_end': 494,
            'v_germline_start': 0,
            'v_germline_end': 295,
            'd_germline_start': 4,
            'd_germline_end': 14,
            'j_germline_start': 2,
            'j_germline_end': 51,
            'junction_sequence_start': 400,
            'junction_sequence_end': 463,
            'v_call': 'IGHVF10-G52*05,IGHVF10-G52*04',
            'd_call': ['IGHD4-23*01','IGHD4-23*01'],
            'j_call': ['IGHJ5*02'],
            'c_call': 'IGHG3*09',
            'mutation_rate': 0.03143418467583497,
            'v_trim_5': 0,
            'v_trim_3': 1,
            'd_trim_5': 4,
            'd_trim_3': 5,
            'j_trim_5': 2,
            'j_trim_3': 0,
            'c_trim_3': 25,
            'corruption_event': 'add',
            'corruption_add_amount': 115,
            'corruption_remove_amount': 0,
            'mutations': {246: 'G>T',
                          267: 'T>C',
                          281: 'G>A',
                          284: 'C>G',
                          391: 'G>A',
                          469: 'G>T'},
            'Ns': {3: 'A>N',
                   50: 'A>N',
                   91: 'A>N',
                   142: 'G>N',
                   221: 'G>N',
                   236: 'C>N',
                   306: 'G>N',
                   336: 'T>N',
                   393: 'G>N',
                   453: 'C>N'},
            'indels': {},
            'productive': False,
            'stop_codon': True,
            'vj_in_frame': False,
            'note': 'Stop codon present.'}
        augmentor.fix_j_position_after_trimming_index_ambiguity(simulated)

        self.assertEqual(simulated['j_trim_3'],0)
        self.assertEqual(simulated['j_trim_5'],0)
        self.assertEqual(simulated['j_sequence_start'],443)
        self.assertEqual(simulated['j_germline_start'],0)

    def test_correct_for_v_end_cut(self):
        # Heavy Chain
        from GenAIRR.simulation import HeavyChainSequenceAugmentor, SequenceAugmentorArguments
        args = SequenceAugmentorArguments(simulate_indels=0.2)
        augmentor = HeavyChainSequenceAugmentor(heavychain_config, args)
        # choose one random v allele
        matches = []
        N = 50
        for _ in range(N): # repeat this N times each time random scenario
            v_allele = random.choice(augmentor.v_alleles)
            v_seq = v_allele.ungapped_seq
            v_seq_end = v_allele.length
            random_trim_end_trim = np.random.randint(1, v_seq_end)
            trimmed_seq = v_seq[:-random_trim_end_trim]
            ambig = []
            for allele in augmentor.v_alleles:
                if trimmed_seq in allele.ungapped_seq:
                    ambig.append(allele.name)

            simulated = {'v_call':[v_allele.name],'v_sequence_end':v_seq_end-random_trim_end_trim,
                         'v_trim_3':random_trim_end_trim}
            augmentor.correct_for_v_end_cut(simulated)
            #print(simulated)
            #print(ambig)
            same_alleles = set(ambig)&set(simulated['v_call'])

            matches.append(len(same_alleles) == len(ambig))


        # Light Chain
        from GenAIRR.simulation import LightChainKappaLambdaSequenceAugmentor, SequenceAugmentorArguments
        args = SequenceAugmentorArguments(simulate_indels=0.2)
        augmentor = LightChainKappaLambdaSequenceAugmentor(kappa_dataconfig=lightchain_kappa_config,
        lambda_dataconfig=lightchain_lambda_config, kappa_args=args,lambda_args=args)
        # choose one random v allele
        for _ in range(N):  # repeat this N times each time random scenario
            v_alleles = augmentor.kappa_augmentor.v_alleles+augmentor.lambda_augmentor.v_alleles
            v_allele = random.choice(v_alleles)
            if 'IGK' in v_allele.name:
                mode = 'kappa'
                v_alleles = augmentor.kappa_augmentor.v_alleles

            else:
                mode = 'lambda'
                v_alleles = augmentor.lambda_augmentor.v_alleles

            v_seq = v_allele.ungapped_seq
            v_seq_end = v_allele.length
            random_trim_end_trim = np.random.randint(1, v_seq_end)
            trimmed_seq = v_seq[:-random_trim_end_trim]
            ambig = []
            for allele in v_alleles:
                if trimmed_seq in allele.ungapped_seq:
                    ambig.append(allele.name)

            simulated = {'v_call': [v_allele.name], 'v_sequence_end': v_seq_end - random_trim_end_trim,
                         'v_trim_3': random_trim_end_trim}

            if mode == 'kappa':
                augmentor.kappa_augmentor.correct_for_v_end_cut(simulated)
            else:
                augmentor.lambda_augmentor.correct_for_v_end_cut(simulated)
            #print(simulated)
            #print(ambig)
            same_alleles = set(ambig) & set(simulated['v_call'])

            matches.append(len(same_alleles) == len(ambig))

        self.assertTrue(sum(matches) == 2*N)

    def test_n_and_removal_ambiguity(self):
        from GenAIRR.simulation import HeavyChainSequenceAugmentor, SequenceAugmentorArguments
        alleles = [j for i in heavychain_config.v_alleles for j in heavychain_config.v_alleles[i]]

        # find single position differing allele pairs
        single_position_diff = []
        for i in alleles:
            for j in alleles:
                if i == j:
                    continue
                else:
                    if i.ungapped_len == j.ungapped_len and sum(
                            [k != z for k, z in zip(i.ungapped_seq, j.ungapped_seq)]) == 1:
                        single_position_diff.append(
                            (i, j, [e for e, (k, z) in enumerate(zip(i.ungapped_seq, j.ungapped_seq)) if k != z][0]))

        if len(single_position_diff) == 0:
            return

        s = 0
        for a,b,p in single_position_diff:
            x = heavychain_config.correction_maps['V_N_AMBIGUITY_CORRECTION_GRAPH'].find_indistinguishable_alleles(
                a.name, [p])
            s = s + (b.name in x and a.name in x)
        self.assertEqual(s,len(single_position_diff))

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
        import base64

        args = SequenceAugmentorArguments(simulate_indels=0.2)

        aug = HeavyChainSequenceAugmentor(heavychain_config, args)
        generated_seqs = []
        for _ in range(100):
            generated_seqs.append(aug.simulate_augmented_sequence())


        self.assertEqual(len(generated_seqs), 100)

    def test_heavychain_sequence_simulator_all_productive(self):
        from GenAIRR.simulation import HeavyChainSequenceAugmentor, SequenceAugmentorArguments
        import base64

        args = SequenceAugmentorArguments(simulate_indels=0.2,productive=True)

        aug = HeavyChainSequenceAugmentor(heavychain_config, args)
        generated_seqs = []
        productive=[]
        for _ in range(100):
            sim = aug.simulate_augmented_sequence()
            generated_seqs.append(sim)
            productive.append(sim['productive'])


        self.assertEqual(len(generated_seqs), 100)
        self.assertEqual(sum(productive), 100)


    def test_tcr_sequence_simulator(self):
        from GenAIRR.TCR.simulation import TCRHeavyChainSequenceAugmentor, SequenceAugmentorArguments
        import base64
        from GenAIRR.data import builtin_tcrb_data_config

        args = SequenceAugmentorArguments(simulate_indels=0.2)

        aug = TCRHeavyChainSequenceAugmentor(builtin_tcrb_data_config(), args)
        generated_seqs = []
        for _ in range(100):
            generated_seqs.append(aug.simulate_augmented_sequence())


        self.assertEqual(len(generated_seqs), 100)

    def test_mutation_rate(self):
        from GenAIRR.simulation import HeavyChainSequenceAugmentor, SequenceAugmentorArguments
        import base64
        args = SequenceAugmentorArguments(simulate_indels=0.2,save_ns_record=True,save_mutations_record=True,
                                          save_corruption_record=True)

        aug = HeavyChainSequenceAugmentor(heavychain_config, args)

        same = 0
        for _ in range(100):
            gseq = aug.simulate_augmented_sequence()
            seq_length = len(gseq['sequence'])
            mutations = len(gseq['mutations'])#len(eval(base64.b64decode(gseq['mutations']).decode('ascii')))
            Ns = len(gseq['Ns'])#len(eval(base64.b64decode(gseq['Ns']).decode('ascii')))
            same += gseq['mutation_rate'] == ((Ns+mutations)/seq_length)
        self.assertEqual(same, 100)

    def test_insert_indels_correct_number(self):
        from GenAIRR.simulation import HeavyChainSequenceAugmentor, SequenceAugmentorArguments
        args = SequenceAugmentorArguments(simulate_indels=0.2)

        augmentor = HeavyChainSequenceAugmentor(heavychain_config, args)
        # Test that the correct number of indels are inserted based on the configured probabilities and limits
        simulated = {'sequence': 'ATCAGCAGTGGTGATTACGATTGGAGCTGGATCCGNCAGCCCCCAGGGAAGGGCCTGGAGTGGATTGGGAACATCTNTTACAGTGGGAGCACCTACTACAACCCGTCCCTCAAGAGTTGAGTTACCATATCATTAGANTNGCCCAAGAACCAGTTCTCCCTGAAACTGAGCTCTGTGACTGCCGCAGACACGGCCGTGTCTTACTGTGTCATCCGGGCCACCGCAGTGTTCGGAGGATCACGGGGTTAAACCTACTGGTACTGCGATCTCTGGGGCCATGGCANCCTGGTCACTGTCTCCTCCG',
                     'v_sequence_start': 0,
                     'v_sequence_end': 211,
                     'd_sequence_start': 221,
                     'd_sequence_end': 228,
                     'j_sequence_start': 251,
                     'j_sequence_end': 304,
                     'junction_sequence_start': -1, #TODO make sure all the test also test for correct movment of  the junction positions
                     'junction_sequence_end': -1,
                     'v_call': 'IGHVF3-G11*06',
                     'd_call': 'IGHD6-19*01',
                     'j_call': 'IGHJ2*01',
                     'mutation_rate': 0.06283912099837222,
                     'v_trim_5': 0,
                     'v_trim_3': 3,
                     'd_trim_5': 6,
                     'd_trim_3': 8,
                     'j_trim_5': 0,
                     'j_trim_3': 0,
                     'corruption_event': 'remove_before_add',
                     'corruption_add_amount': 2,
                     'corruption_remove_amount': 84,
                     'mutations': {18: 'T>G',
                      20: 'C>T',
                      69: 'T>A',
                      117: 'C>T',
                      132: 'G>T',
                      138: 'A>T',
                      141: 'T>C',
                      164: 'G>A',
                      199: 'A>C',
                      208: 'C>T',
                      214: 'T>G',
                      221: 'A>C',
                      237: 'A>T',
                      240: 'A>C',
                      247: 'C>A',
                      262: 'T>G',
                      277: 'G>A',
                      302: 'A>C'},
                     'Ns': {35: 'C>N', 76: 'A>N', 137: 'C>N', 139: 'C>N', 283: 'N>N'},
                     'indels': {}}
        augmentor.max_indels = 5  # Set a limit on the max number of indels for testing
        augmentor.insertion_proba = 0.5
        augmentor.deletion_proba = 0.5
        augmentor.insert_indels(simulated)
        self.assertTrue(0 <= len(simulated['indels']) <= 5)

    def test_simulate_sequence_structure(self):
        from GenAIRR.simulation import HeavyChainSequenceAugmentor, SequenceAugmentorArguments
        args = SequenceAugmentorArguments(simulate_indels=0.2)

        augmentor = HeavyChainSequenceAugmentor(heavychain_config, args)
        # Test that simulate_sequence returns a dictionary with the expected structure and keys
        simulated = augmentor.simulate_sequence()
        expected_keys = ['sequence', 'v_sequence_start', 'v_sequence_end', 'd_sequence_start', 'd_sequence_end',
                         'j_sequence_start', 'j_sequence_end', 'v_call', 'd_call', 'j_call', 'mutation_rate',
                         'v_trim_5', 'v_trim_3', 'd_trim_5', 'd_trim_3', 'j_trim_5', 'j_trim_3', 'corruption_event',
                         'corruption_add_amount', 'corruption_remove_amount', 'mutations', 'Ns', 'indels']
        self.assertTrue(all(key in simulated for key in expected_keys))

    def test_apply_deletion_updates_positions(self):
        from GenAIRR.simulation import HeavyChainSequenceAugmentor, SequenceAugmentorArguments
        args = SequenceAugmentorArguments(simulate_indels=0.2)

        augmentor = HeavyChainSequenceAugmentor(heavychain_config, args)

        simulated = {
            'sequence': 'ATGCGTACGATCG',
            'v_sequence_start': 1,
            'v_sequence_end': 3,
            'd_sequence_start': 4,
            'd_sequence_end': 7,
            'j_sequence_start': 8,
            'j_sequence_end': 12,
            'junction_sequence_start':-1,
            'junction_sequence_end': -1,
            'Ns': {5: 'N'},
            'mutations': {6: 'A'},
            'indels':{}
        }
        augmentor.apply_deletion(simulated, 2)  # Deleting at position 2
        # Test if the start/end positions are updated correctly
        self.assertEqual(simulated['v_sequence_end'], 2)
        self.assertEqual(simulated['d_sequence_start'], 3)
        self.assertEqual(simulated['d_sequence_end'], 6)
        self.assertEqual(simulated['j_sequence_start'], 7)
        self.assertEqual(simulated['j_sequence_end'], 11)
        # Test if N's and Mutations positions are updated correctly
        self.assertNotIn(5, simulated['Ns'])
        self.assertIn(4, simulated['Ns'])
        self.assertNotIn(6, simulated['mutations'])
        self.assertIn(5, simulated['mutations'])

    def test_apply_insertion_updates_positions(self):
        from GenAIRR.simulation import HeavyChainSequenceAugmentor, SequenceAugmentorArguments
        args = SequenceAugmentorArguments(simulate_indels=0.2)

        augmentor = HeavyChainSequenceAugmentor(heavychain_config, args)
        simulated = {
            'sequence': 'ATGCGTACGATCG',
            'v_sequence_start': 1,
            'v_sequence_end': 3,
            'd_sequence_start': 4,
            'd_sequence_end': 7,
            'j_sequence_start': 8,
            'j_sequence_end': 12,
            'junction_sequence_start': -1,
            'junction_sequence_end': -1,
            'Ns': {5: 'N'},
            'mutations': {6: 'A'},
            'indels': {}

        }
        augmentor.apply_insertion(simulated, 2)  # Inserting at position 2
        # Test if the start/end positions are updated correctly
        self.assertEqual(simulated['v_sequence_end'], 4)
        self.assertEqual(simulated['d_sequence_start'], 5)
        self.assertEqual(simulated['d_sequence_end'], 8)
        self.assertEqual(simulated['j_sequence_start'], 9)
        self.assertEqual(simulated['j_sequence_end'], 13)
        # Test if N's and Mutations positions are updated correctly
        self.assertNotIn(5, simulated['Ns'])
        self.assertIn(6, simulated['Ns'])
        self.assertNotIn(6, simulated['mutations'])
        self.assertIn(7, simulated['mutations'])

    def test_multiple_insertions(self):
        from GenAIRR.simulation import HeavyChainSequenceAugmentor, SequenceAugmentorArguments
        args = SequenceAugmentorArguments(simulate_indels=0.2)

        augmentor = HeavyChainSequenceAugmentor(heavychain_config, args)
        original_sequence = 'ATGCGTACGATCG'
        simulated = {
            'sequence': original_sequence,
            'v_sequence_start': 1,
            'v_sequence_end': 3,
            'd_sequence_start': 4,
            'd_sequence_end': 7,
            'j_sequence_start': 8,
            'j_sequence_end': 12,
            'junction_sequence_start': -1,
            'junction_sequence_end': -1,
            'Ns': {5: 'N'},
            'mutations': {6: 'A'},
            'indels':{}
        }
        # ATGCGTACGATCG ACTCGCCGTACGATCG
        # Positions where nucleotides will be inserted
        insert_positions = [2, 4, 1]
        inserted_bases = []
        for idx in range(3):
            augmentor.apply_insertion(simulated,insert_positions[idx])
            inserted_bases.append(simulated['sequence'][insert_positions [idx]])
            for update_idx in range(idx+1, idx + (3 - idx)):
                if insert_positions[update_idx] >= insert_positions[idx]:
                    insert_positions[update_idx] += 1


        self.assertEqual(simulated['v_sequence_end'], 5)
        self.assertEqual(simulated['d_sequence_start'], 7)
        self.assertEqual(simulated['d_sequence_end'], 10)
        self.assertEqual(simulated['j_sequence_start'], 11)
        self.assertEqual(simulated['j_sequence_end'], 15)

        self.assertNotIn(5, simulated['Ns'])
        self.assertNotIn(6, simulated['mutations'])
        self.assertIn(8, simulated['Ns'])
        self.assertIn(9, simulated['mutations'])

        expected_sequence = 'ATGCGTACGATCG'
        expected_sequence = expected_sequence[:2]+inserted_bases[0]+expected_sequence[2:]
        expected_sequence = expected_sequence[:5]+inserted_bases[1]+expected_sequence[5:]
        expected_sequence = expected_sequence[:1]+inserted_bases[2]+expected_sequence[1:]

        self.assertEqual(simulated['sequence'], expected_sequence)

    def test_multiple_deletions(self):
        from GenAIRR.simulation import HeavyChainSequenceAugmentor, SequenceAugmentorArguments
        args = SequenceAugmentorArguments(simulate_indels=0.2)

        augmentor = HeavyChainSequenceAugmentor(heavychain_config, args)
        simulated = {
            'sequence': 'ATGCGTACGATCG',
            'v_sequence_start': 1,
            'junction_sequence_start':-1,
            'junction_sequence_end':-1,
            'v_sequence_end': 3,
            'd_sequence_start': 4,
            'd_sequence_end': 7,
            'j_sequence_start': 8,
            'j_sequence_end': 12,
            'Ns': {5: 'N', 10: 'N'},
            'mutations': {6: 'A', 12: 'G'},
            'indels': {}
        }
        # Positions where nucleotides will be deleted
        delete_positions = [2, 4, 9]  # Choose positions carefully to avoid index out of range errors

        # Apply multiple deletions
        for pos in sorted(delete_positions, reverse=True):  # Reverse order to maintain correct indexes during deletion
            augmentor.apply_deletion(simulated, pos)

        # After deletions, the end positions should have been shifted left by the number of deletions before them
        self.assertEqual(simulated['v_sequence_end'], 3 - sum(pos < 3 for pos in delete_positions))
        self.assertEqual(simulated['d_sequence_start'], 4 - sum(pos < 4 for pos in delete_positions))
        self.assertEqual(simulated['d_sequence_end'], 7 - sum(pos < 7 for pos in delete_positions))
        self.assertEqual(simulated['j_sequence_start'], 8 - sum(pos < 8 for pos in delete_positions))
        self.assertEqual(simulated['j_sequence_end'], 12 - len(delete_positions))  # All deletions affect the j_sequence_end

        # Test if N's and Mutations positions are updated correctly
        # Their new positions should account for the total shift caused by deletions before their original positions
        self.assertNotIn(5, simulated['Ns'])  # Original position 5 should be deleted
        self.assertNotIn(6, simulated['mutations'])  # Original position 6 should be shifted left
        self.assertIn(3, simulated['Ns'])  # New position after deletion at 4
        self.assertIn(9, simulated['mutations'])  # New position after deletions at 2 and 4

        # Ensure the sequence is updated correctly
        expected_sequence = 'ATCTACGTCG'  # Expected sequence after deletions at positions 2, 4, and 9
        self.assertEqual(simulated['sequence'], expected_sequence)

    def test_random_dataconfig_generator(self):
        dcg = RandomDataConfigGenerator(convert_to_asc=False)
        random_dataconfig = dcg.make_dataconfig_from_existing_reference_files(v_reference_path='./IGHV.fasta',
                                                                              d_reference_path='./IGHD.fasta',
                                                                              j_reference_path='./IGHJ.fasta',
                                                                              c_reference_path='./IGHC.fasta')
        for gene in ['V','D','J']:
            self.assertGreater(len(random_dataconfig.gene_use_dict[gene]),0)

        for allele in ['V','D','J']:
            self.assertGreaterEqual(len(random_dataconfig.trim_dicts[allele + '_5']),0)
            self.assertGreaterEqual(len(random_dataconfig.trim_dicts[allele + '_3']),0)

        for np in ['NP1','NP2']:
            self.assertGreaterEqual(len(random_dataconfig.NP_lengths[np]),0)
            self.assertGreaterEqual(len(random_dataconfig.NP_first_bases[np]),0)
            self.assertGreaterEqual(len(random_dataconfig.NP_transitions[np]),0)

    def test_custom_dataconfig_generator(self):
        dcg = CustomDataConfigGenerator(convert_to_asc=False)
        random_dataconfig = dcg.make_dataconfig_from_existing_reference_files(v_reference_path='./IGHV.fasta',
                                                                              d_reference_path='./IGHD.fasta',
                                                                              j_reference_path='./IGHJ.fasta',
                                                                              custom_data='./inference_sample.csv')
        for gene in ['V', 'D', 'J']:
            self.assertGreater(len(random_dataconfig.gene_use_dict[gene]), 0)

        for allele in ['V', 'D', 'J']:
            self.assertGreaterEqual(len(random_dataconfig.trim_dicts[allele + '_5']), 0)
            self.assertGreaterEqual(len(random_dataconfig.trim_dicts[allele + '_3']), 0)

        for np in ['NP1', 'NP2']:
            self.assertGreaterEqual(len(random_dataconfig.NP_lengths[np]), 0)
            self.assertGreaterEqual(len(random_dataconfig.NP_first_bases[np]), 0)
            self.assertGreaterEqual(len(random_dataconfig.NP_transitions[np]), 0)


    #TODO: Fix This Test
    # def test_productive_check(self):
    #     # test the assesment of productive sequence
    #     from GenAIRR.simulation import HeavyChainSequenceAugmentor, SequenceAugmentorArguments
    #     args = SequenceAugmentorArguments(productive=True)
    #     augmentor = HeavyChainSequenceAugmentor(heavychain_config, args)
    #     simulated = {'sequence': 'CAGGTNACTTTGAAGGAGTCTGGTCACTCTGCTGGTGAAACACACAGAGACCGTCACGTTGACCTGCAACGNCTCTGGCTNCTCACTCAGCAAGGCTAGAATGGGTGTGACCTGGCTCCGTCAGTCCCAAGGGAAGNGCCTGGAGTGACTTTTACACGTTTTTTCGAATGACGAATCATAAAGGAGCACGNCTCTGAAGTGAAGTCTCACCATCTCCAAGGACACCTCCAAGAGCCAGGTGGTCGTAACCATGAACAACATGGACCCTGTGGACACAGCCACATATTAGTATGCCGGCTGCCTTCGAGGACTCNCTCATCAGTGAATACTTCCAACATCGGGGCCAGGGCACCCTGGTCACCGTNTCCTCAA',
    #                     'v_sequence_start': 0,
    #                     'v_sequence_end': 301,
    #                     'j_sequence_start': 321,
    #                     'j_sequence_end': 372,
    #                     'v_call': ['IGHVF1-G1*02'],
    #                     'j_call': ['IGHJ1*01'],
    #                     'v_trim_5': 0,
    #                     'v_trim_3': 1,
    #                     'j_trim_5': 1,
    #                     'j_trim_3': 0,
    #                     'corruption_event': 'no-corruption',
    #                     'corruption_add_amount': 0,
    #                     'corruption_remove_amount': 0,
    #                      'productive': False,
    #                      'stop_codon': False,
    #                      'vj_in_frame': False,
    #                      'note': ''
    #                  }
    #     augmentor.fix_productive_call_after_corruption_indel(simulated)
    #     self.assertEqual(simulated['productive'], False) # simulated sequence is not productive
    #     self.assertNotEqual(simulated['stop_codon'], False) # simulated sequence has a stop codon
    #     self.assertEqual(simulated['vj_in_frame'], False) # simulated sequence is not in frame
    #     self.assertEqual(simulated['cdr3_sequence_start'], 291) # simulated sequence cdr3 starts at 291
    #     self.assertNotEqual(simulated['cdr3_sequence_end'], 400) # simulated sequence cdr3 ends at 338
    #
if __name__ == '__main__':
    unittest.main()

