import unittest

from GenAIRR.alleles import VAllele
from GenAIRR.generateDataConfig import RandomDataConfigGenerator, CustomDataConfigGenerator
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

        print([i['indels'] for i in generated_seqs])
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
                                                                              j_reference_path='./IGHJ.fasta')
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


if __name__ == '__main__':
    unittest.main()
