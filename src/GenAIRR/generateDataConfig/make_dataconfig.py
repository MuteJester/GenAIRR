import pandas as pd

from ..utilities import DataConfig
from ..utilities.data_utilities import create_allele_dict
from ..utilities.asc_utilities import create_asc_germline_set, hamming_distance
from ..utilities import AlleleNComparer
import pickle
from collections import defaultdict
from functools import partial
def normalize_and_filter_convert_to_dict(obj):
    if isinstance(obj, defaultdict):
        nested_dict = {k: normalize_and_filter_convert_to_dict(v) for k, v in obj.items()}
        if all(isinstance(val, float) for val in nested_dict.values()):
            # Filter out keys that are not 'A', 'T', 'G', or 'C'
            filtered_dict = {k: v for k, v in nested_dict.items() if k in ['A', 'T', 'G', 'C']}
            total = sum(filtered_dict.values())
            return {k: v / total for k, v in filtered_dict.items()}
        return nested_dict
    return obj

class RandomDataConfigGenerator:
    def __init__(self, convert_to_asc=True):
        self.convert_to_asc = convert_to_asc
        self.v_asc_table = None
        self.has_d = False
        # initialize an empty dataconfig object this class will populate
        self.dataconfig = DataConfig()
        # initialize auxiliary list for allele tracking
        self.alleles = ['V', 'J']

    def _get_reference_pointers(self):
        # aux list
        pointer_to_reference = {'V': self.dataconfig.v_alleles, 'J': self.dataconfig.j_alleles}
        if self.has_d:
            pointer_to_reference['D'] = self.dataconfig.d_alleles
        return pointer_to_reference

    def generate_decaying_probabilities(self, n, base=0.8):
        probabilities = {}
        total_sum = sum(base ** i for i in range(n + 1))
        for i in range(n + 1):
            probabilities[i] = (base ** i) / total_sum
        return probabilities

    def _load_alleles(self, v_alleles, j_alleles, d_alleles=None):
        self.dataconfig.v_alleles = v_alleles
        self.dataconfig.d_alleles = d_alleles
        self.dataconfig.j_alleles = j_alleles

    def _load_random_gene_usage(self):

        pointer_to_reference = self._get_reference_pointers()

        gene_use_dict = {'V': dict(), 'J': dict()}
        if self.has_d:
            gene_use_dict['D'] = dict()

        for allele in self.alleles:
            # get all alleles from reference per gene
            current_alleles = list(pointer_to_reference[allele])
            n = len(current_alleles)
            # uniformly attribute probability to each allele
            random_gene_usage = {i: (1 / n) for i in current_alleles}
            gene_use_dict[allele] = random_gene_usage

        self.dataconfig.gene_use_dict = gene_use_dict

    def _load_random_trimming_proportions(self, max_trim=50):
        pointer_to_reference = self._get_reference_pointers()
        trim_dicts = dict()
        for allele in self.alleles:
            # get only the families
            allele_families = sorted(list(set([i.split('-')[0] for i in pointer_to_reference[allele]])))
            # each family gets up to max_trim trimming amount options with decaying probabilities
            probabilities = self.generate_decaying_probabilities(max_trim)
            trim_dict = {i: probabilities for i in allele_families}

            trim_dicts[allele + '_5'] = trim_dict
            trim_dicts[allele + '_3'] = trim_dict

        self.dataconfig.trim_dicts = trim_dicts

    def _load_random_np_lengths(self, max_size=50):
        NP_lengths = dict()
        np_regions = ['NP1']
        if self.has_d:
            np_regions.append('NP2')

        for np_region in np_regions:
            probabilities = self.generate_decaying_probabilities(max_size)
            NP_lengths[np_region] = probabilities

        self.dataconfig.NP_lengths = NP_lengths

    def _load_random_np_first_base_use(self):
        NP_first_bases = dict()
        np_regions = ['NP1']
        if self.has_d:
            np_regions.append('NP2')

        for np_region in np_regions:
            base_probability = 1 / 4
            NP_first_bases[np_region] = {i: base_probability for i in ['A', 'T', 'C', 'G']}

        self.dataconfig.NP_first_bases = NP_first_bases

    def _load_random_np_transition_probabilities(self, max_size=50):
        NP_transitions = dict()
        nucleotides = ['A', 'T', 'C', 'G']
        np_regions = ['NP1']
        if self.has_d:
            np_regions.append('NP2')

        for np_region in np_regions:
            NP_transitions[np_region] = dict()
            for position in range(max_size):
                NP_transitions[np_region][position] = dict()
                for base_at_position in nucleotides:
                    # uniform transition probabilities
                    actions = {next_base: 1 / 4 for next_base in nucleotides}
                    NP_transitions[np_region][position][base_at_position] = actions

        self.dataconfig.NP_transitions = NP_transitions

    def _derive_3_prime_correction_map(self, target_alleles):
        t_dict = {i.name: i for j in target_alleles for i in target_alleles[j]}
        trim_map = dict()
        for t_allele in t_dict.values():
            trim_map[t_allele.name] = dict()
            seq_length = len(t_allele.ungapped_seq)
            for trim_3 in range(seq_length + 1):
                # Trim from the right (3' end)
                trimmed = t_allele.ungapped_seq[:seq_length - trim_3] if trim_3 > 0 else t_allele.ungapped_seq
                trim_map[t_allele.name][seq_length - trim_3] = []
                for v_c_allele in t_dict.values():
                    # Check if the trimmed sequence is a substring of the v_c_allele sequence
                    if trimmed in v_c_allele.ungapped_seq:
                        trim_map[t_allele.name][seq_length - trim_3].append(v_c_allele.name)
        r_allele = list(t_dict.values())[0]
        allele_ = str(r_allele.type).split('.')[1]  # returns "V" , "D" or "J"
        self.dataconfig.correction_maps[allele_ + '_3_TRIM_SIMILARITY_MAP'] = trim_map

    def _derive_5_prime_correction_map(self, target_alleles):
        t_dict = {i.name: i for j in target_alleles for i in target_alleles[j]}

        trim_map = dict()
        for t_allele in t_dict.values():
            trim_map[t_allele.name] = dict()
            for trim_5 in range(len(t_allele.ungapped_seq) + 1):
                trimmed = t_allele.ungapped_seq[trim_5:] if trim_5 > 0 else t_allele.ungapped_seq
                trim_map[t_allele.name][trim_5] = []
                for v_c_allele in t_dict.values():
                    # Check if the trimmed sequence is a substring of the v_c_allele sequence
                    if trimmed in v_c_allele.ungapped_seq:
                        trim_map[t_allele.name][trim_5].append(v_c_allele.name)
        r_allele = list(t_dict.values())[0]
        allele_ = str(r_allele.type).split('.')[1]  # returns "V" , "D" or "J"
        self.dataconfig.correction_maps[allele_ + '_5_TRIM_SIMILARITY_MAP'] = trim_map

    def _derive_5_and_3_prime_correction_map(self, target_alleles):
        t_dict = {i.name: i for j in target_alleles for i in target_alleles[j]}
        t_list = [i for j in target_alleles for i in target_alleles[j]]
        trim_map = dict()
        for t_allele in t_list:
            trim_map[t_allele.name] = dict()
            for trim_5 in range(len(t_allele.ungapped_seq) + 1):
                for trim_3 in range(len(t_allele.ungapped_seq) - trim_5 + 1):
                    # Correctly handle the trimming for t_allele
                    trimmed = t_allele.ungapped_seq[trim_5:] if trim_5 > 0 else t_allele.ungapped_seq
                    trimmed = trimmed[:-trim_3] if trim_3 > 0 else trimmed

                    trim_map[t_allele.name][(trim_5, trim_3)] = []
                    for d_c_allele in t_list:
                        # Check if the trimmed sequence is a substring of the d_c_allele sequence
                        if trimmed in d_c_allele.ungapped_seq:
                            trim_map[t_allele.name][(trim_5, trim_3)].append(d_c_allele.name)

        r_allele = list(t_dict.values())[0]
        allele_ = str(r_allele.type).split('.')[1]  # returns "V" , "D" or "J"
        self.dataconfig.correction_maps[allele_ + '_5_3_TRIM_SIMILARITY_MAP'] = trim_map

    def _derive_n_ambiguity_map(self, target_alleles):
        t_dict = {i.name: i for j in target_alleles for i in target_alleles[j]}
        comparer = AlleleNComparer()
        for v in t_dict:
            comparer.add_allele(v, t_dict[v].ungapped_seq.upper())

        r_allele = list(t_dict.values())[0]
        allele_ = str(r_allele.type).split('.')[1]  # returns "V" , "D" or "J"
        self.dataconfig.correction_maps[allele_ + '_N_AMBIGUITY_CORRECTION_GRAPH'] = comparer

    def make_dataconfig_from_existing_reference_files(self, v_reference_path, j_reference_path, d_reference_path=None):

        # update d flag
        self.has_d = d_reference_path is not None
        user_d_reference = None
        if self.has_d:
            # add D to aux list to calculate properties for D allele as well
            self.alleles.append('D')

        # 1. read fasta references
        if self.convert_to_asc:
            # ASC logic goes here to resulting variables should be of the following foramt:
            user_v_reference, v_asc_table = create_asc_germline_set(v_reference_path, segment="V")
            # save asc table so reverse transformation will be available to the user
            self.dataconfig.asc_tables['V'] = v_asc_table

            user_j_reference = create_allele_dict(j_reference_path)
            if self.has_d:
                user_d_reference = create_allele_dict(d_reference_path)
        else:

            user_v_reference = create_allele_dict(v_reference_path)
            if d_reference_path is not None:
                user_d_reference = create_allele_dict(d_reference_path)
            user_j_reference = create_allele_dict(j_reference_path)

        print('=' * 50)
        # 2. Fill in Data Config

        # LOAD ALLELES
        self._load_alleles(v_alleles=user_v_reference, d_alleles=user_d_reference, j_alleles=user_j_reference)
        print('Alleles Mounted to DataConfig!...')
        # RANDOM GENE USAGE
        self._load_random_gene_usage()
        print('Random Gene Usage Mounted to DataConfig!...')

        # TRIMMING PROPORTIONS
        self._load_random_trimming_proportions()
        print('Random Trimming Proportions Mounted to DataConfig!...')

        # N REGIONS LENGTHS
        self._load_random_np_lengths()
        print('Random NP Region Lengths Mounted to DataConfig!...')

        # N REGIONS  FIRST BASE USAGE
        self._load_random_np_first_base_use()
        print('Random NP Initial States Mounted to DataConfig!...')
        # N REGIONS MARKOV TRANSITION MATRICES
        self._load_random_np_transition_probabilities()
        print('Random NP Markov Chain Mounted to DataConfig!...')

        # 3. Fill in Data Config correction maps
        self._derive_n_ambiguity_map(self.dataconfig.v_alleles)
        print('V Ns Ambiguity Map Mounted to DataConfig!...')

        self._derive_3_prime_correction_map(self.dataconfig.v_alleles)
        print('V 3 Prime Ambiguity Map Mounted to DataConfig!...')
        self._derive_5_prime_correction_map(self.dataconfig.v_alleles)
        print('V 5 Prime Ambiguity Map Mounted to DataConfig!...')
        self._derive_3_prime_correction_map(self.dataconfig.j_alleles)
        print('J 3 Prime Ambiguity Map Mounted to DataConfig!...')
        self._derive_5_prime_correction_map(self.dataconfig.j_alleles)
        print('J 5 Prime Ambiguity Map Mounted to DataConfig!...')
        if self.has_d:
            self._derive_5_and_3_prime_correction_map(self.dataconfig.d_alleles)
            print('D (5,3) Prime Ambiguity Map Mounted to DataConfig!...')

        print('=' * 50)

        return self.dataconfig


class CustomDataConfigGenerator:
    def __init__(self, convert_to_asc=True):
        self.convert_to_asc = convert_to_asc
        self.v_asc_table = None
        self.has_d = False
        # initialize an empty dataconfig object this class will populate
        self.dataconfig = DataConfig()
        # initialize auxiliary list for allele tracking
        self.alleles = ['V', 'J']

    def _get_reference_pointers(self):
        # aux list
        pointer_to_reference = {'V': self.dataconfig.v_alleles, 'J': self.dataconfig.j_alleles}
        if self.has_d:
            pointer_to_reference['D'] = self.dataconfig.d_alleles
        return pointer_to_reference

    def generate_decaying_probabilities(self, n, base=0.8):
        probabilities = {}
        total_sum = sum(base ** i for i in range(n + 1))
        for i in range(n + 1):
            probabilities[i] = (base ** i) / total_sum
        return probabilities

    def _load_alleles(self, v_alleles, j_alleles, d_alleles=None):
        self.dataconfig.v_alleles = v_alleles
        self.dataconfig.d_alleles = d_alleles
        self.dataconfig.j_alleles = j_alleles

    def _derive_gene_usage(self, data):
        gene_use_dict = {'V': (data['v_call'].apply(lambda x: x.split('*')[0]).value_counts() / len(data)).to_dict(),
                         'J': (data['j_call'].apply(lambda x: x.split('*')[0]).value_counts() / len(data)).to_dict()
                         }

        if 'd_call' in data:
            gene_use_dict['D'] = (data['v_call'].apply(lambda x: x.split('*')[0]).value_counts() / len(data)).to_dict(),

        self.dataconfig.gene_use_dict = gene_use_dict

    def _derive_trimming_proportions(self, data):

        triming_dict = dict()

        for gene in ['v', 'd', 'j']:
            alleles = getattr(self.dataconfig, f'{gene}_alleles')
            families = [i.name.split('-')[0] for j in alleles for i in alleles[j]]
            families = list(set(families))
            for trim in ['5', '3']:
                triming_dict[f'{gene.upper()}_{trim}'] = dict()
                for fam in families:
                    samples = data[data[f'{gene}_call'].str.contains(fam)]
                    trim_values = (samples[f'{gene}_trim_{trim}'].value_counts() / len(samples)).to_dict()
                    triming_dict[f'{gene.upper()}_{trim}'][fam] = trim_values

        self.dataconfig.trim_dicts = triming_dict

    def _derive_np_lengths(self, data):
        np_lengths = {"NP1": {}, 'NP2': {}}
        np1_lengths = data['d_sequence_start'] - data['v_sequence_end']
        np1_lengths = (np1_lengths.value_counts() / len(np1_lengths)).to_dict()

        np2_lengths = data['j_sequence_start'] - data['d_sequence_end']
        np2_lengths = (np2_lengths.value_counts() / len(np2_lengths)).to_dict()

        np_lengths['NP1'] = np1_lengths
        np_lengths['NP2'] = np2_lengths

        self.dataconfig.NP_lengths = np_lengths

    def _derive_np_first_base_use(self, data):
        NP_first_bases = {"NP1": {}, 'NP2': {}}
        NP_first_bases['NP1'] = (
                    data.apply(lambda x: x['sequence'][x['v_sequence_end'] + 1], axis=1).value_counts() / len(
                data)).to_dict()
        NP_first_bases['NP2'] = (
                    data.apply(lambda x: x['sequence'][x['d_sequence_end'] + 1], axis=1).value_counts() / len(
                data)).to_dict()

        if 'N' in NP_first_bases['NP1']:
            NP_first_bases['NP1'].pop('N')
        if 'N' in NP_first_bases['NP2']:
            NP_first_bases['NP2'].pop('N')

        self.dataconfig.NP_first_bases = NP_first_bases

    def _derive_np_transition_probabilities(self, data):
        NP_transitions = {'NP1': {}, 'NP2': {}}

        NP1 = data.apply(lambda x: x['sequence'][x['v_sequence_end']:x['d_sequence_start']], axis=1)
        NP2 = data.apply(lambda x: x['sequence'][x['d_sequence_end']:x['j_sequence_start']], axis=1)

        agg_dict = defaultdict(partial(defaultdict, float))
        for np_seq in NP1:
            for pos in range(len(np_seq) - 1):
                if np_seq[pos] not in agg_dict[pos]:
                    agg_dict[pos][np_seq[pos]] = defaultdict(float)
                agg_dict[pos][np_seq[pos]][np_seq[pos + 1]] += 1
        NP_transitions['NP1'] = normalize_and_filter_convert_to_dict(agg_dict)

        agg_dict = defaultdict(partial(defaultdict, float))
        for np_seq in NP2:
            for pos in range(len(np_seq) - 1):
                if np_seq[pos] not in agg_dict[pos]:
                    agg_dict[pos][np_seq[pos]] = defaultdict(float)
                agg_dict[pos][np_seq[pos]][np_seq[pos + 1]] += 1
        NP_transitions['NP2'] = normalize_and_filter_convert_to_dict(agg_dict)

        self.dataconfig.NP_transitions = NP_transitions

    def _derive_3_prime_correction_map(self, target_alleles):
        t_dict = {i.name: i for j in target_alleles for i in target_alleles[j]}
        trim_map = dict()
        for t_allele in t_dict.values():
            trim_map[t_allele.name] = dict()
            seq_length = len(t_allele.ungapped_seq)
            for trim_3 in range(seq_length + 1):
                # Trim from the right (3' end)
                trimmed = t_allele.ungapped_seq[:seq_length - trim_3] if trim_3 > 0 else t_allele.ungapped_seq
                trim_map[t_allele.name][seq_length - trim_3] = []
                for v_c_allele in t_dict.values():
                    # Check if the trimmed sequence is a substring of the v_c_allele sequence
                    if trimmed in v_c_allele.ungapped_seq:
                        trim_map[t_allele.name][seq_length - trim_3].append(v_c_allele.name)
        r_allele = list(t_dict.values())[0]
        allele_ = str(r_allele.type).split('.')[1]  # returns "V" , "D" or "J"
        self.dataconfig.correction_maps[allele_ + '_3_TRIM_SIMILARITY_MAP'] = trim_map

    def _derive_5_prime_correction_map(self, target_alleles):
        t_dict = {i.name: i for j in target_alleles for i in target_alleles[j]}

        trim_map = dict()
        for t_allele in t_dict.values():
            trim_map[t_allele.name] = dict()
            for trim_5 in range(len(t_allele.ungapped_seq) + 1):
                trimmed = t_allele.ungapped_seq[trim_5:] if trim_5 > 0 else t_allele.ungapped_seq
                trim_map[t_allele.name][trim_5] = []
                for v_c_allele in t_dict.values():
                    # Check if the trimmed sequence is a substring of the v_c_allele sequence
                    if trimmed in v_c_allele.ungapped_seq:
                        trim_map[t_allele.name][trim_5].append(v_c_allele.name)
        r_allele = list(t_dict.values())[0]
        allele_ = str(r_allele.type).split('.')[1]  # returns "V" , "D" or "J"
        self.dataconfig.correction_maps[allele_ + '_5_TRIM_SIMILARITY_MAP'] = trim_map

    def _derive_5_and_3_prime_correction_map(self, target_alleles):
        t_dict = {i.name: i for j in target_alleles for i in target_alleles[j]}
        t_list = [i for j in target_alleles for i in target_alleles[j]]
        trim_map = dict()
        for t_allele in t_list:
            trim_map[t_allele.name] = dict()
            for trim_5 in range(len(t_allele.ungapped_seq) + 1):
                for trim_3 in range(len(t_allele.ungapped_seq) - trim_5 + 1):
                    # Correctly handle the trimming for t_allele
                    trimmed = t_allele.ungapped_seq[trim_5:] if trim_5 > 0 else t_allele.ungapped_seq
                    trimmed = trimmed[:-trim_3] if trim_3 > 0 else trimmed

                    trim_map[t_allele.name][(trim_5, trim_3)] = []
                    for d_c_allele in t_list:
                        # Check if the trimmed sequence is a substring of the d_c_allele sequence
                        if trimmed in d_c_allele.ungapped_seq:
                            trim_map[t_allele.name][(trim_5, trim_3)].append(d_c_allele.name)

        r_allele = list(t_dict.values())[0]
        allele_ = str(r_allele.type).split('.')[1]  # returns "V" , "D" or "J"
        self.dataconfig.correction_maps[allele_ + '_5_3_TRIM_SIMILARITY_MAP'] = trim_map

    def _derive_n_ambiguity_map(self, target_alleles):
        t_dict = {i.name: i for j in target_alleles for i in target_alleles[j]}
        comparer = AlleleNComparer()
        for v in t_dict:
            comparer.add_allele(v, t_dict[v].ungapped_seq.upper())

        r_allele = list(t_dict.values())[0]
        allele_ = str(r_allele.type).split('.')[1]  # returns "V" , "D" or "J"
        self.dataconfig.correction_maps[allele_ + '_N_AMBIGUITY_CORRECTION_GRAPH'] = comparer

    def make_dataconfig_from_existing_reference_files(self, v_reference_path, j_reference_path,
                                                      custom_data, d_reference_path=None,
                                                      ):
        # load data
        data = pd.read_csv(custom_data)
        # update d flag
        self.has_d = d_reference_path is not None
        user_d_reference = None
        if self.has_d:
            # add D to aux list to calculate properties for D allele as well
            self.alleles.append('D')

        # 1. read fasta references
        if self.convert_to_asc:
            # ASC logic goes here to resulting variables should be of the following foramt:
            user_v_reference, v_asc_table = create_asc_germline_set(v_reference_path, segment="V")
            # save asc table so reverse transformation will be available to the user
            self.dataconfig.asc_tables['V'] = v_asc_table

            user_j_reference = create_allele_dict(j_reference_path)
            if self.has_d:
                user_d_reference = create_allele_dict(d_reference_path)
        else:

            user_v_reference = create_allele_dict(v_reference_path)
            if d_reference_path is not None:
                user_d_reference = create_allele_dict(d_reference_path)
            user_j_reference = create_allele_dict(j_reference_path)

        print('=' * 50)
        # 2. Fill in Data Config

        # LOAD ALLELES
        self._load_alleles(v_alleles=user_v_reference, d_alleles=user_d_reference, j_alleles=user_j_reference)
        print('Alleles Mounted to DataConfig!...')
        # RANDOM GENE USAGE
        self._derive_gene_usage(data)
        print('Gene Usage Mounted to DataConfig!...')

        # TRIMMING PROPORTIONS
        self._derive_trimming_proportions(data)
        print('Trimming Proportions Mounted to DataConfig!...')

        # N REGIONS LENGTHS
        self._derive_np_lengths(data)
        print('NP Region Lengths Mounted to DataConfig!...')

        # N REGIONS  FIRST BASE USAGE
        self._derive_np_first_base_use(data)
        print('NP Initial States Mounted to DataConfig!...')
        # N REGIONS MARKOV TRANSITION MATRICES
        self._derive_np_transition_probabilities(data)
        print('NP Markov Chain Mounted to DataConfig!...')

        # 3. Fill in Data Config correction maps
        self._derive_n_ambiguity_map(self.dataconfig.v_alleles)
        print('V Ns Ambiguity Map Mounted to DataConfig!...')

        self._derive_3_prime_correction_map(self.dataconfig.v_alleles)
        print('V 3 Prime Ambiguity Map Mounted to DataConfig!...')
        self._derive_5_prime_correction_map(self.dataconfig.v_alleles)
        print('V 5 Prime Ambiguity Map Mounted to DataConfig!...')
        self._derive_3_prime_correction_map(self.dataconfig.j_alleles)
        print('J 3 Prime Ambiguity Map Mounted to DataConfig!...')
        self._derive_5_prime_correction_map(self.dataconfig.j_alleles)
        print('J 5 Prime Ambiguity Map Mounted to DataConfig!...')
        if self.has_d:
            self._derive_5_and_3_prime_correction_map(self.dataconfig.d_alleles)
            print('D (5,3) Prime Ambiguity Map Mounted to DataConfig!...')

        print('=' * 50)

        return self.dataconfig


class AdaptiveDataConfigGenerator:
    def __init__(self, path_to_base_data_config, convert_to_asc=True):
        self.convert_to_asc = convert_to_asc

        # here we load the dataconfig from which we will estimate the closest matches
        with open(path_to_base_data_config, 'rb') as h:
            self.base_dataconfig = pickle.load(h)

        self.v_asc_table = None
        self.has_d = False
        # initialize an empty dataconfig object this class will populate
        self.dataconfig = DataConfig()
        # initialize auxiliary list for allele tracking
        self.alleles = ['V', 'J']

    def _get_reference_pointers(self):
        # aux list
        pointer_to_reference = {'V': self.dataconfig.v_alleles, 'J': self.dataconfig.j_alleles}
        if self.has_d:
            pointer_to_reference['D'] = self.dataconfig.d_alleles
        return pointer_to_reference

    def _get_base_reference_pointers(self):
        # aux list
        pointer_to_reference = {'V': self.base_dataconfig.v_alleles, 'J': self.base_dataconfig.j_alleles}
        if self.has_d:
            pointer_to_reference['D'] = self.base_dataconfig.d_alleles
        return pointer_to_reference

    def generate_decaying_probabilities(self, n, base=0.8):
        probabilities = {}
        total_sum = sum(base ** i for i in range(n + 1))
        for i in range(n + 1):
            probabilities[i] = (base ** i) / total_sum
        return probabilities

    def _load_alleles(self, v_alleles, j_alleles, d_alleles=None):
        self.dataconfig.v_alleles = v_alleles
        self.dataconfig.d_alleles = d_alleles
        self.dataconfig.j_alleles = j_alleles

    def _closest_match(self):
        
        pointer_to_reference = self._get_reference_pointers()
        pointer_to_base_reference = self._get_base_reference_pointers()

        gene_match_dict = {'V': dict(), 'J': dict()}
        if self.has_d:
            gene_match_dict['D'] = dict()

        for allele in self.alleles:
            # get all alleles from reference per gene
            current_alleles = {j.name:j for i in pointer_to_reference[allele] for j in pointer_to_reference[allele][i]}
            base_alleles = {j.name:j for i in pointer_to_base_reference[allele] for j in pointer_to_base_reference[allele][i]}
        
            match_alleles = {i: None for i in list(pointer_to_reference[allele])}    
            
            # check if the alleles are the same
            for name, current_allele in current_alleles.items():
                seq = current_allele.gapped_seq
                for base_name, base_allele in base_alleles.items():
                    if base_allele.gapped_seq == seq:
                        if match_alleles[current_allele.gene] is None:
                            match_alleles[current_allele.gene] = [base_allele.gene]    
                        else:
                            match_alleles[current_allele.gene].append(base_allele.gene)
                        break
            
            # get the alleles that are not matched
            not_matched = [i for i in match_alleles if match_alleles[i] is None]
            
            if not_matched:
                for current_allele in not_matched:
                    current_allele = current_alleles[name]
                    min_distance = float('inf')
                    for base_name, base_allele in base_alleles.items():
                        s1 = current_allele.gapped_seq
                        s2 = base_allele.gapped_seq
                        max_length = max(len(s1), len(s2))
                        s1 = s1.ljust(max_length, 'N') 
                        s2 = s2.ljust(max_length, 'N') 
                        distance = hamming_distance(s1, s2) 
                        if distance < min_distance: 
                            min_distance = distance 
                            if match_alleles[current_allele.gene] is None:
                                match_alleles[current_allele.gene] = [base_allele.gene]    
                            else:
                                match_alleles[current_allele.gene].append(base_allele.gene)
            gene_match_dict[allele] = match_alleles
        self.dataconfig.closest_matches = gene_match_dict
        
    def _match_closest_gene_usage(self):
        # TO DO: This is the Random Implementation change this to populate
        # self.dataconfig.gene_use_dict with an updated version with the new alleles based on closest ASC
        
        # it assumes that self.base_dataconfig.gene_use_dict genes are the same as in self.base_dataconfig.v_alleles
        pointer_to_reference = self._get_reference_pointers()
        
        gene_use_dict = {'V': dict(), 'J': dict()}
        if self.has_d:
            gene_use_dict['D'] = dict()

        for allele in self.alleles:
            # get all alleles from reference per gene
            current_alleles = list(pointer_to_reference[allele])
            closest_match = self.dataconfig.closest_matches[allele]
            n = len(current_alleles)
            # fills up in case that self.dataconfig.gene_use_dict and self.dataconfig.v_alleles are not the same
            closest_gene_usage = {i: (1 / n) for i in current_alleles}
            for i in current_alleles:
                unique_gene = list(set(closest_match[i]))
                if all(gene in self.base_dataconfig.gene_use_dict[allele].keys() for gene in unique_gene):
                    if len(unique_gene) == 1:
                        closest_gene_usage[i] = self.base_dataconfig.gene_use_dict[allele][unique_gene[0]]
                    else:
                        closest_gene_usage[i] = sum(self.base_dataconfig.gene_use_dict[allele][gene] for gene in unique_gene) / len(unique_gene)

            gene_use_dict[allele] = closest_gene_usage

        self.dataconfig.gene_use_dict = gene_use_dict

    def _match_closest_trimming_proportions(self, max_trim=50):

        # TO DO: This is the Random Implementation change this to populate
        # self.dataconfig.trim_dicts with an updated version with the new alleles based on closest ASC

        # it assumes that self.base_dataconfig.trim_dict genes are the same as in self.base_dataconfig.v_alleles
        # getting the probabilites for 'V_3', 'D_5', 'D_3', 'J_5'
        
        pointer_to_reference = self._get_reference_pointers()
        trim_dicts = dict()
        for allele in self.alleles:
            # get only the families
            allele_families = sorted(list(set([i.split('-')[0] for i in pointer_to_reference[allele]])))
            current_families = dict()
            for i in list(pointer_to_reference[allele]):
                f = i.split('-')[0]
                closest_f = list(set([j.split('-')[0] for j in self.dataconfig.closest_matches[i]]))
                if f not in current_families:
                    current_families[f] = closest_f
                else:
                    current_families[f] = current_families[f] + closest_f
            current_families = {i: list(set(current_families[i])) for i in current_families}
            

            # each family gets up to max_trim trimming amount options with decaying probabilities
            probabilities = self.generate_decaying_probabilities(max_trim)
            trim_dict = {i: probabilities for i in allele_families}

            if allele == 'V':
                # get the values for V_3
                trim_dict_closest = dict()
                for k,v in current_families.items():
                    trim_dict_closest[k] = self.base_dataconfig.trim_dicts[allele + '_3'][v[0]]
                trim_dicts[allele + '_3'] = trim_dict_closest  
                trim_dicts[allele + '_5'] = trim_dict
            
            if allele == 'J':
                trim_dict_closest = dict()
                for k,v in current_families.items():
                    trim_dict_closest[k] = self.base_dataconfig.trim_dicts[allele + '_5'][v[0]]
                trim_dicts[allele + '_5'] = trim_dict_closest  
                trim_dicts[allele + '_3'] = trim_dict
            
            if allele == 'D':
                trim_dict_closest = dict()
                for k,v in current_families.items():
                    trim_dict_closest[k] = self.base_dataconfig.trim_dicts[allele + '_5'][v[0]]
                trim_dicts[allele + '_5'] = trim_dict_closest
                
                trim_dict_closest = dict()
                for k,v in current_families.items():
                    trim_dict_closest[k] = self.base_dataconfig.trim_dicts[allele + '_3'][v[0]]
                trim_dicts[allele + '_3'] = trim_dict
            
            

        self.dataconfig.trim_dicts = trim_dicts

    def _match_closest_np_lengths(self):

        # TO DO: This is the Random Implementation change this to populate
        # self.dataconfig.NP_lengths with an updated version with the new alleles based on closest ASC

        # the data is not depended on the closest match. thus taking the base dataconfig
        
        NP_lengths = dict()
        np_regions = ['NP1']
        if self.has_d:
            np_regions.append('NP2')

        for np_region in np_regions:
            NP_lengths[np_region] = self.base_dataconfig.NP_lengths[np_region]
            
        self.dataconfig.NP_lengths = NP_lengths

    def _match_closest_np_first_base_use(self):

        # TO DO: This is the Random Implementation change this to populate
        # self.dataconfig.NP_first_bases with an updated version with the new alleles based on closest ASC
        # the data is not depended on the closest match. thus taking the base dataconfig

        NP_first_bases = dict()
        np_regions = ['NP1']
        if self.has_d:
            np_regions.append('NP2')

        for np_region in np_regions:
            NP_first_bases[np_region] = self.base_dataconfig.NP_first_bases[np_region]

        self.dataconfig.NP_first_bases = NP_first_bases

    def _match_closest_np_transition_probabilities(self):

        # TO DO: This is the Random Implementation change this to populate
        # self.dataconfig.NP_transitions with an updated version with the new alleles based on closest ASC
        # the data is not depended on the closest match. thus taking the base dataconfig

        NP_transitions = dict()
        nucleotides = ['A', 'T', 'C', 'G']
        np_regions = ['NP1']
        if self.has_d:
            np_regions.append('NP2')

        for np_region in np_regions:
            NP_transitions[np_region] = dict()
            for position in range(max_size):
                NP_transitions[np_region][position] = dict()
                for base_at_position in nucleotides:
                    # uniform transition probabilities
                    NP_transitions[np_region][position][base_at_position] = self.base_dataconfig.NP_transitions[np_region][position][base_at_position]

        self.dataconfig.NP_transitions = NP_transitions

    def _derive_3_prime_correction_map(self, target_alleles):
        t_dict = {i.name: i for j in target_alleles for i in target_alleles[j]}
        trim_map = dict()
        for t_allele in t_dict.values():
            trim_map[t_allele.name] = dict()
            seq_length = len(t_allele.ungapped_seq)
            for trim_3 in range(seq_length + 1):
                # Trim from the right (3' end)
                trimmed = t_allele.ungapped_seq[:seq_length - trim_3] if trim_3 > 0 else t_allele.ungapped_seq
                trim_map[t_allele.name][seq_length - trim_3] = []
                for v_c_allele in t_dict.values():
                    # Check if the trimmed sequence is a substring of the v_c_allele sequence
                    if trimmed in v_c_allele.ungapped_seq:
                        trim_map[t_allele.name][seq_length - trim_3].append(v_c_allele.name)
        r_allele = list(t_dict.values())[0]
        allele_ = str(r_allele.type).split('.')[1]  # returns "V" , "D" or "J"
        self.dataconfig.correction_maps[allele_ + '_3_TRIM_SIMILARITY_MAP'] = trim_map

    def _derive_5_prime_correction_map(self, target_alleles):
        t_dict = {i.name: i for j in target_alleles for i in target_alleles[j]}

        trim_map = dict()
        for t_allele in t_dict.values():
            trim_map[t_allele.name] = dict()
            for trim_5 in range(len(t_allele.ungapped_seq) + 1):
                trimmed = t_allele.ungapped_seq[trim_5:] if trim_5 > 0 else t_allele.ungapped_seq
                trim_map[t_allele.name][trim_5] = []
                for v_c_allele in t_dict.values():
                    # Check if the trimmed sequence is a substring of the v_c_allele sequence
                    if trimmed in v_c_allele.ungapped_seq:
                        trim_map[t_allele.name][trim_5].append(v_c_allele.name)
        r_allele = list(t_dict.values())[0]
        allele_ = str(r_allele.type).split('.')[1]  # returns "V" , "D" or "J"
        self.dataconfig.correction_maps[allele_ + '_5_TRIM_SIMILARITY_MAP'] = trim_map

    def _derive_5_and_3_prime_correction_map(self, target_alleles):
        t_dict = {i.name: i for j in target_alleles for i in target_alleles[j]}
        t_list = [i for j in target_alleles for i in target_alleles[j]]
        trim_map = dict()
        for t_allele in t_list:
            trim_map[t_allele.name] = dict()
            for trim_5 in range(len(t_allele.ungapped_seq) + 1):
                for trim_3 in range(len(t_allele.ungapped_seq) - trim_5 + 1):
                    # Correctly handle the trimming for t_allele
                    trimmed = t_allele.ungapped_seq[trim_5:] if trim_5 > 0 else t_allele.ungapped_seq
                    trimmed = trimmed[:-trim_3] if trim_3 > 0 else trimmed

                    trim_map[t_allele.name][(trim_5, trim_3)] = []
                    for d_c_allele in t_list:
                        # Check if the trimmed sequence is a substring of the d_c_allele sequence
                        if trimmed in d_c_allele.ungapped_seq:
                            trim_map[t_allele.name][(trim_5, trim_3)].append(d_c_allele.name)

        r_allele = list(t_dict.values())[0]
        allele_ = str(r_allele.type).split('.')[1]  # returns "V" , "D" or "J"
        self.dataconfig.correction_maps[allele_ + '_5_3_TRIM_SIMILARITY_MAP'] = trim_map

    def _derive_n_ambiguity_map(self, target_alleles):
        t_dict = {i.name: i for j in target_alleles for i in target_alleles[j]}
        comparer = AlleleNComparer()
        for v in t_dict:
            comparer.add_allele(v, t_dict[v].ungapped_seq.upper())

        r_allele = list(t_dict.values())[0]
        allele_ = str(r_allele.type).split('.')[1]  # returns "V" , "D" or "J"
        self.dataconfig.correction_maps[allele_ + '_N_AMBIGUITY_CORRECTION_GRAPH'] = comparer

    def make_dataconfig_from_existing_reference_files(self, v_reference_path, j_reference_path, d_reference_path=None):

        # update d flag
        self.has_d = d_reference_path is not None
        user_d_reference = None
        if self.has_d:
            # add D to aux list to calculate properties for D allele as well
            self.alleles.append('D')

        # 1. read fasta references
        if self.convert_to_asc:
            # ASC logic goes here to resulting variables should be of the following foramt:
            user_v_reference, v_asc_table = create_asc_germline_set(v_reference_path, segment="V")
            # save asc table so reverse transformation will be available to the user
            self.dataconfig.asc_tables['V'] = v_asc_table

            user_j_reference = create_allele_dict(j_reference_path)
            if self.has_d:
                user_d_reference = create_allele_dict(d_reference_path)
        else:
            user_v_reference = create_allele_dict(v_reference_path)
            if d_reference_path is not None:
                user_d_reference = create_allele_dict(d_reference_path)
            user_j_reference = create_allele_dict(j_reference_path)

        print('=' * 50)
        # 2. Fill in Data Config

        # LOAD ALLELES
        self._load_alleles(v_alleles=user_v_reference, d_alleles=user_d_reference, j_alleles=user_j_reference)
        print('Alleles Mounted to DataConfig!...')
        # RANDOM GENE USAGE
        self._match_closest_gene_usage()
        print('Random Gene Usage Mounted to DataConfig!...')

        # TRIMMING PROPORTIONS
        self._match_closest_trimming_proportions()
        print('Random Trimming Proportions Mounted to DataConfig!...')

        # N REGIONS LENGTHS
        self._match_closest_np_lengths()
        print('Random NP Region Lengths Mounted to DataConfig!...')

        # N REGIONS  FIRST BASE USAGE
        self._match_closest_np_first_base_use()
        print('Random NP Initial States Mounted to DataConfig!...')
        # N REGIONS MARKOV TRANSITION MATRICES
        self._match_closest_np_transition_probabilities()
        print('Random NP Markov Chain Mounted to DataConfig!...')

        # ======================================================================= #
        # 3. Fill in Data Config correction maps
        self._derive_n_ambiguity_map(self.dataconfig.v_alleles)
        print('V Ns Ambiguity Map Mounted to DataConfig!...')

        self._derive_3_prime_correction_map(self.dataconfig.v_alleles)
        print('V 3 Prime Ambiguity Map Mounted to DataConfig!...')
        self._derive_5_prime_correction_map(self.dataconfig.v_alleles)
        print('V 5 Prime Ambiguity Map Mounted to DataConfig!...')
        self._derive_3_prime_correction_map(self.dataconfig.j_alleles)
        print('J 3 Prime Ambiguity Map Mounted to DataConfig!...')
        self._derive_5_prime_correction_map(self.dataconfig.j_alleles)
        print('J 5 Prime Ambiguity Map Mounted to DataConfig!...')
        if self.has_d:
            self._derive_5_and_3_prime_correction_map(self.dataconfig.d_alleles)
            print('D (5,3) Prime Ambiguity Map Mounted to DataConfig!...')

        print('=' * 50)

        return self.dataconfig
