from collections import defaultdict
from functools import partial
from ....utilities.misc import normalize_and_filter_convert_to_dict


class NPMarkovParameterBuilder:
    """
    Generates or derives NP-region related parameters:
    - Length distributions
    - First base probabilities
    - Markov transition matrices
    """

    def __init__(self, dataconfig, has_d=False):
        """
        Args:
            dataconfig (DataConfig): Object to populate.
            has_d (bool): Whether D alleles are used (affects whether NP2 is used).
        """
        self.dataconfig = dataconfig
        self.has_d = has_d
        self.np_regions = ['NP1', 'NP2'] if has_d else ['NP1']

    @staticmethod
    def _generate_decaying_probs(n, base=0.8):
        total = sum(base ** i for i in range(n + 1))
        return {i: (base ** i) / total for i in range(n + 1)}

    def generate_all_random(self, max_length=50):
        """
        Generate all NP-region parameters using random assumptions.
        """
        self.set_random_lengths(max_length)
        self.set_random_first_bases()
        self.set_random_transitions(max_length)

    def set_random_lengths(self, max_size=50):
        """
        Set NP region lengths with decaying probabilities.
        """
        self.dataconfig.NP_lengths = {
            region: self._generate_decaying_probs(max_size)
            for region in self.np_regions
        }

    def set_random_first_bases(self):
        """
        Assign uniform base distributions to first base in each NP region.
        """
        bases = ['A', 'T', 'C', 'G']
        base_probs = {b: 0.25 for b in bases}
        self.dataconfig.NP_first_bases = {
            region: base_probs for region in self.np_regions
        }

    def set_random_transitions(self, max_size=50):
        """
        Assign uniform transition probabilities for NP regions.
        """
        bases = ['A', 'T', 'C', 'G']
        self.dataconfig.NP_transitions = {
            region: {
                pos: {
                    b: {next_b: 0.25 for next_b in bases} for b in bases
                } for pos in range(max_size)
            } for region in self.np_regions
        }

    def derive_all_from_data(self, data):
        """
        Derive all NP-region parameters directly from observed sequence data.

        Args:
            data (pd.DataFrame): AIRR-style annotated dataframe with V, D, J calls and positions.
        """
        self.set_lengths_from_data(data)
        self.set_first_bases_from_data(data)
        self.set_transitions_from_data(data)

    def set_lengths_from_data(self, data):
        """
        Derive NP1 and NP2 length distributions from real AIRR data.
        """
        lengths = {}

        np1_lengths = data['d_sequence_start'] - data['v_sequence_end']
        np2_lengths = data['j_sequence_start'] - data['d_sequence_end']

        lengths['NP1'] = (np1_lengths.value_counts() / len(np1_lengths)).to_dict()
        if self.has_d:
            lengths['NP2'] = (np2_lengths.value_counts() / len(np2_lengths)).to_dict()

        self.dataconfig.NP_lengths = lengths

    def set_first_bases_from_data(self, data):
        """
        Derive first base probabilities in NP1 and NP2 from real sequence data.
        """
        np_first_bases = {}

        # NP1 starts right after V
        np1_first = data.apply(lambda x: x['sequence'][x['v_sequence_end'] + 1], axis=1)
        np_first_bases['NP1'] = (np1_first.value_counts() / len(np1_first)).to_dict()

        if self.has_d:
            # NP2 starts right after D
            np2_first = data.apply(lambda x: x['sequence'][x['d_sequence_end'] + 1], axis=1)
            np_first_bases['NP2'] = (np2_first.value_counts() / len(np2_first)).to_dict()

        # Remove ambiguous 'N'
        for region in np_first_bases:
            np_first_bases[region].pop('N', None)

        self.dataconfig.NP_first_bases = np_first_bases

    def set_transitions_from_data(self, data):
        """
        Derive transition probabilities for NP1 and NP2 from real sequence data.
        """
        bases = ['A', 'T', 'C', 'G']
        transitions = {}

        def extract_sequences(region):
            if region == 'NP1':
                return data.apply(lambda x: x['sequence'][x['v_sequence_end']:x['d_sequence_start']], axis=1)
            elif region == 'NP2':
                return data.apply(lambda x: x['sequence'][x['d_sequence_end']:x['j_sequence_start']], axis=1)

        for region in self.np_regions:
            agg_dict = defaultdict(partial(defaultdict, float))
            sequences = extract_sequences(region)

            for seq in sequences:
                for pos in range(len(seq) - 1):
                    curr, next_ = seq[pos], seq[pos + 1]
                    if curr not in bases or next_ not in bases:
                        continue
                    if curr not in agg_dict[pos]:
                        agg_dict[pos][curr] = defaultdict(float)
                    agg_dict[pos][curr][next_] += 1

            transitions[region] = normalize_and_filter_convert_to_dict(agg_dict)

        self.dataconfig.NP_transitions = transitions
