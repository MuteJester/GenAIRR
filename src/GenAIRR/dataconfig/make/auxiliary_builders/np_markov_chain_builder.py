from collections import Counter, defaultdict
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
            data: List of dicts (AIRR-style rows) with keys like
                  'v_sequence_end', 'd_sequence_start', 'd_sequence_end',
                  'j_sequence_start', 'sequence'.
        """
        self.set_lengths_from_data(data)
        self.set_first_bases_from_data(data)
        self.set_transitions_from_data(data)

    @staticmethod
    def _counter_to_prob_dict(counter):
        """Convert a Counter to a {value: probability} dict."""
        total = sum(counter.values())
        if total == 0:
            return {}
        return {k: v / total for k, v in counter.items()}

    def set_lengths_from_data(self, data):
        """
        Derive NP1 and NP2 length distributions from real AIRR data.
        """
        np1_counts = Counter(
            int(row['d_sequence_start']) - int(row['v_sequence_end'])
            for row in data
        )
        lengths = {'NP1': self._counter_to_prob_dict(np1_counts)}

        if self.has_d:
            np2_counts = Counter(
                int(row['j_sequence_start']) - int(row['d_sequence_end'])
                for row in data
            )
            lengths['NP2'] = self._counter_to_prob_dict(np2_counts)

        self.dataconfig.NP_lengths = lengths

    def set_first_bases_from_data(self, data):
        """
        Derive first base probabilities in NP1 and NP2 from real sequence data.
        """
        np_first_bases = {}

        # NP1 starts right after V
        np1_first = Counter(
            row['sequence'][int(row['v_sequence_end']) + 1]
            for row in data
            if int(row['v_sequence_end']) + 1 < len(row['sequence'])
        )
        np1_first.pop('N', None)
        np_first_bases['NP1'] = self._counter_to_prob_dict(np1_first)

        if self.has_d:
            # NP2 starts right after D
            np2_first = Counter(
                row['sequence'][int(row['d_sequence_end']) + 1]
                for row in data
                if int(row['d_sequence_end']) + 1 < len(row['sequence'])
            )
            np2_first.pop('N', None)
            np_first_bases['NP2'] = self._counter_to_prob_dict(np2_first)

        self.dataconfig.NP_first_bases = np_first_bases

    def set_transitions_from_data(self, data):
        """
        Derive transition probabilities for NP1 and NP2 from real sequence data.
        """
        bases = ['A', 'T', 'C', 'G']
        transitions = {}

        def extract_sequences(region):
            if region == 'NP1':
                return [
                    row['sequence'][int(row['v_sequence_end']):int(row['d_sequence_start'])]
                    for row in data
                ]
            elif region == 'NP2':
                return [
                    row['sequence'][int(row['d_sequence_end']):int(row['j_sequence_start'])]
                    for row in data
                ]

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
