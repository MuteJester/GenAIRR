import random
from enum import Enum, auto
from ..mutation import MutationModel
from ..sequence import NP_Region
from ..sequence.sequence import BaseSequence
from ..utilities import translate
from GenAIRR.dataconfig.data_config import DataConfig


class LightChainType(Enum):
    KAPPA = auto()
    LAMBDA = auto()


class LightChainSequence(BaseSequence):
    """
    Represents a light chain sequence in an immunoglobulin, including V and J segments and an NP region.

    This class extends `BaseSequence` to include functionality specific to light chain sequences, such as
    sequence simulation, mutation, and functionality checks based on junction properties.

    Args:
        alleles (list): A list of `Allele` instances for the V and J gene segments.
        dataconfig (DataConfig): A `DataConfig` instance containing configuration data for sequence simulation.

    Attributes:
        Inherits all attributes from `BaseSequence` and adds/modifies specific ones for light chain sequences.
    """

    def __init__(self, alleles, dataconfig: DataConfig):
        """Initializes a `LightChainSequence` with the given alleles and simulates the sequence."""
        super().__init__(alleles)
        self.simulate_sequence(dataconfig)

    def simulate_sequence(self, dataconfig: DataConfig):
        self.simulate_NP_regions(dataconfig)
        self.simulate_trimmed_sequences(dataconfig)
        self.assemble_sequence()
        self.calculate_junction_properties()
        self.update_metadata()
        self.check_functionality()

    def simulate_NP_regions(self, dataconfig: DataConfig):
        self.NP1_region = NP_Region.create_np_region(dataconfig.NP_lengths, dataconfig.NP_transitions, "NP1", dataconfig.NP_first_bases)
        self.NP1_length = len(self.NP1_region)

    def simulate_trimmed_sequences(self, dataconfig: DataConfig):
        self.v_trimmed_seq, self.v_trim_5, self.v_trim_3 = self.v_allele.get_trimmed(dataconfig.trim_dicts)
        self.j_trimmed_seq, self.j_trim_5, self.j_trim_3 = self.j_allele.get_trimmed(dataconfig.trim_dicts)
        if self.c_allele is not None:
            self.c_trimmed_seq, self.c_trim_5, self.c_trim_3 = self.c_allele.get_trimmed(dataconfig.trim_dicts)


    def assemble_sequence(self):
        self.ungapped_seq = (
            self.v_trimmed_seq
            + self.NP1_region
            + self.j_trimmed_seq

        )
        if self.c_allele is not None:
            self.ungapped_seq += self.c_trimmed_seq

        self.ungapped_seq = self.ungapped_seq.upper()

    def calculate_junction_properties(self):
        self.junction_length = self.get_junction_length()
        self.junction = self.ungapped_seq[self.v_allele.anchor:self.v_allele.anchor + self.junction_length].upper()

    def update_metadata(self):
        self.v_seq_start = 0
        self.v_seq_end = self.v_allele.ungapped_len - self.v_trim_3
        self.j_seq_start = self.v_seq_end + self.NP1_length
        self.j_seq_end = self.j_seq_start + self.j_allele.ungapped_len - self.j_trim_5
        self.v_germline_start = 0
        self.v_germline_end = self.v_allele.ungapped_len - self.v_trim_3
        self.j_germline_start = self.j_trim_5
        self.j_germline_end = self.j_allele.ungapped_len - self.j_trim_3
        self.junction_start = self.v_allele.anchor
        self.junction_end = self.v_allele.anchor + self.junction_length

    def check_functionality(self,sequence=None):
        self.functional = False
        self.stop_codon = self.check_stops(self.ungapped_seq if sequence is None else sequence)
        self.vj_in_frame = (
            (self.junction_end % 3) == 0
            and (self.junction_start % 3 == 0)
            and (self.junction_length % 3 == 0)
            and not self.stop_codon
        )
        self.note = ''
        if (self.junction_length % 3) == 0 and not self.stop_codon:
            self.junction_aa = translate(self.junction)
            if self.junction_aa.startswith("C"):
                if self.junction_aa.endswith("F") or self.junction_aa.endswith("W"):
                    self.functional = True
                else:
                    self.note += 'J anchor (W/F) not present.'
            else:
                self.note += 'V second C not present.'

    def mutate(self, mutation_model: MutationModel):
        mutated_sequence, mutations, mutation_rate = mutation_model.apply_mutation(self)
        self.mutated_seq = mutated_sequence
        self.mutations = mutations
        self.mutation_freq = mutation_rate
        self.mutation_count = len(mutations)
        self.junction = self.mutated_seq[self.v_allele.anchor:self.v_allele.anchor + self.junction_length].upper()
        self.check_functionality(self.mutated_seq)

    def get_junction_length(self):
        return (self.v_allele.length - (self.v_allele.anchor - 1) - self.v_trim_3 +
                self.NP1_length + (self.j_allele.anchor + 2) - self.j_trim_5)

    def check_stops(self, seq):
        stops = ["TAG", "TAA", "TGA"]
        for x in range(0, len(seq), 3):
            if seq[x:x + 3] in stops:
                return True
        return False

    @classmethod
    def create_random(cls, dataconfig: DataConfig, specific_v=None, specific_j=None):
        if specific_v is None:
            random_v_allele = random.choice([i for j in dataconfig.v_alleles for i in dataconfig.v_alleles[j]])
        else:
            random_v_allele = specific_v

        if specific_j is None:
            random_j_allele = random.choice([i for j in dataconfig.j_alleles for i in dataconfig.j_alleles[j]])
        else:
            random_j_allele = specific_j

        random_c_allele = random.choice([i for j in dataconfig.c_alleles for i in dataconfig.c_alleles[j]])

        return cls([random_v_allele, random_j_allele], dataconfig) # removed C allele will be fixed later

    def __repr__(self):
        total_length = len(self.ungapped_seq)
        proportional_length = lambda start, end: int((end - start) / total_length * 100)

        v_length = proportional_length(self.v_seq_start, self.v_seq_end)
        j_length = proportional_length(self.j_seq_start, self.j_seq_end)
        c_length = proportional_length(self.j_seq_end, total_length)

        v_part = f"{self.v_seq_start}|{'-' * v_length}V({self.v_allele.name})|{self.v_seq_end}"
        j_part = f"{self.j_seq_start}|{'-' * j_length}J({self.j_allele.name})|{self.j_seq_end}"
        c_part = f"{self.j_seq_end}|{'-' * c_length}C({self.c_allele.name})|{total_length}"

        return f"{v_part}|{j_part}|{c_part}"
