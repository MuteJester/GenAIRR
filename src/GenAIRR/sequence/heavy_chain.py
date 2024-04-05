import random

from ..mutation import MutationModel
from ..sequence import NP_Region
from ..sequence.sequence import BaseSequence
from ..utilities import translate
from ..utilities.data_config import DataConfig


class HeavyChainSequence(BaseSequence):
    """
        Represents a heavy chain sequence in an immunoglobulin, including V, D, and J segments,
        along with NP1 and NP2 regions.

        This class extends `BaseSequence` to include functionality specific to heavy chain sequences, such as
        sequence simulation with D segment and two NP regions, mutation, and functionality checks based on junction properties.

        Args:
            alleles (list): A list of `Allele` instances for the V, D, and J gene segments.
            dataconfig (DataConfig): A `DataConfig` instance containing configuration data for sequence simulation.
        """
    def __init__(self, alleles, dataconfig: DataConfig):
        """
        Initializes a `HeavyChainSequence` with the given alleles and simulates the sequence using the provided data configuration.
        """
        super().__init__(alleles)
        self.simulate_sequence(dataconfig)

    def simulate_sequence(self, dataconfig: DataConfig):
        """
        Simulates the heavy chain sequence by trimming V, D, J alleles, adding NP1 and NP2 regions,
        and assembling the final sequence.

        Args:
            dataconfig (DataConfig): Configuration data for sequence simulation, including trimming dictionaries,
            NP lengths, and transition probabilities.
        """
        v_allele = self.v_allele
        d_allele = self.d_allele
        j_allele = self.j_allele


        self.NP1_region = NP_Region.create_np_region(dataconfig.NP_lengths, dataconfig.NP_transitions, "NP1",
                                                     dataconfig.NP_first_bases)
        self.NP2_region = NP_Region.create_np_region(dataconfig.NP_lengths, dataconfig.NP_transitions, "NP2",
                                                     dataconfig.NP_first_bases)

        self.NP1_length = len(self.NP1_region)
        self.NP2_length = len(self.NP2_region)

        v_trimmed_seq, v_trim_5, v_trim_3 = v_allele.get_trimmed(dataconfig.trim_dicts)
        d_trimmed_seq, d_trim_5, d_trim_3 = d_allele.get_trimmed(dataconfig.trim_dicts)
        j_trimmed_seq, j_trim_5, j_trim_3 = j_allele.get_trimmed(dataconfig.trim_dicts)

        nuc_seq = (
                v_trimmed_seq
                + self.NP1_region
                + d_trimmed_seq
                + self.NP2_region
                + j_trimmed_seq
        )

        # log trims
        self.v_trim_5 = v_trim_5
        self.v_trim_3 = v_trim_3
        self.d_trim_5 = d_trim_5
        self.d_trim_3 = d_trim_3
        self.j_trim_5 = j_trim_5
        self.j_trim_3 = j_trim_3

        self.ungapped_seq = nuc_seq.upper()

        
        self.junction_length = self.get_junction_length()
        
        self.junction = self.ungapped_seq[self.v_allele.anchor:
                                          self.v_allele.anchor + self.junction_length].upper()
        
        self.update_metadata()
        self._is_functional(self.ungapped_seq)

    def update_metadata(self):
        """
        Updates the metadata for the heavy chain sequence, including the start and end positions
        of the V, D, and J segments, and the Junction.
        """
        self.v_seq_start = 0
        self.v_seq_end = self.v_allele.ungapped_len - self.v_trim_3
        self.d_seq_start = self.v_seq_end + self.NP1_length
        self.d_seq_end = self.d_seq_start + self.d_allele.ungapped_len - self.d_trim_3 - self.d_trim_5
        self.j_seq_start = self.d_seq_end + self.NP2_length
        self.j_seq_end = self.j_seq_start + self.j_allele.ungapped_len - self.j_trim_5
        self.v_germline_start = 0
        self.v_germline_end = self.v_allele.ungapped_len - self.v_trim_3
        self.d_germline_start = self.d_trim_5
        self.d_germline_end = self.d_allele.ungapped_len - self.d_trim_3
        self.j_germline_start = self.j_trim_5
        self.j_germline_end = self.j_allele.ungapped_len - self.j_trim_3
        self.junction_start = self.v_allele.anchor
        self.junction_end = self.v_allele.anchor + self.junction_length
        

    def _is_functional(self, sequence):
        """
       Evaluates whether the heavy chain sequence is functional based on the presence of in-frame junctions
       without stop codons and specific amino acid motifs at the junction boundaries.

       Args:
           sequence (str): The nucleotide sequence to evaluate for functionality.
       """
        self.functional = False
        self.stop_codon = self.check_stops(sequence)
        self.vj_in_frame = (self.junction_end % 3) == 0 and (self.junction_start % 3 == 0) and (self.junction_length % 3 == 0) and self.stop_codon is False
        self.note = ''
        if (self.junction_length % 3) == 0 and self.stop_codon is False:
            self.junction_aa = translate(self.junction)
            if self.junction_aa.startswith("C"):
                if self.junction_aa.endswith("F") or self.junction_aa.endswith("W"):
                    self.functional = True
                else:
                    self.note += 'J anchor (W/F) not present.'
            else:
                self.note += 'V second C not present.'

    def mutate(self, mutation_model: MutationModel):
        """
        Applies mutations to the heavy chain sequence using the given mutation model.

        Args:
            mutation_model (MutationModel): The mutation model to use for applying mutations to the sequence.

        Updates:
            The method updates the mutated sequence, mutations, mutation frequency, mutation count,
            junction sequence, and checks for functionality post-mutation.
        """
        mutated_sequence, mutations, mutation_rate = mutation_model.apply_mutation(self)
        self.mutated_seq = mutated_sequence
        self.mutations = mutations
        self.mutation_freq = mutation_rate
        self.mutation_count = len(mutations)
        self.junction = self.mutated_seq[self.v_allele.anchor:
                                          self.v_allele.anchor + self.junction_length].upper()
        # mutation metadata updates
        self._is_functional(self.mutated_seq)

    def get_junction_length(self):
        """
       Calculates the length of the junction region in the heavy chain sequence, taking into account
       the trimmed V, D, J segments and the NP1, NP2 regions.

       Returns:
           int: The total length of the junction region, including CDR3 and flanking conserved regions.
       """
        junction_length = self.v_allele.length - (
                self.v_allele.anchor - 1) - self.v_trim_3 + self.NP1_length + self.d_allele.length - \
                          self.d_trim_5 - self.d_trim_3 + \
                          self.NP2_length + (self.j_allele.anchor + 2) - self.j_trim_5
        return junction_length

    def check_stops(self, seq):
        """
        Checks the given sequence for the presence of stop codons.

        Args:
            seq (str): The nucleotide sequence to check for stop codons.

        Returns:
            bool: True if a stop codon is found, False otherwise.
        """
        stops = ["TAG", "TAA", "TGA"]
        for x in range(0, len(seq), 3):
            if seq[x:x + 3] in stops:
                return True
        return False

    @classmethod
    def create_random(cls, dataconfig: DataConfig):
        """
        Creates a random instance of `HeavyChainSequence` with randomly selected V, D, and J alleles from the given
        `DataConfig`.

        Args:
            dataconfig (DataConfig): A `DataConfig` instance providing allele choices and other configuration data.

        Returns:
            HeavyChainSequence: A new `HeavyChainSequence` instance with randomly chosen alleles and simulated sequence.
        """
        random_v_allele = random.choice([i for j in dataconfig.v_alleles for i in dataconfig.v_alleles[j]])
        random_d_allele = random.choice([i for j in dataconfig.d_alleles for i in dataconfig.d_alleles[j]])
        random_j_allele = random.choice([i for j in dataconfig.j_alleles for i in dataconfig.j_alleles[j]])
        return cls([random_v_allele, random_d_allele, random_j_allele], dataconfig)

    def __repr__(self):
        """
        Provides a textual representation of the heavy chain sequence, showing the V, D, and J segment positions
        and allele names within the overall sequence context.

        Returns:
            str: A string representation of the heavy chain sequence, including V, D, and J segments.
        """
        total_length = self.j_seq_end  # Assuming j_seq_end is the end of the sequence
        proportional_length = lambda start, end: int((end - start) / total_length * 100)  # Example scale factor

        v_length = proportional_length(self.v_seq_start, self.v_seq_end)
        d_length = proportional_length(self.d_seq_start, self.d_seq_end) if self.d_seq_start != self.d_seq_end else 0
        j_length = proportional_length(self.j_seq_start, self.j_seq_end)

        # Construct ASCII drawing
        v_part = f"{self.v_seq_start}|{'-' * v_length}V({self.v_allele.name})|{self.v_seq_end}"
        d_part = f"{self.d_seq_start}|{'-' * d_length}D({self.d_allele.name})|{self.d_seq_end}" if d_length > 0 else ""
        j_part = f"{self.j_seq_start}|{'-' * j_length}J({self.j_allele.name})|{self.j_seq_end}"

        return f"{v_part}{'|' if d_length > 0 else ''}{d_part}{'|' if d_length > 0 else ''}{j_part}"
