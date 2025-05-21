import random
from ...mutation import MutationModel
from ..sequence import NP_Region
from ...sequence.sequence import BaseSequence
from ...utilities import translate
from GenAIRR.dataconfig.data_config import DataConfig


class TCRHeavyChainSequence(BaseSequence):
    """
        Represents a TCR heavy chain sequence in an immunoglobulin, including V, D, and J segments,
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
        self.simulate_trimmed_sequences(dataconfig)
        self.simulate_NP_regions(dataconfig)
        self.assemble_sequence()
        self.calculate_junction_properties()
        self.update_metadata()
        self.check_functionality()

    def simulate_NP_regions(self, dataconfig: DataConfig):
        self.NP1_region = NP_Region.create_np_region(dataconfig.NP_lengths, dataconfig.NP_transitions, "NP1",
                                                     dataconfig.NP_first_bases)
        self.NP2_region = NP_Region.create_np_region(dataconfig.NP_lengths, dataconfig.NP_transitions, "NP2",
                                                     dataconfig.NP_first_bases)
        self.NP1_length = len(self.NP1_region)
        self.NP2_length = len(self.NP2_region)

    def simulate_trimmed_sequences(self, dataconfig: DataConfig):
        self.v_trimmed_seq, self.v_trim_5, self.v_trim_3 = self.v_allele.get_trimmed(dataconfig.trim_dicts)
        self.d_trimmed_seq, self.d_trim_5, self.d_trim_3 = self.d_allele.get_trimmed(dataconfig.trim_dicts)

        self.j_trimmed_seq, self.j_trim_5, self.j_trim_3 = self.j_allele.get_trimmed(dataconfig.trim_dicts)
        #self.c_trimmed_seq, self.c_trim_5, self.c_trim_3 = self.c_allele.get_trimmed(dataconfig.trim_dicts)

    def assemble_sequence(self):
        self.ungapped_seq = (
                self.v_trimmed_seq
                + self.NP1_region
                + self.d_trimmed_seq
                + self.NP2_region
                + self.j_trimmed_seq
                #+ self.c_trimmed_seq
        ).upper()
    def calculate_junction_properties(self):
        self.junction_length = self.get_junction_length()
        self.junction = self.ungapped_seq[self.v_allele.anchor:self.v_allele.anchor + self.junction_length].upper()

    def update_metadata(self):
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
        self.check_functionality(self.mutated_seq)

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
    def create_random(cls, dataconfig: DataConfig,specific_v=None,specific_d=None,specific_j=None):
        """
        Creates a random instance of `HeavyChainSequence` with randomly selected V, D, and J alleles from the given
        `DataConfig`.

        Args:
            dataconfig (DataConfig): A `DataConfig` instance providing allele choices and other configuration data.

        Returns:
            HeavyChainSequence: A new `HeavyChainSequence` instance with randomly chosen alleles and simulated sequence.
        """


        if specific_v is None:
            random_v_allele = random.choice([i for j in dataconfig.v_alleles for i in dataconfig.v_alleles[j]])
        else:
            random_v_allele = specific_v



        if specific_d is None:
            random_d_allele = random.choice([i for j in dataconfig.d_alleles for i in dataconfig.d_alleles[j]])
        else:
            random_d_allele = specific_d

        d_name = random_d_allele.name
        group = d_name.split('*')[0].replace('TRBD','')



        if specific_j is None:
            all_j_alleles = [i for j in dataconfig.j_alleles for i in dataconfig.j_alleles[j]]
            filtered_j_alleles = list(filter(lambda x: 'TRBJ' + group in x.name, all_j_alleles))

            random_j_allele = random.choice(filtered_j_alleles)
        else:
            random_j_allele = specific_j

        # deprecated on 2025-09-04
        #all_c_alleles = [i for j in dataconfig.c_alleles for i in dataconfig.c_alleles[j]]
        #filtered_c_alleles = list(filter(lambda x: 'TRBC' + group in x.name, all_c_alleles))

        #random_c_allele = random.choice(filtered_c_alleles)

        return cls([random_v_allele, random_d_allele, random_j_allele], dataconfig)

    def __repr__(self):
        """
        Provides a textual representation of the heavy chain sequence, showing the V, D, and J segment positions
        and allele names within the overall sequence context.

        Returns:
            str: A string representation of the heavy chain sequence, including V, D, and J segments.
        """
        total_length = len(self.ungapped_seq)  # Assuming j_seq_end is the end of the sequence
        proportional_length = lambda start, end: int((end - start) / total_length * 100)  # Example scale factor

        v_length = proportional_length(self.v_seq_start, self.v_seq_end)
        d_length = proportional_length(self.d_seq_start, self.d_seq_end) if self.d_seq_start != self.d_seq_end else 0
        j_length = proportional_length(self.j_seq_start, self.j_seq_end)
        c_length = proportional_length(self.j_seq_end, total_length)

        # Construct ASCII drawing
        v_part = f"{self.v_seq_start}|{'-' * v_length}V({self.v_allele.name})|{self.v_seq_end}"
        d_part = f"{self.d_seq_start}|{'-' * d_length}D({self.d_allele.name})|{self.d_seq_end}" if d_length > 0 else ""
        j_part = f"{self.j_seq_start}|{'-' * j_length}J({self.j_allele.name})|{self.j_seq_end}"
        #c_part = f"{self.j_seq_end}|{'-' * c_length}C({self.c_allele.name})|{total_length}"

        return f"{v_part}{'|' if d_length > 0 else ''}{d_part}{'|' if d_length > 0 else ''}{j_part}"#|{c_part}"
