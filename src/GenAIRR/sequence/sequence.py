from ..alleles import AlleleTypes
from abc import ABC, abstractmethod

from GenAIRR.dataconfig.data_config import DataConfig


class BaseSequence(ABC):
    """
    Abstract base class representing a recombined Immunoglobulin (Ig) sequence,
    consisting of variable (V), diversity (D), and joining (J) gene segments.

    This class provides a template for creating and manipulating Ig sequences, including
    features like (NP) regions, junctions, and mutations.

    Attributes:
        v_allele (Allele): The V gene allele.
        d_allele (Allele): The D gene allele, if present.
        j_allele (Allele): The J gene allele.
        NP1_region (str): The nucleotide sequence of the NP1 region between V and D genes.
        NP1_length (int): The length of the NP1 region.
        NP2_region (str): The nucleotide sequence of the NP2 region between D and J genes.
        NP2_length (int): The length of the NP2 region.
        ungapped_seq (str): The complete ungapped nucleotide sequence.
        gapped_seq (str): The complete nucleotide sequence with gaps.
        mutated_seq (str): The ungapped nucleotide sequence after applying mutations.
        gapped_mutated_seq (str): The gapped nucleotide sequence after applying mutations.
        junction (str): The nucleotide sequence of the junction region (CDR3 plus anchors).
        v_seq_start (int): The start position of the V segment in the sequence.
        d_seq_start (int): The start position of the D segment in the sequence, if present.
        j_seq_start (int): The start position of the J segment in the sequence.
        v_seq_end (int): The end position of the V segment in the sequence.
        d_seq_end (int): The end position of the D segment in the sequence, if present.
        j_seq_end (int): The end position of the J segment in the sequence.
        mutations (str): A string representation of mutation events.
        mutation_count (int): The total count of mutations.
        mutation_freq (float): The frequency of mutations in the sequence.

    Args:
        alleles (List[Allele]): A list of  V, D (optional), and J gene alleles.
    """

    def __init__(self, alleles):
        """
        Initializes a BaseSequence instance with given alleles for V, D (optional), and J segments.
        """
        self.v_allele = next((allele for allele in alleles if allele.type == AlleleTypes.V), None)
        assert self.v_allele is not None  # Must Have V Allele!
        self.d_allele = next((allele for allele in alleles if allele.type == AlleleTypes.D), None)
        self.j_allele = next((allele for allele in alleles if allele.type == AlleleTypes.J), None)
        assert self.j_allele is not None  # Must Have J Allele!
        self.c_allele = next((allele for allele in alleles if allele.type == AlleleTypes.C), None)
        self.NP1_region = ""
        self.NP2_region = ""
        self.NP1_length = 0
        self.NP2_length = 0
        self.junction = ""
        self.v_seq_start = 0
        self.d_seq_start = 0
        self.j_seq_start = 0
        self.v_seq_end = 0
        self.d_seq_end = 0
        self.j_seq_end = 0
        self.mutations = ""
        self.mutation_count = 0
        self.mutation_freq = 0
        self.functional = False
        self.seq_stop_codon = False
        self.ungapped_seq = ""
        self.mutated_seq = None
        self.gapped_seq = ""
        self.gapped_mutated_seq = None

    @abstractmethod
    def simulate_sequence(self, dataconfig: DataConfig):
        """
        Abstract method to simulate the recombined nucleotide sequence including V(D)J recombination,
        trimming, and NP region addition based on provided data configurations.

        Args:
            dataconfig (DataConfig): An instance of DataConfig containing all necessary parameters
            and configurations for sequence simulation.

        Returns:
            str: The simulated nucleotide sequence, optionally including gaps.
        """
        pass

    @abstractmethod
    def get_junction_length(self):
        """
        Abstract method to calculate the length of the junction region of the sequence.
        The junction region typically includes the CDR3 region along with conserved residues
        from V and J gene segments.

        Returns:
            int: The length of the junction region in nucleotides.
        """
        pass



