
class DataConfig:
    """
        Configuration class for storing data related to sequence generation, allele usage, trimming, and mutation rates.

        This class encapsulates various dictionaries and settings used in the simulation and analysis of immunoglobulin sequences, including information about gene usage, trimming preferences, non-polymorphic (NP) region transitions, k-mer dictionaries, and allele-specific corrections.

        Attributes:
            family_use_dict (dict): A dictionary specifying the usage frequencies of gene families.
            gene_use_dict (dict): A dictionary specifying the usage frequencies of individual genes.
            trim_dicts (dict): A dictionary containing information on how to trim each gene segment (V, D, J).
            NP_transitions (dict): A dictionary detailing the transition probabilities for bases in NP regions.
            NP_first_bases (dict): A dictionary specifying the probabilities of the first base in NP regions.
            NP_lengths (dict): A dictionary defining the distribution of lengths for NP regions.
            mut_rate_per_seq (dict): A dictionary specifying the mutation rate per sequence for simulating mutations.
            kmer_dicts (dict): A dictionary of k-mer dictionaries used for various analyses or simulations.
            v_alleles, d_alleles, j_alleles: family seperated dictionary of lists for storing allele information for V, D, and J gene segments, respectively.
            correction_maps (dict): A dictionary for storing maps used for correction or adjustment of sequences or simulation parameters.
            asc_tables (dict): A dictionary for storing allele sequence cluster (ASC) tables, which group alleles based on sequence similarity and other criteria.
    """
    def __init__(self,):

        # Config Variables
        self.family_use_dict = {}
        self.gene_use_dict = {}
        self.trim_dicts = {}
        self.NP_transitions = {}
        self.NP_first_bases = {}
        self.NP_lengths = {}
        self.mut_rate_per_seq = {}
        self.kmer_dicts = {}
        self.v_alleles = None
        self.d_alleles = None
        self.j_alleles = None

        self.correction_maps = dict()
        self.asc_tables = dict()
