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

    def __init__(self):

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
        self.c_alleles = None
        self.name = None

        self.correction_maps = dict()
        self.asc_tables = dict()


    def _unfold_alleles(self,gene):
        """
        Unfolds the alleles for a given gene family into a flat list.
        """
        # check if there an attribite self.{gene}_alleles
        if not hasattr(self, f"{gene}_alleles"):
            raise ValueError(f"Gene family '{gene}' not found in DataConfig.")
        alleles = getattr(self, f"{gene}_alleles")
        if alleles is None:
            return []
        else:
            alleles =[i for j in alleles for i in alleles[j]]
        return alleles
    def _count_alleles(self,gene):
        """
        Counts the number of alleles for a given gene family.
        """
        return len(self._unfold_alleles(gene))


    @property
    def number_of_v_alleles(self):
        """
        Returns the number of V alleles in the data configuration.
        """
        return self._count_alleles('v')
    @property
    def number_of_d_alleles(self):
        """
        Returns the number of D alleles in the data configuration.
        """
        return self._count_alleles('d')
    @property
    def number_of_j_alleles(self):
        """
        Returns the number of J alleles in the data configuration.
        """
        return self._count_alleles('j')

    def allele_list(self,gene):
        """
        Returns a list of alleles for a given gene family.
        """
        return self._unfold_alleles(gene)

    def __repr__(self):
        fmst = f"<{self.name} - Data Config>"
        if self.v_alleles is not None:
            fmst += f'-<{self._count_alleles("v")} V Alleles>'
        if self.d_alleles is not None:
            fmst += f'-<{self._count_alleles("d")} D Alleles>'
        if self.j_alleles is not None:
            fmst += f'-<{self._count_alleles("j")} J Alleles>'
        if self.c_alleles is not None:
            fmst += f'-<{self._count_alleles("c")} C Alleles>'
        return fmst
