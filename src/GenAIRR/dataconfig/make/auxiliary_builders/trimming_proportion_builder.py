class TrimmingProbabilityGenerator:
    """
    Class responsible for generating and assigning random trimming proportions
    to V, D, J, and C alleles using decaying probability distributions.
    """

    def __init__(self, dataconfig, allele_types, reference_getter):
        """
        Initialize the generator.

        Args:
            dataconfig (DataConfig): The configuration object to populate.
            allele_types (list): List of allele types to process (e.g., ['V', 'J', 'C', 'D']).
            reference_getter (Callable): A function returning the reference dicts per allele type.
        """
        self.dataconfig = dataconfig
        self.allele_types = allele_types
        self.get_reference_pointers = reference_getter

    @staticmethod
    def generate_decaying_probabilities(n, base=0.8):
        """
        Generate a dictionary of exponentially decaying probabilities.

        Args:
            n (int): Maximum trimming amount.
            base (float): Decay base (default 0.8).

        Returns:
            dict: Mapping from trimming amount (0..n) to probability.
        """
        total = sum(base ** i for i in range(n + 1))
        return {i: (base ** i) / total for i in range(n + 1)}

    def load_trimming_proportions(self, max_trim=50):
        """
        Generate trimming proportions for each allele family and gene and store
        them in the dataconfig object.

        Args:
            max_trim (int): Maximum trimming amount.
        """
        pointer_to_reference = self.get_reference_pointers()
        trim_dicts = {}

        for allele_type in self.allele_types:
            alleles = {
                allele for gene in pointer_to_reference[allele_type]
                for allele in pointer_to_reference[allele_type][gene]
            }

            probabilities = self.generate_decaying_probabilities(max_trim)

            family_gene_map = {}
            for allele in alleles:
                if allele.family not in family_gene_map:
                    family_gene_map[allele.family] = {}
                family_gene_map[allele.family][allele.gene] = probabilities

            # Store for both 5' and 3' trims
            trim_dicts[f"{allele_type}_5"] = family_gene_map
            trim_dicts[f"{allele_type}_3"] = family_gene_map

        self.dataconfig.trim_dicts = trim_dicts
