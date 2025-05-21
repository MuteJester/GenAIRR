from GenAIRR.dataconfig.make.base_dataconfig_builder import BaseDataConfigGenerator
from GenAIRR.dataconfig.make.auxiliary_builders.correction_map_builder import CorrectionMapBuilder
from GenAIRR.dataconfig.make.auxiliary_builders.trimming_proportion_builder import TrimmingProbabilityGenerator
from GenAIRR.dataconfig.make.auxiliary_builders.np_markov_chain_builder import NPMarkovParameterBuilder
from GenAIRR.utilities.AlleleNComparer import AlleleNComparer
from GenAIRR.utilities.data_utilities import create_allele_dict
from GenAIRR.utilities.asc_utilities import create_asc_germline_set


class RandomDataConfigBuilder(BaseDataConfigGenerator):
    """
    Generator that creates a DataConfig object with random (uniform/decaying) distributions
    for gene usage, trimming, and NP region parameters.
    """

    def _load_random_gene_usage(self):
        """
        Assign uniform probabilities to gene usage across all known alleles.
        """
        pointer_to_reference = self._get_reference_pointers()

        gene_use_dict = {allele: {} for allele in self.alleles}
        if self.has_d:
            gene_use_dict['D'] = {}

        for allele_type in gene_use_dict:
            current_alleles = list(pointer_to_reference[allele_type])
            n = len(current_alleles)
            gene_use_dict[allele_type] = {name: 1 / n for name in current_alleles}

        self.dataconfig.gene_use_dict = gene_use_dict

    def make(self, v_reference_path, j_reference_path, c_reference_path, d_reference_path=None):
        """
        Main method to construct and return a DataConfig with all randomly simulated components.
        """
        # Read and load reference sequences
        v, j, c, d = self._read_reference_files(v_reference_path, j_reference_path, c_reference_path, d_reference_path)
        self._load_alleles(v_alleles=v, j_alleles=j, c_alleles=c, d_alleles=d)
        print('Alleles Mounted to DataConfig!...')

        # Random gene usage and trimming proportions
        self._load_random_gene_usage()
        print('Random Gene Usage Mounted to DataConfig!...')

        self._load_trimming_probs()
        print('Random Trimming Proportions Mounted to DataConfig!...')

        np_generator = NPMarkovParameterBuilder(self.dataconfig, has_d=self.has_d)
        np_generator.generate_all_random(max_length=50)
        print('Random NP Parameters Mounted to DataConfig!...')

        # Derive correction maps
        self.correction_builder.build_n_ambiguity_map(self.dataconfig.v_alleles)
        print('V Ns Ambiguity Map Mounted to DataConfig!...')

        self.correction_builder.build_5_prime_trim_map(self.dataconfig.v_alleles)
        self.correction_builder.build_3_prime_trim_map(self.dataconfig.v_alleles)
        self.correction_builder.build_5_prime_trim_map(self.dataconfig.j_alleles)
        self.correction_builder.build_3_prime_trim_map(self.dataconfig.j_alleles)

        if self.has_d and self.dataconfig.d_alleles:
            self.correction_builder.build_5_and_3_prime_trim_map(self.dataconfig.d_alleles)

        print('Allele Trim Ambiguity Maps Mounted to DataConfig!...')
        print('=' * 50)

        return self.dataconfig
