import logging

from GenAIRR.dataconfig.make.base_dataconfig_builder import BaseDataConfigGenerator
from GenAIRR.dataconfig.make.auxiliary_builders.trimming_proportion_builder import TrimmingProbabilityGenerator
from GenAIRR.dataconfig.make.auxiliary_builders.np_markov_chain_builder import NPMarkovParameterBuilder
from GenAIRR.utilities.data_utilities import create_allele_dict
from GenAIRR.utilities.asc_utilities import create_asc_germline_set

logger = logging.getLogger(__name__)


class RandomDataConfigBuilder(BaseDataConfigGenerator):
    """
    Generator that creates a DataConfig object with random (uniform/decaying) distributions
    for gene usage, trimming, and NP region parameters.
    """

    def __init__(self, convert_to_asc=False, *, species=None, chain_type=None, reference_set=None):
        super().__init__(
            convert_to_asc,
            species=species, chain_type=chain_type, reference_set=reference_set,
        )

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

    def make(self, v_reference_path, j_reference_path, c_reference_path=None,
             d_reference_path=None, *, v_anchor_finder=None, j_anchor_finder=None,
             keep_anchorless=False):
        """
        Construct and return a DataConfig with all randomly simulated components.

        Args:
            v_reference_path: V gene FASTA file path or preloaded dict.
            j_reference_path: J gene FASTA file path or preloaded dict.
            c_reference_path: C gene FASTA file path or preloaded dict (optional).
            d_reference_path: D gene FASTA file path or preloaded dict (optional).
            v_anchor_finder: Custom V anchor callable (see ``create_allele_dict``).
            j_anchor_finder: Custom J anchor callable (see ``create_allele_dict``).
            keep_anchorless: If True, keep V/J alleles without anchors
                (they can be sampled but never produce productive sequences).

        Returns:
            DataConfig with random distributions and validated metadata.
        """
        # Read and load reference sequences
        v, j, c, d = self._read_reference_files(
            v_reference_path, j_reference_path, c_reference_path, d_reference_path,
            v_anchor_finder=v_anchor_finder, j_anchor_finder=j_anchor_finder,
            keep_anchorless=keep_anchorless,
        )
        self._load_alleles(v_alleles=v, j_alleles=j, c_alleles=c, d_alleles=d)
        logger.info('Alleles mounted to DataConfig')

        # Random gene usage and trimming proportions
        self._load_random_gene_usage()
        logger.info('Random gene usage mounted to DataConfig')

        self._load_trimming_probs()
        logger.info('Random trimming proportions mounted to DataConfig')

        np_generator = NPMarkovParameterBuilder(self.dataconfig, has_d=self.has_d)
        np_generator.generate_all_random(max_length=50)
        logger.info('Random NP parameters mounted to DataConfig')

        # Correction maps are computed lazily during graph compilation
        self._finalize()
        logger.info('DataConfig build complete')

        return self.dataconfig
