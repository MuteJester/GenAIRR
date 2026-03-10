import csv
import logging
from collections import Counter

from GenAIRR.dataconfig.make.base_dataconfig_builder import BaseDataConfigGenerator
from GenAIRR.dataconfig.make.auxiliary_builders.trimming_proportion_builder import TrimmingProbabilityGenerator
from GenAIRR.dataconfig.make.auxiliary_builders.np_markov_chain_builder import NPMarkovParameterBuilder
from GenAIRR.utilities.data_utilities import create_allele_dict
from GenAIRR.utilities.asc_utilities import create_asc_germline_set

logger = logging.getLogger(__name__)


def _load_airr_data(custom_data):
    """Load AIRR data from a CSV path or list of dicts."""
    if isinstance(custom_data, str):
        with open(custom_data, newline='') as f:
            return list(csv.DictReader(f))
    elif isinstance(custom_data, list):
        return custom_data
    else:
        # Backwards compat: accept pandas DataFrame
        try:
            return custom_data.to_dict(orient='records')
        except AttributeError:
            raise TypeError(
                "custom_data must be a CSV file path, a list of dicts, "
                "or a pandas DataFrame"
            )


def _counter_to_prob_dict(counter):
    """Convert a Counter to a {value: probability} dict."""
    total = sum(counter.values())
    if total == 0:
        return {}
    return {k: v / total for k, v in counter.items()}


class CustomDataConfigBuilder(BaseDataConfigGenerator):
    def __init__(self, convert_to_asc=True, *, species=None, chain_type=None, reference_set=None):
        super().__init__(
            convert_to_asc,
            species=species, chain_type=chain_type, reference_set=reference_set,
        )

    def _derive_gene_usage(self, data):
        v_counts = Counter(row['v_call'].split('*')[0] for row in data)
        j_counts = Counter(row['j_call'].split('*')[0] for row in data)

        gene_use_dict = {
            'V': _counter_to_prob_dict(v_counts),
            'J': _counter_to_prob_dict(j_counts),
        }

        # Check if d_call column exists in the data
        if data and 'd_call' in data[0]:
            d_counts = Counter(row['d_call'].split('*')[0] for row in data)
            gene_use_dict['D'] = _counter_to_prob_dict(d_counts)

        self.dataconfig.gene_use_dict = gene_use_dict

    def _derive_trimming_proportions(self, data):
        for gene in ['v', 'd', 'j', 'c']:
            alleles = self.dataconfig.allele_list(gene)

            for trim_side in ['5', '3']:
                col_call = f'{gene}_call'
                col_trim = f'{gene}_trim_{trim_side}'

                for allele in alleles:
                    # Check if the call column exists in the data
                    has_col = data and col_call in data[0]

                    if not has_col:
                        trim_values = TrimmingProbabilityGenerator.generate_decaying_probabilities(50)
                    else:
                        # Filter rows where the call column contains the allele gene name
                        samples = [
                            row for row in data
                            if allele.gene in str(row.get(col_call, ''))
                        ]
                        if len(samples) > 100:
                            trim_counts = Counter(
                                int(float(row[col_trim])) for row in samples
                                if col_trim in row and row[col_trim] not in (None, '', 'NA')
                            )
                            trim_values = _counter_to_prob_dict(trim_counts)
                        else:
                            trim_values = TrimmingProbabilityGenerator.generate_decaying_probabilities(50)

                    gene_trim_dict = self.dataconfig.trim_dicts.setdefault(f'{gene.upper()}_{trim_side}', {})
                    family_dict = gene_trim_dict.setdefault(allele.family, {})
                    family_dict[allele.gene] = trim_values

    def make(self, v_reference_path, j_reference_path, custom_data, c_reference_path=None,
             d_reference_path=None, *, v_anchor_finder=None, j_anchor_finder=None,
             keep_anchorless=False):
        """
        Construct a DataConfig derived from empirical AIRR data.

        Args:
            v_reference_path: V gene FASTA file path or preloaded dict.
            j_reference_path: J gene FASTA file path or preloaded dict.
            custom_data: AIRR CSV path, list of dicts, or pandas DataFrame.
            c_reference_path: C gene FASTA file path or preloaded dict (optional).
            d_reference_path: D gene FASTA file path or preloaded dict (optional).
            v_anchor_finder: Custom V anchor callable (see ``create_allele_dict``).
            j_anchor_finder: Custom J anchor callable (see ``create_allele_dict``).
            keep_anchorless: If True, keep V/J alleles without anchors
                (they can be sampled but never produce productive sequences).

        Returns:
            DataConfig with empirically derived distributions and validated metadata.
        """
        data = _load_airr_data(custom_data)
        self.has_d = d_reference_path is not None
        if self.has_d and 'D' not in self.alleles:
            self.alleles.append('D')

        # Load references
        v, j, c, d = self._read_reference_files(
            v_reference_path, j_reference_path, c_reference_path, d_reference_path,
            v_anchor_finder=v_anchor_finder, j_anchor_finder=j_anchor_finder,
            keep_anchorless=keep_anchorless,
        )
        self._load_alleles(v_alleles=v, j_alleles=j, c_alleles=c, d_alleles=d)
        logger.info('Alleles mounted to DataConfig')

        # Gene usage
        self._derive_gene_usage(data)
        logger.info('Gene usage mounted to DataConfig')

        # Trimming proportions
        self._load_trimming_probs()
        self._derive_trimming_proportions(data)
        logger.info('Trimming proportions mounted to DataConfig')

        # NP parameters
        npgen = NPMarkovParameterBuilder(self.dataconfig, self.has_d)
        npgen.derive_all_from_data(data)
        logger.info('NP parameters mounted to DataConfig')

        # Correction maps are computed lazily during graph compilation
        self._finalize()
        logger.info('DataConfig build complete')

        return self.dataconfig
