import pandas as pd
from GenAIRR.dataconfig.make.base_dataconfig_builder import BaseDataConfigGenerator
from GenAIRR.dataconfig.make.auxiliary_builders.correction_map_builder import CorrectionMapBuilder
from GenAIRR.dataconfig.make.auxiliary_builders.trimming_proportion_builder import TrimmingProbabilityGenerator
from GenAIRR.dataconfig.make.auxiliary_builders.np_markov_chain_builder import NPMarkovParameterBuilder
from GenAIRR.utilities.data_utilities import create_allele_dict
from GenAIRR.utilities.asc_utilities import create_asc_germline_set


class CustomDataConfigBuilder(BaseDataConfigGenerator):
    def __init__(self, convert_to_asc=True):
        super().__init__(convert_to_asc)
        self.correction_builder = CorrectionMapBuilder(self.dataconfig)

    def _derive_gene_usage(self, data):
        gene_use_dict = {
            'V': (data['v_call'].apply(lambda x: x.split('*')[0]).value_counts() / len(data)).to_dict(),
            'J': (data['j_call'].apply(lambda x: x.split('*')[0]).value_counts() / len(data)).to_dict()
        }
        if 'd_call' in data:
            gene_use_dict['D'] = (data['d_call'].apply(lambda x: x.split('*')[0]).value_counts() / len(data)).to_dict()
        self.dataconfig.gene_use_dict = gene_use_dict

    def _derive_trimming_proportions(self, data):
        for gene in ['v', 'd', 'j', 'c']:
            alleles = self.dataconfig.allele_list(gene)

            for trim_side in ['5', '3']:
                col_call = f'{gene}_call'
                col_trim = f'{gene}_trim_{trim_side}'

                for allele in alleles:
                    # Fallback to decaying probabilities if column missing or too few samples
                    if col_call not in data.columns:
                        trim_values = TrimmingProbabilityGenerator.generate_decaying_probabilities(50)
                    else:
                        samples = data[data[col_call].astype(str).str.contains(allele.gene)]
                        if len(samples) > 100:
                            trim_values = (samples[col_trim].value_counts() / len(samples)).to_dict()
                        else:
                            trim_values = TrimmingProbabilityGenerator.generate_decaying_probabilities(50)

                    gene_trim_dict = self.dataconfig.trim_dicts.setdefault(f'{gene.upper()}_{trim_side}', {})
                    family_dict = gene_trim_dict.setdefault(allele.family, {})
                    family_dict[allele.gene] = trim_values

    def make(self, v_reference_path, j_reference_path,custom_data,c_reference_path=None,
             d_reference_path=None):
        data = custom_data if type(custom_data) == pd.DataFrame else pd.read_csv(custom_data)
        self.has_d = d_reference_path is not None
        if self.has_d and 'D' not in self.alleles:
            self.alleles.append('D')

        # Load references
        v, j, c, d = self._read_reference_files(v_reference_path, j_reference_path, c_reference_path, d_reference_path)
        self._load_alleles(v_alleles=v, j_alleles=j, c_alleles=c, d_alleles=d)
        print('Alleles Mounted to DataConfig!...')

        # Gene usage
        self._derive_gene_usage(data)
        print('Gene Usage Mounted to DataConfig!...')

        # Trimming proportions
        self._load_trimming_probs()
        self._derive_trimming_proportions(data)
        print('Trimming Proportions Mounted to DataConfig!...')

        # NP parameters
        npgen = NPMarkovParameterBuilder(self.dataconfig, self.has_d)
        npgen.derive_all_from_data(data)
        print('NP Parameters Mounted to DataConfig!...')

        # Correction maps
        self.correction_builder.build_n_ambiguity_map(self.dataconfig.v_alleles)
        print('V Ns Ambiguity Map Mounted to DataConfig!...')

        self.correction_builder.build_3_prime_trim_map(self.dataconfig.v_alleles)
        self.correction_builder.build_5_prime_trim_map(self.dataconfig.v_alleles)
        self.correction_builder.build_3_prime_trim_map(self.dataconfig.j_alleles)
        self.correction_builder.build_5_prime_trim_map(self.dataconfig.j_alleles)
        if self.has_d and self.dataconfig.d_alleles:
            self.correction_builder.build_5_and_3_prime_trim_map(self.dataconfig.d_alleles)
            print('D (5,3) Prime Ambiguity Map Mounted to DataConfig!...')

        print('=' * 50)
        return self.dataconfig
