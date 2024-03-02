import io
import pickle
from importlib import resources

from GenAIRR.utilities import DataConfig

#
# class RenameUnpickler(pickle.Unpickler):
#     def find_class(self, module, name):
#         renamed_module = module
#         if module.startswith("SequenceSimulation"):
#             # Update the prefix to your new module path
#             renamed_module = module.replace("SequenceSimulation", "GenAIRR")
#
#         return super(RenameUnpickler, self).find_class(renamed_module, name)
#
#
# def renamed_load(file_obj):
#     return RenameUnpickler(file_obj).load()
#
#
# def renamed_loads(pickled_bytes):
#     file_obj = io.BytesIO(pickled_bytes)
#     return renamed_load(file_obj)
#
#
#
# with resources.path('GenAIRR.data', 'LightChain_KAPPA_DataConfigV2.pkl') as data_path:
#     with open(data_path, 'rb') as h:
#         config = renamed_load(h)
#
# new_config = DataConfig()
#
# new_config.family_use_dict = config.family_use_dict
# new_config.gene_use_dict = config.gene_use_dict
# new_config.trim_dicts = config.trim_dicts
# new_config.NP_transitions = config.NP_transitions
# new_config.NP_first_bases = config.NP_first_bases
# new_config.NP_lengths = config.NP_lengths
# new_config.mut_rate_per_seq = config.mut_rate_per_seq
# new_config.kmer_dicts = config.kmer_dicts
# new_config.v_alleles = config.v_alleles
# new_config.d_alleles = config.d_alleles
# new_config.j_alleles = config.j_alleles
# new_config.correction_maps = config.correction_maps
# new_config.asc_tables = config.asc_tables
#
# print(dir(new_config))
# with resources.path('GenAIRR.data', 'LightChain_KAPPA_DataConfigV2.pkl') as data_path:
#     with open(data_path, 'wb') as h:
#         pickle.dump(new_config,h)


