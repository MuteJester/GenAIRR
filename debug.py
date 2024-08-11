from GenAIRR.utilities.data_config import DataConfig
from GenAIRR.mutation import S5F, Uniform
from GenAIRR.simulation import HeavyChainSequenceAugmentor, SequenceAugmentorArguments
from GenAIRR.simulation import LightChainKappaLambdaSequenceAugmentor
from GenAIRR.simulation import HeavyChainSequenceAugmentor, SequenceAugmentorArguments
from GenAIRR.data import builtin_kappa_chain_data_config, builtin_lambda_chain_data_config,builtin_heavy_chain_data_config

kappa_config = builtin_kappa_chain_data_config()
lambda_config = builtin_lambda_chain_data_config()

# Create initial CSV with headers
# args = SequenceAugmentorArguments(mutation_model=Uniform, custom_mutation_model_path=None, simulate_indels=0.5,
#                                   corrupt_proba=0, productive=True,min_mutation_rate=0,max_mutation_rate=0)

# simulator = LightChainKappaLambdaSequenceAugmentor(lambda_args=args, kappa_args=args,
#                                                    lambda_dataconfig=lambda_config, kappa_dataconfig=kappa_config)
dc = builtin_heavy_chain_data_config()
ag = HeavyChainSequenceAugmentor(dc,SequenceAugmentorArguments(productive=True,corrupt_proba=0,insertion_proba=0))
# simulator = HeavyChainSequenceAugmentor(dataconfig=dc,args=args)
# d_alleles = [i.ungapped_seq.upper() for j in dc.d_alleles for i in dc.d_alleles[j]]
#
# from tqdm.auto import tqdm
# for i in tqdm(range(100000)):
#     simulated = simulator.simulate_augmented_sequence()
#     n_d_calls = len(simulated['d_call'].split(','))
#     if (n_d_calls) > 10:
#         print('\n=================================')
#         print(simulated)
#         print(' Print D Length: ',simulated['d_sequence_end']-simulated['d_sequence_start'])
#         print(f'Trims: 5: {simulated["d_trim_5"]}  | 3: {simulated["d_trim_3"]}')
#         sub_string = simulated['sequence'][simulated['d_sequence_start']:simulated['d_sequence_end']]
#         print('Sustring in: ',sum([sub_string in i for i in d_alleles]))
#
#         print('=================================\n')

v_alleles = {i.name:i for j in dc.v_alleles for i in dc.v_alleles[j]}
d_alleles = {i.name:i for j in dc.d_alleles for i in dc.d_alleles[j]}
j_alleles = {i.name:i for j in dc.j_alleles for i in dc.j_alleles[j]}
c_alleles = {i.name:i for j in dc.c_alleles for i in dc.c_alleles[j]}



sim = ag.simulate_augmented_sequence(specific_v=v_alleles['IGHVF9-G32*02'])
sim
