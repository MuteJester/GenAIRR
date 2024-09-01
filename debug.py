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
ag = HeavyChainSequenceAugmentor(dc,SequenceAugmentorArguments(productive=False))
v_alleles = {i.name:i for j in dc.v_alleles for i in dc.v_alleles[j]}
d_alleles = {i.name:i for j in dc.d_alleles for i in dc.d_alleles[j]}
j_alleles = {i.name:i for j in dc.j_alleles for i in dc.j_alleles[j]}
c_alleles = {i.name:i for j in dc.c_alleles for i in dc.c_alleles[j]}


for _ in range(1000):
    sim = ag.simulate_augmented_sequence()

