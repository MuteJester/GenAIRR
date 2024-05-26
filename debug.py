from GenAIRR.utilities.data_config import DataConfig
from GenAIRR.mutation import S5F, Uniform
from GenAIRR.simulation import HeavyChainSequenceAugmentor, SequenceAugmentorArguments
from GenAIRR.simulation import LightChainKappaLambdaSequenceAugmentor
from GenAIRR.simulation import HeavyChainSequenceAugmentor, SequenceAugmentorArguments
from GenAIRR.data import builtin_kappa_chain_data_config, builtin_lambda_chain_data_config

kappa_config = builtin_kappa_chain_data_config()
lambda_config = builtin_lambda_chain_data_config()

# Create initial CSV with headers
args = SequenceAugmentorArguments(mutation_model=S5F, custom_mutation_model_path=None, simulate_indels=0.5,
                                  corrupt_proba=0, productive=True)

simulator = LightChainKappaLambdaSequenceAugmentor(lambda_args=args, kappa_args=args,
                                                   lambda_dataconfig=lambda_config, kappa_dataconfig=kappa_config)

from tqdm.auto import tqdm
for i in tqdm(range(1000)):
    simulator.simulate_augmented_sequence()