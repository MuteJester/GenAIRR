import random
from dataclasses import dataclass
from itertools import product
from multiprocessing import cpu_count, Pool

import pandas as pd
from tqdm.auto import tqdm

from GenAIRR.alleles import AlleleTypes
from GenAIRR.pipeline import AugmentationPipeline
from GenAIRR.steps import AugmentationStep
from GenAIRR.dataconfig import DataConfig

from GenAIRR.data import builtin_heavy_chain_data_config, builtin_kappa_chain_data_config
from GenAIRR.pipeline import CHAIN_TYPE_BCR_HEAVY
from GenAIRR.steps import (
    SimulateSequence, FixVPositionAfterTrimmingIndexAmbiguity, FixDPositionAfterTrimmingIndexAmbiguity,
    FixJPositionAfterTrimmingIndexAmbiguity, CorrectForVEndCut, CorrectForDTrims, CorruptSequenceBeginning,
    InsertNs, InsertIndels, ShortDValidation, DistillMutationRate
)
from GenAIRR.mutation import S5F, Uniform


def process_combinations_worker(args):
    combinations, _alleles, dataconfig, chain_type, pipeline_config = args

    AugmentationStep.set_dataconfig(dataconfig, chain_type=chain_type)
    pipeline = AugmentationPipeline(pipeline_config)

    def mount_genes(v=None, d=None, j=None):
        pipeline[0].specific_v = v
        pipeline[0].specific_d = d
        pipeline[0].specific_j = j

    def simulate(v, d, j):
        def specific():
            mount_genes(v, d, j)
            gen = pipeline.execute().get_dict()
            mount_genes()
            return gen

        if d != "Short-D":
            gen = specific()
            t = 0
            while "Short-D" in gen["d_call"] and t < 10:
                gen = specific()
                t += 1
            return gen
        else:
            mount_genes(v=v, j=j)
            gen = pipeline.execute().get_dict()
            while "Short-D" not in gen["d_call"]:
                gen = pipeline.execute().get_dict()
            mount_genes()
            return gen

    dataset = []
    for combination in combinations:
        if len(combination) == 2:
            if type(combination[0]) == str or combination[0].type == AlleleTypes.D:
                D, J = combination
                v, d, j = random.choice(_alleles), D, J
            else:
                V, D = combination
                v, d, j = V, D, random.choice(_alleles)
        else:
            v, d, j = combination
        dataset.append(simulate(v, d, j))
    return dataset


class UnbiasedDatasetSimulator:
    def __init__(self, k, dataconfig, chain_type, pipeline=None):
        self.K = k
        self.pipeline = None
        self.dataconfig = dataconfig
        self.chain_type = chain_type
        AugmentationStep.set_dataconfig(dataconfig, chain_type=chain_type)
        self._set_pipeline(pipeline)

    def _set_pipeline(self, pipeline):
        if pipeline:
            self.pipeline = pipeline
        else:
            self.pipeline = AugmentationPipeline([
                SimulateSequence(mutation_model=Uniform(min_mutation_rate=0.003, max_mutation_rate=0.01), productive=True),
                FixVPositionAfterTrimmingIndexAmbiguity(),
                FixDPositionAfterTrimmingIndexAmbiguity(),
                FixJPositionAfterTrimmingIndexAmbiguity(),
                CorrectForVEndCut(),
                CorrectForDTrims(),
                CorruptSequenceBeginning(
                    corruption_probability=0.7,
                    corrupt_events_proba=[0.4, 0.4, 0.2],
                    max_sequence_length=576,
                    nucleotide_add_coefficient=110,
                    nucleotide_remove_coefficient=110,
                    nucleotide_add_after_remove_coefficient=20,
                    random_sequence_add_proba=1,
                    single_base_stream_proba=0,
                    duplicate_leading_proba=0,
                    random_allele_proba=0
                ),
                InsertNs(n_ratio=0.02, proba=0.5),
                ShortDValidation(short_d_length=5),
                InsertIndels(indel_probability=0.5, max_indels=5, insertion_proba=0.5, deletion_proba=0.5),
                DistillMutationRate()
            ])

    def simulate(self):
        V_alleles = [i for j in self.dataconfig.v_alleles for i in self.dataconfig.v_alleles[j]]
        D_alleles = [i for j in self.dataconfig.d_alleles for i in self.dataconfig.d_alleles[j]]
        J_alleles = [i for j in self.dataconfig.j_alleles for i in self.dataconfig.j_alleles[j]]

        if self.pipeline[0].productive:
            def check_stops(seq):
                stops = ["TAG", "TAA", "TGA"]
                return any(seq[x:x + 3] in stops for x in range(0, len(seq), 3))

            V_alleles = [i for i in V_alleles if not check_stops(i.ungapped_seq.upper())]

        D_alleles.append("Short-D")
        V_alleles *= self.K
        D_alleles *= self.K
        J_alleles *= self.K

        random.shuffle(V_alleles)
        random.shuffle(D_alleles)
        random.shuffle(J_alleles)

        V_D_combinations = list(product(V_alleles, D_alleles))
        D_J_combinations = list(product(D_alleles, J_alleles))
        V_D_J_combinations = list(product(V_alleles, D_alleles, J_alleles))

        num_cores = cpu_count()
        chunk_size_VD = len(V_D_combinations) // num_cores
        chunk_size_DJ = len(D_J_combinations) // num_cores
        chunk_size_VDJ = len(V_D_J_combinations) // num_cores

        pipeline_config = self.pipeline.steps

        args_list = [
            (V_D_combinations[i:i + chunk_size_VD], J_alleles, self.dataconfig, self.chain_type, pipeline_config)
            for i in range(0, len(V_D_combinations), chunk_size_VD)
        ] + [
            (D_J_combinations[i:i + chunk_size_DJ], V_alleles, self.dataconfig, self.chain_type, pipeline_config)
            for i in range(0, len(D_J_combinations), chunk_size_DJ)
        ] + [
            (V_D_J_combinations[i:i + chunk_size_VDJ], None, self.dataconfig, self.chain_type, pipeline_config)
            for i in range(0, len(V_D_J_combinations), chunk_size_VDJ)
        ]

        with Pool(num_cores) as pool:
            results = list(
                tqdm(pool.imap(process_combinations_worker, args_list), total=len(args_list), desc='Generating dataset')
            )

        dataset = [item for sublist in results for item in sublist]
        return pd.DataFrame(dataset)
#


#Sample Use Case
# if __name__ == '__main__':
#     from GenAIRR.data import builtin_tcrb_data_config
#     from GenAIRR.pipeline.pipeline_parameters import CHAIN_TYPE_TCR_BETA
#
#     dataconfig = builtin_tcrb_data_config()
#     ub_simulator = UnbiasedDatasetSimulator(k=2, dataconfig=dataconfig, chain_type=CHAIN_TYPE_TCR_BETA)
#     rdf = ub_simulator.simulate()
#     for allele in ['v','d','j']:
#         rdf[f'{allele}_call'] = rdf[f'{allele}_call'].str.join(',')
#     rdf.to_csv('path')
#     print(rdf)
