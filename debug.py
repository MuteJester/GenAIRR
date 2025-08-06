import random
import logging
import time
import os
from dataclasses import dataclass
from itertools import product
from multiprocessing import cpu_count, Pool, Manager
from typing import List, Optional, Tuple, Any, Dict
import math

import pandas as pd
from tqdm.auto import tqdm

from GenAIRR.alleles import AlleleTypes
from GenAIRR.pipeline import AugmentationPipeline
from GenAIRR.steps import AugmentationStep
from GenAIRR.dataconfig import DataConfig
from GenAIRR.steps import (
    SimulateSequence, FixVPositionAfterTrimmingIndexAmbiguity, FixDPositionAfterTrimmingIndexAmbiguity,
    FixJPositionAfterTrimmingIndexAmbiguity, CorrectForVEndCut, CorrectForDTrims, CorruptSequenceBeginning,
    InsertNs, InsertIndels, ShortDValidation, DistillMutationRate
)
import pickle
from GenAIRR.mutation import S5F, Uniform
# Load data configuration
dataconfig_path = "C:/Users/tomas/Downloads/HUMAN_IGL_OGRDB.pkl"

with open(dataconfig_path, 'rb') as h:
    dataconfig = pickle.load(h)

AugmentationStep.set_dataconfig(dataconfig)
# Initialize simulator
is_productive = True  # Set to True for productive sequences
SimulateSequence.MAX_GENERATION_ATTEMPTS = 25

pipeline=AugmentationPipeline([
        SimulateSequence(S5F(min_mutation_rate=0.003, max_mutation_rate=0.2), is_productive),
        # FixVPositionAfterTrimmingIndexAmbiguity(),
        # FixJPositionAfterTrimmingIndexAmbiguity(),
        # CorrectForVEndCut(),
        # CorruptSequenceBeginning(0.7, [0.4, 0.4, 0.2], 576, 210, 310, 50),
        # InsertNs(0.02, 0.5),
        # InsertIndels(0.5, 5, 0.5, 0.5),
        # DistillMutationRate()
    ])

TOTAL = 900
prod_count = 0
results = []
for i in tqdm(range(TOTAL), desc="Simulating sequences"):
    if prod_count >= 100:
        break
    try:
        # Execute the pipeline to simulate a sequence
        seq = pipeline.execute()
        results.append(seq.get_dict())
        # Check if the sequence is productive
        if seq.productive():
            prod_count += 1
            print(f"Productive Sequence {prod_count}: {seq.sequence}")
    except Exception as e:
        print(f"Error during simulation: {e}")

#save as dataframe
#
pd.DataFrame(results).to_csv("C:/Users/tomas/Downloads/test_results.csv", index=False)