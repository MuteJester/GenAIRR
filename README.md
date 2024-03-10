# GenAIRR: AIRR Sequence Simulator

GenAIRR is a Python module designed to generate synthetic Adaptive Immune Receptor Repertoire (AIRR) sequences for the purpose of benchmarking alignment algorithms and conducting sequence analysis in a non-biased manner.


- **Realistic Sequence Simulation**: Generate heavy and light immunoglobulin chain sequences with extensive customization options.
- **Advanced Mutation and Augmentation**: Introduce mutations and augment sequences to closely mimic the natural diversity and sequencing artifacts.
- **Precision in Allele-Specific Corrections**: Utilize sophisticated correction maps to accurately handle allele-specific trimming events and ambiguities.
- **Indel Simulation Capability**: Reflect the intricacies of sequencing data by simulating insertions and deletions within sequences.

# Quick Start Guide to GenAIRR

Welcome to the Quick Start Guide for GenAIRR, a Python module designed for generating synthetic Adaptive Immune Receptor Repertoire (AIRR) sequences. This guide will walk you through the basic usage of GenAIRR, including setting up your environment, simulating heavy and light chain sequences, and customizing your simulations.


## Installation

Before you begin, ensure that you have Python 3.x installed on your system. GenAIRR can be installed using pip, Python's package installer. Execute the following command in your terminal:



```python
import pandas as pd
# Install GenAIRR using pip
#!pip install GenAIRR
```

## Setting Up

To start using GenAIRR, you need to import the necessary classes from the module. We'll also set up a `DataConfig` object to specify our configuration.



```python
# Importing GenAIRR classes
from GenAIRR.simulation import HeavyChainSequenceAugmentor, LightChainSequenceAugmentor, SequenceAugmentorArguments
from GenAIRR.utilities import DataConfig
from GenAIRR.data import builtin_heavy_chain_data_config,builtin_kappa_chain_data_config,builtin_lambda_chain_data_config
# Initialize DataConfig with the path to your configuration
#data_config = DataConfig('/path/to/your/config')
# Or Use one of Our Builtin Data Configs
data_config_builtin = builtin_heavy_chain_data_config()


# Set up augmentation arguments (if you have specific requirements)
args = SequenceAugmentorArguments()

```

## Simulating Heavy Chain Sequences

Let's simulate a heavy chain sequence using `HeavyChainSequenceAugmentor`. This example demonstrates a simple simulation with default settings.



```python
# Initialize the HeavyChainSequenceAugmentor
heavy_augmentor = HeavyChainSequenceAugmentor(data_config_builtin, args)

# Simulate a heavy chain sequence
heavy_sequence = heavy_augmentor.simulate_augmented_sequence

# Print the simulated heavy chain sequence
print("Simulated Heavy Chain Sequence:", heavy_sequence)

```

    Simulated Heavy Chain Sequence: <bound method HeavyChainSequenceAugmentor.simulate_augmented_sequence of <GenAIRR.simulation.heavy_chain_sequence_augmentor.HeavyChainSequenceAugmentor object at 0x000001FD56378D90>>
    

## Customizing Simulations

GenAIRR allows for extensive customization to closely mimic the natural diversity of immune sequences. Below is an example of how to customize mutation rates and indel simulations.



```python
# Customize augmentation arguments
custom_args = SequenceAugmentorArguments(min_mutation_rate=0.01, max_mutation_rate=0.05, simulate_indels=True, max_indels=3,
                                         corrupt_proba=0.7,save_ns_record=True,save_mutations_record=True)

# Use custom arguments to simulate a heavy chain sequence
custom_heavy_augmentor = HeavyChainSequenceAugmentor(data_config_builtin, custom_args)
custom_heavy_sequence = custom_heavy_augmentor.simulate_augmented_sequence()

# Print the customized heavy chain sequence
print("Customized Simulated Heavy Chain Sequence:", custom_heavy_sequence)

```

    Customized Simulated Heavy Chain Sequence: {'sequence': 'GTGTTGGAGTACGAACGCGGAGTTCTGTTGTGAATTGGGCGGTGAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTCCAGCCTGGGGGGCCCCTGNGACTCTCCTGTGCAGCCTCTGGANTCACCTTTAGTAGCTATTGGNTGAGGTGNGTCCGCCAGGCTCCAGGGAAGGGACTGGAGTGGGTGGCCAACATAAAACAAGATGGAAGTGAGAAATACTATGTNGACTCTGTGAAGGGCCGATTCACCATCTCCAGAGACAACGCCAAGAACTCACTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACNCGGCNGTGTATTACTGTGCGAGAGTCCGACAGGAGCAGCCAAATCGTCTCTTCGGCTACTCAGGGACCCTTTCTGGTTNGACCCCTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAG', 'v_sequence_start': 43, 'v_sequence_end': 338, 'd_sequence_start': 347, 'd_sequence_end': 353, 'j_sequence_start': 386, 'j_sequence_end': 433, 'v_call': 'IGHVF10-G49*03,IGHVF10-G49*04', 'd_call': 'IGHD6-13*01,IGHD6-25*01,IGHD6-6*01', 'j_call': 'IGHJ5*02', 'mutation_rate': 0.02771362586605081, 'v_trim_5': 0, 'v_trim_3': 1, 'd_trim_5': 6, 'd_trim_3': 6, 'j_trim_5': 4, 'j_trim_3': 0, 'corruption_event': 'add', 'corruption_add_amount': 43, 'corruption_remove_amount': 0, 'mutations': b'ezkxOiAnVD5DJywgMTQ3OiAnQz5HJywgMTc0OiAnRz5BJywgMTk4OiAnRz5BJ30=', 'Ns': b'ezk3OiAnQT5OJywgMTIxOiAnVD5OJywgMTQyOiAnQT5OJywgMTUwOiAnRz5OJywgMjI1OiAnRz5OJywgMzEzOiAnQT5OJywgMzE4OiAnVD5OJywgMzkyOiAnQz5OJ30=', 'indels': {}}
    

## Generating Naïve Sequences

In immunogenetics, a naïve sequence refers to an antibody sequence that has not undergone the process of somatic hypermutation. GenAIRR allows you to simulate such naïve sequences using the `HeavyChainSequence` class. Let's start by generating a naïve heavy chain sequence.



```python
from GenAIRR.sequence import HeavyChainSequence

# Create a naive heavy chain sequence
naive_heavy_sequence = HeavyChainSequence.create_random(data_config_builtin)

# Access the generated naive sequence
naive_sequence = naive_heavy_sequence

print("Naïve Heavy Chain Sequence:", naive_sequence)
print('Ungapped Sequence: ')
print(naive_sequence.ungapped_seq)

```

    Naïve Heavy Chain Sequence: 0|-----------------------------------------------------------------------------V(IGHVF3-G8*01)|294|296|----D(IGHD2-8*02)|312|332|------------J(IGHJ2*01)|381
    Ungapped Sequence: 
    CAGGTGCAGCTGCAGGAGTCGGGCCCAGGACTGGTGAAGCCTCCGGGGACCCTGTCCCTCACCTGCGCTGTCTCTGGTGGCTCCATCAGCAGTAGTAACTGGTGGAGTTGGGTCCGCCAGCCCCCAGGGAAGGGGCTGGAGTGGATTGGGGAAATCTATCATAGTGGGAGCACCAACTACAACCCGTCCCTCAAGAGTCGAGTCACCATATCAGTAGACAAGTCCAAGAACCAGTTCTCCCTGAAGCTGAGCTCTGTGACCGCCGCGGACACGGCCGTGTATTGCTGTGCGAGAAAACTGGTGGTGTATGCTCCTATCTCCCACATAGGGCTTGGTACTTCGATCTCTGGGGCCGTGGCACCCTGGTCACTGTCTCCTCAG
    

## Applying Mutations

To mimic the natural diversity and evolution of immune sequences, GenAIRR supports the simulation of mutations through various models. Here, we demonstrate how to apply mutations to a naïve sequence using the `S5F` and `Uniform` mutation models from the mutations submodule.


### Using the S5F Mutation Model

The `S5F` model is a sophisticated mutation model that considers context-dependent mutation probabilities. It's particularly useful for simulating realistic somatic hypermutations.



```python
from GenAIRR.mutation import S5F

# Initialize the S5F mutation model with custom mutation rates
s5f_model = S5F(min_mutation_rate=0.01, max_mutation_rate=0.05)

# Apply mutations to the naive sequence using the S5F model
s5f_mutated_sequence, mutations, mutation_rate = s5f_model.apply_mutation(naive_heavy_sequence)

print("S5F Mutated Heavy Chain Sequence:", s5f_mutated_sequence)
print("S5F Mutation Details:", mutations)
print("S5F Mutation Rate:", mutation_rate)

```

    S5F Mutated Heavy Chain Sequence: CAGGTGCAGCTGCAGGAGTCGGGCCCAGGACTGGTGAAGCCTCCGGGGACCCTGTCCCTCACCTGCGCTGTCTCTGCTGGCTCCATCAGCAGTAGTAACTGGTGGAGTTGGGTCCGCCAGCCCCCAGGGAAGGGGCTGGAGTGGATTGGGGAAATCTATCATAGTGGGAGCACCAACTACAACCCGTCCCTCTAGAGTCGAGTCACCATATCAGTAGACAAGTCCAAGAACCAGTTCTCCCTGAAGCAGAGCTCTGTGACCGCCGCGGACTCGGCCGTGTATTGCTGTGCGAGAAAACTGGTGGTGTATGCTCCTATCTCCCACATAGGGCTTGGTACTTCGATCTCTGGGGCCGTGGCACCCTGGTCACTGTCTCCTCAG
    S5F Mutation Details: {270: 'A>T', 192: 'A>T', 247: 'T>A', 76: 'G>C'}
    S5F Mutation Rate: 0.011222406361310347
    

### Using the Uniform Mutation Model

The `Uniform` mutation model applies mutations at a uniform rate across the sequence, providing a simpler alternative to the context-dependent models.



```python
from GenAIRR.mutation import Uniform

# Initialize the Uniform mutation model with custom mutation rates
uniform_model = Uniform(min_mutation_rate=0.01, max_mutation_rate=0.05)

# Apply mutations to the naive sequence using the Uniform model
uniform_mutated_sequence, mutations, mutation_rate = uniform_model.apply_mutation(naive_heavy_sequence)

print("Uniform Mutated Heavy Chain Sequence:", uniform_mutated_sequence)
print("Uniform Mutation Details:", mutations)
print("Uniform Mutation Rate:", mutation_rate)

```

    Uniform Mutated Heavy Chain Sequence: CAGGTGCACCTGCAGGAGTCGGGCCGAGGAGTGGTGAAGCCTCCGGGGACCCTGTCCCTCACCTGCGCTGTCTCTGGTGGCTCCATCAGCAGTAGTTACTTGTGGAGTTGGGTCCGCCAGCCACCAGGGAAGGGGCTGGAGTGGATTGGGGAAATCTATCATAGTGGGAGCACCAACTACAACCCGTCCCTCAAGAGTCGAGTCACCATATCAGTAGACAAGTCCAAGAACCAGTTCTCCCTGAAGCTGAGCTCTGTGACCGCCGCGGACACGGCCGTGTATTGCTGTGCGAGAAAACTGGTGGTGTATGCTCCTATCTCCCACATAGGGCTTGGTACTTCGATCTGTGGGGCCGTGGCACCCTGGTCACTGTCTCCTCAG
    Uniform Mutation Details: {122: 'C>A', 8: 'G>C', 100: 'G>T', 96: 'A>T', 25: 'C>G', 346: 'C>G', 30: 'C>G'}
    Uniform Mutation Rate: 0.019802269687583134
    

## Common Use Cases

GenAIRR is a versatile tool designed to meet a broad range of needs in immunogenetics research. This section provides examples and explanations for some common use cases, including generating multiple sequences, simulating specific allele combinations, and more.


### Generating Many Sequences

One common requirement is to generate a large dataset of synthetic AIRR sequences for analysis or benchmarking. Below is an example of how to generate multiple sequences using GenAIRR in a loop.



```python
num_sequences = 5  # Number of sequences to generate

heavy_sequences = []
for _ in range(num_sequences):
    # Simulate a heavy chain sequence
    heavy_sequence = heavy_augmentor.simulate_augmented_sequence()
    heavy_sequences.append(heavy_sequence)

# Display the generated sequences
for i, seq in enumerate(heavy_sequences, start=1):
    print(f"Heavy Chain Sequence {i}: {seq}")

```

    Heavy Chain Sequence 1: {'sequence': 'TTGGNAAGCCAGGCCCTGGAGTGACTTTCACACACTGATCGGTGCGANGCTGAACTCCACAACGCCTCCCTCAAAACCTGACCCACCACGTCCAGGGACCCGTCCGGTAGTCACATGGTCCTGACACTGTCGAACATGGACCCTGTGGACACAGTCACACATTACTGTGCACCGATNCCCCCCCCTACGANGATTCCGGCCGGGCCCTGGCTAATCCAATCACTTGTTGGAGGTCTGGGGCAAAGGGACCACGGCCACCGACTCNTAAG', 'v_sequence_start': 9, 'v_sequence_end': 176, 'd_sequence_start': 185, 'd_sequence_end': 199, 'j_sequence_start': 207, 'j_sequence_end': 269, 'v_call': 'IGHVF1-G3*06,IGHVF1-G3*05,IGHVF1-G3*04', 'd_call': 'IGHD3-10*03', 'j_call': 'IGHJ6*03', 'mutation_rate': 0.241635687732342, 'v_trim_5': 0, 'v_trim_3': 2, 'd_trim_5': 4, 'd_trim_3': 13, 'j_trim_5': 2, 'j_trim_3': 0, 'corruption_event': 'remove_before_add', 'corruption_add_amount': 9, 'corruption_remove_amount': 132, 'indels': {}}
    Heavy Chain Sequence 2: {'sequence': 'CAGCTGCAGTTGCAGGAGTCGGGCCCNGGACTGGTGAAGCCTTTGGAGGCCCAGTCCCTCTCGTACCCTGTCTCTGGTGACTCCATCAGCAATAGTGGTTACTCCTGGGGCTGAATCCGTCCCCCCNCAGGGAAGGGGCTGGAGTGGATNGCGACTATANATTATAGGGGCAGCTCCTGCTACAACCCGTCCCTCAAGAGTCGAGTCACCATCTCCACAGACACGTCCAAGAAGCAGGTCTCCCTGATGCTGAGCTCTATGACCGCCGCANACACGACTGTNTATTACTGTGCGAGAGTCATGGTTCTGATGTTTTGGAGCAACTGGTTCGACCCCTGGGACCAGGGAAGCCTGGTCACCCTCTCCTCAN', 'v_sequence_start': 0, 'v_sequence_end': 297, 'd_sequence_start': 308, 'd_sequence_end': 317, 'j_sequence_start': 320, 'j_sequence_end': 370, 'v_call': 'IGHVF3-G10*06', 'd_call': 'IGHD3-9*01', 'j_call': 'IGHJ5*02', 'mutation_rate': 0.11621621621621622, 'v_trim_5': 0, 'v_trim_3': 2, 'd_trim_5': 8, 'd_trim_3': 15, 'j_trim_5': 1, 'j_trim_3': 0, 'corruption_event': 'no-corruption', 'corruption_add_amount': 0, 'corruption_remove_amount': 0, 'indels': {}}
    Heavy Chain Sequence 3: {'sequence': 'CAGAAGAGACTGGTGCAGTCTGGGGTTGACATGAAGACGACTGGGTCGTAATTGAAACTTTCACGAAAGACTTCTGAATACACTCGCACANACCGCTATCTGCACTGGGTCCGACAGGCCCCCAGACGGGCGTTTGAGTGGGTGGGGNGGATCACGCCTTTCAGTGGTAACACCCACTACGTGCAGACGTCCCAGGACAGAGTCCCCATTACCAGGNACAAGTNTACGAGTCCAGCCTATATAGAACTGAACACCCTNAAATGCGAGGACACAGACATATATTAATGCGCANGATCCACGGGAACCCCAGCNGAGAACTGGTACTTCGATCTTTGGGGCCGTGGCCCCCTGATCACCGTCTACTCTG', 'v_sequence_start': 0, 'v_sequence_end': 295, 'd_sequence_start': 295, 'd_sequence_end': 305, 'j_sequence_start': 316, 'j_sequence_end': 367, 'v_call': 'IGHVF6-G20*02', 'd_call': 'IGHD4-11*01,IGHD4-4*01', 'j_call': 'IGHJ2*01', 'mutation_rate': 0.1989100817438692, 'v_trim_5': 0, 'v_trim_3': 1, 'd_trim_5': 3, 'd_trim_3': 3, 'j_trim_5': 3, 'j_trim_3': 0, 'corruption_event': 'no-corruption', 'corruption_add_amount': 0, 'corruption_remove_amount': 0, 'indels': {}}
    Heavy Chain Sequence 4: {'sequence': 'CAGTTTCAGCTGGTGCCGTCTGGAGCTGAGGTGAAGAAGNCTGNGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCTGANTACACCTTCACCAGGTATGATATNAGCTNGGTTCGACAGGCCCCTGGACAAGGGCTTGAGTGGGTGGGATGGATCAGCGCTTACAAGGGTAACACAAACTATGAACAGAAGCTCCAGGGCAGAGTCACCATGACCACTGACACATCCACGAGCACAGCCTACATAGAGCTGAGGAGTCTGAGATCTGACGACACGGCCGTGTATCACTGTGCGAGAATCGGCGGCAGGGACGAGTCCGCAGATATCTCGCATCCCTATTGCTACTCCGGTATGGACGTCTGGGGCCAAGNNACCACGGTCACCGTCTCCTCAG', 'v_sequence_start': 0, 'v_sequence_end': 294, 'd_sequence_start': 316, 'd_sequence_end': 323, 'j_sequence_start': 332, 'j_sequence_end': 391, 'v_call': 'IGHVF6-G25*02', 'd_call': 'IGHD5-18*01,IGHD5-5*01', 'j_call': 'IGHJ6*02', 'mutation_rate': 0.0639386189258312, 'v_trim_5': 0, 'v_trim_3': 2, 'd_trim_5': 8, 'd_trim_3': 6, 'j_trim_5': 4, 'j_trim_3': 0, 'corruption_event': 'no-corruption', 'corruption_add_amount': 0, 'corruption_remove_amount': 0, 'indels': {}}
    Heavy Chain Sequence 5: {'sequence': 'GAGGTGCAACTGCTGCAGACCGGGTCAGACTTGATACAGCCAGGGGAGTCCCTCANACTGCCCTGTGCAGCCTCTGGATTCACCTGGTGAANGNATGCCGTGAATTGGGGCCGGCGGCCTCCAGGGATGGGACTTGATTGGGTCTCAGTTCTNAGTGCTAGTGGTGAGAGAACNTTCTCCATAGACTCCATGAAGGGCCGGGTCACCACCTCCAGGGTCAATTGCAAGAGTACGCTGTATCTGAAAATGAAGGGCCTGAGAGCCGAGGACGCGGCTGTTTATTATTGAGCGAGAGAGGCCTTAGGGTCGGATTACTACTCCTTTTACATGGACGTCTGGGGCACAGGGACCGCGGNCACCGTCTCGTCAC', 'v_sequence_start': 0, 'v_sequence_end': 296, 'd_sequence_start': 302, 'd_sequence_end': 306, 'j_sequence_start': 310, 'j_sequence_end': 370, 'v_call': 'IGHVF10-G41*02', 'd_call': 'Short-D', 'j_call': 'IGHJ6*03', 'mutation_rate': 0.1972972972972973, 'v_trim_5': 0, 'v_trim_3': 0, 'd_trim_5': 7, 'd_trim_3': 6, 'j_trim_5': 4, 'j_trim_3': 0, 'corruption_event': 'no-corruption', 'corruption_add_amount': 0, 'corruption_remove_amount': 0, 'indels': {}}
    


```python
import pandas as pd
pd.DataFrame(heavy_sequences)
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>sequence</th>
      <th>v_sequence_start</th>
      <th>v_sequence_end</th>
      <th>d_sequence_start</th>
      <th>d_sequence_end</th>
      <th>j_sequence_start</th>
      <th>j_sequence_end</th>
      <th>v_call</th>
      <th>d_call</th>
      <th>j_call</th>
      <th>...</th>
      <th>v_trim_5</th>
      <th>v_trim_3</th>
      <th>d_trim_5</th>
      <th>d_trim_3</th>
      <th>j_trim_5</th>
      <th>j_trim_3</th>
      <th>corruption_event</th>
      <th>corruption_add_amount</th>
      <th>corruption_remove_amount</th>
      <th>indels</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>TTGGNAAGCCAGGCCCTGGAGTGACTTTCACACACTGATCGGTGCG...</td>
      <td>9</td>
      <td>176</td>
      <td>185</td>
      <td>199</td>
      <td>207</td>
      <td>269</td>
      <td>IGHVF1-G3*06,IGHVF1-G3*05,IGHVF1-G3*04</td>
      <td>IGHD3-10*03</td>
      <td>IGHJ6*03</td>
      <td>...</td>
      <td>0</td>
      <td>2</td>
      <td>4</td>
      <td>13</td>
      <td>2</td>
      <td>0</td>
      <td>remove_before_add</td>
      <td>9</td>
      <td>132</td>
      <td>{}</td>
    </tr>
    <tr>
      <th>1</th>
      <td>CAGCTGCAGTTGCAGGAGTCGGGCCCNGGACTGGTGAAGCCTTTGG...</td>
      <td>0</td>
      <td>297</td>
      <td>308</td>
      <td>317</td>
      <td>320</td>
      <td>370</td>
      <td>IGHVF3-G10*06</td>
      <td>IGHD3-9*01</td>
      <td>IGHJ5*02</td>
      <td>...</td>
      <td>0</td>
      <td>2</td>
      <td>8</td>
      <td>15</td>
      <td>1</td>
      <td>0</td>
      <td>no-corruption</td>
      <td>0</td>
      <td>0</td>
      <td>{}</td>
    </tr>
    <tr>
      <th>2</th>
      <td>CAGAAGAGACTGGTGCAGTCTGGGGTTGACATGAAGACGACTGGGT...</td>
      <td>0</td>
      <td>295</td>
      <td>295</td>
      <td>305</td>
      <td>316</td>
      <td>367</td>
      <td>IGHVF6-G20*02</td>
      <td>IGHD4-11*01,IGHD4-4*01</td>
      <td>IGHJ2*01</td>
      <td>...</td>
      <td>0</td>
      <td>1</td>
      <td>3</td>
      <td>3</td>
      <td>3</td>
      <td>0</td>
      <td>no-corruption</td>
      <td>0</td>
      <td>0</td>
      <td>{}</td>
    </tr>
    <tr>
      <th>3</th>
      <td>CAGTTTCAGCTGGTGCCGTCTGGAGCTGAGGTGAAGAAGNCTGNGG...</td>
      <td>0</td>
      <td>294</td>
      <td>316</td>
      <td>323</td>
      <td>332</td>
      <td>391</td>
      <td>IGHVF6-G25*02</td>
      <td>IGHD5-18*01,IGHD5-5*01</td>
      <td>IGHJ6*02</td>
      <td>...</td>
      <td>0</td>
      <td>2</td>
      <td>8</td>
      <td>6</td>
      <td>4</td>
      <td>0</td>
      <td>no-corruption</td>
      <td>0</td>
      <td>0</td>
      <td>{}</td>
    </tr>
    <tr>
      <th>4</th>
      <td>GAGGTGCAACTGCTGCAGACCGGGTCAGACTTGATACAGCCAGGGG...</td>
      <td>0</td>
      <td>296</td>
      <td>302</td>
      <td>306</td>
      <td>310</td>
      <td>370</td>
      <td>IGHVF10-G41*02</td>
      <td>Short-D</td>
      <td>IGHJ6*03</td>
      <td>...</td>
      <td>0</td>
      <td>0</td>
      <td>7</td>
      <td>6</td>
      <td>4</td>
      <td>0</td>
      <td>no-corruption</td>
      <td>0</td>
      <td>0</td>
      <td>{}</td>
    </tr>
  </tbody>
</table>
<p>5 rows × 21 columns</p>
</div>



### Generating a Specific Allele Combination Sequence

In some cases, you might want to simulate sequences with specific V, D, and J allele combinations. Here's how to specify alleles for your simulations.



```python
# Define your specific alleles
v_allele = 'IGHVF6-G21*01'
d_allele = 'IGHD5-18*01'
j_allele = 'IGHJ6*03'

# Extract the allele objects from data_config
v_allele = next((allele for family in data_config_builtin.v_alleles.values() for allele in family if allele.name == v_allele), None)
d_allele = next((allele for family in data_config_builtin.d_alleles.values() for allele in family if allele.name == d_allele), None)
j_allele = next((allele for family in data_config_builtin.j_alleles.values() for allele in family if allele.name == j_allele), None)

# Check if all alleles were found
if not v_allele or not d_allele or not j_allele:
    raise ValueError("One or more specified alleles could not be found in the data config.")


# Generate a sequence with the specified allele combination
specific_allele_sequence = HeavyChainSequence([v_allele, d_allele, j_allele], data_config_builtin)
specific_allele_sequence.mutate(s5f_model)



print("Specific Allele Combination Sequence:", specific_allele_sequence.mutated_seq)

```

    Specific Allele Combination Sequence: CAGGTGCAGTTGGTGCAGTCTGGGACTGAGTTGAAGACGCCTGGGTCCTCGGTGAAGGTCTCCTGCAAGGCTTCTAGAGGCACCTTCAGCAGCTCTGCTATCAGCTGGGTGCGACAGGCCCCTGGACAAGGGCTTGAGTGGATGGGAGGGATCATCCCTATCTTTGGTACAGCAAACTACGCACAGAAGTTCCAGGGCAGAGTCACGATTACCGCGGATAAATCCACGAGCACAGCCTACATGGAGCTGAGCAGCCTGAGATCTGAGGATACGGCCGTGTATTACTGTGCGAGAGAGGATGGGTCCGGATCCCACCCCATTTACTATTACTACTACTACATGGACGTCTGGGGCAAAGGGACCACGGTCACCGTCTCCTCAG
    

### Simulating Sequences with Custom Mutation Rates

Adjusting mutation rates allows for the simulation of sequences at various stages of affinity maturation. Here's how to customize mutation rates in your simulations.



```python
# Customize augmentation arguments with your desired mutation rates
custom_args = SequenceAugmentorArguments(min_mutation_rate=0.15, max_mutation_rate=0.3)

# Initialize the augmentor with custom arguments
custom_augmentor = HeavyChainSequenceAugmentor(data_config_builtin, custom_args)

# Generate a sequence with the custom mutation rates
custom_mutation_sequence = custom_augmentor.simulate_augmented_sequence()

print("Custom Mutation Rate Sequence:", custom_mutation_sequence)

```

    Custom Mutation Rate Sequence: {'sequence': 'GGTTGGAGCTCATTGGGAGCTNCTATTCTAGTGGGACTACCTAGTACAACCTGTCCCTCAAGAATCGCGTCACCATATCAGTCGACACGTCCAAGAATCANTCCTCCCTGGAGCTGAGCTCCGTGACCGCAGCGGACACGGCCGTGCCTNGTTGNGCGGGAAAGTTGAATATAGTGGCTAACTCTGCCTTTTGCTCTCTGGGGCCAGGGGACAGTGGCCACTGTTTTTTCAG', 'v_sequence_start': 0, 'v_sequence_end': 161, 'd_sequence_start': 165, 'd_sequence_end': 180, 'j_sequence_start': 186, 'j_sequence_end': 232, 'v_call': 'IGHVF3-G10*04', 'd_call': 'IGHD5-12*01,IGHD5-18*02', 'j_call': 'IGHJ3*02', 'mutation_rate': 0.15517241379310345, 'v_trim_5': 0, 'v_trim_3': 1, 'd_trim_5': 4, 'd_trim_3': 8, 'j_trim_5': 4, 'j_trim_3': 0, 'corruption_event': 'remove', 'corruption_add_amount': 0, 'corruption_remove_amount': 136, 'indels': {}}
    

### Generating Naïve vs. Mutated Sequence Pairs

Comparing naïve and mutated versions of the same sequence can be useful for studying somatic hypermutation effects. Here's how to generate such pairs with GenAIRR.



```python
# Generate a naive sequence
sequence_object = HeavyChainSequence.create_random(data_config_builtin)
sequence_object.mutate(s5f_model)

print("Naïve Sequence:", sequence_object.ungapped_seq)
print("Mutated Sequence:", sequence_object.mutated_seq)

```

    Naïve Sequence: CAGGTGCAGCTGCAGGAGTCGGGCCCAGGACTGGTGAAGCCTTCGGACACCCTGTCCCTCACCTGCGCTGTCTCTGGTTACTCCATCAGCAGTAGTAACTGGTGGGGCTGGATCCGGCAGCCCCCAGGGAAGGGACTGGAGTGGATTGGGTACATCTATTATAGTGGGAGCATCTACTACAACCCGTCCCTCAAGAGTCGAGTCACCATGTCAGTAGACACGTCCAAGAACCAGTTCTCCCTGAAGCTGAGCTCTGTGACCGCCGTGGACACGGCCGTGTATTACTGTGCGAGAAAGCCACTCGGTCACACTACGGTGGTAACTCATGCTTTTGATATCTGGGGCCAAGGGACAATGGTCACCGTCTCTTCAG
    Mutated Sequence: CAGGTGCAGCTGCAGGAGTCGGGCCCAGGACTGGTGAAGCCTTCGGACACCCTGTCCCTCACCTGCGCTGTCTCTGGTTACTCCATCAGCAGTAGTAACTGGTGGGGCTGGATCCGGCAGCCCCCAGGGAAGGGACTGGAGTGGATTGGGAACATCCATTATAGTGGGAGCATCTACTACAACCCGTCCCTCAAGAGTCGAGTCACCATGTCACTAGACACGTCCAAGAACCAGTTCTCCCTGAAACTGAGCTCTGTGGCCGCCGTGGACACGGCCGTGTATTACTGTGCGAGAACGCCACTCGGTCACACTACGGTGGTAATTCATGCTTTTGATATCTGGGGCCAAGGGACAATGGTCACCGTCTCTTCAG
    



## Conclusion

This section highlighted some common use cases for GenAIRR, demonstrating its flexibility in simulating AIRR sequences for various research purposes. Whether you need large datasets, specific allele combinations, custom mutation rates, or comparative analyses of naïve and mutated sequences, GenAIRR provides the necessary tools to achieve your objectives.

