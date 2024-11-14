# Creating a Custom `DataConfig` Instance

GenAIRR support for various species is till under development, specific use cases or unique data requirements may necessitate creating a custom `DataConfig` instance. This guide will help you understand how to build a custom `DataConfig` instance to simulate sequences tailored to your reference files and empirical distributions.

Before diving into the calculations for each parameter, let’s first examine the structure of the `DataConfig` class and the attributes it contains:

```python
def __init__(self):
    # Config Variables
    self.trim_dicts = {}
    self.NP_transitions = {}
    self.NP_first_bases = {}
    self.NP_lengths = {}
    self.kmer_dicts = {}
    self.v_alleles = None
    self.d_alleles = None
    self.j_alleles = None
    self.c_alleles = None
    self.name = None

    self.correction_maps = dict()
    self.asc_tables = dict()
```

This constructor outlines the essential attributes in a `DataConfig` instance, which you will need to populate with your specific data in the format discussed below.

**NOTE** - After elaborating on all attributes there is a guide on how to use an automation
we wrote for your convenience that will manage the creation of all these parts for you.
Scroll down to `Automated DataConfig Creation` for details.

---

## Attribute Breakdown

### `trim_dicts`
- **Description**: A nested dictionary holding the probability distributions for different trimming lengths of alleles during V(D)J recombination.
- **Structure**:
  - The top level should have keys representing allele types and the side being trimmed: `['V_3', 'D_5', 'D_3', 'J_5', 'C_3']`. For example, `V_3` represents the 3' end trimming of the V allele.
  - The second level should have keys for each allele family or ASC (Allele Cluster), such as `'IGHVF1'`, `'IGHVF2'`, etc., for a BCR heavy chain.
  - The third level should include specific gene names under each family/ASC, such as `'IGHVF10-G33'`, `'IGHVF10-G34'`, etc.

**Example**: Accessing `trim_dicts` in your custom `DataConfig` instance:

```python
dataconfig.trim_dicts['V_3']['some_family/asc']['some_gene_of_the_same_family']
```

This should return a dictionary mapping trimming lengths to their likelihoods:

```python
defaultdict(float, {
    0: 0.47,
    1: 0.26,
    2: 0.095,
    3: 0.1,
    4: 0.026,
    5: 0.017,
    6: 0.0091,
    7: 0.004,
    8: 0.004,
    9: 0.0015,
    ...
})
```

#### Explanation of Each Level
1. **Top Level**: Specifies which part of the allele is being trimmed (e.g., `V_3` for the 3' end of V alleles).
2. **Second Level**: Contains keys for each allele family or ASC, grouping genes under a common classification.
3. **Third Level**: Lists the specific genes within each family or ASC. The values are dictionaries mapping trimming lengths (keys) to their likelihoods (values).

#### How to Populate `trim_dicts`
- **Data Collection**: Derive the likelihood of different trimming lengths from empirical data or computational predictions.
- **Data Formatting**: Ensure the values are stored as dictionaries with keys ranging from 0 to the maximum trimming length and values as the probabilities.
- **Valid Likelihood Distribution**: Ensure the values are stored at the final dictionary sum to 1.

---


### `NP_transitions`

The `NP_transitions` attribute defines a first-order Markov chain that guides GenAIRR in simulating the NP region segments of sequences. It is structured as a nested dictionary with transition probabilities to model nucleotide occurrences at different positions in the NP regions.

**Key Structure**:
- The `NP_transitions` dictionary must have two main keys: `'NP1'` and `'NP2'`. If your simulation does not include a D allele (e.g., for VJ recombination), only the `'NP1'` key is relevant, as it represents the VJ NP region. For simulations involving D alleles (VDJ recombination), both `'NP1'` and `'NP2'` should be populated.
  
**Nested Key Hierarchy**:
1. **First Level**: Position index in the NP region (`0`, `1`, `2`, etc.).
2. **Second Level**: The nucleotide present at that position (`'A'`, `'T'`, `'C'`, `'G'`).
3. **Third Level**: A dictionary mapping each nucleotide (`'A'`, `'T'`, `'C'`, `'G'`) to the likelihood of it appearing as the next nucleotide in the sequence.

**Example**:
Accessing transition probabilities in your custom `DataConfig` instance:

```python
dataconfig.NP_transitions['NP1'][0]['T']
```

This might return:

```python
{'A': 0.169,
 'C': 0.390,
 'G': 0.233,
 'T': 0.208}
```

**Interpretation**:
This example can be understood as follows: when simulating the NP1 region, if the current nucleotide at position `0` is `'T'`, the dictionary returned provides the likelihoods for the next nucleotide being `'A'`, `'C'`, `'G'`, or `'T'`.

**Guidelines**:
- Ensure that the probabilities for each set of transitions sum to 1 to maintain consistency in the Markov chain.
- Populate this structure based on empirical data or estimated probabilities to simulate realistic nucleotide transitions.

---

### `NP_first_bases`

The `NP_first_bases` attribute is a component of the NP region simulation data, guiding the selection of the starting nucleotide in the NP regions. It is structured similarly to the `NP_transitions` attribute and provides the initial probabilities for the first nucleotide at the beginning of each NP region.

**Key Structure**:
- `NP_first_bases` is a dictionary with two main keys: `'NP1'` and `'NP2'`. If your simulation does not involve a D allele (e.g., VJ recombination), only the `'NP1'` key is relevant. For VDJ recombination simulations, populate both `'NP1'` and `'NP2'`.

**Nested Key Hierarchy**:
- Each main key (`'NP1'`, `'NP2'`) maps to a dictionary that holds the probabilities for the first nucleotide in that NP region. The sub-dictionary should include the nucleotides `'A'`, `'T'`, `'C'`, and `'G'` as keys, each associated with the likelihood of it being the starting base.

**Example**:
Accessing `NP_first_bases` in your custom `DataConfig` instance:

```python
dataconfig.NP_first_bases
```

This might return:

```python
{
  'NP1': {'A': 0.1117,
          'C': 0.2461,
          'G': 0.2846,
          'T': 0.3576},
  'NP2': {'A': 0.1673,
          'C': 0.3713,
          'G': 0.2751,
          'T': 0.1864}
}
```

**Interpretation**:
When simulating the NP regions, GenAIRR will refer to this attribute to sample the first nucleotide of the NP region. For example, if the simulation starts generating the NP1 region, it will use the probabilities under `'NP1'` to pick the initial nucleotide.

**Guidelines**:
- Ensure that the sum of probabilities for each sub-dictionary (e.g., for `'NP1'` or `'NP2'`) is equal to 1 to maintain consistency.
- Populate these probabilities based on empirical data or estimated distributions to reflect realistic biological nucleotide distributions.

---

### `NP_lengths`

The `NP_lengths` attribute is the final component related to NP region simulation in GenAIRR. It provides the probability distribution for the lengths of the NP regions, guiding how long each region will be during sequence simulation.

**Key Structure**:
- Like the other NP region attributes, `NP_lengths` is a dictionary with two main keys: `'NP1'` and `'NP2'`. Populate both keys for VDJ recombination, but only `'NP1'` is needed for VJ recombination.
- Each main key maps to a sub-dictionary representing the distribution of potential NP region lengths.

**Nested Key Hierarchy**:
- Under each main key (`'NP1'`, `'NP2'`), the sub-dictionary contains integer keys representing possible NP region lengths (e.g., 0, 1, 2, etc.). Each key maps to a value that indicates the likelihood of the NP region being that length.

**Example**:
Accessing `NP_lengths` in your custom `DataConfig` instance:

```python
dataconfig.NP_lengths['NP1']
```

This might return:

```python
{
  0: 0.052,
  1: 0.040,
  2: 0.056,
  3: 0.070,
  4: 0.078,
  5: 0.079,
  6: 0.077,
  7: 0.074,
  ...
}
```

**Interpretation**:
When generating an NP region, GenAIRR will use this distribution to sample the length of the NP region. For example, if simulating the NP1 region, it will use the probabilities under `'NP1'` to determine its length.

**Guidelines**:
- Ensure the length distribution starts with a key of 0, representing the possibility of an NP region length of 0.
- The values (likelihoods) in each sub-dictionary should sum to 1 to maintain a valid probability distribution.

**Usage Note**:
Properly defining the `NP_lengths` attribute helps maintain realistic simulations by ensuring that the NP regions are generated with appropriate lengths that match empirical or expected biological distributions.

---

### `v_alleles`, `d_alleles`, `j_alleles`, `c_alleles`

In GenAIRR, each allele from the reference set is encapsulated within specific allele classes. The base class `Allele` is extended by `VAllele`, `DAllele`, `JAllele`, and `CAllele` classes, each designed to hold relevant properties and methods for V, D, J, and constant (C) alleles, respectively.

**Purpose**:
The `v_alleles`, `d_alleles`, `j_alleles`, and `c_alleles` attributes in the `DataConfig` instance are dictionaries representing these alleles at the gene level. They help guide the sequence simulation process by defining the allele repertoire used in recombination.

**Structure**:
- Each attribute (`v_alleles`, `d_alleles`, `j_alleles`, `c_alleles`) is a dictionary where:
  - **Key**: The gene name (e.g., `'IGHVF1-G1'`).
  - **Value**: A list of allele objects for that gene in the relevant `Allele` format (e.g., `VAllele` for `v_alleles`).

**Example**:
Accessing `dataconfig.v_alleles` might look like:

```python
dataconfig.v_alleles['IGHVF1-G1']
```

Which would return:

```python
[
    VAllele(name='IGHVF1-G1*01', sequence='...gcatggatac', length=301, ungapped_len=301, family='IGHVF1', gene='IGHVF1-G1'),
    VAllele(name='IGHVF1-G1*02', sequence='...gcacggatac', length=301, ungapped_len=301, family='IGHVF1', gene='IGHVF1-G1'),
    VAllele(name='IGHVF1-G1*03', sequence='...gcacggatac', length=301, ungapped_len=301, family='IGHVF1', gene='IGHVF1-G1')
]
```

**Explanation**:
- Each key represents the gene-level name, derived from the reference set.
- Each list under a key contains instances of the corresponding allele class (`VAllele`, `DAllele`, `JAllele`, `CAllele`), encapsulating details like the full sequence, name, length, ungapped length, family, and gene.

**How to Populate**:
Ensure that when creating a custom `DataConfig`, you parse the reference file and map each allele to its appropriate `Allele` subclass before storing it in the respective attribute.

**Benefits**:
Having this structured format allows GenAIRR to access and utilize the full range of alleles for simulating sequences, while also maintaining specific allele-level details necessary for accurate simulations and analyses.


---

### `name`
A human-readable identifier for the DataConfig.

---

Here's the refined version for the `correction_maps` attribute:

---

### `correction_maps`

The `correction_maps` attribute in `DataConfig` holds various structures used by GenAIRR to resolve ambiguities in sequence and allele identification. These maps help identify which alleles become indistinguishable under certain modifications, such as trimming or the insertion of "N" bases. Each `DataConfig` instance should have the following keys in the `correction_maps` attribute:

- `'V_N_AMBIGUITY_CORRECTION_GRAPH'`
- `'V_3_TRIM_SIMILARITY_MAP'`
- `'V_5_TRIM_SIMILARITY_MAP'`
- `'J_3_TRIM_SIMILARITY_MAP'`
- `'J_5_TRIM_SIMILARITY_MAP'`
- `'D_5_3_TRIM_SIMILARITY_MAP'`

Each key maps to a data structure that helps in adjusting ground truth information based on sequence modifications. Here's an overview of how each map is used:

#### Trim Similarity Maps (`V_3_TRIM_SIMILARITY_MAP`, `V_5_TRIM_SIMILARITY_MAP`, `J_3_TRIM_SIMILARITY_MAP`, `J_5_TRIM_SIMILARITY_MAP`)

These maps help determine which alleles are indistinguishable after a specific amount of trimming. For example, `V_3_TRIM_SIMILARITY_MAP` holds data for the V alleles trimmed at the 3' end. If an allele is trimmed by a certain number of bases, the map lists other alleles that become indistinguishable due to this modification.

**Structure**:
- The first-level key is the name of an allele from the reference set.
- The second-level key is the number of bases trimmed from the specified end (5' or 3').

**Example**:
```python
dataconfig.correction_maps['V_5_TRIM_SIMILARITY_MAP']['IGHVF1-G1*01'][10]
```
This would return a list of alleles that are indistinguishable from `IGHVF1-G1*01` when it has been trimmed by 10 bases at the 5' end.

#### D Allele Trim Similarity Map (`D_5_3_TRIM_SIMILARITY_MAP`)

This map is slightly different due to the unique nature of D alleles. The key for D alleles is a tuple representing both the 5' and 3' trimming amounts.

**Example**:
```python
dataconfig.correction_maps['D_5_3_TRIM_SIMILARITY_MAP']['IGHD1-1*01'][(5, 10)]
```
This would return a list of alleles that are indistinguishable from `IGHD1-1*01` when it has been trimmed by 5 bases at the 5' end and 10 bases at the 3' end.

#### `V_N_AMBIGUITY_CORRECTION_GRAPH`

This graph structure is used to handle ambiguities that arise when "N" bases are inserted at critical positions in the sequence. GenAIRR can identify crucial positions within alleles that are essential for distinguishing between them. If "N" bases replace these positions, the affected alleles become indistinguishable, and the correction graph updates the ground truth list to include these ambiguous alleles.

**Example Scenario**:
Suppose two V alleles differ only at positions 150, 160, and 170. If "N" bases are inserted at these positions during simulation, distinguishing between these alleles becomes impossible. The `V_N_AMBIGUITY_CORRECTION_GRAPH` detects such scenarios and ensures the simulated data reflects this ambiguity by updating the list of V alleles in the container.

**Usage**:
These correction maps ensure that simulated sequences accurately reflect potential ambiguities that may occur during biological processes or sequence modifications, allowing for more realistic and comprehensive simulation outcomes.

---

## Automated DataConfig Creation

Manually calculating all of the above attributes can be a complex and time-consuming process. While there may be situations where manual configuration is unavoidable, GenAIRR offers a convenient solution to automate most of this work: the `CustomDataConfigGenerator` class. This specialized class calculates all required attributes for a `DataConfig` instance based on provided reference FASTA files and a dataset.

### Example Usage

To create a custom `DataConfig` instance, follow this example:

```python
from GenAIRR.generateDataConfig import CustomDataConfigGenerator

# Initialize the DataConfig generator
dcg = CustomDataConfigGenerator(convert_to_asc=False)

# Generate the DataConfig instance using reference files and a custom data CSV
custom_dataconfig = dcg.make_dataconfig_from_existing_reference_files(
    v_reference_path='./IGHV.fasta',
    d_reference_path='./IGHD.fasta',
    j_reference_path='./IGHJ.fasta',
    c_reference_path='./IGHC.fasta',
    custom_data='./inference_sample.csv'
)
```

This will automatically derive all necessary attributes for a `DataConfig` instance. Note that the `custom_data` CSV file should contain the following columns:

```python
['sequence', 'v_sequence_start', 'v_sequence_end', 'd_sequence_start',
 'd_sequence_end', 'j_sequence_start', 'j_sequence_end', 'v_call',
 'd_call', 'j_call', 'mutation_rate', 'v_trim_5', 'v_trim_3', 'd_trim_5',
 'd_trim_3', 'j_trim_5', 'j_trim_3', 'corruption_event',
 'corruption_add_amount', 'corruption_remove_amount']
```

Ensure that the `v_call`, `d_call`, and `j_call` entries use names that match those in the provided V, D, and J reference FASTA files.

### Saving and Reusing Your DataConfig

Once the `DataConfig` instance is created (with detailed logs provided during the process), you can save it as a pickle file for future use. This allows you to easily reload the `DataConfig` without having to generate it again.

**Saving the generated DataConfig:**
```python
import pickle
with open('your/save/path/my_DataConfig.pkl', 'wb') as file:
    pickle.dump(custom_dataconfig, file)
```

**Loading a pre-calculated DataConfig:**
```python
import pickle
with open('your/save/path/my_DataConfig.pkl', 'rb') as file:
    custom_dataconfig = pickle.load(file)
```

By using the automated `CustomDataConfigGenerator`, you can save significant time and ensure that your custom `DataConfig` instance is consistent and accurate, making it easy to integrate with GenAIRR’s full range of functionalities.