# Migration Guide: Updating to Current GenAIRR Syntax

This guide helps users migrate from older GenAIRR syntax to the current API. The changes improve usability and consistency across the library.

## Key Changes Summary

### 1. Data Configuration Loading

**Old Syntax:**
```python
from GenAIRR.data import builtin_heavy_chain_data_config
from GenAIRR.parameters import ChainType, CHAIN_TYPE_INFO

data_cfg = builtin_heavy_chain_data_config()
AugmentationStep.set_dataconfig(config=data_cfg, chain_type=ChainType.BCR_HEAVY)
```

**New Syntax:**
```python
from GenAIRR.data import HUMAN_IGH_OGRDB

AugmentationStep.set_dataconfig(HUMAN_IGH_OGRDB)
```

### 2. SimulateSequence Step Parameters

**Old Syntax:**
```python
SimulateSequence(mutation_model=S5F(min_mutation_rate=0.003, max_mutation_rate=0.25), productive=True)
```

**New Syntax:**
```python
SimulateSequence(S5F(min_mutation_rate=0.003, max_mutation_rate=0.25), True)
```

### 3. Pipeline Step Parameters

**Old Syntax:**
```python
CorruptSequenceBeginning(
    corruption_probability=0.7,
    corrupt_events_proba=[0.4, 0.4, 0.2],
    max_sequence_length=576,
    nucleotide_add_coefficient=210,
    nucleotide_remove_coefficient=310,
    nucleotide_add_after_remove_coefficient=50,
    random_sequence_add_proba=1
)
InsertNs(n_ratio=0.02, proba=0.5)
ShortDValidation(short_d_length=5)
InsertIndels(indel_probability=0.5, max_indels=5, insertion_proba=0.5, deletion_proba=0.5)
```

**New Syntax:**
```python
CorruptSequenceBeginning(0.7, [0.4, 0.4, 0.2], 576, 210, 310, 50)
InsertNs(0.02, 0.5)
ShortDValidation()
InsertIndels(0.5, 5, 0.5, 0.5)
```

### 4. Available Data Configurations

**Current Built-in Configurations:**
- `HUMAN_IGH_OGRDB` - Human heavy chain immunoglobulin data
- `HUMAN_IGK_OGRDB` - Human kappa light chain immunoglobulin data  
- `HUMAN_IGL_OGRDB` - Human lambda light chain immunoglobulin data
- `HUMAN_TCRB_IMGT` - Human T-cell receptor beta chain data

### 5. Accessing Alleles

**Old Syntax:**
```python
specific_v = data_cfg.allele_list('v')[0]
specific_d = data_cfg.allele_list('d')[0]
specific_j = data_cfg.allele_list('j')[0]
```

**New Syntax:**
```python
specific_v = HUMAN_IGH_OGRDB.v_alleles['IGHV1-2*02'][0]
specific_d = HUMAN_IGH_OGRDB.d_alleles['IGHD3-10*01'][0]
specific_j = HUMAN_IGH_OGRDB.j_alleles['IGHJ4*02'][0]
```

## Complete Migration Example

### Before (Old Syntax):
```python
from GenAIRR.pipeline import AugmentationPipeline
from GenAIRR.parameters import ChainType, CHAIN_TYPE_INFO
from GenAIRR.steps import SimulateSequence, FixVPositionAfterTrimmingIndexAmbiguity
from GenAIRR.mutation import S5F
from GenAIRR.data import builtin_heavy_chain_data_config
from GenAIRR.steps.StepBase import AugmentationStep

data_cfg = builtin_heavy_chain_data_config()
AugmentationStep.set_dataconfig(config=data_cfg, chain_type=ChainType.BCR_HEAVY)

pipeline = AugmentationPipeline([
    SimulateSequence(mutation_model=S5F(min_mutation_rate=0.003, max_mutation_rate=0.25), productive=True),
    FixVPositionAfterTrimmingIndexAmbiguity(),
    CorruptSequenceBeginning(
        corruption_probability=0.7,
        corrupt_events_proba=[0.4, 0.4, 0.2],
        max_sequence_length=576,
        nucleotide_add_coefficient=210,
        nucleotide_remove_coefficient=310,
        nucleotide_add_after_remove_coefficient=50,
        random_sequence_add_proba=1
    ),
    InsertNs(n_ratio=0.02, proba=0.5),
    InsertIndels(indel_probability=0.5, max_indels=5, insertion_proba=0.5, deletion_proba=0.5),
])
```

### After (New Syntax):
```python
from GenAIRR.pipeline import AugmentationPipeline
from GenAIRR.steps import SimulateSequence, FixVPositionAfterTrimmingIndexAmbiguity
from GenAIRR.mutation import S5F
from GenAIRR.data import HUMAN_IGH_OGRDB
from GenAIRR.steps.StepBase import AugmentationStep

AugmentationStep.set_dataconfig(HUMAN_IGH_OGRDB)

pipeline = AugmentationPipeline([
    SimulateSequence(S5F(min_mutation_rate=0.003, max_mutation_rate=0.25), True),
    FixVPositionAfterTrimmingIndexAmbiguity(),
    CorruptSequenceBeginning(0.7, [0.4, 0.4, 0.2], 576, 210, 310, 50),
    InsertNs(0.02, 0.5),
    InsertIndels(0.5, 5, 0.5, 0.5),
])
```

## Benefits of the New Syntax

1. **Simpler Imports**: Fewer imports needed, more direct access to data configurations
2. **Cleaner Code**: Positional arguments reduce verbosity 
3. **Better Performance**: Direct access to pre-loaded configurations
4. **Consistency**: Uniform parameter patterns across all steps
5. **Type Safety**: More predictable parameter types and order

## Troubleshooting Common Issues

### Issue: `AttributeError: module has no attribute 'builtin_heavy_chain_data_config'`
**Solution**: Replace with `from GenAIRR.data import HUMAN_IGH_OGRDB`

### Issue: `TypeError: takes 3 positional arguments but 4 were given`
**Solution**: Remove named parameters and use positional arguments

### Issue: `AttributeError: 'DataConfig' object has no attribute 'allele_list'`
**Solution**: Use the new allele access pattern: `config.v_alleles['family_name'][index]`

## Need Help?

If you encounter issues during migration:
1. Check the [test files](../tests/) for current usage examples
2. Review the [tutorials](./tutorials/) for updated code patterns
3. Open an issue on the [GitHub repository](https://github.com/MuteJester/GenAIRR/issues)
