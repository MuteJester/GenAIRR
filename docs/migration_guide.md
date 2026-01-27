# Migration Guide: Updating to Current GenAIRR Syntax

This guide helps users migrate from older GenAIRR syntax to the current API. The changes improve usability and consistency across the library.

!!! warning "Action required"
    If you are using `AugmentationStep.set_dataconfig()` or positional arguments for step constructors, your code will still work but emits deprecation warnings. Update to the new syntax shown below to avoid breakage in future releases.

## Key Changes Summary

### 1. Data Configuration Loading (Recommended: Pipeline Config)

The recommended approach is to pass the config directly to the Pipeline:

**Old Syntax (Deprecated):**
```python
from GenAIRR.data import HUMAN_IGH_OGRDB
from GenAIRR.steps.StepBase import AugmentationStep

AugmentationStep.set_dataconfig(HUMAN_IGH_OGRDB)  # Global state - deprecated!
```

**New Syntax (Recommended):**
```python
from GenAIRR import Pipeline, steps, HUMAN_IGH_OGRDB, S5F

# Config passed directly to Pipeline - no global state
pipeline = Pipeline(
    config=HUMAN_IGH_OGRDB,
    steps=[...]
)
```

**Or use the convenience function for simple cases:**
```python
from GenAIRR import simulate, HUMAN_IGH_OGRDB, S5F

result = simulate(HUMAN_IGH_OGRDB, S5F(0.003, 0.25))
```

### 1b. Legacy Data Configuration (Very Old Syntax)

**Very Old Syntax:**
```python
from GenAIRR.data import builtin_heavy_chain_data_config
from GenAIRR.parameters import ChainType, CHAIN_TYPE_INFO

data_cfg = builtin_heavy_chain_data_config()
AugmentationStep.set_dataconfig(config=data_cfg, chain_type=ChainType.BCR_HEAVY)
```

**Intermediate Syntax (Still works, but deprecated):**
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
SimulateSequence(S5F(min_mutation_rate=0.003, max_mutation_rate=0.25), productive=True)
```

**Note:** The `productive` parameter is now a keyword argument for clarity.

### 3. Pipeline Step Parameters (Keyword-Only with Defaults)

Steps now use **keyword-only arguments** with sensible defaults. You can use defaults or customize only what you need.

**Old Syntax (Positional - No longer supported):**
```python
CorruptSequenceBeginning(0.7, [0.4, 0.4, 0.2], 576, 210, 310, 50)
InsertNs(0.02, 0.5)
InsertIndels(0.5, 5, 0.5, 0.5)
```

**New Syntax (Keyword-Only with Defaults):**
```python
# Use sensible defaults - simplest usage
CorruptSequenceBeginning()
InsertNs()
InsertIndels()

# Or customize specific parameters as needed
CorruptSequenceBeginning(probability=0.9)
EnforceSequenceLength(max_length=400)  # Separate step for length enforcement
InsertNs(n_ratio=0.05, probability=0.8)
InsertIndels(probability=0.3, max_indels=3)
```

**Available Parameters:**

`CorruptSequenceBeginning`:
- `probability` (default: 0.7) - Probability of corruption
- `event_weights` (default: (0.4, 0.4, 0.2)) - Weights for [add, remove, remove_then_add]
- `nucleotide_add_coefficient` (default: 210)
- `nucleotide_remove_coefficient` (default: 310)
- `nucleotide_add_after_remove_coefficient` (default: 50)

`EnforceSequenceLength`:
- `max_length` (default: 576) - Maximum sequence length (simulates sequencing read limits)

`InsertNs`:
- `n_ratio` (default: 0.02) - Ratio of N's to insert
- `probability` (default: 0.5) - Probability step is applied

`InsertIndels`:
- `probability` (default: 0.5) - Probability indels are inserted
- `max_indels` (default: 5) - Maximum number of indels
- `insertion_probability` (default: 0.5) - Weight for insertions
- `deletion_probability` (default: 0.5) - Weight for deletions

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

### Before (Very Old Syntax):
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

### Intermediate (Deprecated - still works with warning):
```python
from GenAIRR.pipeline import AugmentationPipeline
from GenAIRR.steps import SimulateSequence, FixVPositionAfterTrimmingIndexAmbiguity
from GenAIRR.mutation import S5F
from GenAIRR.data import HUMAN_IGH_OGRDB
from GenAIRR.steps.StepBase import AugmentationStep

AugmentationStep.set_dataconfig(HUMAN_IGH_OGRDB)  # Deprecated!

pipeline = AugmentationPipeline([
    SimulateSequence(S5F(min_mutation_rate=0.003, max_mutation_rate=0.25), True),
    FixVPositionAfterTrimmingIndexAmbiguity(),
    CorruptSequenceBeginning(0.7, [0.4, 0.4, 0.2], 576, 210, 310, 50),
    InsertNs(0.02, 0.5),
    InsertIndels(0.5, 5, 0.5, 0.5),
])
```

### After (Recommended - Current Syntax):
```python
from GenAIRR import Pipeline, steps, HUMAN_IGH_OGRDB, S5F

# Config passed directly to Pipeline - clean, explicit, no global state
# Steps use keyword-only arguments with sensible defaults
pipeline = Pipeline(
    config=HUMAN_IGH_OGRDB,
    steps=[
        steps.SimulateSequence(S5F(min_mutation_rate=0.003, max_mutation_rate=0.25), productive=True),
        steps.FixVPositionAfterTrimmingIndexAmbiguity(),
        steps.CorruptSequenceBeginning(),    # 5' corruption (add/remove nucleotides)
        steps.EnforceSequenceLength(),       # enforce max read length (default: 576)
        steps.InsertNs(),                    # uses defaults, or customize: n_ratio=0.02
        steps.InsertIndels(),                # uses defaults, or customize: probability=0.5
    ]
)
```

### Simplest (For Common Use Cases):
```python
from GenAIRR import simulate, HUMAN_IGH_OGRDB, S5F

# One-liner for quick simulations
result = simulate(HUMAN_IGH_OGRDB, S5F(0.003, 0.25))

# Generate multiple sequences
results = simulate(HUMAN_IGH_OGRDB, S5F(0.003, 0.25), n=100)
```

## Benefits of the New Syntax

1. **No Global State**: Config is passed explicitly to Pipeline, avoiding hidden dependencies
2. **Simpler Imports**: Single import line `from GenAIRR import ...` covers most use cases
3. **Cleaner Code**: Positional arguments reduce verbosity
4. **Concurrent-Safe**: Multiple pipelines can use different configs simultaneously
5. **Easier Debugging**: Config ownership is clear - always check the Pipeline
6. **Convenience Function**: `simulate()` for common one-liner use cases

## Troubleshooting Common Issues

### Issue: `DeprecationWarning: Using class-level AugmentationStep.set_dataconfig() is deprecated`
**Solution**: Pass config directly to the Pipeline constructor:
```python
# Instead of:
AugmentationStep.set_dataconfig(HUMAN_IGH_OGRDB)
pipeline = AugmentationPipeline([...])

# Use:
pipeline = Pipeline(config=HUMAN_IGH_OGRDB, steps=[...])
```

### Issue: `ValueError: No DataConfig provided`
**Solution**: Make sure to pass `config=` to the Pipeline:
```python
pipeline = Pipeline(config=HUMAN_IGH_OGRDB, steps=[...])
```

### Issue: `AttributeError: module has no attribute 'builtin_heavy_chain_data_config'`
**Solution**: Use the new imports: `from GenAIRR import HUMAN_IGH_OGRDB`

### Issue: `TypeError: takes 3 positional arguments but 4 were given`
**Solution**: Remove named parameters and use positional arguments

### Issue: `AttributeError: 'DataConfig' object has no attribute 'allele_list'`
**Solution**: Use the new allele access pattern: `config.v_alleles['family_name'][index]`

## Need Help?

If you encounter issues during migration:
1. Check the [step-by-step tutorial](step_by_step_tutorial.md) for current usage examples
2. Review the [API reference](api_reference.md) for updated code patterns
3. Open an issue on the [GitHub repository](https://github.com/MuteJester/GenAIRR/issues)
