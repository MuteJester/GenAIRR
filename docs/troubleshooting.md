# Troubleshooting Guide

Common issues and their solutions when using GenAIRR.

## Common Errors

### 1. "AttributeError" or "TypeError" when creating pipelines

!!! note "Most common cause"
    This almost always means you are using old API syntax. See the [Migration Guide](migration_guide.md) for a full before/after comparison.

**Problem**: Using the old API or wrong constructor syntax.

**Solution**: Pass `config` to the Pipeline constructor and use keyword arguments for steps:
```python
from GenAIRR import Pipeline, steps, HUMAN_IGH_OGRDB, S5F

pipeline = Pipeline(
    config=HUMAN_IGH_OGRDB,
    steps=[
        steps.SimulateSequence(S5F(min_mutation_rate=0.003, max_mutation_rate=0.25), productive=True),
        # ... more steps
    ]
)
```

### 2. "No functional sequences generated"

**Problem**: All generated sequences have stop codons or are out-of-frame.

**Solutions**:
- Set `productive=True` in SimulateSequence
- Lower your mutation rates
- Check your data configuration has functional alleles

```python
# Use productive sequences
steps.SimulateSequence(S5F(min_mutation_rate=0.003, max_mutation_rate=0.25), productive=True)
```

### 3. "Empty mutations dictionary"

**Problem**: Using Uniform(0, 0) or very low mutation rates.

**Solution**: Increase mutation rates or use S5F for realistic mutations:
```python
# Instead of this (no mutations)
steps.SimulateSequence(Uniform(min_mutation_rate=0, max_mutation_rate=0), productive=True)

# Use this for naive sequences with potential for mutations
steps.SimulateSequence(S5F(min_mutation_rate=0.001, max_mutation_rate=0.01), productive=True)
```

### 4. "Allele not found in data config"

**Problem**: Trying to use specific alleles that don't exist in your DataConfig.

**Solution**: Check available alleles first:
```python
# List available V allele families
print("V allele families:", list(HUMAN_IGH_OGRDB.v_alleles.keys())[:5])

# Access specific allele properly (using family name)
v_allele = HUMAN_IGH_OGRDB.v_alleles['IGHVF1-G1'][0]  # First allele in family
d_allele = HUMAN_IGH_OGRDB.d_alleles['IGHD1-1'][0]
j_allele = HUMAN_IGH_OGRDB.j_alleles['IGHJ1'][0]

# View all alleles in a family
for allele in HUMAN_IGH_OGRDB.v_alleles['IGHVF1-G1']:
    print(allele.name)
```

### 5. "Pipeline execution is very slow"

**Potential causes and solutions**:
- **High mutation rates**: Lower max_mutation_rate
- **Complex pipelines**: Remove unnecessary steps
- **Productive=True with high mutation rates**: The library keeps retrying until finding functional sequences

```python
# Faster execution
steps.SimulateSequence(S5F(min_mutation_rate=0.003, max_mutation_rate=0.05), productive=True)
```

## Performance Optimization

### Memory Usage
- Generate sequences in batches rather than all at once
- Use `get_dict()` only when needed (it creates copies)

### Speed Tips

!!! tip "Quick wins for faster generation"
    - Use built-in data configs (they're pre-optimized and lazily loaded)
    - Avoid very high mutation rates (>0.3) unless necessary
    - Use `EnforceSequenceLength` to cap sequence length
    - Set `productive=False` if you don't need the productivity guarantee

## Debugging Tips

### Check Your Pipeline
```python
# Print pipeline structure
for i, step in enumerate(pipeline.steps):
    print(f"Step {i}: {type(step).__name__}")
```

### Inspect Container State
```python
# After execution
container = pipeline.execute()
print("Sequence length:", len(container.sequence))
print("Mutations count:", len(container.mutations))
print("Is productive:", container.productive)
```

### Validate Data Config
```python
# Check if data config is loaded
print("V alleles count:", len(HUMAN_IGH_OGRDB.v_alleles))
print("Chain type:", HUMAN_IGH_OGRDB.metadata.chain_type)
```

## Getting Help

1. Check this troubleshooting guide first
2. Review the parameter reference for correct usage
3. Look at the examples in documentation
4. Check if your issue is covered in existing GitHub issues
5. Create a minimal reproducible example when reporting bugs
