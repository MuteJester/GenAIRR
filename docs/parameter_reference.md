# Parameter Reference Guide

This document provides detailed explanations of key parameters used throughout GenAIRR.

## Pipeline Step Parameters

### SimulateSequence
- **mutation_model**: The mutation model to apply (S5F, Uniform, or custom)
- **productive**: If `True`, ensures sequences are in-frame and functional
- **specific_v/d/j**: Forces use of specific alleles instead of random selection

### CorruptSequenceBeginning
- **corrupt_proba** (0.7): Probability of applying corruption to a sequence
- **corrupt_events_proba** ([0.4, 0.4, 0.2]): Probabilities for [add, remove, both] events
- **max_sequence_length** (576): Maximum allowed sequence length
- **nucleotide_add_coefficient** (210): Controls distribution of added nucleotides
- **nucleotide_remove_coefficient** (310): Controls distribution of removed nucleotides
- **nucleotide_add_after_remove_coefficient** (50): Controls nucleotides added after removal

### InsertNs
- **n_ratio** (0.02): Fraction of sequence positions that can become 'N'
- **n_proba** (0.5): Probability of actually inserting an 'N' at eligible positions

### InsertIndels
- **indel_proba** (0.5): Probability of inserting indels
- **max_indels** (5): Maximum number of indels to insert
- **deletion_proba** (0.5): Probability of deletion vs insertion
- **insertion_proba** (0.5): Probability of insertion vs deletion

## Mutation Model Parameters

### S5F (Context-Dependent Model)
- **min_mutation_rate**: Minimum mutation frequency (e.g., 0.003 = 0.3%)
- **max_mutation_rate**: Maximum mutation frequency (e.g., 0.25 = 25%)
- **custom_model**: Path to custom mutation model file (optional)
- **productive**: If `True`, avoids mutations that create stop codons and preserves CDR3 anchors

### Uniform (Simple Model)
- **min_mutation_rate**: Minimum mutation frequency
- **max_mutation_rate**: Maximum mutation frequency
- **productive**: If `True`, avoids mutations that create stop codons and preserves CDR3 anchors

## Data Configuration

### Built-in Configs
- **HUMAN_IGH_OGRDB**: Heavy chain immunoglobulin data (OGRDB)
- **HUMAN_IGH_EXTENDED**: Extended heavy chain immunoglobulin data
- **HUMAN_IGK_OGRDB**: Kappa light chain data (OGRDB)
- **HUMAN_IGL_OGRDB**: Lambda light chain data (OGRDB)
- **HUMAN_TCRB_IMGT**: T-cell receptor beta data (IMGT)

## Performance Tips

- Use lower mutation rates (0.003-0.03) for realistic sequences
- Set `productive=True` for functional sequences only
- Adjust `max_sequence_length` based on your sequencing platform
- Use specific alleles only when you need controlled experiments
