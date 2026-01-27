"""
Mutation Models for GenAIRR

This module provides mutation models for simulating somatic hypermutation
in immunoglobulin sequences.

Available Models:
    - Uniform: Applies mutations uniformly at random positions
    - S5F: Context-dependent 5-mer mutation model

Core Classes:
    - MutationModel: Abstract base class for all mutation models
    - MutationContext: Input data for mutation models
    - MutationResult: Output from mutation models

Components:
    - StopCodonChecker: Utility for stop codon detection
    - AnchorChecker: Utility for anchor residue preservation
"""

from .mutation_model import MutationModel
from .context import MutationContext, ChainType
from .result import MutationResult
from .uniform import Uniform
from .s5f import S5F

__all__ = [
    # Abstract base
    "MutationModel",
    # Context and result
    "MutationContext",
    "MutationResult",
    "ChainType",
    # Nucleotide-level models
    "Uniform",
    "S5F",
]
