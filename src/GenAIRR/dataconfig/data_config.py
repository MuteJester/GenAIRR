import pickle
from dataclasses import dataclass, field
from typing import Optional, Dict, List, Any

from GenAIRR.dataconfig.config_info import ConfigInfo
import copy
from GenAIRR.alleles.allele import Allele

@dataclass
class DataConfig:
    """
    Configuration class for storing data related to sequence generation, allele usage, trimming, and mutation rates.

    This class encapsulates various dictionaries and settings used in the simulation and analysis of immunoglobulin sequences.
    """
    # --- Attributes ---
    name: Optional[str] = None
    metadata:ConfigInfo = None

    # Gene usage frequencies
    family_use_dict: Dict[str, float] = field(default_factory=dict)
    gene_use_dict: Dict[str, float] = field(default_factory=dict)

    # Allele dictionaries (grouped by family)
    v_alleles: Optional[Dict[str, List[Allele]]] = None
    d_alleles: Optional[Dict[str, List[Allele]]] = None
    j_alleles: Optional[Dict[str, List[Allele]]] = None
    c_alleles: Optional[Dict[str, List[Allele]]] = None

    # Trimming and NP region parameters
    trim_dicts: Dict[str, Any] = field(default_factory=dict)
    NP_transitions: Dict[str, Any] = field(default_factory=dict)
    NP_first_bases: Dict[str, Any] = field(default_factory=dict)
    NP_lengths: Dict[str, Any] = field(default_factory=dict)

    # Simulation and analysis parameters
    mut_rate_per_seq: Dict[str, float] = field(default_factory=dict)
    kmer_dicts: Dict[str, Any] = field(default_factory=dict)
    correction_maps: Dict[str, Any] = field(default_factory=dict)
    asc_tables: Dict[str, Any] = field(default_factory=dict)

    def _unfold_alleles(self, gene_segment: str) -> List[str]:
        """Unfolds the alleles for a given gene segment (v, d, j, c) into a flat list."""
        alleles_dict = getattr(self, f"{gene_segment}_alleles")
        if alleles_dict is None:
            return []
        # This comprehension is clearer: iterate through the lists of alleles and flatten them.
        return [allele for allele_list in alleles_dict.values() for allele in allele_list]

    def _count_alleles(self, gene_segment: str) -> int:
        """Counts the number of alleles for a given gene segment."""
        return len(self._unfold_alleles(gene_segment))

    # --- Public Properties ---
    @property
    def number_of_v_alleles(self) -> int:
        return self._count_alleles('v')

    @property
    def number_of_d_alleles(self) -> int:
        return self._count_alleles('d')

    @property
    def number_of_j_alleles(self) -> int:
        return self._count_alleles('j')

    @property
    def number_of_c_alleles(self) -> int:
        return self._count_alleles('c')

    # --- Public Methods ---
    def allele_list(self, gene_segment: str) -> List[str]:
        """Returns a flattened list of all alleles for a given gene segment."""
        return self._unfold_alleles(gene_segment)

    def copy(self):
        """
        Creates a deep, independent copy of this DataConfig object.

        Returns:
            DataConfig: A new DataConfig object with all attributes and nested
                        data structures duplicated.
        """
        return copy.deepcopy(self)

    def __repr__(self) -> str:
        parts = [f"<{self.name or 'Unnamed'} - Data Config>"]
        # Check each allele type before adding it to the representation
        if self.v_alleles is not None:
            parts.append(f"<{self.number_of_v_alleles} V Alleles>")
        if self.d_alleles is not None:
            parts.append(f"<{self.number_of_d_alleles} D Alleles>")
        if self.j_alleles is not None:
            parts.append(f"<{self.number_of_j_alleles} J Alleles>")
        if self.c_alleles is not None:
            parts.append(f"<{self.number_of_c_alleles} C Alleles>")
        return "-".join(parts)