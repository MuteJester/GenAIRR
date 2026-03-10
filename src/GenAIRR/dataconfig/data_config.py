from dataclasses import dataclass, field
from typing import Optional, Dict, List, Any

from GenAIRR.dataconfig.config_info import ConfigInfo
from GenAIRR.dataconfig.enums import ChainType
import copy
from GenAIRR.alleles.allele import Allele


DEFAULT_P_NUCLEOTIDE_LENGTH_PROBS = {0: 0.50, 1: 0.25, 2: 0.15, 3: 0.07, 4: 0.03}


class DataConfigError(Exception):
    """Raised when a DataConfig fails validation."""


@dataclass
class DataConfig:
    """
    Configuration class for storing data related to sequence generation, allele usage, trimming, and mutation rates.

    This class encapsulates various dictionaries and settings used in the simulation and analysis of immunoglobulin sequences.
    """
    # --- Attributes ---
    name: Optional[str] = None
    metadata: Optional[ConfigInfo] = None

    # Gene usage frequencies
    gene_use_dict: Dict[str, Any] = field(default_factory=dict)

    # Allele dictionaries (grouped by gene name → list of allele objects)
    v_alleles: Optional[Dict[str, List[Allele]]] = None
    d_alleles: Optional[Dict[str, List[Allele]]] = None
    j_alleles: Optional[Dict[str, List[Allele]]] = None
    c_alleles: Optional[Dict[str, List[Allele]]] = None

    # Trimming and NP region parameters
    trim_dicts: Dict[str, Any] = field(default_factory=dict)
    NP_transitions: Dict[str, Any] = field(default_factory=dict)
    NP_first_bases: Dict[str, Any] = field(default_factory=dict)
    NP_lengths: Dict[str, Any] = field(default_factory=dict)

    # Correction maps and ASC tables
    correction_maps: Dict[str, Any] = field(default_factory=dict)
    asc_tables: Dict[str, Any] = field(default_factory=dict)

    # P-nucleotide length distribution (0-4 bp, geometric decay)
    p_nucleotide_length_probs: Dict[int, float] = field(
        default_factory=lambda: dict(DEFAULT_P_NUCLEOTIDE_LENGTH_PROBS)
    )

    # D-J pairing map: {D_gene_name: [compatible_J_gene_name, ...]}
    # When None, no D-J family pairing constraint is applied.
    dj_pairing_map: Optional[Dict[str, list]] = None

    def __getattr__(self, name):
        # Backward compat for pickled DataConfigs missing new fields
        if name == 'p_nucleotide_length_probs':
            return dict(DEFAULT_P_NUCLEOTIDE_LENGTH_PROBS)
        if name == 'dj_pairing_map':
            return None
        raise AttributeError(f"'{type(self).__name__}' object has no attribute '{name}'")

    def validate(self):
        """
        Validate that this DataConfig has the minimum required fields for simulation.

        Raises:
            DataConfigError: If any required field is missing or malformed.
        """
        errors = []

        # Metadata
        if self.metadata is not None and not isinstance(self.metadata, ConfigInfo):
            errors.append(f"metadata must be a ConfigInfo instance, got {type(self.metadata).__name__}")

        # Alleles
        if not self.v_alleles or len(self.v_alleles) == 0:
            errors.append("v_alleles is required and must be non-empty")
        if not self.j_alleles or len(self.j_alleles) == 0:
            errors.append("j_alleles is required and must be non-empty")

        # D alleles: required when metadata says has_d
        if self.metadata and self.metadata.has_d:
            if not self.d_alleles or len(self.d_alleles) == 0:
                errors.append("d_alleles is required for chains with D segments (metadata.has_d=True)")
        # Consistency: has_d=False but d_alleles provided
        if self.metadata and not self.metadata.has_d:
            if self.d_alleles and len(self.d_alleles) > 0:
                errors.append(
                    "d_alleles is non-empty but metadata.has_d=False. "
                    "Either remove d_alleles or set has_d=True."
                )

        # Gene usage
        if not self.gene_use_dict:
            errors.append("gene_use_dict is required and must be non-empty")
        else:
            for key in ("V", "J"):
                if key not in self.gene_use_dict:
                    errors.append(f"gene_use_dict must contain '{key}' key")

        # Trim dicts
        if not self.trim_dicts:
            errors.append("trim_dicts is required and must be non-empty")
        else:
            for key in ("V_3", "J_5"):
                if key not in self.trim_dicts:
                    errors.append(f"trim_dicts must contain '{key}' key")
            if self.metadata and self.metadata.has_d:
                for key in ("D_5", "D_3"):
                    if key not in self.trim_dicts:
                        errors.append(f"trim_dicts must contain '{key}' key for chains with D segments")

        if errors:
            raise DataConfigError(
                f"DataConfig '{self.name or 'Unnamed'}' failed validation:\n  - " + "\n  - ".join(errors)
            )

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