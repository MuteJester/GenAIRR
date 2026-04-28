from dataclasses import dataclass, field
from typing import Optional, Dict, List, Any

from GenAIRR.dataconfig.config_info import ConfigInfo
from GenAIRR.dataconfig.enums import ChainType
import copy
import hashlib
import pickle
from GenAIRR.alleles.allele import Allele


DEFAULT_P_NUCLEOTIDE_LENGTH_PROBS = {0: 0.50, 1: 0.25, 2: 0.15, 3: 0.07, 4: 0.03}

# Bump on any breaking change to DataConfig fields. Old pickles without
# this version (or with an older one) refuse to load and prompt the user
# to re-migrate. See _check_schema() / verify_integrity() below.
SCHEMA_VERSION = 1


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

    # Schema version + integrity checksum (T2-1). Populated by the migration
    # script for all builtins; pickles without these fields fail load-time
    # verification with a clear migration hint. schema_sha256 is sha256 over
    # pickle.dumps(self, protocol=4) with schema_sha256 set to "" — see
    # compute_checksum() below.
    schema_version: int = SCHEMA_VERSION
    schema_sha256: str = ""

    # Build-time provenance (T2-8 closure step 3). Populated by
    # `RandomDataConfigBuilder.make_from_reference` with a BuildReport
    # describing per-allele anchor resolution outcomes, rejected
    # alleles + reasons, and bucket counts. ``None`` for legacy
    # builtins built before this field existed.
    build_report: Optional[Any] = None

    def __getattr__(self, name):
        # Backward-compat shim for pickled DataConfigs missing post-v1
        # fields. Note: schema_version / schema_sha256 fall through to
        # the dataclass class-level defaults (SCHEMA_VERSION and "")
        # when the instance __dict__ doesn't have them, so they don't
        # need entries here — verify_integrity() then sees version=1
        # but checksum="", which trips the missing-checksum branch and
        # produces a clean migration error.
        if name == 'p_nucleotide_length_probs':
            return dict(DEFAULT_P_NUCLEOTIDE_LENGTH_PROBS)
        if name == 'dj_pairing_map':
            return None
        if name == 'build_report':
            return None
        raise AttributeError(f"'{type(self).__name__}' object has no attribute '{name}'")

    def compute_checksum(self) -> str:
        """Compute the canonical sha256 of this DataConfig's contents.

        The hash is taken over ``pickle.dumps(self, protocol=4)`` with
        ``schema_sha256`` temporarily zeroed (so the field cannot
        include itself in its own hash) and ``build_report`` removed
        from ``__dict__`` entirely (so diagnostic provenance — which
        has no semantic effect on simulation — doesn't invalidate
        the integrity check, AND legacy pickles built before
        ``build_report`` existed continue to verify against their
        original stored checksums).
        """
        saved_sha = self.schema_sha256
        # KEY-LEVEL pop: legacy pickles don't have 'build_report' in
        # __dict__ at all. Setting self.build_report = None would
        # ADD the key, changing the pickled bytes. Removing it from
        # __dict__ instead keeps the wire shape identical for both
        # legacy and new pickles.
        had_report = 'build_report' in self.__dict__
        saved_report = self.__dict__.get('build_report')

        self.schema_sha256 = ""
        if had_report:
            del self.__dict__['build_report']
        try:
            blob = pickle.dumps(self, protocol=4)
            return hashlib.sha256(blob).hexdigest()
        finally:
            self.schema_sha256 = saved_sha
            if had_report:
                self.__dict__['build_report'] = saved_report

    def verify_integrity(self) -> None:
        """Validate schema_version and schema_sha256.

        Raises ``DataConfigError`` if the pickle predates the current
        schema version or if its stored checksum does not match the
        recomputed checksum (corrupted/tampered file).
        """
        version = getattr(self, "schema_version", 0)
        if version != SCHEMA_VERSION:
            raise DataConfigError(
                f"DataConfig '{self.name or 'Unnamed'}' has schema_version="
                f"{version}, but the current code expects "
                f"{SCHEMA_VERSION}. Re-migrate this pickle (see "
                f".private/scripts/migrate_pkl_to_v1.py) or rebuild it "
                f"from the IMGT source via dataconfig builders."
            )
        stored = getattr(self, "schema_sha256", "") or ""
        if not stored:
            raise DataConfigError(
                f"DataConfig '{self.name or 'Unnamed'}' is missing "
                f"schema_sha256. Re-migrate this pickle to populate the "
                f"checksum field."
            )
        actual = self.compute_checksum()
        if actual != stored:
            raise DataConfigError(
                f"DataConfig '{self.name or 'Unnamed'}' checksum mismatch: "
                f"stored={stored[:16]}…, actual={actual[:16]}…. "
                f"The pickle has been modified or corrupted since it was "
                f"saved. Reinstall GenAIRR or rebuild this config from "
                f"its IMGT source."
            )

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