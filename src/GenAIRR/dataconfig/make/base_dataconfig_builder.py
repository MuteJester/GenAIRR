import logging
import re
from collections import defaultdict
from datetime import date

from GenAIRR.dataconfig.make.auxiliary_builders.trimming_proportion_builder import TrimmingProbabilityGenerator
from ..data_config import DataConfig
from ..config_info import ConfigInfo
from ..enums import ChainType
from ...utilities.data_utilities import create_allele_dict
from ...utilities.asc_utilities import create_asc_germline_set

logger = logging.getLogger(__name__)


class BaseDataConfigGenerator:
    """
    Base class for DataConfig generation, containing common functionality
    between CustomDataConfigBuilder and RandomDataConfigBuilder.
    """

    def __init__(self, convert_to_asc=False, *, species=None, chain_type=None, reference_set=None):
        """
        Initialize the base generator.

        Args:
            convert_to_asc: Whether to convert to ASC (Allele Specific Correction) format.
            species: Optional Species enum value (e.g., ``Species.HUMAN``).
                When both *species* and *chain_type* are provided, the builder
                auto-sets ``ConfigInfo`` metadata on the resulting DataConfig.
            chain_type: Optional ChainType enum value (e.g., ``ChainType.BCR_HEAVY``).
            reference_set: Optional reference set name (e.g., ``"OGRDB 2024-01"``).
                Defaults to ``"user-provided"`` when metadata is auto-set.
        """
        self.convert_to_asc = convert_to_asc
        self.v_asc_table = None
        self.has_d = False
        # initialize an empty dataconfig object this class will populate
        self.dataconfig = DataConfig()
        # initialize auxiliary list for allele tracking
        self.alleles = ['V', 'J', 'C']

        # Metadata params (set ConfigInfo in _finalize if both species+chain_type provided)
        self._species = species
        self._chain_type = chain_type
        self._reference_set = reference_set

    def _get_reference_pointers(self):
        """
        Get references to the allele dictionaries.

        Returns:
            Dictionary mapping gene types to their allele collections
        """
        pointer_to_reference = {
            'V': self.dataconfig.v_alleles,
            'J': self.dataconfig.j_alleles,
        }
        if self.dataconfig.c_alleles:
            pointer_to_reference['C'] = self.dataconfig.c_alleles
        if self.has_d:
            pointer_to_reference['D'] = self.dataconfig.d_alleles
        return pointer_to_reference

    def _load_alleles(self, v_alleles, j_alleles, c_alleles=None, d_alleles=None):
        """
        Load alleles into the dataconfig object.

        Args:
            v_alleles: Dictionary of V gene alleles
            j_alleles: Dictionary of J gene alleles
            c_alleles: Dictionary of C gene alleles (optional)
            d_alleles: Dictionary of D gene alleles (optional)
        """
        self.dataconfig.v_alleles = v_alleles
        self.dataconfig.d_alleles = d_alleles
        self.dataconfig.j_alleles = j_alleles
        self.dataconfig.c_alleles = c_alleles
        # Update alleles tracking list based on what's actually loaded
        if c_alleles and 'C' not in self.alleles:
            self.alleles.append('C')
        elif not c_alleles and 'C' in self.alleles:
            self.alleles.remove('C')

    def _load_trimming_probs(self, max_trim=50):
        if self.has_d and 'D' not in self.alleles:
            self.alleles.append('D')
        generator = TrimmingProbabilityGenerator(
            self.dataconfig,
            allele_types=self.alleles,
            reference_getter=self._get_reference_pointers
        )
        generator.load_trimming_proportions(max_trim=max_trim)

    def load_reference(self, ref, segment=None, *, v_anchor_finder=None, j_anchor_finder=None,
                       keep_anchorless=False):
        """Load a reference from a file path or pre-loaded dict.

        Args:
            ref: FASTA file path (str) or pre-loaded allele dict.
            segment: Segment type ('V', 'D', 'J', 'C') passed through
                to ``create_allele_dict`` as ``segment_type``.
            v_anchor_finder: Custom V anchor callable (forwarded to
                ``create_allele_dict``).
            j_anchor_finder: Custom J anchor callable (forwarded to
                ``create_allele_dict``).
            keep_anchorless: If True, keep V/J alleles without anchors
                (forwarded to ``create_allele_dict``).
        """
        if isinstance(ref, str):
            if self.convert_to_asc and segment == "V":
                ref_data, asc_table = create_asc_germline_set(ref, segment="V")
                self.dataconfig.asc_tables['V'] = asc_table
                return ref_data
            return create_allele_dict(
                ref,
                segment_type=segment,
                v_anchor_finder=v_anchor_finder,
                j_anchor_finder=j_anchor_finder,
                keep_anchorless=keep_anchorless,
            )
        return ref

    def _read_reference_files(self, v_reference_path, j_reference_path,
                              c_reference_path, d_reference_path=None,
                              *, v_anchor_finder=None, j_anchor_finder=None,
                              keep_anchorless=False):
        """
        Read reference files and convert them to the appropriate format.

        Args:
            v_reference_path (str or dict): V gene reference file path or preloaded dict
            j_reference_path (str or dict): J gene reference file path or preloaded dict
            c_reference_path (str or dict): C gene reference file path or preloaded dict
            d_reference_path (str or dict, optional): D gene reference file path or preloaded dict
            v_anchor_finder: Custom V anchor callable
            j_anchor_finder: Custom J anchor callable
            keep_anchorless: If True, keep V/J alleles without anchors

        Returns:
            Tuple of (v_reference, j_reference, c_reference, d_reference)
        """

        # Handle optional D gene usage
        self.has_d = d_reference_path is not None
        if self.has_d and 'D' not in self.alleles:
            self.alleles.append('D')

        # Load all references with explicit segment types
        user_v_reference = self.load_reference(
            v_reference_path, segment="V",
            v_anchor_finder=v_anchor_finder,
            keep_anchorless=keep_anchorless,
        )
        user_j_reference = self.load_reference(
            j_reference_path, segment="J",
            j_anchor_finder=j_anchor_finder,
            keep_anchorless=keep_anchorless,
        )
        user_c_reference = (
            self.load_reference(c_reference_path, segment="C")
            if c_reference_path is not None else None
        )
        user_d_reference = (
            self.load_reference(d_reference_path, segment="D")
            if self.has_d else None
        )

        return user_v_reference, user_j_reference, user_c_reference, user_d_reference

    # ─────────────────────────────────────────────────────────
    # Metadata and finalization
    # ─────────────────────────────────────────────────────────

    def _set_metadata(self):
        """Set ConfigInfo metadata on the DataConfig if species and chain_type were provided."""
        if self._species is not None and self._chain_type is not None:
            self.dataconfig.metadata = ConfigInfo(
                species=self._species,
                chain_type=self._chain_type,
                reference_set=self._reference_set or "user-provided",
                last_updated=date.today(),
                has_d=self.has_d,
            )
            logger.info("Metadata set: %s", self.dataconfig.metadata)
        elif self._species is not None or self._chain_type is not None:
            logger.warning(
                "Both species and chain_type must be provided to auto-set metadata. "
                "Got species=%r, chain_type=%r. Metadata not set.",
                self._species, self._chain_type,
            )

    def _auto_generate_dj_pairing_map(self):
        """Auto-generate D-J pairing map for TCR-beta from allele names.

        Groups D and J alleles by the numeric family prefix:
        TRBD1 -> TRBJ1-*, TRBD2 -> TRBJ2-*, etc.

        Only runs when ``chain_type == ChainType.TCR_BETA`` and D alleles
        are present.
        """
        if self._chain_type != ChainType.TCR_BETA or not self.has_d:
            return
        if not self.dataconfig.d_alleles:
            return

        dj_map = {}
        for d_gene in self.dataconfig.d_alleles:
            match = re.search(r'(\d+)', d_gene)
            if match:
                num = match.group(1)
                compatible = [
                    j for j in self.dataconfig.j_alleles
                    if re.search(rf'J{num}', j)
                ]
                if compatible:
                    dj_map[d_gene] = compatible

        if dj_map:
            self.dataconfig.dj_pairing_map = dj_map
            logger.info("Auto-generated D-J pairing map: %s", dj_map)

    def _finalize(self):
        """Set metadata, generate D-J pairing map if applicable, and validate."""
        self._set_metadata()
        if self._chain_type is not None:
            self._auto_generate_dj_pairing_map()
        self.dataconfig.validate()
