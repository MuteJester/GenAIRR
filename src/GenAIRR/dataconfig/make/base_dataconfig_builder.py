from collections import defaultdict

from GenAIRR.dataconfig.make.auxiliary_builders.correction_map_builder import CorrectionMapBuilder
from GenAIRR.dataconfig.make.auxiliary_builders.trimming_proportion_builder import TrimmingProbabilityGenerator
from ..data_config import DataConfig  # Import your DataConfig class
from ...utilities.data_utilities import create_allele_dict
from ...utilities.asc_utilities import create_asc_germline_set





class BaseDataConfigGenerator:
    """
    Base class for DataConfig generation, containing common functionality
    between CustomDataConfigBuilder and RandomDataConfigBuilder.
    """

    def __init__(self, convert_to_asc=False):
        """
        Initialize the base generator.

        Args:
            convert_to_asc: Whether to convert to ASC (Allele Specific Correction) format
        """
        self.convert_to_asc = convert_to_asc
        self.v_asc_table = None
        self.has_d = False
        # initialize an empty dataconfig object this class will populate
        self.dataconfig = DataConfig()
        self.correction_builder = CorrectionMapBuilder(self.dataconfig)
        # initialize auxiliary list for allele tracking
        self.alleles = ['V', 'J', 'C']

    def _get_reference_pointers(self):
        """
        Get references to the allele dictionaries.

        Returns:
            Dictionary mapping gene types to their allele collections
        """
        # aux list
        pointer_to_reference = {
            'V': self.dataconfig.v_alleles,
            'J': self.dataconfig.j_alleles,
            'C': self.dataconfig.c_alleles
        }
        if self.has_d:
            pointer_to_reference['D'] = self.dataconfig.d_alleles
        return pointer_to_reference

    def _load_alleles(self, v_alleles, j_alleles, c_alleles, d_alleles=None):
        """
        Load alleles into the dataconfig object.

        Args:
            v_alleles: Dictionary of V gene alleles
            j_alleles: Dictionary of J gene alleles
            c_alleles: Dictionary of C gene alleles
            d_alleles: Dictionary of D gene alleles (optional)
        """
        self.dataconfig.v_alleles = v_alleles
        self.dataconfig.d_alleles = d_alleles
        self.dataconfig.j_alleles = j_alleles
        self.dataconfig.c_alleles = c_alleles

    def _load_trimming_probs(self, max_trim=50):
        if self.has_d and 'D' not in self.alleles:
            self.alleles.append('D')
        generator = TrimmingProbabilityGenerator(
            self.dataconfig,
            allele_types=self.alleles,
            reference_getter=self._get_reference_pointers
        )
        generator.load_trimming_proportions(max_trim=max_trim)

    def load_reference(self,ref, segment=None):
        if isinstance(ref, str):
            if self.convert_to_asc and segment == "V":
                ref_data, asc_table = create_asc_germline_set(ref, segment="V")
                self.dataconfig.asc_tables['V'] = asc_table
                return ref_data
            return create_allele_dict(ref)
        return ref

    def _read_reference_files(self, v_reference_path, j_reference_path, c_reference_path, d_reference_path=None):
        """
        Read reference files and convert them to the appropriate format.

        Args:
            v_reference_path (str or dict): V gene reference file path or preloaded dict
            j_reference_path (str or dict): J gene reference file path or preloaded dict
            c_reference_path (str or dict): C gene reference file path or preloaded dict
            d_reference_path (str or dict, optional): D gene reference file path or preloaded dict

        Returns:
            Tuple of (v_reference, j_reference, c_reference, d_reference)
        """

        # Handle optional D gene usage
        self.has_d = d_reference_path is not None
        if self.has_d and 'D' not in self.alleles:
            self.alleles.append('D')

        # Load all references
        user_v_reference = self.load_reference(v_reference_path, segment="V")
        user_j_reference = self.load_reference(j_reference_path)
        user_c_reference = self.load_reference(c_reference_path)
        user_d_reference = self.load_reference(d_reference_path) if self.has_d else None

        return user_v_reference, user_j_reference, user_c_reference, user_d_reference
