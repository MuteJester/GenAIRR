class CorrectionMapBuilder:
    """
    Helper class for building various correction maps used in immunoglobulin sequence alignment.

    This class encapsulates logic for computing:
        - 5' trimming similarity maps
        - 3' trimming similarity maps
        - Combined 5' and 3' trimming maps
        - N-ambiguity correction graphs

    Correction maps are stored in the `correction_maps` attribute of the provided DataConfig instance.
    """

    def __init__(self, dataconfig):
        """
        Initialize the CorrectionMapBuilder.

        Args:
            dataconfig (DataConfig): The DataConfig object where correction maps will be stored.
        """
        self.dataconfig = dataconfig

    def _flatten_alleles(self, target_alleles):
        """
        Flatten nested allele dictionary to a flat name → allele object mapping.

        Args:
            target_alleles (dict): Dictionary mapping family/gene → list of allele objects.

        Returns:
            dict: Flattened dictionary mapping allele names to allele objects.
        """
        return {allele.name: allele for group in target_alleles.values() for allele in group}

    def build_5_prime_trim_map(self, target_alleles):
        """
        Build a 5' trimming correction map.

        For each allele, this map lists all alleles that contain a trimmed version of its
        sequence (trimmed from the 5' end) as a substring.

        Args:
            target_alleles (dict): Dictionary mapping family/gene → list of allele objects.
        """
        allele_dict = self._flatten_alleles(target_alleles)
        trim_map = {
            allele.name: {
                trim_5: [
                    other.name for other in allele_dict.values()
                    if allele.ungapped_seq[trim_5:] in other.ungapped_seq
                ]
                for trim_5 in range(len(allele.ungapped_seq) + 1)
            }
            for allele in allele_dict.values()
        }
        key = self._get_correction_key(allele_dict, "_5_TRIM_SIMILARITY_MAP")
        self.dataconfig.correction_maps[key] = trim_map

    def build_3_prime_trim_map(self, target_alleles):
        """
        Build a 3' trimming correction map.

        For each allele, this map lists all alleles that contain a trimmed version of its
        sequence (trimmed from the 3' end) as a substring.

        Args:
            target_alleles (dict): Dictionary mapping family/gene → list of allele objects.
        """
        allele_dict = self._flatten_alleles(target_alleles)
        trim_map = {
            allele.name: {
                trim_3: [
                    other.name for other in allele_dict.values()
                    if allele.ungapped_seq[:len(allele.ungapped_seq) - trim_3] in other.ungapped_seq
                ]
                for trim_3 in range(len(allele.ungapped_seq) + 1)
            }
            for allele in allele_dict.values()
        }
        key = self._get_correction_key(allele_dict, "_3_TRIM_SIMILARITY_MAP")
        self.dataconfig.correction_maps[key] = trim_map

    def build_5_and_3_prime_trim_map(self, target_alleles):
        """
        Build a combined 5' and 3' trimming correction map.

        For each allele, this map lists all alleles that contain a doubly-trimmed version
        of its sequence (trimmed from both ends) as a substring.

        Args:
            target_alleles (dict): Dictionary mapping family/gene → list of allele objects.
        """
        allele_dict = self._flatten_alleles(target_alleles)
        allele_list = list(allele_dict.values())
        trim_map = {}

        for allele in allele_list:
            seq = allele.ungapped_seq
            allele_map = {}
            for trim_5 in range(len(seq) + 1):
                for trim_3 in range(len(seq) - trim_5 + 1):
                    trimmed = seq[trim_5:len(seq) - trim_3]
                    matches = [
                        other.name for other in allele_list
                        if trimmed in other.ungapped_seq
                    ]
                    allele_map[(trim_5, trim_3)] = matches
            trim_map[allele.name] = allele_map

        key = self._get_correction_key(allele_dict, "_5_3_TRIM_SIMILARITY_MAP")
        self.dataconfig.correction_maps[key] = trim_map

    def build_n_ambiguity_map(self, target_alleles):
        """
        Build an N-ambiguity correction graph.

        This graph encodes similarity relationships between alleles where ambiguous base calls ('N') may
        match multiple alleles. It uses an external AlleleNComparer class to build the graph.

        Args:
            target_alleles (dict): Dictionary mapping family/gene → list of allele objects.
        """
        from ....utilities import AlleleNComparer
        allele_dict = self._flatten_alleles(target_alleles)
        comparer = AlleleNComparer()
        for name, allele in allele_dict.items():
            comparer.add_allele(name, allele.ungapped_seq.upper())
        key = self._get_correction_key(allele_dict, "_N_AMBIGUITY_CORRECTION_GRAPH")
        self.dataconfig.correction_maps[key] = comparer

    def _get_correction_key(self, allele_dict, suffix):
        """
        Generate a correction map key using the allele type and given suffix.

        Args:
            allele_dict (dict): Dictionary mapping allele name → allele object.
            suffix (str): Suffix to append to the key (e.g., '_5_TRIM_SIMILARITY_MAP').

        Returns:
            str: A full key like 'V_5_TRIM_SIMILARITY_MAP'
        """
        allele_type = str(next(iter(allele_dict.values())).type).split('.')[1]
        return f"{allele_type}{suffix}"
