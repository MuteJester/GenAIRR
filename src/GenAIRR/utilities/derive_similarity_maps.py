class AmbiguityAlleleTrimMapDeriver:
    """
     A class used to derive ambiguity maps for V, D, and J alleles based on sequence trimming.

     The class provides methods to generate maps that describe which alleles remain ambiguous
     (i.e., indistinguishable) after various degrees of trimming at the 5' and/or 3' ends.
     """
    def __init__(self):
        """
               Initializes the AmbiguityAlleleTrimMapDeriver object.

               This constructor does not require any parameters and does not perform any
               specific initialization tasks. It simply prepares an instance of the class.
           """
        pass

    def derive_v_trim_ambiguity_map(self, dataconfig):
        """
                Derives ambiguity maps for V alleles based on 5' and 3' trimming.

                This method generates two dictionaries: one for 5' trimming and one for 3' trimming.
                Each dictionary maps V allele names to another dictionary that maps the number of
                trimmed nucleotides to a list of allele names that contain the trimmed sequence.

                Parameters
                ----------
                dataconfig : object
                    A configuration object containing information about the V alleles.
                    The object should have a structure where `dataconfig.v_alleles` is a dictionary
                    with allele names as keys and objects representing those alleles as values.

                Returns
                -------
                v_5_trim_map : dict
                    A dictionary where the keys are V allele names, and the values are dictionaries
                    mapping the number of 5' trimmed nucleotides to a list of V allele names that
                    contain the trimmed sequence.

                v_3_trim_map : dict
                    A dictionary where the keys are V allele names, and the values are dictionaries
                    mapping the number of 3' trimmed nucleotides to a list of V allele names that
                    contain the trimmed sequence.

                Example
                -------
                >>> v_5_trim_map, v_3_trim_map = deriver.derive_v_trim_ambiguity_map(dataconfig)
                >>> print(v_5_trim_map['IGHV1-2'][2])
                ['IGHV1-2', 'IGHV1-3']
            """
        v_dict = {i.name:i for j in dataconfig.v_alleles for i in dataconfig.v_alleles[j]}

        v_5_trim_map = dict()
        for v_allele in v_dict.values():
            v_5_trim_map[v_allele.name] = dict()
            for trim_5 in range(len(v_allele.ungapped_seq) + 1):
                trimmed = v_allele.ungapped_seq[trim_5:] if trim_5 > 0 else v_allele.ungapped_seq
                v_5_trim_map[v_allele.name][trim_5] = []
                for v_c_allele in v_dict.values():
                    # Check if the trimmed sequence is a substring of the v_c_allele sequence
                    if trimmed in v_c_allele.ungapped_seq:
                        v_5_trim_map[v_allele.name][trim_5].append(v_c_allele.name)

        v_3_trim_map = dict()
        for v_allele in v_dict.values():
            v_3_trim_map[v_allele.name] = dict()
            for trim_3 in range(len(v_allele.ungapped_seq) + 1):
                trimmed = v_allele.ungapped_seq[:-trim_3] if trim_3 > 0 else v_allele.ungapped_seq
                v_3_trim_map[v_allele.name][trim_3] = []
                for v_c_allele in v_dict.values():
                    # Check if the trimmed sequence is a substring of the v_c_allele sequence
                    if trimmed in v_c_allele.ungapped_seq:
                        v_3_trim_map[v_allele.name][trim_3].append(v_c_allele.name)

        return v_5_trim_map, v_3_trim_map

    def derive_j_trim_ambiguity_map(self, dataconfig):
        """
                Derives ambiguity maps for J alleles based on 5' and 3' trimming.

                This method generates two dictionaries: one for 5' trimming and one for 3' trimming.
                Each dictionary maps J allele names to another dictionary that maps the number of
                trimmed nucleotides to a list of allele names that contain the trimmed sequence.

                Parameters
                ----------
                dataconfig : object
                    A configuration object containing information about the J alleles.
                    The object should have a structure where `dataconfig.j_alleles` is a dictionary
                    with allele names as keys and objects representing those alleles as values.

                Returns
                -------
                j_5_trim_map : dict
                    A dictionary where the keys are J allele names, and the values are dictionaries
                    mapping the number of 5' trimmed nucleotides to a list of J allele names that
                    contain the trimmed sequence.

                j_3_trim_map : dict
                    A dictionary where the keys are J allele names, and the values are dictionaries
                    mapping the number of 3' trimmed nucleotides to a list of J allele names that
                    contain the trimmed sequence.

                Example
                -------
                >>> j_5_trim_map, j_3_trim_map = deriver.derive_j_trim_ambiguity_map(dataconfig)
                >>> print(j_5_trim_map['IGHJ4'][3])
                ['IGHJ4', 'IGHJ6']
            """
        j_dict = {i.name: i for j in dataconfig.j_alleles for i in dataconfig.j_alleles[j]}

        j_5_trim_map = dict()
        for j_allele in j_dict.values():
            j_5_trim_map[j_allele.name] = dict()
            for trim_5 in range(len(j_allele.ungapped_seq) + 1):
                trimmed = j_allele.ungapped_seq[trim_5:] if trim_5 > 0 else j_allele.ungapped_seq
                j_5_trim_map[j_allele.name][trim_5] = []
                for j_c_allele in j_dict.values():
                    # Check if the trimmed sequence is a substring of the j_c_allele sequence
                    if trimmed in j_c_allele.ungapped_seq:
                        j_5_trim_map[j_allele.name][trim_5].append(j_c_allele.name)

        j_3_trim_map = dict()
        for j_allele in j_dict.values():
            j_3_trim_map[j_allele.name] = dict()
            for trim_3 in range(len(j_allele.ungapped_seq) + 1):
                trimmed = j_allele.ungapped_seq[:-trim_3] if trim_3 > 0 else j_allele.ungapped_seq
                j_3_trim_map[j_allele.name][trim_3] = []
                for j_c_allele in j_dict.values():
                    # Check if the trimmed sequence is a substring of the j_c_allele sequence
                    if trimmed in j_c_allele.ungapped_seq:
                        j_3_trim_map[j_allele.name][trim_3].append(j_c_allele.name)

        return j_5_trim_map, j_3_trim_map

    def derive_d_trim_ambiguity_map(self, dataconfig):
        """
                Derives ambiguity maps for D alleles based on both 5' and 3' trimming.

                This method generates a dictionary where each D allele name maps to another dictionary
                that maps a tuple of (5' trimmed nucleotides, 3' trimmed nucleotides) to a list of
                allele names that contain the trimmed sequence.

                Parameters
                ----------
                dataconfig : object
                    A configuration object containing information about the D alleles.
                    The object should have a structure where `dataconfig.d_alleles` is a dictionary
                    with allele names as keys and objects representing those alleles as values.

                Returns
                -------
                d_5_3_trim_map : dict
                    A dictionary where the keys are D allele names, and the values are dictionaries
                    mapping a tuple of (5' trimmed nucleotides, 3' trimmed nucleotides) to a list
                    of D allele names that contain the trimmed sequence.

                Example
                -------
                >>> d_5_3_trim_map = deriver.derive_d_trim_ambiguity_map(dataconfig)
                >>> print(d_5_3_trim_map['IGHD1-1'][(2, 1)])
                ['IGHD1-1', 'IGHD2-2']
                """
        d_alleles = [i for j in dataconfig.d_alleles for i in dataconfig.d_alleles[j]]
        d_5_3_trim_map = dict()
        for d_allele in d_alleles:
            d_5_3_trim_map[d_allele.name] = dict()
            for trim_5 in range(len(d_allele.ungapped_seq) + 1):
                for trim_3 in range(len(d_allele.ungapped_seq) - trim_5 + 1):
                    # Correctly handle the trimming for d_allele
                    trimmed = d_allele.ungapped_seq[trim_5:] if trim_5 > 0 else d_allele.ungapped_seq
                    trimmed = trimmed[:-trim_3] if trim_3 > 0 else trimmed

                    d_5_3_trim_map[d_allele.name][(trim_5, trim_3)] = []
                    for d_c_allele in d_alleles:
                        # Check if the trimmed sequence is a substring of the d_c_allele sequence
                        if trimmed in d_c_allele.ungapped_seq:
                            d_5_3_trim_map[d_allele.name][(trim_5, trim_3)].append(d_c_allele.name)
        return d_5_3_trim_map

