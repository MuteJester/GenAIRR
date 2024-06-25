# TODO: Implement Proper Logic for this file

trim_map = dict()
for d_allele in d_alleles:
    trim_map[d_allele.name] = dict()
    for trim_5 in range(len(d_allele.ungapped_seq) + 1):
        for trim_3 in range(len(d_allele.ungapped_seq) - trim_5 + 1):
            # Correctly handle the trimming for d_allele
            trimmed = d_allele.ungapped_seq[trim_5:] if trim_5 > 0 else d_allele.ungapped_seq
            trimmed = trimmed[:-trim_3] if trim_3 > 0 else trimmed

            trim_map[d_allele.name][(trim_5, trim_3)] = []
            for d_c_allele in d_alleles:
                # Check if the trimmed sequence is a substring of the d_c_allele sequence
                if trimmed in d_c_allele.ungapped_seq:
                    trim_map[d_allele.name][(trim_5, trim_3)].append(d_c_allele.name)


trim_map = dict()
for v_allele in v_dict.values():
    trim_map[v_allele.name] = dict()
    for trim_5 in range(len(v_allele.ungapped_seq) + 1):
        trimmed = v_allele.ungapped_seq[trim_5:] if trim_5 > 0 else v_allele.ungapped_seq
        trim_map[v_allele.name][trim_5] = []
        for v_c_allele in v_dict.values():
            # Check if the trimmed sequence is a substring of the v_c_allele sequence
            if trimmed in v_c_allele.ungapped_seq:
                trim_map[v_allele.name][trim_5].append(v_c_allele.name)
