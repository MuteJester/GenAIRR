from collections import defaultdict
from ..utilities.data_utilities import create_allele_dict
from ..alleles.allele import VAllele


def _require_scipy():
    """Lazy import for scipy (only needed during DataConfig creation)."""
    try:
        from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
        from scipy.spatial.distance import squareform
        return dendrogram, linkage, fcluster, squareform
    except ImportError:
        raise ImportError(
            "scipy is required for allele clustering (DataConfig creation). "
            "Install it with: pip install scipy"
        )


def _require_numpy():
    """Lazy import for numpy (only needed during DataConfig creation)."""
    try:
        import numpy as np
        return np
    except ImportError:
        raise ImportError(
            "numpy is required for allele clustering (DataConfig creation). "
            "Install it with: pip install numpy"
        )

def allele_diff(reference_allele, sample_allele, position_threshold=0, snps=True):
    """
        Identifies differences between two alleles, optionally considering only positions beyond a threshold and reporting single nucleotide polymorphisms (SNPs).

        Args:
            reference_allele (str): The sequence of the reference allele.
            sample_allele (str): The sequence of the sample allele to compare against the reference.
            position_threshold (int, optional): The position from which to start comparing the sequences. Defaults to 0, which means comparison starts from the beginning.
            snps (bool, optional): If True, differences are reported as SNPs in the format 'RefPositionSample'. If False, only the positions of differences are reported. Defaults to True.

        Returns:
            list: A list of SNPs or position indices where the two alleles differ, considering the specified threshold and SNP reporting preference.
        """
    germs = [reference_allele, sample_allele]
    max_length = max(len(germ) for germ in germs)
    for i in range(len(germs)):
        germs[i] += '.' * (max_length - len(germs[i]))
    def setdiff_mat(x):
        unique_chars = set(x)
        filter_chars = {'.', 'N', '-'}
        return len(unique_chars - filter_chars)
    idx_strings = []
    for i in range(max_length):
        column_chars = [germ[i] for germ in germs]
        diff_count = setdiff_mat(column_chars)
        if diff_count > 1 and i >= (position_threshold-1):
            if snps:
                concatenated_str = column_chars[0].upper() + str(i+1) + column_chars[1].upper()
                idx_strings.append(concatenated_str)
            else:
                idx_strings.append(i+1)
    return idx_strings

def hamming_distance(s1, s2):
    """
        Calculates the Hamming distance between two sequences, ignoring positions with 'N' and gaps represented by '.'.

        Args:
            s1 (str): The first sequence for comparison.
            s2 (str): The second sequence for comparison.

        Returns:
            float: The Hamming distance between the two sequences, normalized by the number of comparable positions (excluding gaps).
    """
    np = _require_numpy()
    s1 = np.array(list(s1))
    s2 = np.array(list(s2))
    mismatches = np.sum((s1 != 'N') & (s2 != 'N') & (s1 != s2))
    gapGapMatches = np.sum((s1 == '.') & (s2 == '.'))
    count = len(s1) - gapGapMatches
    distance = mismatches / count if count > 0 else 0
    return distance

def asc_distance(germline_set, trim_3prime_side=318):
    """
        Calculates the pairwise Hamming distances between alleles in a set, optionally trimming sequences to a specified length before comparison.

        Args:
            germline_set (dict): A dictionary of allele sequences with allele names as keys.
            trim_3prime_side (int, optional): The position up to which sequences should be considered. If None, the entire length is used. Defaults to 318.

        Returns:
            ndarray: A square matrix of Hamming distances between all pairs of alleles in the set.
    """
    if trim_3prime_side is not None:
        germline_set = {allele: seq[:trim_3prime_side] for allele, seq in germline_set.items()}
    max_length = max(len(s) for s in germline_set.values())
    germline_set_padded = {allele: sequence.ljust(max_length, 'N') for allele, sequence in germline_set.items()}

    np = _require_numpy()
    sequences_array = np.array([list(seq) for seq in germline_set_padded.values()])

    num_sequences = len(sequences_array)
    germline_distance_matrix = np.zeros((num_sequences, num_sequences))

    for i in range(num_sequences):
        s1 = sequences_array[i]
        for j in range(num_sequences):
            s2 = sequences_array[j]
            mismatches = np.sum((s1 != 'N') & (s2 != 'N') & (s1 != s2))
            gapGapMatches = np.sum((s1 == '.') & (s2 == '.'))
            count = len(s1) - gapGapMatches
            germline_distance_matrix[i, j] = mismatches / count if count > 0 else 0

    return germline_distance_matrix


def _factorize(arr):
    """Assign sequential integer labels to values in order of first appearance.

    Returns a list of 1-based integer labels (same length as input).
    """
    seen = {}
    result = []
    counter = 0
    for val in arr:
        if val not in seen:
            counter += 1
            seen[val] = counter
        result.append(seen[val])
    return result


def asc_clust(germline_distance, germline_set, family_threshold=75, allele_cluster_threshold=95, trim_3prime_side=318):
    """
        Performs hierarchical clustering on alleles based on their Hamming distances and classifies them into families and allele clusters.

        Args:
            germline_distance (ndarray): A square matrix of pairwise Hamming distances between alleles.
            germline_set (dict): A dictionary of allele sequences with allele names as keys.
            family_threshold (int, optional): The distance threshold for classifying alleles into families. Defaults to 75.
            allele_cluster_threshold (int, optional): The distance threshold for classifying alleles within families into clusters. Defaults to 95.
            trim_3prime_side (int, optional): The position up to which sequences should be considered for clustering. Defaults to 318.

        Returns:
            list[dict]: A list of dicts with keys 'Family', 'Allele_Cluster', 'Allele',
                        'duplicated_allele', 'diff_pos_past_trim', 'new_allele'.
    """
    dendrogram, linkage, fcluster, squareform = _require_scipy()
    condensed_distance = squareform(germline_distance)
    germline_cluster = linkage(condensed_distance, method='complete')
    dendro = dendrogram(germline_cluster, labels=list(germline_set.keys()), no_plot=True)
    labels_order = dendro['ivl']
    leaf_order = dendro['leaves']

    labels = list(germline_set.keys())
    segment = labels[0][0:4]

    thresh = [1-family_threshold/100, 1-allele_cluster_threshold/100]
    family_cluster = fcluster(germline_cluster, t=thresh[0], criterion='distance')
    family_clusters_in_order = family_cluster[leaf_order]
    family_clusters_renumbered = _factorize(family_clusters_in_order)

    allele_cluster = fcluster(germline_cluster, t=thresh[1], criterion='distance')
    allele_clusters_in_order = allele_cluster[leaf_order]
    allele_clusters_renumbered = _factorize(allele_clusters_in_order)

    # Build table as list of dicts
    table = []
    for i in range(len(labels_order)):
        table.append({
            'Family': family_clusters_renumbered[i],
            'Allele_Cluster': allele_clusters_renumbered[i],
            'Allele': labels_order[i],
            'duplicated_allele': None,
            'diff_pos_past_trim': None,
        })

    # Build name->index lookup for fast access
    allele_to_idx = {row['Allele']: i for i, row in enumerate(table)}

    for i, row in enumerate(germline_cluster):
        cluster_1, cluster_2, distance, _ = row
        if distance == 0:
            # check which allele is longer, this should be the reference allele
            if len(germline_set[labels[int(cluster_1)]]) > len(germline_set[labels[int(cluster_2)]]):
                cluster_1_label = labels[int(cluster_1)]
                cluster_2_label = labels[int(cluster_2)]
            else:
                cluster_2_label = labels[int(cluster_1)]
                cluster_1_label = labels[int(cluster_2)]

            if trim_3prime_side is not None:
                snps = allele_diff(germline_set[cluster_1_label], germline_set[cluster_2_label], trim_3prime_side)
                if snps:
                    idx2 = allele_to_idx[cluster_2_label]
                    idx1 = allele_to_idx[cluster_1_label]
                    table[idx2]['diff_pos_past_trim'] = '_'.join(map(str, snps))
                    table[idx1]['duplicated_allele'] = cluster_1_label + "," + cluster_2_label
                    table[idx2]['duplicated_allele'] = cluster_1_label + "," + cluster_2_label
                else:
                    idx1 = allele_to_idx[cluster_1_label]
                    if table[idx1]['duplicated_allele'] is not None:
                        table[idx1]['duplicated_allele'] += cluster_2_label
                    else:
                        table[idx1]['duplicated_allele'] = cluster_2_label
                        idx2 = allele_to_idx[cluster_2_label]
                        table[idx2]['duplicated_allele'] = 'remove'
            else:
                idx1 = allele_to_idx[cluster_1_label]
                if table[idx1]['duplicated_allele'] is not None:
                    table[idx1]['duplicated_allele'] += cluster_2_label
                else:
                    table[idx1]['duplicated_allele'] = cluster_2_label
                    idx2 = allele_to_idx[cluster_2_label]
                    table[idx2]['duplicated_allele'] = 'remove'

    # Filter out 'remove' rows
    filtered = [row for row in table if row['duplicated_allele'] != 'remove']

    # Assign allele_idx per Allele_Cluster group (sequential 1-based)
    cluster_counters = defaultdict(int)
    for row in filtered:
        cluster_counters[row['Allele_Cluster']] += 1
        row['allele_idx'] = cluster_counters[row['Allele_Cluster']]

    # For rows sharing the same (Allele_Cluster, Family, duplicated_allele),
    # set allele_index to the minimum allele_idx in that group
    group_min = defaultdict(lambda: float('inf'))
    for row in filtered:
        key = (row['Allele_Cluster'], row['Family'], row['duplicated_allele'])
        group_min[key] = min(group_min[key], row['allele_idx'])

    for row in filtered:
        key = (row['Allele_Cluster'], row['Family'], row['duplicated_allele'])
        row['allele_index'] = group_min[key]

    # Re-rank allele_index within each Allele_Cluster to be sequential
    cluster_seen = defaultdict(list)
    for row in filtered:
        cluster_seen[row['Allele_Cluster']].append(row)

    for cluster_rows in cluster_seen.values():
        # Get unique allele_index values in order of appearance
        seen_indices = []
        for r in cluster_rows:
            if r['allele_index'] not in seen_indices:
                seen_indices.append(r['allele_index'])
        index_map = {old: new + 1 for new, old in enumerate(seen_indices)}
        for r in cluster_rows:
            r['allele_index'] = index_map[r['allele_index']]

    # Generate new allele names
    for row in filtered:
        idx = int(row['allele_index'])
        idx_str = str(idx).zfill(2) if idx >= 1 else str(idx)
        row['new_allele'] = f"{segment}F{row['Family']}-G{row['Allele_Cluster']}*{idx_str}"
        if row['diff_pos_past_trim'] is not None:
            row['new_allele'] += f"_{row['diff_pos_past_trim']}"

    return filtered


def asc_dict(allele_cluster_table, germline_set):
        """
        Converts an allele cluster table into a dictionary of VAllele objects organized by allele cluster sequence (ASC).

        Args:
            allele_cluster_table (list[dict]): A list of dicts containing allele family and cluster assignments.
            germline_set (dict): A dictionary of original allele sequences with allele names as keys.

        Returns:
            defaultdict: A dictionary where keys are ASC identifiers and values are lists of VAllele objects belonging to each ASC.
        """
        germline_set_asc = defaultdict(list)

        # T2-8: previously this passed no anchor_override, so VAllele's
        # legacy `_find_anchor` (the buggy `rfind+3` heuristic) was the
        # only anchor source for ASC-built configs — which silently
        # returned 2 on rfind miss. Now we resolve via the C-side
        # AnchorResolver and pass an explicit, validated anchor.
        from .._native._anchor import (
            LoadedAlleleRecord, Locus, Segment, FunctionalStatus,
            AnchorConfidence, resolve_anchor, locus_from_gene_name,
        )

        for row in allele_cluster_table:
            allele = row['Allele']
            new_allele = row['new_allele']
            asc = new_allele.split("*")[0]
            seq = germline_set[allele]
            ungapped = seq.replace(".", "")
            ungapped_length = len(ungapped)

            rec = LoadedAlleleRecord(
                name=new_allele, aliases=(),
                segment=Segment.V,
                locus=locus_from_gene_name(new_allele),
                species=None,
                sequence=ungapped, gapped_sequence=seq,
                gap_convention_imgt=True,
                functional_status=FunctionalStatus.UNKNOWN,
                explicit_anchor=-1,
                source="asc_dict",
            )
            result = resolve_anchor(rec, segment=Segment.V, locus=rec.locus)
            anchor = (result.position
                      if result.confidence != AnchorConfidence.REJECTED
                      else 0)  # 0 sentinel preserves legacy "anchorless" behavior

            germline_set_asc[asc].append(
                VAllele(new_allele, seq, ungapped_length,
                        anchor_override=anchor)
            )

        return germline_set_asc

def create_asc_germline_set(user_reference, segment="V", trim_3prime_side=318, family_threshold=75, allele_cluster_threshold=95):
    """
        Creates an ASC (Allele Sequence Cluster) germline set from a user-provided reference, applying trimming and clustering thresholds.

        Args:
            user_reference (dict or path): A dictionary of allele sequences or a path to a fasta file containing alleles.
            segment (str, optional): The segment type to consider (e.g., "V"). Defaults to "V".
            trim_3prime_side (int, optional): The position up to which sequences should be considered for clustering. Defaults to 318.
            family_threshold (int, optional): The distance threshold for classifying alleles into families. Defaults to 75.
            allele_cluster_threshold (int, optional): The distance threshold for classifying alleles within families into clusters. Defaults to 95.

        Returns:
            tuple: A tuple containing the ASC germline set dictionary and the allele cluster table.
        """
    germline_set_user = defaultdict(dict)
    for key, value in create_allele_dict(user_reference).items():
        for val in value:
            if segment in val.name:
                germline_set_user[val.name] = val.gapped_seq

    germline_distance = asc_distance(germline_set_user)

    asc_table = asc_clust(germline_distance, germline_set_user, family_threshold, allele_cluster_threshold, trim_3prime_side)

    germline_set_asc = asc_dict(asc_table, germline_set_user)

    return germline_set_asc, asc_table

def test_asc(fasta):
    """
        Tests the ASC germline set creation process using a fasta file or dictionary of alleles.

        Args:
            fasta (dict or path): A dictionary of allele sequences or a path to a fasta file containing alleles.

        Returns:
            defaultdict: A dictionary representing the ASC germline set.
    """
    germline_set = create_allele_dict(fasta)
    germline_set_asc = create_asc_germline_set(germline_set)
    return germline_set_asc

def profile(fasta):
    """
       Profiles the performance of the ASC germline set creation process using a fasta file or dictionary of alleles.

       Args:
           fasta (dict or path): A dictionary of allele sequences or a path to a fasta file containing alleles.
   """
    from cProfile import Profile
    germline_set = create_allele_dict(fasta)
    cProfile.run('create_asc_germline_set(germline_set)')
