import random
from collections import defaultdict

STOP_CODONS = {"TAG", "TAA", "TGA"}

def weighted_choice(choices):
    """
        Selects a key from a dictionary of choices based on their associated weights.

        Args:
            choices (dict): A dictionary where keys are the choices to select from and values are the weights of each choice.

        Returns:
            Any: A randomly selected key from the dictionary, where the probability of each key is weighted according to its associated value.
    """
    return random.choices(population=list(choices.keys()), weights=choices.values(), k=1)[0]


def weighted_choice_zero_break(choices):
    """
    Selects a key from a dictionary of choices based on their associated weights, returning zero if the choices dictionary is empty.

    Args:
        choices (dict): A dictionary where keys are the choices to select from and values are the weights of each choice.

    Returns:
        Any: A randomly selected key from the dictionary, where the probability of each key is weighted according to its associated value, or zero if the dictionary is empty.
    """
    if len(choices) < 1:
        return 0
    else:
        return weighted_choice(choices)


def translate(seq):
    """
    Translates a nucleotide sequence into a corresponding amino acid sequence using the standard genetic code.

    This function translates codons (three nucleotides) into amino acids and handles incomplete codons by adding a '.' symbol. Special symbols ('.', '-', 'N') in the sequence are handled as follows: '.' is replaced with '_', '-' with 'X', and 'N' with 'X'. Stop codons are represented by '*'.

    Args:
        seq (str): The nucleotide sequence to be translated. The sequence should be a string of 'A', 'T', 'C', 'G', 'N', '.', or '-' characters.

    Returns:
        str: The translated amino acid sequence. Amino acids are represented by their one-letter codes. Incomplete codons, gaps, and unknown nucleotides are represented by '.', 'X', and 'X', respectively.
    """

    table = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
        'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W'
    }
    protein = ""
    for i in range(0, len(seq), 3):
        codon = seq[i:i + 3].upper()
        if len(codon) < 3:
            protein += '.'
        elif "." in codon:
            protein += '_'
        elif "-" in codon:
            protein += 'X'
        elif "N" in codon:
            protein += 'X'
        else:
            protein += table[codon]
    return protein

def check_stops(seq,return_pos=False):
        """
        Checks the given sequence for the presence of stop codons.

        Args:
            seq (str): The nucleotide sequence to check for stop codons.

        Returns:
            bool: True if a stop codon is found, False otherwise.
        """
        for x in range(0, len(seq), 3):
            if seq[x:x + 3] in STOP_CODONS:
                if return_pos:
                    return True,x
                else:
                    return True
        if return_pos:
            return False,-1
        else:
            return False

def normalize_and_filter_convert_to_dict(obj):
    if isinstance(obj, defaultdict):
        nested_dict = {k: normalize_and_filter_convert_to_dict(v) for k, v in obj.items()}
        if all(isinstance(val, float) for val in nested_dict.values()):
            # Filter out keys that are not 'A', 'T', 'G', or 'C'
            filtered_dict = {k: v for k, v in nested_dict.items() if k in ['A', 'T', 'G', 'C']}
            total = sum(filtered_dict.values())
            return {k: v / total for k, v in filtered_dict.items()}
        return nested_dict
    return obj