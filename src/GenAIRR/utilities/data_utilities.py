import csv
import re
from collections import defaultdict
from ..alleles.allele import VAllele, JAllele, DAllele,CAllele
from ..utilities import parse_fasta

def create_allele_dict(fasta):
    """
    Constructs a dictionary that maps immunoglobulin gene types (V, D, or J) to their corresponding alleles,
    based on a provided IMGT-formatted reference FASTA file.

    Parameters:
    - fasta (file): A FASTA file that includes reference alleles for V, D, or J immunoglobulin genes,
      with sequences gapped as per IMGT guidelines. The headers of these sequences should comply with
      IMGT formatting standards.

    Returns:
    - allele_dict (dict): Returns a dictionary where keys are the gene types (V, D, or J) and the values
      are lists containing instances of the Allele class for each allele. This dictionary includes only
      alleles identified as functional or as open reading frames, indicated by "F" or "ORF" in their
      sequence headers, formatted as {Gene: [Allele, Allele, ...]}.
    """

    allele_dict = defaultdict(list)
    allele_class_map = {'V':VAllele,'D':DAllele,'J':JAllele,'C':CAllele}
    with open(fasta) as f:
        for header, seq in parse_fasta(f):
            
            header = header.replace(">","") # clean the fasta sequence tag
            seq = seq.lower() # make sure the sequences is in lower cases
            
            allele = header.split("|") # leave the IMGT option, but allow OGRDB reference
            if len(allele) == 1:
                allele = allele[0]
            else:
                allele = allele[1]
                
            segment = None

            if "-" in allele:
                family = allele.split("-")[0]
            else:
                family = allele.split("*")[0]
            if "D" in family:
                seq = seq.replace(".", "").lower()
                segment = 'D'
            if "V" in family:
                cys = seq[309:312]
                if cys != "tgt" and cys != "tgc":  # if allele doesn't have expected V anchor then don't use it
                    continue
                segment = 'V'
            if "J" in family:  # if allele doesn't have expected J anchor then don't use it
                anchor = None
                for frame in range(3):
                    motif = re.compile('(ttt|ttc|tgg)(ggt|ggc|gga|ggg)[a-z]{3}(ggt|ggc|gga|ggg)')
                    match = motif.search(seq[frame:])
                    if match:
                        if match.span()[0] % 3 == 0:
                            anchor = match.span()[0] + frame - 1
                if anchor is None:
                    continue
                segment = 'J'
            if "C" in family:
                segment = 'C'

            gene = allele.split("*")[0]
            
            coding = header.split("|") # extract the coding tag from IMGT references, OGRDB only include functional sequences 
            if len(coding) == 1:
                coding = "F"
            else:
                coding = coding[3]  
                
            ungapped_length = len(seq.replace(".","")) # get the length of the sequence
            
            if "partial" not in header and (coding == "F" or coding == "ORF"):
                if segment is not None:
                    allele_dict[gene].append(allele_class_map[segment](allele, seq, ungapped_length))

    return dict(allele_dict)