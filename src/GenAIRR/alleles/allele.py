"""
Note:
    Adapted from https://github.com/Cowanlab/airrship
"""
import re
import random
from abc import ABC, abstractmethod
from ..utilities import weighted_choice_zero_break
from enum import Enum, auto


class AlleleTypes(Enum):
    """Enumeration for allele types.

        Attributes:
            V: Represents the V allele type.
            D: Represents the D allele type.
            J: Represents the J allele type.
            C: Represents the Constant allele type.
        """
    V = auto()
    D = auto()
    J = auto()
    C = auto() # constant


class Allele(ABC):
    """Abstract base class for alleles.

        This class provides a template for alleles, defining common properties and methods
        that all allele types share.

        Attributes:
            name (str): The name of the allele.
            length (int): The length of the gapped sequence.
            gapped_seq (str): The gapped nucleotide sequence of the allele.
            ungapped_seq (str): The ungapped nucleotide sequence derived from the gapped sequence.
            ungapped_len (int): The length of the ungapped sequence.
            family (str): The gene family derived from the allele name.
            gene (str): The gene name derived from the allele name.
            anchor (int): The anchor position within the sequence, specific to allele type.
        """

    def __init__(self, name, gapped_sequence, length):
        """Initializes an Allele instance.

            Args:
                name (str): The name of the allele.
                gapped_sequence (str): The gapped nucleotide sequence.
                length (int): The length of the gapped sequence.
        """
        self.name = name
        self.length = int(length)
        self.gapped_seq = gapped_sequence
        self.ungapped_seq = gapped_sequence.replace(".", "")
        self.ungapped_len = len(self.ungapped_seq)
        self.family = self.name.split("-")[0] if "-" in self.name else self.name.split("*")[0]
        self.gene = self.name.split("*")[0]
        self.anchor = None
        self._find_anchor()

    @abstractmethod
    def _find_anchor(self):
        """Finds the anchor position within the allele sequence.

           This method must be implemented by subclasses based on allele-specific criteria.
        """
        pass

    def __repr__(self):
        """Represents the Allele instance as a string.

               Returns:
                   str: A string representation of the Allele instance.
       """
        return (f"{self.__class__.__name__}(name='{self.name}', "
                f"sequence='...{self.ungapped_seq[-10:]}', "
                f"length={self.length}, "
                f"ungapped_len={self.ungapped_len}, "
                f"family='{self.family}', "
                f"gene='{self.gene}')")

    @abstractmethod
    def _get_trim_length(self, trim_dicts):
        """Generates the length by which the allele's ungapped sequence should be trimmed.

        This method must be implemented by subclasses to provide allele-specific trimming logic.

        Args:
            trim_dicts (dict): A dictionary of dictionaries of trimming length
                proportions by gene family for each segment (V, D, or J).

        Returns:
            tuple: A tuple of integers representing the trim lengths at the 5' and 3' ends, respectively.
        """

        pass

    @abstractmethod
    def get_trimmed(self, trim_dict):
        """Returns the trimmed sequence along with the trim lengths.

                This method must be implemented by subclasses to apply allele-specific trimming and return the result.

                Args:
                    trim_dict (dict): A dictionary specifying trimming lengths for the allele.

                Returns:
                    tuple: The trimmed sequence, and the trim lengths at the 5' and 3' ends, respectively.
        """
        pass


class VAllele(Allele):
    """Represents a V allele with specific trimming and analysis behaviors.

        This class extends the Allele base class, providing implementations for V allele-specific methods.

        Attributes:
            type (AlleleTypes): The type of the allele, fixed to AlleleTypes.V for V alleles.
    """
    type = AlleleTypes.V

    def _find_anchor(self):
        """Finds the anchor position within a V allele sequence.

            The anchor is identified based on the presence of a specific motif within the gapped sequence.
        """
        cys_wider = self.gapped_seq[306:315]
        self.anchor = self.ungapped_seq.rfind(cys_wider) + 3

    def _get_trim_length(self, trim_dicts):
        """Determines the trim lengths for the V allele's sequence based on provided trimming dictionaries.

        Args:
            trim_dicts (dict): A dictionary specifying trimming lengths for different allele types.

        Returns:
            tuple: A tuple of integers representing the trim lengths at the 5' and 3' ends, respectively.
        """
        trim_3 = 0  # set to 0 - J will never be trimmed at 3'
        trim_5 = 0  # set to 0 - V will never be trimmed at 5'

        trim_3_dict = trim_dicts["V_3"]
        # choose trim length/prob dict by gene family
        if self.family in trim_3_dict:
            prob_dict = trim_3_dict[self.family][self.gene]
        else:
            random_fam = random.choice(list(trim_3_dict.values()))
            prob_dict = random.choice(list(random_fam.values()))


        # prevent entire allele or anchor from being removed
        valid_trim_amounts = filter(lambda amount: (amount < self.length) or \
                                                   (amount < (self.length - self.anchor - 1)), prob_dict)

        prob_dict = {amount: prob_dict[amount] for amount in valid_trim_amounts}

        trim_3 = weighted_choice_zero_break(prob_dict)

        return int(trim_5), int(trim_3)  # make sure type is not float

    def get_trimmed(self, trim_dict):
        """Returns the trimmed sequence of the V allele along with the 5' and 3' trim lengths.

        Args:
            trim_dict (dict): A dictionary specifying trimming lengths for the V allele.

        Returns:
            tuple: The trimmed sequence, and the trim lengths at the 5' and 3' ends, respectively.
        """
        sequence = self.ungapped_seq
        trim_5, trim_3 = self._get_trim_length(trim_dict)
        trimmed_seq = sequence[:-trim_3 if trim_3 > 0 else None]
        return trimmed_seq, trim_5, trim_3


class DAllele(Allele):
    """Represents a D allele with specific trimming and analysis behaviors.

        This class extends the Allele base class, providing implementations for D allele-specific methods.

        Attributes:
            type (AlleleTypes): The type of the allele, fixed to AlleleTypes.D for D alleles.
    """
    type = AlleleTypes.D

    def _find_anchor(self):
        """Finds the anchor position within a D allele sequence.

        D alleles might not have a specific anchor finding logic, so this method can be
                implemented based on D allele-specific criteria if needed.
        """
        pass

    def _get_trim_length(self, trim_dicts):
        """Determines the trim lengths for the D allele's sequence based on provided trimming dictionaries.

           Args:
               trim_dicts (dict): A dictionary specifying trimming lengths for different allele types.

           Returns:
               tuple: A tuple of integers representing the trim lengths at the 5' and 3' ends, respectively.
       """
        trim_3 = 0  # set to 0 - J will never be trimmed at 3'
        trim_5 = 0  # set to 0 - V will never be trimmed at 5'

        trim_5_dict = trim_dicts["D_5"]
        if self.family in trim_5_dict:
            prob_5_dict = trim_5_dict[self.family][self.gene]
        else:
            random_fam = random.choice(list(trim_5_dict.values()))
            prob_5_dict = random.choice(list(random_fam.values()))

        valid_d5_trim_amounts = filter(lambda amount: amount + trim_5 < self.length, prob_5_dict)
        valid_d5_trim_amounts = {amount: prob_5_dict[amount] for amount in valid_d5_trim_amounts}
        trim_5 = weighted_choice_zero_break(valid_d5_trim_amounts)

        trim_3_dict = trim_dicts["D_3"]
        if self.family in trim_3_dict:
            prob_3_dict = trim_3_dict[self.family][self.gene]
        else:
            random_fam = random.choice(list(trim_3_dict.values()))
            prob_3_dict = random.choice(list(random_fam.values()))



        valid_d3_trim_amounts = filter(lambda amount: amount + trim_5 < self.length, prob_3_dict)

        prob_3_dict = {amount: prob_3_dict[amount] for amount in valid_d3_trim_amounts}

        trim_3 = weighted_choice_zero_break(prob_3_dict)

        return int(trim_5), int(trim_3)  # make sure type is not float

    def get_trimmed(self, trim_dict):
        """Returns the trimmed sequence of the D allele along with the 5' and 3' trim lengths.

           Args:
               trim_dict (dict): A dictionary specifying trimming lengths for the D allele.

           Returns:
               tuple: The trimmed sequence, and the trim lengths at the 5' and 3' ends, respectively.
       """
        sequence = self.ungapped_seq
        trim_5, trim_3 = self._get_trim_length(trim_dict)
        trimmed_seq = sequence[trim_5:len(sequence) - trim_3 if trim_3 > 0 else None]
        return trimmed_seq, trim_5, trim_3


class JAllele(Allele):
    """Represents a J allele with specific trimming and analysis behaviors.

    This class extends the Allele base class, providing implementations for J allele-specific methods.

    Attributes:
        type (AlleleTypes): The type of the allele, fixed to AlleleTypes.J for J alleles.
    """
    type = AlleleTypes.J

    def _find_anchor(self):
        """Finds the anchor position within a J allele sequence.

            The anchor is identified based on a specific motif within the ungapped sequence,
            which is characteristic of J alleles.

            The motif is defined by a regular expression that searches for specific nucleotide patterns.
        """
        # motif = re.compile('(ttt|ttc|tgg)(ggt|ggc|gga|ggg)[a-z]{3}(ggt|ggc|gga|ggg)')
        # match = motif.search(self.ungapped_seq)
        # self.anchor = match.span()[0] if match else None
        ### ---------------------------------------------------------------------- ###
        ### testing the influence of including the correct frame for the J anchor. ###
        ### same concept but we will look for the motif based on the frame. Span%3 ###
        self.anchor = None
        self.frame = None
        for frame in range(3):
            motif = re.compile('(ttt|ttc|tgg)(ggt|ggc|gga|ggg)[a-z]{3}(ggt|ggc|gga|ggg)')
            match = motif.search(self.ungapped_seq[frame:])
            if match:
                if match.span()[0] % 3 == 0:
                    self.anchor = match.span()[0] + frame # correct for the end of the cdr3
                    self.frame = frame # retain the frame of the J


    def _get_trim_length(self, trim_dicts):
        """Determines the trim lengths for the J allele's sequence based on provided trimming dictionaries.

        Args:
            trim_dicts (dict): A dictionary specifying trimming lengths for different allele types.

        Returns:
            tuple: A tuple of integers representing the trim lengths at the 5' end, as the 3' end is not trimmed for J alleles.
        """
        trim_3 = 0  # set to 0 - J will never be trimmed at 3'
        trim_5 = 0  # set to 0 - V will never be trimmed at 5'

        trim_5_dict = trim_dicts["J_5"]
        if self.family in trim_5_dict:
            prob_dict = trim_5_dict[self.family][self.gene]
        else:
            random_fam = random.choice(list(trim_5_dict.values()))
            prob_dict = random.choice(list(random_fam.values()))

        valid_5_trims = filter(lambda t5: (t5 < self.length) or (t5 < self.anchor), prob_dict)
        prob_dict = {amount: prob_dict[amount] for amount in valid_5_trims}
        trim_5 = weighted_choice_zero_break(prob_dict)
        return int(trim_5), int(trim_3)  # make sure type is not float

    def get_trimmed(self, trim_dict):
        """Returns the trimmed sequence of the J allele along with the 5' trim length.

        Args:
            trim_dict (dict): A dictionary specifying trimming lengths for the J allele.

        Returns:
            tuple: The trimmed sequence and the trim length at the 5' end.
        """
        trim_5, trim_3 = self._get_trim_length(trim_dict)
        sequence = self.ungapped_seq
        return sequence[trim_5:], trim_5, trim_3


class CAllele(Allele):
    """Represents a Constant allele with specific trimming and analysis behaviors.

    This class extends the Allele base class, providing implementations for C allele-specific methods.

    Attributes:
        type (AlleleTypes): The type of the allele, fixed to AlleleTypes.C for Constant alleles.
    """
    type = AlleleTypes.C

    def _find_anchor(self):
        """Finds the anchor position within the allele sequence.

           This method must be implemented by subclasses based on allele-specific criteria.
        """
        pass

    def _get_trim_length(self, trim_dicts):
        """Determines the trim lengths for the J allele's sequence based on provided trimming dictionaries.

        Args:
            trim_dicts (dict): A dictionary specifying trimming lengths for different allele types.

        Returns:
            tuple: A tuple of integers representing the trim lengths at the 5' end, as the 3' end is not trimmed for J alleles.
        """
        trim_3 = 0  # set to 0 - J will never be trimmed at 3'
        trim_5 = 0  # set to 0 - V will never be trimmed at 5'

        trim_3_dict = trim_dicts["C_3"]
        if self.family in trim_3_dict:
            prob_dict = trim_3_dict[self.family][self.gene]
        else:
            random_fam = random.choice(list(trim_3_dict.values()))
            prob_dict = random.choice(list(random_fam.values()))

        valid_3_trims = filter(lambda amount: amount < self.length, prob_dict)
        prob_dict = {amount: prob_dict[amount] for amount in valid_3_trims}
        trim_3 = weighted_choice_zero_break(prob_dict)
        return int(trim_5), int(trim_3)  # make sure type is not float

    def get_trimmed(self, trim_dict):
        """Returns the trimmed sequence of the J allele along with the 5' trim length.

        Args:
            trim_dict (dict): A dictionary specifying trimming lengths for the J allele.

        Returns:
            tuple: The trimmed sequence and the trim length at the 5' end.
        """
        trim_5, trim_3 = self._get_trim_length(trim_dict)
        sequence = self.ungapped_seq
        trimmed_seq = sequence[:-trim_3 if trim_3 > 0 else None]
        return trimmed_seq, trim_5, trim_3
