import pickle

class AlleleNComparer:
    """
        A class to compare alleles and find indistinguishable alleles based on specific positions.

        This class allows for the addition of alleles with their sequences, comparison of alleles to identify
        indistinguishable ones based on specified positions, and saving/loading of allele data.

        Attributes:
            alleles (dict): A dictionary storing the sequences of each allele, with allele IDs as keys and sequences as values.
        """
    def __init__(self):
        """
         Initializes the AlleleNComparer with an empty dictionary to store allele sequences.
       """
        self.alleles = {}  # Stores the sequences of each allele

    def add_allele(self, allele_id, sequence):
        """
            Adds an allele and its sequence to the alleles dictionary.

            Args:
                allele_id (str): The unique identifier for the allele.
                sequence (str): The nucleotide sequence of the allele.
        """
        self.alleles[allele_id] = sequence

    def find_indistinguishable_alleles(self, allele_id, n_positions):
        """
            Finds alleles that are indistinguishable from a specified allele, ignoring specified positions.

            Args:
                allele_id (str): The ID of the target allele to compare against.
                n_positions (list or set): Positions to be ignored during the comparison.

            Returns:
                set: A set of allele IDs that are indistinguishable from the specified allele, excluding the specified positions.
                str: An error message if the specified allele ID is not found in the alleles dictionary.
        """
        if allele_id not in self.alleles:
            return f"Allele {allele_id} not found."

        target_sequence = self.alleles[allele_id]
        indistinguishable_alleles = set()

        for other_allele_id, sequence in self.alleles.items():
            is_indistinguishable = True
            for pos in range(min(len(target_sequence), len(sequence))):
                if pos in n_positions:
                    continue
                if target_sequence[pos] != sequence[pos]:
                    is_indistinguishable = False
                    break

            if is_indistinguishable:
                indistinguishable_alleles.add(other_allele_id)

        return indistinguishable_alleles

    def save(self, filename):
        """
        Saves the alleles dictionary to a file using pickle for serialization.

        Args:
            filename (str): The name of the file where the alleles dictionary will be saved.
        """
        with open(filename, 'wb') as file:
            pickle.dump(self.alleles, file)

    def load(self, filename):
        """
        Loads an alleles dictionary from a file using pickle for deserialization.

        Args:
            filename (str): The name of the file from which the alleles dictionary will be loaded.

        Updates:
            The method updates the alleles dictionary with the contents loaded from the file.
        """
        with open(filename, 'rb') as file:
            self.alleles = pickle.load(file)

