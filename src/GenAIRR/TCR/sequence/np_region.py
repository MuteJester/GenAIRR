from ...utilities import weighted_choice


class NP_Region:
    """
    Represents an NP region in an immunoglobulin sequence, where the base
    transitions are determined by a first-order Markov chain with position-specific transition matrices.

    Attributes:
        transition_probs (dict): A dictionary where keys are position indices and values are
            transition matrices (dict) indicating the probability of transitioning from one base to another.
        first_base (str): The initial base of the NP region, used as the starting point for the Markov chain.
        length (int): The total length of the NP region to be generated.

    Args:
        transition_probs (dict): Position-specific transition probabilities for base changes.
        first_base (str): The starting base of the NP region sequence.
        length (int): The desired length of the NP region sequence.
    """

    def __init__(self, transition_probs, first_base, length):
        """Initializes an NP_Region instance with specified transition probabilities, initial base, and length."""

        self.transition_probs = transition_probs
        self.bases = ["A", "C", "G", "T"]
        self.first_base = first_base
        self.length = length

    def next_base(self, current_base, position):
        """
        Determines the next base in the NP region sequence based on the current base and position
        using the provided transition probabilities.

        Args:
            current_base (str): The current base in the sequence.
            position (int): The current position in the NP region, used to select the appropriate transition matrix.

        Returns:
            str: The next base in the sequence, chosen based on the transition probabilities.
        """

        next_base_options = list(self.transition_probs[position][current_base])
        base = weighted_choice(
            {next_base: self.transition_probs[position][current_base][next_base] for next_base in next_base_options})
        return base

    def validate_next_base(self,current_base,position):
        """
        Validates whether transition data is available for the current base at the given position,
        ensuring that the Markov chain can continue.

        Args:
            current_base (str): The current base in the sequence.
            position (int): The current position in the NP region.

        Returns:
            bool: True if transition data is available and the sequence generation can continue, False otherwise.
        """
        return True if current_base in self.transition_probs[position] else False


    def generate_np_seq(self):
        """
        Generates an NP region sequence using a first-order Markov chain based on the provided transition
        probabilities and initial base.

        Returns:
            str: The generated NP region sequence. If the sequence generation is halted due to lack of data,
            the sequence is returned as is and its length is updated to reflect the actual generated length.
        """
        sequence = ""
        current_base = self.first_base
        sequence += current_base
        for i in range(self.length - 1):

            if self.validate_next_base(current_base,i): # valid position
                next_base = self.next_base(current_base, i)
                sequence += next_base
                current_base = next_base
            else: #not way to continue halt and update metadata!
                self.length = len(sequence)
                return sequence.lower()

        return sequence.lower()

    @classmethod
    def create_np_region(cls, NP_lengths, NP_transitions, which_NP, first_base_dict):
        """
        Class method to generate an NP region sequence based on specified lengths, transitions,
        the type of NP region (NP1 or NP2), and initial base probabilities.

        Args:
            NP_lengths (dict): A dictionary mapping NP region types (e.g., 'NP1', 'NP2') to
                their length distributions.
            NP_transitions (dict): A dictionary mapping NP region types to their position-specific
                transition probabilities.
            which_NP (str): Specifies the type of NP region ('NP1' or 'NP2') for which the sequence is generated.
            first_base_dict (dict): A dictionary mapping NP region types to the distribution of possible
                initial bases.

        Returns:
            str: The generated NP region sequence. Returns an empty string if the specified length is zero.
        """
        length = weighted_choice(NP_lengths[which_NP])
        if length > 0:
            first_base = weighted_choice(first_base_dict[which_NP])
            np_region = cls(NP_transitions[which_NP], first_base, length)
            return np_region.generate_np_seq()
        return ""
