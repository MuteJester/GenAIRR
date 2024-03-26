import random
from ..mutation.mutation_model import MutationModel
import pickle


class Nucleotide:
    """Represents a single nucleotide in a DNA sequence.

   Attributes:
       value (str): The nucleotide character (A, T, C, G).
       adjacent (list): A list of references to adjacent `FiveMER` objects containing this nucleotide.

   Args:
       value (str): The nucleotide character.
   """
    def __init__(self, value):
        """Initialize a Nucleotide with a given value."""
        self.value = value
        self.adjacent = []  # References to adjacent nucleotides

    def update_value(self, new_value, update_callback=None):
        """Updates the nucleotide's value and optionally triggers a callback for each adjacent FiveMER.

        Args:
            new_value (str): The new nucleotide value.
            update_callback (function, optional): A callback function to be called for each adjacent FiveMER.
        """
        self.value = new_value
        if update_callback:
            for five_mer in self.adjacent:
                update_callback(five_mer)

    def __repr__(self):
        """String representation of the Nucleotide object."""
        return self.value


class FiveMER:
    """Represents a 5-mer, a sequence of 5 nucleotides, in a DNA sequence.

    Attributes:
        nucleotides (list): A list of `Nucleotide` objects that make up the 5-mer.
        sequence (str): The nucleotide sequence of the 5-mer.
        position (int): The position of the 5-mer within a larger sequence.
        likelihood (float): The likelihood or probability of this 5-mer being mutated.
        modified (bool): Indicates whether the 5-mer has been modified.

    Args:
        nucleotides (list): A list of `Nucleotide` objects.
    """
    def __init__(self, nucleotides):
        """Initialize a FiveMER with a list of Nucleotide objects."""
        self.nucleotides = nucleotides  # List of Nucleotide objects
        for nuc in self.nucleotides:
            nuc.adjacent.append(self)
        self.sequence = ''.join([nuc.value for nuc in nucleotides])
        self.position = None
        self.likelihood = 0
        self.modified = False

    def update_sequence(self, mutability=None):
        """Updates the 5-mer sequence and optionally its mutability likelihood.

        Args:
            mutability (dict, optional): A dictionary mapping 5-mer sequences to their mutability likelihoods.
        """
        self.sequence = ''.join([nuc.value for nuc in self.nucleotides])
        if mutability is not None:
            self.likelihood = self.likelihood if self.sequence not in mutability else mutability[self.sequence]

    def change_center(self, new_value, mutability=None):
        """Changes the central nucleotide of the 5-mer and updates sequences and likelihoods.

        Args:
            new_value (str): The new value for the central nucleotide.
            mutability (dict, optional): A dictionary mapping 5-mer sequences to their mutability likelihoods.
        """
        center_nucleotide = self.nucleotides[2]  # Assuming 0-indexed, 2 is the center
        # Pass the update callback to update_value
        center_nucleotide.update_value(new_value, lambda fm: fm.update_sequence(mutability))
        self.modified = True

    def __repr__(self):
        """String representation of the FiveMER object."""
        return ''.join([i.value for i in self.nucleotides])

    def __eq__(self, other):
        """Defines equality comparison for FiveMER objects with other FiveMER objects or strings."""
        if isinstance(other, FiveMER):
            return self.sequence == other.sequence
        elif isinstance(other, str):
            return self.sequence == other
        else:
            raise TypeError("Unsupported comparison between FiveMER and {}".format(type(other)))

    @staticmethod
    def create_five_mers(dna_sequence, mutability=None):
        """Creates a list of FiveMER objects from a DNA sequence.

        Args:
            dna_sequence (str): The DNA sequence from which to create 5-mers.
            mutability (dict, optional): A dictionary mapping 5-mer sequences to their mutability likelihoods.

        Returns:
            list: A list of FiveMER objects.
        """
        # Step 1: Pad the sequence
        padded_sequence = 'NN' + dna_sequence + 'NN'

        # Step 2: Create Nucleotide objects for the padded sequence
        nucleotides = [Nucleotide(nuc) for nuc in padded_sequence]
        # Step 3: Group nucleotides into FiveMER objects
        five_mers = []
        for i in range(len(dna_sequence)):  # Iterate based on the original sequence length
            five_mer_nucleotides = nucleotides[i:i + 5]
            five_mer = FiveMER(five_mer_nucleotides)
            five_mer.position = i  # Position of the original second nucleotide (central in unpadded)

            # if mutability map was supplied
            if mutability is not None:
                str_five_mer = five_mer.sequence
                five_mer.likelihood = 0 if str_five_mer not in mutability else mutability[str_five_mer]

            five_mers.append(five_mer)

        return five_mers

    @staticmethod
    def five_mers_to_dna(five_mers):
        """Converts a list of FiveMER objects back into a DNA sequence string.

        Args:
            five_mers (list): A list of FiveMER objects.

        Returns:
            str: The DNA sequence reconstructed from the FiveMER objects.
        """
        if not five_mers:
            return ""

        # Start with the central nucleotide of the first FiveMER (skipping initial padding)
        dna_sequence = five_mers[0].nucleotides[2].value

        # Iterate over the remaining FiveMERs and append only the last nucleotide of each
        for five_mer in five_mers[1:-1]:  # Skip the last FiveMER with padding
            dna_sequence += five_mer.nucleotides[2].value

        # Add the central nucleotide of the last FiveMER (which is the actual last nucleotide of the DNA)
        dna_sequence += five_mers[-1].nucleotides[2].value

        return dna_sequence


class S5F(MutationModel):
    """Implements the S5F mutation model, a specific model for simulating mutations in DNA sequences.

    This class extends `MutationModel` and provides an implementation for applying mutations based on the S5F model.

    Attributes:
        min_mutation_rate (float): The minimum mutation rate.
        max_mutation_rate (float): The maximum mutation rate.
        targeting (dict): A dictionary mapping targeting information.
        substitution (dict): A dictionary mapping substitution probabilities.
        mutability (dict): A dictionary mapping mutability scores.
        bases (set): A set of valid nucleotide bases.
        loaded_metadata (bool): Indicates whether the mutation metadata has been loaded.
        custom_model (str, optional): Path to a custom mutation model file.

    Args:
        min_mutation_rate (float): The minimum mutation rate.
        max_mutation_rate (float): The maximum mutation rate.
        custom_model (str, optional): Path to a custom mutation model file.
        productive (bool): Whether to ensure the sequence is productive (No stop codons and mutation in cdr3 anchors), defaulting to false.
    """
    def __init__(self, min_mutation_rate=0, max_mutation_rate=0, custom_model=None, productive=False):
        """Initialize an S5F mutation model with specified parameters."""
        self.targeting = None
        self.substitution = None
        self.mutability = None
        self.max_mutation_rate = max_mutation_rate
        self.min_mutation_rate = min_mutation_rate
        self.bases = {'A', 'T', 'C', 'G'}
        self.loaded_metadata = False
        self.custom_model = custom_model
        self.productive = productive

    def load_metadata(self, sequence):
        """Loads mutation model metadata based on the sequence type.

        Args:
            sequence (Sequence): The sequence object, used to determine the appropriate metadata to load.

        Raises:
            ValueError: If the sequence type is unsupported.
        """
        from ..sequence import HeavyChainSequence
        from ..sequence.light_chain import LightChainSequence
        from importlib import resources

        if self.custom_model is None:
            if type(sequence) == HeavyChainSequence:
                with resources.path('GenAIRR.data', 'HH_S5F_META.pkl') as data_path:
                    with open(data_path, 'rb') as h:
                        self.mutability, self.substitution, self.targeting = pickle.load(h)
            elif type(sequence) == LightChainSequence:
                with resources.path('GenAIRR.data', 'HKL_S5F_META.pkl') as data_path:
                    with open(data_path, 'rb') as h:
                        self.mutability, self.substitution, self.targeting = pickle.load(h)
            else:
                raise ValueError('Unsupported Sequence Type')
        else:
            with open(self.custom_model, 'rb') as h:
                self.mutability, self.substitution, self.targeting = pickle.load(h)

    def apply_mutation(self, sequence_object):
        """Applies mutations to a given sequence object based on the S5F mutation model.

        Args:
            sequence_object (Sequence): The sequence object to which mutations will be applied.

        Returns:
            tuple: A tuple containing the mutated sequence, a dictionary of mutations, and the mutation rate.
        """
        # 1. Load the Likelihoods File
        if not self.loaded_metadata:
            self.load_metadata(sequence_object)
            self.loaded_metadata = True

        # 1.1 Sample Mutation Rate
        mutation_rate = random.uniform(self.min_mutation_rate, self.max_mutation_rate)
        target_number_of_mutations = int(mutation_rate * len(sequence_object.ungapped_seq))

        # have a copy of the original neucliotdies
        naive_sequence = list(sequence_object.ungapped_seq)

        # Log mutations
        mutations = dict()


        # 2. Extract 5-Mers
        fiver_mers = FiveMER.create_five_mers(sequence_object.ungapped_seq, self.mutability)
                            
        # add a failsafe to insure while loop does not get locked
        patience = 0

        # duplicate the loop for unproductive to not repeat if.
        # productive
        if self.productive:
            # if productive do not include the positions of the anchors. 
            reading_frame = [[2, 1, 0][idx % 3] for idx, element in enumerate(fiver_mers)]
            v_anchor = sequence_object.junction_start
            j_anchor = sequence_object.junction_end
            restricted_positions = {v_anchor, v_anchor+1, v_anchor+2, j_anchor-3, j_anchor-2, j_anchor-1}
            
            while len(mutations) < target_number_of_mutations:
                # 3. Mutability, Weighted Choice of Position Based on 5-Mer Likelihoods
                sampled_position, chosen_index = self.weighted_choice(fiver_mers)  # likelihoods are normalized here
                # re-sample the position if it's in the anchor
                in_anchor = chosen_index in restricted_positions
                while in_anchor:
                    sampled_position, chosen_index = self.weighted_choice(fiver_mers)  # likelihoods are normalized here
                    in_anchor = chosen_index in restricted_positions
                # 4. Substitution
                substitutions = self.substitution[sampled_position.sequence].dropna()  # drop Nan's - N's and Same Base
                mutable_bases = substitutions.index
                bases_likelihoods = substitutions.values
                mutation_to_apply = random.choices(mutable_bases, weights=bases_likelihoods, k=1)[0]
                # if mutation creates a stop codon and productive is true. generate a new mutations.
                # under the assumption that the reading frame is 0. 
                stop, codon = self._is_stop(sampled_position.nucleotides, mutation_to_apply, reading_frame[chosen_index])
                while stop:
                    sampled_position, chosen_index = self.weighted_choice(fiver_mers)
                    in_anchor = chosen_index in restricted_positions
                    while in_anchor:
                        sampled_position, chosen_index = self.weighted_choice(fiver_mers)  # likelihoods are normalized here
                        in_anchor = chosen_index in restricted_positions
                    substitutions = self.substitution[sampled_position.sequence].dropna()
                    mutable_bases = substitutions.index
                    bases_likelihoods = substitutions.values
                    mutation_to_apply = random.choices(mutable_bases, weights=bases_likelihoods, k=1)[0]
                    stop, codon = self._is_stop(sampled_position.nucleotides, mutation_to_apply, reading_frame[chosen_index])
                # log
                if sampled_position.position not in mutations:
                    mutations[sampled_position.position] = f'{sampled_position.sequence[2]}>{mutation_to_apply}'
                else:
                    mutations[sampled_position.position] += f'>{mutation_to_apply}'
                # if mutation reverted previous mutation back to naive state, drop that record from the log
                if mutation_to_apply == naive_sequence[sampled_position.position]:
                    mutations.pop(sampled_position.position)
                # 5. Apply Mutation
                # This will also update all relevant 5-MERS and their likelihood with pointer like logic
                sampled_position.change_center(mutation_to_apply, self.mutability)
                patience += 1
                # Patience logic
                if patience > (target_number_of_mutations * 30) and target_number_of_mutations > 1:
                    patience = 0
                    # restart process
                    mutation_rate = random.uniform(self.min_mutation_rate, self.max_mutation_rate)
                    target_number_of_mutations = int(mutation_rate * len(sequence_object.ungapped_seq))
                    mutations = dict()
                    # 2. Extract 5-Mers
                    fiver_mers = FiveMER.create_five_mers(sequence_object.ungapped_seq, self.mutability)
        else:
            ## normal
            while len(mutations) < target_number_of_mutations:
                # 3. Mutability, Weighted Choice of Position Based on 5-Mer Likelihoods
                sampled_position, chosen_index = self.weighted_choice(fiver_mers)  # likelihoods are normalized here
                
                # 4. Substitution
                substitutions = self.substitution[sampled_position.sequence].dropna()  # drop Nan's - N's and Same Base
                mutable_bases = substitutions.index
                bases_likelihoods = substitutions.values
                mutation_to_apply = random.choices(mutable_bases, weights=bases_likelihoods, k=1)[0]
        
                # log
                if sampled_position.position not in mutations:
                    mutations[sampled_position.position] = f'{sampled_position.sequence[2]}>{mutation_to_apply}'
                else:
                    mutations[sampled_position.position] += f'>{mutation_to_apply}'

                # if mutation reverted previous mutation back to naive state, drop that record from the log
                if mutation_to_apply == naive_sequence[sampled_position.position]:
                    mutations.pop(sampled_position.position)

                # 5. Apply Mutation
                # This will also update all relevant 5-MERS and their likelihood with pointer like logic
                sampled_position.change_center(mutation_to_apply, self.mutability)

                patience += 1
                # Patience logic
                if patience > (target_number_of_mutations * 30) and target_number_of_mutations > 1:
                    patience = 0
                    # restart process
                    mutation_rate = random.uniform(self.min_mutation_rate, self.max_mutation_rate)
                    target_number_of_mutations = int(mutation_rate * len(sequence_object.ungapped_seq))
                    mutations = dict()
                    # 2. Extract 5-Mers
                    fiver_mers = FiveMER.create_five_mers(sequence_object.ungapped_seq, self.mutability)
                        
        mutated_sequence = FiveMER.five_mers_to_dna(fiver_mers)
        return mutated_sequence, mutations, mutation_rate

    def _mutate_base(self, base):
        """Randomly selects a different base as a mutation of the given base.

            Args:
                base (str): The original nucleotide base.

            Returns:
                str: The mutated nucleotide base.
            """
        return random.choice(list(self.bases - {base}))

    @staticmethod
    def weighted_choice(five_mers):
        """Selects a FiveMER object from a list, weighted by their likelihoods.

        Args:
            five_mers (list): A list of FiveMER objects.

        Returns:
            FiveMER: A randomly selected FiveMER object, weighted by likelihood.
        """
        weights = [fm.likelihood if fm.likelihood == fm.likelihood else 0 for fm in five_mers]
        # Choose an index instead of the object
        chosen_index = random.choices(range(len(five_mers)), weights, k=1)[0]
        return five_mers[chosen_index], chosen_index

    def _is_stop(self, nucleotides, new_base, reading_frame):
        """Check if the mutation introduce a stop codon

        Args:
            nucleotides (_type_): The nucleotides of the five mer
            new_base (_type_): The new mutated base
        """
        # replace center
        nucs = [str(nuc) for nuc in nucleotides]
        codon = nucs[0:2] + [new_base] + nucs[3:]
        codon = ''.join(codon[reading_frame:reading_frame+3])
        stop = codon in ["TAG", "TAA", "TGA"]
        return stop, codon