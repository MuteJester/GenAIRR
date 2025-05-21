import random

import pandas as pd
from ..utilities.misc import check_stops, STOP_CODONS
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
            self.likelihood = mutability.get(self.sequence, self.likelihood)

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
        # Step 1: Pad the sequence and precompute some constants
        padded_sequence = f'NN{dna_sequence}NN'
        sequence_length = len(dna_sequence)

        # Step 2: Create FiveMER objects and populate the list
        five_mers = []
        for i in range(sequence_length):  # Iterate based on the original sequence length
            # Extract the 5-mer string directly
            five_mer_nucleotides = padded_sequence[i:i + 5]

            # Create the FiveMER object
            five_mer = FiveMER([Nucleotide(nuc) for nuc in five_mer_nucleotides])
            five_mer.position = i  # Position of the original second nucleotide (central in unpadded)

            # Apply mutability likelihood if available
            if mutability is not None:
                five_mer.likelihood = mutability.get(five_mer_nucleotides, 0)

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
        self.stop_codons = {"TAG", "TAA", "TGA"}

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
                with resources.path('GenAIRR.data.mutation_model_parameters', 'HH_S5F_META.pkl') as data_path:
                    with open(data_path, 'rb') as h:
                        self.mutability, self.substitution, self.targeting = pickle.load(h)
            elif type(sequence) == LightChainSequence:
                with resources.path('GenAIRR.data.mutation_model_parameters', 'HKL_S5F_META.pkl') as data_path:
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
            # check if substitution map is a dataframe convert it to a nested dict to make access more efficient
            self.substitution = self.substitution.to_dict(orient='dict')
            # Remove keys with NaN values in each dictionary
            self.substitution = {
                outer_key: {inner_key: inner_value for inner_key, inner_value in outer_dict.items() if
                            not pd.isna(inner_value)}
                for outer_key, outer_dict in self.substitution.items()
            }
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
            restricted_positions = {v_anchor: 'v',
                                    v_anchor + 1: 'v',
                                    v_anchor + 2: 'v',
                                    j_anchor - 3: 'j',
                                    j_anchor - 2: 'j',
                                    j_anchor - 1: 'j'}
            while len(mutations) < target_number_of_mutations:
                sampled_position, mutation_to_apply = self._productive_constrained_mutation(fiver_mers, reading_frame,
                                                                                            restricted_positions)
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
                #debug
                #prev_postion = sampled_position.sequence[2]
                sampled_position.change_center(mutation_to_apply, self.mutability)
                # if check_stops(FiveMER.five_mers_to_dna(fiver_mers)):
                #     sampled_position.change_center(prev_postion, self.mutability)
                #     prev_postion = None

                patience += 1
                # Patience logic
                if check_stops(FiveMER.five_mers_to_dna(fiver_mers)):
                    mutation_rate = random.uniform(self.min_mutation_rate, self.max_mutation_rate)
                    target_number_of_mutations = int(mutation_rate * len(sequence_object.ungapped_seq))
                    mutations = dict()
                    # 2. Extract 5-Mers
                    fiver_mers = FiveMER.create_five_mers(sequence_object.ungapped_seq, self.mutability)
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
                substitutions = self.substitution[sampled_position.sequence]
                mutable_bases = list(substitutions.keys())
                bases_likelihoods = list(substitutions.values())
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
        total_weight = 0
        cumulative_weights = []

        for fm in five_mers:
            weight = fm.likelihood if fm.likelihood == fm.likelihood else 0
            total_weight += weight
            cumulative_weights.append(total_weight)

        r = random.uniform(0, total_weight)
        for i, cumulative_weight in enumerate(cumulative_weights):
            if r < cumulative_weight:
                return five_mers[i], i

    def _is_stop_codon(self, nucleotides, new_base, reading_frame):
        """Check if the mutation introduces a stop codon in the current, previous, or next reading frame.

        Args:
            nucleotides (list): The nucleotides of the five mer.
            new_base (str): The new mutated base.
            reading_frame (int): The current reading frame.

        Returns:
            tuple: A tuple containing a boolean indicating if a stop codon is present and the codon sequence.
        """
        # Replace center nucleotide with the new base
        nucs = [str(nuc) for nuc in nucleotides]
        nucs[2] = new_base

        # Check the current reading frame
        codon = ''.join(nucs[reading_frame:reading_frame + 3])
        tests = []
        if codon in STOP_CODONS:
            tests.append(True)

        # Check the previous reading frame if possible
        if reading_frame > 0:
            prev_codon = ''.join(nucs[reading_frame - 1:reading_frame + 2])
            if prev_codon in STOP_CODONS:
                tests.append(True)

        # Check the next reading frame if possible
        if reading_frame < 2:
            next_codon = ''.join(nucs[reading_frame + 1:reading_frame + 4])
            if next_codon in STOP_CODONS:
                tests.append(True)

        return any(tests), codon

    def _test_stop_codon_formation(self, fivemer):
        stop_codon_nucleotides = []

        # Check each possible nucleotide replacement
        for nucleotide in self.bases:
            # Replace the center nucleotide
            new_sequence = list(fivemer.sequence)
            new_sequence[2] = nucleotide  # Center position is index 2
            new_sequence = ''.join(new_sequence)

            # Check if any 3-mer in the new sequence is a stop codon
            for i in range(len(new_sequence) - 2):
                if new_sequence[i:i + 3] in self.stop_codons:
                    stop_codon_nucleotides.append(nucleotide)
                    break
        return stop_codon_nucleotides

    def _reweight_substitution_likelihoods(self, likelihood_dict, stop_codon_nucleotides):
        # Remove stop codon nucleotides from the likelihood dictionary
        copy_of_likelihood_dict = likelihood_dict.copy()

        for stop_codon_nucleotide in stop_codon_nucleotides:
            copy_of_likelihood_dict.pop(stop_codon_nucleotide, None)

        # Normalize the remaining likelihoods
        total_likelihood = sum(copy_of_likelihood_dict.values())
        for base in copy_of_likelihood_dict:
            copy_of_likelihood_dict[base] /= total_likelihood

        return copy_of_likelihood_dict

    def _productive_constrained_mutation(self, fiver_mers, reading_frame, restricted_positions):
        counter = 0
        while counter < 1000:
            counter += 1
            sampled_position, chosen_index = self.weighted_choice(fiver_mers)
            stop_codon_potential_nucleotides = self._test_stop_codon_formation(sampled_position)
            substitutions = self.substitution[sampled_position.sequence]
            substitutions_modified = self._reweight_substitution_likelihoods(substitutions,
                                                                             stop_codon_potential_nucleotides)
            mutable_bases = list(substitutions_modified.keys())
            bases_likelihoods = list(substitutions_modified.values())

            if len(bases_likelihoods) == 0:  # no valid mutation will avoid stop codon formation
                continue

            mutation_to_apply = random.choices(mutable_bases, weights=bases_likelihoods, k=1)[0]
            stop, codon = self._is_stop_codon(sampled_position.nucleotides, mutation_to_apply,
                                              reading_frame[chosen_index])
            if stop:
                continue
            elif chosen_index in restricted_positions.keys():
                tag = restricted_positions[chosen_index]
                _reading_frame = reading_frame[chosen_index]
                nucs = [str(nuc) for nuc in sampled_position.nucleotides]
                codon = nucs[0:2] + [mutation_to_apply] + nucs[3:]
                codon = ''.join(codon[_reading_frame:_reading_frame + 3])
                pattern = {'TGC', 'TGT'} if tag == 'v' else {'TTT', 'TTC', 'TGG'}
                if codon in pattern:
                    return sampled_position, mutation_to_apply
                else:
                    continue
            return sampled_position, mutation_to_apply
        raise AssertionError("Maximum iteration depth exceeded")

    def _productive_recursive(self, counter, fiver_mers, reading_frame, restricted_positions):
        if counter >= 1000:
            raise RecursionError("Maximum recursion depth exceeded")
        counter += 1
        sampled_position, chosen_index = self.weighted_choice(fiver_mers)
        stop_codon_potential_nucleotides = self._test_stop_codon_formation(sampled_position)
        substitutions = self.substitution[sampled_position.sequence]
        substitutions_modified = self._reweight_substitution_likelihoods(substitutions,
                                                                         stop_codon_potential_nucleotides)
        mutable_bases = list(substitutions_modified.keys())
        bases_likelihoods = list(substitutions_modified.values())
        if len(bases_likelihoods) == 0:  # no valid mutation will avoid stop codon formation
            return self._productive_recursive(counter, fiver_mers, reading_frame, restricted_positions)
        mutation_to_apply = random.choices(mutable_bases, weights=bases_likelihoods, k=1)[0]
        stop, codon = self._is_stop_codon(sampled_position.nucleotides, mutation_to_apply, reading_frame[chosen_index])
        if stop:
            return self._productive_recursive(counter, fiver_mers, reading_frame, restricted_positions)
        elif chosen_index in restricted_positions.keys():
            tag = restricted_positions[chosen_index]
            _reading_frame = reading_frame[chosen_index]
            nucs = [str(nuc) for nuc in sampled_position.nucleotides]
            codon = nucs[0:2] + [mutation_to_apply] + nucs[3:]
            codon = ''.join(codon[_reading_frame:_reading_frame + 3])
            pattern = {'TGC', 'TGT'} if tag == 'v' else {'TTT', 'TTC', 'TGG'}
            if codon in pattern:
                return sampled_position, mutation_to_apply
            else:
                return self._productive_recursive(counter, fiver_mers, reading_frame, restricted_positions)
        return sampled_position, mutation_to_apply
