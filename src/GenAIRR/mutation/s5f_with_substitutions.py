import random

import numpy as np

from .s5f import FiveMER
from ..mutation.mutation_model import MutationModel
import pickle

class S5F_w_Substitutions(MutationModel):
    """Implements the S5F mutation model with Uniform Substituions, a specific model for simulating mutations in DNA sequences.

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

    def __init__(self, min_mutation_rate=0, max_mutation_rate=0, custom_model=None, productive=False,
                 substitution_probability=0):
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
        self.substitution_probability = substitution_probability

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

    def mutable_positions(self, sequence):
        """Identifies mutable positions in the sequence, excluding (NP) regions.

        Args:
            sequence (Sequence): The sequence object containing V, D, and J region start and end positions.
            ignor_anchors (bool): Whether to ignor the anchors positions, defaulting to false.
        Returns:
            list: A list of positions that are eligible for mutation, combining positions from V, D, and J regions.
        """

        positions_to_mutate = []

        # add v region positions
        positions_to_mutate += list(range(sequence.v_seq_start, sequence.v_seq_end))
        positions_to_mutate += list(range(sequence.d_seq_start, sequence.d_seq_end))
        positions_to_mutate += list(range(sequence.j_seq_start, sequence.j_seq_end))

        return positions_to_mutate
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
            restricted_positions = {v_anchor: 'v',
                                    v_anchor + 1: 'v',
                                    v_anchor + 2: 'v',
                                    j_anchor - 3: 'j',
                                    j_anchor - 2: 'j',
                                    j_anchor - 1: 'j'}
            while len(mutations) < target_number_of_mutations:
                sampled_position, mutation_to_apply = self._productive_recursive(0, fiver_mers, reading_frame,
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

                # if true perform random substitution instead of s5f mutation:
                if np.random.binomial(1,self.substitution_probability,size=1).item():
                    sampled_position, chosen_index = self.weighted_choice(fiver_mers)
                    base = sampled_position[2].nucleotides.value
                    mutation_to_apply = random.choice(list(self.bases - {base}))

                    # log
                    if sampled_position.position not in mutations:
                        mutations[sampled_position.position] = f'S|{sampled_position.sequence[2]}>{mutation_to_apply}'
                    else:
                        mutations[sampled_position.position] += f'>{mutation_to_apply}'

                    # if mutation reverted previous mutation back to naive state, drop that record from the log
                    if mutation_to_apply == naive_sequence[sampled_position.position]:
                        mutations.pop(sampled_position.position)

                else:
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

    def _is_stop_codon(self, nucleotides, new_base, reading_frame):
        """Check if the mutation introduce a stop codon

        Args:
            nucleotides (_type_): The nucleotides of the five mer
            new_base (_type_): The new mutated base
        """
        # replace center
        nucs = [str(nuc) for nuc in nucleotides]
        codon = nucs[0:2] + [new_base] + nucs[3:]
        codon = ''.join(codon[reading_frame:reading_frame + 3])
        stop = codon in ["TAG", "TAA", "TGA"]
        return stop, codon

    def _productive_recursive(self, counter, fiver_mers, reading_frame, restricted_positions):
        if counter >= 1000:
            raise RecursionError("Maximum recursion depth exceeded")
        counter += 1
        sampled_position, chosen_index = self.weighted_choice(fiver_mers)
        substitutions = self.substitution[sampled_position.sequence].dropna()
        mutable_bases = substitutions.index
        bases_likelihoods = substitutions.values
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


