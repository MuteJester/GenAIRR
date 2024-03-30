import random
from ..mutation.mutation_model import MutationModel


class Uniform(MutationModel):
    """Implements a uniform mutation model as a subclass of MutationModel.

    This mutation model applies mutations across the sequence at a uniform rate, ignoring non-polymorphic (NP) regions.

    Attributes:
        min_mutation_rate (float): The minimum mutation rate to be considered when applying mutations.
        max_mutation_rate (float): The maximum mutation rate to be considered when applying mutations.
        bases (set): A set containing the nucleotide bases ('A', 'T', 'C', 'G') used in mutations.

    Args:
        min_mutation_rate (float): The minimum mutation rate, defaulting to 0.
        max_mutation_rate (float): The maximum mutation rate, defaulting to 0.
        productive (bool): Wheter to ensure the sequence is productive (No stop codons and mutation in cdr3 anchors), defaulting to false.
    """
    def __init__(self, min_mutation_rate=0, max_mutation_rate=0, productive=False):
        """Initializes the Uniform mutation model with specified minimum and maximum mutation rates."""
        self.max_mutation_rate = max_mutation_rate
        self.min_mutation_rate = min_mutation_rate
        self.bases = {'A','T','C','G'}
        self.productive = productive


    def mutable_positions(self,sequence):
        """Identifies mutable positions in the sequence, excluding (NP) regions.

        Args:
            sequence (Sequence): The sequence object containing V, D, and J region start and end positions.
            ignor_anchors (bool): Whether to ignor the anchors positions, defaulting to false.
        Returns:
            list: A list of positions that are eligible for mutation, combining positions from V, D, and J regions.
        """

        positions_to_mutate = []

        # add v region positions
        positions_to_mutate += list(range(sequence.v_seq_start,sequence.v_seq_end))
        positions_to_mutate += list(range(sequence.d_seq_start,sequence.d_seq_end))
        positions_to_mutate += list(range(sequence.j_seq_start,sequence.j_seq_end))
    
        return positions_to_mutate

    def apply_mutation(self, sequence_object):
        """Applies mutations to the given sequence object based on the uniform mutation model.

        The method randomly selects positions to mutate based on the calculated mutation rate and mutates them to a different base.

        Args:
            sequence_object (Sequence): The sequence object to which mutations will be applied. It must have an 'ungapped_seq' attribute.

        Returns:
            tuple: A tuple containing the mutated sequence, a dictionary of mutations with positions as keys and mutations as values, and the mutation rate.
        """
        sequence = sequence_object.ungapped_seq
        mutation_rate = random.uniform(self.min_mutation_rate, self.max_mutation_rate)
        number_of_mutations = int(mutation_rate*len(sequence))
        possible_positions_to_mutate = self.mutable_positions(sequence_object)
        positions_to_mutate = random.sample(possible_positions_to_mutate,k=number_of_mutations)

        # log mutations
        mutations = dict()
        mutated_sequence = list(sequence)
        if self.productive:
            v_anchor = sequence_object.junction_start
            j_anchor = sequence_object.junction_end
            restricted_positions = {v_anchor:'v', 
                                    v_anchor+1:'v', 
                                    v_anchor+2:'v', 
                                    j_anchor-3:'j', 
                                    j_anchor-2:'j', 
                                    j_anchor-1:'j'}
            for position in positions_to_mutate:
                position, new_base = self._productive_recursive(0, position, mutated_sequence, restricted_positions, possible_positions_to_mutate)
                mutations[position] = f'{mutated_sequence[position]}>{new_base}'
                mutated_sequence[position] = new_base
        else:
            for position in positions_to_mutate:
                new_base = self._mutate_base(mutated_sequence[position])
                mutations[position] = f'{mutated_sequence[position]}>{new_base}'
                mutated_sequence[position] = new_base

        mutated_sequence = ''.join(mutated_sequence)
        return mutated_sequence, mutations, mutation_rate

    def _mutate_base(self, base):
        """Mutates a given nucleotide base to a different base.

        Args:
            base (str): The original nucleotide base to be mutated.

        Returns:
            str: A new nucleotide base, different from the original.
        """
        return random.choice(list(self.bases - {base}))
    
    def _is_stop_codon(self, mutated_sequence, position, new_base):
        """Asses if the mutation in the creates a stop codon

        Args:
            sequence (str): baseline mutated sequence
            position (int): mutated position
            new_base (str): mutated nucleotide 
        
        Returns:
            bool: Whether the codon is stop.
        """
        if(position>=(len(mutated_sequence)-1)):
            return False
        positions = [(position - (position % 3)) + i for i in range(3)]
        codon = ''.join( [new_base if pos == position else mutated_sequence[pos] for pos in positions])
        return codon in ["TAG", "TAA", "TGA"]
        
    def _productive_recursive(self, counter, position, mutated_sequence, restricted_positions, positions_to_mutate):
        if counter >= 100:
            raise RecursionError("Maximum recursion depth exceeded")
        new_base = self._mutate_base(mutated_sequence[position])
        stop = self._is_stop_codon(mutated_sequence, position, new_base)
        if stop:
            position = random.sample(positions_to_mutate,k=1)[0]
            return self._productive_recursive(counter, position, mutated_sequence, restricted_positions, positions_to_mutate = positions_to_mutate)
        elif position in restricted_positions.keys():
            tag = restricted_positions[position]
            codon = ''.join([new_base if k ==position else mutated_sequence[k] for k, v in restricted_positions.items() if v == tag])
            pattern = {'TGC', 'TGT'} if tag == 'v' else {'TTT', 'TTC', 'TGG'}
            if codon in pattern:
                return position, new_base
            else:
                position = random.sample(positions_to_mutate,k=1)[0]
                return self._productive_recursive(counter, position, mutated_sequence, restricted_positions, positions_to_mutate = positions_to_mutate)
        return position, new_base