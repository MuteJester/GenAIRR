import pickle
import random
import numpy as np
import scipy.stats as st
from ..utilities import AlleleNComparer
from ..sequence import LightChainType
from ..sequence import LightChainSequence
from ..simulation import SequenceAugmentorArguments
from ..simulation.sequence_augmentor_base import SequenceAugmentorBase
from ..utilities.data_config import DataConfig

class LightChainSequenceAugmentor(SequenceAugmentorBase):
    """
        A class for augmenting light chain immunoglobulin sequences, extending the capabilities of the SequenceAugmentorBase class.

        This class specifically handles the augmentation of light chain sequences, taking into account the unique aspects of kappa and lambda light chains.

        Attributes:
            alleles_used (list): Specifies the alleles ('v', 'j') used in light chain sequences.

        Args:
            dataconfig (DataConfig): The data configuration object containing settings for sequence simulation.
            args (SequenceAugmentorArguments): The arguments object containing simulation parameters.
    """
    alleles_used = ['v', 'j']

    def __init__(self, dataconfig: DataConfig, args: SequenceAugmentorArguments = SequenceAugmentorArguments()):
        super().__init__(dataconfig, args)

        self.nucleotide_add_distribution = st.beta(2, 3)
        self.nucleotide_remove_distribution = st.beta(2, 3)
        self.nucleotide_add_after_remove_distribution = st.beta(1, 3)
        self.chain_type = LightChainType.KAPPA if 'IGKV' in list(self.dataconfig.v_alleles)[
            0] else LightChainType.LAMBDA
        # Loading Routines
        self.load_correction_maps()

    # Loading Routines
    def load_correction_maps(self):
        """
               Loads allele correction maps for V gene segment start and end positions to account for similarities between different alleles after trimming events.

               These maps help in adjusting allele choices based on the amount removed from the start or end of a given allele, identifying equally likely options.
       """
        from importlib import resources
        """
        This will load the V start and V end maps that tell us given the amount removed from the start or end
        of a given allele what option are equally likely
        :return:
        """

        self.v_start_allele_correction_map = self.dataconfig.correction_maps['V_5_TRIM_SIMILARITY_MAP']
        self.max_v_start_correction_map_value = max(
        self.v_start_allele_correction_map[list(self.v_start_allele_correction_map)[0]])

        self.v_end_allele_correction_map = self.dataconfig.correction_maps['V_3_TRIM_SIMILARITY_MAP']
        self.max_v_end_correction_map_value = max(
        self.v_end_allele_correction_map[list(self.v_end_allele_correction_map)[0]])

        self.v_n_ambiguity_comparer = self.dataconfig.correction_maps['V_N_AMBIGUITY_CORRECTION_GRAPH']

    def attribute_position(self, simulated, position):
        """
        Determines the region (V, NP1, J) of a given position in a simulated light chain sequence.

        Args:
            simulated (dict): A dictionary containing the simulated sequence and metadata.
            position (int): The position in the sequence to be attributed.

        Returns:
            str: The region name ('v', 'np1', 'j') corresponding to the given position.

        Raises:
            ValueError: If the position does not fall within defined sequence regions.
        """
        if simulated['v_sequence_start'] <= position < simulated['v_sequence_end']:
            return 'v'
        elif simulated['v_sequence_end'] <= position < simulated['j_sequence_start']:
            return 'np1'
        elif simulated['j_sequence_start'] <= position < simulated['j_sequence_end']:
            return 'j'
        else:
            raise ValueError('Undefined Position!')


    # Sequence Simulation
    def simulate_sequence(self):
        """
            Simulates a light chain sequence, applying mutations according to the configured mutation model.

                This method generates a random light chain sequence, mutates it, and compiles relevant metadata into a dictionary.

                Returns:
                    dict: A dictionary containing the simulated sequence, mutation rate, trimming details, and allele information.
        """
        # Sample Sequence
        gen = LightChainSequence.create_random(self.dataconfig)

        gen.mutate(self.mutation_model)

        data = {
            "sequence": gen.mutated_seq,
            "v_sequence_start": gen.v_seq_start,
            "v_sequence_end": gen.v_seq_end,
            "j_sequence_start": gen.j_seq_start,
            "j_sequence_end": gen.j_seq_end,
            "v_allele": [gen.v_allele.name],
            "j_allele": [gen.j_allele.name],
            'mutation_rate': gen.mutation_freq,
            'v_trim_5': gen.v_trim_5,
            'v_trim_3': gen.v_trim_3,
            'j_trim_5': gen.j_trim_5,
            'j_trim_3': gen.j_trim_3,
            'type': self.chain_type,
            'corruption_event': 'no-corruption',
            'corruption_add_amount': 0,
            'corruption_remove_amount': 0,
            'mutations': {pos: gen.mutations[pos] for pos in sorted(gen.mutations)},  # sort the mutations by position
            "Ns": dict(),
        }
        return data

    def fix_v_position_after_trimming_index_ambiguity(self, simulation):
        """
                Adjusts the V gene segment end position in the simulated sequence to account for potential overlaps between the generated junction and the reference after a trimming event.

                Args:
                    simulation (dict): A dictionary containing the simulated sequence and metadata, including V allele positions and trimming details.
        """
        # Extract Current V Metadata
        v_start, v_end = simulation['v_sequence_start'], simulation['v_sequence_end']
        v_allele_remainder = simulation['sequence'][v_start:v_end]
        v_allele_ref = self.v_dict[simulation['v_allele'][0]]

        # Get the junction inserted after trimming to the sequence
        junction_3 = simulation['sequence'][v_end:simulation['j_sequence_start']]

        # Get the trimming lengths
        v_trim_3 = simulation['v_trim_3']

        # Get the trimmed off sections from the reference
        trimmed_3 = v_allele_ref[len(v_allele_ref) - v_trim_3:]

        # check for overlap between generated junction and reference in the 3' trim
        for a, b in zip(trimmed_3, junction_3):
            # in case the current poistion in the junction matches the reference exapnd the d segment
            if a == b:
                v_end += 1
            else:  # if the continuous streak is broken or non-existent break!
                break

        simulation['v_sequence_end'] = v_end

    def fix_j_position_after_trimming_index_ambiguity(self, simulation):
        """
                Corrects the J gene segment start position in the simulated sequence to address possible overlaps between the generated junction and the reference following a trimming event.

                Args:
                    simulation (dict): A dictionary containing the simulated sequence and metadata, including J allele positions and trimming information.
        """
        # Extract Current J Metadata
        j_start, j_end = simulation['j_sequence_start'], simulation['j_sequence_end']
        j_allele_remainder = simulation['sequence'][j_start:j_end]
        j_allele_ref = self.j_dict[simulation['j_allele'][0]]

        # Get the junction inserted after trimming to the sequence
        junction_5 = simulation['sequence'][simulation['v_sequence_end']:j_start]

        # Get the trimming lengths
        j_trim_5 = simulation['j_trim_5']

        # Get the trimmed off sections from the reference
        trimmed_5 = j_allele_ref[:j_trim_5]

        # check for overlap between generated junction and reference in the 5' trim
        for a, b in zip(trimmed_5[::-1], junction_5[::-1]):
            # in case the current poistion in the junction matches the reference exapnd the d segment
            if a == b:
                j_start -= 1
            else:  # if the continuous streak is broken or non-existent break!
                break

        simulation['j_sequence_start'] = j_start

    def remove_event(self, simulated):
        """
                Simulates the removal of a portion of the sequence from the beginning for a light chain sequence, updating metadata to reflect the changes.

                Args:
                    simulated (dict): A dictionary containing the simulated sequence and metadata.
        """
        # Update Simulation Metadata
        simulated['corruption_event'] = 'remove'
        # Sample how much to remove
        v_length = simulated['v_sequence_end'] - simulated['v_sequence_start']
        amount_to_remove = self._sample_nucleotide_remove_distribution(v_length)
        # remove from the start of the sequence the sampled amount
        simulated['sequence'] = simulated['sequence'][amount_to_remove:]
        # Update Simulation Metadata
        simulated['corruption_remove_amount'] = amount_to_remove
        # Update mutation log, remove the mutations that were removed with the remove event
        # and while updating the mutations log correct the position of the mutations accordingly
        simulated['mutations'] = {i - amount_to_remove: j for i, j in simulated['mutations'].items() if
                                  i > amount_to_remove}
        # Adjust mutation rate
        self.correct_mutation_rate(simulated)

        # Adjust Start/End Position Accordingly
        simulated['v_sequence_start'] = 0
        simulated['v_sequence_end'] -= amount_to_remove
        simulated['j_sequence_start'] -= amount_to_remove
        simulated['j_sequence_end'] -= amount_to_remove

        # Correction - Add All V Alleles That Cant be Distinguished Based on the Amount Cut from the V Allele
        self.correct_for_v_start_cut(simulated)

    def add_event(self, simulated, amount=None):
        """
                Simulates the addition of bases to the beginning of a light chain sequence, updating metadata to reflect the changes.

                Args:
                    simulated (dict): A dictionary containing the simulated sequence and metadata.
                    amount (int, optional): The number of bases to add. If not specified, the amount is sampled from a configured distribution.
        """
        # Update Simulation Metadata
        simulated['corruption_event'] = 'add'
        # Sample the Amount to Add by default, if a spesific value is given via the amount variable,use it.
        amount_to_add = self._sample_nucleotide_add_distribution() if amount is None else amount
        # Sample the method by which addition will be made
        method = self._sample_corruption_add_method()
        # Modify the sequence
        modified_sequence = method(amount_to_add, simulated['sequence'])
        # Validate the modified sequence, make sure we didn't over add pass our max sequence size
        modified_sequence, amount_to_add = self.validate_sequence_length_after_addition(modified_sequence,
                                                                                        amount_to_add)
        # Update Simulation Sequence
        simulated['sequence'] = modified_sequence

        # Update Simulation Metadata
        simulated['corruption_add_amount'] = amount_to_add
        # Update mutation log positions
        simulated['mutations'] = {i + amount_to_add: j for i, j in simulated['mutations'].items()}

        # Adjust Start/End Position Accordingly
        simulated['v_sequence_start'] += amount_to_add
        simulated['v_sequence_end'] += amount_to_add
        simulated['j_sequence_start'] += amount_to_add
        simulated['j_sequence_end'] += amount_to_add

    def remove_before_add_event(self, simulated):
        """
                Simulates a sequence modification event where a portion is first removed from the start of the sequence, and then bases are added, updating metadata accordingly.

                Args:
                    simulated (dict): A dictionary containing the simulated sequence and metadata.
        """
        # ----- REMOVE PART -----#
        # Sample how much to remove
        self.remove_event(simulated)

        # Update Simulation Metadata
        simulated['corruption_event'] = 'remove_before_add'

        # ----- ADD PART -----#
        # Sample how much to add after removal occurred
        amount_to_add = self._sample_nucleotide_add_after_remove_distribution()
        self.add_event(simulated, amount=amount_to_add)

    def simulate_augmented_sequence(self):
        """
            Simulates and augments a light chain sequence, incorporating sequence corrections, potential corruption at the beginning, and insertion of 'N' bases.

                Returns:
                    dict: A dictionary containing the augmented sequence and associated metadata.
        """
        # 1. Simulate a Sequence from AIRRship
        simulated = self.simulate_sequence()

        # 1.1 Correction - Correct Start/End Positions Based on Generated Junctions
        self.fix_v_position_after_trimming_index_ambiguity(simulated)
        self.fix_j_position_after_trimming_index_ambiguity(simulated)

        # 1.2 Correction - Add All V Alleles That Cant be Distinguished Based on the Amount Cut from the V Allele
        self.correct_for_v_end_cut(simulated)

        # 2. Corrupt Begging of Sequence ( V Start )
        if self.perform_corruption():
            # Inside this method, based on the corruption event we will also adjust the respective ground truth v
            # alleles
            self.corrupt_sequence_beginning(simulated)

        # 3. Add N's
        self.insert_Ns(simulated)

        self.process_before_return(simulated)

        return simulated


class LightChainKappaLambdaSequenceAugmentor:
    """
        A class for augmenting both kappa and lambda light chain immunoglobulin sequences.

        This class manages two separate augmentors for kappa and lambda chains, allowing for the simulation of sequences based on a specified kappa to lambda ratio.

        Args:
            lambda_dataconfig (DataConfig): The data configuration object for lambda chain sequence simulation.
            kappa_dataconfig (DataConfig): The data configuration object for kappa chain sequence simulation.
            lambda_args (SequenceAugmentorArguments): The arguments object containing simulation parameters for lambda chains.
            kappa_args (SequenceAugmentorArguments): The arguments object containing simulation parameters for kappa chains.
    """
    def __init__(self, lambda_dataconfig: DataConfig, kappa_dataconfig: DataConfig,
                 lambda_args: SequenceAugmentorArguments = SequenceAugmentorArguments(),
                 kappa_args: SequenceAugmentorArguments = SequenceAugmentorArguments()):
        self.lambda_augmentor = LightChainSequenceAugmentor(lambda_dataconfig, lambda_args)
        self.kappa_augmentor = LightChainSequenceAugmentor(kappa_dataconfig, kappa_args)
        self.kappa_lambda_ratio = lambda_args.kappa_lambda_ratio

    def simulate_augmented_sequence(self):
        """
                Simulates and augments either a kappa or lambda light chain sequence based on a predefined kappa to lambda ratio.

                Randomly chooses between kappa and lambda chain types to simulate a sequence accordingly.

                Returns:
                    dict: A dictionary containing the augmented sequence and associated metadata for the chosen chain type.
        """
        chain_type = np.random.choice([LightChainType.KAPPA, LightChainType.LAMBDA], size=1, p=[self.kappa_lambda_ratio,
                                                                                                1 - self.kappa_lambda_ratio]).item()
        if chain_type == LightChainType.KAPPA:
            return self.kappa_augmentor.simulate_augmented_sequence()
        else:
            return self.lambda_augmentor.simulate_augmented_sequence()

    @property
    def columns(self):
        """
                Provides a list of column names expected in the output of the simulate_augmented_sequence method.

                The column names are derived from the kappa chain augmentor, assuming both kappa and lambda augmentors produce the same set of columns.

                Returns:
                    list: A list of strings representing the column names in the augmented sequence output.
        """
        return self.kappa_augmentor.columns