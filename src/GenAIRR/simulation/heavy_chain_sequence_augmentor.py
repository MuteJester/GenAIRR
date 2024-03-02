import pickle
import random
import numpy as np
import scipy.stats as st
from ..utilities import AlleleNComparer
from ..sequence import HeavyChainSequence
from ..simulation import SequenceAugmentorArguments
from ..simulation.sequence_augmentor_base import SequenceAugmentorBase
from ..utilities.data_config import DataConfig
import base64


class HeavyChainSequenceAugmentor(SequenceAugmentorBase):
    """
        Augmentor class for simulating and augmenting heavy chain immunoglobulin sequences.

        This class extends the SequenceAugmentorBase to include functionalities specific to heavy chain sequences,
        such as handling V, D, and J gene segments, and applying heavy chain-specific corrections and augmentations.

        Attributes:
            alleles_used (list): List of gene segment types used in heavy chain sequences, typically including 'v', 'd', and 'j'.

        Args:
            dataconfig (DataConfig): Configuration object containing allele information and simulation parameters.
            args (SequenceAugmentorArguments): Object containing arguments for sequence simulation and augmentation.
    """
    alleles_used = ['v', 'd', 'j']

    def __init__(self, dataconfig: DataConfig, args: SequenceAugmentorArguments = SequenceAugmentorArguments()):
        super().__init__(dataconfig, args)

        self.nucleotide_add_distribution = st.beta(2, 3)
        self.nucleotide_remove_distribution = st.beta(2, 3)
        self.nucleotide_add_after_remove_distribution = st.beta(1, 3)

        self.short_d_length = args.short_d_length
        # Class Misc
        self.d_alleles = sorted([i for j in self.dataconfig.d_alleles for i in self.dataconfig.d_alleles[j]],
                                key=lambda x: x.name)
        self.d_dict = {i.name: i.ungapped_seq.upper() for i in self.d_alleles}

        # Loading Routines
        self.load_correction_maps()

    # Loading Routines
    def load_correction_maps(self):
        """
                Loads correction maps for V and D gene segments from the configuration object.

                These maps are used to adjust allele choices based on trimming events and to resolve ambiguities caused by similar alleles.
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

        self.d_trim_correction_map = self.dataconfig.correction_maps['D_5_3_TRIM_SIMILARITY_MAP']

        self.v_n_ambiguity_comparer = self.dataconfig.correction_maps['V_N_AMBIGUITY_CORRECTION_GRAPH']

    def attribute_position(self, simulated, position):
        """
        Determines the gene segment (V, D, J, or NP regions) associated with a given position in a simulated heavy chain sequence.

        Args:
            simulated (dict): Dictionary containing the simulated sequence and its metadata.
            position (int): The position within the sequence to be attributed.

        Returns:
            str: The gene segment or region ('v', 'd', 'j', 'np1', 'np2') associated with the given position.

        Raises:
            ValueError: If the position does not fall within the defined segments or regions.
        """
        if simulated['v_sequence_start'] <= position < simulated['v_sequence_end']:
            return 'v'
        elif simulated['v_sequence_end'] <= position < simulated['d_sequence_start']:
            return 'np1'
        elif simulated['d_sequence_start'] <= position < simulated['d_sequence_end']:
            return 'd'
        elif simulated['d_sequence_end'] <= position < simulated['j_sequence_start']:
            return 'np2'
        elif simulated['j_sequence_start'] <= position < simulated['j_sequence_end']:
            return 'j'
        else:
            raise ValueError('Undefined Position!')

    def apply_deletion(self, simulated, position):
        pass

        sequence = simulated['sequence']
        nucleotides_list = list(sequence)
        region = self.attribute_position(simulated, position)

    def correct_for_d_trims(self, simulated):
        """
                Adjusts the D allele choices in the simulated sequence based on 5' and 3' trimming information.

                Args:
                    simulated (dict): Dictionary containing the simulated sequence and its metadata, including D allele trims.
        """
        # Get the 5' and 3' trims of the d allele in the simulated sequence
        trim_5 = simulated['d_trim_5']
        trim_3 = simulated['d_trim_3']
        # infer the precalculated map what alleles should be the ground truth for this sequence based on the trim
        simulated['d_allele'] = list(self.d_trim_correction_map[simulated['d_allele'][0]][(trim_5, trim_3)])

    def short_d_validation(self, simulated):
        """
                Validates the length of the D gene segment in the simulated sequence. If the segment is shorter than a predefined threshold, it is labeled as "Short-D".

                Args:
                    simulated (dict): Dictionary containing the simulated sequence and its metadata.
        """
        d_length = simulated['d_sequence_end'] - simulated['d_sequence_start']
        if d_length < self.short_d_length:
            simulated['d_allele'] = ['Short-D']

    def fix_d_position_after_trimming_index_ambiguity(self, simulation):
        """
                Corrects the start and end positions of the D gene segment in the simulated sequence to resolve ambiguities caused by trimming events.

                Args:
                    simulation (dict): Dictionary containing the simulated sequence and its metadata, including D allele positions and trimming details.
        """
        # Extract Current D Metadata
        d_start, d_end = simulation['d_sequence_start'], simulation['d_sequence_end']
        d_allele_remainder = simulation['sequence'][d_start:d_end]
        d_allele_ref = self.d_dict[simulation['d_allele'][0]]

        # Get the junction inserted after trimming to the sequence
        junction_5 = simulation['sequence'][simulation['v_sequence_end']:d_start]
        junction_3 = simulation['sequence'][d_end:simulation['j_sequence_start']]

        # Get the trimming lengths
        d_trim_5 = simulation['d_trim_5']
        d_trim_3 = simulation['d_trim_3']

        # Get the trimmed off sections from the reference
        trimmed_5 = d_allele_ref[:d_trim_5]
        trimmed_3 = d_allele_ref[len(d_allele_ref) - d_trim_3:]

        # check for overlap between generated junction and reference in the 5' trim
        for a, b in zip(trimmed_5[::-1], junction_5[::-1]):
            # in case the current poistion in the junction matches the reference exapnd the d segment
            if a == b:
                d_start -= 1
            else:  # if the continuous streak is broken or non-existent break!
                break

        # check for overlap between generated junction and reference in the 3' trim
        for a, b in zip(trimmed_3, junction_3):
            # in case the current poistion in the junction matches the reference exapnd the d segment
            if a == b:
                d_end += 1
            else:  # if the continious streak is broken or non existant break!
                break

        simulation['d_sequence_start'] = d_start
        simulation['d_sequence_end'] = d_end

    # Sequence Simulation
    def simulate_sequence(self):
        """
                Simulates a heavy chain sequence, applying mutations and generating relevant metadata.

                Returns:
                    dict: A dictionary containing the simulated sequence, its metadata including gene segment positions, alleles, mutation rate, and trimming information.
        """
        gen = HeavyChainSequence.create_random(self.dataconfig)
        gen.mutate(self.mutation_model)

        data = {
            "sequence": gen.mutated_seq,
            "v_sequence_start": gen.v_seq_start,
            "v_sequence_end": gen.v_seq_end,
            "d_sequence_start": gen.d_seq_start,
            "d_sequence_end": gen.d_seq_end,
            "j_sequence_start": gen.j_seq_start,
            "j_sequence_end": gen.j_seq_end,
            "v_allele": [gen.v_allele.name],
            "d_allele": [gen.d_allele.name],
            "j_allele": [gen.j_allele.name],
            'mutation_rate': gen.mutation_freq,
            'v_trim_5': gen.v_trim_5,
            'v_trim_3': gen.v_trim_3,
            'd_trim_5': gen.d_trim_5,
            'd_trim_3': gen.d_trim_3,
            'j_trim_5': gen.j_trim_5,
            'j_trim_3': gen.j_trim_3,
            'corruption_event': 'no-corruption',
            'corruption_add_amount': 0,
            'corruption_remove_amount': 0,
            'mutations': {pos: gen.mutations[pos] for pos in sorted(gen.mutations)},  # sort the mutations by position
            "Ns": dict(),
        }
        return data

    def simulate_augmented_sequence(self):
        """
                Generates an augmented heavy chain sequence, applying sequence corrections, potential corruption events at the beginning, and insertion of 'N' bases.

                Returns:
                    dict: A dictionary containing the augmented sequence and associated metadata, ready for further analysis or training processes.
        """
        # 1. Simulate a Sequence from AIRRship
        simulated = self.simulate_sequence()

        # 1.1 Correction - Correct Start/End Positions Based on Generated Junctions
        self.fix_v_position_after_trimming_index_ambiguity(simulated)
        self.fix_d_position_after_trimming_index_ambiguity(simulated)
        self.fix_j_position_after_trimming_index_ambiguity(simulated)

        # 1.2 Correction - Add All V Alleles That Cant be Distinguished Based on the Amount Cut from the V Allele
        self.correct_for_v_end_cut(simulated)

        # 1.3 Correction - Add All D Alleles That Cant be Distinguished Based on the 5' and 3' Trims
        self.correct_for_d_trims(simulated)

        # 2. Corrupt Begging of Sequence ( V Start )
        if self.perform_corruption():
            # Inside this method, based on the corruption event we will also adjust the respective ground truth v
            # alleles
            self.corrupt_sequence_beginning(simulated)

        # 3. Add N's
        self.insert_Ns(simulated)

        # 4. Adjust D Allele , if the simulated length is smaller than short_d_length
        # class property, change the d allele to the "Short-D" label
        self.short_d_validation(simulated)


        return simulated
