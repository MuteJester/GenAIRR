import random
import numpy as np
import scipy.stats as st
from ..utilities import translate
from ..sequence import LightChainType
from ..sequence import LightChainSequence
from ..simulation import SequenceAugmentorArguments
from ..simulation.sequence_augmentor_base import SequenceAugmentorBase
from GenAIRR.dataconfig.data_config import DataConfig

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
    alleles_used = ['v', 'j','c']

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

    def apply_deletion(self, simulated, position):
        sequence = simulated['sequence']
        nucleotides_list = list(sequence)

        deleted = nucleotides_list[position]
        # print('-=-'*20)
        # print(f'Deleting {deleted} From Position: {position}')
        # print('-=-'*20)

        after_deletion = nucleotides_list[:position] + nucleotides_list[position + 1:]

        # update start/end positions
        for reg in ['v_sequence_start', 'v_sequence_end', 'j_sequence_start', 'j_sequence_end']:
            if simulated[reg] > position:
                simulated[reg] = simulated[reg] - 1

        simulated['indels'][position] = 'D > ' + deleted

        # correct N's and Mutations
        corrected_mutations = {}
        for pos in simulated['mutations']:
            if pos > position:
                corrected_mutations[max(position, pos - 1)] = simulated['mutations'][pos]
            else:
                corrected_mutations[pos] = simulated['mutations'][pos]

        simulated['mutations'] = corrected_mutations

        corrected_Ns = {}
        for pos in simulated['Ns']:
            if pos >= position:
                corrected_Ns[max(position, pos - 1)] = simulated['Ns'][pos]
            else:
                corrected_Ns[pos] = simulated['Ns'][pos]

        simulated['Ns'] = corrected_Ns
        # Update log positions for deletions
        updated_log = {}
        for log_pos, log_entry in simulated['indels'].items():
            if log_pos >= position:
                updated_log[log_pos - 1] = log_entry
            elif log_pos < position:
                updated_log[log_pos] = log_entry

        simulated['sequence'] = ''.join(after_deletion)

    def apply_insertion(self, simulated, position):
        sequence = simulated['sequence']
        nucleotides_list = list(sequence)
        random_base = random.choice(['A', 'T', 'C', 'G'])

        # print('-=-'*20)
        # print(f'Inserting {random_base} To Position: {position}')
        # print('-=-'*20)
        after_insertion = nucleotides_list[:position] + [random_base] + nucleotides_list[position:]
        # update start/end positions
        for reg in ['v_sequence_start', 'v_sequence_end', 'j_sequence_start', 'j_sequence_end']:
            if simulated[reg] >= position:
                simulated[reg] += 1

        simulated['indels'][position] = 'I < ' + random_base
        # correct N's and Mutations
        corrected_mutations = {}
        for pos in simulated['mutations']:
            if pos >= position:
                corrected_mutations[pos + 1] = simulated['mutations'][pos]
            else:
                corrected_mutations[pos] = simulated['mutations'][pos]

        simulated['mutations'] = corrected_mutations

        corrected_Ns = {}
        for pos in simulated['Ns']:
            if pos >= position:
                corrected_Ns[pos + 1] = simulated['Ns'][pos]
            else:
                corrected_Ns[pos] = simulated['Ns'][pos]

        simulated['Ns'] = corrected_Ns

        # Update log positions for insertions
        updated_log = {}
        for log_pos, log_entry in simulated['indels'].items():
            if log_pos > position:
                updated_log[log_pos + 1] = log_entry
            else:
                updated_log[log_pos] = log_entry
        simulated['indels'] = updated_log

        simulated['sequence'] = ''.join(after_insertion)

    def valid_indel_positions(self, simulated):
        all_positions = set(range(0, simulated['j_sequence_end']))
        # remove np regions
        np1_positions = set(range(simulated['v_sequence_end'], simulated['j_sequence_start']))
        all_positions -= np1_positions

        # remove ns positions
        all_positions -= set(simulated['Ns'].keys())
        # remove mutations position
        all_positions -= set(simulated['mutations'].keys())

        return all_positions

    def insert_indels(self, simulated):
        # get valid position for indels excluding np regions, n's and mutated positions
        valid_positions = list(self.valid_indel_positions(simulated))
        num_indels = np.random.randint(1, self.max_indels, size=1).item()
        num_indels = min(num_indels, len(valid_positions))
        random.shuffle(valid_positions)
        n_valid_positions = len(valid_positions)

        for idx in range(num_indels):
            indel_position = valid_positions[idx]

            # choose action 1 = insertion -1 = deletion
            action = np.random.choice([1, -1], size=1, p=[self.insertion_proba, self.deletion_proba]).item()
            if action == 1:  # insertion case
                # print('###'*30)
                # print('Before Insertion')
                # print("V segment: Start - {v_start}, End - {v_end}\n"
                #       "D segment: Start - {d_start}, End - {d_end}\n"
                #       "J segment: Start - {j_start}, End - {j_end}".format(
                #           v_start=simulated['v_sequence_start'],
                #           v_end=simulated['v_sequence_end'],
                #           d_start=simulated['d_sequence_start'],
                #           d_end=simulated['d_sequence_end'],
                #           j_start=simulated['j_sequence_start'],
                #           j_end=simulated['j_sequence_end']))
                # print(simulated['Ns'])
                # print(simulated['mutations'])

                self.apply_insertion(simulated, indel_position)
                # print('After Insertion')
                # print("V segment: Start - {v_start}, End - {v_end}\n"
                #       "D segment: Start - {d_start}, End - {d_end}\n"
                #       "J segment: Start - {j_start}, End - {j_end}".format(
                #           v_start=simulated['v_sequence_start'],
                #           v_end=simulated['v_sequence_end'],
                #           d_start=simulated['d_sequence_start'],
                #           d_end=simulated['d_sequence_end'],
                #           j_start=simulated['j_sequence_start'],
                #           j_end=simulated['j_sequence_end']))
                # print(simulated['Ns'])
                # print(simulated['mutations'])
                for update_idx in range(idx, idx + (num_indels - idx)):
                    if valid_positions[update_idx] >= indel_position:
                        valid_positions[update_idx] += 1
            else:  # deletion case

                # print('###'*30)
                # print('Before Deletion')
                # print("V segment: Start - {v_start}, End - {v_end}\n"
                #       "D segment: Start - {d_start}, End - {d_end}\n"
                #       "J segment: Start - {j_start}, End - {j_end}".format(
                #           v_start=simulated['v_sequence_start'],
                #           v_end=simulated['v_sequence_end'],
                #           d_start=simulated['d_sequence_start'],
                #           d_end=simulated['d_sequence_end'],
                #           j_start=simulated['j_sequence_start'],
                #           j_end=simulated['j_sequence_end']))
                # print(simulated['Ns'])
                # print(simulated['mutations'])
                self.apply_deletion(simulated, indel_position)
                #                 print('After Deletion')

                #                 print("V segment: Start - {v_start}, End - {v_end}\n"
                #                       "D segment: Start - {d_start}, End - {d_end}\n"
                #                       "J segment: Start - {j_start}, End - {j_end}".format(
                #                           v_start=simulated['v_sequence_start'],
                #                           v_end=simulated['v_sequence_end'],
                #                           d_start=simulated['d_sequence_start'],
                #                           d_end=simulated['d_sequence_end'],
                #                           j_start=simulated['j_sequence_start'],
                #                           j_end=simulated['j_sequence_end']))
                #                 print(simulated['Ns'])
                #                 print(simulated['mutations'])

                for update_idx in range(idx, idx + (num_indels - idx)):
                    if valid_positions[update_idx] > indel_position:
                        valid_positions[update_idx] -= 1

    # Sequence Simulation
    def simulate_sequence(self, specific_v=None, specific_j=None):
        """
            Simulates a light chain sequence, applying mutations according to the configured mutation model.

                This method generates a random light chain sequence, mutates it, and compiles relevant metadata into a dictionary.

                Returns:
                    dict: A dictionary containing the simulated sequence, mutation rate, trimming details, and allele information.
        """
        # Sample Sequence
        gen = LightChainSequence.create_random(self.dataconfig,specific_v=specific_v,specific_j=specific_j)
        if self.productive:
            while not gen.functional:
                gen = LightChainSequence.create_random(self.dataconfig)
        
        gen.mutate(self.mutation_model)
        
        data = {
            "sequence": gen.mutated_seq,
            "v_sequence_start": gen.v_seq_start,
            "v_sequence_end": gen.v_seq_end,
            "j_sequence_start": gen.j_seq_start,
            "j_sequence_end": gen.j_seq_end,
            "v_germline_start": gen.v_germline_start,
            "v_germline_end": gen.v_germline_end,
            "j_germline_start": gen.j_germline_start,
            "j_germline_end": gen.j_germline_end,
            "junction_sequence_start": gen.junction_start,
            "junction_sequence_end": gen.junction_end,
            'v_call': [gen.v_allele.name],
            'j_call': [gen.j_allele.name],
            'c_call': [gen.c_allele.name],
            'mutation_rate': gen.mutation_freq,
            'v_trim_5': gen.v_trim_5,
            'v_trim_3': gen.v_trim_3,
            'j_trim_5': gen.j_trim_5,
            'j_trim_3': gen.j_trim_3,
            'c_trim_3': gen.c_trim_3,
            'type': self.chain_type,
            'corruption_event': 'no-corruption',
            'corruption_add_amount': 0,
            'corruption_remove_amount': 0,
            'corruption_removed_section': '',
            'corruption_added_section': '',
            'mutations': {pos: gen.mutations[pos] for pos in sorted(gen.mutations)},  # sort the mutations by position
            "Ns": dict(),
            'indels': dict(),
            'productive': gen.functional,
            'stop_codon': gen.stop_codon,
            'vj_in_frame': gen.vj_in_frame,
            'note': gen.note
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
        v_germline_start, v_germline_end = simulation['v_germline_start'], simulation['v_germline_end']
        v_allele_remainder = simulation['sequence'][v_start:v_end]
        v_allele_ref = self.v_dict[simulation['v_call'][0]]

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
                v_germline_end += 1
            else:  # if the continuous streak is broken or non-existent break!
                break

        simulation['v_sequence_end'] = v_end
        simulation['v_germline_end'] = v_germline_end

    def fix_j_position_after_trimming_index_ambiguity(self, simulation):
        """
                Corrects the J gene segment start position in the simulated sequence to address possible overlaps between the generated junction and the reference following a trimming event.

                Args:
                    simulation (dict): A dictionary containing the simulated sequence and metadata, including J allele positions and trimming information.
        """
        # Extract Current J Metadata
        j_start, j_end = simulation['j_sequence_start'], simulation['j_sequence_end']
        j_germline_start, j_germline_end = simulation['j_germline_start'], simulation['j_germline_end']
        j_allele_remainder = simulation['sequence'][j_start:j_end]
        j_allele_ref = self.j_dict[simulation['j_call'][0]]

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
                j_germline_start -= 1
            else:  # if the continuous streak is broken or non-existent break!
                break

        simulation['j_sequence_start'] = j_start
        simulation['j_germline_start'] = j_germline_start

    def distill_mutation_rate(self, simulated):
        # remove mutations in NP region
        np_positions = list(range(simulated['v_sequence_end'], simulated['j_sequence_start'] + 1))
        mutations_in_np_regions = list(set(np_positions) & set(list(simulated['mutations'].keys())))
        for pos in mutations_in_np_regions:
            simulated['mutations'].pop(pos)

        distilled_mutation_rate = (len(simulated['mutations']) + len(simulated['Ns'])) / len(simulated['sequence'])
        simulated['mutation_rate'] = distilled_mutation_rate

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
        simulated['junction_sequence_end'] -= amount_to_remove
        simulated['junction_sequence_end'] -= amount_to_remove
        simulated['v_germline_start'] += amount_to_remove

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
        modified_sequence = method(amount_to_add, simulated)
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
        simulated['junction_sequence_end'] += amount_to_add
        simulated['junction_sequence_end'] += amount_to_add
        
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

    def fix_productive_call_after_corruption_indel(self,simulated):
        sequence = simulated['sequence']
        functional = False
        stop_codon = False
        note = simulated['note']
        # stop codon
        stops = ["TAG", "TAA", "TGA"]
        for x in range(simulated['junction_sequence_end'], simulated['v_sequence_start'], -3):
            if sequence[x-3:x] in stops:
                stop_codon = True
        if not stop_codon:
            for x in range(simulated['junction_sequence_end'], len(sequence), 3):
                if sequence[x:x+3] in stops:
                    stop_codon = True
        # vj in frame
        from_j_to_start = (simulated['junction_sequence_end'] - simulated['v_sequence_start']) % 3 == 0
        junction_length = (simulated['junction_sequence_end'] - simulated['junction_sequence_start']) % 3 == 0
        from_v_to_start = (simulated['junction_sequence_start'] - simulated['v_sequence_start']) % 3 == 0
        vj_in_frame = from_j_to_start and junction_length and from_v_to_start and stop_codon == False
        junction = sequence[simulated['junction_sequence_start']:simulated['junction_sequence_end']]
        # prooductivity
        if junction_length and stop_codon == False:
            junction_aa = translate(junction)
            if junction_aa.startswith("C"):
                if junction_aa.endswith("F") or junction_aa.endswith("W"):
                    functional = True
                else:
                    functional = False
                    note += 'J anchor (W/F) not present.'
            else:
                functional = False
                note += 'V second C not present.'
        else:
            if stop_codon == True:
                note += 'Stop codon present.'
            if junction_length == False:
                note += 'Junction length not divisible by 3.'
            functional = False
        # update the values
        simulated['productive'] = functional
        simulated['stop_codon'] = stop_codon
        simulated['vj_in_frame'] = vj_in_frame
        simulated['note'] = note
        
    def simulate_augmented_sequence(self, specific_v=None, specific_j=None):
        """
            Simulates and augments a light chain sequence, incorporating sequence corrections, potential corruption at the beginning, and insertion of 'N' bases.

                Returns:
                    dict: A dictionary containing the augmented sequence and associated metadata.
        """
        # 1. Simulate a Naive Sequence
        simulated = self.simulate_sequence(specific_v=specific_v, specific_j=specific_j)

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
            self.fix_productive_call_after_corruption_indel(simulated)

        # 3. Add N's
        if bool(np.random.binomial(1, self.n_proba)):
            self.insert_Ns(simulated)

        # Insert Indels:
        if bool(np.random.binomial(1,self.simulate_indels)):
            self.insert_indels(simulated)
            self.fix_productive_call_after_corruption_indel(simulated)
            
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

    def simulate_augmented_sequence(self, specific_v=None, specific_j=None, chain_type=None):
        """
        Simulates and augments either a kappa or lambda light chain sequence based on a predefined kappa to lambda ratio.

        Parameters:
            specific_v (Optional[allele object]): The specific V allele to use for simulation. If provided, `chain_type` must also be specified.
            specific_j (Optional[allele object]): The specific J allele to use for simulation. If provided, `chain_type` must also be specified.
            chain_type (Optional[LightChainType]): The type of light chain to simulate. Must be either LightChainType.KAPPA or LightChainType.LAMBDA.
                                                   If provided, the method will use this value and not sample randomly.

        If `specific_v` or `specific_j` is provided, `chain_type` must also be specified. If only one of `specific_v` or `specific_j` is provided, the
        remaining allele will be randomly generated. If neither `specific_v` nor `specific_j` is provided, a random chain type will be chosen unless
        `chain_type` is specified.

        Returns:
            dict: A dictionary containing the augmented sequence and associated metadata for the chosen chain type.

        Raises:
            ValueError: If `specific_v` or `specific_j` is provided without specifying `chain_type`.
        """

        if (specific_v or specific_j) and not chain_type:
            raise ValueError("If `specific_v` or `specific_j` is provided, `chain_type` must also be specified.")

        if chain_type:
            if chain_type == LightChainType.KAPPA:
                return self.kappa_augmentor.simulate_augmented_sequence(specific_v=specific_v, specific_j=specific_j)
            elif chain_type == LightChainType.LAMBDA:
                return self.lambda_augmentor.simulate_augmented_sequence(specific_v=specific_v, specific_j=specific_j)
            else:
                raise ValueError("`chain_type` must be either LightChainType.KAPPA or LightChainType.LAMBDA.")
        else:
            chain_type = np.random.choice([LightChainType.KAPPA, LightChainType.LAMBDA], size=1,
                                          p=[self.kappa_lambda_ratio, 1 - self.kappa_lambda_ratio]).item()
            if chain_type == LightChainType.KAPPA:
                return self.kappa_augmentor.simulate_augmented_sequence(specific_v,specific_j)
            else:
                return self.lambda_augmentor.simulate_augmented_sequence(specific_v,specific_j)

    @property
    def columns(self):
        """
                Provides a list of column names expected in the output of the simulate_augmented_sequence method.

                The column names are derived from the kappa chain augmentor, assuming both kappa and lambda augmentors produce the same set of columns.

                Returns:
                    list: A list of strings representing the column names in the augmented sequence output.
        """
        return self.kappa_augmentor.columns
