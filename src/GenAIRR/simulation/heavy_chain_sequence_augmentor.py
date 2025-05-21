import random
import numpy as np
import scipy.stats as st
from ..utilities import translate
from ..sequence import HeavyChainSequence
from ..simulation import SequenceAugmentorArguments
from ..simulation.sequence_augmentor_base import SequenceAugmentorBase
from GenAIRR.dataconfig.data_config import DataConfig


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
    alleles_used = ['v', 'd', 'j','c']

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

        # list alleles with 

    # Loading Routines
    def load_correction_maps(self):
        """
                Loads correction maps for V and D gene segments from the configuration object.

                These maps are used to adjust allele choices based on trimming events and to resolve ambiguities caused by similar alleles.
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
        sequence = simulated['sequence']
        nucleotides_list = list(sequence)

        deleted = nucleotides_list[position]

        after_deletion = nucleotides_list[:position] + nucleotides_list[position + 1:]

        # update start/end positions
        for reg in ['v_sequence_start', 'v_sequence_end', 'j_sequence_start', 'j_sequence_end', 'd_sequence_start',
                    'd_sequence_end', 'junction_sequence_start' , 'junction_sequence_end']:
            if simulated[reg] > position:
                simulated[reg] = simulated[reg] - 1

        simulated['indels'][position] = 'D > ' + deleted

        # correct N's and Mutations
        corrected_mutations = {}
        for pos in simulated['mutations']:
            if pos > position:
                corrected_mutations[pos-1] = simulated['mutations'][pos]
            else:
                corrected_mutations[pos] = simulated['mutations'][pos]

        simulated['mutations'] = corrected_mutations

        corrected_Ns = {}
        for pos in simulated['Ns']:
            if pos >= position:
                corrected_Ns[pos-1] = simulated['Ns'][pos]
            else:
                corrected_Ns[pos] = simulated['Ns'][pos]

        simulated['Ns'] = corrected_Ns


        simulated['sequence'] = ''.join(after_deletion)
        # Update log positions for deletions
        updated_log = {}
        for log_pos, log_entry in simulated['indels'].items():
            if log_pos >= position:
                updated_log[log_pos - 1] = log_entry
            elif log_pos < position:
                updated_log[log_pos] = log_entry

        simulated['indels'] = updated_log

    def apply_insertion(self, simulated, position):
        sequence = simulated['sequence']
        nucleotides_list = list(sequence)
        random_base = random.choice(['A', 'T', 'C', 'G'])


        after_insertion = nucleotides_list[:position] + [random_base] + nucleotides_list[position:]
        # update start/end positions
        for reg in ['v_sequence_start', 'v_sequence_end', 'j_sequence_start', 'j_sequence_end', 'd_sequence_start',
                    'd_sequence_end', 'junction_sequence_start' , 'junction_sequence_end']:
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

        simulated['sequence'] = ''.join(after_insertion)

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
        np1_positions = set(range(simulated['v_sequence_end'], simulated['d_sequence_start']))
        np2_positions = set(range(simulated['d_sequence_end'], simulated['j_sequence_start']))
        all_positions -= np1_positions
        all_positions -= np2_positions

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

                self.apply_insertion(simulated, indel_position)

                for update_idx in range(idx, idx + (num_indels - idx)):
                    if valid_positions[update_idx] >= indel_position:
                        valid_positions[update_idx] += 1
            else:  # deletion case

                self.apply_deletion(simulated, indel_position)

                for update_idx in range(idx, idx + (num_indels - idx)):
                    if valid_positions[update_idx] > indel_position:
                        valid_positions[update_idx] -= 1

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
        #simulated['d_call'] = list(self.d_trim_correction_map[simulated['d_call'][0]][(trim_5, trim_3)])
        # retaining the order such that the selected d allele is the first.
        sampled_d = simulated['d_call'][0]
        simulated['d_call'] = [sampled_d] + list(set(self.d_trim_correction_map[sampled_d][(trim_5, trim_3)]) - set([sampled_d]))

    def short_d_validation(self, simulated):
        """
                Validates the length of the D gene segment in the simulated sequence. If the segment is shorter than a predefined threshold, it is labeled as "Short-D".

                Args:
                    simulated (dict): Dictionary containing the simulated sequence and its metadata.
        """
        d_length = simulated['d_sequence_end'] - simulated['d_sequence_start']
        if d_length < self.short_d_length:
            simulated['d_call'] = ['Short-D']

    def fix_d_position_after_trimming_index_ambiguity(self, simulation):
        """
                Corrects the start and end positions of the D gene segment in the simulated sequence to resolve ambiguities caused by trimming events.

                Args:
                    simulation (dict): Dictionary containing the simulated sequence and its metadata, including D allele positions and trimming details.
        """
        # Extract Current D Metadata
        d_start, d_end = simulation['d_sequence_start'], simulation['d_sequence_end']
        d_germline_start, d_germline_end = simulation['d_germline_start'], simulation['d_germline_end']
        d_allele_remainder = simulation['sequence'][d_start:d_end]
        d_allele_ref = self.d_dict[simulation['d_call'][0]]

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
                d_germline_start -= 1
                d_trim_5 -= 1
            else:  # if the continuous streak is broken or non-existent break!
                break

        # check for overlap between generated junction and reference in the 3' trim
        for a, b in zip(trimmed_3, junction_3):
            # in case the current poistion in the junction matches the reference exapnd the d segment
            if a == b:
                d_end += 1
                d_germline_end += 1
                d_trim_3 -= 1

            else:  # if the continious streak is broken or non existant break!
                break

        simulation['d_sequence_start'] = d_start
        simulation['d_sequence_end'] = d_end
        simulation['d_germline_start'] = d_germline_start
        simulation['d_germline_end'] = d_germline_end
        simulation['d_trim_5'] = d_trim_5
        simulation['d_trim_3'] = d_trim_3

    def distill_mutation_rate(self,simulated):
        # remove mutations in NP region
        np_positions = list(range(simulated['v_sequence_end'],simulated['d_sequence_start']+1))
        np_positions += list(range(simulated['d_sequence_end'],simulated['j_sequence_start']+1))
        mutations_in_np_regions = list(set(np_positions)&set(list(simulated['mutations'].keys())))
        for pos in mutations_in_np_regions:
            simulated['mutations'].pop(pos)

        distilled_mutation_rate = (len(simulated['mutations'])+len(simulated['Ns']))/len(simulated['sequence'])
        simulated['mutation_rate'] = distilled_mutation_rate

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
    
    # Sequence Simulation
    def simulate_sequence(self,specific_v=None,specific_d=None,specific_j=None):
        """
                Simulates a heavy chain sequence, applying mutations and generating relevant metadata.
                specific_v,specific_d,specific_j are argument you can pass alleles object to in order to force simulation
                of a sequence with specific v,d,j parameters.
                If only 1 of the arguments is passed,the remaining will be randomly generated.

                Returns:
                    dict: A dictionary containing the simulated sequence, its metadata including gene segment positions, alleles, mutation rate, and trimming information.
        """
        # if productive is true. ensure the base sequence is productive. This will keep the VJ in frame prior to indels.
        gen = HeavyChainSequence.create_random(self.dataconfig,specific_v=specific_v,
                                               specific_d=specific_d,specific_j=specific_j)
        if self.productive:
            while not gen.functional:
                gen = HeavyChainSequence.create_random(self.dataconfig,specific_v=specific_v,
                                                       specific_d=specific_d,specific_j=specific_j)
        
        gen.mutate(self.mutation_model)
        
        data = {
            "sequence": gen.mutated_seq,
            "v_sequence_start": gen.v_seq_start,
            "v_sequence_end": gen.v_seq_end,
            "d_sequence_start": gen.d_seq_start,
            "d_sequence_end": gen.d_seq_end,
            "j_sequence_start": gen.j_seq_start,
            "j_sequence_end": gen.j_seq_end,
            "v_germline_start": gen.v_germline_start,
            "v_germline_end": gen.v_germline_end,
            "d_germline_start": gen.d_germline_start,
            "d_germline_end": gen.d_germline_end,
            "j_germline_start": gen.j_germline_start,
            "j_germline_end": gen.j_germline_end,
            "junction_sequence_start": gen.junction_start,
            "junction_sequence_end": gen.junction_end,
            "v_call": [gen.v_allele.name],
            "d_call": [gen.d_allele.name],
            "j_call": [gen.j_allele.name],
            "c_call": [gen.c_allele.name],
            'mutation_rate': gen.mutation_freq,
            'v_trim_5': gen.v_trim_5,
            'v_trim_3': gen.v_trim_3,
            'd_trim_5': gen.d_trim_5,
            'd_trim_3': gen.d_trim_3,
            'j_trim_5': gen.j_trim_5,
            'j_trim_3': gen.j_trim_3,
            'c_trim_3': gen.c_trim_3,
            'corruption_event': 'no-corruption',
            'corruption_add_amount': 0,
            'corruption_remove_amount': 0,
            'corruption_removed_section':'',
            'corruption_added_section':'',
            'mutations': {pos: gen.mutations[pos] for pos in sorted(gen.mutations)},  # sort the mutations by position
            "Ns": dict(),
            'indels': dict(),
            'productive': gen.functional,
            'stop_codon': gen.stop_codon,
            'vj_in_frame': gen.vj_in_frame,
            'note': gen.note
        }
        return data

    def simulate_augmented_sequence(self,specific_v=None,specific_d=None,specific_j=None):
        """
            Generates an augmented heavy chain sequence, applying sequence corrections, potential corruption events at the beginning,indels, and insertion of 'N' bases.

            Parameters:
                specific_v (Optional[allele object]): The specific V allele to use for simulation. If provided, will override random generation for V allele.
                specific_d (Optional[allele object]): The specific D allele to use for simulation. If provided, will override random generation for D allele.
                specific_j (Optional[allele object]): The specific J allele to use for simulation. If provided, will override random generation for J allele.

            If only one or two of `specific_v`, `specific_d`, or `specific_j` are provided, the remaining alleles will be randomly generated.

            Returns:
                dict: A dictionary containing the augmented sequence and associated metadata, ready for further analysis or training processes.
         """
        # 1. Simulate a Naive Sequence
        simulated = self.simulate_sequence(specific_v=specific_v,specific_d=specific_d,specific_j=specific_j)

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
            self.fix_productive_call_after_corruption_indel(simulated)

        # 3. Add N's
        if bool(np.random.binomial(1, self.n_proba)):
            self.insert_Ns(simulated)

        # 4. Adjust D Allele , if the simulated length is smaller than short_d_length
        # class property, change the d allele to the "Short-D" label
        self.short_d_validation(simulated)

        # Insert Indels with probability = simulate_indels :
        if bool(np.random.binomial(1,self.simulate_indels)):
            self.insert_indels(simulated)
            self.fix_productive_call_after_corruption_indel(simulated)

        self.distill_mutation_rate(simulated)
        
        # add sequence productivey assesment
        #self.productive_cdr3_and_stop_codon(simulated)
        
        self.process_before_return(simulated)

        return simulated
