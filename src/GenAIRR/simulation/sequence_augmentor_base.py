from abc import ABC, abstractmethod
from GenAIRR.dataconfig.data_config import DataConfig
from dataclasses import dataclass, field
import enum
from ..mutation import MutationModel, S5F
import scipy.stats as st
import random
import numpy as np


# Corruption Events
class Event(enum.Enum):
    Remove = 1
    Add = 2
    Remove_Before_Add = 3


# Arg Class
@dataclass
class SequenceAugmentorArguments:
    """
        Configuration arguments for sequence augmentation processes, specifying mutation rates, indel simulation parameters, and other augmentation settings.

        Attributes:
            min_mutation_rate (float): The minimum mutation rate for simulating sequence mutations, defaulting to 0.003.
            max_mutation_rate (float): The maximum mutation rate for simulating sequence mutations, defaulting to 0.25.
            simulate_indels (float): The probability of simulating a random amount of indels.
            max_indels (int): The maximum number of indels to simulate, defaulting to 5.
            deletion_proba (float): The probability of simulating a deletion event, defaulting to 0.5.
            insertion_proba (float): The probability of simulating an insertion event, defaulting to 0.5.
            n_ratio (float): The ratio of 'N' bases to introduce as noise into sequences, defaulting to 0.02.
            n_proba (float): The probability of simulating 'N'.
            max_sequence_length (int): The maximum length of sequences to simulate, defaulting to 512.
            mutation_model (MutationModel): The mutation model to use for sequence mutation simulation, defaulting to S5F.
            custom_mutation_model_path (str): The path to a custom mutation model file, if any, defaulting to None.
            nucleotide_add_coefficient (float): The coefficient for the distribution of nucleotides added, defaulting to 210.
            nucleotide_remove_coefficient (float): The coefficient for the distribution of nucleotides removed, defaulting to 310.
            nucleotide_add_after_remove_coefficient (float): The coefficient for the distribution of nucleotides added after removal, defaulting to 50.
            random_sequence_add_proba (float): The probability of adding a random sequence during augmentation, defaulting to 1.
            single_base_stream_proba (float): The probability of adding a single base stream during augmentation, defaulting to 0.
            duplicate_leading_proba (float): The probability of duplicating the leading base during augmentation, defaulting to 0.
            random_allele_proba (float): The probability of adding a random allele during augmentation, defaulting to 0.
            corrupt_proba (float): The probability of corrupting the sequence from the start, defaulting to 0.7.
            corrupt_events_proba (list): The probability for each of the possible corruption events: Add, Remove, Remove and Add, defaulting to [0.5, 0.5, 0] respectively.
            short_d_length (int): The minimum length required from the D allele to not be tagged as "Short-D", defaulting to 5.
            kappa_lambda_ratio (float): The ratio of kappa to lambda light chains, defaulting to 0.5.
            save_mutations_record (bool): Whether to save the record of mutations in the sequence, defaulting to False.
            save_ns_record (bool): Whether to save the record of 'N' bases in the sequence, defaulting to False.
            productive (bool): Whether to generate a productive sequence (VJ in frame and no stop codons), defaulting to False.
            substitution_probability (float) : if this value is above 0, when applying mutation with the given min and max mutation rate, there will be a probability equals to this value that instead of a mutation a random substitution will be applied
        """
    min_mutation_rate: float = 0.003
    max_mutation_rate: float = 0.25
    simulate_indels: float = 0.2
    max_indels:int = 5
    deletion_proba: float = 0.5
    insertion_proba: float = 0.5
    n_ratio: float = 0.02
    n_proba: float = 1.0
    max_sequence_length: int = 576
    mutation_model: MutationModel = S5F
    custom_mutation_model_path: str = None
    nucleotide_add_coefficient: float = 210
    nucleotide_remove_coefficient: float = 310
    nucleotide_add_after_remove_coefficient: float = 50
    random_sequence_add_proba: float = 1
    single_base_stream_proba: float = 0
    duplicate_leading_proba: float = 0
    random_allele_proba: float = 0
    corrupt_proba: float = 0.7
    corrupt_events_proba: list = field(default_factory=lambda: [0.5,0.5,0])
    short_d_length: int = 5
    kappa_lambda_ratio: float = 0.5
    save_mutations_record: bool = True
    save_ns_record: bool = True
    save_corruption_record: bool = False
    productive: bool = False
    substitution_probability : float = 0

class SequenceAugmentorBase(ABC):
    """
       An abstract base class for sequence augmentors, providing methods and attributes to simulate and augment immunoglobulin sequences.

       This class includes functionalities for simulating mutations, introducing noise (such as indels and base substitutions), and managing allele information for V, D, and J gene segments.

       Attributes:
           dataconfig (DataConfig): Configuration data for sequence simulation.
           min_mutation_rate (float): The minimum mutation rate for simulating sequence mutations.
           max_mutation_rate (float): The maximum mutation rate for simulating sequence mutations.
           mutation_model (MutationModel): The model used for simulating mutations within sequences.
           n_ratio (float): The ratio of 'N' bases to introduce as noise into sequences.
           corrupt_proba (float): The probability of introducing corruption events at the beginning of sequences.
           simulate_indels (bool): Flag indicating whether to simulate insertions and deletions.
           max_indels (int): The maximum number of indels to simulate.
           deletion_proba (float): The probability of simulating a deletion event.
           insertion_proba (float): The probability of simulating an insertion event.
           save_mutations_record (bool): Flag indicating whether to save the record of mutations applied to sequences.
           save_ns_record (bool): Flag indicating whether to save the record of 'N' bases introduced as noise.
           productive (bool): Wheter to ensure the sequence is productive (No stop codons and mutation in cdr3 anchors), defaulting to false.
            
       Args:
           dataconfig (DataConfig): The data configuration object containing settings for sequence simulation.
           args (SequenceAugmentorArguments): The arguments object containing simulation parameters.
       """
    alleles_used = None

    def __init__(self, dataconfig: DataConfig, args: SequenceAugmentorArguments = SequenceAugmentorArguments()):

        # Sequence Generation Parameters
        self.dataconfig = dataconfig
        self.min_mutation_rate = args.min_mutation_rate
        self.max_mutation_rate = args.max_mutation_rate
        if args.custom_mutation_model_path is not None:
            self.mutation_model = args.mutation_model(self.min_mutation_rate, self.max_mutation_rate,
                                                      custom_model=args.custom_mutation_model_path, 
                                                      productive=args.productive)
        else:
            self.mutation_model = args.mutation_model(self.min_mutation_rate, self.max_mutation_rate, 
                                                      productive=args.productive)

        # Noising Parameters
        self.n_ratio = args.n_ratio
        self.n_proba = args.n_proba
        self.corrupt_proba = args.corrupt_proba
        self.corrupt_events_proba = args.corrupt_events_proba
        self.nucleotide_add_distribution = st.beta(2, 3)
        self.nucleotide_remove_distribution = st.beta(2, 3)
        self.nucleotide_add_after_remove_distribution = st.beta(1, 3)

        self.nucleotide_add_coef = args.nucleotide_add_coefficient
        self.nucleotide_remove_coef = args.nucleotide_remove_coefficient
        self.nucleotide_add_after_remove_coef = args.nucleotide_add_after_remove_coefficient

        self.random_sequence_add_proba = args.random_sequence_add_proba
        self.single_base_stream_proba = args.single_base_stream_proba
        self.duplicate_leading_proba = args.duplicate_leading_proba
        self.random_allele_proba = args.random_allele_proba

        self.simulate_indels = args.simulate_indels
        self.max_indels = args.max_indels
        self.deletion_proba = args.deletion_proba
        self.insertion_proba = args.insertion_proba
        # Class Misc
        self.chain_type = None
        self.save_mutations_record = args.save_mutations_record
        self.save_corruption_record = args.save_corruption_record
        self.save_ns_record = args.save_ns_record
        self.v_alleles = sorted([i for j in self.dataconfig.v_alleles for i in self.dataconfig.v_alleles[j]],
                                key=lambda x: x.name)

        self.j_alleles = sorted([i for j in self.dataconfig.j_alleles for i in self.dataconfig.j_alleles[j]],
                                key=lambda x: x.name)
        self.v_dict = {i.name: i.ungapped_seq.upper() for i in self.v_alleles}
        self.j_dict = {i.name: i.ungapped_seq.upper() for i in self.j_alleles}
        self.max_v_length = max(map(lambda x: len(x.ungapped_seq), self.v_alleles))
        self.max_sequence_length = args.max_sequence_length
        self.productive = args.productive
        # productivity parameters
        self.v_anchor = {i.name: i.anchor for i in self.v_alleles}
        self.j_anchor = {i.name: i.anchor for i in self.j_alleles}
        self.j_frame = {i.name: i.frame for i in self.j_alleles}
        # Loading Routines
        self.load_correction_maps()

    # Loading Routines
    @abstractmethod
    def load_correction_maps(self):
        """
            Abstract method to load allele correction maps, which are used to adjust allele choices based on observed mutations and indels.

            Implementations should load correction maps specific to the sequence type being augmented (e.g., heavy or light chain sequences).
        """
        pass

    def get_original_index(self, corrupted_index, simulated):
        """
        Calculates the original index in the uncorrupted sequence corresponding to a given index in the corrupted sequence.

        Args:
            corrupted_index (int): The index in the corrupted sequence.
            simulated (dict): A dictionary containing simulation metadata, including details of corruption events and amounts.

        Returns:
            int or None: The original index in the uncorrupted sequence, or None if the index cannot be mapped back due to the nature of the corruption.
        """
        corruption_event = simulated['corruption_event']
        amount_added = simulated['corruption_add_amount']
        amount_removed = simulated['corruption_remove_amount']

        if corruption_event == 'add':
            # If characters were added at the beginning, the original index is shifted to the right
            return corrupted_index - amount_added if corrupted_index >= amount_added else None
        elif corruption_event == 'remove':
            # If characters were removed from the beginning, the original index is shifted to the left
            return corrupted_index + amount_removed
        elif corruption_event == 'remove_before_add':
            # Combined effect of removal and addition
            adjusted_index = corrupted_index + amount_removed
            return adjusted_index - amount_added if adjusted_index >= amount_added else None
        else:
            # If no corruption or unknown corruption event, assume the index is unchanged
            return corrupted_index

    def get_allele_spesific_n_positions(self, simulated, n_positions, allele):
        """
                Introduces 'N' bases into the simulated sequence based on the configured 'N' ratio, simulating sequencing errors or low-quality bases.

                Args:
                    simulated (dict): A dictionary containing the simulated sequence and metadata.
        """
        return [i for i in n_positions if
                simulated[f'{allele}_sequence_start'] <= i <= simulated[f'{allele}_sequence_end']]

    @abstractmethod
    def attribute_position(self, simulated, position):
        """
        This method will take a simulated sequence and a position and return the name of the region it is in
        :param simulated:
        :param position:
        :return:
        """
        pass

    # Noise Introducing Methods
    def insert_Ns(self, simulated):
        """
        Introduces 'N' bases into the simulated sequence based on the configured 'N' ratio, simulating sequencing errors or low-quality bases.

                Args:
                    simulated (dict): A dictionary containing the simulated sequence and metadata.
        """
        sequence = simulated['sequence']
        # Calculate how many Ns we should insert
        num_replacements = int(len(sequence) * self.n_ratio)
        nucleotides_list = list(sequence)

        for _ in range(num_replacements):
            # Get random position in the sequence
            index = random.randint(0, len(nucleotides_list) - 1)
            # Log the N insertion event in the sequence metadata
            simulated['Ns'][index] = nucleotides_list[index] + '>' + 'N'
            # Make the N insertion
            nucleotides_list[index] = "N"

        # Concatenate the list back into a string
        simulated['sequence'] = ''.join(nucleotides_list)

        # Sort the N's insertion log
        simulated['Ns'] = {pos: simulated['Ns'][pos] for pos in sorted(simulated['Ns'])}

        # Get All the N Poistion Inserted to the V Allele
        v_allele_n_positions = self.get_allele_spesific_n_positions(simulated, list(simulated['Ns']), 'v')
        # Get Original Positions For Correction
        v_allele_n_positions_original = [self.get_original_index(i, simulated) for i in v_allele_n_positions]

        indistinguishable_v_alleles = self.v_n_ambiguity_comparer.find_indistinguishable_alleles(
            simulated['v_call'][0], v_allele_n_positions_original)
        simulated['v_call'] = simulated['v_call'] + list(
            set(indistinguishable_v_alleles) - set(simulated['v_call']))

        # Adjust mutation rate to account for added "N's" by adding the ratio of N's added to the mutation rate
        # because we randomly sample position we might get a ratio of N's less than what was defined by the n_ratio
        # property, thus we recalculate the actual inserted ration of N's
        # Also make sure we only count the N's the were inserted in the sequence and not in the added slack
        # as we might consider all the slack as noise in general
        Ns_in_pure_sequence = [i for i in simulated['Ns'] if
                               simulated['v_sequence_start'] <= i <= simulated['j_sequence_end']]
        pure_sequence_length = simulated['j_sequence_end'] - simulated['v_sequence_start']

        # exception handling
        if pure_sequence_length == 0:
            # handle division by zero
            simulated_n_ratio = 0
        else:
            simulated_n_ratio = len(Ns_in_pure_sequence) / pure_sequence_length
        simulated['mutation_rate'] += simulated_n_ratio

    def remove_event(self, simulated, autocorrect=True):
        """
        Simulates the removal of a portion of the sequence from the beginning, adjusting metadata accordingly.

                Args:
                    simulated (dict): A dictionary containing the simulated sequence and metadata.
                    autocorrect (bool): Whether to automatically correct allele choices based on the removal event.
        """
        # Update Simulation Metadata
        simulated['corruption_event'] = 'remove'
        # Sample how much to remove
        v_length = simulated['v_sequence_end'] - simulated['v_sequence_start']
        amount_to_remove = self._sample_nucleotide_remove_distribution(v_length)
        # remove from the start of the sequence the sampled amount
        # log removal
        simulated['corruption_removed_section'] = simulated['sequence'][:amount_to_remove]
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
        simulated['d_sequence_start'] -= amount_to_remove
        simulated['d_sequence_end'] -= amount_to_remove
        simulated['j_sequence_start'] -= amount_to_remove
        simulated['j_sequence_end'] -= amount_to_remove
        simulated['junction_sequence_start'] -= amount_to_remove
        simulated['junction_sequence_end'] -= amount_to_remove
        simulated['v_germline_start'] += amount_to_remove

        # Correction - Add All V Alleles That Cant be Distinguished Based on the Amount Cut from the V Allele
        if autocorrect:
            self.correct_for_v_start_cut(simulated)

    def add_event(self, simulated, amount=None):
        """
        Simulates the addition of bases to the beginning of the sequence, adjusting metadata accordingly.

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
        modified_sequence = method(amount_to_add,simulated)
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
        simulated['d_sequence_start'] += amount_to_add
        simulated['d_sequence_end'] += amount_to_add
        simulated['j_sequence_start'] += amount_to_add
        simulated['j_sequence_end'] += amount_to_add
        simulated['junction_sequence_start'] += amount_to_add
        simulated['junction_sequence_end'] += amount_to_add

    def remove_before_add_event(self, simulated):
        # ----- REMOVE PART -----#
        # Sample how much to remove
        self.remove_event(simulated, autocorrect=False)  # Dont Correct For Removed Section YET! First Lets Add

        # ----- ADD PART -----#
        # Sample how much to add after removal occurred
        amount_to_add = self._sample_nucleotide_add_after_remove_distribution()
        self.add_event(simulated, amount=amount_to_add)

        # Update Simulation Metadata
        simulated['corruption_event'] = 'remove_before_add'

        # Check If The Addition Recreated Some of the Removed Section by Chance
        self.correct_for_v_start_add(simulated)
        # Correct for the removed section of the V now that we have check for the reconstruction by chance
        self.correct_for_v_start_cut(simulated)

    def corrupt_sequence_beginning(self, simulated):
        # Sample Corruption Event
        event = self.sample_random_event()
        if event is Event.Remove:
            self.remove_event(simulated)
        elif event is Event.Add:
            self.add_event(simulated)
        elif event is Event.Remove_Before_Add:
            self.remove_before_add_event(simulated)
        else:
            raise ValueError(f'Unknown Corruption Event {event}')

    # Different V Start "Add" Event Scenarios
    @staticmethod
    def random_nucleotides(amount, simulated):
        random_seq = ''.join(random.choices(['A', 'T', 'C', 'G'], k=amount))
        simulated['corruption_added_section'] = random_seq
        return random_seq + simulated['sequence']

    @staticmethod
    def duplicate_leading(amount, simulated):
        sequence = simulated['sequence']
        cap = amount if amount < len(sequence) else len(sequence) - 1
        simulated['corruption_added_section'] = sequence[:cap]
        return sequence[:cap] + sequence

    def random_allele_section(self, amount, simulated):
        random_allele = random.choice(self.v_alleles).ungapped_seq.upper()
        cap = amount if amount < len(random_allele) else len(random_allele) - 1
        simulated['corruption_added_section'] = random_allele[:cap]
        return random_allele[:cap] + simulated['sequence']

    @staticmethod
    def single_base_stream(amount, simulated):
        random_base = random.choice(['A', 'T', 'G', 'C', 'N']) * amount
        simulated['corruption_added_section'] = random_base
        return random_base + simulated['sequence']

    # Correction Functions
    def correct_for_v_end_cut(self, simulated):
        """
        Corrects the V allele choices based on the trimming of the sequence's end, ensuring the chosen alleles remain consistent with the modified sequence.

                Args:
                    simulated (dict): A dictionary containing the simulated sequence and metadata, including the current V allele and its trimmed length.
        """
        sampled_v = simulated['v_call'][0]
        valid_length = min(simulated['v_trim_3'], self.max_v_end_correction_map_value)
        equivalent_alleles = self.v_end_allele_correction_map[sampled_v][valid_length]
        simulated['v_call'] = simulated['v_call'] + list(set(equivalent_alleles) - set(simulated['v_call']))

    def correct_for_v_start_cut(self, simulated):
        """
       Adjusts the V allele choices in response to the sequence's start being trimmed, maintaining consistency between the chosen alleles and the altered sequence.

               Args:
                   simulated (dict): A dictionary containing the simulated sequence and metadata, particularly details about the V allele and the amount trimmed from the start.
       """
        removed = simulated['corruption_remove_amount']
        equivalent_alleles = self.v_start_allele_correction_map[simulated['v_call'][0]][
            min(removed, self.max_v_start_correction_map_value)]
        simulated['v_call'] = simulated['v_call'] + list(set(equivalent_alleles) - set(simulated['v_call']))

    def correct_for_v_start_add(self, simulated):
        """
        Evaluates and corrects the V allele start position after an addition event, ensuring the V allele start position remains accurate despite sequence modifications.

        Args:
            simulated (dict): A dictionary containing the simulated sequence and metadata, including information about addition events and their impact on the V allele start position.
        """

        if simulated['corruption_event'] == 'remove_before_add':
            amount_added = simulated['corruption_add_amount']
            amount_removed = simulated['corruption_remove_amount']

            add_section = simulated['sequence'][:amount_added]
            removed_section = self.v_dict[simulated['v_call'][0]][:amount_removed]

            min_length = min(len(add_section), len(removed_section))
            to_adjust = 0
            for i in range(1, min_length + 1):
                # Compare characters from the end
                if add_section[-i] == removed_section[-i]:
                    to_adjust += 1
                else:  # Mismatch found, halt the iteration
                    break

            # Adjust V Start
            simulated['v_sequence_start'] -= to_adjust
            simulated['v_germline_start'] -= to_adjust
            # Adjust Removed Amount
            simulated['corruption_remove_amount'] -= to_adjust

    def correct_mutation_rate(self, simulated):
        """
        Recalculates and updates the mutation rate in the simulated sequence metadata to reflect the actual mutation rate after sequence modifications such as insertions or deletions.

        Args:
            simulated (dict): A dictionary containing the simulated sequence and metadata, including the current mutation log and sequence length.
        """

        # recalculate the mutation ratio based on the current mutation log
        n_mutations = len(simulated['mutations'])
        # the actual sequence length without any noise at the start or at the end
        sequence_length = simulated['j_sequence_end'] - simulated['v_sequence_start']
        mutation_ratio = n_mutations / sequence_length
        # update the mutation_rate in the simulation log
        simulated['mutation_rate'] = mutation_ratio

    def validate_sequence_length_after_addition(self, sequence, added):
        """
                Ensures that the sequence length does not exceed the maximum allowed length after an addition event, trimming excess bases if necessary.

                Args:
                    sequence (str): The sequence after the addition event.
                    added (int): The number of bases added to the sequence.

                Returns:
                    tuple: A tuple containing the validated sequence and the adjusted number of bases added.
        """
        if len(sequence) > self.max_sequence_length:
            # Calculate how much did we over add to based on our maximum sequence length
            slack = len(sequence) - self.max_sequence_length
            # remove the slack from the begging and update the added variable to match the final added amount
            return sequence[slack:], added - slack
        else:
            return sequence, added

    def fix_v_position_after_trimming_index_ambiguity(self, simulation):
        """
                Resolves any ambiguities in V allele positions resulting from sequence trimming, ensuring accurate representation of V allele boundaries.

                Args:
                    simulation (dict): A dictionary containing the simulated sequence and metadata, including V allele positions and trimming details.
        """
        # Extract Current V Metadata
        v_start, v_end = simulation['v_sequence_start'], simulation['v_sequence_end']
        v_germline_start, v_germline_end = simulation['v_germline_start'], simulation['v_germline_end']
        v_allele_remainder = simulation['sequence'][v_start:v_end]
        v_allele_ref = self.v_dict[simulation['v_call'][0]]

        # Get the junction inserted after trimming to the sequence
        junction_3 = simulation['sequence'][v_end:simulation['d_sequence_start']]

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
                v_trim_3 -= 1
            else:  # if the continuous streak is broken or non-existent break!
                break

        simulation['v_sequence_end'] = v_end
        simulation['v_germline_end'] = v_germline_end
        simulation['v_trim_3'] = v_trim_3

    def fix_j_position_after_trimming_index_ambiguity(self, simulation):
        """
            Addresses ambiguities in J allele positions due to sequence trimming, ensuring J allele boundaries are accurately represented.

                Args:
                    simulation (dict): A dictionary containing the simulated sequence and metadata, including J allele positions and trimming details.
        """
        # Extract Current J Metadata
        j_start, j_end = simulation['j_sequence_start'], simulation['j_sequence_end']
        j_germline_start, j_germline_end = simulation['j_germline_start'], simulation['j_germline_end']
        j_allele_remainder = simulation['sequence'][j_start:j_end]
        j_allele_ref = self.j_dict[simulation['j_call'][0]]

        # Get the junction inserted after trimming to the sequence
        junction_5 = simulation['sequence'][simulation['d_sequence_end']:j_start]

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
                j_trim_5 -= 1
            else:  # if the continuous streak is broken or non-existent break!
                break

        simulation['j_sequence_start'] = j_start
        simulation['j_germline_start'] = j_germline_start
        simulation['j_trim_5'] = j_trim_5

    @abstractmethod
    def distill_mutation_rate(self):
        pass

    # Sequence Simulation
    @abstractmethod
    def simulate_sequence(self):
        """
        Abstract method to simulate an immunoglobulin sequence, incorporating V(D)J recombination, mutations, and potentially other forms of sequence variation.

        Implementations should generate a simulated sequence along with relevant metadata, such as the positions of gene segments and mutations.

        Returns:
            dict: A dictionary containing the simulated sequence and associated metadata.
        """
        pass

    # Interface
    def process_before_return(self, simulated):
        """
        Performs final processing on the simulated sequence before returning it, including formatting allele lists and encoding mutation and 'N' base records.

        Args:
            simulated (dict): A dictionary containing the simulated sequence and metadata.
        """
        # Convert Allele From List to Comma Seperated String
        for gene in self.alleles_used:
            simulated[f'{gene}_call'] = ','.join(simulated[f'{gene}_call'])

        # convert mutations and N log to base64 for proper tabular preservation
        if not self.save_ns_record:
        #     simulated['Ns'] = base64.b64encode(str(simulated['Ns']).encode('ascii'))
        # else:
            simulated.pop('Ns')

        if not self.save_mutations_record:
        #     simulated['mutations'] = base64.b64encode(str(simulated['mutations']).encode('ascii'))
        # else:
            simulated.pop('mutations')

        if not self.save_corruption_record:
            simulated.pop('corruption_removed_section')
            simulated.pop('corruption_added_section')

    @abstractmethod
    def simulate_augmented_sequence(self):
        """
                Abstract method to simulate and augment an immunoglobulin sequence, incorporating additional forms of sequence variation beyond basic V(D)J recombination and mutations.

                Implementations might include more complex forms of noise or variation, such as indels or base modifications, and should adjust the sequence and metadata accordingly.

                Returns:
                    dict: A dictionary containing the augmented sequence and associated metadata.
        """
        pass

    # Misc Methods
    def _sample_nucleotide_add_distribution(self):
        """
            Samples from a predefined distribution to determine the number of nucleotides to add during a sequence augmentation event.

                Returns:
                    int: The number of nucleotides to be added to the sequence.
        """
        sample = (self.nucleotide_add_coef * self.nucleotide_add_distribution.rvs(size=1)).astype(int)
        return sample.item()

    def _sample_nucleotide_remove_distribution(self, v_length):
        """
                Samples from a predefined distribution to decide the number of nucleotides to remove, with a maximum limit based on the V segment length.

                Args:
                    v_length (int): The length of the V segment in the sequence.

                Returns:
                    int: The number of nucleotides to be removed from the sequence.
        """
        # Sample amount based on predefined distribution
        sample = (self.nucleotide_remove_coef * self.nucleotide_remove_distribution.rvs(size=1)).astype(int).item()
        # make sure that no matter how much we get from sampling the predefined distribution we wont get a value
        # larger than the total v length we have in our sequence
        sample = min(sample, v_length)
        return sample

    def _sample_nucleotide_add_after_remove_distribution(self):
        """
                Samples from a predefined distribution to determine the number of nucleotides to add after a removal event during sequence augmentation.

                Returns:
                    int: The number of nucleotides to be added to the sequence following a removal event.
        """
        sample = (self.nucleotide_add_after_remove_coef * self.nucleotide_add_after_remove_distribution.rvs(
            size=1)).astype(int)
        return sample.item()

    def perform_corruption(self):
        """
                Determines whether a corruption event should be performed on the sequence based on the configured probability.

                Returns:
                    bool: True if a corruption event should be performed, False otherwise.
        """
        return bool(np.random.binomial(1, self.corrupt_proba))

    #@staticmethod
    def sample_random_event(self):
        """
                Randomly selects a corruption event type (e.g., add, remove, remove before add) based on predefined probabilities.

                Returns:
                    Event: An enumerated value representing the selected corruption event type.
        """
        return np.random.choice([Event.Remove, Event.Add, Event.Remove_Before_Add], size=1, p=self.corrupt_events_proba).item()

    def _sample_corruption_add_method(self):
        """
               Randomly selects a method for adding nucleotides to the sequence during an augmentation event based on predefined probabilities.

               Returns:
                   function: A reference to the function that will be used to add nucleotides to the sequence.
       """
        method = random.choices([self.random_nucleotides, self.single_base_stream, self.duplicate_leading,
                                 self.random_allele_section],
                                weights=[self.random_sequence_add_proba, self.single_base_stream_proba,
                                         self.duplicate_leading_proba, self.random_allele_proba], k=1)[0]
        return method

    @property
    def columns(self):
        """
               Returns a list of column names that would be present in the output of the simulate_augmented_sequence method.

               Returns:
                   list: A list of strings representing the column names.
       """
        return list(self.simulate_augmented_sequence())
