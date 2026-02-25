import random

from ..container.SimulationContainer import SimulationContainer
from ..pipeline.plot_parameters import CORRUPTION_STEP_BOX_COLOR
from ..simulation import Event
from ..steps.StepBase import AugmentationStep
from ..utilities import translate
from ..dataconfig import DataConfig


class CorruptSequenceBeginning(AugmentationStep):
    """
    Simulates 5' end corruption of sequences, modeling sequencing artifacts.

    This step can add, remove, or remove-then-add nucleotides at the beginning
    of sequences to simulate real-world sequencing issues like PCR artifacts
    and primer binding.

    Note:
        This step handles corruption events only. For enforcing sequencing
        platform read length limits, use EnforceSequenceLength after this step.

    Args:
        probability: Probability that corruption occurs (0.0-1.0). Default: 0.7
        event_weights: Weights for [add, remove, remove_then_add] events. Default: (0.4, 0.4, 0.2)
        nucleotide_add_coefficient: Scale factor for nucleotide addition distribution. Default: 210
        nucleotide_remove_coefficient: Scale factor for nucleotide removal distribution. Default: 310
        nucleotide_add_after_remove_coefficient: Scale factor for add-after-remove. Default: 50
        random_sequence_add_probability: Weight for random nucleotide addition method. Default: 1.0
        single_base_stream_probability: Weight for single base repeat method. Default: 0.0
        duplicate_leading_probability: Weight for duplicating leading bases method. Default: 0.0
        random_allele_probability: Weight for random allele section method. Default: 0.0

    Example:
        # Use defaults
        step = CorruptSequenceBeginning()

        # Customize specific parameters
        step = CorruptSequenceBeginning(probability=0.9, nucleotide_add_coefficient=150)

        # For length enforcement, add EnforceSequenceLength after this step:
        # pipeline = AugmentationPipeline([
        #     ...,
        #     CorruptSequenceBeginning(),
        #     EnforceSequenceLength(max_length=576),
        #     ...
        # ])
    """

    def __init__(
        self,
        *,
        probability: float = 0.7,
        event_weights: tuple = (0.4, 0.4, 0.2),
        nucleotide_add_coefficient: int = 210,
        nucleotide_remove_coefficient: int = 310,
        nucleotide_add_after_remove_coefficient: int = 50,
        random_sequence_add_probability: float = 1.0,
        single_base_stream_probability: float = 0.0,
        duplicate_leading_probability: float = 0.0,
        random_allele_probability: float = 0.0,
    ):
        super().__init__()

        self.corrupt_events_proba = list(event_weights)
        self.random_sequence_add_proba = random_sequence_add_probability
        self.single_base_stream_proba = single_base_stream_probability
        self.duplicate_leading_proba = duplicate_leading_probability
        self.random_allele_proba = random_allele_probability
        self.corruption_probability = probability

        self.nucleotide_add_beta_params = (2, 3)
        self.nucleotide_remove_beta_params = (2, 3)
        self.nucleotide_add_after_remove_beta_params = (1, 3)

        self.nucleotide_add_coef = nucleotide_add_coefficient
        self.nucleotide_remove_coef = nucleotide_remove_coefficient
        self.nucleotide_add_after_remove_coef = nucleotide_add_after_remove_coefficient

        # Lazy-initialized attributes (populated on first access after config is bound)
        self._v_start_allele_correction_map = None
        self._max_v_start_correction_map_value = None
        self._v_alleles = None
        self._j_alleles = None
        self._v_dict = None

    def _ensure_config_loaded(self):
        """Initialize dataconfig-dependent attributes on first use."""
        if self._v_start_allele_correction_map is None:
            self._v_start_allele_correction_map = self.dataconfig.correction_maps['V_5_TRIM_SIMILARITY_MAP']
            self._max_v_start_correction_map_value = max(
                self._v_start_allele_correction_map[list(self._v_start_allele_correction_map)[0]])

            self._v_alleles = sorted(
                [i for j in self.dataconfig.v_alleles for i in self.dataconfig.v_alleles[j]],
                key=lambda x: x.name
            )

            self._j_alleles = sorted(
                [i for j in self.dataconfig.j_alleles for i in self.dataconfig.j_alleles[j]],
                key=lambda x: x.name
            )
            self._v_dict = {i.name: i.ungapped_seq.upper() for i in self._v_alleles}

    @property
    def v_start_allele_correction_map(self):
        self._ensure_config_loaded()
        return self._v_start_allele_correction_map

    @property
    def max_v_start_correction_map_value(self):
        self._ensure_config_loaded()
        return self._max_v_start_correction_map_value

    @property
    def v_alleles(self):
        self._ensure_config_loaded()
        return self._v_alleles

    @property
    def j_alleles(self):
        self._ensure_config_loaded()
        return self._j_alleles

    @property
    def v_dict(self):
        self._ensure_config_loaded()
        return self._v_dict

    @staticmethod
    def random_nucleotides(amount, container: SimulationContainer):
        random_seq = ''.join(random.choices(['A', 'T', 'C', 'G'], k=amount))
        container.corruption_added_section = random_seq
        return random_seq + container.sequence

    @staticmethod
    def duplicate_leading(amount, container: SimulationContainer):
        sequence = container.sequence
        cap = amount if amount < len(sequence) else len(sequence) - 1
        container.corruption_added_section = sequence[:cap]
        return sequence[:cap] + sequence

    def random_allele_section(self, amount, container: SimulationContainer):
        random_allele = random.choice(self.v_alleles).ungapped_seq.upper()
        cap = amount if amount < len(random_allele) else len(random_allele) - 1
        container.corruption_added_section = random_allele[:cap]
        return random_allele[:cap] + container.sequence

    @staticmethod
    def single_base_stream(amount, container: SimulationContainer):
        random_base = random.choice(['A', 'T', 'G', 'C', 'N']) * amount
        container.corruption_added_section = random_base
        return random_base + container.sequence

    def _sample_nucleotide_add_distribution(self):
        """
            Samples from a predefined distribution to determine the number of nucleotides to add during a sequence augmentation event.

                Returns:
                    int: The number of nucleotides to be added to the sequence.
        """
        sample = int(self.nucleotide_add_coef * random.betavariate(*self.nucleotide_add_beta_params))
        return max(1, sample)

    def _sample_nucleotide_remove_distribution(self, v_length):
        """
                Samples from a predefined distribution to decide the number of nucleotides to remove, with a maximum limit based on the V segment length.

                Args:
                    v_length (int): The length of the V segment in the sequence.

                Returns:
                    int: The number of nucleotides to be removed from the sequence.
        """
        # Sample amount based on predefined distribution
        sample = int(self.nucleotide_remove_coef * random.betavariate(*self.nucleotide_remove_beta_params))
        # make sure that no matter how much we get from sampling the predefined distribution we wont get a value
        # larger than the total v length we have in our sequence
        sample = min(sample, v_length)
        return max(1, sample)

    def _sample_nucleotide_add_after_remove_distribution(self):
        """
                Samples from a predefined distribution to determine the number of nucleotides to add after a removal event during sequence augmentation.

                Returns:
                    int: The number of nucleotides to be added to the sequence following a removal event.
        """
        sample = int(self.nucleotide_add_after_remove_coef * random.betavariate(*self.nucleotide_add_after_remove_beta_params))
        return max(1, sample)

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

    def correct_mutation_rate(self, container: SimulationContainer):
        """
        Recalculates and updates the mutation rate in the simulated sequence metadata to reflect the actual mutation rate after sequence modifications such as insertions or deletions.

        Args:
            container (dict): A dictionary containing the simulated sequence and metadata, including the current mutation log and sequence length.
        """

        # recalculate the mutation ratio based on the current mutation log
        n_mutations = len(container.mutations)
        mutation_ratio = n_mutations / len(container.sequence)
        # update the mutation_rate in the simulation log
        container.mutation_rate = mutation_ratio

    def correct_for_v_start_cut(self, container: SimulationContainer):
        """
       Adjusts the V allele choices in response to the sequence's start being trimmed, maintaining consistency between the chosen alleles and the altered sequence.

               Args:
                   container (dict): A dictionary containing the simulated sequence and metadata, particularly details about the V allele and the amount trimmed from the start.
       """
        removed = container.corruption_remove_amount
        equivalent_alleles = self.v_start_allele_correction_map[container.v_call[0]][
            min(removed, self.max_v_start_correction_map_value)]
        container.v_call = container.v_call + list(set(equivalent_alleles) - set(container.v_call))

    def correct_for_v_start_add(self, container: SimulationContainer):
        """
        Evaluates and corrects the V allele start position after an addition event, ensuring the V allele start position remains accurate despite sequence modifications.

        Args:
            container (SimulationContainer): A dictionary containing the simulated sequence and metadata, including information about addition events and their impact on the V allele start position.
        """

        if container.corruption_event == 'remove_before_add':
            amount_added = container.corruption_add_amount
            amount_removed = container.corruption_remove_amount

            add_section = container.sequence[:amount_added]
            removed_section = self.v_dict[container.v_call[0]][:amount_removed]

            min_length = min(len(add_section), len(removed_section))
            to_adjust = 0
            for i in range(1, min_length + 1):
                # Compare characters from the end
                if add_section[-i] == removed_section[-i]:
                    to_adjust += 1
                else:  # Mismatch found, halt the iteration
                    break

            # Adjust V Start
            container.v_sequence_start -= to_adjust
            container.v_germline_start -= to_adjust
            # Adjust Removed Amount
            container.corruption_remove_amount -= to_adjust

    @property
    def perform_corruption(self):
        """
                Determines whether a corruption event should be performed on the sequence based on the configured probability.

                Returns:
                    bool: True if a corruption event should be performed, False otherwise.
        """
        return random.random() < self.corruption_probability

    def sample_random_event(self):
        """
                Randomly selects a corruption event type (e.g., add, remove, remove before add) based on predefined probabilities.

                Returns:
                    Event: An enumerated value representing the selected corruption event type.
        """
        return random.choices(
            [Event.Remove, Event.Add, Event.Remove_Before_Add],
            weights=self.corrupt_events_proba, k=1
        )[0]

    def remove_event(self, container, autocorrect=True):
        """
        Simulates the removal of a portion of the sequence from the beginning, adjusting metadata accordingly.

                Args:
                    simulated (dict): A dictionary containing the simulated sequence and metadata.
                    autocorrect (bool): Whether to automatically correct allele choices based on the removal event.
        """
        # Update Simulation Metadata
        container.corruption_event = 'remove'
        # Sample how much to remove
        v_length = container.v_sequence_end - container.v_sequence_start
        amount_to_remove = self._sample_nucleotide_remove_distribution(v_length)
        # remove from the start of the sequence the sampled amount
        # log removal
        container.corruption_removed_section = container.sequence[:amount_to_remove]
        container.sequence = container.sequence[amount_to_remove:]
        # Update Simulation Metadata
        container.corruption_remove_amount = amount_to_remove
        # Update mutation log, remove the mutations that were removed with the remove event
        # and while updating the mutations log correct the position of the mutations accordingly
        container.mutations = {i - amount_to_remove: j for i, j in container.mutations.items() if
                               i >= amount_to_remove}

        # Adjust Start/End Position Accordingly
        container.v_sequence_start = 0
        container.v_sequence_end -= amount_to_remove
        container.d_sequence_start -= amount_to_remove
        container.d_sequence_end -= amount_to_remove
        container.j_sequence_start -= amount_to_remove
        container.j_sequence_end -= amount_to_remove
        container.junction_sequence_start -= amount_to_remove
        container.junction_sequence_end -= amount_to_remove
        container.v_germline_start += amount_to_remove

        # Adjust mutation rate (after positions are updated so denominator is correct)
        self.correct_mutation_rate(container)

        # Correction - Add All V Alleles That Cant be Distinguished Based on the Amount Cut from the V Allele
        if autocorrect:
            self.correct_for_v_start_cut(container)

    def add_event(self, container: SimulationContainer, amount=None):
        """
        Simulates the addition of bases to the beginning of the sequence, adjusting metadata accordingly.

                Args:
                    simulated (dict): A dictionary containing the simulated sequence and metadata.
                    amount (int, optional): The number of bases to add. If not specified, the amount is sampled from a configured distribution.
        """
        # Update Simulation Metadata
        container.corruption_event = 'add'
        # Sample the Amount to Add by default, if a specific value is given via the amount variable, use it.
        amount_to_add = self._sample_nucleotide_add_distribution() if amount is None else amount
        # Sample the method by which addition will be made
        method = self._sample_corruption_add_method()
        # Modify the sequence
        modified_sequence = method(amount_to_add, container)
        # Update Simulation Sequence
        container.sequence = modified_sequence

        # Update Simulation Metadata
        container.corruption_add_amount = amount_to_add
        # Update mutation log positions
        container.mutations = {i + amount_to_add: j for i, j in container.mutations.items()}

        # Adjust Start/End Position Accordingly
        container.v_sequence_start += amount_to_add
        container.v_sequence_end += amount_to_add
        container.d_sequence_start += amount_to_add
        container.d_sequence_end += amount_to_add
        container.j_sequence_start += amount_to_add
        container.j_sequence_end += amount_to_add
        container.junction_sequence_start += amount_to_add
        container.junction_sequence_end += amount_to_add

    def remove_before_add_event(self, container: SimulationContainer):
        # ----- REMOVE PART -----#
        # Sample how much to remove
        self.remove_event(container, autocorrect=False)  # Dont Correct For Removed Section YET! First Lets Add

        # ----- ADD PART -----#
        # Sample how much to add after removal occurred
        amount_to_add = self._sample_nucleotide_add_after_remove_distribution()
        self.add_event(container, amount=amount_to_add)

        # Update Simulation Metadata
        container.corruption_event = 'remove_before_add'

        # Check If The Addition Recreated Some of the Removed Section by Chance
        self.correct_for_v_start_add(container)
        # Correct for the removed section of the V now that we have check for the reconstruction by chance
        self.correct_for_v_start_cut(container)

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

    def fix_productive_call_after_corruption_indel(self, container: SimulationContainer):
        sequence = container.sequence
        functional = False
        stop_codon = False
        note = ''
        # stop codon
        stops = ["TAG", "TAA", "TGA"]
        for x in range(container.junction_sequence_end, container.v_sequence_start, -3):
            if sequence[x - 3:x] in stops:
                stop_codon = True
        if not stop_codon:
            for x in range(container.junction_sequence_end, len(sequence), 3):
                if sequence[x:x + 3] in stops:
                    stop_codon = True
        # vj in frame
        from_j_to_start = (container.junction_sequence_end - container.v_sequence_start) % 3 == 0
        junction_length = (container.junction_sequence_end - container.junction_sequence_start) % 3 == 0
        from_v_to_start = (container.junction_sequence_start - container.v_sequence_start) % 3 == 0
        vj_in_frame = from_j_to_start and junction_length and from_v_to_start and stop_codon == False
        junction = sequence[container.junction_sequence_start:container.junction_sequence_end]
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
        container.productive = functional
        container.stop_codon = stop_codon
        container.vj_in_frame = vj_in_frame
        container.note = note

    def apply(self, container: SimulationContainer) -> None:
        if self.perform_corruption:
            self.corrupt_sequence_beginning(container)
            self.fix_productive_call_after_corruption_indel(container)

    def get_graph_node(self):
        """Generates a detailed GraphViz node representation with constructor details in an HTML-like format."""
        step_name = "Corrupt Sequence 5' End"

        # Constructing an HTML-like label using a table for detailed formatting
        label = f"""
        <TABLE BORDER="0" CELLBORDER="1" CELLSPACING="0" CELLPADDING="4">
        <TR><TD COLSPAN="2" BGCOLOR="lightsteelblue"><B>{step_name}</B></TD></TR>
        <TR><TD ALIGN="LEFT"><B>Corruption Probability</B></TD><TD ALIGN="LEFT">{self.corruption_probability}</TD></TR>
        <TR><TD ALIGN="LEFT"><B>Add Event Probability</B></TD><TD ALIGN="LEFT">{self.corrupt_events_proba[0]}</TD></TR>
        <TR><TD ALIGN="LEFT"><B>Remove Event Probability</B></TD><TD ALIGN="LEFT">{self.corrupt_events_proba[1]}</TD></TR>
        <TR><TD ALIGN="LEFT"><B>Remove Than Add Event Probability</B></TD><TD ALIGN="LEFT">{self.corrupt_events_proba[2]}</TD></TR>
        <TR><TD ALIGN="LEFT"><B>Add Type: |Random Sequence| Probability</B></TD><TD ALIGN="LEFT">{self.random_sequence_add_proba}</TD></TR>
        <TR><TD ALIGN="LEFT"><B>Add Type: |Single Base Stream | Probability</B></TD><TD ALIGN="LEFT">{self.single_base_stream_proba}</TD></TR>
        <TR><TD ALIGN="LEFT"><B>Add Type: |Duplicate Leading| Probability</B></TD><TD ALIGN="LEFT">{self.duplicate_leading_proba}</TD></TR>
        <TR><TD ALIGN="LEFT"><B>Add Type: |Random Allele Section| Probability</B></TD><TD ALIGN="LEFT">{self.random_allele_proba}</TD></TR>
        </TABLE>
        
        """


        return label, 'box', "filled,rounded", CORRUPTION_STEP_BOX_COLOR, "Helvetica", "black"