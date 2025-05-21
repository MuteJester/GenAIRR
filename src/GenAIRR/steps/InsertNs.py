import random

import numpy as np

from ..container.SimulationContainer import SimulationContainer
from ..pipeline.plot_parameters import CORRUPTION_STEP_BOX_COLOR
from ..steps.StepBase import AugmentationStep
from ..dataconfig import DataConfig


class InsertNs(AugmentationStep):
    def __init__(self, n_ratio,proba):
        super().__init__()
        """

        :param n_ratio: this will control how many N's to add to the sequence, n_ratio*sequence_length is the number of
        random N's that will be added to the sequence when this step is applied
        :param proba: the probability that this step is applied
        """

        self.v_n_ambiguity_comparer = self.dataconfig.correction_maps['V_N_AMBIGUITY_CORRECTION_GRAPH']
        self.n_ratio = n_ratio
        self.n_proba = proba

    def get_allele_spesific_n_positions(self, container, allele):
        """
                Introduces 'N' bases into the simulated sequence based on the configured 'N' ratio, simulating sequencing errors or low-quality bases.

                Args:
                    simulated (dict): A dictionary containing the simulated sequence and metadata.
        """
        return [i for i in list(container['Ns']) if
                container[f'{allele}_sequence_start'] <= i <= container[f'{allele}_sequence_end']]

    def insert_Ns(self, container: SimulationContainer):
        """
        Introduces 'N' bases into the simulated sequence based on the configured 'N' ratio, simulating sequencing errors or low-quality bases.

                Args:
                    simulated (dict): A dictionary containing the simulated sequence and metadata.
        """
        sequence = container['sequence']
        # Calculate how many Ns we should insert
        num_replacements = int(len(sequence) * self.n_ratio)
        nucleotides_list = list(sequence)

        for _ in range(num_replacements):
            # Get random position in the sequence
            index = random.randint(0, len(nucleotides_list) - 1)
            # Log the N insertion event in the sequence metadata
            container.add_N_insertion(index, nucleotides_list[index])
            # Make the N insertion
            nucleotides_list[index] = "N"

        # Concatenate the list back into a string
        container.sequence = ''.join(nucleotides_list)

        # Sort the N's insertion log
        container.Ns = {pos: container['Ns'][pos] for pos in sorted(container['Ns'])}

        # Get All the N Poistion Inserted to the V Allele
        v_allele_n_positions = self.get_allele_spesific_n_positions(container, 'v')
        # Get Original Positions For Correction
        v_allele_n_positions_original = [container.get_original_index(i) for i in v_allele_n_positions]

        indistinguishable_v_alleles = self.v_n_ambiguity_comparer.find_indistinguishable_alleles(
            container['v_call'][0], v_allele_n_positions_original)
        container['v_call'] = container['v_call'] + list(
            set(indistinguishable_v_alleles) - set(container['v_call']))

        # Adjust mutation rate to account for added "N's" by adding the ratio of N's added to the mutation rate
        # because we randomly sample position we might get a ratio of N's less than what was defined by the n_ratio
        # property, thus we recalculate the actual inserted ration of N's
        # Also make sure we only count the N's the were inserted in the sequence and not in the added slack
        # as we might consider all the slack as noise in general
        Ns_in_pure_sequence = [i for i in container['Ns'] if
                               container['v_sequence_start'] <= i <= container['j_sequence_end']]
        pure_sequence_length = container['j_sequence_end'] - container['v_sequence_start']

        # exception handling
        if pure_sequence_length == 0:
            # handle division by zero
            simulated_n_ratio = 0
        else:
            simulated_n_ratio = len(Ns_in_pure_sequence) / pure_sequence_length
        container['mutation_rate'] += simulated_n_ratio

    @property
    def simulate_Ns(self):
        return bool(np.random.binomial(1, self.n_proba))
    def apply(self, container: SimulationContainer) -> None:
        if self.simulate_Ns:
            self.insert_Ns(container)

    def get_graph_node(self):
        """Generates a detailed GraphViz node representation with constructor details in an HTML-like format."""
        step_name = "Insert Ns"

        # Constructing an HTML-like label using a table for detailed formatting
        label = f"""
        <TABLE BORDER="0" CELLBORDER="1" CELLSPACING="0" CELLPADDING="4">
        <TR><TD COLSPAN="2" BGCOLOR="lightsteelblue"><B>{step_name}</B></TD></TR>
        <TR><TD ALIGN="LEFT"><B>N Ratio</B></TD><TD ALIGN="LEFT">{self.n_ratio}</TD></TR>
        <TR><TD ALIGN="LEFT"><B>Probability</B></TD><TD ALIGN="LEFT">{self.n_proba}</TD></TR>
        </TABLE>
        
        """

        return label, 'box', "filled,rounded", CORRUPTION_STEP_BOX_COLOR, "Helvetica", "black"
