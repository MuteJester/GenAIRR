import random

import numpy as np

from GenAIRR.container.SimulationContainer import SimulationContainer
from GenAIRR.pipeline.plot_parameters import CORRUPTION_STEP_BOX_COLOR
from GenAIRR.steps.StepBase import AugmentationStep
from GenAIRR.utilities import translate
from GenAIRR.dataconfig import DataConfig

class InsertIndels(AugmentationStep):
    def __init__(self, indel_probability,max_indels,insertion_proba,deletion_proba):
        super().__init__()
        self.deletion_proba = deletion_proba
        self.insertion_proba = insertion_proba
        self.indel_probability = indel_probability
        self.max_indels = max_indels

    def valid_indel_positions(self, container: SimulationContainer):
        all_positions = set(range(0, container.j_sequence_end))
        # remove np regions
        np1_positions = set(range(container.v_sequence_end, container.d_sequence_start))
        np2_positions = set(range(container.d_sequence_end, container.j_sequence_start))
        all_positions -= np1_positions
        all_positions -= np2_positions

        # remove ns positions
        all_positions -= set(container.Ns.keys())
        # remove mutations position
        all_positions -= set(container.mutations.keys())

        return all_positions

    def apply_deletion(self, container: SimulationContainer, position):
        sequence = container.sequence
        nucleotides_list = list(sequence)

        deleted = nucleotides_list[position]

        after_deletion = nucleotides_list[:position] + nucleotides_list[position + 1:]

        # update start/end positions
        for reg in ['v_sequence_start', 'v_sequence_end', 'j_sequence_start', 'j_sequence_end', 'd_sequence_start',
                    'd_sequence_end', 'junction_sequence_start' , 'junction_sequence_end']:
            if container[reg] > position:
                container[reg] = container[reg] - 1

        container.indels[position] = 'D > ' + deleted

        # correct N's and Mutations
        corrected_mutations = {}
        for pos in container.mutations:
            if pos > position:
                corrected_mutations[pos-1] = container.mutations[pos]
            else:
                corrected_mutations[pos] = container.mutations[pos]

        container.mutations = corrected_mutations

        corrected_Ns = {}
        for pos in container.Ns:
            if pos >= position:
                corrected_Ns[pos-1] = container.Ns[pos]
            else:
                corrected_Ns[pos] = container.Ns[pos]

        container.Ns = corrected_Ns


        container.sequence = ''.join(after_deletion)
        # Update log positions for deletions
        updated_log = {}
        for log_pos, log_entry in container.indels.items():
            if log_pos >= position:
                updated_log[log_pos - 1] = log_entry
            elif log_pos < position:
                updated_log[log_pos] = log_entry

        container.indels = updated_log

    def apply_insertion(self, container: SimulationContainer, position):
        sequence = container.sequence
        nucleotides_list = list(sequence)
        random_base = random.choice(['A', 'T', 'C', 'G'])


        after_insertion = nucleotides_list[:position] + [random_base] + nucleotides_list[position:]
        # update start/end positions
        for reg in ['v_sequence_start', 'v_sequence_end', 'j_sequence_start', 'j_sequence_end', 'd_sequence_start',
                    'd_sequence_end', 'junction_sequence_start' , 'junction_sequence_end']:
            if container[reg] >= position:
                container[reg] += 1

        container.indels[position] = 'I < ' + random_base
        # correct N's and Mutations
        corrected_mutations = {}
        for pos in container.mutations:
            if pos >= position:
                corrected_mutations[pos + 1] = container.mutations[pos]
            else:
                corrected_mutations[pos] = container.mutations[pos]

        container.mutations = corrected_mutations

        corrected_Ns = {}
        for pos in container.Ns:
            if pos >= position:
                corrected_Ns[pos + 1] = container.Ns[pos]
            else:
                corrected_Ns[pos] = container.Ns[pos]

        container.Ns = corrected_Ns

        container.sequence = ''.join(after_insertion)

        # Update log positions for insertions
        updated_log = {}
        for log_pos, log_entry in container.indels.items():
            if log_pos > position:
                updated_log[log_pos + 1] = log_entry
            else:
                updated_log[log_pos] = log_entry
        container.indels = updated_log

        container.sequence = ''.join(after_insertion)

    def insert_indels(self, container: SimulationContainer):
        # get valid position for indels excluding np regions, n's and mutated positions
        valid_positions = list(self.valid_indel_positions(container))
        num_indels = np.random.randint(1, self.max_indels, size=1).item()
        num_indels = min(num_indels, len(valid_positions))
        random.shuffle(valid_positions)
        n_valid_positions = len(valid_positions)

        for idx in range(num_indels):
            indel_position = valid_positions[idx]

            # choose action 1 = insertion -1 = deletion
            action = np.random.choice([1, -1], size=1, p=[self.insertion_proba, self.deletion_proba]).item()
            if action == 1:  # insertion case

                self.apply_insertion(container, indel_position)

                for update_idx in range(idx, idx + (num_indels - idx)):
                    if valid_positions[update_idx] >= indel_position:
                        valid_positions[update_idx] += 1
            else:  # deletion case

                self.apply_deletion(container, indel_position)

                for update_idx in range(idx, idx + (num_indels - idx)):
                    if valid_positions[update_idx] > indel_position:
                        valid_positions[update_idx] -= 1

    def fix_productive_call_after_corruption_indel(self,container: SimulationContainer):
        sequence = container.sequence
        functional = False
        stop_codon = False
        note = container.note
        # stop codon
        stops = ["TAG", "TAA", "TGA"]
        for x in range(container.junction_sequence_end, container.v_sequence_start, -3):
            if sequence[x-3:x] in stops:
                stop_codon = True
        if not stop_codon:
            for x in range(container.junction_sequence_end, len(sequence), 3):
                if sequence[x:x+3] in stops:
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
        # Insert indels if the condition is met
        if bool(np.random.binomial(1, self.indel_probability)):
            self.insert_indels(container)
            self.fix_productive_call_after_corruption_indel(container)

    def get_graph_node(self):
        """Generates a detailed GraphViz node representation with constructor details in an HTML-like format."""
        step_name = "Insert Indels"

        # Constructing an HTML-like label using a table for detailed formatting
        label = f"""
        <TABLE BORDER="0" CELLBORDER="1" CELLSPACING="0" CELLPADDING="4">
        <TR><TD COLSPAN="2" BGCOLOR="lightsteelblue"><B>{step_name}</B></TD></TR>
        <TR><TD ALIGN="LEFT"><B>Insertion Probability</B></TD><TD ALIGN="LEFT">{self.indel_probability}</TD></TR>
        <TR><TD ALIGN="LEFT"><B>Max Indels</B></TD><TD ALIGN="LEFT">{self.max_indels}</TD></TR>
        <TR><TD ALIGN="LEFT"><B>Insertion Proba</B></TD><TD ALIGN="LEFT">{self.insertion_proba}</TD></TR>
        <TR><TD ALIGN="LEFT"><B>Deletion Proba</B></TD><TD ALIGN="LEFT">{self.deletion_proba}</TD></TR>
        </TABLE>
        
        """

        return label, 'box', "filled,rounded", CORRUPTION_STEP_BOX_COLOR, "Helvetica", "black"
