from GenAIRR.container.SimulationContainer import SimulationContainer
from GenAIRR.pipeline.plot_parameters import CORRECTION_STEP_HEADING_COLOR
from GenAIRR.steps.StepBase import AugmentationStep
from GenAIRR.dataconfig import DataConfig


class FixDPositionAfterTrimmingIndexAmbiguity(AugmentationStep):

    def __init__(self):
        super().__init__()
        self.d_alleles = sorted([i for j in self.dataconfig.d_alleles for i in self.dataconfig.d_alleles[j]],
                                key=lambda x: x.name)
        self.d_dict = {i.name: i.ungapped_seq.upper() for i in self.d_alleles}

    def fix_d_position_after_trimming_index_ambiguity(self, container):
        """
                Corrects the start and end positions of the D gene segment in the simulated sequence to resolve ambiguities caused by trimming events.

                Args:
                    simulation (dict): Dictionary containing the simulated sequence and its metadata, including D allele positions and trimming details.
        """
        # Extract Current D Metadata
        d_start, d_end = container.d_sequence_start, container.d_sequence_end
        d_germline_start, d_germline_end = container.d_germline_start, container.d_germline_end
        d_allele_remainder = container.sequence[d_start:d_end]
        d_allele_ref = self.d_dict[container.d_call[0]]

        # Get the junction inserted after trimming to the sequence
        junction_5 = container.sequence[container.v_sequence_end:d_start]
        junction_3 = container.sequence[d_end:container.j_sequence_start]

        # Get the trimming lengths
        d_trim_5 = container.d_trim_5
        d_trim_3 = container.d_trim_3

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

        container.d_sequence_start = d_start
        container.d_sequence_end = d_end
        container.d_germline_start = d_germline_start
        container.d_germline_end = d_germline_end
        container.d_trim_5 = d_trim_5
        container.d_trim_3 = d_trim_3
    def apply(self, container: SimulationContainer) -> None:
        # Implement the logic to correct D position after trimming index ambiguity
        self.fix_d_position_after_trimming_index_ambiguity(container)

    def get_graph_node(self):
        """Generates a detailed GraphViz node representation with constructor details in an HTML-like format."""
        step_name = "Fix D Position After Trimming Ambiguity"

        # Constructing an HTML-like label using a table for detailed formatting
        label = f"""
        <TABLE BORDER="0" CELLBORDER="1" CELLSPACING="0" CELLPADDING="4">
        <TR><TD COLSPAN="2" BGCOLOR="lightsteelblue"><B>{step_name}</B></TD></TR>
        <TR><TD ALIGN="LEFT"><B>Description</B></TD><TD ALIGN="LEFT">Corrects D segment start/end positions</TD></TR>
        </TABLE>
        
        """

        return label, 'box', "filled,rounded", CORRECTION_STEP_HEADING_COLOR, "Helvetica", "black"
