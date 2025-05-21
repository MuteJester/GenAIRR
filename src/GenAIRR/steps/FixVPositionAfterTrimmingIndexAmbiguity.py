from GenAIRR.container.SimulationContainer import SimulationContainer
from GenAIRR.pipeline.plot_parameters import CORRECTION_STEP_HEADING_COLOR
from GenAIRR.steps.StepBase import AugmentationStep
from GenAIRR.dataconfig import DataConfig



class FixVPositionAfterTrimmingIndexAmbiguity(AugmentationStep):
    def __init__(self):
        super().__init__()
        self.v_alleles = sorted([i for j in self.dataconfig.v_alleles for i in self.dataconfig.v_alleles[j]],
                                key=lambda x: x.name)

        self.j_alleles = sorted([i for j in self.dataconfig.j_alleles for i in self.dataconfig.j_alleles[j]],
                                key=lambda x: x.name)
        self.v_dict = {i.name: i.ungapped_seq.upper() for i in self.v_alleles}

    def fix_v_position_after_trimming_index_ambiguity(self, container):
        """
                Resolves any ambiguities in V allele positions resulting from sequence trimming, ensuring accurate representation of V allele boundaries.

                Args:
                    container (SimulationContainer): An object containing the simulated sequence and metadata, including V allele positions and trimming details.
        """
        # Extract Current V Metadata
        v_start, v_end = container.v_sequence_start, container.v_sequence_end
        v_germline_start, v_germline_end = container.v_germline_start, container.v_germline_end
        v_allele_remainder = container.sequence[v_start:v_end]
        v_allele_ref = self.v_dict[container.v_call[0]]

        # Get the junction inserted after trimming to the sequence
        junction_3 = container.sequence[v_end:container.d_sequence_start]

        # Get the trimming lengths
        v_trim_3 = container.v_trim_3

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

        container.v_sequence_end = v_end
        container.v_germline_end = v_germline_end
        container.v_trim_3 = v_trim_3

    def apply(self, container: SimulationContainer) -> None:
        self.fix_v_position_after_trimming_index_ambiguity(container)

    def get_graph_node(self):
        """Generates a detailed GraphViz node representation with constructor details in an HTML-like format."""
        step_name = "Fix V Position After Trimming Ambiguity"

        # Constructing an HTML-like label using a table for detailed formatting
        label = f"""
        <TABLE BORDER="0" CELLBORDER="1" CELLSPACING="0" CELLPADDING="4">
        <TR><TD COLSPAN="2" BGCOLOR="lightsteelblue"><B>{step_name}</B></TD></TR>
        <TR><TD ALIGN="LEFT"><B>Description</B></TD><TD ALIGN="LEFT">Corrects V segment end positions</TD></TR>
        </TABLE>
        
        """

        return label, 'box', "filled,rounded", CORRECTION_STEP_HEADING_COLOR, "Helvetica", "black"
