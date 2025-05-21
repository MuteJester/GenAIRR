from ..container.SimulationContainer import SimulationContainer
from ..pipeline.plot_parameters import CORRECTION_STEP_HEADING_COLOR
from ..steps.StepBase import AugmentationStep
from ..dataconfig import DataConfig

class FixJPositionAfterTrimmingIndexAmbiguity(AugmentationStep):
    def __init__(self):
        super().__init__()
        self.j_alleles = sorted([i for j in self.dataconfig.j_alleles for i in self.dataconfig.j_alleles[j]],
                                key=lambda x: x.name)
        self.j_dict = {i.name: i.ungapped_seq.upper() for i in self.j_alleles}
        self.j_anchor = {i.name: i.anchor for i in self.j_alleles}
        self.j_frame = {i.name: i.frame for i in self.j_alleles}

    def fix_j_position_after_trimming_index_ambiguity(self, container):
        """
            Addresses ambiguities in J allele positions due to sequence trimming, ensuring J allele boundaries are accurately represented.

                Args:
                    simulation (dict): A dictionary containing the simulated sequence and metadata, including J allele positions and trimming details.
        """
        # Extract Current J Metadata
        j_start, j_end = container.j_sequence_start, container.j_sequence_end
        j_germline_start, j_germline_end = container.j_germline_start, container.j_germline_end
        j_allele_remainder = container.sequence[j_start:j_end]
        j_allele_ref = self.j_dict[container.j_call[0]]

        # Get the junction inserted after trimming to the sequence
        junction_5 = container['sequence'][container['d_sequence_end']:j_start]

        # Get the trimming lengths
        j_trim_5 = container.j_trim_5

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

        container.j_sequence_start = j_start
        container.j_germline_start = j_germline_start
        container.j_trim_5 = j_trim_5
    def apply(self, container: SimulationContainer) -> None:
        self.fix_j_position_after_trimming_index_ambiguity(container)

    def get_graph_node(self):
        """Generates a detailed GraphViz node representation with constructor details in an HTML-like format."""
        step_name = "Fix J Position After Trimming Ambiguity"

        # Constructing an HTML-like label using a table for detailed formatting
        label = f"""
        <TABLE BORDER="0" CELLBORDER="1" CELLSPACING="0" CELLPADDING="4">
        <TR><TD COLSPAN="2" BGCOLOR="lightsteelblue"><B>{step_name}</B></TD></TR>
        <TR><TD ALIGN="LEFT"><B>Description</B></TD><TD ALIGN="LEFT">Corrects J segment start positions</TD></TR>
        </TABLE>
        
        """

        return label, 'box', "filled,rounded", CORRECTION_STEP_HEADING_COLOR, "Helvetica", "black"
