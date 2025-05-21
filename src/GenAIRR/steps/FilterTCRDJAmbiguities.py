from ..container.SimulationContainer import SimulationContainer
from ..pipeline.plot_parameters import CORRECTION_STEP_HEADING_COLOR
from ..steps.StepBase import AugmentationStep
from ..dataconfig import DataConfig
import re

class FilterTCRDJAmbiguities(AugmentationStep):

    def __init__(self):
        super().__init__()

    def correct_for_dj_ambiguities(self, container):
        """
         During the simulation of a sequence, a D allele is selected, and the subsequent J allele selection ensures that
          the family number of the J gene is equal to or greater than that of the D gene. However, after introducing
          noise, such as indels, and applying corrections for D trims, ambiguities may arise where the D allele becomes
          indistinguishable from other D alleles in different families. This can lead to logical inconsistencies.
           By leveraging the family number of the J allele, these ambiguities can be resolved by filtering out
           conflicting entries from the D allele list.
        """
        d_calls = container.d_call
        j_calls = container.j_call


        valid_d_families = []
        #extracting family numbers from J calls
        for j_call in j_calls:
            j_gene_num = int(re.search(r'\d+', j_call).group())
            valid_d_families.append(j_gene_num)

        d_calls = [d_calls[0]] + [call for call in d_calls[1:] if int(re.search(r'\d+', call).group()) in valid_d_families]

        container.d_call = d_calls



    def apply(self, container: SimulationContainer) -> None:
        self.correct_for_dj_ambiguities(container)

    def get_graph_node(self):
        """Generates a detailed GraphViz node representation with constructor details in an HTML-like format."""
        step_name = "Filter TCR DJ Ambiguities"

        # Constructing an HTML-like label using a table for detailed formatting
        label = f"""
        <TABLE BORDER="0" CELLBORDER="1" CELLSPACING="0" CELLPADDING="4">
        <TR><TD COLSPAN="2" BGCOLOR="lightsteelblue"><B>{step_name}</B></TD></TR>
        <TR><TD ALIGN="LEFT"><B>D Trim Correction Map</B></TD><TD ALIGN="LEFT">Loaded</TD></TR>
        </TABLE>

        """

        return label, 'box', "filled,rounded", CORRECTION_STEP_HEADING_COLOR, "Helvetica", "black"
