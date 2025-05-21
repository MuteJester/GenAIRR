from GenAIRR.container.SimulationContainer import SimulationContainer
from GenAIRR.pipeline.plot_parameters import CORRECTION_STEP_HEADING_COLOR
from GenAIRR.steps.StepBase import AugmentationStep
from GenAIRR.dataconfig import DataConfig

class DistillMutationRate(AugmentationStep):

    def __init__(self):
        super().__init__()

    def distill_mutation_rate(self,container: SimulationContainer):
        # remove mutations in NP region
        np_positions = list(range(container.v_sequence_end,container.d_sequence_start+1))
        np_positions += list(range(container.d_sequence_end,container.j_sequence_start+1))
        mutations_in_np_regions = list(set(np_positions)&set(list(container.mutations.keys())))
        for pos in mutations_in_np_regions:
            container.mutations.pop(pos)

        distilled_mutation_rate = (len(container.mutations)+len(container.Ns))/len(container.sequence)
        container.mutation_rate = distilled_mutation_rate
    def apply(self, container: SimulationContainer) -> None:
        self.distill_mutation_rate(container)

    def get_graph_node(self):
        """Generates a detailed GraphViz node representation with constructor details in an HTML-like format."""
        step_name = "Distill Mutation Rate"

        # Constructing an HTML-like label using a table for detailed formatting
        label = f"""
        <TABLE BORDER="0" CELLBORDER="1" CELLSPACING="0" CELLPADDING="4">
        <TR><TD COLSPAN="2" BGCOLOR="lightsteelblue"><B>{step_name}</B></TD></TR>
        <TR><TD ALIGN="LEFT"><B>Description</B></TD><TD ALIGN="LEFT">Distills mutation rate by excluding NP regions</TD></TR>
        </TABLE>
        
        """

        return label, 'box', "filled,rounded", CORRECTION_STEP_HEADING_COLOR, "Helvetica", "black"
