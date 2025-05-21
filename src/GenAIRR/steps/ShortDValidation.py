from GenAIRR.container.SimulationContainer import SimulationContainer
from GenAIRR.pipeline.plot_parameters import VALIDATION_STEP_BOX_COLOR
from GenAIRR.steps.StepBase import AugmentationStep
from GenAIRR.dataconfig import DataConfig


class ShortDValidation(AugmentationStep):
    def __init__(self,short_d_length = 5):
        super().__init__()
        self.short_d_length = short_d_length

    def short_d_validation(self, container: SimulationContainer):
        """
                Validates the length of the D gene segment in the simulated sequence. If the segment is shorter than a predefined threshold, it is labeled as "Short-D".

                Args:
                    simulated (SimulationContainer):  containing the simulated sequence and its metadata.
        """
        d_length = container.d_sequence_end - container.d_sequence_start
        if d_length < self.short_d_length:
            container.d_call = ['Short-D']
    def apply(self, container: SimulationContainer) -> None:
        self.short_d_validation(container)

    def get_graph_node(self):
        """Generates a detailed GraphViz node representation with constructor details in an HTML-like format."""
        step_name = "Short D Validation"

        # Constructing an HTML-like label using a table for detailed formatting
        label = f"""
        <TABLE BORDER="0" CELLBORDER="1" CELLSPACING="0" CELLPADDING="4">
        <TR><TD COLSPAN="2" BGCOLOR="lightsteelblue"><B>{step_name}</B></TD></TR>
        <TR><TD ALIGN="LEFT"><B>Short D Length Threshold</B></TD><TD ALIGN="LEFT">{self.short_d_length}</TD></TR>
        </TABLE>
        
        """

        return label, 'box', "filled,rounded", VALIDATION_STEP_BOX_COLOR, "Helvetica", "black"
