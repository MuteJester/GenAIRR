from GenAIRR.container.SimulationContainer import SimulationContainer
from GenAIRR.steps.StepBase import AugmentationStep
from GenAIRR.utilities import DataConfig

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
        """Returns a string representation of the step for GraphViz with relevant information."""
        return (
            f'"ShortDValidation" [label="ShortDValidation\\n'
            f'Short D Length Threshold: {self.short_d_length}"]'
        )