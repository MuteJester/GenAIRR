from ..container.SimulationContainer import SimulationContainer
from ..steps.StepBase import AugmentationStep
from ..utilities import DataConfig

class CorrectForDTrims(AugmentationStep):

    def __init__(self):
        super().__init__()
        self.d_trim_correction_map = self.dataconfig.correction_maps['D_5_3_TRIM_SIMILARITY_MAP']

    def correct_for_d_trims(self, container):
        """
                Adjusts the D allele choices in the simulated sequence based on 5' and 3' trimming information.

                Args:
                    simulated (dict): Dictionary containing the simulated sequence and its metadata, including D allele trims.
        """
        # Get the 5' and 3' trims of the d allele in the simulated sequence
        trim_5 = container.d_trim_5
        trim_3 = container.d_trim_3
        # infer the precalculated map what alleles should be the ground truth for this sequence based on the trim
        #simulated['d_call'] = list(self.d_trim_correction_map[simulated['d_call'][0]][(trim_5, trim_3)])
        # retaining the order such that the selected d allele is the first.
        sampled_d = container.d_call[0]
        container.d_call = [sampled_d] + list(set(self.d_trim_correction_map[sampled_d][(trim_5, trim_3)]) - {
            sampled_d})

    def apply(self, container: SimulationContainer) -> None:
        self.correct_for_d_trims(container)

    def get_graph_node(self):
        """Returns a string representation of the step for GraphViz with relevant information."""
        return f'"CorrectForDTrims" [label="CorrectForDTrims\\nD trim correction map used"]'