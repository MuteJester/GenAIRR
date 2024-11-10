from ..container.SimulationContainer import SimulationContainer
from ..steps.StepBase import AugmentationStep
from ..utilities import DataConfig


class CorrectForVEndCut(AugmentationStep):
    def __init__(self):
        super().__init__()
        self.v_end_allele_correction_map = self.dataconfig.correction_maps['V_3_TRIM_SIMILARITY_MAP']
        self.max_v_end_correction_map_value = max(
            self.v_end_allele_correction_map[list(self.v_end_allele_correction_map)[0]])

    def correct_for_v_end_cut(self, container):
        """
        Corrects the V allele choices based on the trimming of the sequence's end, ensuring the chosen alleles remain consistent with the modified sequence.

                Args:
                    simulated (dict): A dictionary containing the simulated sequence and metadata, including the current V allele and its trimmed length.
        """
        sampled_v = container.v_call[0]
        valid_length = min(container.v_trim_3, self.max_v_end_correction_map_value)
        equivalent_alleles = self.v_end_allele_correction_map[sampled_v][valid_length]
        container.v_call = container.v_call + list(set(equivalent_alleles) - set(container.v_call))

    def apply(self, container: SimulationContainer) -> None:
        self.correct_for_v_end_cut(container)

    def get_graph_node(self):
        """Returns a string representation of the step for GraphViz with relevant information."""
        return f'"CorrectForVEndCut" [label="CorrectForVEndCut\\nMax V end correction: {self.max_v_end_correction_map_value}"]'