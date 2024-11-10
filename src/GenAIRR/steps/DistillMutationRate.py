from GenAIRR.container.SimulationContainer import SimulationContainer
from GenAIRR.steps.StepBase import AugmentationStep
from GenAIRR.utilities import DataConfig

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
        """Returns a string representation of the step for GraphViz with relevant information."""
        return (
            f'"DistillMutationRate" [label="DistillMutationRate\\n'
            f'Distills mutation rate\\nExcludes NP region mutations"]'
        )