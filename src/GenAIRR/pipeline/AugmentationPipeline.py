from ..container.SimulationContainer import SimulationContainer
from ..utilities import DataConfig
from graphviz import Digraph


class AugmentationPipeline:
    def __init__(self, steps: list,dataconfig: DataConfig):
        self.steps = steps
        self.dataconfig = dataconfig
    def execute(self, container: SimulationContainer) -> None:
        for step in self.steps:
            step.apply(container,self.dataconfig)

    def plot(self, filename='pipeline_diagram'):
        dot = Digraph(comment='Augmentation Pipeline')

        # Add nodes for each step
        for i, step in enumerate(self.steps):
            step_node = step.get_graph_node()
            dot.node(f'step_{i}', step_node)

            # Connect each step to the next
            if i > 0:
                dot.edge(f'step_{i - 1}', f'step_{i}')

        # Render the diagram
        dot.render(filename, format='png', cleanup=True)
        print(f'Pipeline diagram saved as {filename}.png')
