from ..container.SimulationContainer import SimulationContainer
from ..steps import SimulateSequence
from ..steps.StepBase import AugmentationStep
from ..utilities import DataConfig
from graphviz import Digraph
from ..container.SimulationContainer import SimulationContainer


class AugmentationPipeline:
    def __init__(self, steps: list):
        self.steps = steps

    def execute(self) -> SimulationContainer:
        if len(self.steps) == 0 or type(self.steps[0]) != SimulateSequence:
            raise Exception("First Step must be an instance of SimulateSequence")

        container = SimulationContainer()
        for step in self.steps:
            step.apply(container)
        return container

    def plot(self, filename='pipeline_diagram'):
        """
        Creates and saves a visually enhanced GraphViz diagram of the pipeline showing each step as a detailed block.

        Args:
            filename (str): The name of the file to save the diagram as (without extension).
        """
        dot = Digraph(comment='Augmentation Pipeline')


        # Define general style for all nodes
        dot.attr('node', shape='box', style='rounded, filled', color='lightblue', fontname='Helvetica')

        # Add nodes for each step with detailed labels and color coding
        for i, step in enumerate(self.steps):
            label,shape, style, color, fontname,fontcolor = step.get_graph_node()
            dot.node(
                f'step_{i}',
                label=f'<{label}>',
                shape=shape,
                style=style,
                color=color,
                fontname=fontname,
                fontcolor=fontcolor
            )

            # Connect each step to the next
            if i > 0:
                dot.edge(f'step_{i - 1}', f'step_{i}', color='black')

        # Render the diagram
        dot.render(filename, format='png', cleanup=True)
        print(f'Pipeline diagram saved as {filename}.png')