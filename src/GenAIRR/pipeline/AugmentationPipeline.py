import warnings
from typing import Optional, List

from ..container.SimulationContainer import SimulationContainer
from ..steps import SimulateSequence
from ..steps.StepBase import AugmentationStep
from ..dataconfig import DataConfig


class AugmentationPipeline:
    """
    Pipeline for executing a sequence of augmentation steps.

    Args:
        steps: List of AugmentationStep instances to execute in order.
        config: DataConfig containing reference data and parameters.
                If not provided, falls back to class-level config (deprecated).

    Example:
        >>> from GenAIRR import Pipeline, steps, HUMAN_IGH_OGRDB, S5F
        >>> pipeline = Pipeline(
        ...     config=HUMAN_IGH_OGRDB,
        ...     steps=[
        ...         steps.SimulateSequence(S5F(0.003, 0.25), productive=True),
        ...         steps.FixVPositionAfterTrimmingIndexAmbiguity(),
        ...     ]
        ... )
        >>> result = pipeline.execute()
    """

    def __init__(self, steps: List[AugmentationStep], config: Optional[DataConfig] = None):
        self.steps = steps
        self._config = config

    @property
    def config(self) -> DataConfig:
        """Returns the pipeline's DataConfig, with fallback to class-level config."""
        if self._config is not None:
            return self._config

        # Fallback to class-level config for backwards compatibility
        if AugmentationStep._class_dataconfig is not None:
            warnings.warn(
                "Using class-level AugmentationStep.set_dataconfig() is deprecated. "
                "Pass config directly to Pipeline: Pipeline(config=your_config, steps=[...])",
                DeprecationWarning,
                stacklevel=3
            )
            return AugmentationStep._class_dataconfig

        raise ValueError(
            "No DataConfig provided. Pass config to Pipeline: "
            "Pipeline(config=HUMAN_IGH_OGRDB, steps=[...])"
        )

    def execute(self) -> SimulationContainer:
        """
        Execute all pipeline steps and return the simulation result.

        Returns:
            SimulationContainer with the simulated sequence and metadata.

        Raises:
            Exception: If the first step is not a SimulateSequence instance.
            ValueError: If no DataConfig is available.
        """
        if len(self.steps) == 0 or type(self.steps[0]) != SimulateSequence:
            raise Exception("First Step must be an instance of SimulateSequence")

        # Get config (either from pipeline or class-level with deprecation warning)
        config = self.config

        # Bind config to all steps before execution
        for step in self.steps:
            step._bind_config(config)

        container = SimulationContainer()
        for step in self.steps:
            step.apply(container)
        return container

    def __getitem__(self, index) -> AugmentationStep:
        """
        Returns the step at the specified index.

        Args:
            index (int): The index of the step to retrieve.

        Returns:
            AugmentationStep: The step at the specified index.
        """
        return self.steps[index]

    def __setitem__(self, key, value):
        self.steps[key] = value

    def plot(self, filename='pipeline_diagram'):
        """
        Creates and saves a visually enhanced GraphViz diagram of the pipeline showing each step as a detailed block.

        Args:
            filename (str): The name of the file to save the diagram as (without extension).
        """
        try:
            from graphviz import Digraph
        except ImportError:
            raise ImportError(
                "graphviz is required for pipeline visualization. "
                "Install it with: pip install graphviz"
            )
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
