import warnings
from abc import ABC, abstractmethod
from typing import Optional

from .. import dataconfig
from ..container.SimulationContainer import SimulationContainer
from ..dataconfig import DataConfig
from ..dataconfig.enums import ChainType


class AugmentationStep(ABC):
    """
    Abstract base class for all augmentation pipeline steps.

    Configuration is provided by the Pipeline during execution via _bind_config().
    The class-level set_dataconfig() method is deprecated but still supported
    for backwards compatibility.
    """

    # Class-level attributes (deprecated, for backwards compatibility)
    _class_dataconfig: Optional[DataConfig] = None
    _class_chain_type: Optional[ChainType] = None

    def __init__(self):
        # Instance-level config (set by pipeline via _bind_config)
        self._bound_config: Optional[DataConfig] = None

    def _bind_config(self, config: DataConfig) -> None:
        """
        Internal method called by Pipeline to bind configuration to this step.

        Args:
            config: The DataConfig to use for this step's execution.
        """
        self._bound_config = config

    @property
    def dataconfig(self) -> DataConfig:
        """
        Returns the DataConfig for this step.

        Priority:
        1. Instance-level config (bound by pipeline)
        2. Class-level config (deprecated, set via set_dataconfig)

        Raises:
            ValueError: If no config is available.
        """
        if self._bound_config is not None:
            return self._bound_config

        if AugmentationStep._class_dataconfig is not None:
            return AugmentationStep._class_dataconfig

        raise ValueError(
            "No DataConfig available. Pass config to Pipeline: "
            "Pipeline(config=HUMAN_IGH_OGRDB, steps=[...])"
        )

    @property
    def chain_type(self) -> ChainType:
        """Returns the chain type from the current DataConfig."""
        return self.dataconfig.metadata.chain_type

    @classmethod
    def set_dataconfig(cls, config: DataConfig) -> None:
        """
        Set the class-level DataConfig (deprecated).

        This method is deprecated. Pass config directly to Pipeline instead:

            # New way (recommended):
            pipeline = Pipeline(config=HUMAN_IGH_OGRDB, steps=[...])

            # Old way (deprecated):
            AugmentationStep.set_dataconfig(HUMAN_IGH_OGRDB)
            pipeline = AugmentationPipeline([...])

        Args:
            config: The DataConfig to use for all steps.
        """
        warnings.warn(
            "AugmentationStep.set_dataconfig() is deprecated. "
            "Pass config directly to Pipeline: Pipeline(config=your_config, steps=[...])",
            DeprecationWarning,
            stacklevel=2
        )
        cls._class_dataconfig = config
        cls._class_chain_type = config.metadata.chain_type

    @abstractmethod
    def apply(self, container: SimulationContainer) -> None:
        """
        Apply the augmentation step to the provided simulated sequence data.

        Args:
            container (SimulationContainer): The instance containing the simulated sequence and its metadata.
        """
        pass
