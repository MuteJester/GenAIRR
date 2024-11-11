from abc import ABC, abstractmethod
from ..container.SimulationContainer import SimulationContainer
from ..utilities import DataConfig


class AugmentationStep(ABC):
    # Class-level attributes
    dataconfig = None
    chain_type = None
    @classmethod
    def set_dataconfig(cls, config: DataConfig,chain_type: int):
        cls.dataconfig = config
        cls.chain_type = chain_type
    @abstractmethod
    def apply(self, container: SimulationContainer) -> None:
        """
        Apply the augmentation step to the provided simulated sequence data.

        Args:
            container (SimulationContainer): The instance containing the simulated sequence and its metadata.
        """
        pass
