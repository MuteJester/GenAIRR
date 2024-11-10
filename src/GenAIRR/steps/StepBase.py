from abc import ABC, abstractmethod
from ..container.SimulationContainer import SimulationContainer
from ..utilities import DataConfig


class AugmentationStep(ABC):
    dataconfig = None  # Class-level attribute to store the shared dataconfig

    @classmethod
    def set_dataconfig(cls, config):
        cls.dataconfig = config
    @abstractmethod
    def apply(self, container: SimulationContainer) -> None:
        """
        Apply the augmentation step to the provided simulated sequence data.

        Args:
            simulated (dict): The dictionary containing the simulated sequence and its metadata.
        """
        pass
