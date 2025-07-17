import os
import pickle
from dataclasses import dataclass
from datetime import date
from .enums import *

@dataclass(frozen=True)
class ConfigInfo:
    """A structured container for data configuration metadata."""
    species: Species
    chain_type: ChainType
    reference_set: str
    last_updated: date
    has_d: bool

    def __repr__(self) -> str:
        """Provides a clean, readable representation of the config info."""
        # Added the 'has_d' flag for a more complete summary
        return (
            f"Species: {self.species.value}, "
            f"Chain: {self.chain_type.value}, "
            f"Has D: {self.has_d}, "
            f"Reference: '{self.reference_set}', "
            f"Updated: {self.last_updated.strftime('%Y-%m-%d')}"
        )
