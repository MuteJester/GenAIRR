import os
import pickle
import logging
from importlib.resources import files

_CONFIG_NAMES = {
    'HUMAN_IGH_OGRDB',
    'HUMAN_IGH_EXTENDED',
    'HUMAN_IGL_OGRDB',
    'HUMAN_IGK_OGRDB',
    'HUMAN_TCRB_IMGT',
}

_logger = logging.getLogger(__name__)
_CACHE = {}  # Cache to store already loaded data


def __getattr__(name: str):
    """
    This function lazy-loads the DataConfig object, which now
    already contains its own metadata.
    """
    if name in _CACHE:
        return _CACHE[name]

    if name not in _CONFIG_NAMES:
        raise AttributeError(f"Module '{__name__}' has no attribute '{name}'")

    filename = f"{name}.pkl"
    resource = files(__package__).joinpath('builtin_dataconfigs', filename)

    _logger.info("Loading '%s' data config...", name)

    try:
        with resource.open('rb') as f:
            # The loaded object is now complete and already has the .metadata attribute
            data_config = pickle.load(f)
    except FileNotFoundError:
        raise FileNotFoundError(f"Could not find data file for '{name}' (resource: {resource})")
    except Exception as e:
        raise ImportError(f"Could not load data config '{name}': {e}") from e

    # Cache the loaded data
    _CACHE[name] = data_config
    globals()[name] = data_config  # Make it a real module attribute

    return data_config


def __dir__() -> list[str]:
    """
    This function tells IDEs and `dir()` what names are available
    for import and autocompletion, hiding internal variables.
    """
    # Change this function to ONLY return the list of config names.
    return list(_CONFIG_NAMES)
