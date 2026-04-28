import os
import pickle
import logging
from importlib.resources import files

# All available DataConfig names — each corresponds to a .pkl file
# in the builtin_dataconfigs directory. All are lazy-loaded on first access.
_CONFIG_NAMES = {
    # ── Original configs (OGRDB / custom) ──
    'HUMAN_IGH_OGRDB',
    'HUMAN_IGH_EXTENDED',
    'HUMAN_IGL_OGRDB',
    'HUMAN_IGK_OGRDB',

    # ── IMGT V-QUEST reference configs ──
    # Human
    'HUMAN_IGH_IMGT',
    'HUMAN_IGK_IMGT',
    'HUMAN_IGL_IMGT',
    'HUMAN_TCRA_IMGT',
    'HUMAN_TCRB_IMGT',
    'HUMAN_TCRD_IMGT',
    'HUMAN_TCRG_IMGT',

    # Mouse
    'MOUSE_IGH_IMGT',
    'MOUSE_IGK_IMGT',
    'MOUSE_IGL_IMGT',
    'MOUSE_TCRA_IMGT',
    'MOUSE_TCRB_IMGT',
    'MOUSE_TCRD_IMGT',
    'MOUSE_TCRG_IMGT',

    # Mouse C57BL/6J
    'MOUSE_C57BL6J_IGH_IMGT',
    'MOUSE_C57BL6J_IGK_IMGT',
    'MOUSE_C57BL6J_IGL_IMGT',
    'MOUSE_C57BL6J_TCRA_IMGT',
    'MOUSE_C57BL6J_TCRB_IMGT',
    'MOUSE_C57BL6J_TCRD_IMGT',
    'MOUSE_C57BL6J_TCRG_IMGT',

    # Rat
    'RAT_IGH_IMGT',
    'RAT_IGK_IMGT',
    'RAT_IGL_IMGT',

    # Rabbit
    'RABBIT_IGH_IMGT',
    'RABBIT_IGK_IMGT',
    'RABBIT_IGL_IMGT',
    'RABBIT_TCRA_IMGT',
    'RABBIT_TCRB_IMGT',
    'RABBIT_TCRD_IMGT',
    'RABBIT_TCRG_IMGT',

    # Rhesus macaque
    'RHESUS_IGH_IMGT',
    'RHESUS_IGK_IMGT',
    'RHESUS_IGL_IMGT',
    'RHESUS_TCRA_IMGT',
    'RHESUS_TCRB_IMGT',
    'RHESUS_TCRD_IMGT',
    'RHESUS_TCRG_IMGT',

    # Cynomolgus macaque
    'CYNOMOLGUS_IGH_IMGT',
    'CYNOMOLGUS_TCRB_IMGT',

    # Gorilla
    'GORILLA_IGH_IMGT',
    'GORILLA_IGK_IMGT',
    'GORILLA_IGL_IMGT',
    'GORILLA_TCRA_IMGT',
    'GORILLA_TCRB_IMGT',
    'GORILLA_TCRD_IMGT',
    'GORILLA_TCRG_IMGT',

    # Cow
    'COW_IGH_IMGT',
    'COW_IGK_IMGT',
    'COW_IGL_IMGT',
    'COW_TCRA_IMGT',
    'COW_TCRB_IMGT',
    'COW_TCRD_IMGT',
    'COW_TCRG_IMGT',

    # Sheep
    'SHEEP_IGH_IMGT',
    'SHEEP_IGK_IMGT',
    'SHEEP_IGL_IMGT',
    'SHEEP_TCRA_IMGT',
    'SHEEP_TCRB_IMGT',
    'SHEEP_TCRD_IMGT',

    # Goat
    'GOAT_IGK_IMGT',
    'GOAT_IGL_IMGT',

    # Pig
    'PIG_IGH_IMGT',
    'PIG_IGK_IMGT',
    'PIG_IGL_IMGT',
    'PIG_TCRB_IMGT',
    'PIG_TCRG_IMGT',

    # Horse
    'HORSE_IGH_IMGT',
    'HORSE_IGK_IMGT',

    # Dog
    'DOG_IGH_IMGT',
    'DOG_IGK_IMGT',
    'DOG_IGL_IMGT',
    'DOG_TCRA_IMGT',
    'DOG_TCRB_IMGT',
    'DOG_TCRD_IMGT',
    'DOG_TCRG_IMGT',

    # Cat
    'CAT_IGK_IMGT',
    'CAT_IGL_IMGT',
    'CAT_TCRA_IMGT',
    'CAT_TCRB_IMGT',
    'CAT_TCRD_IMGT',
    'CAT_TCRG_IMGT',

    # Camelids
    'DROMEDARY_IGK_IMGT',
    'DROMEDARY_TCRB_IMGT',
    'DROMEDARY_TCRG_IMGT',
    'ALPACA_IGH_IMGT',

    # Mustelids
    'FERRET_IGH_IMGT',
    'FERRET_IGK_IMGT',
    'FERRET_IGL_IMGT',
    'FERRET_TCRA_IMGT',
    'FERRET_TCRB_IMGT',
    'FERRET_TCRD_IMGT',
    'FERRET_TCRG_IMGT',

    # Platypus
    'PLATYPUS_IGH_IMGT',

    # Chicken
    'CHICKEN_IGH_IMGT',
    'CHICKEN_IGL_IMGT',

    # Fish
    'ZEBRAFISH_IGH_IMGT',
    'ZEBRAFISH_TCRA_IMGT',
    'ZEBRAFISH_TCRD_IMGT',
    'TROUT_IGH_IMGT',
    'TROUT_TCRB_IMGT',
    'SALMON_IGH_IMGT',
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

    # T2-1: refuse to ship a stale-schema or tampered builtin. A mismatch
    # here means the wheel itself is broken (someone shipped without
    # running the migration) or the file was modified post-install.
    try:
        data_config.verify_integrity()
    except Exception as e:
        raise ImportError(
            f"Builtin DataConfig '{name}' failed integrity check: {e}"
        ) from e

    # Cache the loaded data
    _CACHE[name] = data_config
    globals()[name] = data_config  # Make it a real module attribute

    return data_config


def list_configs():
    """Return a sorted list of all available DataConfig names.

    Example::

        from GenAIRR.data import list_configs
        print(list_configs())
        # ['ALPACA_IGH_IMGT', 'CAT_IGK_IMGT', ...]
    """
    return sorted(_CONFIG_NAMES)


def __dir__() -> list[str]:
    """
    This function tells IDEs and ``dir()`` what names are available
    for import and autocompletion.
    """
    return sorted(_CONFIG_NAMES | {'list_configs'})
