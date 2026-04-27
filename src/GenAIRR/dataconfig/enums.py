from enum import Enum,auto

from enum import Enum


class Productivity(Enum):
    """Filter mode for which V(D)J rearrangements appear in the output.

    PRODUCTIVE_ONLY:
        The simulator retries the rearrangement phase (up to
        ``max_productive_attempts``, default 25) until ``rec.productive``
        is true. Every returned record has ``productive == True`` after
        recombination. **Caveat (T2-16):** SHM is applied AFTER the retry
        boundary, so heavy mutation rates can flip individual records to
        ``productive=False`` in the final AIRR output.

    NON_PRODUCTIVE_ONLY:
        Symmetric to PRODUCTIVE_ONLY — retries until ``rec.productive``
        is false. Useful for studying out-of-frame statistics, training
        an aligner on negative examples, or tuning productive filters.

    PRODUCTIVE_MIXED (default):
        No filtering — emits whatever recombination produces. Pure V(D)J
        biology yields ~22.8% productive on human IGH (T1-15), matching
        theoretical 33% in-frame × ~68% stop-free rates.
    """
    PRODUCTIVE_ONLY     = "productive_only"
    NON_PRODUCTIVE_ONLY = "non_productive_only"
    PRODUCTIVE_MIXED    = "mixed"


class ChainType(Enum):
    BCR_HEAVY = "BCR_HEAVY"
    BCR_LIGHT_KAPPA = "BCR_LIGHT_KAPPA"
    BCR_LIGHT_LAMBDA = "BCR_LIGHT_LAMBDA"
    TCR_ALPHA = "TCR_ALPHA"
    TCR_BETA = "TCR_BETA"
    TCR_GAMMA = "TCR_GAMMA"
    TCR_DELTA = "TCR_DELTA"

    @property
    def has_d(self) -> bool:
        """Returns True if the chain type utilizes a D segment."""
        # A set provides a fast lookup of members that have a D segment.
        d_segment_chains = {
            ChainType.BCR_HEAVY,
            ChainType.TCR_BETA,
            ChainType.TCR_DELTA
        }
        return self in d_segment_chains

class Species(Enum):
    # Common Mammalian Models
    HUMAN = "Human"
    MOUSE = "Mouse"
    MOUSE_C57BL6J = "Mouse C57BL/6J"
    RAT = "Rat"
    RABBIT = "Rabbit"
    GUINEA_PIG = "Guinea Pig"

    # Primates
    RHESUS_MACAQUE = "Rhesus Macaque"
    CYNOMOLGUS_MACAQUE = "Cynomolgus Macaque"
    MARMOSET = "Marmoset"
    CHIMPANZEE = "Chimpanzee"
    GORILLA = "Gorilla"
    SUMATRAN_ORANGUTAN = "Sumatran Orangutan"
    BORNEAN_ORANGUTAN = "Bornean Orangutan"
    OWL_MONKEY = "Owl Monkey"
    RING_TAILED_LEMUR = "Ring-tailed Lemur"

    # Agricultural and Domestic Animals
    PIG = "Pig"
    COW = "Cow"
    SHEEP = "Sheep"
    GOAT = "Goat"
    HORSE = "Horse"
    DOG = "Dog"
    CAT = "Cat"

    # Camelids
    LLAMA = "Llama"
    ALPACA = "Alpaca"
    DROMEDARY_CAMEL = "Dromedary Camel"

    # Mustelids
    FERRET = "Ferret"
    AMERICAN_MINK = "American Mink"

    # Other Mammals
    NAKED_MOLE_RAT = "Naked Mole-rat"
    PLATYPUS = "Platypus"

    # Birds
    CHICKEN = "Chicken"
    DUCK = "Duck"
    TURKEY = "Turkey"

    # Fish and Aquatic Vertebrates
    ZEBRAFISH = "Zebrafish"
    TROUT = "Trout"
    SALMON = "Salmon"
    CATFISH = "Catfish"
    ATLANTIC_COD = "Atlantic Cod"
    SHARK = "Shark"

    # Other
    BAT = "Bat"

