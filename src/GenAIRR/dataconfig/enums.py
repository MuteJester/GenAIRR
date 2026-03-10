from enum import Enum,auto

from enum import Enum

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

