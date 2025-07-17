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
    RAT = "Rat"
    RABBIT = "Rabbit"
    GUINEA_PIG = "Guinea Pig"

    # Primates
    RHESUS_MACAQUE = "Rhesus Macaque"
    CYNOMOLGUS_MACAQUE = "Cynomolgus Macaque"
    MARMOSET = "Marmoset"

    # Agricultural and Domestic Animals
    PIG = "Pig"
    COW = "Cow"  # or Bovine
    SHEEP = "Sheep"
    GOAT = "Goat"
    HORSE = "Horse"
    DOG = "Dog"
    CAT = "Cat"

    # Camelids (for Heavy-Chain-Only Antibodies)
    LLAMA = "Llama"
    ALPACA = "Alpaca"
    DROMEDARY_CAMEL = "Dromedary Camel"

    # Birds
    CHICKEN = "Chicken"
    DUCK = "Duck"
    TURKEY = "Turkey"

    # Fish and Aquatic Vertebrates
    ZEBRAFISH = "Zebrafish"
    TROUT = "Trout"  # e.g., Rainbow Trout
    SALMON = "Salmon"  # e.g., Atlantic Salmon
    CATFISH = "Catfish"  # e.g., Channel Catfish
    SHARK = "Shark"  # For IgNAR Receptors

    # Other
    FERRET = "Ferret"  # Important model for influenza
    BAT = "Bat"  # Increasing interest in viral immunology

