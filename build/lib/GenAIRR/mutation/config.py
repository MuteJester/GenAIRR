from enum import Enum, auto

class MutationModels(Enum):
    """Enumeration for specifying mutation models.

    This enumeration defines the types of mutation models that can be used in simulations or analyses.

    Attributes:
        UNIFORM: Represents a uniform mutation model where mutations are assumed to occur at a constant rate across all sites.
        S5F: Represents a specific, more complex 5-mer based mutation model.
    """
    UNIFORM = auto()
    S5F = auto()

class MutationConfig:
    """Configuration for mutation simulations or analyses.

    This class encapsulates the configuration parameters for mutations, including the mutation rate and the model used to simulate mutations.

    Attributes:
        mutation_rate (float): The rate at which mutations occur. This value is model-dependent and can represent different concepts depending on the model, such as the probability of mutation per site per .
        model (MutationModels): The mutation model to be used. Defaults to `MutationModels.UNIFORM`, indicating a uniform mutation rate across all sites.

    Args:
        mutation_rate (float): The rate of mutations. This could be, for example, the probability of mutation per site per generation.
        model (MutationModels, optional): The model of mutation to apply. Defaults to `MutationModels.UNIFORM`.
    """
    def __init__(self,mutation_rate,model=MutationModels.UNIFORM):
        self.mutation_rate = mutation_rate
        self.model = model