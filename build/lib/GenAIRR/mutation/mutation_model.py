from abc import ABC, abstractmethod


class MutationModel(ABC):
    """Abstract base class for mutation models.

   This class serves as a template for defining various mutation models. Subclasses are expected to implement the `apply_mutation` method, which applies the mutation logic to a given sequence.

   Methods:
       apply_mutation(sequence): Abstract method to apply mutation to a sequence. Must be implemented by subclasses.
   """
    @abstractmethod
    def apply_mutation(self, sequence):
        """Applies mutation to a given sequence according to the model.

            This is an abstract method and must be implemented by subclasses to define the specific mutation behavior.

            Args:
                sequence (sequence_object): The original sequence (DNA) to which mutations will be applied.

            Returns:
                str: The mutated sequence after applying the mutation model.

            Raises:
                NotImplementedError: If the subclass does not implement this method.
        """

        pass
