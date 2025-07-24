from ..container.SimulationContainer import SimulationContainer
from ..pipeline.plot_parameters import SIMULATION_STEP_BOX_COLOR
from ..steps.StepBase import AugmentationStep
from ..dataconfig import DataConfig
from ..sequence import HeavyChainSequence
from ..dataconfig.enums import ChainType


class SimulateSequence(AugmentationStep):
    MAX_GENERATION_ATTEMPTS = 25
    def __init__(self, mutation_model, productive=False, specific_v=None, specific_d=None, specific_j=None):
        """
        Initializes the step for simulating a heavy chain sequence.

        Args:
            mutation_model: The mutation model to apply during sequence simulation.
            productive (bool): If True, ensures that the simulated sequence is productive.
        """
        super().__init__()
        self.mutation_model = mutation_model
        self.productive = productive
        self.specific_v = specific_v
        self.specific_d = specific_d
        self.specific_j = specific_j
        self.sequence_constructor_instance = self._get_sequence_class()

    def _get_sequence_class(self):
        """
        Returns the appropriate sequence class based on the chain type.

        Args:
            chain_type (int): The type of chain to be simulated.

        Returns:
            class: The sequence class corresponding to the chain type.
        """
        if self.chain_type == ChainType.BCR_HEAVY:
            from ..sequence import HeavyChainSequence
            return HeavyChainSequence
        elif self.chain_type == ChainType.BCR_LIGHT_LAMBDA or self.chain_type == ChainType.BCR_LIGHT_KAPPA:
            from ..sequence import LightChainSequence
            return LightChainSequence
        elif self.chain_type == ChainType.TCR_BETA:
            from ..TCR.sequence.heavy_chain import TCRHeavyChainSequence
            return TCRHeavyChainSequence
        else:
            raise ValueError(f"Unknown chain type: {self.chain_type}")

    def simulate_sequence(self, container: SimulationContainer):
        """
        Simulates a heavy chain sequence and populates the SimulationContainer with relevant data.

        Args:
            container (SimulationContainer): The container to store the simulated sequence data.
            """

        gen_args = {
            'specific_v': self.specific_v,
            'specific_j': self.specific_j
        }

        # Only add specific_d if the chain type has a D segment
        if self.chain_type.has_d:
            gen_args['specific_d'] = self.specific_d

        # Create the sequence using the constructor instance
        gen = self.sequence_constructor_instance.create_random(self.dataconfig, **gen_args)

        # Ensure productivity if specified
        if self.productive and not gen.functional:
            for _ in range(self.MAX_GENERATION_ATTEMPTS):
                gen = self.sequence_constructor_instance.create_random(self.dataconfig, **gen_args)
                if gen.functional:
                    break  #exit the loop
            else:
                # This block now prints a warning instead of raising an error
                print(
                    f"Warning: Failed to generate a productive sequence after {self.MAX_GENERATION_ATTEMPTS} attempts. "
                    f"Proceeding with the last non-functional sequence. "
                    f"Parameters: V={self.specific_v}, J={self.specific_j}" + "D={}".format(self.specific_d if self.chain_type.has_d else "N/A")
                )


        container.from_instance(gen)
        # Apply mutation -
        gen.mutate(self.mutation_model)

        # Populate the container with generated sequence data
        container.sequence = gen.mutated_seq
        container.v_sequence_start = gen.v_seq_start
        container.v_sequence_end = gen.v_seq_end
        container.d_sequence_start = gen.d_seq_start
        container.d_sequence_end = gen.d_seq_end
        container.j_sequence_start = gen.j_seq_start
        container.j_sequence_end = gen.j_seq_end
        container.v_germline_start = gen.v_germline_start
        container.v_germline_end = gen.v_germline_end
        container.d_germline_start = gen.d_germline_start if hasattr(gen,'d_germline_start') else None
        container.d_germline_end = gen.d_germline_end if hasattr(gen,'d_germline_end') else None
        container.j_germline_start = gen.j_germline_start
        container.j_germline_end = gen.j_germline_end
        container.junction_sequence_start = gen.junction_start
        container.junction_sequence_end = gen.junction_end
        container.v_call = [gen.v_allele.name]
        container.d_call = [] if not self.chain_type.has_d else [gen.d_allele.name]
        container.j_call = [gen.j_allele.name]
        container.c_call = [gen.c_allele.name] if getattr(gen, 'c_allele', None) else [None]
        container.mutation_rate = gen.mutation_freq
        container.mutations = {pos: gen.mutations[pos] for pos in sorted(gen.mutations)}
        container.productive = gen.functional
        container.stop_codon = gen.stop_codon
        container.vj_in_frame = gen.vj_in_frame
        container.note = gen.note

    def apply(self, container: SimulationContainer) -> None:
        """
        Applies the simulation step to generate a heavy chain sequence.
        """
        self.simulate_sequence(container)

    def get_graph_node(self):
        """Generates a detailed GraphViz node representation with constructor details in an HTML-like format."""
        step_name = "Simulate Sequence"
        mutation_model_name = self.mutation_model.__class__.__name__
        # Constructing an HTML-like label using a table for detailed formatting
        label = f"""
        <TABLE BORDER="0" CELLBORDER="1" CELLSPACING="0" CELLPADDING="4">
        <TR><TD COLSPAN="2" BGCOLOR="lightsteelblue"><B>{step_name}</B></TD></TR>
        <TR><TD ALIGN="LEFT"><B>Chain Type</B></TD><TD ALIGN="LEFT">{ChainType[self.chain_type].name}</TD></TR>
        <TR><TD ALIGN="LEFT"><B>Mutation Model</B></TD><TD ALIGN="LEFT">{mutation_model_name}</TD></TR>
        <TR><TD ALIGN="LEFT"><B>Productive</B></TD><TD ALIGN="LEFT">{self.productive}</TD></TR>
        <TR><TD ALIGN="LEFT"><B>Specific V Allele</B></TD><TD ALIGN="LEFT">{self.specific_v}</TD></TR>
        <TR><TD ALIGN="LEFT"><B>Specific D Allele</B></TD><TD ALIGN="LEFT">{self.specific_d}</TD></TR>
        <TR><TD ALIGN="LEFT"><B>Specific J Allele</B></TD><TD ALIGN="LEFT">{self.specific_j}</TD></TR>

        </TABLE>
        
        """

        return label, 'box', "filled,rounded", SIMULATION_STEP_BOX_COLOR, "Helvetica", "black"
