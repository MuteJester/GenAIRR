from ..container.SimulationContainer import SimulationContainer
from ..pipeline.plot_parameters import SIMULATION_STEP_BOX_COLOR
from ..steps.StepBase import AugmentationStep
from ..utilities import DataConfig
from ..sequence import HeavyChainSequence


class SimulateSequence(AugmentationStep):
    def __init__(self,mutation_model, productive=False, specific_v=None, specific_d=None, specific_j=None):
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

    def simulate_sequence(self, container: SimulationContainer):
        """
        Simulates a heavy chain sequence and populates the SimulationContainer with relevant data.

        Args:
            container (SimulationContainer): The container to store the simulated sequence data.
        """
        gen = HeavyChainSequence.create_random(self.dataconfig, specific_v=self.specific_v,
                                               specific_d=self.specific_d, specific_j=self.specific_j)

        # Ensure productivity if specified
        if self.productive:
            while not gen.functional:
                gen = HeavyChainSequence.create_random(self.dataconfig, specific_v=self.specific_v,
                                                       specific_d=self.specific_d, specific_j=self.specific_j)

        container.from_instance(gen)
        # Apply mutation
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
        container.d_germline_start = gen.d_germline_start
        container.d_germline_end = gen.d_germline_end
        container.j_germline_start = gen.j_germline_start
        container.j_germline_end = gen.j_germline_end
        container.junction_sequence_start = gen.junction_start
        container.junction_sequence_end = gen.junction_end
        container.v_call = [gen.v_allele.name]
        container.d_call = [gen.d_allele.name]
        container.j_call = [gen.j_allele.name]
        container.c_call = [gen.c_allele.name]
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
        <TR><TD ALIGN="LEFT"><B>Mutation Model</B></TD><TD ALIGN="LEFT">{mutation_model_name}</TD></TR>
        <TR><TD ALIGN="LEFT"><B>Productive</B></TD><TD ALIGN="LEFT">{self.productive}</TD></TR>
        <TR><TD ALIGN="LEFT"><B>Specific V Allele</B></TD><TD ALIGN="LEFT">{self.specific_v}</TD></TR>
        <TR><TD ALIGN="LEFT"><B>Specific D Allele</B></TD><TD ALIGN="LEFT">{self.specific_d}</TD></TR>
        <TR><TD ALIGN="LEFT"><B>Specific J Allele</B></TD><TD ALIGN="LEFT">{self.specific_j}</TD></TR>

        </TABLE>
        
        """

        return label, 'box', "filled,rounded", SIMULATION_STEP_BOX_COLOR, "Helvetica", "black"