from ..container.SimulationContainer import SimulationContainer
from ..pipeline.plot_parameters import CORRUPTION_STEP_BOX_COLOR
from ..steps.StepBase import AugmentationStep


class EnforceSequenceLength(AugmentationStep):
    """
    Enforces a maximum sequence length by trimming from the 5' end.

    This step simulates sequencing platform read length limits. When a biological
    sequence exceeds the maximum read length of a sequencing platform (e.g.,
    Illumina 300bp or 500bp reads), this step truncates from the 5' end to fit
    the platform constraints.

    This step should typically be placed AFTER CorruptSequenceBeginning in the
    pipeline, as it handles a separate concern (platform limits) from corruption
    events (sequencing artifacts).

    Args:
        max_length: Maximum allowed sequence length. Default: 576

    Example:
        # Use default (576bp - typical for many pipelines)
        step = EnforceSequenceLength()

        # Simulate Illumina 300bp reads
        step = EnforceSequenceLength(max_length=300)

        # Simulate longer reads
        step = EnforceSequenceLength(max_length=500)
    """

    def __init__(
        self,
        *,
        max_length: int = 576,
    ):
        super().__init__()
        self.max_length = max_length

    def apply(self, container: SimulationContainer) -> None:
        """
        Applies the sequence length enforcement.

        If the sequence exceeds max_length, trims from the 5' end and updates
        all position metadata accordingly.
        """
        current_length = len(container.sequence)

        if current_length <= self.max_length:
            return  # Nothing to do

        # Calculate how much to trim from the 5' end
        trim_amount = current_length - self.max_length

        # Trim the sequence
        container.sequence = container.sequence[trim_amount:]

        # Update mutation positions - remove mutations that were trimmed away,
        # adjust positions for remaining mutations
        container.mutations = {
            pos - trim_amount: val
            for pos, val in container.mutations.items()
            if pos >= trim_amount
        }

        # Update Ns positions - remove trimmed entries, shift remaining
        container.Ns = {
            pos - trim_amount: val
            for pos, val in container.Ns.items()
            if pos >= trim_amount
        }

        # Update indel positions - remove trimmed entries, shift remaining
        container.indels = {
            pos - trim_amount: val
            for pos, val in container.indels.items()
            if pos >= trim_amount
        }

        # Save original v_sequence_start before shifting
        original_v_start = container.v_sequence_start

        # Adjust all position metadata
        container.v_sequence_start = max(0, container.v_sequence_start - trim_amount)
        container.v_sequence_end = max(0, container.v_sequence_end - trim_amount)
        container.d_sequence_start = max(0, container.d_sequence_start - trim_amount)
        container.d_sequence_end = max(0, container.d_sequence_end - trim_amount)
        container.j_sequence_start = max(0, container.j_sequence_start - trim_amount)
        container.j_sequence_end = max(0, container.j_sequence_end - trim_amount)
        container.junction_sequence_start = max(0, container.junction_sequence_start - trim_amount)
        container.junction_sequence_end = max(0, container.junction_sequence_end - trim_amount)

        # Adjust germline positions only if we trimmed into the V region
        if trim_amount > original_v_start:
            container.v_germline_start = container.v_germline_start + (trim_amount - original_v_start)

        # Track what we trimmed (useful for debugging/analysis)
        if not hasattr(container, 'length_enforcement_trim') or container.length_enforcement_trim is None:
            container.length_enforcement_trim = trim_amount
        else:
            container.length_enforcement_trim += trim_amount

    def get_graph_node(self):
        """Generates a detailed GraphViz node representation with constructor details in an HTML-like format."""
        step_name = "Enforce Sequence Length"

        label = f"""
        <TABLE BORDER="0" CELLBORDER="1" CELLSPACING="0" CELLPADDING="4">
        <TR><TD COLSPAN="2" BGCOLOR="lightsteelblue"><B>{step_name}</B></TD></TR>
        <TR><TD ALIGN="LEFT"><B>Max Length</B></TD><TD ALIGN="LEFT">{self.max_length}</TD></TR>
        <TR><TD ALIGN="LEFT"><B>Description</B></TD><TD ALIGN="LEFT">Trims 5' end to fit sequencing read length</TD></TR>
        </TABLE>

        """

        return label, 'box', "filled,rounded", CORRUPTION_STEP_BOX_COLOR, "Helvetica", "black"
