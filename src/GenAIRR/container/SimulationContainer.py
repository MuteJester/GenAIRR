class SimulationContainer:
    def __init__(self, sequence_instance=None):
        """
        Initializes the SimulationContainer. If a sequence_instance is not provided,
        all values are set to None by default.

        Args:
            sequence_instance (BaseSequence, optional): An instance of a class inherited from BaseSequence.
        """
        if sequence_instance:
            # Extract sequence and related properties from the sequence instance
            self.from_instance(sequence_instance)
            self.sequence_instance = sequence_instance
        else:
            # Initialize all attributes to None if no sequence_instance is provided
            self.sequence = 'None'
            self.sequence_instance = None
            self.v_call = []
            self.d_call = []
            self.j_call = []
            self.c_call = []

            # Position metadata
            self.v_sequence_start = 0
            self.v_sequence_end = 0
            self.d_sequence_start = 0
            self.d_sequence_end = 0
            self.j_sequence_start = 0
            self.j_sequence_end = 0

            # Germline positions
            self.v_germline_start = 0
            self.v_germline_end = 0
            self.d_germline_start = 0
            self.d_germline_end = 0
            self.j_germline_start = 0
            self.j_germline_end = 0

            # Junction positions
            self.junction_sequence_start = 0
            self.junction_sequence_end = 0

            # Mutation and trimming information
            self.mutation_rate = 0
            self.mutations = dict()
            self.indels = dict()
            self.Ns = dict()

            # Trimming data
            self.v_trim_5 = 0
            self.v_trim_3 = 0
            self.d_trim_5 = 0
            self.d_trim_3 = 0
            self.j_trim_5 = 0
            self.j_trim_3 = 0
            self.c_trim_3 = 0

            # Productivity and assessment flags
            self.productive = None
            self.stop_codon = None
            self.vj_in_frame = None
            self.note = ''

            # Corruption metadata
            self.corruption_event = 'no-corruption'
            self.corruption_add_amount = 0
            self.corruption_remove_amount = 0
            self.corruption_removed_section = ''
            self.corruption_added_section = ''

    def from_instance(self,sequence_instance):
        self.sequence = sequence_instance.mutated_seq
        self.v_call = [sequence_instance.v_allele.name]
        self.d_call = [] if sequence_instance.d_allele is None else [sequence_instance.d_allele.name]
        self.j_call = [sequence_instance.j_allele.name]
        self.c_call = [sequence_instance.c_allele.name] if getattr(sequence_instance, "c_allele", None) else [None]

        # Position metadata
        self.v_sequence_start = sequence_instance.v_seq_start
        self.v_sequence_end = sequence_instance.v_seq_end
        self.d_sequence_start = sequence_instance.d_seq_start
        self.d_sequence_end = sequence_instance.d_seq_end
        self.j_sequence_start = sequence_instance.j_seq_start
        self.j_sequence_end = sequence_instance.j_seq_end

        # Germline positions
        self.v_germline_start = sequence_instance.v_germline_start
        self.v_germline_end = sequence_instance.v_germline_end
        self.d_germline_start = sequence_instance.d_germline_start if hasattr(sequence_instance,'d_germline_start') else None
        self.d_germline_end = sequence_instance.d_germline_end if hasattr(sequence_instance,'d_germline_end') else None
        self.j_germline_start = sequence_instance.j_germline_start
        self.j_germline_end = sequence_instance.j_germline_end

        # Junction positions
        self.junction_sequence_start = sequence_instance.junction_start
        self.junction_sequence_end = sequence_instance.junction_end

        # Mutation and trimming information
        self.mutation_rate = sequence_instance.mutation_freq
        self.mutations = {pos: sequence_instance.mutations[pos] for pos in sorted(sequence_instance.mutations)}
        self.indels = {}  # Initialize empty; can be updated during augmentation
        self.Ns = {}  # Initialize empty; can be updated during augmentation

        # Trimming data
        self.v_trim_5 = sequence_instance.v_trim_5
        self.v_trim_3 = sequence_instance.v_trim_3
        self.d_trim_5 = sequence_instance.d_trim_5 if hasattr(sequence_instance,'d_trim_5') else None
        self.d_trim_3 = sequence_instance.d_trim_3 if hasattr(sequence_instance,'d_trim_3') else None
        self.j_trim_5 = sequence_instance.j_trim_5
        self.j_trim_3 = sequence_instance.j_trim_3
        self.c_trim_3 = sequence_instance.c_trim_3 if getattr(sequence_instance, "c_trim_3", None) else None

        # Productivity and assessment flags
        self.productive = sequence_instance.functional
        self.stop_codon = sequence_instance.stop_codon
        self.vj_in_frame = sequence_instance.vj_in_frame
        self.note = sequence_instance.note

        # Corruption metadata
        self.corruption_event = 'no-corruption'
        self.corruption_add_amount = 0
        self.corruption_remove_amount = 0
        self.corruption_removed_section = ''
        self.corruption_added_section = ''

    def add_mutation(self, position, mutation):
        """Add a mutation log entry."""
        self.mutations[position] = mutation

    def add_indel(self, position, indel):
        """Add an indel log entry."""
        self.indels[position] = indel

    def add_N_insertion(self, position, original_base):
        """Add an 'N' insertion log entry."""
        self.Ns[position] = f"{original_base} > N"

    def update_sequence(self, new_sequence):
        """Update the sequence and adjust related attributes as needed."""
        self.sequence = new_sequence

    def shift_positions(self, shift_amount):
        """Shift all start and end positions by a given amount (used for additions/removals)."""
        self.v_sequence_start += shift_amount
        self.v_sequence_end += shift_amount
        self.d_sequence_start += shift_amount
        self.d_sequence_end += shift_amount
        self.j_sequence_start += shift_amount
        self.j_sequence_end += shift_amount
        self.junction_sequence_start += shift_amount
        self.junction_sequence_end += shift_amount

    def get_original_index(self, corrupted_index):
        """
        Calculates the original index in the uncorrupted sequence corresponding to a given index in the corrupted sequence.

        Args:
            corrupted_index (int): The index in the corrupted sequence.
            simulated (dict): A dictionary containing simulation metadata, including details of corruption events and amounts.

        Returns:
            int or None: The original index in the uncorrupted sequence, or None if the index cannot be mapped back due to the nature of the corruption.
        """
        corruption_event = self.corruption_event
        amount_added = self.corruption_add_amount
        amount_removed = self.corruption_remove_amount

        if corruption_event == 'add':
            # If characters were added at the beginning, the original index is shifted to the right
            return corrupted_index - amount_added if corrupted_index >= amount_added else None
        elif corruption_event == 'remove':
            # If characters were removed from the beginning, the original index is shifted to the left
            return corrupted_index + amount_removed
        elif corruption_event == 'remove_before_add':
            # Combined effect of removal and addition
            adjusted_index = corrupted_index + amount_removed
            return adjusted_index - amount_added if adjusted_index >= amount_added else None
        else:
            # If no corruption or unknown corruption event, assume the index is unchanged
            return corrupted_index

    def update_from_dict(self, updates):
        """
        Updates the attributes of the container from a dictionary.

        Args:
            updates (dict): A dictionary where keys match attribute names of the container.
        """
        for key, value in updates.items():
            if hasattr(self, key):
                setattr(self, key, value)
            else:
                raise KeyError(f"'{key}' not found in the SimulationContainer.")


    def get_dict(self):
        """Return a summary dictionary similar."""
        return {
            "sequence": self.sequence,
            "v_call": self.v_call,
            "d_call": self.d_call,
            "j_call": self.j_call,
            "c_call": self.c_call,
            "v_sequence_start": self.v_sequence_start,
            "v_sequence_end": self.v_sequence_end,
            "d_sequence_start": self.d_sequence_start,
            "d_sequence_end": self.d_sequence_end,
            "j_sequence_start": self.j_sequence_start,
            "j_sequence_end": self.j_sequence_end,
            "v_germline_start": self.v_germline_start,
            "v_germline_end": self.v_germline_end,
            "d_germline_start": self.d_germline_start,
            "d_germline_end": self.d_germline_end,
            "j_germline_start": self.j_germline_start,
            "j_germline_end": self.j_germline_end,
            "junction_sequence_start": self.junction_sequence_start,
            "junction_sequence_end": self.junction_sequence_end,
            "mutation_rate": self.mutation_rate,
            "mutations": self.mutations,
            "indels": self.indels,
            "Ns": self.Ns,
            "v_trim_5": self.v_trim_5,
            "v_trim_3": self.v_trim_3,
            "d_trim_5": self.d_trim_5,
            "d_trim_3": self.d_trim_3,
            "j_trim_5": self.j_trim_5,
            "j_trim_3": self.j_trim_3,
            "c_trim_3": self.c_trim_3,
            "productive": self.productive,
            "stop_codon": self.stop_codon,
            "vj_in_frame": self.vj_in_frame,
            "note": self.note,
            "corruption_event": self.corruption_event,
            "corruption_add_amount": self.corruption_add_amount,
            "corruption_remove_amount": self.corruption_remove_amount,
            "corruption_removed_section": self.corruption_removed_section,
            "corruption_added_section": self.corruption_added_section,
        }

    def __getitem__(self, key):
        """Allow dictionary-like access to attributes."""
        try:
            return getattr(self, key)
        except AttributeError:
            raise KeyError(f"'{key}' not found in the SimulationContainer.")

    def __setitem__(self, key, value):
        """Allow setting attributes using dictionary-like access."""
        if hasattr(self, key):
            setattr(self, key, value)
        else:
            raise KeyError(f"'{key}' not found in the SimulationContainer.")
