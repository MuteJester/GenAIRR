"""
S5F Mutation Model

Context-dependent 5-mer mutation model that uses position-specific mutation
probabilities based on the surrounding sequence context.
"""

import random
import pickle
from typing import Dict, Optional, Tuple, List, Set

import pandas as pd

from .mutation_model import MutationModel
from .context import MutationContext, ChainType
from .result import MutationResult
from .components import StopCodonChecker, AnchorChecker


class Nucleotide:
    """
    Represents a single nucleotide in a DNA sequence.

    Attributes:
        value: The nucleotide character (A, T, C, G).
        adjacent: List of FiveMER objects containing this nucleotide.
    """

    def __init__(self, value: str):
        """Initialize a Nucleotide with a given value."""
        self.value = value
        self.adjacent: List["FiveMER"] = []

    def update_value(self, new_value: str, update_callback=None):
        """
        Update the nucleotide's value and trigger callbacks for adjacent 5-mers.

        Args:
            new_value: The new nucleotide value.
            update_callback: Optional callback for each adjacent FiveMER.
        """
        self.value = new_value
        if update_callback:
            for five_mer in self.adjacent:
                update_callback(five_mer)

    def __repr__(self) -> str:
        return self.value


class FiveMER:
    """
    Represents a 5-mer (sequence of 5 nucleotides) in a DNA sequence.

    The S5F model uses 5-mer context to determine mutation probabilities.
    Each 5-mer has a likelihood (mutability) based on the sequence context.

    Attributes:
        nucleotides: List of Nucleotide objects making up the 5-mer.
        sequence: The nucleotide sequence string.
        position: Position within the larger sequence.
        likelihood: Probability of this 5-mer being mutated.
        modified: Whether this 5-mer has been modified.
    """

    def __init__(self, nucleotides: List[Nucleotide]):
        """Initialize a FiveMER with a list of Nucleotide objects."""
        self.nucleotides = nucleotides
        for nuc in self.nucleotides:
            nuc.adjacent.append(self)
        self.sequence = "".join([nuc.value for nuc in nucleotides])
        self.position: Optional[int] = None
        self.likelihood: float = 0.0
        self.modified: bool = False

    def update_sequence(self, mutability: Optional[Dict] = None):
        """Update the 5-mer sequence and optionally its mutability likelihood."""
        self.sequence = "".join([nuc.value for nuc in self.nucleotides])
        if mutability is not None:
            self.likelihood = mutability.get(self.sequence, self.likelihood)

    def change_center(self, new_value: str, mutability: Optional[Dict] = None):
        """
        Change the central nucleotide and update sequences and likelihoods.

        Args:
            new_value: New value for the central nucleotide.
            mutability: Dictionary mapping 5-mer sequences to mutability likelihoods.
        """
        center_nucleotide = self.nucleotides[2]  # 0-indexed, 2 is center
        center_nucleotide.update_value(
            new_value, lambda fm: fm.update_sequence(mutability)
        )
        self.modified = True

    def __repr__(self) -> str:
        return "".join([i.value for i in self.nucleotides])

    def __eq__(self, other) -> bool:
        if isinstance(other, FiveMER):
            return self.sequence == other.sequence
        elif isinstance(other, str):
            return self.sequence == other
        raise TypeError(f"Unsupported comparison between FiveMER and {type(other)}")

    @staticmethod
    def create_five_mers(
        dna_sequence: str, mutability: Optional[Dict] = None
    ) -> List["FiveMER"]:
        """
        Create a list of FiveMER objects from a DNA sequence.

        Args:
            dna_sequence: The DNA sequence to create 5-mers from.
            mutability: Dictionary mapping 5-mer sequences to mutability likelihoods.

        Returns:
            List of FiveMER objects.
        """
        padded_sequence = f"NN{dna_sequence}NN"
        sequence_length = len(dna_sequence)

        five_mers = []
        for i in range(sequence_length):
            five_mer_nucleotides = padded_sequence[i : i + 5]
            five_mer = FiveMER([Nucleotide(nuc) for nuc in five_mer_nucleotides])
            five_mer.position = i

            if mutability is not None:
                five_mer.likelihood = mutability.get(five_mer_nucleotides, 0)

            five_mers.append(five_mer)

        return five_mers

    @staticmethod
    def five_mers_to_dna(five_mers: List["FiveMER"]) -> str:
        """Convert a list of FiveMER objects back into a DNA sequence string."""
        if not five_mers:
            return ""

        dna_sequence = five_mers[0].nucleotides[2].value

        for five_mer in five_mers[1:-1]:
            dna_sequence += five_mer.nucleotides[2].value

        dna_sequence += five_mers[-1].nucleotides[2].value

        return dna_sequence


class S5F(MutationModel):
    """
    S5F (Somatic 5-mer Frequency) mutation model.

    This context-dependent mutation model uses 5-mer sequence context to
    determine mutation probabilities. Different metadata files are used
    for heavy chains vs. light chains.

    Args:
        min_mutation_rate: Minimum mutation rate (0.0 to 1.0). Default: 0.0
        max_mutation_rate: Maximum mutation rate (0.0 to 1.0). Default: 0.0
        custom_model: Path to custom mutation model pickle file. Default: None
        productive: If True, avoid stop codons and preserve anchors. Default: False

    Example:
        >>> from GenAIRR import S5F, set_seed
        >>> set_seed(42)
        >>> model = S5F(min_mutation_rate=0.003, max_mutation_rate=0.25)

    Note:
        The S5F model caches metadata per chain type (heavy vs light) to avoid
        the stateful bug where reusing the same instance for different chain
        types would use incorrect metadata.
    """

    NUCLEOTIDES: Set[str] = {"A", "T", "C", "G"}

    # Class-level cache for metadata (shared across instances)
    _metadata_cache: Dict[str, Tuple[Dict, Dict, Dict]] = {}

    def __init__(
        self,
        min_mutation_rate: float = 0.0,
        max_mutation_rate: float = 0.0,
        custom_model: Optional[str] = None,
        productive: bool = False,
    ):
        """Initialize the S5F mutation model."""
        super().__init__(
            min_mutation_rate=min_mutation_rate,
            max_mutation_rate=max_mutation_rate,
            productive=productive,
        )
        self.bases = self.NUCLEOTIDES  # Legacy attribute
        self.custom_model = custom_model

        # Shared components
        self._stop_checker = StopCodonChecker()
        self._anchor_checker = AnchorChecker()

        # Per-chain-type metadata (loaded lazily)
        self._mutability: Optional[Dict] = None
        self._substitution: Optional[Dict] = None
        self._targeting: Optional[Dict] = None
        self._current_chain_type: Optional[str] = None

    @property
    def name(self) -> str:
        """Human-readable name of this mutation model."""
        return "S5F"

    @property
    def supports_productive_mode(self) -> bool:
        """This model supports productive mode."""
        return True

    @property
    def is_context_dependent(self) -> bool:
        """S5F uses 5-mer context for mutation probabilities."""
        return True

    @property
    def requires_external_data(self) -> bool:
        """S5F requires loading external pickle files."""
        return True

    def _get_chain_category(self, chain_type: ChainType) -> str:
        """Map chain type to metadata category (heavy or light)."""
        heavy_types = {ChainType.BCR_HEAVY, ChainType.TCR_BETA, ChainType.TCR_DELTA}
        return "heavy" if chain_type in heavy_types else "light"

    def _load_metadata_for_chain(self, chain_category: str) -> Tuple[Dict, Dict, Dict]:
        """
        Load mutation model metadata for a chain category.

        Args:
            chain_category: Either "heavy" or "light".

        Returns:
            Tuple of (mutability, substitution, targeting) dictionaries.
        """
        cache_key = self.custom_model or chain_category

        if cache_key in S5F._metadata_cache:
            return S5F._metadata_cache[cache_key]

        from importlib import resources

        if self.custom_model is not None:
            with open(self.custom_model, "rb") as h:
                mutability, substitution, targeting = pickle.load(h)
        else:
            if chain_category == "heavy":
                pkg = "GenAIRR.data.mutation_model_parameters"
                filename = "HH_S5F_META.pkl"
            else:
                pkg = "GenAIRR.data.mutation_model_parameters"
                filename = "HKL_S5F_META.pkl"

            with resources.path(pkg, filename) as data_path:
                with open(data_path, "rb") as h:
                    mutability, substitution, targeting = pickle.load(h)

        # Convert substitution DataFrame to dict and clean NaN values
        if hasattr(substitution, "to_dict"):
            substitution = substitution.to_dict(orient="dict")
            substitution = {
                outer_key: {
                    inner_key: inner_value
                    for inner_key, inner_value in outer_dict.items()
                    if not pd.isna(inner_value)
                }
                for outer_key, outer_dict in substitution.items()
            }

        # Cache for future use
        S5F._metadata_cache[cache_key] = (mutability, substitution, targeting)

        return mutability, substitution, targeting

    def _ensure_metadata_loaded(self, chain_type: ChainType):
        """Ensure metadata is loaded for the given chain type."""
        chain_category = self._get_chain_category(chain_type)

        if self._current_chain_type != chain_category:
            self._mutability, self._substitution, self._targeting = (
                self._load_metadata_for_chain(chain_category)
            )
            self._current_chain_type = chain_category

    # Termination constants
    MAX_CONSECUTIVE_FAILURES = 500  # Exit if stuck for this many attempts
    TOTAL_ATTEMPT_MULTIPLIER = 100  # Hard cap: target * this multiplier

    def apply(self, context: MutationContext) -> MutationResult:
        """
        Apply S5F mutations to a sequence.

        Args:
            context: MutationContext containing sequence and region information.

        Returns:
            MutationResult with the mutated sequence and mutation details.
        """
        # Load appropriate metadata for chain type
        self._ensure_metadata_loaded(context.chain_type)

        sequence = context.sequence
        mutation_rate = random.uniform(self.min_mutation_rate, self.max_mutation_rate)
        target_mutations = int(mutation_rate * len(sequence))

        if target_mutations == 0:
            return MutationResult(
                mutated_sequence=sequence,
                mutations={},
                mutation_rate=mutation_rate,
                mutation_count=0,
            )

        # Create 5-mers
        five_mers = FiveMER.create_five_mers(sequence, self._mutability)
        naive_sequence = list(sequence)

        if self.productive:
            mutations, five_mers = self._apply_productive_mutations(
                context,
                five_mers,
                naive_sequence,
                target_mutations,
            )
        else:
            mutations, five_mers = self._apply_standard_mutations(
                five_mers,
                naive_sequence,
                target_mutations,
            )

        mutated_sequence = FiveMER.five_mers_to_dna(five_mers)

        # Return actual achieved rate (may differ from target if we hit limits)
        actual_rate = len(mutations) / len(sequence) if len(sequence) > 0 else 0.0

        return MutationResult(
            mutated_sequence=mutated_sequence,
            mutations=mutations,
            mutation_rate=actual_rate,
            mutation_count=len(mutations),
        )

    def _apply_standard_mutations(
        self,
        five_mers: List[FiveMER],
        naive_sequence: List[str],
        target_mutations: int,
    ) -> Tuple[Dict, List[FiveMER]]:
        """
        Apply mutations without productive constraints.

        Uses fixed target with termination guarantees:
        - Exits if consecutive failures exceed MAX_CONSECUTIVE_FAILURES
        - Exits if total attempts exceed target * TOTAL_ATTEMPT_MULTIPLIER
        """
        mutations = {}
        consecutive_failures = 0
        total_attempts = 0
        max_total_attempts = max(target_mutations * self.TOTAL_ATTEMPT_MULTIPLIER, 1000)

        while len(mutations) < target_mutations:
            total_attempts += 1

            # Hard safety cap
            if total_attempts >= max_total_attempts:
                break

            sampled_position, chosen_index = self._weighted_choice(five_mers)

            if sampled_position is None:
                consecutive_failures += 1
                if consecutive_failures >= self.MAX_CONSECUTIVE_FAILURES:
                    break
                continue

            # Get substitution probabilities for this 5-mer
            substitutions = self._substitution.get(sampled_position.sequence, {})
            if not substitutions:
                consecutive_failures += 1
                if consecutive_failures >= self.MAX_CONSECUTIVE_FAILURES:
                    break
                continue

            mutable_bases = list(substitutions.keys())
            bases_likelihoods = list(substitutions.values())
            mutation_to_apply = random.choices(
                mutable_bases, weights=bases_likelihoods, k=1
            )[0]

            # Log mutation
            pos = sampled_position.position
            if pos not in mutations:
                mutations[pos] = f"{sampled_position.sequence[2]}>{mutation_to_apply}"
            else:
                mutations[pos] += f">{mutation_to_apply}"

            # Remove if reverted to naive
            if mutation_to_apply == naive_sequence[pos]:
                mutations.pop(pos, None)

            # Apply mutation
            sampled_position.change_center(mutation_to_apply, self._mutability)

            # Reset consecutive failures on successful mutation application
            consecutive_failures = 0

        return mutations, five_mers

    def _apply_productive_mutations(
        self,
        context: MutationContext,
        five_mers: List[FiveMER],
        naive_sequence: List[str],
        target_mutations: int,
    ) -> Tuple[Dict, List[FiveMER]]:
        """
        Apply mutations while ensuring sequence remains productive.

        Uses fixed target with termination guarantees:
        - Exits if consecutive failures exceed MAX_CONSECUTIVE_FAILURES
        - Exits if total attempts exceed target * TOTAL_ATTEMPT_MULTIPLIER
        - Post-mutation verification ensures no stop codons slip through
        """
        mutations = {}
        consecutive_failures = 0
        total_attempts = 0
        max_total_attempts = max(target_mutations * self.TOTAL_ATTEMPT_MULTIPLIER, 1000)

        reading_frame = [[2, 1, 0][idx % 3] for idx in range(len(five_mers))]
        restricted_positions = self._anchor_checker.get_restricted_positions(
            context.junction_start, context.junction_end
        )

        while len(mutations) < target_mutations:
            total_attempts += 1

            # Hard safety cap
            if total_attempts >= max_total_attempts:
                break

            # Try to find a valid productive mutation
            result = self._try_productive_mutation(
                five_mers, reading_frame, restricted_positions
            )

            if result is None:
                consecutive_failures += 1
                if consecutive_failures >= self.MAX_CONSECUTIVE_FAILURES:
                    break
                continue

            sampled_position, mutation_to_apply = result
            pos = sampled_position.position

            # Save the original base in case we need to undo
            original_base = sampled_position.sequence[2]

            # Apply mutation
            sampled_position.change_center(mutation_to_apply, self._mutability)

            # Post-mutation verification: check if a stop codon was created
            current_sequence = FiveMER.five_mers_to_dna(five_mers)
            if self._stop_checker.has_stop_codon(current_sequence):
                # Undo the mutation
                sampled_position.change_center(original_base, self._mutability)
                consecutive_failures += 1
                if consecutive_failures >= self.MAX_CONSECUTIVE_FAILURES:
                    break
                continue

            # Mutation verified safe - log it
            if pos not in mutations:
                mutations[pos] = f"{original_base}>{mutation_to_apply}"
            else:
                mutations[pos] += f">{mutation_to_apply}"

            # Remove if reverted to naive
            if mutation_to_apply == naive_sequence[pos]:
                mutations.pop(pos, None)

            # Reset consecutive failures on successful mutation
            consecutive_failures = 0

        return mutations, five_mers

    def _weighted_choice(
        self, five_mers: List[FiveMER]
    ) -> Tuple[Optional[FiveMER], int]:
        """Select a FiveMER weighted by likelihood."""
        total_weight = 0.0
        cumulative_weights = []

        for fm in five_mers:
            weight = fm.likelihood if fm.likelihood == fm.likelihood else 0
            total_weight += weight
            cumulative_weights.append(total_weight)

        if total_weight == 0:
            return None, -1

        r = random.uniform(0, total_weight)
        for i, cumulative_weight in enumerate(cumulative_weights):
            if r < cumulative_weight:
                return five_mers[i], i

        return five_mers[-1], len(five_mers) - 1

    def _try_productive_mutation(
        self,
        five_mers: List[FiveMER],
        reading_frame: List[int],
        restricted_positions: Dict[int, str],
    ) -> Optional[Tuple[FiveMER, str]]:
        """
        Try to find a single valid mutation that doesn't create stop codons or destroy anchors.

        Returns:
            Tuple of (FiveMER, new_base) if successful, None if no valid mutation found.
        """
        sampled_position, chosen_index = self._weighted_choice(five_mers)

        if sampled_position is None:
            return None

        # Get bases that would form stop codons
        stop_forming = self._stop_checker.get_stop_forming_bases(
            sampled_position.sequence
        )

        # Get substitution probabilities
        substitutions = self._substitution.get(sampled_position.sequence, {})

        # Remove stop-forming bases
        filtered_subs = {
            base: prob
            for base, prob in substitutions.items()
            if base not in stop_forming
        }

        if not filtered_subs:
            return None

        # Normalize
        total = sum(filtered_subs.values())
        filtered_subs = {base: prob / total for base, prob in filtered_subs.items()}

        mutable_bases = list(filtered_subs.keys())
        bases_likelihoods = list(filtered_subs.values())

        mutation_to_apply = random.choices(
            mutable_bases, weights=bases_likelihoods, k=1
        )[0]

        # Check stop codon in reading frames
        nucs = [str(nuc) for nuc in sampled_position.nucleotides]
        rf = reading_frame[chosen_index]
        is_stop, _ = self._stop_checker.would_create_stop_any_frame(
            nucs, 2, mutation_to_apply, rf
        )

        if is_stop:
            return None

        # Check anchor preservation
        if chosen_index in restricted_positions:
            anchor_type = restricted_positions[chosen_index]
            nucs_copy = nucs[0:2] + [mutation_to_apply] + nucs[3:]
            codon = "".join(nucs_copy[rf : rf + 3])

            if not self._anchor_checker.is_valid_anchor_mutation(codon, anchor_type):
                return None

        return sampled_position, mutation_to_apply

    # =========================================================================
    # Legacy interface (for backwards compatibility)
    # =========================================================================

    def load_metadata(self, sequence_object):
        """
        Load metadata based on sequence object type (legacy interface).

        Deprecated: Metadata is now loaded automatically based on chain type.
        """
        from ..sequence import HeavyChainSequence
        from ..sequence.light_chain import LightChainSequence

        if self.custom_model is not None:
            chain_category = "custom"
        elif isinstance(sequence_object, HeavyChainSequence):
            chain_category = "heavy"
        elif isinstance(sequence_object, LightChainSequence):
            chain_category = "light"
        else:
            raise ValueError("Unsupported Sequence Type")

        self._mutability, self._substitution, self._targeting = (
            self._load_metadata_for_chain(chain_category)
        )
        self._current_chain_type = chain_category

    def apply_mutation(self, sequence_object):
        """
        Apply mutations to a sequence object (legacy interface).

        This method is maintained for backwards compatibility.
        New code should use the `apply` method with MutationContext instead.
        """
        # Create context from sequence object
        context = MutationContext.from_sequence_object(sequence_object)

        # Call the new interface
        result = self.apply(context)

        # Return in legacy format
        return result.mutated_sequence, result.mutations, result.mutation_rate

    # Keep static method for backwards compatibility
    @staticmethod
    def weighted_choice(five_mers):
        """Legacy static method - use instance method _weighted_choice instead."""
        total_weight = 0
        cumulative_weights = []

        for fm in five_mers:
            weight = fm.likelihood if fm.likelihood == fm.likelihood else 0
            total_weight += weight
            cumulative_weights.append(total_weight)

        if total_weight == 0:
            return None, -1

        r = random.uniform(0, total_weight)
        for i, cumulative_weight in enumerate(cumulative_weights):
            if r < cumulative_weight:
                return five_mers[i], i

        return five_mers[-1], len(five_mers) - 1
