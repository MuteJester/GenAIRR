from __future__ import annotations

from typing import Any, TYPE_CHECKING

# Imported for the quoted return annotations below only. At runtime this
# module is imported BY the class it annotates, so a real import would be
# circular. Under TYPE_CHECKING a type checker and an IDE resolve the
# return type instead of seeing an undefined name.
#
# This does NOT make typing.get_type_hints() work on these signatures: the
# import does not exist at runtime, so it still raises NameError, as it did
# before. Fixing that needs the name importable at runtime, which is the
# circularity above, or an explicit namespace at the call site.
if TYPE_CHECKING:
    from ..experiment import Experiment


class _ConstraintsMixin:
    __slots__ = ()

    def with_metadata(self, **fields: Any) -> "Experiment":
        """Attach sample-level metadata to every AIRR record.

        Standard AIRR Repertoire fields like ``sample_id``,
        ``donor``, ``repertoire_id``, and ``cell_id`` are commonly
        used by multi-sample analysis tools — Change-O, scirpy,
        immcantation pipelines all expect them populated. Custom
        keys are also accepted and pass through unchanged.

        Generated sequence/annotation fields are reserved and cannot
        be overwritten. Subsequent ``with_metadata()`` calls **merge** with the
        prior set (per-key replacement). Pass ``key=None`` to clear
        a single key. Values are stored as-is and serialised via
        the standard CSV/TSV writers; non-string values are
        converted with ``str()`` on output.

        Example::

            (Experiment.on("human_igh")
             .with_metadata(sample_id="P1", donor="D001",
                            repertoire_id="R-001-IGH")
             .recombine().mutate(count=8))
        """
        from .._result_export import _DEFAULT_COLUMN_ORDER

        reserved = set(_DEFAULT_COLUMN_ORDER) | {
            "coding_end_processing",
            "duplex_repair", "sequence_length",
            "productivity_policy", "productivity_scope", "productivity_reasons",
            "p_v_3_length", "p_d_5_length", "p_d_3_length", "p_j_5_length",
            "truth_v_call", "truth_d_call", "truth_j_call", "clone_id",
            "parent_id", "duplicate_count", "lineage_node_id", "lineage_parent_id",
            "lineage_generation", "lineage_abundance", "lineage_affinity",
            "lineage_mutations_from_parent",
            "r1_sequence", "r2_sequence",
        }
        collisions = reserved.intersection(fields)
        if collisions:
            raise ValueError(
                f"with_metadata: reserved generated fields cannot be overwritten: {sorted(collisions)}"
            )
        for key, value in fields.items():
            if not isinstance(key, str):
                raise TypeError(
                    f"with_metadata: keys must be strings, got "
                    f"{type(key).__name__}"
                )
            if value is None:
                self._metadata.pop(key, None)
            else:
                self._metadata[key] = value
        return self

    def productive_only(self) -> "Experiment":
        """Require accepted final observed coding spans under the versioned policy.

        ``canonical_observed_vdj_v1`` requires original V/J landmark
        provenance, V Cys and J Phe/Trp residue classes, joining in the
        retained V coding frame, a consistent J-tail frame, and definitely
        no stop in complete observed coding codons. It predicts sequence
        acceptance, not expression, folding, or receptor function. The
        returned ``productivity_policy`` and ``productivity_scope`` identify
        this interpretation. Finite duplex repair explicitly selects
        ``duplex_observed_vdj_v1``, which also accepts certified in-place
        restored landmarks and evaluates their actual current codons. Reference admission overrides do not change it.

        Unknown landmarks or phase and ambiguity that could hide a stop
        cannot satisfy this bundle. Terminal clipping can remain accepted
        when the required landmarks and coding evidence remain observed.
        Every successful emission has ``productive is True`` in strict and
        permissive mode; a final unresolved/failed assessment raises instead
        of emitting a mislabeled record. There is no retry/post-filter loop.

        Sampling narrows each pass's support before drawing. Fixed and mixed
        P lengths and original-end processing use feasible future completion
        support. Draw weights are locally renormalized; this is not a claim
        of a globally productivity-conditioned joint probability model.
        Impossible initial/P-aware support raises even in permissive mode.
        Other passes may use documented no-op behavior if that preserves an
        accepted final record. Strict mode reports their empty supports.
        Replay consumes recorded choices without replacement draws and must
        satisfy the same final assessment.

        This is an order-independent constraint declaration for BCR and TCR
        chains. Exact preservation of the reference anchor residue is a
        separate native ``AnchorPreserved`` contract; the default J class
        permits Phe/Trp changes at the same original landmark.

        Example::

            result = (
                Experiment.on("human_igh")
                .recombine()
                .mutate(count=(5, 15))
                .productive_only()
                .run_records(seed=42)
            )
        """
        if "productive" not in self._contracts:
            self._contracts.append("productive")
        return self
