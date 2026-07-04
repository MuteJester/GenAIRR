from __future__ import annotations

from typing import Any


class _ConstraintsMixin:
    __slots__ = ()

    def with_metadata(self, **fields: Any) -> "Experiment":
        """Attach sample-level metadata to every AIRR record.

        Standard AIRR Repertoire fields like ``sample_id``,
        ``donor``, ``repertoire_id``, and ``cell_id`` are commonly
        used by multi-sample analysis tools — Change-O, scirpy,
        immcantation pipelines all expect them populated. Custom
        keys are also accepted and pass through unchanged.

        Subsequent ``with_metadata()`` calls **merge** with the
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
        """Require every emitted record to be a productive sequence.

        Attaches the canonical productive-sequence contract bundle to
        the experiment: junction frame in-register, no stop codons in
        the junction, and V/J anchor amino acids preserved. The bundle
        is enforced during recombination, mutation, and corruption
        passes by narrowing each pass's action support before
        sampling. If a constrained support is empty, permissive
        execution uses the pass's explicit no-op / sentinel behavior;
        use ``strict=True`` at run time to raise instead.

        **Failure surfaces** (see
        ``docs/productive_failure_mode_audit.md`` for the full matrix):

        - *Compile-time precondition* — when a sampling distribution
          is statically impossible under the bundle (e.g. every NP1
          length violates frame), :meth:`compile` raises
          ``ValueError`` regardless of the ``strict`` flag.
        - *Runtime fresh strict* (``run(..., strict=True)``) — when
          dynamic state makes a sampler's admissible support empty,
          raises :class:`GenAIRR._engine.StrictSamplingError` with
          structured ``(pass_name, address, reason)`` args.
        - *Runtime fresh permissive* (``strict=False``, default) —
          the pass records its declared sentinel (indel ``site=-1``,
          NP length ``0``, NP base ``N``, trim ``0``) or skips the
          slot; the record continues.
        - *Trace replay* (``replay_from_trace_file``) — consumes
          recorded values verbatim, does not re-evaluate
          admissibility. A permissive-sentinel trace replays cleanly
          even with ``strict=True``.

        **Order-independent.** This method is a constraint declaration,
        not a pipeline step — it can be called anywhere in the chain
        and the result is identical. Convention is to place it last
        (right before ``run_records()``) so the constraint reads as a
        post-hoc requirement on the emitted records.

        TCR refdata accepts the call but raises ``ValueError`` at
        :meth:`compile` time because TCRs don't have somatic
        hypermutation and the productive bundle's anchor checks
        assume BCR semantics. Catch this early at the builder if it
        matters to you.

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
