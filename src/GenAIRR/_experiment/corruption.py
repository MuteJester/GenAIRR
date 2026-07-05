from __future__ import annotations

from typing import Iterable, Optional, Tuple, Union

from .._normalize import _normalize_count
from .._pipeline_ir import (
    _CORRUPT_KIND_3PRIME_LOSS,
    _CORRUPT_KIND_5PRIME_LOSS,
    _CORRUPT_KIND_CONTAMINANT,
    _CORRUPT_KIND_INDEL,
    _CORRUPT_KIND_NS,
    _CORRUPT_KIND_PCR,
    _CORRUPT_KIND_QUALITY,
    _CORRUPT_KIND_REV_COMP,
    _CorruptStep,
)


class _CorruptionMixin:
    __slots__ = ()

    def pcr_amplify(
        self,
        *,
        count: Optional[Union[int, Tuple[int, int], Iterable[Tuple[int, float]]]] = None,
        rate: Optional[float] = None,
    ) -> "Experiment":
        """Append a PCR-error step modelling substitution errors
        introduced during PCR amplification.

        **Specify intensity with exactly one of ``rate`` or ``count``.**

        ``rate`` is the per-base PCR error probability (e.g.
        ``1e-4`` per base per cycle, depending on polymerase
        fidelity). At execute time the engine draws
        ``count ~ Poisson(rate × pool_len)`` against each record's
        current sequence length — matching how PCR error is
        universally reported in the literature. This is the
        canonical biology-default form.

        ``count`` is the legacy explicit count distribution
        (fixed int, ``(low, high)`` range, or empirical
        ``(count, weight)`` list). Useful for benchmark scripts.

        Each error samples a uniform position in the assembled pool
        and replaces the base with a uniform A/C/G/T draw. Trace
        addresses: ``corrupt.pcr.{count, error_site[i], error_base[i]}``.

        Passing both ``count`` and ``rate`` raises ``ValueError``.
        Passing neither raises ``ValueError``.
        """
        return self._append_count_or_rate_step(
            kind=_CORRUPT_KIND_PCR,
            label="pcr_amplify",
            count=count,
            rate=rate,
        )

    def sequencing_errors(
        self,
        *,
        count: Optional[Union[int, Tuple[int, int], Iterable[Tuple[int, float]]]] = None,
        rate: Optional[float] = None,
    ) -> "Experiment":
        """Append a sequencing-quality-error step modelling base-call
        errors during sequencing readout.

        **Specify intensity with exactly one of ``rate`` or ``count``.**

        ``rate`` is the per-base sequencing error probability (e.g.
        ``1e-3`` for a Q30 base, ``1e-2`` for Q20). Drawn as
        ``count ~ Poisson(rate × pool_len)`` per record — the
        canonical Phred-quality framing in immunoseq.

        ``count`` is the legacy explicit count distribution.

        Same shape as :meth:`pcr_amplify` but each substitution
        writes the destination base in **lowercase** to mark the
        position as corrupted (the sequencing-error convention).
        """
        return self._append_count_or_rate_step(
            kind=_CORRUPT_KIND_QUALITY,
            label="sequencing_errors",
            count=count,
            rate=rate,
        )

    def _append_count_or_rate_step(
        self,
        *,
        kind: str,
        label: str,
        count: Optional[Union[int, Tuple[int, int], Iterable[Tuple[int, float]]]],
        rate: Optional[float],
    ) -> "Experiment":
        """Shared validator + step appender for the corruption methods
        that accept both ``count`` and ``rate`` (PCR errors, sequencing
        errors). See :meth:`mutate` for the same shape on the
        mutation side."""
        if count is not None and rate is not None:
            raise ValueError(
                f"{label}(): pass exactly one of `rate` or `count`, not both. "
                f"`rate` is the canonical biology default (e.g. rate=1e-4 "
                f"for per-base PCR error); `count` is the explicit per-record "
                f"count for benchmark / deterministic-count workflows."
            )
        if count is None and rate is None:
            raise ValueError(
                f"{label}(): pass exactly one of `rate` or `count`."
            )
        if rate is not None:
            if not isinstance(rate, (int, float)) or isinstance(rate, bool):
                raise TypeError(
                    f"rate must be a finite float in [0.0, 1.0], got "
                    f"{type(rate).__name__}"
                )
            if not (0.0 <= float(rate) <= 1.0):
                raise ValueError(
                    f"rate must be in [0.0, 1.0] (got {rate!r}); rate is a "
                    f"per-base error probability, not an absolute count."
                )
            self._steps.append(
                _CorruptStep(
                    kind=kind,
                    rate=float(rate),
                )
            )
            return self
        pairs = _normalize_count(count)
        self._steps.append(
            _CorruptStep(
                kind=kind,
                count_pairs=pairs,
            )
        )
        return self

    def polymerase_indels(
        self,
        *,
        count: Union[int, Tuple[int, int], Iterable[Tuple[int, float]]],
        insertion_prob: float = 0.5,
    ) -> "Experiment":
        """Append a polymerase-slippage indel step (PCR-stage artifact).

        Models the small insertions / deletions that arise from
        polymerase slippage during PCR amplification — single-base
        indels at low rate, typically <1 % per read at standard
        polymerase fidelity. Larger structural indels are rare and
        not modelled here.

        ``count`` is the per-simulation distribution over the total
        number of indel events. Each event independently chooses
        insertion (probability ``insertion_prob``) vs. deletion.
        Insertions sample a uniform A/C/G/T base. Indels change
        pool length and shift downstream region ranges accordingly.

        Raises ``ValueError`` if ``insertion_prob`` is outside
        ``[0.0, 1.0]`` or non-finite.
        """
        if (
            insertion_prob != insertion_prob  # NaN check
            or not (0.0 <= insertion_prob <= 1.0)
        ):
            raise ValueError(
                f"insertion_prob must be a finite number in [0.0, 1.0], "
                f"got {insertion_prob}"
            )
        pairs = _normalize_count(count)
        self._steps.append(
            _CorruptStep(
                kind=_CORRUPT_KIND_INDEL,
                count_pairs=pairs,
                insertion_prob=float(insertion_prob),
                apply_prob=0.0,
            )
        )
        return self

    def end_loss_5prime(
        self,
        *,
        length: Union[int, Tuple[int, int], Iterable[Tuple[int, float]]],
    ) -> "Experiment":
        """Append a 5'-end loss corruption step.

        Drops bases from the start of the assembled sequence. The
        underlying engine pass is `EndLossPass` — this models
        **observation-stage** read-end / primer-region loss (the
        AIRR record exposes it as ``end_loss_5_length``, the trace
        address is ``corrupt.end_loss.5``). It is distinct from
        recombination-stage trimming (`v_trim_5` etc.), which writes
        allele-instance metadata rather than deleting pool bytes.

        ``length`` accepts the same shapes as ``count`` on other
        corrupt ops:

        - ``length=10`` — strip exactly 10 bases.
        - ``length=(0, 20)`` — strip a uniform integer in ``[0, 20]``.
        - ``length=[(5, 1.0), (10, 1.0), (15, 1.0)]`` — empirical
          (length, weight) distribution.

        The actual loss is clamped to the pool length (so a sample
        larger than the sequence drops the whole pool). The loss is
        permanent for downstream passes — subsequent corruption
        operates on the shorter pool.

        :meth:`primer_trim_5prime` is a backwards-compatible alias.
        """
        pairs = _normalize_count(length)
        self._steps.append(
            _CorruptStep(
                kind=_CORRUPT_KIND_5PRIME_LOSS,
                count_pairs=pairs,
                insertion_prob=0.0,
                apply_prob=0.0,
            )
        )
        return self

    def end_loss_3prime(
        self,
        *,
        length: Union[int, Tuple[int, int], Iterable[Tuple[int, float]]],
    ) -> "Experiment":
        """Append a 3'-end loss corruption step.

        Drops bases from the end of the assembled sequence. Same
        observation-stage semantics as :meth:`end_loss_5prime`
        (trace address ``corrupt.end_loss.3``, AIRR field
        ``end_loss_3_length``). Same ``length`` shapes.

        :meth:`primer_trim_3prime` is a backwards-compatible alias.
        """
        pairs = _normalize_count(length)
        self._steps.append(
            _CorruptStep(
                kind=_CORRUPT_KIND_3PRIME_LOSS,
                count_pairs=pairs,
                insertion_prob=0.0,
                apply_prob=0.0,
            )
        )
        return self

    def primer_trim_5prime(
        self,
        *,
        length: Union[int, Tuple[int, int], Iterable[Tuple[int, float]]],
    ) -> "Experiment":
        """Backwards-compatible alias for :meth:`end_loss_5prime`.

        The DSL name `primer_trim_5prime` predates the engine's
        cleaner vocabulary — the same step is now exposed as
        :meth:`end_loss_5prime`. Behaviour is identical (same trace
        address ``corrupt.end_loss.5``, same AIRR field
        ``end_loss_5_length``, byte-identical records for the same
        seed). New code should prefer the `end_loss_*` form; the
        alias stays so existing scripts keep working.
        """
        return self.end_loss_5prime(length=length)

    def primer_trim_3prime(
        self,
        *,
        length: Union[int, Tuple[int, int], Iterable[Tuple[int, float]]],
    ) -> "Experiment":
        """Backwards-compatible alias for :meth:`end_loss_3prime`.

        See :meth:`primer_trim_5prime` for the rationale.
        """
        return self.end_loss_3prime(length=length)

    def ambiguous_base_calls(
        self,
        *,
        count: Union[int, Tuple[int, int], Iterable[Tuple[int, float]]],
    ) -> "Experiment":
        """Append an N-substitution corruption step.

        With the per-simulation count drawn from ``count``, replace
        that many uniform-random pool positions with the ambiguous
        base ``N``. Models the low-quality positions that real
        sequencers emit when the base caller cannot commit.

        ``count`` accepts the same shapes as other corruption ops:
        - ``count=10`` — fixed.
        - ``count=(0, 20)`` — uniform integer range.
        - ``count=[(5, 1.0), (10, 1.0)]`` — empirical (count, weight).

        Sampled sites are with-replacement, so ``count`` is the
        upper bound on resulting `N` count (collisions reduce it).
        """
        pairs = _normalize_count(count)
        self._steps.append(
            _CorruptStep(
                kind=_CORRUPT_KIND_NS,
                count_pairs=pairs,
                insertion_prob=0.0,
                apply_prob=0.0,
            )
        )
        return self

    def random_strand_orientation(self, *, prob: float = 0.5) -> "Experiment":
        """Append a reverse-complement corruption step.

        With probability ``prob`` the AIRR record-builder reverse-
        complements the ``sequence``, ``np1``, ``np2``, and
        ``junction`` fields and flips the corresponding pool-position
        coords (``*_sequence_start/end``, ``junction_start/end``).
        Alignment / germline / CIGAR / identity fields stay in
        forward orientation per the AIRR Rearrangement spec.
        ``prob=0.0`` is a no-op (coin flip recorded but never fires);
        ``prob=1.0`` always flips. The biological model is the ~50%
        of antisense reads in real immune-seq libraries.

        Raises ``ValueError`` if ``prob`` is outside ``[0.0, 1.0]``
        or non-finite.
        """
        if prob != prob or not (0.0 <= prob <= 1.0):
            raise ValueError(
                f"prob must be a finite number in [0.0, 1.0], got {prob}"
            )
        self._steps.append(
            _CorruptStep(
                kind=_CORRUPT_KIND_REV_COMP,
                count_pairs=(),
                insertion_prob=0.0,
                apply_prob=float(prob),
            )
        )
        return self

    def contaminate(self, *, prob: float) -> "Experiment":
        """Append a contaminant-replacement corruption step.

        With probability ``prob`` the entire assembled pool is
        overwritten with uniform A/C/G/T bases. Models primer
        dimers, bacterial DNA, or any non-receptor sequence in the
        library. ``prob=0.0`` is a no-op (coin flip recorded but
        never fires); ``prob=1.0`` always contaminates.

        Raises ``ValueError`` if ``prob`` is outside ``[0.0, 1.0]``
        or non-finite.
        """
        if prob != prob or not (0.0 <= prob <= 1.0):
            raise ValueError(
                f"prob must be a finite number in [0.0, 1.0], got {prob}"
            )
        self._steps.append(
            _CorruptStep(
                kind=_CORRUPT_KIND_CONTAMINANT,
                count_pairs=(),
                insertion_prob=0.0,
                apply_prob=float(prob),
            )
        )
        return self
