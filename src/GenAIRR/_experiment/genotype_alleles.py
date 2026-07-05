from __future__ import annotations

from typing import Dict, Iterable, List, Optional, Tuple, Union

from ._common import _UNSET


class _GenotypeAllelesMixin:
    __slots__ = ()

    def with_genotype(self, genotype) -> "Experiment":
        """Attach a single-subject diploid genotype.

        With a genotype attached, V(D)J recombination becomes
        haplotype-phased: V, D and J of each rearrangement are drawn from
        a single chromosome, honouring the genotype's allele
        presence/absence, zygosity, and copy-number/deletion. With no
        genotype, the flat (uniform / usage-weighted) path runs unchanged.

        Mutually exclusive with :meth:`restrict_alleles` and the
        ``recombine(*_allele_weights=...)`` kwargs — the genotype owns
        allele presence and within-gene expression.

        Raises ``ValueError`` if the genotype was built against a
        different cartridge (content-hash mismatch), or if allele locks /
        explicit allele weights are already set.
        """
        live_hash = self._refdata.content_hash()
        if genotype._source_hash != live_hash:
            raise ValueError(
                "genotype was built against a different cartridge (content hash "
                f"{genotype._source_hash!r} != experiment {live_hash!r})"
            )
        if any(v is not None for v in self._locks.values()):
            raise ValueError(
                "with_genotype() and restrict_alleles() are mutually exclusive"
            )
        if self._user_allele_weights_set:
            raise ValueError(
                "with_genotype() and recombine(*_allele_weights=...) are mutually "
                "exclusive: the genotype owns allele expression"
            )
        # Snapshot the (mutable) builder so later edits to ``genotype``
        # cannot desync the compiled engine genotype from
        # ``result.genotypes`` (review #8).
        self._genotype = genotype._snapshot()
        return self

    def restrict_alleles(
        self,
        *,
        v: "_LockInput" = _UNSET,
        d: "_LockInput" = _UNSET,
        j: "_LockInput" = _UNSET,
    ) -> "Experiment":
        """Narrow allele sampling to a named subset, per segment.

        Sampling stays uniform — this restricts the *support* of the
        sampling distribution to the listed allele names, not pins a
        single allele. If you supply one name, the result is
        effectively deterministic (uniform over one element); for any
        list of N names, the recombination pass samples uniformly
        over those N alleles.

        Each kwarg accepts:
        - a single allele name string (e.g. ``"IGHV1-2*02"``),
        - a list / tuple of allele names (sample uniformly among them),
        - ``None`` to clear a previously-set restriction for that segment,
        - omitted (default) — the existing restriction for that segment
          is left unchanged.

        The restrictions are applied to the next ``recombine()`` step's
        ``push_sample_allele`` calls at compile time. Calling
        ``restrict_alleles()`` more than once overlays the new values
        onto the previous restrictions (per-segment), so
        ``.restrict_alleles(v="A").restrict_alleles(d="B")`` restricts
        both V and D.

        Raises:
        - ``ValueError`` if an allele name doesn't exist in the
          configured refdata, or if a restriction is set for ``D`` on
          a VJ chain.
        - ``TypeError`` if an unexpected input shape is passed.
        """
        if self._genotype is not None:
            raise ValueError(
                "restrict_alleles() and with_genotype() are mutually exclusive: "
                "a genotype already owns allele presence and within-gene expression"
            )
        for segment, value in (("V", v), ("D", d), ("J", j)):
            if value is _UNSET:
                continue
            if value is None:
                self._locks[segment] = None
                continue
            ids = self._resolve_lock(segment, value)
            self._locks[segment] = ids
        return self

    def _resolve_lock(
        self,
        segment: str,
        value: Union[str, Iterable[str]],
    ) -> Tuple[int, ...]:
        """Resolve allele-name(s) → tuple of allele IDs against this
        experiment's refdata. Raises on unknown names, on D locks
        for VJ chains, or on non-string inputs."""
        if segment == "D" and self._refdata.chain_type != "vdj":
            raise ValueError(
                f"cannot lock D alleles on a {self._refdata.chain_type!r} chain"
            )

        if isinstance(value, str):
            names: Tuple[str, ...] = (value,)
        else:
            try:
                names = tuple(value)
            except TypeError as exc:
                raise TypeError(
                    f"restrict_alleles(): {segment} lock must be a name string or an "
                    f"iterable of name strings; got {type(value).__name__}"
                ) from exc
            if not all(isinstance(n, str) for n in names):
                raise TypeError(
                    f"restrict_alleles(): {segment} lock entries must all be strings"
                )
            if not names:
                raise ValueError(
                    f"restrict_alleles(): {segment} lock list must be non-empty "
                    "(pass None to clear instead)"
                )

        index = self._allele_name_index(segment)
        ids: List[int] = []
        seen: set = set()
        for name in names:
            if name not in index:
                raise ValueError(
                    f"restrict_alleles(): no {segment} allele named {name!r} in refdata "
                    "(check spelling against `list_alleles` or the loaded config)"
                )
            allele_id = index[name]
            if allele_id in seen:
                raise ValueError(
                    f"restrict_alleles(): duplicate {segment} allele {name!r} in lock list"
                )
            seen.add(allele_id)
            ids.append(allele_id)
        return tuple(ids)

    def _allele_name_index(self, segment: str) -> Dict[str, int]:
        """Build a name → allele_id map for the given segment by
        scanning the refdata pool. Cheap enough to do per-call given
        typical pool sizes (≤ a few hundred) and the rarity of
        ``.restrict_alleles()``."""
        if segment == "V":
            n = self._refdata.v_pool_size()
            getter = self._refdata.v_allele
        elif segment == "D":
            n = self._refdata.d_pool_size()
            getter = self._refdata.d_allele
        elif segment == "J":
            n = self._refdata.j_pool_size()
            getter = self._refdata.j_allele
        else:  # pragma: no cover — guarded above.
            raise ValueError(f"unsupported segment {segment!r}")
        return {getter(i).name: i for i in range(n)}

    def _resolve_allele_weights(
        self,
        segment: str,
        user_weights: Optional[Dict[str, float]],
    ) -> Optional[Tuple[float, ...]]:
        """Build a dense pool-aligned weight vector from the
        user-supplied ``{name: weight}`` dict. Listed alleles get the
        supplied weight; everything else gets ``1.0``. Returns
        ``None`` when no weights were supplied (preserving the
        upstream uniform-default behavior).

        Raises ``ValueError`` for unknown allele names, non-positive
        weights, or D weights on a VJ chain.
        """
        if user_weights is None:
            return None
        if not user_weights:
            raise ValueError(
                f"recombine: {segment.lower()}_allele_weights must contain "
                "at least one (name, weight) entry"
            )
        if segment == "D" and self._refdata.chain_type != "vdj":
            raise ValueError(
                f"recombine: cannot weight D alleles on a "
                f"{self._refdata.chain_type!r} chain"
            )

        index = self._allele_name_index(segment)
        if not index:
            return None  # No alleles for this segment (e.g. VJ + D).

        # Dense vector indexed by allele_id; default 1.0.
        pool_size = max(index.values()) + 1
        weights: List[float] = [1.0] * pool_size
        for name, w in user_weights.items():
            if not isinstance(name, str):
                raise TypeError(
                    f"recombine: {segment.lower()}_allele_weights keys must "
                    f"be allele-name strings, got {type(name).__name__}"
                )
            if isinstance(w, bool) or not isinstance(w, (int, float)):
                raise TypeError(
                    f"recombine: {segment.lower()}_allele_weights['{name}'] "
                    f"must be numeric, got {type(w).__name__}"
                )
            wf = float(w)
            if not (wf > 0.0 and wf == wf and wf != float("inf")):
                raise ValueError(
                    f"recombine: {segment.lower()}_allele_weights['{name}'] "
                    f"must be a finite positive number, got {wf}"
                )
            if name not in index:
                raise ValueError(
                    f"recombine: no {segment} allele named {name!r} in "
                    "refdata (check spelling against `list_alleles` or the "
                    "loaded config)"
                )
            weights[index[name]] = wf
        return tuple(weights)
