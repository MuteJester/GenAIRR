"""Input-normalization helpers for :class:`GenAIRR.experiment.Experiment`.

The fluent DSL accepts user-friendly shapes for distributions and
counts (a single int, a ``(low, high)`` range tuple, or a list of
``(value, weight)`` pairs); the underlying typed IR
(:mod:`._pipeline_ir`) expects the canonical ``tuple of (int, float)``
shape. These functions are the boundary.

This module is intentionally pure: no engine handles, no
dataconfig parsing. It depends only on the typed-IR module for
the NP-length default.
"""
from __future__ import annotations

from typing import Iterable, List, Optional, Tuple, Union

from ._pipeline_ir import _DEFAULT_NP_LENGTHS


def _normalize_count(
    count: Union[int, Tuple[int, int], Iterable[Tuple[int, float]]],
) -> Tuple[Tuple[int, float], ...]:
    """Normalise the user-friendly ``count=`` argument to the
    list-of-pairs shape ``PyPassPlan.push_mutate_*`` expects.

    Accepts:
    - ``count=15`` — fixed count.
    - ``count=(5, 25)`` — uniform integer in ``[5, 25]``
      (high inclusive, matches the Pythonic ``range``-with-stop
      conventions; both endpoints are valid sample values).
    - ``count=[(5, 1.0), (10, 1.0), (15, 1.0)]`` — explicit
      empirical distribution.
    """
    # Single int: fixed count. Reject bools (which are int subclasses)
    # so a stray True/False raises TypeError instead of silently
    # becoming ``count=1`` / ``count=0``.
    if isinstance(count, bool) or not isinstance(count, (int, tuple, list)):
        raise TypeError(
            f"count: expected int, (low, high) tuple, or list of "
            f"(count, weight) pairs, got {type(count).__name__}"
        )
    if isinstance(count, int):
        if count < 0:
            raise ValueError(f"count must be non-negative, got {count}")
        return ((int(count), 1.0),)
    # Tuple form: must be (low, high) of two ints OR a single
    # (count, weight) pair we promote to a 1-element list.
    if isinstance(count, tuple) and len(count) == 2:
        a, b = count[0], count[1]
        if (
            isinstance(a, int)
            and isinstance(b, int)
            and not isinstance(a, bool)
            and not isinstance(b, bool)
        ):
            low, high = a, b
            if low < 0 or high < low:
                raise ValueError(
                    f"count range must satisfy 0 <= low <= high, got ({low}, {high})"
                )
            return tuple((c, 1.0) for c in range(low, high + 1))
    # Otherwise: empirical (count, weight) list/tuple.
    pairs: List[Tuple[int, float]] = []
    for pair in count:
        if not (isinstance(pair, (tuple, list)) and len(pair) == 2):
            raise TypeError(
                f"count entries must be (count, weight) pairs, got {pair!r}"
            )
        c, w = pair
        if isinstance(c, bool) or not isinstance(c, int) or c < 0:
            raise ValueError(
                f"count entry must have a non-negative int count, got {c!r}"
            )
        pairs.append((int(c), float(w)))
    if not pairs:
        raise ValueError("count distribution must contain at least one entry")
    return tuple(pairs)


def _to_immutable_pairs(
    pairs: Optional[List[Tuple[int, float]]],
) -> Optional[Tuple[Tuple[int, float], ...]]:
    """Convert an optional ``list[(int, float)]`` (the shape returned
    by :mod:`._dataconfig_extract`) into the hashable tuple-of-tuples
    form that :class:`_RecombineStep`'s frozen dataclass requires.
    ``None`` passes through unchanged so the step can omit the trim
    pass entirely when no empirical distribution is available.
    """
    if not pairs:
        return None
    return tuple((int(c), float(w)) for c, w in pairs)


def _normalize_lengths(
    lengths: Optional[Iterable[Tuple[int, float]]],
) -> Tuple[Tuple[int, float], ...]:
    """Convert a user-supplied length iterable into a hashable tuple of
    ``(int, float)`` pairs. ``None`` falls back to the module default;
    empty iterables raise ``ValueError`` here so the error surfaces at
    builder time."""
    if lengths is None:
        return tuple(_DEFAULT_NP_LENGTHS)
    pairs: List[Tuple[int, float]] = []
    for pair in lengths:
        length, weight = pair
        pairs.append((int(length), float(weight)))
    if not pairs:
        raise ValueError(
            "length distribution must contain at least one (length, weight) pair"
        )
    return tuple(pairs)
