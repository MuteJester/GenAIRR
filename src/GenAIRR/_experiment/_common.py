from __future__ import annotations

from typing import Iterable, Union


# Sentinel for `using(...)` arguments that the caller did not pass at
# all. We can't use ``None`` for "unchanged" because ``None`` already
# means "clear the lock."
class _Unset:
    __slots__ = ()

    def __repr__(self) -> str:
        return "<UNSET>"


_UNSET: _Unset = _Unset()


# Inputs accepted by :meth:`Experiment.restrict_alleles` per segment. ``None``
# clears any prior lock; a single name is a one-allele lock; an
# iterable of names is a multi-allele subset; ``_UNSET`` (the default)
# leaves the existing lock unchanged.
_LockInput = Union[str, Iterable[str], None, _Unset]
