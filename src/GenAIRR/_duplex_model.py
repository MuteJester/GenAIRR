"""Immutable, reference-bound authored repair choices; no inferred enzyme rates."""
from __future__ import annotations

import json
import math
from dataclasses import dataclass
from typing import Tuple

from . import _engine

MODEL = "finite_authored_duplex_history_v1"
CONDITIONING_LAW = "finite_authored_history_local_conditioning_v1"


def _work_budget(value, name):
    if type(value) is not int or not 0 < value < 2**64:
        raise ValueError(f"{name} must be a positive unsigned-64-bit integer")
    return value


@dataclass(frozen=True)
class DuplexHistory:
    """One explicitly named history and its weight within its source context.

    ``script`` supplies complete canonical V/D/J substrates, orientation and
    ordered physical operations. Calls identify the actual generating alleles.
    Different histories with identical final DNA retain their separate weights.
    """

    id: str
    script: object
    v_call: str
    d_call: str
    j_call: str
    weight: float = 1.0

    def __post_init__(self):
        from .experimental import CoupledDuplexRepair

        if not isinstance(self.id, str) or not self.id:
            raise ValueError("history id must be a nonempty string")
        if not isinstance(self.script, CoupledDuplexRepair):
            raise TypeError("history script must be a CoupledDuplexRepair")
        for name in ("v_call", "d_call", "j_call"):
            if not isinstance(getattr(self, name), str) or not getattr(self, name):
                raise ValueError(f"{name} must be a nonempty generating allele name")
        if isinstance(self.weight, bool):
            raise ValueError("history weight must be a finite nonnegative number")
        try:
            weight = float(self.weight)
        except (TypeError, ValueError, OverflowError):
            raise ValueError("history weight must be a finite nonnegative number") from None
        if not math.isfinite(weight) or weight < 0:
            raise ValueError("history weight must be a finite nonnegative number")
        object.__setattr__(self, "weight", weight if weight else 0.0)

    def to_dict(self):
        return {"id": self.id, "weight": self.weight, "v_call": self.v_call,
                "d_call": self.d_call, "j_call": self.j_call,
                "script": self.script.to_dict()}

    @classmethod
    def from_dict(cls, value):
        from .experimental import CoupledDuplexRepair

        fields = {"id", "weight", "v_call", "d_call", "j_call", "script"}
        if not isinstance(value, dict) or set(value) != fields:
            raise ValueError("unsupported duplex history row schema")
        return cls(value["id"], CoupledDuplexRepair.from_dict(value["script"]),
                   value["v_call"], value["d_call"], value["j_call"], value["weight"])


@dataclass(frozen=True)
class DuplexRepairModel:
    """A finite authored model pinned to exact reference content/evidence.

    Create with :meth:`bind`. Compilation validates all reachable positive
    histories before applying productivity conditioning. These weights describe
    choices among authored histories, not measured biochemical rates.
    ``max_preflight_operations`` bounds reachable-context checks, authored
    operations and inspected source weight entries, including distinct
    productivity-filtered flat supports and natural genotype masses.
    ``max_preflight_bases`` bounds the sum of
    active scripts' base budgets reserved before compilation; removal refunds
    no physical allocation.
    """

    reference_content_hash: str
    rows: Tuple[DuplexHistory, ...]
    reference_ends_are_coding_ends: bool = False
    max_preflight_operations: int = 1_000_000
    max_preflight_bases: int = 10_000_000

    def __post_init__(self):
        if not isinstance(self.reference_content_hash, str) or not self.reference_content_hash:
            raise ValueError("reference_content_hash must be a nonempty content identity")
        if type(self.reference_ends_are_coding_ends) is not bool:
            raise ValueError("reference_ends_are_coding_ends must be a boolean")
        try:
            rows = tuple(self.rows)
        except TypeError:
            raise TypeError("rows must be an iterable of DuplexHistory values") from None
        if not rows or any(not isinstance(row, DuplexHistory) for row in rows):
            raise ValueError("rows must contain at least one DuplexHistory")
        if len(rows) > 100_000:
            raise ValueError("a duplex model cannot exceed 100000 history rows")
        if len({row.id for row in rows}) != len(rows):
            raise ValueError("duplex history ids must be unique")
        if not any(row.weight > 0 for row in rows):
            raise ValueError("a duplex model needs positive history support")
        object.__setattr__(self, "rows", rows)
        _work_budget(self.max_preflight_operations, "max_preflight_operations")
        _work_budget(self.max_preflight_bases, "max_preflight_bases")

    @classmethod
    def bind(cls, refdata, histories, *, reference_ends_are_coding_ends=False,
             max_preflight_operations=1_000_000, max_preflight_bases=10_000_000):
        """Pin histories to a native reference; compile performs physical checks.

        Unknown supplied coding boundaries require an explicit assertion.
        Known partial junction-facing boundaries remain unsupported. All
        selectable allele/orientation contexts must have positive support;
        use explicit experiment allele restrictions to narrow that domain.
        """
        if not isinstance(refdata, _engine.RefDataConfig):
            raise TypeError("bind requires a GenAIRR RefDataConfig")
        return cls(refdata.content_hash(), tuple(histories), reference_ends_are_coding_ends,
                   max_preflight_operations, max_preflight_bases)

    def to_dict(self):
        return {"model": MODEL, "reference_content_hash": self.reference_content_hash,
                "reference_ends_are_coding_ends": self.reference_ends_are_coding_ends,
                "max_preflight_operations": self.max_preflight_operations,
                "max_preflight_bases": self.max_preflight_bases,
                "rows": [row.to_dict() for row in self.rows]}

    def to_json(self):
        """Canonical native-serializable choices, with authored order retained."""
        return json.dumps(self.to_dict(), sort_keys=True, allow_nan=False, separators=(",", ":"))

    @classmethod
    def from_dict(cls, value):
        fields = {"model", "reference_content_hash", "reference_ends_are_coding_ends",
                  "max_preflight_operations", "max_preflight_bases", "rows"}
        if not isinstance(value, dict) or set(value) != fields or value["model"] != MODEL:
            raise ValueError("unsupported finite duplex model schema")
        if not isinstance(value["rows"], list):
            raise ValueError("model rows must be a JSON array")
        return cls(value["reference_content_hash"], tuple(DuplexHistory.from_dict(row) for row in value["rows"]),
                   value["reference_ends_are_coding_ends"], value["max_preflight_operations"], value["max_preflight_bases"])
