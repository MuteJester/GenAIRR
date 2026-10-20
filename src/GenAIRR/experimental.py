"""Experimental, explicitly authored duplex repair inside GenAIRR.

The local kernel resolves a pair of complete coding-end hairpins. It is
independent of the AIRR experiment pipeline until the richer ancestry can be
projected without losing joint strand origins. Histories carry no fitted rates.
"""
from __future__ import annotations

import json
from dataclasses import dataclass

from . import _engine


def _nonnegative(value, name):
    if type(value) is not int or value < 0 or value >= 2**63:
        raise ValueError(f"{name} must be a nonnegative signed-64-bit integer")
    return value


def _side(value):
    if value not in ("left", "right"):
        raise ValueError("side must be 'left' or 'right'")
    return value


def _strand(value):
    if value not in ("top", "bottom"):
        raise ValueError("strand must be 'top' or 'bottom'")
    return value


@dataclass(frozen=True)
class DuplexRepair:
    """Immutable script for a bounded two-ended duplex repair experiment.

    Top strands follow the coding axis 5'->3'; bottom strands follow that axis
    3'->5'. ``add_n`` sequences always follow polymerization order 5'->3'.
    ``right_start`` is the coordinate of the right original coding origin on
    the left molecule's axis. Only fully paired, doubly ligated products run.
    """

    left: str
    right: str
    operations_json: str = "[]"
    max_operations: int = 256
    max_bases: int = 100_000

    def __post_init__(self):
        for name in ("left", "right"):
            seq = getattr(self, name)
            if not isinstance(seq, str) or not seq or any(b not in "ACGTacgt" for b in seq):
                raise ValueError("coding ends require nonempty resolved A/C/G/T strings")
            object.__setattr__(self, name, seq.upper())
        _nonnegative(self.max_operations, "max_operations")
        _nonnegative(self.max_bases, "max_bases")
        if self.max_bases > 1_000_000:
            raise ValueError("max_bases cannot exceed 1000000")
        try:
            operations = json.loads(self.operations_json)
        except (TypeError, ValueError):
            raise ValueError("operations_json must contain a JSON array") from None
        if not isinstance(operations, list) or any(not isinstance(x, dict) for x in operations):
            raise ValueError("operations_json must contain an array of operations")
        object.__setattr__(self, "operations_json", json.dumps(operations, sort_keys=True, allow_nan=False))

    @property
    def operations(self):
        return json.loads(self.operations_json)

    def _append(self, operation):
        from dataclasses import replace
        ops = self.operations + [operation]
        if len(ops) > self.max_operations:
            raise ValueError("repair operation budget exceeded")
        return replace(self, operations_json=json.dumps(ops, allow_nan=False))

    def open(self, side, *, offset):
        """Nick an original hairpin; offset produces a 2*offset 3' overhang."""
        return self._append({"operation": "open", "side": _side(side),
                             "offset": _nonnegative(offset, "offset")})

    def resect(self, side, strand, *, count):
        return self._append({"operation": "resect", "side": _side(side),
                             "strand": _strand(strand), "count": _nonnegative(count, "count")})

    def add_n(self, side, sequence):
        return self._append({"operation": "add_n", "side": _side(side), "sequence": sequence})

    def synapse(self, *, right_start):
        if type(right_start) is not int or not -(2**63) <= right_start < 2**63:
            raise ValueError("right_start must be a signed-64-bit integer")
        return self._append({"operation": "synapse", "right_start": right_start})

    def fill(self, strand):
        return self._append({"operation": "fill", "strand": _strand(strand)})

    def ligate(self):
        return self._append({"operation": "ligate"})

    def run(self):
        """Return the resolved sequence, both strand ancestries and full history."""
        return json.loads(_engine.duplex_repair(self.left, self.right, self.operations_json,
                                                self.max_operations, self.max_bases))

    def to_dict(self):
        return {"model": "experimental_duplex_repair_v1", "left": self.left,
                "right": self.right, "operations": self.operations,
                "max_operations": self.max_operations, "max_bases": self.max_bases}

    @classmethod
    def from_dict(cls, value):
        if not isinstance(value, dict) or set(value) != {"model", "left", "right", "operations", "max_operations", "max_bases"} or value["model"] != "experimental_duplex_repair_v1":
            raise ValueError("unsupported duplex repair script schema")
        return cls(value["left"], value["right"], json.dumps(value["operations"], allow_nan=False),
                   value["max_operations"], value["max_bases"])
