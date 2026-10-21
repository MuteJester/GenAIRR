"""Experimental, explicitly authored duplex repair inside GenAIRR.

The kernels resolve authored local joins and coupled successive V-D/D-J joins.
Standalone reference-bound projection preserves both strand origins and exposes
nullable AIRR metadata. Finite reference-bound histories can also be selected
by Experiment.recombine(duplex_repair=model). Histories carry no fitted rates
and start from supplied coding-end hairpins.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, replace

from . import _engine

from ._duplex_projection_export import export_duplex_projection
from ._duplex_model import DuplexHistory, DuplexRepairModel
from ._family_replay import FamilyReplayBundle, FamilyReplayError, FamilyReplayLimits

__all__ = ["DuplexHistory", "DuplexRepairModel", "DuplexRepair", "DuplexJoin", "CoupledDuplexRepair", "export_duplex_projection", "project_duplex_simulation",
           "FamilyReplayBundle", "FamilyReplayError", "FamilyReplayLimits"]


def project_duplex_simulation(simulation, refdata):
    """Project a duplex-backed persistent simulation after audited edits.

    The founder duplex is preserved separately from the observed query.
    This projector returns the full nullable envelope. Normal Experiment
    records use the same nullable metadata; ``export_duplex_projection`` can
    also export this standalone envelope directly.
    """
    if not isinstance(simulation, _engine.Simulation):
        raise TypeError("project_duplex_simulation requires a GenAIRR Simulation")
    if not isinstance(refdata, _engine.RefDataConfig):
        raise TypeError("project_duplex_simulation requires a GenAIRR RefDataConfig")
    return json.loads(_engine.project_duplex_simulation(simulation, refdata))


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


def _validate_dna_fields(script, names):
    for name in names:
        sequence = getattr(script, name)
        if not isinstance(sequence, str) or not sequence or any(b not in "ACGTacgt" for b in sequence):
            raise ValueError("coding ends require nonempty resolved A/C/G/T strings")
        object.__setattr__(script, name, sequence.upper())


def _validate_work_limits(script):
    _nonnegative(script.max_operations, "max_operations")
    _nonnegative(script.max_bases, "max_bases")
    if script.max_bases > 1_000_000:
        raise ValueError("max_bases cannot exceed 1000000")


def _validate_operations(script):
    try:
        operations = json.loads(script.operations_json)
    except (TypeError, ValueError):
        raise ValueError("operations_json must contain a JSON array") from None
    if not isinstance(operations, list) or any(not isinstance(x, dict) for x in operations):
        raise ValueError("operations_json must contain an array of operations")
    object.__setattr__(script, "operations_json", json.dumps(operations, sort_keys=True, allow_nan=False))


class _RepairScriptOperations:
    @property
    def operations(self):
        return json.loads(self.operations_json)

    def _append(self, operation):
        ops = self.operations + [operation]
        limit = getattr(self, "max_operations", None)
        if limit is not None and len(ops) > limit:
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


@dataclass(frozen=True)
class DuplexRepair(_RepairScriptOperations):
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
        _validate_dna_fields(self, ("left", "right"))
        _validate_work_limits(self)
        _validate_operations(self)

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


@dataclass(frozen=True)
class DuplexJoin(_RepairScriptOperations):
    """Authored operations for one ``vd`` or ``dj`` join, without new substrates.

    ``left``/``right`` refer to the exposed ends at that join. At ``vd`` these
    are V's right tip and oriented D's left tip; at ``dj`` they are oriented D's
    right tip and J's left tip. ``right_start`` locates the original oriented D
    origin relative to V for ``vd``, or J's origin relative to oriented D for
    ``dj``. Resection never changes those origin anchors.
    """

    name: str
    operations_json: str = "[]"

    def __post_init__(self):
        if self.name not in ("vd", "dj"):
            raise ValueError("join name must be 'vd' or 'dj'")
        _validate_operations(self)

    def to_dict(self):
        return {"join": self.name, "operations": self.operations}

    @classmethod
    def from_dict(cls, value):
        if not isinstance(value, dict) or set(value) != {"join", "operations"}:
            raise ValueError("unsupported duplex join script schema")
        return cls(value["join"], json.dumps(value["operations"], allow_nan=False))


@dataclass(frozen=True)
class CoupledDuplexRepair:
    """Two successive joins on one shared V/D/J two-strand substrate state.

    All references are supplied in their canonical orientation. ``d_inverted``
    reverses D's placement, retaining each original nucleotide's canonical
    strand and coordinate. The two joins may be authored in either order.
    Only a complete, paired, connected VDJ product is returned.

    Budgets apply cumulatively to both joins. Initial bases on both strands
    and every synthesized base count toward ``max_bases``; resection refunds
    none. These authored histories carry no calibrated probabilities.
    """

    v: str
    d: str
    j: str
    d_inverted: bool = False
    joins_json: str = "[]"
    max_operations: int = 256
    max_bases: int = 100_000

    def __post_init__(self):
        _validate_dna_fields(self, ("v", "d", "j"))
        _validate_work_limits(self)
        if type(self.d_inverted) is not bool:
            raise ValueError("d_inverted must be a boolean")
        try:
            joins = json.loads(self.joins_json)
        except (TypeError, ValueError):
            raise ValueError("joins_json must contain a JSON array") from None
        if not isinstance(joins, list):
            raise ValueError("joins_json must contain a JSON array")
        scripts = [DuplexJoin.from_dict(value) for value in joins]
        if len(scripts) > 2 or len({script.name for script in scripts}) != len(scripts):
            raise ValueError("provide each of the vd and dj joins at most once")
        if sum(len(script.operations) for script in scripts) > self.max_operations:
            raise ValueError("repair operation budget exceeded across both joins")
        object.__setattr__(self, "joins_json", json.dumps([s.to_dict() for s in scripts],
                                                        sort_keys=True, allow_nan=False))

    @property
    def joins(self):
        """Independent JSON copy of the authored joins in execution order."""
        return json.loads(self.joins_json)

    def join(self, script):
        if not isinstance(script, DuplexJoin):
            raise ValueError("join() requires a DuplexJoin script")
        return replace(self, joins_json=json.dumps(self.joins + [script.to_dict()], allow_nan=False))

    def run(self):
        """Resolve both joins while retaining the global nucleotide ledger."""
        return json.loads(_engine.coupled_duplex_repair(
            self.v, self.d, self.j, self.joins_json, d_inverted=self.d_inverted,
            max_operations=self.max_operations, max_bases=self.max_bases))

    def project(self, refdata, *, v_call, d_call, j_call,
                reference_ends_are_coding_ends=False):
        """Return a bound AIRR record, certified paths and paired ancestry.

        Calls select exact generating alleles from ``RefDataConfig``; the
        substrates must match their whole canonical sequences. Unknown coding
        boundaries require the explicit assertion. Known partial junction ends
        cannot be overridden. Valid overlapping paths yield a partial record
        and capability reasons, with no ancestry or history silently removed.
        """
        if not isinstance(refdata, _engine.RefDataConfig):
            raise TypeError("project() requires a GenAIRR RefDataConfig")
        if type(reference_ends_are_coding_ends) is not bool:
            raise ValueError("reference_ends_are_coding_ends must be a boolean")
        for name, value in (("v_call", v_call), ("d_call", d_call), ("j_call", j_call)):
            if not isinstance(value, str) or not value:
                raise ValueError(f"{name} must be a nonempty generating allele name")
        return json.loads(_engine.coupled_duplex_projection(
            self.v, self.d, self.j, self.joins_json, refdata,
            v_call=v_call, d_call=d_call, j_call=j_call,
            d_inverted=self.d_inverted,
            reference_ends_are_coding_ends=reference_ends_are_coding_ends,
            max_operations=self.max_operations, max_bases=self.max_bases))

    def seed(self, refdata, *, v_call, d_call, j_call,
             reference_ends_are_coding_ends=False):
        """Install an admissible founder as a persistent native Simulation.

        Overlapping/nonmonotone layouts remain available through ``project``
        but cannot enter the current single-segment execution representation.
        Audited native continuation passes preserve founder ancestry; source
        rewrites and unsupported legacy output channels are explicitly gated.
        """
        if not isinstance(refdata, _engine.RefDataConfig):
            raise TypeError("seed() requires a GenAIRR RefDataConfig")
        if type(reference_ends_are_coding_ends) is not bool:
            raise ValueError("reference_ends_are_coding_ends must be a boolean")
        for name, value in (("v_call", v_call), ("d_call", d_call), ("j_call", j_call)):
            if not isinstance(value, str) or not value:
                raise ValueError(f"{name} must be a nonempty generating allele name")
        return _engine.coupled_duplex_seed(
            self.v, self.d, self.j, self.joins_json, refdata,
            v_call=v_call, d_call=d_call, j_call=j_call,
            d_inverted=self.d_inverted, reference_ends_are_coding_ends=reference_ends_are_coding_ends,
            max_operations=self.max_operations, max_bases=self.max_bases)

    def to_dict(self):
        return {"model": "experimental_coupled_duplex_repair_v1", "v": self.v,
                "d": self.d, "j": self.j, "d_inverted": self.d_inverted,
                "joins": self.joins, "max_operations": self.max_operations,
                "max_bases": self.max_bases}

    @classmethod
    def from_dict(cls, value):
        fields = {"model", "v", "d", "j", "d_inverted", "joins", "max_operations", "max_bases"}
        if not isinstance(value, dict) or set(value) != fields or value["model"] != "experimental_coupled_duplex_repair_v1":
            raise ValueError("unsupported coupled duplex repair script schema")
        return cls(value["v"], value["d"], value["j"], value["d_inverted"],
                   json.dumps(value["joins"], allow_nan=False),
                   value["max_operations"], value["max_bases"])
