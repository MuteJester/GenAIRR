"""Cartridge model estimators, split one module per estimator (behavior-
preserving). ReferenceCartridgeBuilder inherits the aggregate mixin."""
from .allele_usage import _AlleleUsageEstimatorMixin
from .trim import _TrimEstimatorMixin
from .np_lengths import _NpLengthEstimatorMixin
from .np_base_model import _NpBaseModelEstimatorMixin
from .p_nucleotide_lengths import _PNucleotideLengthEstimatorMixin


class _CartridgeEstimators(
    _AlleleUsageEstimatorMixin,
    _TrimEstimatorMixin,
    _NpLengthEstimatorMixin,
    _NpBaseModelEstimatorMixin,
    _PNucleotideLengthEstimatorMixin,
):
    __slots__ = ()


__all__ = ["_CartridgeEstimators"]
