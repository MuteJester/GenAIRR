"""Regressions for explicit overtrimming and portable bundled SHM models."""

import hashlib
import struct
import subprocess
import sys

import pytest

import GenAIRR as ga
from GenAIRR._s5f_loader import load_builtin_s5f_kernel

from ._p_productive_support import synthetic_cartridge


@pytest.mark.parametrize("strict", [False, True])
def test_explicit_combined_d_overtrim_is_a_python_error(strict):
    cfg = synthetic_cartridge({}, d="CCC")
    exp = ga.Experiment.on(cfg).recombine().trim(d_5=[(2, 1)], d_3=[(2, 1)])
    with pytest.raises((ValueError, ga.StrictSamplingError), match="trim_exceeds_allele_length"):
        exp.run_records(n=1, seed=0, strict=strict)


@pytest.mark.parametrize("name,expected", [
    ("hh_s5f", "0de250c6b58ffe44c0ce0b3a55b133dd78cd67359c54ea39a82ddeac4f5a9dde"),
    ("hh_s5f_60", "98968f4e8cba274b7d33f407f8f2914f89a6692926082b36358951b1c04f3a6c"),
    ("hh_s5f_opposite", "fea14f87e9e11383d343870c36477aa50b7e3046038064370d65740f4fd66545"),
    ("hkl_s5f", "b3ab7ca8b54f335931617fb6de811ca4221497a8c4316bb886e0ab682cae2b91"),
])
def test_builtin_s5f_values_preserve_original_float64_bits(name, expected):
    # Digests captured from the original bundled artifacts before conversion.
    mutability, substitution = load_builtin_s5f_kernel(name)
    packed = struct.pack("<5120d", *(mutability + substitution))
    assert hashlib.sha256(packed).hexdigest() == expected


def test_all_builtin_s5f_models_work_without_pandas_or_numpy():
    script = """
import importlib.abc
import sys

class BlockOptionalPackages(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname.split(".")[0] in {"pandas", "numpy"}:
            raise ModuleNotFoundError("Blocked optional package: " + fullname)

sys.meta_path.insert(0, BlockOptionalPackages())
import GenAIRR as ga
from GenAIRR._s5f_loader import available_s5f_kernels, load_builtin_s5f_kernel

for name in available_s5f_kernels():
    mu, sub = load_builtin_s5f_kernel(name)
    assert len(mu) == 1024 and len(sub) == 4096
    result = ga.Experiment.on("human_igh").recombine().mutate(
        model="s5f", s5f_model=name, count=2,
    ).run_records(n=2, seed=42)
    assert len(result) == 2
    assert all(record["n_mutations"] == 2 for record in result)
assert not any(name.split(".")[0] in {"pandas", "numpy"} for name in sys.modules)
"""
    completed = subprocess.run([sys.executable, "-c", script], capture_output=True, text=True)
    assert completed.returncode == 0, completed.stderr


def test_describe_names_d_inversion_and_its_probability():
    description = ga.Experiment.on("human_igh").recombine().invert_d(prob=0.25).describe()
    assert "D-segment inversion: 25%" in description
    assert "reverse-complemented D allele" in description
    assert "Unknown step" not in description
