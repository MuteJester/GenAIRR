"""
GenAIRR C backend — ctypes bindings to libgenairr.so.

This module provides CSimulator, a thin Python wrapper around the
C simulation engine. The C library handles all rearrangement, mutation,
corruption, and AIRR serialization.

Usage::

    from GenAIRR._native import CSimulator

    sim = CSimulator("/path/to/config.gdc")
    sim.set_feature("mutate", True)
    sim.set_param("min_mutation_rate", 0.01)
    results = sim.simulate(n=1000)  # list[dict]
"""

from __future__ import annotations

import csv
import ctypes
import os
import sys
import tempfile
from pathlib import Path
from typing import Optional

# ── Library loading ──────────────────────────────────────────────

_LIB: Optional[ctypes.CDLL] = None

# Platform-specific shared library name
if sys.platform == "darwin":
    _LIB_NAME = "libgenairr.dylib"
elif sys.platform == "win32":
    _LIB_NAME = "genairr.dll"
else:
    _LIB_NAME = "libgenairr.so"


def _find_library() -> str:
    """Find the GenAIRR shared library for the current platform."""
    candidates = []

    here = Path(__file__).parent  # src/GenAIRR/_native/

    # 1. Next to this file (bundled in wheel or editable install)
    candidates.append(here / _LIB_NAME)

    # 2. Co-located C source build directory (development mode)
    candidates.append(here / "csrc" / "build" / _LIB_NAME)

    # 3. MSVC multi-config build directories (Debug/Release subdirs)
    if sys.platform == "win32":
        candidates.append(here / "csrc" / "build" / "Release" / _LIB_NAME)
        candidates.append(here / "csrc" / "build" / "Debug" / _LIB_NAME)
        # 3b. Setuptools may place the DLL alongside the package root
        candidates.append(here.parent / _LIB_NAME)  # src/GenAIRR/
        candidates.append(here.parent.parent / _LIB_NAME)  # src/

    # 4. System paths (Unix only)
    if sys.platform != "win32":
        candidates.append(Path("/usr/local/lib") / _LIB_NAME)
        candidates.append(Path("/usr/lib") / _LIB_NAME)

    for path in candidates:
        if path.exists():
            return str(path)

    searched = "\n".join(f"  - {p}" for p in candidates)
    raise ImportError(
        f"Cannot find {_LIB_NAME}. Searched:\n{searched}\n"
        "Install with:\n"
        "  pip install -e .      (builds C backend automatically)\n"
        "Or build manually:\n"
        "  cd src/GenAIRR/_native/csrc && mkdir -p build && cd build\n"
        "  cmake .. -DCMAKE_BUILD_TYPE=Release && cmake --build ."
    )


def _load_lib() -> ctypes.CDLL:
    """Load the shared library and set up function signatures."""
    global _LIB
    if _LIB is not None:
        return _LIB

    lib = ctypes.CDLL(_find_library())

    # ── genairr_create ──
    lib.genairr_create.argtypes = [ctypes.c_char_p]
    lib.genairr_create.restype = ctypes.c_void_p

    # ── genairr_create_from_memory ──
    lib.genairr_create_from_memory.argtypes = [ctypes.c_void_p, ctypes.c_int]
    lib.genairr_create_from_memory.restype = ctypes.c_void_p

    # ── genairr_destroy ──
    lib.genairr_destroy.argtypes = [ctypes.c_void_p]
    lib.genairr_destroy.restype = None

    # ── genairr_set_feature ──
    lib.genairr_set_feature.argtypes = [ctypes.c_void_p, ctypes.c_char_p, ctypes.c_int]
    lib.genairr_set_feature.restype = ctypes.c_int

    # ── genairr_set_param_f64 ──
    lib.genairr_set_param_f64.argtypes = [ctypes.c_void_p, ctypes.c_char_p, ctypes.c_double]
    lib.genairr_set_param_f64.restype = ctypes.c_int

    # ── genairr_set_param_i32 ──
    lib.genairr_set_param_i32.argtypes = [ctypes.c_void_p, ctypes.c_char_p, ctypes.c_int]
    lib.genairr_set_param_i32.restype = ctypes.c_int

    # ── genairr_set_seed ──
    lib.genairr_set_seed.argtypes = [ctypes.c_void_p, ctypes.c_uint64]
    lib.genairr_set_seed.restype = None

    # ── genairr_simulate ──
    lib.genairr_simulate.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_char_p]
    lib.genairr_simulate.restype = ctypes.c_int

    # ── genairr_version ──
    lib.genairr_version.argtypes = []
    lib.genairr_version.restype = ctypes.c_char_p

    # ── genairr_lock_allele ──
    lib.genairr_lock_allele.argtypes = [ctypes.c_void_p, ctypes.c_char_p, ctypes.c_char_p]
    lib.genairr_lock_allele.restype = ctypes.c_int

    # ── genairr_clear_locks ──
    lib.genairr_clear_locks.argtypes = [ctypes.c_void_p, ctypes.c_char_p]
    lib.genairr_clear_locks.restype = None

    # ── genairr_simulate_one ──
    lib.genairr_simulate_one.argtypes = [ctypes.c_void_p, ctypes.c_char_p, ctypes.c_int]
    lib.genairr_simulate_one.restype = ctypes.c_int

    # ── genairr_get_header ──
    lib.genairr_get_header.argtypes = [ctypes.c_char_p, ctypes.c_int]
    lib.genairr_get_header.restype = ctypes.c_int

    # ── genairr_set_trace ──
    lib.genairr_set_trace.argtypes = [ctypes.c_void_p, ctypes.c_int]
    lib.genairr_set_trace.restype = None

    # ── genairr_get_trace ──
    lib.genairr_get_trace.argtypes = [ctypes.c_void_p, ctypes.c_char_p, ctypes.c_int]
    lib.genairr_get_trace.restype = ctypes.c_int

    # ── genairr_clear_trace ──
    lib.genairr_clear_trace.argtypes = [ctypes.c_void_p]
    lib.genairr_clear_trace.restype = None

    # ── genairr_set_hooks ──
    lib.genairr_set_hooks.argtypes = [ctypes.c_void_p, ctypes.c_uint32]
    lib.genairr_set_hooks.restype = None

    # ── genairr_get_snapshot_count ──
    lib.genairr_get_snapshot_count.argtypes = [ctypes.c_void_p]
    lib.genairr_get_snapshot_count.restype = ctypes.c_int

    # ── genairr_get_snapshot_name ──
    lib.genairr_get_snapshot_name.argtypes = [ctypes.c_void_p, ctypes.c_int]
    lib.genairr_get_snapshot_name.restype = ctypes.c_char_p

    # ── genairr_dump_snapshot ──
    lib.genairr_dump_snapshot.argtypes = [
        ctypes.c_void_p, ctypes.c_int, ctypes.c_char_p, ctypes.c_int
    ]
    lib.genairr_dump_snapshot.restype = ctypes.c_int

    # ── genairr_dump_snapshot_codon_rail ──
    lib.genairr_dump_snapshot_codon_rail.argtypes = [
        ctypes.c_void_p, ctypes.c_int, ctypes.c_char_p, ctypes.c_int
    ]
    lib.genairr_dump_snapshot_codon_rail.restype = ctypes.c_int

    # ── ABI version sanity check (T1-12) ────────────────────────
    # The Python package and the loaded C library must report the
    # same version string. A mismatch typically means a stale
    # `libgenairr.so` was found on the search path (e.g. an old
    # editable install left behind, or a system-wide build of a
    # different version). Loading a mismatched binary yields silent
    # ABI corruption — fail loudly at import instead.
    _check_abi_version(lib)

    _LIB = lib
    return lib


def _check_abi_version(lib: ctypes.CDLL) -> None:
    """Assert that the loaded C library version matches Python's.

    Reads __version__ from the parent GenAIRR package (already
    populated by importlib.metadata at package import time) and
    compares it to the C library's genairr_version() return value.
    """
    try:
        from GenAIRR import __version__ as _py_version
    except ImportError:
        return  # parent package not initialized — skip check
    c_version = lib.genairr_version().decode("ascii", errors="replace")
    if c_version != _py_version:
        raise RuntimeError(
            f"GenAIRR C library version mismatch: loaded "
            f"libgenairr reports {c_version!r}, but Python package "
            f"is {_py_version!r}. Reinstall to resolve: "
            f"`pip install --force-reinstall --no-cache-dir GenAIRR`."
        )


def get_version() -> str:
    """Return the C library version string."""
    lib = _load_lib()
    return lib.genairr_version().decode()


# ── AIRR column order (cached from C) ───────────────────────────

_AIRR_COLUMNS: Optional[list] = None


def _get_airr_columns() -> list:
    """Get the AIRR TSV column names from the C library (cached)."""
    global _AIRR_COLUMNS
    if _AIRR_COLUMNS is not None:
        return _AIRR_COLUMNS
    lib = _load_lib()
    buf = ctypes.create_string_buffer(4096)
    rc = lib.genairr_get_header(buf, 4096)
    if rc < 0:
        raise RuntimeError("Failed to get AIRR header from C library")
    _AIRR_COLUMNS = buf.value.decode().split('\t')
    return _AIRR_COLUMNS


# ── Boolean AIRR fields (TSV "T"/"F" → Python bool) ─────────────

_BOOL_FIELDS = {
    "productive", "stop_codon", "vj_in_frame",
    "is_reverse_complement", "is_contaminant",
    "d_inverted", "receptor_revised",
}

_INT_FIELDS = {
    "v_sequence_start", "v_sequence_end",
    "v_germline_start", "v_germline_end",
    "d_sequence_start", "d_sequence_end",
    "d_germline_start", "d_germline_end",
    "j_sequence_start", "j_sequence_end",
    "j_germline_start", "j_germline_end",
    "junction_start", "junction_end", "junction_length",
    "v_trim_5", "v_trim_3", "d_trim_5", "d_trim_3", "j_trim_5", "j_trim_3",
    "np1_length", "np2_length",
    "n_mutations", "n_insertions", "n_deletions",
    "n_sequencing_errors", "n_pcr_errors",
    "revision_footprint_length",
    "sequence_length",
}

_FLOAT_FIELDS = {
    "mutation_rate",
}


def _parse_row(row: dict) -> dict:
    """Convert TSV string values to proper Python types."""
    for key in _BOOL_FIELDS:
        if key in row:
            row[key] = row[key] == "T"
    for key in _INT_FIELDS:
        if key in row:
            try:
                row[key] = int(row[key])
            except (ValueError, TypeError):
                row[key] = 0
    for key in _FLOAT_FIELDS:
        if key in row:
            try:
                row[key] = float(row[key])
            except (ValueError, TypeError):
                row[key] = 0.0
    return row


# ── CSimulator ───────────────────────────────────────────────────

class CSimulator:
    """
    Python wrapper around the C simulation engine.

    Manages the lifecycle of a GenAIRRSimulator handle and provides
    Pythonic methods for configuring features and running simulations.
    """

    def __init__(self, gdc_path: Optional[str] = None,
                 gdc_bytes: Optional[bytes] = None):
        """
        Create a simulator from a .gdc file or in-memory bytes.

        Args:
            gdc_path: Path to a .gdc config file.
            gdc_bytes: Raw .gdc bytes (alternative to path).
        """
        self._lib = _load_lib()

        if gdc_path is not None:
            self._handle = self._lib.genairr_create(gdc_path.encode())
        elif gdc_bytes is not None:
            buf = ctypes.create_string_buffer(gdc_bytes)
            self._handle = self._lib.genairr_create_from_memory(buf, len(gdc_bytes))
        else:
            raise ValueError("Either gdc_path or gdc_bytes must be provided")

        if not self._handle:
            raise RuntimeError(
                f"Failed to create simulator from "
                f"{'gdc_path=' + repr(gdc_path) if gdc_path else 'gdc_bytes'}"
            )

    def __del__(self):
        self.destroy()

    def destroy(self):
        """Free the C simulator handle."""
        if hasattr(self, '_handle') and self._handle:
            self._lib.genairr_destroy(self._handle)
            self._handle = None

    def set_feature(self, name: str, enabled: bool) -> None:
        """Enable or disable a named feature."""
        rc = self._lib.genairr_set_feature(
            self._handle, name.encode(), 1 if enabled else 0
        )
        if rc != 0:
            raise ValueError(f"Unknown feature: {name!r}")

    def set_param(self, name: str, value) -> None:
        """Set a named parameter (auto-detects int vs float)."""
        if isinstance(value, float):
            rc = self._lib.genairr_set_param_f64(
                self._handle, name.encode(), value
            )
        elif isinstance(value, int):
            rc = self._lib.genairr_set_param_i32(
                self._handle, name.encode(), value
            )
        else:
            raise TypeError(f"Parameter value must be int or float, got {type(value)}")

        if rc != 0:
            raise ValueError(f"Unknown parameter: {name!r}")

    def set_seed(self, seed: int) -> None:
        """Set the random seed (0 = time-based)."""
        self._lib.genairr_set_seed(self._handle, seed)

    def lock_allele(self, segment: str, name: str) -> None:
        """Lock sampling of a segment to a specific allele by name.

        Call multiple times to allow multiple alleles for one segment.

        Args:
            segment: "v", "d", "j", or "c".
            name: Allele name (e.g. "IGHV1-2*01").
        """
        rc = self._lib.genairr_lock_allele(
            self._handle, segment.encode(), name.encode()
        )
        if rc == -1:
            raise ValueError(
                f"Unknown segment {segment!r} or allele {name!r} not found in pool"
            )
        if rc == -2:
            raise ValueError(f"Too many alleles locked for segment {segment!r} (max 64)")

    def clear_locks(self, segment: str = None) -> None:
        """Clear locked alleles for a segment (or all segments if None)."""
        seg = segment.encode() if segment else None
        self._lib.genairr_clear_locks(self._handle, seg)

    def simulate(self, n: int = 1) -> list[dict]:
        """
        Run N simulations and return AIRR-format dicts.

        Returns:
            List of dicts, one per simulated sequence. Keys match AIRR fields.
        """
        if not self._handle:
            raise RuntimeError("Simulator has been destroyed")

        # Write TSV to a temp file, read back as Python dicts
        with tempfile.NamedTemporaryFile(
            mode='w', suffix='.tsv', delete=False
        ) as tmp:
            tmp_path = tmp.name

        try:
            rc = self._lib.genairr_simulate(
                self._handle, n, tmp_path.encode()
            )
            if rc < 0:
                raise RuntimeError("C simulation failed")

            # Parse TSV output
            with open(tmp_path, 'r') as f:
                reader = csv.DictReader(f, delimiter='\t')
                results = [_parse_row(row) for row in reader]

            return results
        finally:
            try:
                os.unlink(tmp_path)
            except OSError:
                pass

    def simulate_one(self) -> dict:
        """
        Simulate a single sequence and return it as an AIRR dict.

        Unlike simulate(), this does not write to a temp file — the C engine
        writes a single TSV row directly to an in-memory buffer. On the first
        call the RNG is seeded; subsequent calls continue the stream.

        Returns:
            Dict with AIRR fields for one simulated sequence.
        """
        if not self._handle:
            raise RuntimeError("Simulator has been destroyed")

        buf = ctypes.create_string_buffer(32768)
        rc = self._lib.genairr_simulate_one(self._handle, buf, 32768)
        if rc < 0:
            raise RuntimeError("simulate_one failed (buffer too small or C error)")

        columns = _get_airr_columns()
        values = buf.value.decode().split('\t')
        row = dict(zip(columns, values))
        return _parse_row(row)

    def set_trace(self, enabled: bool) -> None:
        """Enable or disable the execution trace log."""
        if not self._handle:
            raise RuntimeError("Simulator has been destroyed")
        self._lib.genairr_set_trace(self._handle, 1 if enabled else 0)

    def get_trace(self) -> str:
        """Return the trace log contents from the last simulate_one() call."""
        if not self._handle:
            raise RuntimeError("Simulator has been destroyed")
        buf = ctypes.create_string_buffer(1 << 20)  # 1 MB
        rc = self._lib.genairr_get_trace(self._handle, buf, 1 << 20)
        if rc < 0:
            return ""
        return buf.value.decode()

    def clear_trace(self) -> None:
        """Clear the trace buffer."""
        if not self._handle:
            raise RuntimeError("Simulator has been destroyed")
        self._lib.genairr_clear_trace(self._handle)

    # ── Hook point constants (mirror C enum) ──────────────────────

    HOOK_POINTS = {
        "post_assembly": 0,
        "post_functionality": 1,
        "post_d_inversion": 2,
        "post_receptor_rev": 3,
        "post_mutation": 4,
        "post_selection": 5,
        "post_corrupt_5": 6,
        "post_corrupt_3": 7,
        "post_indels": 8,
        "post_ns": 9,
        "post_pcr": 10,
        "post_quality": 11,
        "final": 12,
    }

    def set_hooks(self, hook_names: list[str]) -> None:
        """Enable snapshot capture at specified hook points."""
        if not self._handle:
            raise RuntimeError("Simulator has been destroyed")
        mask = 0
        for name in hook_names:
            if name not in self.HOOK_POINTS:
                raise ValueError(
                    f"Unknown hook: {name!r}. Valid: {list(self.HOOK_POINTS)}"
                )
            mask |= (1 << self.HOOK_POINTS[name])
        self._lib.genairr_set_hooks(self._handle, mask)

    def get_snapshots(self) -> list[dict]:
        """Get all captured snapshots from last simulate_one().

        Returns:
            List of dicts with keys: hook (str), nodes (list[dict]).
            Each node dict has: pos, cur, germ, seg, flags, gp, ph, aa, prod.
        """
        import json

        if not self._handle:
            raise RuntimeError("Simulator has been destroyed")

        count = self._lib.genairr_get_snapshot_count(self._handle)
        results = []
        buf_size = 1 << 20  # 1 MB
        buf = ctypes.create_string_buffer(buf_size)

        for i in range(count):
            name_bytes = self._lib.genairr_get_snapshot_name(self._handle, i)
            name = name_bytes.decode() if name_bytes else "unknown"

            rc = self._lib.genairr_dump_snapshot(self._handle, i, buf, buf_size)
            if rc > 0:
                nodes = json.loads(buf.value.decode())
            else:
                nodes = []

            results.append({"hook": name, "nodes": nodes})

        return results

    def get_snapshot_codon_rail(self, index: int) -> list[dict]:
        """Get the codon rail from a specific snapshot.

        Returns:
            List of dicts with keys: idx, bases, aa, is_stop, seg.
        """
        import json

        if not self._handle:
            raise RuntimeError("Simulator has been destroyed")

        buf_size = 1 << 20
        buf = ctypes.create_string_buffer(buf_size)
        rc = self._lib.genairr_dump_snapshot_codon_rail(
            self._handle, index, buf, buf_size
        )
        if rc > 0:
            return json.loads(buf.value.decode())
        return []

    def simulate_one_hooked(self, hooks: list[str]) -> tuple[dict, list[dict], str]:
        """Simulate one sequence with hooks. Returns (airr_record, snapshots, trace).

        Args:
            hooks: List of hook point names (e.g. ["post_assembly", "post_mutation", "final"]).

        Returns:
            Tuple of (airr_record, snapshots, trace_log).
        """
        self.set_hooks(hooks)
        self.set_trace(True)
        rec = self.simulate_one()
        snapshots = self.get_snapshots()
        trace = self.get_trace()
        self.set_trace(False)
        self.set_hooks([])
        return rec, snapshots, trace

    def simulate_to_file(self, n: int, output_path: str) -> int:
        """
        Run N simulations and write AIRR TSV directly to a file.

        More efficient than simulate() for large runs — no Python
        parsing overhead.

        Returns:
            Number of sequences written.
        """
        if not self._handle:
            raise RuntimeError("Simulator has been destroyed")

        rc = self._lib.genairr_simulate(
            self._handle, n, output_path.encode()
        )
        if rc < 0:
            raise RuntimeError("C simulation failed")
        return rc
