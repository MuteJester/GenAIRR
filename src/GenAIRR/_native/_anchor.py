"""
GenAIRR anchor resolution — Python bindings to the C subsystem.

The C side (see csrc/include/genairr/anchor.h, ref_record.h, ref_loader.h)
implements:

  * ``AnchorResolver`` — strategy stack: explicit → IMGT-gapped → motif
  * ``ReferenceLoader`` — vtable-based reader producing
    ``LoadedAlleleRecord`` values one at a time
  * ``BuildReport`` — Builder accumulator for per-allele resolution
    outcomes (Phase 2/3 will wire this through)

This module is the thin ctypes layer that mirrors those types and
exposes them as Pythonic dataclasses, an iterator/context-manager
``ReferenceLoader`` class, and a ``resolve_anchor()`` function.

OWNERSHIP NOTE
--------------
The C side returns string fields in ``LoadedAlleleRecord`` as raw
``const char*`` pointers into loader-owned memory; those pointers are
valid only until the next ``next()`` call or ``close()``. The Python
``LoadedAlleleRecord`` dataclass copies the strings into Python ``str``
on construction, so the returned objects are safe to retain across
loader calls.
"""
from __future__ import annotations

import ctypes
from dataclasses import dataclass, field
from enum import IntEnum
from typing import Optional

from . import _load_lib


# ── Enums (mirror the C side; values must stay in lock-step) ───────


class Locus(IntEnum):
    UNKNOWN = 0
    IGH = 1
    IGK = 2
    IGL = 3
    TRA = 4
    TRB = 5
    TRD = 6
    TRG = 7


class Segment(IntEnum):
    """Mirrors the C `Segment` enum (defined in csrc/include/genairr/types.h).

    Note: the simulator's runtime enum has additional values (NP1, NP2,
    UMI, ADAPTER) that are never used for reference alleles. We expose
    only the meaningful ones; SEG_COUNT (=8) is the "unknown / not yet
    classified" sentinel."""
    V = 0
    NP1 = 1
    D = 2
    NP2 = 3
    J = 4
    C = 5
    UMI = 6
    ADAPTER = 7
    UNKNOWN = 8       # alias for SEG_COUNT


class FunctionalStatus(IntEnum):
    UNKNOWN = 0
    F = 1
    ORF = 2
    PSEUDO = 3
    PARTIAL = 4


class AnchorResidue(IntEnum):
    UNKNOWN = 0
    CYS = 1
    TRP = 2
    PHE = 3


class AnchorConfidence(IntEnum):
    REJECTED = 0
    BEST_GUESS = 1
    ALTERNATIVE = 2
    CANONICAL = 3


class AnchorMethod(IntEnum):
    NONE = 0
    OVERRIDE = 1
    CUSTOM = 2
    EXPLICIT = 3
    IMGT_GAPPED = 4
    MOTIF_SEARCH = 5


# ── ctypes Structures (mirror C struct layouts) ───────────────────


_REF_MAX_ALIASES = 8


class _CAnchorResult(ctypes.Structure):
    """Layout-mirror of `AnchorResult` from anchor.h.

    Order and types must match the C struct exactly. `reason` is a
    ``const char*`` that the C side either leaves as NULL or sets to
    a static string literal — we never own it on the Python side."""
    _fields_ = [
        ("position",   ctypes.c_int),
        ("codon",      ctypes.c_char * 4),
        ("residue",    ctypes.c_int),
        ("confidence", ctypes.c_int),
        ("method",     ctypes.c_int),
        ("reason",     ctypes.c_char_p),
    ]


class _CLoadedAlleleRecord(ctypes.Structure):
    """Layout-mirror of `LoadedAlleleRecord` from ref_record.h."""
    _fields_ = [
        ("name",                 ctypes.c_char_p),
        ("aliases",              ctypes.c_char_p * _REF_MAX_ALIASES),
        ("n_aliases",            ctypes.c_int),
        ("segment",              ctypes.c_int),
        ("locus",                ctypes.c_int),
        ("species",              ctypes.c_char_p),
        ("sequence",             ctypes.c_char_p),
        ("sequence_length",      ctypes.c_int),
        ("gapped_sequence",      ctypes.c_char_p),
        ("gapped_length",        ctypes.c_int),
        ("gap_convention_imgt",  ctypes.c_bool),
        ("functional_status",    ctypes.c_int),
        ("explicit_anchor",      ctypes.c_int),
        ("source",               ctypes.c_char_p),
    ]


class _CAnchorResolverConfig(ctypes.Structure):
    """Layout-mirror of `AnchorResolverConfig`. The `custom_finder`
    function pointer field is exposed as ``c_void_p`` so we can leave
    it NULL — calling Python from C requires CFUNCTYPE marshalling
    that we deliberately defer to Phase 3."""
    _fields_ = [
        ("locus",             ctypes.c_int),
        ("allow_trp_v",       ctypes.c_bool),
        ("strict",            ctypes.c_bool),
        ("custom_finder",     ctypes.c_void_p),
        ("custom_user_data",  ctypes.c_void_p),
    ]


class _CLoaderHints(ctypes.Structure):
    """Layout-mirror of `LoaderHints` from ref_loader.h."""
    _fields_ = [
        ("locus_hint",   ctypes.c_int),
        ("segment_hint", ctypes.c_int),
    ]


# ── Wire ctypes function signatures (called once at module import) ──


def _bind_signatures(lib: ctypes.CDLL) -> None:
    """Set argtypes / restype on every C function we call.

    Required: ctypes can't infer signatures from the compiled symbol;
    without this it defaults to `int` arguments and pointers, which
    silently corrupts wide types like ``c_char_p`` or struct returns."""

    # ── anchor_resolver_config_init(AnchorResolverConfig *cfg) ──
    lib.anchor_resolver_config_init.argtypes = [
        ctypes.POINTER(_CAnchorResolverConfig),
    ]
    lib.anchor_resolver_config_init.restype = None

    # ── AnchorResult anchor_resolve_v(const cfg*, const rec*) ──
    lib.anchor_resolve_v.argtypes = [
        ctypes.POINTER(_CAnchorResolverConfig),
        ctypes.POINTER(_CLoadedAlleleRecord),
    ]
    lib.anchor_resolve_v.restype = _CAnchorResult

    lib.anchor_resolve_j.argtypes = [
        ctypes.POINTER(_CAnchorResolverConfig),
        ctypes.POINTER(_CLoadedAlleleRecord),
    ]
    lib.anchor_resolve_j.restype = _CAnchorResult

    # ── loaded_allele_record_init(rec*) ──
    lib.loaded_allele_record_init.argtypes = [
        ctypes.POINTER(_CLoadedAlleleRecord),
    ]
    lib.loaded_allele_record_init.restype = None

    # ── ReferenceLoader factories + flat wrappers ──
    lib.imgt_vquest_loader_open.argtypes = [
        ctypes.c_char_p,                            # fasta_path
        ctypes.c_int,                               # segment_hint
        ctypes.POINTER(ctypes.c_char_p),            # err_msg out-param
    ]
    lib.imgt_vquest_loader_open.restype = ctypes.c_void_p   # opaque ReferenceLoader*

    # Phase 2 loaders
    lib.plain_fasta_loader_open.argtypes = [
        ctypes.c_char_p, ctypes.c_int, ctypes.c_int,
        ctypes.POINTER(ctypes.c_char_p),
    ]
    lib.plain_fasta_loader_open.restype = ctypes.c_void_p

    lib.airrc_germline_loader_open.argtypes = [
        ctypes.c_char_p, ctypes.POINTER(ctypes.c_char_p),
    ]
    lib.airrc_germline_loader_open.restype = ctypes.c_void_p

    lib.ogrdb_loader_open.argtypes = [
        ctypes.c_char_p,                # fasta_path
        ctypes.c_char_p,                # json_sidecar_path (NULL OK)
        ctypes.c_int, ctypes.c_int,
        ctypes.POINTER(ctypes.c_char_p),
    ]
    lib.ogrdb_loader_open.restype = ctypes.c_void_p

    # Phase 3 loader
    lib.igblast_loader_open.argtypes = [
        ctypes.c_char_p,                # fasta_path
        ctypes.c_char_p,                # aux/ndm sidecar (NULL OK)
        ctypes.c_int, ctypes.c_int,
        ctypes.POINTER(ctypes.c_char_p),
    ]
    lib.igblast_loader_open.restype = ctypes.c_void_p

    # Auto-dispatcher — takes a LoaderHints* by pointer.
    lib.reference_loader_open_auto.argtypes = [
        ctypes.c_char_p,
        ctypes.POINTER(_CLoaderHints),
        ctypes.POINTER(ctypes.c_char_p),
    ]
    lib.reference_loader_open_auto.restype = ctypes.c_void_p

    lib.reference_loader_next.argtypes = [
        ctypes.c_void_p,
        ctypes.POINTER(_CLoadedAlleleRecord),
        ctypes.POINTER(ctypes.c_char_p),
    ]
    lib.reference_loader_next.restype = ctypes.c_int

    lib.reference_loader_close.argtypes = [ctypes.c_void_p]
    lib.reference_loader_close.restype = None

    # ── Locus / segment inference ──
    lib.locus_from_gene_name.argtypes = [ctypes.c_char_p]
    lib.locus_from_gene_name.restype = ctypes.c_int

    lib.segment_from_gene_name.argtypes = [ctypes.c_char_p]
    lib.segment_from_gene_name.restype = ctypes.c_int

    lib.locus_from_filename.argtypes = [ctypes.c_char_p]
    lib.locus_from_filename.restype = ctypes.c_int

    # ── Enum-name lookups (useful for diagnostics) ──
    for fn in (lib.anchor_residue_name, lib.anchor_confidence_name,
               lib.anchor_method_name, lib.locus_name, lib.segment_name,
               lib.functional_status_name):
        fn.argtypes = [ctypes.c_int]
        fn.restype = ctypes.c_char_p


_lib = _load_lib()
_bind_signatures(_lib)


# ── Pythonic dataclasses ───────────────────────────────────────────


@dataclass(frozen=True)
class AnchorResult:
    """Outcome of a single anchor-resolution attempt.

    Mirrors the C `AnchorResult` but converts the ctypes-friendly enum
    integers into Python IntEnums and copies the codon / reason out
    so the dataclass survives the C string going stale.
    """
    position: int
    codon: str
    residue: AnchorResidue
    confidence: AnchorConfidence
    method: AnchorMethod
    reason: Optional[str]

    @classmethod
    def _from_c(cls, c: _CAnchorResult) -> "AnchorResult":
        return cls(
            position=int(c.position),
            codon=c.codon.decode("ascii"),
            residue=AnchorResidue(int(c.residue)),
            confidence=AnchorConfidence(int(c.confidence)),
            method=AnchorMethod(int(c.method)),
            reason=c.reason.decode("ascii") if c.reason else None,
        )

    @property
    def accepted(self) -> bool:
        """True if any strategy yielded a usable anchor."""
        return self.confidence != AnchorConfidence.REJECTED


@dataclass(frozen=True)
class LoadedAlleleRecord:
    """Format-agnostic input record for the anchor resolver.

    All string fields are Python ``str`` (copied out of loader memory
    on construction) so instances are safe to retain across loader
    iterations.
    """
    name: str
    aliases: tuple[str, ...]
    segment: Segment
    locus: Locus
    species: Optional[str]
    sequence: str
    gapped_sequence: Optional[str]
    gap_convention_imgt: bool
    functional_status: FunctionalStatus
    explicit_anchor: int
    source: str

    @classmethod
    def _from_c(cls, c: _CLoadedAlleleRecord) -> "LoadedAlleleRecord":
        aliases = tuple(
            c.aliases[i].decode("ascii")
            for i in range(c.n_aliases)
            if c.aliases[i] is not None
        )
        # ctypes c_char_p decodes bytes/None automatically; we just
        # need to handle the optional-None cases.
        return cls(
            name=(c.name or b"").decode("ascii"),
            aliases=aliases,
            segment=Segment(int(c.segment)),
            locus=Locus(int(c.locus)),
            species=c.species.decode("ascii") if c.species else None,
            sequence=(c.sequence or b"").decode("ascii"),
            gapped_sequence=(c.gapped_sequence.decode("ascii")
                             if c.gapped_sequence else None),
            gap_convention_imgt=bool(c.gap_convention_imgt),
            functional_status=FunctionalStatus(int(c.functional_status)),
            explicit_anchor=int(c.explicit_anchor),
            source=(c.source or b"").decode("ascii"),
        )

    def _to_c(self) -> _CLoadedAlleleRecord:
        """Render back into a C struct so a Python-constructed record
        can be passed to the resolver. The bytes objects must outlive
        the C call — we keep them on the returned struct's `_keepalive`
        list, which is a small dynamic attribute holding refs."""
        c = _CLoadedAlleleRecord()
        keepalive = []

        def _bytes_or_none(s: Optional[str]) -> Optional[bytes]:
            if s is None:
                return None
            b = s.encode("ascii")
            keepalive.append(b)
            return b

        name_b = _bytes_or_none(self.name) or b""
        c.name = name_b
        for i in range(_REF_MAX_ALIASES):
            if i < len(self.aliases):
                c.aliases[i] = _bytes_or_none(self.aliases[i])
            else:
                c.aliases[i] = None
        c.n_aliases = len(self.aliases)
        c.segment = int(self.segment)
        c.locus = int(self.locus)
        c.species = _bytes_or_none(self.species)
        c.sequence = _bytes_or_none(self.sequence)
        c.sequence_length = len(self.sequence)
        c.gapped_sequence = _bytes_or_none(self.gapped_sequence)
        c.gapped_length = len(self.gapped_sequence) if self.gapped_sequence else 0
        c.gap_convention_imgt = bool(self.gap_convention_imgt)
        c.functional_status = int(self.functional_status)
        c.explicit_anchor = int(self.explicit_anchor)
        c.source = _bytes_or_none(self.source)

        # Attach keepalive as a plain attribute. ctypes Structures allow
        # arbitrary attributes; this prevents premature GC of the bytes
        # backing the c_char_p fields.
        c._keepalive = keepalive  # type: ignore[attr-defined]
        return c


# ── ReferenceLoader (Iterator + Context Manager) ──────────────────


class ReferenceLoader:
    """Iterator + context manager wrapping a C ``ReferenceLoader *``.

    Today only the IMGT V-QUEST FASTA backend is wired; OGRDB / AIRR-C /
    IgBLAST / VDJbase loaders land in Phase 2.

    Usage::

        with ReferenceLoader.open_imgt("path/to/IGHV.fasta") as loader:
            for rec in loader:
                ...
    """

    def __init__(self, handle: int):
        if handle == 0 or handle is None:
            raise ValueError("ReferenceLoader handle is NULL")
        self._handle: Optional[int] = handle

    @classmethod
    def open_imgt(cls,
                  fasta_path: str,
                  *,
                  segment_hint: Segment = Segment.UNKNOWN) -> "ReferenceLoader":
        """Open an IMGT V-QUEST FASTA reference (gapped, pipe headers)."""
        err = ctypes.c_char_p()
        handle = _lib.imgt_vquest_loader_open(
            fasta_path.encode("utf-8"),
            ctypes.c_int(int(segment_hint)),
            ctypes.byref(err),
        )
        if not handle:
            msg = err.value.decode("ascii") if err.value else "open failed"
            raise OSError(f"imgt_vquest_loader_open({fasta_path!r}): {msg}")
        return cls(handle)

    @classmethod
    def open_plain(cls,
                   fasta_path: str,
                   *,
                   locus_hint: Locus = Locus.UNKNOWN,
                   segment_hint: Segment = Segment.UNKNOWN) -> "ReferenceLoader":
        """Open a bare FASTA (custom lab references; no metadata).

        Locus hint is used by the J-anchor resolver to set canonical
        vs alternative confidence; pass ``Locus.UNKNOWN`` if locus
        is genuinely unknown (results downgrade to ``BEST_GUESS``).
        """
        err = ctypes.c_char_p()
        handle = _lib.plain_fasta_loader_open(
            fasta_path.encode("utf-8"),
            ctypes.c_int(int(locus_hint)),
            ctypes.c_int(int(segment_hint)),
            ctypes.byref(err),
        )
        if not handle:
            msg = err.value.decode("ascii") if err.value else "open failed"
            raise OSError(f"plain_fasta_loader_open({fasta_path!r}): {msg}")
        return cls(handle)

    @classmethod
    def open_airrc(cls, json_path: str) -> "ReferenceLoader":
        """Open an AIRR-Community GermlineSet JSON reference.

        AIRR-C is the modern source-of-truth format with explicit
        anchor coordinates (V cdr3_start, J cdr3_end). Coordinates
        are 1-based against the unaligned sequence; the loader
        converts to our 0-based internal convention.
        """
        err = ctypes.c_char_p()
        handle = _lib.airrc_germline_loader_open(
            json_path.encode("utf-8"),
            ctypes.byref(err),
        )
        if not handle:
            msg = err.value.decode("ascii") if err.value else "open failed"
            raise OSError(f"airrc_germline_loader_open({json_path!r}): {msg}")
        return cls(handle)

    @classmethod
    def open_ogrdb(cls,
                   fasta_path: str,
                   json_sidecar_path: Optional[str] = None,
                   *,
                   locus_hint: Locus = Locus.UNKNOWN,
                   segment_hint: Segment = Segment.UNKNOWN) -> "ReferenceLoader":
        """Open an OGRDB FASTA + JSON sidecar combo.

        If ``json_sidecar_path`` is None, falls back to FASTA-only
        loading (anchors come from motif resolution; functional
        status is FUNC_UNKNOWN for every record).
        """
        err = ctypes.c_char_p()
        sidecar = json_sidecar_path.encode("utf-8") if json_sidecar_path else None
        handle = _lib.ogrdb_loader_open(
            fasta_path.encode("utf-8"),
            sidecar,
            ctypes.c_int(int(locus_hint)),
            ctypes.c_int(int(segment_hint)),
            ctypes.byref(err),
        )
        if not handle:
            msg = err.value.decode("ascii") if err.value else "open failed"
            raise OSError(f"ogrdb_loader_open({fasta_path!r}): {msg}")
        return cls(handle)

    @classmethod
    def open_igblast(cls,
                     fasta_path: str,
                     aux_or_ndm_path: Optional[str] = None,
                     *,
                     locus_hint: Locus = Locus.UNKNOWN,
                     segment_hint: Segment = Segment.UNKNOWN) -> "ReferenceLoader":
        """Open an IgBLAST custom-DB FASTA + sidecar.

        ``aux_or_ndm_path`` may be either:
          - a ``_gl.aux`` file (J anchors, 0-based)
          - a ``.ndm.imgt`` file (V coordinates, 1-based)
          - or ``None`` (motif fallback for anchors)
        The loader sniffs the path extension to dispatch.
        """
        err = ctypes.c_char_p()
        sidecar = aux_or_ndm_path.encode("utf-8") if aux_or_ndm_path else None
        handle = _lib.igblast_loader_open(
            fasta_path.encode("utf-8"),
            sidecar,
            ctypes.c_int(int(locus_hint)),
            ctypes.c_int(int(segment_hint)),
            ctypes.byref(err),
        )
        if not handle:
            msg = err.value.decode("ascii") if err.value else "open failed"
            raise OSError(f"igblast_loader_open({fasta_path!r}): {msg}")
        return cls(handle)

    @classmethod
    def open_auto(cls,
                  path: str,
                  *,
                  locus_hint: Locus = Locus.UNKNOWN,
                  segment_hint: Segment = Segment.UNKNOWN) -> "ReferenceLoader":
        """Auto-detect the reference format from path + file content.

        Decision rules (in order):
          1. ``.json`` extension → AIRR-C
          2. First non-whitespace char is ``{`` or ``[`` → AIRR-C
          3. FASTA with sibling ``<base>.json`` → OGRDB
          4. FASTA with ``|``-delimited header → IMGT V-QUEST
          5. Otherwise → plain FASTA

        Raises ``OSError`` if no rule matches or the file can't be
        opened.
        """
        err = ctypes.c_char_p()
        hints = _CLoaderHints(
            locus_hint=int(locus_hint),
            segment_hint=int(segment_hint),
        )
        handle = _lib.reference_loader_open_auto(
            path.encode("utf-8"),
            ctypes.byref(hints),
            ctypes.byref(err),
        )
        if not handle:
            msg = err.value.decode("ascii") if err.value else "open failed"
            raise OSError(f"reference_loader_open_auto({path!r}): {msg}")
        return cls(handle)

    def __iter__(self):
        return self

    def __next__(self) -> LoadedAlleleRecord:
        if self._handle is None:
            raise RuntimeError("ReferenceLoader is closed")
        rec = _CLoadedAlleleRecord()
        err = ctypes.c_char_p()
        rc = _lib.reference_loader_next(
            ctypes.c_void_p(self._handle),
            ctypes.byref(rec),
            ctypes.byref(err),
        )
        if rc == 1:
            raise StopIteration
        if rc < 0:
            msg = err.value.decode("ascii") if err.value else "loader error"
            raise RuntimeError(f"ReferenceLoader.next(): {msg}")
        return LoadedAlleleRecord._from_c(rec)

    def close(self) -> None:
        """Release the C handle. Idempotent and safe to call after
        iteration has ended."""
        if self._handle is not None:
            _lib.reference_loader_close(ctypes.c_void_p(self._handle))
            self._handle = None

    def __enter__(self) -> "ReferenceLoader":
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        self.close()

    def __del__(self):
        # Match T2-5 pattern: best-effort cleanup at GC time.
        try:
            self.close()
        except Exception:
            pass


# ── Resolver entry point ──────────────────────────────────────────


def resolve_anchor(record: LoadedAlleleRecord,
                   *,
                   segment: Segment,
                   locus: Locus = Locus.UNKNOWN,
                   allow_trp_v: bool = True,
                   strict: bool = False) -> AnchorResult:
    """Resolve the V or J anchor for an allele record.

    Args:
        record: the LoadedAlleleRecord to resolve. Either freshly
            yielded by a ReferenceLoader or hand-constructed.
        segment: ``Segment.V`` or ``Segment.J`` — selects the resolver
            entry point.
        locus: locus hint. Required for canonical-confidence J results
            (IGH/TRB/TRD expect Trp, others expect Phe). When unset, J
            resolution returns ``CONF_BEST_GUESS`` and the residue
            field reveals what was found.
        allow_trp_v: when ``True`` (the default), Trp at a V anchor is
            accepted as ``CONF_ALTERNATIVE`` (biology — some TCRG and
            rare V genes use Trp). Set ``False`` to reject Trp-anchor
            V alleles outright.
        strict: when ``True``, motif-fallback ``CONF_BEST_GUESS``
            results are downgraded to ``CONF_REJECTED`` — useful for
            production builds where every allele must have an
            authoritative anchor.

    Returns:
        AnchorResult with the position, codon, residue, confidence,
        method, and reason (if rejected).
    """
    cfg = _CAnchorResolverConfig()
    _lib.anchor_resolver_config_init(ctypes.byref(cfg))
    cfg.locus = int(locus)
    cfg.allow_trp_v = bool(allow_trp_v)
    cfg.strict = bool(strict)

    c_rec = record._to_c()

    if segment == Segment.V:
        c_result = _lib.anchor_resolve_v(ctypes.byref(cfg), ctypes.byref(c_rec))
    elif segment == Segment.J:
        c_result = _lib.anchor_resolve_j(ctypes.byref(cfg), ctypes.byref(c_rec))
    else:
        raise ValueError(
            f"resolve_anchor: segment must be V or J, got {segment!r}")

    return AnchorResult._from_c(c_result)


# ── Helper: locus / segment from name (ctypes-thin) ────────────────


def locus_from_gene_name(name: str) -> Locus:
    """``IGHV1-2*01`` → ``Locus.IGH``. Returns ``Locus.UNKNOWN`` on no
    match. Backed by the C ``locus_from_gene_name``."""
    return Locus(int(_lib.locus_from_gene_name(name.encode("ascii"))))


def segment_from_gene_name(name: str) -> Segment:
    """``IGHV1-2*01`` → ``Segment.V``."""
    return Segment(int(_lib.segment_from_gene_name(name.encode("ascii"))))


def locus_from_filename(path: str) -> Locus:
    """Heuristic locus inference from the basename of ``path``."""
    return Locus(int(_lib.locus_from_filename(path.encode("utf-8"))))


__all__ = [
    "Locus",
    "Segment",
    "FunctionalStatus",
    "AnchorResidue",
    "AnchorConfidence",
    "AnchorMethod",
    "AnchorResult",
    "LoadedAlleleRecord",
    "ReferenceLoader",
    "resolve_anchor",
    "locus_from_gene_name",
    "segment_from_gene_name",
    "locus_from_filename",
]
