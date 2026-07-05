"""Tests for `tools/build_imgt_configs.py` (the IMGT cartridge build tool).

No network: `build_cartridge` is exercised with inline FASTA fixtures, so
the download/probe path is not touched. Verifies the port to
`ReferenceCartridgeBuilder` produces valid, compilable structural
cartridges and that the CLI refuses to run without an explicit
`--output-dir`.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

import GenAIRR as ga
from GenAIRR.dataconfig.enums import ChainType, Species

_REPO_ROOT = Path(__file__).resolve().parent.parent
_TOOL_PATH = _REPO_ROOT / "tools" / "build_imgt_configs.py"


def _load_tool():
    """Load the standalone tool module (it lives outside the package)."""
    spec = importlib.util.spec_from_file_location("build_imgt_configs", _TOOL_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


tool = _load_tool()

# Tiny synthetic FASTA — enough length for the recombination passes; the
# builder degrades gracefully without native anchors (compile under
# allow_curatable_refdata, per the v1 builder workflow).
_V_FASTA = (
    ">IGHV1-MOCK*01\n"
    "GAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTACAGCCTGGGGGGTCCCTGAGACTC"
    "TCCTGTGCAGCCTCTGGATTCACCTTCAGTAGCTATGCCATGAGCTGGGTCCGCCAGGCT\n"
    ">IGHV2-MOCK*01\n"
    "CAGGTCAACTTAAGGGAGTCTGGTCCTGCGCTGGTGAAACCCACACAGACCCTCACACTG"
    "ACCTGCACCTTCTCTGGGTTCTCACTCAGCACTAGTGGAATGTGTGTGAGCTGGATCCGT\n"
)
_D_FASTA = (
    ">IGHD1-MOCK*01\n"
    "GGGTATAGCAGCAGCTGGTAC\n"
    ">IGHD2-MOCK*01\n"
    "AGGATATTGTAGTGGTGGTAGCTGCTACTCC\n"
)
_J_FASTA = (
    ">IGHJ1-MOCK*01\n"
    "TACTACTACGGTATGGACGTCTGGGGCCAAGGGACCACGGTCACCGTCTCCTCAG\n"
    ">IGHJ2-MOCK*01\n"
    "ACTACTGGTACTTCGATCTCTGGGGCCGTGGCACCCTGGTCACTGTCTCCTCAG\n"
)


def test_tool_module_imports_and_carries_imgt_maps() -> None:
    """The module loads and exposes the IMGT layout tables."""
    assert tool.IMGT_BASE.startswith("https://")
    assert "Homo_sapiens" in tool.SPECIES_MAP
    assert tool.SPECIES_MAP["Homo_sapiens"][0] is Species.HUMAN
    assert set(tool.LOCUS_DEFS) == {"IGH", "IGK", "IGL", "TRA", "TRB", "TRG", "TRD"}


def test_build_cartridge_vdj_produces_valid_dataconfig() -> None:
    """A VDJ (heavy) build yields a structural DataConfig with the
    expected catalogue and identity."""
    cfg = tool.build_cartridge(
        Species.HUMAN, ChainType.BCR_HEAVY, "IGH", "MOCK",
        v_path=_V_FASTA, j_path=_J_FASTA, d_path=_D_FASTA,
    )
    assert isinstance(cfg, ga.DataConfig)
    assert cfg.name == "MOCK_IGH_IMGT"
    assert len(cfg.v_alleles) == 2
    assert len(cfg.d_alleles) == 2
    assert len(cfg.j_alleles) == 2
    assert cfg.metadata is not None
    assert cfg.metadata.has_d is True
    # Structural only: no data-derived typed empirical models.
    assert cfg.reference_models is None


def test_build_cartridge_vj_light_chain_has_no_d() -> None:
    """A VJ (light) build passes no D FASTA and yields a D-less config."""
    cfg = tool.build_cartridge(
        Species.HUMAN, ChainType.BCR_LIGHT_KAPPA, "IGK", "MOCK",
        v_path=_V_FASTA, j_path=_J_FASTA, d_path=None,
    )
    assert isinstance(cfg, ga.DataConfig)
    assert cfg.name == "MOCK_IGK_IMGT"
    assert cfg.metadata.has_d is False
    assert not cfg.d_alleles


def test_build_cartridge_missing_d_for_d_chain_returns_none() -> None:
    """A D-bearing chain with no D FASTA is skipped (returns None), not
    an exception — matches the tool's per-locus fault tolerance."""
    cfg = tool.build_cartridge(
        Species.HUMAN, ChainType.BCR_HEAVY, "IGH", "MOCK",
        v_path=_V_FASTA, j_path=_J_FASTA, d_path=None,
    )
    assert cfg is None


def test_built_cartridge_compiles_and_runs_through_experiment() -> None:
    """The built structural cartridge is drop-in for Experiment.on(cfg)."""
    cfg = tool.build_cartridge(
        Species.HUMAN, ChainType.BCR_HEAVY, "IGH", "MOCK",
        v_path=_V_FASTA, j_path=_J_FASTA, d_path=_D_FASTA,
    )
    result = (
        ga.Experiment.on(cfg)
        .recombine()
        .allow_curatable_refdata()
        .run_records(n=3, seed=1)
    )
    assert len(result) == 3


def test_cli_requires_output_dir() -> None:
    """--output-dir is mandatory; the parser rejects its absence."""
    with pytest.raises(SystemExit):
        tool._parser().parse_args(["--dry-run"])
    ns = tool._parser().parse_args(["--output-dir", "/tmp/x", "--dry-run"])
    assert ns.output_dir == "/tmp/x"
    assert ns.dry_run is True
