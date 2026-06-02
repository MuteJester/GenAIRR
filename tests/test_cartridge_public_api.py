"""Public API contract for the reference cartridge.

Confirms the cartridge authoring types are reachable from the
top-level ``GenAIRR`` namespace and that a minimal custom cartridge
flows end-to-end through compile. Companion to
``docs/reference_cartridge.md`` — keep in sync.
"""
from __future__ import annotations

from pathlib import Path

import pytest


# ──────────────────────────────────────────────────────────────────
# 1. Discoverability — importable from the top-level package.
# ──────────────────────────────────────────────────────────────────


def test_cartridge_specs_are_importable_from_top_level() -> None:
    """A user pasting the README snippet must succeed without
    reaching into private submodules."""
    from GenAIRR import (  # noqa: F401
        AnchorRuleSpec,
        EmpiricalDistributionSpec,
        ReferenceEmpiricalModels,
        ReferenceRulesSpec,
    )


def test_cartridge_specs_appear_in_all() -> None:
    import GenAIRR

    for name in (
        "AnchorRuleSpec",
        "EmpiricalDistributionSpec",
        "ReferenceEmpiricalModels",
        "ReferenceRulesSpec",
    ):
        assert name in GenAIRR.__all__, f"{name!r} missing from GenAIRR.__all__"
        assert hasattr(GenAIRR, name)


def test_cartridge_docstrings_mention_cartridge_doc() -> None:
    """The authoring types' docstrings should point to the
    cartridge guide so editors and `help(...)` surfaces it."""
    from GenAIRR import (
        AnchorRuleSpec,
        ReferenceEmpiricalModels,
        ReferenceRulesSpec,
    )

    for cls in (ReferenceRulesSpec, ReferenceEmpiricalModels, AnchorRuleSpec):
        assert cls.__doc__ is not None
        assert "docs/reference_cartridge.md" in cls.__doc__, (
            f"{cls.__name__} docstring should reference docs/reference_cartridge.md"
        )


def test_experiment_curation_methods_documented() -> None:
    """``curate_refdata`` + ``allow_curatable_refdata`` should
    surface the cartridge model and the relationship between the
    two opt-ins in their docstrings."""
    from GenAIRR import Experiment

    cura = Experiment.curate_refdata.__doc__
    allow = Experiment.allow_curatable_refdata.__doc__
    assert cura is not None and "Curation" in cura
    assert cura is not None and "validation" in cura
    assert allow is not None and "curate_refdata" in allow
    assert allow is not None and "allow_curatable_refdata" in allow


# ──────────────────────────────────────────────────────────────────
# 2. End-to-end — minimal cartridge with rules + models compiles.
# ──────────────────────────────────────────────────────────────────


def _minimal_human_igh_with_specs():
    """Bundled human_igh with explicit rules + models attached.
    Built from the bundled config so we don't have to hand-author
    Allele subclasses for a smoke test.
    """
    import GenAIRR as ga
    from GenAIRR import (
        AnchorRuleSpec,
        EmpiricalDistributionSpec,
        ReferenceEmpiricalModels,
        ReferenceRulesSpec,
    )
    from GenAIRR._refdata_resolver import _resolve_config_name

    cfg = _resolve_config_name("human_igh").copy()
    cfg.reference_rules = ReferenceRulesSpec(
        v_anchor=AnchorRuleSpec(expected_aa=["C"]),
        j_anchor=AnchorRuleSpec(expected_aa=["W"]),
    )
    cfg.reference_models = ReferenceEmpiricalModels(
        np_lengths={
            "NP1": EmpiricalDistributionSpec([(0, 1.0), (3, 4.0), (6, 2.0)]),
            "NP2": EmpiricalDistributionSpec([(0, 1.0), (2, 3.0), (4, 2.0)]),
        },
    )
    return cfg, ga


def test_minimal_custom_cartridge_compiles_end_to_end() -> None:
    cfg, ga = _minimal_human_igh_with_specs()
    compiled = ga.Experiment.on(cfg).recombine().compile()
    assert compiled is not None


def test_minimal_custom_cartridge_runs_end_to_end() -> None:
    cfg, ga = _minimal_human_igh_with_specs()
    result = ga.Experiment.on(cfg).recombine().run_records(n=5, seed=42)
    assert len(result.records) == 5


# ──────────────────────────────────────────────────────────────────
# 3. Doctest-style cartridge guide snippet — keep README in sync.
# ──────────────────────────────────────────────────────────────────


def test_readme_quickstart_snippet_uses_public_imports() -> None:
    """The snippet documented in ``docs/reference_cartridge.md``
    imports everything from the top-level ``GenAIRR`` namespace.
    This regression-guards that all of those names stay exported."""
    snippet = """
from GenAIRR import (
    AnchorRuleSpec,
    EmpiricalDistributionSpec,
    ReferenceEmpiricalModels,
    ReferenceRulesSpec,
    DataConfig,
    RefDataConfig,
    Experiment,
)
"""
    exec(compile(snippet, "<README snippet>", "exec"), {})


def test_cartridge_doc_exists() -> None:
    """The cartridge guide is the canonical reference for the
    authoring surface; a missing file means the API docstrings
    point at nothing."""
    docs_dir = Path(__file__).resolve().parent.parent / "docs"
    if not docs_dir.is_dir():
        import pytest
        pytest.skip("docs/ is contributor-only; not present in this checkout")

    doc = Path(__file__).resolve().parent.parent / "docs" / "reference_cartridge.md"
    assert doc.is_file(), f"docs/reference_cartridge.md is missing at {doc}"
    text = doc.read_text(encoding="utf-8")
    # Sanity-check the four planes are all named in the doc.
    for section in ("Identity", "Catalogue", "Rules", "Empirical models", "Curation"):
        assert section in text, f"reference_cartridge.md missing section {section!r}"
