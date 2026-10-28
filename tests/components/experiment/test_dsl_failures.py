"""DSL misuse that must fail loudly rather than silently.
"""
from __future__ import annotations
import warnings
import pytest
from GenAIRR import CompiledExperiment, Experiment
import GenAIRR as ga  # noqa: E402
from tests.support.fixtures.experiment_dsl import _vdj_refdata, _vj_refdata

pytestmark = pytest.mark.behavior

def test_recombine_rejects_np2_lengths_on_vj_chain():
    """VJ chains have no NP2 region; passing np2_lengths was previously
    silently ignored, which hid genuine user mistakes (binding the wrong
    refdata, copy-pasted VDJ snippet, etc.). Phase 1 of the v2.0 DSL
    cleanup turns this into an explicit ValueError."""
    exp = Experiment.on(_vj_refdata()).allow_curatable_refdata()
    with pytest.raises(ValueError, match="np2_lengths.*VDJ"):
        exp.recombine(np2_lengths=[(0, 1.0), (3, 1.0)])


def test_recombine_accepts_np2_lengths_on_vdj_chain():
    """The np2_lengths kwarg is still valid for VDJ; the new check
    only fires on VJ."""
    exp = Experiment.on(_vdj_refdata()).allow_curatable_refdata()
    # Should not raise.
    compiled = exp.recombine(np2_lengths=[(0, 1.0), (3, 1.0)]).trim(enabled=False).compile(allow_curatable_refdata=True)
    assert compiled.refdata.chain_type == "vdj"


def test_recombine_accepts_no_np2_on_vj_silently():
    """Not passing np2_lengths on a VJ chain is fine — that's the
    normal case. We only complain when the user explicitly supplies a
    nonsensical kwarg."""
    exp = Experiment.on(_vj_refdata()).allow_curatable_refdata()
    # Should not raise and should not warn about np2.
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        # Trim warning is expected on raw RefData; suppress that
        # specifically so we only catch unexpected ones.
        warnings.filterwarnings("ignore", message=".*trim.*", category=UserWarning)
        warnings.filterwarnings("ignore", message=".*NP-length.*", category=UserWarning)
        exp.recombine().compile(allow_curatable_refdata=True)


def test_recombine_warns_on_raw_refdata_uniform_np_fallback():
    """Raw RefDataConfig has no DataConfig, so NP-length distributions
    fall back to a uniform synthetic default. That used to be silent;
    now we surface it so users notice they're not getting empirical
    biology."""
    exp = Experiment.on(_vj_refdata()).allow_curatable_refdata()
    with pytest.warns(UserWarning, match="uniform"):
        exp.recombine().trim(enabled=False)


def test_recombine_warns_on_raw_refdata_trim_noop():
    """trim is default-on but on raw RefDataConfig there's no DataConfig
    with trim distributions, so it's a silent no-op. v2.0 surfaces a
    warning at compile() time when the user hasn't explicitly called
    :meth:`Experiment.trim`. Calling :meth:`trim` (with any setting)
    is the documented way to suppress the warning."""
    exp = Experiment.on(_vj_refdata()).allow_curatable_refdata().recombine(
        np1_lengths=[(0, 1.0), (3, 1.0)]
    )
    with pytest.warns(UserWarning, match="trim"):
        exp.compile(allow_curatable_refdata=True)


def test_recombine_warn_silenced_by_explicit_trim_disable():
    """Calling .trim(enabled=False) explicitly silences the raw-RefDataConfig
    trim-noop warning at compile() time."""
    exp = (
        Experiment.on(_vj_refdata()).allow_curatable_refdata()
        .recombine(np1_lengths=[(0, 1.0), (3, 1.0)])
        .trim(enabled=False)
    )
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        warnings.filterwarnings("ignore", message=".*NP-length.*", category=UserWarning)
        exp.compile(allow_curatable_refdata=True)


def test_recombine_warning_silenced_when_caller_supplies_lengths_and_trim_false():
    """User passes everything explicitly → no warnings."""
    exp = Experiment.on(_vdj_refdata()).allow_curatable_refdata()
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        exp.recombine(
            np1_lengths=[(0, 1.0), (3, 1.0)],
            np2_lengths=[(0, 1.0), (3, 1.0)],
        ).trim(enabled=False)


def test_recombine_no_warnings_on_dataconfig_path():
    """Bound to a real DataConfig (via the builtin string alias) →
    empirical NP + empirical trim available → no fallback warnings."""
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        ga.Experiment.on("human_igh").recombine()
