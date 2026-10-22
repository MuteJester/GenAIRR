"""Running experiments on a portable ``Reference`` and binding results to it."""
import builtins
import hashlib
import io
import json
import pathlib
import warnings

import pytest

import GenAIRR as ga
from GenAIRR import Experiment, Reference, ReferenceArtifactError, SimulationResult
from GenAIRR.experimental import FamilyReplayError

# A human IGHV1-2*02-like V (TGT Cys codon at 285), two IGHD-like Ds and two IGHJ-like Js
# (TGG-GGC-CAG/A-GG motif at 14 and 29).  Software fixtures, not evidence.
V1 = ("CAGGTGCAGCTGGTGCAGTCTGGGGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCTGGATACACCTTCACC"
      "GGCTACTATATGCACTGGGTGCGACAGGCCCCTGGACAAGGGCTTGAGTGGATGGGATGGATCAACCCTAACAGTGGTGGCACAAACTAT"
      "GCACAGAAGTTTCAGGGCAGGGTCACCATGACCAGGGACACGTCCATCAGCACAGCCTACATGGAGCTGAGCAGGCTGAGATCTGACGAC"
      "ACGGCCGTGTATTACTGTGCGAGAGA")
V2 = "G" + V1[1:]
D1 = "GTATTACTATGGTTCGGGGAGTTATTATAAC"
D2 = "AGGATATTGTAGTAGTACCAGCTGCTATGCC"
J1 = "ACTACTTTGACTACTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAG"
J2 = "ATTACTACTACTACTACGGTATGGACGTCTGGGGCCAAGGGACCACGGTCACCGTCTCCTCAG"

LINEAGE = dict(n_clones=3, max_generations=6, n_max=200, n_sample=20, rate=0.05, lambda_base=1.6)
DEFERRED = "not yet supported for portable references"


def _digest(seq):
    return "sha256:" + hashlib.sha256(seq.encode("ascii")).hexdigest()


def _assertion(value):
    return {"value": value, "support": "assertion", "source": None, "method": "test fixture"}


def _entry(segment, label, seq, **annotations):
    base = {"segment": segment, "label": label, "sequence_sha256": _digest(seq), "gene": None,
            "functionality": None, "anchor": None, "coding_phase": None,
            "boundary_5": None, "boundary_3": None, "notes": []}
    base.update(annotations)
    return base


def build_reference(directory, *, name="lab_igh", version="1", v2_anchor=True):
    directory.mkdir(parents=True, exist_ok=True)
    files = {"v": f">IGHV1-2*02\n{V1}\n>IGHV1-3*01\n{V2}\n",
             "d": f">IGHD3-10*01\n{D1}\n>IGHD2-2*01\n{D2}\n",
             "j": f">IGHJ4*02\n{J1}\n>IGHJ6*02\n{J2}\n"}
    paths = {}
    for segment, text in files.items():
        paths[segment] = directory / f"{segment}.fasta"
        paths[segment].write_bytes(text.encode("utf-8"))
    sidecar = {
        "schema": "genairr_reference_annotations", "schema_version": 1,
        "reference": {"name": name, "version": version, "locus": "IGH", "chain_type": "vdj",
                      "taxon": {"scientific_name": "Homo sapiens", "taxonomy_id": 9606},
                      "build_date": None},
        "conventions": {"coordinates": "zero_based_half_open", "v_alignment": "none"},
        "alleles": [_entry("V", "IGHV1-2*02", V1, anchor=_assertion(285)),
                    _entry("V", "IGHV1-3*01", V2, anchor=_assertion(285) if v2_anchor else None),
                    _entry("J", "IGHJ4*02", J1, anchor=_assertion(14)),
                    _entry("J", "IGHJ6*02", J2, anchor=_assertion(29))],
        "relationships": [], "suggestions": [],
    }
    annotations = directory / "annotations.json"
    annotations.write_bytes(json.dumps(sidecar).encode("utf-8"))
    return Reference.from_fasta(v=paths["v"], d=paths["d"], j=paths["j"], annotations=annotations)


@pytest.fixture
def ref(tmp_path):
    return build_reference(tmp_path / "ref")


@pytest.fixture
def other_ref(tmp_path):
    return build_reference(tmp_path / "other", version="2")


def _lineage(e):
    return e.recombine().clonal_lineage(**LINEAGE)


def _expand(e):
    with pytest.warns(DeprecationWarning):
        return e.recombine().expand_clones(n_clones=2, per_clone=3)


BUILDERS = [
    pytest.param(lambda e: e.recombine(), {"n": 20}, id="plain"),
    pytest.param(_expand, {}, id="expand_clones"),
    pytest.param(lambda e: e.recombine().clonal_repertoire(n_clones=5, max_size=10), {}, id="clonal_repertoire"),
    pytest.param(_lineage, {}, id="clonal_lineage"),
]


def test_run_records_binds_hash_without_warnings(ref):
    exp = Experiment.on(ref)
    assert exp._reference is ref
    assert exp._refdata is ref.native_refdata()
    assert exp._dataconfig is None
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        result = exp.recombine().run_records(n=50, seed=1, validate_records=True)
    assert len(result) == 50
    assert result.reference_content_hash == ref.execution_hash


@pytest.mark.parametrize("builder,run_kwargs", BUILDERS)
def test_every_compiled_path_binds_reference(ref, builder, run_kwargs):
    result = builder(Experiment.on(ref)).run_records(seed=3, **run_kwargs)
    assert len(result) > 0
    assert result.reference_content_hash == ref.execution_hash
    assert result.validate_records(ref).ok


@pytest.mark.parametrize("builder,run_kwargs", BUILDERS)
def test_every_compiled_path_binds_reference_when_run_directly(ref, builder, run_kwargs):
    compiled = builder(Experiment.on(ref)).compile()
    assert compiled.refdata is ref.native_refdata()
    result = compiled.run_records(seed=4, **run_kwargs)
    assert result.reference_content_hash == ref.execution_hash
    assert result.validate_records(ref).ok


def test_compiled_experiment_run_records_binds_reference(ref):
    result = Experiment.on(ref).recombine().compile().run_records(n=10, seed=2)
    assert len(result) == 10
    assert result.reference_content_hash == ref.execution_hash
    assert result.validate_records(ref).ok


def test_generating_calls_are_view_labels(ref):
    labels = {allele.label for allele in ref.alleles}
    result = Experiment.on(ref).recombine().run_records(n=60, seed=5)
    seen = set()
    for record in result:
        assert record["v_call"] and record["j_call"]
        for field in ("v_call", "d_call", "j_call"):
            for name in record[field].split(","):
                if name:
                    assert name in labels, (field, name)
                    seen.add(name)
    assert {"IGHV1-2*02", "IGHV1-3*01", "IGHJ4*02", "IGHJ6*02"} <= seen


def test_validate_records_accepts_reference_and_rejects_other(ref, other_ref):
    result = Experiment.on(ref).recombine().run_records(n=10, seed=6)
    assert result.validate_records(ref).ok
    assert result.validate_records(ref.native_refdata()).ok
    assert other_ref.execution_hash != ref.execution_hash
    for wrong in (other_ref, other_ref.native_refdata(), ref.native_refdata().curated("raw")):
        with pytest.raises(ValueError, match="records were generated against a different reference"):
            result.validate_records(wrong)


def test_validate_families_with_parents_accepts_reference_and_rejects_other(ref, other_ref):
    result = _expand(Experiment.on(ref)).run_records(seed=8, expose_provenance=True)
    assert result.parents is not None
    assert result.validate_families_with_parents(ref).ok
    with pytest.raises(ValueError, match="records were generated against a different reference"):
        result.validate_families_with_parents(other_ref)


def test_builtin_result_is_bound_and_rejects_a_reference(ref):
    exp = Experiment.on("human_igh").recombine()
    result = exp.run_records(n=3, seed=1)
    assert result.reference_content_hash == exp.refdata.content_hash()
    with pytest.raises(ValueError, match="records were generated against a different reference"):
        result.validate_records(ref)


def test_records_only_results_are_unbound():
    assert SimulationResult([{"sequence_id": "a"}]).reference_content_hash is None


def test_explicit_defaults_drive_np_and_trim(ref):
    result = Experiment.on(ref).recombine().run_records(n=200, seed=9)
    for record in result:
        assert 0 <= record["np1_length"] <= 6
        assert 0 <= record["np2_length"] <= 6
        for field in ("v_trim_3", "d_trim_5", "d_trim_3", "j_trim_5"):
            assert record[field] == 0, field
    assert len({record["np1_length"] for record in result}) > 1

    # Each default can still be overridden per run.
    custom = (Experiment.on(ref).recombine(np1_lengths=[(9, 1.0)], np2_lengths=[(8, 1.0)])
              .trim(v_3=[(2, 1.0)]).run_records(n=20, seed=9))
    assert {(r["np1_length"], r["np2_length"], r["v_trim_3"]) for r in custom} == {(9, 8, 2)}
    assert custom.reference_content_hash == ref.execution_hash


def test_curate_refdata_on_reference_rejected_at_the_call(ref):
    exp = Experiment.on(ref).recombine()
    for policy in ("raw", "functional_anchors_only", "functional_status"):
        with pytest.raises(ReferenceArtifactError) as excinfo:
            exp.curate_refdata(policy)
        assert excinfo.value.kind == "capability"
        message = str(excinfo.value)
        assert "cannot be curated in place" in message
        assert "allele selection is not yet supported for portable references" in message
    assert exp.refdata is ref.native_refdata()
    assert exp.run_records(n=1, seed=0).reference_content_hash == ref.execution_hash


def test_replaced_native_reference_rejected_at_compile(ref):
    exp = Experiment.on(ref).recombine()
    exp._refdata = ref.native_refdata().curated("raw")
    with pytest.raises(ReferenceArtifactError) as excinfo:
        exp.compile()
    assert excinfo.value.kind == "integrity"
    assert ref.execution_hash in str(excinfo.value)
    with pytest.raises(ReferenceArtifactError):
        exp.run_records(n=1, seed=0)


def test_instance_curatable_flag_overrides_the_admission_default(tmp_path):
    admitted = build_reference(tmp_path / "admitted", v2_anchor=False)
    assert "MissingAnchor" in {f["kind"] for f in admitted.admission["accepted_curatable_findings"]}
    # Not set: the reference's admitted findings select the lenient mode.
    Experiment.on(admitted).recombine().compile()
    # Explicitly set to False: strict validation rejects the admitted finding.
    strict = Experiment.on(admitted).allow_curatable_refdata(False).recombine()
    with pytest.raises(ValueError, match="V allele 1 has no anchor"):
        strict.compile()
    with pytest.raises(ValueError, match="V allele 1 has no anchor"):
        strict.run_records(n=1, seed=0)
    # The compile() argument still wins over both.
    strict.compile(allow_curatable_refdata=True)
    with pytest.raises(ValueError, match="V allele 1 has no anchor"):
        Experiment.on(admitted).recombine().compile(allow_curatable_refdata=False)
    Experiment.on(admitted).allow_curatable_refdata(True).recombine().compile()


def test_genotype_and_cohort_rejected_with_capability_error(ref):
    genotype = ga.Genotype.from_dataconfig(ga.data.HUMAN_IGH_OGRDB)

    with pytest.raises(ReferenceArtifactError, match=DEFERRED) as excinfo:
        Experiment.on(ref).recombine().with_genotype(genotype)
    assert excinfo.value.kind == "capability"

    with pytest.raises(ReferenceArtifactError, match=DEFERRED) as excinfo:
        Experiment.on(ref).recombine().run_cohort([genotype])
    assert excinfo.value.kind == "capability"

    # A genotype that reaches compilation by any route is still rejected there.
    exp = Experiment.on(ref).recombine()
    exp._genotype = genotype._snapshot()
    with pytest.raises(ReferenceArtifactError, match=DEFERRED) as excinfo:
        exp.compile()
    assert excinfo.value.kind == "capability"


def test_family_replay_capture_with_reference_unsupported(ref, monkeypatch):
    import GenAIRR._s5f_loader as loader

    monkeypatch.setattr(loader, "load_builtin_s5f_kernel",
                        lambda *a, **k: pytest.fail("sample preparation reached"))
    with pytest.raises(FamilyReplayError, match="reference_spec_v1") as excinfo:
        _lineage(Experiment.on(ref)).run_records(seed=7, capture_family_replay=True)
    assert DEFERRED in str(excinfo.value)
    with pytest.raises(FamilyReplayError, match="reference_spec_v1") as excinfo:
        _lineage(Experiment.on(ref)).compile().run_records(seed=7, capture_family_replay=True)
    assert DEFERRED in str(excinfo.value)


def test_family_replay_restored_result_binds_reconstructed_reference():
    from GenAIRR._family_replay import _bundle_from_result
    from tests.duplex_experiment_cases import history, model, references

    refs = references()
    exp = (Experiment.on(refs).recombine(duplex_repair=model(refs, [history(restore=True)]))
           .clonal_lineage(n_clones=1, max_generations=2, n_max=8, n_sample=6, rate=.05, lambda_base=3))
    result = exp.run_records(seed=7, capture_family_replay=True)
    assert result.reference_content_hash == refs.content_hash()
    bundle = _bundle_from_result(result)
    restored = bundle.replay_records()
    assert restored.reference_content_hash == bundle.refdata.content_hash() == refs.content_hash()
    assert restored.validate_records(bundle.refdata).ok


def test_experiment_on_path_string_points_to_reference_read(ref, tmp_path, monkeypatch):
    path = tmp_path / "lab.genairr-ref.json"
    ref.write(path)

    def no_reading(*args, **kwargs):
        raise AssertionError("the artifact file was read")

    monkeypatch.setattr(builtins, "open", no_reading)
    monkeypatch.setattr(io, "open", no_reading)
    monkeypatch.setattr(pathlib.Path, "read_text", no_reading)
    monkeypatch.setattr(pathlib.Path, "read_bytes", no_reading)
    monkeypatch.setattr(pathlib.Path, "open", no_reading)
    for source in (str(path), "missing.genairr-ref.json", "x.json"):
        with pytest.raises(TypeError, match=r"GenAIRR\.Reference\.read\(path\)"):
            Experiment.on(source)
    with pytest.raises(TypeError, match=r"GenAIRR\.Reference\.read\(path\)"):
        Experiment.on(path)
    monkeypatch.undo()
    assert Experiment.on(Reference.read(path))._reference.execution_hash == ref.execution_hash


def test_existing_inputs_unchanged_smoke():
    by_name = Experiment.on("human_igh")
    by_config = Experiment.on(ga.data.HUMAN_IGH_OGRDB)
    by_refdata = Experiment.on(by_name.refdata)
    assert by_name._reference is None and by_name._dataconfig is ga.data.HUMAN_IGH_OGRDB
    assert by_config._reference is None and by_config._dataconfig is ga.data.HUMAN_IGH_OGRDB
    assert by_refdata._reference is None and by_refdata._dataconfig is None
    assert by_refdata.refdata is by_name.refdata
    for exp in (by_name, by_config):
        result = exp.recombine().run_records(n=3, seed=1, validate_records=True)
        assert result.reference_content_hash == exp.refdata.content_hash()
    with pytest.warns(UserWarning, match="raw RefDataConfig"):
        by_refdata.recombine()
    with pytest.warns(UserWarning, match="no trim distributions"):
        result = by_refdata.run_records(n=3, seed=1, validate_records=True)
    assert result.reference_content_hash == by_refdata.refdata.content_hash()
    with pytest.raises(ValueError, match="Unknown config name"):
        Experiment.on("not_a_config")


# ── describe() ───────────────────────────────────────────────────────

BUILTIN_DESCRIBE = (
    "Experiment on Human heavy-chain BCR\n"
    "  1. V(D)J recombination: sample V/D/J alleles; custom V allele weights; empirical exonuclease "
    "trim (V3', D5', D3', J5'); insert NP1 (0–40 weighted bases, empirical) and NP2 (0–40 weighted "
    "bases, empirical)\n"
    "  2. Somatic hypermutation (S5F context model, human heavy-chain (HH_S5F)): 5–15 mutations/record"
)
RAW_DESCRIBE = ("Experiment on raw RefDataConfig (VJ)\n"
                "  1. V(D)J recombination: sample V/J alleles; insert NP1 (0–6 bases, uniform fallback)")
CLONAL_DESCRIBE = (
    "Experiment on Human heavy-chain BCR\n"
    "  1. V(D)J recombination: sample V/D/J alleles; empirical exonuclease trim (V3', D5', D3', J5'); "
    "insert NP1 (3 bases, empirical) and NP2 (0–6 bases, uniform fallback)\n"
    "  ── Clonal expansion: 2 lineages × 3 reads/clone ──\n"
    "      (steps above run once per clone; steps below run once per descendant)\n"
    "  2. Somatic hypermutation (S5F context model, human heavy-chain (HH_S5F)): 2 mutations/record"
)
PROFILE = "genairr.unfitted_defaults.v1"
REFERENCE_DESCRIBE = (
    "Experiment on portable reference 'lab_igh' (VDJ)\n"
    "  1. V(D)J recombination: sample V/D/J alleles; uniform V/D/J allele usage (unfitted defaults "
    "profile genairr.unfitted_defaults.v1); no exonuclease trim or P-nucleotides (unfitted defaults "
    "profile); insert NP1 (0–6 bases, unfitted defaults profile) and NP2 (0–6 bases, unfitted "
    "defaults profile)"
)


def test_builtin_describe_output_unchanged():
    exp = (Experiment.on("human_igh").recombine(v_allele_weights={"IGHVF1-G1*01": 5.0})
           .mutate(count=(5, 15)))
    assert exp.describe() == BUILTIN_DESCRIBE
    assert exp.compile().describe() == BUILTIN_DESCRIBE
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        raw = Experiment.on(Experiment.on("human_igk").refdata).recombine().trim(enabled=False)
        clonal = (Experiment.on(ga.data.HUMAN_IGH_OGRDB)
                  .recombine(np1_lengths=[(3, 1.0)], np2_lengths=[(i, 1.0) for i in range(7)])
                  .expand_clones(n_clones=2, per_clone=3).mutate(count=2))
    assert raw.describe() == RAW_DESCRIBE
    assert clonal.compile().describe() == CLONAL_DESCRIBE


def test_reference_describe_names_reference_and_profile(ref):
    exp = Experiment.on(ref).recombine()
    assert exp.describe() == REFERENCE_DESCRIBE
    texts = [exp.describe(), exp.compile().describe(), _expand(Experiment.on(ref)).compile().describe(),
             _expand(Experiment.on(ref)).describe()]
    for text in texts:
        assert "portable reference 'lab_igh'" in text
        assert PROFILE in text
        for wrong in ("raw RefDataConfig", "fallback", "custom"):
            assert wrong not in text, (wrong, text)


def test_reference_describe_marks_user_overrides(ref):
    exp = (Experiment.on(ref)
           .recombine(v_allele_weights={"IGHV1-2*02": 3.0}, np1_lengths=[(9, 1.0)])
           .trim(v_3=[(2, 1.0)]))
    text = exp.describe()
    assert "custom V allele weights" in text
    assert f"uniform D/J allele usage (unfitted defaults profile {PROFILE})" in text
    assert "NP1 (9 bases, user-specified)" in text
    assert "NP2 (0–6 bases, unfitted defaults profile)" in text
    assert "user-specified exonuclease trim (V3'); no P-nucleotides (unfitted defaults profile)" in text
    for wrong in ("raw RefDataConfig", "fallback", "empirical"):
        assert wrong not in text
    assert exp.compile().describe() == text

    # A user override that happens to equal the default is still the user's.
    same = Experiment.on(ref).recombine(np1_lengths=[(i, 1.0) for i in range(7)]).describe()
    assert "NP1 (0–6 bases, user-specified)" in same
    disabled = Experiment.on(ref).recombine().trim(enabled=False).describe()
    assert "no exonuclease trim (user-specified); no P-nucleotides (unfitted defaults profile)" in disabled


def test_unknown_config_name_points_to_reference_read(monkeypatch):
    def no_reading(*args, **kwargs):
        raise AssertionError("a file was read")

    monkeypatch.setattr(builtins, "open", no_reading)
    monkeypatch.setattr(io, "open", no_reading)
    monkeypatch.setattr(pathlib.Path, "read_text", no_reading)
    monkeypatch.setattr(pathlib.Path, "read_bytes", no_reading)
    monkeypatch.setattr(pathlib.Path, "open", no_reading)
    with pytest.raises(ValueError, match="Unknown config name") as excinfo:
        Experiment.on("refs/lab.genairr-ref")
    assert "To load a portable reference file, use GenAIRR.Reference.read(path)." in str(excinfo.value)
