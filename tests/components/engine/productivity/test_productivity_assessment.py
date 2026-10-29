"""Observed-span policy regressions independent of production translation helpers."""
import csv
import hashlib
import json

import pytest
import GenAIRR as ga
from GenAIRR import _engine as ge

pytestmark = pytest.mark.behavior


def refdata(v=b"AAACCCTGT", j=b"TTTGGG", *, v_anchor=6, j_anchor=0):
    cfg = ge.RefDataConfig.vj()
    cfg.add_v_allele("V*01", "V", v, anchor=v_anchor)
    cfg.add_j_allele("J*01", "J", j, anchor=j_anchor)
    cfg.set_allowed_bases(bases=list("ACGTRYSWKMBDHVN"))
    return cfg


def experiment(cfg, n=0, original_ends=False):
    exp = ga.Experiment.on(cfg).allow_curatable_refdata().recombine(np1_lengths=[(n, 1)])
    if original_ends:
        return exp.process_coding_ends(reference_ends_are_coding_ends=True)
    return exp.trim(enabled=False)


@pytest.mark.parametrize("original_ends", [False, True])
def test_reference_admission_override_does_not_change_productivity_policy(original_ends):
    cfg = refdata(v=b"AAACCCGGG")
    cfg.set_v_anchor_rule(expected_aa=["G"])
    assert cfg.validate() == []
    result = experiment(cfg, original_ends=original_ends).run_records(n=1, seed=3)
    rec = result[0]
    assert rec["productive"] is None
    assert rec["stop_codon"] is None and rec["vj_in_frame"] is None
    assert rec["junction"] == rec["junction_aa"] == ""
    assert rec["productivity_policy"] == "canonical_observed_vdj_v1"
    assert rec["productivity_scope"] == "observed_coding_span"
    assert "reference_anchor_unresolved" in rec["productivity_reasons"]
    assert result.validate_records(cfg).ok
    for strict in (False, True):
        with pytest.raises(ValueError):
            experiment(cfg, original_ends=original_ends).productive_only().run_records(n=1, seed=0, strict=strict)


@pytest.mark.parametrize("v,productive,stop", [
    (b"AAANNNTGT", None, None),
    (b"AAAGCNTGT", True, False),
    (b"AAATARTGT", False, True),
    (b"TAANNNTGT", False, True),
])
def test_iupac_resolutions_and_false_over_unknown_precedence(v, productive, stop):
    cfg = refdata(v=v)
    result = experiment(cfg).run_records(n=1, seed=2)
    assert result[0]["productive"] is productive
    assert result[0]["stop_codon"] is stop
    assert result.validate_records(cfg).ok


def test_independent_phase_does_not_invent_an_unknown_reference_landmark():
    v = b"AAATAAGGG"
    cfg = refdata(v=v)
    evidence = dict(schema_version=1, sequence_sha256="sha256:" + hashlib.sha256(v).hexdigest(),
                    coding_phase=0, anchor_status="invalid", boundary_5="unknown", boundary_3="unknown",
                    source="explicit test annotation", method="fixture", notes=[])
    cfg.set_allele_evidence("V", 0, json.dumps(evidence))
    result = experiment(cfg).run_records(n=1, seed=2)
    rec = result[0]
    assert rec["productive"] is False  # definite stop dominates unresolved landmark
    assert rec["stop_codon"] is True and rec["vj_in_frame"] is None
    assert rec["junction_start"] is None
    assert result.validate_records(cfg).ok


@pytest.mark.parametrize("loss,expected", [(1, True), (9, None), (100, None)])
def test_partial_observation_is_assessed_on_available_coding_evidence(loss, expected):
    cfg = refdata()
    result = experiment(cfg).end_loss_5prime(length=loss).run_records(n=1, seed=4)
    assert result[0]["productive"] is expected
    assert result.validate_records(cfg).ok


@pytest.mark.parametrize("original_ends", [False, True])
def test_ambiguous_boundary_conditions_np_support_before_drawing(original_ends):
    # The crossing codon is T + N-addition + N (reference ambiguity).
    # TCN/TTN are definitely sense; TAN/TGN can hide a stop.
    cfg = refdata(v=b"AAACCCTGTT", j=b"NTTTGGG", j_anchor=1)
    result = experiment(cfg, 1, original_ends).productive_only().run_records(n=40, seed=10)
    assert {r["np1"] for r in result} == {"C", "T"}
    assert all(r["productive"] is True and r["stop_codon"] is False for r in result)
    assert result.validate_records(cfg).ok


@pytest.mark.parametrize("original_ends", [False, True])
@pytest.mark.parametrize("strict", [False, True])
def test_definitely_canonical_ambiguous_landmarks_compile_generate_and_replay(original_ends, strict):
    # TGY has only Cys realizations; TTY has only Phe realizations. Neither
    # ambiguity changes the accepted class or admits a possible stop.
    cfg = refdata(v=b"AAACCCTGY", j=b"TTYGGG")
    compiled = experiment(cfg, original_ends=original_ends).productive_only().compile()
    original = compiled.simulator.run(seed=17, strict=strict)
    rec = original.airr_record(cfg)
    assert rec["productive"] is True
    assert rec["stop_codon"] is False and rec["vj_in_frame"] is True
    assert rec["junction"] == "TGYTTY"
    assert original.validate_record(cfg) == []
    trace = compiled.simulator.trace_file_from(original, seed=17)
    replayed = compiled.simulator.replay_from_trace_file(trace, strict=strict)
    assert replayed.airr_record(cfg) == rec
    assert replayed.validate_record(cfg) == []


def test_nullable_flags_diagnostics_export_and_tamper_detection(tmp_path):
    cfg = refdata(v=b"AAANNNTGT")
    result = experiment(cfg).run_records(n=1, seed=0)
    path = tmp_path / "assessment.tsv"
    result.to_tsv(str(path))
    row = next(csv.DictReader(path.open(), delimiter="\t"))
    assert row["productive"] == row["stop_codon"] == ""
    assert json.loads(row["productivity_reasons"]) == result[0]["productivity_reasons"]
    result[0]["stop_codon"] = False
    report = result.validate_records(cfg)
    assert not report.ok
    assert any(i.get("field") == "stop_codon" for f in report.failures for i in f["issues"])


@pytest.mark.parametrize("field", ["productivity_policy", "productivity_scope", "productivity_reasons"])
def test_assessment_extensions_are_reserved(field):
    with pytest.raises(ValueError, match="generated fields"):
        ga.Experiment.on("human_igh").with_metadata(**{field: "spoofed"})


@pytest.mark.parametrize("strict", [False, True])
@pytest.mark.parametrize("read_loss", [False, True])
def test_lineage_divisions_inherit_productive_constraints(strict, read_loss):
    exp = (ga.Experiment.on("human_igh").recombine().productive_only()
           .clonal_lineage(n_clones=1, max_generations=2, n_max=50, n_sample=20, rate=.15))
    if read_loss:
        exp.end_loss_5prime(length=1)
    result = exp.run_records(seed=9, strict=strict)
    assert result.records
    assert all(r["productive"] is True for r in result)
    assert any(r["n_mutations"] > 0 for r in result)
    assert result.validate_records(exp.refdata).ok


@pytest.mark.parametrize("alias, config_name", [("human_igh", "HUMAN_IGH_OGRDB"), ("human_igl", "HUMAN_IGL_OGRDB")])
def test_junction_spans_whatever_trimming_left_of_each_conserved_codon(alias, config_name):
    """The junction runs from the first surviving base of the V Cys codon to the last surviving base of
    the J Trp/Phe codon. Expected bounds come from each record's own germline and sequence coordinates and
    the reference anchors; trimming into a codon narrows the junction, and only a codon trimmed away
    entirely leaves it undefined."""
    from GenAIRR import data
    config = getattr(data, config_name)
    anchors = {(segment, allele.name): allele.anchor
               for segment, pool in (("V", config.v_alleles), ("J", config.j_alleles))
               for alleles in pool.values() for allele in alleles}
    result = ga.Experiment.on(alias).recombine().run_records(n=5000, seed=11, validate_records=True)
    shapes = {"intact": 0, "partly trimmed": 0, "codon lost": 0}
    for rec in result:
        v, j = anchors[("V", rec["v_call"].split(",")[0])], anchors[("J", rec["j_call"].split(",")[0])]
        v_kept = [g for g in range(v, v + 3) if rec["v_germline_start"] <= g < rec["v_germline_end"]]
        j_kept = [g for g in range(j, j + 3) if rec["j_germline_start"] <= g < rec["j_germline_end"]]
        if not (v_kept and j_kept):
            shapes["codon lost"] += 1
            assert (rec["junction_start"], rec["junction_end"], rec["junction"]) == (None, None, ""), rec["sequence_id"]
            continue
        shapes["partly trimmed" if len(v_kept) + len(j_kept) < 6 else "intact"] += 1
        start = rec["v_sequence_start"] + v_kept[0] - rec["v_germline_start"]
        end = rec["j_sequence_start"] + j_kept[-1] - rec["j_germline_start"] + 1
        assert (rec["junction_start"], rec["junction_end"]) == (start, end), rec["sequence_id"]
        assert rec["junction"] == rec["sequence"][start:end]
    assert shapes["partly trimmed"] > 0 and shapes["codon lost"] > 0, shapes
