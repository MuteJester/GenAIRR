"""Literal annotation oracles independent of the native placement walker."""
import csv
import hashlib
import json
from itertools import product

import pytest

import GenAIRR as ga
from GenAIRR import RefDataConfig
from GenAIRR.experimental import CoupledDuplexRepair, DuplexJoin, export_duplex_projection


def reference(v, d, j, *, va=None, ja=None):
    r = RefDataConfig.vdj()
    r.add_v_allele("V*01", "V", v.encode(), anchor=va)
    r.add_d_allele("D*01", "D", d.encode())
    r.add_j_allele("J*01", "J", j.encode(), anchor=ja)
    return r


def evidence(r, segment, sequence, *, status="verified", phase=0, five="complete", three="complete"):
    r.set_allele_evidence(segment, 0, json.dumps({
        "schema_version": 1, "sequence_sha256": "sha256:" + hashlib.sha256(sequence.encode()).hexdigest(),
        "anchor_status": status, "coding_phase": phase, "boundary_5": five, "boundary_3": three,
        "source": "independent literal fixture", "method": "authored_coordinate", "notes": [],
    }))


def opened(name, left=0, right=0):
    return DuplexJoin(name).open("left", offset=left).open("right", offset=right)


def blunt(name, at):
    return opened(name).synapse(right_start=at).ligate()


def authored(v, d, j, *, order=("vd", "dj"), inverted=False, scripts=None):
    scripts = scripts or {"vd": blunt("vd", len(v)), "dj": blunt("dj", len(d))}
    x = CoupledDuplexRepair(v, d, j, d_inverted=inverted)
    for name in order:
        x = x.join(scripts[name])
    return x


def project(script, r=None, **kwargs):
    r = r or reference(script.v, script.d, script.j)
    return script.project(r, v_call=kwargs.pop("v_call", "V*01"), d_call="D*01", j_call="J*01",
                          reference_ends_are_coding_ends=True, **kwargs)


def assert_path(p, segment, q, g, cigar):
    path = p["paths"][segment]
    assert path["available"] is True
    assert (path["query_start"], path["query_end"]) == q
    assert (path["reference_start"], path["reference_end"]) == g
    assert path["cigar"] == cigar
    assert len(path["sequence_alignment"]) == len(path["germline_alignment"])


@pytest.mark.parametrize("order", [("vd", "dj"), ("dj", "vd")])
@pytest.mark.parametrize("inverted", [False, True])
def test_blunt_body_paths_are_exact_and_chemical_ledger_unchanged(order, inverted):
    script = authored("ACG", "AAT", "TCC", order=order, inverted=inverted)
    p = project(script)
    rec = p["record"]
    assert rec["sequence"] == ("ACGATTTCC" if inverted else "ACGAATTCC")
    assert rec["duplex_repair"]["product"] == script.run()
    for i, name in enumerate(("v", "d", "j")):
        assert_path(p, name, (3*i, 3*i+3), (0, 3), "3M")
        assert rec[f"{name}_trim_5"] == rec[f"{name}_trim_3"] == 0
    assert rec["np1"] == rec["np2"] == ""
    assert rec["np1_length"] == rec["np2_length"] == 0
    assert rec["sequence_alignment"] == rec["sequence"]
    assert rec["n_indels"] == rec["n_mutations"] == 0
    assert rec["p_v_3_length"] is None
    assert p["capabilities"]["experiment_integration"] is False


@pytest.mark.parametrize("order", [("vd", "dj"), ("dj", "vd")])
def test_restored_d_is_body_but_remains_synthesized(order):
    scripts = {"vd": opened("vd").resect("right", "top", count=1).synapse(right_start=3).fill("top").ligate(),
               "dj": blunt("dj", 3)}
    p = project(authored("ACG", "CAT", "TAA", order=order, scripts=scripts))
    assert_path(p, "d", (3, 6), (0, 3), "3M")
    rec = p["record"]
    assert rec["d_trim_5"] == 0
    truth = rec["duplex_repair"]
    top = truth["product"]["columns"][3]["top"]
    assert top["kind"] == "templated_synthesis" and top["original"] is None
    assert top["template_id"] == 7  # canonical D-bottom position 0
    assert truth["evidence"]["placements"][str(top["id"])]["role"] == "restored_body"
    assert truth["counts"]["physical"]["resected_nucleotides"] == 1
    assert truth["counts"]["physical"]["templated_synthesis"] == 1
    assert rec["n_d_indels"] == 0


def test_restored_v_landmark_has_new_policy_and_actual_junction():
    script = authored("AAATGT", "CCC", "TTTAAA", scripts={
        "vd": opened("vd").resect("left", "top", count=1).synapse(right_start=6).fill("top").ligate(),
        "dj": blunt("dj", 3)})
    p = project(script, reference(script.v, script.d, script.j, va=3, ja=0))
    rec = p["record"]
    assert rec["productive"] is True and rec["vj_in_frame"] is True and rec["stop_codon"] is False
    assert rec["junction"] == "TGTCCCTTT" and rec["junction_aa"] == "CPF"
    assert (rec["junction_start"], rec["junction_end"]) == (3, 12)
    assert rec["productivity_policy"] == "duplex_observed_vdj_v1"
    assert_path(p, "v", (0, 6), (0, 6), "6M")


def test_nonproductive_record_retains_junction_translation():
    script = authored("AAATGT", "TAA", "TTTAAA")
    rec = project(script, reference(script.v, script.d, script.j, va=3, ja=0))["record"]
    assert rec["productive"] is False and rec["stop_codon"] is True
    assert rec["junction"] == "TGTTAATTT"
    assert rec["junction_aa"] == "C*F"


def test_verified_noncanonical_landmark_is_locatable_without_phase():
    script = authored("AAACCCGGG", "CCC", "TTTAAA")
    r = reference(script.v, script.d, script.j, va=6, ja=0)
    evidence(r, "V", script.v, phase=None)
    rec = project(script, r)["record"]
    assert rec["productive"] is False
    assert rec["junction_start"] == 6 and rec["junction_end"] == 15
    assert rec["junction_aa"] == "GPF"
    assert rec["vj_in_frame"] is None and rec["stop_codon"] is None
    assert "anchor_residue_rejected" in rec["productivity_reasons"]


def test_identical_sequence_cannot_restore_a_removed_anchor_by_motif_matching():
    r = reference("AAATGT", "CCC", "TTTAAA", va=3, ja=0)
    natural = project(authored("AAATGT", "CCC", "TTTAAA"), r)["record"]
    script = authored("AAATGT", "CCC", "TTTAAA", scripts={
        "vd": opened("vd").resect("left", "top", count=3).resect("left", "bottom", count=3)
        .add_n("left", "TGT").synapse(right_start=6).fill("bottom").ligate(),
        "dj": blunt("dj", 3)})
    p = project(script, r)
    rec = p["record"]
    assert rec["sequence"] == natural["sequence"]
    assert rec["sequence_id"] != natural["sequence_id"]
    assert natural["productive"] is True and rec["productive"] is False
    assert rec["v_trim_3"] == 3 and rec["v_cigar"] == "3M"
    assert rec["np1"] == "TGT" and rec["np1_length"] == 3
    assert rec["junction"] is None and "anchor_removed" in rec["productivity_reasons"]
    assert rec["duplex_repair"]["counts"]["physical"]["untemplated_additions"] == 3


@pytest.mark.parametrize("order", [("vd", "dj"), ("dj", "vd")])
def test_paired_n_additions_count_query_columns_separately_from_physical_events(order):
    script = authored("AAA", "CCC", "GGG", order=order, scripts={
        "vd": opened("vd").add_n("left", "G").add_n("right", "C").synapse(right_start=4).ligate(),
        "dj": opened("dj").add_n("left", "A").add_n("right", "T").synapse(right_start=4).ligate()})
    p = project(script)
    rec = p["record"]
    assert rec["sequence"] == "AAAGCCCAGGG"
    assert rec["np1"] == "G" and rec["np2"] == "A"
    assert rec["np1_length"] == rec["np2_length"] == 1
    counts = rec["duplex_repair"]["counts"]
    assert counts["physical"]["untemplated_additions"] == 4
    assert counts["founder_query"]["placement_roles"]["n_addition"] == 2
    assert counts["founder_query"]["timepoint"] == "after_both_joins"


def test_cross_source_body_overlap_is_preserved_without_forced_owner():
    script = authored("ACG", "GTA", "CCC", scripts={
        "vd": opened("vd").resect("left", "bottom", count=1).resect("right", "top", count=1).synapse(right_start=2).ligate(),
        "dj": blunt("dj", 3)})
    p = project(script)
    assert p["record"]["sequence"] == "ACGTACCC"
    assert_path(p, "v", (0, 3), (0, 3), "3M")
    assert_path(p, "d", (2, 5), (0, 3), "3M")
    assert p["record"]["sequence_alignment"] is None and p["record"]["germline_alignment"] is None
    assert p["record"]["v_alignment_start"] is None
    assert p["record"]["np1_length"] == 0
    assert p["capabilities"]["simulation_seed"] is False
    overlap = p["record"]["duplex_repair"]["overlaps"][0]
    assert overlap["query_position"] == 2
    assert {(c["segment"], c["reference_position"]) for c in overlap["claims"]} == {("v", 2), ("d", 0)}


def test_top_p_pairs_with_d_body_without_relabeling_physical_origin():
    scripts = {n: opened(n, 1, 1).synapse(right_start=at).ligate() for n, at in (("vd", 3), ("dj", 2))}
    p = project(authored("AAG", "CA", "TAA", scripts=scripts))
    assert p["record"]["sequence"] == "AAGCATAA"
    assert_path(p, "d", (3, 5), (0, 2), "2M")
    truth = p["record"]["duplex_repair"]
    top = truth["product"]["columns"][3]["top"]
    assert top["original"]["segment"] == "v"
    assert truth["evidence"]["placements"][str(top["id"])]["role"] == "hairpin_transfer"
    assert p["paths"]["d"]["body"][0]["support_ids"] == [7]
    assert truth["counts"]["physical"]["hairpin_transfer_events"] == 4


def test_inverted_d_tip_events_keep_canonical_and_assembled_ends_distinct():
    scripts = {n: opened(n, 1, 1).synapse(right_start=at).ligate() for n, at in (("vd", 3), ("dj", 2))}
    p = project(authored("AAG", "CG", "CAA", inverted=True, scripts=scripts))
    events = p["record"]["duplex_repair"]["evidence"]["events"]
    tips = [event["tip"] for event in events if event["event"] == "transferred" and event["tip"]["segment"] == "d"]
    assert {(tip["end"], tip["canonical_end"]) for tip in tips} == {("five", "three"), ("three", "five")}


def test_moved_first_join_synthesis_loses_body_claim_and_removed_template_chain_survives():
    script = authored("ATG", "CAA", "TAA", order=("dj", "vd"), scripts={
        "dj": opened("dj").resect("left", "top", count=2).synapse(right_start=3).fill("top").ligate(),
        "vd": opened("vd", 2, 2).synapse(right_start=3).ligate()})
    truth = project(script)["record"]["duplex_repair"]
    moved = truth["product"]["columns"][1]["bottom"]
    placement = truth["evidence"]["placements"][str(moved["id"])]
    assert moved["kind"] == "templated_synthesis"
    assert placement["body"] is None and placement["role"] == "hairpin_transfer"
    script = authored("ACG", "CAA", "TAA", order=("dj", "vd"), scripts={
        "dj": opened("dj").resect("left", "top", count=2).synapse(right_start=3).fill("top").ligate(),
        "vd": opened("vd").resect("right", "bottom", count=2).synapse(right_start=3).fill("bottom").ligate()})
    p = project(script)
    assert_path(p, "d", (3, 6), (0, 3), "3M")
    truth = p["record"]["duplex_repair"]
    top, bottom = truth["product"]["columns"][4]["top"], truth["product"]["columns"][4]["bottom"]
    assert top["template_id"] in {b["id"] for b in truth["product"]["removed"]}
    assert bottom["template_id"] == top["id"]


def test_reference_binding_requires_explicit_allele_identity_and_boundary_contract():
    script = authored("AAATGT", "CCC", "TTTAAA")
    r = reference(script.v, script.d, script.j, va=3, ja=0)
    r.add_v_allele("V*02", "V", script.v.encode(), anchor=3)
    p = project(script, r, v_call="V*02")
    assert p["record"]["v_call"] == "V*02"
    assert p["record"]["duplex_repair"]["bindings"]["v"]["allele_id"] == 1
    with pytest.raises(ValueError, match="unknown generating"):
        project(script, r, v_call="missing")
    with pytest.raises(ValueError, match="unknown coding end"):
        script.project(r, v_call="V*01", d_call="D*01", j_call="J*01")
    evidence(r, "V", script.v, three="partial")
    with pytest.raises(ValueError, match="known partial"):
        project(script, r)


def test_declared_complete_boundaries_and_replay_are_explicit():
    script = authored("AAATGT", "CCC", "TTTAAA")
    r = reference(script.v, script.d, script.j, va=3, ja=0)
    for segment, sequence in (("V", script.v), ("D", script.d), ("J", script.j)):
        evidence(r, segment, sequence, status="unknown" if segment == "D" else "verified",
                 phase=None if segment == "D" else 0)
    p = script.project(r, v_call="V*01", d_call="D*01", j_call="J*01")
    assert p["record"]["productive"] is True
    assert p["record"]["duplex_repair"]["bindings"]["v"]["boundary_contract"] == "declared_complete"
    replay = CoupledDuplexRepair.from_dict(json.loads(json.dumps(script.to_dict())))
    assert replay.project(r, v_call="V*01", d_call="D*01", j_call="J*01") == p
    changed = reference("AAATGC", script.d, script.j, va=3, ja=0)
    with pytest.raises(ValueError, match="entire canonical"):
        project(script, changed)


@pytest.mark.parametrize("estimator", ["estimate_np_length_distributions", "estimate_np_base_model",
                                     "estimate_p_nucleotide_lengths", "estimate_trim_distributions"])
@pytest.mark.parametrize("as_tsv", [False, True])
def test_legacy_estimators_reject_duplex_provenance_before_changing_models(estimator, as_tsv, tmp_path):
    builder = ga.ReferenceCartridgeBuilder.from_fasta(
        v_fasta=">IGKV1-MOCK*01\nAAACCCTGT\n", j_fasta=">IGKJ1-MOCK*01\nTGGAAA\n",
        chain_type="BCR_LIGHT_KAPPA")
    row = project(authored("AAA", "CCC", "GGG"))["record"]
    source = [row]
    if as_tsv:
        target = tmp_path / "duplex.tsv"
        export_duplex_projection(project(authored("AAA", "CCC", "GGG")), target, format="tsv")
        source = target
    stages = list(builder.report().stages)
    models = builder._reference_models
    with pytest.raises(ValueError, match="duplex repair provenance"):
        getattr(builder, estimator)(source)
    assert builder._reference_models is models
    assert builder.report().stages == stages


def test_retained_outcome_fitting_cannot_reinterpret_duplex_history():
    from GenAIRR.coding_end_models import CodingEndOutcomeModel
    row = project(authored("AAA", "CCC", "GGG"))["record"]
    with pytest.raises(ValueError, match="history-aware fitting"):
        CodingEndOutcomeModel.fit_generating_records(
            [row], chain_type="vdj", reference_hash="literal-fixture", source="authored histories",
            training_sampling_law="unconditioned_generating_choices", independent_generating_events=True)


def test_strict_export_preserves_nested_coordinates_nulls_and_asymmetric_inverted_d(tmp_path):
    script = authored("ACG", "AATG", "CCC", inverted=True, scripts={
        "vd": blunt("vd", 3),
        "dj": opened("dj").resect("left", "top", count=1).resect("left", "bottom", count=1).synapse(right_start=3).ligate()})
    p = project(script)
    assert_path(p, "d", (3, 6), (1, 4), "3M")
    assert p["paths"]["d"]["oriented_reference_start"] == 0
    before = json.dumps(p, sort_keys=True)
    export_duplex_projection(p, tmp_path / "duplex.json")
    assert json.loads((tmp_path / "duplex.json").read_text()) == p
    export_duplex_projection(p, tmp_path / "duplex.tsv", format="tsv")
    row = next(csv.DictReader((tmp_path / "duplex.tsv").open(), delimiter="\t"))
    assert row["d_germline_start"] == "2" and row["d_germline_end"] == "4"
    assert row["d_sequence_start"] == "4" and row["d_cigar"] == "3S3M"
    assert row["p_d_5_length"] == ""
    assert json.loads(row["duplex_repair"]) == p["record"]["duplex_repair"]
    assert json.dumps(p, sort_keys=True) == before


def test_small_blunt_sweep_paths_consume_exact_canonical_references():
    # 256 controls, with a literal concatenation and coordinate oracle.
    rc = lambda s: s.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    for v, d, j in product(("AA", "CG", "GC", "TT"), repeat=3):
        for inverted, order in product((False, True), (("vd", "dj"), ("dj", "vd"))):
            p = project(authored(v, d, j, inverted=inverted, order=order))
            assert p["record"]["sequence"] == v + (rc(d) if inverted else d) + j
            for i, name in enumerate(("v", "d", "j")):
                assert_path(p, name, (2*i, 2*i+2), (0, 2), "2M")
