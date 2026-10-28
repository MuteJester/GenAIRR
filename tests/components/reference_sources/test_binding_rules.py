"""What a locked source binds: ids, mappings, limits and labels.

These are the rules that make a reference reproducible -- a content-derived id
belongs to one lock, a mapping must match the values actually observed, and
identical bytes must seal identically.
"""
from __future__ import annotations
import copy
import json
import pytest
from GenAIRR._artifact_io import canonical_json, strict_loads, typed_digest
from GenAIRR._reference_sources import model as m
from tests.support.fixtures.reference_source_models import ADAPTER, CLAIM_ID, FINDING_ID, FIXED, GENERATOR, IGH, INV_RECORDS, Inv, L, ROLE, ROOT, RUNTIME, SELECTED, SOURCES, _RECORDS, _blocks, _capabilities, _claim, _fasta, _find, _finding, _gapped_index, _load_json, _nested, _pointer, _raises, _recipe, _records, _reseal, _sorted, _source_set, decision, inv, inventory, lock_dict, manifest_dict, recipe_dict, report, report_d, report_decisions, report_dict, selection_entry

def test_content_ids_change_with_the_lock_even_for_identical_content(inventory):
    other_lock = m.SourceLock.from_dict(lock_dict(origin="local_file"))
    other = Inv(other_lock).build()
    assert [r.locator for r in other.records] == [r.locator for r in inventory.records]
    for name in ("claims", "mappings", "findings"):
        mine, theirs = getattr(inventory, name), getattr(other, name)
        assert len(mine) == len(theirs)
        assert not {x.id for x in mine} & {x.id for x in theirs}, name
        assert {x.source_lock_hash for x in theirs} == {other_lock.commitment}
    body = {k: v for k, v in inventory.claims[0].to_dict().items() if k != "id"}
    assert body["source_lock_hash"] == inventory.source_lock_hash
    assert inventory.claims[0].id == "claim:" + typed_digest(body, "genairr_reference_source_claim_id_v1")[7:]


def test_items_from_another_lock_cannot_join_an_inventory(inv, inventory):
    other = Inv(m.SourceLock.from_dict(lock_dict(origin="local_file")))
    for name in ("claims", "findings"):
        with _raises("digest_mismatch") as caught:
            inv.build(**{name: list(getattr(inv, name))[1:] + [getattr(other, name)[0]]})
        assert caught.error.path.startswith(f"/{name}/") and caught.error.path.endswith("/source_lock_hash")


def test_decisions_go_stale_when_the_recipe_is_pointed_at_another_lock(inventory, inv):
    decisions = report_decisions(inventory, inv)
    lock = m.SourceLock.from_dict(lock_dict())
    m.check_build_inputs(lock, inventory, _recipe(inventory, decisions=decisions))
    other_lock = m.SourceLock.from_dict(lock_dict(origin="local_file"))
    other_inventory = Inv(other_lock).build()
    repointed = _recipe(inventory, decisions=decisions,
                        edit=lambda d: d.update(source_lock_hash=other_lock.commitment))
    with _raises("stale_decision", "/decisions/0/findings/0"):
        m.check_build_inputs(other_lock, other_inventory, repointed)


def _replace_mapping(inv, index, **changes):
    """Rebuild the inventory with mapping ``index`` changed (claims follow its new id)."""
    import dataclasses
    old = inv.mappings[index]
    new = dataclasses.replace(old, **changes)
    claims = [dataclasses.replace(c, mappings=tuple(sorted(new.id if x == old.id else x for x in c.mappings)))
              for c in inv.claims]
    mappings = list(inv.mappings)
    mappings[index] = new
    return inv.build(mappings=mappings, claims=claims)


def test_alignment_blocks_are_the_exact_runs_of_ungapped_positions():
    gapped = _gapped_index("IGHV1-58*03", "gapped")[1]
    assert m.alignment_blocks(gapped) == _blocks(gapped)
    assert m.alignment_blocks("AC..G.T") == (m.Block((0, 2), (0, 2)), m.Block((2, 3), (4, 5)), m.Block((3, 4), (6, 7)))
    assert m.IMGT_GAP == "."


def test_fasta_record_observations_have_a_frozen_value_shape():
    loc = L("gapped", "/3")
    m.Observation(locator=loc, record=None, presence="present", value={"header": "IGHV1-2*02", "sequence": "CAG..GTG"},
                  interpretation="interpreted")
    for value in ("CAG..GTG", {"header": "IGHV1-2*02"}, {"header": "IGHV1-2*02", "sequence": 7},
                  {"header": None, "sequence": "CAG"}, {"header": "x", "sequence": "CAG", "description": ""}):
        with _raises("malformed", "/value"):
            m.Observation(locator=loc, record=None, presence="present", value=value, interpretation="interpreted")
    with _raises("malformed", "/presence"):
        m.Observation(locator=loc, record=None, presence="absent", value=None, interpretation="interpreted")


def test_gapped_mapping_output_digest_must_match_the_observed_gapped_record(inv):
    with _raises("digest_mismatch") as caught:
        _replace_mapping(inv, 0, output_sha256="sha256:" + "a" * 64)
    assert caught.error.path.startswith("/mappings/") and caught.error.path.endswith("/output_sha256")


def test_exact_gapped_blocks_must_reproduce_the_candidate(inv):
    n = inv.mappings[0].blocks[-1].input[1]
    with _raises("malformed") as caught:
        _replace_mapping(inv, 0, blocks=(m.Block(input=(0, n), output=(1, n + 1)),))
    assert caught.error.path.endswith("/blocks")
    with _raises("malformed") as caught:
        _replace_mapping(inv, 0, blocks=(), validation="mismatch")
    assert caught.error.path.endswith("/validation")


def test_coding_interval_mapping_is_verified_against_the_record_field(inv):
    index, raw = _find("IGHJ4*02")
    rec, cand = inv.rec["airr", "IGHJ4*02"], inv.cand["airr", "IGHJ4*02"]
    start, end = raw["gene_start"] - 1, raw["gene_end"]
    seq_loc = L("airr", rec.pointer + "/sequence")
    observation = m.Observation(locator=seq_loc, record=rec, presence="present", value=raw["sequence"],
                                interpretation="interpreted")

    def mapping(**changes):
        fields = dict(source_lock_hash=inv.lock.commitment, kind="coding_interval", input=seq_loc, output=cand,
                      input_sha256=m.text_digest(raw["sequence"]), output_sha256=m.sequence_digest(raw["coding_sequence"]),
                      input_space="record_sequence", output_space="coding_sequence",
                      blocks=(m.Block((start, end), (0, end - start)),), transform="extract_interval",
                      validation="exact")
        fields.update(changes)
        return inv.build(observations=inv.observations + [observation], mappings=inv.mappings + [m.Mapping(**fields)])

    mapping()
    with _raises("digest_mismatch") as caught:
        mapping(input_sha256="sha256:" + "b" * 64)
    assert caught.error.path.endswith("/input_sha256")
    with _raises("malformed") as caught:
        mapping(blocks=(m.Block((start + 1, end + 1), (0, end - start)),))
    assert caught.error.path.endswith("/blocks")


def test_record_fields_are_observed_with_their_record(inv):
    rec = inv.rec["airr", "IGHV1-69*01"]
    orphan = m.Observation(locator=L("airr", rec.pointer + "/curation"), record=None, presence="present", value="x",
                           interpretation="preserved")
    with _raises("malformed") as caught:
        inv.build(observations=inv.observations + [orphan])
    assert caught.error.path.startswith("/observations/") and caught.error.path.endswith("/record")


def test_mention_locators_are_observed_fields(inv):
    ex69 = inv.rec["airr_ex", "IGHV1-69*01"]
    unobserved = m.MentionedIdentity(asset="airr_ex", segment="V", label="IGHV1-69D*01", mentions=(
        m.Mention(record=ex69, locator=L("airr_ex", ex69.pointer + "/paralogs/1"), relationship="paralog"),))
    with _raises("malformed", "/mentioned_identities/0/mentions/0/locator"):
        inv.build(mentioned_identities=[unobserved])


def test_alignment_derived_claims_rest_on_exact_mappings(inv):
    cand, rec = inv.cand["airr", "IGHV1-69*01"], inv.rec["airr", "IGHV1-69*01"]
    gindex, gapped = _gapped_index("IGHV1-69*01", "gapped")
    last = max(i for i, char in enumerate(gapped) if char != ".")
    gapped = gapped[:last] + ("A" if gapped[last] != "A" else "C") + gapped[last + 1:]   # one base differs
    gloc = L("gapped", f"/{gindex}")
    observation = m.Observation(locator=gloc, record=None, presence="present",
                                value={"header": "IGHV1-69*01", "sequence": gapped}, interpretation="interpreted")
    mismatch = m.Mapping(source_lock_hash=inv.lock.commitment, kind="gapped_alignment", input=cand, output=gloc,
                         input_sha256=m.sequence_digest(_find("IGHV1-69*01")[1]["coding_sequence"]),
                         output_sha256=m.text_digest(gapped), input_space="coding_sequence",
                         output_space="imgt_gapped", blocks=(), transform="remove_imgt_gaps", validation="mismatch")
    fields = dict(observations=inv.observations + [observation], mappings=inv.mappings + [mismatch])
    inv.build(**fields)
    cdr3 = inv.obs_at(rec, "v_gene_delineations/0/cdr3_start")
    claim = inv.claim(cand, "anchor", 285, [cdr3, gloc], system="coding_sequence", derivation="imgt_alignment",
                      mappings=[mismatch.id])
    with _raises("malformed") as caught:
        inv.build(claims=inv.claims + [claim], **fields)
    assert caught.error.path.endswith("/mappings/0") and "mismatch" in caught.error.message


def _shared_bytes_lock(bytes_delta=0):
    d = lock_dict()
    (gapped,) = [a for a in d["assets"] if a["format"] == "gapped"]
    for asset in d["assets"]:
        if asset["format"] == "gapped_ex":
            asset["sha256"], asset["bytes"] = gapped["sha256"], gapped["bytes"] + bytes_delta
    return _reseal(d)


def test_one_digest_has_one_size_in_a_lock():
    m.SourceLock.from_dict(_shared_bytes_lock())
    with _raises("malformed", "/assets") as caught:
        m.SourceLock.from_dict(_shared_bytes_lock(bytes_delta=1))
    assert "gapped" in caught.error.message


def test_manifest_holds_one_blob_for_assets_with_identical_bytes(report):
    lock = m.SourceLock.from_dict(_shared_bytes_lock())
    d = manifest_dict(report)
    d["hashes"]["source_lock"] = lock.commitment
    names = {"cache/blobs/sha256/" + a.blob_name: a for a in lock.assets}
    d["members"] = [x for x in d["members"] if not x["name"].startswith("cache/")]
    d["members"] += [{"name": n, "sha256": a.sha256, "bytes": a.bytes} for n, a in names.items()]
    d["members"].sort(key=lambda x: x["name"])
    assert len(names) == 3
    m.check_manifest(m.BuildManifest.from_dict(d), lock)
    shared = "cache/blobs/sha256/" + lock.asset("gapped").blob_name
    d["members"] = [x for x in d["members"] if x["name"] != shared]
    with _raises("missing_asset", "/members") as caught:
        m.check_manifest(m.BuildManifest.from_dict(d), lock)
    assert "gapped" in caught.error.message


def test_withholding_a_claim_the_runtime_cannot_hold_is_a_projection_error(inventory, inv):
    ex = inv.cand["airr_ex", "IGHV1-58*03"]
    extension = _claim(inventory, "extension", ex)
    selected = tuple(label for label in SELECTED if label != "IGHV1-58*03")

    def edit(d):
        d["selection"].insert(0, selection_entry("IGHV1-58*03", fmt="airr_ex"))
    entry = decision("withhold_claim", claims=[extension.id], sequence_sha256=inventory.candidate(ex).sha256,
                     reason="do not trust the extension")
    recipe = _recipe(inventory, selected=selected, decisions=[entry], edit=edit)
    with _raises("projection", "/decisions/0/claims/0"):
        m.check_decisions(recipe, inventory)
    with _raises("projection", "/decisions/0/claims/0"):
        m.build_input_hash(m.SourceLock.from_dict(lock_dict()), inventory, recipe)


def test_record_claim_decisions_carry_the_selected_sequence_digest(inv):
    rec = inv.rec["airr", "IGHV1-69*01"]
    other_seq = "ACGT" * 30
    extra = m.SequenceCandidate(locator=L("airr", rec.pointer + "/sequence"), record=rec, sequence=other_seq,
                                sha256=m.sequence_digest(other_seq))
    inventory = inv.build(candidates=list(inv.candidates) + [extra])
    functional = _claim(inventory, "functionality", rec)
    with _raises("stale_decision", "/decisions/0/sequence_sha256"):
        m.check_decisions(_recipe(inventory, decisions=[
            decision("withhold_claim", claims=[functional.id], sequence_sha256=extra.sha256)]), inventory)
    selected_sha = selection_entry("IGHV1-69*01")["sequence_sha256"]
    m.check_decisions(_recipe(inventory, decisions=[
        decision("withhold_claim", claims=[functional.id], sequence_sha256=selected_sha)]), inventory)


def test_independent_phase_is_never_anchor_derived(inv):
    jcand = inv.cand["airr", "IGHJ4*02"]
    frame = inv.obs_at(inv.rec["airr", "IGHJ4*02"], "j_codon_frame")
    with _raises("malformed", "/derivation"):
        inv.claim(jcand, "phase_independence", "independent", [frame], derivation="anchor_derived")
    inv.claim(jcand, "phase_independence", "independent", [frame])


def test_a_record_of_a_selected_identity_is_an_alternative_or_duplicate(report_d):
    d = copy.deepcopy(report_d)
    (index,) = [i for i, r in enumerate(d["records"]) if r["exclusion"] == "duplicate_occurrence"][:1]
    d["records"][index].update(exclusion="outside_explicit_selection", related=None)
    with _raises("malformed", f"/records/{index}/exclusion"):
        m.BuildReport.from_dict(d)


def test_one_finding_is_resolved_by_one_selection_decision():
    d = recipe_dict()
    d["decisions"] = [decision(reason="first"), decision(reason="second")]
    with _raises("malformed", "/decisions"):
        m.Recipe.from_dict(d)


def test_over_long_names_are_resource_limits(report):
    d = lock_dict()
    d["assets"][0]["asset_id"] = "a" * 65
    with _raises("resource_limit", "/assets/0/asset_id"):
        m.SourceLock.from_dict(_reseal(d))
    d = manifest_dict(report)
    d["members"][0]["name"] = "x" * 256
    with _raises("resource_limit", "/members/0/name") as caught:
        m.BuildManifest.from_dict(d)
    assert "256" in caught.error.message and "255" in caught.error.message


def test_seal_locates_a_field_that_is_not_json():
    lock = m.SourceLock.from_dict(lock_dict())
    with _raises("malformed", "/release"):
        m.SourceLock.seal(provider="ogrdb", set_id=lock.set_id, release=object(), identity=lock.identity,
                          adapter=lock.adapter, assets=lock.assets, missing_assets=lock.missing_assets, terms=lock.terms)


def test_selection_checks_scale_linearly():
    import random
    import time
    lock = m.SourceLock.from_dict(lock_dict())
    rng = random.Random(7)
    records, candidates, selection = [], [], []
    for i in range(20_000):
        loc = L("airr", _pointer(i))
        label = f"IGHV9-{i}*01"
        records.append(m.SourceRecord(locator=loc, status="parsed", problem=None, label=label, provider_id=None,
                                      provider_ref=None, segment="V", taxonomy_id=9606, locus="IGH"))
        sequence = "".join(rng.choice("ACGT") for _ in range(24))
        candidates.append(m.SequenceCandidate(locator=L("airr", loc.pointer + "/coding_sequence"), record=loc,
                                              sequence=sequence, sha256=m.sequence_digest(sequence)))
        selection.append({"segment": "V", "label": label, "record": loc.to_dict(),
                          "sequence_sha256": m.sequence_digest(sequence), "reason": "x"})
    inventory = m.SourceInventory.assemble(source_lock_hash=lock.commitment, adapter=lock.adapter, records=records,
                                           candidates=candidates, mentioned_identities=[], observations=[],
                                           support=[], claims=[], mappings=[], findings=[])
    d = recipe_dict()
    d["selection"] = selection
    recipe = m.Recipe.from_dict(d)
    started = time.perf_counter()
    m.check_selection(recipe, inventory)
    m.check_decisions(recipe, inventory)
    assert len(m.selected_subjects(recipe, inventory)) == 40_000
    assert time.perf_counter() - started < 2.0


@pytest.mark.parametrize("set_id", [IGH, "9606.IGKappa_VJ", "9606.IGLambda_VJ"])
def test_pinned_gapped_headers_are_exactly_record_labels(set_id):
    for fmt, gapped_fmt in (("airr", "gapped"), ("airr_ex", "gapped_ex")):
        labels = {r["label"] for r in _records(fmt, set_id)}
        headers = [name for name, _ in _fasta(gapped_fmt, set_id)]
        assert headers and all(h in labels for h in headers)


def test_gapped_mapping_header_must_equal_the_subject_label(inv):
    cand = inv.cand["airr", "IGHV1-58*03"]
    # The real gapped IGHV1-58*01 record: its header differs only by the allele suffix.
    gindex, gapped = _gapped_index("IGHV1-58*01", "gapped")
    gloc = L("gapped", f"/{gindex}")
    sibling = m.Observation(locator=gloc, record=None, presence="present",
                            value={"header": "IGHV1-58*01", "sequence": gapped}, interpretation="interpreted")
    mapping = m.Mapping(source_lock_hash=inv.lock.commitment, kind="gapped_alignment", input=cand, output=gloc,
                        input_sha256=inv.build().candidate(cand).sha256, output_sha256=m.text_digest(gapped),
                        input_space="coding_sequence", output_space="imgt_gapped", blocks=(),
                        transform="remove_imgt_gaps", validation="mismatch")
    with _raises("malformed") as caught:
        inv.build(observations=inv.observations + [sibling], mappings=[mapping, inv.mappings[1]],
                  claims=[c for c in inv.claims if not c.mappings])
    assert caught.error.path.startswith("/mappings/") and caught.error.path.endswith("/output")
    assert "IGHV1-58*01" in caught.error.message and "IGHV1-58*03" in caught.error.message
    # The exact own-label alignment is accepted; relabelling its header is not.
    own = inv.mappings[0].output
    relabelled = []
    for observation in inv.observations:
        if observation.locator == own:
            value = observation.value.get()
            observation = m.Observation(locator=own, record=None, presence="present",
                                        value={"header": "IGHV1-58*01", "sequence": value["sequence"]},
                                        interpretation="interpreted")
        relabelled.append(observation)
    with _raises("malformed") as caught:
        inv.build(observations=relabelled)
    assert caught.error.path.endswith("/output")


def _value_roots(value, depth=0, path=()):
    """Depths of every preserved-value root (observation value, claim value, finding observations)."""
    if isinstance(value, dict):
        if "interpretation" in value and "presence" in value:
            yield depth + 1, path + ("value",)
        if "predicate" in value and "derivation" in value:
            yield depth + 1, path + ("value",)
        if "code" in value and "severity" in value:
            yield depth + 1, path + ("observations",)
        for key, item in value.items():
            yield from _value_roots(item, depth + 1, path + (key,))
    elif isinstance(value, list):
        for i, item in enumerate(value):
            yield from _value_roots(item, depth + 1, path + (i,))


def test_embedded_value_depth_leaves_room_for_the_deepest_embedding(report):
    deepest = max(depth for depth, _ in _value_roots(report.to_dict()))
    assert deepest == m.VALUE_EMBEDDING_DEPTH == 4
    assert m.MAX_EMBEDDED_VALUE_DEPTH == m.MAX_JSON_DEPTH - m.VALUE_EMBEDDING_DEPTH == 60
    assert m.MAX_JSON_DEPTH == 64


def _deep_report(inv, depth):
    deep = m.Observation(locator=L("airr", "/Info/deep"), record=None, presence="present", value=_nested(depth),
                         interpretation="preserved")
    inventory = inv.build(observations=inv.observations + [deep])
    return report_dict(inventory, inv)


def test_a_value_at_the_embedded_depth_limit_survives_a_report(inv):
    d = _deep_report(inv, m.MAX_EMBEDDED_VALUE_DEPTH)
    report = m.BuildReport.from_dict(d)
    assert m.report_hash(report) == typed_digest(d, "genairr_reference_source_build_report_v1")
    text = canonical_json(report.to_dict())
    again = m.BuildReport.from_json(text)
    assert again == report and canonical_json(again.to_dict()) == text
    assert m.BuildReport.from_dict(report.to_dict()) == report


def test_a_value_one_level_deeper_fails_at_construction(inv):
    with _raises("resource_limit") as caught:
        m.Observation(locator=L("airr", "/Info/deep"), record=None, presence="present",
                      value=_nested(m.MAX_EMBEDDED_VALUE_DEPTH + 1), interpretation="preserved")
    assert caught.error.path.startswith("/value/0") and str(m.MAX_EMBEDDED_VALUE_DEPTH) in caught.error.message
    for build in (lambda v: m.Claim(source_lock_hash=inv.lock.commitment, subject=inv.cand["airr", "IGHJ4*02"],
                                    predicate="extension", value=v, coordinate_system=None, support="none",
                                    derivation="variant_comparison", status="candidate",
                                    evidence=(L("airr", inv.cand["airr", "IGHJ4*02"].pointer),), mappings=()),
                  lambda v: m.Finding(source_lock_hash=inv.lock.commitment, code="deep", severity="info",
                                      category="informational", subjects=(), locators=(L("airr", "/Info/version"),),
                                      observations={"value": v}, capabilities=())):
        with _raises("resource_limit"):
            build(_nested(m.MAX_EMBEDDED_VALUE_DEPTH + 1))


def test_raw_provider_payloads_still_parse_to_depth_64():
    from GenAIRR._reference_sources._fields import parse_json
    assert parse_json(json.dumps(_nested(m.MAX_JSON_DEPTH)), max_bytes=1 << 20) == _nested(m.MAX_JSON_DEPTH)
    with _raises("resource_limit"):
        parse_json(json.dumps(_nested(m.MAX_JSON_DEPTH + 1)), max_bytes=1 << 20)
