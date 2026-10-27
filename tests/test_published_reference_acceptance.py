"""End-to-end acceptance of published germline sources and the portable references built from them.

``ACCEPTANCE`` maps every acceptance criterion, A01 to A30, to the tests that exercise it, here
and in the source, projection, dossier, packaging and runtime test files.  A test checks that
the table is closed and that every test it names exists.  The tests defined here cover what the
others do not:

- an IGK VJ reference built through the public API from the pinned IGK release 5 exports, with a
  recipe made by a stated mechanical rule: record accounting, usage weights, capabilities, a
  capability failure, a moved snapshot and record-validated simulation;
- ordinary generation, somatic hypermutation and polymerase indels on the packaged human IGH
  starter, checked record by record against the exact selected source alleles, which are read from
  the raw OGRDB ``airr`` export with the standard library;
- runtime unknowns, and the workflows a portable reference still rejects, on the packaged starter.

Expected values come from the independently authored fixtures under
``tests/data/reference_sources/`` (``expectations.json``, ``interpretation_evidence.json``, the
selection fixture and the reviewed starter recipe), from raw export bytes read with the standard
library, or from literals; never from the code under test.  Nothing here touches the network:
sockets, name lookups and the fetch transport are refused for every test and fixture.
"""
from __future__ import annotations

import collections
import contextlib
import dataclasses
import hashlib
import http.client
import importlib
import json
import pathlib
import re
import shutil
import socket
from typing import NamedTuple, Optional

import pytest

import GenAIRR
from GenAIRR import Experiment, ReferenceCapabilityError
from GenAIRR._reference_sources import api
from GenAIRR._reference_sources.acquisition import HttpsTransport
from GenAIRR.experimental import FamilyReplayError
from GenAIRR.reference_sources import OgrdbSource, ReferenceBuild, ReferenceSourceError
from tests.test_ogrdb_source_inventory import EXPECT, IGH, IGK, RECORDS, SELECTION, raw, raw_records

TESTS = pathlib.Path(__file__).resolve().parent
REPO = TESTS.parent
DATA = TESTS / "data" / "reference_sources"
EVIDENCE = json.loads((DATA / "interpretation_evidence.json").read_bytes().decode("utf-8"))
STARTER_ID = "human_igh_ogrdb_10_observed_v1"
STARTER_RECIPE = json.loads((DATA / f"{STARTER_ID}.recipe.json").read_bytes().decode("utf-8"))
PROFILE = "ogrdb_api_v2_human_ig_observed_v1"
ALL = ("airr", "airr_ex", "gapped", "gapped_ex")
SUFFIX = {"airr": "json", "airr_ex": "json", "gapped": "fasta", "gapped_ex": "fasta"}
IGK_KEY = "9606.IGKappa_VJ/5"
DEFERRED = "not yet supported for portable references"
LINEAGE = dict(n_clones=3, max_generations=6, n_max=200, n_sample=20, rate=0.05, lambda_base=1.6)
CODONS = {a + b + c: aa for (a, b, c), aa in zip(
    ((a, b, c) for a in "TCAG" for b in "TCAG" for c in "TCAG"),
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")}
V_REGIONS = ("n_fwr1_mutations", "n_cdr1_mutations", "n_fwr2_mutations", "n_cdr2_mutations", "n_fwr3_mutations")
P_LENGTHS = ("p_v_3_length", "p_d_5_length", "p_d_3_length", "p_j_5_length")

# ── acceptance criteria and the tests that exercise them ─────────────

_HERE = "tests/test_published_reference_acceptance.py::"
_ACQ = "tests/test_reference_source_acquisition.py::"
_API = "tests/test_reference_source_api.py::"
_CACHE = "tests/test_reference_source_cache.py::"
_DOSSIER = "tests/test_reference_source_dossier.py::"
_FIXTURES = "tests/test_reference_source_fixtures.py::"
_MODELS = "tests/test_reference_source_models.py::"
_PROJECTION = "tests/test_reference_source_projection.py::"
_COORDINATES = "tests/test_ogrdb_source_coordinates.py::"
_EVIDENCE = "tests/test_ogrdb_source_evidence.py::"
_INVENTORY = "tests/test_ogrdb_source_inventory.py::"
_RESOURCES = "tests/test_published_reference_resources.py::"
_EXPERIMENT = "tests/test_reference_experiment.py::"

ACCEPTANCE = {
    "A01": ("Two clean acquisitions of identical bytes have equal source locks; receipt times never change a "
            "build identity", (
                _ACQ + "test_independent_clean_acquisitions_of_the_same_bytes_have_equal_lock_identity",
                _MODELS + "test_receipts_with_different_times_describe_the_same_lock",
                _MODELS + "test_receipts_enter_no_digest",
                _DOSSIER + "test_builds_with_different_receipts_have_identical_identities")),
    "A02": ("A moved or renamed local snapshot gives the same lock and build, and no canonical output names a "
            "location", (
                _API + "test_a_moved_and_renamed_snapshot_has_the_same_lock",
                _ACQ + "test_identical_local_bytes_moved_between_directories_yield_an_identical_lock",
                _MODELS + "test_lock_assets_cannot_carry_acquisition_or_filesystem_details",
                _PROJECTION + "test_two_temporary_roots_give_identical_files_runtime_hashes_and_report",
                _DOSSIER + "test_a_clean_cache_build_and_a_moved_offline_dossier_agree",
                _HERE + "test_a_moved_and_renamed_igk_snapshot_builds_the_same_reference_and_names_no_location",
                _HERE + "test_the_location_check_skips_a_home_directory_at_the_filesystem_root")),
    "A03": ("An offline cache miss or corrupt blob is a typed missing/digest error with no transport call", (
        _API + "test_an_offline_cache_miss_or_corrupt_blob_is_a_typed_error",
        _ACQ + "test_a_missing_blob_names_its_asset_and_digest_and_never_downloads",
        _ACQ + "test_a_corrupt_cache_blob_fails_an_offline_read",
        _ACQ + "test_offline_reads_never_touch_the_network",
        _ACQ + "test_offline_operations_take_no_transport",
        _CACHE + "test_a_missing_blob_names_the_asset_and_expected_digest",
        _CACHE + "test_a_corrupt_cache_blob_is_rejected_on_every_read")),
    "A04": ("New bytes served under a pinned release are drift: the pin is rejected, nothing is refreshed or "
            "overwritten", (
                _API + "test_changed_bytes_under_the_same_release_are_drift_and_change_nothing",
                _ACQ + "test_changing_upstream_bytes_are_reported_as_drift_and_change_nothing",
                _ACQ + "test_failures_never_replace_a_valid_cache_blob_or_lock",
                _CACHE + "test_a_conflicting_existing_file_is_an_error_and_is_never_overwritten")),
    "A05": ("An unsupported profile, release or schema is an explicit unsupported result with no motif-based "
            "fallback", (
                _API + "test_fetch_rejects_implicit_or_unsupported_requests_before_any_transport_call",
                _API + "test_a_lock_of_another_schema_version_is_unsupported",
                _API + "test_build_errors_keep_recipe_paths_and_name_the_recipe_file",
                _API + "test_an_out_of_scope_release_is_inspected_and_reported_as_such",
                _MODELS + "test_lock_rejects_unknown_semantic_profiles",
                _MODELS + "test_recipe_rejects_unknown_schema_and_profiles",
                _MODELS + "test_locks_outside_the_profile_scope_are_valid_but_cannot_build",
                _MODELS + "test_drifted_bytes_for_a_validated_release_are_out_of_scope",
                _COORDINATES + "test_profile_lookup_is_explicit_and_never_guessed",
                _COORDINATES + "test_the_airr_prose_reading_is_never_a_fallback",
                _EVIDENCE + "test_outside_the_validated_scope_no_claim_is_accepted",
                _DOSSIER + "test_a_dossier_outside_this_versions_validated_scope_is_unsupported")),
    "A06": ("The IGL endpoint pair keeps 97 versus 98 explicit records and the exact suffix and variant "
            "observations", (
                _INVENTORY + "test_igl_release4_explicit_records_and_mentioned_identities",
                _INVENTORY + "test_igl_release4_ext_suffix_on_every_extended_ref_and_no_plain_ref",
                _INVENTORY + "test_igl_release4_six_coding_differences_match_the_coding_end_table",
                _INVENTORY + "test_pinned_pairs_differ_only_in_provider_ref_and_coding_sequence",
                _FIXTURES + "test_set_level_counts_and_variant_differences")),
    "A07": ("The five historical two-base IGL cases have their own code, bases and context, and keep every note "
            "and genomic-support entry", (
                _EVIDENCE + "test_the_five_igl_coding_ends_before_a_declared_v_rss",
                _EVIDENCE + "test_igl_extensions_two_base_class_and_the_separate_disagreeing_extension",
                _EVIDENCE + "test_only_genomic_support_entries_are_listed_as_genomic_support",
                _COORDINATES + "test_note_quotes_keep_the_providers_own_confirmations",
                _FIXTURES + "test_igl_coding_end_observations")),
    "A08": ("IGLV3-9*01 has a separate context-disagreement code, and CACCC is never relocated", (
        _EVIDENCE + "test_igl_extensions_two_base_class_and_the_separate_disagreeing_extension",
        _EVIDENCE + "test_the_five_igl_coding_ends_before_a_declared_v_rss",
        _INVENTORY + "test_coding_sequence_differences_are_the_expected_extensions")),
    "A09": ("V/RSS interval accounting: 82 adjacent, five two-base, one missing RSS", (
        _EVIDENCE + "test_igl_v_rss_gaps_name_the_two_base_intervals_and_the_missing_rss",
        _EVIDENCE + "test_the_gap_between_coding_end_and_v_rss_start_is_counted_per_export")),
    "A10": ("IGLJ2/IGLJ3: explicit and mentioned identities are kept apart, with no cloned support or extra alias "
            "sampling weight", (
                _INVENTORY + "test_j2_and_j3_keep_their_own_records_ids_and_support",
                _INVENTORY + "test_igl_release4_j3_is_absent_from_the_extended_export_and_named_there_by_j2",
                _PROJECTION + "test_a_paralog_name_without_a_selected_record_stays_a_named_entity_without_weight")),
    "A11": ("J coordinate and phase conventions match independent hand-checked values for both conventions, all "
            "phase classes and nonzero source offsets", (
                _COORDINATES + "test_hand_checked_j_literals",
                _COORDINATES + "test_j_examples_cover_all_three_phase_classes_and_nonzero_context_offsets",
                _COORDINATES + "test_the_airr_prose_reading_is_a_different_convention",
                _COORDINATES + "test_j_landmark_follows_the_producer_convention",
                _COORDINATES + "test_every_j_occurrence_gets_a_source_stated_anchor_and_phase_value",
                _FIXTURES + "test_j_examples")),
    "A12": ("V delineations follow the proven convention, unsupported conventions stay unresolved, and "
            "noncanonical residues are never repaired by motif", (
                _COORDINATES + "test_hand_checked_v_literals",
                _COORDINATES + "test_v_landmark_is_bound_to_its_own_gapped_record",
                _COORDINATES + "test_every_v_occurrence_is_bound_to_its_own_gapped_record_and_accepted",
                _COORDINATES + "test_a_noncanonical_v_landmark_is_never_moved_to_a_nearby_cysteine",
                _COORDINATES + "test_noncanonical_landmark_codons_are_findings_and_never_relocated",
                _COORDINATES + "test_two_imgt_delineations_are_not_chosen_between",
                _COORDINATES + "test_profile_keeps_the_unresolved_conventions_unresolved",
                _FIXTURES + "test_v_examples")),
    "A13": ("A matching, mismatching or missing gapped asset gives an exact map or an explicitly unresolved "
            "affirmed span; a gapped column is never used as an ungapped index", (
                _COORDINATES + "test_a_reviewed_note_maps_through_an_exact_gapped_alignment",
                _COORDINATES + "test_a_reviewed_note_stays_unresolved_without_a_matching_alignment_or_sequence",
                _COORDINATES + "test_a_gapped_record_that_does_not_ungap_to_the_coding_sequence_is_not_bound",
                _COORDINATES + "test_without_gapped_exports_no_v_anchor_is_accepted",
                _EVIDENCE + "test_a_reviewed_note_would_give_an_affirmed_span_through_the_exact_alignment",
                _EVIDENCE + "test_every_span_limiting_note_is_reported_unresolved")),
    "A14": ("Equal DNA under distinct genes stays distinct executable identities where selected", (
        _INVENTORY + "test_identical_dna_under_distinct_labels_stays_distinct_identities",
        _MODELS + "test_recipe_selects_distinct_identities_with_identical_dna",
        _EVIDENCE + "test_identical_dna_identities_cite_only_their_own_records",
        _PROJECTION + "test_identical_dna_under_distinct_labels_stays_two_identities_with_their_own_evidence",
        _RESOURCES + "test_identical_dna_reasons_name_the_other_identity_and_keep_each_ones_own_evidence",
        _HERE + "test_an_igk_vj_build_from_a_mechanical_recipe_accounts_for_every_source_record",
        _HERE + "test_ordinary_igh_generation_uses_the_exact_selected_source_alleles")),
    "A15": ("The same label with different DNA needs an exact record choice; the first file never wins", (
        _PROJECTION + "test_a_selection_binds_an_exact_record_and_digest_never_the_first_label_match",
        _PROJECTION + "test_one_executable_label_cannot_select_two_records",
        _MODELS + "test_recipe_rejects_two_alternatives_for_one_executable_label",
        _MODELS + "test_check_selection_binds_exact_record_and_digest",
        _RESOURCES + "test_the_recipe_selects_the_airr_record_of_each_record_choice_explicitly")),
    "A16": ("Missing phase or end evidence gives honest capability blockers and null outputs, never an "
            "inference from functional or full-length tags", (
                _PROJECTION + "test_the_starter_reports_honest_blockers_for_coding_ends_it_does_not_establish",
                _PROJECTION + "test_withheld_functionality_says_the_provider_value_cannot_be_classified",
                _PROJECTION + "test_one_capability_error_names_every_unmet_requirement",
                _EVIDENCE + "test_no_coding_boundary_is_ever_claimed_complete",
                _EVIDENCE + "test_tags_functionality_rss_adjacency_and_exact_slices_set_no_boundary",
                _EVIDENCE + "test_provider_functionality_stays_raw_and_unclassified",
                _RESOURCES + "test_the_bundled_starter_states_its_capabilities_truthfully",
                _HERE + "test_the_igk_build_states_its_capabilities_and_keeps_unknowns_null",
                _HERE + "test_an_igk_recipe_requiring_original_coding_ends_fails_with_a_capability_error",
                _HERE + "test_the_bundled_starter_keeps_unknowns_unknown_and_names_its_blockers")),
    "A17": ("An exact explicit selection accounts for every exclusion, in deterministic order, with the selected "
            "DNA unchanged", (
                _PROJECTION + "test_every_source_record_has_exactly_one_outcome_and_the_counts_reconcile",
                _PROJECTION + "test_starter_fasta_holds_exactly_the_selected_labels_and_supplied_dna",
                _PROJECTION + "test_starter_runtime_alleles_are_exactly_the_selected_records",
                _API + "test_the_report_text_accounts_for_every_source_record",
                _HERE + "test_an_igk_vj_build_from_a_mechanical_recipe_accounts_for_every_source_record")),
    "A18": ("A duplicate alias or record never increases an allele's usage weight", (
        _INVENTORY + "test_duplicate_occurrences_and_same_label_alternatives_are_counted_apart",
        _INVENTORY + "test_declared_aliases_are_names_of_their_own_record_not_identities",
        _PROJECTION + "test_starter_runtime_alleles_are_exactly_the_selected_records",
        _PROJECTION + "test_starter_relationships_are_the_selected_records_own_aliases_and_paralogs",
        _HERE + "test_an_igk_vj_build_from_a_mechanical_recipe_accounts_for_every_source_record")),
    "A19": ("A sequence-bound decision goes stale after a lock, digest or subject change", (
        _MODELS + "test_stale_decisions_are_rejected",
        _MODELS + "test_decisions_go_stale_when_the_recipe_is_pointed_at_another_lock",
        _MODELS + "test_decisions_on_unselected_sequences_are_stale",
        _PROJECTION + "test_stale_decisions_are_rejected")),
    "A20": ("The runtime projection round-trips through the unchanged reader with no provider-only fields; the "
            "loss report names the missing subregions; unsupported phase/landmark pairs are refused or withheld", (
                _PROJECTION + "test_the_runtime_reader_accepts_the_generated_files_unchanged",
                _PROJECTION + "test_the_sidecar_is_canonical_lf_json_in_the_v1_shape_without_provider_fields",
                _PROJECTION + "test_projection_losses_name_what_the_runtime_reference_does_not_hold",
                _PROJECTION + "test_the_subregion_loss_names_its_consequence_for_simulated_v_mutation_counts",
                _PROJECTION + "test_a_frameshifted_selected_allele_without_a_withholding_decision_is_refused",
                _PROJECTION + "test_a_frameshifted_allele_whose_anchor_and_phase_are_withheld_builds_with_null_"
                              "annotations",
                _API + "test_the_report_text_names_every_projection_loss",
                _HERE + "test_the_bundled_starter_keeps_unknowns_unknown_and_names_its_blockers",
                _HERE + "test_shm_is_recorded_against_the_generating_allele_and_never_against_a_recall")),
    "A21": ("A rehashed but semantically tampered dossier is rejected by regeneration, including altered "
            "provenance strength and hidden exclusions", (
                _DOSSIER + "test_an_altered_mapping_is_rejected",
                _DOSSIER + "test_an_altered_claim_support_is_rejected",
                _DOSSIER + "test_a_hidden_exclusion_is_rejected",
                _DOSSIER + "test_a_selected_label_changed_everywhere_but_the_source_is_rejected",
                _DOSSIER + "test_altered_native_output_is_rejected",
                _DOSSIER + "test_a_report_claiming_another_runtime_identity_is_rejected",
                _DOSSIER + "test_a_runtime_body_resealed_under_another_generator_is_rejected",
                _DOSSIER + "test_a_recipe_that_is_not_the_reports_recipe_is_rejected",
                _MODELS + "test_tampered_report_runtime_body_and_report_hash_are_both_caught")),
    "A22": ("A failure or race while writing leaves no valid-looking partial manifest, keeps complete prior files "
            "and cleans up nothing foreign", (
                _DOSSIER + "test_a_failed_write_leaves_an_incomplete_dossier_and_removes_nothing",
                _DOSSIER + "test_a_second_writer_fails_before_writing_anything",
                _DOSSIER + "test_a_manifest_published_by_another_writer_is_never_claimed",
                _DOSSIER + "test_a_failure_after_the_manifest_is_published_says_the_dossier_is_complete",
                _DOSSIER + "test_a_leftover_temporary_file_is_named_as_such",
                _CACHE + "test_a_concurrent_conflicting_writer_is_an_error_and_its_file_survives",
                _ACQ + "test_concurrent_fetches_into_one_cache_agree")),
    "A23": ("The 65,535/65,536-nt boundary and an invalid alphabet fail in Python before native entry", (
        _PROJECTION + "test_a_selected_coding_sequence_of_65535_bases_builds",
        _PROJECTION + "test_a_selected_coding_sequence_the_runtime_cannot_hold_is_refused_before_native_entry",
        _PROJECTION + "test_an_unselected_record_the_runtime_cannot_hold_is_a_structural_exclusion")),
    "A24": ("LF/CRLF/Unicode and Windows defaults: generated files are exact LF/UTF-8 bytes, and raw source "
            "digests truthfully distinguish source byte changes", (
                _PROJECTION + "test_crlf_source_bytes_change_lock_and_citations_but_not_dna_labels_or_annotation_"
                              "values",
                _PROJECTION + "test_unicode_labels_are_kept_literally_in_utf8_files_and_the_runtime",
                _API + "test_write_lock_writes_canonical_lf_text_and_never_replaces_silently",
                _FIXTURES + "test_generated_text_is_utf8_with_lf_endings",
                _INVENTORY + "test_fasta_line_terminators_are_removed_and_everything_else_is_kept",
                _MODELS + "test_lock_from_json_accepts_a_byte_order_mark_and_reads_files",
                _DOSSIER + "test_generated_members_must_be_in_canonical_form",
                _RESOURCES + "test_the_packaged_archive_and_its_sources_are_checked_out_byte_for_byte_on_every_"
                             "platform")),
    "A25": ("A moved dossier checks, rebuilds from its own sources and simulates in a new offline process with "
            "no network attempt", (
                _DOSSIER + "test_a_moved_dossier_checks_rebuilds_and_simulates_in_a_fresh_offline_process",
                _DOSSIER + "test_a_clean_cache_build_and_a_moved_offline_dossier_agree",
                _RESOURCES + "test_the_starter_loads_and_simulates_in_a_fresh_process_without_optional_packages_or_"
                             "network")),
    "A26": ("Human IGH and IGK VJ: actual calls, landmarks and alignments are checked against the selected source "
            "alleles, and record validation passes", (
                _RESOURCES + "test_the_bundled_starter_simulates_bound_to_its_execution_hash",
                _HERE + "test_an_igk_vj_build_from_a_mechanical_recipe_accounts_for_every_source_record",
                _HERE + "test_igk_simulation_is_validated_against_the_selected_source_alleles",
                _HERE + "test_ordinary_igh_generation_uses_the_exact_selected_source_alleles")),
    "A27": ("Mutation and indel truth: generating alleles stay authoritative, recalls are separate, and no "
            "alignment edit is lost", (
                _HERE + "test_shm_is_recorded_against_the_generating_allele_and_never_against_a_recall",
                _HERE + "test_polymerase_indels_keep_every_alignment_edit_against_the_generating_allele")),
    # The 106 bundled cartridges' native content hashes are not pinned in this suite: the comparison needs
    # the tree from before published sources were added, so it is made outside the suite.
    "A28": ("Compatibility: the existing golden reference artifact loads unchanged, the human_igh alias is "
            "unchanged, and unsupported integrations still reject", (
                "tests/test_reference_first_use_acceptance.py::test_golden_igh_artifact",
                _RESOURCES + "test_loading_the_starter_leaves_the_human_igh_alias_unchanged",
                _EXPERIMENT + "test_existing_inputs_unchanged_smoke",
                _EXPERIMENT + "test_genotype_and_cohort_rejected_with_capability_error",
                _EXPERIMENT + "test_curate_refdata_on_reference_rejected_at_the_call",
                _EXPERIMENT + "test_family_replay_capture_with_reference_unsupported",
                "tests/test_reference_capabilities.py::test_unsupported_capabilities_unavailable_with_note",
                _HERE + "test_workflows_portable_references_do_not_support_still_reject_the_bundled_starter")),
    # Only a fresh process with the optional packages blocked runs in this suite; installing the built
    # wheel into a fresh environment and running the starter there is a release check.
    "A29": ("The packaged starter works without a checkout, cache, network or optional dependencies", (
        _RESOURCES + "test_the_starter_loads_and_simulates_in_a_fresh_process_without_optional_packages_or_network",
        _RESOURCES + "test_the_starter_loads_from_a_read_only_package_location",
        _RESOURCES + "test_the_starter_needs_no_source_cache_or_network",
        _RESOURCES + "test_the_packaged_archive_is_included_in_the_wheel_and_the_sdist")),
    "A30": ("A dossier of an older generator checks under a supported contract; an offline rebuild records its "
            "own version with new artifact, report and manifest hashes and an equal execution hash; byte "
            "equality is claimed only for the same version", (
                _DOSSIER + "test_an_older_generators_dossier_checks_and_rebuilds_under_the_installed_one",
                _RESOURCES + "test_rebuilding_from_the_pinned_bytes_and_the_recipe_reproduces_the_packaged_starter",
                _RESOURCES + "test_the_packaged_dossier_was_generated_by_the_version_this_checkout_declares",
                _RESOURCES + "test_the_generator_writes_one_deterministic_stored_zip_of_the_dossier_members_"
                             "offline")),
}


def _test_function(node: str):
    """The test function a ``path::name`` node id names, or ``None``; found by importing its module."""
    path, separator, name = node.partition("::")
    if not separator or not path.endswith(".py") or "/" not in path or not (REPO / path).is_file():
        return None
    module = importlib.import_module(path[:-len(".py")].replace("/", "."))
    function = getattr(module, name.split("[", 1)[0], None)
    return function if name.startswith("test_") and callable(function) else None


def test_every_acceptance_criterion_is_mapped_to_tests_that_exist():
    assert list(ACCEPTANCE) == [f"A{n:02d}" for n in range(1, 31)]
    missing = []
    for row, (scenario, nodes) in ACCEPTANCE.items():
        assert scenario.strip() and nodes, row
        assert len(set(nodes)) == len(nodes), row
        missing.extend((row, node) for node in nodes if _test_function(node) is None)
    assert missing == []
    named_here = {node[len(_HERE):] for _, nodes in ACCEPTANCE.values() for node in nodes if node.startswith(_HERE)}
    defined_here = {name for name, value in globals().items() if name.startswith("test_") and callable(value)}
    assert defined_here - {"test_every_acceptance_criterion_is_mapped_to_tests_that_exist"} == named_here


# ── offline guard ────────────────────────────────────────────────────


def _refuse(*args, **kwargs):
    pytest.fail("an acceptance check tried to reach the network")


@contextlib.contextmanager
def no_network():
    """Every connection, name lookup and the fetch transport fail the test while this is active."""
    with pytest.MonkeyPatch.context() as patch:
        for owner, name in ((socket, "create_connection"), (socket, "getaddrinfo"), (socket.socket, "connect"),
                            (http.client.HTTPConnection, "connect"), (http.client.HTTPSConnection, "connect"),
                            (HttpsTransport, "get"), (api, "_transport")):
            patch.setattr(owner, name, _refuse)
        yield


@pytest.fixture(autouse=True)
def offline():
    with no_network():
        yield


# ── the selected source alleles, read from the raw exports ───────────


@dataclasses.dataclass(frozen=True)
class SourceAllele:
    """A selected allele as its raw 'airr' record states it."""
    segment: str
    label: str
    sequence: str
    anchor: Optional[int]          # zero-based landmark codon start in the coding sequence (V and J)


def sha(sequence: str) -> str:
    return "sha256:" + hashlib.sha256(sequence.encode("ascii")).hexdigest()


def gapped_landmarks(source) -> dict:
    """label -> number of residues before gapped column 310 (IMGT codon 104) in the raw 'gapped' export."""
    found = {}
    for chunk in raw(source, "gapped").decode("ascii").split(">")[1:]:
        header, _, body = chunk.partition("\n")
        found[header.strip()] = sum(char != "." for char in body.replace("\n", "")[:309])
    return found


def source_allele(source, index: int, landmarks: dict) -> SourceAllele:
    """The allele of 'airr' record ``index``; anchors follow the fixture's stated rules
    (V: cdr3_start - 4 of the IMGT delineation, bound to the label's own gapped record; J: j_cdr3_end - 1)."""
    record = raw_records(source, "airr")[index]
    rules = EXPECT["rules_under_test"]
    assert (rules["v_anchor_coding0"], rules["j_anchor_coding0"]) == ("cdr3_start - 4", "j_cdr3_end - 1")
    anchor = None
    if record["sequence_type"] == "V":
        (imgt,) = [d for d in record["v_gene_delineations"] if d["delineation_scheme"] == "IMGT"]
        anchor = imgt["cdr3_start"] - 4
        assert anchor == landmarks[record["label"]], record["label"]
    elif record["sequence_type"] == "J":
        anchor = record["j_cdr3_end"] - 1
    return SourceAllele(record["sequence_type"], record["label"], record["coding_sequence"], anchor)


def starter_alleles() -> dict:
    """label -> SourceAllele for every allele the reviewed starter recipe selects, from its own record."""
    landmarks = gapped_landmarks(IGH)
    alleles = {}
    for entry in STARTER_RECIPE["selection"]:
        pointer = entry["record"]["pointer"]
        assert entry["record"]["asset"] == "airr" and pointer.startswith(RECORDS + "/")
        allele = source_allele(IGH, int(pointer[len(RECORDS) + 1:]), landmarks)
        assert (allele.segment, allele.label, sha(allele.sequence)) == (
            entry["segment"], entry["label"], entry["sequence_sha256"])
        alleles[allele.label] = allele
    assert [(a.segment, a.label) for a in alleles.values()] == [
        (segment, label) for segment in ("V", "D", "J") for label in SELECTION["selected"][segment]]
    return alleles


def translate(sequence: str) -> str:
    return "".join(CODONS.get(sequence[i:i + 3], "X") for i in range(0, len(sequence) - 2, 3))


# ── record truth against the generating source alleles ───────────────


class Walk(NamedTuple):
    mapped: dict           # germline position -> sequence position, for aligned (M) bases
    differs: list          # sequence positions of aligned bases that differ from the source base
    edits: collections.Counter    # CIGAR operation -> count
    inserted: list         # sequence positions of inserted (I) bases
    deleted: list          # germline positions of deleted (D) bases


def walk(record, segment: str, allele: SourceAllele) -> Walk:
    """Walk the segment's CIGAR over the record's alignment strings against the generating allele's raw
    source sequence."""
    s = segment.lower()
    cigar = record[f"{s}_cigar"]
    assert re.fullmatch(r"(?:[0-9]+[MID])+", cigar), cigar
    sequence, observed, germline = record["sequence"], record["sequence_alignment"], record["germline_alignment"]
    start, end = record[f"{s}_alignment_start"], record[f"{s}_alignment_end"]
    query, position, column = record[f"{s}_sequence_start"], record[f"{s}_germline_start"], start
    assert position == record[f"{s}_trim_5"] == 0 and record[f"{s}_trim_3"] == 0     # trimming is off
    mapped, differs, edits, inserted, deleted = {}, [], collections.Counter(), [], []
    for count, op in re.findall(r"([0-9]+)([MID])", cigar):
        for _ in range(int(count)):
            assert germline[column] == (allele.sequence[position] if op in "MD" else "-"), (cigar, column)
            assert observed[column] == (sequence[query] if op in "MI" else "-"), (cigar, column)
            if op == "M":
                mapped[position] = query
                if sequence[query] != allele.sequence[position]:
                    differs.append(query)
            elif op == "I":
                inserted.append(query)
            else:
                deleted.append(position)
            edits[op] += 1
            query += op in "MI"
            position += op in "MD"
            column += 1
    assert (column, query) == (end, record[f"{s}_sequence_end"])
    assert position == record[f"{s}_germline_end"] == len(allele.sequence)
    assert germline[start:end].replace("-", "") == allele.sequence
    return Walk(mapped, differs, edits, inserted, deleted)


def check_calls(record, alleles: dict, segments: str, outcome, native) -> None:
    """Generating calls are single selected labels naming the alleles the engine sampled; recalls are
    separate fields naming only selected labels of the same segment."""
    simulation = outcome.final_simulation()
    for s in "vdj":
        if s not in segments:
            assert record[f"{s}_call"] == record[f"recalled_{s}_call"] == ""
            continue
        call = record[f"{s}_call"]
        assert call in alleles and alleles[call].segment == s.upper(), call
        assert getattr(native, f"{s}_allele")(getattr(simulation, f"{s}_allele_id")()).name == call
        recalled = record[f"recalled_{s}_call"].split(",")
        assert recalled and all(name in alleles and alleles[name].segment == s.upper() for name in recalled), recalled


def check_junction(record, mapped: dict, alleles: dict) -> str:
    """The junction runs from the observed first base of the V landmark codon to the last base of an intact
    J landmark codon.  With trimming off, only a deleted landmark base leaves it undefined."""
    v, j = alleles[record["v_call"]], alleles[record["j_call"]]
    v_first, j_codon = mapped["V"].get(v.anchor), [mapped["J"].get(j.anchor + k) for k in range(3)]
    if not record["junction"]:
        assert v_first is None or j_codon[0] is None, record["sequence_id"]
        assert record["junction_start"] is None and "anchor_removed" in record["productivity_reasons"]
        return "undefined"
    assert v_first is not None and j_codon[0] is not None, record["sequence_id"]
    assert record["junction_start"] == v_first
    intact = None not in j_codon and j_codon[2] - j_codon[0] == 2
    if intact:
        assert record["junction_end"] == j_codon[2] + 1
    # The reported junction end when an indel lies inside the J landmark codon is not asserted here.
    assert record["junction"] == record["sequence"][record["junction_start"]:record["junction_end"]]
    assert record["junction_length"] == len(record["junction"])
    if record["junction_length"] % 3 == 0:
        assert record["junction_aa"] == translate(record["junction"]), record["sequence_id"]
    # An out-of-frame junction's amino acids are left as the record states them.
    return "intact" if intact else "landmark_codon_edited"


def check_productivity_without_indels(record, alleles: dict) -> None:
    """canonical_observed_vdj_v1 recomputed from the sequence: the junction joins V and J in frame, no complete
    codon of the V-to-J coding frame is a stop, the V landmark reads Cys and the J landmark Phe or Trp."""
    assert record["junction"] and record["junction_length"] is not None, record["sequence_id"]   # both landmarks kept
    v = alleles[record["v_call"]]
    first = record["v_sequence_start"] + v.anchor % 3                      # V coding phase = anchor % 3
    stop = "*" in translate(record["sequence"][first:record["j_sequence_end"]])
    in_frame = record["junction_length"] % 3 == 0
    residues = translate(record["junction"])
    assert (record["vj_in_frame"], record["stop_codon"]) == (in_frame, stop), record["sequence_id"]
    expected = in_frame and not stop and residues[:1] == "C" and residues[-1:] in ("F", "W")
    assert record["productive"] is expected, (record["sequence_id"], record["productivity_reasons"])
    assert (record["productivity_reasons"] == ["accepted"]) is expected


def check_record(record, alleles: dict, segments: str, outcome, native, *, indels: bool) -> dict:
    """Every truth check shared by the scenarios; returns the walks by segment."""
    check_calls(record, alleles, segments, outcome, native)
    assert record["sequence_alignment"].replace("-", "") == record["sequence"]
    assert len(record["sequence_alignment"]) == len(record["germline_alignment"])
    assert all(record[name] == 0 for name in P_LENGTHS)                    # P additions are off by default
    assert record["d_inverted"] is False
    assert (record["productivity_policy"], record["productivity_scope"]) == (
        "canonical_observed_vdj_v1", "observed_coding_span")
    assert record["productive"] in (True, False)                           # never unknown: landmarks are known
    first, last = segments[0], segments[-1]                                # no end loss in any scenario here
    assert (record[f"{first}_sequence_start"], record[f"{last}_sequence_end"]) == (0, len(record["sequence"]))
    assert record["sequence_length"] == len(record["sequence"])
    walks = {s.upper(): walk(record, s, alleles[record[f"{s}_call"]]) for s in segments}
    for s in segments:
        edits = walks[s.upper()].edits
        assert record[f"n_{s}_indels"] == edits["I"] + edits["D"]          # every edit is in the alignment
    assert record["n_indels"] >= sum(record[f"n_{s}_indels"] for s in segments)
    assert record["n_v_unannotated_mutations"] == record["n_v_mutations"]
    assert all(record[name] == 0 for name in V_REGIONS)                    # no V subregions are projected
    if not indels:
        assert record["n_indels"] == 0
        check_productivity_without_indels(record, alleles)
    elif record["productive"]:
        residues = translate(record["junction"])
        assert record["junction_length"] % 3 == 0 and residues[0] == "C" and residues[-1] in "FW"
        assert "*" not in residues
    return walks


def check_np(record, segments: str) -> None:
    """Without indels, each NP region is exactly the bases between its segments, 0-6 of them (the default)."""
    order = [s for s in "vdj" if s in segments]
    for i, (a, b) in enumerate(zip(order, order[1:]), start=1):
        observed = record["sequence"][record[f"{a}_sequence_end"]:record[f"{b}_sequence_start"]]
        assert observed == record[f"np{i}"] and len(observed) == record[f"np{i}_length"], record["sequence_id"]
        assert 0 <= len(observed) <= 6 and set(observed) <= set("ACGT")


def shm_events(outcome) -> list:
    """(sequence position, new base) for each somatic hypermutation event, in order, from the trace."""
    found = collections.defaultdict(dict)
    for choice in outcome.trace().prefix_query("mutate.s5f."):
        match = re.fullmatch(r"mutate\.s5f\.(site|base)\[([0-9]+)\]", choice.address)
        if match:
            found[int(match.group(2))][match.group(1)] = choice.value
    return [(found[i]["site"], found[i]["base"].decode("ascii")) for i in sorted(found)]


def check_shm(record, outcome, walks: dict, alleles: dict) -> bool:
    """Without indels, every mutation event the trace records is counted in the segment it hit, the last event
    at a site sets its base, and the bases that differ from the generating allele's source sequence are exactly
    the hit sites whose final base is not the source base.  Returns whether the V recall names another allele;
    if it does, that allele's DNA differs and the germline alignment is still the generating allele's."""
    events = shm_events(outcome)
    assert len(events) == record["n_mutations"] > 0
    regions = {s: (record[f"{s}_sequence_start"], record[f"{s}_sequence_end"]) for s in "vdj"}

    def segment_of(site):
        return next((s for s, (a, b) in regions.items() if a <= site < b), "np")

    hits = collections.Counter(segment_of(site) for site, _ in events)
    assert [record[f"n_{s}_mutations"] for s in ("v", "d", "j", "np")] == [hits[s] for s in ("v", "d", "j", "np")]
    final = dict(events)
    assert all(record["sequence"][site] == base for site, base in final.items())
    expected = set()
    for site, base in final.items():
        s = segment_of(site)
        if s != "np" and base != alleles[record[f"{s}_call"]].sequence[site - regions[s][0]]:
            expected.add(site)
    assert {q for w in walks.values() for q in w.differs} == expected
    call, recalled = record["v_call"], record["recalled_v_call"].split(",")
    if call in recalled:
        return False
    generating = alleles[call].sequence
    assert all(alleles[name].sequence != generating for name in recalled)
    assert record["germline_alignment"][record["v_alignment_start"]:record["v_alignment_end"]] == generating
    return True


def indel_events(outcome) -> tuple:
    """(drawn count, [(kind "I" or "D", sequence position when the event happened, inserted base or None)]) for
    the polymerase indel step, in order, from the trace."""
    count, found = None, collections.defaultdict(dict)
    for choice in outcome.trace().prefix_query("corrupt.indel."):
        if choice.address == "corrupt.indel.count":
            count = choice.value
            continue
        match = re.fullmatch(r"corrupt\.indel\.(kind|site|base)\[([0-9]+)\]", choice.address)
        assert match, choice.address
        found[int(match.group(2))][match.group(1)] = choice.value
    events = []
    for i in sorted(found):
        event = found[i]
        assert ("base" in event) is event["kind"]                         # an insertion draws its base
        events.append(("I", event["site"], event["base"].decode("ascii")) if event["kind"]
                      else ("D", event["site"], None))
    return count, events


def check_indels(record, outcome, walks: dict, segments: str, low: int, high: int) -> None:
    """The trace's indel events, replayed in order over the bases present before the indel step, are exactly
    the record's alignment edits: each inserted base sits at an I column (with its drawn base), each deleted
    segment base is a D column at its germline position, every aligned base is the base that preceded the
    indel step at that germline position, and the number of events is within the configured count range."""
    count, events = indel_events(outcome)
    assert low <= count <= high and count == len(events) == record["n_indels"], (count, record["n_indels"])
    insertions = sum(kind == "I" for kind, _, _ in events)
    before = len(record["sequence"]) - insertions + (len(events) - insertions)
    tokens, removed = list(range(before)), []               # positions before the indel step; inserted bases
    for kind, site, base in events:
        if kind == "I":
            tokens.insert(site, ("inserted", base))
        else:
            removed.append(tokens.pop(site))
    assert len(tokens) == len(record["sequence"])
    new = [q for q, token in enumerate(tokens) if isinstance(token, tuple)]
    assert all(record["sequence"][q] == tokens[q][1] for q in new)
    offsets = {}
    for s in segments:
        w = walks[s.upper()]
        start, end = record[f"{s}_sequence_start"], record[f"{s}_sequence_end"]
        assert w.inserted == [q for q in new if start <= q < end], s
        assert all(isinstance(tokens[q], int) for q in w.mapped.values()), s
        layout = {tokens[q] - g for g, q in w.mapped.items()}
        assert len(layout) == 1, s                                      # one layout before the indel step
        (offsets[s],) = layout
        length = record[f"{s}_germline_end"]
        assert w.deleted == sorted(i - offsets[s] for i in removed
                                   if isinstance(i, int) and 0 <= i - offsets[s] < length), s
    assert offsets[segments[0]] == 0 and offsets[segments[-1]] + record[f"{segments[-1]}_germline_end"] == before


# ── the IGK VJ build: pinned IGK release 5, a mechanical recipe ──────


def igk_exports(directory=None, names=None) -> dict:
    """The pinned IGK release 5 exports, or copies of them in ``directory`` under ``names``."""
    exports = {fmt: DATA / "ogrdb" / IGK[0] / str(IGK[1]) / f"{fmt}.{SUFFIX[fmt]}" for fmt in ALL}
    if directory is None:
        return exports
    directory.mkdir(parents=True)
    copies = {}
    for fmt, name in zip(ALL, names):
        copies[fmt] = directory / name
        shutil.copyfile(exports[fmt], copies[fmt])
    return copies


def igk_rule() -> tuple:
    """The recipe rule, applied to the raw exports with the standard library.

    Every V and J record of the 'airr' export is selected from its own record, V before J and in export
    order within a segment, unless
      1. its label has more than one record in that export (no first-record choice);
      2. its label has an 'airr_ex' record with another coding sequence (the exports disagree, and choosing
         either would need an explicit record choice);
      3. its landmark codon does not encode the residue canonical_productive needs, which the recipe requires:
         Cys at cdr3_start - 4 for V (IMGT 104), Phe at j_cdr3_end - 1 for an IGK J (IMGT 118).
    The provider's 'functional' field is deliberately not a selection criterion: the profile cannot classify it
    and the runtime claims no functionality, so alleles the provider marks false are selected like any other.
    Returns (selection entries, {excluded label: reason})."""
    records, extended = raw_records(IGK, "airr"), raw_records(IGK, "airr_ex")
    landmarks = gapped_landmarks(IGK)
    per_label = collections.Counter(record["label"] for record in records)
    extended_dna = collections.defaultdict(set)
    for record in extended:
        extended_dna[record["label"]].add(record["coding_sequence"])
    selected, excluded = {"V": [], "J": []}, {}
    for index, record in enumerate(records):
        allele = source_allele(IGK, index, landmarks)
        codon = allele.sequence[allele.anchor:allele.anchor + 3]
        if per_label[allele.label] != 1:
            excluded[allele.label] = "more than one 'airr' record"
        elif extended_dna[allele.label] - {allele.sequence}:
            excluded[allele.label] = "the 'airr_ex' coding sequence differs"
        elif CODONS[codon] != {"V": "C", "J": "F"}[allele.segment]:
            excluded[allele.label] = f"landmark codon {codon}"
        else:
            selected[allele.segment].append(
                {"segment": allele.segment, "label": allele.label,
                 "record": {"asset": "airr", "pointer": f"{RECORDS}/{index}"},
                 "sequence_sha256": sha(allele.sequence),
                 "reason": "selected by the test rule: one 'airr' record, both exports agree, canonical landmark"})
    return selected["V"] + selected["J"], excluded


IGK_REQUIRES = (("observed_span_simulation", "assertion"), ("source_annotated_junction", "source"),
                ("canonical_productive", "assertion"))


def igk_recipe(lock, selection, requires=IGK_REQUIRES) -> dict:
    identity = lock.identity
    return {"artifact": "genairr_reference_recipe", "schema_version": 1, "source_lock_hash": lock.commitment,
            "adapter": lock.adapter.to_dict(), "projection": {"id": "observed_span", "version": 1},
            "reference": {"name": "human IGK test reference (OGRDB IGKappa_VJ release 5)", "version": "1",
                          "taxon": identity.taxon.to_dict(), "locus": identity.locus,
                          "chain_type": identity.chain_type, "build_date": "2026-10-28"},
            "selection": list(selection), "decisions": [],
            "requires": [{"capability": c, "evidence": e} for c, e in requires],
            "defaults_profile": "genairr.unfitted_defaults.v1"}


@pytest.fixture(scope="module")
def igk(tmp_path_factory):
    """(source, build, recipe document, selection, excluded) for the IGK VJ test reference."""
    with no_network():
        source = OgrdbSource.from_files(set_id=IGK[0], release=IGK[1], profile=PROFILE, exports=igk_exports(),
                                        cache=tmp_path_factory.mktemp("igk") / "cache")
        selection, excluded = igk_rule()
        document = igk_recipe(source.lock, selection)
        return source, source.build(document), document, selection, excluded


def test_an_igk_vj_build_from_a_mechanical_recipe_accounts_for_every_source_record(igk):
    _, build, _, selection, excluded = igk
    stated = EXPECT["sets"][IGK_KEY]
    # The labels the rule excludes are the fixtures' disagreeing export pair and non-Cys V landmarks.
    non_cys = EVIDENCE["conventions"]["v_landmark"]["residues_other_than_cys_at_imgt_104"]["IGKappa_VJ.5"]
    assert non_cys == {"IGKV2-29*01": "TGA", "IGKV1D-37*01": "GGT"}
    assert stated["airr_ex_coding_suffix_beyond_airr"] == {"IGKV1D-8*04": "TCC"}
    assert excluded == {"IGKV1D-8*04": "the 'airr_ex' coding sequence differs",
                        **{label: f"landmark codon {codon}" for label, codon in non_cys.items()}}
    assert [entry["segment"] for entry in selection].count("V") == stated["airr"]["V"] - len(excluded)
    assert [entry["segment"] for entry in selection].count("J") == stated["airr"]["J"]

    reference, report = build.reference, build.report
    records = raw_records(IGK, "airr")
    assert [(a.segment, a.label, a.sequence_sha256) for a in reference.alleles] == [
        (e["segment"], e["label"], e["sequence_sha256"]) for e in selection]
    usage = json.loads(reference.to_json())["models"]["allele_usage"]
    assert usage["V"] == [e["label"] for e in selection if e["segment"] == "V"]   # one weight per identity
    assert usage["J"] == [e["label"] for e in selection if e["segment"] == "J"]
    assert "D" not in usage
    alleles = {a.label: a for a in reference.alleles}
    for group in stated["airr_identical_coding_groups"]:                          # equal DNA, distinct identities
        assert len({alleles[label].sequence_sha256 for label in group}) == 1 and len(group) == 2

    index = {record["label"]: i for i, record in enumerate(records)}
    chosen = {e["label"] for e in selection}
    expected = {}
    for fmt in ("airr", "airr_ex"):
        for i, record in enumerate(raw_records(IGK, fmt)):
            label, at = record["label"], f"{fmt}#{RECORDS}/{i}"
            if label in chosen and fmt == "airr":
                expected[at] = ("selected", f"{record['sequence_type']}:{label}", None, None)
            elif label in chosen:
                assert record["coding_sequence"] == records[index[label]]["coding_sequence"]
                expected[at] = ("excluded", None, "duplicate_occurrence", f"airr#{RECORDS}/{index[label]}")
            elif label in non_cys:                       # canonical_productive is required
                expected[at] = ("excluded", None, "missing_required_evidence", None)
            else:                                        # the two exports disagree and neither is chosen
                expected[at] = ("excluded", None, "unresolved_conflicting_sequence", None)
    assert len(expected) == stated["airr"]["records"] + stated["airr_ex"]["records"]
    assert [o.record.id for o in report.records] == list(expected)
    for outcome in report.records:
        related = None if outcome.related is None else outcome.related.id
        assert (outcome.outcome, outcome.allele, outcome.exclusion, related) == expected[outcome.record.id]
        label = report.inventory.record(outcome.record).label
        if outcome.exclusion == "missing_required_evidence":
            codon = non_cys[label]
            assert f"IMGT 104 codon {codon} encodes {CODONS[codon]}, not C" in outcome.detail
        if outcome.exclusion == "unresolved_conflicting_sequence":
            assert "extends the airr one at the 3' end by 3 nt (TCC)" in outcome.detail
    counts = collections.Counter(o.outcome for o in report.records)
    assert counts == {"selected": len(selection), "excluded": len(expected) - len(selection)}
    for export in report.counts.exports:
        assert (export.record_occurrences.V, export.record_occurrences.J) == (
            stated[export.asset]["V"], stated[export.asset]["J"])


def test_the_igk_build_states_its_capabilities_and_keeps_unknowns_null(igk):
    _, build, _, selection, _ = igk
    landmarks = gapped_landmarks(IGK)
    records = raw_records(IGK, "airr")
    # The provider marks three V false; the rule ignores that field, so the two with a Cys landmark are selected.
    provider_false = {record["label"] for record in records if record["functional"] is False}
    assert provider_false == {"IGKV1-13*01", "IGKV2-29*01", "IGKV3D-15*02"}
    alleles = {a.label: a for a in build.reference.alleles}
    assert provider_false & set(alleles) == {"IGKV1-13*01", "IGKV3D-15*02"}
    assert all(alleles[label].functionality is None for label in ("IGKV1-13*01", "IGKV3D-15*02"))
    pointer = {e["label"]: int(e["record"]["pointer"].rsplit("/", 1)[1]) for e in selection}
    for capabilities in (build.reference.capabilities(), build.report.capabilities):
        status = capabilities.status
        assert (status("source_annotated_junction").available, status("source_annotated_junction").evidence) == (
            True, "source")
        assert (status("canonical_productive").available, status("canonical_productive").evidence) == (
            True, "assertion")
        assert capabilities.phase_support == "anchor_derived"
        for name, count in (("original_coding_ends", len(selection)), ("complete_reference_span", 2 * len(selection))):
            assert not status(name).available and len(status(name).unmet) == count, name   # V 3' + J 5'; both ends
            assert all('"complete"' in u.requirement and u.hint for u in status(name).unmet)
        for name in ("genotype_constrained", "duplex_family_replay"):
            assert (status(name).available, status(name).note) == (False, DEFERRED)
    assert [r.to_dict() for r in build.report.requirements] == [
        {"capability": c, "evidence": e, "satisfied": True} for c, e in IGK_REQUIRES]
    for allele in build.reference.alleles:
        assert (allele.functionality, allele.boundary_5, allele.boundary_3) == (None, None, None)
        source = source_allele(IGK, pointer[allele.label], landmarks)
        assert source.label == records[pointer[allele.label]]["label"] == allele.label
        assert (allele.anchor.value, allele.anchor.support) == (source.anchor, "source")
        assert (allele.coding_phase.value, allele.coding_phase.derivation) == (source.anchor % 3, "anchor_derived")


def test_an_igk_recipe_requiring_original_coding_ends_fails_with_a_capability_error(igk):
    source, _, _, selection, _ = igk
    requires = IGK_REQUIRES + (("original_coding_ends", "assertion"),)
    with pytest.raises(ReferenceSourceError) as info:
        source.build(igk_recipe(source.lock, selection, requires))
    assert (info.value.kind, info.value.path) == ("capability", "/requires/3")
    assert "1 of 4 required capabilities are not met" in info.value.message
    assert "original_coding_ends" in info.value.message
    assert all(f"{e['segment']}:{e['label']}" in info.value.message for e in selection)


def location_fragments(*places) -> set:
    """The text forms of ``places`` (paths in both separators, or names) that would reveal a location.  A
    fragment shorter than 8 characters, or a bare root or drive such as '/' or 'C:\\', matches unrelated bytes
    and reveals nothing, so it is skipped."""
    fragments = set()
    for place in places:
        for form in ({str(place), place.as_posix()} if isinstance(place, pathlib.PurePath) else {place}):
            if len(form) >= 8 and form != pathlib.PurePath(form).anchor:
                fragments.add(form)
    return fragments


def assert_names_no_location(files: dict, *places) -> None:
    for name, data in files.items():
        for fragment in location_fragments(*places):
            assert fragment.encode("utf-8") not in data, (name, fragment)


def written_members(build, directory: pathlib.Path) -> dict:
    """The manifest members and the manifest of ``build`` written as a dossier into ``directory``."""
    build.write_dossier(directory)
    return {name: (directory / name).read_bytes()
            for name in [member.name for member in build.manifest.members] + ["manifest.json"]}


def test_a_moved_and_renamed_igk_snapshot_builds_the_same_reference_and_names_no_location(igk, tmp_path):
    source, build, document, _, _ = igk
    names = [f"renamed-export-{n}.bin" for n in range(len(ALL))]
    moved = OgrdbSource.from_files(set_id=IGK[0], release=IGK[1], profile=PROFILE,
                                   exports=igk_exports(tmp_path / "moved snapshot" / "elsewhere", names),
                                   cache=tmp_path / "another cache")
    assert moved.lock == source.lock
    again = moved.build(document)
    assert again.manifest == build.manifest
    assert again.report.to_dict() == build.report.to_dict()
    assert again.reference.artifact_hash == build.reference.artifact_hash
    first, second = written_members(build, tmp_path / "first"), written_members(again, tmp_path / "second")
    assert first == second
    assert_names_no_location(second, tmp_path, REPO, DATA / "ogrdb", pathlib.Path.home(), "moved snapshot",
                             "another cache", "renamed-export-", "airr.json")


def test_the_location_check_skips_a_home_directory_at_the_filesystem_root(igk, tmp_path, monkeypatch):
    _, build, _, _, _ = igk
    root = pathlib.Path(tmp_path.anchor)                       # '/' here; a drive root such as 'C:\\' on Windows
    monkeypatch.setenv("HOME", str(root))
    monkeypatch.setenv("USERPROFILE", str(root))
    assert pathlib.Path.home() == root
    assert location_fragments(pathlib.Path.home(), "/", "C:\\", "short") == set()
    assert_names_no_location(written_members(build, tmp_path / "dossier"), tmp_path, REPO, pathlib.Path.home())


def test_igk_simulation_is_validated_against_the_selected_source_alleles(igk):
    _, build, _, selection, _ = igk
    reference = build.reference
    landmarks = gapped_landmarks(IGK)
    alleles = {e["label"]: source_allele(IGK, int(e["record"]["pointer"].rsplit("/", 1)[1]), landmarks)
               for e in selection}
    result = Experiment.on(reference).recombine().run_records(n=200, seed=17, validate_records=True)
    assert result.reference_content_hash == reference.execution_hash
    assert result.validate_records(reference).ok
    assert len(result) == 200
    shapes = collections.Counter()
    for outcome, record in zip(result.outcomes, result):
        walks = check_record(record, alleles, "vj", outcome, reference.native_refdata(), indels=False)
        assert all(not w.differs for w in walks.values())                     # no mutation step: exact source DNA
        v, j = alleles[record["v_call"]], alleles[record["j_call"]]
        assert CODONS[v.sequence[v.anchor:v.anchor + 3]] == "C" and CODONS[j.sequence[j.anchor:j.anchor + 3]] == "F"
        shapes[check_junction(record, {s: walks[s].mapped for s in "VJ"}, alleles)] += 1
        check_np(record, "vj")
        assert record["v_call"] in record["recalled_v_call"].split(",")
        assert record["j_call"] in record["recalled_j_call"].split(",")
    assert shapes == {"intact": 200}
    assert sum(record["productive"] for record in result) > 0
    seen = {record["v_call"] for record in result} | {record["j_call"] for record in result}
    assert len(seen & {e["label"] for e in selection if e["segment"] == "J"}) == 7


# ── the packaged human IGH starter: generation truth ─────────────────


@pytest.fixture(scope="module")
def starter():
    with no_network():
        return ReferenceBuild.bundled(STARTER_ID)


@pytest.fixture(scope="module")
def igh_alleles():
    return starter_alleles()


def identical_dna_siblings() -> dict:
    siblings = {}
    for group in SELECTION["identical_coding_sequence_groups_selected"]:
        for label in group:
            siblings[label] = set(group)
    return siblings


def test_ordinary_igh_generation_uses_the_exact_selected_source_alleles(starter, igh_alleles):
    reference = starter.reference
    assert [(a.segment, a.label, a.sequence_sha256) for a in reference.alleles] == [
        (a.segment, a.label, sha(a.sequence)) for a in igh_alleles.values()]
    assert all(a.sequence[a.anchor:a.anchor + 3] == {"V": "TGT", "J": "TGG"}[a.segment]
               for a in igh_alleles.values() if a.segment != "D")
    result = Experiment.on(reference).recombine().run_records(n=300, seed=101, validate_records=True)
    assert result.reference_content_hash == reference.execution_hash
    assert result.validate_records(reference).ok
    siblings, shapes, listed = identical_dna_siblings(), collections.Counter(), collections.Counter()
    for outcome, record in zip(result.outcomes, result):
        walks = check_record(record, igh_alleles, "vdj", outcome, reference.native_refdata(), indels=False)
        assert all(not w.differs for w in walks.values())
        shapes[check_junction(record, {s: walks[s].mapped for s in "VJ"}, igh_alleles)] += 1
        check_np(record, "vdj")
        for s in "vdj":
            call, recalled = record[f"{s}_call"], set(record[f"recalled_{s}_call"].split(","))
            assert call in recalled                                   # unmutated: the source allele is recalled
            if call in siblings:                                      # and so is each identical-DNA identity
                assert siblings[call] <= recalled, (call, recalled)
                listed[s] += 1
    assert shapes == {"intact": 300}
    assert listed["v"] > 0 and listed["d"] > 0
    assert 0 < sum(record["productive"] for record in result) < 300


def test_shm_is_recorded_against_the_generating_allele_and_never_against_a_recall(starter, igh_alleles):
    reference = starter.reference
    result = Experiment.on(reference).recombine().mutate(rate=0.1).run_records(n=300, seed=202,
                                                                                validate_records=True)
    assert result.reference_content_hash == reference.execution_hash
    assert result.validate_records(reference).ok
    recalled_elsewhere = 0
    for outcome, record in zip(result.outcomes, result):
        walks = check_record(record, igh_alleles, "vdj", outcome, reference.native_refdata(), indels=False)
        assert check_junction(record, {s: walks[s].mapped for s in "VJ"}, igh_alleles) == "intact"
        recalled_elsewhere += check_shm(record, outcome, walks, igh_alleles)
    assert recalled_elsewhere > 0, "exercise records whose V recall is another allele"


def test_polymerase_indels_keep_every_alignment_edit_against_the_generating_allele(starter, igh_alleles):
    reference = starter.reference
    result = (Experiment.on(reference).recombine().mutate(count=(5, 15)).polymerase_indels(count=(1, 3))
              .run_records(n=300, seed=303, validate_records=True))
    assert result.reference_content_hash == reference.execution_hash
    assert result.validate_records(reference).ok
    edits, shapes = collections.Counter(), collections.Counter()
    for outcome, record in zip(result.outcomes, result):
        walks = check_record(record, igh_alleles, "vdj", outcome, reference.native_refdata(), indels=True)
        check_indels(record, outcome, walks, "vdj", 1, 3)
        for s in "VDJ":
            assert len(walks[s].differs) <= record[f"n_{s.lower()}_mutations"]  # only SHM changes an aligned base
            edits.update({(s, op): n for op, n in walks[s].edits.items() if op in "ID"})
        shapes[check_junction(record, {s: walks[s].mapped for s in "VJ"}, igh_alleles)] += 1
    assert sum(record["n_indels"] for record in result) > 300
    assert all(edits[(s, op)] > 0 for s in "VJ" for op in "ID"), edits
    assert shapes["intact"] >= 280 and set(shapes) <= {"intact", "landmark_codon_edited", "undefined"}


# ── unknowns and unsupported workflows on the packaged starter ───────


def test_the_bundled_starter_keeps_unknowns_unknown_and_names_its_blockers(starter, igh_alleles):
    reference = starter.reference
    for allele in reference.alleles:
        assert (allele.functionality, allele.boundary_5, allele.boundary_3) == (None, None, None), allele.label
        source = igh_alleles[allele.label]
        if allele.segment == "D":
            assert (allele.anchor, allele.coding_phase) == (None, None)
        else:
            assert (allele.anchor.value, allele.anchor.support) == (source.anchor, "source")
            assert (allele.coding_phase.value, allele.coding_phase.derivation) == (source.anchor % 3, "anchor_derived")
    capabilities = reference.capabilities()
    for name in ("original_coding_ends", "complete_reference_span"):
        status = capabilities.status(name)
        assert not status.available and status.unmet
        assert all('"complete"' in u.requirement and '"value": "complete"' in u.hint for u in status.unmet)
    for name in ("genotype_constrained", "duplex_family_replay"):
        assert (capabilities.status(name).available, capabilities.status(name).note) == (False, DEFERRED)
    reference.require("observed_span_simulation", "source_annotated_junction", evidence="source")
    reference.require("canonical_productive")
    for name, evidence in (("canonical_productive", "source"), ("original_coding_ends", "assertion"),
                           ("complete_reference_span", "assertion"), ("genotype_constrained", "assertion")):
        with pytest.raises(ReferenceCapabilityError) as info:
            reference.require(name, evidence=evidence)
        assert info.value.kind == "capability" and name in str(info.value)


def test_workflows_portable_references_do_not_support_still_reject_the_bundled_starter(starter, monkeypatch):
    import GenAIRR._s5f_loader as loader
    reference = starter.reference
    genotype = GenAIRR.Genotype.from_dataconfig(GenAIRR.data.HUMAN_IGH_OGRDB)
    for attempt in (lambda: Experiment.on(reference).recombine().with_genotype(genotype),
                    lambda: Experiment.on(reference).recombine().run_cohort([genotype])):
        with pytest.raises(ReferenceCapabilityError, match=DEFERRED) as info:
            attempt()
        assert info.value.kind == "capability"
    compiled = Experiment.on(reference).recombine()
    compiled._genotype = genotype._snapshot()                      # a genotype reaching compilation by any route
    with pytest.raises(ReferenceCapabilityError, match=DEFERRED):
        compiled.compile()
    experiment = Experiment.on(reference).recombine()
    for policy in ("raw", "functional_anchors_only", "functional_status"):
        with pytest.raises(ReferenceCapabilityError) as info:
            experiment.curate_refdata(policy)
        assert info.value.kind == "capability"
        assert "allele selection is not yet supported for portable references" in str(info.value)
    assert experiment.refdata is reference.native_refdata()
    monkeypatch.setattr(loader, "load_builtin_s5f_kernel", lambda *a, **k: pytest.fail("sample preparation reached"))
    for run in (lambda e: e.run_records(seed=7, capture_family_replay=True),
                lambda e: e.compile().run_records(seed=7, capture_family_replay=True)):
        with pytest.raises(FamilyReplayError, match="reference_spec_v1") as info:
            run(Experiment.on(reference).recombine().clonal_lineage(**LINEAGE))
        assert DEFERRED in str(info.value)
