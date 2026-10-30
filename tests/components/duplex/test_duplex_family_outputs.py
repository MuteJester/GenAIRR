"""Duplex founder truth survives family forks, export and branch replay."""
from collections import Counter
from copy import deepcopy
import csv
import json

import pytest

from GenAIRR import Experiment, SimulationResult
from tests.support.cases.duplex import references, model, history

pytestmark = pytest.mark.behavior


def _base():
    refs = references()
    return refs, Experiment.on(refs).recombine(duplex_repair=model(refs, [
        history("first"), history("same-dna-second", weight=1),
    ]))


def _fixed(*, strict=False):
    refs, exp = _base()
    with pytest.warns(DeprecationWarning):
        exp.expand_clones(n_clones=2, per_clone=8)
    result = (exp.mutate(model="uniform", count=2)
              .end_loss_5prime(length=1)
              .random_strand_orientation(prob=0.5)
              .run_records(seed=17, strict=strict, validate_records=True))
    return refs, result


@pytest.mark.parametrize("strict", [False, True])
def test_fixed_forks_keep_founder_evidence_and_independent_current_maps(strict):
    refs, result = _fixed(strict=strict)
    assert len(result.parents) == 2 and len(result) == 16
    assert result.validate_families_with_parents(refs).ok
    for clone in range(2):
        descendants = [r for r in result if r["clone_id"] == clone]
        founder = result.parents[clone].airr_record(refs)["duplex_repair"]
        for record in descendants:
            duplex = record["duplex_repair"]
            assert record["parent_id"] == clone
            assert duplex["product"] == founder["product"]
            assert duplex["evidence"] == founder["evidence"]
            assert duplex["selection"] == founder["selection"]
            # Mutation runs AFTER the fork here, so each descendant draws
            # its own two attempts and they legitimately differ -- which is
            # also why the sequences below differ. A uniform draw landing on
            # the base already there is not counted, so on this 15-base
            # molecule the realised count is anywhere from zero to two; two
            # of the sixteen descendants realise none, so the band is the
            # strongest statement available per record. The recoverable count
            # is pinned exactly instead, against the alignment pair.
            assert record["n_mutations"] <= 2
            assert record["n_observable_mutations"] == sum(
                1 for observed, germline in zip(record["sequence_alignment"],
                                                record["germline_alignment"])
                if observed.upper() != germline.upper())
            assert record["end_loss_5_length"] == 1
            assert duplex["current_paths"]["v"]["reference_start"] == 1
        assert len({r["sequence"] for r in descendants}) > 1
    # Python output dictionaries are independent as well as native maps.
    other = deepcopy(result[1]["duplex_repair"])
    result[0]["duplex_repair"]["observation"]["links"].clear()
    assert result[1]["duplex_repair"] == other
    assert not result.validate_records(refs).ok


@pytest.mark.parametrize("tamper", ["ledger", "selection", "missing", "integer-type"])
def test_record_only_and_parent_validators_detect_changed_founder(tamper):
    refs, result = _fixed()
    record = result[0]
    if tamper == "ledger":
        record["duplex_repair"]["product"]["ledger"][0]["base"] = "G"
    elif tamper == "selection":
        record["duplex_repair"]["selection"]["history_id"] = "forged"
    elif tamper == "missing":
        del record["duplex_repair"]
    else:
        record["duplex_repair"]["selection"]["row_index"] = False
    records_only = SimulationResult(result.records)
    assert not records_only.validate_families().ok
    assert not result.validate_families_with_parents(refs).ok


def test_parent_check_detects_unanimously_forged_or_dropped_founder():
    refs, result = _fixed()
    for record in result:
        if record["clone_id"] == 0:
            del record["duplex_repair"]
    # An absent entire family payload needs the actual parent as evidence.
    assert not result.validate_families_with_parents(refs).ok


def test_repertoire_keeps_equal_dna_with_distinct_observation_histories():
    refs, exp = _base()
    compiled = (exp.clonal_repertoire(
        n_clones=1, size_distribution="lognormal", mu=4, sigma=0,
        max_size=40,
    ).polymerase_indels(count=1, insertion_prob=0).compile())
    result = compiled.run_records(seed=9, validate_records=True)
    assert len(result.parents) == 1
    assert sum(r["duplicate_count"] for r in result) == 40
    assert result.validate_families_with_parents(refs).ok
    # Deleting different founder bases in the CCC tract produces identical DNA
    # with different removed physical ancestry. Keep those records separate.
    by_dna = {}
    for record in result:
        key = json.dumps({k: v for k, v in record.items()
                          if k not in {"sequence_id", "duplicate_count"}},
                         sort_keys=True)
        assert key not in by_dna.get(record["sequence"], set())
        by_dna.setdefault(record["sequence"], set()).add(key)
    assert any(len(keys) > 1 for keys in by_dna.values())


def test_repertoire_no_edit_shortcut_and_nullable_export(tmp_path):
    refs, exp = _base()
    result = exp.clonal_repertoire(
        n_clones=2, size_distribution="lognormal", mu=2, sigma=0, max_size=8,
    ).run_records(seed=2, validate_records=True)
    assert len(result) == 2 and len(result.parents) == 2
    assert all(r["duplicate_count"] == 7 for r in result)
    path = tmp_path / "families.tsv"
    result.to_tsv(path)
    with path.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert all(r["p_v_3_length"] == "" for r in rows)
    for original, exported in zip(result, rows):
        assert json.loads(exported["duplex_repair"]) == original["duplex_repair"]
        assert int(exported["parent_id"]) == original["parent_id"]
    assert result.validate_families_with_parents(refs).ok


@pytest.mark.parametrize("corrupt", [False, True])
@pytest.mark.parametrize("strict", [False, True])
def test_lineage_output_retains_families_founders_and_exact_node_replay(corrupt, strict):
    refs, exp = _base()
    exp.clonal_lineage(n_clones=2, max_generations=3, n_max=30,
                       n_sample=12, rate=0.05, lambda_base=2.5)
    if corrupt:
        exp.end_loss_3prime(length=1).random_strand_orientation(prob=1)
    result = exp.run_records(seed=8, strict=strict, validate_records=True)
    assert len(result.parents) == len(result.families) == len(result.lineage_trees) == 2
    assert result.validate_families_with_parents(refs).ok
    assert Counter(r["clone_id"] for r in result).keys() == {0, 1}
    for record in result:
        family = result.families[record["parent_id"]]
        node = result.lineage_trees[record["parent_id"]].nodes()[record["lineage_node_id"]]
        assert record["lineage_mutations_from_parent"] == node.mutations_from_parent
        replayed = family.replay_node(record["lineage_node_id"], refs)
        replay_record = replayed.airr_record(refs)
        assert replayed.validate_record(refs) == []
        assert replay_record["duplex_repair"]["product"] == record["duplex_repair"]["product"]
        assert replay_record["n_mutations"] == record["n_mutations"]
        if not corrupt:
            assert replay_record["sequence"] == record["sequence"]
            assert replay_record["duplex_repair"] == record["duplex_repair"]
        else:
            assert record["end_loss_3_length"] == 1
            assert record["rev_comp"] is True
    pytest.importorskip("pandas")
    assert str(result.to_dataframe()["lineage_mutations_from_parent"].dtype) == "Int64"


def test_extinct_duplex_lineage_retains_empty_tree_and_founder():
    refs, exp = _base()
    result = exp.clonal_lineage(
        n_clones=2, max_generations=1, n_max=5, n_sample=2,
        lambda_base=0, allow_extinction=True,
    ).run_records(seed=3, validate_records=True)
    assert len(result) == 0
    assert len(result.parents) == len(result.families) == len(result.lineage_trees) == 2
    assert result.validate_records(refs).ok


def test_equal_dna_from_distinct_lineage_cells_keeps_separate_branch_truth():
    refs, exp = _base()
    result = exp.clonal_lineage(
        n_clones=1, max_generations=2, n_max=30, n_sample=40,
        lambda_base=4, rate=0,
    ).run_records(seed=4, validate_records=True)
    assert len(result) > 1
    assert {r["sequence"] for r in result} == {"AAATGTCCCTTTAAA"}
    assert len({r["lineage_node_id"] for r in result}) == len(result)
    assert sum(r["duplicate_count"] for r in result) == 40
    assert result.validate_families_with_parents(refs).ok


def test_survival_retry_keeps_the_once_sampled_duplex_founder(monkeypatch):
    from GenAIRR import _engine

    refs, exp = _base()
    compiled = exp.clonal_lineage(
        n_clones=1, max_generations=2, n_max=20, n_sample=4,
        lambda_base=4, rate=0,
    ).compile()
    native_pre = compiled._pre
    source_seeds, growth_founders, growth_seeds = [], [], []

    class CountingSource:
        def pass_names(self):
            return native_pre.pass_names()

        def run(self, *, seed, strict):
            source_seeds.append(seed)
            return native_pre.run(seed=seed, strict=strict)

    class ForcedExtinction:
        def observed_outcomes(self):
            return []

    native_growth = _engine.simulate_family_outcomes

    def grow(founder, *args, **kwargs):
        growth_founders.append(founder.airr_record(refs)["duplex_repair"]["selection"])
        growth_seeds.append(args[9])
        if len(growth_founders) == 1:
            return ForcedExtinction()
        return native_growth(founder, *args, **kwargs)

    compiled._pre = CountingSource()
    monkeypatch.setattr(_engine, "simulate_family_outcomes", grow)
    result = compiled.run_records(seed=6, validate_records=True)
    assert source_seeds == [6]
    assert len(growth_seeds) >= 2 and growth_seeds[:2] == [6, 6 + 7919]
    assert all(history == growth_founders[0] for history in growth_founders)
    assert result.parents[0].airr_record(refs)["duplex_repair"]["selection"] == growth_founders[0]


@pytest.mark.parametrize("field,replacement", [
    ("lineage_node_id", False),
    ("lineage_parent_id", False),
    ("lineage_generation", "1"),
    ("lineage_abundance", 0),
    ("lineage_affinity", 0),  # int is not the generated float value.
    ("lineage_mutations_from_parent", False),
    ("duplicate_count", 0),
])
def test_lineage_actual_record_metadata_is_checked_against_native_node(field, replacement):
    refs, exp = _base()
    result = exp.clonal_lineage(
        n_clones=1, max_generations=1, n_max=10, n_sample=4,
        lambda_base=4, rate=0,
    ).run_records(seed=6)
    result[0][field] = replacement
    report = result.validate_records(refs)
    assert not report.ok
    assert any(issue.get("field") == field and issue["kind"] == "LineageMetadataMismatch"
               for failure in report.failures for issue in failure["issues"])


def test_equal_dna_does_not_allow_reassigning_a_record_to_another_clone_tree():
    refs = references()
    result = Experiment.on(refs).recombine(duplex_repair=model(refs)).clonal_lineage(
        n_clones=2, max_generations=1, n_max=10, n_sample=4,
        lambda_base=4, rate=0,
    ).run_records(seed=6, validate_records=True)
    assert {record["sequence"] for record in result} == {"AAATGTCCCTTTAAA"}
    first = next(record for record in result if record["clone_id"] == 0)
    other = next(record for record in result if record["clone_id"] == 1)
    for field in ("clone_id", "parent_id", "lineage_node_id"):
        first[field] = other[field]
    report = result.validate_records(refs)
    assert not report.ok
    assert {issue.get("field") for failure in report.failures for issue in failure["issues"]
            if issue["kind"] == "LineageMetadataMismatch"} >= {"clone_id", "parent_id"}


def test_duplex_family_stamps_custom_metadata():
    refs, exp = _base()
    result = (exp.clonal_repertoire(n_clones=1, max_size=1)
              .with_metadata(sample_id="control")
              .run_records(seed=0, validate_records=True))
    assert result[0]["sample_id"] == "control"
    assert result.validate_records(refs).ok


def _mutated_duplex(**passes):
    """A duplex family whose descendants carry hypermutation and nothing else."""
    refs, exp = _base()
    exp = exp.mutate(model="uniform", count=3)
    for name, kwargs in passes.items():
        exp = getattr(exp, name)(**kwargs)
    return refs, exp.run_records(n=6, seed=17)


def test_duplex_records_report_the_recoverable_mutation_count():
    """`n_observable_mutations` is in the default column order, so a duplex or
    clonal-duplex export that does not produce it writes a blank column and
    says nothing about why. The duplex plane has a mechanism for a field it
    cannot produce -- Null plus an availability reason -- and a field that
    bypasses it leaves a consumer with no way to tell "none" from "unknown".

    The expected count is derived from the alignment pair rather than from the
    duplex edit ledger the implementation reads, so the two sides do not share
    an input.
    """
    _, records = _mutated_duplex()
    total = 0
    for index, record in enumerate(records):
        assert "n_observable_mutations" in record, (
            f"record {index} omits n_observable_mutations entirely")
        expected = sum(1 for observed, germline
                       in zip(record["sequence_alignment"],
                              record["germline_alignment"])
                       if observed.upper() != germline.upper())
        assert record["n_observable_mutations"] == expected, (
            f"record {index} reports {record['n_observable_mutations']} "
            f"recoverable where the alignment shows {expected} mismatches")
        assert record["n_observable_mutations"] <= record["n_mutations"]
        total += record["n_observable_mutations"]
    assert total > 0, "nothing was mutated, so this proves nothing"


def test_duplex_recoverable_count_ignores_library_substitutions():
    """A duplex record's recoverable count is attributable to hypermutation.
    The duplex ledger attributes every write to its pass, so a PCR
    misincorporation is distinguishable from a mutation and must not be
    counted as one."""
    _, clean = _mutated_duplex()
    _, noisy = _mutated_duplex(pcr_amplify={"count": 4})
    assert [r["n_observable_mutations"] for r in noisy] == \
           [r["n_observable_mutations"] for r in clean], (
        "a PCR pass moved the recoverable hypermutation count")
    assert sum(r["n_pcr_errors"] for r in noisy) > 0, (
        "fixture error -- no PCR substitution was applied")
