"""Independent literal ancestry checks across saved/reloaded lineage bundles.

Uses only GenAIRR and the standard library. These short authored controls
establish replay/model consistency, not measured biochemical frequencies.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import itertools
import json
from pathlib import Path
import tempfile

from GenAIRR import Experiment, _engine
from GenAIRR.experimental import FamilyReplayBundle, FamilyReplayLimits
from GenAIRR._family_replay import _content_decode
from .duplex_family_sweep import source, check_record


def _packed(value):
    return _engine._family_replay_pack(value)


def run():
    counts = Counter()
    plans = []
    controls = [(order, inverted, strict, productive, edited, "neutral")
                for order, inverted, strict, productive, edited in itertools.product(
                    (("vd", "dj"), ("dj", "vd")), (False, True), (False, True),
                    (False, True), (False, True))]
    controls += [(('vd', 'dj'), inverted, strict, True, edited, target)
                 for inverted, strict, edited, target in itertools.product(
                     (False, True), (False, True), (False, True), ("explicit", "automatic"))]
    controls = [profile + ("ordinary",) for profile in controls]
    controls += [(('vd', 'dj'), False, strict, True, True, "neutral", growth)
                 for strict, growth in itertools.product((False, True), ("retry", "extinct"))]
    limits = FamilyReplayLimits()
    extension = Path(_engine.__file__)
    binary_hash = hashlib.sha256(extension.read_bytes()).hexdigest()
    with tempfile.TemporaryDirectory(prefix="genairr-durable-sweep-") as directory:
        path = Path(directory) / "family.json"
        for index, (order, inverted, strict, productive, edited, target, growth) in enumerate(controls):
            ref, model, expected = source(inverted, order)
            experiment = Experiment.on(ref).recombine(duplex_repair=model)
            if inverted:
                experiment.invert_d(prob=1)
            options = {"selection_strength": 0.0}
            if target != "neutral":
                options = {"selection_strength": 0.75, "mature_substitutions": 3}
                if target == "explicit":
                    options["target_aa"] = "KCSFK" if inverted else "KCTFK"
            experiment.clonal_lineage(n_clones=2, max_generations=3, n_max=18,
                                      n_sample=6, rate=0.12,
                                      lambda_base={"ordinary": 2.5, "retry": 0.8, "extinct": 0}[growth],
                                      allow_extinction=growth == "extinct", **options)
            if edited:
                (experiment.polymerase_indels(count=2).pcr_amplify(count=1)
                 .sequencing_errors(count=1).end_loss_5prime(length=1)
                 .end_loss_3prime(length=1).random_strand_orientation(prob=1))
            if productive:
                experiment.productive_only()
            experiment.with_metadata(unfitted_uniform=True,
                                     control={"negative_zero": -0.0, "literal_hex": "8000000000000000"})
            original = experiment.run_records(seed=23, strict=strict, validate_records=True,
                                               capture_family_replay=True)
            uncaptured = experiment.run_records(seed=23, strict=strict, validate_records=True)
            assert _packed(original.records) == _packed(uncaptured.records)
            original.save_family_replay(path)
            loaded = FamilyReplayBundle.read_from(path)
            replayed = loaded.replay_records()
            assert _packed(original.records) == _packed(replayed.records)
            assert replayed.validate_records(loaded.refdata).ok
            assert replayed.validate_families_with_parents(loaded.refdata).ok
            assert [tree.to_newick() for tree in replayed.lineage_trees] == [tree.to_newick() for tree in original.lineage_trees]
            assert [tree.to_node_table_tsv() for tree in replayed.lineage_trees] == [tree.to_node_table_tsv() for tree in original.lineage_trees]
            assert sum(record["duplicate_count"] for record in replayed) == (0 if growth == "extinct" else 12)
            for record, outcome in zip(replayed, replayed.outcomes, strict=True):
                founder = replayed.parents[record["parent_id"]].airr_record(loaded.refdata)["duplex_repair"]
                check_record(record, outcome, expected, inverted, founder)
                if productive:
                    assert record["productive"] is True
                counts["literal_record_checks"] += 1
            for clone, (old, new, tree) in enumerate(zip(original.families, replayed.families,
                                                        replayed.lineage_trees, strict=True)):
                for node in tree.nodes():
                    a = old.replay_node(node.id, ref).airr_record(ref, sequence_id="same")
                    b = new.replay_node(node.id, loaded.refdata).airr_record(loaded.refdata, sequence_id="same")
                    assert _packed(a) == _packed(b)
                    check_record(b, new.replay_node(node.id, loaded.refdata), expected, inverted)
                    counts["node_checks_including_roots_unsampled"] += 1
            wire = json.loads(loaded.to_json())
            payload = _content_decode(wire["root"], wire["content"], limits,
                                      _engine.FamilyReplayBudget(limits._json()))
            attempt_count = sum(len(clone["attempts"]) for clone in payload["clones"])
            counts["attempts_including_rejected"] += attempt_count
            counts["rejected_attempts"] += attempt_count - 2
            if growth == "retry":
                assert attempt_count > 2
            if growth == "extinct":
                assert len(replayed) == 0
                counts["accepted_extinct_families"] += len(replayed.families)
            for clone in payload["clones"]:
                for attempt in clone["attempts"]:
                    counts["growth_parent_proposals"] += len(attempt["growth"])
                    counts["ordered_sample_choices"] += len(attempt["samples"])
            counts["plans"] += 1
            counts["founders"] += len(original.parents)
            counts["artifact_bytes"] += path.stat().st_size
            plans.append({"join_order": order, "d_inverted": inverted, "strict": strict,
                          "productive_only": productive, "observation_suffix": edited,
                          "target": target, "records": len(replayed),
                          "growth_control": growth,
                          "attempts": attempt_count, "bytes": path.stat().st_size})
            path.unlink()
    assert hashlib.sha256(extension.read_bytes()).hexdigest() == binary_hash
    return {"schema": "durable_family_sweep_v1", "extension_sha256": binary_hash,
            "counts": dict(counts), "errors": 0, "plans": plans,
            "checks": ["capture-on/off seeded identity", "saved and loaded exact typed records",
                       "literal repair/query/edit/link/counter reconstruction", "canonical inverted D",
                       "actual cumulative mutations including reversions", "all accepted nodes including unsampled cells",
                       "tree topology, affinity and abundance", "complete attempted-family and ordered sampling receipts"],
            "scope": "52 short restored-landmark controls including retries and accepted extinction; equal authored weights and declared neutral/selection controls; no measured repair rates or full-length performance claim"}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.write_text(json.dumps(run(), indent=2, allow_nan=False) + "\n")
