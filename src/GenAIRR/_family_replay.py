"""Self-contained, source-backed duplex lineage replay bundles.

The native codecs and replay kernels own scientific semantics. This module
captures the Python family boundary, seals actual output and publishes bounded
JSON artifacts without deserializing arbitrary Python or simulation objects.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
import math
import os
from pathlib import Path
from typing import Any

from GenAIRR import _engine
from GenAIRR._artifact_io import ArtifactIOError, write_atomic


_ARTIFACT = "genairr_family_replay"
_WORKFLOW = "duplex_lineage_v1"
_CAPABILITIES = ["typed_content_dag_v1", "native_family_attempt_v1", "sealed_output_v1"]


class FamilyReplayError(ValueError):
    """A structured failure to capture, publish or replay a family bundle."""

    def __init__(self, category: str, path: str, message: str):
        self.category, self.path = category, path
        super().__init__(f"family_replay[{category}] at {path}: {message}")


@dataclass(frozen=True)
class FamilyReplayLimits:
    """Computational ceilings for capture/loading; these are not biological priors."""

    max_file_bytes: int = 256 * 1024 * 1024
    max_depth: int = 64
    max_nodes: int = 1_000_000
    max_choices: int = 5_000_000
    max_bases: int = 50_000_000
    max_work: int = 10_000_000
    max_attempts: int = 10_000
    max_metadata_bytes: int = 16 * 1024 * 1024
    max_string_bytes: int = 16 * 1024 * 1024
    max_components: int = 100_000
    max_history_entries: int = 5_000_000
    max_retained_bytes: int = 256 * 1024 * 1024
    max_values: int = 5_000_000

    def __post_init__(self):
        for name, value in asdict(self).items():
            if type(value) is not int or not 0 < value <= 2**63 - 1:
                raise FamilyReplayError("resource", f"limits.{name}", "expected a positive signed 64-bit integer")

    def _json(self):
        return _json(asdict(self))


def _limits(value):
    if value is None:
        return FamilyReplayLimits()
    if not isinstance(value, FamilyReplayLimits):
        raise TypeError("family_replay_limits must be FamilyReplayLimits")
    return value


def _json(value):
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _native(path, function, *args, category="invalid_evidence", **kwargs):
    try:
        return function(*args, **kwargs)
    except FamilyReplayError:
        raise
    except Exception as exc:
        message = str(exc)
        if any(word in message.lower() for word in ("resource", "budget exceeded", "budget exhausted", "exhaust", "limit exceeded")):
            category = "resource"
        raise FamilyReplayError(category, path, message) from exc


def _parse(text, limits, path):
    if not isinstance(text, str):
        raise FamilyReplayError("format", path, "expected JSON text")
    normalized = _native(path, _engine._family_replay_parse_json, text, limits._json(), category="format")
    return json.loads(normalized)


def _digest(value, domain):
    return _native("commitment", _engine._family_replay_digest, _json(value), domain)


def _pack(value, path, limits=None):
    return _native(path, _engine._family_replay_pack, value,
                   None if limits is None else limits._json(), category="metadata")


def _unpack(value, path, limits=None):
    return _native(path, _engine._family_replay_unpack, value,
                   None if limits is None else limits._json(), category="format")


def _charge(budget, kind, count, path):
    _native(path, getattr(budget, f"charge_{kind}"), count, category="resource")


def _value_count(value):
    """Count dynamic output values before the lossless codec copies them."""
    pending, count = [value], 0
    while pending:
        item = pending.pop()
        count += 1
        if isinstance(item, dict):
            pending.extend(item.values())
        elif isinstance(item, list):
            pending.extend(item)
    return count


def _truth_aliases(record):
    """Compatibility aliases of the authoritative generating calls."""
    for segment in ("v", "d", "j"):
        record[f"truth_{segment}_call"] = record[f"{segment}_call"]


class _OutputFence:
    """Bound cumulative dynamic rows before the next row is retained."""

    def __init__(self, limits, budget, annotation_text, annotation_count):
        self.limits, self.budget = limits, budget
        self.annotation_bytes = len(annotation_text.encode())
        self.annotation_count = annotation_count
        self.output_bytes = 32

    def record(self, record):
        _charge(self.budget, "work", _value_count(record) + self.annotation_count, "records")
        # One tagged annotation object overestimates the cost of merging its
        # entries into a row by only its object envelope and a separator.
        self.output_bytes += (len(_pack(record, "records", self.limits).encode())
                              + self.annotation_bytes + 1)
        if self.output_bytes > min(self.limits.max_file_bytes, self.limits.max_retained_bytes):
            raise FamilyReplayError("resource", "records", "cumulative output byte limit exceeded")


def _require_keys(value, required, path):
    if not isinstance(value, dict) or value.keys() != set(required):
        raise FamilyReplayError("format", path, f"expected exactly {sorted(required)}")


def _index(value, path):
    if type(value) is not int or value < 0 or value > 2**32 - 1:
        raise FamilyReplayError("format", path, "expected an unsigned 32-bit integer")
    return value


def _seed(value, path):
    if (type(value) is not str or not value or any(c not in "0123456789" for c in value)
            or (len(value) > 1 and value.startswith("0")) or len(value) > 20):
        raise FamilyReplayError("format", path, "expected a canonical u64 decimal string")
    result = int(value)
    if result > 2**64 - 1:
        raise FamilyReplayError("format", path, "seed exceeds u64")
    return result


def _stage_policy(receipt, strict, seed, path):
    if not isinstance(receipt, dict) or receipt.get("strict") is not strict:
        raise FamilyReplayError("plan_mismatch", path, "stage execution policy differs from the captured run")
    if _seed(receipt.get("seed_label"), path + ".seed_label") != seed:
        raise FamilyReplayError("invalid_scope", path + ".seed_label", "seed scheduling label mismatch")


def _reference_policy(refdata, path):
    if refdata.chain_type != "vdj" or (refdata.v_pool_size() and refdata.v_allele(0).name.upper().startswith("TR")):
        raise FamilyReplayError("unsupported", path, "durable clonal_lineage replay requires BCR VDJ reference data")


def _stage_shapes(source, post):
    names = source.pass_names()
    prefix = ["sample_allele.v", "sample_allele.d", "sample_allele.j"]
    if names not in (prefix + ["duplex.recombine"], prefix + ["invert_d", "duplex.recombine"]):
        raise FamilyReplayError("unsupported", "source.plan", "source stages must match the public duplex clonal_lineage prefix")
    allowed = {"corrupt.pcr", "corrupt.quality", "corrupt.indel", "corrupt.ns",
               "corrupt.end_loss.5", "corrupt.end_loss.3", "corrupt.rev_comp"}
    if post is not None and (not post.pass_names() or any(name not in allowed for name in post.pass_names())):
        raise FamilyReplayError("unsupported", "observations.plan", "only supported per-cell read passes may follow clonal_lineage")


def _provenance():
    from GenAIRR import __version__

    extension = getattr(_engine, "__file__", None)
    identity = None
    if extension is not None and Path(extension).is_file():
        with open(extension, "rb") as handle:
            identity = hashlib.file_digest(handle, "sha256").hexdigest() if hasattr(hashlib, "file_digest") else hashlib.sha256(handle.read()).hexdigest()
    return {"package_version": __version__, "engine_version": _engine.version(),
            "extension_sha256": identity}


def _content_encode(value, limits):
    """Intern structured JSON values; literal user dictionaries are never refs."""
    table = {}
    values = 0

    def visit(item, depth):
        nonlocal values
        values += 1
        if values > limits.max_values or depth > limits.max_depth:
            raise FamilyReplayError("resource", "content", "value/depth limit exceeded")
        if isinstance(item, dict):
            if any(type(key) is not str for key in item):
                raise FamilyReplayError("format", "content", "object keys must be strings")
            node = {"kind": "object", "items": [[key, visit(child, depth + 1)]
                    for key, child in item.items()]}
        elif isinstance(item, list):
            node = {"kind": "array", "items": [visit(child, depth + 1) for child in item]}
        else:
            if item is not None and type(item) not in (bool, int, str):
                raise FamilyReplayError("format", "content", "native values must be typed before interning")
            return {"literal": item}
        identity = _digest(node, "genairr_family_content_v1")
        table.setdefault(identity, node)
        if len(table) > limits.max_components:
            raise FamilyReplayError("resource", "content", "component limit exceeded")
        return {"ref": identity}

    root = visit(value, 0)
    return root, [{"id": identity, "value": table[identity]} for identity in sorted(table)]


def _content_decode(root, entries, limits, budget):
    if not isinstance(entries, list) or len(entries) > limits.max_components:
        raise FamilyReplayError("resource", "content", "invalid component count")
    table = {}
    for i, entry in enumerate(entries):
        path = f"content[{i}]"
        _require_keys(entry, ("id", "value"), path)
        identity, node = entry["id"], entry["value"]
        if type(identity) is not str or identity in table:
            raise FamilyReplayError("format", path, "duplicate or invalid content ID")
        _require_keys(node, ("kind", "items"), path + ".value")
        if node["kind"] not in ("object", "array") or not isinstance(node["items"], list):
            raise FamilyReplayError("format", path, "unknown content variant")
        _charge(budget, "work", len(node["items"]) + 1, path)
        if identity != _digest(node, "genairr_family_content_v1"):
            raise FamilyReplayError("integrity", path, "content commitment mismatch")
        table[identity] = node
    active, visited, costs, children = set(), set(), {}, {}

    def token_cost(token):
        if not isinstance(token, dict) or len(token) != 1:
            raise FamilyReplayError("format", "content", "invalid reference/literal token")
        if "literal" in token:
            literal = token["literal"]
            if literal is not None and type(literal) not in (bool, int, str):
                raise FamilyReplayError("format", "content", "invalid literal type")
            size = len(_json(literal).encode())
            if type(literal) is str and len(literal.encode()) > limits.max_string_bytes:
                raise FamilyReplayError("resource", "content", "string byte limit exceeded")
            return 1, size, 0
        if "ref" not in token or type(token["ref"]) is not str or token["ref"] not in table:
            raise FamilyReplayError("format", "content", "unresolved content reference")
        return None

    # Validate all closed tokens before following references. The explicit DFS
    # stack also supports caller depth ceilings above Python's recursion limit.
    token_cost(root)
    for identity, node in table.items():
        if node["kind"] == "array":
            descendants = [(child, 0) for child in node["items"]]
        else:
            descendants, keys = [], set()
            for pair in node["items"]:
                if not isinstance(pair, list) or len(pair) != 2 or type(pair[0]) is not str:
                    raise FamilyReplayError("format", "content", "invalid object entry")
                key = pair[0]
                if key in keys:
                    raise FamilyReplayError("format", "content", "duplicate object key")
                keys.add(key)
                if len(key.encode()) > limits.max_string_bytes:
                    raise FamilyReplayError("resource", "content", "key byte limit exceeded")
                descendants.append((pair[1], len(_json(key).encode()) + 1))
        for child, _ in descendants:
            token_cost(child)
        children[identity] = descendants

    pending = [(root["ref"], False)] if "ref" in root else []
    while pending:
        identity, finishing = pending.pop()
        if not finishing:
            visited.add(identity)
            if identity in costs:
                continue
            if identity in active:
                raise FamilyReplayError("format", "content", "cyclic content reference")
            active.add(identity)
            if len(active) > limits.max_depth + 1:
                raise FamilyReplayError("resource", "content", "expanded depth limit exceeded")
            pending.append((identity, True))
            pending.extend((child["ref"], False) for child, _ in reversed(children[identity]) if "ref" in child)
            continue
        count, size, depth = 1, 2 + max(0, len(children[identity]) - 1), 0
        for child, key_bytes in children[identity]:
            child_count, child_bytes, child_depth = (costs[child["ref"]] if "ref" in child else token_cost(child))
            count += child_count
            size += child_bytes + key_bytes
            depth = max(depth, child_depth + 1)
            if count > limits.max_values or size > limits.max_retained_bytes or depth > limits.max_depth:
                raise FamilyReplayError("resource", "content", "expanded value/byte/depth limit exceeded")
        active.remove(identity)
        costs[identity] = count, size, depth

    expanded_values, _, _ = costs[root["ref"]] if "ref" in root else token_cost(root)
    if visited != table.keys():
        raise FamilyReplayError("format", "content", "surplus unreachable components")
    _charge(budget, "work", expanded_values, "content")
    cache = {}

    def materialized(token):
        return cache[token["ref"]] if "ref" in token else token["literal"]

    # Costs are inserted in postorder, so every child is ready before its parent.
    for identity in costs:
        node = table[identity]
        cache[identity] = ([materialized(child) for child in node["items"]]
                           if node["kind"] == "array" else
                           {key: materialized(child) for key, child in node["items"]})
    return materialized(root)


class _CaptureRequest:
    def __init__(self, annotations, limits=None, *, deferred=False):
        self.limits = _limits(limits)
        self.annotations = _pack(annotations, "annotations", self.limits)
        if len(self.annotations.encode()) > self.limits.max_metadata_bytes:
            raise FamilyReplayError("resource", "annotations", "metadata byte limit exceeded")
        self.deferred = deferred
        self.value_count = _value_count(annotations)

    def annotation_values(self):
        return _unpack(self.annotations, "annotations", self.limits)


class _FamilyCapture:
    def __init__(self, compiled, request, strict, seed):
        self.request, self.limits = request, request.limits
        self.budget = _engine.FamilyReplayBudget(self.limits._json())
        _charge(self.budget, "work", request.value_count, "annotations")
        if "duplex.recombine" not in compiled._pre.pass_names():
            raise FamilyReplayError("unsupported", "capture.workflow", "only source-backed duplex clonal_lineage is supported")
        _reference_policy(compiled._refdata, "capture.reference")
        _stage_shapes(compiled._pre, compiled._post)
        step = compiled._step
        if type(seed) is not int or not 0 <= seed <= 2**64 - 1:
            raise FamilyReplayError("format", "capture.seed", "expected an unsigned 64-bit integer")
        last_clone_seed = seed + (step.n_clones - 1) * 1_000_000
        if last_clone_seed + (0 if step.allow_extinction else 49 * 7919) > 2**64 - 1:
            raise FamilyReplayError("resource", "capture.seed", "clone/retry seed schedule exceeds u64")
        self.seed, self.n_clones = seed, step.n_clones
        bound = step.lambda_base * (1 + step.selection_strength)
        if not math.isfinite(bound) or bound < 0 or bound > 64:
            raise FamilyReplayError("unsupported", "capture.growth", "captured effective offspring rate bound must not exceed 64")
        self.source_spec = _parse(_native("source.plan", compiled._pre.family_replay_spec, budget=self.budget), self.limits, "source.plan")
        self.post_spec = (_parse(_native("observations.plan", compiled._post.family_replay_spec, budget=self.budget), self.limits, "observations.plan")
                          if compiled._post is not None else None)
        self.profile = _parse(_engine._family_replay_numeric_profile(), self.limits, "numeric_profile")
        self.provenance = _native("provenance", _provenance, category="io")
        self.strict = strict
        self.clones = []
        self.family_config = None
        self.sealed = None
        self.expose_provenance = False
        self.output_fence = _OutputFence(self.limits, self.budget, request.annotations, request.value_count)

    def source(self, simulator, founder, clone_id):
        receipt = _parse(_native(f"clones[{clone_id}].source", simulator.family_trace_receipt,
                                 founder, self.strict, budget=self.budget), self.limits, "source.receipt")
        clone = {"clone_id": clone_id, "parent_id": len(self.clones), "source": receipt,
                 "attempts": [], "accepted_attempt": None, "observations": []}
        self.clones.append(clone)
        return clone

    def attempt(self, clone, family):
        path = f"clones[{clone['clone_id']}].attempts[{len(clone['attempts'])}]"
        receipt = _parse(_native(path, family.capture_receipt), self.limits, path)
        config = dict(receipt["config"])
        config.pop("seed_label")
        if self.family_config is None:
            self.family_config = config
        elif _json(config) != _json(self.family_config):
            raise FamilyReplayError("plan_mismatch", path + ".config", "family law changed between attempts")
        clone["attempts"].append(receipt)

    def observation(self, clone, node_id, simulator=None, outcome=None):
        path = f"clones[{clone['clone_id']}].observations[{node_id}]"
        receipt = (None if simulator is None else _parse(
            _native(path, simulator.family_trace_receipt, outcome, self.strict, budget=self.budget), self.limits, path))
        clone["observations"].append({"node_id": node_id, "receipt": receipt})

    def record(self, record):
        """Fence cumulative wrapper output before retaining another row."""
        self.output_fence.record(record)

    def seal(self, result):
        if self.sealed is not None:
            raise FamilyReplayError("sealing", "result", "capture is already sealed")
        for report in (result.validate_records(result._family_replay_refdata, _family_replay_budget=self.budget),
                       result.validate_families(),
                       result.validate_families_with_parents(result._family_replay_refdata, _family_replay_budget=self.budget)):
            if not report.ok:
                raise FamilyReplayError("record_mismatch", "result", str(report.failures[:1]))
        annotations = self.request.annotation_values()
        for i, record in enumerate(result):
            for key, value in annotations.items():
                if key not in record or _pack(record[key], f"records[{i}].{key}", self.limits) != _pack(value, f"annotations.{key}", self.limits):
                    raise FamilyReplayError("record_mismatch", f"records[{i}].{key}", "declared annotation was not applied exactly")
        packed_records = _pack(result.records, "records", self.limits)
        # These native objects expose read-only scientific evidence. Pin their
        # identities separately from the public, mutable container lists.
        self.evidence = tuple(tuple(items) for items in
                              (result.parents, result.outcomes, result.families, result.lineage_trees))
        self.sealed = {
            "workflow": _WORKFLOW,
            "numeric_profile": self.profile,
            "provenance": self.provenance,
            "source_plan": self.source_spec,
            "post_plan": self.post_spec,
            "execution": {"strict": self.strict, "expose_provenance": self.expose_provenance,
                          "allow_extinction": self.allow_extinction, "retry_bound": 50,
                          "retry_stride": 7919, "n_clones": self.n_clones,
                          "seed_label": str(self.seed), "seed_schedule": "clone_1000000_attempt_7919_node_plus_one_v1",
                          "family_config": self.family_config},
            "annotations": json.loads(self.request.annotations),
            "records": json.loads(packed_records),
            "clones": self.clones,
        }
        self.sealed = _json(self.sealed)
        result._family_replay_seal = self.sealed


def _capture_request(annotations, limits=None, *, deferred=False):
    return _CaptureRequest(annotations, limits, deferred=deferred)


def _start_capture(compiled, request, strict, expose_provenance, seed):
    if type(strict) is not bool:
        raise FamilyReplayError("format", "capture.strict", "expected Boolean")
    if type(expose_provenance) is not bool:
        raise FamilyReplayError("format", "capture.expose_provenance", "expected Boolean")
    capture = _FamilyCapture(compiled, request, strict, seed)
    capture.allow_extinction = compiled._step.allow_extinction
    capture.expose_provenance = expose_provenance
    return capture


def _finish_capture(result):
    capture = result._family_replay_capture
    if capture is None:
        raise FamilyReplayError("missing_capture", "result", "run with capture_family_replay=True")
    capture.seal(result)


def _bundle_from_result(result):
    capture = result._family_replay_capture
    if capture is None or result._family_replay_seal is None:
        raise FamilyReplayError("missing_capture", "result", "run with capture_family_replay=True")
    if capture.sealed != result._family_replay_seal:
        raise FamilyReplayError("sealing", "result", "the immutable capture envelope was replaced")
    current = (result.parents, result.outcomes, result.families, result.lineage_trees)
    if any(items is None or len(items) != len(expected)
           or any(item is not original for item, original in zip(items, expected))
           for items, expected in zip(current, capture.evidence)):
        raise FamilyReplayError("record_mismatch", "result.evidence", "the captured native evidence was changed")
    payload = json.loads(result._family_replay_seal)
    if _pack(result.records, "records", capture.limits) != _json(payload["records"]):
        raise FamilyReplayError("record_mismatch", "result.records", "the sealed records or annotations changed")
    budget = _engine.FamilyReplayBudget(capture.limits._json())
    for report in (result.validate_records(result._family_replay_refdata, _family_replay_budget=budget),
                   result.validate_families(),
                   result.validate_families_with_parents(result._family_replay_refdata, _family_replay_budget=budget)):
        if not report.ok:
            raise FamilyReplayError("record_mismatch", "result", str(report.failures[:1]))
    return FamilyReplayBundle._from_payload(payload, capture.limits)


class FamilyReplayBundle:
    """A bounded, self-contained choice replay artifact for duplex lineages."""

    SCHEMA_VERSION = 1

    def __init__(self, text, limits):
        self._limits = _limits(limits)
        self._wire = _parse(text, self._limits, "bundle")
        _require_keys(self._wire, ("artifact", "schema_version", "capabilities", "content", "root", "commitment"), "bundle")
        if (self._wire["artifact"] != _ARTIFACT or type(self._wire["schema_version"]) is not int
                or self._wire["schema_version"] != 1):
            raise FamilyReplayError("unsupported_schema", "bundle", "unsupported artifact or schema version")
        if self._wire["capabilities"] != _CAPABILITIES:
            raise FamilyReplayError("unsupported_capability", "bundle.capabilities", "unknown or missing required capabilities")
        committed = {key: value for key, value in self._wire.items() if key != "commitment"}
        if self._wire["commitment"] != _digest(committed, "genairr_family_bundle_v1"):
            raise FamilyReplayError("integrity", "bundle.commitment", "root commitment mismatch")
        self._text = _json(self._wire)
        self._verified = None

    @classmethod
    def _from_payload(cls, payload, limits):
        root, content = _content_encode(payload, limits)
        wire = {"artifact": _ARTIFACT, "schema_version": 1, "capabilities": _CAPABILITIES,
                "content": content, "root": root}
        wire["commitment"] = _digest(wire, "genairr_family_bundle_v1")
        text = _json(wire)
        if len(text.encode()) > limits.max_file_bytes:
            raise FamilyReplayError("resource", "bundle", "file byte limit exceeded")
        return cls(text, limits)

    @classmethod
    def from_json(cls, text, *, limits=None):
        return cls(text, _limits(limits))

    @classmethod
    def read_from(cls, path, *, limits=None):
        limits = _limits(limits)
        try:
            with open(path, "rb") as handle:
                data = handle.read(limits.max_file_bytes + 1)
            if len(data) > limits.max_file_bytes:
                raise FamilyReplayError("resource", "bundle", "file byte limit exceeded")
            return cls(data.decode("utf-8"), limits)
        except FamilyReplayError:
            raise
        except (OSError, UnicodeError) as exc:
            raise FamilyReplayError("io", str(path), str(exc)) from exc

    def to_json(self):
        return self._text

    def write_to(self, path, *, replace=False):
        """Publish atomically, refusing to replace an existing file by default."""
        if type(replace) is not bool:
            raise TypeError("replace must be Boolean")
        try:
            write_atomic(path, self._text.encode("utf-8"), replace=replace)
        except ArtifactIOError as exc:
            raise FamilyReplayError("io", str(path), str(exc)) from exc

    def replay_records(self, *, limits=None):
        """Validate every source, attempt and observed read before returning truth."""
        selected_limits = _limits(limits) if limits is not None else self._limits
        self._verified = None
        # Reapply caller limits to the frozen bytes. Parsed private dictionaries
        # and earlier successful returned results are never replay authority.
        wire = _parse(self._text, selected_limits, "bundle")
        budget = _engine.FamilyReplayBudget(selected_limits._json())
        payload = _content_decode(wire["root"], wire["content"], selected_limits, budget)
        _require_keys(payload, ("workflow", "numeric_profile", "provenance", "source_plan", "post_plan", "execution",
                                "annotations", "records", "clones"), "execution")
        _require_keys(payload["provenance"], ("package_version", "engine_version", "extension_sha256"), "provenance")
        for field in ("package_version", "engine_version"):
            if type(payload["provenance"][field]) is not str:
                raise FamilyReplayError("format", f"provenance.{field}", "expected string")
        extension_hash = payload["provenance"]["extension_sha256"]
        if extension_hash is not None and (type(extension_hash) is not str or len(extension_hash) != 64
                                           or any(c not in "0123456789abcdef" for c in extension_hash)):
            raise FamilyReplayError("format", "provenance.extension_sha256", "expected optional SHA256")
        if payload["workflow"] != _WORKFLOW:
            raise FamilyReplayError("unsupported", "execution.workflow", "only duplex_lineage_v1 is supported")
        profile = _parse(_engine._family_replay_numeric_profile(), selected_limits, "numeric_profile")
        if payload["numeric_profile"] != profile:
            raise FamilyReplayError("numerical_incompatibility", "numeric_profile", "this numeric execution profile is not qualified")
        source = _native("source.plan", _engine.CompiledSimulator.from_family_replay_spec, _json(payload["source_plan"]), budget=budget)
        post = (_native("observations.plan", _engine.CompiledSimulator.from_family_replay_spec, _json(payload["post_plan"]), budget=budget)
                if payload["post_plan"] is not None else None)
        refdata = _native("reference", source.family_replay_refdata)
        _reference_policy(refdata, "reference")
        _stage_shapes(source, post)
        if payload["post_plan"] is not None:
            for field in ("reference", "contracts", "productivity_policy"):
                if payload["post_plan"][field] != payload["source_plan"][field]:
                    raise FamilyReplayError("plan_mismatch", f"observations.plan.{field}", "source and read contexts differ")
        execution = payload["execution"]
        _require_keys(execution, ("strict", "expose_provenance", "allow_extinction", "retry_bound", "retry_stride",
                                  "n_clones", "seed_label", "seed_schedule", "family_config"), "execution.settings")
        if (any(type(execution[key]) is not bool for key in ("strict", "expose_provenance", "allow_extinction"))
                or execution["retry_bound"] != 50 or execution["retry_stride"] != 7919
                or execution["seed_schedule"] != "clone_1000000_attempt_7919_node_plus_one_v1"):
            raise FamilyReplayError("format", "execution.settings", "unsupported execution policy")
        run_seed = _seed(execution["seed_label"], "execution.seed_label")
        n_clones = _index(execution["n_clones"], "execution.n_clones")
        if n_clones == 0:
            raise FamilyReplayError("format", "execution.n_clones", "expected positive clone count")
        config = execution["family_config"]
        if not isinstance(config, dict) or "seed_label" in config or config.get("strict") is not execution["strict"]:
            raise FamilyReplayError("plan_mismatch", "execution.family_config", "invalid shared family law")
        expected_contracts = (None if payload["source_plan"]["contracts"] is None else source.active_contracts())
        if config.get("contracts") != expected_contracts:
            raise FamilyReplayError("plan_mismatch", "execution.family_config.contracts", "division/source contracts differ")
        if len(_json(payload["annotations"]).encode()) > selected_limits.max_metadata_bytes:
            raise FamilyReplayError("resource", "annotations", "metadata byte limit exceeded")
        annotations = _unpack(_json(payload["annotations"]), "annotations", selected_limits)
        expected_records = _unpack(_json(payload["records"]), "records", selected_limits)
        if not isinstance(annotations, dict) or not isinstance(expected_records, list):
            raise FamilyReplayError("format", "records", "invalid records or annotation payload")
        output_fence = _OutputFence(selected_limits, budget, _json(payload["annotations"]), _value_count(annotations))
        from .result import SimulationResultWithLineages
        records, outcomes, parents, families, trees, origins = [], [], [], [], [], []
        if not isinstance(payload["clones"], list) or len(payload["clones"]) != n_clones:
            raise FamilyReplayError("format", "clones", "clone count does not match execution specification")
        for clone_index, clone in enumerate(payload["clones"]):
            path = f"clones[{clone_index}]"
            _require_keys(clone, ("clone_id", "parent_id", "source", "attempts", "accepted_attempt", "observations"), path)
            if _index(clone["clone_id"], path) != clone_index or _index(clone["parent_id"], path) != len(parents):
                raise FamilyReplayError("invalid_scope", path, "clone/parent addressing mismatch")
            clone_seed = run_seed + clone_index * 1_000_000
            _stage_policy(clone["source"], execution["strict"], clone_seed, path + ".source")
            founder = _native(path + ".source", source.replay_family_trace, None, _json(clone["source"]), budget=budget)
            attempts = clone["attempts"]
            accepted = _index(clone["accepted_attempt"], path + ".accepted_attempt")
            if (not isinstance(attempts, list) or not attempts or len(attempts) > execution["retry_bound"]
                    or accepted != len(attempts) - 1):
                raise FamilyReplayError("invalid_order", path + ".attempts", "accepted attempt must follow all recorded failed attempts")
            family = None
            for attempt_index, receipt in enumerate(attempts):
                stage_path = f"{path}.attempts[{attempt_index}]"
                if not isinstance(receipt, dict) or not isinstance(receipt.get("config"), dict):
                    raise FamilyReplayError("format", stage_path, "missing native family configuration")
                attempt_config = dict(receipt["config"])
                attempt_seed = _seed(attempt_config.pop("seed_label", None), stage_path + ".config.seed_label")
                if _json(attempt_config) != _json(config) or attempt_seed != clone_seed + attempt_index * 7919:
                    raise FamilyReplayError("plan_mismatch", stage_path + ".config", "attempt law or scheduling differs")
                candidate = _native(stage_path, _engine.replay_family_outcomes, founder, refdata, _json(receipt), budget=budget)
                survived = _native(stage_path, _engine._family_replay_has_observed, candidate, budget=budget)
                eligible = survived or execution["allow_extinction"]
                if eligible != (attempt_index == accepted):
                    raise FamilyReplayError("invalid_choice", stage_path, "survival/first-accepted attempt policy mismatch")
                if eligible:
                    family = candidate
            parents.append(founder)
            families.append(family)
            tree = family.tree()
            trees.append(tree)
            observed = [node for node in tree.nodes() if node.observed]
            if len(records) + len(observed) > len(expected_records):
                raise FamilyReplayError("record_mismatch", "records", "observed cells exceed the sealed row count")
            observations = clone["observations"]
            if not isinstance(observations, list) or len(observations) != len(observed):
                raise FamilyReplayError("invalid_order", path + ".observations", "one observation per observed cell is required")
            native_nodes = _native(path + ".nodes", _engine._family_replay_node_outcomes, family, budget=budget)
            for node, observation in zip(observed, observations):
                stage_path = f"{path}.observations[{node.id}]"
                _require_keys(observation, ("node_id", "receipt"), stage_path)
                if _index(observation["node_id"], stage_path) != node.id:
                    raise FamilyReplayError("invalid_scope", stage_path, "observed cell order mismatch")
                outcome = native_nodes[node.id]
                if post is None:
                    if observation["receipt"] is not None:
                        raise FamilyReplayError("invalid_choice", stage_path, "pristine observations cannot execute passes")
                else:
                    if observation["receipt"] is None:
                        raise FamilyReplayError("invalid_choice", stage_path, "missing post-read stage")
                    _stage_policy(observation["receipt"], execution["strict"], clone_seed + 1 + node.id, stage_path)
                    initial = _native(stage_path, _engine._family_replay_final_simulation, outcome, budget=budget)
                    read = _native(stage_path, post.replay_family_trace, initial,
                                   _json(observation["receipt"]), budget=budget)
                    outcome = _native(stage_path, _engine._family_replay_merge_outcomes, outcome, read, budget=budget)
                record = _native(stage_path, _engine._family_replay_airr_record, outcome, refdata,
                                 f"clone{clone_index}_node{node.id}", budget=budget)
                record.update({"clone_id": clone_index, "parent_id": len(parents) - 1,
                               "lineage_node_id": node.id, "lineage_parent_id": node.parent_id if node.parent_id is not None else -1,
                               "lineage_generation": node.generation, "lineage_abundance": node.abundance,
                               "duplicate_count": node.abundance, "lineage_affinity": node.affinity,
                               "lineage_mutations_from_parent": node.mutations_from_parent})
                if execution["expose_provenance"]:
                    _truth_aliases(record)
                output_fence.record(record)
                for key, value in annotations.items():
                    if key in record:
                        raise FamilyReplayError("metadata", f"annotations.{key}", "annotation overwrites a generated field")
                    record[key] = value
                records.append(record)
                outcomes.append(outcome)
                origins.append((len(parents) - 1, clone_index, node.id))
        if _pack(records, "records", selected_limits) != _pack(expected_records, "records.expected", selected_limits):
            raise FamilyReplayError("record_mismatch", "records", "replayed records differ in fields, values, types or order")
        result = SimulationResultWithLineages(records, outcomes=outcomes, parents=parents, lineage_trees=trees,
                                             families=families, lineage_record_origins=origins)
        for report in (result.validate_records(refdata, _family_replay_budget=budget), result.validate_families(),
                       result.validate_families_with_parents(refdata, _family_replay_budget=budget)):
            if not report.ok:
                raise FamilyReplayError("record_mismatch", "records", str(report.failures[:1]))
        self._verified = (tuple(families), refdata, budget)
        result._family_replay_refdata = refdata
        return result

    @property
    def refdata(self):
        """The reconstructed reference after full semantic replay has passed."""
        if self._verified is None:
            raise FamilyReplayError("not_verified", "reference", "call replay_records before requesting native reference data")
        return self._verified[1]

    def replay_node(self, *, clone_id, node_id, limits=None):
        clone_id = _index(clone_id, "clone_id")
        node_id = _index(node_id, "node_id")
        # Validate the whole artifact, including read stages, before publishing a
        # cellular node. Never reuse mutable result dictionaries as evidence.
        self.replay_records(limits=limits)
        families, refdata, budget = self._verified
        if clone_id >= len(families):
            raise FamilyReplayError("invalid_scope", "clone_id", "clone is outside the bundle")
        return _native(f"clones[{clone_id}].nodes[{node_id}]", families[clone_id].replay_node,
                       node_id, refdata, budget=budget)
