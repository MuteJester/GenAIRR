"""Finite retained coding-end models and supervised, in-house estimation.

Training consumes recorded generating choices. It never infers a hairpin
opening or P ancestry from a palindromic match in an observed junction.
"""
from __future__ import annotations

import json
import math
from collections import Counter
from dataclasses import dataclass
from typing import Iterable, Mapping, Optional, Tuple

_VERSION = "coding_end_outcome_model_v1"
_PROCESSING = "original_coding_end_v1"
_LAW = "unconditioned_generating_choices"


def _integer(value, name):
    if type(value) is not int or not -(2**63) <= value < 2**63:
        raise ValueError(f"{name} must be a signed 64-bit integer")
    return value


def _rows(values, segment):
    combined = {}
    for row in values:
        if not isinstance(row, (tuple, list)) or len(row) != 3:
            raise ValueError(f"{segment} outcomes require (five, three, weight) rows")
        five, three = (_integer(row[i], f"{segment} outcome") for i in (0, 1))
        if (segment == "V" and five != 0) or (segment == "J" and three != 0):
            raise ValueError(f"{segment} supports processing only at its junction-facing end")
        if isinstance(row[2], bool) or not isinstance(row[2], (int, float)):
            raise ValueError(f"{segment} weights must be finite positive numbers")
        weight = float(row[2])
        if not math.isfinite(weight) or weight <= 0:
            raise ValueError(f"{segment} weights must be finite positive numbers")
        combined.setdefault((five, three), []).append(weight)
    if not combined:
        raise ValueError(f"{segment} outcome support must be nonempty")
    try:
        sums = {key: math.fsum(weights) for key, weights in combined.items()}
        total = math.fsum(sums.values())
    except OverflowError:
        raise ValueError(f"{segment} weight sum must be finite") from None
    if not math.isfinite(total):
        raise ValueError(f"{segment} weight sum must be finite")
    return tuple((a, b, w) for (a, b), w in sorted(sums.items()))


def _choices(record, segments):
    if record.get("duplex_repair") not in (None, ""):
        raise ValueError("duplex repair histories require a history-aware fitting model")
    metadata = record.get("coding_end_processing")
    if isinstance(metadata, str):
        try:
            metadata = json.loads(metadata)
        except (TypeError, ValueError):
            raise ValueError("coding_end_processing must be valid JSON") from None
    if not isinstance(metadata, Mapping) or metadata.get("model") != _PROCESSING:
        raise ValueError("training requires original_coding_end_v1 generating provenance")
    if metadata.get("recombination_complete") is not True:
        raise ValueError("training requires completed recombination provenance")
    if metadata.get("boundary_contract") != "provided_allele_coding_ends":
        raise ValueError("training requires the original coding-boundary contract")
    found = {}
    choices = metadata.get("choices")
    if not isinstance(choices, (tuple, list)):
        raise ValueError("generating choices must be an array")
    for choice in choices:
        if not isinstance(choice, Mapping):
            raise ValueError("malformed generating choice")
        segment = choice.get("segment")
        if segment not in segments or segment in found:
            raise ValueError("each training row requires exactly one choice per segment")
        allele = choice.get("original_allele")
        if not isinstance(allele, str) or not allele:
            raise ValueError("generating choices require an original allele identity")
        five = _integer(choice.get("five_prime_outcome"), "five_prime_outcome")
        three = _integer(choice.get("three_prime_outcome"), "three_prime_outcome")
        _rows(((five, three, 1.0),), segment)
        found[segment] = (five, three)
    if set(found) != set(segments):
        raise ValueError("each training row requires exactly one choice per segment")
    return found


@dataclass(frozen=True)
class CodingEndOutcomeModel:
    """A portable finite model of canonical signed retained outcomes.

    ``v`` and ``j`` contain joint ``(five, three, weight)`` rows with their
    non-junction-facing coordinate fixed at zero. Optional ``d`` rows preserve
    correlation between the two D ends. Weights need not be normalized.
    ``reference_hashes`` restrict use to declared compatible native references.
    """

    v: Tuple[Tuple[int, int, float], ...]
    j: Tuple[Tuple[int, int, float], ...]
    d: Optional[Tuple[Tuple[int, int, float], ...]] = None
    reference_hashes: Tuple[str, ...] = ()
    provenance_json: str = "{}"

    def __post_init__(self):
        if self.v is None or self.j is None:
            raise ValueError("V and J outcome support must be nonempty")
        for name in ("v", "j", "d"):
            rows = getattr(self, name)
            if rows is not None:
                object.__setattr__(self, name, _rows(rows, name.upper()))
        hashes = tuple(self.reference_hashes)
        if any(not isinstance(h, str) or not h.startswith("sha256:") or
               len(h) != 71 or any(c not in "0123456789abcdef" for c in h[7:]) for h in hashes):
            raise ValueError("reference_hashes must contain native sha256 content identities")
        object.__setattr__(self, "reference_hashes", tuple(sorted(set(hashes))))
        try:
            provenance = json.loads(self.provenance_json)
            if not isinstance(provenance, dict):
                raise ValueError()
            canonical = json.dumps(provenance, sort_keys=True, allow_nan=False)
        except (TypeError, ValueError):
            raise ValueError("provenance_json must encode a finite JSON object") from None
        object.__setattr__(self, "provenance_json", canonical)

    @property
    def provenance(self):
        """Return an independent copy of the authoring or calibration report."""
        return json.loads(self.provenance_json)

    def to_dict(self):
        return {"schema": _VERSION, "processing_model": _PROCESSING,
                "v": [list(r) for r in self.v], "j": [list(r) for r in self.j],
                "d": None if self.d is None else [list(r) for r in self.d],
                "reference_hashes": list(self.reference_hashes),
                "provenance": self.provenance}

    @classmethod
    def from_dict(cls, value):
        if not isinstance(value, Mapping) or value.get("schema") != _VERSION or value.get("processing_model") != _PROCESSING:
            raise ValueError("unsupported coding-end model schema")
        allowed = {"schema", "processing_model", "v", "d", "j", "reference_hashes", "provenance"}
        if set(value) != allowed:
            raise ValueError("coding-end model fields are missing or unrecognized")
        return cls(v=value["v"], j=value["j"], d=value["d"],
                   reference_hashes=value["reference_hashes"],
                   provenance_json=json.dumps(value["provenance"], allow_nan=False))

    @classmethod
    def fit_generating_records(cls, records: Iterable[Mapping], *, chain_type: str,
                               reference_hash: str, source: str,
                               training_sampling_law: str, independent_generating_events: bool = False,
                               pseudocount: float = 0.0,
                               support: Optional["CodingEndOutcomeModel"] = None):
        """Estimate pooled retained outcomes from genuine generating histories.

        The required sampling-law declaration is an input contract. Productive
        or other outcome-selected samples require an ascertainment likelihood
        and are rejected here. The result estimates observed generating-choice
        frequencies, not microscopic enzyme rates or unique biological histories.
        With varying allele feasibility, these pooled frequencies need not equal
        the pre-filter authored weights. ``support`` supplies explicit prior
        support; pseudocounts never invent unobserved outcomes implicitly.
        """
        if chain_type not in ("vj", "vdj"):
            raise ValueError("chain_type must be 'vj' or 'vdj'")
        if training_sampling_law != _LAW:
            raise ValueError(f"supervised fitting requires {_LAW!r}; selected samples need an ascertainment likelihood")
        if independent_generating_events is not True:
            raise ValueError("declare independent_generating_events=True after selecting independent recombination events")
        if not isinstance(source, str) or not source.strip():
            raise ValueError("source must describe the training data")
        if isinstance(pseudocount, bool) or not isinstance(pseudocount, (int, float)) or not math.isfinite(pseudocount) or pseudocount < 0:
            raise ValueError("pseudocount must be a finite nonnegative number")
        segments = ("V", "D", "J") if chain_type == "vdj" else ("V", "J")
        if support is not None:
            if not isinstance(support, cls) or (support.d is not None) != (chain_type == "vdj"):
                raise ValueError("support must be a model with the same chain type")
            if support.reference_hashes and reference_hash not in support.reference_hashes:
                raise ValueError("training reference is incompatible with model support")
        counters = {s: Counter() for s in segments}
        n = 0
        for row in records:
            if not isinstance(row, Mapping):
                raise ValueError("training records must be mappings")
            if any(row.get(key) is not None for key in ("clone_id", "family_id", "lineage_node_id", "node_id")):
                raise ValueError("clonal/lineage rows repeat generating events; supply independently selected founders")
            if row.get("duplicate_count", 1) not in (None, 1):
                raise ValueError("collapsed observations require independent-event accounting before estimation")
            choices = _choices(row, segments)
            for segment, pair in choices.items():
                counters[segment][pair] += 1
            n += 1
        if not n:
            raise ValueError("training data must contain at least one complete generating history")
        outputs, report = {}, {}
        for segment, counts in counters.items():
            domain = set(counts) if support is None else {(a, b) for a, b, _ in getattr(support, segment.lower())}
            if not set(counts) <= domain:
                raise ValueError(f"{segment} training outcome lies outside declared support")
            total = n + pseudocount * len(domain)
            if not math.isfinite(total):
                raise ValueError("prior support mass must remain finite")
            rows, intervals = [], []
            for five, three in sorted(domain):
                count = counts[(five, three)]
                weight = (count + pseudocount) / total
                if weight == 0 and (count > 0 or pseudocount > 0):
                    raise ValueError("outcome probability underflow; rescale the prior")
                if weight > 0:
                    rows.append((five, three, weight))
                # Marginal Wilson interval on raw observed proportions; this is
                # not a simultaneous interval and excludes prior pseudocounts.
                z2, p = 1.959963984540054**2, count / n
                center = (p + z2 / (2*n)) / (1 + z2/n)
                half = math.sqrt(z2 * (p*(1-p)/n + z2/(4*n*n))) / (1 + z2/n)
                intervals.append({"five": five, "three": three, "count": count,
                                  "marginal_95_interval": [max(0.0, center-half), min(1.0, center+half)]})
            outputs[segment.lower()] = tuple(rows)
            report[segment] = intervals
        provenance = {"estimator": "supervised_retained_outcomes_v1", "source": source,
                      "training_sampling_law": training_sampling_law, "records": n,
                      "conditioning": "pooled_by_segment", "pseudocount": pseudocount,
                      "independent_generating_events": True,
                      "uncertainty": "marginal Wilson intervals on observed counts",
                      "counts": report,
                      "limitations": ["not biochemical rate calibration",
                                      "pooled frequencies reflect allele feasibility",
                                      "sampling-law and source declarations supplied by caller"]}
        return cls(**outputs, reference_hashes=(reference_hash,),
                   provenance_json=json.dumps(provenance, allow_nan=False))

    def score_generating_records(self, records: Iterable[Mapping]):
        """Score held-out generating outcomes under the segment-factorized model.

        This evaluates the authored outcome table; it is not a read likelihood
        and includes neither allele-feasibility renormalization nor selection.
        """
        segments = ("V", "D", "J") if self.d is not None else ("V", "J")
        log_probabilities = {}
        for segment in segments:
            rows = getattr(self, segment.lower())
            total = math.fsum(w for _, _, w in rows)
            log_probabilities[segment] = {(a, b): math.log(w)-math.log(total) for a, b, w in rows}
        score, n, unseen = 0.0, 0, 0
        for row in records:
            choices = _choices(row, segments)
            for segment, pair in choices.items():
                log_p = log_probabilities[segment].get(pair)
                if log_p is None:
                    unseen += 1
                else:
                    score += log_p
            n += 1
        if not n:
            raise ValueError("scoring requires at least one generating history")
        return {"records": n, "log_likelihood": -math.inf if unseen else score,
                "unseen_outcomes": unseen, "scope": "generating_outcome_table"}

    def fit_compatible_histories(self, observations, *, source: str,
                                 training_sampling_law: str,
                                 independent_observations: bool = False,
                                 complete_scenarios: bool = False,
                                 max_iterations: int = 200, tolerance: float = 1e-9):
        """Fit a finite mixture by marginalizing explicitly compatible histories.

        Each observation is a list of ``{"outcomes": {"V": (five, three),
        "D": (five, three), "J": (five, three)}, "likelihood": positive}``
        scenarios. The supplied likelihood is P(observation | history), with
        other modeled processes already marginalized. Scenario enumeration,
        allele-feasibility normalization and ascertainment must be accounted for
        by the caller's forward model; generic AIRR annotations do not provide
        them. This method estimates finite retained-outcome probabilities, not
        microscopic histories or enzyme rates. Identifiability is not assumed.
        """
        if training_sampling_law != _LAW:
            raise ValueError("selected samples require an explicit ascertainment likelihood")
        if independent_observations is not True or complete_scenarios is not True:
            raise ValueError("declare independent_observations=True and complete_scenarios=True for the supplied forward model")
        if not isinstance(source, str) or not source.strip():
            raise ValueError("source must describe observations and the forward likelihood")
        if type(max_iterations) is not int or not 1 <= max_iterations <= 10_000:
            raise ValueError("max_iterations must be an integer in [1, 10000]")
        if isinstance(tolerance, bool) or not isinstance(tolerance, (int, float)) or not math.isfinite(tolerance) or tolerance <= 0:
            raise ValueError("tolerance must be finite and positive")
        segments = ("V", "D", "J") if self.d is not None else ("V", "J")
        domains = {s: [(a, b) for a, b, _ in getattr(self, s.lower())] for s in segments}
        indices = {s: {pair: i for i, pair in enumerate(domain)} for s, domain in domains.items()}
        probabilities = {}
        for s in segments:
            rows = getattr(self, s.lower())
            total = math.fsum(w for _, _, w in rows)
            probabilities[s] = [w/total for _, _, w in rows]
            if any(p == 0 for p in probabilities[s]):
                raise ValueError("initial probability underflow; rescale model weights")
        data, scenario_count = [], 0
        for observation in observations:
            scenarios = []
            for scenario in observation:
                if not isinstance(scenario, Mapping) or set(scenario) != {"outcomes", "likelihood"}:
                    raise ValueError("each compatible history requires outcomes and likelihood")
                outcomes = scenario["outcomes"]
                if not isinstance(outcomes, Mapping) or set(outcomes) != set(segments):
                    raise ValueError("compatible history has a different chain type")
                path = []
                for s in segments:
                    pair = outcomes[s]
                    if not isinstance(pair, (list, tuple)) or len(pair) != 2:
                        raise ValueError("history outcomes must be (five, three) pairs")
                    pair = tuple(_integer(v, "history outcome") for v in pair)
                    if pair not in indices[s]:
                        raise ValueError("compatible history lies outside declared finite support")
                    path.append(indices[s][pair])
                weight = scenario["likelihood"]
                if isinstance(weight, bool) or not isinstance(weight, (int, float)) or not math.isfinite(weight) or weight <= 0:
                    raise ValueError("history likelihoods must be finite and positive")
                scenarios.append((tuple(path), math.log(weight)))
                scenario_count += 1
                if scenario_count > 1_000_000:
                    raise ValueError("compatible-history scenario budget exceeded")
            if not scenarios:
                raise ValueError("each observation needs at least one compatible history")
            data.append(scenarios)
        if not data:
            raise ValueError("fitting requires at least one observation")

        def expectation(params):
            counts = {s: [0.0]*len(domains[s]) for s in segments}
            log_params = {s: [math.log(p) if p > 0 else -math.inf for p in params[s]] for s in segments}
            log_likelihood, entropy = 0.0, 0.0
            for scenarios in data:
                logs = [w + math.fsum(log_params[s][path[i]] for i, s in enumerate(segments)) for path, w in scenarios]
                peak = max(logs)
                if not math.isfinite(peak):
                    raise ValueError("observation has zero probability under the fitted support")
                normalizer = peak + math.log(math.fsum(math.exp(w-peak) for w in logs))
                log_likelihood += normalizer
                for (path, _), log_p in zip(scenarios, logs):
                    posterior = math.exp(log_p-normalizer)
                    if posterior:
                        entropy -= posterior*(log_p-normalizer)
                    for i, s in enumerate(segments):
                        counts[s][path[i]] += posterior
            return counts, log_likelihood, entropy/len(data)

        trajectory, converged = [], False
        for _ in range(max_iterations):
            counts, ll, entropy = expectation(probabilities)
            trajectory.append(ll)
            if len(trajectory) > 1 and abs(trajectory[-1]-trajectory[-2]) <= tolerance*(1+abs(trajectory[-2])):
                converged = True
                break
            probabilities = {s: [c/len(data) for c in counts[s]] for s in segments}
        _, final_ll, entropy = expectation(probabilities)
        if final_ll != trajectory[-1]:
            trajectory.append(final_ll)
        outputs = {s.lower(): tuple((a,b,p) for (a,b),p in zip(domains[s], probabilities[s]) if p > 0)
                   for s in segments}
        report = {"estimator": "finite_compatible_history_em_v1", "source": source,
                  "training_sampling_law": training_sampling_law, "observations": len(data),
                  "independent_observations": True, "complete_scenarios": True,
                  "scenarios": scenario_count, "converged": converged,
                  "log_likelihood_trajectory": trajectory,
                  "mean_history_posterior_entropy": entropy,
                  "identifiability": "not established; equivalent histories and parameter confounding may remain",
                  "uncertainty_intervals": None,
                  "limitations": ["caller supplies complete compatible scenarios and forward likelihoods",
                                  "no automatic inference from generic AIRR records",
                                  "no biochemical rate calibration"]}
        return type(self)(**outputs, reference_hashes=self.reference_hashes,
                          provenance_json=json.dumps(report, allow_nan=False))
