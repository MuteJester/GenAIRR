#!/usr/bin/env python3
"""
export_dataconfig.py — Convert a Python DataConfig + S5F model data
to a C source file with static arrays.

Usage:
    python3 scripts/export_dataconfig.py HUMAN_IGH_IMGT

Generates:
    data/<config_name>.c   — Static data arrays
    data/<config_name>.h   — Header with loader declarations
"""

import sys
import os
import pickle
import math
from importlib import resources
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

from GenAIRR.data import _CONFIG_NAMES
from GenAIRR.dataconfig.enums import ChainType


def load_dataconfig(name):
    """Load a DataConfig by name."""
    pkg = "GenAIRR.data.builtin_dataconfigs"
    filename = f"{name}.pkl"
    with resources.path(pkg, filename) as path:
        with open(path, "rb") as f:
            return pickle.load(f)


def load_s5f_model(chain_category):
    """Load S5F mutability and substitution data."""
    pkg = "GenAIRR.data.mutation_model_parameters"
    filename = "HH_S5F_META.pkl" if chain_category == "heavy" else "HKL_S5F_META.pkl"
    with resources.path(pkg, filename) as path:
        with open(path, "rb") as f:
            mutability, substitution, targeting = pickle.load(f)

    # Clean NaN values from mutability
    mutability = {
        k: (v if isinstance(v, (int, float)) and not math.isnan(v) else 0.0)
        for k, v in mutability.items()
    }

    # Convert substitution DataFrame if needed
    if hasattr(substitution, "to_dict"):
        substitution = substitution.to_dict(orient="dict")
        substitution = {
            outer_key: {
                inner_key: inner_value
                for inner_key, inner_value in outer_dict.items()
                if isinstance(inner_value, (int, float)) and not math.isnan(inner_value)
            }
            for outer_key, outer_dict in substitution.items()
        }

    return mutability, substitution


def chain_type_to_c(ct):
    """Map Python ChainType to C enum name."""
    mapping = {
        ChainType.BCR_HEAVY: "CHAIN_IGH",
        ChainType.BCR_LIGHT_KAPPA: "CHAIN_IGK",
        ChainType.BCR_LIGHT_LAMBDA: "CHAIN_IGL",
        ChainType.TCR_ALPHA: "CHAIN_TCRA",
        ChainType.TCR_BETA: "CHAIN_TCRB",
        ChainType.TCR_DELTA: "CHAIN_TCRD",
        ChainType.TCR_GAMMA: "CHAIN_TCRG",
    }
    return mapping.get(ct, "CHAIN_IGH")


def segment_to_c(seg_str):
    """Map segment string to C enum."""
    return {"V": "SEG_V", "D": "SEG_D", "J": "SEG_J", "C": "SEG_C"}[seg_str]


BASE_MAP = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4,
            'a': 0, 'c': 1, 'g': 2, 't': 3}


def encode_kmer(kmer):
    """Encode a 5-mer string to integer key (base-5)."""
    return (BASE_MAP.get(kmer[0], 4) * 625 +
            BASE_MAP.get(kmer[1], 4) * 125 +
            BASE_MAP.get(kmer[2], 4) * 25 +
            BASE_MAP.get(kmer[3], 4) * 5 +
            BASE_MAP.get(kmer[4], 4))


def flatten_alleles(allele_dict, seg_type):
    """Flatten gene→allele_list dict into a flat list of allele dicts."""
    result = []
    for gene_name, allele_list in allele_dict.items():
        for allele in allele_list:
            result.append({
                "name": allele.name,
                "gene": allele.gene,
                "family": allele.family,
                "seq": allele.ungapped_seq.upper(),
                "length": allele.ungapped_len,
                "anchor": allele.anchor if allele.anchor is not None else 0,
                "segment": seg_type,
            })
    return result


def aggregate_trim_dist(trim_dict):
    """Aggregate per-gene trim distributions into a single distribution.

    Takes the average across all genes, normalized.
    """
    all_amounts = {}
    count = 0
    for fam in trim_dict.values():
        for gene in fam.values():
            for amount, prob in gene.items():
                amount = int(amount)
                all_amounts[amount] = all_amounts.get(amount, 0.0) + prob
            count += 1

    if count == 0:
        return []

    max_amount = max(all_amounts.keys()) if all_amounts else 0
    probs = []
    for i in range(max_amount + 1):
        probs.append(all_amounts.get(i, 0.0) / count)

    # Normalize
    total = sum(probs)
    if total > 0:
        probs = [p / total for p in probs]

    return probs


def write_c_files(config_name, dc, s5f_mutability, s5f_substitution, out_dir):
    """Generate .c and .h files with static data."""
    c_chain = chain_type_to_c(dc.metadata.chain_type)
    has_d = dc.metadata.has_d
    prefix = config_name.lower()

    # Flatten alleles
    v_alleles = flatten_alleles(dc.v_alleles, "V")
    d_alleles = flatten_alleles(dc.d_alleles, "D") if dc.d_alleles else []
    j_alleles = flatten_alleles(dc.j_alleles, "J")

    # Aggregate trim distributions
    v_trim_3 = aggregate_trim_dist(dc.trim_dicts.get("V_3", {}))
    j_trim_5 = aggregate_trim_dist(dc.trim_dicts.get("J_5", {}))
    d_trim_5 = aggregate_trim_dist(dc.trim_dicts.get("D_5", {})) if has_d else []
    d_trim_3 = aggregate_trim_dist(dc.trim_dicts.get("D_3", {})) if has_d else []

    # NP lengths
    np1_lengths = dc.NP_lengths.get("NP1", {})
    np2_lengths = dc.NP_lengths.get("NP2", {})
    np1_max = max(np1_lengths.keys()) if np1_lengths else 0
    np2_max = max(np2_lengths.keys()) if np2_lengths else 0

    # Gene usage
    v_usage = dc.gene_use_dict.get("V", {})
    j_usage = dc.gene_use_dict.get("J", {})
    d_usage = dc.gene_use_dict.get("D", {}) if has_d else {}

    # S5F: build flat arrays
    # mutability: kmer_string -> float
    mutability_flat = [0.0] * 3125
    for kmer, val in s5f_mutability.items():
        if len(kmer) == 5:
            key = encode_kmer(kmer)
            if 0 <= key < 3125:
                mutability_flat[key] = float(val)

    # substitution: kmer_string -> {base: weight}
    # We store as: for each kmer key, up to 4 (base, weight) pairs
    sub_data = []  # list of (key, bases_str, weights_list)
    for kmer, subs in s5f_substitution.items():
        if len(kmer) != 5 or not subs:
            continue
        key = encode_kmer(kmer)
        if key < 0 or key >= 3125:
            continue
        bases = ""
        weights = []
        for base, weight in sorted(subs.items()):
            if base in "ACGT" and isinstance(weight, (int, float)):
                bases += base
                weights.append(float(weight))
        if bases:
            sub_data.append((key, bases, weights))

    # ── Write .h file ────────────────────────────────────────────
    h_path = os.path.join(out_dir, f"{prefix}.h")
    with open(h_path, "w") as h:
        guard = f"GENAIRR_DATA_{config_name}_H"
        h.write(f"/* Auto-generated from {config_name} DataConfig. Do not edit. */\n\n")
        h.write(f"#ifndef {guard}\n#define {guard}\n\n")
        h.write('#include "genairr/genairr.h"\n\n')
        h.write(f"/* Load the {config_name} data into a SimConfig. */\n")
        h.write(f"void {prefix}_load_config(SimConfig *cfg);\n\n")
        h.write(f"/* Load the {config_name} S5F model data. */\n")
        h.write(f"void {prefix}_load_s5f(S5FModel *model);\n\n")
        h.write(f"#endif /* {guard} */\n")

    # ── Write .c file ────────────────────────────────────────────
    c_path = os.path.join(out_dir, f"{prefix}.c")
    with open(c_path, "w") as c:
        c.write(f"/* Auto-generated from {config_name} DataConfig. Do not edit. */\n\n")
        c.write(f'#include "{prefix}.h"\n')
        c.write('#include <string.h>\n')
        c.write('#include <stdlib.h>\n\n')

        # ── Allele data ──────────────────────────────────────────
        def write_alleles(f, alleles, var_name):
            f.write(f"static const int {var_name}_count = {len(alleles)};\n")
            if not alleles:
                f.write(f"static const Allele *{var_name} = NULL; /* empty */\n\n")
                return
            f.write(f"static const Allele {var_name}[] = {{\n")
            for a in alleles:
                # Escape any quotes in name
                name = a["name"].replace('"', '\\"')
                gene = a["gene"].replace('"', '\\"')
                family = a["family"].replace('"', '\\"')
                seg = segment_to_c(a["segment"])
                f.write(f'    {{ "{name}", "{gene}", "{family}",\n')
                f.write(f'      "{a["seq"]}",\n')
                f.write(f'      {a["length"]}, {a["anchor"]}, {seg} }},\n')
            f.write("};\n\n")

        write_alleles(c, v_alleles, "v_alleles")
        write_alleles(c, d_alleles, "d_alleles")
        write_alleles(c, j_alleles, "j_alleles")

        # ── Trim distributions ───────────────────────────────────
        def write_trim(f, probs, var_name):
            if not probs:
                f.write(f"static const int {var_name}_len = 0;\n")
                f.write(f"static const double *{var_name} = NULL;\n\n")
                return
            f.write(f"static const int {var_name}_len = {len(probs)};\n")
            f.write(f"static const double {var_name}[] = {{\n    ")
            for i, p in enumerate(probs):
                f.write(f"{p:.10f}")
                if i < len(probs) - 1:
                    f.write(", ")
                if (i + 1) % 8 == 0:
                    f.write("\n    ")
            f.write("\n};\n\n")

        write_trim(c, v_trim_3, "v_trim_3_probs")
        write_trim(c, j_trim_5, "j_trim_5_probs")
        write_trim(c, d_trim_5, "d_trim_5_probs")
        write_trim(c, d_trim_3, "d_trim_3_probs")

        # ── NP length distributions ──────────────────────────────
        def write_np_lengths(f, np_dict, var_name):
            max_len = max(np_dict.keys()) if np_dict else 0
            probs = [np_dict.get(i, 0.0) for i in range(max_len + 1)]
            total = sum(probs)
            if total > 0:
                probs = [p / total for p in probs]
            write_trim(f, probs, var_name)

        write_np_lengths(c, np1_lengths, "np1_length_probs")
        write_np_lengths(c, np2_lengths, "np2_length_probs")

        # ── S5F mutability (flat 3125 array) ─────────────────────
        c.write("static const double s5f_mutability_data[3125] = {\n    ")
        for i, val in enumerate(mutability_flat):
            c.write(f"{val:.10e}")
            if i < 3124:
                c.write(", ")
            if (i + 1) % 5 == 0:
                c.write("\n    ")
        c.write("\n};\n\n")

        # ── S5F substitution data ────────────────────────────────
        c.write(f"/* {len(sub_data)} 5-mers with substitution data */\n")
        c.write(f"static const int s5f_sub_count = {len(sub_data)};\n\n")

        c.write("static const int s5f_sub_keys[] = {\n    ")
        for i, (key, _, _) in enumerate(sub_data):
            c.write(f"{key}")
            if i < len(sub_data) - 1:
                c.write(", ")
            if (i + 1) % 16 == 0:
                c.write("\n    ")
        c.write("\n};\n\n")

        # Bases as strings (each up to 4 chars)
        c.write("static const char *s5f_sub_bases[] = {\n    ")
        for i, (_, bases, _) in enumerate(sub_data):
            c.write(f'"{bases}"')
            if i < len(sub_data) - 1:
                c.write(", ")
            if (i + 1) % 10 == 0:
                c.write("\n    ")
        c.write("\n};\n\n")

        # Counts
        c.write("static const int s5f_sub_counts[] = {\n    ")
        for i, (_, bases, _) in enumerate(sub_data):
            c.write(f"{len(bases)}")
            if i < len(sub_data) - 1:
                c.write(", ")
            if (i + 1) % 16 == 0:
                c.write("\n    ")
        c.write("\n};\n\n")

        # Weights (flattened, 4 per entry, padded with 0)
        c.write("static const double s5f_sub_weights[] = {\n    ")
        idx = 0
        for _, _, weights in sub_data:
            padded = weights + [0.0] * (4 - len(weights))
            for j, w in enumerate(padded):
                c.write(f"{w:.10e}")
                idx += 1
                if idx < len(sub_data) * 4:
                    c.write(", ")
            c.write("\n    ")
        c.write("\n};\n\n")

        # ── Load config function ─────────────────────────────────
        c.write(f"void {prefix}_load_config(SimConfig *cfg) {{\n")
        c.write(f"    sim_config_init(cfg, {c_chain});\n\n")

        # V alleles
        c.write("    for (int i = 0; i < v_alleles_count; i++)\n")
        c.write("        allele_pool_add(&cfg->v_alleles, &v_alleles[i]);\n\n")

        # D alleles
        if d_alleles:
            c.write("    for (int i = 0; i < d_alleles_count; i++)\n")
            c.write("        allele_pool_add(&cfg->d_alleles, &d_alleles[i]);\n\n")

        # J alleles
        c.write("    for (int i = 0; i < j_alleles_count; i++)\n")
        c.write("        allele_pool_add(&cfg->j_alleles, &j_alleles[i]);\n\n")

        # Trim distributions
        c.write("    /* Trim distributions */\n")
        for var, field in [("v_trim_3_probs", "v_trim_3"),
                           ("j_trim_5_probs", "j_trim_5"),
                           ("d_trim_5_probs", "d_trim_5"),
                           ("d_trim_3_probs", "d_trim_3")]:
            c.write(f"    if ({var}_len > 0) {{\n")
            c.write(f"        cfg->{field}.max_trim = {var}_len;\n")
            c.write(f"        cfg->{field}.probs = malloc({var}_len * sizeof(double));\n")
            c.write(f"        memcpy(cfg->{field}.probs, {var}, {var}_len * sizeof(double));\n")
            c.write(f"    }}\n")

        # NP parameters
        c.write(f"\n    cfg->np1_length_max = {np1_max};\n")
        c.write(f"    cfg->np2_length_max = {np2_max};\n")

        c.write("}\n\n")

        # ── Load S5F function ────────────────────────────────────
        c.write(f"void {prefix}_load_s5f(S5FModel *model) {{\n")
        c.write("    /* Copy mutability table */\n")
        c.write("    memcpy(model->mutability, s5f_mutability_data, sizeof(s5f_mutability_data));\n\n")
        c.write("    /* Load substitution data */\n")
        c.write("    for (int i = 0; i < s5f_sub_count; i++) {\n")
        c.write("        int key = s5f_sub_keys[i];\n")
        c.write("        s5f_set_substitution(model, key,\n")
        c.write("                              s5f_sub_bases[i],\n")
        c.write("                              &s5f_sub_weights[i * 4],\n")
        c.write("                              s5f_sub_counts[i]);\n")
        c.write("    }\n")
        c.write("}\n")

    print(f"Generated {h_path}")
    print(f"Generated {c_path}")
    print(f"  V alleles: {len(v_alleles)}")
    print(f"  D alleles: {len(d_alleles)}")
    print(f"  J alleles: {len(j_alleles)}")
    print(f"  V_3 trim dist: {len(v_trim_3)} entries")
    print(f"  J_5 trim dist: {len(j_trim_5)} entries")
    print(f"  D_5 trim dist: {len(d_trim_5)} entries")
    print(f"  D_3 trim dist: {len(d_trim_3)} entries")
    print(f"  S5F mutability: {sum(1 for v in mutability_flat if v > 0)} nonzero entries")
    print(f"  S5F substitution: {len(sub_data)} 5-mers")


def main():
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} CONFIG_NAME")
        print(f"Available: {', '.join(sorted(_CONFIG_NAMES))}")
        sys.exit(1)

    config_name = sys.argv[1].upper()
    if config_name not in _CONFIG_NAMES:
        print(f"Unknown config: {config_name}")
        print(f"Available: {', '.join(sorted(_CONFIG_NAMES))}")
        sys.exit(1)

    out_dir = os.path.join(os.path.dirname(__file__), "..", "data")
    os.makedirs(out_dir, exist_ok=True)

    print(f"Loading {config_name}...")
    dc = load_dataconfig(config_name)

    # Determine chain category for S5F
    chain_type = dc.metadata.chain_type
    if chain_type in (ChainType.BCR_HEAVY,):
        chain_category = "heavy"
    else:
        chain_category = "light"

    print(f"Loading S5F model ({chain_category})...")
    s5f_mut, s5f_sub = load_s5f_model(chain_category)

    print("Generating C source files...")
    write_c_files(config_name, dc, s5f_mut, s5f_sub, out_dir)
    print("Done.")


if __name__ == "__main__":
    main()
