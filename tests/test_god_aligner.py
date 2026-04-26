"""
God-Aligner Validation Tests.

These tests verify that the metadata of every simulated sequence is
exactly what an ideal aligner -- one that flawlessly extracts all
possible information from the sequence -- would report.

The core principle:
    metadata = f(sequence, germline_reference)

Tested invariants (10 checks):
  1. Segment boundary accuracy (V/D/J vs germline)
  2. Mutation accuracy (pos:from>to annotations)
  3. NP region accuracy (np1_region, np2_region extraction)
  4. Junction accuracy (junction_nt, junction_aa, junction_length)
  5. Trimming consistency (trim fields vs germline positions)
  6. Productivity oracle (vj_in_frame, stop_codon, productive)
  7. Germline alignment oracle (germline_alignment vs allele reference)
  8. Sequence reconstruction (rebuild from V+NP1+D+NP2+J + mutations)
  9. Coordinate consistency (non-overlapping, in-bounds)
 10. Ambiguous D calls (all listed alleles match equally)

Multi-config: human_igh, human_igk, human_tcrb, mouse_igh, rabbit_igh
"""

import pytest
from GenAIRR import Experiment
from GenAIRR.protocol import _resolve_config


# ====================================================================
# Helpers
# ====================================================================

CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}


def _flatten_alleles(allele_dict):
    """Flatten {gene: [Allele, ...]} to {name: ungapped_seq}."""
    out = {}
    if allele_dict is None:
        return out
    for _gene, alleles in allele_dict.items():
        for allele in alleles:
            out[allele.name] = allele.ungapped_seq
    return out


def _parse_mutations(mutation_str):
    """Parse 'pos:X>Y,pos:X>Y,...' into {pos: (from_base, to_base)}."""
    result = {}
    if not mutation_str:
        return result
    for entry in mutation_str.split(","):
        entry = entry.strip()
        if not entry or ":" not in entry:
            continue
        pos_str, change = entry.split(":", 1)
        try:
            pos = int(pos_str)
        except ValueError:
            continue
        if ">" not in change:
            continue
        parts = change.split(">", 1)
        if len(parts) == 2:
            result[pos] = (parts[0], parts[1])
    return result


def _translate_sequence(seq):
    """Translate nucleotide sequence from position 0, every 3 bases."""
    aas = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i + 3].upper()
        aas.append(CODON_TABLE.get(codon, 'X'))
    return ''.join(aas)


def _has_stop_codon_in_reading_frame(seq):
    """Check if sequence contains a stop codon in reading frame."""
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i + 3].upper()
        if CODON_TABLE.get(codon) == '*':
            return True
    return False


def _score_allele(seq, allele_seq, seq_start, seq_end, germ_start, germ_end):
    """Score how well an allele matches the sequence in the given region.

    Returns n_matches or None if invalid range.
    """
    if germ_start < 0 or germ_end <= germ_start:
        return None
    region_len = min(seq_end - seq_start, germ_end - germ_start)
    if region_len <= 0:
        return None
    n_match = 0
    for i in range(region_len):
        si = seq_start + i
        gi = germ_start + i
        s = seq[si].upper() if si < len(seq) else 'N'
        if gi < len(allele_seq):
            a = allele_seq[gi].upper()
        else:
            continue
        if s == 'N' or a == 'N':
            continue
        if s == a:
            n_match += 1
    return n_match


def _has_d_segment(config_name):
    """Return True if the config has a D gene segment."""
    return config_name in ("human_igh", "mouse_igh", "rabbit_igh",
                           "human_tcrb")


# ====================================================================
# Config definitions
# ====================================================================

CONFIGS = [
    "human_igh",
    "human_igk",
    "human_tcrb",
    "mouse_igh",
    "rabbit_igh",
]


# ====================================================================
# Fixtures -- module-scoped for efficiency
# ====================================================================

@pytest.fixture(scope="module")
def all_configs():
    """Load all DataConfig objects and flatten alleles."""
    configs = {}
    for name in CONFIGS:
        dc = _resolve_config(name)
        configs[name] = {
            "dc": dc,
            "v": _flatten_alleles(dc.v_alleles),
            "d": _flatten_alleles(dc.d_alleles),
            "j": _flatten_alleles(dc.j_alleles),
            "has_d": _has_d_segment(name),
        }
    return configs


@pytest.fixture(scope="module")
def mutated_records(all_configs):  # noqa: ARG001
    """100 mutated sequences per config."""
    from GenAIRR.ops import rate
    results = {}
    for name in CONFIGS:
        results[name] = (Experiment.on(name)
                         .mutate(rate(0.02, 0.08))
                         .run(n=100, seed=42))
    return results


@pytest.fixture(scope="module")
def noisy_records():
    """50 sequences with mutation + corruption + indels + Ns for heavy chains."""
    from GenAIRR.ops import rate, with_5prime_loss, with_indels, with_ns
    results = {}
    for name in ("human_igh", "human_tcrb", "mouse_igh"):
        results[name] = (Experiment.on(name)
                         .mutate(rate(0.02, 0.08))
                         .sequence(with_5prime_loss())
                         .observe(with_indels(0.01), with_ns(0.005))
                         .run(n=50, seed=42))
    return results


# ====================================================================
# 1. Segment Boundary Accuracy
# ====================================================================

class TestSegmentBoundaryAccuracy:
    """Extract sequence[v_sequence_start:v_sequence_end] and align to
    the V germline allele. Every non-mutated position must match."""

    def test_v_segment_matches_germline(self, mutated_records, all_configs):
        """V segment of the sequence matches the V allele at non-mutated positions.

        Uses `v_call_true` (the truly-sampled allele) rather than `v_call`
        because the latter is god-aligner-derived and can drift to a
        different allele under heavy SHM, breaking AIRR self-consistency
        with the recorded `mutations` field (which is computed against
        the truly-sampled allele's germline)."""
        for config_name in CONFIGS:
            v_alleles = all_configs[config_name]["v"]
            for i, rec in enumerate(mutated_records[config_name]):
                seq = rec["sequence"]
                v_start = rec["v_sequence_start"]
                v_end = rec["v_sequence_end"]
                g_start = rec["v_germline_start"]
                v_call = rec.get("v_call_true") or rec["v_call"].split(",")[0]
                if v_call not in v_alleles or v_end <= v_start:
                    continue
                allele_seq = v_alleles[v_call]
                mutations = _parse_mutations(rec.get("mutations", ""))
                for offset in range(v_end - v_start):
                    spos = v_start + offset
                    gpos = g_start + offset
                    if gpos >= len(allele_seq) or spos >= len(seq):
                        break
                    if spos in mutations:
                        continue  # mutated position, skip
                    s = seq[spos].upper()
                    a = allele_seq[gpos].upper()
                    if s == 'N' or a == 'N':
                        continue
                    assert s == a, (
                        f"[{config_name} seq {i}] V segment mismatch at "
                        f"seq[{spos}]='{s}' vs allele[{gpos}]='{a}' "
                        f"for {v_call}")

    def test_j_segment_matches_germline(self, mutated_records, all_configs):
        """J segment of the sequence matches the J allele at non-mutated positions."""
        for config_name in CONFIGS:
            j_alleles = all_configs[config_name]["j"]
            for i, rec in enumerate(mutated_records[config_name]):
                seq = rec["sequence"]
                j_start = rec["j_sequence_start"]
                j_end = rec["j_sequence_end"]
                g_start = rec["j_germline_start"]
                j_call = rec["j_call"].split(",")[0]
                if j_call not in j_alleles or j_end <= j_start:
                    continue
                allele_seq = j_alleles[j_call]
                mutations = _parse_mutations(rec.get("mutations", ""))
                for offset in range(j_end - j_start):
                    spos = j_start + offset
                    gpos = g_start + offset
                    if gpos >= len(allele_seq) or spos >= len(seq):
                        break
                    if spos in mutations:
                        continue
                    s = seq[spos].upper()
                    a = allele_seq[gpos].upper()
                    if s == 'N' or a == 'N':
                        continue
                    assert s == a, (
                        f"[{config_name} seq {i}] J segment mismatch at "
                        f"seq[{spos}]='{s}' vs allele[{gpos}]='{a}' "
                        f"for {j_call}")

    def test_d_segment_matches_germline(self, mutated_records, all_configs):
        """D segment matches when present (heavy chains, TCR beta).

        Uses d_call_true (truly-sampled allele) rather than d_call,
        which is god-aligner-derived and can drift to a different
        allele under heavy SHM (recorded mutations are computed
        against the truly-sampled allele's germline)."""
        for config_name in CONFIGS:
            if not all_configs[config_name]["has_d"]:
                continue
            d_alleles = all_configs[config_name]["d"]
            for i, rec in enumerate(mutated_records[config_name]):
                d_start = rec["d_sequence_start"]
                d_end = rec["d_sequence_end"]
                if d_start == d_end or not rec["d_call"]:
                    continue
                # Prefer the truly-sampled allele (d_call_true) for the
                # self-consistency check; fall back to the first listed
                # d_call when d_call_true is absent (older snapshots).
                d_call = rec.get("d_call_true") or rec["d_call"].split(",")[0]
                if d_call == "Short-D" or d_call not in d_alleles:
                    continue
                seq = rec["sequence"]
                g_start = rec["d_germline_start"]
                allele_seq = d_alleles[d_call]
                mutations = _parse_mutations(rec.get("mutations", ""))
                for offset in range(d_end - d_start):
                    spos = d_start + offset
                    gpos = g_start + offset
                    if gpos >= len(allele_seq) or spos >= len(seq):
                        break
                    if spos in mutations:
                        continue
                    s = seq[spos].upper()
                    a = allele_seq[gpos].upper()
                    if s == 'N' or a == 'N':
                        continue
                    assert s == a, (
                        f"[{config_name} seq {i}] D segment mismatch at "
                        f"seq[{spos}]='{s}' vs allele[{gpos}]='{a}' "
                        f"for {d_call}")

    def test_segment_span_equals_germline_span(self, mutated_records):
        """sequence span == germline span for each segment."""
        for config_name in CONFIGS:
            for i, rec in enumerate(mutated_records[config_name]):
                v_seq_span = rec["v_sequence_end"] - rec["v_sequence_start"]
                v_germ_span = rec["v_germline_end"] - rec["v_germline_start"]
                assert v_seq_span == v_germ_span, (
                    f"[{config_name} seq {i}] V seq span {v_seq_span} "
                    f"!= V germ span {v_germ_span}")

                if rec["d_sequence_start"] != rec["d_sequence_end"]:
                    d_seq_span = rec["d_sequence_end"] - rec["d_sequence_start"]
                    d_germ_span = rec["d_germline_end"] - rec["d_germline_start"]
                    assert d_seq_span == d_germ_span, (
                        f"[{config_name} seq {i}] D seq span {d_seq_span} "
                        f"!= D germ span {d_germ_span}")

                j_seq_span = rec["j_sequence_end"] - rec["j_sequence_start"]
                j_germ_span = rec["j_germline_end"] - rec["j_germline_start"]
                assert j_seq_span == j_germ_span, (
                    f"[{config_name} seq {i}] J seq span {j_seq_span} "
                    f"!= J germ span {j_germ_span}")


# ====================================================================
# 2. Mutation Accuracy
# ====================================================================

class TestMutationAccuracy:
    """Parse the mutations string and verify every annotation."""

    def test_mutation_to_base_matches_sequence(self, mutated_records):
        """For each mutation pos:X>Y, sequence[pos] == Y."""
        for config_name in CONFIGS:
            for i, rec in enumerate(mutated_records[config_name]):
                seq = rec["sequence"]
                mutations = _parse_mutations(rec.get("mutations", ""))
                for pos, (_, to_base) in mutations.items():
                    assert pos < len(seq), (
                        f"[{config_name} seq {i}] mutation pos {pos} "
                        f"out of bounds (len={len(seq)})")
                    assert seq[pos] == to_base, (
                        f"[{config_name} seq {i}] mutation says pos {pos} "
                        f"is '{to_base}' but sequence has '{seq[pos]}'")

    def test_mutation_from_base_matches_germline(self, mutated_records):
        """For each mutation pos:X>Y, germline_alignment[pos] == X."""
        for config_name in CONFIGS:
            for i, rec in enumerate(mutated_records[config_name]):
                germ = rec["germline_alignment"]
                mutations = _parse_mutations(rec.get("mutations", ""))
                for pos, (from_base, _to_base) in mutations.items():
                    if pos >= len(germ):
                        continue
                    assert germ[pos] == from_base, (
                        f"[{config_name} seq {i}] mutation says germline at "
                        f"{pos} is '{from_base}' but germline_alignment "
                        f"has '{germ[pos]}'")

    def test_n_mutations_equals_annotation_count(self, mutated_records):
        """n_mutations == len(parsed mutations)."""
        for config_name in CONFIGS:
            for i, rec in enumerate(mutated_records[config_name]):
                mutations = _parse_mutations(rec.get("mutations", ""))
                assert rec["n_mutations"] == len(mutations), (
                    f"[{config_name} seq {i}] n_mutations={rec['n_mutations']} "
                    f"but {len(mutations)} annotations parsed")

    def test_mutation_rate_approximate(self, mutated_records):
        """mutation_rate approximately equals n_mutations / germline_length."""
        for config_name in CONFIGS:
            for i, rec in enumerate(mutated_records[config_name]):
                n_mut = rec["n_mutations"]
                rate = rec["mutation_rate"]
                v_span = rec["v_sequence_end"] - rec["v_sequence_start"]
                j_span = rec["j_sequence_end"] - rec["j_sequence_start"]
                d_span = rec["d_sequence_end"] - rec["d_sequence_start"]
                total = v_span + j_span + d_span
                if total == 0:
                    continue
                expected = n_mut / total
                assert abs(rate - expected) < 0.05, (
                    f"[{config_name} seq {i}] stored rate={rate:.4f}, "
                    f"computed={expected:.4f} ({n_mut}/{total})")

    def test_no_mutation_at_matching_positions(self, mutated_records):
        """If sequence[pos] == germline[pos], no mutation at that pos."""
        for config_name in CONFIGS:
            for i, rec in enumerate(mutated_records[config_name]):
                seq = rec["sequence"]
                germ = rec["germline_alignment"]
                mutations = _parse_mutations(rec.get("mutations", ""))
                for pos in mutations:
                    if pos >= len(seq) or pos >= len(germ):
                        continue
                    assert seq[pos] != germ[pos], (
                        f"[{config_name} seq {i}] mutation at pos {pos} "
                        f"but seq='{seq[pos]}' == germline='{germ[pos]}'")


# ====================================================================
# 3. NP Region Accuracy
# ====================================================================

class TestNPRegionAccuracy:
    """Verify NP region fields match the sequence between segments."""

    def test_np1_region_matches_sequence(self, mutated_records):
        """sequence[v_end:d_start] == np1_region (heavy chains).
        For light chains: sequence[v_end:j_start] == np1_region."""
        for config_name in CONFIGS:
            has_d = _has_d_segment(config_name)
            for i, rec in enumerate(mutated_records[config_name]):
                v_end = rec["v_sequence_end"]
                if has_d:
                    np1_end = rec["d_sequence_start"]
                else:
                    np1_end = rec["j_sequence_start"]
                seq = rec["sequence"]
                np1 = rec["np1_region"]
                extracted = seq[v_end:np1_end]
                assert extracted.upper() == np1.upper(), (
                    f"[{config_name} seq {i}] NP1: "
                    f"seq[{v_end}:{np1_end}]='{extracted}' != "
                    f"np1_region='{np1}'")

    def test_np1_length_field(self, mutated_records):
        """np1_length == len(np1_region)."""
        for config_name in CONFIGS:
            for i, rec in enumerate(mutated_records[config_name]):
                assert rec["np1_length"] == len(rec["np1_region"]), (
                    f"[{config_name} seq {i}] np1_length={rec['np1_length']} "
                    f"!= len(np1_region)={len(rec['np1_region'])}")

    def test_np2_region_matches_sequence(self, mutated_records):
        """sequence[d_end:j_start] == np2_region (heavy chains only)."""
        for config_name in CONFIGS:
            if not _has_d_segment(config_name):
                continue
            for i, rec in enumerate(mutated_records[config_name]):
                d_end = rec["d_sequence_end"]
                j_start = rec["j_sequence_start"]
                seq = rec["sequence"]
                np2 = rec["np2_region"]
                extracted = seq[d_end:j_start]
                assert extracted.upper() == np2.upper(), (
                    f"[{config_name} seq {i}] NP2: "
                    f"seq[{d_end}:{j_start}]='{extracted}' != "
                    f"np2_region='{np2}'")

    def test_np2_length_field(self, mutated_records):
        """np2_length == len(np2_region)."""
        for config_name in CONFIGS:
            for i, rec in enumerate(mutated_records[config_name]):
                assert rec["np2_length"] == len(rec["np2_region"]), (
                    f"[{config_name} seq {i}] np2_length={rec['np2_length']} "
                    f"!= len(np2_region)={len(rec['np2_region'])}")

    def test_vj_chain_np2_empty(self, mutated_records):
        """Light chains (no D) should have empty NP2."""
        for config_name in CONFIGS:
            if _has_d_segment(config_name):
                continue
            for i, rec in enumerate(mutated_records[config_name]):
                assert rec["np2_region"] == "", (
                    f"[{config_name} seq {i}] VJ chain has non-empty "
                    f"np2_region='{rec['np2_region']}'")
                assert rec["np2_length"] == 0


# ====================================================================
# 4. Junction Accuracy
# ====================================================================

class TestJunctionAccuracy:
    """Verify junction_nt, junction_aa, junction_length consistency."""

    def test_junction_nt_matches_sequence_slice(self, mutated_records):
        """junction_nt == sequence[junction_start:junction_end]."""
        for config_name in CONFIGS:
            for i, rec in enumerate(mutated_records[config_name]):
                js = rec["junction_start"]
                je = rec["junction_end"]
                seq = rec["sequence"]
                junction_nt = rec.get("junction_nt", "")
                if not junction_nt:
                    continue
                extracted = seq[js:je]
                assert extracted == junction_nt, (
                    f"[{config_name} seq {i}] junction_nt mismatch: "
                    f"seq[{js}:{je}]='{extracted[:30]}...' != "
                    f"junction_nt='{junction_nt[:30]}...'")

    def test_junction_length_field(self, mutated_records):
        """junction_length == junction_end - junction_start."""
        for config_name in CONFIGS:
            for i, rec in enumerate(mutated_records[config_name]):
                jlen = rec["junction_end"] - rec["junction_start"]
                assert rec["junction_length"] == jlen, (
                    f"[{config_name} seq {i}] junction_length="
                    f"{rec['junction_length']} != "
                    f"{rec['junction_end']}-{rec['junction_start']}={jlen}")

    def test_junction_aa_matches_translation(self, mutated_records):
        """Independently translate junction_nt and compare to junction_aa."""
        for config_name in CONFIGS:
            for i, rec in enumerate(mutated_records[config_name]):
                junction_nt = rec.get("junction_nt", "")
                junction_aa = rec.get("junction_aa", "")
                if not junction_nt or not junction_aa:
                    continue
                expected_aa = _translate_sequence(junction_nt)
                assert junction_aa == expected_aa, (
                    f"[{config_name} seq {i}] junction_aa mismatch:\n"
                    f"  stored:   '{junction_aa}'\n"
                    f"  computed: '{expected_aa}'\n"
                    f"  junction_nt: '{junction_nt}'")

    def test_junction_aa_length(self, mutated_records):
        """len(junction_aa) == junction_length // 3."""
        for config_name in CONFIGS:
            for i, rec in enumerate(mutated_records[config_name]):
                junction_aa = rec.get("junction_aa", "")
                jlen = rec["junction_length"]
                if junction_aa:
                    expected = jlen // 3
                    assert len(junction_aa) == expected, (
                        f"[{config_name} seq {i}] junction_aa length "
                        f"{len(junction_aa)} != junction_length//3={expected}")


# ====================================================================
# 5. Trimming Consistency
# ====================================================================

class TestTrimmingConsistency:
    """Verify trim fields are consistent with germline positions.

    Note: Boundary correction (NP reconstruction / D-J ambiguity resolution)
    can adjust germline boundaries after the initial trim. The trim fields
    record the *original* exonuclease trim, while germline_start/end reflect
    the *corrected* boundaries. Therefore:
      - v_germline_end <= allele_len - v_trim_3 (correction can only shrink or keep)
      - For most sequences they are equal; a small fraction differ by 1-3 bp
    """

    def test_v_germline_end_within_allele_bounds(self, mutated_records, all_configs):
        """v_germline_end must not exceed allele length."""
        for config_name in CONFIGS:
            v_alleles = all_configs[config_name]["v"]
            for i, rec in enumerate(mutated_records[config_name]):
                v_call = rec["v_call"].split(",")[0]
                if v_call not in v_alleles:
                    continue
                allele_len = len(v_alleles[v_call])
                assert rec["v_germline_end"] <= allele_len, (
                    f"[{config_name} seq {i}] v_germline_end="
                    f"{rec['v_germline_end']} > allele_len={allele_len}")

    def test_v_germline_start_is_zero(self, mutated_records):
        """V is never trimmed at 5', so v_germline_start should be 0."""
        for config_name in CONFIGS:
            for i, rec in enumerate(mutated_records[config_name]):
                assert rec["v_germline_start"] == 0, (
                    f"[{config_name} seq {i}] v_germline_start="
                    f"{rec['v_germline_start']} (expected 0)")

    def test_v_trim_5_is_zero(self, mutated_records):
        """V genes are never trimmed at 5'."""
        for config_name in CONFIGS:
            for i, rec in enumerate(mutated_records[config_name]):
                assert rec["v_trim_5"] == 0, (
                    f"[{config_name} seq {i}] v_trim_5={rec['v_trim_5']}")

    def test_j_germline_start_equals_trim_5(self, mutated_records):
        """j_germline_start should equal j_trim_5."""
        for config_name in CONFIGS:
            for i, rec in enumerate(mutated_records[config_name]):
                assert rec["j_germline_start"] == rec["j_trim_5"], (
                    f"[{config_name} seq {i}] j_germline_start="
                    f"{rec['j_germline_start']} != j_trim_5={rec['j_trim_5']}")

    def test_j_germline_end_within_allele_bounds(self, mutated_records, all_configs):
        """j_germline_end must not exceed allele length."""
        for config_name in CONFIGS:
            j_alleles = all_configs[config_name]["j"]
            for i, rec in enumerate(mutated_records[config_name]):
                j_call = rec["j_call"].split(",")[0]
                if j_call not in j_alleles:
                    continue
                allele_len = len(j_alleles[j_call])
                assert rec["j_germline_end"] <= allele_len, (
                    f"[{config_name} seq {i}] j_germline_end="
                    f"{rec['j_germline_end']} > allele_len={allele_len}")

    def test_j_trim_3_is_zero(self, mutated_records):
        """J genes are never trimmed at 3'."""
        for config_name in CONFIGS:
            for i, rec in enumerate(mutated_records[config_name]):
                assert rec["j_trim_3"] == 0, (
                    f"[{config_name} seq {i}] j_trim_3={rec['j_trim_3']}")

    def test_d_germline_start_equals_trim_5(self, mutated_records):
        """d_germline_start should equal d_trim_5."""
        for config_name in CONFIGS:
            if not _has_d_segment(config_name):
                continue
            for i, rec in enumerate(mutated_records[config_name]):
                if rec["d_sequence_start"] == rec["d_sequence_end"]:
                    continue
                assert rec["d_germline_start"] == rec["d_trim_5"], (
                    f"[{config_name} seq {i}] d_germline_start="
                    f"{rec['d_germline_start']} != d_trim_5={rec['d_trim_5']}")

    def test_d_germline_end_within_allele_bounds(self, mutated_records, all_configs):
        """d_germline_end must not exceed allele length."""
        for config_name in CONFIGS:
            if not all_configs[config_name]["has_d"]:
                continue
            d_alleles = all_configs[config_name]["d"]
            for i, rec in enumerate(mutated_records[config_name]):
                if rec["d_sequence_start"] == rec["d_sequence_end"]:
                    continue
                d_call = rec["d_call"].split(",")[0]
                if d_call == "Short-D" or d_call not in d_alleles:
                    continue
                allele_len = len(d_alleles[d_call])
                assert rec["d_germline_end"] <= allele_len, (
                    f"[{config_name} seq {i}] d_germline_end="
                    f"{rec['d_germline_end']} > allele_len={allele_len}")

    def test_trims_nonnegative(self, mutated_records):
        """All trim values must be >= 0."""
        for config_name in CONFIGS:
            for i, rec in enumerate(mutated_records[config_name]):
                for field in ("v_trim_3", "v_trim_5", "j_trim_3", "j_trim_5",
                              "d_trim_3", "d_trim_5"):
                    assert rec[field] >= 0, (
                        f"[{config_name} seq {i}] {field}={rec[field]} < 0")


# ====================================================================
# 6. Productivity Oracle
# ====================================================================

class TestProductivityOracle:
    """Independently verify productivity metadata."""

    def test_stop_codon_matches_sequence_scan(self, mutated_records):
        """Scan sequence in reading frame for stop codons."""
        for config_name in CONFIGS:
            for i, rec in enumerate(mutated_records[config_name]):
                seq = rec["sequence"].upper()
                has_stop = _has_stop_codon_in_reading_frame(seq)
                assert rec["stop_codon"] == has_stop, (
                    f"[{config_name} seq {i}] Python scan stop_codon="
                    f"{has_stop}, stored={rec['stop_codon']}")

    def test_productive_implies_no_stop_and_in_frame(self, mutated_records):
        """productive == True requires stop_codon=False AND vj_in_frame=True."""
        for config_name in CONFIGS:
            for i, rec in enumerate(mutated_records[config_name]):
                if rec["productive"]:
                    assert not rec["stop_codon"], (
                        f"[{config_name} seq {i}] productive but stop_codon")
                    assert rec["vj_in_frame"], (
                        f"[{config_name} seq {i}] productive but not in frame")

    def test_stop_codon_implies_not_productive(self, mutated_records):
        """stop_codon == True requires productive == False."""
        for config_name in CONFIGS:
            for i, rec in enumerate(mutated_records[config_name]):
                if rec["stop_codon"]:
                    assert not rec["productive"], (
                        f"[{config_name} seq {i}] stop_codon but productive")

    def test_vj_in_frame_geometry(self, mutated_records):
        """vj_in_frame=True requires junction geometrically in frame."""
        for config_name in CONFIGS:
            for i, rec in enumerate(mutated_records[config_name]):
                js = rec["junction_start"]
                je = rec["junction_end"]
                jlen = je - js
                geom = (js % 3 == 0) and (je % 3 == 0) and (jlen % 3 == 0)
                if rec["vj_in_frame"]:
                    assert geom, (
                        f"[{config_name} seq {i}] vj_in_frame=True but "
                        f"junction [{js}:{je}] len {jlen} not in frame")
                if not geom:
                    assert not rec["vj_in_frame"], (
                        f"[{config_name} seq {i}] junction out of frame "
                        f"but vj_in_frame=True")

    def test_productive_anchors(self, mutated_records):
        """Productive sequences: junction_aa starts with C, ends with W/F."""
        for config_name in CONFIGS:
            for i, rec in enumerate(mutated_records[config_name]):
                if not rec["productive"]:
                    continue
                jaa = rec.get("junction_aa", "")
                if not jaa:
                    continue
                assert jaa[0] == 'C', (
                    f"[{config_name} seq {i}] productive but "
                    f"junction_aa[0]='{jaa[0]}', expected 'C'")
                assert jaa[-1] in ('W', 'F'), (
                    f"[{config_name} seq {i}] productive but "
                    f"junction_aa[-1]='{jaa[-1]}', expected W/F")


# ====================================================================
# 7. Germline Alignment Oracle
# ====================================================================

class TestGermlineAlignmentOracle:
    """Validate germline_alignment against allele reference."""

    def test_germline_alignment_length(self, mutated_records):
        """germline_alignment length == sequence length."""
        for config_name in CONFIGS:
            for i, rec in enumerate(mutated_records[config_name]):
                assert len(rec["germline_alignment"]) == len(rec["sequence"]), (
                    f"[{config_name} seq {i}] germline_alignment length "
                    f"{len(rec['germline_alignment'])} != "
                    f"sequence length {len(rec['sequence'])}")

    def test_germline_has_bases_in_vdj(self, mutated_records):
        """Within V/D/J spans, germline must have real bases (ACGT).

        Some IMGT alleles contain IUPAC ambiguity codes (N, Y, etc.)
        at specific positions. These are faithfully reproduced in the
        germline alignment and are acceptable.
        """
        # IUPAC nucleotide codes that can appear in IMGT alleles
        iupac_bases = set('ACGTRYSWKMBDHVN')
        for config_name in CONFIGS:
            for i, rec in enumerate(mutated_records[config_name]):
                germ = rec["germline_alignment"].upper()
                spans = [
                    ("V", rec["v_sequence_start"], rec["v_sequence_end"]),
                    ("J", rec["j_sequence_start"], rec["j_sequence_end"]),
                ]
                d_s, d_e = rec["d_sequence_start"], rec["d_sequence_end"]
                if d_s != d_e:
                    spans.append(("D", d_s, d_e))
                for seg, start, end in spans:
                    for pos in range(start, min(end, len(germ))):
                        assert germ[pos] in iupac_bases, (
                            f"[{config_name} seq {i}] germline[{pos}]="
                            f"'{germ[pos]}' in {seg} span, expected "
                            f"nucleotide base")

    def test_germline_has_n_in_np_regions(self, mutated_records):
        """Between V and D (NP1), germline must have N."""
        for config_name in CONFIGS:
            has_d = _has_d_segment(config_name)
            for i, rec in enumerate(mutated_records[config_name]):
                germ = rec["germline_alignment"]
                v_end = rec["v_sequence_end"]
                if has_d:
                    np1_end = rec["d_sequence_start"]
                else:
                    np1_end = rec["j_sequence_start"]
                for pos in range(v_end, min(np1_end, len(germ))):
                    assert germ[pos] == 'N', (
                        f"[{config_name} seq {i}] germline[{pos}]='{germ[pos]}'"
                        f" in NP1 region, expected 'N'")

    def test_v_germline_matches_allele_reference(self, mutated_records, all_configs):
        """germline_alignment in V region must match the V allele reference.

        Note: v_call reflects the god-aligner's best match to the *mutated*
        sequence, which may differ from the original allele. The germline
        alignment stores the original allele's bases. We therefore check
        unmutated sequences only, or skip positions with mutations.
        """
        for config_name in CONFIGS:
            v_alleles = all_configs[config_name]["v"]
            for i, rec in enumerate(mutated_records[config_name]):
                germ = rec["germline_alignment"]
                # Use v_call_true (truly-sampled allele) rather than the
                # god-aligner-derived v_call, which can drift away from
                # the recorded mutations under heavy SHM.
                v_call = rec.get("v_call_true") or rec["v_call"].split(",")[0]
                if v_call not in v_alleles:
                    continue
                allele_seq = v_alleles[v_call]
                v_start = rec["v_sequence_start"]
                v_end = rec["v_sequence_end"]
                g_start = rec["v_germline_start"]
                mutations = _parse_mutations(rec.get("mutations", ""))
                for offset in range(v_end - v_start):
                    spos = v_start + offset
                    gpos = g_start + offset
                    if gpos >= len(allele_seq) or spos >= len(germ):
                        break
                    if spos in mutations:
                        continue  # mutation may have shifted the call
                    g = germ[spos].upper()
                    a = allele_seq[gpos].upper()
                    if g == 'N' or a == 'N':
                        continue
                    assert g == a, (
                        f"[{config_name} seq {i}] germline[{spos}]='{g}' "
                        f"!= allele[{gpos}]='{a}' for {v_call}")


# ====================================================================
# 8. Sequence Reconstruction
# ====================================================================

class TestSequenceReconstruction:
    """Reconstruct the full sequence from germline + NP + mutations."""

    def test_reconstruct_unmutated(self, all_configs):
        """For unmutated sequences, rebuild from V+NP1+D+NP2+J and compare."""
        for config_name in CONFIGS:
            records = Experiment.on(config_name).run(n=50, seed=42)
            v_alleles = all_configs[config_name]["v"]
            d_alleles = all_configs[config_name]["d"]
            j_alleles = all_configs[config_name]["j"]
            has_d = all_configs[config_name]["has_d"]

            for i, rec in enumerate(records):
                seq = rec["sequence"]
                v_call = rec["v_call"].split(",")[0]
                j_call = rec["j_call"].split(",")[0]
                if v_call not in v_alleles or j_call not in j_alleles:
                    continue

                # Build V portion
                v_allele = v_alleles[v_call]
                v_part = v_allele[rec["v_germline_start"]:rec["v_germline_end"]]

                # Build NP1
                np1 = rec["np1_region"]

                # Build D portion
                if has_d and rec["d_sequence_start"] != rec["d_sequence_end"]:
                    d_call = rec["d_call"].split(",")[0]
                    if d_call != "Short-D" and d_call in d_alleles:
                        d_allele = d_alleles[d_call]
                        d_part = d_allele[rec["d_germline_start"]:rec["d_germline_end"]]
                    else:
                        # Short-D: extract from sequence directly
                        d_part = seq[rec["d_sequence_start"]:rec["d_sequence_end"]]
                else:
                    d_part = ""

                # Build NP2
                np2 = rec["np2_region"]

                # Build J portion
                j_allele = j_alleles[j_call]
                j_part = j_allele[rec["j_germline_start"]:rec["j_germline_end"]]

                # Reconstruct
                reconstructed = v_part + np1 + d_part + np2 + j_part

                # Compare case-insensitively
                assert reconstructed.upper() == seq.upper(), (
                    f"[{config_name} seq {i}] Reconstruction mismatch:\n"
                    f"  reconstructed({len(reconstructed)}): "
                    f"'{reconstructed[:50]}...'\n"
                    f"  sequence({len(seq)}): '{seq[:50]}...'")

    def test_reconstruct_from_germline_alignment(self, mutated_records):
        """Rebuild sequence from germline_alignment + mutations.

        The germline_alignment stores the original allele bases at V/D/J
        positions and 'N' at NP regions. Applying the mutations to the
        germline_alignment should reproduce the sequence exactly.

        Note: allele calls may differ from the original allele after
        mutation shifts the best match. We use germline_alignment
        (which always reflects the original allele) rather than
        re-deriving from the call.
        """
        for config_name in CONFIGS:
            for i, rec in enumerate(mutated_records[config_name]):
                seq = rec["sequence"]
                germ = rec["germline_alignment"]

                if len(germ) != len(seq):
                    continue  # skip if lengths differ (shouldn't happen)

                # Start with germline, replace NP regions with actual sequence
                reconstructed = list(germ)

                # NP1: germline has 'N', sequence has the actual NP bases
                has_d = _has_d_segment(config_name)
                v_end = rec["v_sequence_end"]
                if has_d:
                    np1_end = rec["d_sequence_start"]
                else:
                    np1_end = rec["j_sequence_start"]
                for pos in range(v_end, np1_end):
                    if pos < len(reconstructed):
                        reconstructed[pos] = seq[pos]

                # NP2: only for chains with D
                if has_d:
                    d_end = rec["d_sequence_end"]
                    j_start = rec["j_sequence_start"]
                    for pos in range(d_end, j_start):
                        if pos < len(reconstructed):
                            reconstructed[pos] = seq[pos]

                # Apply mutations (they replace germline bases)
                mutations = _parse_mutations(rec.get("mutations", ""))
                for pos, (_from_base, to_base) in mutations.items():
                    if pos < len(reconstructed):
                        reconstructed[pos] = to_base

                reconstructed = ''.join(reconstructed)
                assert reconstructed.upper() == seq.upper(), (
                    f"[{config_name} seq {i}] Reconstruction mismatch:\n"
                    f"  first diff at: "
                    f"{next((j for j in range(min(len(reconstructed), len(seq))) if reconstructed[j].upper() != seq[j].upper()), 'N/A')}")


# ====================================================================
# 9. Coordinate Consistency
# ====================================================================

class TestCoordinateConsistency:
    """Verify all coordinates are self-consistent."""

    def test_segments_non_overlapping(self, mutated_records):
        """V ends before D starts, D ends before J starts."""
        for config_name in CONFIGS:
            for i, rec in enumerate(mutated_records[config_name]):
                d_start = rec["d_sequence_start"]
                d_end = rec["d_sequence_end"]
                if d_start != d_end:
                    assert rec["v_sequence_end"] <= d_start, (
                        f"[{config_name} seq {i}] V end "
                        f"({rec['v_sequence_end']}) > D start ({d_start})")
                    assert d_end <= rec["j_sequence_start"], (
                        f"[{config_name} seq {i}] D end ({d_end}) > "
                        f"J start ({rec['j_sequence_start']})")
                assert rec["v_sequence_end"] <= rec["j_sequence_start"], (
                    f"[{config_name} seq {i}] V end "
                    f"({rec['v_sequence_end']}) > J start "
                    f"({rec['j_sequence_start']})")

    def test_all_positions_in_bounds(self, mutated_records):
        """All coordinates within [0, sequence_length]."""
        for config_name in CONFIGS:
            for i, rec in enumerate(mutated_records[config_name]):
                slen = len(rec["sequence"])
                for field in ("v_sequence_start", "v_sequence_end",
                              "d_sequence_start", "d_sequence_end",
                              "j_sequence_start", "j_sequence_end",
                              "junction_start", "junction_end"):
                    val = rec[field]
                    assert 0 <= val <= slen, (
                        f"[{config_name} seq {i}] {field}={val} "
                        f"out of [0, {slen}]")

    def test_start_before_end(self, mutated_records):
        """start <= end for all segment pairs."""
        for config_name in CONFIGS:
            for _i, rec in enumerate(mutated_records[config_name]):
                assert rec["v_sequence_start"] <= rec["v_sequence_end"]
                assert rec["j_sequence_start"] <= rec["j_sequence_end"]
                assert rec["junction_start"] <= rec["junction_end"]
                # junction_start == junction_end == 0 is valid when
                # anchors are missing (anchorless alleles from IMGT)
                if rec["junction_start"] != 0 or rec["junction_end"] != 0:
                    assert rec["junction_start"] < rec["junction_end"]
                d_s, d_e = rec["d_sequence_start"], rec["d_sequence_end"]
                if d_s != d_e:
                    assert d_s <= d_e

    def test_sequence_length_field(self, mutated_records):
        """sequence_length == len(sequence)."""
        for config_name in CONFIGS:
            for i, rec in enumerate(mutated_records[config_name]):
                assert rec["sequence_length"] == len(rec["sequence"]), (
                    f"[{config_name} seq {i}] sequence_length="
                    f"{rec['sequence_length']} != len={len(rec['sequence'])}")

    def test_v_starts_at_zero_without_corruption(self, mutated_records):
        """Without 5' corruption, V should start at position 0."""
        for config_name in CONFIGS:
            for i, rec in enumerate(mutated_records[config_name]):
                assert rec["v_sequence_start"] == 0, (
                    f"[{config_name} seq {i}] v_sequence_start="
                    f"{rec['v_sequence_start']} (expected 0, no corruption)")

    def test_valid_nucleotides_only(self, mutated_records):
        """Sequences contain only valid nucleotides."""
        for config_name in CONFIGS:
            for i, rec in enumerate(mutated_records[config_name]):
                invalid = set(rec["sequence"].upper()) - set('ATCGN')
                assert not invalid, (
                    f"[{config_name} seq {i}] invalid chars: {invalid}")


# ====================================================================
# 10. Ambiguous D Calls
# ====================================================================

class TestAmbiguousDCalls:
    """When d_call contains commas, all listed alleles must match equally."""

    def test_all_ambiguous_d_alleles_match_equally(self, mutated_records, all_configs):
        """All D alleles in an ambiguous call must produce the same score."""
        for config_name in CONFIGS:
            if not all_configs[config_name]["has_d"]:
                continue
            d_alleles = all_configs[config_name]["d"]
            for i, rec in enumerate(mutated_records[config_name]):
                d_call = rec["d_call"]
                if not d_call or "," not in d_call or d_call == "Short-D":
                    continue
                calls = d_call.split(",")
                seq = rec["sequence"]
                d_start = rec["d_sequence_start"]
                d_end = rec["d_sequence_end"]
                g_start = rec["d_germline_start"]
                g_end = rec["d_germline_end"]

                scores = []
                for name in calls:
                    name = name.strip()
                    if name not in d_alleles:
                        continue
                    s = _score_allele(seq, d_alleles[name],
                                      d_start, d_end, g_start, g_end)
                    if s is not None:
                        scores.append((name, s))

                if len(scores) > 1:
                    first_score = scores[0][1]
                    for name, score in scores[1:]:
                        assert score == first_score, (
                            f"[{config_name} seq {i}] Ambiguous D call: "
                            f"{scores[0][0]}={first_score} but "
                            f"{name}={score}")


# ====================================================================
# 11. V/D/J Allele Call God-Aligner
# ====================================================================

class TestAlleleCallGodAligner:
    """Every called allele must have the best possible match score."""

    def test_v_call_is_best_match(self, mutated_records, all_configs):
        """Called V alleles must have the god-aligner best score."""
        for config_name in CONFIGS:
            v_alleles = all_configs[config_name]["v"]
            n_checked = 0
            for i, rec in enumerate(mutated_records[config_name]):
                seq = rec["sequence"]
                v_start = rec["v_sequence_start"]
                v_end = rec["v_sequence_end"]
                g_start = rec["v_germline_start"]
                g_end = rec["v_germline_end"]
                if v_end <= v_start:
                    continue
                scores = {}
                for name, allele_seq in v_alleles.items():
                    s = _score_allele(seq, allele_seq, v_start, v_end,
                                      g_start, g_end)
                    if s is not None:
                        scores[name] = s
                if not scores:
                    continue
                best = max(scores.values())
                actual = set(rec["v_call"].split(","))
                n_checked += 1
                for called in actual:
                    if called in scores:
                        assert scores[called] == best, (
                            f"[{config_name} seq {i}] V call {called} "
                            f"scores {scores[called]} but best is {best}")
            assert n_checked > 0, f"No V calls checked for {config_name}"

    def test_j_call_is_best_match(self, mutated_records, all_configs):
        """Called J alleles must have the god-aligner best score."""
        for config_name in CONFIGS:
            j_alleles = all_configs[config_name]["j"]
            n_checked = 0
            for i, rec in enumerate(mutated_records[config_name]):
                seq = rec["sequence"]
                j_start = rec["j_sequence_start"]
                j_end = rec["j_sequence_end"]
                g_start = rec["j_germline_start"]
                g_end = rec["j_germline_end"]
                if j_end <= j_start:
                    continue
                scores = {}
                for name, allele_seq in j_alleles.items():
                    s = _score_allele(seq, allele_seq, j_start, j_end,
                                      g_start, g_end)
                    if s is not None:
                        scores[name] = s
                if not scores:
                    continue
                best = max(scores.values())
                actual = set(rec["j_call"].split(","))
                n_checked += 1
                for called in actual:
                    if called in scores:
                        assert scores[called] == best, (
                            f"[{config_name} seq {i}] J call {called} "
                            f"scores {scores[called]} but best is {best}")
            assert n_checked > 0, f"No J calls checked for {config_name}"

    def test_d_call_is_best_match_or_short_d(self, mutated_records, all_configs):
        """D call should be Short-D or match god-aligner best score."""
        for config_name in CONFIGS:
            if not all_configs[config_name]["has_d"]:
                continue
            d_alleles = all_configs[config_name]["d"]
            n_checked = 0
            for i, rec in enumerate(mutated_records[config_name]):
                d_call = rec["d_call"]
                if d_call == "Short-D" or not d_call:
                    continue
                if rec["d_sequence_start"] == rec["d_sequence_end"]:
                    continue
                seq = rec["sequence"]
                d_start = rec["d_sequence_start"]
                d_end = rec["d_sequence_end"]
                g_start = rec["d_germline_start"]
                g_end = rec["d_germline_end"]
                scores = {}
                for name, allele_seq in d_alleles.items():
                    s = _score_allele(seq, allele_seq, d_start, d_end,
                                      g_start, g_end)
                    if s is not None:
                        scores[name] = s
                if not scores:
                    continue
                best = max(scores.values())
                actual = set(d_call.split(","))
                n_checked += 1
                for called in actual:
                    if called in scores:
                        assert scores[called] == best, (
                            f"[{config_name} seq {i}] D call {called} "
                            f"scores {scores[called]} but best is {best}")
            assert n_checked > 0, f"No D calls checked for {config_name}"


# ====================================================================
# 12. Noisy Sequences -- Metadata Survives Corruption
# ====================================================================

class TestNoisySequenceIntegrity:
    """With 5' corruption + indels + Ns, basic invariants still hold."""

    def test_positions_valid(self, noisy_records):
        """Segment positions remain in bounds after corruption."""
        for config_name, records in noisy_records.items():
            for i, rec in enumerate(records):
                slen = len(rec["sequence"])
                for field in ("v_sequence_start", "v_sequence_end",
                              "j_sequence_start", "j_sequence_end",
                              "junction_start", "junction_end"):
                    assert 0 <= rec[field] <= slen, (
                        f"[{config_name} seq {i}] {field}={rec[field]} "
                        f"out of [0, {slen}]")

    def test_germline_alignment_length(self, noisy_records):
        """germline_alignment must match sequence length."""
        for config_name, records in noisy_records.items():
            for i, rec in enumerate(records):
                assert len(rec["germline_alignment"]) == len(rec["sequence"]), (
                    f"[{config_name} seq {i}] germline length "
                    f"{len(rec['germline_alignment'])} != "
                    f"seq length {len(rec['sequence'])}")

    def test_mutation_positions_valid(self, noisy_records):
        """Mutation positions must be within sequence bounds."""
        for config_name, records in noisy_records.items():
            for i, rec in enumerate(records):
                mutations = _parse_mutations(rec.get("mutations", ""))
                for pos in mutations:
                    assert 0 <= pos < len(rec["sequence"]), (
                        f"[{config_name} seq {i}] mutation at {pos} "
                        f"out of bounds (len={len(rec['sequence'])})")

    def test_v_call_nonempty(self, noisy_records):
        """v_call should always be non-empty for non-contaminants."""
        for config_name, records in noisy_records.items():
            for i, rec in enumerate(records):
                if not rec.get("is_contaminant"):
                    assert rec["v_call"], (
                        f"[{config_name} seq {i}] v_call empty")

    def test_junction_nt_matches(self, noisy_records):
        """junction_nt should still match sequence slice after corruption."""
        for config_name, records in noisy_records.items():
            for i, rec in enumerate(records):
                jnt = rec.get("junction_nt", "")
                if not jnt:
                    continue
                js = rec["junction_start"]
                je = rec["junction_end"]
                extracted = rec["sequence"][js:je]
                assert extracted == jnt, (
                    f"[{config_name} seq {i}] junction_nt mismatch after "
                    f"corruption")


# ====================================================================
# 13. VJ Chains (No D Segment)
# ====================================================================

class TestVJChainConsistency:
    """Light chains (IGK) should have no D segment."""

    def test_no_d_segment(self, mutated_records):
        """d_call empty, d positions zero for VJ chains."""
        for config_name in CONFIGS:
            if _has_d_segment(config_name):
                continue
            for i, rec in enumerate(mutated_records[config_name]):
                assert rec["d_call"] == "", (
                    f"[{config_name} seq {i}] d_call='{rec['d_call']}'")
                assert rec["d_sequence_start"] == 0
                assert rec["d_sequence_end"] == 0

    def test_j_follows_v_via_np1(self, mutated_records):
        """For VJ chains, gap between V and J equals NP1 length."""
        for config_name in CONFIGS:
            if _has_d_segment(config_name):
                continue
            for i, rec in enumerate(mutated_records[config_name]):
                gap = rec["j_sequence_start"] - rec["v_sequence_end"]
                assert gap >= 0, (
                    f"[{config_name} seq {i}] J starts before V ends")
                assert gap == rec["np1_length"], (
                    f"[{config_name} seq {i}] gap={gap} != "
                    f"np1_length={rec['np1_length']}")


# ====================================================================
# 14. Mutation Count Positional Oracle
# ====================================================================

class TestMutationCountOracle:
    """Count mutations by scanning seq vs germline, compare to stored count."""

    def test_scan_count_matches_stored(self, mutated_records):
        """Walk V/D/J regions, count seq!=germline positions."""
        for config_name in CONFIGS:
            for i, rec in enumerate(mutated_records[config_name]):
                seq = rec["sequence"].upper()
                germ = rec["germline_alignment"].upper()
                if len(seq) != len(germ):
                    continue
                regions = [
                    (rec["v_sequence_start"], rec["v_sequence_end"]),
                    (rec["j_sequence_start"], rec["j_sequence_end"]),
                ]
                d_s, d_e = rec["d_sequence_start"], rec["d_sequence_end"]
                if d_s != d_e:
                    regions.append((d_s, d_e))
                n_mut = 0
                for start, end in regions:
                    for pos in range(start, min(end, len(seq), len(germ))):
                        s, g = seq[pos], germ[pos]
                        if s == 'N' or g == 'N':
                            continue
                        if s != g:
                            n_mut += 1

                # Account for same-base mutations (SHM targeted but same base)
                n_same = 0
                mutations = _parse_mutations(rec.get("mutations", ""))
                for pos, (fb, tb) in mutations.items():
                    if fb.upper() == tb.upper():
                        n_same += 1
                expected = rec["n_mutations"] - n_same
                assert n_mut == expected, (
                    f"[{config_name} seq {i}] scan={n_mut}, "
                    f"stored={rec['n_mutations']}, same_base={n_same}")


# ====================================================================
# 15. Full Oracle -- All Checks per Sequence
# ====================================================================

class TestFullOracle:
    """The ultimate test: re-derive ALL metadata for every sequence
    and report all discrepancies at once."""

    @pytest.fixture(scope="class")
    def oracle_data(self, all_configs):  # noqa: ARG001
        """200 sequences per config with moderate SHM."""
        from GenAIRR.ops import rate
        results = {}
        for name in CONFIGS:
            results[name] = (Experiment.on(name)
                             .mutate(rate(0.03, 0.10))
                             .run(n=200, seed=42))
        return results

    def test_full_oracle(self, oracle_data, all_configs):
        """Re-derive all metadata for every sequence across all configs."""
        for config_name in CONFIGS:
            v_alleles = all_configs[config_name]["v"]
            j_alleles = all_configs[config_name]["j"]
            has_d = all_configs[config_name]["has_d"]

            for i, rec in enumerate(oracle_data[config_name]):
                seq = rec["sequence"]
                germ = rec["germline_alignment"]
                errors = []

                # -- sequence_length --
                if rec["sequence_length"] != len(seq):
                    errors.append(
                        f"sequence_length: {rec['sequence_length']} != {len(seq)}")

                # -- germline_alignment length --
                if len(germ) != len(seq):
                    errors.append(
                        f"germline_alignment length: {len(germ)} != {len(seq)}")

                # -- junction_nt --
                jnt = rec.get("junction_nt", "")
                js, je = rec["junction_start"], rec["junction_end"]
                if jnt and seq[js:je] != jnt:
                    errors.append(
                        f"junction_nt != seq[{js}:{je}]")

                # -- junction_length --
                if rec["junction_length"] != je - js:
                    errors.append(
                        f"junction_length: {rec['junction_length']} != {je - js}")

                # -- junction_aa translation --
                jaa = rec.get("junction_aa", "")
                if jnt:
                    expected_aa = _translate_sequence(jnt)
                    if jaa != expected_aa:
                        errors.append(
                            f"junction_aa: '{jaa}' != computed '{expected_aa}'")

                # -- stop_codon --
                has_stop = _has_stop_codon_in_reading_frame(seq.upper())
                if rec["stop_codon"] != has_stop:
                    errors.append(
                        f"stop_codon: stored={rec['stop_codon']}, "
                        f"scanned={has_stop}")

                # -- productivity consistency --
                if rec["stop_codon"] and rec["productive"]:
                    errors.append("productive=True but stop_codon=True")
                if rec["productive"] and not rec["vj_in_frame"]:
                    errors.append("productive=True but vj_in_frame=False")

                # -- vj_in_frame geometry --
                jlen = je - js
                geom = (js % 3 == 0) and (je % 3 == 0) and (jlen % 3 == 0)
                if rec["vj_in_frame"] and not geom:
                    errors.append(
                        f"vj_in_frame=True but junction not in frame "
                        f"[{js}:{je}] len={jlen}")

                # -- NP regions --
                v_end = rec["v_sequence_end"]
                if has_d:
                    np1_end = rec["d_sequence_start"]
                else:
                    np1_end = rec["j_sequence_start"]
                np1_extracted = seq[v_end:np1_end]
                if np1_extracted.upper() != rec["np1_region"].upper():
                    errors.append(
                        f"np1_region: '{rec['np1_region']}' != "
                        f"seq[{v_end}:{np1_end}]='{np1_extracted}'")
                if rec["np1_length"] != len(rec["np1_region"]):
                    errors.append(
                        f"np1_length: {rec['np1_length']} != "
                        f"{len(rec['np1_region'])}")

                if has_d:
                    d_end = rec["d_sequence_end"]
                    j_start = rec["j_sequence_start"]
                    np2_extracted = seq[d_end:j_start]
                    if np2_extracted.upper() != rec["np2_region"].upper():
                        errors.append(
                            f"np2_region: '{rec['np2_region']}' != "
                            f"seq[{d_end}:{j_start}]='{np2_extracted}'")

                # -- V allele call oracle --
                v_start = rec["v_sequence_start"]
                v_end_pos = rec["v_sequence_end"]
                g_start = rec["v_germline_start"]
                g_end = rec["v_germline_end"]
                if v_end_pos > v_start and g_end > g_start:
                    v_scores = {}
                    for name, aseq in v_alleles.items():
                        s = _score_allele(seq, aseq, v_start, v_end_pos,
                                          g_start, g_end)
                        if s is not None:
                            v_scores[name] = s
                    if v_scores:
                        best = max(v_scores.values())
                        for called in rec["v_call"].split(","):
                            if called in v_scores and v_scores[called] != best:
                                errors.append(
                                    f"v_call: {called} scores "
                                    f"{v_scores[called]} but best is {best}")

                # -- J allele call oracle --
                j_start = rec["j_sequence_start"]
                j_end = rec["j_sequence_end"]
                jg_start = rec["j_germline_start"]
                jg_end = rec["j_germline_end"]
                if j_end > j_start and jg_end > jg_start:
                    j_scores = {}
                    for name, aseq in j_alleles.items():
                        s = _score_allele(seq, aseq, j_start, j_end,
                                          jg_start, jg_end)
                        if s is not None:
                            j_scores[name] = s
                    if j_scores:
                        best = max(j_scores.values())
                        for called in rec["j_call"].split(","):
                            if called in j_scores and j_scores[called] != best:
                                errors.append(
                                    f"j_call: {called} scores "
                                    f"{j_scores[called]} but best is {best}")

                # -- Report --
                assert not errors, (
                    f"[{config_name} seq {i}] oracle found "
                    f"{len(errors)} discrepancies:\n" +
                    "\n".join(f"  - {e}" for e in errors))
