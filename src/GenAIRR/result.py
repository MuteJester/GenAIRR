"""
SimulationResult: wrapper around a list of simulation result dicts.

Implements ``collections.abc.Sequence`` so it is backward-compatible
with code that expects a plain list, and ``pd.DataFrame(result)``
works unchanged. Adds convenience export methods.
"""

from __future__ import annotations

import csv
from collections.abc import Sequence
from typing import Any, Dict, Iterator, List, Optional, Union


class SimulationResult(Sequence):
    """
    Container for batch simulation results.

    Wraps a ``list[dict]`` returned by ``SimulationGraph.simulate()``
    and adds export methods. Fully backward-compatible — ``pd.DataFrame(result)``
    still works, as does indexing, slicing, and iteration.

    Example::

        result = graph.simulate(1000)
        df = result.to_dataframe()
        result.to_csv("output.csv")
        result.to_fasta("output.fasta")
    """

    __slots__ = ("_records",)

    def __init__(self, records: List[Dict[str, Any]]) -> None:
        self._records = records

    # --- Sequence protocol ---

    def __getitem__(self, index: Union[int, slice]) -> Any:
        if isinstance(index, slice):
            return SimulationResult(self._records[index])
        return self._records[index]

    def __len__(self) -> int:
        return len(self._records)

    def __iter__(self) -> Iterator[Dict[str, Any]]:
        return iter(self._records)

    def __bool__(self) -> bool:
        return bool(self._records)

    # --- Export methods ---

    def to_dataframe(self) -> Any:
        """Convert to a ``pandas.DataFrame``.

        Raises:
            ImportError: If pandas is not installed.
        """
        try:
            import pandas as pd
        except ImportError:
            raise ImportError(
                "pandas is required for to_dataframe(). "
                "Install it with: pip install pandas"
            )
        return pd.DataFrame(self._records)

    def to_csv(self, path: str, delimiter: str = ",") -> None:
        """Write results to a CSV file using stdlib csv.

        Always creates a file. If there are no records, writes a
        header-only CSV.

        Args:
            path: Output file path.
            delimiter: Field delimiter (default comma).
        """
        if self._records:
            fieldnames = list(self._records[0].keys())
        else:
            fieldnames = ["sequence"]
        with open(path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter=delimiter)
            writer.writeheader()
            writer.writerows(self._records)

    def to_fasta(self, path: str) -> None:
        """Write sequences to a FASTA file with descriptive headers.

        Each header includes v_call, j_call, and mutation_rate.

        Args:
            path: Output file path.
        """
        with open(path, "w") as f:
            for i, rec in enumerate(self._records):
                v = rec.get("v_call", "")
                j = rec.get("j_call", "")
                mr = rec.get("mutation_rate", 0)
                header = f">seq_{i} v_call={v} j_call={j} mutation_rate={mr:.4f}"
                f.write(header + "\n")
                seq = rec.get("sequence", "")
                # Wrap at 80 characters
                for start in range(0, len(seq), 80):
                    f.write(seq[start:start + 80] + "\n")

    # --- repr ---

    def __repr__(self) -> str:
        n = len(self._records)
        fields = len(self._records[0]) if self._records else 0
        return f"SimulationResult({n} sequences, {fields} fields)"


# ── Live execution trace narration ──────────────────────────────────


def narrate(experiment, *, seed: int = 42, color: bool = True) -> str:
    """
    Simulate one sequence with full execution tracing and return a
    human-readable narrative of every internal decision.

    Unlike ``narrate_from_record()`` which reconstructs the story from
    the final AIRR fields, this function hooks into the C simulation engine
    and records every step as it happens: allele selection, trimming amounts,
    NP generation, assembly, productivity assessment, mutation events
    (position, context, from→to), corruption, boundary ambiguity resolution,
    and AIRR serialization.

    Args:
        experiment: An ``Experiment`` instance (before ``.run()``).
        seed: Random seed for reproducibility (default 42).
        color: Use ANSI color codes (default True; set False for plain text).

    Returns:
        A multi-line string with the formatted trace narrative.

    Example::

        from GenAIRR import Experiment, narrate
        exp = Experiment.on("human_igh").somatic_hypermutation(min_rate=0.02, max_rate=0.08).using_s5f()
        print(narrate(exp, seed=42))
    """
    # Compile a one-shot simulator from the Experiment
    sim = experiment.compile(seed=seed)
    csim = sim._sim

    # Enable tracing, run one sequence, read trace
    csim.set_trace(True)
    record = csim.simulate_one()
    trace_text = csim.get_trace()
    csim.set_trace(False)

    # ANSI codes
    if color:
        DIM = "\033[2m"
        BOLD = "\033[1m"
        RST = "\033[0m"
        GRN = "\033[32m"
        RED = "\033[31m"
        YEL = "\033[33m"
        CYN = "\033[36m"
        MAG = "\033[35m"
        BLU = "\033[34m"
    else:
        DIM = BOLD = RST = GRN = RED = YEL = CYN = MAG = BLU = ""

    lines: List[str] = []
    lines.append(f"{BOLD}{CYN}{'═' * 70}{RST}")
    lines.append(f"{BOLD}{CYN}  GENAIRR EXECUTION TRACE  —  Live simulation narrative{RST}")
    lines.append(f"{BOLD}{CYN}{'═' * 70}{RST}")

    # Group trace lines by phase
    current_phase = None
    phase_map = {
        "sample_v": "Allele Selection",
        "sample_d": "Allele Selection",
        "sample_j": "Allele Selection",
        "sample_c": "Allele Selection",
        "trim_v": "Exonuclease Trimming",
        "trim_d": "Exonuclease Trimming",
        "trim_j": "Exonuclease Trimming",
        "assemble": "Sequence Assembly",
        "assess": "Productivity Assessment",
        "s5f": "Somatic Hypermutation (S5F)",
        "uniform": "Somatic Hypermutation (Uniform)",
        "csr": "Class Switch Recombination",
        "selection": "Antigen Selection Pressure",
        "d_inversion": "D Gene Inversion",
        "receptor_revision": "Receptor Revision",
        "corrupt_5'": "5' End Corruption",
        "corrupt_3'": "3' End Corruption",
        "quality_errors": "Sequencing Quality Errors",
        "pcr": "PCR Amplification Errors",
        "indels": "Indel Introduction",
        "insert_ns": "Ambiguous N Insertion",
        "contaminant": "Contaminant Spike-in",
        "umi": "UMI Barcode",
        "primer_mask": "Primer Masking (FR1)",
        "reverse_complement": "Reverse Complement",
        "paired_end": "Paired-End Merge Profile",
        "long_read": "Long-Read Homopolymer Errors",
        "trim_to_length": "Sequence Length Enforcement",
        "airr": "AIRR Serialization",
        "pipeline": "Pipeline Control",
        "begin": "Simulation Start",
        "end": "Simulation End",
        "serialize": "Output Serialization",
    }

    for trace_line in trace_text.strip().split('\n'):
        if not trace_line:
            continue

        # Extract tag: [tag] rest
        tag = ""
        msg = trace_line
        if trace_line.startswith('['):
            bracket_end = trace_line.find(']')
            if bracket_end > 0:
                tag = trace_line[1:bracket_end]
                msg = trace_line[bracket_end + 2:] if bracket_end + 2 < len(trace_line) else ""

        # Determine phase from tag
        phase_name = phase_map.get(tag, None)
        if phase_name and phase_name != current_phase:
            current_phase = phase_name
            lines.append(f"\n{BOLD}{CYN}{'─' * 60}{RST}")
            lines.append(f"{BOLD}{CYN}  {phase_name}{RST}")
            lines.append(f"{BOLD}{CYN}{'─' * 60}{RST}")

        # Format the message with color based on content
        if "REJECTED" in msg:
            lines.append(f"  {RED}✗{RST} {DIM}{msg}{RST}")
        elif "mutate pos=" in msg:
            lines.append(f"  {MAG}↻{RST} {msg}")
        elif "allele" in tag or "sample" in tag:
            lines.append(f"  {GRN}▸{RST} {msg}")
        elif "trim" in tag:
            lines.append(f"  {YEL}✂{RST} {msg}")
        elif "ambiguity" in msg.lower():
            lines.append(f"  {BLU}⇌{RST} {msg}")
        elif "productive" in msg.lower() or "stop_codon" in msg.lower():
            if "yes" in msg:
                lines.append(f"  {GRN}✓{RST} {msg}")
            else:
                lines.append(f"  {RED}✗{RST} {msg}")
        elif "contaminant" in msg.lower():
            lines.append(f"  {RED}⚠{RST} {msg}")
        elif "done:" in msg or "total:" in msg:
            lines.append(f"  {BOLD}{msg}{RST}")
        else:
            lines.append(f"  {DIM}▪{RST} {msg}")

    # Summary footer
    lines.append(f"\n{BOLD}{'═' * 70}{RST}")
    prod = record.get('productive', False)
    n_mut = record.get('n_mutations', 0)
    seq_len = record.get('sequence_length', 0)
    verdict = f"{GRN}PRODUCTIVE{RST}" if prod else f"{RED}NON-PRODUCTIVE{RST}"
    lines.append(f"{BOLD}  Result: {seq_len}bp, {n_mut} mutations, {verdict}{RST}")
    lines.append(f"  V: {record.get('v_call', '')}, "
                 f"D: {record.get('d_call', '')}, "
                 f"J: {record.get('j_call', '')}")
    lines.append(f"{BOLD}{'═' * 70}{RST}")

    return '\n'.join(lines)


def narrate_from_record(record: Dict[str, Any], *, color: bool = True) -> str:
    """
    Generate a human-readable narrative of how a sequence was simulated.

    Reconstructs the biological "life story" from the AIRR record fields:
    allele selection, trimming, N/P addition, assembly, productivity
    assessment, somatic hypermutation, and sequencing artifacts.

    Args:
        record: An AIRR dict (from ``stream().get()`` or ``simulate()``).
        color:  Use ANSI color codes (default True; set False for plain text).

    Returns:
        A multi-line string with the formatted narrative.

    Example::

        sim = Experiment.on("human_igh").mutate(rate(0.02, 0.08)).compile(seed=42)
        seq = sim.stream().get()
        print(narrate(seq))
    """
    r = record

    # ANSI codes (disabled if color=False)
    if color:
        DIM = "\033[2m"
        BOLD = "\033[1m"
        RST = "\033[0m"
        GRN = "\033[32m"
        RED = "\033[31m"
        YEL = "\033[33m"
        CYN = "\033[36m"
        MAG = "\033[35m"
    else:
        DIM = BOLD = RST = GRN = RED = YEL = CYN = MAG = ""

    lines: List[str] = []

    def phase(title: str) -> None:
        lines.append(f"\n{BOLD}{CYN}{'─' * 60}{RST}")
        lines.append(f"{BOLD}{CYN}  {title}{RST}")
        lines.append(f"{BOLD}{CYN}{'─' * 60}{RST}")

    def step(icon: str, msg: str) -> None:
        lines.append(f"  {icon} {msg}")

    # ── Phase 1: V(D)J Recombination ────────────────────────────
    phase("Phase 1: V(D)J Recombination")
    step("\u2022", f"Selected {BOLD}V{RST} allele: {GRN}{r['v_call']}{RST}")
    step("\u2022", f"Selected {BOLD}D{RST} allele: {GRN}{r['d_call']}{RST}")
    step("\u2022", f"Selected {BOLD}J{RST} allele: {GRN}{r['j_call']}{RST}")
    if r.get('c_call'):
        step("\u2022", f"Selected {BOLD}C{RST} allele: {GRN}{r['c_call']}{RST} (isotype)")

    # ── Phase 2: Exonuclease Trimming ───────────────────────────
    phase("Phase 2: Exonuclease Trimming")
    vt5, vt3 = r['v_trim_5'], r['v_trim_3']
    dt5, dt3 = r['d_trim_5'], r['d_trim_3']
    jt5, jt3 = r['j_trim_5'], r['j_trim_3']

    if vt5 == 0 and vt3 == 0:
        step(" ", f"V gene: {DIM}no trimming{RST}")
    else:
        if vt5 > 0:
            step(" ", f"V gene: trimmed {BOLD}{vt5}bp{RST} from 5\u2032 end")
        if vt3 > 0:
            step(" ", f"V gene: trimmed {BOLD}{vt3}bp{RST} from 3\u2032 end")

    if dt5 > 0:
        step(" ", f"D gene: trimmed {BOLD}{dt5}bp{RST} from 5\u2032 end")
    if dt3 > 0:
        step(" ", f"D gene: trimmed {BOLD}{dt3}bp{RST} from 3\u2032 end")
    d_remaining = r['d_germline_end'] - r['d_germline_start']
    step(" ", f"{DIM}D gene: {d_remaining}bp remaining after trimming{RST}")

    if jt5 > 0:
        step(" ", f"J gene: trimmed {BOLD}{jt5}bp{RST} from 5\u2032 end")
    if jt3 > 0:
        step(" ", f"J gene: trimmed {BOLD}{jt3}bp{RST} from 3\u2032 end")

    # ── Phase 3: TdT N/P-Nucleotide Addition ────────────────────
    phase("Phase 3: TdT N/P-Nucleotide Addition")
    np1 = r['np1_region']
    np2 = r['np2_region']
    if np1:
        step(" ", f"NP1 (V\u2192D junction): inserted {BOLD}{len(np1)}bp{RST}: {YEL}{np1}{RST}")
    else:
        step(" ", f"NP1 (V\u2192D junction): {DIM}no nucleotides added{RST}")
    if np2:
        step(" ", f"NP2 (D\u2192J junction): inserted {BOLD}{len(np2)}bp{RST}: {YEL}{np2}{RST}")
    else:
        step(" ", f"NP2 (D\u2192J junction): {DIM}no nucleotides added{RST}")

    # ── Phase 4: Sequence Assembly ──────────────────────────────
    phase("Phase 4: Sequence Assembly")
    vs, ve = r['v_sequence_start'], r['v_sequence_end']
    ds, de = r['d_sequence_start'], r['d_sequence_end']
    js, je = r['j_sequence_start'], r['j_sequence_end']
    step(" ", f"Assembled: V({ve - vs}bp) + NP1({r['np1_length']}bp) "
              f"+ D({de - ds}bp) + NP2({r['np2_length']}bp) + J({je - js}bp)")
    step(" ", f"Total assembled length: {BOLD}{r['sequence_length']}bp{RST}")

    jnc_s, jnc_e = r['junction_start'], r['junction_end']
    jnc_len = jnc_e - jnc_s
    step(" ", f"Junction (CDR3): positions [{jnc_s}:{jnc_e}] = {jnc_len}bp")
    step(" ", f"  nt: {r['junction_nt']}")
    step(" ", f"  aa: {r['junction_aa']}")

    # ── Phase 5: Productivity Assessment ────────────────────────
    phase("Phase 5: Productivity Assessment")
    in_frame = r['vj_in_frame']
    has_stop = r['stop_codon']
    productive = r['productive']

    if in_frame:
        step(" ", f"VJ reading frame: {GRN}IN FRAME{RST} "
                  f"(junction {jnc_len}bp, divisible by 3)")
    else:
        step(" ", f"VJ reading frame: {RED}OUT OF FRAME{RST} "
                  f"(junction {jnc_len}bp mod 3 = {jnc_len % 3})")

    jaa = r.get('junction_aa', '')
    if jaa:
        first_aa, last_aa = jaa[0], jaa[-1]
        cys_ok = first_aa == 'C'
        anc_ok = last_aa in ('W', 'F')
        step(" ", f"Conserved Cys (V anchor): "
                  f"{GRN + 'PRESENT' + RST if cys_ok else RED + 'MISSING' + RST} "
                  f"('{first_aa}')")
        step(" ", f"Conserved W/F (J anchor): "
                  f"{GRN + 'PRESENT' + RST if anc_ok else RED + 'MISSING' + RST} "
                  f"('{last_aa}')")

    if has_stop:
        stop_idx = jaa.find('*') if jaa else -1
        if stop_idx >= 0:
            step(" ", f"Stop codon: {RED}FOUND{RST} "
                      f"(junction AA position {stop_idx})")
        else:
            step(" ", f"Stop codon: {RED}FOUND{RST}")
    else:
        step(" ", f"Stop codon: {GRN}NONE{RST}")

    if productive:
        step(" ", f"Verdict: {BOLD}{GRN}PRODUCTIVE{RST}")
    else:
        note = r.get('note', '')
        step(" ", f"Verdict: {BOLD}{RED}NON-PRODUCTIVE{RST}"
                  + (f" \u2014 {note}" if note else ""))

    # ── Phase 6: Somatic Hypermutation ──────────────────────────
    phase("Phase 6: Somatic Hypermutation")
    n_mut = r['n_mutations']
    mut_rate = r['mutation_rate']

    if n_mut == 0:
        step(" ", f"{DIM}No somatic mutations{RST}")
    else:
        step(" ", f"Applied {BOLD}{n_mut} mutations{RST} "
                  f"(rate = {mut_rate:.4f} = {mut_rate * 100:.2f}%)")

        purines = set('AGag')
        pyrimidines = set('CTct')
        v_n = d_n = j_n = 0
        ti = tv = 0
        mut_lines: List[str] = []

        muts = [m for m in r['mutations'].split(',') if m]
        for m in muts:
            pos_s, change = m.split(':')
            pos = int(pos_s)
            orig, new = change[0], change[2]

            if pos < ve:
                region = "V"
                v_n += 1
            elif pos < ds:
                region = "NP1"
            elif pos < de:
                region = "D"
                d_n += 1
            elif pos < js:
                region = "NP2"
            else:
                region = "J"
                j_n += 1

            if (orig in purines and new in purines) or \
               (orig in pyrimidines and new in pyrimidines):
                ti += 1
            else:
                tv += 1

            mut_lines.append(
                f"    pos {pos:3d} ({region:3s}): {orig} \u2192 {MAG}{new}{RST}"
            )

        step(" ", f"  V region: {v_n} mutations")
        if de > ds:
            step(" ", f"  D region: {d_n} mutations")
        step(" ", f"  J region: {j_n} mutations")
        if tv > 0:
            step(" ", f"  Transitions: {ti}, Transversions: {tv} "
                      f"(Ti/Tv = {ti / tv:.1f})")
        else:
            step(" ", f"  All transitions ({ti})")

        for ml in mut_lines:
            lines.append(ml)

    # ── Phase 7: Sequencing Artifacts ───────────────────────────
    n_seq_err = r.get('n_sequencing_errors', 0)
    n_pcr_err = r.get('n_pcr_errors', 0)
    is_rc = r.get('is_reverse_complement', False)
    is_contam = r.get('is_contaminant', False)

    phase("Phase 7: Sequencing Artifacts")
    if n_pcr_err > 0:
        step(" ", f"PCR errors: {BOLD}{n_pcr_err}{RST}")
        pcr_str = r.get('pcr_errors', '')
        if pcr_str:
            for m in pcr_str.split(','):
                if m:
                    lines.append(f"    {m}")
    if n_seq_err > 0:
        step(" ", f"Sequencing errors: {BOLD}{n_seq_err}{RST}")
        seq_str = r.get('sequencing_errors', '')
        if seq_str:
            for m in seq_str.split(','):
                if m:
                    lines.append(f"    {m}")
    if is_rc:
        step(" ", f"Reverse complemented: {YEL}YES{RST}")
    if is_contam:
        step(" ", f"Contaminant: {RED}YES{RST}")
    if not (n_pcr_err or n_seq_err or is_rc or is_contam):
        step(" ", f"{DIM}No sequencing artifacts{RST}")

    # ── Summary ─────────────────────────────────────────────────
    lines.append(f"\n{BOLD}{'=' * 60}{RST}")
    verdict = f"{GRN}PRODUCTIVE{RST}" if productive else f"{RED}NON-PRODUCTIVE{RST}"
    lines.append(f"{BOLD}  Final: {r['sequence_length']}bp, "
                 f"{n_mut} mutations, {verdict}{RST}")
    lines.append(f"{BOLD}{'=' * 60}{RST}")

    return '\n'.join(lines)
