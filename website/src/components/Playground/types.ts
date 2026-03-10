/* ───────────────────────────────────────────────────────────
 *  types.ts — Typed AIRR result interface and parser.
 *  Maps raw Pyodide JSON output to the SequenceRecord interface
 *  used by DissectionBay, ResultsTable, etc.
 * ─────────────────────────────────────────────────────────── */

// ── Segment types ─────────────────────────────────────────

export type SegmentId = 'V' | 'NP1' | 'D' | 'NP2' | 'J';

export interface SegmentMeta {
  id: SegmentId;
  label: string;
  fullLabel: string;
  color: string;
  start: number;
  end: number;
}

export const SEGMENT_COLORS: Record<SegmentId, string> = {
  V:   '#4A90D9',
  NP1: '#5DC48C',
  D:   '#E8685A',
  NP2: '#5DC48C',
  J:   '#E8A838',
};

// ── Sequence record ───────────────────────────────────────

export interface SequenceRecord {
  id: number;
  sequence_id: string;
  sequence: string;
  v_call: string;
  d_call: string;
  j_call: string;
  junction_aa: string;
  productive: boolean;
  mutation_count: number;
  cdr3_length: number;
  mutation_rate: number;
  mutations: Record<number, string>;
  // Region positions
  v_sequence_start: number;
  v_sequence_end: number;
  d_sequence_start: number;
  d_sequence_end: number;
  j_sequence_start: number;
  j_sequence_end: number;
  junction_start: number;
  junction_end: number;
  // NP provenance
  np1_region: string;
  np2_region: string;
  np1_p_prefix: string;
  np1_n_region: string;
  np1_p_suffix: string;
  np2_p_prefix: string;
  np2_n_region: string;
  np2_p_suffix: string;
  // Artifacts
  corruption_5prime: number;
  corruption_3prime: number;
  d_inverted: boolean;
  v_trim_3: number;
  d_trim_5: number;
  d_trim_3: number;
  j_trim_5: number;
}

// ── Parser ────────────────────────────────────────────────

function parseMutations(raw: any): Record<number, string> {
  if (!raw) return {};
  // May come as a JSON string or already parsed dict
  if (typeof raw === 'string') {
    try {
      const parsed = JSON.parse(raw);
      const out: Record<number, string> = {};
      for (const [k, v] of Object.entries(parsed)) {
        out[Number(k)] = String(v);
      }
      return out;
    } catch {
      return {};
    }
  }
  if (typeof raw === 'object') {
    const out: Record<number, string> = {};
    for (const [k, v] of Object.entries(raw)) {
      out[Number(k)] = String(v);
    }
    return out;
  }
  return {};
}

export function parseAIRRResults(
  raw: Record<string, any>[],
): SequenceRecord[] {
  return raw.map((r, i) => {
    const mutations = parseMutations(r.mutations);
    const junctionAA = r.junction_aa ?? r.junction ?? '';
    return {
      id: i + 1,
      sequence_id: r.sequence_id ?? `SEQ_${String(i + 1).padStart(4, '0')}`,
      sequence: r.sequence ?? '',
      v_call: r.v_call ?? '',
      d_call: r.d_call ?? '',
      j_call: r.j_call ?? '',
      junction_aa: junctionAA,
      productive: Boolean(r.productive),
      mutation_count: Object.keys(mutations).length,
      cdr3_length: junctionAA.length,
      mutation_rate: Number(r.mutation_rate ?? 0),
      mutations,
      v_sequence_start: Number(r.v_sequence_start ?? 0),
      v_sequence_end: Number(r.v_sequence_end ?? 0),
      d_sequence_start: Number(r.d_sequence_start ?? 0),
      d_sequence_end: Number(r.d_sequence_end ?? 0),
      j_sequence_start: Number(r.j_sequence_start ?? 0),
      j_sequence_end: Number(r.j_sequence_end ?? 0),
      junction_start: Number(r.junction_sequence_start ?? r.junction_start ?? 0),
      junction_end: Number(r.junction_sequence_end ?? r.junction_end ?? 0),
      np1_region: r.np1_region ?? '',
      np2_region: r.np2_region ?? '',
      np1_p_prefix: r.np1_p_prefix ?? '',
      np1_n_region: r.np1_n_region ?? '',
      np1_p_suffix: r.np1_p_suffix ?? '',
      np2_p_prefix: r.np2_p_prefix ?? '',
      np2_n_region: r.np2_n_region ?? '',
      np2_p_suffix: r.np2_p_suffix ?? '',
      corruption_5prime: Number(r.corruption_5prime ?? 0),
      corruption_3prime: Number(r.corruption_3prime ?? 0),
      d_inverted: Boolean(r.d_inverted),
      v_trim_3: Number(r.v_trim_3 ?? 0),
      d_trim_5: Number(r.d_trim_5 ?? 0),
      d_trim_3: Number(r.d_trim_3 ?? 0),
      j_trim_5: Number(r.j_trim_5 ?? 0),
    };
  });
}
