/* ───────────────────────────────────────────────────────────
 *  ops.ts — Single source of truth for the 20 GenAIRR ops.
 *  Drives the catalog, pipeline validation, param editors,
 *  and Python code generation.
 * ─────────────────────────────────────────────────────────── */

// ── Types ──────────────────────────────────────────────────

export type ParamType = 'number' | 'boolean' | 'select' | 'range';
export type Category = 'rearrangement' | 'mutation' | 'biology' | 'artifacts';

export interface OpParam {
  name: string;
  label: string;
  type: ParamType;
  default: number | boolean | string;
  min?: number;
  max?: number;
  step?: number;
  options?: { label: string; value: string }[];
  description: string;
  /** Return false to hide this param based on sibling values. */
  condition?: (params: Record<string, any>) => boolean;
}

export interface OpDef {
  id: string;
  label: string;
  category: Category;
  description: string;
  pythonClass: string;
  /** Extra Python imports needed (mutation models, etc.) */
  extraImports?: string[];
  params: OpParam[];
  /** Op IDs that must already be in the pipeline. */
  requires?: string[];
  /** Dependency names this op satisfies for downstream ops. */
  satisfies?: string[];
  /** Mutually exclusive op IDs. */
  excludes?: string[];
  /** Build the Python constructor string from current params. */
  toPython: (params: Record<string, any>) => string;
}

// ── Category colours (reused in UI) ───────────────────────

export const CATEGORY_META: Record<
  Category,
  { label: string; color: string; darkColor: string }
> = {
  rearrangement: { label: 'Rearrangement', color: '#1a6fac', darkColor: '#5ba8d9' },
  mutation:      { label: 'Mutation & Selection', color: '#00b4a0', darkColor: '#00d4be' },
  biology:       { label: 'Biological Events', color: '#8b5cc4', darkColor: '#b38ee6' },
  artifacts:     { label: 'Sequencing Artifacts', color: '#e8863a', darkColor: '#f0a060' },
};

// ── Helper: format a float without trailing zeros ─────────

function f(v: number): string {
  return Number(v).toString();
}

function optStr(v: string | undefined): string {
  return v ? `"${v}"` : 'None';
}

// ── Op Definitions ────────────────────────────────────────

export const OP_DEFS: OpDef[] = [
  // ═══════════════════════════════════════════════════════
  //  REARRANGEMENT
  // ═══════════════════════════════════════════════════════
  {
    id: 'Rearrange',
    label: 'Rearrange',
    category: 'rearrangement',
    description: 'V(D)J recombination — the foundation of every protocol.',
    pythonClass: 'Rearrange',
    params: [
      { name: 'v_allele', label: 'Force V allele', type: 'select', default: '', options: [{ label: 'Random (default)', value: '' }], description: 'Leave empty for weighted random selection.' },
      { name: 'd_allele', label: 'Force D allele', type: 'select', default: '', options: [{ label: 'Random (default)', value: '' }], description: 'Leave empty for weighted random selection.' },
      { name: 'j_allele', label: 'Force J allele', type: 'select', default: '', options: [{ label: 'Random (default)', value: '' }], description: 'Leave empty for weighted random selection.' },
    ],
    satisfies: ['Rearrange'],
    excludes: ['SimulateAllelePolymorphism'],
    toPython: (p) => {
      const args: string[] = [];
      if (p.v_allele) args.push(`v_allele="${p.v_allele}"`);
      if (p.d_allele) args.push(`d_allele="${p.d_allele}"`);
      if (p.j_allele) args.push(`j_allele="${p.j_allele}"`);
      return `Rearrange(${args.join(', ')})`;
    },
  },

  {
    id: 'SimulateAllelePolymorphism',
    label: 'Allele Polymorphism',
    category: 'rearrangement',
    description: 'Generate novel germline allele variants for benchmarking.',
    pythonClass: 'SimulateAllelePolymorphism',
    params: [
      { name: 'n_novel_v', label: 'Novel V alleles', type: 'number', default: 5, min: 0, max: 20, step: 1, description: 'Number of novel V alleles to create.' },
      { name: 'n_novel_d', label: 'Novel D alleles', type: 'number', default: 1, min: 0, max: 10, step: 1, description: 'Number of novel D alleles to create.' },
      { name: 'n_novel_j', label: 'Novel J alleles', type: 'number', default: 1, min: 0, max: 10, step: 1, description: 'Number of novel J alleles to create.' },
      { name: 'min_snps', label: 'Min SNPs', type: 'number', default: 1, min: 1, max: 10, step: 1, description: 'Minimum SNPs per novel allele.' },
      { name: 'max_snps_v', label: 'Max SNPs (V)', type: 'number', default: 5, min: 1, max: 20, step: 1, description: 'Maximum SNPs for V alleles.' },
      { name: 'max_snps_d', label: 'Max SNPs (D)', type: 'number', default: 3, min: 1, max: 10, step: 1, description: 'Maximum SNPs for D alleles.' },
      { name: 'max_snps_j', label: 'Max SNPs (J)', type: 'number', default: 2, min: 1, max: 10, step: 1, description: 'Maximum SNPs for J alleles.' },
      { name: 'protect_anchors', label: 'Protect anchors', type: 'boolean', default: true, description: 'Avoid SNPs at conserved anchor positions.' },
      { name: 'protect_codons', label: 'Protect codons', type: 'boolean', default: true, description: 'Reject stop-codon-creating variants.' },
    ],
    satisfies: ['Rearrange'],
    excludes: ['Rearrange'],
    toPython: (p) =>
      `SimulateAllelePolymorphism(n_novel_v=${p.n_novel_v}, n_novel_d=${p.n_novel_d}, n_novel_j=${p.n_novel_j}, min_snps=${p.min_snps}, max_snps_v=${p.max_snps_v}, max_snps_d=${p.max_snps_d}, max_snps_j=${p.max_snps_j}, protect_anchors=${p.protect_anchors ? 'True' : 'False'}, protect_codons=${p.protect_codons ? 'True' : 'False'})`,
  },

  // ═══════════════════════════════════════════════════════
  //  MUTATION & SELECTION
  // ═══════════════════════════════════════════════════════
  {
    id: 'Mutate',
    label: 'Somatic Hypermutation',
    category: 'mutation',
    description: 'Apply SHM using S5F or Uniform mutation model.',
    pythonClass: 'Mutate',
    extraImports: ['from GenAIRR.mutation import S5F, Uniform'],
    params: [
      { name: 'model', label: 'Mutation model', type: 'select', default: 'S5F', options: [{ label: 'S5F (context-dependent)', value: 'S5F' }, { label: 'Uniform', value: 'Uniform' }], description: 'S5F uses 5-mer context; Uniform is position-independent.' },
      { name: 'min_mutation_rate', label: 'Min mutation rate', type: 'range', default: 0.01, min: 0, max: 0.3, step: 0.001, description: 'Lower bound of per-sequence mutation rate.' },
      { name: 'max_mutation_rate', label: 'Max mutation rate', type: 'range', default: 0.05, min: 0, max: 0.3, step: 0.001, description: 'Upper bound of per-sequence mutation rate.' },
    ],
    requires: ['Rearrange'],
    satisfies: ['Mutate'],
    excludes: ['SimulateCSR'],
    toPython: (p) =>
      `Mutate(${p.model}(${f(p.min_mutation_rate)}, ${f(p.max_mutation_rate)}))`,
  },

  {
    id: 'SimulateCSR',
    label: 'Class Switch Recombination',
    category: 'mutation',
    description: 'Isotype-aware SHM — mutation rates vary by C allele isotype.',
    pythonClass: 'SimulateCSR',
    extraImports: ['from GenAIRR.mutation import S5F, Uniform'],
    params: [
      { name: 'model', label: 'Mutation model', type: 'select', default: 'S5F', options: [{ label: 'S5F (context-dependent)', value: 'S5F' }, { label: 'Uniform', value: 'Uniform' }], description: 'Base mutation model (rates adjusted per isotype).' },
      { name: 'min_mutation_rate', label: 'Fallback min rate', type: 'range', default: 0.01, min: 0, max: 0.3, step: 0.001, description: 'Default min rate for unknown isotypes.' },
      { name: 'max_mutation_rate', label: 'Fallback max rate', type: 'range', default: 0.05, min: 0, max: 0.3, step: 0.001, description: 'Default max rate for unknown isotypes.' },
    ],
    requires: ['Rearrange'],
    satisfies: ['Mutate'],
    excludes: ['Mutate'],
    toPython: (p) =>
      `SimulateCSR(mutation_model=${p.model}(${f(p.min_mutation_rate)}, ${f(p.max_mutation_rate)}))`,
  },

  {
    id: 'SelectionPressure',
    label: 'Selection Pressure',
    category: 'mutation',
    description: 'Antigen-driven selection — filter replacement mutations by IMGT region.',
    pythonClass: 'SelectionPressure',
    params: [
      { name: 'strength', label: 'Strength', type: 'range', default: 0.5, min: 0, max: 1, step: 0.05, description: '0 = no selection, 1 = maximum selection pressure.' },
      { name: 'cdr_r_acceptance', label: 'CDR R acceptance', type: 'range', default: 0.85, min: 0, max: 1, step: 0.05, description: 'Probability of keeping replacement mutations in CDR regions at full strength.' },
      { name: 'fwr_r_acceptance', label: 'FWR R acceptance', type: 'range', default: 0.40, min: 0, max: 1, step: 0.05, description: 'Probability of keeping replacement mutations in FWR regions at full strength.' },
    ],
    requires: ['Mutate'],
    toPython: (p) =>
      `SelectionPressure(strength=${f(p.strength)}, cdr_r_acceptance=${f(p.cdr_r_acceptance)}, fwr_r_acceptance=${f(p.fwr_r_acceptance)})`,
  },

  // ═══════════════════════════════════════════════════════
  //  BIOLOGICAL EVENTS
  // ═══════════════════════════════════════════════════════
  {
    id: 'SimulateDGeneInversion',
    label: 'D-Gene Inversion',
    category: 'biology',
    description: 'Reverse-complement the D segment (~15% of rearrangements).',
    pythonClass: 'SimulateDGeneInversion',
    params: [
      { name: 'probability', label: 'Probability', type: 'range', default: 0.15, min: 0, max: 1, step: 0.05, description: 'Per-sequence probability of D inversion.' },
    ],
    requires: ['Rearrange'],
    toPython: (p) => `SimulateDGeneInversion(probability=${f(p.probability)})`,
  },

  {
    id: 'SimulateReceptorRevision',
    label: 'Receptor Revision',
    category: 'biology',
    description: 'V-gene replacement with a junction footprint of the original V.',
    pythonClass: 'SimulateReceptorRevision',
    params: [
      { name: 'probability', label: 'Probability', type: 'range', default: 0.05, min: 0, max: 1, step: 0.01, description: 'Per-sequence probability of receptor revision.' },
      { name: 'footprint_min', label: 'Min footprint (nt)', type: 'number', default: 5, min: 1, max: 30, step: 1, description: 'Minimum footprint length from original V.' },
      { name: 'footprint_max', label: 'Max footprint (nt)', type: 'number', default: 20, min: 1, max: 50, step: 1, description: 'Maximum footprint length from original V.' },
    ],
    requires: ['Rearrange'],
    toPython: (p) =>
      `SimulateReceptorRevision(probability=${f(p.probability)}, footprint_min=${p.footprint_min}, footprint_max=${p.footprint_max})`,
  },

  // ═══════════════════════════════════════════════════════
  //  SEQUENCING ARTIFACTS
  // ═══════════════════════════════════════════════════════
  {
    id: 'Corrupt5Prime',
    label: "5' Corruption",
    category: 'artifacts',
    description: "Add/remove nucleotides at the 5' end of the sequence.",
    pythonClass: 'Corrupt5Prime',
    params: [
      { name: 'probability', label: 'Probability', type: 'range', default: 0.7, min: 0, max: 1, step: 0.05, description: 'Per-sequence probability of corruption.' },
    ],
    requires: ['Rearrange'],
    toPython: (p) => `Corrupt5Prime(probability=${f(p.probability)})`,
  },

  {
    id: 'Corrupt3Prime',
    label: "3' Corruption",
    category: 'artifacts',
    description: "Add/remove nucleotides at the 3' end of the sequence.",
    pythonClass: 'Corrupt3Prime',
    params: [
      { name: 'probability', label: 'Probability', type: 'range', default: 0.7, min: 0, max: 1, step: 0.05, description: 'Per-sequence probability of corruption.' },
    ],
    requires: ['Rearrange'],
    toPython: (p) => `Corrupt3Prime(probability=${f(p.probability)})`,
  },

  {
    id: 'CorruptQuality',
    label: 'Quality Errors (Illumina SE)',
    category: 'artifacts',
    description: 'Position-dependent sequencing errors with linear error profile.',
    pythonClass: 'CorruptQuality',
    params: [
      { name: 'base_error_rate', label: "Base error rate (5')", type: 'range', default: 0.001, min: 0, max: 0.05, step: 0.001, description: "Error probability at the 5' end (~Q30)." },
      { name: 'peak_error_rate', label: "Peak error rate (3')", type: 'range', default: 0.02, min: 0, max: 0.1, step: 0.001, description: "Error probability at the 3' end (~Q17)." },
      { name: 'probability', label: 'Probability', type: 'range', default: 1.0, min: 0, max: 1, step: 0.05, description: 'Per-sequence probability of applying errors.' },
    ],
    requires: ['Rearrange'],
    excludes: ['SimulatePairedEnd', 'SkewBaseComposition'],
    toPython: (p) =>
      `CorruptQuality(base_error_rate=${f(p.base_error_rate)}, peak_error_rate=${f(p.peak_error_rate)}, probability=${f(p.probability)})`,
  },

  {
    id: 'SimulatePairedEnd',
    label: 'Paired-End Merge',
    category: 'artifacts',
    description: 'Paired-end merge artifacts with bathtub error profile.',
    pythonClass: 'SimulatePairedEnd',
    params: [
      { name: 'read_length', label: 'Read length', type: 'number', default: 300, min: 50, max: 600, step: 50, description: 'Cycles per read (common: 150, 250, 300).' },
      { name: 'base_error_rate', label: 'Base error rate', type: 'range', default: 0.001, min: 0, max: 0.05, step: 0.001, description: 'Error probability at cycle 1.' },
      { name: 'peak_error_rate', label: 'Peak error rate', type: 'range', default: 0.02, min: 0, max: 0.1, step: 0.001, description: 'Error probability at last cycle.' },
      { name: 'probability', label: 'Probability', type: 'range', default: 1.0, min: 0, max: 1, step: 0.05, description: 'Per-sequence probability.' },
    ],
    requires: ['Rearrange'],
    excludes: ['CorruptQuality', 'SkewBaseComposition'],
    toPython: (p) =>
      `SimulatePairedEnd(read_length=${p.read_length}, base_error_rate=${f(p.base_error_rate)}, peak_error_rate=${f(p.peak_error_rate)}, probability=${f(p.probability)})`,
  },

  {
    id: 'SkewBaseComposition',
    label: 'Long-Read Errors (Nanopore)',
    category: 'artifacts',
    description: 'Homopolymer-targeted indels for Nanopore/PacBio profiles.',
    pythonClass: 'SkewBaseComposition',
    params: [
      { name: 'error_rate', label: 'Error rate', type: 'range', default: 0.03, min: 0, max: 0.15, step: 0.005, description: 'Per-run-base error rate.' },
      { name: 'min_run_length', label: 'Min run length', type: 'number', default: 3, min: 2, max: 10, step: 1, description: 'Minimum homopolymer length to target.' },
      { name: 'insertion_bias', label: 'Insertion bias', type: 'range', default: 0.6, min: 0, max: 1, step: 0.05, description: 'Fraction of errors that are insertions.' },
      { name: 'probability', label: 'Probability', type: 'range', default: 1.0, min: 0, max: 1, step: 0.05, description: 'Per-sequence probability.' },
    ],
    requires: ['Rearrange'],
    excludes: ['CorruptQuality', 'SimulatePairedEnd'],
    toPython: (p) =>
      `SkewBaseComposition(error_rate=${f(p.error_rate)}, min_run_length=${p.min_run_length}, insertion_bias=${f(p.insertion_bias)}, probability=${f(p.probability)})`,
  },

  {
    id: 'PCRAmplification',
    label: 'PCR Amplification',
    category: 'artifacts',
    description: 'Polymerase substitution errors accumulated over PCR cycles.',
    pythonClass: 'PCRAmplification',
    params: [
      { name: 'error_rate', label: 'Error rate per cycle', type: 'range', default: 0.0001, min: 0, max: 0.001, step: 0.00001, description: 'Per-base per-cycle error rate (Taq: ~1e-4).' },
      { name: 'n_cycles', label: 'PCR cycles', type: 'number', default: 30, min: 1, max: 50, step: 1, description: 'Number of amplification cycles.' },
      { name: 'probability', label: 'Probability', type: 'range', default: 1.0, min: 0, max: 1, step: 0.05, description: 'Per-sequence probability.' },
    ],
    requires: ['Rearrange'],
    toPython: (p) =>
      `PCRAmplification(error_rate=${f(p.error_rate)}, n_cycles=${p.n_cycles}, probability=${f(p.probability)})`,
  },

  {
    id: 'SimulateUMI',
    label: 'UMI Barcode',
    category: 'artifacts',
    description: 'Prepend a random UMI barcode and shift annotations.',
    pythonClass: 'SimulateUMI',
    params: [
      { name: 'umi_length', label: 'UMI length (nt)', type: 'number', default: 12, min: 4, max: 32, step: 1, description: 'Barcode length (10x Genomics: 12).' },
      { name: 'probability', label: 'Probability', type: 'range', default: 1.0, min: 0, max: 1, step: 0.05, description: 'Per-sequence probability.' },
    ],
    requires: ['Rearrange'],
    toPython: (p) =>
      `SimulateUMI(umi_length=${p.umi_length}, probability=${f(p.probability)})`,
  },

  {
    id: 'ReverseComplement',
    label: 'Reverse Complement',
    category: 'artifacts',
    description: 'Reverse-complement the sequence (~50% of NGS reads).',
    pythonClass: 'ReverseComplement',
    params: [
      { name: 'probability', label: 'Probability', type: 'range', default: 0.5, min: 0, max: 1, step: 0.05, description: 'Per-sequence probability.' },
    ],
    requires: ['Rearrange'],
    toPython: (p) => `ReverseComplement(probability=${f(p.probability)})`,
  },

  {
    id: 'SpikeContaminants',
    label: 'Spike Contaminants',
    category: 'artifacts',
    description: 'Replace sequences with phiX or random contaminant sequences.',
    pythonClass: 'SpikeContaminants',
    params: [
      { name: 'probability', label: 'Contamination rate', type: 'range', default: 0.01, min: 0, max: 0.1, step: 0.005, description: 'Per-sequence contamination probability.' },
      { name: 'contaminant_type', label: 'Type', type: 'select', default: 'random', options: [{ label: 'Random nucleotides', value: 'random' }, { label: 'phiX174', value: 'phix' }], description: 'Contaminant source.' },
    ],
    requires: ['Rearrange'],
    toPython: (p) =>
      `SpikeContaminants(probability=${f(p.probability)}, contaminant_type="${p.contaminant_type}")`,
  },

  {
    id: 'PrimerMask',
    label: 'Primer Mask',
    category: 'artifacts',
    description: 'Overwrite FR1 region with germline primer sequence.',
    pythonClass: 'PrimerMask',
    params: [
      { name: 'probability', label: 'Probability', type: 'range', default: 1.0, min: 0, max: 1, step: 0.05, description: 'Per-sequence probability.' },
    ],
    requires: ['Rearrange'],
    toPython: (p) => `PrimerMask(probability=${f(p.probability)})`,
  },

  {
    id: 'TrimToLength',
    label: 'Trim to Length',
    category: 'artifacts',
    description: "Enforce maximum sequence length by trimming from 5' end.",
    pythonClass: 'TrimToLength',
    params: [
      { name: 'max_length', label: 'Max length (nt)', type: 'number', default: 576, min: 50, max: 2000, step: 1, description: 'Maximum allowed sequence length.' },
    ],
    requires: ['Rearrange'],
    toPython: (p) => `TrimToLength(max_length=${p.max_length})`,
  },

  {
    id: 'InsertNs',
    label: 'Insert N-bases',
    category: 'artifacts',
    description: "Insert ambiguous 'N' bases by replacing random nucleotides.",
    pythonClass: 'InsertNs',
    params: [
      { name: 'n_ratio', label: 'N ratio', type: 'range', default: 0.02, min: 0, max: 0.1, step: 0.005, description: 'Fraction of positions replaced with N.' },
      { name: 'probability', label: 'Probability', type: 'range', default: 0.5, min: 0, max: 1, step: 0.05, description: 'Per-sequence probability.' },
    ],
    requires: ['Rearrange'],
    toPython: (p) =>
      `InsertNs(n_ratio=${f(p.n_ratio)}, probability=${f(p.probability)})`,
  },

  {
    id: 'InsertIndels',
    label: 'Insert Indels',
    category: 'artifacts',
    description: 'Insert insertions and deletions while avoiding NP regions.',
    pythonClass: 'InsertIndels',
    params: [
      { name: 'probability', label: 'Probability', type: 'range', default: 0.5, min: 0, max: 1, step: 0.05, description: 'Per-sequence probability.' },
      { name: 'max_indels', label: 'Max indels', type: 'number', default: 5, min: 1, max: 20, step: 1, description: 'Maximum indels per sequence.' },
    ],
    requires: ['Rearrange'],
    toPython: (p) =>
      `InsertIndels(probability=${f(p.probability)}, max_indels=${p.max_indels})`,
  },
];

// ── Lookup helpers ────────────────────────────────────────

const _byId = new Map(OP_DEFS.map((d) => [d.id, d]));

export function getOpDef(id: string): OpDef | undefined {
  return _byId.get(id);
}

export function getOpsByCategory(cat: Category): OpDef[] {
  return OP_DEFS.filter((d) => d.category === cat);
}

/** What dependency names does the current pipeline satisfy? */
export function satisfiedDeps(pipelineOpIds: string[]): Set<string> {
  const deps = new Set<string>();
  for (const id of pipelineOpIds) {
    const def = _byId.get(id);
    if (def?.satisfies) def.satisfies.forEach((s) => deps.add(s));
  }
  return deps;
}

/** Can this op be added given the current pipeline? Returns null if OK, or a reason string. */
export function canAddOp(
  opId: string,
  pipelineOpIds: string[],
): string | null {
  const def = _byId.get(opId);
  if (!def) return 'Unknown op.';

  // Check requires
  if (def.requires) {
    const sat = satisfiedDeps(pipelineOpIds);
    for (const req of def.requires) {
      if (!sat.has(req)) {
        // Build a human-readable message listing ops that satisfy this dependency
        const satisfiers = OP_DEFS.filter((d) => d.satisfies?.includes(req)).map((d) => d.label);
        const names = satisfiers.length > 0 ? satisfiers.join(' or ') : req;
        return `Requires ${names} first.`;
      }
    }
  }

  // Check excludes — skip for rearrangement ops that replace each other
  // (addOp in index.tsx swaps the rearrangement slot rather than adding alongside)
  if (def.excludes) {
    for (const exc of def.excludes) {
      if (pipelineOpIds.includes(exc)) {
        const excDef = _byId.get(exc);
        if (def.category === 'rearrangement' && excDef?.category === 'rearrangement') {
          continue; // rearrangement ops replace, not conflict
        }
        return `Conflicts with ${excDef?.label ?? exc}.`;
      }
    }
  }

  return null;
}

/** Build default params for an op. */
export function defaultParams(opId: string): Record<string, any> {
  const def = _byId.get(opId);
  if (!def) return {};
  const out: Record<string, any> = {};
  for (const p of def.params) {
    out[p.name] = p.default;
  }
  return out;
}
