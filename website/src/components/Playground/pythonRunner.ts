/* ───────────────────────────────────────────────────────────
 *  pythonRunner.ts — Convert pipeline state into executable
 *  Python code for Pyodide.
 * ─────────────────────────────────────────────────────────── */

import { getOpDef, type OpDef } from './ops';

export interface PipelineOp {
  instanceId: string;
  opId: string;
  params: Record<string, any>;
}

export interface RunConfig {
  pipeline: PipelineOp[];
  useCustomConfig: boolean;
  n: number;
  seed: number | null;
  productive: boolean;
}

/**
 * Generate the Python script that Pyodide will execute.
 * Returns the script as a string. The last expression evaluates
 * to a JSON string of results.
 */
export function buildPythonScript(cfg: RunConfig): string {
  const defs = cfg.pipeline
    .map((op) => getOpDef(op.opId))
    .filter((d): d is OpDef => d != null);

  // Collect unique imports
  const opClasses = new Set(defs.map((d) => d.pythonClass));
  const extraImports = new Set<string>();
  for (const d of defs) {
    d.extraImports?.forEach((i) => extraImports.add(i));
  }

  const opImportLine = `from GenAIRR.ops import ${[...opClasses].join(', ')}`;
  const allImports = [
    opImportLine,
    ...extraImports,
    'from GenAIRR.protocol import Protocol',
    'import json',
  ];

  // Build op constructor lines
  const opLines = cfg.pipeline.map((op) => {
    const def = getOpDef(op.opId);
    if (!def) return `# Unknown op: ${op.opId}`;
    return `    ${def.toPython(op.params)},`;
  });

  // Config resolution
  const configLine = cfg.useCustomConfig
    ? '_custom_config'
    : "getattr(__import__('GenAIRR.data', fromlist=['HUMAN_IGH_IMGT']), 'HUMAN_IGH_IMGT')";

  const seedArg = cfg.seed !== null ? `seed=${cfg.seed}, ` : '';
  const productiveArg = cfg.productive ? 'productive=True, ' : '';

  const script = `
${allImports.join('\n')}

# Build protocol
protocol = Protocol([
${opLines.join('\n')}
])

# Compile and run
config = ${configLine}
graph = protocol.compile(config=config, ${productiveArg}${seedArg})
result = graph.simulate(n=${cfg.n}, airr=True)

# Extract results as JSON
_results = [dict(r) for r in result]
json.dumps(_results)
`.trim();

  return script;
}

/**
 * Generate the "display" version of the Python code — the
 * copy-pasteable snippet shown to the user. Same logic but
 * with nicer formatting and comments.
 */
export function buildDisplayScript(cfg: RunConfig): string {
  const defs = cfg.pipeline
    .map((op) => getOpDef(op.opId))
    .filter((d): d is OpDef => d != null);

  const opClasses = new Set(defs.map((d) => d.pythonClass));
  const extraImports = new Set<string>();
  for (const d of defs) {
    d.extraImports?.forEach((i) => extraImports.add(i));
  }

  const opImportLine = `from GenAIRR.ops import ${[...opClasses].join(', ')}`;
  const allImports = [opImportLine, ...extraImports, 'from GenAIRR.protocol import Protocol'];

  const opLines = cfg.pipeline.map((op) => {
    const def = getOpDef(op.opId);
    if (!def) return `    # Unknown op: ${op.opId}`;
    return `    ${def.toPython(op.params)},`;
  });

  const configStr = cfg.useCustomConfig
    ? 'my_config  # your uploaded DataConfig'
    : 'from GenAIRR.data import HUMAN_IGH_IMGT\nconfig = HUMAN_IGH_IMGT';

  const compileArgs: string[] = ['config=config'];
  if (cfg.productive) compileArgs.push('productive=True');
  if (cfg.seed !== null) compileArgs.push(`seed=${cfg.seed}`);

  let configBlock: string;
  if (cfg.useCustomConfig) {
    configBlock = `import pickle\nwith open("my_config.pkl", "rb") as f:\n    config = pickle.load(f)`;
  } else {
    configBlock = `from GenAIRR.data import HUMAN_IGH_IMGT\nconfig = HUMAN_IGH_IMGT`;
  }

  return `${allImports.join('\n')}

${configBlock}

protocol = Protocol([
${opLines.join('\n')}
])

graph = protocol.compile(${compileArgs.join(', ')})
result = graph.simulate(n=${cfg.n}, airr=True)
df = result.to_dataframe()
print(df.head())`;
}
