import React from 'react';
import CodeBlock from '@theme/CodeBlock';
import styles from './WikiPage.module.css';

const whatYouGetCode = `result = simulate("human_igh", n=3)
record = result[0]  # dict`;

const airrCode = `result = simulate("human_igh", n=3, airr=False)`;

const airrExample = `mutations type: str
mutations sample: {"65": "C>A", "91": "G>A", "193": "A>G", "196": "G>T"}
v_call type: str
v_call: IGHVF6-G21*12`;

const recordExample = `{
  "sequence": "CAGGTCCAGCTGGTGCAGTCTGGGGCTGAGGTGAAGAAGCCTGGGTCCTCAGTGAAGGTCTCCTGAAAGGCTTCT...",
  "v_call": ["IGHVF6-G21*12"],
  "d_call": ["IGHD2-21*02"],
  "j_call": ["IGHJ4*02"],
  "v_sequence_start": 0,
  "v_sequence_end": 296,
  "junction_sequence_start": 285,
  "junction_sequence_end": 332,
  "mutation_rate": 0.010899182561307902,
  "mutations": {65: "C>A", 91: "G>A", 193: "A>G", 196: "G>T"},
  "productive": false
}`;

const histCode = `from GenAIRR import simulate

res = simulate("human_igh", n=200, seed=13, airr=False)
vals = [r["mutation_rate"] for r in res]
print(f"min={min(vals):.4f} max={max(vals):.4f} mean={sum(vals)/len(vals):.4f}")`;

const histOutput = `min=0.0079 max=0.0495 mean=0.0289`;

export default function InterpretingResultsPage(): JSX.Element {
  return (
    <article className={styles.article}>
      <p className={styles.lead}>
        GenAIRR returns rich, ground-truth annotations. This page explains what you get and how to read it.
      </p>

      <aside className={styles.infobox}>
        <div className={styles.infoboxTitle}>Result facts</div>
        <div className={styles.infoboxTable}>
          <div className={styles.infoboxLabel}>Single run</div>
          <div>SimulationContainer</div>
          <div className={styles.infoboxLabel}>Batch run</div>
          <div>SimulationResult</div>
          <div className={styles.infoboxLabel}>Indexing</div>
          <div>0-based</div>
          <div className={styles.infoboxLabel}>AIRR</div>
          <div>Default</div>
        </div>
      </aside>

      <h2>What you get</h2>
      <ul>
        <li><code>graph.execute()</code> returns a single <code>SimulationContainer</code>.</li>
        <li><code>simulate(...)</code> and <code>graph.simulate(n)</code> return a <code>SimulationResult</code>.</li>
      </ul>
      <CodeBlock language="python">{whatYouGetCode}</CodeBlock>

      <h2>AIRR serialization</h2>
      <p>
        By default, batch results use AIRR-friendly serialization:
      </p>
      <ul>
        <li><code>v_call</code>, <code>d_call</code>, <code>j_call</code>, <code>c_call</code> are comma-separated strings</li>
        <li>dict fields like <code>mutations</code>, <code>indels</code>, <code>Ns</code> are JSON strings</li>
      </ul>
      <p>To keep raw Python types, use <code>airr=False</code>:</p>
      <CodeBlock language="python">{airrCode}</CodeBlock>
      <CodeBlock language="text">{airrExample}</CodeBlock>

      <h2>Key fields</h2>
      <p>Positions are 0-based indices into the final <code>sequence</code>.</p>
      <table>
        <thead>
          <tr>
            <th>Field</th>
            <th>Meaning</th>
          </tr>
        </thead>
        <tbody>
          <tr><td><code>sequence</code></td><td>Final nucleotide sequence.</td></tr>
          <tr><td><code>sequence_aa</code></td><td>Amino-acid translation of the full sequence.</td></tr>
          <tr><td><code>v_call</code>, <code>d_call</code>, <code>j_call</code></td><td>Allele calls. The first entry is the true allele; later entries are indistinguishable alternatives.</td></tr>
          <tr><td><code>v_sequence_start/end</code></td><td>Segment boundaries in the final sequence.</td></tr>
          <tr><td><code>junction_sequence_start/end</code></td><td>Junction boundaries.</td></tr>
          <tr><td><code>junction_nt</code>, <code>junction_aa</code></td><td>Junction sequence in nt and aa.</td></tr>
          <tr><td><code>mutation_rate</code></td><td>Fraction of mutated positions.</td></tr>
          <tr><td><code>mutations</code></td><td>Dict of <code>{`{position: "A>T"}`}</code> changes.</td></tr>
          <tr><td><code>indels</code>, <code>Ns</code>, <code>sequencing_errors</code>, <code>pcr_errors</code></td><td>Artifact metadata.</td></tr>
          <tr><td><code>productive</code>, <code>stop_codon</code>, <code>vj_in_frame</code>, <code>note</code></td><td>Functionality assessment.</td></tr>
        </tbody>
      </table>

      <h2>Example record (shortened, airr=False, seed=7)</h2>
      <CodeBlock language="python">{recordExample}</CodeBlock>

      <div className={styles.note}>
        If <code>v_call</code> includes multiple entries, those alleles are indistinguishable given trimming and mutations.
        <br />
        <code>Short-D</code> in <code>d_call</code> means the D segment is too short to identify reliably.
      </div>

      <h2>Export helpers</h2>
      <CodeBlock language="python">{`result.to_dataframe()
result.to_csv("output.csv")
result.to_fasta("output.fasta")`}</CodeBlock>

      <h2>Example plot: mutation rate distribution</h2>
      <figure className={styles.figure}>
        <img src="/plots/mutation_rate_hist.svg" alt="Mutation rate histogram" />
        <figcaption className={styles.figcaption}>
          Mutation-rate variability across 200 sequences (seed=13).
        </figcaption>
      </figure>
      <CodeBlock language="python">{histCode}</CodeBlock>
      <CodeBlock language="text">{histOutput}</CodeBlock>
    </article>
  );
}
