import React from 'react';
import Link from '@docusaurus/Link';
import CodeBlock from '@theme/CodeBlock';
import styles from './QuickStartPage.module.css';

const installCmd = `pip install GenAIRR`;

const verifyCode = `import GenAIRR
print(GenAIRR.__version__)`;

const simulateCode = `from GenAIRR import simulate

result = simulate("human_igh", n=3, seed=7)
print(result[0]["sequence"])
print(f"Length: {len(result[0]['sequence'])} nt")`;

const simulateOutput = `CAGGTCCAGCTGGTGCAGTCTGGGGCTGAGGTGAAGAAGCCTGGGTCCTCAGTGAAGGTCTCCTGAAAGGCTTCTGGAGG...
Length: 367 nt`;

const unmutatedCode = `from GenAIRR import simulate

result = simulate("human_igh", n=100, mutated=False)`;

const listConfigsCode = `from GenAIRR import list_configs
print(list_configs()[:10])`;

const useConfigCode = `from GenAIRR import HUMAN_IGH_OGRDB
result = simulate(HUMAN_IGH_OGRDB, n=10)`;

const protocolCode = `from GenAIRR import Protocol, Rearrange, Mutate, S5F, HUMAN_IGH_OGRDB

protocol = Protocol([
    Rearrange(),
    Mutate(S5F(0.01, 0.05)),
])

graph = protocol.compile(config=HUMAN_IGH_OGRDB, productive=True, seed=42)
result = graph.simulate(n=1000, progress=True)`;

export default function QuickStartPage(): JSX.Element {
  return (
    <div className={styles.page}>
      <section className={styles.hero}>
        <p className={styles.heroLead}>
          GenAIRR simulates immune receptor sequences with ground-truth annotations.
          This page gets you from install to your first sequences.
        </p>
        <div className={styles.heroActions}>
          <Link className={styles.primaryCta} to="/docs/getting-started/choosing-config">
            Choose a Config
          </Link>
          <Link className={styles.secondaryCta} to="/docs/getting-started/interpreting-results">
            Interpret Results
          </Link>
        </div>
      </section>

      <section className={styles.section}>
        <h2>Install</h2>
        <div className={styles.twoCol}>
          <div className={styles.card}>
            <div className={styles.cardTitle}>Install</div>
            <CodeBlock language="bash">{installCmd}</CodeBlock>
            <div className={styles.note}>Python 3.9 to 3.12 supported.</div>
          </div>
          <div className={styles.card}>
            <div className={styles.cardTitle}>Verify</div>
            <CodeBlock language="python">{verifyCode}</CodeBlock>
          </div>
        </div>
      </section>

      <section className={styles.section}>
        <h2>Simulate in 30 seconds</h2>
        <div className={styles.twoCol}>
          <div className={styles.card}>
            <div className={styles.cardTitle}>Code</div>
            <CodeBlock language="python">{simulateCode}</CodeBlock>
          </div>
          <div className={styles.card}>
            <div className={styles.cardTitle}>Output (seed=7)</div>
            <CodeBlock language="text">{simulateOutput}</CodeBlock>
            <div className={styles.note}>Sequence output is truncated for readability.</div>
          </div>
        </div>
        <div className={styles.noteWide}>
          <strong>Tip:</strong> `simulate` returns a `SimulationResult` (list-like). Use `result[0]` to access a record or `result.to_dataframe()` for pandas.
        </div>
      </section>

      <section className={styles.section}>
        <h2>Common next steps</h2>
        <div className={styles.pathGrid}>
          <div className={styles.card}>
            <div className={styles.cardTitle}>Unmutated sequences</div>
            <CodeBlock language="python">{unmutatedCode}</CodeBlock>
          </div>
          <div className={styles.card}>
            <div className={styles.cardTitle}>List available configs</div>
            <CodeBlock language="python">{listConfigsCode}</CodeBlock>
          </div>
          <div className={styles.card}>
            <div className={styles.cardTitle}>Use a specific config</div>
            <CodeBlock language="python">{useConfigCode}</CodeBlock>
          </div>
        </div>
      </section>

      <section className={styles.section}>
        <h2>When you need more control</h2>
        <p className={styles.sectionLead}>
          Use the Protocol API to assemble explicit ops, compile once, then run fast batches.
        </p>
        <CodeBlock language="python">{protocolCode}</CodeBlock>
        <div className={styles.nextRow}>
          <Link to="/docs/getting-started/protocol-api">Using the Protocol API</Link>
          <Link to="/docs/cookbook/basics/protocols-101">Protocols 101</Link>
        </div>
      </section>
    </div>
  );
}
