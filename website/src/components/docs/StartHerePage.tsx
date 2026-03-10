import React from 'react';
import Link from '@docusaurus/Link';
import useBaseUrl from '@docusaurus/useBaseUrl';
import CodeBlock from '@theme/CodeBlock';
import styles from './WikiPage.module.css';

const quickCode = `from GenAIRR import simulate

result = simulate("human_igh", n=100, seed=7)
print(result[0]["sequence"])`;

const protocolCode = `from GenAIRR import Protocol, Rearrange, Mutate, S5F, HUMAN_IGH_OGRDB

protocol = Protocol([Rearrange(), Mutate(S5F(0.01, 0.05))])
graph = protocol.compile(config=HUMAN_IGH_OGRDB, seed=42)
result = graph.simulate(n=1000)`;

export default function StartHerePage(): JSX.Element {
  const flowSrc = useBaseUrl('/diagrams/start-here-flow.svg');
  return (
    <article className={styles.article}>
      <div className={styles.shortdesc}>Getting started guide for GenAIRR.</div>
      <div className={styles.hatnote}>
        For a full API summary, see the <Link to="/docs/reference/api-overview">API Overview</Link>.
      </div>
      <p className={styles.lead}>
        <strong>Start Here</strong> is the recommended entry point for GenAIRR. It explains the two
        main usage paths and links you to the pages you should read first.
      </p>

      <aside className={styles.infobox}>
        <table className={styles.infoboxTable}>
          <tbody>
            <tr>
              <th colSpan={2}>Start Here</th>
            </tr>
            <tr>
              <td className={styles.infoboxLabel}>Audience</td>
              <td className={styles.infoboxValue}>New users</td>
            </tr>
            <tr>
              <td className={styles.infoboxLabel}>Goal</td>
              <td className={styles.infoboxValue}>First successful run</td>
            </tr>
            <tr>
              <td className={styles.infoboxLabel}>Paths</td>
              <td className={styles.infoboxValue}>Quick / Protocol</td>
            </tr>
            <tr>
              <td className={styles.infoboxLabel}>Next</td>
              <td className={styles.infoboxValue}>
                <Link to="/docs/getting-started/quick-start">Quick Start</Link>
              </td>
            </tr>
            <tr>
              <td className={styles.infoboxLabel}>Reference</td>
              <td className={styles.infoboxValue}>
                <Link to="/docs/reference/api-overview">API Overview</Link>
              </td>
            </tr>
          </tbody>
        </table>
      </aside>

      <div className={styles.toc}>
        <div className={styles.tocTitle}>Contents</div>
        <ol className={styles.tocList}>
          <li><a href="#overview">Overview</a></li>
          <li><a href="#flow">Getting started flow</a></li>
          <li><a href="#path-a">Path A: Quick simulation</a></li>
          <li><a href="#path-b">Path B: Build a protocol</a></li>
          <li><a href="#learn-first">What to learn first</a></li>
          <li><a href="#reading-order">Suggested reading order</a></li>
        </ol>
      </div>

      <h2 id="overview">Overview</h2>
      <p>
        Most users start with the one-call API (<code>simulate</code>) to generate sequences fast.
        If you need full control, use the Protocol API to assemble explicit ops, compile once,
        and simulate large batches efficiently.
      </p>

      <h2 id="flow">Getting started flow</h2>
      <p>
        Use this diagram as the map for how the docs are organized and where to go next.
      </p>
      <figure className={styles.figure}>
        <img src={flowSrc} alt="Getting started flow" />
        <figcaption className={styles.figcaption}>
          Overview of the getting-started flow.
        </figcaption>
      </figure>

      <h2 id="path-a">Path A: Quick simulation</h2>
      <p>Use the one-call API when you want sequences fast.</p>
      <CodeBlock language="python">{quickCode}</CodeBlock>
      <ul>
        <li><Link to="/docs/getting-started/quick-start">Quick Start</Link></li>
        <li><Link to="/docs/getting-started/interpreting-results">Interpreting Results</Link></li>
      </ul>

      <h2 id="path-b">Path B: Build a protocol</h2>
      <p>Use the Protocol API when you need full control over the pipeline.</p>
      <CodeBlock language="python">{protocolCode}</CodeBlock>
      <ul>
        <li><Link to="/docs/getting-started/protocol-api">Using the Protocol API</Link></li>
        <li><Link to="/docs/cookbook/basics/protocols-101">Protocols 101</Link></li>
      </ul>

      <h2 id="learn-first">What you should learn first</h2>
      <ul>
        <li>How to choose a DataConfig and chain type</li>
        <li>How to interpret output fields and AIRR serialization</li>
        <li>How auto-injected corrections affect the final annotations</li>
      </ul>

      <h2 id="reading-order">Suggested reading order</h2>
      <ol>
        <li><Link to="/docs/getting-started/quick-start">Quick Start</Link></li>
        <li><Link to="/docs/getting-started/choosing-config">Choosing a Config</Link></li>
        <li><Link to="/docs/getting-started/interpreting-results">Interpreting Results</Link></li>
        <li><Link to="/docs/getting-started/protocol-api">Using the Protocol API</Link></li>
        <li><Link to="/docs/architecture/signals-corrections">Signals and Corrections</Link></li>
      </ol>

      <h2>See also</h2>
      <ul>
        <li><Link to="/docs/reference/parameters">Parameter Reference</Link></li>
        <li><Link to="/docs/reference/ops-reference">Ops Reference</Link></li>
        <li><Link to="/docs/reference/species-coverage">Species Coverage</Link></li>
      </ul>
    </article>
  );
}
