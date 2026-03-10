import React from 'react';
import Link from '@docusaurus/Link';
import CodeBlock from '@theme/CodeBlock';
import styles from './ChoosingConfigPage.module.css';

const aliasCode = `from GenAIRR import simulate

result = simulate("human_igh", n=100)`;

const directCode = `from GenAIRR import HUMAN_IGH_OGRDB
result = simulate(HUMAN_IGH_OGRDB, n=100)`;

const listConfigsCode = `from GenAIRR import list_configs
print(list_configs()[:10])`;

const listConfigsOutput = `['ALPACA_IGH_IMGT', 'CAT_IGK_IMGT', 'CAT_IGL_IMGT', 'CAT_TCRA_IMGT',
 'CAT_TCRB_IMGT', 'CAT_TCRD_IMGT', 'CAT_TCRG_IMGT', 'CHICKEN_IGH_IMGT',
 'CHICKEN_IGL_IMGT', 'COW_IGH_IMGT']`;

const metadataCode = `from GenAIRR import HUMAN_IGH_OGRDB
print(HUMAN_IGH_OGRDB.metadata)`;

const metadataOutput = `Species: Human, Chain: BCR_HEAVY, Has D: True, Reference: 'OGRDB V8', Updated: 2025-02-01`;

const dSegmentCode = `from GenAIRR import simulate

igh = simulate("human_igh", n=1, seed=3)
igk = simulate("human_igk", n=1, seed=3)
print("IGH d_call:", igh[0]["d_call"])
print("IGK d_call:", igk[0]["d_call"])`;

const dSegmentOutput = `IGH d_call: IGHD2-21*01,IGHD2-21*02
IGK d_call:`;

export default function ChoosingConfigPage(): JSX.Element {
  return (
    <div className={styles.page}>
      <section className={styles.hero}>
        <p className={styles.heroLead}>
          A DataConfig contains germline alleles, trimming distributions, NP models, and metadata
          for a specific species and chain type. Picking the right config is the first real
          decision in GenAIRR.
        </p>
      </section>

      <section className={styles.section}>
        <h2>Two ways to choose</h2>
        <div className={styles.twoCol}>
          <div className={styles.card}>
            <div className={styles.cardTitle}>Use a string alias</div>
            <CodeBlock language="python">{aliasCode}</CodeBlock>
            <div className={styles.note}>Aliases are case-insensitive. Hyphens become underscores.</div>
          </div>
          <div className={styles.card}>
            <div className={styles.cardTitle}>Import a config directly</div>
            <CodeBlock language="python">{directCode}</CodeBlock>
          </div>
        </div>
      </section>

      <section className={styles.section}>
        <h2>OGRDB vs IMGT</h2>
        <div className={styles.keyGrid}>
          <div className={styles.keyItem}>OGRDB: curated for human BCR (IGH, IGK, IGL)</div>
          <div className={styles.keyItem}>IMGT: broad coverage across species and TCR chains</div>
        </div>
        <div className={styles.noteWide}>
          Default aliases map human BCR chains to OGRDB and human TCR chains to IMGT.
        </div>
      </section>

      <section className={styles.section}>
        <h2>Common aliases</h2>
        <div className={styles.aliasTable}>
          <div className={styles.aliasRow}><code>human_igh</code><span>HUMAN_IGH_OGRDB</span></div>
          <div className={styles.aliasRow}><code>human_igk</code><span>HUMAN_IGK_OGRDB</span></div>
          <div className={styles.aliasRow}><code>human_igl</code><span>HUMAN_IGL_OGRDB</span></div>
          <div className={styles.aliasRow}><code>human_tcrb</code><span>HUMAN_TCRB_IMGT</span></div>
          <div className={styles.aliasRow}><code>mouse_igh</code><span>MOUSE_IGH_IMGT</span></div>
          <div className={styles.aliasRow}><code>rabbit_igh</code><span>RABBIT_IGH_IMGT</span></div>
        </div>
        <div className={styles.noteWide}>
          See full mapping in <Link to="/docs/reference/config-aliases">Config Aliases</Link>.
        </div>
      </section>

      <section className={styles.section}>
        <h2>Chain types and D segments</h2>
        <div className={styles.chainGrid}>
          <div className={styles.chainCard}>
            <div className={styles.chainTitle}>VDJ (has D segment)</div>
            <div>IGH, TCRB, TCRD</div>
          </div>
          <div className={styles.chainCard}>
            <div className={styles.chainTitle}>VJ (no D segment)</div>
            <div>IGK, IGL, TCRA, TCRG</div>
          </div>
        </div>
        <div className={styles.twoCol}>
          <div className={styles.card}>
            <div className={styles.cardTitle}>Example</div>
            <CodeBlock language="python">{dSegmentCode}</CodeBlock>
          </div>
          <div className={styles.card}>
            <div className={styles.cardTitle}>Output (seed=3)</div>
            <CodeBlock language="text">{dSegmentOutput}</CodeBlock>
          </div>
        </div>
      </section>

      <section className={styles.section}>
        <h2>List all configs</h2>
        <div className={styles.twoCol}>
          <div className={styles.card}>
            <div className={styles.cardTitle}>Code</div>
            <CodeBlock language="python">{listConfigsCode}</CodeBlock>
          </div>
          <div className={styles.card}>
            <div className={styles.cardTitle}>Output (first 10)</div>
            <CodeBlock language="text">{listConfigsOutput}</CodeBlock>
          </div>
        </div>
        <div className={styles.noteWide}>
          For detailed coverage by species, see <Link to="/docs/reference/species-coverage">Species Coverage</Link>.
        </div>
      </section>

      <section className={styles.section}>
        <h2>Inspect config metadata</h2>
        <div className={styles.twoCol}>
          <div className={styles.card}>
            <div className={styles.cardTitle}>Code</div>
            <CodeBlock language="python">{metadataCode}</CodeBlock>
          </div>
          <div className={styles.card}>
            <div className={styles.cardTitle}>Output</div>
            <CodeBlock language="text">{metadataOutput}</CodeBlock>
          </div>
        </div>
      </section>
    </div>
  );
}
