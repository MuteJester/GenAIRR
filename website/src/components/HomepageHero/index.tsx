import React from 'react';
import Link from '@docusaurus/Link';
import {Github, ArrowRight} from 'lucide-react';
import styles from './styles.module.css';

export default function HomepageHero(): React.JSX.Element {
  return (
    <section className={styles.hero}>
      <div className={styles.heroInner}>
        <div className={styles.badge}>Open-Source AIRR Simulation Framework</div>

        <h1 className={styles.heroTitle}>
          Simulate immune receptor
          <br />
          sequences{' '}
          <span className={styles.highlight}>at scale</span>
        </h1>

        <p className={styles.heroSubtitle}>
          GenAIRR is a high-performance Python + C engine for simulating V(D)J
          recombination, somatic hypermutation, and sequencing artifacts —
          across 23 species with full ground-truth annotations.
        </p>

        <div className={styles.heroCta}>
          <Link
            className={`${styles.heroBtn} ${styles.heroBtnPrimary}`}
            to="/docs/getting-started/quick-start"
          >
            Get Started <ArrowRight size={15} />
          </Link>
          <Link
            className={`${styles.heroBtn} ${styles.heroBtnSecondary}`}
            href="https://github.com/MuteJester/GenAIRR"
          >
            <Github size={16} /> GitHub
          </Link>
        </div>

        <div className={styles.codeWindow}>
          <div className={styles.codeHeader}>
            <div className={styles.dots}>
              <span className={styles.dot} data-color="red" />
              <span className={styles.dot} data-color="yellow" />
              <span className={styles.dot} data-color="green" />
            </div>
            <span className={styles.codeLabel}>Python</span>
          </div>
          <pre className={styles.codeBody}>
            <code>
{`\
`}<span className={styles.keyword}>from</span>{` GenAIRR `}<span className={styles.keyword}>import</span>{` Experiment
`}<span className={styles.keyword}>from</span>{` GenAIRR.ops `}<span className={styles.keyword}>import</span>{` (
    rate, model, with_antigen_selection,
    with_5prime_loss, with_3prime_loss,
    with_indels, with_ns,
)

result = (
    `}<span className={styles.func}>Experiment</span>{`.`}<span className={styles.func}>on</span>{`(`}<span className={styles.string}>"human_igh"</span>{`)

    `}<span className={styles.comment}># Somatic hypermutation with selection pressure</span>{`
    .`}<span className={styles.func}>mutate</span>{`(
        `}<span className={styles.func}>model</span>{`(`}<span className={styles.string}>"s5f"</span>{`),
        `}<span className={styles.func}>rate</span>{`(`}<span className={styles.number}>0.01</span>{`, `}<span className={styles.number}>0.05</span>{`),
        `}<span className={styles.func}>with_antigen_selection</span>{`(`}<span className={styles.number}>0.5</span>{`),
    )

    `}<span className={styles.comment}># Sequencing artifacts</span>{`
    .`}<span className={styles.func}>sequence</span>{`(
        `}<span className={styles.func}>with_5prime_loss</span>{`(`}<span className={styles.attr}>min_remove</span>{`=`}<span className={styles.number}>5</span>{`, `}<span className={styles.attr}>max_remove</span>{`=`}<span className={styles.number}>30</span>{`),
        `}<span className={styles.func}>with_3prime_loss</span>{`(`}<span className={styles.attr}>min_remove</span>{`=`}<span className={styles.number}>5</span>{`, `}<span className={styles.attr}>max_remove</span>{`=`}<span className={styles.number}>20</span>{`),
    )

    `}<span className={styles.comment}># Post-sequencing noise</span>{`
    .`}<span className={styles.func}>observe</span>{`(
        `}<span className={styles.func}>with_indels</span>{`(`}<span className={styles.attr}>prob</span>{`=`}<span className={styles.number}>0.005</span>{`),
        `}<span className={styles.func}>with_ns</span>{`(`}<span className={styles.attr}>prob</span>{`=`}<span className={styles.number}>0.005</span>{`),
    )

    .`}<span className={styles.func}>run</span>{`(`}<span className={styles.attr}>n</span>{`=`}<span className={styles.number}>1000</span>{`, `}<span className={styles.attr}>seed</span>{`=`}<span className={styles.number}>42</span>{`)
)`}
            </code>
          </pre>
        </div>

        <div className={styles.stats}>
          <div className={styles.stat}>
            <span className={styles.statValue}>23</span>
            <span className={styles.statLabel}>Species</span>
          </div>
          <div className={styles.statDivider} />
          <div className={styles.stat}>
            <span className={styles.statValue}>106</span>
            <span className={styles.statLabel}>Built-in Configs</span>
          </div>
          <div className={styles.statDivider} />
          <div className={styles.stat}>
            <span className={styles.statValue}>15K+</span>
            <span className={styles.statLabel}>Seqs / Second</span>
          </div>
          <div className={styles.statDivider} />
          <div className={styles.stat}>
            <span className={styles.statValue}>47</span>
            <span className={styles.statLabel}>Output Fields</span>
          </div>
        </div>
      </div>
    </section>
  );
}
