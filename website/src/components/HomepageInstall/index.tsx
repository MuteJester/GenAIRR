import React from 'react';
import {Copy} from 'lucide-react';
import styles from './styles.module.css';

export default function HomepageInstall(): React.JSX.Element {
  const [copied, setCopied] = React.useState(false);

  const handleCopy = () => {
    navigator.clipboard.writeText('pip install GenAIRR');
    setCopied(true);
    setTimeout(() => setCopied(false), 2000);
  };

  return (
    <>
      <section className={`${styles.section} ${styles.sectionAlt}`}>
        <div className={styles.sectionInner}>
          <div className={styles.installBlock}>
            <h2 className={styles.sectionTitle}>Quick Install</h2>
            <div className={styles.installCode}>
              <code>pip install GenAIRR</code>
              <button
                className={styles.copyBtn}
                onClick={handleCopy}
                title="Copy to clipboard"
                aria-label="Copy install command"
              >
                {copied ? (
                  <span className={styles.copiedLabel}>Copied!</span>
                ) : (
                  <Copy size={14} />
                )}
              </button>
            </div>
            <p className={styles.installNote}>
              Python 3.9+ &middot; Pre-built wheels for Linux, macOS, Windows &middot; Zero mandatory dependencies
            </p>
          </div>
        </div>
      </section>

      <section className={`${styles.section} ${styles.cite}`}>
        <div className={styles.sectionInner}>
          <p className={styles.citeText}>
            If you use GenAIRR in your research, please cite:{' '}
            <strong>GenAIRR</strong> &mdash;{' '}
            <em>Briefings in Bioinformatics</em>, 2024.{' '}
            <a href="https://academic.oup.com/bib/article/25/6/bbae556/7863770">
              DOI: 10.1093/bib/bbae556
            </a>
          </p>
        </div>
      </section>
    </>
  );
}
