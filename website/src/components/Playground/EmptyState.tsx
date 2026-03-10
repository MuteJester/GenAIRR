/* ───────────────────────────────────────────────────────────
 *  EmptyState.tsx — Empty results placeholder with a
 *  "Try This" protocol hint. Swiss Lab aesthetic.
 * ─────────────────────────────────────────────────────────── */

import React from 'react';
import { FlaskConical, ArrowRight } from 'lucide-react';
import styles from './styles.module.css';

export default function EmptyState() {
  return (
    <div className={styles.emptyState}>
      <FlaskConical size={32} style={{ color: 'var(--pg-border)' }} strokeWidth={1.5} />

      <div style={{ display: 'flex', flexDirection: 'column', alignItems: 'center', gap: '0.5rem', textAlign: 'center' }}>
        <p className={styles.emptyStateTitle}>No results yet</p>
        <p className={styles.emptyStateDesc}>
          Build your protocol above by adding operations, then click Run Simulation
          to generate synthetic antibody sequences.
        </p>
      </div>

      {/* Sample hint */}
      <div className={styles.emptyStateHint}>
        <p className={styles.emptyStateHintLabel}>Try This</p>
        <div className={styles.emptyStateFlow}>
          <span className={styles.emptyStateChip}>Rearrange</span>
          <ArrowRight size={12} style={{ color: 'var(--pg-muted-fg)' }} />
          <span className={styles.emptyStateChip}>Somatic Hypermutation</span>
          <ArrowRight size={12} style={{ color: 'var(--pg-muted-fg)' }} />
          <span className={styles.emptyStateChip}>Quality Errors (Illumina SE)</span>
        </div>
        <p className={styles.emptyStateHintText}>
          A typical 3-step protocol: rearrange V(D)J, introduce SHM, then add sequencing artifacts.
        </p>
      </div>
    </div>
  );
}
