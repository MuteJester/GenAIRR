import React from 'react';
import {Shuffle, Dna, FlaskConical, Radio, Eye, ChevronDown} from 'lucide-react';
import styles from './styles.module.css';

interface PhaseItem {
  color: string;
  icon: React.ReactNode;
  method: string;
  title: string;
  description: string;
  ops: string[];
}

const phases: PhaseItem[] = [
  {
    color: '#1a6fac',
    icon: <Shuffle size={18} />,
    method: '.recombine()',
    title: 'V(D)J Recombination',
    description: 'Control biological events during gene rearrangement.',
    ops: ['with_d_inversion', 'with_receptor_revision'],
  },
  {
    color: '#0d9488',
    icon: <Dna size={18} />,
    method: '.mutate()',
    title: 'Somatic Hypermutation',
    description: 'Apply context-aware S5F mutation with optional CSR and selection pressure.',
    ops: ['model', 'rate', 'with_isotype_rates', 'with_antigen_selection'],
  },
  {
    color: '#7c3aed',
    icon: <FlaskConical size={18} />,
    method: '.prepare()',
    title: 'Library Preparation',
    description: 'Simulate wet-lab sample processing steps.',
    ops: ['with_primer_mask', 'with_umi', 'with_pcr'],
  },
  {
    color: '#ea580c',
    icon: <Radio size={18} />,
    method: '.sequence()',
    title: 'Sequencing Artifacts',
    description: "Model instrument-level degradation and error profiles.",
    ops: ['with_5prime_loss', 'with_3prime_loss', 'with_quality_profile'],
  },
  {
    color: '#b45309',
    icon: <Eye size={18} />,
    method: '.observe()',
    title: 'Post-Sequencing Noise',
    description: 'Add indels, N-bases, contaminants, and other observational artifacts.',
    ops: ['with_indels', 'with_ns', 'with_contaminants', 'with_reverse_complement'],
  },
];

function PhaseCard({color, icon, method, title, description, ops}: PhaseItem) {
  return (
    <div className={styles.phaseCard} style={{'--phase-color': color} as React.CSSProperties}>
      <div className={styles.phaseHeader}>
        <div className={styles.phaseIcon} style={{background: `${color}15`, color}}>
          {icon}
        </div>
        <code className={styles.phaseMethod} style={{color}}>{method}</code>
      </div>
      <h3 className={styles.phaseTitle}>{title}</h3>
      <p className={styles.phaseDesc}>{description}</p>
      <div className={styles.opsRow}>
        {ops.map((op) => (
          <span key={op} className={styles.opTag}>{op}()</span>
        ))}
      </div>
    </div>
  );
}

export default function HomepageHowItWorks(): React.JSX.Element {
  return (
    <section className={styles.section}>
      <div className={styles.sectionInner}>
        <h2 className={styles.sectionTitle}>The Experiment DSL</h2>
        <p className={styles.sectionSub}>
          Model every stage of a sequencing experiment as a composable pipeline.
          Include only the phases you need — each one is optional.
        </p>

        <div className={styles.pipelineContainer}>
          {/* Left: the code */}
          <div className={styles.pipelineCode}>
            <div className={styles.miniCodeWindow}>
              <div className={styles.miniCodeHeader}>
                <div className={styles.miniDots}>
                  <span className={styles.miniDot} data-color="red" />
                  <span className={styles.miniDot} data-color="yellow" />
                  <span className={styles.miniDot} data-color="green" />
                </div>
              </div>
              <pre className={styles.miniCodeBody}>
                <code>{`Experiment.on("human_igh")

    .recombine(...)

    .mutate(...)

    .prepare(...)

    .sequence(...)

    .observe(...)

    .run(n=1000, seed=42)`}</code>
              </pre>
            </div>
          </div>

          {/* Right: phase cards */}
          <div className={styles.pipelinePhases}>
            {phases.map((phase, idx) => (
              <React.Fragment key={idx}>
                <PhaseCard {...phase} />
                {idx < phases.length - 1 && (
                  <div className={styles.phaseConnector}>
                    <ChevronDown size={16} />
                  </div>
                )}
              </React.Fragment>
            ))}
          </div>
        </div>
      </div>
    </section>
  );
}
