import React from 'react';
import {
  Target,
  Dna,
  Layers,
  Globe,
  Zap,
  AudioWaveform,
} from 'lucide-react';
import styles from './styles.module.css';

interface FeatureItem {
  icon: React.ReactNode;
  title: string;
  description: string;
}

const featureList: FeatureItem[] = [
  {
    icon: <Target size={20} />,
    title: 'Ground-Truth Annotations',
    description:
      '47 AIRR-format fields per record. Every V/D/J call, mutation position, and segment boundary is absolute truth — not statistically inferred.',
  },
  {
    icon: <Dna size={20} />,
    title: 'Realistic Biology',
    description:
      'D-gene inversion, receptor revision, class switch recombination, and antigen selection pressure — model the processes that shape real repertoires.',
  },
  {
    icon: <Layers size={20} />,
    title: 'Composable Experiment DSL',
    description:
      'Five biological phases from recombination to sequencing noise. A compiled C pipeline with automatic signal corrections handles the complexity for you.',
  },
  {
    icon: <AudioWaveform size={20} />,
    title: 'Full Wet-Lab Simulation',
    description:
      "5'/3' degradation, PCR amplification, quality profiles, indels, N-insertions, UMIs, and primer masking — simulate every stage of a sequencing experiment.",
  },
  {
    icon: <Globe size={20} />,
    title: '23 Species, 106 Configs',
    description:
      'Human, mouse, rabbit, cow, and 19 more. Pre-built germline references from IMGT and OGRDB — ready to simulate out of the box.',
  },
  {
    icon: <Zap size={20} />,
    title: 'C Engine, Zero Dependencies',
    description:
      '15,000\u201330,000 sequences/second on a single core. Pre-built wheels for Linux, macOS, and Windows with no mandatory runtime dependencies.',
  },
];

function Feature({icon, title, description}: FeatureItem) {
  return (
    <div className={styles.featureCard}>
      <div className={styles.featureIcon}>{icon}</div>
      <div>
        <h3 className={styles.featureTitle}>{title}</h3>
        <p className={styles.featureDesc}>{description}</p>
      </div>
    </div>
  );
}

export default function HomepageFeatures(): React.JSX.Element {
  return (
    <section className={styles.section}>
      <div className={styles.sectionInner}>
        <h2 className={styles.sectionTitle}>Built for Immunoinformatics</h2>
        <p className={styles.sectionSub}>
          Everything you need to generate realistic adaptive immune receptor
          repertoires for benchmarking, training, and analysis.
        </p>
        <div className={styles.features}>
          {featureList.map((props, idx) => (
            <Feature key={idx} {...props} />
          ))}
        </div>
      </div>
    </section>
  );
}
