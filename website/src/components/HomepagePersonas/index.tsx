import React from 'react';
import Link from '@docusaurus/Link';
import {FlaskConical, BrainCircuit, Code2} from 'lucide-react';
import styles from './styles.module.css';

interface PersonaItem {
  variant: 'bio' | 'ml' | 'dev';
  icon: React.ReactNode;
  title: string;
  description: string;
  link: string;
}

const personas: PersonaItem[] = [
  {
    variant: 'bio',
    icon: <FlaskConical size={22} />,
    title: 'Bioinformatician',
    description:
      'Generate benchmarking datasets with realistic V(D)J recombination, SHM, and sequencing artifacts across 23 species.',
    link: '/docs/getting-started/quick-start',
  },
  {
    variant: 'ml',
    icon: <BrainCircuit size={22} />,
    title: 'ML Engineer',
    description:
      'Create training data at scale with 47 ground-truth labels, deterministic seeding, and DataFrame export.',
    link: '/docs/getting-started/quick-start',
  },
  {
    variant: 'dev',
    icon: <Code2 size={22} />,
    title: 'Developer',
    description:
      'Build custom DataConfigs from your own germline databases and integrate GenAIRR into analysis pipelines.',
    link: '/docs/reference/',
  },
];

function PersonaCard({variant, icon, title, description, link}: PersonaItem) {
  return (
    <Link
      to={link}
      className={`${styles.pathCard} ${styles[`pathCard_${variant}`]}`}
    >
      <div className={`${styles.pathIcon} ${styles[`pathIcon_${variant}`]}`}>
        {icon}
      </div>
      <h3 className={styles.pathTitle}>{title}</h3>
      <p className={styles.pathDesc}>{description}</p>
      <span className={styles.pathLink}>Learn more &rarr;</span>
    </Link>
  );
}

export default function HomepagePersonas(): React.JSX.Element {
  return (
    <section className={`${styles.section} ${styles.sectionAlt}`}>
      <div className={styles.sectionInner}>
        <h2 className={styles.sectionTitle}>Choose Your Path</h2>
        <div className={styles.paths}>
          {personas.map((props, idx) => (
            <PersonaCard key={idx} {...props} />
          ))}
        </div>
      </div>
    </section>
  );
}
