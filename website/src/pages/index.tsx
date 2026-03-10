import React from 'react';
import Layout from '@theme/Layout';
import HomepageHero from '@site/src/components/HomepageHero';
import HomepageFeatures from '@site/src/components/HomepageFeatures';
import HomepageHowItWorks from '@site/src/components/HomepageHowItWorks';
import HomepagePersonas from '@site/src/components/HomepagePersonas';
import HomepageInstall from '@site/src/components/HomepageInstall';

export default function Home(): React.JSX.Element {
  return (
    <Layout
      title="Simulate Immune Receptor Sequences at Scale"
      description="GenAIRR — Modular Immunoglobulin Sequence Simulation Framework"
    >
      <HomepageHero />
      <main>
        <HomepageFeatures />
        <HomepageHowItWorks />
        <HomepagePersonas />
        <HomepageInstall />
      </main>
    </Layout>
  );
}
