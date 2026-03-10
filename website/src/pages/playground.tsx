import React from 'react';
import Layout from '@theme/Layout';
import Playground from '@site/src/components/Playground';

export default function PlaygroundPage(): React.JSX.Element {
  return (
    <Layout
      title="Playground"
      description="Build GenAIRR protocols and simulate sequences in your browser"
    >
      <Playground />
    </Layout>
  );
}
