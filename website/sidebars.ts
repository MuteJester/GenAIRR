import type {SidebarsConfig} from '@docusaurus/plugin-content-docs';

const sidebars: SidebarsConfig = {
  startHere: [
    'intro',
    {
      type: 'category',
      label: 'Start Here',
      collapsible: false,
      items: [
        'getting-started/quick-start',
        'getting-started/choosing-config',
        'getting-started/interpreting-results',
      ],
    },
  ],
  guides: [
    'guides/index',
    {
      type: 'category',
      label: 'Basics',
      collapsed: false,
      items: [
        'guides/basics/experiment-dsl',
        'guides/basics/chain-types',
        'guides/basics/export',
      ],
    },
    {
      type: 'category',
      label: 'Simulation Options',
      collapsed: false,
      items: [
        'guides/options/shm',
        'guides/options/artifacts',
        'guides/options/biology',
        'guides/options/productive',
        'guides/options/reproducibility',
      ],
    },
  ],
  concepts: [
    'concepts/index',
    'concepts/aseq-linked-list',
    'concepts/simulation-pipeline',
    'concepts/metadata-accuracy',
  ],
  showcase: [
    'showcase/index',
    'showcase/visualize-sequence',
    'showcase/narrate',
  ],
};

export default sidebars;
