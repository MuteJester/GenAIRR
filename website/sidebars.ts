import type {SidebarsConfig} from '@docusaurus/plugin-content-docs';

const sidebars: SidebarsConfig = {
  startHere: [
    'intro',
    {
      type: 'category',
      label: 'Getting Started',
      collapsed: false,
      items: [
        'getting-started/quick-start',
        'getting-started/choosing-config',
        'getting-started/interpreting-results',
      ],
    },
    {
      type: 'category',
      label: 'Core Concepts',
      items: [
        'concepts/index',
        'concepts/persistent-ir',
        'concepts/simulation-pipeline',
        'concepts/sampling-mechanics',
        'concepts/metadata-accuracy',
      ],
    },
  ],
  guides: [
    'guides/index',
    {
      type: 'category',
      label: 'Simulation Basics',
      collapsed: false,
      items: [
        'guides/basics/experiment-dsl',
        'guides/basics/chain-types',
        'guides/basics/export',
      ],
    },
    {
      type: 'category',
      label: 'Biological Variation',
      collapsed: false,
      items: [
        'guides/options/shm',
        'guides/options/biology',
        'guides/options/clonal-families',
      ],
    },
    {
      type: 'category',
      label: 'Technical Realism',
      collapsed: false,
      items: [
        'guides/options/artifacts',
        'guides/options/reproducibility',
      ],
    },
    {
      type: 'category',
      label: 'Safety & Constraints',
      collapsed: false,
      items: [
        'guides/options/productive',
      ],
    },
  ],
  showcase: [
    'showcase/index',
    'showcase/visualize-sequence',
    'showcase/narrate',
  ],
  reference: [
    'reference/index',
  ],
  help: [
    'help/index',
  ],
};

export default sidebars;
