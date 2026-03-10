import {themes as prismThemes} from 'prism-react-renderer';
import type {Config} from '@docusaurus/types';
import type * as Preset from '@docusaurus/preset-classic';

const config: Config = {
  title: 'GenAIRR Documentation',
  tagline: 'Modular Immunoglobulin Sequence Simulation',
  favicon: 'img/favicon.ico',


  url: 'https://mutejester.github.io',
  baseUrl: '/GenAIRR/',

  organizationName: 'MuteJester',
  projectName: 'GenAIRR',
  trailingSlash: false,

  onBrokenLinks: 'warn',

  i18n: {
    defaultLocale: 'en',
    locales: ['en'],
  },

  markdown: {
    mermaid: true,
  },
  themes: ['@docusaurus/theme-mermaid'],

  presets: [
    [
      'classic',
      {
        docs: {
          sidebarPath: './sidebars.ts',
          editUrl:
            'https://github.com/MuteJester/GenAIRR/tree/master/website/',
        },
        blog: false,
        theme: {
          customCss: './src/css/custom.css',
        },
      } satisfies Preset.Options,
    ],
  ],

  themeConfig: {
    colorMode: {
      defaultMode: 'light',
      disableSwitch: false,
      respectPrefersColorScheme: true,
    },
    navbar: {
      title: 'GenAIRR',
      logo: {
        alt: 'GenAIRR Logo',
        src: 'img/logo.svg',
      },
      items: [
        {
          type: 'docSidebar',
          sidebarId: 'startHere',
          position: 'left',
          label: 'Start Here',
        },
        {
          type: 'docSidebar',
          sidebarId: 'guides',
          position: 'left',
          label: 'Guides',
        },
        {
          type: 'docSidebar',
          sidebarId: 'concepts',
          position: 'left',
          label: 'Concepts',
        },
        {
          type: 'docSidebar',
          sidebarId: 'showcase',
          position: 'left',
          label: 'Showcase',
        },
{
          href: 'https://github.com/MuteJester/GenAIRR',
          label: 'GitHub',
          position: 'right',
        },
      ],
    },
    footer: {
      style: 'dark',
      links: [
        {
          title: 'Docs',
          items: [
            {
              label: 'Quick Start',
              to: '/docs/getting-started/quick-start',
            },
            {
              label: 'Docs Home',
              to: '/docs/',
            },
          ],
        },
        {
          title: 'Community',
          items: [
            {
              label: 'GitHub',
              href: 'https://github.com/MuteJester/GenAIRR',
            },
            {
              label: 'Issues',
              href: 'https://github.com/MuteJester/GenAIRR/issues',
            },
          ],
        },
      ],
      copyright: `If you use GenAIRR in your research, please cite: GenAIRR — <em>Briefings in Bioinformatics</em>, 2024. <a href="https://academic.oup.com/bib/article/25/6/bbae556/7863770" style="color: var(--ifm-footer-link-color)">DOI: 10.1093/bib/bbae556</a>`,
    },
    prism: {
      theme: prismThemes.github,
      darkTheme: prismThemes.dracula,
      additionalLanguages: ['bash', 'python'],
    },
  } satisfies Preset.ThemeConfig,
};

export default config;
