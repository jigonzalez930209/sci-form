import { defineConfig } from 'vitepress'
import mathjax3 from 'markdown-it-mathjax3'

const customElements = [
  'mjx-container', 'mjx-assistive-mml', 'math', 'maction', 'maligngroup',
  'malignmark', 'menclose', 'merror', 'mfenced', 'mfrac', 'mi', 'mlongdiv',
  'mmultiscripts', 'mn', 'mo', 'mover', 'mpadded', 'mphantom', 'mroot',
  'mrow', 'ms', 'mscarries', 'mscarry', 'msgroup', 'mstack', 'msline',
  'mspace', 'msqrt', 'msrow', 'mstyle', 'msub', 'msup', 'msubsup',
  'mtable', 'mtd', 'mtext', 'mtr', 'munder', 'munderover', 'semantics',
  'annotation', 'annotation-xml',
]

export default defineConfig({
  base: '/sci-form/',
  title: 'sci-form',
  description: 'High-performance 3D molecular conformer generation from SMILES',
  
  head: [
    ['link', { rel: 'icon', type: 'image/svg+xml', href: '/logo.svg' }],
  ],

  markdown: {
    config: (md) => {
      md.use(mathjax3)
    },
  },

  vue: {
    template: {
      compilerOptions: {
        isCustomElement: (tag) => customElements.includes(tag),
      },
    },
  },

  themeConfig: {
    logo: '/logo.svg',
    
    nav: [
      { text: 'Guide', link: '/guide/getting-started' },
      { text: 'Algorithm', link: '/algorithm/overview' },
      { text: 'API', link: '/api/rust' },
      { text: 'Benchmarks', link: '/benchmarks' },
    ],

    sidebar: {
      '/guide/': [
        {
          text: 'Introduction',
          items: [
            { text: 'What is sci-form?', link: '/guide/what-is-sci-form' },
            { text: 'Getting Started', link: '/guide/getting-started' },
          ],
        },
        {
          text: 'Platforms',
          items: [
            { text: 'Rust', link: '/guide/rust' },
            { text: 'Python', link: '/guide/python' },
            { text: 'TypeScript / JS', link: '/guide/typescript' },
            { text: 'CLI', link: '/guide/cli' },
          ],
        },
      ],
      '/algorithm/': [
        {
          text: 'Algorithm',
          items: [
            { text: 'Pipeline Overview', link: '/algorithm/overview' },
            { text: 'SMILES Parsing', link: '/algorithm/smiles-parsing' },
            { text: 'Distance Geometry', link: '/algorithm/distance-geometry' },
            { text: 'Bounds Matrix', link: '/algorithm/bounds-matrix' },
            { text: 'Embedding', link: '/algorithm/embedding' },
            { text: 'Force Fields', link: '/algorithm/force-fields' },
            { text: 'ETKDG Refinement', link: '/algorithm/etkdg-refinement' },
            { text: 'SMARTS Matching', link: '/algorithm/smarts-matching' },
            { text: 'Validation', link: '/algorithm/validation' },
          ],
        },
      ],
      '/api/': [
        {
          text: 'API Reference',
          items: [
            { text: 'Rust', link: '/api/rust' },
            { text: 'Python', link: '/api/python' },
            { text: 'TypeScript / JS', link: '/api/typescript' },
            { text: 'CLI', link: '/api/cli' },
          ],
        },
      ],
    },

    socialLinks: [
      { icon: 'github', link: 'https://github.com/jigonzalez930209/sci-form' },
    ],

    search: {
      provider: 'local',
    },

    footer: {
      message: 'Released under the MIT License.',
      copyright: 'Copyright © 2024-2026 sci-form contributors',
    },

    outline: {
      level: [2, 3],
    },
  },
})
