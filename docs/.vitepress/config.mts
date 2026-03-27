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
  
  ignoreDeadLinks: [
    // Ignore links to source code files outside /docs/
    /^\.\.\/\.\.\//,
  ],
  
  head: [
    ['link', { rel: 'icon', type: 'image/svg+xml', href: '/sci-form/logo.svg' }],
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
    logo: '/sci-form/logo.svg',
    
    nav: [
      { text: 'Guide', link: '/guide/getting-started' },
      { text: 'Algorithm', link: '/algorithm/overview' },
      { text: 'Experimental', link: '/experimental/overview' },
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
            { text: 'Metal EHT Validity', link: '/guide/metal-eht-validity' },
          ],
        },
      ],
      '/algorithm/': [
        {
          text: 'Conformer Generation',
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
        {
          text: 'Properties & Analysis',
          items: [
            { text: 'Extended Hückel Theory', link: '/algorithm/extended-huckel-theory' },
            { text: 'HF-3c Quantum Engine', link: '/algorithm/hf3c-quantum-engine' },
            { text: 'Machine Learning Potentials', link: '/algorithm/machine-learning-potentials' },
            { text: 'Electrostatic Potential', link: '/algorithm/electrostatic-potential' },
            { text: 'Density of States', link: '/algorithm/density-of-states' },
            { text: 'Population Analysis', link: '/algorithm/population-analysis' },
            { text: 'Dipole Moments', link: '/algorithm/dipole-moments' },
            { text: 'Strain Energy', link: '/algorithm/strain-energy' },
            { text: 'Molecular Alignment', link: '/algorithm/molecular-alignment' },
          ],
        },
        {
          text: 'Spectroscopy (Track D)',
          items: [
            { text: 'UV-Vis Spectroscopy', link: '/algorithm/uvvis-spectroscopy' },
            { text: 'IR Spectroscopy', link: '/algorithm/ir-spectroscopy' },
            { text: 'NMR Spectroscopy', link: '/algorithm/nmr-spectroscopy' },
          ],
        },
        {
          text: 'Materials',
          items: [
            { text: 'Materials Assembly', link: '/algorithm/materials-assembly' },
            { text: 'Web Transport', link: '/algorithm/web-transport' },
            { text: 'WebGPU Validation', link: '/algorithm/webgpu-validation' },
          ],
        },
      ],
      '/experimental/': [
        {
          text: 'Experimental Modules',
          items: [
            { text: 'Overview', link: '/experimental/overview' },
          ],
        },
        {
          text: 'Group A — Advanced Geometry',
          items: [
            { text: 'E1: CGA', link: '/experimental/cga' },
            { text: 'E2: RandNLA', link: '/experimental/randnla' },
            { text: 'E3: Riemannian', link: '/experimental/riemannian' },
          ],
        },
        {
          text: 'Group B — Electronic Structure',
          items: [
            { text: 'E4: KPM', link: '/experimental/kpm' },
            { text: 'E5: EEQ', link: '/experimental/eeq' },
            { text: 'E6: ALPB', link: '/experimental/alpb' },
          ],
        },
        {
          text: 'Group C — Precision Methods',
          items: [
            { text: 'E7: D4', link: '/experimental/d4' },
            { text: 'E8: SDR', link: '/experimental/sdr' },
            { text: 'E9: MBH', link: '/experimental/mbh' },
            { text: 'E10: CPM', link: '/experimental/cpm' },
            { text: 'E11: GSM', link: '/experimental/gsm' },
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
