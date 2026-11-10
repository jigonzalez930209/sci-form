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

      // Mermaid diagrams are rendered on the client by theme/Mermaid.vue.
      // Hand the source over as base64 so quotes and line breaks inside the
      // diagram cannot break the component attribute.
      const fence = md.renderer.rules.fence
      md.renderer.rules.fence = (tokens, idx, options, env, self) => {
        const token = tokens[idx]
        if (token.info.trim().toLowerCase() === 'mermaid') {
          const code = Buffer.from(token.content, 'utf8').toString('base64')
          return `<Mermaid code="${code}" />\n`
        }
        return fence ? fence(tokens, idx, options, env, self) : self.renderToken(tokens, idx, options)
      }
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
      { text: 'Experimental', link: '/experimental/overview' },
      { text: 'Validation', link: '/validation/scientific-validation' },
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
            { text: 'Live examples (run in the browser)', link: '/guide/live-examples' },
            { text: 'CLI', link: '/guide/cli' },
            { text: 'Metal EHT Validity', link: '/guide/metal-eht-validity' },
          ],
        },
        {
          text: 'Examples',
          items: [
            { text: '3D Reaction Viewer', link: '/guide/reaction-viewer' },
          ],
        },
        {
          text: 'Advanced',
          items: [
            { text: 'Reaction Backends', link: '/guide/reaction-backends' },
            { text: 'GPU backend portability', link: '/guide/gpu-portability' },
          ],
        },
        {
          text: 'Upgrading',
          items: [
            { text: 'Migrating to 2.0', link: '/guide/migration-2.0' },
            { text: 'Migrating to 1.0 (historical)', link: '/guide/migration-1.0' },
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
            { text: 'PM3 & xTB: Surfaces and Gradients', link: '/algorithm/semiempirical-gradients' },
            { text: 'Kohn–Sham DFT', link: '/algorithm/ks-dft' },
            { text: 'Machine Learning Potentials', link: '/algorithm/machine-learning-potentials' },
            { text: 'QM/MM', link: '/algorithm/qmmm' },
            { text: 'Electrostatic Potential', link: '/algorithm/electrostatic-potential' },
            { text: 'Density of States', link: '/algorithm/density-of-states' },
            { text: 'Population Analysis', link: '/algorithm/population-analysis' },
            { text: 'Dipole Moments', link: '/algorithm/dipole-moments' },
            { text: 'Redox & pKa', link: '/algorithm/redox-pka' },
            { text: 'Enhanced Sampling', link: '/algorithm/enhanced-sampling' },
            { text: 'Alchemical Free Energy', link: '/algorithm/free-energy' },
            { text: 'Kinetic Monte Carlo', link: '/algorithm/kmc' },
            { text: 'Active Learning', link: '/algorithm/active-learning' },
            { text: 'Strain Energy', link: '/algorithm/strain-energy' },
            { text: 'Molecular Alignment', link: '/algorithm/molecular-alignment' },
          ],
        },
        {
          text: 'Spectroscopy (Track D)',
          items: [
            { text: 'UV-Vis Spectroscopy', link: '/algorithm/uvvis-spectroscopy' },
            { text: 'Fluorescence', link: '/algorithm/fluorescence' },
            { text: 'IR Spectroscopy', link: '/algorithm/ir-spectroscopy' },
            { text: 'NMR Spectroscopy', link: '/algorithm/nmr-spectroscopy' },
            { text: 'Raman Spectroscopy', link: '/algorithm/raman' },
            { text: 'EPR Spectroscopy', link: '/algorithm/epr' },
          ],
        },
        {
          text: 'Materials',
          items: [
            { text: 'Materials Assembly', link: '/algorithm/materials-assembly' },
            { text: 'GCMC Adsorption', link: '/algorithm/gcmc' },
            { text: 'Web Transport', link: '/algorithm/web-transport' },
            { text: 'WebGPU Validation', link: '/algorithm/webgpu-validation' },
          ],
        },
        {
          text: 'Roadmaps',
          items: [
            { text: 'Quantum Validation Roadmap', link: '/algorithm/quantum-validation-roadmap' },
          ],
        },
      ],
      '/validation/': [
        {
          text: 'Validation',
          items: [
            {
              text: 'Scientific validation',
              link: '/validation/scientific-validation',
            },
            {
              text: 'Full matrix (129 rows)',
              link: '/validation/full-matrix',
            },
            {
              text: 'Cross-library plan',
              link: '/validation/cross-library-validation.es',
            },
          ],
        },
      ],
      '/experimental/': [
        {
          text: 'Experimental Modules',
          items: [
            { text: 'Overview', link: '/experimental/overview' },
            { text: 'Alpha tier', link: '/experimental/alpha' },
            { text: 'Beta tier', link: '/experimental/beta' },
            { text: 'Electrical double layer', link: '/experimental/edl' },
          ],
        },
        {
          text: 'Group A — Advanced Geometry',
          items: [
            { text: 'Conformal geometric algebra', link: '/experimental/cga' },
            { text: 'Randomized eigensolver', link: '/experimental/randnla' },
            { text: 'Riemannian geometry', link: '/experimental/riemannian' },
          ],
        },
        {
          text: 'Group B — Electronic Structure',
          items: [
            { text: 'Kernel Polynomial Method', link: '/experimental/kpm' },
            { text: 'Electronegativity equalization', link: '/experimental/eeq' },
            { text: 'Analytical solvation', link: '/experimental/alpb' },
          ],
        },
        {
          text: 'Group C — Precision Methods',
          items: [
            { text: 'Dispersion correction', link: '/experimental/d4' },
            { text: 'Semidefinite relaxation', link: '/experimental/sdr' },
            { text: 'Mobile Block Hessian', link: '/experimental/mbh' },
            { text: 'Constant Potential Method', link: '/experimental/cpm' },
            { text: 'Growing String Method', link: '/experimental/gsm' },
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
            { text: 'Complete Supported Reference', link: '/api/complete-supported-reference' },
          ],
        },
      ],
      '/': [
        {
          text: 'Roadmaps',
          items: [
            { text: 'v2.x Feature Roadmap', link: '/v2-roadmap' },
            { text: 'Dynamics & Reactivity', link: '/dynamics_reactivity_roadmap' },
            { text: 'NMR — Detailed Plan', link: '/nmr_detailed_roadmap' },
            { text: 'Spectroscopy (Track D)', link: '/spectroscopy_roadmap' },
          ],
        },
        {
          text: 'Releasing',
          items: [
            { text: 'v2.0.0 Release Checklist', link: '/release-checklist-2.0' },
            { text: 'v1.0.0 Release Checklist', link: '/release-checklist-1.0' },
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
