import DefaultTheme from 'vitepress/theme'
import type { Theme } from 'vitepress'
import SvgDiagram from './SvgDiagram.vue'
import LiveExample from './LiveExample.vue'
import LiveExamples from './LiveExamples.vue'
import Mermaid from './Mermaid.vue'
import './custom.css'

export default {
  extends: DefaultTheme,
  enhanceApp({ app }) {
    app.component('SvgDiagram', SvgDiagram)
    app.component('LiveExample', LiveExample)
    app.component('LiveExamples', LiveExamples)
    app.component('Mermaid', Mermaid)
  },
} satisfies Theme
