import DefaultTheme from 'vitepress/theme'
import type { Theme } from 'vitepress'
import SvgDiagram from './SvgDiagram.vue'
import './custom.css'

export default {
  extends: DefaultTheme,
  enhanceApp({ app }) {
    app.component('SvgDiagram', SvgDiagram)
  },
} satisfies Theme
