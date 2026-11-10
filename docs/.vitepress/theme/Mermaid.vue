<script setup lang="ts">
import { onMounted, ref, watch } from 'vue'
import { useData } from 'vitepress'

const props = defineProps<{ code: string }>()

const { isDark } = useData()
const svg = ref('')
const failed = ref(false)
let seq = 0

function decode(value: string): string {
  const bytes = Uint8Array.from(atob(value), (char) => char.charCodeAt(0))
  return new TextDecoder().decode(bytes)
}

function themeVariables(dark: boolean) {
  return dark
    ? {
        background: 'transparent',
        primaryColor: '#1b2b4a',
        primaryTextColor: '#e6edf7',
        primaryBorderColor: '#3b82f6',
        secondaryColor: '#123035',
        secondaryTextColor: '#e6edf7',
        secondaryBorderColor: '#2dd4bf',
        tertiaryColor: '#172033',
        tertiaryTextColor: '#e6edf7',
        tertiaryBorderColor: '#475569',
        lineColor: '#94a3b8',
        textColor: '#e6edf7',
        edgeLabelBackground: '#1b1b1f',
        clusterBkg: 'rgba(96,165,250,.08)',
        clusterBorder: '#475569',
        fontFamily: 'ui-sans-serif, system-ui, -apple-system, "Segoe UI", Roboto, sans-serif',
        fontSize: '14px',
      }
    : {
        background: 'transparent',
        primaryColor: '#eef4ff',
        primaryTextColor: '#0f172a',
        primaryBorderColor: '#93b4f7',
        secondaryColor: '#e8faf6',
        secondaryTextColor: '#0f172a',
        secondaryBorderColor: '#5eead4',
        tertiaryColor: '#f8fafc',
        tertiaryTextColor: '#0f172a',
        tertiaryBorderColor: '#cbd5e1',
        lineColor: '#64748b',
        textColor: '#0f172a',
        edgeLabelBackground: '#ffffff',
        clusterBkg: 'rgba(148,163,184,.10)',
        clusterBorder: '#cbd5e1',
        fontFamily: 'ui-sans-serif, system-ui, -apple-system, "Segoe UI", Roboto, sans-serif',
        fontSize: '14px',
      }
}

async function render() {
  try {
    const mermaid = (await import('mermaid')).default
    mermaid.initialize({
      startOnLoad: false,
      securityLevel: 'strict',
      theme: 'base',
      themeVariables: themeVariables(isDark.value),
      flowchart: { htmlLabels: true, curve: 'basis', useMaxWidth: true },
    })
    const { svg: markup } = await mermaid.render(`sf-mermaid-${seq++}`, decode(props.code))
    svg.value = markup
    failed.value = false
  } catch (error) {
    failed.value = true
    console.error('sci-form docs: could not render a Mermaid diagram', error)
  }
}

onMounted(render)
watch(isDark, render)
</script>

<template>
  <div class="sf-mermaid" :class="{ 'sf-mermaid-failed': failed }">
    <div v-if="svg" class="sf-mermaid-canvas" v-html="svg" />
    <pre v-else class="sf-mermaid-source">{{ decode(code) }}</pre>
  </div>
</template>
