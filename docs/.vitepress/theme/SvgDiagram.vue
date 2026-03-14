<script setup lang="ts">
import { ref, onMounted } from 'vue'

const props = defineProps<{
  src: string
  alt?: string
}>()

const svgHtml = ref('')
const loaded = ref(false)

onMounted(async () => {
  try {
    // Handle VitePress base path (e.g. /sci-form/ in production)
    const base = import.meta.env.BASE_URL.replace(/\/$/, '')
    const url = base + props.src
    const res = await fetch(url)
    if (!res.ok) throw new Error(`HTTP ${res.status}`)
    const text = await res.text()
    // Strip inline style from <svg> root — our wrapper CSS controls sizing
    svgHtml.value = text.replace(
      /(<svg\b(?:[^>](?!style=))*)\s+style="[^"]*"/,
      '$1'
    )
    loaded.value = true
  } catch (e) {
    console.warn('[SvgDiagram] failed to load', props.src, e)
  }
})
</script>

<template>
  <figure class="svg-diagram-wrap" :aria-label="alt">
    <!-- Skeleton while loading -->
    <div v-if="!loaded" class="svg-skeleton" role="img" :aria-label="alt" />
    <!-- Inline SVG: scales perfectly, no img aspect-ratio bug -->
    <div
      v-else
      class="svg-inline"
      role="img"
      :aria-label="alt"
      v-html="svgHtml"
    />
  </figure>
</template>
