<script setup lang="ts">
// The whole catalogue: one runnable example per public function, by topic.
import { computed, ref } from 'vue'
import catalog from '../../examples/catalog.json'
import LiveExample from './LiveExample.vue'

const query = ref('')
const groups = computed(() => {
  const q = query.value.trim().toLowerCase()
  const out: { name: string; items: any[] }[] = []
  for (const e of catalog as any[]) {
    if (q && !(e.id + ' ' + e.title + ' ' + e.about).toLowerCase().includes(q)) continue
    let g = out.find((x) => x.name === e.group)
    if (!g) out.push((g = { name: e.group, items: [] }))
    g.items.push(e)
  }
  return out
})
const total = (catalog as any[]).length
</script>

<template>
  <div class="le">
    <input v-model="query" class="le-search" type="search" :placeholder="`Search ${total} examples (e.g. pm3, nmr, solvation)…`" />
    <nav class="le-index">
      <a v-for="g in groups" :key="g.name" :href="'#' + g.name.replace(/[^a-z0-9]+/gi, '-').toLowerCase()">{{ g.name }} ({{ g.items.length }})</a>
    </nav>
    <section v-for="g in groups" :key="g.name">
      <h2 :id="g.name.replace(/[^a-z0-9]+/gi, '-').toLowerCase()">{{ g.name }}</h2>
      <details v-for="e in g.items" :key="e.id" class="le-item">
        <summary><strong>{{ e.title }}</strong> <code>{{ e.id }}</code></summary>
        <LiveExample :id="e.id" />
      </details>
    </section>
    <p v-if="!groups.length"><em>No example matches “{{ query }}”.</em></p>
  </div>
</template>

<style scoped>
.le-search { width: 100%; padding: 9px 12px; border-radius: 8px; border: 1px solid var(--vp-c-divider); background: var(--vp-c-bg-soft); color: var(--vp-c-text-1); font: inherit; }
.le-index { display: flex; flex-wrap: wrap; gap: 6px 14px; margin: 12px 0 4px; font-size: 13px; }
.le-item { margin: 6px 0; }
.le-item summary { cursor: pointer; padding: 6px 2px; }
.le-item summary code { font-size: 12px; color: var(--vp-c-text-2); margin-left: 8px; }
</style>
