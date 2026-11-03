<script setup lang="ts">
// One live example: shows the code, runs it in the visitor's browser against
// the compiled sci-form WASM package and prints the result.
import { computed, onMounted, ref } from 'vue'
import { withBase } from 'vitepress'
import catalog from '../../examples/catalog.json'
import { makeScope, runExample } from '../../examples/runtime.mjs'

const props = defineProps<{ id: string; open?: boolean }>()
const ex = computed(() => (catalog as any[]).find((e) => e.id === props.id))

type State = 'idle' | 'loading' | 'done' | 'error'
const state = ref<State>('idle')
const output = ref('')
const ms = ref(0)
const full = ref(false)
const MAX_LINES = 40

let scopePromise: Promise<any> | null = null
function getScope() {
  if (!scopePromise) {
    scopePromise = (async () => {
      const url = withBase('/wasm/sci_form_wasm.js')
      const mod: any = await import(/* @vite-ignore */ url)
      await mod.default(withBase('/wasm/sci_form_wasm_bg.wasm'))
      return makeScope(mod)
    })()
  }
  return scopePromise
}

const shown = computed(() => {
  if (full.value) return output.value
  const lines = output.value.split('\n')
  return lines.length > MAX_LINES
    ? lines.slice(0, MAX_LINES).join('\n') + `\n… ${lines.length - MAX_LINES} more lines`
    : output.value
})
const truncated = computed(() => output.value.split('\n').length > MAX_LINES)

async function run() {
  if (!ex.value) return
  state.value = 'loading'
  full.value = false
  try {
    const scope = await getScope()
    await new Promise((r) => setTimeout(r, 0)) // let the UI paint "running"
    const r = await runExample(scope, ex.value.code)
    const v: any = r.value
    if (v && typeof v === 'object' && !Array.isArray(v) && v.error) throw new Error(String(v.error))
    output.value = typeof v === 'string' ? v : JSON.stringify(v, null, 2)
    ms.value = r.ms
    state.value = 'done'
  } catch (e: any) {
    output.value = String(e?.message ?? e)
    state.value = 'error'
  }
}

onMounted(() => {
  if (props.open) run()
})
</script>

<template>
  <div v-if="ex" class="live">
    <div class="live-head">
      <strong>{{ ex.title }}</strong>
      <code class="live-id">{{ ex.id }}</code>
    </div>
    <p class="live-about">{{ ex.about }}</p>
    <pre class="live-code"><code>{{ ex.code }}</code></pre>
    <div class="live-bar">
      <button class="live-run" :disabled="state === 'loading'" @click="run">
        {{ state === 'loading' ? 'Running…' : state === 'done' ? '▶ Run again' : '▶ Run in your browser' }}
      </button>
      <span v-if="state === 'done'" class="live-ms">{{ ms.toFixed(0) }} ms</span>
    </div>
    <pre v-if="state === 'done'" class="live-out"><code>{{ shown }}</code></pre>
    <button v-if="state === 'done' && truncated" class="live-more" @click="full = !full">
      {{ full ? 'Show less' : 'Show all' }}
    </button>
    <pre v-if="state === 'error'" class="live-err"><code>{{ output }}</code></pre>
    <p v-if="state === 'error' && /import|fetch|wasm|Failed to load/i.test(output)" class="live-note">
      The compiled WebAssembly package is only present on the published site (it is built together with each release).
    </p>
  </div>
  <p v-else><em>Unknown example “{{ id }}”.</em></p>
</template>

<style scoped>
.live { border: 1px solid var(--vp-c-divider); border-radius: 10px; padding: 14px 16px; margin: 16px 0; background: var(--vp-c-bg-soft); }
.live-head { display: flex; gap: 10px; align-items: baseline; flex-wrap: wrap; }
.live-id { font-size: 12px; color: var(--vp-c-text-2); }
.live-about { margin: 6px 0 10px; color: var(--vp-c-text-2); }
.live-code, .live-out, .live-err { margin: 8px 0; padding: 10px 12px; border-radius: 8px; overflow-x: auto; font-size: 12.5px; line-height: 1.5; background: var(--vp-c-bg); border: 1px solid var(--vp-c-divider); }
.live-out { max-height: 420px; overflow: auto; }
.live-err { border-color: var(--vp-c-danger-1); color: var(--vp-c-danger-1); }
.live-bar { display: flex; gap: 12px; align-items: center; }
.live-run, .live-more { font: inherit; font-size: 13px; padding: 5px 12px; border-radius: 8px; border: 1px solid var(--vp-c-brand-1); background: var(--vp-c-brand-soft); color: var(--vp-c-brand-1); cursor: pointer; }
.live-run:disabled { opacity: .6; cursor: progress; }
.live-ms { font-size: 12px; color: var(--vp-c-text-2); }
.live-note { font-size: 12.5px; color: var(--vp-c-text-2); }
</style>
