import { defineConfig } from 'vite';
import { fileURLToPath } from 'node:url';

const wasmRoot = fileURLToPath(new URL('../../', import.meta.url));

export default defineConfig({
  worker: {
    format: 'es',
  },
  server: {
    headers: {
      'Cross-Origin-Opener-Policy': 'same-origin',
      'Cross-Origin-Embedder-Policy': 'require-corp',
    },
    fs: {
      allow: [wasmRoot],
    },
  },
  preview: {
    headers: {
      'Cross-Origin-Opener-Policy': 'same-origin',
      'Cross-Origin-Embedder-Policy': 'require-corp',
    },
  },
});