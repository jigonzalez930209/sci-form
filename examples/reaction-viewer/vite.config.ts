import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'

export default defineConfig({
  plugins: [react()],
  base: '/sci-form/reaction-viewer/',
  build: {
    outDir: '../../docs/public/reaction-viewer',
    emptyOutDir: true,
  },
})
