# Vite WebGPU Hybrid Example

This example demonstrates the browser path with both WebGPU and threaded CPU enabled.

## 1. Build the WASM package

From [crates/wasm/build.sh](/home/lestad/github/sci-form/crates/wasm/build.sh):

```bash
cd /home/lestad/github/sci-form/crates/wasm
./build.sh --web-only --web-features "parallel experimental-gpu"
```

That creates [crates/wasm/pkg](/home/lestad/github/sci-form/crates/wasm/pkg), which this Vite app imports via `../../../pkg/` from [src/main.js](src/main.js).

## 2. Run the Vite app

```bash
cd /home/lestad/github/sci-form/crates/wasm/examples/vite-webgpu-hybrid
npm install
npm run dev
```

The Vite config already sets:

- `Cross-Origin-Opener-Policy: same-origin`
- `Cross-Origin-Embedder-Policy: require-corp`

Those headers are required for `SharedArrayBuffer`, which in turn is required by `initThreadPool()`.

## 3. What it exercises

- `initThreadPool()` for browser CPU parallelism
- `init_webgpu()` for WebGPU device initialization
- `compute_esp_grid_accelerated(..., "hybrid")`
- `eht_orbital_grid_accelerated(..., "gpu")`
- `eht_density_grid_accelerated(..., "auto")`