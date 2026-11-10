# GPU backend portability in sci-form

This page summarizes how portable the current GPU backend is, which systems it
can run on, and what its real limits are.

## Short answer

The backend is designed to be cross-platform and multi-vendor. It does not
depend on CUDA, ROCm or a particular GPU brand: it uses `wgpu`, which abstracts
the operating system's native backend and selects a compatible adapter
automatically.

In principle, that means it can run on:

- macOS with Metal-compatible GPUs
- Windows with Direct3D 12 or Vulkan-compatible GPUs and drivers
- Linux with Vulkan-compatible GPUs and drivers

It can run on hardware from:

- NVIDIA
- AMD
- Intel
- Apple Silicon

That does not mean every card works without exception. Real compatibility
depends on four things:

- the crate is compiled with the `experimental-gpu` feature;
- the operating system exposes a backend supported by `wgpu`;
- the GPU driver is installed correctly and is recent enough;
- the adapter supports the compute and buffer limits these kernels need.

## How backend selection works today

The GPU context uses `GpuContext::best_available()` and tries to create a native
`wgpu` device.

The current flow is:

1. If `experimental-gpu` is not enabled, the code never attempts a GPU.
2. If the feature is enabled, `wgpu` asks for a high-performance adapter.
3. If an adapter is found and the device can be created, the GPU is used.
4. If no adapter is found or initialization fails, the code falls back to the
   CPU.

In other words, the behaviour is not "GPU or error"; it is "GPU when available,
CPU otherwise".

## By operating system

### macOS

On macOS, `wgpu` normally uses Metal. That includes:

- Macs with Apple Silicon;
- Intel Macs with Metal-compatible GPUs.

The GPU does not need to be NVIDIA, AMD or Intel in particular. What matters is
working Metal support.

### Windows

On Windows, the usual path is Direct3D 12; some systems can also use Vulkan
depending on the environment and drivers.

In practice this works on NVIDIA, AMD and Intel hardware as long as the driver
supports the backend `wgpu` ends up using.

### Linux

On Linux, the important path for this project is Vulkan. Again, the GPU brand is
not the deciding factor; the graphics stack is:

- Mesa or a proprietary driver;
- Vulkan support;
- device and permission visibility from the process environment.

## What "cross-platform" means here

It means the shader code and the compute pipeline are written against a
portable abstraction, not against a single vendor's proprietary API. The same
kernels can be reused on different machines without maintaining a CUDA version
for NVIDIA, a Metal version for Apple and a DirectX version for Windows.

## What it does not mean

It does not mean that any old GPU, any virtual GPU or any broken driver will
work. It also does not mean that every adapter gives the same performance. For
example:

- an Intel integrated GPU can run the kernel but be slower than a discrete GPU;
- an Apple GPU can run the kernel well, with a different performance profile
  from an NVIDIA card;
- an old driver can accept the adapter but fail while creating buffers or
  compiling complex shaders.

## Current state in sci-form

The native GPU backend is designed for desktop. Initialization is conditional
on:

- the `experimental-gpu` feature;
- a native target, not `wasm32`.

That means this path is not enabled today for WebAssembly in the browser with
the same infrastructure.

## Current state in WASM / browser

There is now a separate WebGPU path for WASM, distinct from the synchronous
native backend. This matters because the browser cannot force the native
runtime model: WebGPU in WASM is asynchronous by nature.

WebGPU support in WASM is therefore a dedicated layer that leaves the existing
native backend almost intact.

### What the first WASM layer supports

- asynchronous WebGPU initialization;
- WebGPU runtime status queries;
- accelerated execution in WASM for experimental Tier 3 paths;
- a mode selector for `cpu`, `gpu`, `hybrid` and `auto`.

In this iteration the accelerated WASM path is connected mainly to calculations
where the O(N²) pattern was already clear:

- EEQ
- D4
- CPM

### What "hybrid" means in WASM

In the browser, "hybrid" does not mean exactly the same thing as in a native
desktop binary. It means taking advantage of the fact that part of an algorithm
runs well on the GPU while another part is still better on the CPU, or that the
CPU can prepare and post-process data while the GPU runs the dominant kernel.

Examples:

- **EEQ:** the GPU builds the damped interaction matrix; the CPU assembles and
  solves the linear system.
- **D4:** the GPU handles the dominant O(N²) term; the CPU keeps the auxiliary
  and corrective parts.
- **CPM:** the GPU builds the Coulomb matrix; the CPU keeps the charge
  iteration.

This is cooperative hybrid execution, not blind duplication of the full
calculation on both sides.

## Real limits of WebGPU in WASM

The WASM/WebGPU support has practical limits that should be explicit:

- it requires a browser with WebGPU available;
- it normally requires a secure context, meaning HTTPS or `localhost`;
- the exposed API is asynchronous;
- performance depends heavily on the browser, driver and operating system;
- some native-backend kernels are not yet exposed in the WASM WebGPU layer.

## Parallel CPU in WASM

WASM can also accelerate the CPU with `wasm-bindgen-rayon`, but there is another
important restriction. Compiling the `parallel` feature is not enough: the web
environment must provide threads and atomics.

In practice that means:

- `SharedArrayBuffer` is available;
- the browser has an appropriate origin isolation policy;
- the WASM target/configuration supports `atomics` and `bulk-memory`.

Without that configuration, the parallel WASM build fails or the runtime cannot
use the thread pool.

## Reproducible path for a complete hybrid browser build

The recommended reproducible path in this repository is not `cargo check` from
the root to validate browser threads. Use the canonical `crates/wasm` package.

### Step 1: build the web package with the official script

```bash
cd crates/wasm
./build.sh --web-only --web-features "parallel experimental-gpu"
```

That script forces the `web` target, keeps the `snippets` directory for
`wasm-bindgen-rayon`, and adds the WASM target features needed when `parallel`
is active.

### Step 2: serve the frontend with COOP and COEP headers

Minimum headers:

- `Cross-Origin-Opener-Policy: same-origin`
- `Cross-Origin-Embedder-Policy: require-corp`

### Step 3: initialize both runtimes in the browser

```ts
await init();
await initThreadPool(navigator.hardwareConcurrency ?? 4);
await init_webgpu();
```

That closes the hybrid path end to end:

- parallel CPU in the browser through `wasm-bindgen-rayon`;
- GPU through WebGPU;
- `hybrid` mode to split work between both for volumetric kernels and
  experimental paths.

**About the browser example:** there is no minimal example in the repository.
The `crates/wasm/examples/vite-webgpu-hybrid/` folder that earlier versions of
this document linked **does not exist in git**, and the npm scripts
`webgpu:benchmark:dev`, `webgpu:benchmark:build` and `webgpu:benchmark:preview`
point at it, so they fail. The same applies to `npm run webgpu:test`, which
needs that folder to build the harness that Playwright drives.

What does exist and is reproducible:

- the WASM package with the `experimental-gpu` feature compiled (step 1 above);
- the `init_webgpu()` and `webgpu_status()` exports;
- the accelerated surfaces that accept `execution_mode`
  (`compute_esp_grid_accelerated`, `eht_orbital_grid_accelerated`,
  `eht_density_grid_accelerated`, `compute_orbital_mesh_accelerated`,
  `compute_eeq_charges_accelerated`, `compute_d4_energy_accelerated`,
  `compute_cpm_charges_accelerated`);
- the report comparator `tests/integration/compare_webgpu_reports.js` with its
  policy in `tests/fixtures/webgpu_regression_policy.json`.

The numerical validation of the WebGPU path that runs in CI today is the Rust
side (GPU↔CPU gates in `src/gpu/`), not a browser run. Details are in
[WebGPU Validation](../algorithm/webgpu-validation.md).

## Current practical summary

- **Native desktop:** main GPU backend, synchronous, with CPU fallback.
- **WASM:** dedicated WebGPU backend, asynchronous, with CPU fallback.
- **Hybrid in WASM:** cooperative and limited by the browser environment.
- **Parallel CPU in WASM:** possible, but depends on atomics, bulk-memory and
  the web environment configuration.

## Practical consequence for users

The correct answer to "can it run on any graphics adapter?" is:

It can run on a wide variety of adapters and brands, but there is no universal
per-model guarantee. Real compatibility depends on the system graphics backend,
the driver and the adapter limits.

The correct answer to "does it work on NVIDIA, Intel or AMD?" is:

Yes, in principle. The design is not tied to one brand and can work with NVIDIA,
Intel, AMD and Apple hardware as long as the system exposes a backend supported
by `wgpu` and device initialization succeeds.

## Safety behaviour

If the GPU is unavailable, `sci-form` should not become unusable. The current
design falls back to the CPU when:

- the `experimental-gpu` feature is missing;
- no adapter is found;
- GPU device creation fails.

That reduces operational risk: the GPU accelerates work, but it is not required
for the code to keep running.

## Operational recommendation

For production or serious validation, assume this rule:

Validate by operating system + driver + GPU family, not just by brand.

A reasonable matrix is:

- macOS + Apple Silicon
- Windows + NVIDIA
- Windows + Intel
- Linux + AMD or NVIDIA with Vulkan

That verifies the real portability of the backend, not just the code's
cross-platform intention.
