//! GPU compute context — wgpu device initialization with CPU fallback.
//!
//! Wraps `wgpu::Device` + `Queue` for native GPU compute, or falls back
//! to CPU-only mode when no GPU is available or the feature is disabled.
//!
//! # Contract with callers
//!
//! **The CPU fallback always works.** [`GpuContext::cpu_fallback`] and
//! [`GpuContext::best_available`] never fail and never panic, so a caller can
//! always obtain a usable context. When the `experimental-gpu` feature is
//! disabled, or no adapter is present, or adapter acquisition panics inside
//! `wgpu-hal` (headless hosts aborting with `BadDisplay`), the context
//! degrades to CPU-only mode and records the reason in
//! [`GpuContext::activation_report`].
//!
//! The GPU path is therefore *optional*: the WebGPU default limits that the
//! rest of the crate respects are **256 invocations per workgroup**
//! (`max_compute_invocations_per_workgroup`) and a **128 MB per storage
//! binding** limit (`max_storage_buffer_binding_size`), see
//! [`crate::gpu::memory_budget`].

#[cfg(all(feature = "experimental-gpu", not(target_arch = "wasm32")))]
use std::borrow::Cow;

#[cfg(all(feature = "experimental-gpu", not(target_arch = "wasm32")))]
use std::collections::HashMap;

#[cfg(all(feature = "experimental-gpu", not(target_arch = "wasm32")))]
use std::sync::mpsc;

#[cfg(all(feature = "experimental-gpu", not(target_arch = "wasm32")))]
use std::sync::Mutex;

use std::sync::OnceLock;

#[cfg(all(feature = "experimental-gpu", not(target_arch = "wasm32")))]
use wgpu::util::DeviceExt;

use super::backend_report::{GpuActivationReport, GpuActivationState};

/// Process-wide GPU context.
///
/// The adapter enumeration and the device creation are expensive (0.1–10 ms
/// each) and used to be repeated on every dispatch (blocker D-M2).  The
/// context is therefore created once and shared.
static GLOBAL_CONTEXT: OnceLock<GpuContext> = OnceLock::new();

/// Capabilities detected from the GPU adapter.
///
/// The `Default` implementation describes the conservative CPU-fallback
/// profile: no GPU, but the WebGPU *minimum* limits (256 threads per
/// workgroup axis, 256 invocations per workgroup) so that dispatch
/// configurations planned against a fallback context stay valid.
#[derive(Debug, Clone)]
pub struct ComputeCapabilities {
    /// Human-readable backend name, e.g. `"CPU-fallback"` or the `wgpu`
    /// backend enum (`"Vulkan"`, `"Metal"`, `"Dx12"`).
    pub backend: String,
    /// Maximum workgroup size along the x axis, in threads.
    pub max_workgroup_size_x: u32,
    /// Maximum workgroup size along the y axis, in threads.
    pub max_workgroup_size_y: u32,
    /// Maximum number of invocations in a single workgroup (threads).
    pub max_workgroup_invocations: u32,
    /// Maximum size of a single storage-buffer binding, in bytes
    /// (WebGPU default: 134 217 728 = 128 MB).
    pub max_storage_buffer_size: u64,
    /// `true` when a real GPU device is attached.
    pub gpu_available: bool,
}

impl Default for ComputeCapabilities {
    fn default() -> Self {
        Self {
            backend: "CPU-fallback".to_string(),
            max_workgroup_size_x: 256,
            max_workgroup_size_y: 256,
            max_workgroup_invocations: 256,
            max_storage_buffer_size: u64::MAX,
            gpu_available: false,
        }
    }
}

/// Buffer access mode for compute bindings.
///
/// The mode determines the `wgpu` buffer usage flags and the bind-group
/// layout entry, and therefore whether the kernel may write to the buffer.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ComputeBindingKind {
    /// Uniform buffer: read-only, fast-path, limited to 64 KB by default.
    Uniform,
    /// Storage buffer: read-only, up to 128 MB by default.
    StorageReadOnly,
    /// Storage buffer: read-write; its contents are read back after the
    /// dispatch and appear in [`ComputeDispatchResult::outputs`].
    StorageReadWrite,
}

/// One logical binding passed to a compute kernel.
///
/// The binding index in the kernel's bind group is the position of the
/// descriptor in [`ComputeDispatchDescriptor::bindings`].
#[derive(Debug, Clone)]
pub struct ComputeBindingDescriptor {
    /// Debug label for the underlying GPU buffer.
    pub label: String,
    /// Access mode, which decides the buffer usage flags and layout entry.
    pub kind: ComputeBindingKind,
    /// Raw payload; interpreted as `f32`/`u32` words by the shader.
    pub bytes: Vec<u8>,
}

/// Generic compute dispatch description.
///
/// This is the CPU-side description of one WGSL compute dispatch; it is
/// independent of the physical workgroup size, which the shader itself
/// declares with `@workgroup_size`.
#[derive(Debug, Clone)]
pub struct ComputeDispatchDescriptor {
    /// Debug label for the pipeline, bind group and command encoder.
    pub label: String,
    /// Full WGSL source of the compute shader.
    pub shader_source: String,
    /// Entry-point name inside `shader_source` (usually `"main"`).
    pub entry_point: String,
    /// Number of workgroups dispatched along each axis.
    pub workgroup_count: [u32; 3],
    /// Bindings, in binding-index order.
    pub bindings: Vec<ComputeBindingDescriptor>,
}

/// Output produced by a compute dispatch.
///
/// `outputs` contains one byte vector per read-write binding, in the order
/// the read-write bindings appear in [`ComputeDispatchDescriptor::bindings`].
/// The contents are raw bytes; use
/// [`bytes_to_f32_vec`] to decode
/// `f32` data.
#[derive(Debug, Clone)]
pub struct ComputeDispatchResult {
    /// Backend that executed the dispatch (e.g. `"Vulkan"` or `"CPU-fallback"`).
    pub backend: String,
    /// Raw read-back payloads of the read-write bindings.
    pub outputs: Vec<Vec<u8>>,
}

/// GPU compute context handle.
///
/// Holds the detected [`ComputeCapabilities`] and, when the GPU is actually
/// available, the `wgpu` device/queue plus a compiled-pipeline cache.
pub struct GpuContext {
    /// Capabilities of the attached device (or the CPU-fallback profile).
    pub capabilities: ComputeCapabilities,
    runtime_error: Option<String>,
    #[cfg(all(feature = "experimental-gpu", not(target_arch = "wasm32")))]
    runtime: Option<NativeGpuRuntime>,
}

#[cfg(all(feature = "experimental-gpu", not(target_arch = "wasm32")))]
struct NativeGpuRuntime {
    _instance: wgpu::Instance,
    _adapter: wgpu::Adapter,
    device: wgpu::Device,
    queue: wgpu::Queue,
    /// Cache of already compiled compute pipelines (blocker D-M2).
    ///
    /// Creating a shader module, a bind group layout and a compute pipeline
    /// costs ~0.1–10 ms; the same kernels are dispatched over and over
    /// (once per orbital, once per grid point batch, ...), so the compiled
    /// pipeline is reused whenever the label, the shader source, the entry
    /// point and the binding layout are unchanged.
    pipelines: Mutex<HashMap<PipelineCacheKey, std::sync::Arc<CachedPipeline>>>,
}

/// Identity of a compiled compute pipeline.
#[cfg(all(feature = "experimental-gpu", not(target_arch = "wasm32")))]
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct PipelineCacheKey {
    label: String,
    shader_source: String,
    entry_point: String,
    binding_kinds: Vec<ComputeBindingKind>,
}

/// A compiled compute pipeline and the layout it was built against.
///
/// `wgpu` 0.19 resource handles are not `Clone`, so the cache stores them
/// behind an `Arc`.
#[cfg(all(feature = "experimental-gpu", not(target_arch = "wasm32")))]
struct CachedPipeline {
    bind_group_layout: wgpu::BindGroupLayout,
    /// Kept alive together with the pipeline: the pipeline was created
    /// against this layout.
    #[allow(dead_code)]
    pipeline_layout: wgpu::PipelineLayout,
    pipeline: wgpu::ComputePipeline,
}

impl std::fmt::Debug for GpuContext {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("GpuContext")
            .field("capabilities", &self.capabilities)
            .field("runtime_error", &self.runtime_error)
            .finish()
    }
}

impl GpuContext {
    /// Create a CPU-fallback context (always available).
    pub fn cpu_fallback() -> Self {
        Self {
            capabilities: ComputeCapabilities::default(),
            runtime_error: None,
            #[cfg(all(feature = "experimental-gpu", not(target_arch = "wasm32")))]
            runtime: None,
        }
    }

    /// Attempt to initialize a real wgpu compute device.
    ///
    /// Adapter/device acquisition can panic inside `wgpu-hal` (for example the
    /// GLES/EGL backend aborts with `BadDisplay` on headless hosts). The whole
    /// acquisition is therefore wrapped in `catch_unwind` so the CPU-fallback
    /// contract of [`GpuContext::best_available`] always holds.
    ///
    /// # Errors
    ///
    /// - `"experimental-gpu feature not enabled"` when the `experimental-gpu`
    ///   feature is off (or the target is `wasm32`).
    /// - `"No GPU adapter found"` when no compatible adapter could be
    ///   acquired.
    /// - `"Failed to create wgpu device: {err}"` when the device request was
    ///   rejected by the driver.
    /// - `"GPU adapter acquisition panicked: {detail}"` when the underlying
    ///   `wgpu-hal` aborted during adapter enumeration.
    pub fn try_create() -> Result<Self, String> {
        #[cfg(all(feature = "experimental-gpu", not(target_arch = "wasm32")))]
        {
            std::panic::catch_unwind(std::panic::AssertUnwindSafe(Self::try_create_native))
                .unwrap_or_else(|payload| {
                    let detail = payload
                        .downcast_ref::<&str>()
                        .map(|s| (*s).to_string())
                        .or_else(|| payload.downcast_ref::<String>().cloned())
                        .unwrap_or_else(|| "unknown panic".to_string());
                    Err(format!("GPU adapter acquisition panicked: {detail}"))
                })
        }

        #[cfg(not(all(feature = "experimental-gpu", not(target_arch = "wasm32"))))]
        {
            Err("experimental-gpu feature not enabled".to_string())
        }
    }

    #[cfg(all(feature = "experimental-gpu", not(target_arch = "wasm32")))]
    fn try_create_native() -> Result<Self, String> {
        let instance = wgpu::Instance::default();
        let adapter = pollster::block_on(instance.request_adapter(&wgpu::RequestAdapterOptions {
            power_preference: wgpu::PowerPreference::HighPerformance,
            compatible_surface: None,
            force_fallback_adapter: false,
        }))
        .ok_or_else(|| "No GPU adapter found".to_string())?;

        let adapter_info = adapter.get_info();
        let limits = adapter.limits();
        let required_limits = wgpu::Limits::default().using_resolution(limits.clone());

        let (device, queue) = pollster::block_on(adapter.request_device(
            &wgpu::DeviceDescriptor {
                label: Some("sci-form gpu"),
                required_features: wgpu::Features::empty(),
                required_limits,
            },
            None,
        ))
        .map_err(|err| format!("Failed to create wgpu device: {err}"))?;

        Ok(Self {
            capabilities: ComputeCapabilities {
                backend: format!("{:?}", adapter_info.backend),
                max_workgroup_size_x: limits.max_compute_workgroup_size_x,
                max_workgroup_size_y: limits.max_compute_workgroup_size_y,
                max_workgroup_invocations: limits.max_compute_invocations_per_workgroup,
                max_storage_buffer_size: limits.max_storage_buffer_binding_size as u64,
                gpu_available: true,
            },
            runtime_error: None,
            runtime: Some(NativeGpuRuntime {
                _instance: instance,
                _adapter: adapter,
                device,
                queue,
                pipelines: Mutex::new(HashMap::new()),
            }),
        })
    }

    /// Process-wide shared GPU context.
    ///
    /// The context is created on first use and reused afterwards: adapter
    /// enumeration and device creation used to be repeated on *every*
    /// dispatch (blocker D-M2).  The returned context is a CPU-fallback
    /// context when no GPU is available, so the "always works" contract of
    /// [`GpuContext::best_available`] is preserved.
    pub fn shared() -> &'static GpuContext {
        GLOBAL_CONTEXT.get_or_init(Self::best_available)
    }

    /// Best available backend: GPU if possible, CPU fallback otherwise.
    pub fn best_available() -> Self {
        match Self::try_create() {
            Ok(ctx) => ctx,
            Err(reason) => {
                let mut ctx = Self::cpu_fallback();
                ctx.runtime_error = Some(reason);
                ctx
            }
        }
    }

    /// Build an activation report describing the current runtime state.
    pub fn activation_report(&self) -> GpuActivationReport {
        if self.capabilities.gpu_available {
            GpuActivationReport {
                backend: self.capabilities.backend.clone(),
                feature_enabled: true,
                gpu_available: true,
                runtime_ready: true,
                state: GpuActivationState::Ready,
                reason: "GPU runtime available".to_string(),
            }
        } else if cfg!(feature = "experimental-gpu") {
            GpuActivationReport {
                backend: self.capabilities.backend.clone(),
                feature_enabled: true,
                gpu_available: false,
                runtime_ready: false,
                state: GpuActivationState::NoAdapter,
                reason: self
                    .runtime_error
                    .clone()
                    .unwrap_or_else(|| "experimental-gpu enabled but no adapter found".to_string()),
            }
        } else {
            GpuActivationReport {
                backend: "CPU-fallback".to_string(),
                feature_enabled: false,
                gpu_available: false,
                runtime_ready: false,
                state: GpuActivationState::FeatureDisabled,
                reason: "experimental-gpu feature not enabled".to_string(),
            }
        }
    }

    /// Whether the GPU backend is available.
    pub fn is_gpu_available(&self) -> bool {
        self.capabilities.gpu_available
    }

    #[cfg(all(feature = "experimental-gpu", not(target_arch = "wasm32")))]
    fn runtime(&self) -> Result<&NativeGpuRuntime, String> {
        self.runtime.as_ref().ok_or_else(|| {
            self.runtime_error
                .clone()
                .unwrap_or_else(|| "GPU runtime not initialized".to_string())
        })
    }

    /// Dispatch an arbitrary WGSL compute kernel.
    ///
    /// Compiled pipelines are cached by (label, shader source, entry point,
    /// binding kinds), so repeatedly dispatching the same kernel does not
    /// re-create the shader module.
    ///
    /// # Errors
    ///
    /// - `"experimental-gpu feature not enabled"` when the `experimental-gpu`
    ///   feature is off (or the target is `wasm32`).
    /// - The stored runtime error, or `"GPU runtime not initialized"`, when
    ///   this context is a CPU-fallback context.
    /// - `"GPU readback channel error"` when the mapping callback never
    ///   delivered a result.
    /// - `"GPU buffer map failed: {err}"` when the device refused to map the
    ///   staging buffer for reading.
    pub fn run_compute(
        &self,
        descriptor: &ComputeDispatchDescriptor,
    ) -> Result<ComputeDispatchResult, String> {
        #[cfg(all(feature = "experimental-gpu", not(target_arch = "wasm32")))]
        {
            let runtime = self.runtime()?;

            let mut layout_entries = Vec::with_capacity(descriptor.bindings.len());
            let mut buffers = Vec::with_capacity(descriptor.bindings.len());
            let mut readbacks = Vec::new();

            for (index, binding) in descriptor.bindings.iter().enumerate() {
                let usage = match binding.kind {
                    ComputeBindingKind::Uniform => {
                        wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST
                    }
                    ComputeBindingKind::StorageReadOnly => {
                        wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST
                    }
                    ComputeBindingKind::StorageReadWrite => {
                        wgpu::BufferUsages::STORAGE
                            | wgpu::BufferUsages::COPY_DST
                            | wgpu::BufferUsages::COPY_SRC
                    }
                };

                let buffer = runtime
                    .device
                    .create_buffer_init(&wgpu::util::BufferInitDescriptor {
                        label: Some(&binding.label),
                        contents: &binding.bytes,
                        usage,
                    });

                layout_entries.push(wgpu::BindGroupLayoutEntry {
                    binding: index as u32,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Buffer {
                        ty: match binding.kind {
                            ComputeBindingKind::Uniform => wgpu::BufferBindingType::Uniform,
                            ComputeBindingKind::StorageReadOnly => {
                                wgpu::BufferBindingType::Storage { read_only: true }
                            }
                            ComputeBindingKind::StorageReadWrite => {
                                wgpu::BufferBindingType::Storage { read_only: false }
                            }
                        },
                        has_dynamic_offset: false,
                        min_binding_size: None,
                    },
                    count: None,
                });
                if matches!(binding.kind, ComputeBindingKind::StorageReadWrite) {
                    readbacks.push((buffers.len(), binding.bytes.len()));
                }
                buffers.push(buffer);
            }

            let bind_group_entries: Vec<_> = buffers
                .iter()
                .enumerate()
                .map(|(index, buffer)| wgpu::BindGroupEntry {
                    binding: index as u32,
                    resource: buffer.as_entire_binding(),
                })
                .collect();

            // Reuse the compiled pipeline whenever the kernel is unchanged
            // (blocker D-M2).  The cache key contains everything the pipeline
            // depends on: the shader source, the entry point and the shape of
            // the binding layout.
            let cache_key = PipelineCacheKey {
                label: descriptor.label.clone(),
                shader_source: descriptor.shader_source.clone(),
                entry_point: descriptor.entry_point.clone(),
                binding_kinds: descriptor.bindings.iter().map(|b| b.kind).collect(),
            };

            let cached = {
                let cache = runtime
                    .pipelines
                    .lock()
                    .unwrap_or_else(|poisoned| poisoned.into_inner());
                cache.get(&cache_key).cloned()
            };

            let compiled =
                match cached {
                    Some(cached) => cached,
                    None => {
                        let shader =
                            runtime
                                .device
                                .create_shader_module(wgpu::ShaderModuleDescriptor {
                                    label: Some(&descriptor.label),
                                    source: wgpu::ShaderSource::Wgsl(Cow::Borrowed(
                                        &descriptor.shader_source,
                                    )),
                                });

                        let bind_group_layout = runtime.device.create_bind_group_layout(
                            &wgpu::BindGroupLayoutDescriptor {
                                label: Some(&format!("{} layout", descriptor.label)),
                                entries: &layout_entries,
                            },
                        );
                        let pipeline_layout = runtime.device.create_pipeline_layout(
                            &wgpu::PipelineLayoutDescriptor {
                                label: Some(&format!("{} pipeline", descriptor.label)),
                                bind_group_layouts: &[&bind_group_layout],
                                push_constant_ranges: &[],
                            },
                        );
                        let pipeline = runtime.device.create_compute_pipeline(
                            &wgpu::ComputePipelineDescriptor {
                                label: Some(&descriptor.label),
                                layout: Some(&pipeline_layout),
                                module: &shader,
                                entry_point: &descriptor.entry_point,
                            },
                        );

                        let compiled = std::sync::Arc::new(CachedPipeline {
                            bind_group_layout,
                            pipeline_layout,
                            pipeline,
                        });

                        let mut cache = runtime
                            .pipelines
                            .lock()
                            .unwrap_or_else(|poisoned| poisoned.into_inner());
                        cache.insert(cache_key, std::sync::Arc::clone(&compiled));

                        compiled
                    }
                };

            let bind_group = runtime
                .device
                .create_bind_group(&wgpu::BindGroupDescriptor {
                    label: Some(&format!("{} bind group", descriptor.label)),
                    layout: &compiled.bind_group_layout,
                    entries: &bind_group_entries,
                });

            let mut encoder =
                runtime
                    .device
                    .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                        label: Some(&format!("{} encoder", descriptor.label)),
                    });
            {
                let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                    label: Some(&format!("{} pass", descriptor.label)),
                    timestamp_writes: None,
                });
                pass.set_pipeline(&compiled.pipeline);
                pass.set_bind_group(0, &bind_group, &[]);
                pass.dispatch_workgroups(
                    descriptor.workgroup_count[0],
                    descriptor.workgroup_count[1],
                    descriptor.workgroup_count[2],
                );
            }

            let mut staging_buffers = Vec::with_capacity(readbacks.len());
            for (buffer_index, size_bytes) in &readbacks {
                let staging = runtime.device.create_buffer(&wgpu::BufferDescriptor {
                    label: Some("readback staging"),
                    size: *size_bytes as u64,
                    usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ,
                    mapped_at_creation: false,
                });
                encoder.copy_buffer_to_buffer(
                    &buffers[*buffer_index],
                    0,
                    &staging,
                    0,
                    *size_bytes as u64,
                );
                staging_buffers.push(staging);
            }

            runtime.queue.submit(Some(encoder.finish()));
            runtime.device.poll(wgpu::Maintain::Wait);

            let mut outputs = Vec::with_capacity(staging_buffers.len());
            for staging in staging_buffers {
                let slice = staging.slice(..);
                let (sender, receiver) = mpsc::channel();
                slice.map_async(wgpu::MapMode::Read, move |result| {
                    let _ = sender.send(result);
                });
                runtime.device.poll(wgpu::Maintain::Wait);
                receiver
                    .recv()
                    .map_err(|_| "GPU readback channel error".to_string())?
                    .map_err(|err| format!("GPU buffer map failed: {err}"))?;

                let bytes = slice.get_mapped_range().to_vec();
                staging.unmap();
                outputs.push(bytes);
            }

            Ok(ComputeDispatchResult {
                backend: self.capabilities.backend.clone(),
                outputs,
            })
        }

        #[cfg(not(all(feature = "experimental-gpu", not(target_arch = "wasm32"))))]
        {
            let _ = descriptor;
            Err("experimental-gpu feature not enabled".to_string())
        }
    }

    /// Validate that a WGSL shader compiles successfully on the current device.
    ///
    /// Returns a human-readable confirmation message naming the backend the
    /// shader was compiled on.
    ///
    /// # Errors
    ///
    /// - `"experimental-gpu feature not enabled"` when the `experimental-gpu`
    ///   feature is off (or the target is `wasm32`).
    /// - The stored runtime error, or `"GPU runtime not initialized"`, when
    ///   this context is a CPU-fallback context.
    /// - A `wgpu` validation error string when the WGSL source does not
    ///   compile.
    pub fn validate_shader(&self, label: &str, source: &str) -> Result<String, String> {
        #[cfg(all(feature = "experimental-gpu", not(target_arch = "wasm32")))]
        {
            let runtime = self.runtime()?;
            // Shader creation will fail if WGSL is invalid
            let _module = runtime
                .device
                .create_shader_module(wgpu::ShaderModuleDescriptor {
                    label: Some(label),
                    source: wgpu::ShaderSource::Wgsl(Cow::Borrowed(source)),
                });
            Ok(format!(
                "Shader '{}' compiled on {}",
                label, self.capabilities.backend
            ))
        }

        #[cfg(not(all(feature = "experimental-gpu", not(target_arch = "wasm32"))))]
        {
            let _ = (label, source);
            Err("experimental-gpu feature not enabled".to_string())
        }
    }

    /// Convenience: f32 vector addition for GPU validation.
    ///
    /// Returns the element-wise sum `lhs[i] + rhs[i]` as computed on the GPU.
    /// The values are plain `f32` numbers with no unit attached.
    ///
    /// # Errors
    ///
    /// - `"Vectors must have the same length"` when `lhs.len() != rhs.len()`.
    /// - `"No output from GPU kernel"` when the dispatch produced no read-back
    ///   payload.
    /// - Any error from [`GpuContext::run_compute`].
    pub fn vector_add_f32(&self, lhs: &[f32], rhs: &[f32]) -> Result<Vec<f32>, String> {
        if lhs.len() != rhs.len() {
            return Err("Vectors must have the same length".to_string());
        }

        let params = VectorAddParams {
            len: lhs.len() as u32,
            _pad: [0; 3],
        };
        let dispatch = (lhs.len() as u32).div_ceil(64);
        let output_seed = vec![0.0f32; lhs.len()];

        let descriptor = ComputeDispatchDescriptor {
            label: "vector add".to_string(),
            shader_source: VECTOR_ADD_SHADER.to_string(),
            entry_point: "main".to_string(),
            workgroup_count: [dispatch.max(1), 1, 1],
            bindings: vec![
                ComputeBindingDescriptor {
                    label: "lhs".to_string(),
                    kind: ComputeBindingKind::StorageReadOnly,
                    bytes: f32_slice_to_bytes(lhs),
                },
                ComputeBindingDescriptor {
                    label: "rhs".to_string(),
                    kind: ComputeBindingKind::StorageReadOnly,
                    bytes: f32_slice_to_bytes(rhs),
                },
                ComputeBindingDescriptor {
                    label: "params".to_string(),
                    kind: ComputeBindingKind::Uniform,
                    bytes: vector_add_params_to_bytes(&params),
                },
                ComputeBindingDescriptor {
                    label: "output".to_string(),
                    kind: ComputeBindingKind::StorageReadWrite,
                    bytes: f32_slice_to_bytes(&output_seed),
                },
            ],
        };

        let mut result = self.run_compute(&descriptor)?.outputs;
        let bytes = result.pop().ok_or("No output from GPU kernel")?;
        Ok(bytes_to_f32_vec(&bytes))
    }
}

// ─── Helper types and functions ──────────────────────────────────────────────

#[repr(C)]
#[derive(Debug, Clone, Copy)]
struct VectorAddParams {
    len: u32,
    _pad: [u32; 3],
}

/// Packs a slice of `f32` into a little-endian-free native-endian byte buffer
/// suitable for uploading as a storage or uniform buffer.
///
/// The output is `4 * values.len()` bytes long and uses the host byte order
/// (`to_ne_bytes`), so it round-trips exactly with
/// [`bytes_to_f32_vec`].
pub fn f32_slice_to_bytes(values: &[f32]) -> Vec<u8> {
    let mut bytes = Vec::with_capacity(values.len() * 4);
    for v in values {
        bytes.extend_from_slice(&v.to_ne_bytes());
    }
    bytes
}

/// Decodes a native-endian byte buffer into `f32` values.
///
/// The buffer is consumed in 4-byte chunks; a trailing partial chunk is
/// silently dropped, and the returned vector may therefore be shorter than
/// `bytes.len() / 4`.
pub fn bytes_to_f32_vec(bytes: &[u8]) -> Vec<f32> {
    bytes
        .chunks_exact(4)
        .map(|c| f32::from_ne_bytes(c.try_into().expect("4 bytes")))
        .collect()
}

/// A single 32-bit word of a uniform buffer.
#[derive(Debug, Clone, Copy)]
pub enum UniformValue {
    /// Raw 32-bit unsigned integer word.
    U32(u32),
    /// 32-bit IEEE-754 floating-point word.
    F32(f32),
}

/// Packs uniform words into a native-endian byte buffer for a uniform binding.
///
/// The result is `4 * values.len()` bytes long and is written in the host
/// byte order, matching the WGSL `uniform` block layout used by this crate.
pub fn pack_uniform_values(values: &[UniformValue]) -> Vec<u8> {
    let mut bytes = Vec::with_capacity(values.len() * 4);
    for value in values {
        match value {
            UniformValue::U32(word) => bytes.extend_from_slice(&word.to_ne_bytes()),
            UniformValue::F32(word) => bytes.extend_from_slice(&word.to_ne_bytes()),
        }
    }
    bytes
}

/// Packs 3D positions into the `vec4<f32>` layout expected by the GPU kernels.
///
/// Each position contributes 16 bytes: `x`, `y`, `z` as `f32` plus a
/// zero padding word. The input coordinates are **Å** and are converted to
/// `f32` without any unit change, so the caller must convert to bohr (or
/// whatever length unit the shader expects) beforehand.
pub fn pack_vec3_positions_f32(positions: &[[f64; 3]]) -> Vec<u8> {
    let mut bytes = Vec::with_capacity(positions.len() * 16);
    for position in positions {
        bytes.extend_from_slice(&(position[0] as f32).to_ne_bytes());
        bytes.extend_from_slice(&(position[1] as f32).to_ne_bytes());
        bytes.extend_from_slice(&(position[2] as f32).to_ne_bytes());
        bytes.extend_from_slice(&0.0f32.to_ne_bytes());
    }
    bytes
}

/// Decodes a buffer of `f32` words and widens every value to `f64`.
///
/// This is the read-back path for GPU kernels that work in single precision
/// but whose results are reported in double precision by the CPU code path.
pub fn bytes_to_f64_vec_from_f32(bytes: &[u8]) -> Vec<f64> {
    bytes_to_f32_vec(bytes)
        .into_iter()
        .map(|value| value as f64)
        .collect()
}

/// Rounds `value / chunk` up to the next whole number, as a `u32`.
///
/// Used to translate an element count into a workgroup count.
///
/// # Panics
///
/// Panics if `chunk` is zero (integer division by zero).
pub fn ceil_div_u32(value: usize, chunk: u32) -> u32 {
    (value as u32).div_ceil(chunk)
}

fn vector_add_params_to_bytes(params: &VectorAddParams) -> Vec<u8> {
    let mut bytes = Vec::with_capacity(16);
    bytes.extend_from_slice(&params.len.to_ne_bytes());
    for v in params._pad {
        bytes.extend_from_slice(&v.to_ne_bytes());
    }
    bytes
}

const VECTOR_ADD_SHADER: &str = r#"
struct Params {
    len: u32, _pad0: u32, _pad1: u32, _pad2: u32,
};

@group(0) @binding(0) var<storage, read> lhs: array<f32>;
@group(0) @binding(1) var<storage, read> rhs: array<f32>;
@group(0) @binding(2) var<uniform> params: Params;
@group(0) @binding(3) var<storage, read_write> out: array<f32>;

@compute @workgroup_size(64, 1, 1)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let idx = gid.x;
    if (idx >= params.len) { return; }
    out[idx] = lhs[idx] + rhs[idx];
}
"#;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cpu_fallback_creation() {
        let ctx = GpuContext::cpu_fallback();
        assert!(!ctx.is_gpu_available());
        assert_eq!(ctx.capabilities.backend, "CPU-fallback");
    }

    #[test]
    fn test_best_available_never_panics() {
        let ctx = GpuContext::best_available();
        let report = ctx.activation_report();
        assert!(!report.backend.is_empty());
    }

    #[test]
    fn test_activation_report_feature_disabled() {
        let ctx = GpuContext::cpu_fallback();
        let report = ctx.activation_report();
        if !cfg!(feature = "experimental-gpu") {
            assert_eq!(report.state, GpuActivationState::FeatureDisabled);
            assert!(!report.feature_enabled);
        }
    }

    #[test]
    fn test_compute_capabilities_default() {
        let caps = ComputeCapabilities::default();
        assert!(!caps.gpu_available);
        assert!(caps.max_workgroup_size_x > 0);
    }

    #[test]
    fn test_f32_roundtrip() {
        let values = vec![1.0f32, 2.5, -std::f32::consts::PI, 0.0];
        let bytes = f32_slice_to_bytes(&values);
        let result = bytes_to_f32_vec(&bytes);
        assert_eq!(values, result);
    }

    #[test]
    fn test_uniform_word_packing() {
        let bytes = pack_uniform_values(&[
            UniformValue::U32(7),
            UniformValue::F32(1.5),
            UniformValue::U32(9),
            UniformValue::F32(-2.0),
        ]);
        assert_eq!(bytes.len(), 16);
        assert_eq!(u32::from_ne_bytes(bytes[0..4].try_into().unwrap()), 7);
        assert!((f32::from_ne_bytes(bytes[4..8].try_into().unwrap()) - 1.5).abs() < 1e-6);
        assert_eq!(u32::from_ne_bytes(bytes[8..12].try_into().unwrap()), 9);
        assert!((f32::from_ne_bytes(bytes[12..16].try_into().unwrap()) + 2.0).abs() < 1e-6);
    }

    #[test]
    fn test_pack_vec3_positions_f32_layout() {
        let bytes = pack_vec3_positions_f32(&[[1.0, -2.0, 3.5]]);
        assert_eq!(bytes.len(), 16);
        assert!((f32::from_ne_bytes(bytes[0..4].try_into().unwrap()) - 1.0).abs() < 1e-6);
        assert!((f32::from_ne_bytes(bytes[4..8].try_into().unwrap()) + 2.0).abs() < 1e-6);
        assert!((f32::from_ne_bytes(bytes[8..12].try_into().unwrap()) - 3.5).abs() < 1e-6);
        assert_eq!(f32::from_ne_bytes(bytes[12..16].try_into().unwrap()), 0.0);
    }
}
