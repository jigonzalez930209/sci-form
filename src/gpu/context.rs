//! GPU compute context — wgpu device initialization with CPU fallback.
//!
//! Wraps `wgpu::Device` + `Queue` for native GPU compute, or falls back
//! to CPU-only mode when no GPU is available or the feature is disabled.

#[cfg(all(feature = "experimental-gpu", not(target_arch = "wasm32")))]
use std::borrow::Cow;

#[cfg(all(feature = "experimental-gpu", not(target_arch = "wasm32")))]
use std::sync::mpsc;

#[cfg(all(feature = "experimental-gpu", not(target_arch = "wasm32")))]
use wgpu::util::DeviceExt;

use super::backend_report::{GpuActivationReport, GpuActivationState};

/// Capabilities detected from the GPU adapter.
#[derive(Debug, Clone)]
pub struct ComputeCapabilities {
    pub backend: String,
    pub max_workgroup_size_x: u32,
    pub max_workgroup_size_y: u32,
    pub max_workgroup_invocations: u32,
    pub max_storage_buffer_size: u64,
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
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ComputeBindingKind {
    Uniform,
    StorageReadOnly,
    StorageReadWrite,
}

/// One logical binding passed to a compute kernel.
#[derive(Debug, Clone)]
pub struct ComputeBindingDescriptor {
    pub label: String,
    pub kind: ComputeBindingKind,
    pub bytes: Vec<u8>,
}

/// Generic compute dispatch description.
#[derive(Debug, Clone)]
pub struct ComputeDispatchDescriptor {
    pub label: String,
    pub shader_source: String,
    pub entry_point: String,
    pub workgroup_count: [u32; 3],
    pub bindings: Vec<ComputeBindingDescriptor>,
}

/// Output produced by a compute dispatch.
#[derive(Debug, Clone)]
pub struct ComputeDispatchResult {
    pub backend: String,
    pub outputs: Vec<Vec<u8>>,
}

/// GPU compute context handle.
pub struct GpuContext {
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
    pub fn try_create() -> Result<Self, String> {
        #[cfg(all(feature = "experimental-gpu", not(target_arch = "wasm32")))]
        {
            let instance = wgpu::Instance::default();
            let adapter =
                pollster::block_on(instance.request_adapter(&wgpu::RequestAdapterOptions {
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
                }),
            })
        }

        #[cfg(not(all(feature = "experimental-gpu", not(target_arch = "wasm32"))))]
        {
            Err("experimental-gpu feature not enabled".to_string())
        }
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
    pub fn run_compute(
        &self,
        descriptor: &ComputeDispatchDescriptor,
    ) -> Result<ComputeDispatchResult, String> {
        #[cfg(all(feature = "experimental-gpu", not(target_arch = "wasm32")))]
        {
            let runtime = self.runtime()?;

            let shader = runtime
                .device
                .create_shader_module(wgpu::ShaderModuleDescriptor {
                    label: Some(&descriptor.label),
                    source: wgpu::ShaderSource::Wgsl(Cow::Borrowed(&descriptor.shader_source)),
                });

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

            let bind_group_layout =
                runtime
                    .device
                    .create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
                        label: Some(&format!("{} layout", descriptor.label)),
                        entries: &layout_entries,
                    });
            let pipeline_layout =
                runtime
                    .device
                    .create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
                        label: Some(&format!("{} pipeline", descriptor.label)),
                        bind_group_layouts: &[&bind_group_layout],
                        push_constant_ranges: &[],
                    });
            let pipeline =
                runtime
                    .device
                    .create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
                        label: Some(&descriptor.label),
                        layout: Some(&pipeline_layout),
                        module: &shader,
                        entry_point: &descriptor.entry_point,
                    });
            let bind_group = runtime
                .device
                .create_bind_group(&wgpu::BindGroupDescriptor {
                    label: Some(&format!("{} bind group", descriptor.label)),
                    layout: &bind_group_layout,
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
                pass.set_pipeline(&pipeline);
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

pub fn f32_slice_to_bytes(values: &[f32]) -> Vec<u8> {
    let mut bytes = Vec::with_capacity(values.len() * 4);
    for v in values {
        bytes.extend_from_slice(&v.to_ne_bytes());
    }
    bytes
}

pub fn bytes_to_f32_vec(bytes: &[u8]) -> Vec<f32> {
    bytes
        .chunks_exact(4)
        .map(|c| f32::from_ne_bytes(c.try_into().expect("4 bytes")))
        .collect()
}

#[derive(Debug, Clone, Copy)]
pub enum UniformValue {
    U32(u32),
    F32(f32),
}

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

pub fn bytes_to_f64_vec_from_f32(bytes: &[u8]) -> Vec<f64> {
    bytes_to_f32_vec(bytes)
        .into_iter()
        .map(|value| value as f64)
        .collect()
}

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
