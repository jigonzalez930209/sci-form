//! WebGPU runtime for sci-form WASM bindings.
//!
//! This module adds an async WebGPU compute path for browsers supporting
//! `wgpu`/WebGPU, while keeping the existing CPU-oriented WASM exports intact.

use std::cell::RefCell;
use std::rc::Rc;

use futures_channel::oneshot;
use sci_form::gpu::context::{
    bytes_to_f32_vec, ComputeBindingKind, ComputeDispatchDescriptor, ComputeDispatchResult,
};
use wasm_bindgen::prelude::*;
use wgpu::util::DeviceExt;

thread_local! {
    static WEBGPU_RUNTIME: RefCell<Option<WebGpuRuntime>> = const { RefCell::new(None) };
}

#[derive(Clone)]
pub(crate) struct WebGpuRuntime {
    pub(crate) device: Rc<wgpu::Device>,
    pub(crate) queue: Rc<wgpu::Queue>,
    pub(crate) backend: String,
    pub(crate) max_storage_buffer_size: u64,
    pub(crate) max_workgroup_size_x: u32,
    pub(crate) max_workgroup_size_y: u32,
    pub(crate) max_workgroup_invocations: u32,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WasmExecutionMode {
    Auto,
    Cpu,
    Gpu,
    Hybrid,
}

pub fn parse_execution_mode(mode: &str) -> Result<WasmExecutionMode, String> {
    match mode.trim().to_ascii_lowercase().as_str() {
        "" | "auto" => Ok(WasmExecutionMode::Auto),
        "cpu" => Ok(WasmExecutionMode::Cpu),
        "gpu" => Ok(WasmExecutionMode::Gpu),
        "hybrid" => Ok(WasmExecutionMode::Hybrid),
        other => Err(format!("unsupported execution mode: {other}")),
    }
}

async fn create_runtime() -> Result<WebGpuRuntime, String> {
    let instance = wgpu::Instance::default();
    let adapter = instance
        .request_adapter(&wgpu::RequestAdapterOptions {
            power_preference: wgpu::PowerPreference::HighPerformance,
            compatible_surface: None,
            force_fallback_adapter: false,
        })
        .await
        .map_err(|_| "No WebGPU adapter found".to_string())?;

    let adapter_info = adapter.get_info();
    let limits = adapter.limits();
    let required_limits = wgpu::Limits::default().using_resolution(limits.clone());

    let (device, queue) = adapter
        .request_device(&wgpu::DeviceDescriptor {
            label: Some("sci-form wasm webgpu"),
            required_features: wgpu::Features::empty(),
            required_limits,
            experimental_features: wgpu::ExperimentalFeatures::disabled(),
            memory_hints: wgpu::MemoryHints::Performance,
            trace: wgpu::Trace::Off,
        })
        .await
        .map_err(|err| format!("Failed to create WebGPU device: {err}"))?;

    Ok(WebGpuRuntime {
        device: Rc::new(device),
        queue: Rc::new(queue),
        backend: format!("{:?}", adapter_info.backend),
        max_storage_buffer_size: limits.max_storage_buffer_binding_size as u64,
        max_workgroup_size_x: limits.max_compute_workgroup_size_x,
        max_workgroup_size_y: limits.max_compute_workgroup_size_y,
        max_workgroup_invocations: limits.max_compute_invocations_per_workgroup,
    })
}

pub(crate) async fn ensure_runtime() -> Result<WebGpuRuntime, String> {
    if let Some(runtime) = WEBGPU_RUNTIME.with(|slot| slot.borrow().clone()) {
        return Ok(runtime);
    }

    let runtime = create_runtime().await?;
    WEBGPU_RUNTIME.with(|slot| {
        *slot.borrow_mut() = Some(runtime.clone());
    });
    Ok(runtime)
}

pub(crate) async fn try_runtime() -> Option<WebGpuRuntime> {
    ensure_runtime().await.ok()
}

pub(crate) async fn run_compute_async(
    descriptor: &ComputeDispatchDescriptor,
) -> Result<ComputeDispatchResult, String> {
    let runtime = ensure_runtime().await?;

    let shader = runtime
        .device
        .create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some(&descriptor.label),
            source: wgpu::ShaderSource::Wgsl(descriptor.shader_source.clone().into()),
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
    let pipeline_layout = runtime
        .device
        .create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some(&format!("{} pipeline", descriptor.label)),
            bind_group_layouts: &[Some(&bind_group_layout)],
            immediate_size: 0,
        });
    let pipeline = runtime
        .device
        .create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
            label: Some(&descriptor.label),
            layout: Some(&pipeline_layout),
            module: &shader,
            entry_point: Some(descriptor.entry_point.as_str()),
            compilation_options: wgpu::PipelineCompilationOptions::default(),
            cache: None,
        });
    let bind_group = runtime
        .device
        .create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some(&format!("{} bind group", descriptor.label)),
            layout: &bind_group_layout,
            entries: &bind_group_entries,
        });

    let mut encoder = runtime
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
        encoder.copy_buffer_to_buffer(&buffers[*buffer_index], 0, &staging, 0, *size_bytes as u64);
        staging_buffers.push(staging);
    }

    runtime.queue.submit(Some(encoder.finish()));
    let _ = runtime.device.poll(wgpu::PollType::Poll);

    let mut outputs = Vec::with_capacity(staging_buffers.len());
    for staging in staging_buffers {
        let slice = staging.slice(..);
        let (sender, receiver) = oneshot::channel();
        slice.map_async(wgpu::MapMode::Read, move |result| {
            let _ = sender.send(result);
        });
        let _ = runtime.device.poll(wgpu::PollType::Poll);
        receiver
            .await
            .map_err(|_| "WebGPU readback channel error".to_string())?
            .map_err(|err| format!("WebGPU buffer map failed: {err}"))?;

        let mapped = slice.get_mapped_range();
        let bytes = mapped.to_vec();
        drop(mapped);
        staging.unmap();
        outputs.push(bytes);
    }

    Ok(ComputeDispatchResult {
        backend: runtime.backend,
        outputs,
    })
}

pub fn current_status_json() -> String {
    if let Some(runtime) = WEBGPU_RUNTIME.with(|slot| slot.borrow().clone()) {
        serde_json::json!({
            "initialized": true,
            "backend": runtime.backend,
            "gpu_available": true,
            "max_storage_buffer_size": runtime.max_storage_buffer_size,
            "max_workgroup_size_x": runtime.max_workgroup_size_x,
            "max_workgroup_size_y": runtime.max_workgroup_size_y,
            "max_workgroup_invocations": runtime.max_workgroup_invocations,
        })
        .to_string()
    } else {
        serde_json::json!({
            "initialized": false,
            "backend": "CPU-fallback",
            "gpu_available": false,
        })
        .to_string()
    }
}

#[wasm_bindgen]
pub async fn init_webgpu() -> String {
    match ensure_runtime().await {
        Ok(runtime) => serde_json::json!({
            "initialized": true,
            "backend": runtime.backend,
            "gpu_available": true,
            "max_storage_buffer_size": runtime.max_storage_buffer_size,
            "max_workgroup_size_x": runtime.max_workgroup_size_x,
            "max_workgroup_size_y": runtime.max_workgroup_size_y,
            "max_workgroup_invocations": runtime.max_workgroup_invocations,
        })
        .to_string(),
        Err(err) => serde_json::json!({
            "initialized": false,
            "gpu_available": false,
            "error": err,
        })
        .to_string(),
    }
}

#[wasm_bindgen]
pub fn webgpu_status() -> String {
    current_status_json()
}

pub(crate) fn bytes_to_f64_vec(bytes: &[u8]) -> Vec<f64> {
    bytes_to_f32_vec(bytes)
        .into_iter()
        .map(|value| value as f64)
        .collect()
}
