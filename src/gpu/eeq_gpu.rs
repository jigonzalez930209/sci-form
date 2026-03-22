//! GPU-assisted EEQ charges: GPU damped Coulomb matrix + CPU linear solve.

use nalgebra::{DMatrix, DVector};

use super::context::{
    bytes_to_f64_vec_from_f32, ceil_div_u32, f32_slice_to_bytes, pack_uniform_values,
    pack_vec3_positions_f32, ComputeBindingDescriptor, ComputeBindingKind,
    ComputeDispatchDescriptor, GpuContext, UniformValue,
};
use crate::charges_eeq::{fractional_coordination, get_eeq_params, EeqChargeResult, EeqConfig};

const GPU_DISPATCH_THRESHOLD: usize = 2;

pub fn compute_eeq_charges_gpu(
    ctx: &GpuContext,
    elements: &[u8],
    positions: &[[f64; 3]],
    config: &EeqConfig,
) -> Result<EeqChargeResult, String> {
    let n = elements.len();
    if positions.len() != n {
        return Err("elements/position length mismatch".to_string());
    }
    if n < GPU_DISPATCH_THRESHOLD {
        return Err("System too small for GPU dispatch".to_string());
    }

    let params: Vec<_> = elements.iter().map(|&z| get_eeq_params(z)).collect();
    let cn = fractional_coordination(elements, positions);
    let gamma = build_eeq_coulomb_gpu(
        ctx,
        &params.iter().map(|p| p.r_eeq).collect::<Vec<_>>(),
        positions,
    )?;

    let dim = n + 1;
    let mut a = DMatrix::zeros(dim, dim);
    let mut b_vec = vec![0.0; dim];

    for i in 0..n {
        a[(i, i)] = params[i].eta + config.regularization;

        for j in (i + 1)..n {
            let gij = gamma[(i, j)];
            a[(i, j)] = gij;
            a[(j, i)] = gij;
        }

        a[(i, n)] = 1.0;
        a[(n, i)] = 1.0;

        let cn_correction = -0.1 * (cn[i] - 2.0);
        b_vec[i] = -(params[i].chi + cn_correction);
    }

    b_vec[n] = config.total_charge;

    let solution = a.lu().solve(&DVector::from_vec(b_vec));
    let charges = match solution {
        Some(sol) => (0..n).map(|i| sol[i]).collect(),
        None => vec![0.0; n],
    };
    let total_charge = charges.iter().sum();

    Ok(EeqChargeResult {
        charges,
        coordination_numbers: cn,
        total_charge,
    })
}

pub fn build_eeq_coulomb_gpu(
    ctx: &GpuContext,
    radii: &[f64],
    positions: &[[f64; 3]],
) -> Result<DMatrix<f64>, String> {
    let n = radii.len();
    if positions.len() != n {
        return Err("radii/position length mismatch".to_string());
    }
    if n < GPU_DISPATCH_THRESHOLD {
        return Err("Matrix too small for GPU dispatch".to_string());
    }

    let descriptor = ComputeDispatchDescriptor {
        label: "eeq coulomb".to_string(),
        shader_source: EEQ_COULOMB_SHADER.to_string(),
        entry_point: "main".to_string(),
        workgroup_count: [ceil_div_u32(n, 16), ceil_div_u32(n, 16), 1],
        bindings: vec![
            ComputeBindingDescriptor {
                label: "positions".to_string(),
                kind: ComputeBindingKind::StorageReadOnly,
                bytes: pack_vec3_positions_f32(positions),
            },
            ComputeBindingDescriptor {
                label: "radii".to_string(),
                kind: ComputeBindingKind::StorageReadOnly,
                bytes: f32_slice_to_bytes(
                    &radii.iter().map(|value| *value as f32).collect::<Vec<_>>(),
                ),
            },
            ComputeBindingDescriptor {
                label: "params".to_string(),
                kind: ComputeBindingKind::Uniform,
                bytes: pack_uniform_values(&[
                    UniformValue::U32(n as u32),
                    UniformValue::U32(0),
                    UniformValue::U32(0),
                    UniformValue::U32(0),
                ]),
            },
            ComputeBindingDescriptor {
                label: "output".to_string(),
                kind: ComputeBindingKind::StorageReadWrite,
                bytes: f32_slice_to_bytes(&vec![0.0f32; n * n]),
            },
        ],
    };

    let mut outputs = ctx.run_compute(&descriptor)?.outputs;
    let bytes = outputs.pop().ok_or("No output from EEQ Coulomb kernel")?;
    if bytes.len() != n * n * 4 {
        return Err(format!(
            "EEQ Coulomb output size mismatch: expected {}, got {}",
            n * n * 4,
            bytes.len()
        ));
    }

    Ok(DMatrix::from_row_slice(
        n,
        n,
        &bytes_to_f64_vec_from_f32(&bytes),
    ))
}

pub const EEQ_COULOMB_SHADER: &str = r#"
struct AtomPos {
    x: f32, y: f32, z: f32, _pad: f32,
};

struct Params {
    n_atoms: u32,
    _pad0: u32,
    _pad1: u32,
    _pad2: u32,
};

@group(0) @binding(0) var<storage, read> positions: array<AtomPos>;
@group(0) @binding(1) var<storage, read> radii: array<f32>;
@group(0) @binding(2) var<uniform> params: Params;
@group(0) @binding(3) var<storage, read_write> output: array<f32>;

fn erf_approx(x: f32) -> f32 {
    let a1 = 0.254829592;
    let a2 = -0.284496736;
    let a3 = 1.421413741;
    let a4 = -1.453152027;
    let a5 = 1.061405429;
    let p = 0.3275911;

    let sign = select(-1.0, 1.0, x >= 0.0);
    let ax = abs(x);
    let t = 1.0 / (1.0 + p * ax);
    let y = 1.0 - (((((a5 * t + a4) * t + a3) * t + a2) * t + a1) * t) * exp(-ax * ax);
    return sign * y;
}

@compute @workgroup_size(16, 16, 1)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    let j = gid.y;
    let n = params.n_atoms;
    if (i >= n || j >= n) {
        return;
    }
    if (i == j) {
        output[i * n + j] = 0.0;
        return;
    }

    let pi = positions[i];
    let pj = positions[j];
    let dx = pi.x - pj.x;
    let dy = pi.y - pj.y;
    let dz = pi.z - pj.z;
    let r = sqrt(dx * dx + dy * dy + dz * dz);
    if (r < 1e-10) {
        output[i * n + j] = 0.0;
        return;
    }

    let ri = radii[i];
    let rj = radii[j];
    let sigma = sqrt(ri * ri + rj * rj);
    let arg = 1.41421356237 * r / sigma;
    output[i * n + j] = erf_approx(arg) / r;
}
"#;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_eeq_gpu_threshold() {
        let ctx = GpuContext::cpu_fallback();
        let result = compute_eeq_charges_gpu(&ctx, &[8], &[[0.0, 0.0, 0.0]], &EeqConfig::default());
        assert!(result.is_err());
    }
}
