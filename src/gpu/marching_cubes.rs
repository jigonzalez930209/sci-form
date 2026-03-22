//! Marching cubes isosurface extraction from 3D scalar fields.
//!
//! Converts a scalar field (orbital ψ or density ρ) into a triangle mesh
//! at a specified isovalue. Reference: Lorensen & Cline, SIGGRAPH 1987.
//!
//! GPU path: dispatches a WGSL compute shader that processes each voxel
//! independently, then reads back the generated triangle vertices.
//! CPU path: triple-nested loop (always available as fallback).

use super::backend_report::OrbitalGridReport;
use super::context::{
    bytes_to_f32_vec, ceil_div_u32, f32_slice_to_bytes, pack_uniform_values,
    ComputeBindingDescriptor, ComputeBindingKind, ComputeDispatchDescriptor, GpuContext,
    UniformValue,
};
use super::orbital_grid::GridParams;

/// Triangle mesh output from marching cubes.
#[derive(Debug, Clone)]
pub struct McOutput {
    /// Triangle vertices (flat: [x0,y0,z0, x1,y1,z1, ...]).
    pub vertices: Vec<f32>,
    /// Triangle normals (flat: [nx0,ny0,nz0, ...]).
    pub normals: Vec<f32>,
    /// Triangle indices.
    pub indices: Vec<u32>,
    /// Number of triangles.
    pub n_triangles: usize,
}

/// CPU marching cubes on a 3D scalar field.
pub fn marching_cubes_cpu(values: &[f64], params: &GridParams, isovalue: f64) -> McOutput {
    let [nx, ny, nz] = params.dimensions;
    let mut vertices = Vec::new();
    let mut normals = Vec::new();
    let mut indices = Vec::new();
    let mut vertex_count = 0u32;

    for ix in 0..nx.saturating_sub(1) {
        for iy in 0..ny.saturating_sub(1) {
            for iz in 0..nz.saturating_sub(1) {
                let corners = [
                    values[params.flat_index(ix, iy, iz)],
                    values[params.flat_index(ix + 1, iy, iz)],
                    values[params.flat_index(ix + 1, iy + 1, iz)],
                    values[params.flat_index(ix, iy + 1, iz)],
                    values[params.flat_index(ix, iy, iz + 1)],
                    values[params.flat_index(ix + 1, iy, iz + 1)],
                    values[params.flat_index(ix + 1, iy + 1, iz + 1)],
                    values[params.flat_index(ix, iy + 1, iz + 1)],
                ];

                let mut cube_index = 0u8;
                for (i, &val) in corners.iter().enumerate() {
                    if val > isovalue {
                        cube_index |= 1 << i;
                    }
                }

                let edge_bits = EDGE_TABLE[cube_index as usize];
                if edge_bits == 0 {
                    continue;
                }

                let corner_positions = [
                    params.point(ix, iy, iz),
                    params.point(ix + 1, iy, iz),
                    params.point(ix + 1, iy + 1, iz),
                    params.point(ix, iy + 1, iz),
                    params.point(ix, iy, iz + 1),
                    params.point(ix + 1, iy, iz + 1),
                    params.point(ix + 1, iy + 1, iz + 1),
                    params.point(ix, iy + 1, iz + 1),
                ];

                let mut edge_vertices = [[0.0f64; 3]; 12];
                for edge in 0..12 {
                    if edge_bits & (1 << edge) != 0 {
                        let (v0, v1) = EDGE_VERTICES[edge];
                        let p0 = corner_positions[v0];
                        let p1 = corner_positions[v1];
                        let val0 = corners[v0];
                        let val1 = corners[v1];

                        let t = if (val1 - val0).abs() > 1e-10 {
                            ((isovalue - val0) / (val1 - val0)).clamp(0.0, 1.0)
                        } else {
                            0.5
                        };

                        edge_vertices[edge] = [
                            p0[0] + t * (p1[0] - p0[0]),
                            p0[1] + t * (p1[1] - p0[1]),
                            p0[2] + t * (p1[2] - p0[2]),
                        ];
                    }
                }

                let tri_entry = &TRI_TABLE[cube_index as usize];
                let mut i = 0;
                while i < tri_entry.len() && tri_entry[i] != -1 {
                    let e0 = tri_entry[i] as usize;
                    let e1 = tri_entry[i + 1] as usize;
                    let e2 = tri_entry[i + 2] as usize;

                    let v0 = edge_vertices[e0];
                    let v1 = edge_vertices[e1];
                    let v2 = edge_vertices[e2];

                    // Face normal via cross product
                    let u = [v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]];
                    let v = [v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]];
                    let mut n = [
                        u[1] * v[2] - u[2] * v[1],
                        u[2] * v[0] - u[0] * v[2],
                        u[0] * v[1] - u[1] * v[0],
                    ];
                    let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
                    if len > 1e-10 {
                        n[0] /= len;
                        n[1] /= len;
                        n[2] /= len;
                    }

                    for vert in &[v0, v1, v2] {
                        vertices.push(vert[0] as f32);
                        vertices.push(vert[1] as f32);
                        vertices.push(vert[2] as f32);
                        normals.push(n[0] as f32);
                        normals.push(n[1] as f32);
                        normals.push(n[2] as f32);
                        indices.push(vertex_count);
                        vertex_count += 1;
                    }

                    i += 3;
                }
            }
        }
    }

    McOutput {
        n_triangles: indices.len() / 3,
        vertices,
        normals,
        indices,
    }
}

/// Smooth vertex normals by averaging at shared positions.
pub fn smooth_mesh_normals(mesh: &mut McOutput) {
    if mesh.vertices.is_empty() {
        return;
    }

    let n_verts = mesh.vertices.len() / 3;
    let tolerance = 1e-4f32;
    let mut smoothed = vec![[0.0f32; 3]; n_verts];

    for i in 0..n_verts {
        let pi = [
            mesh.vertices[i * 3],
            mesh.vertices[i * 3 + 1],
            mesh.vertices[i * 3 + 2],
        ];
        let mut nx = 0.0f32;
        let mut ny = 0.0f32;
        let mut nz = 0.0f32;
        let mut count = 0u32;

        for j in 0..n_verts {
            let pj = [
                mesh.vertices[j * 3],
                mesh.vertices[j * 3 + 1],
                mesh.vertices[j * 3 + 2],
            ];
            let dx = pi[0] - pj[0];
            let dy = pi[1] - pj[1];
            let dz = pi[2] - pj[2];
            if dx * dx + dy * dy + dz * dz < tolerance * tolerance {
                nx += mesh.normals[j * 3];
                ny += mesh.normals[j * 3 + 1];
                nz += mesh.normals[j * 3 + 2];
                count += 1;
            }
        }

        if count > 0 {
            let len = (nx * nx + ny * ny + nz * nz).sqrt();
            if len > 1e-10 {
                smoothed[i] = [nx / len, ny / len, nz / len];
            } else {
                smoothed[i] = [
                    mesh.normals[i * 3],
                    mesh.normals[i * 3 + 1],
                    mesh.normals[i * 3 + 2],
                ];
            }
        }
    }

    for i in 0..n_verts {
        mesh.normals[i * 3] = smoothed[i][0];
        mesh.normals[i * 3 + 1] = smoothed[i][1];
        mesh.normals[i * 3 + 2] = smoothed[i][2];
    }
}

/// Marching cubes with automatic GPU/CPU backend selection.
///
/// Uses GPU when available, falls back to CPU otherwise.
pub fn marching_cubes_with_report(
    values: &[f64],
    params: &GridParams,
    isovalue: f64,
) -> (McOutput, OrbitalGridReport) {
    let ctx = GpuContext::best_available();
    if ctx.is_gpu_available() {
        if let Ok(mesh) = marching_cubes_gpu(&ctx, values, params, isovalue) {
            return (
                mesh,
                OrbitalGridReport {
                    backend: ctx.capabilities.backend.clone(),
                    used_gpu: true,
                    attempted_gpu: true,
                    n_points: params.n_points(),
                    note: format!("GPU marching cubes on {}", ctx.capabilities.backend),
                },
            );
        }
    }
    let mesh = marching_cubes_cpu(values, params, isovalue);
    (
        mesh,
        OrbitalGridReport {
            backend: "CPU".to_string(),
            used_gpu: false,
            attempted_gpu: ctx.is_gpu_available(),
            n_points: params.n_points(),
            note: "CPU marching cubes".to_string(),
        },
    )
}

/// GPU dispatch for marching cubes isosurface extraction.
///
/// Each workgroup processes one voxel. The shader classifies corners,
/// looks up the edge/triangle tables, interpolates edge vertices, and
/// writes triangle data into an append-style output buffer.
pub fn marching_cubes_gpu(
    ctx: &GpuContext,
    values: &[f64],
    params: &GridParams,
    isovalue: f64,
) -> Result<McOutput, String> {
    let [nx, ny, nz] = params.dimensions;
    let vx = nx.saturating_sub(1).max(1);
    let vy = ny.saturating_sub(1).max(1);
    let vz = nz.saturating_sub(1).max(1);
    let n_voxels = vx * vy * vz;

    // Each voxel can produce at most 5 triangles (15 vertices × 6 floats each)
    let max_floats_per_voxel = 5 * 3 * 6; // 5 triangles × 3 verts × (pos3 + norm3)
    let max_output = n_voxels * max_floats_per_voxel;

    let values_f32: Vec<f32> = values.iter().map(|&v| v as f32).collect();

    let params_bytes = pack_uniform_values(&[
        UniformValue::F32(params.origin[0] as f32),
        UniformValue::F32(params.origin[1] as f32),
        UniformValue::F32(params.origin[2] as f32),
        UniformValue::F32(params.spacing as f32),
        UniformValue::U32(nx as u32),
        UniformValue::U32(ny as u32),
        UniformValue::U32(nz as u32),
        UniformValue::F32(isovalue as f32),
    ]);

    // Pack edge table as u32 array
    let edge_table_u32: Vec<u32> = EDGE_TABLE.iter().map(|&v| v as u32).collect();
    let edge_table_bytes: Vec<u8> = edge_table_u32
        .iter()
        .flat_map(|v| v.to_le_bytes())
        .collect();

    // Pack tri table as i32 array (flattened 256 × 16)
    let mut tri_table_i32: Vec<i32> = Vec::with_capacity(256 * 16);
    for row in &TRI_TABLE {
        for &val in row.iter() {
            tri_table_i32.push(val as i32);
        }
    }
    let tri_table_bytes: Vec<u8> = tri_table_i32.iter().flat_map(|v| v.to_le_bytes()).collect();

    // Atomic counter buffer (1 u32 for triangle count)
    let counter_bytes: Vec<u8> = vec![0u8; 4];

    // Output buffer for vertex+normal data
    let output_bytes = f32_slice_to_bytes(&vec![0.0f32; max_output]);

    let descriptor = ComputeDispatchDescriptor {
        label: "marching cubes".to_string(),
        shader_source: MARCHING_CUBES_SHADER.to_string(),
        entry_point: "main".to_string(),
        workgroup_count: [
            ceil_div_u32(vx, 4),
            ceil_div_u32(vy, 4),
            ceil_div_u32(vz, 4),
        ],
        bindings: vec![
            ComputeBindingDescriptor {
                label: "scalar_field".to_string(),
                kind: ComputeBindingKind::StorageReadOnly,
                bytes: f32_slice_to_bytes(&values_f32),
            },
            ComputeBindingDescriptor {
                label: "edge_table".to_string(),
                kind: ComputeBindingKind::StorageReadOnly,
                bytes: edge_table_bytes,
            },
            ComputeBindingDescriptor {
                label: "tri_table".to_string(),
                kind: ComputeBindingKind::StorageReadOnly,
                bytes: tri_table_bytes,
            },
            ComputeBindingDescriptor {
                label: "params".to_string(),
                kind: ComputeBindingKind::Uniform,
                bytes: params_bytes,
            },
            ComputeBindingDescriptor {
                label: "counter".to_string(),
                kind: ComputeBindingKind::StorageReadWrite,
                bytes: counter_bytes,
            },
            ComputeBindingDescriptor {
                label: "output".to_string(),
                kind: ComputeBindingKind::StorageReadWrite,
                bytes: output_bytes,
            },
        ],
    };

    let result = ctx.run_compute(&descriptor)?;
    if result.outputs.len() < 2 {
        return Err("Missing outputs from marching cubes kernel".to_string());
    }

    // First ReadWrite output is the counter, second is the vertex data
    let counter_out = &result.outputs[0];
    let tri_count = if counter_out.len() >= 4 {
        u32::from_le_bytes([
            counter_out[0],
            counter_out[1],
            counter_out[2],
            counter_out[3],
        ])
    } else {
        0
    };

    let vertex_data = bytes_to_f32_vec(&result.outputs[1]);
    let n_triangles = tri_count as usize;
    let n_verts = n_triangles * 3;
    let n_floats = n_verts * 6; // 3 pos + 3 normal per vertex

    let actual_floats = n_floats.min(vertex_data.len());
    let actual_verts = actual_floats / 6;
    let actual_tris = actual_verts / 3;

    let mut vertices = Vec::with_capacity(actual_verts * 3);
    let mut normals = Vec::with_capacity(actual_verts * 3);
    let mut indices = Vec::with_capacity(actual_verts);

    for v in 0..actual_verts {
        let base = v * 6;
        if base + 5 < vertex_data.len() {
            vertices.push(vertex_data[base]);
            vertices.push(vertex_data[base + 1]);
            vertices.push(vertex_data[base + 2]);
            normals.push(vertex_data[base + 3]);
            normals.push(vertex_data[base + 4]);
            normals.push(vertex_data[base + 5]);
            indices.push(v as u32);
        }
    }

    Ok(McOutput {
        n_triangles: actual_tris,
        vertices,
        normals,
        indices,
    })
}

/// WGSL compute shader for marching cubes isosurface extraction.
///
/// Each invocation processes one voxel, classifying corners against the isovalue,
/// looking up the edge and triangle tables, and emitting interpolated triangle
/// vertices with face normals into an append-style output buffer.
const MARCHING_CUBES_SHADER: &str = r#"
struct McParams {
    origin_x: f32, origin_y: f32, origin_z: f32,
    spacing: f32,
    dims_x: u32, dims_y: u32, dims_z: u32,
    isovalue: f32,
};

@group(0) @binding(0) var<storage, read> scalar_field: array<f32>;
@group(0) @binding(1) var<storage, read> edge_table: array<u32>;
@group(0) @binding(2) var<storage, read> tri_table: array<i32>;
@group(0) @binding(3) var<uniform> params: McParams;
@group(0) @binding(4) var<storage, read_write> tri_counter: array<atomic<u32>>;
@group(0) @binding(5) var<storage, read_write> output: array<f32>;

fn field_index(ix: u32, iy: u32, iz: u32) -> u32 {
    return ix * params.dims_y * params.dims_z + iy * params.dims_z + iz;
}

fn grid_point(ix: u32, iy: u32, iz: u32) -> vec3<f32> {
    return vec3<f32>(
        params.origin_x + f32(ix) * params.spacing,
        params.origin_y + f32(iy) * params.spacing,
        params.origin_z + f32(iz) * params.spacing,
    );
}

fn interp_vertex(p0: vec3<f32>, p1: vec3<f32>, v0: f32, v1: f32, iso: f32) -> vec3<f32> {
    let dv = v1 - v0;
    var t: f32 = 0.5;
    if (abs(dv) > 1e-10) {
        t = clamp((iso - v0) / dv, 0.0, 1.0);
    }
    return mix(p0, p1, vec3<f32>(t));
}

@compute @workgroup_size(4, 4, 4)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let ix = gid.x;
    let iy = gid.y;
    let iz = gid.z;

    let vx = params.dims_x - 1u;
    let vy = params.dims_y - 1u;
    let vz = params.dims_z - 1u;

    if (ix >= vx || iy >= vy || iz >= vz) {
        return;
    }

    // Sample 8 corners
    let c0 = scalar_field[field_index(ix,     iy,     iz)];
    let c1 = scalar_field[field_index(ix + 1, iy,     iz)];
    let c2 = scalar_field[field_index(ix + 1, iy + 1, iz)];
    let c3 = scalar_field[field_index(ix,     iy + 1, iz)];
    let c4 = scalar_field[field_index(ix,     iy,     iz + 1)];
    let c5 = scalar_field[field_index(ix + 1, iy,     iz + 1)];
    let c6 = scalar_field[field_index(ix + 1, iy + 1, iz + 1)];
    let c7 = scalar_field[field_index(ix,     iy + 1, iz + 1)];

    let iso = params.isovalue;

    var cube_index: u32 = 0u;
    if (c0 > iso) { cube_index |= 1u; }
    if (c1 > iso) { cube_index |= 2u; }
    if (c2 > iso) { cube_index |= 4u; }
    if (c3 > iso) { cube_index |= 8u; }
    if (c4 > iso) { cube_index |= 16u; }
    if (c5 > iso) { cube_index |= 32u; }
    if (c6 > iso) { cube_index |= 64u; }
    if (c7 > iso) { cube_index |= 128u; }

    let edge_bits = edge_table[cube_index];
    if (edge_bits == 0u) {
        return;
    }

    // Corner positions
    let p0 = grid_point(ix,     iy,     iz);
    let p1 = grid_point(ix + 1, iy,     iz);
    let p2 = grid_point(ix + 1, iy + 1, iz);
    let p3 = grid_point(ix,     iy + 1, iz);
    let p4 = grid_point(ix,     iy,     iz + 1);
    let p5 = grid_point(ix + 1, iy,     iz + 1);
    let p6 = grid_point(ix + 1, iy + 1, iz + 1);
    let p7 = grid_point(ix,     iy + 1, iz + 1);

    // Interpolate edge vertices
    var ev: array<vec3<f32>, 12>;
    if ((edge_bits &    1u) != 0u) { ev[0]  = interp_vertex(p0, p1, c0, c1, iso); }
    if ((edge_bits &    2u) != 0u) { ev[1]  = interp_vertex(p1, p2, c1, c2, iso); }
    if ((edge_bits &    4u) != 0u) { ev[2]  = interp_vertex(p2, p3, c2, c3, iso); }
    if ((edge_bits &    8u) != 0u) { ev[3]  = interp_vertex(p3, p0, c3, c0, iso); }
    if ((edge_bits &   16u) != 0u) { ev[4]  = interp_vertex(p4, p5, c4, c5, iso); }
    if ((edge_bits &   32u) != 0u) { ev[5]  = interp_vertex(p5, p6, c5, c6, iso); }
    if ((edge_bits &   64u) != 0u) { ev[6]  = interp_vertex(p6, p7, c6, c7, iso); }
    if ((edge_bits &  128u) != 0u) { ev[7]  = interp_vertex(p7, p4, c7, c4, iso); }
    if ((edge_bits &  256u) != 0u) { ev[8]  = interp_vertex(p0, p4, c0, c4, iso); }
    if ((edge_bits &  512u) != 0u) { ev[9]  = interp_vertex(p1, p5, c1, c5, iso); }
    if ((edge_bits & 1024u) != 0u) { ev[10] = interp_vertex(p2, p6, c2, c6, iso); }
    if ((edge_bits & 2048u) != 0u) { ev[11] = interp_vertex(p3, p7, c3, c7, iso); }

    // Emit triangles from tri_table
    let tri_base = cube_index * 16u;
    var i: u32 = 0u;
    loop {
        if (i >= 16u) { break; }
        let e0_i = tri_table[tri_base + i];
        if (e0_i < 0) { break; }
        let e1_i = tri_table[tri_base + i + 1u];
        let e2_i = tri_table[tri_base + i + 2u];

        let v0 = ev[u32(e0_i)];
        let v1 = ev[u32(e1_i)];
        let v2 = ev[u32(e2_i)];

        // Face normal
        let u_vec = v1 - v0;
        let v_vec = v2 - v0;
        var n = cross(u_vec, v_vec);
        let len = length(n);
        if (len > 1e-10) {
            n = n / len;
        }

        // Atomically allocate 1 triangle slot
        let tri_idx = atomicAdd(&tri_counter[0], 1u);
        let out_base = tri_idx * 18u; // 3 verts × 6 floats (pos3 + normal3)

        // Vertex 0
        output[out_base +  0u] = v0.x;
        output[out_base +  1u] = v0.y;
        output[out_base +  2u] = v0.z;
        output[out_base +  3u] = n.x;
        output[out_base +  4u] = n.y;
        output[out_base +  5u] = n.z;
        // Vertex 1
        output[out_base +  6u] = v1.x;
        output[out_base +  7u] = v1.y;
        output[out_base +  8u] = v1.z;
        output[out_base +  9u] = n.x;
        output[out_base + 10u] = n.y;
        output[out_base + 11u] = n.z;
        // Vertex 2
        output[out_base + 12u] = v2.x;
        output[out_base + 13u] = v2.y;
        output[out_base + 14u] = v2.z;
        output[out_base + 15u] = n.x;
        output[out_base + 16u] = n.y;
        output[out_base + 17u] = n.z;

        i = i + 3u;
    }
}
"#;

// ─── Lookup tables ───────────────────────────────────────────────────────────

const EDGE_VERTICES: [(usize, usize); 12] = [
    (0, 1),
    (1, 2),
    (2, 3),
    (3, 0),
    (4, 5),
    (5, 6),
    (6, 7),
    (7, 4),
    (0, 4),
    (1, 5),
    (2, 6),
    (3, 7),
];

#[rustfmt::skip]
const EDGE_TABLE: [u16; 256] = [
    0x000, 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
    0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
    0x190, 0x099, 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
    0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
    0x230, 0x339, 0x033, 0x13a, 0x636, 0x73f, 0x435, 0x53c,
    0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
    0x3a0, 0x2a9, 0x1a3, 0x0aa, 0x7a6, 0x6af, 0x5a5, 0x4ac,
    0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
    0x460, 0x569, 0x663, 0x76a, 0x066, 0x16f, 0x265, 0x36c,
    0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
    0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0x0ff, 0x3f5, 0x2fc,
    0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
    0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x055, 0x15c,
    0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
    0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0x0cc,
    0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
    0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,
    0x0cc, 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
    0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,
    0x15c, 0x055, 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
    0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
    0x2fc, 0x3f5, 0x0ff, 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
    0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,
    0x36c, 0x265, 0x16f, 0x066, 0x76a, 0x663, 0x569, 0x460,
    0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
    0x4ac, 0x5a5, 0x6af, 0x7a6, 0x0aa, 0x1a3, 0x2a9, 0x3a0,
    0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c,
    0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x033, 0x339, 0x230,
    0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,
    0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x099, 0x190,
    0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
    0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x000,
];

/// Triangle table: edge triples for each cube configuration, terminated by -1.
const TRI_TABLE: [[i8; 16]; 256] = {
    let mut t = [[-1i8; 16]; 256];
    t[1] = [0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1];
    t[2] = [0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1];
    t[3] = [1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1];
    t[4] = [1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1];
    t[5] = [0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1];
    t[6] = [9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1];
    t[7] = [2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1];
    t[8] = [3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1];
    t[9] = [0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1];
    t[10] = [1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1];
    t[11] = [1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1];
    t[12] = [3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1];
    t[13] = [0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1];
    t[14] = [3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1];
    t[15] = [9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1];
    t
};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_marching_cubes_sphere() {
        let params = GridParams {
            origin: [-2.0, -2.0, -2.0],
            spacing: 0.5,
            dimensions: [9, 9, 9],
        };

        let mut values = vec![0.0; params.n_points()];
        let [nx, ny, nz] = params.dimensions;

        for ix in 0..nx {
            for iy in 0..ny {
                for iz in 0..nz {
                    let r = params.point(ix, iy, iz);
                    let r2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
                    values[params.flat_index(ix, iy, iz)] = 1.0 - r2;
                }
            }
        }

        let mesh = marching_cubes_cpu(&values, &params, 0.0);
        assert!(mesh.n_triangles > 0);
        assert_eq!(mesh.vertices.len(), mesh.n_triangles * 9);
        assert_eq!(mesh.normals.len(), mesh.n_triangles * 9);
    }

    #[test]
    fn test_marching_cubes_empty() {
        let params = GridParams {
            origin: [0.0, 0.0, 0.0],
            spacing: 1.0,
            dimensions: [3, 3, 3],
        };
        let values = vec![-1.0; params.n_points()];
        let mesh = marching_cubes_cpu(&values, &params, 0.0);
        assert_eq!(mesh.n_triangles, 0);
    }

    #[test]
    fn test_marching_cubes_all_inside() {
        let params = GridParams {
            origin: [0.0, 0.0, 0.0],
            spacing: 1.0,
            dimensions: [3, 3, 3],
        };
        let values = vec![1.0; params.n_points()];
        let mesh = marching_cubes_cpu(&values, &params, 0.0);
        assert_eq!(mesh.n_triangles, 0);
    }

    #[test]
    fn test_smooth_normals_empty() {
        let mut mesh = McOutput {
            vertices: vec![],
            normals: vec![],
            indices: vec![],
            n_triangles: 0,
        };
        smooth_mesh_normals(&mut mesh); // Should not panic
    }
}
