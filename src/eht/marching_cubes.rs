//! Marching Cubes isosurface extraction from volumetric grids.
//!
//! Given a VolumetricGrid and an isovalue, produces a triangle mesh
//! (vertices, normals, indices) suitable for WebGL/WebGPU rendering.

use super::volume::VolumetricGrid;
use serde::{Deserialize, Serialize};

/// A triangle mesh extracted from an isosurface.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IsosurfaceMesh {
    /// Flat array of vertex positions [x0,y0,z0, x1,y1,z1, ...] in Ångström.
    pub vertices: Vec<f32>,
    /// Flat array of vertex normals [nx0,ny0,nz0, ...].
    pub normals: Vec<f32>,
    /// Triangle indices (every 3 = one triangle).
    pub indices: Vec<u32>,
    /// Number of triangles.
    pub num_triangles: usize,
    /// The isovalue used.
    pub isovalue: f32,
}

impl IsosurfaceMesh {
    pub fn num_vertices(&self) -> usize {
        self.vertices.len() / 3
    }
}

/// Extract an isosurface from a volumetric grid using Marching Cubes.
pub fn marching_cubes(grid: &VolumetricGrid, isovalue: f32) -> IsosurfaceMesh {
    let [nx, ny, nz] = grid.dims;
    let mut vertices = Vec::new();
    let mut normals = Vec::new();
    let mut indices = Vec::new();

    if nx < 2 || ny < 2 || nz < 2 {
        return IsosurfaceMesh {
            vertices,
            normals,
            indices,
            num_triangles: 0,
            isovalue,
        };
    }

    for ix in 0..nx - 1 {
        for iy in 0..ny - 1 {
            for iz in 0..nz - 1 {
                process_cube(
                    grid,
                    ix,
                    iy,
                    iz,
                    isovalue,
                    &mut vertices,
                    &mut normals,
                    &mut indices,
                );
            }
        }
    }

    let num_triangles = indices.len() / 3;
    IsosurfaceMesh {
        vertices,
        normals,
        indices,
        num_triangles,
        isovalue,
    }
}

fn process_cube(
    grid: &VolumetricGrid,
    ix: usize,
    iy: usize,
    iz: usize,
    isovalue: f32,
    vertices: &mut Vec<f32>,
    normals: &mut Vec<f32>,
    indices: &mut Vec<u32>,
) {
    // 8 corners of the cube
    let corners = [
        (ix, iy, iz),
        (ix + 1, iy, iz),
        (ix + 1, iy + 1, iz),
        (ix, iy + 1, iz),
        (ix, iy, iz + 1),
        (ix + 1, iy, iz + 1),
        (ix + 1, iy + 1, iz + 1),
        (ix, iy + 1, iz + 1),
    ];

    let vals: [f32; 8] = [
        grid.values[grid.index(corners[0].0, corners[0].1, corners[0].2)],
        grid.values[grid.index(corners[1].0, corners[1].1, corners[1].2)],
        grid.values[grid.index(corners[2].0, corners[2].1, corners[2].2)],
        grid.values[grid.index(corners[3].0, corners[3].1, corners[3].2)],
        grid.values[grid.index(corners[4].0, corners[4].1, corners[4].2)],
        grid.values[grid.index(corners[5].0, corners[5].1, corners[5].2)],
        grid.values[grid.index(corners[6].0, corners[6].1, corners[6].2)],
        grid.values[grid.index(corners[7].0, corners[7].1, corners[7].2)],
    ];

    // Compute cube index from corner classifications
    let mut cube_idx = 0u8;
    for i in 0..8 {
        if vals[i] >= isovalue {
            cube_idx |= 1 << i;
        }
    }

    if cube_idx == 0 || cube_idx == 255 {
        return; // Entirely inside or outside
    }

    let edge_flags = EDGE_TABLE[cube_idx as usize];
    if edge_flags == 0 {
        return;
    }

    let positions: [[f64; 3]; 8] =
        core::array::from_fn(|i| grid.point_position(corners[i].0, corners[i].1, corners[i].2));

    // Interpolate edge vertices
    let mut edge_verts = [[0.0f64; 3]; 12];
    let edges = [
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

    for (ei, &(a, b)) in edges.iter().enumerate() {
        if edge_flags & (1 << ei) != 0 {
            edge_verts[ei] =
                interpolate_vertex(&positions[a], &positions[b], vals[a], vals[b], isovalue);
        }
    }

    // Emit triangles
    let tri = &TRI_TABLE[cube_idx as usize];
    let mut i = 0;
    while i < 16 && tri[i] != -1 {
        let base_idx = (vertices.len() / 3) as u32;

        for k in 0..3 {
            let edge = tri[i + k] as usize;
            let v = edge_verts[edge];
            vertices.push(v[0] as f32);
            vertices.push(v[1] as f32);
            vertices.push(v[2] as f32);

            // Approximate normal from gradient at vertex position
            let n = estimate_normal(grid, &v, isovalue);
            normals.push(n[0] as f32);
            normals.push(n[1] as f32);
            normals.push(n[2] as f32);
        }

        indices.push(base_idx);
        indices.push(base_idx + 1);
        indices.push(base_idx + 2);

        i += 3;
    }
}

fn interpolate_vertex(p1: &[f64; 3], p2: &[f64; 3], v1: f32, v2: f32, iso: f32) -> [f64; 3] {
    let diff = v2 - v1;
    let t = if diff.abs() < 1e-10 {
        0.5
    } else {
        ((iso - v1) / diff) as f64
    };
    [
        p1[0] + t * (p2[0] - p1[0]),
        p1[1] + t * (p2[1] - p1[1]),
        p1[2] + t * (p2[2] - p1[2]),
    ]
}

fn estimate_normal(grid: &VolumetricGrid, point: &[f64; 3], _iso: f32) -> [f64; 3] {
    let h = grid.spacing;
    let sample = |p: &[f64; 3]| -> f64 {
        // Nearest-neighbor sampling
        let ix = ((p[0] - grid.origin[0]) / h).round() as isize;
        let iy = ((p[1] - grid.origin[1]) / h).round() as isize;
        let iz = ((p[2] - grid.origin[2]) / h).round() as isize;
        if ix < 0
            || iy < 0
            || iz < 0
            || ix >= grid.dims[0] as isize
            || iy >= grid.dims[1] as isize
            || iz >= grid.dims[2] as isize
        {
            return 0.0;
        }
        grid.values[grid.index(ix as usize, iy as usize, iz as usize)] as f64
    };

    let gx =
        sample(&[point[0] + h, point[1], point[2]]) - sample(&[point[0] - h, point[1], point[2]]);
    let gy =
        sample(&[point[0], point[1] + h, point[2]]) - sample(&[point[0], point[1] - h, point[2]]);
    let gz =
        sample(&[point[0], point[1], point[2] + h]) - sample(&[point[0], point[1], point[2] - h]);

    let len = (gx * gx + gy * gy + gz * gz).sqrt();
    if len < 1e-12 {
        [0.0, 0.0, 1.0]
    } else {
        [-gx / len, -gy / len, -gz / len]
    }
}

/// Dual-phase isosurface result for orbital visualization.
///
/// Extracts both positive and negative lobes of an orbital (or ESP, etc.)
/// at the given absolute isovalue, producing two separate meshes.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DualPhaseMesh {
    /// Positive lobe (values >= +isovalue).
    pub positive: IsosurfaceMesh,
    /// Negative lobe (values <= -isovalue).
    pub negative: IsosurfaceMesh,
}

/// Extract dual-phase isosurfaces (positive and negative lobes).
///
/// Useful for orbital visualization where positive/negative phases
/// should be rendered in different colors.
pub fn marching_cubes_dual(grid: &VolumetricGrid, isovalue: f32) -> DualPhaseMesh {
    let positive = marching_cubes(grid, isovalue.abs());
    let negative = marching_cubes(grid, -isovalue.abs());
    DualPhaseMesh { positive, negative }
}

/// Mesh simplification: weld duplicate vertices and remove degenerate triangles.
///
/// Vertices closer than `weld_distance` are merged. Degenerate triangles
/// (with zero area or repeated vertex indices) are removed.
pub fn simplify_mesh(mesh: &IsosurfaceMesh, weld_distance: f32) -> IsosurfaceMesh {
    let n_verts = mesh.vertices.len() / 3;
    if n_verts == 0 {
        return mesh.clone();
    }

    let weld_dist_sq = weld_distance * weld_distance;

    // Build vertex map: for each vertex, find the canonical (first seen) vertex index
    let mut vertex_map: Vec<usize> = vec![0; n_verts];
    let mut new_vertices: Vec<[f32; 3]> = Vec::new();
    let mut new_normals: Vec<[f32; 3]> = Vec::new();
    let mut old_to_new: Vec<usize> = vec![0; n_verts];

    for i in 0..n_verts {
        let vx = mesh.vertices[3 * i];
        let vy = mesh.vertices[3 * i + 1];
        let vz = mesh.vertices[3 * i + 2];

        // Check if this vertex matches an existing one
        let mut found = None;
        for (j, v) in new_vertices.iter().enumerate() {
            let dx = vx - v[0];
            let dy = vy - v[1];
            let dz = vz - v[2];
            if dx * dx + dy * dy + dz * dz < weld_dist_sq {
                found = Some(j);
                break;
            }
        }

        if let Some(j) = found {
            old_to_new[i] = j;
        } else {
            old_to_new[i] = new_vertices.len();
            new_vertices.push([vx, vy, vz]);
            new_normals.push([
                mesh.normals[3 * i],
                mesh.normals[3 * i + 1],
                mesh.normals[3 * i + 2],
            ]);
        }
        vertex_map[i] = old_to_new[i];
    }

    // Rebuild indices with mapped vertices, removing degenerate triangles
    let mut new_indices = Vec::new();
    let n_tris = mesh.indices.len() / 3;
    for t in 0..n_tris {
        let i0 = vertex_map[mesh.indices[3 * t] as usize] as u32;
        let i1 = vertex_map[mesh.indices[3 * t + 1] as usize] as u32;
        let i2 = vertex_map[mesh.indices[3 * t + 2] as usize] as u32;

        // Skip degenerate triangles
        if i0 == i1 || i1 == i2 || i0 == i2 {
            continue;
        }

        // Check for zero-area triangle
        let v0 = new_vertices[i0 as usize];
        let v1 = new_vertices[i1 as usize];
        let v2 = new_vertices[i2 as usize];
        let e1 = [v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]];
        let e2 = [v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]];
        let cx = e1[1] * e2[2] - e1[2] * e2[1];
        let cy = e1[2] * e2[0] - e1[0] * e2[2];
        let cz = e1[0] * e2[1] - e1[1] * e2[0];
        let area_sq = cx * cx + cy * cy + cz * cz;
        if area_sq < 1e-12 {
            continue;
        }

        new_indices.push(i0);
        new_indices.push(i1);
        new_indices.push(i2);
    }

    let flat_verts: Vec<f32> = new_vertices
        .iter()
        .flat_map(|v| v.iter().copied())
        .collect();
    let flat_norms: Vec<f32> = new_normals.iter().flat_map(|n| n.iter().copied()).collect();

    IsosurfaceMesh {
        vertices: flat_verts,
        normals: flat_norms,
        num_triangles: new_indices.len() / 3,
        indices: new_indices,
        isovalue: mesh.isovalue,
    }
}

/// Compute angle-weighted vertex normals for smooth shading.
///
/// For each vertex, averages the face normals of adjacent triangles,
/// weighted by the angle at that vertex. This produces smoother shading
/// than per-face normals.
pub fn compute_angle_weighted_normals(mesh: &mut IsosurfaceMesh) {
    let n_verts = mesh.vertices.len() / 3;
    if n_verts == 0 {
        return;
    }

    let mut normals = vec![[0.0f32; 3]; n_verts];

    let n_tris = mesh.indices.len() / 3;
    for t in 0..n_tris {
        let idx = [
            mesh.indices[3 * t] as usize,
            mesh.indices[3 * t + 1] as usize,
            mesh.indices[3 * t + 2] as usize,
        ];

        let v: [[f32; 3]; 3] = [
            [
                mesh.vertices[3 * idx[0]],
                mesh.vertices[3 * idx[0] + 1],
                mesh.vertices[3 * idx[0] + 2],
            ],
            [
                mesh.vertices[3 * idx[1]],
                mesh.vertices[3 * idx[1] + 1],
                mesh.vertices[3 * idx[1] + 2],
            ],
            [
                mesh.vertices[3 * idx[2]],
                mesh.vertices[3 * idx[2] + 1],
                mesh.vertices[3 * idx[2] + 2],
            ],
        ];

        // Face normal
        let e1 = [v[1][0] - v[0][0], v[1][1] - v[0][1], v[1][2] - v[0][2]];
        let e2 = [v[2][0] - v[0][0], v[2][1] - v[0][1], v[2][2] - v[0][2]];
        let face_normal = [
            e1[1] * e2[2] - e1[2] * e2[1],
            e1[2] * e2[0] - e1[0] * e2[2],
            e1[0] * e2[1] - e1[1] * e2[0],
        ];

        // For each vertex in the triangle, compute the angle at that vertex
        for vi in 0..3 {
            let i0 = vi;
            let i1 = (vi + 1) % 3;
            let i2 = (vi + 2) % 3;

            let a = [
                v[i1][0] - v[i0][0],
                v[i1][1] - v[i0][1],
                v[i1][2] - v[i0][2],
            ];
            let b = [
                v[i2][0] - v[i0][0],
                v[i2][1] - v[i0][1],
                v[i2][2] - v[i0][2],
            ];

            let dot = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
            let la = (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt();
            let lb = (b[0] * b[0] + b[1] * b[1] + b[2] * b[2]).sqrt();

            let cos_angle = if la > 1e-8 && lb > 1e-8 {
                (dot / (la * lb)).clamp(-1.0, 1.0)
            } else {
                0.0
            };
            let angle = cos_angle.acos();

            // Weight the face normal by the angle at this vertex
            normals[idx[vi]][0] += angle * face_normal[0];
            normals[idx[vi]][1] += angle * face_normal[1];
            normals[idx[vi]][2] += angle * face_normal[2];
        }
    }

    // Normalize all vertex normals
    for i in 0..n_verts {
        let n = &mut normals[i];
        let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        if len > 1e-8 {
            n[0] /= len;
            n[1] /= len;
            n[2] /= len;
        } else {
            n[0] = 0.0;
            n[1] = 0.0;
            n[2] = 1.0;
        }
        mesh.normals[3 * i] = n[0];
        mesh.normals[3 * i + 1] = n[1];
        mesh.normals[3 * i + 2] = n[2];
    }
}

/// Ensure consistent outward-facing normal orientation.
///
/// Flips normals that point inward (based on the overall winding of the mesh).
pub fn flip_normals_outward(mesh: &mut IsosurfaceMesh) {
    let n_verts = mesh.vertices.len() / 3;
    if n_verts == 0 {
        return;
    }

    // Compute centroid of all vertices
    let mut cx = 0.0f32;
    let mut cy = 0.0f32;
    let mut cz = 0.0f32;
    for i in 0..n_verts {
        cx += mesh.vertices[3 * i];
        cy += mesh.vertices[3 * i + 1];
        cz += mesh.vertices[3 * i + 2];
    }
    cx /= n_verts as f32;
    cy /= n_verts as f32;
    cz /= n_verts as f32;

    // For each vertex, check if normal points away from centroid
    for i in 0..n_verts {
        let vx = mesh.vertices[3 * i] - cx;
        let vy = mesh.vertices[3 * i + 1] - cy;
        let vz = mesh.vertices[3 * i + 2] - cz;

        let dot =
            vx * mesh.normals[3 * i] + vy * mesh.normals[3 * i + 1] + vz * mesh.normals[3 * i + 2];

        // If normal points inward (dot < 0), flip it
        if dot < 0.0 {
            mesh.normals[3 * i] = -mesh.normals[3 * i];
            mesh.normals[3 * i + 1] = -mesh.normals[3 * i + 1];
            mesh.normals[3 * i + 2] = -mesh.normals[3 * i + 2];
        }
    }
}

/// Export mesh data in WASM-ready format: interleaved position+normal Float32Arrays.
///
/// Returns a flat array of [x, y, z, nx, ny, nz, x, y, z, nx, ny, nz, ...].
pub fn mesh_to_interleaved(mesh: &IsosurfaceMesh) -> Vec<f32> {
    let n_verts = mesh.vertices.len() / 3;
    let mut interleaved = Vec::with_capacity(n_verts * 6);
    for i in 0..n_verts {
        interleaved.push(mesh.vertices[3 * i]);
        interleaved.push(mesh.vertices[3 * i + 1]);
        interleaved.push(mesh.vertices[3 * i + 2]);
        interleaved.push(mesh.normals[3 * i]);
        interleaved.push(mesh.normals[3 * i + 1]);
        interleaved.push(mesh.normals[3 * i + 2]);
    }
    interleaved
}

// ─── Marching Cubes lookup tables ────────────────────────────────────────────
// Standard MC edge and triangle tables (Lorensen & Cline, 1987).

#[rustfmt::skip]
static EDGE_TABLE: [u16; 256] = [
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

#[rustfmt::skip]
static TRI_TABLE: [[i8; 16]; 256] = [
    [-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 0, 8, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 0, 1, 9,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 1, 8, 3, 9, 8, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 1, 2,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 0, 8, 3, 1, 2,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 9, 2,10, 0, 2, 9,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 2, 8, 3, 2,10, 8,10, 9, 8,-1,-1,-1,-1,-1,-1,-1],
    [ 3,11, 2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 0,11, 2, 8,11, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 1, 9, 0, 2, 3,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 1,11, 2, 1, 9,11, 9, 8,11,-1,-1,-1,-1,-1,-1,-1],
    [ 3,10, 1,11,10, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 0,10, 1, 0, 8,10, 8,11,10,-1,-1,-1,-1,-1,-1,-1],
    [ 3, 9, 0, 3,11, 9,11,10, 9,-1,-1,-1,-1,-1,-1,-1],
    [ 9, 8,10,10, 8,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 4, 7, 8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 4, 3, 0, 7, 3, 4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 0, 1, 9, 8, 4, 7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 4, 1, 9, 4, 7, 1, 7, 3, 1,-1,-1,-1,-1,-1,-1,-1],
    [ 1, 2,10, 8, 4, 7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 3, 4, 7, 3, 0, 4, 1, 2,10,-1,-1,-1,-1,-1,-1,-1],
    [ 9, 2,10, 9, 0, 2, 8, 4, 7,-1,-1,-1,-1,-1,-1,-1],
    [ 2,10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4,-1,-1,-1,-1],
    [ 8, 4, 7, 3,11, 2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [11, 4, 7,11, 2, 4, 2, 0, 4,-1,-1,-1,-1,-1,-1,-1],
    [ 9, 0, 1, 8, 4, 7, 2, 3,11,-1,-1,-1,-1,-1,-1,-1],
    [ 4, 7,11, 9, 4,11, 9,11, 2, 9, 2, 1,-1,-1,-1,-1],
    [ 3,10, 1, 3,11,10, 7, 8, 4,-1,-1,-1,-1,-1,-1,-1],
    [ 1,11,10, 1, 4,11, 1, 0, 4, 7,11, 4,-1,-1,-1,-1],
    [ 4, 7, 8, 9, 0,11, 9,11,10,11, 0, 3,-1,-1,-1,-1],
    [ 4, 7,11, 4,11, 9, 9,11,10,-1,-1,-1,-1,-1,-1,-1],
    [ 9, 5, 4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 9, 5, 4, 0, 8, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 0, 5, 4, 1, 5, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 8, 5, 4, 8, 3, 5, 3, 1, 5,-1,-1,-1,-1,-1,-1,-1],
    [ 1, 2,10, 9, 5, 4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 3, 0, 8, 1, 2,10, 4, 9, 5,-1,-1,-1,-1,-1,-1,-1],
    [ 5, 2,10, 5, 4, 2, 4, 0, 2,-1,-1,-1,-1,-1,-1,-1],
    [ 2,10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8,-1,-1,-1,-1],
    [ 9, 5, 4, 2, 3,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 0,11, 2, 0, 8,11, 4, 9, 5,-1,-1,-1,-1,-1,-1,-1],
    [ 0, 5, 4, 0, 1, 5, 2, 3,11,-1,-1,-1,-1,-1,-1,-1],
    [ 2, 1, 5, 2, 5, 8, 2, 8,11, 4, 8, 5,-1,-1,-1,-1],
    [10, 3,11,10, 1, 3, 9, 5, 4,-1,-1,-1,-1,-1,-1,-1],
    [ 4, 9, 5, 0, 8, 1, 8,10, 1, 8,11,10,-1,-1,-1,-1],
    [ 5, 4, 0, 5, 0,11, 5,11,10,11, 0, 3,-1,-1,-1,-1],
    [ 5, 4, 8, 5, 8,10,10, 8,11,-1,-1,-1,-1,-1,-1,-1],
    [ 9, 7, 8, 5, 7, 9,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 9, 3, 0, 9, 5, 3, 5, 7, 3,-1,-1,-1,-1,-1,-1,-1],
    [ 0, 7, 8, 0, 1, 7, 1, 5, 7,-1,-1,-1,-1,-1,-1,-1],
    [ 1, 5, 3, 3, 5, 7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 9, 7, 8, 9, 5, 7,10, 1, 2,-1,-1,-1,-1,-1,-1,-1],
    [10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3,-1,-1,-1,-1],
    [ 8, 0, 2, 8, 2, 5, 8, 5, 7,10, 5, 2,-1,-1,-1,-1],
    [ 2,10, 5, 2, 5, 3, 3, 5, 7,-1,-1,-1,-1,-1,-1,-1],
    [ 7, 9, 5, 7, 8, 9, 3,11, 2,-1,-1,-1,-1,-1,-1,-1],
    [ 9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7,11,-1,-1,-1,-1],
    [ 2, 3,11, 0, 1, 8, 1, 7, 8, 1, 5, 7,-1,-1,-1,-1],
    [11, 2, 1,11, 1, 7, 7, 1, 5,-1,-1,-1,-1,-1,-1,-1],
    [ 9, 5, 8, 8, 5, 7,10, 1, 3,10, 3,11,-1,-1,-1,-1],
    [ 5, 7, 0, 5, 0, 9, 7,11, 0, 1, 0,10,11,10, 0,-1],
    [11,10, 0,11, 0, 3,10, 5, 0, 8, 0, 7, 5, 7, 0,-1],
    [11,10, 5, 7,11, 5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [10, 6, 5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 0, 8, 3, 5,10, 6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 9, 0, 1, 5,10, 6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 1, 8, 3, 1, 9, 8, 5,10, 6,-1,-1,-1,-1,-1,-1,-1],
    [ 1, 6, 5, 2, 6, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 1, 6, 5, 1, 2, 6, 3, 0, 8,-1,-1,-1,-1,-1,-1,-1],
    [ 9, 6, 5, 9, 0, 6, 0, 2, 6,-1,-1,-1,-1,-1,-1,-1],
    [ 5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8,-1,-1,-1,-1],
    [ 2, 3,11,10, 6, 5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [11, 0, 8,11, 2, 0,10, 6, 5,-1,-1,-1,-1,-1,-1,-1],
    [ 0, 1, 9, 2, 3,11, 5,10, 6,-1,-1,-1,-1,-1,-1,-1],
    [ 5,10, 6, 1, 9, 2, 9,11, 2, 9, 8,11,-1,-1,-1,-1],
    [ 6, 3,11, 6, 5, 3, 5, 1, 3,-1,-1,-1,-1,-1,-1,-1],
    [ 0, 8,11, 0,11, 5, 0, 5, 1, 5,11, 6,-1,-1,-1,-1],
    [ 3,11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9,-1,-1,-1,-1],
    [ 6, 5, 9, 6, 9,11,11, 9, 8,-1,-1,-1,-1,-1,-1,-1],
    [ 5,10, 6, 4, 7, 8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 4, 3, 0, 4, 7, 3, 6, 5,10,-1,-1,-1,-1,-1,-1,-1],
    [ 1, 9, 0, 5,10, 6, 8, 4, 7,-1,-1,-1,-1,-1,-1,-1],
    [10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4,-1,-1,-1,-1],
    [ 6, 1, 2, 6, 5, 1, 4, 7, 8,-1,-1,-1,-1,-1,-1,-1],
    [ 1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7,-1,-1,-1,-1],
    [ 8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6,-1,-1,-1,-1],
    [ 7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9,-1],
    [ 3,11, 2, 7, 8, 4,10, 6, 5,-1,-1,-1,-1,-1,-1,-1],
    [ 5,10, 6, 4, 7, 2, 4, 2, 0, 2, 7,11,-1,-1,-1,-1],
    [ 0, 1, 9, 4, 7, 8, 2, 3,11, 5,10, 6,-1,-1,-1,-1],
    [ 9, 2, 1, 9,11, 2, 9, 4,11, 7,11, 4, 5,10, 6,-1],
    [ 8, 4, 7, 3,11, 5, 3, 5, 1, 5,11, 6,-1,-1,-1,-1],
    [ 5, 1,11, 5,11, 6, 1, 0,11, 7,11, 4, 0, 4,11,-1],
    [ 0, 5, 9, 0, 6, 5, 0, 3, 6,11, 6, 3, 8, 4, 7,-1],
    [ 6, 5, 9, 6, 9,11, 4, 7, 9, 7,11, 9,-1,-1,-1,-1],
    [10, 4, 9, 6, 4,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 4,10, 6, 4, 9,10, 0, 8, 3,-1,-1,-1,-1,-1,-1,-1],
    [10, 0, 1,10, 6, 0, 6, 4, 0,-1,-1,-1,-1,-1,-1,-1],
    [ 8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1,10,-1,-1,-1,-1],
    [ 1, 4, 9, 1, 2, 4, 2, 6, 4,-1,-1,-1,-1,-1,-1,-1],
    [ 3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4,-1,-1,-1,-1],
    [ 0, 2, 4, 4, 2, 6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 8, 3, 2, 8, 2, 4, 4, 2, 6,-1,-1,-1,-1,-1,-1,-1],
    [10, 4, 9,10, 6, 4,11, 2, 3,-1,-1,-1,-1,-1,-1,-1],
    [ 0, 8, 2, 2, 8,11, 4, 9,10, 4,10, 6,-1,-1,-1,-1],
    [ 3,11, 2, 0, 1, 6, 0, 6, 4, 6, 1,10,-1,-1,-1,-1],
    [ 6, 4, 1, 6, 1,10, 4, 8, 1, 2, 1,11, 8,11, 1,-1],
    [ 9, 6, 4, 9, 3, 6, 9, 1, 3,11, 6, 3,-1,-1,-1,-1],
    [ 8,11, 1, 8, 1, 0,11, 6, 1, 9, 1, 4, 6, 4, 1,-1],
    [ 3,11, 6, 3, 6, 0, 0, 6, 4,-1,-1,-1,-1,-1,-1,-1],
    [ 6, 4, 8,11, 6, 8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 7,10, 6, 7, 8,10, 8, 9,10,-1,-1,-1,-1,-1,-1,-1],
    [ 0, 7, 3, 0,10, 7, 0, 9,10, 6, 7,10,-1,-1,-1,-1],
    [10, 6, 7, 1,10, 7, 1, 7, 8, 1, 8, 0,-1,-1,-1,-1],
    [10, 6, 7,10, 7, 1, 1, 7, 3,-1,-1,-1,-1,-1,-1,-1],
    [ 1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7,-1,-1,-1,-1],
    [ 2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9,-1],
    [ 7, 8, 0, 7, 0, 6, 6, 0, 2,-1,-1,-1,-1,-1,-1,-1],
    [ 7, 3, 2, 6, 7, 2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 2, 3,11,10, 6, 8,10, 8, 9, 8, 6, 7,-1,-1,-1,-1],
    [ 2, 0, 7, 2, 7,11, 0, 9, 7, 6, 7,10, 9,10, 7,-1],
    [ 1, 8, 0, 1, 7, 8, 1,10, 7, 6, 7,10, 2, 3,11,-1],
    [11, 2, 1,11, 1, 7,10, 6, 1, 6, 7, 1,-1,-1,-1,-1],
    [ 8, 9, 6, 8, 6, 7, 9, 1, 6,11, 6, 3, 1, 3, 6,-1],
    [ 0, 9, 1,11, 6, 7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 7, 8, 0, 7, 0, 6, 3,11, 0,11, 6, 0,-1,-1,-1,-1],
    [ 7,11, 6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 7, 6,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 3, 0, 8,11, 7, 6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 0, 1, 9,11, 7, 6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 8, 1, 9, 8, 3, 1,11, 7, 6,-1,-1,-1,-1,-1,-1,-1],
    [10, 1, 2, 6,11, 7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 1, 2,10, 3, 0, 8, 6,11, 7,-1,-1,-1,-1,-1,-1,-1],
    [ 2, 9, 0, 2,10, 9, 6,11, 7,-1,-1,-1,-1,-1,-1,-1],
    [ 6,11, 7, 2,10, 3,10, 8, 3,10, 9, 8,-1,-1,-1,-1],
    [ 7, 2, 3, 6, 2, 7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 7, 0, 8, 7, 6, 0, 6, 2, 0,-1,-1,-1,-1,-1,-1,-1],
    [ 2, 7, 6, 2, 3, 7, 0, 1, 9,-1,-1,-1,-1,-1,-1,-1],
    [ 1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6,-1,-1,-1,-1],
    [10, 7, 6,10, 1, 7, 1, 3, 7,-1,-1,-1,-1,-1,-1,-1],
    [10, 7, 6, 1, 7,10, 1, 8, 7, 1, 0, 8,-1,-1,-1,-1],
    [ 0, 3, 7, 0, 7,10, 0,10, 9, 6,10, 7,-1,-1,-1,-1],
    [ 7, 6,10, 7,10, 8, 8,10, 9,-1,-1,-1,-1,-1,-1,-1],
    [ 6, 8, 4,11, 8, 6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 3, 6,11, 3, 0, 6, 0, 4, 6,-1,-1,-1,-1,-1,-1,-1],
    [ 8, 6,11, 8, 4, 6, 9, 0, 1,-1,-1,-1,-1,-1,-1,-1],
    [ 9, 4, 6, 9, 6, 3, 9, 3, 1,11, 3, 6,-1,-1,-1,-1],
    [ 6, 8, 4, 6,11, 8, 2,10, 1,-1,-1,-1,-1,-1,-1,-1],
    [ 1, 2,10, 3, 0,11, 0, 6,11, 0, 4, 6,-1,-1,-1,-1],
    [ 4,11, 8, 4, 6,11, 0, 2, 9, 2,10, 9,-1,-1,-1,-1],
    [10, 9, 3,10, 3, 2, 9, 4, 3,11, 3, 6, 4, 6, 3,-1],
    [ 8, 2, 3, 8, 4, 2, 4, 6, 2,-1,-1,-1,-1,-1,-1,-1],
    [ 0, 4, 2, 4, 6, 2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8,-1,-1,-1,-1],
    [ 1, 9, 4, 1, 4, 2, 2, 4, 6,-1,-1,-1,-1,-1,-1,-1],
    [ 8, 1, 3, 8, 6, 1, 8, 4, 6, 6,10, 1,-1,-1,-1,-1],
    [10, 1, 0,10, 0, 6, 6, 0, 4,-1,-1,-1,-1,-1,-1,-1],
    [ 4, 6, 3, 4, 3, 8, 6,10, 3, 0, 3, 9,10, 9, 3,-1],
    [10, 9, 4, 6,10, 4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 4, 9, 5, 7, 6,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 0, 8, 3, 4, 9, 5,11, 7, 6,-1,-1,-1,-1,-1,-1,-1],
    [ 5, 0, 1, 5, 4, 0, 7, 6,11,-1,-1,-1,-1,-1,-1,-1],
    [11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5,-1,-1,-1,-1],
    [ 9, 5, 4,10, 1, 2, 7, 6,11,-1,-1,-1,-1,-1,-1,-1],
    [ 6,11, 7, 1, 2,10, 0, 8, 3, 4, 9, 5,-1,-1,-1,-1],
    [ 7, 6,11, 5, 4,10, 4, 2,10, 4, 0, 2,-1,-1,-1,-1],
    [ 3, 4, 8, 3, 5, 4, 3, 2, 5,10, 5, 2,11, 7, 6,-1],
    [ 7, 2, 3, 7, 6, 2, 5, 4, 9,-1,-1,-1,-1,-1,-1,-1],
    [ 9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7,-1,-1,-1,-1],
    [ 3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0,-1,-1,-1,-1],
    [ 6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8,-1],
    [ 9, 5, 4,10, 1, 6, 1, 7, 6, 1, 3, 7,-1,-1,-1,-1],
    [ 1, 6,10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4,-1],
    [ 4, 0,10, 4,10, 5, 0, 3,10, 6,10, 7, 3, 7,10,-1],
    [ 7, 6,10, 7,10, 8, 5, 4,10, 4, 8,10,-1,-1,-1,-1],
    [ 6, 9, 5, 6,11, 9,11, 8, 9,-1,-1,-1,-1,-1,-1,-1],
    [ 3, 6,11, 0, 6, 3, 0, 5, 6, 0, 9, 5,-1,-1,-1,-1],
    [ 0,11, 8, 0, 5,11, 0, 1, 5, 5, 6,11,-1,-1,-1,-1],
    [ 6,11, 3, 6, 3, 5, 5, 3, 1,-1,-1,-1,-1,-1,-1,-1],
    [ 1, 2,10, 9, 5,11, 9,11, 8,11, 5, 6,-1,-1,-1,-1],
    [ 0,11, 3, 0, 6,11, 0, 9, 6, 5, 6, 9, 1, 2,10,-1],
    [11, 8, 5,11, 5, 6, 8, 0, 5,10, 5, 2, 0, 2, 5,-1],
    [ 6,11, 3, 6, 3, 5, 2,10, 3,10, 5, 3,-1,-1,-1,-1],
    [ 5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2,-1,-1,-1,-1],
    [ 9, 5, 6, 9, 6, 0, 0, 6, 2,-1,-1,-1,-1,-1,-1,-1],
    [ 1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8,-1],
    [ 1, 5, 6, 2, 1, 6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 1, 3, 6, 1, 6,10, 3, 8, 6, 5, 6, 9, 8, 9, 6,-1],
    [10, 1, 0,10, 0, 6, 9, 5, 0, 5, 6, 0,-1,-1,-1,-1],
    [ 0, 3, 8, 5, 6,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [10, 5, 6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [11, 5,10, 7, 5,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [11, 5,10,11, 7, 5, 8, 3, 0,-1,-1,-1,-1,-1,-1,-1],
    [ 5,11, 7, 5,10,11, 1, 9, 0,-1,-1,-1,-1,-1,-1,-1],
    [10, 7, 5,10,11, 7, 9, 8, 1, 8, 3, 1,-1,-1,-1,-1],
    [11, 1, 2,11, 7, 1, 7, 5, 1,-1,-1,-1,-1,-1,-1,-1],
    [ 0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2,11,-1,-1,-1,-1],
    [ 9, 7, 5, 9, 2, 7, 9, 0, 2, 2,11, 7,-1,-1,-1,-1],
    [ 7, 5, 2, 7, 2,11, 5, 9, 2, 3, 2, 8, 9, 8, 2,-1],
    [ 2, 5,10, 2, 3, 5, 3, 7, 5,-1,-1,-1,-1,-1,-1,-1],
    [ 8, 2, 0, 8, 5, 2, 8, 7, 5,10, 2, 5,-1,-1,-1,-1],
    [ 9, 0, 1, 5,10, 3, 5, 3, 7, 3,10, 2,-1,-1,-1,-1],
    [ 9, 8, 2, 9, 2, 1, 8, 7, 2,10, 2, 5, 7, 5, 2,-1],
    [ 1, 3, 5, 3, 7, 5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 0, 8, 7, 0, 7, 1, 1, 7, 5,-1,-1,-1,-1,-1,-1,-1],
    [ 9, 0, 3, 9, 3, 5, 5, 3, 7,-1,-1,-1,-1,-1,-1,-1],
    [ 9, 8, 7, 5, 9, 7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 5, 8, 4, 5,10, 8,10,11, 8,-1,-1,-1,-1,-1,-1,-1],
    [ 5, 0, 4, 5,11, 0, 5,10,11,11, 3, 0,-1,-1,-1,-1],
    [ 0, 1, 9, 8, 4,10, 8,10,11,10, 4, 5,-1,-1,-1,-1],
    [10,11, 4,10, 4, 5,11, 3, 4, 9, 4, 1, 3, 1, 4,-1],
    [ 2, 5, 1, 2, 8, 5, 2,11, 8, 4, 5, 8,-1,-1,-1,-1],
    [ 0, 4,11, 0,11, 3, 4, 5,11, 2,11, 1, 5, 1,11,-1],
    [ 0, 2, 5, 0, 5, 9, 2,11, 5, 4, 5, 8,11, 8, 5,-1],
    [ 9, 4, 5, 2,11, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 2, 5,10, 3, 5, 2, 3, 4, 5, 3, 8, 4,-1,-1,-1,-1],
    [ 5,10, 2, 5, 2, 4, 4, 2, 0,-1,-1,-1,-1,-1,-1,-1],
    [ 3,10, 2, 3, 5,10, 3, 8, 5, 4, 5, 8, 0, 1, 9,-1],
    [ 5,10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2,-1,-1,-1,-1],
    [ 8, 4, 5, 8, 5, 3, 3, 5, 1,-1,-1,-1,-1,-1,-1,-1],
    [ 0, 4, 5, 1, 0, 5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5,-1,-1,-1,-1],
    [ 9, 4, 5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 4,11, 7, 4, 9,11, 9,10,11,-1,-1,-1,-1,-1,-1,-1],
    [ 0, 8, 3, 4, 9, 7, 9,11, 7, 9,10,11,-1,-1,-1,-1],
    [ 1,10,11, 1,11, 4, 1, 4, 0, 7, 4,11,-1,-1,-1,-1],
    [ 3, 1, 4, 3, 4, 8, 1,10, 4, 7, 4,11,10,11, 4,-1],
    [ 4,11, 7, 9,11, 4, 9, 2,11, 9, 1, 2,-1,-1,-1,-1],
    [ 9, 7, 4, 9,11, 7, 9, 1,11, 2,11, 1, 0, 8, 3,-1],
    [11, 7, 4,11, 4, 2, 2, 4, 0,-1,-1,-1,-1,-1,-1,-1],
    [11, 7, 4,11, 4, 2, 8, 3, 4, 3, 2, 4,-1,-1,-1,-1],
    [ 2, 9,10, 2, 7, 9, 2, 3, 7, 7, 4, 9,-1,-1,-1,-1],
    [ 9,10, 7, 9, 7, 4,10, 2, 7, 8, 7, 0, 2, 0, 7,-1],
    [ 3, 7,10, 3,10, 2, 7, 4,10, 1,10, 0, 4, 0,10,-1],
    [ 1,10, 2, 8, 7, 4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 4, 9, 1, 4, 1, 7, 7, 1, 3,-1,-1,-1,-1,-1,-1,-1],
    [ 4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1,-1,-1,-1,-1],
    [ 4, 0, 3, 7, 4, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 4, 8, 7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 9,10, 8,10,11, 8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 3, 0, 9, 3, 9,11,11, 9,10,-1,-1,-1,-1,-1,-1,-1],
    [ 0, 1,10, 0,10, 8, 8,10,11,-1,-1,-1,-1,-1,-1,-1],
    [ 3, 1,10,11, 3,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 1, 2,11, 1,11, 9, 9,11, 8,-1,-1,-1,-1,-1,-1,-1],
    [ 3, 0, 9, 3, 9,11, 1, 2, 9, 2,11, 9,-1,-1,-1,-1],
    [ 0, 2,11, 8, 0,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 3, 2,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 2, 3, 8, 2, 8,10,10, 8, 9,-1,-1,-1,-1,-1,-1,-1],
    [ 9,10, 2, 0, 9, 2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 2, 3, 8, 2, 8,10, 0, 1, 8, 1,10, 8,-1,-1,-1,-1],
    [ 1,10, 2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 1, 3, 8, 9, 1, 8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 0, 9, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [ 0, 3, 8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    [-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
];

#[cfg(test)]
mod tests {
    use super::super::volume::VolumetricGrid;
    use super::*;

    fn sphere_grid(radius: f64, spacing: f64) -> VolumetricGrid {
        let extent = radius + 1.0;
        let n = (2.0 * extent / spacing).ceil() as usize + 1;
        let origin = [-extent, -extent, -extent];
        let mut values = vec![0.0f32; n * n * n];
        for ix in 0..n {
            for iy in 0..n {
                for iz in 0..n {
                    let x = origin[0] + ix as f64 * spacing;
                    let y = origin[1] + iy as f64 * spacing;
                    let z = origin[2] + iz as f64 * spacing;
                    let r = (x * x + y * y + z * z).sqrt();
                    values[ix * n * n + iy * n + iz] = (radius - r) as f32;
                }
            }
        }
        VolumetricGrid {
            origin,
            spacing,
            dims: [n, n, n],
            values,
        }
    }

    #[test]
    fn test_marching_cubes_sphere() {
        let grid = sphere_grid(2.0, 0.3);
        let mesh = marching_cubes(&grid, 0.0);
        assert!(
            mesh.num_triangles > 50,
            "sphere should produce many triangles, got {}",
            mesh.num_triangles,
        );
        // All indices should be valid
        let nv = mesh.num_vertices() as u32;
        for &idx in &mesh.indices {
            assert!(idx < nv, "index {} out of bound (nv={})", idx, nv);
        }
    }

    #[test]
    fn test_marching_cubes_empty() {
        let grid = VolumetricGrid {
            origin: [0.0; 3],
            spacing: 1.0,
            dims: [3, 3, 3],
            values: vec![0.0f32; 27],
        };
        let mesh = marching_cubes(&grid, 1.0);
        assert_eq!(mesh.num_triangles, 0);
    }

    #[test]
    fn test_normals_same_count_as_vertices() {
        let grid = sphere_grid(2.0, 0.3);
        let mesh = marching_cubes(&grid, 0.0);
        assert_eq!(mesh.vertices.len(), mesh.normals.len());
    }

    #[test]
    fn test_dual_phase_isosurface() {
        // Scalar field: positive inside sphere -> positive lobe,
        // negative outside -> negative lobe (with offset)
        let grid = sphere_grid(2.0, 0.3);
        let dual = marching_cubes_dual(&grid, 0.5);
        assert!(
            dual.positive.num_triangles > 0,
            "positive lobe should have triangles"
        );
        // The negative phase at -0.5 also produces a surface
        assert!(
            dual.negative.num_triangles > 0,
            "negative lobe should have triangles"
        );
    }

    #[test]
    fn test_simplify_mesh_welds_vertices() {
        let grid = sphere_grid(2.0, 0.3);
        let mesh = marching_cubes(&grid, 0.0);
        let before = mesh.num_vertices();
        let simplified = simplify_mesh(&mesh, 0.01);
        // Welding should reduce vertex count
        assert!(
            simplified.num_vertices() <= before,
            "simplified should have <= vertices: {} vs {}",
            simplified.num_vertices(),
            before
        );
        assert!(simplified.num_triangles > 0, "should still have triangles");
    }

    #[test]
    fn test_angle_weighted_normals_sphere() {
        let grid = sphere_grid(2.0, 0.3);
        let mut mesh = marching_cubes(&grid, 0.0);
        compute_angle_weighted_normals(&mut mesh);

        // For a sphere centered at origin, normals should point radially outward
        // Check cosine similarity between vertex position (= radial direction) and normal
        let n_verts = mesh.num_vertices();
        let mut cos_sum = 0.0f64;
        let mut count = 0;
        for i in 0..n_verts {
            let vx = mesh.vertices[3 * i] as f64;
            let vy = mesh.vertices[3 * i + 1] as f64;
            let vz = mesh.vertices[3 * i + 2] as f64;
            let vlen = (vx * vx + vy * vy + vz * vz).sqrt();
            if vlen < 1e-6 {
                continue;
            }
            let nx = mesh.normals[3 * i] as f64;
            let ny = mesh.normals[3 * i + 1] as f64;
            let nz = mesh.normals[3 * i + 2] as f64;
            let nlen = (nx * nx + ny * ny + nz * nz).sqrt();
            if nlen < 1e-6 {
                continue;
            }
            let cos = (vx * nx + vy * ny + vz * nz) / (vlen * nlen);
            cos_sum += cos.abs();
            count += 1;
        }
        let avg_cos = cos_sum / count as f64;
        assert!(
            avg_cos > 0.90,
            "angle-weighted sphere normals should align with radial direction, got avg cos = {:.3}",
            avg_cos
        );
    }

    #[test]
    fn test_flip_normals_outward_sphere() {
        let grid = sphere_grid(2.0, 0.3);
        let mut mesh = marching_cubes(&grid, 0.0);
        flip_normals_outward(&mut mesh);

        // After flipping, all normals should point away from centroid (≈ origin)
        let n_verts = mesh.num_vertices();
        let mut inward = 0;
        for i in 0..n_verts {
            let dot = mesh.vertices[3 * i] * mesh.normals[3 * i]
                + mesh.vertices[3 * i + 1] * mesh.normals[3 * i + 1]
                + mesh.vertices[3 * i + 2] * mesh.normals[3 * i + 2];
            if dot < 0.0 {
                inward += 1;
            }
        }
        assert_eq!(inward, 0, "no normals should point inward after flip");
    }

    #[test]
    fn test_mesh_to_interleaved() {
        let grid = sphere_grid(2.0, 0.3);
        let mesh = marching_cubes(&grid, 0.0);
        let interleaved = mesh_to_interleaved(&mesh);
        let n_verts = mesh.num_vertices();
        assert_eq!(interleaved.len(), n_verts * 6, "6 floats per vertex");

        // Verify interleaving
        for i in 0..n_verts.min(10) {
            assert_eq!(interleaved[6 * i], mesh.vertices[3 * i]);
            assert_eq!(interleaved[6 * i + 1], mesh.vertices[3 * i + 1]);
            assert_eq!(interleaved[6 * i + 2], mesh.vertices[3 * i + 2]);
            assert_eq!(interleaved[6 * i + 3], mesh.normals[3 * i]);
            assert_eq!(interleaved[6 * i + 4], mesh.normals[3 * i + 1]);
            assert_eq!(interleaved[6 * i + 5], mesh.normals[3 * i + 2]);
        }
    }
}
