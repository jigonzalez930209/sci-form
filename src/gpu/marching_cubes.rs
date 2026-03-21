//! Marching cubes isosurface extraction from 3D scalar fields.
//!
//! Converts a scalar field (orbital ψ or density ρ) into a triangle mesh
//! at a specified isovalue. Reference: Lorensen & Cline, SIGGRAPH 1987.

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

// ─── Lookup tables ───────────────────────────────────────────────────────────

const EDGE_VERTICES: [(usize, usize); 12] = [
    (0, 1), (1, 2), (2, 3), (3, 0),
    (4, 5), (5, 6), (6, 7), (7, 4),
    (0, 4), (1, 5), (2, 6), (3, 7),
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
    t[1]  = [0,8,3, -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1];
    t[2]  = [0,1,9, -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1];
    t[3]  = [1,8,3, 9,8,1, -1,-1,-1,-1,-1,-1,-1,-1,-1,-1];
    t[4]  = [1,2,10, -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1];
    t[5]  = [0,8,3, 1,2,10, -1,-1,-1,-1,-1,-1,-1,-1,-1,-1];
    t[6]  = [9,2,10, 0,2,9, -1,-1,-1,-1,-1,-1,-1,-1,-1,-1];
    t[7]  = [2,8,3, 2,10,8, 10,9,8, -1,-1,-1,-1,-1,-1,-1];
    t[8]  = [3,11,2, -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1];
    t[9]  = [0,11,2, 8,11,0, -1,-1,-1,-1,-1,-1,-1,-1,-1,-1];
    t[10] = [1,9,0, 2,3,11, -1,-1,-1,-1,-1,-1,-1,-1,-1,-1];
    t[11] = [1,11,2, 1,9,11, 9,8,11, -1,-1,-1,-1,-1,-1,-1];
    t[12] = [3,10,1, 11,10,3, -1,-1,-1,-1,-1,-1,-1,-1,-1,-1];
    t[13] = [0,10,1, 0,8,10, 8,11,10, -1,-1,-1,-1,-1,-1,-1];
    t[14] = [3,9,0, 3,11,9, 11,10,9, -1,-1,-1,-1,-1,-1,-1];
    t[15] = [9,8,10, 10,8,11, -1,-1,-1,-1,-1,-1,-1,-1,-1,-1];
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
