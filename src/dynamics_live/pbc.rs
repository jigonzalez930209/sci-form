//! Periodic boundary conditions for molecular dynamics.
//!
//! Implements minimum-image convention and coordinate wrapping for
//! orthorhombic and triclinic simulation boxes.

use serde::{Deserialize, Serialize};

/// Simulation box for periodic boundary conditions.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SimulationBox {
    /// Lattice vectors as rows: [[ax,ay,az],[bx,by,bz],[cx,cy,cz]].
    pub lattice: [[f64; 3]; 3],
    /// Inverse lattice matrix for fractional coordinate conversion.
    pub inv_lattice: [[f64; 3]; 3],
    /// Whether PBC is active.
    pub periodic: bool,
}

impl SimulationBox {
    /// Create an orthorhombic box with side lengths (Å).
    pub fn orthorhombic(lx: f64, ly: f64, lz: f64) -> Self {
        let lattice = [[lx, 0.0, 0.0], [0.0, ly, 0.0], [0.0, 0.0, lz]];
        let inv_lattice = [
            [1.0 / lx, 0.0, 0.0],
            [0.0, 1.0 / ly, 0.0],
            [0.0, 0.0, 1.0 / lz],
        ];
        Self {
            lattice,
            inv_lattice,
            periodic: true,
        }
    }

    /// Create a general triclinic box from lattice vectors.
    pub fn triclinic(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> Self {
        let lattice = [a, b, c];
        let inv_lattice = invert_3x3(&lattice);
        Self {
            lattice,
            inv_lattice,
            periodic: true,
        }
    }

    /// No periodic boundaries (isolated molecule).
    pub fn none() -> Self {
        Self {
            lattice: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            inv_lattice: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            periodic: false,
        }
    }

    /// Wrap a coordinate into the primary cell using fractional coordinates.
    pub fn wrap(&self, pos: &mut [f64; 3]) {
        if !self.periodic {
            return;
        }
        // Convert to fractional
        let frac = self.to_fractional(pos);
        let wrapped = [
            frac[0] - frac[0].floor(),
            frac[1] - frac[1].floor(),
            frac[2] - frac[2].floor(),
        ];
        // Convert back to Cartesian
        *pos = self.to_cartesian(&wrapped);
    }

    /// Minimum-image displacement vector from r1 to r2.
    ///
    /// The closest periodic image is located by searching neighbouring
    /// lattice images. Rounding the fractional components alone is only
    /// exact for orthogonal boxes (there the components are independent),
    /// so skewed boxes need the explicit image search.
    pub fn minimum_image(&self, r1: &[f64; 3], r2: &[f64; 3]) -> [f64; 3] {
        if !self.periodic {
            return [r2[0] - r1[0], r2[1] - r1[1], r2[2] - r1[2]];
        }

        let dr = [r2[0] - r1[0], r2[1] - r1[1], r2[2] - r1[2]];
        let frac = self.to_fractional(&dr);
        // Fractional rounding: the exact minimum image for orthogonal boxes.
        let base = [frac[0].round(), frac[1].round(), frac[2].round()];
        let mut best = [frac[0] - base[0], frac[1] - base[1], frac[2] - base[2]];
        if !self.is_orthogonal() {
            // Skewed box: search the ±1-cell neighbourhood, then refine by
            // stepping while a neighbouring image is closer.
            let mut best_d2 = self.frac_dist2(best);
            let mut shift = base;
            for i in -1..=1 {
                for j in -1..=1 {
                    for k in -1..=1 {
                        let candidate_shift =
                            [base[0] + i as f64, base[1] + j as f64, base[2] + k as f64];
                        let candidate = [
                            frac[0] - candidate_shift[0],
                            frac[1] - candidate_shift[1],
                            frac[2] - candidate_shift[2],
                        ];
                        let d2 = self.frac_dist2(candidate);
                        if d2 < best_d2 {
                            best_d2 = d2;
                            best = candidate;
                            shift = candidate_shift;
                        }
                    }
                }
            }
            for _ in 0..64 {
                let mut improved = false;
                for axis in 0..3 {
                    for step in [-1.0f64, 1.0] {
                        let mut candidate_shift = shift;
                        candidate_shift[axis] += step;
                        let candidate = [
                            frac[0] - candidate_shift[0],
                            frac[1] - candidate_shift[1],
                            frac[2] - candidate_shift[2],
                        ];
                        let d2 = self.frac_dist2(candidate);
                        if d2 < best_d2 {
                            best_d2 = d2;
                            best = candidate;
                            shift = candidate_shift;
                            improved = true;
                        }
                    }
                }
                if !improved {
                    break;
                }
            }
        }
        self.to_cartesian(&best)
    }

    /// Length of the shortest lattice vector (Å).
    ///
    /// Used to decide whether a cell list with edge `r` can be searched with
    /// only the 27 neighbouring cells: that is valid when the shortest box
    /// vector is at least `2 r` long.
    pub fn shortest_vector(&self) -> f64 {
        let l = &self.lattice;
        let norm = |v: &[f64; 3]| (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        norm(&l[0]).min(norm(&l[1])).min(norm(&l[2]))
    }

    /// Squared Cartesian length of a fractional displacement.
    fn frac_dist2(&self, frac: [f64; 3]) -> f64 {
        let d = self.to_cartesian(&frac);
        d[0] * d[0] + d[1] * d[1] + d[2] * d[2]
    }

    /// Whether the lattice vectors are mutually orthogonal.
    fn is_orthogonal(&self) -> bool {
        let l = self.lattice;
        for i in 0..3 {
            for j in (i + 1)..3 {
                let dot = l[i][0] * l[j][0] + l[i][1] * l[j][1] + l[i][2] * l[j][2];
                let scale = norm3(l[i]) * norm3(l[j]);
                if dot.abs() > 1e-12 * scale {
                    return false;
                }
            }
        }
        true
    }

    /// Minimum-image distance squared.
    pub fn distance_sq(&self, r1: &[f64; 3], r2: &[f64; 3]) -> f64 {
        let dr = self.minimum_image(r1, r2);
        dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]
    }

    /// Convert Cartesian → fractional coordinates.
    fn to_fractional(&self, cart: &[f64; 3]) -> [f64; 3] {
        let m = &self.inv_lattice;
        [
            m[0][0] * cart[0] + m[0][1] * cart[1] + m[0][2] * cart[2],
            m[1][0] * cart[0] + m[1][1] * cart[1] + m[1][2] * cart[2],
            m[2][0] * cart[0] + m[2][1] * cart[1] + m[2][2] * cart[2],
        ]
    }

    /// Convert fractional → Cartesian coordinates.
    fn to_cartesian(&self, frac: &[f64; 3]) -> [f64; 3] {
        let m = &self.lattice;
        [
            m[0][0] * frac[0] + m[1][0] * frac[1] + m[2][0] * frac[2],
            m[0][1] * frac[0] + m[1][1] * frac[1] + m[2][1] * frac[2],
            m[0][2] * frac[0] + m[1][2] * frac[1] + m[2][2] * frac[2],
        ]
    }

    /// Box volume in ų.
    pub fn volume(&self) -> f64 {
        let a = self.lattice[0];
        let b = self.lattice[1];
        let c = self.lattice[2];
        // Volume = |a · (b × c)|
        let cross = [
            b[1] * c[2] - b[2] * c[1],
            b[2] * c[0] - b[0] * c[2],
            b[0] * c[1] - b[1] * c[0],
        ];
        (a[0] * cross[0] + a[1] * cross[1] + a[2] * cross[2]).abs()
    }
}

/// Euclidean length of a 3-vector.
fn norm3(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

/// Invert a 3×3 matrix.
fn invert_3x3(m: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    let inv_det = 1.0 / det;
    [
        [
            (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * inv_det,
            (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * inv_det,
            (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * inv_det,
        ],
        [
            (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * inv_det,
            (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * inv_det,
            (m[0][2] * m[1][0] - m[0][0] * m[1][2]) * inv_det,
        ],
        [
            (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * inv_det,
            (m[0][1] * m[2][0] - m[0][0] * m[2][1]) * inv_det,
            (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * inv_det,
        ],
    ]
}

/// Wrap all atom coordinates in a flat buffer into the primary cell.
pub fn wrap_all_positions(positions_flat: &mut [f64], sim_box: &SimulationBox) {
    if !sim_box.periodic {
        return;
    }
    let n = positions_flat.len() / 3;
    for i in 0..n {
        let mut pos = [
            positions_flat[3 * i],
            positions_flat[3 * i + 1],
            positions_flat[3 * i + 2],
        ];
        sim_box.wrap(&mut pos);
        positions_flat[3 * i] = pos[0];
        positions_flat[3 * i + 1] = pos[1];
        positions_flat[3 * i + 2] = pos[2];
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn orthorhombic_wrap() {
        let b = SimulationBox::orthorhombic(10.0, 10.0, 10.0);
        let mut pos = [11.5, -1.0, 25.3];
        b.wrap(&mut pos);
        assert!((pos[0] - 1.5).abs() < 1e-10);
        assert!((pos[1] - 9.0).abs() < 1e-10);
        assert!((pos[2] - 5.3).abs() < 1e-10);
    }

    #[test]
    fn minimum_image_distance() {
        let b = SimulationBox::orthorhombic(10.0, 10.0, 10.0);
        let r1 = [1.0, 1.0, 1.0];
        let r2 = [9.0, 1.0, 1.0];
        let dr = b.minimum_image(&r1, &r2);
        assert!((dr[0] - (-2.0)).abs() < 1e-10);
    }
}
