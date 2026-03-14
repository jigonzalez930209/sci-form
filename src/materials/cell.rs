//! Unit cell definition with lattice vectors and periodic boundary conditions.
//!
//! A crystallographic unit cell is defined by three lattice vectors **a**, **b**, **c**
//! (or equivalently six parameters: a, b, c, α, β, γ).
//!
//! Fractional ↔ Cartesian conversions:
//!   r_cart = M · r_frac    where M = [a | b | c] column matrix
//!   r_frac = M⁻¹ · r_cart

use serde::{Deserialize, Serialize};

/// A 3D periodic unit cell defined by lattice vectors.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct UnitCell {
    /// Lattice vectors as rows: \[\[ax,ay,az\], \[bx,by,bz\], \[cx,cy,cz\]\].
    pub lattice: [[f64; 3]; 3],
}

/// Cell parameters in crystallographic notation (a, b, c, α, β, γ).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CellParameters {
    /// Lattice lengths (Å).
    pub a: f64,
    pub b: f64,
    pub c: f64,
    /// Lattice angles (degrees).
    pub alpha: f64,
    pub beta: f64,
    pub gamma: f64,
}

impl UnitCell {
    /// Create a unit cell from three lattice vectors (row vectors).
    pub fn new(lattice: [[f64; 3]; 3]) -> Self {
        Self { lattice }
    }

    /// Create a cubic unit cell with edge length `a`.
    pub fn cubic(a: f64) -> Self {
        Self {
            lattice: [[a, 0.0, 0.0], [0.0, a, 0.0], [0.0, 0.0, a]],
        }
    }

    /// Create a cell from crystallographic parameters (a, b, c, α, β, γ in degrees).
    pub fn from_parameters(params: &CellParameters) -> Self {
        let alpha = params.alpha.to_radians();
        let beta = params.beta.to_radians();
        let gamma = params.gamma.to_radians();

        let cos_alpha = alpha.cos();
        let cos_beta = beta.cos();
        let cos_gamma = gamma.cos();
        let sin_gamma = gamma.sin();

        // Standard crystallographic convention: a along x, b in xy-plane
        let ax = params.a;
        let bx = params.b * cos_gamma;
        let by = params.b * sin_gamma;
        let cx = params.c * cos_beta;
        let cy = params.c * (cos_alpha - cos_beta * cos_gamma) / sin_gamma;
        let cz = (params.c * params.c - cx * cx - cy * cy).sqrt();

        Self {
            lattice: [[ax, 0.0, 0.0], [bx, by, 0.0], [cx, cy, cz]],
        }
    }

    /// Extract crystallographic parameters from lattice vectors.
    pub fn parameters(&self) -> CellParameters {
        let a_vec = self.lattice[0];
        let b_vec = self.lattice[1];
        let c_vec = self.lattice[2];

        let a = norm3(a_vec);
        let b = norm3(b_vec);
        let c = norm3(c_vec);

        let alpha = angle_between(b_vec, c_vec).to_degrees();
        let beta = angle_between(a_vec, c_vec).to_degrees();
        let gamma = angle_between(a_vec, b_vec).to_degrees();

        CellParameters {
            a,
            b,
            c,
            alpha,
            beta,
            gamma,
        }
    }

    /// Unit cell volume: |a · (b × c)|.
    pub fn volume(&self) -> f64 {
        let a = self.lattice[0];
        let b = self.lattice[1];
        let c = self.lattice[2];
        // b × c
        let bxc = cross3(b, c);
        dot3(a, bxc).abs()
    }

    /// Convert fractional coordinates to Cartesian (Å).
    /// r_cart = f0 * a + f1 * b + f2 * c
    pub fn frac_to_cart(&self, frac: [f64; 3]) -> [f64; 3] {
        let a = self.lattice[0];
        let b = self.lattice[1];
        let c = self.lattice[2];
        [
            frac[0] * a[0] + frac[1] * b[0] + frac[2] * c[0],
            frac[0] * a[1] + frac[1] * b[1] + frac[2] * c[1],
            frac[0] * a[2] + frac[1] * b[2] + frac[2] * c[2],
        ]
    }

    /// Convert Cartesian coordinates (Å) to fractional.
    /// frac_to_cart computes: r = f[0]*a + f[1]*b + f[2]*c = M^T · f
    /// So: f = (M^T)^{-1} · r = (M^{-1})^T · r
    pub fn cart_to_frac(&self, cart: [f64; 3]) -> [f64; 3] {
        let inv = self.inverse_matrix();
        // We need (M^{-1})^T applied to cart, so use columns of inv as rows
        [
            inv[0][0] * cart[0] + inv[1][0] * cart[1] + inv[2][0] * cart[2],
            inv[0][1] * cart[0] + inv[1][1] * cart[1] + inv[2][1] * cart[2],
            inv[0][2] * cart[0] + inv[1][2] * cart[1] + inv[2][2] * cart[2],
        ]
    }

    /// Wrap fractional coordinates into [0, 1) — periodic boundary conditions.
    pub fn wrap_frac(frac: [f64; 3]) -> [f64; 3] {
        [
            frac[0] - frac[0].floor(),
            frac[1] - frac[1].floor(),
            frac[2] - frac[2].floor(),
        ]
    }

    /// Wrap Cartesian coordinates into the unit cell.
    pub fn wrap_cart(&self, cart: [f64; 3]) -> [f64; 3] {
        let frac = self.cart_to_frac(cart);
        let wrapped = Self::wrap_frac(frac);
        self.frac_to_cart(wrapped)
    }

    /// Minimum-image distance between two Cartesian points under PBC.
    pub fn minimum_image_distance(&self, a: [f64; 3], b: [f64; 3]) -> f64 {
        let fa = self.cart_to_frac(a);
        let fb = self.cart_to_frac(b);
        let mut df = [fa[0] - fb[0], fa[1] - fb[1], fa[2] - fb[2]];
        // Apply minimum image convention
        for d in &mut df {
            *d -= d.round();
        }
        let dc = self.frac_to_cart(df);
        norm3(dc)
    }

    /// Build a supercell by replicating (na × nb × nc) times.
    /// Returns new lattice and list of fractional offsets for each image.
    pub fn supercell(&self, na: usize, nb: usize, nc: usize) -> (UnitCell, Vec<[f64; 3]>) {
        let new_lattice = [
            [
                self.lattice[0][0] * na as f64,
                self.lattice[0][1] * na as f64,
                self.lattice[0][2] * na as f64,
            ],
            [
                self.lattice[1][0] * nb as f64,
                self.lattice[1][1] * nb as f64,
                self.lattice[1][2] * nb as f64,
            ],
            [
                self.lattice[2][0] * nc as f64,
                self.lattice[2][1] * nc as f64,
                self.lattice[2][2] * nc as f64,
            ],
        ];

        let mut offsets = Vec::with_capacity(na * nb * nc);
        for ia in 0..na {
            for ib in 0..nb {
                for ic in 0..nc {
                    offsets.push([ia as f64, ib as f64, ic as f64]);
                }
            }
        }

        (UnitCell::new(new_lattice), offsets)
    }

    /// Inverse of the 3×3 matrix whose rows are the lattice vectors.
    fn inverse_matrix(&self) -> [[f64; 3]; 3] {
        let m = self.lattice;
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
}

fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn norm3(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

fn angle_between(a: [f64; 3], b: [f64; 3]) -> f64 {
    let cos_angle = dot3(a, b) / (norm3(a) * norm3(b));
    cos_angle.clamp(-1.0, 1.0).acos()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cubic_cell_volume() {
        let cell = UnitCell::cubic(10.0);
        assert!((cell.volume() - 1000.0).abs() < 1e-10);
    }

    #[test]
    fn test_cubic_parameters() {
        let cell = UnitCell::cubic(5.0);
        let p = cell.parameters();
        assert!((p.a - 5.0).abs() < 1e-10);
        assert!((p.b - 5.0).abs() < 1e-10);
        assert!((p.c - 5.0).abs() < 1e-10);
        assert!((p.alpha - 90.0).abs() < 1e-10);
        assert!((p.beta - 90.0).abs() < 1e-10);
        assert!((p.gamma - 90.0).abs() < 1e-10);
    }

    #[test]
    fn test_frac_cart_roundtrip() {
        let cell = UnitCell::from_parameters(&CellParameters {
            a: 10.0, b: 12.0, c: 8.0,
            alpha: 90.0, beta: 90.0, gamma: 120.0,
        });
        let frac = [0.3, 0.4, 0.7];
        let cart = cell.frac_to_cart(frac);
        let back = cell.cart_to_frac(cart);
        for i in 0..3 {
            assert!((frac[i] - back[i]).abs() < 1e-10,
                "Roundtrip failed at {i}: {:.6} vs {:.6}", frac[i], back[i]);
        }
    }

    #[test]
    fn test_wrap_frac() {
        let wrapped = UnitCell::wrap_frac([1.3, -0.2, 2.7]);
        assert!((wrapped[0] - 0.3).abs() < 1e-10);
        assert!((wrapped[1] - 0.8).abs() < 1e-10);
        assert!((wrapped[2] - 0.7).abs() < 1e-10);
    }

    #[test]
    fn test_minimum_image_distance_cubic() {
        let cell = UnitCell::cubic(10.0);
        // Two points: one at (1,0,0) and one at (9,0,0)
        // Under PBC, minimum distance should be 2.0 Å (not 8.0)
        let dist = cell.minimum_image_distance([1.0, 0.0, 0.0], [9.0, 0.0, 0.0]);
        assert!((dist - 2.0).abs() < 1e-10,
            "Minimum image distance should be 2.0, got {dist:.6}");
    }

    #[test]
    fn test_supercell_offsets() {
        let cell = UnitCell::cubic(5.0);
        let (super_cell, offsets) = cell.supercell(2, 2, 2);
        assert_eq!(offsets.len(), 8);
        assert!((super_cell.lattice[0][0] - 10.0).abs() < 1e-10);
        assert!((super_cell.volume() - 5.0f64.powi(3) * 8.0).abs() < 1e-6);
    }

    #[test]
    fn test_from_parameters_orthorhombic() {
        let p = CellParameters { a: 5.0, b: 7.0, c: 9.0, alpha: 90.0, beta: 90.0, gamma: 90.0 };
        let cell = UnitCell::from_parameters(&p);
        assert!((cell.volume() - 315.0).abs() < 1e-6);
        let back = cell.parameters();
        assert!((back.a - 5.0).abs() < 1e-6);
        assert!((back.b - 7.0).abs() < 1e-6);
        assert!((back.c - 9.0).abs() < 1e-6);
    }
}
