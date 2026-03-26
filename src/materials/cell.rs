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
    ///
    /// Returns a degenerate cell if angles do not produce positive volume.
    pub fn from_parameters(params: &CellParameters) -> Self {
        let alpha = params.alpha.to_radians();
        let beta = params.beta.to_radians();
        let gamma = params.gamma.to_radians();

        let cos_alpha = alpha.cos();
        let cos_beta = beta.cos();
        let cos_gamma = gamma.cos();
        let sin_gamma = gamma.sin();

        // Validate: angles must produce positive volume.
        // The volume factor is: 1 - cos²α - cos²β - cos²γ + 2·cosα·cosβ·cosγ > 0
        let vol_factor = 1.0 - cos_alpha * cos_alpha - cos_beta * cos_beta - cos_gamma * cos_gamma
            + 2.0 * cos_alpha * cos_beta * cos_gamma;
        if vol_factor <= 0.0 {
            eprintln!(
                "Warning: cell angles α={:.1}° β={:.1}° γ={:.1}° produce non-positive volume",
                params.alpha, params.beta, params.gamma
            );
        }

        // γ must be in (0°, 180°) for a valid cell
        if sin_gamma.abs() < 1e-12 {
            return Self {
                lattice: [
                    [params.a, 0.0, 0.0],
                    [0.0, params.b, 0.0],
                    [0.0, 0.0, params.c],
                ],
            };
        }

        // Standard crystallographic convention: a along x, b in xy-plane
        let ax = params.a;
        let bx = params.b * cos_gamma;
        let by = params.b * sin_gamma;
        let cx = params.c * cos_beta;
        let cy = params.c * (cos_alpha - cos_beta * cos_gamma) / sin_gamma;
        let cz_sq = params.c * params.c - cx * cx - cy * cy;

        // Negative cz² means angles are incompatible — clamp to avoid NaN
        let cz = if cz_sq > 0.0 { cz_sq.sqrt() } else { 0.0 };

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
    pub(crate) fn inverse_matrix(&self) -> [[f64; 3]; 3] {
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

/// Result of powder XRD simulation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PowderXrdResult {
    /// 2θ angles (degrees) for each reflection.
    pub two_theta: Vec<f64>,
    /// Relative intensities (0–100).
    pub intensities: Vec<f64>,
    /// Miller indices (h, k, l) for each reflection.
    pub miller_indices: Vec<[i32; 3]>,
    /// d-spacings (Å) for each reflection.
    pub d_spacings: Vec<f64>,
}

/// Simulate powder X-ray diffraction pattern from unit cell and atoms.
///
/// Uses Bragg's law: $n\lambda = 2d\sin\theta$ with Cu-Kα radiation (λ=1.5406 Å).
///
/// `cell`: unit cell definition.
/// `elements`: atomic numbers of all atoms.
/// `frac_coords`: fractional coordinates of all atoms.
/// `two_theta_max`: maximum 2θ angle in degrees (default 90°).
pub fn simulate_powder_xrd(
    cell: &UnitCell,
    elements: &[u8],
    frac_coords: &[[f64; 3]],
    two_theta_max: f64,
) -> PowderXrdResult {
    let lambda = 1.5406; // Cu-Kα wavelength in Å

    // Max h, k, l to consider
    let params = cell.parameters();
    let d_min = lambda / (2.0 * (two_theta_max.to_radians() / 2.0).sin());
    let h_max = (params.a / d_min).ceil() as i32 + 1;
    let k_max = (params.b / d_min).ceil() as i32 + 1;
    let l_max = (params.c / d_min).ceil() as i32 + 1;

    let mut reflections = Vec::new();

    for h in -h_max..=h_max {
        for k in -k_max..=k_max {
            for l in -l_max..=l_max {
                if h == 0 && k == 0 && l == 0 {
                    continue;
                }

                // d-spacing from reciprocal lattice
                let g = [
                    h as f64 * cell.inverse_matrix()[0][0]
                        + k as f64 * cell.inverse_matrix()[1][0]
                        + l as f64 * cell.inverse_matrix()[2][0],
                    h as f64 * cell.inverse_matrix()[0][1]
                        + k as f64 * cell.inverse_matrix()[1][1]
                        + l as f64 * cell.inverse_matrix()[2][1],
                    h as f64 * cell.inverse_matrix()[0][2]
                        + k as f64 * cell.inverse_matrix()[1][2]
                        + l as f64 * cell.inverse_matrix()[2][2],
                ];
                let g_len = norm3(g);
                let d = 1.0 / g_len;

                // Bragg: 2d sin(θ) = λ → sin(θ) = λ/(2d)
                let sin_theta = lambda / (2.0 * d);
                if !(0.0..=1.0).contains(&sin_theta) {
                    continue;
                }
                let two_theta_val = 2.0 * sin_theta.asin().to_degrees();
                if !(1.0..=two_theta_max).contains(&two_theta_val) {
                    continue;
                }

                // Structure factor (simplified — uses atomic number as scattering factor)
                let mut f_real = 0.0f64;
                let mut f_imag = 0.0f64;
                for (i, &elem) in elements.iter().enumerate() {
                    let f_atom = elem as f64; // Approximation: f ≈ Z for low angles
                    let phase = 2.0
                        * std::f64::consts::PI
                        * (h as f64 * frac_coords[i][0]
                            + k as f64 * frac_coords[i][1]
                            + l as f64 * frac_coords[i][2]);
                    f_real += f_atom * phase.cos();
                    f_imag += f_atom * phase.sin();
                }
                let intensity = f_real * f_real + f_imag * f_imag;

                if intensity > 1e-6 {
                    reflections.push((two_theta_val, intensity, [h, k, l], d));
                }
            }
        }
    }

    // Merge equivalent reflections (within 0.01° of 2θ)
    reflections.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

    let mut merged_theta: Vec<f64> = Vec::new();
    let mut merged_intensity = Vec::new();
    let mut merged_hkl = Vec::new();
    let mut merged_d = Vec::new();

    for (tt, intens, hkl, d) in &reflections {
        if let Some(last_tt) = merged_theta.last() {
            if (*tt - *last_tt).abs() < 0.01 {
                let idx = merged_intensity.len() - 1;
                merged_intensity[idx] += intens;
                continue;
            }
        }
        merged_theta.push(*tt);
        merged_intensity.push(*intens);
        merged_hkl.push(*hkl);
        merged_d.push(*d);
    }

    // Normalize intensities to 0–100
    let max_i = merged_intensity
        .iter()
        .cloned()
        .fold(0.0f64, f64::max)
        .max(1e-10);
    for i in merged_intensity.iter_mut() {
        *i = *i / max_i * 100.0;
    }

    PowderXrdResult {
        two_theta: merged_theta,
        intensities: merged_intensity,
        miller_indices: merged_hkl,
        d_spacings: merged_d,
    }
}

/// A symmetry operation (rotation + translation in fractional coords).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SymmetryOperation {
    /// 3×3 rotation matrix (integer in fractional space).
    pub rotation: [[i32; 3]; 3],
    /// Translation vector in fractional coordinates.
    pub translation: [f64; 3],
    /// Human-readable label (e.g., "x,y,z" or "-x,-y,z").
    pub label: String,
}

/// Space group definition with symmetry operations.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SpaceGroup {
    /// International Tables number (1–230).
    pub number: u16,
    /// Hermann-Mauguin symbol.
    pub symbol: String,
    /// Crystal system.
    pub crystal_system: String,
    /// Symmetry operations (general positions).
    pub operations: Vec<SymmetryOperation>,
}

/// Get symmetry operations for common space groups.
pub fn get_space_group(number: u16) -> Option<SpaceGroup> {
    match number {
        1 => Some(SpaceGroup {
            number: 1,
            symbol: "P1".to_string(),
            crystal_system: "triclinic".to_string(),
            operations: vec![SymmetryOperation {
                rotation: [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                translation: [0.0, 0.0, 0.0],
                label: "x,y,z".to_string(),
            }],
        }),
        2 => Some(SpaceGroup {
            number: 2,
            symbol: "P-1".to_string(),
            crystal_system: "triclinic".to_string(),
            operations: vec![
                SymmetryOperation {
                    rotation: [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                    translation: [0.0, 0.0, 0.0],
                    label: "x,y,z".to_string(),
                },
                SymmetryOperation {
                    rotation: [[-1, 0, 0], [0, -1, 0], [0, 0, -1]],
                    translation: [0.0, 0.0, 0.0],
                    label: "-x,-y,-z".to_string(),
                },
            ],
        }),
        14 => Some(SpaceGroup {
            number: 14,
            symbol: "P2_1/c".to_string(),
            crystal_system: "monoclinic".to_string(),
            operations: vec![
                SymmetryOperation {
                    rotation: [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                    translation: [0.0, 0.0, 0.0],
                    label: "x,y,z".to_string(),
                },
                SymmetryOperation {
                    rotation: [[-1, 0, 0], [0, 1, 0], [0, 0, -1]],
                    translation: [0.0, 0.5, 0.5],
                    label: "-x,y+1/2,-z+1/2".to_string(),
                },
                SymmetryOperation {
                    rotation: [[-1, 0, 0], [0, -1, 0], [0, 0, -1]],
                    translation: [0.0, 0.0, 0.0],
                    label: "-x,-y,-z".to_string(),
                },
                SymmetryOperation {
                    rotation: [[1, 0, 0], [0, -1, 0], [0, 0, 1]],
                    translation: [0.0, 0.5, 0.5],
                    label: "x,-y+1/2,z+1/2".to_string(),
                },
            ],
        }),
        225 => Some(SpaceGroup {
            number: 225,
            symbol: "Fm-3m".to_string(),
            crystal_system: "cubic".to_string(),
            operations: vec![
                SymmetryOperation {
                    rotation: [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                    translation: [0.0, 0.0, 0.0],
                    label: "x,y,z".to_string(),
                },
                SymmetryOperation {
                    rotation: [[-1, 0, 0], [0, -1, 0], [0, 0, 1]],
                    translation: [0.0, 0.0, 0.0],
                    label: "-x,-y,z".to_string(),
                },
                SymmetryOperation {
                    rotation: [[-1, 0, 0], [0, 1, 0], [0, 0, -1]],
                    translation: [0.0, 0.0, 0.0],
                    label: "-x,y,-z".to_string(),
                },
                SymmetryOperation {
                    rotation: [[1, 0, 0], [0, -1, 0], [0, 0, -1]],
                    translation: [0.0, 0.0, 0.0],
                    label: "x,-y,-z".to_string(),
                },
            ],
        }),
        _ => None,
    }
}

/// Apply a symmetry operation to a fractional coordinate.
pub fn apply_symmetry(op: &SymmetryOperation, frac: [f64; 3]) -> [f64; 3] {
    let r = op.rotation;
    let t = op.translation;
    [
        r[0][0] as f64 * frac[0] + r[0][1] as f64 * frac[1] + r[0][2] as f64 * frac[2] + t[0],
        r[1][0] as f64 * frac[0] + r[1][1] as f64 * frac[1] + r[1][2] as f64 * frac[2] + t[1],
        r[2][0] as f64 * frac[0] + r[2][1] as f64 * frac[1] + r[2][2] as f64 * frac[2] + t[2],
    ]
}

/// Generate all symmetry-equivalent positions from a set of asymmetric unit atoms.
pub fn expand_by_symmetry(
    space_group: &SpaceGroup,
    frac_coords: &[[f64; 3]],
    elements: &[u8],
) -> (Vec<[f64; 3]>, Vec<u8>) {
    let mut all_coords = Vec::new();
    let mut all_elements = Vec::new();

    for (i, &fc) in frac_coords.iter().enumerate() {
        for op in &space_group.operations {
            let new_fc = apply_symmetry(op, fc);
            let wrapped = UnitCell::wrap_frac(new_fc);

            // Check for duplicate (within tolerance)
            let is_dup = all_coords.iter().any(|existing: &[f64; 3]| {
                let dx = (existing[0] - wrapped[0]).abs();
                let dy = (existing[1] - wrapped[1]).abs();
                let dz = (existing[2] - wrapped[2]).abs();
                // Handle PBC wrap
                let dx = dx.min(1.0 - dx);
                let dy = dy.min(1.0 - dy);
                let dz = dz.min(1.0 - dz);
                dx < 0.01 && dy < 0.01 && dz < 0.01
            });

            if !is_dup {
                all_coords.push(wrapped);
                all_elements.push(elements[i]);
            }
        }
    }

    (all_coords, all_elements)
}

/// Result of geometric pore analysis.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PorosityResult {
    /// Geometric pore volume (Å³).
    pub pore_volume: f64,
    /// Porosity fraction (void / total volume).
    pub porosity: f64,
    /// Largest cavity diameter (Å) — largest sphere that fits in a pore.
    pub largest_cavity_diameter: f64,
    /// Pore-limiting diameter (Å) — largest sphere that can pass through.
    pub pore_limiting_diameter: f64,
}

/// Compute geometric porosity using grid-based cavity detection.
///
/// Probes the unit cell with a grid and marks points as "void" if they are
/// farther than any atom's van der Waals radius + probe radius.
pub fn compute_porosity(
    cell: &UnitCell,
    elements: &[u8],
    frac_coords: &[[f64; 3]],
    probe_radius: f64,
    grid_spacing: f64,
) -> PorosityResult {
    let params = cell.parameters();
    let nx = (params.a / grid_spacing).ceil() as usize;
    let ny = (params.b / grid_spacing).ceil() as usize;
    let nz = (params.c / grid_spacing).ceil() as usize;

    let total_points = nx * ny * nz;
    let mut void_count = 0usize;
    let mut max_void_dist = 0.0f64;
    let mut min_void_dist = f64::INFINITY;

    // Precompute Cartesian atom positions
    let atom_positions: Vec<[f64; 3]> = frac_coords
        .iter()
        .map(|&fc| cell.frac_to_cart(fc))
        .collect();
    let vdw_radii: Vec<f64> = elements.iter().map(|&z| vdw_radius(z)).collect();

    for ix in 0..nx {
        for iy in 0..ny {
            for iz in 0..nz {
                let frac = [
                    (ix as f64 + 0.5) / nx as f64,
                    (iy as f64 + 0.5) / ny as f64,
                    (iz as f64 + 0.5) / nz as f64,
                ];
                let cart = cell.frac_to_cart(frac);

                // Find minimum distance to any atom surface
                let mut min_surface_dist = f64::INFINITY;
                for (j, &pos) in atom_positions.iter().enumerate() {
                    let d = cell.minimum_image_distance(cart, pos);
                    let surface_d = d - vdw_radii[j];
                    min_surface_dist = min_surface_dist.min(surface_d);
                }

                if min_surface_dist > probe_radius {
                    void_count += 1;
                    max_void_dist = max_void_dist.max(min_surface_dist);
                    min_void_dist = min_void_dist.min(min_surface_dist);
                }
            }
        }
    }

    let vol = cell.volume();
    let voxel_vol = vol / total_points as f64;
    let pore_volume = void_count as f64 * voxel_vol;
    let porosity = void_count as f64 / total_points as f64;

    PorosityResult {
        pore_volume,
        porosity,
        largest_cavity_diameter: if max_void_dist > 0.0 {
            2.0 * max_void_dist
        } else {
            0.0
        },
        pore_limiting_diameter: if min_void_dist < f64::INFINITY && void_count > 0 {
            2.0 * min_void_dist
        } else {
            0.0
        },
    }
}

/// Van der Waals radius for common elements (Å).
fn vdw_radius(z: u8) -> f64 {
    match z {
        1 => 1.20,
        6 => 1.70,
        7 => 1.55,
        8 => 1.52,
        9 => 1.47,
        15 => 1.80,
        16 => 1.80,
        17 => 1.75,
        22 => 2.00, // Ti
        26 => 2.00, // Fe
        29 => 1.40, // Cu
        30 => 1.39, // Zn
        35 => 1.85,
        53 => 1.98,
        _ => 1.80,
    }
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
            a: 10.0,
            b: 12.0,
            c: 8.0,
            alpha: 90.0,
            beta: 90.0,
            gamma: 120.0,
        });
        let frac = [0.3, 0.4, 0.7];
        let cart = cell.frac_to_cart(frac);
        let back = cell.cart_to_frac(cart);
        for i in 0..3 {
            assert!(
                (frac[i] - back[i]).abs() < 1e-10,
                "Roundtrip failed at {i}: {:.6} vs {:.6}",
                frac[i],
                back[i]
            );
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
        assert!(
            (dist - 2.0).abs() < 1e-10,
            "Minimum image distance should be 2.0, got {dist:.6}"
        );
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
        let p = CellParameters {
            a: 5.0,
            b: 7.0,
            c: 9.0,
            alpha: 90.0,
            beta: 90.0,
            gamma: 90.0,
        };
        let cell = UnitCell::from_parameters(&p);
        assert!((cell.volume() - 315.0).abs() < 1e-6);
        let back = cell.parameters();
        assert!((back.a - 5.0).abs() < 1e-6);
        assert!((back.b - 7.0).abs() < 1e-6);
        assert!((back.c - 9.0).abs() < 1e-6);
    }
}
