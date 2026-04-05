use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde::{Deserialize, Serialize};

const AMU_ANGFS2_TO_KCAL_MOL: f64 = 2_390.057_361_533_49;
const R_GAS_KCAL_MOLK: f64 = 0.001_987_204_258_640_83;
const EV_TO_KCAL_MOL: f64 = 23.060_5;
const HARTREE_TO_KCAL_MOL: f64 = 627.509_474_063_1;

/// Force backend for molecular dynamics.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum MdBackend {
    /// UFF force field (fast, approximate).
    Uff,
    /// PM3 semi-empirical (slower, more accurate).
    Pm3,
    /// GFN-xTB tight-binding (moderate speed, good for metals).
    Xtb,
}

/// Energy backend for NEB path calculations.
///
/// Supports all methods that can provide energy and gradients (analytical
/// or numerical) for NEB image relaxation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum NebBackend {
    /// UFF force field — fast exploratory paths.
    Uff,
    /// MMFF94 force field — better organic coverage than UFF.
    Mmff94,
    /// PM3 semi-empirical — analytical gradients.
    Pm3,
    /// GFN0-xTB tight-binding — analytical gradients.
    Xtb,
    /// GFN1-xTB — numerical gradients (energy-only + finite differences).
    Gfn1,
    /// GFN2-xTB — numerical gradients (energy-only + finite differences).
    Gfn2,
    /// HF-3c minimal-basis Hartree-Fock — numerical gradients (expensive).
    Hf3c,
}

impl NebBackend {
    /// Parse a method string into a `NebBackend`.
    pub fn from_method(s: &str) -> Result<Self, String> {
        match s.trim().to_ascii_lowercase().as_str() {
            "uff" => Ok(Self::Uff),
            "mmff94" | "mmff" => Ok(Self::Mmff94),
            "pm3" => Ok(Self::Pm3),
            "xtb" | "gfn0" | "gfn0-xtb" | "gfn0_xtb" => Ok(Self::Xtb),
            "gfn1" | "gfn1-xtb" | "gfn1_xtb" => Ok(Self::Gfn1),
            "gfn2" | "gfn2-xtb" | "gfn2_xtb" => Ok(Self::Gfn2),
            "hf3c" | "hf-3c" => Ok(Self::Hf3c),
            other => Err(format!(
                "Unknown NEB backend '{}'. Expected: uff, mmff94, pm3, xtb, gfn1, gfn2, hf3c",
                other
            )),
        }
    }

    /// method label for human display.
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Uff => "uff",
            Self::Mmff94 => "mmff94",
            Self::Pm3 => "pm3",
            Self::Xtb => "xtb",
            Self::Gfn1 => "gfn1",
            Self::Gfn2 => "gfn2",
            Self::Hf3c => "hf3c",
        }
    }
}

/// One trajectory frame for molecular-dynamics sampling.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MdFrame {
    /// Integration step index.
    pub step: usize,
    /// Elapsed simulation time in femtoseconds.
    pub time_fs: f64,
    /// Flat xyz coordinates in angstroms.
    pub coords: Vec<f64>,
    /// Potential energy from UFF in kcal/mol.
    pub potential_energy_kcal_mol: f64,
    /// Kinetic energy in kcal/mol.
    pub kinetic_energy_kcal_mol: f64,
    /// Instantaneous temperature estimate in K.
    pub temperature_k: f64,
}

/// Full trajectory output for exploratory molecular dynamics.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MdTrajectory {
    /// Stored MD frames.
    pub frames: Vec<MdFrame>,
    /// Timestep in femtoseconds.
    pub dt_fs: f64,
    /// Notes and caveats for interpretation.
    pub notes: Vec<String>,
    /// Energy conservation drift: (E_total_final - E_total_initial) / |E_total_initial| * 100%.
    /// Only meaningful for NVE (no thermostat) simulations.
    pub energy_drift_percent: Option<f64>,
}

/// One image (node) on a simplified NEB path.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NebImage {
    /// Image index from reactant (0) to product (n-1).
    pub index: usize,
    /// Flat xyz coordinates in angstroms.
    pub coords: Vec<f64>,
    /// UFF potential energy in kcal/mol.
    pub potential_energy_kcal_mol: f64,
}

/// Simplified NEB output for low-cost pathway exploration.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NebPathResult {
    /// Ordered path images from reactant to product.
    pub images: Vec<NebImage>,
    /// Notes and caveats.
    pub notes: Vec<String>,
}

pub fn atomic_mass_amu(z: u8) -> f64 {
    match z {
        1 => 1.008,
        5 => 10.81,
        6 => 12.011,
        7 => 14.007,
        8 => 15.999,
        9 => 18.998,
        14 => 28.085,
        15 => 30.974,
        16 => 32.06,
        17 => 35.45,
        35 => 79.904,
        53 => 126.904,
        26 => 55.845,
        46 => 106.42,
        78 => 195.084,
        _ => 12.0,
    }
}

fn sample_standard_normal(rng: &mut StdRng) -> f64 {
    let u1 = (1.0 - rng.gen::<f64>()).max(1e-12);
    let u2 = rng.gen::<f64>();
    (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
}

fn kinetic_energy_and_temperature(
    velocities: &[f64],
    masses_amu: &[f64],
    n_atoms: usize,
) -> (f64, f64) {
    let mut ke = 0.0;
    for i in 0..n_atoms {
        let vx = velocities[3 * i];
        let vy = velocities[3 * i + 1];
        let vz = velocities[3 * i + 2];
        let v2 = vx * vx + vy * vy + vz * vz;
        ke += 0.5 * masses_amu[i] * v2 * AMU_ANGFS2_TO_KCAL_MOL;
    }
    let dof = (3 * n_atoms).saturating_sub(6).max(1) as f64;
    let t = 2.0 * ke / (dof * R_GAS_KCAL_MOLK);
    (ke, t)
}

/// Run short exploratory molecular dynamics using Velocity Verlet and optional Berendsen NVT.
pub fn simulate_velocity_verlet_uff(
    smiles: &str,
    coords: &[f64],
    n_steps: usize,
    dt_fs: f64,
    seed: u64,
    target_temp_and_tau: Option<(f64, f64)>,
) -> Result<MdTrajectory, String> {
    if n_steps == 0 {
        return Err("n_steps must be > 0".to_string());
    }
    if dt_fs <= 0.0 {
        return Err("dt_fs must be > 0".to_string());
    }

    let mol = crate::graph::Molecule::from_smiles(smiles)?;
    let n_atoms = mol.graph.node_count();
    if coords.len() != n_atoms * 3 {
        return Err(format!(
            "coords length {} != 3 * atoms {}",
            coords.len(),
            n_atoms
        ));
    }

    let masses_amu: Vec<f64> = (0..n_atoms)
        .map(petgraph::graph::NodeIndex::new)
        .map(|idx| atomic_mass_amu(mol.graph[idx].element))
        .collect();

    let ff = crate::forcefield::builder::build_uff_force_field(&mol);
    let mut x = coords.to_vec();
    let mut grad = vec![0.0f64; n_atoms * 3];
    let mut potential = ff.compute_system_energy_and_gradients(&x, &mut grad);

    let mut rng = StdRng::seed_from_u64(seed);
    let mut v = vec![0.0f64; n_atoms * 3];

    if let Some((target_temp_k, _tau_fs)) = target_temp_and_tau {
        for i in 0..n_atoms {
            let sigma = ((R_GAS_KCAL_MOLK * target_temp_k)
                / (masses_amu[i] * AMU_ANGFS2_TO_KCAL_MOL))
                .sqrt();
            v[3 * i] = sigma * sample_standard_normal(&mut rng);
            v[3 * i + 1] = sigma * sample_standard_normal(&mut rng);
            v[3 * i + 2] = sigma * sample_standard_normal(&mut rng);
        }
    }

    let (ke0, t0) = kinetic_energy_and_temperature(&v, &masses_amu, n_atoms);
    let mut frames = vec![MdFrame {
        step: 0,
        time_fs: 0.0,
        coords: x.clone(),
        potential_energy_kcal_mol: potential,
        kinetic_energy_kcal_mol: ke0,
        temperature_k: t0,
    }];

    for step in 1..=n_steps {
        let mut a = vec![0.0f64; n_atoms * 3];
        for i in 0..n_atoms {
            let m = masses_amu[i];
            a[3 * i] = -grad[3 * i] / (m * AMU_ANGFS2_TO_KCAL_MOL);
            a[3 * i + 1] = -grad[3 * i + 1] / (m * AMU_ANGFS2_TO_KCAL_MOL);
            a[3 * i + 2] = -grad[3 * i + 2] / (m * AMU_ANGFS2_TO_KCAL_MOL);
        }

        for i in 0..(n_atoms * 3) {
            x[i] += v[i] * dt_fs + 0.5 * a[i] * dt_fs * dt_fs;
        }

        potential = ff.compute_system_energy_and_gradients(&x, &mut grad);

        let mut a_new = vec![0.0f64; n_atoms * 3];
        for i in 0..n_atoms {
            let m = masses_amu[i];
            a_new[3 * i] = -grad[3 * i] / (m * AMU_ANGFS2_TO_KCAL_MOL);
            a_new[3 * i + 1] = -grad[3 * i + 1] / (m * AMU_ANGFS2_TO_KCAL_MOL);
            a_new[3 * i + 2] = -grad[3 * i + 2] / (m * AMU_ANGFS2_TO_KCAL_MOL);
        }

        for i in 0..(n_atoms * 3) {
            v[i] += 0.5 * (a[i] + a_new[i]) * dt_fs;
        }

        let (mut ke, mut temp_k) = kinetic_energy_and_temperature(&v, &masses_amu, n_atoms);
        if let Some((target_temp_k, tau_fs)) = target_temp_and_tau {
            let tau = tau_fs.max(1e-6);
            let lambda = (1.0 + (dt_fs / tau) * (target_temp_k / temp_k.max(1e-6) - 1.0)).sqrt();
            let lambda = lambda.clamp(0.5, 2.0);
            for vi in &mut v {
                *vi *= lambda;
            }
            let kt = kinetic_energy_and_temperature(&v, &masses_amu, n_atoms);
            ke = kt.0;
            temp_k = kt.1;
        }

        if !x.iter().all(|v| v.is_finite()) || !potential.is_finite() || !ke.is_finite() {
            return Err(format!(
                "MD diverged at step {} (non-finite coordinates/energy)",
                step
            ));
        }

        frames.push(MdFrame {
            step,
            time_fs: step as f64 * dt_fs,
            coords: x.clone(),
            potential_energy_kcal_mol: potential,
            kinetic_energy_kcal_mol: ke,
            temperature_k: temp_k,
        });
    }

    let mut notes = vec![
        "Velocity Verlet integration over UFF force-field gradients for short exploratory trajectories."
            .to_string(),
    ];
    if target_temp_and_tau.is_some() {
        notes.push(
            "Berendsen thermostat rescaling applied for approximate constant-temperature sampling."
                .to_string(),
        );
    } else {
        notes.push(
            "No thermostat applied (NVE-like propagation with current numerical approximations)."
                .to_string(),
        );
    }

    Ok(MdTrajectory {
        energy_drift_percent: if target_temp_and_tau.is_none() {
            let e0 = frames[0].potential_energy_kcal_mol + frames[0].kinetic_energy_kcal_mol;
            let ef = frames
                .last()
                .map(|f| f.potential_energy_kcal_mol + f.kinetic_energy_kcal_mol)
                .unwrap_or(e0);
            if e0.abs() > 1e-10 {
                Some((ef - e0) / e0.abs() * 100.0)
            } else {
                None
            }
        } else {
            None
        },
        frames,
        dt_fs,
        notes,
    })
}

/// Build a simplified NEB-like path using linear interpolation and spring-coupled relaxation.
pub fn compute_simplified_neb_path(
    smiles: &str,
    start_coords: &[f64],
    end_coords: &[f64],
    n_images: usize,
    n_iter: usize,
    spring_k: f64,
    step_size: f64,
) -> Result<NebPathResult, String> {
    if n_images < 2 {
        return Err("n_images must be >= 2".to_string());
    }
    if n_iter == 0 {
        return Err("n_iter must be > 0".to_string());
    }
    if step_size <= 0.0 {
        return Err("step_size must be > 0".to_string());
    }

    let mol = crate::graph::Molecule::from_smiles(smiles)?;
    let n_atoms = mol.graph.node_count();
    let n_xyz = n_atoms * 3;
    if start_coords.len() != n_xyz || end_coords.len() != n_xyz {
        return Err(format!(
            "start/end coords must each have length {} (3 * n_atoms)",
            n_xyz
        ));
    }

    let ff = crate::forcefield::builder::build_uff_force_field(&mol);

    let mut images = vec![vec![0.0f64; n_xyz]; n_images];
    for (img_idx, img) in images.iter_mut().enumerate() {
        let t = img_idx as f64 / (n_images - 1) as f64;
        for k in 0..n_xyz {
            img[k] = (1.0 - t) * start_coords[k] + t * end_coords[k];
        }
    }

    for _ in 0..n_iter {
        let prev = images.clone();

        #[cfg(feature = "parallel")]
        {
            use rayon::prelude::*;
            let updated: Vec<(usize, Vec<f64>)> = (1..(n_images - 1))
                .into_par_iter()
                .map(|i| {
                    let ff_local = crate::forcefield::builder::build_uff_force_field(&mol);
                    let mut grad = vec![0.0f64; n_xyz];
                    let _ = ff_local.compute_system_energy_and_gradients(&prev[i], &mut grad);
                    let mut new_img = prev[i].clone();
                    for k in 0..n_xyz {
                        let spring_force =
                            spring_k * (prev[i + 1][k] - 2.0 * prev[i][k] + prev[i - 1][k]);
                        let total_force = -grad[k] + spring_force;
                        new_img[k] = prev[i][k] + step_size * total_force;
                    }
                    (i, new_img)
                })
                .collect();
            for (i, img) in updated {
                images[i] = img;
            }
        }

        #[cfg(not(feature = "parallel"))]
        {
            for i in 1..(n_images - 1) {
                let mut grad = vec![0.0f64; n_xyz];
                let _ = ff.compute_system_energy_and_gradients(&prev[i], &mut grad);
                for k in 0..n_xyz {
                    let spring_force =
                        spring_k * (prev[i + 1][k] - 2.0 * prev[i][k] + prev[i - 1][k]);
                    let total_force = -grad[k] + spring_force;
                    images[i][k] = prev[i][k] + step_size * total_force;
                }
            }
        }
    }

    let mut out_images = Vec::with_capacity(n_images);
    for (i, coords) in images.into_iter().enumerate() {
        let mut grad = vec![0.0f64; n_xyz];
        let e = ff.compute_system_energy_and_gradients(&coords, &mut grad);
        out_images.push(NebImage {
            index: i,
            coords,
            potential_energy_kcal_mol: e,
        });
    }

    Ok(NebPathResult {
        images: out_images,
        notes: vec![
            "Simplified NEB: linear interpolation + spring-coupled UFF gradient relaxation on internal images."
                .to_string(),
            "This is a low-cost exploratory path tool and not a full climbing-image / tangent-projected NEB implementation."
                .to_string(),
        ],
    })
}

// ─── Configurable NEB backend dispatch ──────────────────────────────────────

/// Compute energy (kcal/mol) for a single point using a backend that only
/// exposes energy (no analytical gradients). Used by the numerical gradient
/// fallback for GFN1, GFN2, and HF-3c.
fn neb_point_energy_kcal(
    backend: NebBackend,
    elements: &[u8],
    coords: &[f64],
) -> Result<f64, String> {
    let positions: Vec<[f64; 3]> = coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();
    match backend {
        NebBackend::Gfn1 => {
            let r = crate::xtb::gfn1::solve_gfn1(elements, &positions)?;
            Ok(r.total_energy * EV_TO_KCAL_MOL)
        }
        NebBackend::Gfn2 => {
            let r = crate::xtb::gfn2::solve_gfn2(elements, &positions)?;
            Ok(r.total_energy * EV_TO_KCAL_MOL)
        }
        NebBackend::Hf3c => {
            let r = crate::hf::solve_hf3c(elements, &positions, &crate::hf::HfConfig::default())?;
            Ok(r.energy * HARTREE_TO_KCAL_MOL)
        }
        _ => unreachable!("neb_point_energy_kcal only for energy-only backends"),
    }
}

/// Compute energy and numerical gradient (central finite differences) for
/// backends without analytical gradients.
fn neb_numerical_gradient(
    backend: NebBackend,
    elements: &[u8],
    coords: &[f64],
    grad: &mut [f64],
) -> Result<f64, String> {
    let delta = 1e-5; // Å
    let e0 = neb_point_energy_kcal(backend, elements, coords)?;
    let mut displaced = coords.to_vec();
    for i in 0..coords.len() {
        displaced[i] = coords[i] + delta;
        let e_plus = neb_point_energy_kcal(backend, elements, &displaced)?;
        displaced[i] = coords[i] - delta;
        let e_minus = neb_point_energy_kcal(backend, elements, &displaced)?;
        displaced[i] = coords[i]; // restore
        grad[i] = (e_plus - e_minus) / (2.0 * delta);
    }
    Ok(e0)
}

/// Compute energy (kcal/mol) and gradients (kcal/mol/Å) for a NEB image
/// using the specified backend.
pub fn neb_energy_and_gradient(
    backend: NebBackend,
    _smiles: &str,
    elements: &[u8],
    mol: &crate::graph::Molecule,
    coords: &[f64],
    grad: &mut [f64],
) -> Result<f64, String> {
    match backend {
        NebBackend::Uff => {
            let ff = crate::forcefield::builder::build_uff_force_field(mol);
            Ok(ff.compute_system_energy_and_gradients(coords, grad))
        }
        NebBackend::Mmff94 => {
            let bonds: Vec<(usize, usize, u8)> = mol
                .graph
                .edge_indices()
                .map(|e| {
                    let (a, b) = mol.graph.edge_endpoints(e).unwrap();
                    let order = match mol.graph[e].order {
                        crate::graph::BondOrder::Single => 1u8,
                        crate::graph::BondOrder::Double => 2,
                        crate::graph::BondOrder::Triple => 3,
                        crate::graph::BondOrder::Aromatic => 2,
                        crate::graph::BondOrder::Unknown => 1,
                    };
                    (a.index(), b.index(), order)
                })
                .collect();
            let terms = crate::forcefield::mmff94::Mmff94Builder::build(elements, &bonds);
            let (energy, g) =
                crate::forcefield::mmff94::Mmff94Builder::total_energy(&terms, coords);
            grad[..g.len()].copy_from_slice(&g);
            Ok(energy)
        }
        NebBackend::Pm3 => {
            let positions: Vec<[f64; 3]> = coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();
            let r = crate::pm3::gradients::compute_pm3_gradient(elements, &positions)?;
            let energy_kcal = r.energy * EV_TO_KCAL_MOL;
            for (a, g) in r.gradients.iter().enumerate() {
                for d in 0..3 {
                    grad[a * 3 + d] = g[d] * EV_TO_KCAL_MOL;
                }
            }
            Ok(energy_kcal)
        }
        NebBackend::Xtb => {
            let positions: Vec<[f64; 3]> = coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();
            let r = crate::xtb::gradients::compute_xtb_gradient(elements, &positions)?;
            let energy_kcal = r.energy * EV_TO_KCAL_MOL;
            for (a, g) in r.gradients.iter().enumerate() {
                for d in 0..3 {
                    grad[a * 3 + d] = g[d] * EV_TO_KCAL_MOL;
                }
            }
            Ok(energy_kcal)
        }
        NebBackend::Gfn1 | NebBackend::Gfn2 | NebBackend::Hf3c => {
            neb_numerical_gradient(backend, elements, coords, grad)
        }
    }
}

/// Build a simplified NEB path with a configurable energy backend.
///
/// This is the multi-method version of [`compute_simplified_neb_path`].
/// Supply `method` as one of: `"uff"`, `"mmff94"`, `"pm3"`, `"xtb"`,
/// `"gfn1"`, `"gfn2"`, `"hf3c"`.
pub fn compute_simplified_neb_path_configurable(
    smiles: &str,
    start_coords: &[f64],
    end_coords: &[f64],
    n_images: usize,
    n_iter: usize,
    spring_k: f64,
    step_size: f64,
    method: &str,
) -> Result<NebPathResult, String> {
    let backend = NebBackend::from_method(method)?;
    if n_images < 2 {
        return Err("n_images must be >= 2".to_string());
    }
    if n_iter == 0 {
        return Err("n_iter must be > 0".to_string());
    }
    if step_size <= 0.0 {
        return Err("step_size must be > 0".to_string());
    }

    let mol = crate::graph::Molecule::from_smiles(smiles)?;
    let n_atoms = mol.graph.node_count();
    let n_xyz = n_atoms * 3;
    if start_coords.len() != n_xyz || end_coords.len() != n_xyz {
        return Err(format!(
            "start/end coords must each have length {} (3 * n_atoms)",
            n_xyz
        ));
    }

    let elements: Vec<u8> = (0..n_atoms)
        .map(|i| mol.graph[petgraph::graph::NodeIndex::new(i)].element)
        .collect();

    // Linear interpolation
    let mut images = vec![vec![0.0f64; n_xyz]; n_images];
    for (img_idx, img) in images.iter_mut().enumerate() {
        let t = img_idx as f64 / (n_images - 1) as f64;
        for k in 0..n_xyz {
            img[k] = (1.0 - t) * start_coords[k] + t * end_coords[k];
        }
    }

    // Spring-coupled NEB relaxation
    for _ in 0..n_iter {
        let prev = images.clone();
        for i in 1..(n_images - 1) {
            let mut grad = vec![0.0f64; n_xyz];
            let _ = neb_energy_and_gradient(backend, smiles, &elements, &mol, &prev[i], &mut grad)?;
            for k in 0..n_xyz {
                let spring_force = spring_k * (prev[i + 1][k] - 2.0 * prev[i][k] + prev[i - 1][k]);
                let total_force = -grad[k] + spring_force;
                images[i][k] = prev[i][k] + step_size * total_force;
            }
        }
    }

    // Evaluate final energies
    let mut out_images = Vec::with_capacity(n_images);
    for (i, coords) in images.into_iter().enumerate() {
        let mut grad = vec![0.0f64; n_xyz];
        let e = neb_energy_and_gradient(backend, smiles, &elements, &mol, &coords, &mut grad)?;
        out_images.push(NebImage {
            index: i,
            coords,
            potential_energy_kcal_mol: e,
        });
    }

    Ok(NebPathResult {
        images: out_images,
        notes: vec![
            format!(
                "Simplified NEB ({}) with spring-coupled gradient relaxation on {} internal images.",
                backend.as_str(),
                n_images.saturating_sub(2),
            ),
            "Low-cost exploratory path; not a full climbing-image / tangent-projected NEB."
                .to_string(),
        ],
    })
}

/// Compute energy-only for any NEB backend (used for single-point comparisons).
///
/// Returns energy in kcal/mol.
pub fn neb_backend_energy_kcal(method: &str, smiles: &str, coords: &[f64]) -> Result<f64, String> {
    let backend = NebBackend::from_method(method)?;
    let mol = crate::graph::Molecule::from_smiles(smiles)?;
    let n_atoms = mol.graph.node_count();
    let n_xyz = n_atoms * 3;
    if coords.len() != n_xyz {
        return Err(format!(
            "coords length {} != 3 * atoms {}",
            coords.len(),
            n_atoms
        ));
    }
    let elements: Vec<u8> = (0..n_atoms)
        .map(|i| mol.graph[petgraph::graph::NodeIndex::new(i)].element)
        .collect();
    let mut grad = vec![0.0f64; n_xyz];
    neb_energy_and_gradient(backend, smiles, &elements, &mol, coords, &mut grad)
}

/// Compute energy and return both energy (kcal/mol) and flat gradient (kcal/mol/Å).
pub fn neb_backend_energy_and_gradient(
    method: &str,
    smiles: &str,
    coords: &[f64],
) -> Result<(f64, Vec<f64>), String> {
    let backend = NebBackend::from_method(method)?;
    let mol = crate::graph::Molecule::from_smiles(smiles)?;
    let n_atoms = mol.graph.node_count();
    let n_xyz = n_atoms * 3;
    if coords.len() != n_xyz {
        return Err(format!(
            "coords length {} != 3 * atoms {}",
            coords.len(),
            n_atoms
        ));
    }
    let elements: Vec<u8> = (0..n_atoms)
        .map(|i| mol.graph[petgraph::graph::NodeIndex::new(i)].element)
        .collect();
    let mut grad = vec![0.0f64; n_xyz];
    let energy = neb_energy_and_gradient(backend, smiles, &elements, &mol, coords, &mut grad)?;
    Ok((energy, grad))
}

/// Compute energy and gradients using the specified backend.
///
/// Returns (energy_kcal_mol, gradients) or an error.
pub fn compute_backend_energy_and_gradients(
    backend: MdBackend,
    smiles: &str,
    elements: &[u8],
    coords: &[f64],
    grad: &mut [f64],
) -> Result<f64, String> {
    match backend {
        MdBackend::Uff => {
            let mol = crate::graph::Molecule::from_smiles(smiles)?;
            let ff = crate::forcefield::builder::build_uff_force_field(&mol);
            let e = ff.compute_system_energy_and_gradients(coords, grad);
            Ok(e)
        }
        MdBackend::Pm3 => {
            let positions: Vec<[f64; 3]> = coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();
            let grad_result = crate::pm3::gradients::compute_pm3_gradient(elements, &positions)?;
            // PM3 gradient returns eV/Å; convert to kcal/mol/Å
            let energy_kcal = grad_result.energy * 23.0605;
            for (a, g) in grad_result.gradients.iter().enumerate() {
                for d in 0..3 {
                    grad[a * 3 + d] = g[d] * 23.0605;
                }
            }
            Ok(energy_kcal)
        }
        MdBackend::Xtb => {
            let positions: Vec<[f64; 3]> = coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();
            let grad_result = crate::xtb::gradients::compute_xtb_gradient(elements, &positions)?;
            let energy_kcal = grad_result.energy * 23.0605;
            for (a, g) in grad_result.gradients.iter().enumerate() {
                for d in 0..3 {
                    grad[a * 3 + d] = g[d] * 23.0605;
                }
            }
            Ok(energy_kcal)
        }
    }
}

/// Run molecular dynamics using Velocity Verlet with a configurable force backend.
///
/// Supports UFF, PM3, and xTB backends. Note that PM3/xTB use numerical gradients
/// and are significantly slower than UFF.
pub fn simulate_velocity_verlet(
    smiles: &str,
    coords: &[f64],
    elements: &[u8],
    n_steps: usize,
    dt_fs: f64,
    seed: u64,
    target_temp_and_tau: Option<(f64, f64)>,
    backend: MdBackend,
) -> Result<MdTrajectory, String> {
    if n_steps == 0 {
        return Err("n_steps must be > 0".to_string());
    }
    let n_atoms = elements.len();
    if coords.len() != n_atoms * 3 {
        return Err(format!(
            "coords length {} != 3*atoms {}",
            coords.len(),
            n_atoms
        ));
    }

    let masses_amu: Vec<f64> = elements.iter().map(|&z| atomic_mass_amu(z)).collect();

    let mut x = coords.to_vec();
    let mut grad = vec![0.0f64; n_atoms * 3];
    let mut potential =
        compute_backend_energy_and_gradients(backend, smiles, elements, &x, &mut grad)?;

    let mut rng = StdRng::seed_from_u64(seed);
    let mut v = vec![0.0f64; n_atoms * 3];

    if let Some((target_temp_k, _)) = target_temp_and_tau {
        for i in 0..n_atoms {
            let sigma = ((R_GAS_KCAL_MOLK * target_temp_k)
                / (masses_amu[i] * AMU_ANGFS2_TO_KCAL_MOL))
                .sqrt();
            v[3 * i] = sigma * sample_standard_normal(&mut rng);
            v[3 * i + 1] = sigma * sample_standard_normal(&mut rng);
            v[3 * i + 2] = sigma * sample_standard_normal(&mut rng);
        }
    }

    let (ke0, t0) = kinetic_energy_and_temperature(&v, &masses_amu, n_atoms);
    let mut frames = vec![MdFrame {
        step: 0,
        time_fs: 0.0,
        coords: x.clone(),
        potential_energy_kcal_mol: potential,
        kinetic_energy_kcal_mol: ke0,
        temperature_k: t0,
    }];

    for step in 1..=n_steps {
        let mut a = vec![0.0f64; n_atoms * 3];
        for i in 0..n_atoms {
            let m = masses_amu[i];
            for k in 0..3 {
                a[3 * i + k] = -grad[3 * i + k] / (m * AMU_ANGFS2_TO_KCAL_MOL);
            }
        }

        for i in 0..(n_atoms * 3) {
            x[i] += v[i] * dt_fs + 0.5 * a[i] * dt_fs * dt_fs;
        }

        potential = compute_backend_energy_and_gradients(backend, smiles, elements, &x, &mut grad)?;

        let mut a_new = vec![0.0f64; n_atoms * 3];
        for i in 0..n_atoms {
            let m = masses_amu[i];
            for k in 0..3 {
                a_new[3 * i + k] = -grad[3 * i + k] / (m * AMU_ANGFS2_TO_KCAL_MOL);
            }
        }

        for i in 0..(n_atoms * 3) {
            v[i] += 0.5 * (a[i] + a_new[i]) * dt_fs;
        }

        let (mut ke, mut temp_k) = kinetic_energy_and_temperature(&v, &masses_amu, n_atoms);
        if let Some((target_temp_k, tau_fs)) = target_temp_and_tau {
            let tau = tau_fs.max(1e-6);
            let lambda = (1.0 + (dt_fs / tau) * (target_temp_k / temp_k.max(1e-6) - 1.0)).sqrt();
            let lambda = lambda.clamp(0.5, 2.0);
            for vi in &mut v {
                *vi *= lambda;
            }
            let kt = kinetic_energy_and_temperature(&v, &masses_amu, n_atoms);
            ke = kt.0;
            temp_k = kt.1;
        }

        if !x.iter().all(|v| v.is_finite()) || !potential.is_finite() {
            return Err(format!("MD diverged at step {}", step));
        }

        frames.push(MdFrame {
            step,
            time_fs: step as f64 * dt_fs,
            coords: x.clone(),
            potential_energy_kcal_mol: potential,
            kinetic_energy_kcal_mol: ke,
            temperature_k: temp_k,
        });
    }

    let energy_drift_percent = if target_temp_and_tau.is_none() {
        let e0 = frames[0].potential_energy_kcal_mol + frames[0].kinetic_energy_kcal_mol;
        let ef = frames
            .last()
            .map(|f| f.potential_energy_kcal_mol + f.kinetic_energy_kcal_mol)
            .unwrap_or(e0);
        if e0.abs() > 1e-10 {
            Some((ef - e0) / e0.abs() * 100.0)
        } else {
            None
        }
    } else {
        None
    };

    Ok(MdTrajectory {
        frames,
        dt_fs,
        notes: vec![format!("Velocity Verlet with {:?} backend.", backend)],
        energy_drift_percent,
    })
}

/// Nosé-Hoover chain thermostat for rigorous NVT sampling.
///
/// Implements a chain of `chain_length` thermostats coupled to velocities.
/// Produces canonical (NVT) ensemble with correct fluctuations.
pub fn simulate_nose_hoover(
    smiles: &str,
    coords: &[f64],
    elements: &[u8],
    n_steps: usize,
    dt_fs: f64,
    target_temp_k: f64,
    thermostat_mass: f64,
    chain_length: usize,
    seed: u64,
    backend: MdBackend,
) -> Result<MdTrajectory, String> {
    if n_steps == 0 {
        return Err("n_steps must be > 0".to_string());
    }
    let n_atoms = elements.len();
    if coords.len() != n_atoms * 3 {
        return Err("coords length mismatch".to_string());
    }

    let masses_amu: Vec<f64> = elements.iter().map(|&z| atomic_mass_amu(z)).collect();
    let dof = (3 * n_atoms).saturating_sub(6).max(1) as f64;
    let target_ke = 0.5 * dof * R_GAS_KCAL_MOLK * target_temp_k;

    let mut x = coords.to_vec();
    let mut grad = vec![0.0f64; n_atoms * 3];
    let mut potential =
        compute_backend_energy_and_gradients(backend, smiles, elements, &x, &mut grad)?;

    let mut rng = StdRng::seed_from_u64(seed);
    let mut v = vec![0.0f64; n_atoms * 3];
    for i in 0..n_atoms {
        let sigma =
            ((R_GAS_KCAL_MOLK * target_temp_k) / (masses_amu[i] * AMU_ANGFS2_TO_KCAL_MOL)).sqrt();
        v[3 * i] = sigma * sample_standard_normal(&mut rng);
        v[3 * i + 1] = sigma * sample_standard_normal(&mut rng);
        v[3 * i + 2] = sigma * sample_standard_normal(&mut rng);
    }

    // Nosé-Hoover chain variables
    let chain_len = chain_length.max(1);
    let q = vec![thermostat_mass; chain_len]; // thermostat masses
    let mut xi = vec![0.0f64; chain_len]; // thermostat positions
    let mut v_xi = vec![0.0f64; chain_len]; // thermostat velocities

    let (ke0, t0) = kinetic_energy_and_temperature(&v, &masses_amu, n_atoms);
    let mut frames = vec![MdFrame {
        step: 0,
        time_fs: 0.0,
        coords: x.clone(),
        potential_energy_kcal_mol: potential,
        kinetic_energy_kcal_mol: ke0,
        temperature_k: t0,
    }];

    let dt2 = dt_fs * 0.5;
    let dt4 = dt_fs * 0.25;

    for step in 1..=n_steps {
        let (ke, _) = kinetic_energy_and_temperature(&v, &masses_amu, n_atoms);

        // Nosé-Hoover chain: propagate thermostat (Yoshida-Suzuki-like)
        // Force on first thermostat
        let mut g_xi = vec![0.0f64; chain_len];
        g_xi[0] = (2.0 * ke - 2.0 * target_ke) / q[0];
        for j in 1..chain_len {
            g_xi[j] =
                (q[j - 1] * v_xi[j - 1] * v_xi[j - 1] - R_GAS_KCAL_MOLK * target_temp_k) / q[j];
        }

        // Update thermostat velocities (outside-in)
        if chain_len > 1 {
            v_xi[chain_len - 1] += g_xi[chain_len - 1] * dt4;
        }
        for j in (0..chain_len.saturating_sub(1)).rev() {
            let exp_factor = (-v_xi[j + 1] * dt4).exp();
            v_xi[j] = v_xi[j] * exp_factor + g_xi[j] * dt4;
        }

        // Scale velocities
        let scale = (-v_xi[0] * dt2).exp();
        for vi in v.iter_mut() {
            *vi *= scale;
        }

        // Update thermostat positions
        for j in 0..chain_len {
            xi[j] += v_xi[j] * dt2;
        }

        // Recompute KE after scaling
        let (ke_scaled, _) = kinetic_energy_and_temperature(&v, &masses_amu, n_atoms);
        g_xi[0] = (2.0 * ke_scaled - 2.0 * target_ke) / q[0];

        // Update thermostat velocities again (inside-out)
        for j in 0..chain_len.saturating_sub(1) {
            let exp_factor = (-v_xi[j + 1] * dt4).exp();
            v_xi[j] = v_xi[j] * exp_factor + g_xi[j] * dt4;
        }
        if chain_len > 1 {
            // Recompute force on last thermostat
            g_xi[chain_len - 1] = (q[chain_len - 2] * v_xi[chain_len - 2] * v_xi[chain_len - 2]
                - R_GAS_KCAL_MOLK * target_temp_k)
                / q[chain_len - 1];
            v_xi[chain_len - 1] += g_xi[chain_len - 1] * dt4;
        }

        // Velocity Verlet: half-step velocity
        for i in 0..n_atoms {
            let m = masses_amu[i];
            for k in 0..3 {
                v[3 * i + k] -= 0.5 * dt_fs * grad[3 * i + k] / (m * AMU_ANGFS2_TO_KCAL_MOL);
            }
        }

        // Position update
        for i in 0..(n_atoms * 3) {
            x[i] += v[i] * dt_fs;
        }

        // New forces
        potential = compute_backend_energy_and_gradients(backend, smiles, elements, &x, &mut grad)?;

        // Velocity Verlet: second half-step velocity
        for i in 0..n_atoms {
            let m = masses_amu[i];
            for k in 0..3 {
                v[3 * i + k] -= 0.5 * dt_fs * grad[3 * i + k] / (m * AMU_ANGFS2_TO_KCAL_MOL);
            }
        }

        let (ke_final, temp_k) = kinetic_energy_and_temperature(&v, &masses_amu, n_atoms);

        if !x.iter().all(|v| v.is_finite()) || !potential.is_finite() {
            return Err(format!("NH-MD diverged at step {}", step));
        }

        frames.push(MdFrame {
            step,
            time_fs: step as f64 * dt_fs,
            coords: x.clone(),
            potential_energy_kcal_mol: potential,
            kinetic_energy_kcal_mol: ke_final,
            temperature_k: temp_k,
        });
    }

    // Temperature analysis
    let temps: Vec<f64> = frames
        .iter()
        .skip(frames.len() / 5)
        .map(|f| f.temperature_k)
        .collect();
    let avg_temp = temps.iter().sum::<f64>() / temps.len() as f64;
    let temp_variance =
        temps.iter().map(|t| (t - avg_temp).powi(2)).sum::<f64>() / temps.len() as f64;

    Ok(MdTrajectory {
        frames,
        dt_fs,
        notes: vec![
            format!(
                "Nosé-Hoover chain thermostat (chain={}) with {:?} backend.",
                chain_len, backend
            ),
            format!(
                "Target T = {:.1} K, avg T = {:.1} K, σ(T) = {:.1} K",
                target_temp_k,
                avg_temp,
                temp_variance.sqrt()
            ),
        ],
        energy_drift_percent: None,
    })
}

/// Export an MD trajectory to XYZ format string.
///
/// Each frame becomes a block with atom count, comment line (step/energy/temp),
/// and atom coordinates.
pub fn trajectory_to_xyz(trajectory: &MdTrajectory, elements: &[u8]) -> String {
    let n_atoms = elements.len();
    let mut output = String::new();

    for frame in &trajectory.frames {
        output.push_str(&format!("{}\n", n_atoms));
        output.push_str(&format!(
            "step={} t={:.1}fs E_pot={:.4} E_kin={:.4} T={:.1}K\n",
            frame.step,
            frame.time_fs,
            frame.potential_energy_kcal_mol,
            frame.kinetic_energy_kcal_mol,
            frame.temperature_k,
        ));
        for i in 0..n_atoms {
            let sym = element_symbol(elements[i]);
            output.push_str(&format!(
                "{} {:.6} {:.6} {:.6}\n",
                sym,
                frame.coords[3 * i],
                frame.coords[3 * i + 1],
                frame.coords[3 * i + 2],
            ));
        }
    }
    output
}

fn element_symbol(z: u8) -> &'static str {
    match z {
        1 => "H",
        2 => "He",
        3 => "Li",
        4 => "Be",
        5 => "B",
        6 => "C",
        7 => "N",
        8 => "O",
        9 => "F",
        10 => "Ne",
        11 => "Na",
        12 => "Mg",
        13 => "Al",
        14 => "Si",
        15 => "P",
        16 => "S",
        17 => "Cl",
        18 => "Ar",
        19 => "K",
        20 => "Ca",
        22 => "Ti",
        24 => "Cr",
        25 => "Mn",
        26 => "Fe",
        27 => "Co",
        28 => "Ni",
        29 => "Cu",
        30 => "Zn",
        35 => "Br",
        44 => "Ru",
        46 => "Pd",
        47 => "Ag",
        53 => "I",
        78 => "Pt",
        79 => "Au",
        _ => "X",
    }
}

// ═══════════════════════════════════════════════════════════════════════════
//  Reaction dynamics: SMILES → embed → complex → NEB → full frame path
// ═══════════════════════════════════════════════════════════════════════════

/// Configuration for a full reaction dynamics computation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReactionDynamicsConfig {
    /// Number of NEB images for the reactive region.
    pub n_neb_images: usize,
    /// NEB spring-coupled relaxation iterations.
    pub neb_iterations: usize,
    /// NEB spring constant (kcal/mol/Å²).
    pub spring_k: f64,
    /// NEB optimisation step size (Å).
    pub step_size: f64,
    /// Number of approach frames (molecules approaching).
    pub n_approach_frames: usize,
    /// Number of departure frames (products separating).
    pub n_departure_frames: usize,
    /// Far distance (Å) at start/end of approach/departure.
    pub far_distance: f64,
    /// Target distance between reactive atoms in the complexes (Å).
    pub reactive_distance: f64,
    /// Random seed for conformer generation.
    pub seed: u64,
}

impl Default for ReactionDynamicsConfig {
    fn default() -> Self {
        Self {
            n_neb_images: 30,
            neb_iterations: 100,
            spring_k: 0.1,
            step_size: 0.01,
            n_approach_frames: 15,
            n_departure_frames: 15,
            far_distance: 8.0,
            reactive_distance: 2.0,
            seed: 42,
        }
    }
}

/// A single frame along the reaction coordinate.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReactionDynamicsFrame {
    /// Frame index (0-based).
    pub index: usize,
    /// Flat xyz coordinates `[x0,y0,z0, x1,y1,z1, ...]` in Å.
    pub coords: Vec<f64>,
    /// Energy at this frame (kcal/mol).
    pub energy_kcal_mol: f64,
    /// Phase label: `"approach"`, `"reaction"`, or `"departure"`.
    pub phase: String,
}

/// Full result of a reaction dynamics computation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReactionDynamicsResult {
    /// All frames ordered along the reaction coordinate.
    pub frames: Vec<ReactionDynamicsFrame>,
    /// Atomic numbers for every atom in each frame (same for all frames).
    pub elements: Vec<u8>,
    /// Frame index of the transition state (highest energy).
    pub ts_frame_index: usize,
    /// Activation energy (kcal/mol) = E(TS) − E(first frame).
    pub activation_energy_kcal_mol: f64,
    /// Reaction energy (kcal/mol) = E(last frame) − E(first frame).
    pub reaction_energy_kcal_mol: f64,
    /// Method used for energy/gradient evaluation.
    pub method: String,
    /// Number of atoms per frame.
    pub n_atoms: usize,
    /// Informational notes.
    pub notes: Vec<String>,
}

/// Centre-of-mass of a flat coordinate array.
fn com_flat(coords: &[f64]) -> [f64; 3] {
    let n = coords.len() / 3;
    if n == 0 {
        return [0.0; 3];
    }
    let mut cx = 0.0;
    let mut cy = 0.0;
    let mut cz = 0.0;
    for i in 0..n {
        cx += coords[i * 3];
        cy += coords[i * 3 + 1];
        cz += coords[i * 3 + 2];
    }
    let nf = n as f64;
    [cx / nf, cy / nf, cz / nf]
}

/// Translate coords so COM is at origin.
fn centre_at_origin(coords: &mut [f64]) {
    let [cx, cy, cz] = com_flat(coords);
    for i in (0..coords.len()).step_by(3) {
        coords[i] -= cx;
        coords[i + 1] -= cy;
        coords[i + 2] -= cz;
    }
}

/// Rodrigues rotation: rotate flat coords so direction `from` aligns with `to`.
fn rotate_to_align(coords: &mut [f64], from: [f64; 3], to: [f64; 3]) {
    let fl = (from[0] * from[0] + from[1] * from[1] + from[2] * from[2]).sqrt();
    let tl = (to[0] * to[0] + to[1] * to[1] + to[2] * to[2]).sqrt();
    if fl < 1e-10 || tl < 1e-10 {
        return;
    }
    let fx = from[0] / fl;
    let fy = from[1] / fl;
    let fz = from[2] / fl;
    let tx = to[0] / tl;
    let ty = to[1] / tl;
    let tz = to[2] / tl;

    let kx = fy * tz - fz * ty;
    let ky = fz * tx - fx * tz;
    let kz = fx * ty - fy * tx;
    let sin_a = (kx * kx + ky * ky + kz * kz).sqrt();
    let cos_a = fx * tx + fy * ty + fz * tz;

    if sin_a < 1e-10 {
        if cos_a > 0.0 {
            return; // already aligned
        }
        // Anti-aligned: reflect through a perpendicular plane
        for i in (0..coords.len()).step_by(3) {
            coords[i] = -coords[i];
        }
        return;
    }

    let nkx = kx / sin_a;
    let nky = ky / sin_a;
    let nkz = kz / sin_a;
    let c = cos_a;
    let s = sin_a;
    let t1 = 1.0 - c;

    for i in (0..coords.len()).step_by(3) {
        let x = coords[i];
        let y = coords[i + 1];
        let z = coords[i + 2];
        let dot = nkx * x + nky * y + nkz * z;
        coords[i] = x * c + (nky * z - nkz * y) * s + nkx * dot * t1;
        coords[i + 1] = y * c + (nkz * x - nkx * z) * s + nky * dot * t1;
        coords[i + 2] = z * c + (nkx * y - nky * x) * s + nkz * dot * t1;
    }
}

/// Build a reactive complex from separately-embedded molecule conformers.
///
/// Each molecule is centred at its own COM, then the first two are oriented so
/// their closest-pair atoms face each other at `reactive_dist` Å apart.
/// Returns `(flat_coords, elements)` with the overall COM at origin.
fn build_reaction_complex(
    conformers: &[crate::ConformerResult],
    reactive_dist: f64,
) -> (Vec<f64>, Vec<u8>) {
    if conformers.is_empty() {
        return (vec![], vec![]);
    }

    // Centre each molecule at its own COM
    let mut mols: Vec<(Vec<f64>, Vec<u8>)> = conformers
        .iter()
        .map(|c| {
            let mut coords = c.coords.clone();
            centre_at_origin(&mut coords);
            (coords, c.elements.clone())
        })
        .collect();

    if mols.len() == 1 {
        let (mut coords, elems) = mols.remove(0);
        centre_at_origin(&mut coords);
        return (coords, elems);
    }

    // Find closest inter-molecular pair between mol 0 and mol 1
    let n0 = mols[0].1.len();
    let n1 = mols[1].1.len();
    let mut best_i = 0usize;
    let mut best_j = 0usize;
    let mut best_d2 = f64::INFINITY;
    // Use a reference offset so we rank "face-to-face" proximity
    let ref_off = 4.0;
    for i in 0..n0 {
        let xi = mols[0].0[i * 3];
        let yi = mols[0].0[i * 3 + 1];
        let zi = mols[0].0[i * 3 + 2];
        for j in 0..n1 {
            let dx = mols[1].0[j * 3] + ref_off - xi;
            let dy = mols[1].0[j * 3 + 1] - yi;
            let dz = mols[1].0[j * 3 + 2] - zi;
            let d2 = dx * dx + dy * dy + dz * dz;
            if d2 < best_d2 {
                best_d2 = d2;
                best_i = i;
                best_j = j;
            }
        }
    }

    // Orient mol 0 reactive atom → +X
    {
        let rx = mols[0].0[best_i * 3];
        let ry = mols[0].0[best_i * 3 + 1];
        let rz = mols[0].0[best_i * 3 + 2];
        let rl = (rx * rx + ry * ry + rz * rz).sqrt();
        if rl > 0.05 {
            rotate_to_align(&mut mols[0].0, [rx / rl, ry / rl, rz / rl], [1.0, 0.0, 0.0]);
        }
    }

    // Orient mol 1 reactive atom → −X
    {
        let rx = mols[1].0[best_j * 3];
        let ry = mols[1].0[best_j * 3 + 1];
        let rz = mols[1].0[best_j * 3 + 2];
        let rl = (rx * rx + ry * ry + rz * rz).sqrt();
        if rl > 0.05 {
            rotate_to_align(
                &mut mols[1].0,
                [rx / rl, ry / rl, rz / rl],
                [-1.0, 0.0, 0.0],
            );
        }
    }

    // Place mol 1 so reactive atoms are `reactive_dist` apart along X
    let ra_x = mols[0].0[best_i * 3];
    let rb_x = mols[1].0[best_j * 3];
    let offset_x = ra_x - rb_x + reactive_dist;

    // Assemble combined coords
    let total_atoms: usize = mols.iter().map(|(_, e)| e.len()).sum();
    let mut all_coords = Vec::with_capacity(total_atoms * 3);
    let mut all_elements = Vec::with_capacity(total_atoms);

    // Mol 0 — no offset
    all_coords.extend_from_slice(&mols[0].0);
    all_elements.extend_from_slice(&mols[0].1);

    // Mol 1 — shifted along X
    for k in 0..n1 {
        all_coords.push(mols[1].0[k * 3] + offset_x);
        all_coords.push(mols[1].0[k * 3 + 1]);
        all_coords.push(mols[1].0[k * 3 + 2]);
    }
    all_elements.extend_from_slice(&mols[1].1);

    // Spectator molecules (m ≥ 2) — stacked further along +X
    let mut extra = offset_x + 4.0;
    for mol in mols.iter().skip(2) {
        for k in 0..mol.1.len() {
            all_coords.push(mol.0[k * 3] + extra);
            all_coords.push(mol.0[k * 3 + 1]);
            all_coords.push(mol.0[k * 3 + 2]);
        }
        all_elements.extend_from_slice(&mol.1);
        extra += 4.0;
    }

    centre_at_origin(&mut all_coords);
    (all_coords, all_elements)
}

/// Greedy atom mapping between two complexes by element and distance.
///
/// Returns `mapping[i]` = index in the product complex that corresponds
/// to reactant atom `i`. Matches atoms of the same element, choosing the
/// closest unmatched pair greedily.
fn map_atoms_greedy(
    r_elements: &[u8],
    r_coords: &[f64],
    p_elements: &[u8],
    p_coords: &[f64],
) -> Result<Vec<usize>, String> {
    let n = r_elements.len();
    if n != p_elements.len() {
        return Err(format!(
            "Atom count mismatch: {} reactant atoms vs {} product atoms",
            n,
            p_elements.len(),
        ));
    }

    // Group indices by element
    let mut r_by_elem: std::collections::HashMap<u8, Vec<usize>> = std::collections::HashMap::new();
    let mut p_by_elem: std::collections::HashMap<u8, Vec<usize>> = std::collections::HashMap::new();
    for i in 0..n {
        r_by_elem.entry(r_elements[i]).or_default().push(i);
        p_by_elem.entry(p_elements[i]).or_default().push(i);
    }

    let mut mapping = vec![0usize; n];
    let mut used_p = vec![false; n];

    for (&elem, r_indices) in &r_by_elem {
        let p_indices = p_by_elem
            .get(&elem)
            .ok_or_else(|| format!("Element Z={} present in reactants but not products", elem,))?;
        if r_indices.len() != p_indices.len() {
            return Err(format!(
                "Element Z={}: {} in reactants vs {} in products",
                elem,
                r_indices.len(),
                p_indices.len(),
            ));
        }

        // Build distance pairs and sort by distance
        let mut pairs: Vec<(usize, usize, f64)> = Vec::new();
        for &ri in r_indices {
            for &pi in p_indices {
                let dx = r_coords[ri * 3] - p_coords[pi * 3];
                let dy = r_coords[ri * 3 + 1] - p_coords[pi * 3 + 1];
                let dz = r_coords[ri * 3 + 2] - p_coords[pi * 3 + 2];
                pairs.push((ri, pi, dx * dx + dy * dy + dz * dz));
            }
        }
        pairs.sort_by(|a, b| a.2.partial_cmp(&b.2).unwrap_or(std::cmp::Ordering::Equal));

        let mut used_r = vec![false; n];
        for (ri, pi, _) in pairs {
            if used_r[ri] || used_p[pi] {
                continue;
            }
            mapping[ri] = pi;
            used_r[ri] = true;
            used_p[pi] = true;
        }
    }

    Ok(mapping)
}

/// Reorder product coords into reactant atom order using `mapping`.
fn reorder_coords(coords: &[f64], mapping: &[usize]) -> Vec<f64> {
    let mut out = vec![0.0; coords.len()];
    for (i, &j) in mapping.iter().enumerate() {
        out[i * 3] = coords[j * 3];
        out[i * 3 + 1] = coords[j * 3 + 1];
        out[i * 3 + 2] = coords[j * 3 + 2];
    }
    out
}

/// Rigid-body slide: move each molecule's atoms toward/away from the global COM.
///
/// `mol_ranges[m] = (start_atom, end_atom)` (exclusive end) per molecule.
/// `alpha` = 0 → atoms at `far_dist` from global COM; `alpha` = 1 → at complex position.
fn slide_molecules(
    complex_coords: &[f64],
    mol_ranges: &[(usize, usize)],
    alpha: f64,
    far_dist: f64,
) -> Vec<f64> {
    let mut out = complex_coords.to_vec();
    let [gcx, gcy, gcz] = com_flat(complex_coords);

    for &(start, end) in mol_ranges {
        let n = end - start;
        if n == 0 {
            continue;
        }
        let mut mx = 0.0;
        let mut my = 0.0;
        let mut mz = 0.0;
        for a in start..end {
            mx += complex_coords[a * 3];
            my += complex_coords[a * 3 + 1];
            mz += complex_coords[a * 3 + 2];
        }
        let nf = n as f64;
        mx /= nf;
        my /= nf;
        mz /= nf;

        let dx = mx - gcx;
        let dy = my - gcy;
        let dz = mz - gcz;
        let d = (dx * dx + dy * dy + dz * dz).sqrt();
        if d < 0.01 {
            continue;
        }
        let target = d + (far_dist - d) * (1.0 - alpha);
        let shift = target - d;
        let nx = dx / d;
        let ny = dy / d;
        let nz = dz / d;
        for a in start..end {
            out[a * 3] += nx * shift;
            out[a * 3 + 1] += ny * shift;
            out[a * 3 + 2] += nz * shift;
        }
    }
    out
}

/// Compute a full reaction dynamics path: embed reactants + products,
/// build oriented complexes, run NEB for the reactive region, and generate
/// approach/departure frames — all energies computed in Rust with the
/// chosen quantum-chemistry method.
///
/// # Arguments
///
/// * `reactant_smiles` — SMILES of each reactant fragment (e.g. `["CC(=O)O", "N"]`).
/// * `product_smiles`  — SMILES of each product fragment (e.g. `["CC(=O)N", "O"]`).
/// * `method`          — NEB backend: `"uff"`, `"mmff94"`, `"pm3"`, `"xtb"`, `"gfn1"`, `"gfn2"`, `"hf3c"`.
/// * `config`          — [`ReactionDynamicsConfig`] with NEB/path parameters.
///
/// Returns [`ReactionDynamicsResult`] with all frames, energies, and TS info.
pub fn compute_reaction_dynamics(
    reactant_smiles: &[&str],
    product_smiles: &[&str],
    method: &str,
    config: &ReactionDynamicsConfig,
) -> Result<ReactionDynamicsResult, String> {
    if reactant_smiles.is_empty() {
        return Err("At least one reactant SMILES is required".into());
    }
    if product_smiles.is_empty() {
        return Err("At least one product SMILES is required".into());
    }

    let backend = NebBackend::from_method(method)?;

    // ── 1. Embed all fragments ─────────────────────────────────────────
    let r_confs: Vec<crate::ConformerResult> = reactant_smiles
        .iter()
        .map(|s| crate::embed(s, config.seed))
        .collect();
    for (i, c) in r_confs.iter().enumerate() {
        if let Some(ref e) = c.error {
            return Err(format!(
                "Failed to embed reactant '{}': {}",
                reactant_smiles[i], e
            ));
        }
    }

    let p_confs: Vec<crate::ConformerResult> = product_smiles
        .iter()
        .map(|s| crate::embed(s, config.seed))
        .collect();
    for (i, c) in p_confs.iter().enumerate() {
        if let Some(ref e) = c.error {
            return Err(format!(
                "Failed to embed product '{}': {}",
                product_smiles[i], e
            ));
        }
    }

    // Validate atom conservation
    let r_total: usize = r_confs.iter().map(|c| c.num_atoms).sum();
    let p_total: usize = p_confs.iter().map(|c| c.num_atoms).sum();
    if r_total != p_total {
        return Err(format!(
            "Atom count mismatch: {} atoms in reactants vs {} in products — \
             atoms must be conserved",
            r_total, p_total,
        ));
    }

    // ── 2. Build oriented reactive complexes ───────────────────────────
    let (r_coords, r_elements) = build_reaction_complex(&r_confs, config.reactive_distance);
    let (p_coords, p_elements) = build_reaction_complex(&p_confs, config.reactive_distance);

    // ── 3. Atom mapping: greedy by element + distance ──────────────────
    let mapping = map_atoms_greedy(&r_elements, &r_coords, &p_elements, &p_coords)?;
    let p_reordered = reorder_coords(&p_coords, &mapping);

    // ── 4. Combined SMILES (dot-separated reactants) for NEB topology ──
    let combined_smiles = reactant_smiles.join(".");

    // ── 5. NEB between reactant and product complexes ──────────────────
    let neb = compute_simplified_neb_path_configurable(
        &combined_smiles,
        &r_coords,
        &p_reordered,
        config.n_neb_images,
        config.neb_iterations,
        config.spring_k,
        config.step_size,
        method,
    )?;

    // ── 6. Compute molecule ranges for approach/departure sliding ──────
    let r_mol_ranges: Vec<(usize, usize)> = {
        let mut ranges = Vec::new();
        let mut off = 0usize;
        for c in &r_confs {
            ranges.push((off, off + c.num_atoms));
            off += c.num_atoms;
        }
        ranges
    };

    // Product molecule ranges (after reordering into reactant atom order)
    let p_mol_assign: Vec<usize> = {
        let mut assign = vec![0usize; p_total];
        let mut off = 0usize;
        for (m, c) in p_confs.iter().enumerate() {
            for a in off..off + c.num_atoms {
                assign[a] = m;
            }
            off += c.num_atoms;
        }
        // Remap into reactant order via mapping
        mapping.iter().map(|&pi| assign[pi]).collect()
    };
    let p_mol_ranges: Vec<(usize, usize)> = {
        let n_mols = p_confs.len();
        let mut ranges = Vec::with_capacity(n_mols);
        for m in 0..n_mols {
            let atoms: Vec<usize> = p_mol_assign
                .iter()
                .enumerate()
                .filter(|(_, &a)| a == m)
                .map(|(i, _)| i)
                .collect();
            if let (Some(&lo), Some(&hi)) = (atoms.first(), atoms.last()) {
                ranges.push((lo, hi + 1));
            }
        }
        ranges
    };

    // ── 7. Build approach frames with energies ─────────────────────────
    let mut frames = Vec::new();
    let mol = crate::graph::Molecule::from_smiles(&combined_smiles)?;
    let n_xyz = r_total * 3;

    let na = config.n_approach_frames;
    for i in 0..na {
        let alpha = if na > 1 {
            i as f64 / (na - 1) as f64
        } else {
            1.0
        };
        let coords = slide_molecules(&r_coords, &r_mol_ranges, alpha, config.far_distance);
        let mut grad = vec![0.0; n_xyz];
        let energy = neb_energy_and_gradient(
            backend,
            &combined_smiles,
            &r_elements,
            &mol,
            &coords,
            &mut grad,
        )
        .unwrap_or(0.0);
        frames.push(ReactionDynamicsFrame {
            index: frames.len(),
            coords,
            energy_kcal_mol: energy,
            phase: "approach".into(),
        });
    }

    // ── 8. Reaction frames from NEB ────────────────────────────────────
    for img in &neb.images {
        frames.push(ReactionDynamicsFrame {
            index: frames.len(),
            coords: img.coords.clone(),
            energy_kcal_mol: img.potential_energy_kcal_mol,
            phase: "reaction".into(),
        });
    }

    // ── 9. Departure frames with energies ──────────────────────────────
    let nd = config.n_departure_frames;
    for i in 0..nd {
        let alpha = if nd > 1 {
            1.0 - i as f64 / (nd - 1) as f64
        } else {
            0.0
        };
        let coords = slide_molecules(&p_reordered, &p_mol_ranges, alpha, config.far_distance);
        let mut grad = vec![0.0; n_xyz];
        let energy = neb_energy_and_gradient(
            backend,
            &combined_smiles,
            &r_elements,
            &mol,
            &coords,
            &mut grad,
        )
        .unwrap_or(0.0);
        frames.push(ReactionDynamicsFrame {
            index: frames.len(),
            coords,
            energy_kcal_mol: energy,
            phase: "departure".into(),
        });
    }

    // ── 10. Find TS and compute energetics ─────────────────────────────
    let ts_idx = frames
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| {
            a.energy_kcal_mol
                .partial_cmp(&b.energy_kcal_mol)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .map(|(i, _)| i)
        .unwrap_or(0);

    let e_first = frames.first().map(|f| f.energy_kcal_mol).unwrap_or(0.0);
    let e_last = frames.last().map(|f| f.energy_kcal_mol).unwrap_or(0.0);
    let e_ts = frames[ts_idx].energy_kcal_mol;

    let mut notes = neb.notes;
    notes.push(format!(
        "Full reaction path: {} approach + {} NEB + {} departure = {} total frames.",
        na,
        neb.images.len(),
        nd,
        frames.len(),
    ));

    Ok(ReactionDynamicsResult {
        frames,
        elements: r_elements,
        ts_frame_index: ts_idx,
        activation_energy_kcal_mol: e_ts - e_first,
        reaction_energy_kcal_mol: e_last - e_first,
        method: method.to_string(),
        n_atoms: r_total,
        notes,
    })
}
