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
