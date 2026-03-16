use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde::{Deserialize, Serialize};

const AMU_ANGFS2_TO_KCAL_MOL: f64 = 2_390.057_361_533_49;
const R_GAS_KCAL_MOLK: f64 = 0.001_987_204_258_640_83;

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

fn atomic_mass_amu(z: u8) -> f64 {
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
        for i in 1..(n_images - 1) {
            let mut grad = vec![0.0f64; n_xyz];
            let _ = ff.compute_system_energy_and_gradients(&prev[i], &mut grad);

            for k in 0..n_xyz {
                let spring_force = spring_k * (prev[i + 1][k] - 2.0 * prev[i][k] + prev[i - 1][k]);
                let total_force = -grad[k] + spring_force;
                images[i][k] = prev[i][k] + step_size * total_force;
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
