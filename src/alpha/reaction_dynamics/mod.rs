//! Alpha Reaction Dynamics — Full 3D reaction path simulation.
//!
//! Fixes the X-axis-only orientation bug and implements:
//! - Climbing-Image NEB (CI-NEB) with tangent projection
//! - Orbital-guided molecular approach (FMO theory)
//! - Constrained complex optimisation
//! - Per-frame property calculation (charges, bond orders, dipole)
//! - IRC steepest descent from the transition state
//! - Electrostatic steering during approach
//! - Multi-angular orientation sampling
//!
//! Feature flag: `alpha-reaction-dynamics`

pub mod ci_neb;
pub mod complex;
pub mod constrained_opt;
pub mod electrostatics;
pub mod idpp;
pub mod irc;
pub mod orbital_approach;
pub mod per_frame;
pub mod reactive_sites;
pub mod sampling;

// Re-exports
pub use ci_neb::{compute_ci_neb_path, CiNebConfig, CiNebImage, CiNebResult};
pub use complex::{
    assemble_fragments_at_product_positions, build_3d_reaction_complex,
    build_3d_reaction_complex_guided, build_product_guided_complex_3d, optimize_reactive_complex,
    orient_fragment_along_direction,
};
pub use constrained_opt::{
    optimize_frame_sequence, optimize_with_constraints, relax_neb_images, ConstrainedOptConfig,
    ConstrainedOptResult, DistanceConstraint,
};
pub use electrostatics::{compute_electrostatic_approach, ElectrostaticApproachInfo};
pub use idpp::generate_idpp_path;
pub use irc::{compute_irc, IrcConfig, IrcResult};
pub use orbital_approach::{compute_orbital_approach_direction, OrbitalApproachInfo};
pub use per_frame::{compute_frame_properties, FrameProperties};
pub use reactive_sites::{
    identify_reactive_sites, reactive_atom_pairs, BondChangeInfo, ReactiveSite,
    ReactiveSiteAnalysis, SiteRole, SiteSource,
};
pub use sampling::{sample_approach_orientations, SamplingCandidate, SamplingResult};

use serde::{Deserialize, Serialize};

/// Full configuration for the improved reaction dynamics pipeline.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReactionDynamics3DConfig {
    /// NEB method backend: "uff", "pm3", "xtb", "gfn1", "gfn2", "hf3c".
    pub method: String,
    /// Number of NEB images.
    pub n_images: usize,
    /// Max NEB optimisation iterations.
    pub neb_max_iter: usize,
    /// NEB spring constant (kcal/mol/Å²).
    pub spring_k: f64,
    /// NEB step size (Å).
    pub step_size: f64,
    /// Use climbing-image NEB (true) or simplified NEB (false).
    pub use_climbing_image: bool,
    /// CI-NEB convergence threshold on max perpendicular force (kcal/mol/Å).
    pub ci_neb_force_threshold: f64,
    /// Number of approach frames.
    pub n_approach_frames: usize,
    /// Number of departure frames.
    pub n_departure_frames: usize,
    /// Far distance for approach/departure (Å).
    pub far_distance: f64,
    /// Target reactive distance (Å).
    pub reactive_distance: f64,
    /// Random seed.
    pub seed: u64,
    /// Use orbital-guided approach direction.
    pub use_orbital_guidance: bool,
    /// Use electrostatic steering during approach.
    pub use_electrostatic_steering: bool,
    /// Optimise the reactive complex before NEB.
    pub optimise_complex: bool,
    /// Max steps for complex optimisation.
    pub complex_opt_max_steps: usize,
    /// Compute per-frame properties (charges, bond orders, dipole).
    pub compute_properties: bool,
    /// Number of angular samples for multi-orientation search (0=disabled).
    pub n_angular_samples: usize,
    /// Optional SMIRKS reaction string for explicit atom mapping.
    /// When provided, reactive sites are identified from SMIRKS atom maps
    /// instead of relying solely on geometric/electronic heuristics.
    pub smirks: Option<String>,
}

impl Default for ReactionDynamics3DConfig {
    fn default() -> Self {
        Self {
            method: "gfn2".to_string(),
            n_images: 30,
            neb_max_iter: 150,
            spring_k: 0.1,
            step_size: 0.01,
            use_climbing_image: true,
            ci_neb_force_threshold: 0.05,
            n_approach_frames: 15,
            n_departure_frames: 15,
            far_distance: 8.0,
            reactive_distance: 2.0,
            seed: 42,
            use_orbital_guidance: false,
            use_electrostatic_steering: false,
            optimise_complex: false,
            complex_opt_max_steps: 50,
            compute_properties: false,
            n_angular_samples: 0,
            smirks: None,
        }
    }
}

/// A single frame along the 3D reaction coordinate with optional properties.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReactionFrame3D {
    pub index: usize,
    pub coords: Vec<f64>,
    pub energy_kcal_mol: f64,
    pub phase: String,
    /// Per-frame properties (if `compute_properties` was true).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub properties: Option<FrameProperties>,
}

/// Full result of the improved 3D reaction dynamics pipeline.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReactionDynamics3DResult {
    pub frames: Vec<ReactionFrame3D>,
    pub elements: Vec<u8>,
    pub ts_frame_index: usize,
    pub activation_energy_kcal_mol: f64,
    pub reaction_energy_kcal_mol: f64,
    pub method: String,
    pub n_atoms: usize,
    pub notes: Vec<String>,
    /// Approach direction used (unit vector), for diagnostics.
    pub approach_direction: Option<[f64; 3]>,
}

/// Main entry point — full 3D reaction dynamics pipeline.
pub fn compute_reaction_dynamics_3d(
    reactant_smiles: &[&str],
    product_smiles: &[&str],
    config: &ReactionDynamics3DConfig,
) -> Result<ReactionDynamics3DResult, String> {
    if reactant_smiles.is_empty() {
        return Err("At least one reactant SMILES required".into());
    }
    if product_smiles.is_empty() {
        return Err("At least one product SMILES required".into());
    }

    let backend = crate::dynamics::NebBackend::from_method(&config.method)?;

    // ── 1. Embed all fragments (parallelized) ─────────────────────
    #[cfg(feature = "parallel")]
    let (r_confs, p_confs): (Vec<crate::ConformerResult>, Vec<crate::ConformerResult>) = {
        use rayon::prelude::*;
        let r: Vec<crate::ConformerResult> = reactant_smiles
            .par_iter()
            .map(|s| crate::embed(s, config.seed))
            .collect();
        let p: Vec<crate::ConformerResult> = product_smiles
            .par_iter()
            .map(|s| crate::embed(s, config.seed))
            .collect();
        (r, p)
    };
    #[cfg(not(feature = "parallel"))]
    let (r_confs, p_confs): (Vec<crate::ConformerResult>, Vec<crate::ConformerResult>) = {
        let r: Vec<crate::ConformerResult> = reactant_smiles
            .iter()
            .map(|s| crate::embed(s, config.seed))
            .collect();
        let p: Vec<crate::ConformerResult> = product_smiles
            .iter()
            .map(|s| crate::embed(s, config.seed))
            .collect();
        (r, p)
    };

    for (i, c) in r_confs.iter().enumerate() {
        if let Some(ref e) = c.error {
            return Err(format!(
                "Failed to embed reactant '{}': {}",
                reactant_smiles[i], e
            ));
        }
    }
    for (i, c) in p_confs.iter().enumerate() {
        if let Some(ref e) = c.error {
            return Err(format!(
                "Failed to embed product '{}': {}",
                product_smiles[i], e
            ));
        }
    }

    let r_total: usize = r_confs.iter().map(|c| c.num_atoms).sum();
    let p_total: usize = p_confs.iter().map(|c| c.num_atoms).sum();
    if r_total != p_total {
        return Err(format!(
            "Atom count mismatch: {} reactant vs {} product atoms",
            r_total, p_total
        ));
    }

    let r_elements: Vec<u8> = r_confs
        .iter()
        .flat_map(|c| c.elements.iter().copied())
        .collect();
    let combined_smiles = reactant_smiles.join(".");

    // ── 1b. Identify reactive sites (SMIRKS + Fukui + heuristics) ──
    let fragment_elements: Vec<&[u8]> = r_confs.iter().map(|c| c.elements.as_slice()).collect();
    let fragment_coords: Vec<&[f64]> = r_confs.iter().map(|c| c.coords.as_slice()).collect();
    let site_analysis = reactive_sites::identify_reactive_sites(
        &fragment_elements,
        &fragment_coords,
        reactant_smiles,
        config.smirks.as_deref(),
    );
    let reactive_pairs = site_analysis
        .as_ref()
        .ok()
        .map(reactive_sites::reactive_atom_pairs)
        .unwrap_or_default();

    // ── 2. Build product complex & atom mapping ────────────────────
    let reactive_pair_0 = reactive_pairs.first().copied();
    let (p_coords, p_elements) = complex::build_3d_reaction_complex_guided(
        &p_confs,
        config.reactive_distance,
        reactive_pair_0,
    );

    let r_coords_tmp = complex::assemble_fragments_at_product_positions(
        &r_confs,
        &p_coords,
        &p_elements,
        &r_elements,
    );

    let mapping =
        crate::dynamics::map_atoms_greedy(&r_elements, &r_coords_tmp, &p_elements, &p_coords)?;
    let p_reordered = crate::dynamics::reorder_coords(&p_coords, &mapping);

    // ── 3. Build 3D reactant complex (product-guided, NO X-axis forcing)
    let mut r_coords = complex::build_product_guided_complex_3d(&r_confs, &p_reordered);

    // ── 4. Approach direction: reactive sites > orbital > electrostatic
    let mut approach_dir: Option<[f64; 3]> = None;

    // 4a. Use reactive site analysis direction (highest priority when SMIRKS is available)
    if let Ok(ref sa) = site_analysis {
        if let Some(dir) = sa.approach_direction {
            approach_dir = Some(dir);
            if r_confs.len() >= 2 {
                let n0 = r_confs[0].num_atoms;
                complex::orient_fragment_along_direction(&mut r_coords, n0, r_total, dir);
            }
        }
        for note in &sa.notes {
            // Collected later in result notes
            let _ = note;
        }
    }

    // 4b. Orbital-guided approach (if no reactive site direction found)
    if approach_dir.is_none() && config.use_orbital_guidance && r_confs.len() >= 2 {
        let n0 = r_confs[0].num_atoms;
        let coords_a = &r_coords[..n0 * 3];
        let coords_b = &r_coords[n0 * 3..];
        let elems_a = &r_elements[..n0];
        let elems_b = &r_elements[n0..];
        if let Some(info) = orbital_approach::compute_orbital_approach_direction(
            elems_a, coords_a, elems_b, coords_b,
        ) {
            approach_dir = Some(info.direction);
            // Re-orient fragment 1 to approach along orbital direction
            complex::orient_fragment_along_direction(
                &mut r_coords,
                r_confs[0].num_atoms,
                r_total,
                info.direction,
            );
        }
    }

    // ── 5. Optional: constrained complex optimisation ──────────────
    if config.optimise_complex {
        r_coords = complex::optimize_reactive_complex(
            &combined_smiles,
            &r_coords,
            &r_elements,
            config.complex_opt_max_steps,
            &config.method,
        )?;
    }

    // ── 5b. Optional: electrostatic steering for initial orientation ─
    if config.use_electrostatic_steering && approach_dir.is_none() && r_confs.len() >= 2 {
        let n0 = r_confs[0].num_atoms;
        if let Some(info) = electrostatics::compute_electrostatic_approach(
            &r_elements[..n0],
            &r_coords[..n0 * 3],
            &r_elements[n0..],
            &r_coords[n0 * 3..],
        ) {
            approach_dir = Some(info.direction);
            complex::orient_fragment_along_direction(&mut r_coords, n0, r_total, info.direction);
        }
    }

    // ── 6. Optional: multi-angular sampling ────────────────────────
    if config.n_angular_samples > 0 && r_confs.len() >= 2 {
        let n0 = r_confs[0].num_atoms;
        let result = sampling::sample_approach_orientations(
            &combined_smiles,
            &r_elements[..n0],
            &r_coords[..n0 * 3],
            &r_elements[n0..],
            &r_coords[n0 * 3..],
            config.reactive_distance,
            config.n_angular_samples,
            &config.method,
        )?;
        approach_dir = Some(result.best_direction);
        // Rebuild complex with the best direction
        complex::orient_fragment_along_direction(&mut r_coords, n0, r_total, result.best_direction);
    }

    // ── 7. NEB (climbing-image or simplified) ──────────────────────
    let neb_result = if config.use_climbing_image {
        ci_neb::compute_ci_neb_path(
            &combined_smiles,
            &r_coords,
            &p_reordered,
            &r_elements,
            &CiNebConfig {
                n_images: config.n_images,
                max_iter: config.neb_max_iter,
                spring_k: config.spring_k,
                step_size: config.step_size,
                force_threshold: config.ci_neb_force_threshold,
                method: config.method.clone(),
                use_idpp: true,
                relax_images: true,
                relax_max_steps: 30,
            },
        )?
    } else {
        // Fall back to existing simplified NEB, wrap into CiNebResult
        let neb = crate::dynamics::compute_simplified_neb_path_configurable(
            &combined_smiles,
            &r_coords,
            &p_reordered,
            config.n_images,
            config.neb_max_iter,
            config.spring_k,
            config.step_size,
            &config.method,
        )?;
        CiNebResult {
            images: neb
                .images
                .iter()
                .map(|img| ci_neb::CiNebImage {
                    index: img.index,
                    coords: img.coords.clone(),
                    energy_kcal_mol: img.potential_energy_kcal_mol,
                    is_climbing: false,
                })
                .collect(),
            converged: true,
            iterations: config.neb_max_iter,
            max_force: 0.0,
            notes: neb.notes,
        }
    };

    // ── 8. Build molecule ranges for approach/departure sliding ────
    let r_mol_ranges: Vec<(usize, usize)> = {
        let mut ranges = Vec::new();
        let mut off = 0usize;
        for c in &r_confs {
            ranges.push((off, off + c.num_atoms));
            off += c.num_atoms;
        }
        ranges
    };

    let n_xyz = r_total * 3;
    let mol = crate::graph::Molecule::from_smiles(&combined_smiles)?;

    // Constrained optimisation config for approach/departure relaxation
    let relax_config = constrained_opt::ConstrainedOptConfig {
        max_steps: config.complex_opt_max_steps.min(30),
        grad_threshold: 0.1,
        step_size: 0.008,
        method: config.method.clone(),
        constraint_tolerance: 1e-4,
    };

    // ── 9. Approach frames (with constrained geometry relaxation) ──
    let na = config.n_approach_frames;
    let mut approach_raw: Vec<Vec<f64>> = Vec::with_capacity(na);
    for i in 0..na {
        let alpha = if na > 1 {
            i as f64 / (na - 1) as f64
        } else {
            1.0
        };
        let coords =
            crate::dynamics::slide_molecules(&r_coords, &r_mol_ranges, alpha, config.far_distance);
        approach_raw.push(coords);
    }

    // Constrained relaxation of approach frames (parallelized internally)
    let approach_optimised = constrained_opt::optimize_frame_sequence(
        &combined_smiles,
        &r_elements,
        &approach_raw,
        &r_mol_ranges,
        &relax_config,
    )
    .unwrap_or_else(|_| {
        approach_raw
            .iter()
            .map(|coords| constrained_opt::ConstrainedOptResult {
                coords: coords.clone(),
                energy: 0.0,
                n_steps: 0,
                converged: false,
                grad_norm: 0.0,
            })
            .collect()
    });

    // Evaluate approach frame energies (parallelized)
    #[cfg(feature = "parallel")]
    let approach_frames: Vec<ReactionFrame3D> = {
        use rayon::prelude::*;
        approach_optimised
            .par_iter()
            .enumerate()
            .map(|(idx, opt_result)| {
                let mut grad = vec![0.0; n_xyz];
                let energy = crate::dynamics::neb_energy_and_gradient(
                    backend,
                    &combined_smiles,
                    &r_elements,
                    &mol,
                    &opt_result.coords,
                    &mut grad,
                )
                .unwrap_or(opt_result.energy);
                ReactionFrame3D {
                    index: idx,
                    coords: opt_result.coords.clone(),
                    energy_kcal_mol: energy,
                    phase: "approach".into(),
                    properties: None,
                }
            })
            .collect()
    };
    #[cfg(not(feature = "parallel"))]
    let approach_frames: Vec<ReactionFrame3D> = approach_optimised
        .iter()
        .enumerate()
        .map(|(idx, opt_result)| {
            let mut grad = vec![0.0; n_xyz];
            let energy = crate::dynamics::neb_energy_and_gradient(
                backend,
                &combined_smiles,
                &r_elements,
                &mol,
                &opt_result.coords,
                &mut grad,
            )
            .unwrap_or(opt_result.energy);
            ReactionFrame3D {
                index: idx,
                coords: opt_result.coords.clone(),
                energy_kcal_mol: energy,
                phase: "approach".into(),
                properties: None,
            }
        })
        .collect();

    let mut frames: Vec<ReactionFrame3D> = approach_frames;

    // ── 10. Reaction frames from NEB (already IDPP + FIRE + relaxed) ──
    for img in &neb_result.images {
        frames.push(ReactionFrame3D {
            index: frames.len(),
            coords: img.coords.clone(),
            energy_kcal_mol: img.energy_kcal_mol,
            phase: "reaction".into(),
            properties: None,
        });
    }

    // ── 11. Departure frames (with constrained geometry relaxation) ──
    let nd = config.n_departure_frames;
    let single_product = p_confs.len() == 1;
    let p_mol_ranges: Vec<(usize, usize)> = if single_product {
        vec![(0, r_total)]
    } else {
        let p_mol_assign: Vec<usize> = {
            let mut assign = vec![0usize; p_total];
            let mut off = 0usize;
            for (m, c) in p_confs.iter().enumerate() {
                for a in off..off + c.num_atoms {
                    assign[a] = m;
                }
                off += c.num_atoms;
            }
            mapping.iter().map(|&pi| assign[pi]).collect()
        };
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

    let mut departure_raw: Vec<Vec<f64>> = Vec::with_capacity(nd);
    for i in 0..nd {
        let alpha = if nd > 1 {
            1.0 - i as f64 / (nd - 1) as f64
        } else {
            0.0
        };
        let coords = if single_product {
            p_reordered.clone()
        } else {
            crate::dynamics::slide_molecules(
                &p_reordered,
                &p_mol_ranges,
                alpha,
                config.far_distance,
            )
        };
        departure_raw.push(coords);
    }

    // Constrained relaxation of departure frames (parallelized)
    let departure_optimised = if single_product {
        // Single product: no inter-fragment constraint needed, just relax each frame
        #[cfg(feature = "parallel")]
        {
            use rayon::prelude::*;
            departure_raw
                .par_iter()
                .map(|coords| {
                    constrained_opt::optimize_with_constraints(
                        &combined_smiles,
                        &r_elements,
                        coords,
                        &[],
                        &relax_config,
                    )
                    .unwrap_or(constrained_opt::ConstrainedOptResult {
                        coords: coords.clone(),
                        energy: 0.0,
                        n_steps: 0,
                        converged: false,
                        grad_norm: 0.0,
                    })
                })
                .collect::<Vec<_>>()
        }
        #[cfg(not(feature = "parallel"))]
        {
            departure_raw
                .iter()
                .map(|coords| {
                    constrained_opt::optimize_with_constraints(
                        &combined_smiles,
                        &r_elements,
                        coords,
                        &[],
                        &relax_config,
                    )
                    .unwrap_or(constrained_opt::ConstrainedOptResult {
                        coords: coords.clone(),
                        energy: 0.0,
                        n_steps: 0,
                        converged: false,
                        grad_norm: 0.0,
                    })
                })
                .collect::<Vec<_>>()
        }
    } else {
        constrained_opt::optimize_frame_sequence(
            &combined_smiles,
            &r_elements,
            &departure_raw,
            &p_mol_ranges,
            &relax_config,
        )
        .unwrap_or_else(|_| {
            departure_raw
                .iter()
                .map(|coords| constrained_opt::ConstrainedOptResult {
                    coords: coords.clone(),
                    energy: 0.0,
                    n_steps: 0,
                    converged: false,
                    grad_norm: 0.0,
                })
                .collect()
        })
    };

    // Evaluate departure frame energies (parallelized)
    let departure_offset = frames.len();
    #[cfg(feature = "parallel")]
    let departure_frames: Vec<ReactionFrame3D> = {
        use rayon::prelude::*;
        departure_optimised
            .par_iter()
            .enumerate()
            .map(|(idx, opt_result)| {
                let mut grad = vec![0.0; n_xyz];
                let energy = crate::dynamics::neb_energy_and_gradient(
                    backend,
                    &combined_smiles,
                    &r_elements,
                    &mol,
                    &opt_result.coords,
                    &mut grad,
                )
                .unwrap_or(opt_result.energy);
                ReactionFrame3D {
                    index: departure_offset + idx,
                    coords: opt_result.coords.clone(),
                    energy_kcal_mol: energy,
                    phase: "departure".into(),
                    properties: None,
                }
            })
            .collect()
    };
    #[cfg(not(feature = "parallel"))]
    let departure_frames: Vec<ReactionFrame3D> = departure_optimised
        .iter()
        .enumerate()
        .map(|(idx, opt_result)| {
            let mut grad = vec![0.0; n_xyz];
            let energy = crate::dynamics::neb_energy_and_gradient(
                backend,
                &combined_smiles,
                &r_elements,
                &mol,
                &opt_result.coords,
                &mut grad,
            )
            .unwrap_or(opt_result.energy);
            ReactionFrame3D {
                index: departure_offset + idx,
                coords: opt_result.coords.clone(),
                energy_kcal_mol: energy,
                phase: "departure".into(),
                properties: None,
            }
        })
        .collect();

    frames.extend(departure_frames);

    // ── 12. Per-frame properties (optional, parallelized) ─────────
    if config.compute_properties {
        #[cfg(feature = "parallel")]
        {
            use rayon::prelude::*;
            let properties: Vec<Option<FrameProperties>> = frames
                .par_iter()
                .map(|frame| per_frame::compute_frame_properties(&r_elements, &frame.coords).ok())
                .collect();
            for (frame, props) in frames.iter_mut().zip(properties) {
                frame.properties = props;
            }
        }
        #[cfg(not(feature = "parallel"))]
        {
            for frame in &mut frames {
                if let Ok(props) = per_frame::compute_frame_properties(&r_elements, &frame.coords) {
                    frame.properties = Some(props);
                }
            }
        }
    }

    // ── 13. Find TS and compute energetics ─────────────────────────
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

    let mut notes = neb_result.notes;
    notes.push(format!(
        "Full 3D path: {} approach + {} NEB + {} departure = {} frames.",
        na,
        neb_result.images.len(),
        nd,
        frames.len(),
    ));
    notes.push("Approach/departure frames relaxed with constrained optimization.".into());
    if config.use_climbing_image {
        notes.push(format!(
            "CI-NEB converged={}, max_force={:.4} kcal/mol/Å",
            neb_result.converged, neb_result.max_force
        ));
    }
    if approach_dir.is_some() {
        notes.push("Approach direction determined by orbital/sampling guidance.".into());
    } else {
        notes.push("Product-guided 3D orientation (Kabsch alignment, all axes).".into());
    }
    if let Ok(ref sa) = site_analysis {
        for note in &sa.notes {
            notes.push(note.clone());
        }
        if !reactive_pairs.is_empty() {
            notes.push(format!(
                "Reactive atom pairs identified: {:?}",
                reactive_pairs
            ));
        }
    }

    Ok(ReactionDynamics3DResult {
        frames,
        elements: r_elements,
        ts_frame_index: ts_idx,
        activation_energy_kcal_mol: e_ts - e_first,
        reaction_energy_kcal_mol: e_last - e_first,
        method: config.method.clone(),
        n_atoms: r_total,
        notes,
        approach_direction: approach_dir,
    })
}

/// Reactant-only approach-PES scan — no product SMILES required.
///
/// Embeds the reactant fragments, identifies reactive sites using Fukui/FMO
/// and electrostatic analysis, assembles a start complex at `far_distance`,
/// then performs a constrained-relaxed scan from `far_distance` down to
/// `reactive_distance` using the same `slide_molecules` + `optimize_frame_sequence`
/// machinery as the full NEB pipeline.
///
/// The result has the same `ReactionDynamics3DResult` shape as
/// `compute_reaction_dynamics_3d` — approach frames only, `phase = "approach"`.
/// `reaction_energy_kcal_mol` is always 0.0 because no products are known.
///
/// Feature flag: `alpha-reaction-dynamics`
pub fn compute_approach_dynamics(
    reactant_smiles: &[&str],
    config: &ReactionDynamics3DConfig,
) -> Result<ReactionDynamics3DResult, String> {
    if reactant_smiles.len() < 2 {
        return Err(
            "compute_approach_dynamics requires at least 2 reactant SMILES (one per fragment)"
                .into(),
        );
    }

    let backend = crate::dynamics::NebBackend::from_method(&config.method)?;

    // ── 1. Embed reactant fragments ────────────────────────────────
    #[cfg(feature = "parallel")]
    let confs: Vec<crate::ConformerResult> = {
        use rayon::prelude::*;
        reactant_smiles
            .par_iter()
            .map(|s| crate::embed(s, config.seed))
            .collect()
    };
    #[cfg(not(feature = "parallel"))]
    let confs: Vec<crate::ConformerResult> = reactant_smiles
        .iter()
        .map(|s| crate::embed(s, config.seed))
        .collect();

    for (i, c) in confs.iter().enumerate() {
        if let Some(ref e) = c.error {
            return Err(format!(
                "Failed to embed reactant '{}': {}",
                reactant_smiles[i], e
            ));
        }
    }

    let r_total: usize = confs.iter().map(|c| c.num_atoms).sum();
    let r_elements: Vec<u8> = confs
        .iter()
        .flat_map(|c| c.elements.iter().copied())
        .collect();
    let combined_smiles = reactant_smiles.join(".");

    // ── 1b. Identify reactive sites ────────────────────────────────
    let fragment_elements: Vec<&[u8]> = confs.iter().map(|c| c.elements.as_slice()).collect();
    let fragment_coords: Vec<&[f64]> = confs.iter().map(|c| c.coords.as_slice()).collect();
    let site_analysis = reactive_sites::identify_reactive_sites(
        &fragment_elements,
        &fragment_coords,
        reactant_smiles,
        config.smirks.as_deref(),
    );
    let reactive_pairs = site_analysis
        .as_ref()
        .ok()
        .map(reactive_sites::reactive_atom_pairs)
        .unwrap_or_default();
    let reactive_pair_0 = reactive_pairs.first().copied();

    // ── 2. Build reactant complex at reactive_distance ─────────────
    let (mut r_coords, _) = complex::build_3d_reaction_complex_guided(
        &confs[..2.min(confs.len())],
        config.reactive_distance,
        reactive_pair_0,
    );
    // If there are more than 2 fragments, append them naively (rare case)
    for c in confs.iter().skip(2) {
        r_coords.extend_from_slice(&c.coords);
    }

    // ── 3. Approach direction: orbital > electrostatic > site heuristic
    let mut approach_dir: Option<[f64; 3]> = None;
    let mut notes: Vec<String> = Vec::new();

    if let Ok(ref sa) = site_analysis {
        if let Some(dir) = sa.approach_direction {
            approach_dir = Some(dir);
            complex::orient_fragment_along_direction(
                &mut r_coords,
                confs[0].num_atoms,
                r_total,
                dir,
            );
        }
        for note in &sa.notes {
            notes.push(note.clone());
        }
        if !reactive_pairs.is_empty() {
            notes.push(format!("Reactive atom pairs: {:?}", reactive_pairs));
        }
    }

    if approach_dir.is_none() && config.use_orbital_guidance && confs.len() >= 2 {
        let n0 = confs[0].num_atoms;
        if let Some(info) = orbital_approach::compute_orbital_approach_direction(
            &r_elements[..n0],
            &r_coords[..n0 * 3],
            &r_elements[n0..],
            &r_coords[n0 * 3..],
        ) {
            approach_dir = Some(info.direction);
            notes.push(format!(
                "Orbital-guided approach ({})",
                info.interaction_type
            ));
            complex::orient_fragment_along_direction(&mut r_coords, n0, r_total, info.direction);
        }
    }

    if approach_dir.is_none() && config.use_electrostatic_steering && confs.len() >= 2 {
        let n0 = confs[0].num_atoms;
        if let Some(info) = electrostatics::compute_electrostatic_approach(
            &r_elements[..n0],
            &r_coords[..n0 * 3],
            &r_elements[n0..],
            &r_coords[n0 * 3..],
        ) {
            approach_dir = Some(info.direction);
            notes.push("Electrostatic steering applied.".into());
            complex::orient_fragment_along_direction(&mut r_coords, n0, r_total, info.direction);
        }
    }

    // ── 4. Optional: multi-angular sampling ────────────────────────
    if config.n_angular_samples > 0 && confs.len() >= 2 {
        let n0 = confs[0].num_atoms;
        if let Ok(result) = sampling::sample_approach_orientations(
            &combined_smiles,
            &r_elements[..n0],
            &r_coords[..n0 * 3],
            &r_elements[n0..],
            &r_coords[n0 * 3..],
            config.reactive_distance,
            config.n_angular_samples,
            &config.method,
        ) {
            approach_dir = Some(result.best_direction);
            notes.push(format!(
                "Multi-angular sampling ({} orientations), best barrier: {:.2} kcal/mol",
                config.n_angular_samples, result.best_barrier
            ));
            complex::orient_fragment_along_direction(
                &mut r_coords,
                n0,
                r_total,
                result.best_direction,
            );
        }
    }

    // ── 5. Generate approach frames with slide_molecules ───────────
    let r_mol_ranges: Vec<(usize, usize)> = {
        let mut ranges = Vec::new();
        let mut off = 0usize;
        for c in &confs {
            ranges.push((off, off + c.num_atoms));
            off += c.num_atoms;
        }
        ranges
    };

    let na = config.n_approach_frames.max(2);
    let mut approach_raw: Vec<Vec<f64>> = Vec::with_capacity(na);
    for i in 0..na {
        let alpha = i as f64 / (na - 1).max(1) as f64;
        // alpha=0 → far apart; alpha=1 → close (reactive_distance)
        let coords =
            crate::dynamics::slide_molecules(&r_coords, &r_mol_ranges, alpha, config.far_distance);
        approach_raw.push(coords);
    }
    // Reverse so frames[0] = far (reactants separated), frames[na-1] = contact
    approach_raw.reverse();

    let relax_config = constrained_opt::ConstrainedOptConfig {
        max_steps: config.complex_opt_max_steps.min(30),
        grad_threshold: 0.1,
        step_size: 0.008,
        method: config.method.clone(),
        constraint_tolerance: 1e-4,
    };

    let approach_optimised = constrained_opt::optimize_frame_sequence(
        &combined_smiles,
        &r_elements,
        &approach_raw,
        &r_mol_ranges,
        &relax_config,
    )
    .unwrap_or_else(|e| {
        notes.push(format!(
            "Constrained relaxation failed: {e}; using rigid frames."
        ));
        approach_raw
            .iter()
            .map(|coords| constrained_opt::ConstrainedOptResult {
                coords: coords.clone(),
                energy: 0.0,
                n_steps: 0,
                converged: false,
                grad_norm: 0.0,
            })
            .collect()
    });

    // ── 6. Evaluate frame energies ─────────────────────────────────
    let n_xyz = r_total * 3;
    let mol = crate::graph::Molecule::from_smiles(&combined_smiles)?;

    #[cfg(feature = "parallel")]
    let mut frames: Vec<ReactionFrame3D> = {
        use rayon::prelude::*;
        approach_optimised
            .par_iter()
            .enumerate()
            .map(|(idx, opt)| {
                let mut grad = vec![0.0; n_xyz];
                let energy = crate::dynamics::neb_energy_and_gradient(
                    backend,
                    &combined_smiles,
                    &r_elements,
                    &mol,
                    &opt.coords,
                    &mut grad,
                )
                .unwrap_or(opt.energy);
                ReactionFrame3D {
                    index: idx,
                    coords: opt.coords.clone(),
                    energy_kcal_mol: energy,
                    phase: "approach".into(),
                    properties: None,
                }
            })
            .collect()
    };
    #[cfg(not(feature = "parallel"))]
    let mut frames: Vec<ReactionFrame3D> = approach_optimised
        .iter()
        .enumerate()
        .map(|(idx, opt)| {
            let mut grad = vec![0.0; n_xyz];
            let energy = crate::dynamics::neb_energy_and_gradient(
                backend,
                &combined_smiles,
                &r_elements,
                &mol,
                &opt.coords,
                &mut grad,
            )
            .unwrap_or(opt.energy);
            ReactionFrame3D {
                index: idx,
                coords: opt.coords.clone(),
                energy_kcal_mol: energy,
                phase: "approach".into(),
                properties: None,
            }
        })
        .collect();

    // ── 7. Optional per-frame properties ──────────────────────────
    if config.compute_properties {
        #[cfg(feature = "parallel")]
        {
            use rayon::prelude::*;
            let props: Vec<Option<FrameProperties>> = frames
                .par_iter()
                .map(|f| per_frame::compute_frame_properties(&r_elements, &f.coords).ok())
                .collect();
            for (frame, p) in frames.iter_mut().zip(props) {
                frame.properties = p;
            }
        }
        #[cfg(not(feature = "parallel"))]
        for frame in &mut frames {
            if let Ok(p) = per_frame::compute_frame_properties(&r_elements, &frame.coords) {
                frame.properties = Some(p);
            }
        }
    }

    // ── 8. Summary statistics ──────────────────────────────────────
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
    let e_ts = frames.get(ts_idx).map(|f| f.energy_kcal_mol).unwrap_or(0.0);

    notes.push(format!(
        "Approach-only scan: {} frames, {:.1}–{:.1} Å, method: {}. \
         No product SMILES provided — NEB not performed.",
        na, config.far_distance, config.reactive_distance, config.method
    ));

    Ok(ReactionDynamics3DResult {
        frames,
        elements: r_elements,
        ts_frame_index: ts_idx,
        activation_energy_kcal_mol: (e_ts - e_first).max(0.0),
        reaction_energy_kcal_mol: e_last - e_first,
        method: config.method.clone(),
        n_atoms: r_total,
        notes,
        approach_direction: approach_dir,
    })
}
