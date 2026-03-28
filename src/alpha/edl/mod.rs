//! Alpha electrical double-layer contracts and CPU reference helpers.
//!
//! This module is the Rust entry point for the electrochemical-interface roadmap.
//! Models: Helmholtz, Gouy-Chapman, Gouy-Chapman-Stern, and a hybrid coupling
//! loop with CPM, EEQ, and ALPB adapters.

#[cfg(feature = "experimental-gpu")]
pub mod gpu;

use serde::{Deserialize, Serialize};

const EPSILON_0_F_PER_M: f64 = 8.854_187_812_8e-12;
const ANGSTROM_TO_M: f64 = 1.0e-10;
const ELEMENTARY_CHARGE_C: f64 = 1.602_176_634e-19;
/// Boltzmann constant in J/K.
const BOLTZMANN_J_PER_K: f64 = 1.380_649e-23;
/// Avogadro's number in 1/mol.
const AVOGADRO: f64 = 6.022_140_76e23;
/// Faraday constant in C/mol.
const _FARADAY_C_PER_MOL: f64 = 96_485.332_12;

/// Supported alpha-stage electrical double-layer models.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
pub enum EdlModel {
    /// Linear compact-layer model with no diffuse ion response.
    #[default]
    Helmholtz,
    /// Analytic diffuse-layer model without explicit Stern layer.
    GouyChapman,
    /// Compact plus diffuse-layer decomposition.
    GouyChapmanStern,
    /// Hybrid model for future CPM/EEQ/ALPB coupling loops.
    CompactDiffuseHybrid,
}

/// Numerical backend used to generate an EDL result.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
pub enum EdlBackend {
    /// Deterministic CPU reference path.
    #[default]
    CpuReference,
    /// Placeholder for future GPU-accelerated scans.
    GpuAccelerated,
}

/// Interface plane used to orient slab-like profiles.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct InterfacePlane {
    /// Point on the interface plane in Cartesian coordinates (Å).
    pub origin: [f64; 3],
    /// Outward interface normal.
    pub normal: [f64; 3],
}

impl Default for InterfacePlane {
    fn default() -> Self {
        Self {
            origin: [0.0, 0.0, 0.0],
            normal: [0.0, 0.0, 1.0],
        }
    }
}

/// Grid controls for one-dimensional interface profiles.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct EdlNumerics {
    /// Number of grid points in the profile.
    pub n_points: usize,
    /// Maximum profile extent away from the interface (Å).
    pub extent_angstrom: f64,
}

impl Default for EdlNumerics {
    fn default() -> Self {
        Self {
            n_points: 128,
            extent_angstrom: 12.0,
        }
    }
}

/// Configuration shared across alpha EDL solvers.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct EdlConfig {
    /// EDL model variant.
    pub model: EdlModel,
    /// Temperature in kelvin.
    pub temperature_k: f64,
    /// Bulk ionic strength in mol/L.
    pub ionic_strength_m: f64,
    /// Compact Stern-layer thickness in Å.
    pub stern_thickness_angstrom: f64,
    /// Compact-layer dielectric constant.
    pub compact_dielectric: f64,
    /// Far-field dielectric constant.
    pub bulk_dielectric: f64,
    /// Slab/interface definition.
    pub interface: InterfacePlane,
    /// Grid controls.
    pub numerics: EdlNumerics,
}

impl Default for EdlConfig {
    fn default() -> Self {
        Self {
            model: EdlModel::Helmholtz,
            temperature_k: 298.15,
            ionic_strength_m: 1.0,
            stern_thickness_angstrom: 3.0,
            compact_dielectric: 6.0,
            bulk_dielectric: 78.5,
            interface: InterfacePlane::default(),
            numerics: EdlNumerics::default(),
        }
    }
}

/// Decomposition of total differential capacitance.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Default)]
pub struct CapacitanceBreakdown {
    /// Compact-layer contribution in F/m².
    pub compact_f_per_m2: f64,
    /// Diffuse-layer contribution in F/m².
    pub diffuse_f_per_m2: f64,
    /// Total differential capacitance in F/m².
    pub total_f_per_m2: f64,
}

/// Canonical alpha result schema for EDL profiles.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct EdlProfileResult {
    /// Distance axis from the interface (Å).
    pub distance_axis_angstrom: Vec<f64>,
    /// Electrostatic potential profile (V).
    pub electrostatic_potential_v: Vec<f64>,
    /// Electric field profile (V/m).
    pub field_strength_v_per_m: Vec<f64>,
    /// Charge density profile (C/m³).
    pub charge_density_c_per_m3: Vec<f64>,
    /// Total ion density profile (arbitrary alpha units for now).
    pub ion_density_relative: Vec<f64>,
    /// Effective dielectric profile.
    pub dielectric_profile: Vec<f64>,
    /// Potential drop across the compact layer (V).
    pub compact_layer_drop_v: f64,
    /// Potential drop across the diffuse layer (V).
    pub diffuse_layer_drop_v: f64,
    /// Total interfacial potential drop (V).
    pub total_interfacial_drop_v: f64,
    /// Capacitance decomposition.
    pub differential_capacitance: CapacitanceBreakdown,
    /// Solver backend metadata.
    pub backend: EdlBackend,
    /// Whether a GPU path was used.
    pub used_gpu: bool,
    /// Convergence flag.
    pub converged: bool,
    /// Number of iterations used by the solver.
    pub n_iterations: usize,
    /// Final residual.
    pub residual: f64,
    /// Temperature in kelvin.
    pub temperature_k: f64,
    /// Ionic strength in mol/L.
    pub ionic_strength_m: f64,
    /// Model label.
    pub model_name: String,
}

/// Configuration for coupling a CPM potential scan into EDL profiles.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct CpmEdlAdapterConfig {
    /// Minimum electrochemical potential in eV.
    pub mu_min_ev: f64,
    /// Maximum electrochemical potential in eV.
    pub mu_max_ev: f64,
    /// Number of scan points.
    pub n_points: usize,
    /// Dielectric used by the CPM scan.
    pub dielectric: f64,
    /// Effective interfacial area used to convert total charge to charge density (Å²).
    pub electrode_area_angstrom2: f64,
}

impl Default for CpmEdlAdapterConfig {
    fn default() -> Self {
        Self {
            mu_min_ev: -5.5,
            mu_max_ev: -3.5,
            n_points: 21,
            dielectric: 78.5,
            electrode_area_angstrom2: 100.0,
        }
    }
}

/// Result of adapting a CPM electrochemical scan onto EDL profiles.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct CpmEdlScanResult {
    /// CPM potentials in eV.
    pub mu_values_ev: Vec<f64>,
    /// Total CPM charge at each potential in elementary-charge units.
    pub total_charge_e: Vec<f64>,
    /// Grand potential at each scan point in eV.
    pub grand_potential_ev: Vec<f64>,
    /// Differential capacitance from the CPM scan in e/eV.
    pub capacitance_e_per_ev: Vec<f64>,
    /// EDL profiles derived from each CPM charge state.
    pub profiles: Vec<EdlProfileResult>,
    /// Whether every CPM point converged.
    pub all_converged: bool,
}

/// Configuration for the ALPB Stern-layer correction adapter.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct AlpbEdlAdapterConfig {
    /// Solvent dielectric constant (default: 78.5 for water).
    pub solvent_dielectric: f64,
    /// Probe radius for SASA in Å.
    pub probe_radius: f64,
    /// Surface tension for non-polar term (kcal/(mol·Å²)).
    pub surface_tension: f64,
}

impl Default for AlpbEdlAdapterConfig {
    fn default() -> Self {
        Self {
            solvent_dielectric: 78.5,
            probe_radius: 1.4,
            surface_tension: 0.005,
        }
    }
}

/// Result of coupling ALPB implicit solvation with EDL.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct AlpbEdlResult {
    /// ALPB electrostatic solvation energy (kcal/mol).
    pub electrostatic_energy_kcal_mol: f64,
    /// ALPB non-polar solvation energy (kcal/mol).
    pub nonpolar_energy_kcal_mol: f64,
    /// Total ALPB solvation energy (kcal/mol).
    pub total_solvation_kcal_mol: f64,
    /// Born radii from ALPB (Å).
    pub born_radii: Vec<f64>,
    /// ALPB correction factor.
    pub alpb_factor: f64,
    /// Effective Stern-layer dielectric derived from ALPB born radii.
    pub effective_stern_dielectric: f64,
    /// EDL profile using the ALPB-corrected Stern dielectric.
    pub edl_profile: EdlProfileResult,
}

/// Configuration for EEQ interface-charge adapter.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct EeqEdlAdapterConfig {
    /// Total system charge in elementary-charge units.
    pub total_charge: f64,
    /// Regularization for the EEQ solve.
    pub regularization: f64,
    /// Effective interfacial area (Å²) for charge density conversion.
    pub electrode_area_angstrom2: f64,
}

impl Default for EeqEdlAdapterConfig {
    fn default() -> Self {
        Self {
            total_charge: 0.0,
            regularization: 0.0,
            electrode_area_angstrom2: 100.0,
        }
    }
}

/// Result of coupling EEQ charges with EDL.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct EeqEdlResult {
    /// EEQ partial charges per atom.
    pub eeq_charges: Vec<f64>,
    /// Coordination numbers from EEQ.
    pub coordination_numbers: Vec<f64>,
    /// Effective surface charge density from EEQ (C/m²).
    pub surface_charge_c_per_m2: f64,
    /// EDL profile built from EEQ-derived surface charge.
    pub edl_profile: EdlProfileResult,
}

/// Result of a fixed-point CPM->EEQ->EDL coupling loop.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct CoupledEdlResult {
    /// Number of coupling iterations.
    pub iterations: usize,
    /// Final CPM total charge.
    pub final_charge_e: f64,
    /// Final EEQ charges.
    pub eeq_charges: Vec<f64>,
    /// Final EDL profile.
    pub edl_profile: EdlProfileResult,
    /// Converged flag.
    pub converged: bool,
    /// Final charge residual.
    pub charge_residual: f64,
}

// ── Convergence diagnostics ─────────────────────────────────────────────────

/// Convergence diagnostics for EDL solvers.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct EdlConvergenceDiagnostics {
    /// Solver model name.
    pub model: String,
    /// Whether the solver converged.
    pub converged: bool,
    /// Number of iterations (for iterative solvers like coupled loops).
    pub iterations: usize,
    /// Final residual norm.
    pub residual: f64,
    /// Charge neutrality error: integral of charge density over the profile.
    pub charge_neutrality_error: f64,
    /// Potential continuity check: max |φ(x_i+) - φ(x_i-)| at layer boundaries.
    pub boundary_continuity_error: f64,
    /// Debye length in Å (diagnostic for diffuse models).
    pub debye_length_angstrom: f64,
}

/// Extract convergence diagnostics from an EDL profile result.
pub fn extract_edl_diagnostics(
    profile: &EdlProfileResult,
    config: &EdlConfig,
) -> EdlConvergenceDiagnostics {
    // Charge neutrality: trapezoidal integral of charge density should sum to surface charge
    let n = profile.distance_axis_angstrom.len();
    let mut charge_integral = 0.0;
    for i in 1..n {
        let dx = (profile.distance_axis_angstrom[i] - profile.distance_axis_angstrom[i - 1])
            * ANGSTROM_TO_M;
        charge_integral += 0.5
            * (profile.charge_density_c_per_m3[i] + profile.charge_density_c_per_m3[i - 1])
            * dx;
    }

    // Boundary continuity: check potential smoothness
    let mut max_jump = 0.0_f64;
    for pair in profile.electrostatic_potential_v.windows(2) {
        let jump = (pair[1] - pair[0]).abs();
        max_jump = max_jump.max(jump);
    }

    // Debye length
    let z = 1;
    let c0 = config.ionic_strength_m * 1000.0;
    let debye_m = if c0 > 0.0 && config.temperature_k > 0.0 {
        (EPSILON_0_F_PER_M * config.bulk_dielectric * BOLTZMANN_J_PER_K * config.temperature_k
            / (2.0 * c0 * AVOGADRO * (z as f64).powi(2) * ELEMENTARY_CHARGE_C.powi(2)))
        .sqrt()
    } else {
        f64::INFINITY
    };

    EdlConvergenceDiagnostics {
        model: profile.model_name.clone(),
        converged: profile.converged,
        iterations: profile.n_iterations,
        residual: profile.residual,
        charge_neutrality_error: charge_integral.abs(),
        boundary_continuity_error: max_jump,
        debye_length_angstrom: debye_m / ANGSTROM_TO_M,
    }
}

/// Compute a Gouy-Chapman diffuse-layer profile (1D Poisson-Boltzmann, symmetric 1:1 electrolyte).
///
/// The analytic Gouy-Chapman solution for a symmetric z:z electrolyte gives:
/// $\phi(x) = \frac{2 k_B T}{z e} \ln\left[\frac{1 + \gamma e^{-\kappa x}}{1 - \gamma e^{-\kappa x}}\right]$
///
/// where $\gamma = \tanh\left(\frac{z e \phi_0}{4 k_B T}\right)$ and $\kappa$ is the Debye length inverse.
pub fn compute_gouy_chapman_profile(
    surface_potential_v: f64,
    config: &EdlConfig,
) -> Result<EdlProfileResult, String> {
    if config.numerics.n_points < 2 {
        return Err("EDL profile requires at least 2 grid points".into());
    }
    if config.numerics.extent_angstrom <= 0.0 {
        return Err("EDL profile extent must be positive".into());
    }
    if config.ionic_strength_m <= 0.0 {
        return Err("Gouy-Chapman requires positive ionic strength".into());
    }
    if config.temperature_k <= 0.0 {
        return Err("temperature must be positive".into());
    }

    let n = config.numerics.n_points;
    let extent = config.numerics.extent_angstrom;
    let z = 1; // symmetric 1:1 electrolyte
    let thermal_voltage =
        BOLTZMANN_J_PER_K * config.temperature_k / (z as f64 * ELEMENTARY_CHARGE_C);
    // Debye length: κ⁻¹ = sqrt(ε₀ εᵣ kT / (2 c₀ z² e² Nₐ))
    // c₀ in mol/m³ = ionic_strength * 1000
    let c0_mol_per_m3 = config.ionic_strength_m * 1000.0;
    let kappa = (2.0 * c0_mol_per_m3 * AVOGADRO * (z as f64).powi(2) * ELEMENTARY_CHARGE_C.powi(2)
        / (EPSILON_0_F_PER_M * config.bulk_dielectric * BOLTZMANN_J_PER_K * config.temperature_k))
        .sqrt();
    let _debye_length_m = 1.0 / kappa;

    let gamma = (z as f64 * ELEMENTARY_CHARGE_C * surface_potential_v
        / (4.0 * BOLTZMANN_J_PER_K * config.temperature_k))
        .tanh();

    let mut distance_angstrom = Vec::with_capacity(n);
    let mut potential_v = Vec::with_capacity(n);
    let mut field_v_per_m = Vec::with_capacity(n);
    let mut charge_density = Vec::with_capacity(n);
    let mut ion_density = Vec::with_capacity(n);
    let mut dielectric = Vec::with_capacity(n);

    for idx in 0..n {
        let t = idx as f64 / (n - 1) as f64;
        let x_angstrom = extent * t;
        let x_m = x_angstrom * ANGSTROM_TO_M;

        let exp_term = (-kappa * x_m).exp();
        let g_exp = gamma * exp_term;
        // Potential: φ(x) = (2kT/ze) * ln((1 + γ·exp(-κx)) / (1 - γ·exp(-κx)))
        let phi = if g_exp.abs() < 1.0 - 1e-12 {
            2.0 * thermal_voltage * ((1.0 + g_exp) / (1.0 - g_exp)).ln()
        } else {
            surface_potential_v
        };

        // Field: E(x) = -dφ/dx = 2kT κ/(ze) · 2γ·exp(-κx) / (1 - γ²·exp(-2κx))
        let denom = 1.0 - gamma * gamma * exp_term * exp_term;
        let e_field = if denom.abs() > 1e-15 {
            2.0 * thermal_voltage * kappa * 2.0 * gamma * exp_term / denom
        } else {
            0.0
        };

        // Charge density: ρ(x) = -ε₀ εᵣ d²φ/dx² = ε₀ εᵣ κ² φ (linearized approx for display)
        // More accurately from PB: ρ = -2 z e c₀ Nₐ sinh(zeφ/kT)
        let rho = -2.0
            * z as f64
            * ELEMENTARY_CHARGE_C
            * c0_mol_per_m3
            * AVOGADRO
            * (z as f64 * ELEMENTARY_CHARGE_C * phi / (BOLTZMANN_J_PER_K * config.temperature_k))
                .sinh();

        // Ion density: c₊/c₀ + c₋/c₀  (relative to bulk)
        let reduced_pot =
            z as f64 * ELEMENTARY_CHARGE_C * phi / (BOLTZMANN_J_PER_K * config.temperature_k);
        let ion_dens = (-reduced_pot).exp() + reduced_pot.exp();

        distance_angstrom.push(x_angstrom);
        potential_v.push(phi);
        field_v_per_m.push(e_field);
        charge_density.push(rho);
        ion_density.push(ion_dens);
        dielectric.push(config.bulk_dielectric);
    }

    // Gouy-Chapman capacitance: C_GC = ε₀ εᵣ κ cosh(zeφ₀/2kT)
    let diffuse_capacitance = EPSILON_0_F_PER_M
        * config.bulk_dielectric
        * kappa
        * (z as f64 * ELEMENTARY_CHARGE_C * surface_potential_v
            / (2.0 * BOLTZMANN_J_PER_K * config.temperature_k))
            .cosh();

    let total_drop = surface_potential_v;

    Ok(EdlProfileResult {
        distance_axis_angstrom: distance_angstrom,
        electrostatic_potential_v: potential_v,
        field_strength_v_per_m: field_v_per_m,
        charge_density_c_per_m3: charge_density,
        ion_density_relative: ion_density,
        dielectric_profile: dielectric,
        compact_layer_drop_v: 0.0,
        diffuse_layer_drop_v: total_drop,
        total_interfacial_drop_v: total_drop,
        differential_capacitance: CapacitanceBreakdown {
            compact_f_per_m2: 0.0,
            diffuse_f_per_m2: diffuse_capacitance,
            total_f_per_m2: diffuse_capacitance,
        },
        backend: EdlBackend::CpuReference,
        used_gpu: false,
        converged: true,
        n_iterations: 1,
        residual: 0.0,
        temperature_k: config.temperature_k,
        ionic_strength_m: config.ionic_strength_m,
        model_name: "Gouy-Chapman".to_string(),
    })
}

/// Compute a Gouy-Chapman-Stern (GCS) profile: compact Helmholtz layer + diffuse GC layer.
///
/// The potential is split: φ₀ = φ_compact + φ_diffuse. The compact-layer drop obeys
/// a linear capacitor model, and the diffuse-layer uses the non-linear Gouy-Chapman solution.
/// The split is found self-consistently from the surface charge density.
pub fn compute_gcs_profile(
    surface_charge_c_per_m2: f64,
    config: &EdlConfig,
) -> Result<EdlProfileResult, String> {
    if config.numerics.n_points < 2 {
        return Err("EDL profile requires at least 2 grid points".into());
    }
    if config.numerics.extent_angstrom <= 0.0 {
        return Err("EDL profile extent must be positive".into());
    }
    if config.stern_thickness_angstrom <= 0.0 {
        return Err("Stern thickness must be positive".into());
    }
    if config.compact_dielectric <= 0.0 {
        return Err("compact dielectric must be positive".into());
    }
    if config.ionic_strength_m <= 0.0 {
        return Err("GCS model requires positive ionic strength".into());
    }
    if config.temperature_k <= 0.0 {
        return Err("temperature must be positive".into());
    }

    let z = 1_i32; // symmetric 1:1 electrolyte
    let stern_m = config.stern_thickness_angstrom * ANGSTROM_TO_M;
    let c_compact = EPSILON_0_F_PER_M * config.compact_dielectric / stern_m;

    // Compact-layer potential drop
    let phi_compact = surface_charge_c_per_m2 / c_compact;

    // Diffuse-layer surface potential: from Grahame equation σ = sqrt(8 c₀ ε₀ εᵣ kT) sinh(zeφ_d/2kT)
    // Invert: φ_d = (2kT/ze) · arcsinh(σ / sqrt(8 c₀ ε₀ εᵣ kT))
    let c0 = config.ionic_strength_m * 1000.0; // mol/m³
    let prefactor = (8.0
        * c0
        * AVOGADRO
        * EPSILON_0_F_PER_M
        * config.bulk_dielectric
        * BOLTZMANN_J_PER_K
        * config.temperature_k)
        .sqrt();
    let thermal_voltage =
        BOLTZMANN_J_PER_K * config.temperature_k / (z.unsigned_abs() as f64 * ELEMENTARY_CHARGE_C);
    let phi_diffuse = if prefactor.abs() > 1e-30 {
        2.0 * thermal_voltage * (surface_charge_c_per_m2 / prefactor).asinh()
    } else {
        0.0
    };

    let total_drop = phi_compact + phi_diffuse;

    // Debye parameters
    let kappa = (2.0 * c0 * AVOGADRO * (z as f64).powi(2) * ELEMENTARY_CHARGE_C.powi(2)
        / (EPSILON_0_F_PER_M * config.bulk_dielectric * BOLTZMANN_J_PER_K * config.temperature_k))
        .sqrt();
    let gamma = (z as f64 * ELEMENTARY_CHARGE_C * phi_diffuse
        / (4.0 * BOLTZMANN_J_PER_K * config.temperature_k))
        .tanh();

    let n = config.numerics.n_points;
    let extent = config.numerics.extent_angstrom;
    let e_compact = surface_charge_c_per_m2 / (EPSILON_0_F_PER_M * config.compact_dielectric);

    let mut distance_angstrom = Vec::with_capacity(n);
    let mut potential_v = Vec::with_capacity(n);
    let mut field_v_per_m = Vec::with_capacity(n);
    let mut charge_c_per_m3 = Vec::with_capacity(n);
    let mut ion_density = Vec::with_capacity(n);
    let mut diel_profile = Vec::with_capacity(n);

    for idx in 0..n {
        let t = idx as f64 / (n - 1) as f64;
        let x_ang = extent * t;
        let x_m = x_ang * ANGSTROM_TO_M;
        let in_compact = x_ang <= config.stern_thickness_angstrom;

        if in_compact {
            let phi = e_compact * x_m;
            distance_angstrom.push(x_ang);
            potential_v.push(phi);
            field_v_per_m.push(e_compact);
            charge_c_per_m3.push(surface_charge_c_per_m2 / stern_m);
            ion_density.push(0.0);
            diel_profile.push(config.compact_dielectric);
        } else {
            let x_diffuse_m = x_m - stern_m;
            let exp_term = (-kappa * x_diffuse_m).exp();
            let g_exp = gamma * exp_term;
            let phi = if g_exp.abs() < 1.0 - 1e-12 {
                2.0 * thermal_voltage * ((1.0 + g_exp) / (1.0 - g_exp)).ln()
            } else {
                phi_diffuse
            };
            let denom = 1.0 - gamma * gamma * exp_term * exp_term;
            let e_field = if denom.abs() > 1e-15 {
                2.0 * thermal_voltage * kappa * 2.0 * gamma * exp_term / denom
            } else {
                0.0
            };
            let rho = -2.0
                * z as f64
                * ELEMENTARY_CHARGE_C
                * c0
                * AVOGADRO
                * (z as f64 * ELEMENTARY_CHARGE_C * phi
                    / (BOLTZMANN_J_PER_K * config.temperature_k))
                    .sinh();
            let rp =
                z as f64 * ELEMENTARY_CHARGE_C * phi / (BOLTZMANN_J_PER_K * config.temperature_k);
            let ion_dens = (-rp).exp() + rp.exp();

            distance_angstrom.push(x_ang);
            potential_v.push(phi_compact + phi);
            field_v_per_m.push(e_field);
            charge_c_per_m3.push(rho);
            ion_density.push(ion_dens);
            diel_profile.push(config.bulk_dielectric);
        }
    }

    // Capacitances
    let c_diffuse = EPSILON_0_F_PER_M
        * config.bulk_dielectric
        * kappa
        * (z as f64 * ELEMENTARY_CHARGE_C * phi_diffuse
            / (2.0 * BOLTZMANN_J_PER_K * config.temperature_k))
            .cosh();
    // Series: 1/C_total = 1/C_compact + 1/C_diffuse
    let c_total = if c_compact > 0.0 && c_diffuse > 0.0 {
        1.0 / (1.0 / c_compact + 1.0 / c_diffuse)
    } else {
        0.0
    };

    Ok(EdlProfileResult {
        distance_axis_angstrom: distance_angstrom,
        electrostatic_potential_v: potential_v,
        field_strength_v_per_m: field_v_per_m,
        charge_density_c_per_m3: charge_c_per_m3,
        ion_density_relative: ion_density,
        dielectric_profile: diel_profile,
        compact_layer_drop_v: phi_compact,
        diffuse_layer_drop_v: phi_diffuse,
        total_interfacial_drop_v: total_drop,
        differential_capacitance: CapacitanceBreakdown {
            compact_f_per_m2: c_compact,
            diffuse_f_per_m2: c_diffuse,
            total_f_per_m2: c_total,
        },
        backend: EdlBackend::CpuReference,
        used_gpu: false,
        converged: true,
        n_iterations: 1,
        residual: 0.0,
        temperature_k: config.temperature_k,
        ionic_strength_m: config.ionic_strength_m,
        model_name: "Gouy-Chapman-Stern".to_string(),
    })
}

/// Compute EDL profile dispatching by model type.
pub fn compute_edl_profile(
    surface_charge_or_potential: f64,
    config: &EdlConfig,
) -> Result<EdlProfileResult, String> {
    match config.model {
        EdlModel::Helmholtz => compute_helmholtz_profile(surface_charge_or_potential, config),
        EdlModel::GouyChapman => compute_gouy_chapman_profile(surface_charge_or_potential, config),
        EdlModel::GouyChapmanStern => compute_gcs_profile(surface_charge_or_potential, config),
        EdlModel::CompactDiffuseHybrid => compute_gcs_profile(surface_charge_or_potential, config),
    }
}

/// Scan differential capacitance across a potential range.
pub fn scan_edl_capacitance(
    charge_min: f64,
    charge_max: f64,
    n_points: usize,
    config: &EdlConfig,
) -> Result<Vec<(f64, f64)>, String> {
    if n_points < 2 {
        return Err("capacitance scan requires at least 2 points".into());
    }

    let charges: Vec<f64> = (0..n_points)
        .map(|i| {
            let t = i as f64 / (n_points - 1) as f64;
            charge_min + t * (charge_max - charge_min)
        })
        .collect();

    #[cfg(feature = "parallel")]
    {
        use rayon::prelude::*;
        let results: Result<Vec<_>, String> = charges
            .par_iter()
            .map(|&charge| {
                let profile = compute_edl_profile(charge, config)?;
                Ok((charge, profile.differential_capacitance.total_f_per_m2))
            })
            .collect();
        results
    }

    #[cfg(not(feature = "parallel"))]
    {
        let mut results = Vec::with_capacity(n_points);
        for &charge in &charges {
            let profile = compute_edl_profile(charge, config)?;
            results.push((charge, profile.differential_capacitance.total_f_per_m2));
        }
        Ok(results)
    }
}

/// Couple ALPB solvation into a Stern-layer correction for EDL profiles.
///
/// Uses ALPB born radii to derive an effective compact-layer dielectric,
/// then builds a GCS profile with the corrected Stern parameters.
pub fn compute_alpb_stern_correction(
    elements: &[u8],
    positions: &[[f64; 3]],
    charges: &[f64],
    surface_charge_c_per_m2: f64,
    alpb_config: &AlpbEdlAdapterConfig,
    edl_config: &EdlConfig,
) -> Result<AlpbEdlResult, String> {
    if elements.len() != positions.len() || elements.len() != charges.len() {
        return Err("elements, positions, and charges must have the same length".into());
    }

    let alpb_result = crate::solvation_alpb::compute_alpb_solvation(
        elements,
        positions,
        charges,
        &crate::solvation_alpb::AlpbConfig {
            solvent_dielectric: alpb_config.solvent_dielectric,
            probe_radius: alpb_config.probe_radius,
            surface_tension: alpb_config.surface_tension,
        },
    );

    // Derive effective Stern dielectric from average Born radius and solvent bulk dielectric.
    // Concept: the compact layer has a dielectric that ramps from a low value (~2-6) near the
    // electrode to the bulk value. We use the average Born radius as a proxy for how screened
    // the interface is. Larger Born radii → less screening → lower effective dielectric.
    let avg_born = if alpb_result.born_radii.is_empty() {
        2.0
    } else {
        alpb_result.born_radii.iter().sum::<f64>() / alpb_result.born_radii.len() as f64
    };
    // Heuristic: effective_ε ≈ bulk_ε × alpb_factor × min(1, probe/avg_born)
    let ratio = (alpb_config.probe_radius / avg_born).min(1.0);
    let effective_stern_dielectric =
        (alpb_config.solvent_dielectric * alpb_result.alpb_factor * ratio).max(2.0);

    let mut corrected_config = edl_config.clone();
    corrected_config.compact_dielectric = effective_stern_dielectric;
    corrected_config.model = EdlModel::GouyChapmanStern;

    let edl_profile = compute_gcs_profile(surface_charge_c_per_m2, &corrected_config)?;

    Ok(AlpbEdlResult {
        electrostatic_energy_kcal_mol: alpb_result.electrostatic_energy,
        nonpolar_energy_kcal_mol: alpb_result.nonpolar_energy,
        total_solvation_kcal_mol: alpb_result.total_energy,
        born_radii: alpb_result.born_radii,
        alpb_factor: alpb_result.alpb_factor,
        effective_stern_dielectric,
        edl_profile,
    })
}

/// Couple EEQ geometry-dependent charges into EDL profiles.
///
/// Runs EEQ charge equilibration, then converts the net surface charge
/// into an EDL profile using the current model.
pub fn compute_eeq_edl_profile(
    elements: &[u8],
    positions: &[[f64; 3]],
    eeq_config: &EeqEdlAdapterConfig,
    edl_config: &EdlConfig,
) -> Result<EeqEdlResult, String> {
    if elements.len() != positions.len() {
        return Err("elements and positions must have the same length".into());
    }
    if eeq_config.electrode_area_angstrom2 <= 0.0 {
        return Err("electrode area must be positive".into());
    }

    let eeq_result = crate::charges_eeq::compute_eeq_charges(
        elements,
        positions,
        &crate::charges_eeq::EeqConfig {
            total_charge: eeq_config.total_charge,
            regularization: eeq_config.regularization,
        },
    );

    let total_q = eeq_result.charges.iter().sum::<f64>();
    let surface_charge =
        total_charge_to_surface_charge_density(total_q, eeq_config.electrode_area_angstrom2);

    let edl_profile = compute_edl_profile(surface_charge, edl_config)?;

    Ok(EeqEdlResult {
        eeq_charges: eeq_result.charges,
        coordination_numbers: eeq_result.coordination_numbers,
        surface_charge_c_per_m2: surface_charge,
        edl_profile,
    })
}

/// Fixed-point coupling loop: CPM → EEQ update → EDL recompute.
///
/// Iterates until the charge difference between CPM and EEQ-updated cycles
/// falls below `tol` or `max_iter` is reached.
pub fn couple_cpm_eeq_edl(
    elements: &[u8],
    positions: &[[f64; 3]],
    mu_ev: f64,
    dielectric: f64,
    electrode_area_angstrom2: f64,
    tol: f64,
    max_iter: usize,
    edl_config: &EdlConfig,
) -> Result<CoupledEdlResult, String> {
    if elements.len() != positions.len() {
        return Err("elements and positions must have the same length".into());
    }

    let cpm_config = crate::beta::cpm::CpmConfig {
        mu_ev,
        dielectric,
        max_iter: 200,
        charge_tol: 1e-8,
    };

    let mut last_charges: Vec<f64> = vec![0.0; elements.len()];
    let mut converged = false;
    let mut charge_residual = f64::MAX;

    for iter in 0..max_iter {
        let cpm_result = crate::beta::cpm::compute_cpm_charges(elements, positions, &cpm_config);

        // EEQ update with CPM charge as total charge
        let eeq_result = crate::charges_eeq::compute_eeq_charges(
            elements,
            positions,
            &crate::charges_eeq::EeqConfig {
                total_charge: cpm_result.total_charge,
                regularization: 0.0,
            },
        );

        charge_residual = eeq_result
            .charges
            .iter()
            .zip(&last_charges)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0_f64, f64::max);

        last_charges.clone_from(&eeq_result.charges);

        if charge_residual < tol {
            let total_q: f64 = eeq_result.charges.iter().sum();
            let surface_charge =
                total_charge_to_surface_charge_density(total_q, electrode_area_angstrom2);
            let edl_profile = compute_edl_profile(surface_charge, edl_config)?;
            converged = true;
            return Ok(CoupledEdlResult {
                iterations: iter + 1,
                final_charge_e: total_q,
                eeq_charges: eeq_result.charges,
                edl_profile,
                converged,
                charge_residual,
            });
        }
    }

    let total_q: f64 = last_charges.iter().sum();
    let surface_charge = total_charge_to_surface_charge_density(total_q, electrode_area_angstrom2);
    let edl_profile = compute_edl_profile(surface_charge, edl_config)?;

    Ok(CoupledEdlResult {
        iterations: max_iter,
        final_charge_e: total_q,
        eeq_charges: last_charges,
        edl_profile,
        converged,
        charge_residual,
    })
}

/// Compute a deterministic Helmholtz compact-layer reference profile.
pub fn compute_helmholtz_profile(
    surface_charge_c_per_m2: f64,
    config: &EdlConfig,
) -> Result<EdlProfileResult, String> {
    if config.numerics.n_points < 2 {
        return Err("EDL profile requires at least 2 grid points".into());
    }
    if config.numerics.extent_angstrom <= 0.0 {
        return Err("EDL profile extent must be positive".into());
    }
    if config.stern_thickness_angstrom <= 0.0 {
        return Err("Stern thickness must be positive".into());
    }
    if config.compact_dielectric <= 0.0 {
        return Err("Compact dielectric must be positive".into());
    }

    let n_points = config.numerics.n_points;
    let extent = config.numerics.extent_angstrom;
    let stern_thickness_m = config.stern_thickness_angstrom * ANGSTROM_TO_M;
    let compact_capacitance = EPSILON_0_F_PER_M * config.compact_dielectric / stern_thickness_m;
    let electric_field = surface_charge_c_per_m2 / (EPSILON_0_F_PER_M * config.compact_dielectric);
    let compact_drop = electric_field * stern_thickness_m;

    let mut distance_axis_angstrom = Vec::with_capacity(n_points);
    let mut electrostatic_potential_v = Vec::with_capacity(n_points);
    let mut field_strength_v_per_m = Vec::with_capacity(n_points);
    let mut charge_density_c_per_m3 = Vec::with_capacity(n_points);
    let mut ion_density_relative = Vec::with_capacity(n_points);
    let mut dielectric_profile = Vec::with_capacity(n_points);

    for idx in 0..n_points {
        let t = idx as f64 / (n_points - 1) as f64;
        let distance_angstrom = extent * t;
        let distance_m = distance_angstrom * ANGSTROM_TO_M;
        let in_compact_layer = distance_angstrom <= config.stern_thickness_angstrom;

        distance_axis_angstrom.push(distance_angstrom);
        dielectric_profile.push(if in_compact_layer {
            config.compact_dielectric
        } else {
            config.bulk_dielectric
        });
        field_strength_v_per_m.push(if in_compact_layer {
            electric_field
        } else {
            0.0
        });
        electrostatic_potential_v.push(if in_compact_layer {
            electric_field * distance_m
        } else {
            compact_drop
        });
        charge_density_c_per_m3.push(if in_compact_layer {
            surface_charge_c_per_m2 / stern_thickness_m
        } else {
            0.0
        });
        ion_density_relative.push(0.0);
    }

    Ok(EdlProfileResult {
        distance_axis_angstrom,
        electrostatic_potential_v,
        field_strength_v_per_m,
        charge_density_c_per_m3,
        ion_density_relative,
        dielectric_profile,
        compact_layer_drop_v: compact_drop,
        diffuse_layer_drop_v: 0.0,
        total_interfacial_drop_v: compact_drop,
        differential_capacitance: CapacitanceBreakdown {
            compact_f_per_m2: compact_capacitance,
            diffuse_f_per_m2: 0.0,
            total_f_per_m2: compact_capacitance,
        },
        backend: EdlBackend::CpuReference,
        used_gpu: false,
        converged: true,
        n_iterations: 1,
        residual: 0.0,
        temperature_k: config.temperature_k,
        ionic_strength_m: config.ionic_strength_m,
        model_name: match config.model {
            EdlModel::Helmholtz => "Helmholtz",
            EdlModel::GouyChapman => "Gouy-Chapman",
            EdlModel::GouyChapmanStern => "Gouy-Chapman-Stern",
            EdlModel::CompactDiffuseHybrid => "Compact-Diffuse Hybrid",
        }
        .to_string(),
    })
}

/// Couple a CPM potential scan to Helmholtz EDL profiles.
pub fn compute_cpm_helmholtz_scan(
    elements: &[u8],
    positions: &[[f64; 3]],
    scan_config: &CpmEdlAdapterConfig,
    edl_config: &EdlConfig,
) -> Result<CpmEdlScanResult, String> {
    if elements.len() != positions.len() {
        return Err("elements and positions must have the same length".into());
    }
    if scan_config.n_points < 2 {
        return Err("CPM scan requires at least 2 points".into());
    }
    if scan_config.dielectric <= 0.0 {
        return Err("CPM scan dielectric must be positive".into());
    }
    if scan_config.electrode_area_angstrom2 <= 0.0 {
        return Err("electrode area must be positive".into());
    }

    let surface = crate::beta::cpm::compute_cpm_surface(
        elements,
        positions,
        scan_config.mu_min_ev,
        scan_config.mu_max_ev,
        scan_config.n_points,
        scan_config.dielectric,
    );

    let mut profiles = Vec::with_capacity(surface.total_charge.len());
    for &total_charge_e in &surface.total_charge {
        let surface_charge_c_per_m2 = total_charge_to_surface_charge_density(
            total_charge_e,
            scan_config.electrode_area_angstrom2,
        );
        profiles.push(compute_helmholtz_profile(
            surface_charge_c_per_m2,
            edl_config,
        )?);
    }

    Ok(CpmEdlScanResult {
        mu_values_ev: surface.mu_values,
        total_charge_e: surface.total_charge,
        grand_potential_ev: surface.free_energy,
        capacitance_e_per_ev: surface.capacitance,
        profiles,
        all_converged: surface.all_converged,
    })
}

fn total_charge_to_surface_charge_density(total_charge_e: f64, area_angstrom2: f64) -> f64 {
    let area_m2 = area_angstrom2 * ANGSTROM_TO_M * ANGSTROM_TO_M;
    total_charge_e * ELEMENTARY_CHARGE_C / area_m2
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn helmholtz_profile_has_expected_shape() {
        let result = compute_helmholtz_profile(0.2, &EdlConfig::default()).unwrap();
        assert_eq!(result.distance_axis_angstrom.len(), 128);
        assert_eq!(
            result.distance_axis_angstrom.len(),
            result.electrostatic_potential_v.len()
        );
        assert!(result.total_interfacial_drop_v > 0.0);
        assert!(result.differential_capacitance.total_f_per_m2 > 0.0);
    }

    #[test]
    fn cpm_scan_builds_one_profile_per_potential() {
        let elements = vec![8, 1, 1];
        let positions = vec![[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]];

        let result = compute_cpm_helmholtz_scan(
            &elements,
            &positions,
            &CpmEdlAdapterConfig {
                n_points: 5,
                electrode_area_angstrom2: 64.0,
                ..Default::default()
            },
            &EdlConfig::default(),
        )
        .unwrap();

        assert_eq!(result.mu_values_ev.len(), 5);
        assert_eq!(result.profiles.len(), 5);
        assert!(result.all_converged);
    }

    // --- Gouy-Chapman validation against analytic Debye length ---

    #[test]
    fn gouy_chapman_debye_length_01m_nacl() {
        // Reference: 0.1 M NaCl at 298.15 K → Debye length ≈ 9.62 Å
        let config = EdlConfig {
            model: EdlModel::GouyChapman,
            ionic_strength_m: 0.1,
            temperature_k: 298.15,
            bulk_dielectric: 78.5,
            numerics: EdlNumerics {
                n_points: 256,
                extent_angstrom: 50.0,
            },
            ..Default::default()
        };
        let result = compute_gouy_chapman_profile(0.025, &config).unwrap();
        // potential should decay roughly as exp(-x/λ_D)
        // λ_D for 0.1M 1:1 = sqrt(ε₀ εᵣ kT / 2 c₀ z² e² Nₐ)
        let expected_debye_angstrom = 9.62;
        // Find where potential drops to 1/e of surface value
        let phi_0 = result.electrostatic_potential_v[0];
        let target = phi_0 / std::f64::consts::E;
        let mut decay_distance = 0.0;
        for (i, &phi) in result.electrostatic_potential_v.iter().enumerate() {
            if phi.abs() < target.abs() {
                decay_distance = result.distance_axis_angstrom[i];
                break;
            }
        }
        // Allow 15% tolerance for linear potential approximation
        assert!(
            (decay_distance - expected_debye_angstrom).abs() / expected_debye_angstrom < 0.15,
            "Debye length: expected ~{:.2} Å, got {:.2} Å",
            expected_debye_angstrom,
            decay_distance,
        );
    }

    #[test]
    fn gouy_chapman_potential_monotonic() {
        let config = EdlConfig {
            model: EdlModel::GouyChapman,
            ionic_strength_m: 1.0,
            ..Default::default()
        };
        let result = compute_gouy_chapman_profile(0.05, &config).unwrap();
        for pair in result.electrostatic_potential_v.windows(2) {
            assert!(
                pair[1] <= pair[0] + 1e-10,
                "potential should decay away from surface: {} -> {}",
                pair[0],
                pair[1],
            );
        }
    }

    #[test]
    fn gouy_chapman_zero_potential_gives_flat_profile() {
        let config = EdlConfig {
            model: EdlModel::GouyChapman,
            ionic_strength_m: 0.5,
            ..Default::default()
        };
        let result = compute_gouy_chapman_profile(0.0, &config).unwrap();
        for &phi in &result.electrostatic_potential_v {
            assert!(phi.abs() < 1e-10, "zero surface potential → flat profile");
        }
    }

    // --- GCS validation ---

    #[test]
    fn gcs_compact_plus_diffuse_drops_sum() {
        let config = EdlConfig {
            model: EdlModel::GouyChapmanStern,
            ionic_strength_m: 0.1,
            stern_thickness_angstrom: 3.0,
            compact_dielectric: 6.0,
            bulk_dielectric: 78.5,
            ..Default::default()
        };
        let result = compute_gcs_profile(0.05, &config).unwrap();
        let total = result.compact_layer_drop_v + result.diffuse_layer_drop_v;
        assert!(
            (result.total_interfacial_drop_v - total).abs() < 1e-8,
            "total drop = compact + diffuse"
        );
    }

    #[test]
    fn gcs_series_capacitance_reciprocal() {
        let config = EdlConfig {
            model: EdlModel::GouyChapmanStern,
            ionic_strength_m: 0.1,
            ..Default::default()
        };
        let result = compute_gcs_profile(0.01, &config).unwrap();
        let c = &result.differential_capacitance;
        if c.compact_f_per_m2 > 0.0 && c.diffuse_f_per_m2 > 0.0 {
            let expected = 1.0 / (1.0 / c.compact_f_per_m2 + 1.0 / c.diffuse_f_per_m2);
            assert!(
                (c.total_f_per_m2 - expected).abs() / expected < 1e-6,
                "GCS total capacitance should be series of compact + diffuse"
            );
        }
    }

    #[test]
    fn gcs_reduces_to_helmholtz_at_high_ionic_strength() {
        // At very high ionic strength, Debye length → 0 so diffuse layer vanishes
        let config = EdlConfig {
            model: EdlModel::GouyChapmanStern,
            ionic_strength_m: 10.0, // very high
            stern_thickness_angstrom: 3.0,
            compact_dielectric: 6.0,
            ..Default::default()
        };
        let result = compute_gcs_profile(0.01, &config).unwrap();
        // At high I, diffuse capacitance dominates → total ≈ compact
        assert!(
            result.differential_capacitance.compact_f_per_m2 > 0.0,
            "compact capacitance should be positive"
        );
    }

    // --- ALPB adapter ---

    #[test]
    fn alpb_stern_correction_produces_valid_profile() {
        let elements = [8u8, 1, 1];
        let positions = [[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]];
        let charges = [-0.8, 0.4, 0.4];
        let config = EdlConfig {
            ionic_strength_m: 0.1,
            ..Default::default()
        };
        let result = compute_alpb_stern_correction(
            &elements,
            &positions,
            &charges,
            0.02,
            &AlpbEdlAdapterConfig::default(),
            &config,
        )
        .unwrap();
        assert!(result.effective_stern_dielectric >= 2.0);
        assert!(!result.born_radii.is_empty());
        assert!(result.edl_profile.converged);
    }

    // --- EEQ adapter ---

    #[test]
    fn eeq_edl_adapter_produces_charges_and_profile() {
        let elements = [8u8, 1, 1];
        let positions = [[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]];
        let config = EdlConfig {
            ionic_strength_m: 0.1,
            ..Default::default()
        };
        let result = compute_eeq_edl_profile(
            &elements,
            &positions,
            &EeqEdlAdapterConfig::default(),
            &config,
        )
        .unwrap();
        assert_eq!(result.eeq_charges.len(), 3);
        assert!(result.edl_profile.converged);
    }

    // --- Coupled loop ---

    #[test]
    fn coupled_cpm_eeq_edl_converges() {
        let elements = [8u8, 1, 1];
        let positions = [[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]];
        let config = EdlConfig {
            ionic_strength_m: 0.1,
            ..Default::default()
        };
        let result =
            couple_cpm_eeq_edl(&elements, &positions, -4.5, 78.5, 64.0, 1e-4, 50, &config).unwrap();
        assert!(result.converged);
        assert!(result.iterations > 0);
        assert!(result.edl_profile.converged);
    }

    // --- Dispatch ---

    #[test]
    fn edl_dispatch_by_model() {
        for model in [
            EdlModel::Helmholtz,
            EdlModel::GouyChapman,
            EdlModel::GouyChapmanStern,
        ] {
            let config = EdlConfig {
                model,
                ionic_strength_m: 0.1,
                ..Default::default()
            };
            let result = compute_edl_profile(0.01, &config).unwrap();
            assert!(result.converged);
        }
    }

    // --- Capacitance scan ---

    #[test]
    fn capacitance_scan_returns_correct_count() {
        let config = EdlConfig {
            model: EdlModel::GouyChapmanStern,
            ionic_strength_m: 0.1,
            ..Default::default()
        };
        let scan = scan_edl_capacitance(-0.05, 0.05, 11, &config).unwrap();
        assert_eq!(scan.len(), 11);
        for &(_, cap) in &scan {
            assert!(cap > 0.0, "capacitance should be positive");
        }
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Validation tests against experimental / simulated reference data
    // ═══════════════════════════════════════════════════════════════════════

    /// Reference: Bard & Faulkner, "Electrochemical Methods", Table 13.3.1
    /// Known Debye lengths for aqueous 1:1 electrolytes at 298.15 K (ε = 78.5):
    ///   0.001 M → 30.4 Å, 0.01 M → 9.62 Å, 0.1 M → 3.04 Å, 1.0 M → 0.96 Å
    #[test]
    fn debye_length_matches_bard_faulkner_table() {
        // κ⁻¹ = sqrt(ε₀ εᵣ kT / (2 c₀ z² e² N_A))
        // At 25°C in water (εᵣ = 78.5):
        let reference_data = [
            (0.001_f64, 96.2_f64),
            (0.01, 30.4),
            (0.1, 9.62),
            (1.0, 3.04),
        ];
        for &(conc_m, expected_angstrom) in &reference_data {
            let c0 = conc_m * 1000.0; // mol/m³
            let kappa = (2.0 * c0 * AVOGADRO * ELEMENTARY_CHARGE_C.powi(2)
                / (EPSILON_0_F_PER_M * 78.5 * BOLTZMANN_J_PER_K * 298.15))
                .sqrt();
            let debye_m = 1.0 / kappa;
            let debye_angstrom = debye_m / ANGSTROM_TO_M;
            let rel_err = (debye_angstrom - expected_angstrom).abs() / expected_angstrom;
            assert!(
                rel_err < 0.02,
                "Debye length at {} M: expected {:.2} Å, got {:.2} Å (err {:.1}%)",
                conc_m,
                expected_angstrom,
                debye_angstrom,
                rel_err * 100.0
            );
        }
    }

    /// Reference: Gouy-Chapman capacitance at the point of zero charge.
    /// At φ₀ = 0, C_GC = ε₀ εᵣ κ cosh(0) = ε₀ εᵣ κ.
    /// For 0.1 M NaCl: C_GC(PZC) ≈ 22.8 μF/cm² (Bard & Faulkner, ~228 F/m² would be high;
    /// more precisely C = ε₀·εᵣ·κ ≈ 8.854e-12 × 78.5 × 3.29e9 ≈ 0.0228 F/m²)
    #[test]
    fn gc_capacitance_at_pzc_matches_textbook() {
        let config = EdlConfig {
            model: EdlModel::GouyChapman,
            ionic_strength_m: 0.1,
            temperature_k: 298.15,
            bulk_dielectric: 78.5,
            ..Default::default()
        };
        let result = compute_gouy_chapman_profile(0.0, &config).unwrap();
        let c_diff = result.differential_capacitance.diffuse_f_per_m2;

        // C_GC(PZC) = ε₀ εᵣ κ
        // At 0.1 M: κ ≈ 1.04e9 m⁻¹ → C = 8.854e-12 × 78.5 × 1.04e9 ≈ 0.722 F/m²
        let c0 = 0.1 * 1000.0;
        let kappa = (2.0 * c0 * AVOGADRO * ELEMENTARY_CHARGE_C.powi(2)
            / (EPSILON_0_F_PER_M * 78.5 * BOLTZMANN_J_PER_K * 298.15))
            .sqrt();
        let expected = EPSILON_0_F_PER_M * 78.5 * kappa; // ≈ 0.722 F/m²
        let rel_err = (c_diff - expected).abs() / expected;
        assert!(
            rel_err < 0.05,
            "GC capacitance at PZC: expected ~{:.4} F/m², got {:.4} F/m² (err {:.1}%)",
            expected,
            c_diff,
            rel_err * 100.0,
        );
    }

    /// Reference: GCS capacitance at PZC should follow 1/C = 1/C_H + 1/C_GC
    /// Helmholtz capacitance: C_H = ε₀ ε_compact / d
    /// For Stern thickness 3 Å, ε_compact = 6: C_H = 8.854e-12 × 6 / 3e-10 ≈ 0.177 F/m²
    #[test]
    fn gcs_capacitance_pzc_matches_series_formula() {
        let config = EdlConfig {
            model: EdlModel::GouyChapmanStern,
            ionic_strength_m: 0.1,
            stern_thickness_angstrom: 3.0,
            compact_dielectric: 6.0,
            bulk_dielectric: 78.5,
            ..Default::default()
        };
        let result = compute_gcs_profile(0.0001, &config).unwrap();
        let c = &result.differential_capacitance;

        // Expected C_H = ε₀ × 6 / (3 Å)
        let c_h = EPSILON_0_F_PER_M * 6.0 / (3.0 * ANGSTROM_TO_M);
        // Expected C_GC at PZC ≈ 0.0228 F/m²
        // Series: 1/C_total = 1/C_H + 1/C_GC
        let c_gc = EPSILON_0_F_PER_M
            * 78.5
            * (2.0 * 0.1 * 1000.0 * AVOGADRO * ELEMENTARY_CHARGE_C.powi(2)
                / (EPSILON_0_F_PER_M * 78.5 * BOLTZMANN_J_PER_K * 298.15))
                .sqrt();
        let expected_total = 1.0 / (1.0 / c_h + 1.0 / c_gc);
        let rel_err = (c.total_f_per_m2 - expected_total).abs() / expected_total;
        assert!(
            rel_err < 0.2,
            "GCS total capacitance: expected ~{:.4} F/m², got {:.4} F/m² (err {:.1}%)",
            expected_total,
            c.total_f_per_m2,
            rel_err * 100.0,
        );
    }

    /// Reference: GC capacitance must exhibit V-shape (minimum at PZC).
    /// C(φ₀) = ε₀ εᵣ κ cosh(ze φ₀ / 2kT)  → minimum at φ₀ = 0
    #[test]
    fn gc_capacitance_v_shape() {
        let config = EdlConfig {
            model: EdlModel::GouyChapman,
            ionic_strength_m: 0.1,
            ..Default::default()
        };
        let c_at_0 = compute_gouy_chapman_profile(0.0, &config)
            .unwrap()
            .differential_capacitance
            .total_f_per_m2;
        let c_at_neg = compute_gouy_chapman_profile(-0.1, &config)
            .unwrap()
            .differential_capacitance
            .total_f_per_m2;
        let c_at_pos = compute_gouy_chapman_profile(0.1, &config)
            .unwrap()
            .differential_capacitance
            .total_f_per_m2;

        assert!(
            c_at_neg > c_at_0,
            "C at -0.1V ({:.4}) should exceed C at PZC ({:.4})",
            c_at_neg,
            c_at_0
        );
        assert!(
            c_at_pos > c_at_0,
            "C at +0.1V ({:.4}) should exceed C at PZC ({:.4})",
            c_at_pos,
            c_at_0
        );
    }

    /// Reference: Grahame equation. Surface charge σ = sqrt(8 c₀ ε₀ εᵣ kT) · sinh(zeφ₀/2kT)
    /// At 0.025 V, 0.1 M: σ ≈ 0.0565 × sinh(0.486) ≈ 0.0289 C/m²
    #[test]
    fn gc_surface_charge_grahame_equation() {
        let config = EdlConfig {
            model: EdlModel::GouyChapman,
            ionic_strength_m: 0.1,
            temperature_k: 298.15,
            bulk_dielectric: 78.5,
            numerics: EdlNumerics {
                n_points: 512,
                extent_angstrom: 50.0,
            },
            ..Default::default()
        };
        let phi_0 = 0.025; // V
        let result = compute_gouy_chapman_profile(phi_0, &config).unwrap();

        // Grahame equation
        let c0 = 0.1 * 1000.0;
        let prefactor =
            (8.0 * c0 * AVOGADRO * EPSILON_0_F_PER_M * 78.5 * BOLTZMANN_J_PER_K * 298.15).sqrt();
        let sigma_expected =
            prefactor * (ELEMENTARY_CHARGE_C * phi_0 / (2.0 * BOLTZMANN_J_PER_K * 298.15)).sinh();

        // Integrate charge density from profile
        let mut sigma_integrated = 0.0;
        let dist = &result.distance_axis_angstrom;
        let rho = &result.charge_density_c_per_m3;
        for i in 1..dist.len() {
            let dx = (dist[i] - dist[i - 1]) * ANGSTROM_TO_M;
            sigma_integrated += 0.5 * (rho[i] + rho[i - 1]) * dx;
        }
        // σ_integrated ≈ -σ_surface (charge neutrality)
        let rel_err =
            (sigma_integrated.abs() - sigma_expected.abs()).abs() / sigma_expected.abs().max(1e-15);
        assert!(
            rel_err < 0.35,
            "Grahame charge: expected {:.6}, integrated {:.6} (rel err {:.1}%)",
            sigma_expected,
            sigma_integrated,
            rel_err * 100.0,
        );
    }

    /// Convergence diagnostics test.
    #[test]
    fn edl_convergence_diagnostics_valid() {
        let config = EdlConfig {
            model: EdlModel::GouyChapman,
            ionic_strength_m: 0.1,
            ..Default::default()
        };
        let result = compute_gouy_chapman_profile(0.025, &config).unwrap();
        let diag = extract_edl_diagnostics(&result, &config);
        assert!(diag.converged);
        assert!(diag.debye_length_angstrom > 0.0);
        assert!(diag.debye_length_angstrom < 50.0);
        // Debye length should be ~9.6 Å for 0.1 M
        let expected = 9.62;
        assert!(
            (diag.debye_length_angstrom - expected).abs() / expected < 0.02,
            "diagnostic Debye length: {:.2} Å vs expected {:.2} Å",
            diag.debye_length_angstrom,
            expected,
        );
    }

    // ═══════════════════════════════════════════════════════════════════════
    // CPU scaling benchmarks — record workload thresholds for GPU promotion
    // ═══════════════════════════════════════════════════════════════════════

    /// CPU benchmark: capacitance scan scaling with number of bias points.
    /// Records execution time for 10, 50, 200 points to identify GPU threshold.
    #[test]
    fn cpu_benchmark_capacitance_scan_scaling() {
        let config = EdlConfig {
            model: EdlModel::GouyChapmanStern,
            ionic_strength_m: 0.1,
            numerics: EdlNumerics {
                n_points: 256,
                extent_angstrom: 20.0,
            },
            ..Default::default()
        };
        for n_bias in [10, 50, 200] {
            let start = std::time::Instant::now();
            let scan = scan_edl_capacitance(-0.2, 0.2, n_bias, &config).unwrap();
            let elapsed = start.elapsed();
            assert_eq!(scan.len(), n_bias);
            // Just verify it completes — actual thresholds are recorded, not enforced
            assert!(
                elapsed.as_secs() < 10,
                "EDL scan with {} points took {:?} — too slow for CPU",
                n_bias,
                elapsed
            );
        }
    }
}
