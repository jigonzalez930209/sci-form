//! Fluorescence from the same sTDA single excitations used for absorption.
//!
//! The emitting state follows Kasha's rule: the lowest singlet, unless that
//! state is dark, in which case the lowest bright state is used and the
//! departure is recorded. The band shape is the spontaneous-emission
//! lineshape, I(E) ∝ n² E³ |μ|² g(E), not a copy of the absorption envelope.
//! The radiative rate is the Einstein A coefficient. Excited-state geometry
//! relaxation is not computed; a Stokes shift is applied only when supplied.

use serde::{Deserialize, Serialize};

use crate::reactivity::{collect_stda_excitations, BroadeningType, StdaExcitation};

const EV_TO_NM: f64 = 1239.841984;
const EV_TO_CM: f64 = 8065.543937;
/// Einstein A prefactor: A = PREFACTOR · n² · ν̃² · f, with ν̃ in cm⁻¹.
const EINSTEIN_A_CM2_PER_S: f64 = 0.667025;

/// Settings for the simulated fluorescence band.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FluorescenceConfig {
    /// Width parameter of the emission line shape, in **eV**.
    ///
    /// For the Gaussian profile this is the standard deviation; for the
    /// Lorentzian profile it is the half-width at half maximum.
    pub sigma: f64,
    /// Lower end of the emission energy window, in **eV**.
    pub e_min: f64,
    /// Upper end of the emission energy window, in **eV**.
    pub e_max: f64,
    /// Number of grid points sampled across the energy window.
    pub n_points: usize,
    /// Extra red-shift in eV applied to the vertical emitting energy.
    /// Zero means vertical emission: this path does not optimize S1.
    pub stokes_ev: f64,
    /// Refractive index of the medium the emission propagates in
    /// (dimensionless); it scales the emission intensity and the radiative
    /// rate with `n²`.
    pub refractive_index: f64,
    /// Oscillator strength below which a state is treated as dark.
    pub dark_threshold: f64,
    /// Line shape used for the emission band.
    pub broadening: BroadeningType,
}

impl Default for FluorescenceConfig {
    fn default() -> Self {
        Self {
            sigma: 0.25,
            e_min: 0.5,
            e_max: 8.0,
            n_points: 400,
            stokes_ev: 0.0,
            refractive_index: 1.0,
            dark_threshold: 1e-3,
            broadening: BroadeningType::Gaussian,
        }
    }
}

/// The simulated fluorescence band of one molecule.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FluorescenceSpectrum {
    /// Emission energy grid, in **eV**, ascending from `e_min` to `e_max`.
    pub energies_ev: Vec<f64>,
    /// Wavelength grid corresponding to [`FluorescenceSpectrum::energies_ev`],
    /// in **nm** (`hc / E`); 0 for energies at or below 0.01 eV.
    pub wavelengths_nm: Vec<f64>,
    /// Relative emission intensity, normalized so the maximum is 1.
    pub intensity: Vec<f64>,
    /// Vertical excitation energy of the emitting state, in **eV**.
    pub absorption_energy_ev: f64,
    /// Emission energy of the emitting state after the Stokes shift, in **eV**.
    pub emission_energy_ev: f64,
    /// Wavelength of the emission maximum, in **nm**.
    pub emission_wavelength_nm: f64,
    /// Oscillator strength of the emitting transition (dimensionless).
    pub oscillator_strength: f64,
    /// Einstein A coefficient in s⁻¹.
    pub radiative_rate_hz: f64,
    /// Radiative lifetime in ns. Infinite when the emitting state is dark.
    pub radiative_lifetime_ns: f64,
    /// True when the emitter is the lowest singlet and it is bright.
    pub kasha: bool,
    /// Index of the occupied molecular orbital the excitation starts from.
    pub from_mo: usize,
    /// Index of the virtual molecular orbital the excitation ends in.
    pub to_mo: usize,
    /// Notes describing the approximations that were made.
    pub notes: Vec<String>,
}

fn lineshape(broadening: BroadeningType, x: f64, center: f64, width: f64) -> f64 {
    let width = width.max(1e-6);
    let dx = x - center;
    match broadening {
        BroadeningType::Gaussian => {
            let norm = 1.0 / (width * (2.0 * std::f64::consts::PI).sqrt());
            norm * (-0.5 * dx * dx / (width * width)).exp()
        }
        BroadeningType::Lorentzian => width / (std::f64::consts::PI * (dx * dx + width * width)),
    }
}

/// A (s⁻¹) = 0.667025 n² ν̃² f, ν̃ in cm⁻¹.
pub fn einstein_a_coefficient(
    oscillator_strength: f64,
    energy_ev: f64,
    refractive_index: f64,
) -> f64 {
    if oscillator_strength <= 0.0 || energy_ev <= 0.0 {
        return 0.0;
    }
    let nu = energy_ev * EV_TO_CM;
    let n = refractive_index.max(1.0);
    EINSTEIN_A_CM2_PER_S * n * n * nu * nu * oscillator_strength
}

fn select_emitter(
    excitations: &[StdaExcitation],
    dark_threshold: f64,
) -> Result<(usize, bool), String> {
    if excitations.is_empty() {
        return Err("No sTDA excitations available for fluorescence".to_string());
    }
    let mut order: Vec<usize> = (0..excitations.len()).collect();
    order.sort_by(|&a, &b| {
        excitations[a]
            .energy_ev
            .partial_cmp(&excitations[b].energy_ev)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let lowest = order[0];
    if excitations[lowest].oscillator_strength >= dark_threshold {
        return Ok((lowest, true));
    }
    for index in order.into_iter().skip(1) {
        if excitations[index].oscillator_strength >= dark_threshold {
            return Ok((index, false));
        }
    }
    Err("Every sTDA singlet in the window is dark; no fluorescence band is emitted".to_string())
}

/// Simulate the fluorescence band of a molecule from its sTDA singlets.
///
/// `elements` and `positions` describe the molecule in the usual convention
/// of this crate: atomic numbers and Cartesian coordinates in **Ångström**.
///
/// The emitting state follows Kasha's rule: the lowest singlet is used when
/// it is bright (oscillator strength ≥ `dark_threshold`); otherwise the
/// lowest bright singlet is used and the departure is recorded in
/// `notes` and flagged through [`FluorescenceSpectrum::kasha`].
///
/// The band shape is the spontaneous-emission lineshape
/// `I(E) ∝ n² E³ |μ|² g(E)`, normalised to a maximum of 1, and the
/// radiative rate is the Einstein A coefficient
/// `A = 0.667025 · n² · ν̃² · f` (ν̃ in cm⁻¹, `A` in s⁻¹). Excited-state
/// geometry relaxation is *not* computed: `stokes_ev` is an optional
/// user-supplied red shift, and 0 means vertical emission.
///
/// # Errors
///
/// Returns `Err` when
///
/// * `config.stokes_ev` is negative,
/// * `config.refractive_index` is smaller than 1,
/// * the sTDA calculation fails (the error of the underlying
///   [`crate::reactivity::collect_stda_excitations`] is propagated unchanged),
/// * the sTDA window contains no singlet at all, or every singlet in it is
///   dark, or
/// * the Stokes shift is so large that the emitting state falls to
///   0.05 eV or below.
///
/// # Example
///
/// ```
/// use sci_form::spectroscopy::fluorescence::{compute_fluorescence_spectrum, FluorescenceConfig};
///
/// // Formaldehyde, H2C=O.
/// let elements = [6u8, 8, 1, 1];
/// let positions = [
///     [0.0, 0.0, 0.0],
///     [1.22, 0.0, 0.0],
///     [-0.55, 0.94, 0.0],
///     [-0.55, -0.94, 0.0],
/// ];
/// let spectrum = compute_fluorescence_spectrum(&elements, &positions, &FluorescenceConfig::default())
///     .expect("formaldehyde has a bright singlet");
/// // Emission energies are in eV, the radiative rate in s^-1.
/// assert!(spectrum.emission_energy_ev > 0.0);
/// assert!(spectrum.radiative_rate_hz > 0.0);
/// ```
pub fn compute_fluorescence_spectrum(
    elements: &[u8],
    positions: &[[f64; 3]],
    config: &FluorescenceConfig,
) -> Result<FluorescenceSpectrum, String> {
    if config.stokes_ev < 0.0 {
        return Err("stokes_ev must be >= 0".to_string());
    }
    if config.refractive_index < 1.0 {
        return Err("refractive_index must be >= 1".to_string());
    }
    let window = (config.e_max + config.stokes_ev + 2.0 * config.sigma).max(config.e_max);
    let excitations = collect_stda_excitations(elements, positions, window)?;
    let (emitter_index, kasha) = select_emitter(&excitations, config.dark_threshold)?;
    let emitter = &excitations[emitter_index];
    let emission_energy = emitter.energy_ev - config.stokes_ev;
    if emission_energy <= 0.05 {
        return Err(format!(
            "Stokes shift {:.3} eV removes the emitting state at {:.3} eV",
            config.stokes_ev, emitter.energy_ev
        ));
    }

    let n_points = config.n_points.max(2);
    let span = (config.e_max - config.e_min).max(1e-6);
    let step = span / (n_points as f64 - 1.0);
    let energies_ev: Vec<f64> = (0..n_points)
        .map(|i| config.e_min + step * i as f64)
        .collect();
    let wavelengths_nm: Vec<f64> = energies_ev
        .iter()
        .map(|&e| if e > 0.01 { EV_TO_NM / e } else { 0.0 })
        .collect();

    // I(E) ∝ n² E³ |μ|² g(E). With f = (2/3) ΔE |μ|², |μ|² ∝ f / ΔE_abs.
    let dipole_weight = emitter.oscillator_strength / emitter.energy_ev.max(1e-6);
    let n2 = config.refractive_index * config.refractive_index;
    let mut intensity: Vec<f64> = energies_ev
        .iter()
        .map(|&e| {
            let e3 = e * e * e;
            n2 * e3 * dipole_weight * lineshape(config.broadening, e, emission_energy, config.sigma)
        })
        .collect();
    let peak = intensity.iter().cloned().fold(0.0_f64, f64::max);
    if peak > 0.0 {
        for value in &mut intensity {
            *value /= peak;
        }
    }

    let rate = einstein_a_coefficient(
        emitter.oscillator_strength,
        emission_energy,
        config.refractive_index,
    );
    let lifetime_ns = if rate > 0.0 {
        1.0e9 / rate
    } else {
        f64::INFINITY
    };

    let mut notes = vec![
        "Fluorescence uses the same EHT sTDA singlets as the absorption spectrum.".to_string(),
        "The band is the spontaneous-emission lineshape I(E) ∝ n² E³ |μ|² g(E), normalized to 1.".to_string(),
        "The radiative rate is the Einstein A coefficient. Quantum yield is not reported: non-radiative rates are not computed.".to_string(),
        "No excited-state geometry optimization is performed. stokes_ev = 0 is vertical emission.".to_string(),
        "Metals missing from the Hoffmann table use hydrogenic EHT parameters (ζ from the ionization potential, f electrons in the core). Open shells keep the SOMO occupied. This is still not TD-DFT.".to_string(),
    ];
    if kasha {
        notes.push("Emitting state is the lowest singlet (Kasha).".to_string());
    } else {
        notes.push(format!(
            "Lowest singlet is dark (f < {}). Emission is taken from the lowest bright singlet; this violates Kasha's rule.",
            config.dark_threshold
        ));
    }
    if config.stokes_ev > 0.0 {
        notes.push(format!(
            "Applied a supplied Stokes shift of {:.3} eV. It is not derived from an S1 optimization.",
            config.stokes_ev
        ));
    }

    Ok(FluorescenceSpectrum {
        energies_ev,
        wavelengths_nm,
        intensity,
        absorption_energy_ev: emitter.energy_ev,
        emission_energy_ev: emission_energy,
        emission_wavelength_nm: EV_TO_NM / emission_energy,
        oscillator_strength: emitter.oscillator_strength,
        radiative_rate_hz: rate,
        radiative_lifetime_ns: lifetime_ns,
        kasha,
        from_mo: emitter.from_mo,
        to_mo: emitter.to_mo,
        notes,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn einstein_a_matches_the_cm_prefactor() {
        let rate = einstein_a_coefficient(1.0, 3.1, 1.0);
        let nu = 3.1 * EV_TO_CM;
        let expected = EINSTEIN_A_CM2_PER_S * nu * nu;
        assert!((rate - expected).abs() / expected < 1e-12);
        let lifetime_ns = 1.0e9 / rate;
        assert!((lifetime_ns - 2.4).abs() < 0.15);
    }

    #[test]
    fn formaldehyde_emits_below_or_at_the_vertical_gap() {
        let elements = [6u8, 8, 1, 1];
        let positions = [
            [0.0, 0.0, 0.0],
            [1.22, 0.0, 0.0],
            [-0.55, 0.94, 0.0],
            [-0.55, -0.94, 0.0],
        ];
        let config = FluorescenceConfig {
            stokes_ev: 0.15,
            ..Default::default()
        };
        let spectrum = compute_fluorescence_spectrum(&elements, &positions, &config)
            .expect("formaldehyde should have a bright singlet");
        assert!(spectrum.emission_energy_ev < spectrum.absorption_energy_ev);
        assert!(
            (spectrum.emission_wavelength_nm - EV_TO_NM / spectrum.emission_energy_ev).abs() < 1e-6
        );
        assert!(spectrum.intensity.iter().any(|v| (*v - 1.0).abs() < 1e-9));
        assert!(spectrum.radiative_rate_hz.is_finite() && spectrum.radiative_rate_hz > 0.0);
        assert!(spectrum.radiative_lifetime_ns.is_finite());
        assert!(spectrum.notes.iter().any(|n| n.contains("Einstein")));
    }

    fn fake_excitation(energy_ev: f64, oscillator_strength: f64) -> StdaExcitation {
        StdaExcitation {
            energy_ev,
            wavelength_nm: EV_TO_NM / energy_ev,
            oscillator_strength,
            from_mo: 0,
            to_mo: 1,
            transition_dipole: 0.0,
        }
    }

    #[test]
    fn all_dark_singlets_return_an_error() {
        let excitations = [fake_excitation(3.0, 1e-6), fake_excitation(4.0, 5e-7)];
        let err = select_emitter(&excitations, 1e-3).unwrap_err();
        assert!(err.contains("dark"));
    }

    #[test]
    fn dark_s1_selects_the_next_bright_singlet() {
        let excitations = [fake_excitation(2.5, 1e-6), fake_excitation(3.5, 0.2)];
        let (index, kasha) = select_emitter(&excitations, 1e-3).unwrap();
        assert!(!kasha);
        assert_eq!(index, 1);
    }

    #[test]
    fn rejects_negative_stokes_and_subunity_index() {
        let elements = [1u8, 1];
        let positions = [[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]];
        let config = FluorescenceConfig {
            stokes_ev: -0.1,
            ..Default::default()
        };
        assert!(compute_fluorescence_spectrum(&elements, &positions, &config).is_err());
        let config = FluorescenceConfig {
            stokes_ev: 0.0,
            refractive_index: 0.9,
            ..Default::default()
        };
        assert!(compute_fluorescence_spectrum(&elements, &positions, &config).is_err());
    }

    #[test]
    fn nickel_carbonyl_hydrogenic_eht_emits_a_band() {
        let elements = [28u8, 6, 8];
        let positions = [[0.0, 0.0, 0.0], [1.85, 0.0, 0.0], [3.00, 0.0, 0.0]];
        let spectrum =
            compute_fluorescence_spectrum(&elements, &positions, &FluorescenceConfig::default())
                .expect("NiCO should produce a screening fluorescence band");
        assert!(spectrum.intensity.iter().any(|v| *v > 0.0));
        assert!(spectrum.notes.iter().any(|n| n.contains("hydrogenic")));
    }
}
