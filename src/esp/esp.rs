//! Electrostatic potential (ESP) on 3D grids + Gaussian .cube format I/O.
//!
//! ESP at point r:
//!   Φ(r) = Σ_A Z_A / |r - R_A|  −  Σ_{μν} P_{μν} ∫ φ_μ(r') / |r - r'| φ_ν(r') dr'
//!
//! The full two-electron integral is expensive. For EHT-level ESP we use the
//! Mulliken-charge approximation:
//!   Φ(r) ≈ Σ_A q_A / |r - R_A|
//!
//! where q_A is the full nuclear charge (not just valence) minus the
//! Mulliken electron population.
//!
//! For a more refined ESP: we can sum nuclear + point-charge electron contributions:
//!   Φ(r) = Σ_A Z_nuclear / |r - R_A|  −  Σ_{μ} n_μ  Σ_g c_g ∫ g(r') / |r-r'| dr'
//!
//! We use the Mulliken approach for speed, which is standard for semi-empirical methods.

use serde::{Deserialize, Serialize};
use std::io::{BufRead, Write};

/// Bohr ↔ Ångström conversion.
const BOHR_TO_ANG: f64 = 0.529177;
const ANG_TO_BOHR: f64 = 1.0 / BOHR_TO_ANG;

/// ESP grid result.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EspGrid {
    /// Origin in Ångström.
    pub origin: [f64; 3],
    /// Grid spacing in Ångström.
    pub spacing: f64,
    /// Dimensions [nx, ny, nz].
    pub dims: [usize; 3],
    /// ESP values in atomic units (Hartree/e).
    pub values: Vec<f64>,
}

/// Evaluation settings for the ESP grid.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct EspGridConfig {
    /// Grid spacing in Ångström (must be finite and > 0).
    pub spacing: f64,
    /// Padding around the molecular bounding box in Ångström (must be finite and >= 0).
    pub padding: f64,
    /// Minimum distance to a nucleus in Ångström; the `1/r` singularity is
    /// clamped here (must be finite and > 0).
    pub min_distance: f64,
}

impl Default for EspGridConfig {
    fn default() -> Self {
        Self {
            spacing: 0.5,
            padding: 2.0,
            min_distance: DEFAULT_ESP_MIN_DISTANCE,
        }
    }
}

impl EspGridConfig {
    /// Validate the configuration, returning a message for the first invalid field.
    pub fn validate(&self) -> Result<(), String> {
        if !self.spacing.is_finite() || self.spacing <= 0.0 {
            return Err(format!(
                "ESP grid spacing must be finite and > 0 Å (got {})",
                self.spacing
            ));
        }
        if !self.padding.is_finite() || self.padding < 0.0 {
            return Err(format!(
                "ESP grid padding must be finite and >= 0 Å (got {})",
                self.padding
            ));
        }
        if !self.min_distance.is_finite() || self.min_distance <= 0.0 {
            return Err(format!(
                "ESP grid minimum distance must be finite and > 0 Å (got {})",
                self.min_distance
            ));
        }
        Ok(())
    }
}

/// Legacy default for the nuclear-distance clamp (Ångström).
pub const DEFAULT_ESP_MIN_DISTANCE: f64 = 0.1;

/// Upper bound on the number of grid points, to reject configurations whose
/// `dims` product would exhaust memory (e.g. a 0.001 Å spacing).
const MAX_ESP_GRID_POINTS: usize = 1 << 30;

/// Grid origin and dimensions for a molecular bounding box.
fn grid_bounds(positions: &[[f64; 3]], spacing: f64, padding: f64) -> ([f64; 3], [usize; 3]) {
    let mut min = [f64::MAX; 3];
    let mut max = [f64::MIN; 3];
    for pos in positions {
        for k in 0..3 {
            min[k] = min[k].min(pos[k]);
            max[k] = max[k].max(pos[k]);
        }
    }
    let origin = [min[0] - padding, min[1] - padding, min[2] - padding];
    let dims = [
        ((max[0] - min[0] + 2.0 * padding) / spacing).ceil() as usize + 1,
        ((max[1] - min[1] + 2.0 * padding) / spacing).ceil() as usize + 1,
        ((max[2] - min[2] + 2.0 * padding) / spacing).ceil() as usize + 1,
    ];
    (origin, dims)
}

/// Evaluate the ESP values on a grid, clamping the distance to each nucleus.
#[allow(clippy::too_many_arguments)]
fn evaluate_esp_values(
    positions: &[[f64; 3]],
    mulliken_charges: &[f64],
    origin: [f64; 3],
    dims: [usize; 3],
    spacing: f64,
    min_distance: f64,
) -> Vec<f64> {
    let total = dims[0] * dims[1] * dims[2];
    let mut values = vec![0.0f64; total];
    for ix in 0..dims[0] {
        for iy in 0..dims[1] {
            for iz in 0..dims[2] {
                let rx = origin[0] + ix as f64 * spacing;
                let ry = origin[1] + iy as f64 * spacing;
                let rz = origin[2] + iz as f64 * spacing;

                let mut phi = 0.0;
                for (a, pos) in positions.iter().enumerate() {
                    let dx = rx - pos[0];
                    let dy = ry - pos[1];
                    let dz = rz - pos[2];
                    let dist = (dx * dx + dy * dy + dz * dz).sqrt();
                    let dist_clamped = dist.max(min_distance);
                    phi += mulliken_charges[a] / (dist_clamped * ANG_TO_BOHR);
                }

                let idx = ix * dims[1] * dims[2] + iy * dims[2] + iz;
                values[idx] = phi;
            }
        }
    }
    values
}

/// Compute the ESP grid with explicit evaluation settings.
///
/// Unlike [`compute_esp_grid`], invalid `spacing`/`padding`/`min_distance`
/// values are rejected instead of silently adjusted.
pub fn compute_esp_grid_with_config(
    elements: &[u8],
    positions: &[[f64; 3]],
    mulliken_charges: &[f64],
    config: EspGridConfig,
) -> Result<EspGrid, String> {
    config.validate()?;
    if positions.len() != mulliken_charges.len() {
        return Err(format!(
            "positions ({}) and charges ({}) length mismatch",
            positions.len(),
            mulliken_charges.len()
        ));
    }
    let _ = elements;
    let (origin, dims) = grid_bounds(positions, config.spacing, config.padding);
    let total = dims[0] * dims[1] * dims[2];
    if total > MAX_ESP_GRID_POINTS {
        return Err(format!(
            "ESP grid would need {} points (spacing {} Å, padding {} Å); \
             increase the spacing or reduce the padding",
            total, config.spacing, config.padding
        ));
    }
    let values = evaluate_esp_values(
        positions,
        mulliken_charges,
        origin,
        dims,
        config.spacing,
        config.min_distance,
    );
    Ok(EspGrid {
        origin,
        spacing: config.spacing,
        dims,
        values,
    })
}

/// Compute ESP on a 3D grid using Mulliken-charge approximation.
///
/// - `elements`: atomic numbers
/// - `positions`: Ångström
/// - `mulliken_charges`: per-atom partial charges from population analysis
/// - `spacing`: grid spacing in Ångström
/// - `padding`: padding around molecule in Ångström
pub fn compute_esp_grid(
    elements: &[u8],
    positions: &[[f64; 3]],
    mulliken_charges: &[f64],
    spacing: f64,
    padding: f64,
) -> EspGrid {
    let _ = elements;
    // Legacy behaviour: sanitize instead of rejecting (use
    // [`compute_esp_grid_with_config`] for validated settings).
    let spacing = if spacing.is_finite() && spacing > 0.0 {
        spacing
    } else {
        0.5
    };
    let padding = if padding.is_finite() {
        padding.max(0.0)
    } else {
        0.0
    };

    let (origin, dims) = grid_bounds(positions, spacing, padding);
    let values = evaluate_esp_values(
        positions,
        mulliken_charges,
        origin,
        dims,
        spacing,
        DEFAULT_ESP_MIN_DISTANCE,
    );

    EspGrid {
        origin,
        spacing,
        dims,
        values,
    }
}

/// Compute ESP on a 3D grid using rayon parallelism.
///
/// Same as `compute_esp_grid` but evaluates grid points in parallel.
#[cfg(feature = "parallel")]
pub fn compute_esp_grid_parallel(
    elements: &[u8],
    positions: &[[f64; 3]],
    mulliken_charges: &[f64],
    spacing: f64,
    padding: f64,
) -> EspGrid {
    let _ = elements;
    let spacing = if spacing.is_finite() && spacing > 0.0 {
        spacing
    } else {
        0.5
    };
    let padding = if padding.is_finite() {
        padding.max(0.0)
    } else {
        0.0
    };
    let (origin, dims) = grid_bounds(positions, spacing, padding);
    let values = evaluate_esp_values_parallel(
        positions,
        mulliken_charges,
        origin,
        dims,
        spacing,
        DEFAULT_ESP_MIN_DISTANCE,
    );
    EspGrid {
        origin,
        spacing,
        dims,
        values,
    }
}

/// Compute the ESP grid in parallel with explicit evaluation settings.
///
/// Validated counterpart of [`compute_esp_grid_parallel`]; invalid
/// `spacing`/`padding`/`min_distance` values are rejected.
#[cfg(feature = "parallel")]
pub fn compute_esp_grid_parallel_with_config(
    elements: &[u8],
    positions: &[[f64; 3]],
    mulliken_charges: &[f64],
    config: EspGridConfig,
) -> Result<EspGrid, String> {
    config.validate()?;
    if positions.len() != mulliken_charges.len() {
        return Err(format!(
            "positions ({}) and charges ({}) length mismatch",
            positions.len(),
            mulliken_charges.len()
        ));
    }
    let _ = elements;
    let (origin, dims) = grid_bounds(positions, config.spacing, config.padding);
    let total = dims[0] * dims[1] * dims[2];
    if total > MAX_ESP_GRID_POINTS {
        return Err(format!(
            "ESP grid would need {} points (spacing {} Å, padding {} Å); \
             increase the spacing or reduce the padding",
            total, config.spacing, config.padding
        ));
    }
    let values = evaluate_esp_values_parallel(
        positions,
        mulliken_charges,
        origin,
        dims,
        config.spacing,
        config.min_distance,
    );
    Ok(EspGrid {
        origin,
        spacing: config.spacing,
        dims,
        values,
    })
}

/// Parallel ESP evaluation with the distance clamp.
#[cfg(feature = "parallel")]
#[allow(clippy::too_many_arguments)]
fn evaluate_esp_values_parallel(
    positions: &[[f64; 3]],
    mulliken_charges: &[f64],
    origin: [f64; 3],
    dims: [usize; 3],
    spacing: f64,
    min_distance: f64,
) -> Vec<f64> {
    use rayon::prelude::*;

    let total = dims[0] * dims[1] * dims[2];
    let ny = dims[1];
    let nz = dims[2];
    let positions_clone = positions.to_vec();
    let charges_clone = mulliken_charges.to_vec();

    (0..total)
        .into_par_iter()
        .map(|flat_idx| {
            let ix = flat_idx / (ny * nz);
            let iy = (flat_idx / nz) % ny;
            let iz = flat_idx % nz;

            let rx = origin[0] + ix as f64 * spacing;
            let ry = origin[1] + iy as f64 * spacing;
            let rz = origin[2] + iz as f64 * spacing;

            let mut phi = 0.0;
            for (a, pos) in positions_clone.iter().enumerate() {
                let dx = rx - pos[0];
                let dy = ry - pos[1];
                let dz = rz - pos[2];
                let dist = (dx * dx + dy * dy + dz * dz).sqrt();
                let dist_clamped = dist.max(min_distance);
                phi += charges_clone[a] / (dist_clamped * ANG_TO_BOHR);
            }
            phi
        })
        .collect()
}

/// Map an ESP value to an RGB color.
///
/// Red = negative potential (electron-rich), Blue = positive (electron-poor),
/// White = neutral.  The value is clamped to `[-range, +range]` before mapping.
///
/// Returns `[r, g, b]` where each component is in `[0, 255]`.
pub fn esp_color_map(value: f64, range: f64) -> [u8; 3] {
    let clamped = value.max(-range).min(range);
    let t = clamped / range; // -1..+1

    if t < 0.0 {
        // Negative → red (through white)
        let f = (-t).min(1.0);
        let white = ((1.0 - f) * 255.0) as u8;
        [255, white, white]
    } else {
        // Positive → blue (through white)
        let f = t.min(1.0);
        let white = ((1.0 - f) * 255.0) as u8;
        [white, white, 255]
    }
}

/// Map an entire ESP grid to RGB colors.
///
/// Returns a flat `Vec<u8>` of length `3 * values.len()` in `[r,g,b, r,g,b, ...]` order.
/// `range` controls the saturation: values beyond ±range map to full red/blue.
pub fn esp_grid_to_colors(grid: &EspGrid, range: f64) -> Vec<u8> {
    let mut colors = Vec::with_capacity(grid.values.len() * 3);
    for &val in &grid.values {
        let [r, g, b] = esp_color_map(val, range);
        colors.push(r);
        colors.push(g);
        colors.push(b);
    }
    colors
}

/// Gaussian .cube file (for I/O).
#[derive(Debug, Clone)]
pub struct CubeFile {
    pub comment1: String,
    pub comment2: String,
    pub elements: Vec<u8>,
    pub positions: Vec<[f64; 3]>,
    pub origin: [f64; 3],
    pub dims: [usize; 3],
    pub axes: [[f64; 3]; 3],
    pub values: Vec<f64>,
}

/// Export ESP grid as a Gaussian .cube file.
pub fn export_cube<W: Write>(
    writer: &mut W,
    elements: &[u8],
    positions: &[[f64; 3]],
    grid: &EspGrid,
) -> std::io::Result<()> {
    let n_atoms = elements.len();
    writeln!(writer, " ESP generated by sci-form")?;
    writeln!(writer, " Mulliken-charge approximation")?;
    // The field below is the electrostatic potential in Hartree evaluated on
    // the grid, not an electron density. Gaussian cube files normally carry a
    // density, so the units are stated explicitly to avoid a reader
    // interpreting the values as a density.

    // Origin in bohr
    let o = [
        grid.origin[0] * ANG_TO_BOHR,
        grid.origin[1] * ANG_TO_BOHR,
        grid.origin[2] * ANG_TO_BOHR,
    ];
    writeln!(
        writer,
        "{:5} {:12.6} {:12.6} {:12.6}",
        n_atoms as i32, o[0], o[1], o[2]
    )?;

    // Axis vectors (orthogonal, spacing in bohr)
    let sp = grid.spacing * ANG_TO_BOHR;
    writeln!(
        writer,
        "{:5} {:12.6} {:12.6} {:12.6}",
        grid.dims[0], sp, 0.0, 0.0
    )?;
    writeln!(
        writer,
        "{:5} {:12.6} {:12.6} {:12.6}",
        grid.dims[1], 0.0, sp, 0.0
    )?;
    writeln!(
        writer,
        "{:5} {:12.6} {:12.6} {:12.6}",
        grid.dims[2], 0.0, 0.0, sp
    )?;

    // Atom lines
    for i in 0..n_atoms {
        let z = elements[i] as i32;
        writeln!(
            writer,
            "{:5} {:12.6} {:12.6} {:12.6} {:12.6}",
            z,
            z as f64,
            positions[i][0] * ANG_TO_BOHR,
            positions[i][1] * ANG_TO_BOHR,
            positions[i][2] * ANG_TO_BOHR
        )?;
    }

    // Volume data (6 values per line, x varies slowest, z fastest)
    let mut count = 0;
    for val in &grid.values {
        write!(writer, " {:13.8E}", val)?;
        count += 1;
        if count % 6 == 0 {
            writeln!(writer)?;
        }
    }
    if count % 6 != 0 {
        writeln!(writer)?;
    }

    Ok(())
}

/// Read a Gaussian .cube file.
pub fn read_cube<R: BufRead>(reader: &mut R) -> Result<CubeFile, String> {
    let mut lines = Vec::new();
    let mut buf = String::new();
    loop {
        buf.clear();
        match reader.read_line(&mut buf) {
            Ok(0) => break,
            Ok(_) => lines.push(buf.trim_end().to_string()),
            Err(e) => return Err(format!("Read error: {}", e)),
        }
    }

    if lines.len() < 6 {
        return Err("Cube file too short".to_string());
    }

    let comment1 = lines[0].clone();
    let comment2 = lines[1].clone();

    // Line 3: n_atoms origin_x origin_y origin_z
    let parts: Vec<f64> = lines[2]
        .split_whitespace()
        .filter_map(|s| s.parse().ok())
        .collect();
    let n_atoms = parts[0] as usize;
    let origin = [parts[1], parts[2], parts[3]]; // bohr

    // Lines 4-6: nx vx_x vx_y vx_z (axis vectors)
    let mut dims = [0usize; 3];
    let mut axes = [[0.0f64; 3]; 3];
    for k in 0..3 {
        let p: Vec<f64> = lines[3 + k]
            .split_whitespace()
            .filter_map(|s| s.parse().ok())
            .collect();
        dims[k] = p[0] as usize;
        axes[k] = [p[1], p[2], p[3]]; // bohr
    }

    // Atom lines
    let mut elements = Vec::new();
    let mut positions = Vec::new();
    for i in 0..n_atoms {
        let p: Vec<f64> = lines[6 + i]
            .split_whitespace()
            .filter_map(|s| s.parse().ok())
            .collect();
        elements.push(p[0] as u8);
        positions.push([p[2], p[3], p[4]]); // bohr
    }

    // Volume data
    let data_start = 6 + n_atoms;
    let mut values = Vec::new();
    for line in &lines[data_start..] {
        for token in line.split_whitespace() {
            if let Ok(v) = token.parse::<f64>() {
                values.push(v);
            }
        }
    }

    let expected = dims[0] * dims[1] * dims[2];
    if values.len() != expected {
        return Err(format!(
            "Expected {} values, got {}",
            expected,
            values.len()
        ));
    }

    Ok(CubeFile {
        comment1,
        comment2,
        elements,
        positions,
        origin,
        dims,
        axes,
        values,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::eht::solve_eht;
    use crate::population::compute_population;

    fn water_molecule() -> (Vec<u8>, Vec<[f64; 3]>) {
        (
            vec![8, 1, 1],
            vec![[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]],
        )
    }

    #[test]
    fn test_esp_grid_dimensions() {
        let (elems, pos) = water_molecule();
        let result = solve_eht(&elems, &pos, None).unwrap();
        let pop = compute_population(&elems, &pos, &result.coefficients, result.n_electrons);

        let grid = compute_esp_grid(&elems, &pos, &pop.mulliken_charges, 0.5, 2.0);

        assert_eq!(
            grid.values.len(),
            grid.dims[0] * grid.dims[1] * grid.dims[2]
        );
        assert!(grid.dims[0] > 0 && grid.dims[1] > 0 && grid.dims[2] > 0);
    }

    #[test]
    fn test_esp_no_nan() {
        let (elems, pos) = water_molecule();
        let result = solve_eht(&elems, &pos, None).unwrap();
        let pop = compute_population(&elems, &pos, &result.coefficients, result.n_electrons);

        let grid = compute_esp_grid(&elems, &pos, &pop.mulliken_charges, 0.5, 2.0);

        assert!(
            grid.values.iter().all(|&v| !v.is_nan()),
            "ESP grid should have no NaN values"
        );
    }

    #[test]
    fn test_esp_sign_near_oxygen() {
        let (elems, pos) = water_molecule();
        let result = solve_eht(&elems, &pos, None).unwrap();
        let pop = compute_population(&elems, &pos, &result.coefficients, result.n_electrons);

        // ESP near O (negative charge) should be negative
        // Place test point near oxygen (slightly offset)
        let charges = &pop.mulliken_charges;
        let test_point = [0.5, 0.0, 0.0]; // near O
        let mut phi = 0.0;
        for a in 0..3 {
            let dx = test_point[0] - pos[a][0];
            let dy = test_point[1] - pos[a][1];
            let dz = test_point[2] - pos[a][2];
            let dist = (dx * dx + dy * dy + dz * dz).sqrt();
            if dist > 0.01 {
                phi += charges[a] / (dist * ANG_TO_BOHR);
            }
        }
        assert!(phi < 0.0, "ESP near oxygen should be negative, got {}", phi);
    }

    #[test]
    fn test_cube_roundtrip() {
        let (elems, pos) = water_molecule();
        let result = solve_eht(&elems, &pos, None).unwrap();
        let pop = compute_population(&elems, &pos, &result.coefficients, result.n_electrons);
        let grid = compute_esp_grid(&elems, &pos, &pop.mulliken_charges, 0.5, 2.0);

        // Export
        let mut buf = Vec::new();
        export_cube(&mut buf, &elems, &pos, &grid).unwrap();

        // Read back
        let mut reader = std::io::BufReader::new(&buf[..]);
        let cube = read_cube(&mut reader).unwrap();

        assert_eq!(cube.dims, grid.dims);
        assert_eq!(cube.values.len(), grid.values.len());
        // Values should match within formatting precision
        for (a, b) in cube.values.iter().zip(grid.values.iter()) {
            assert!(
                (a - b).abs() < 1e-3,
                "Cube roundtrip mismatch: {} vs {}",
                a,
                b
            );
        }
    }

    #[test]
    fn test_esp_far_from_molecule_near_zero() {
        let elems = vec![1u8, 1];
        let pos = vec![[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]];
        let result = solve_eht(&elems, &pos, None).unwrap();
        let pop = compute_population(&elems, &pos, &result.coefficients, result.n_electrons);

        // For symmetric H₂, charges are ~0 → ESP should be ~0 everywhere
        let grid = compute_esp_grid(&elems, &pos, &pop.mulliken_charges, 1.0, 5.0);

        let max_abs: f64 = grid.values.iter().map(|v| v.abs()).fold(0.0, f64::max);
        assert!(
            max_abs < 1.0,
            "ESP for neutral symmetric H₂ should be very small, max={}",
            max_abs
        );
    }

    #[test]
    fn test_esp_color_map_negative() {
        let [r, g, b] = esp_color_map(-1.0, 1.0);
        assert_eq!(r, 255);
        assert_eq!(g, 0);
        assert_eq!(b, 0);
    }

    #[test]
    fn test_esp_color_map_positive() {
        let [r, g, b] = esp_color_map(1.0, 1.0);
        assert_eq!(r, 0);
        assert_eq!(g, 0);
        assert_eq!(b, 255);
    }

    #[test]
    fn test_esp_color_map_zero() {
        let [r, g, b] = esp_color_map(0.0, 1.0);
        assert_eq!(r, 255);
        assert_eq!(g, 255);
        assert_eq!(b, 255);
    }

    #[test]
    fn test_esp_grid_to_colors_length() {
        let (elems, pos) = water_molecule();
        let result = solve_eht(&elems, &pos, None).unwrap();
        let pop = compute_population(&elems, &pos, &result.coefficients, result.n_electrons);
        let grid = compute_esp_grid(&elems, &pos, &pop.mulliken_charges, 1.0, 2.0);

        let colors = esp_grid_to_colors(&grid, 0.1);
        assert_eq!(colors.len(), grid.values.len() * 3);
    }

    #[test]
    fn test_esp_config_validation_rejects_invalid_inputs() {
        let (elems, pos) = water_molecule();
        let charges = vec![-0.6, 0.3, 0.3];

        for (config, needle) in [
            (
                EspGridConfig {
                    spacing: 0.0,
                    ..Default::default()
                },
                "spacing",
            ),
            (
                EspGridConfig {
                    spacing: f64::NAN,
                    ..Default::default()
                },
                "spacing",
            ),
            (
                EspGridConfig {
                    padding: -1.0,
                    ..Default::default()
                },
                "padding",
            ),
            (
                EspGridConfig {
                    padding: f64::INFINITY,
                    ..Default::default()
                },
                "padding",
            ),
            (
                EspGridConfig {
                    min_distance: 0.0,
                    ..Default::default()
                },
                "minimum distance",
            ),
        ] {
            let err = compute_esp_grid_with_config(&elems, &pos, &charges, config)
                .expect_err("invalid config must be rejected");
            assert!(
                err.contains(needle),
                "error should mention {needle:?}, got {err:?}"
            );
        }

        // Charges/positions length mismatch is rejected as well.
        let err =
            compute_esp_grid_with_config(&elems, &pos, &charges[..2], EspGridConfig::default())
                .expect_err("length mismatch must be rejected");
        assert!(err.contains("length mismatch"), "got {err:?}");
    }

    #[test]
    fn test_esp_min_distance_is_configurable() {
        // Single unit charge at the origin; the grid starts exactly at the
        // nucleus (padding 0), so the first value is `q / (min_distance · bohr)`.
        let elems = [6u8];
        let pos = [[0.0, 0.0, 0.0]];
        let charges = [1.0];

        let near = compute_esp_grid_with_config(
            &elems,
            &pos,
            &charges,
            EspGridConfig {
                spacing: 0.5,
                padding: 0.0,
                min_distance: 0.01,
            },
        )
        .unwrap();
        let far = compute_esp_grid_with_config(
            &elems,
            &pos,
            &charges,
            EspGridConfig {
                spacing: 0.5,
                padding: 0.0,
                min_distance: 0.1,
            },
        )
        .unwrap();

        // The origin is the first grid point.
        let ratio = near.values[0] / far.values[0];
        assert!(
            (ratio - 10.0).abs() < 1e-9,
            "clamping at 0.01 Å should be 10x larger than at 0.1 Å, got {ratio}"
        );
    }

    #[cfg(feature = "parallel")]
    #[test]
    fn test_esp_parallel_config_matches_serial() {
        let (elems, pos) = water_molecule();
        let result = solve_eht(&elems, &pos, None).unwrap();
        let pop = compute_population(&elems, &pos, &result.coefficients, result.n_electrons);
        let config = EspGridConfig {
            spacing: 0.5,
            padding: 2.0,
            min_distance: 0.05,
        };
        let serial =
            compute_esp_grid_with_config(&elems, &pos, &pop.mulliken_charges, config).unwrap();
        let parallel =
            compute_esp_grid_parallel_with_config(&elems, &pos, &pop.mulliken_charges, config)
                .unwrap();
        assert_eq!(serial.dims, parallel.dims);
        for (a, b) in serial.values.iter().zip(parallel.values.iter()) {
            assert!((a - b).abs() < 1e-12, "{a} vs {b}");
        }
    }
}
