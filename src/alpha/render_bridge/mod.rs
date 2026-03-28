//! Alpha render-bridge payload contracts.
//!
//! This module keeps visualization-facing payloads decoupled from solver
//! internals so Python, WASM, and charting frontends can share a stable format.
//! Supports EDL profiles, periodic DOS, band structures, capacitance scans,
//! temperature sweeps, microkinetic traces, and GSM path trajectories.
//!
//! Transport decisions (Sprint 1):
//! - Small results: JSON via `serde_json` (< 64 KiB).
//! - Large profile/DOS arrays: Arrow-style `RecordBatch` via `crate::transport::arrow`.
//! - WASM typed arrays: `Float64Array` / `Float32Array` via flat contiguous buffers.

#[cfg(feature = "experimental-gpu")]
pub mod gpu;

use serde::{Deserialize, Serialize};

use crate::transport::arrow::RecordBatch;

#[cfg(feature = "alpha-edl")]
use crate::alpha::edl::{CpmEdlScanResult, EdlProfileResult};
#[cfg(feature = "alpha-kinetics")]
use crate::alpha::kinetics::{ElementaryRateResult, MicrokineticTrace};
#[cfg(feature = "alpha-periodic-linear")]
use crate::alpha::periodic_linear::{BandStructureAdapterResult, PeriodicKpmDosResult};

/// One XY series payload for charts.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Default)]
pub struct ChartSeries {
    /// Stable series identifier.
    pub series_id: String,
    /// Human-readable label.
    pub label: String,
    /// X values.
    pub x: Vec<f64>,
    /// Y values.
    pub y: Vec<f64>,
    /// Suggested unit label for the x axis.
    pub x_unit: String,
    /// Suggested unit label for the y axis.
    pub y_unit: String,
}

/// Multi-series chart payload.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Default)]
pub struct ChartPayload {
    /// Chart title.
    pub title: String,
    /// Included series.
    pub series: Vec<ChartSeries>,
}

/// One 3D polyline frame, typically for reaction trajectories.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Default)]
pub struct TrajectoryFrame {
    /// Time-like coordinate.
    pub frame_time: f64,
    /// Flattened xyz coordinates.
    pub coords: Vec<f64>,
}

/// Generic bridge payload for render clients.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Default)]
pub struct RenderPayload {
    /// Scalar chart data.
    pub charts: Vec<ChartPayload>,
    /// Optional geometry frames.
    pub trajectory: Vec<TrajectoryFrame>,
}

#[cfg(feature = "alpha-edl")]
/// Convert an EDL profile into a standard chart payload.
pub fn edl_profile_chart(profile: &EdlProfileResult) -> ChartPayload {
    ChartPayload {
        title: format!("EDL profile ({})", profile.model_name),
        series: vec![
            ChartSeries {
                series_id: "potential".into(),
                label: "Potential".into(),
                x: profile.distance_axis_angstrom.clone(),
                y: profile.electrostatic_potential_v.clone(),
                x_unit: "A".into(),
                y_unit: "V".into(),
            },
            ChartSeries {
                series_id: "charge_density".into(),
                label: "Charge density".into(),
                x: profile.distance_axis_angstrom.clone(),
                y: profile.charge_density_c_per_m3.clone(),
                x_unit: "A".into(),
                y_unit: "C/m^3".into(),
            },
        ],
    }
}

#[cfg(feature = "alpha-periodic-linear")]
/// Convert a periodic DOS result into a standard chart payload.
pub fn periodic_dos_chart(result: &PeriodicKpmDosResult) -> ChartPayload {
    ChartPayload {
        title: "Periodic DOS".into(),
        series: vec![ChartSeries {
            series_id: "dos".into(),
            label: "Density of states".into(),
            x: result.energies_ev.clone(),
            y: result.total_dos.clone(),
            x_unit: "eV".into(),
            y_unit: "a.u.".into(),
        }],
    }
}

#[cfg(feature = "alpha-kinetics")]
/// Convert a microkinetic trace into one series per species.
pub fn microkinetic_population_chart(trace: &MicrokineticTrace) -> ChartPayload {
    let Some(first_frame) = trace.frames.first() else {
        return ChartPayload {
            title: "Microkinetic populations".into(),
            series: Vec::new(),
        };
    };

    let mut series = Vec::with_capacity(first_frame.species.len());
    for species in &first_frame.species {
        let mut x = Vec::with_capacity(trace.frames.len());
        let mut y = Vec::with_capacity(trace.frames.len());
        for frame in &trace.frames {
            if let Some(population) = frame
                .species
                .iter()
                .find(|population| population.species_id == species.species_id)
            {
                x.push(frame.time_s);
                y.push(population.population);
            }
        }
        series.push(ChartSeries {
            series_id: species.species_id.clone(),
            label: species.species_id.clone(),
            x,
            y,
            x_unit: "s".into(),
            y_unit: "arb.".into(),
        });
    }

    ChartPayload {
        title: "Microkinetic populations".into(),
        series,
    }
}

#[cfg(feature = "alpha-edl")]
/// Convert a CPM capacitance scan into a chart payload.
pub fn capacitance_scan_chart(scan: &CpmEdlScanResult) -> ChartPayload {
    ChartPayload {
        title: "CPM Capacitance Scan".into(),
        series: vec![
            ChartSeries {
                series_id: "charge_vs_potential".into(),
                label: "Q(μ)".into(),
                x: scan.mu_values_ev.clone(),
                y: scan.total_charge_e.clone(),
                x_unit: "eV".into(),
                y_unit: "e".into(),
            },
            ChartSeries {
                series_id: "capacitance".into(),
                label: "C(μ)".into(),
                x: scan.mu_values_ev.clone(),
                y: scan.capacitance_e_per_ev.clone(),
                x_unit: "eV".into(),
                y_unit: "e/eV".into(),
            },
        ],
    }
}

#[cfg(feature = "alpha-kinetics")]
/// Convert HTST temperature sweep results into an Arrhenius chart.
pub fn arrhenius_chart(results: &[ElementaryRateResult]) -> ChartPayload {
    let x: Vec<f64> = results
        .iter()
        .map(|r| 1000.0 / r.state.temperature_k) // 1000/T for Arrhenius
        .collect();
    let y: Vec<f64> = results.iter().map(|r| r.forward_rate_s_inv.ln()).collect();

    ChartPayload {
        title: "Arrhenius Plot".into(),
        series: vec![ChartSeries {
            series_id: "arrhenius".into(),
            label: "ln(k) vs 1000/T".into(),
            x,
            y,
            x_unit: "1000/K".into(),
            y_unit: "ln(s⁻¹)".into(),
        }],
    }
}

#[cfg(feature = "alpha-periodic-linear")]
/// Convert a band-structure result into a multi-series chart (one series per band).
pub fn band_structure_chart(bs: &BandStructureAdapterResult) -> ChartPayload {
    let mut series = Vec::with_capacity(bs.n_bands.min(20)); // limit to 20 bands for display
    let x: Vec<f64> = (0..bs.n_kpoints).map(|i| i as f64).collect();
    for band_idx in 0..bs.n_bands.min(20) {
        let y: Vec<f64> = bs
            .bands
            .iter()
            .map(|bands_at_k| {
                if band_idx < bands_at_k.len() {
                    bands_at_k[band_idx]
                } else {
                    0.0
                }
            })
            .collect();
        series.push(ChartSeries {
            series_id: format!("band_{}", band_idx),
            label: format!("Band {}", band_idx + 1),
            x: x.clone(),
            y,
            x_unit: "k-index".into(),
            y_unit: "eV".into(),
        });
    }

    ChartPayload {
        title: "Band Structure".into(),
        series,
    }
}

#[cfg(feature = "alpha-kinetics")]
/// Convert GSM path coordinates into trajectory frames for 3D visualization.
pub fn gsm_path_trajectory(
    path_coords: &[Vec<f64>],
    path_energies: &[f64],
) -> Vec<TrajectoryFrame> {
    path_coords
        .iter()
        .zip(path_energies)
        .map(|(coords, &energy)| TrajectoryFrame {
            frame_time: energy, // Use energy as the "time" coordinate for path visualization
            coords: coords.clone(),
        })
        .collect()
}

/// Build a full render payload from all available alpha results.
#[cfg(all(
    feature = "alpha-edl",
    feature = "alpha-periodic-linear",
    feature = "alpha-kinetics"
))]
pub fn build_full_render_payload(
    edl_profile: Option<&EdlProfileResult>,
    periodic_dos: Option<&PeriodicKpmDosResult>,
    band_structure: Option<&BandStructureAdapterResult>,
    microkinetic_trace: Option<&MicrokineticTrace>,
    gsm_path: Option<(&[Vec<f64>], &[f64])>,
) -> RenderPayload {
    let mut charts = Vec::new();
    let mut trajectory = Vec::new();

    if let Some(profile) = edl_profile {
        charts.push(edl_profile_chart(profile));
    }
    if let Some(dos) = periodic_dos {
        charts.push(periodic_dos_chart(dos));
    }
    if let Some(bs) = band_structure {
        charts.push(band_structure_chart(bs));
    }
    if let Some(trace) = microkinetic_trace {
        charts.push(microkinetic_population_chart(trace));
    }
    if let Some((coords, energies)) = gsm_path {
        trajectory = gsm_path_trajectory(coords, energies);
    }

    RenderPayload { charts, trajectory }
}

// ── Transport pack helpers ──────────────────────────────────────────────────

#[cfg(feature = "alpha-edl")]
/// Pack an EDL profile into an Arrow-style `RecordBatch` for typed-array transport.
pub fn pack_edl_profile(profile: &EdlProfileResult) -> RecordBatch {
    let n = profile.distance_axis_angstrom.len();
    let mut batch = RecordBatch::new();
    batch.add_float64(
        "distance_angstrom",
        profile.distance_axis_angstrom.clone(),
        vec![n],
    );
    batch.add_float64(
        "potential_v",
        profile.electrostatic_potential_v.clone(),
        vec![n],
    );
    batch.add_float64(
        "charge_density_c_per_m3",
        profile.charge_density_c_per_m3.clone(),
        vec![n],
    );
    batch
}

#[cfg(feature = "alpha-periodic-linear")]
/// Pack periodic KPM DOS into a `RecordBatch`.
pub fn pack_periodic_dos(result: &PeriodicKpmDosResult) -> RecordBatch {
    let n = result.energies_ev.len();
    let mut batch = RecordBatch::new();
    batch.add_float64("energies_ev", result.energies_ev.clone(), vec![n]);
    batch.add_float64("total_dos", result.total_dos.clone(), vec![n]);
    batch
}

#[cfg(feature = "alpha-kinetics")]
/// Pack an Arrhenius temperature sweep into a `RecordBatch`.
pub fn pack_arrhenius(results: &[ElementaryRateResult]) -> RecordBatch {
    let n = results.len();
    let mut batch = RecordBatch::new();
    batch.add_float64(
        "temperature_k",
        results.iter().map(|r| r.state.temperature_k).collect(),
        vec![n],
    );
    batch.add_float64(
        "forward_rate_s_inv",
        results.iter().map(|r| r.forward_rate_s_inv).collect(),
        vec![n],
    );
    batch.add_float64(
        "equilibrium_constant",
        results.iter().map(|r| r.equilibrium_constant).collect(),
        vec![n],
    );
    batch
}

#[cfg(feature = "alpha-kinetics")]
/// Pack a microkinetic trace into a `RecordBatch` (time column + one per species).
pub fn pack_microkinetic_trace(trace: &MicrokineticTrace) -> RecordBatch {
    let n = trace.frames.len();
    let mut batch = RecordBatch::new();
    batch.add_float64(
        "time_s",
        trace.frames.iter().map(|f| f.time_s).collect(),
        vec![n],
    );
    if let Some(first) = trace.frames.first() {
        for sp in &first.species {
            let col: Vec<f64> = trace
                .frames
                .iter()
                .map(|f| {
                    f.species
                        .iter()
                        .find(|s| s.species_id == sp.species_id)
                        .map_or(0.0, |s| s.population)
                })
                .collect();
            batch.add_float64(&sp.species_id, col, vec![n]);
        }
    }
    batch
}

/// Pack a generic `ChartPayload` into a `RecordBatch`.
pub fn pack_chart_payload(chart: &ChartPayload) -> RecordBatch {
    let mut batch = RecordBatch::new();
    for (i, s) in chart.series.iter().enumerate() {
        let n = s.x.len();
        batch.add_float64(&format!("{}_x", s.series_id), s.x.clone(), vec![n]);
        batch.add_float64(&format!("{}_y", s.series_id), s.y.clone(), vec![n]);
        let _ = i; // series index for diagnostics
    }
    batch
}

/// Pack multiple chart payloads into `RecordBatch` in parallel using rayon.
#[cfg(feature = "parallel")]
pub fn pack_charts_parallel(charts: &[ChartPayload]) -> Vec<RecordBatch> {
    use rayon::prelude::*;
    charts.par_iter().map(pack_chart_payload).collect()
}

/// Pack multiple chart payloads into `RecordBatch` sequentially (non-parallel fallback).
#[cfg(not(feature = "parallel"))]
pub fn pack_charts_parallel(charts: &[ChartPayload]) -> Vec<RecordBatch> {
    charts.iter().map(pack_chart_payload).collect()
}

/// Verify lossless round-trip: `ChartPayload` → `RecordBatch` → reconstructed series.
/// Returns `true` when all values survive within `tol`.
pub fn verify_chart_round_trip(chart: &ChartPayload, tol: f64) -> bool {
    let batch = pack_chart_payload(chart);
    for s in &chart.series {
        let x_col = batch
            .float_columns
            .iter()
            .find(|c| c.name == format!("{}_x", s.series_id));
        let y_col = batch
            .float_columns
            .iter()
            .find(|c| c.name == format!("{}_y", s.series_id));
        match (x_col, y_col) {
            (Some(xc), Some(yc)) => {
                if xc.values.len() != s.x.len() || yc.values.len() != s.y.len() {
                    return false;
                }
                for (a, b) in xc.values.iter().zip(&s.x) {
                    if (a - b).abs() > tol {
                        return false;
                    }
                }
                for (a, b) in yc.values.iter().zip(&s.y) {
                    if (a - b).abs() > tol {
                        return false;
                    }
                }
            }
            _ => return false,
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;
    #[cfg(feature = "alpha-kinetics")]
    use crate::alpha::kinetics::{MicrokineticFrame, SpeciesPopulation};

    #[cfg(feature = "alpha-kinetics")]
    #[test]
    fn microkinetic_chart_emits_one_series_per_species() {
        let chart = microkinetic_population_chart(&MicrokineticTrace {
            rates: Vec::new(),
            frames: vec![
                MicrokineticFrame {
                    time_s: 0.0,
                    species: vec![
                        SpeciesPopulation {
                            species_id: "A*".into(),
                            population: 0.8,
                        },
                        SpeciesPopulation {
                            species_id: "B*".into(),
                            population: 0.2,
                        },
                    ],
                },
                MicrokineticFrame {
                    time_s: 1.0,
                    species: vec![
                        SpeciesPopulation {
                            species_id: "A*".into(),
                            population: 0.6,
                        },
                        SpeciesPopulation {
                            species_id: "B*".into(),
                            population: 0.4,
                        },
                    ],
                },
            ],
        });

        assert_eq!(chart.series.len(), 2);
        assert_eq!(chart.series[0].x.len(), 2);
    }

    #[cfg(feature = "alpha-kinetics")]
    #[test]
    fn arrhenius_chart_produces_linear_plot() {
        use crate::alpha::kinetics::{ElementaryRateResult, ThermodynamicState};

        let results: Vec<ElementaryRateResult> = (300..=600)
            .step_by(100)
            .map(|t| {
                let t = t as f64;
                ElementaryRateResult {
                    step_id: "test".into(),
                    forward_rate_s_inv: 1e13 * (-0.5 / (8.617e-5 * t)).exp(),
                    reverse_rate_s_inv: 1.0,
                    equilibrium_constant: 1.0,
                    state: ThermodynamicState {
                        temperature_k: t,
                        pressure_bar: 1.0,
                    },
                }
            })
            .collect();

        let chart = arrhenius_chart(&results);
        assert_eq!(chart.series.len(), 1);
        assert_eq!(chart.series[0].x.len(), 4);
        // Check linearity: ln(k) vs 1/T should be roughly linear
        let x = &chart.series[0].x;
        let y = &chart.series[0].y;
        // Slope should be approximately constant
        if x.len() >= 3 {
            let slope1 = (y[1] - y[0]) / (x[1] - x[0]);
            let slope2 = (y[2] - y[1]) / (x[2] - x[1]);
            assert!(
                (slope1 - slope2).abs() / slope1.abs() < 0.1,
                "Arrhenius plot should be linear"
            );
        }
    }

    #[cfg(feature = "alpha-kinetics")]
    #[test]
    fn gsm_path_trajectory_correct_frame_count() {
        let coords = vec![
            vec![0.0, 0.0, 0.0],
            vec![1.0, 0.0, 0.0],
            vec![2.0, 0.0, 0.0],
        ];
        let energies = vec![0.0, 1.5, 0.5];
        let frames = gsm_path_trajectory(&coords, &energies);
        assert_eq!(frames.len(), 3);
        assert!((frames[1].frame_time - 1.5).abs() < 1e-10);
    }

    // ── Transport round-trip tests ──────────────────────────────────────────

    #[test]
    fn chart_payload_round_trip_lossless() {
        let chart = ChartPayload {
            title: "test".into(),
            series: vec![
                ChartSeries {
                    series_id: "s1".into(),
                    label: "Series 1".into(),
                    x: vec![0.0, 1.0, 2.0],
                    y: vec![3.0, 4.0, 5.0],
                    x_unit: "eV".into(),
                    y_unit: "a.u.".into(),
                },
                ChartSeries {
                    series_id: "s2".into(),
                    label: "Series 2".into(),
                    x: vec![10.0, 20.0],
                    y: vec![0.5, 1.5],
                    x_unit: "Å".into(),
                    y_unit: "V".into(),
                },
            ],
        };
        assert!(verify_chart_round_trip(&chart, 0.0));
    }

    #[test]
    fn pack_chart_payload_preserves_column_count() {
        let chart = ChartPayload {
            title: "multi".into(),
            series: vec![
                ChartSeries {
                    series_id: "a".into(),
                    label: "A".into(),
                    x: vec![1.0],
                    y: vec![2.0],
                    x_unit: "".into(),
                    y_unit: "".into(),
                },
                ChartSeries {
                    series_id: "b".into(),
                    label: "B".into(),
                    x: vec![3.0],
                    y: vec![4.0],
                    x_unit: "".into(),
                    y_unit: "".into(),
                },
            ],
        };
        let batch = pack_chart_payload(&chart);
        // 2 series × 2 columns (x + y) = 4
        assert_eq!(batch.float_columns.len(), 4);
    }

    #[cfg(feature = "alpha-edl")]
    #[test]
    fn edl_profile_packs_to_record_batch() {
        use crate::alpha::edl::{CapacitanceBreakdown, EdlBackend, EdlProfileResult};
        let profile = EdlProfileResult {
            distance_axis_angstrom: vec![0.0, 1.0, 2.0],
            electrostatic_potential_v: vec![0.5, 0.3, 0.1],
            field_strength_v_per_m: vec![0.0; 3],
            charge_density_c_per_m3: vec![1.0, 0.5, 0.1],
            ion_density_relative: vec![0.0; 3],
            dielectric_profile: vec![78.5; 3],
            compact_layer_drop_v: 0.3,
            diffuse_layer_drop_v: 0.2,
            total_interfacial_drop_v: 0.5,
            differential_capacitance: CapacitanceBreakdown::default(),
            model_name: "test".into(),
            backend: EdlBackend::CpuReference,
            used_gpu: false,
            converged: true,
            n_iterations: 1,
            residual: 0.0,
            temperature_k: 298.15,
            ionic_strength_m: 1.0,
        };
        let batch = pack_edl_profile(&profile);
        assert_eq!(batch.float_columns.len(), 3);
        assert_eq!(batch.float_columns[0].values.len(), 3);
    }

    #[cfg(feature = "alpha-kinetics")]
    #[test]
    fn arrhenius_packs_to_record_batch() {
        use crate::alpha::kinetics::{ElementaryRateResult, ThermodynamicState};

        let results: Vec<ElementaryRateResult> = vec![
            ElementaryRateResult {
                step_id: "test".into(),
                forward_rate_s_inv: 1.0e10,
                reverse_rate_s_inv: 1.0,
                equilibrium_constant: 1.0e10,
                state: ThermodynamicState {
                    temperature_k: 300.0,
                    pressure_bar: 1.0,
                },
            },
            ElementaryRateResult {
                step_id: "test".into(),
                forward_rate_s_inv: 1.0e12,
                reverse_rate_s_inv: 1.0,
                equilibrium_constant: 1.0e12,
                state: ThermodynamicState {
                    temperature_k: 600.0,
                    pressure_bar: 1.0,
                },
            },
        ];
        let batch = pack_arrhenius(&results);
        assert_eq!(batch.float_columns.len(), 3);
        assert_eq!(batch.float_columns[0].values.len(), 2);
    }

    #[cfg(feature = "alpha-kinetics")]
    #[test]
    fn microkinetic_trace_packs_to_record_batch() {
        let trace = MicrokineticTrace {
            rates: Vec::new(),
            frames: vec![
                MicrokineticFrame {
                    time_s: 0.0,
                    species: vec![SpeciesPopulation {
                        species_id: "A".into(),
                        population: 1.0,
                    }],
                },
                MicrokineticFrame {
                    time_s: 1.0,
                    species: vec![SpeciesPopulation {
                        species_id: "A".into(),
                        population: 0.5,
                    }],
                },
            ],
        };
        let batch = pack_microkinetic_trace(&trace);
        // time_s + A = 2 columns
        assert_eq!(batch.float_columns.len(), 2);
    }

    // ── JSON round-trip for all payload types ───────────────────────────────

    #[test]
    fn render_payload_json_round_trip() {
        let payload = RenderPayload {
            charts: vec![ChartPayload {
                title: "test".into(),
                series: vec![ChartSeries {
                    series_id: "s1".into(),
                    label: "S1".into(),
                    x: vec![1.0, 2.0, 3.0],
                    y: vec![4.0, 5.0, 6.0],
                    x_unit: "eV".into(),
                    y_unit: "a.u.".into(),
                }],
            }],
            trajectory: vec![TrajectoryFrame {
                frame_time: 0.5,
                coords: vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
            }],
        };
        let json = serde_json::to_string(&payload).unwrap();
        let recovered: RenderPayload = serde_json::from_str(&json).unwrap();
        assert_eq!(payload, recovered);
    }

    #[test]
    fn pack_charts_parallel_matches_sequential() {
        let charts: Vec<ChartPayload> = (0..8)
            .map(|i| ChartPayload {
                title: format!("chart_{}", i),
                series: vec![ChartSeries {
                    series_id: format!("s{}", i),
                    label: format!("S{}", i),
                    x: vec![i as f64; 10],
                    y: vec![i as f64 * 2.0; 10],
                    x_unit: "".into(),
                    y_unit: "".into(),
                }],
            })
            .collect();
        let parallel_result = pack_charts_parallel(&charts);
        let sequential_result: Vec<RecordBatch> = charts.iter().map(pack_chart_payload).collect();
        assert_eq!(parallel_result.len(), sequential_result.len());
        for (p, s) in parallel_result.iter().zip(&sequential_result) {
            assert_eq!(p.float_columns.len(), s.float_columns.len());
        }
    }
}
