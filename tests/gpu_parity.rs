//! Real GPU-vs-CPU parity gate (native wgpu backend).
//!
//! These tests exercise the actual `wgpu` device. On a machine without a
//! compatible adapter they are skipped with a message; set
//! `SCIFORM_REQUIRE_GPU=1` to turn a missing adapter into a failure (useful
//! for a dedicated GPU runner).
//!
//! The thresholds are the ones used by the candidate benchmarks in
//! `tests/regression/test_gpu_candidates.rs`, promoted here to hard
//! assertions so a kernel regression cannot hide behind the benchmark
//! harness that turns panics into `SKIPPED` lines.
//!
//! Systems are deliberately small: the goal is numerical parity, not a
//! throughput comparison.

#![cfg(feature = "experimental-gpu")]

use sci_form::gpu::context::GpuContext;

fn gpu_context(label: &str) -> Option<GpuContext> {
    let ctx = GpuContext::best_available();
    if ctx.is_gpu_available() {
        eprintln!("{label}: {}", ctx.capabilities.backend);
        return Some(ctx);
    }

    let message = format!(
        "GPU unavailable for {label}: {}",
        ctx.activation_report().reason
    );
    if std::env::var_os("SCIFORM_REQUIRE_GPU").is_some() {
        panic!("{message}");
    }
    eprintln!("{message}; skipping");
    None
}

fn water() -> (Vec<u8>, Vec<[f64; 3]>) {
    (
        vec![8, 1, 1],
        vec![
            [0.000, 0.000, 0.117],
            [0.000, 0.757, -0.469],
            [0.000, -0.757, -0.469],
        ],
    )
}

fn benzene() -> (Vec<u8>, Vec<[f64; 3]>) {
    let conformer = sci_form::embed("c1ccccc1", 42);
    assert!(conformer.error.is_none(), "embed failed");
    let positions = conformer
        .coords
        .chunks_exact(3)
        .map(|p| [p[0], p[1], p[2]])
        .collect();
    (conformer.elements, positions)
}

fn to_bohr(positions: &[[f64; 3]]) -> Vec<[f64; 3]> {
    const ANGSTROM_TO_BOHR: f64 = 1.889_725_988_6;
    positions
        .iter()
        .map(|p| {
            [
                p[0] * ANGSTROM_TO_BOHR,
                p[1] * ANGSTROM_TO_BOHR,
                p[2] * ANGSTROM_TO_BOHR,
            ]
        })
        .collect()
}

fn max_abs_diff(lhs: &[f64], rhs: &[f64]) -> f64 {
    assert_eq!(lhs.len(), rhs.len(), "length mismatch");
    lhs.iter()
        .zip(rhs)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0, f64::max)
}

#[test]
fn vector_add_matches_cpu() {
    let Some(ctx) = gpu_context("vector add") else {
        return;
    };
    let lhs: Vec<f32> = (0..10_000).map(|i| (i as f32) * 0.25).collect();
    let rhs: Vec<f32> = (0..10_000).map(|i| 1.0 - (i as f32) * 0.125).collect();
    let gpu = ctx.vector_add_f32(&lhs, &rhs).expect("GPU vector add");
    let cpu: Vec<f32> = lhs.iter().zip(&rhs).map(|(a, b)| a + b).collect();
    let max_diff = cpu
        .iter()
        .zip(&gpu)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0f32, f32::max);
    assert!(max_diff <= 1e-5, "vector add mismatch: {max_diff:e}");
}

#[test]
fn eri_matches_cpu() {
    let Some(ctx) = gpu_context("ERI") else {
        return;
    };
    // The GPU kernel is s-only; use a 20-atom hydrogen chain so every basis
    // function is an s function. p/d bases are covered by the rejection test
    // in src/gpu/two_electron_gpu.rs.
    let elements = vec![1u8; 20];
    let positions_bohr: Vec<[f64; 3]> = (0..20).map(|i| [i as f64 * 1.7, 0.0, 0.0]).collect();
    let basis = sci_form::scf::basis::BasisSet::sto3g(&elements, &positions_bohr);

    let cpu = sci_form::scf::two_electron::TwoElectronIntegrals::compute_parallel(&basis);
    let gpu =
        sci_form::gpu::two_electron_gpu::compute_eris_gpu(&ctx, &basis).expect("GPU ERI dispatch");
    let n = basis.n_basis;
    let mut worst = 0.0f64;
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                for l in 0..n {
                    worst = worst.max((cpu.get(i, j, k, l) - gpu.get(i, j, k, l)).abs());
                }
            }
        }
    }
    eprintln!("ERI H20 max |Δ| = {worst:e} a.u.");
    assert!(worst < 5e-3, "ERI mismatch: {worst:e}");
}

#[test]
fn fock_matches_cpu() {
    let Some(ctx) = gpu_context("Fock build") else {
        return;
    };
    let (elements, positions) = water();
    let system = sci_form::scf::types::MolecularSystem::from_angstrom(&elements, &positions, 0, 1);
    let basis =
        sci_form::scf::basis::BasisSet::sto3g(&system.atomic_numbers, &system.positions_bohr);
    let eris = sci_form::scf::two_electron::TwoElectronIntegrals::compute_parallel(&basis);
    let n = basis.n_basis;

    // Deterministic symmetric H_core and density trial matrices.
    let h_core: Vec<f64> = (0..n * n)
        .map(|idx| {
            let i = idx / n;
            let j = idx % n;
            if i == j {
                -1.0 - 0.1 * i as f64
            } else {
                0.02 / (1.0 + (i as f64 - j as f64).abs())
            }
        })
        .collect();
    let density_seed: Vec<f64> = (0..n * n)
        .map(|idx| {
            let i = idx / n;
            let j = idx % n;
            if i == j {
                0.5
            } else {
                0.01 * ((i + 2 * j) as f64).sin()
            }
        })
        .collect();
    let mut density = density_seed.clone();
    for i in 0..n {
        for j in 0..i {
            let value = 0.5 * (density_seed[i * n + j] + density_seed[j * n + i]);
            density[i * n + j] = value;
            density[j * n + i] = value;
        }
    }

    let eri_flat: Vec<f64> = (0..n * n * n * n)
        .map(|idx| {
            let l = idx % n;
            let k = (idx / n) % n;
            let j = (idx / (n * n)) % n;
            let i = idx / (n * n * n);
            eris.get(i, j, k, l)
        })
        .collect();
    let gpu = sci_form::gpu::fock_build_gpu::build_fock_gpu(&ctx, &h_core, &density, &eri_flat, n)
        .expect("GPU Fock dispatch");

    let mut cpu = h_core.clone();
    for i in 0..n {
        for j in 0..n {
            let mut g = 0.0;
            for k in 0..n {
                for l in 0..n {
                    let p = density[k * n + l];
                    g += p * (eris.get(i, j, k, l) - 0.5 * eris.get(i, k, j, l));
                }
            }
            cpu[i * n + j] += g;
        }
    }

    let worst = max_abs_diff(&cpu, &gpu);
    eprintln!("Fock H2O max |Δ| = {worst:e}");
    assert!(worst < 5e-3, "Fock mismatch: {worst:e}");
}

#[test]
fn gamma_matrix_matches_cpu() {
    let Some(ctx) = gpu_context("gamma matrix") else {
        return;
    };
    let (elements, positions) = water();
    let pos_bohr = to_bohr(&positions);
    let eta: Vec<f64> = elements
        .iter()
        .map(|&z| match z {
            1 => 0.4720,
            8 => 0.4991,
            _ => 0.3,
        })
        .collect();

    let cpu = sci_form::scf::fock_matrix::build_gamma_matrix(&eta, &pos_bohr);
    let gpu = sci_form::gpu::gamma_matrix_gpu::build_gamma_gpu(&ctx, &eta, &pos_bohr)
        .expect("GPU gamma dispatch");
    let n = cpu.nrows();
    let mut worst = 0.0f64;
    for i in 0..n {
        for j in 0..n {
            worst = worst.max((cpu[(i, j)] - gpu[(i, j)]).abs());
        }
    }
    eprintln!("Gamma H2O max |Δ| = {worst:e}");
    assert!(worst < 1e-5, "gamma matrix mismatch: {worst:e}");
}

#[test]
fn orbital_and_density_grids_match_cpu() {
    let Some(_ctx) = gpu_context("orbital/density grids") else {
        return;
    };
    let (elements, positions) = water();
    let system = sci_form::scf::types::MolecularSystem::from_angstrom(&elements, &positions, 0, 1);
    let config = sci_form::scf::scf_loop::ScfConfig {
        max_iterations: 20,
        energy_threshold: 1e-5,
        density_threshold: 1e-4,
        diis_size: 4,
        level_shift: 0.0,
        damping: 0.15,
        use_parallel_eri: true,
        parallel_threshold: 0,
    };
    let scf = sci_form::scf::scf_loop::run_scf(&system, &config);
    let basis =
        sci_form::scf::basis::BasisSet::sto3g(&system.atomic_numbers, &system.positions_bohr);
    let params =
        sci_form::gpu::orbital_grid::GridParams::from_molecule(&system.positions_bohr, 0.5, 3.0);

    let homo = scf.n_electrons / 2 - 1;
    let cpu_orbital = sci_form::gpu::orbital_grid::evaluate_orbital_cpu(
        &basis,
        &scf.mo_coefficients,
        homo,
        &params,
    );
    let (gpu_orbital, report) = sci_form::gpu::orbital_grid::evaluate_orbital_with_report(
        &basis,
        &scf.mo_coefficients,
        homo,
        &params,
    );
    assert!(
        report.used_gpu,
        "orbital grid did not use GPU: {}",
        report.note
    );
    let orbital_diff = max_abs_diff(&cpu_orbital.values, &gpu_orbital.values);
    eprintln!("Orbital grid H2O max |Δ| = {orbital_diff:e}");
    assert!(
        orbital_diff < 5e-4,
        "orbital grid mismatch: {orbital_diff:e}"
    );

    let density = sci_form::scf::density_matrix::build_density_matrix(
        &scf.mo_coefficients,
        scf.n_electrons / 2,
    );
    let cpu_density = sci_form::gpu::orbital_grid::evaluate_density_cpu(&basis, &density, &params);
    let (gpu_density, density_report) =
        sci_form::gpu::density_grid_gpu::evaluate_density_with_report(&basis, &density, &params);
    assert!(
        density_report.used_gpu,
        "density grid did not use GPU: {}",
        density_report.note
    );
    let density_diff = max_abs_diff(&cpu_density, &gpu_density);
    eprintln!("Density grid H2O max |Δ| = {density_diff:e}");
    assert!(
        density_diff < 1e-3,
        "density grid mismatch: {density_diff:e}"
    );
}

#[test]
fn esp_grid_matches_cpu() {
    let Some(_ctx) = gpu_context("ESP grid") else {
        return;
    };
    let (elements, positions) = water();
    let charges = [-0.8, 0.4, 0.4];
    let config = sci_form::esp::EspGridConfig {
        spacing: 0.7,
        padding: 2.0,
        ..Default::default()
    };

    let cpu = sci_form::esp::compute_esp_grid_with_config(&elements, &positions, &charges, config)
        .expect("CPU ESP grid");
    let (gpu, report) = sci_form::gpu::esp_grid_gpu::compute_esp_grid_with_report(
        &elements, &positions, &charges, config,
    );
    assert!(report.used_gpu, "ESP grid did not use GPU: {}", report.note);
    assert_eq!(cpu.values.len(), gpu.values.len());
    let worst = max_abs_diff(&cpu.values, &gpu.values);
    eprintln!("ESP grid H2O max |Δ| = {worst:e} Ha/e");
    assert!(worst < 1e-4, "ESP grid mismatch: {worst:e}");
}

#[test]
fn eeq_charges_match_cpu() {
    let Some(ctx) = gpu_context("EEQ charges") else {
        return;
    };
    let (elements, positions) = benzene();
    let config = sci_form::charges_eeq::EeqConfig::default();
    let cpu = sci_form::charges_eeq::compute_eeq_charges(&elements, &positions, &config);
    let gpu = sci_form::gpu::eeq_gpu::compute_eeq_charges_gpu(&ctx, &elements, &positions, &config)
        .expect("GPU EEQ dispatch");
    let charge_diff = max_abs_diff(&cpu.charges, &gpu.charges);
    let total_diff = (cpu.total_charge - gpu.total_charge).abs();
    eprintln!("EEQ benzene max |Δq| = {charge_diff:e}, |ΔQ| = {total_diff:e}");
    assert!(charge_diff < 1e-3, "EEQ charge mismatch: {charge_diff:e}");
    assert!(
        total_diff < 1e-6,
        "EEQ total charge mismatch: {total_diff:e}"
    );
}

#[test]
fn cpm_charges_match_cpu() {
    let Some(ctx) = gpu_context("CPM charges") else {
        return;
    };
    // Two small electrodes with the same potential; the 2-site closed form is
    // not used here, the CPU path is the reference.
    let positions = vec![
        [0.0, 0.0, 0.0],
        [1.4, 0.0, 0.0],
        [0.0, 0.0, 6.0],
        [1.4, 0.0, 6.0],
    ];
    let elements = vec![6, 6, 6, 6];
    let config = sci_form::cpm::CpmConfig::default();
    let cpu = sci_form::cpm::compute_cpm_charges(&elements, &positions, &config);
    let gpu = sci_form::gpu::cpm_gpu::compute_cpm_charges_gpu(&ctx, &elements, &positions, &config)
        .expect("GPU CPM dispatch");
    let charge_diff = max_abs_diff(&cpu.charges, &gpu.charges);
    eprintln!("CPM 2-site max |Δq| = {charge_diff:e}");
    assert!(charge_diff < 1e-5, "CPM charge mismatch: {charge_diff:e}");
}
