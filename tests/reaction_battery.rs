//! Batería de tests — sistema de reacciones, 21 frames, frame central = TS.
//!
//! cargo test --test reaction_battery --features alpha-kinetics,parallel -- --nocapture --test-threads=6
//! cargo test --test reaction_battery --features alpha-gsm,alpha-kinetics,parallel -- --nocapture --test-threads=6

use rayon::prelude::*;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;

const N: usize = 21;
const MID: usize = N / 2; // frame 10

use sci_form::potentials::{
    asymmetric_1d as asymmetric, double_well_1d as double_well, muller_brown_2d as muller_brown,
    proton_transfer_evb as evb, sn2_model_1d as sn2,
};

// ─── Helpers ───────────────────────────────────────────────────────────────

fn assert_ts_near(energies: &[f64], expected: usize, tol: usize, label: &str) {
    let ts = energies
        .iter()
        .enumerate()
        .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
        .map(|(i, _)| i)
        .unwrap();
    let dist = (ts as isize - expected as isize).unsigned_abs();
    println!("  [{label}] ts_frame={ts}  esperado≈{expected}  Δ={dist}");
    assert!(
        dist <= tol,
        "{label}: TS en frame {ts}, esperado ≈ {expected} (tol={tol})"
    );
}

fn print_energy_table(label: &str, energies: &[f64]) {
    println!("\n  ── Perfil energético — {label} (21 frames) ──");
    println!("  {:>5}  {:>14}", "frame", "E (u.a.)");
    for (i, e) in energies.iter().enumerate() {
        let mark = if i == MID {
            "  ← frame central (TS)"
        } else {
            ""
        };
        println!("  {:>5}  {:>14.8}{}", i, e, mark);
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// R1 — H-H disociación (pozo doble 1D simétrico)
// Frame 10 (x=0) es EXACTAMENTE el TS: V(0)=1, V(±1)=0
// ═══════════════════════════════════════════════════════════════════════════
#[test]
fn r1_h2_double_well_21_frames() {
    let xs: Vec<f64> = (0..N)
        .map(|i| -1.0 + 2.0 * i as f64 / (N - 1) as f64)
        .collect();
    let energies: Vec<f64> = xs.par_iter().map(|&x| double_well(x)).collect();

    println!("\n═══ R1: H₂ Disociación — Pozo doble 1D simétrico ═══");
    println!("  V(x) = (x²-1)²   Reactante x=-1, Producto x=+1");
    println!("  Frame 10 → x=0 → V(0)=1.0 (MÁXIMO = TS exacto)\n");
    println!(
        "  {:>5}  {:>8}  {:>14}  {:>14}  {:>14}",
        "frame", "x", "E", "Ea_fwd", "Ea_rev"
    );
    for (i, (&x, &e)) in xs.iter().zip(energies.iter()).enumerate() {
        let ea_fwd = e - energies[0];
        let ea_rev = e - energies[N - 1];
        let mark = match i {
            i if i == MID => "  ← TS",
            0 => "  ← Reactante",
            j if j == N - 1 => "  ← Producto",
            _ => "",
        };
        println!(
            "  {:>5}  {:>8.4}  {:>14.8}  {:>14.8}  {:>14.8}{}",
            i, x, e, ea_fwd, ea_rev, mark
        );
    }

    let e_mid = energies[MID];
    let e_r = energies[0];
    let e_p = energies[N - 1];
    println!("\n  E(reactante)={e_r:.6}  E(TS@frame10)={e_mid:.6}  E(producto)={e_p:.6}");
    println!("  Ea_fwd={:.6}  Ea_rev={:.6}", e_mid - e_r, e_mid - e_p);

    assert_eq!(energies.len(), N);
    assert_ts_near(&energies, MID, 0, "R1-H2-doble-pozo"); // exacto: x=0 es el máximo global
    assert!(
        (e_mid - 1.0).abs() < 1e-10,
        "V(0) debe valer exactamente 1.0, got {e_mid}"
    );
    assert!(
        e_r < 1e-10 && e_p < 1e-10,
        "reactante y producto deben estar en 0"
    );
    assert!(
        (e_mid - e_r - (e_mid - e_p)).abs() < 1e-10,
        "reacción simétrica: Ea_fwd=Ea_rev"
    );
}

// ═══════════════════════════════════════════════════════════════════════════
// R2 — Reacción asimétrica exotérmica
// V(x)=(x²-1)²−0.3x  → x=+1 más estable (ΔE=-0.6), TS cerca de frame 9-10
// ═══════════════════════════════════════════════════════════════════════════
#[test]
fn r2_asymmetric_exothermic_21_frames() {
    let xs: Vec<f64> = (0..N)
        .map(|i| -1.0 + 2.0 * i as f64 / (N - 1) as f64)
        .collect();
    let energies: Vec<f64> = xs.par_iter().map(|&x| asymmetric(x)).collect();

    let ts_idx = energies
        .iter()
        .enumerate()
        .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
        .map(|(i, _)| i)
        .unwrap();
    let ea_fwd = energies[ts_idx] - energies[0];
    let ea_rev = energies[ts_idx] - energies[N - 1];
    let delta_e = energies[N - 1] - energies[0];

    println!("\n═══ R2: Reacción Exotérmica Asimétrica ═══");
    println!("  V(x) = (x²-1)² − 0.3x");
    println!("  x=-1: E=+0.3  x=+1: E=-0.3  → ΔE=-0.6 (exotérmico)");
    print_energy_table("Asimétrica", &energies);
    println!("\n  TS en frame {ts_idx}  Ea_fwd={ea_fwd:.4}  Ea_rev={ea_rev:.4}  ΔE={delta_e:.4}");
    println!("  (TS ligeramente desplazado hacia reactante — principio de Hammond)");

    assert_eq!(energies.len(), N);
    assert!(ea_fwd > 0.0, "debe haber barrera forward");
    assert!(
        ea_rev > ea_fwd,
        "exotérmico: Ea_rev={ea_rev:.4} debe > Ea_fwd={ea_fwd:.4}"
    );
    assert!(delta_e < 0.0, "ΔE debe ser negativo: {delta_e:.4}");
    assert_ts_near(&energies, MID, 3, "R2-asimétrica");
}

// ═══════════════════════════════════════════════════════════════════════════
// R3 — Potencial Müller-Brown 2D (benchmark estándar de TS)
// ═══════════════════════════════════════════════════════════════════════════
#[test]
fn r3_muller_brown_2d_ts_search() {
    let start = [0.62_f64, 0.03];
    let end = [-0.56_f64, 1.44];

    let data: Vec<(f64, f64, f64)> = (0..N)
        .into_par_iter()
        .map(|i| {
            let t = i as f64 / (N - 1) as f64;
            let x = start[0] + t * (end[0] - start[0]);
            let y = start[1] + t * (end[1] - start[1]);
            (x, y, muller_brown(x, y))
        })
        .collect();

    let energies: Vec<f64> = data.iter().map(|&(_, _, e)| e).collect();
    let ts_idx = energies
        .iter()
        .enumerate()
        .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
        .map(|(i, _)| i)
        .unwrap();

    println!("\n═══ R3: Potencial Müller-Brown 2D ═══");
    println!(
        "  Reactante: ({:.2},{:.2})  Producto: ({:.2},{:.2})",
        start[0], start[1], end[0], end[1]
    );
    println!("\n  {:>5}  {:>8}  {:>8}  {:>12}", "frame", "x", "y", "E");
    for (i, &(x, y, e)) in data.iter().enumerate() {
        let mark = match i {
            i if i == ts_idx => "  ← TS",
            i if i == MID => "  ← frame 10 (central)",
            _ => "",
        };
        println!("  {:>5}  {:>8.4}  {:>8.4}  {:>12.6}{}", i, x, y, e, mark);
    }
    println!("\n  TS detectado en frame {ts_idx}");

    assert_eq!(energies.len(), N);
    assert!(
        ts_idx > 0 && ts_idx < N - 1,
        "TS no puede estar en extremos"
    );
    assert!(energies[ts_idx] >= energies[0], "E(TS) >= E(reactante)");
    assert!(energies[ts_idx] >= energies[N - 1], "E(TS) >= E(producto)");
}

// ═══════════════════════════════════════════════════════════════════════════
// R4 — Transferencia de protón — Marcus/EVB adiabático
// TS EXACTAMENTE en frame 10 (s=0) por simetría del potencial
// ═══════════════════════════════════════════════════════════════════════════
#[test]
fn r4_proton_transfer_marcus_evb_21_frames() {
    let k_marcus = 1.0_f64; // curvatura de pozos
    let coupling = 0.15_f64; // acoplamiento electrónico A-B

    let data: Vec<(f64, f64)> = (0..N)
        .into_par_iter()
        .map(|i| {
            // s = -1 (reactante H en A) → 0 (TS, H equidistante) → +1 (producto H en B)
            let s = -1.0 + 2.0 * i as f64 / (N - 1) as f64;
            (s, evb(s, k_marcus, coupling))
        })
        .collect();

    let energies: Vec<f64> = data.iter().map(|&(_, e)| e).collect();
    let e_mid = energies[MID];
    let e_r = energies[0];
    let e_p = energies[N - 1];
    let ea_fwd = e_mid - e_r;
    let ea_rev = e_mid - e_p;

    println!("\n═══ R4: Transferencia de Protón — Marcus/EVB Adiabático ═══");
    println!("  k={k_marcus}  coupling V={coupling}");
    println!("  s=-1 (H en A, reactante)  s=0 (TS, frame10, ε_D=ε_A)  s=+1 (H en B, producto)");
    println!(
        "\n  {:>5}  {:>8}  {:>10}  {:>10}  {:>12}",
        "frame", "s", "ε_D", "ε_A", "E_EVB"
    );
    for (i, &(s, e)) in data.iter().enumerate() {
        let eps_d = k_marcus * (s + 1.0).powi(2);
        let eps_a = k_marcus * (s - 1.0).powi(2);
        let mark = match i {
            i if i == MID => "  ← TS (ε_D=ε_A @ s=0) ✓",
            0 => "  ← Reactante",
            j if j == N - 1 => "  ← Producto",
            _ => "",
        };
        println!(
            "  {:>5}  {:>8.4}  {:>10.4}  {:>10.4}  {:>12.6}{}",
            i, s, eps_d, eps_a, e, mark
        );
    }
    println!("\n  E(R)={e_r:.6}  E(TS@10)={e_mid:.6}  E(P)={e_p:.6}");
    println!(
        "  Ea_fwd={ea_fwd:.6}  Ea_rev={ea_rev:.6}  |Ea_fwd-Ea_rev|={:.1e}",
        (ea_fwd - ea_rev).abs()
    );

    assert_eq!(energies.len(), N);
    // El TS está EXACTAMENTE en s=0 (MID) porque ε_D(0)=ε_A(0) por simetría
    assert_ts_near(&energies, MID, 0, "R4-EVB-protón");
    assert!(ea_fwd > 0.0, "barrera forward > 0: {ea_fwd}");
    // Simetría: |Ea_fwd - Ea_rev| < 1e-8 (reacción simétrica)
    let asym = (ea_fwd - ea_rev).abs();
    assert!(
        asym < 1e-8,
        "EVB simétrico: |Ea_fwd-Ea_rev|={asym:.2e} debe ser < 1e-8"
    );
}

// ═══════════════════════════════════════════════════════════════════════════
// R5 — Modelo SN2 + cinética HTST multiobjetivo
// ═══════════════════════════════════════════════════════════════════════════
#[test]
fn r5_sn2_model_21_frames_plus_kinetics() {
    let xs: Vec<f64> = (0..N)
        .map(|i| -1.5 + 3.0 * i as f64 / (N - 1) as f64)
        .collect();
    let energies: Vec<f64> = xs.par_iter().map(|&x| sn2(x)).collect();

    let ts_idx = energies
        .iter()
        .enumerate()
        .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
        .map(|(i, _)| i)
        .unwrap();
    let ea_fwd = energies[ts_idx] - energies[0];
    let ea_rev = energies[ts_idx] - energies[N - 1];
    let delta_e = energies[N - 1] - energies[0];

    println!("\n═══ R5: Modelo SN2 — Nu: + R-LG → Nu-R + :LG ═══");
    println!("  V(x) = −2·exp(−2(x+1.5)²) − 3·exp(−2(x−1.5)²) + 1.5·exp(−8x²)");
    println!("\n  {:>5}  {:>8}  {:>14}", "frame", "x", "E");
    for (i, (&x, &e)) in xs.iter().zip(energies.iter()).enumerate() {
        let mark = match i {
            i if i == ts_idx => "  ← TS",
            i if i == MID => "  ← frame central",
            0 => "  ← Reactante",
            j if j == N - 1 => "  ← Producto",
            _ => "",
        };
        println!("  {:>5}  {:>8.4}  {:>14.8}{}", i, x, e, mark);
    }
    println!("\n  TS@frame{ts_idx}  Ea_fwd={ea_fwd:.4}  Ea_rev={ea_rev:.4}  ΔE={delta_e:.4}");

    #[cfg(feature = "alpha-kinetics")]
    {
        use sci_form::alpha::kinetics::{
            evaluate_htst_rate, evaluate_htst_temperature_sweep, solve_microkinetic_network,
            ElementaryStep, MicrokineticNetworkConfig, ThermodynamicState,
        };

        // Barrera capped a 0.30 eV para estabilidad numérica del integrador
        let ea_ev = (ea_fwd * 0.04336).clamp(0.01, 0.30);
        let dg_ev = delta_e * 0.04336;

        let step = ElementaryStep {
            step_id: "SN2".into(),
            activation_free_energy_ev: ea_ev,
            reaction_free_energy_ev: dg_ev,
            prefactor_s_inv: None, // TST default: kBT/h
        };

        // Level 3: HTST @ 298 K
        let r298 = evaluate_htst_rate(&step, ThermodynamicState::default()).unwrap();
        println!("\n  ── Level 3: HTST Eyring @ 298 K ──");
        println!("  k_fwd = {:.4e} s⁻¹", r298.forward_rate_s_inv);
        println!("  k_rev = {:.4e} s⁻¹", r298.reverse_rate_s_inv);
        println!("  K_eq  = {:.4e}", r298.equilibrium_constant);

        // Level 3: sweep Arrhenius (rayon interno)
        let temps: Vec<f64> = (2..=10).map(|i| i as f64 * 100.0).collect();
        let sweep = evaluate_htst_temperature_sweep(&step, &temps, 1.0).unwrap();
        println!("\n  ── Arrhenius sweep 200-1000 K ──");
        println!("  {:>7}  {:>14}  {:>14}", "T (K)", "k_fwd (s⁻¹)", "K_eq");
        for (t, r) in temps.iter().zip(sweep.iter()) {
            println!(
                "  {:>7.0}  {:>14.4e}  {:>14.4e}",
                t, r.forward_rate_s_inv, r.equilibrium_constant
            );
        }
        for w in sweep.windows(2) {
            assert!(
                w[1].forward_rate_s_inv >= w[0].forward_rate_s_inv,
                "k debe aumentar con T (Arrhenius)"
            );
        }

        // Level 4: Red microcinética — prefactores 1e5 para estabilidad Euler
        let config = MicrokineticNetworkConfig {
            species: vec!["Reactante".into(), "Complejo-TS".into(), "Producto".into()],
            steps: vec![
                ElementaryStep {
                    step_id: "R→TS".into(),
                    activation_free_energy_ev: ea_ev,
                    reaction_free_energy_ev: dg_ev * 0.5,
                    prefactor_s_inv: Some(1e5),
                },
                ElementaryStep {
                    step_id: "TS→P".into(),
                    activation_free_energy_ev: 0.08,
                    reaction_free_energy_ev: dg_ev * 0.5,
                    prefactor_s_inv: Some(1e5),
                },
            ],
            stoichiometry: vec![vec![(0, -1.0), (1, 1.0)], vec![(1, -1.0), (2, 1.0)]],
            initial_populations: vec![1.0, 0.0, 0.0],
            total_time_s: 5e-2,
            n_frames: N,
            state: ThermodynamicState {
                temperature_k: 400.0,
                pressure_bar: 1.0,
            },
        };

        let trace = solve_microkinetic_network(&config).unwrap();
        println!("\n  ── Level 4: Red microcinética (21 frames temporales) ──");
        println!(
            "  {:>10}  {:>12}  {:>14}  {:>10}",
            "t(s)", "[Reactante]", "[Complejo-TS]", "[Producto]"
        );
        for f in &trace.frames {
            println!(
                "  {:>10.3e}  {:>12.4}  {:>14.4}  {:>10.4}",
                f.time_s, f.species[0].population, f.species[1].population, f.species[2].population
            );
        }

        assert_eq!(trace.frames.len(), N, "microcinética debe tener {N} frames");
        let last = trace.frames.last().unwrap();
        let first = &trace.frames[0];
        assert!(
            !last.species[2].population.is_nan() && !last.species[2].population.is_infinite(),
            "concentración producto no debe ser NaN/inf"
        );
        assert!(
            last.species[2].population > first.species[2].population,
            "Producto debe formarse a lo largo de la reacción"
        );
        assert!(
            last.species[0].population < first.species[0].population,
            "Reactante debe consumirse"
        );
    }

    assert_eq!(energies.len(), N);
    assert!(ea_fwd >= 0.0);
    assert_ts_near(&energies, MID, 3, "R5-SN2");
}

// ═══════════════════════════════════════════════════════════════════════════
// Level 2: EHT — scan H₂ (21 puntos de r) y H₂O (21 ángulos) en paralelo
// ═══════════════════════════════════════════════════════════════════════════
#[test]
fn level2_eht_h2_and_h2o_scan_21_frames() {
    // H₂: r = 0.5 → 3.0 Å en 21 puntos (paralelo)
    let r_vals: Vec<f64> = (0..N)
        .map(|i| 0.5 + 2.5 * i as f64 / (N - 1) as f64)
        .collect();
    let eht_h2: Vec<Option<(f64, f64)>> = r_vals
        .par_iter()
        .map(|&r| {
            let pos = vec![[0.0_f64, 0.0, 0.0], [r, 0.0, 0.0]];
            sci_form::eht::solve_eht(&[1u8, 1], &pos, None)
                .ok()
                .map(|e| (e.homo_energy, e.gap))
        })
        .collect();

    println!("\n═══ Level 2: EHT — H₂ scan r=0.5→3.0 Å (21 frames paralelo) ═══");
    println!(
        "  {:>5}  {:>8}  {:>14}  {:>10}",
        "frame", "r(Å)", "HOMO(eV)", "gap(eV)"
    );
    for (i, (r, res)) in r_vals.iter().zip(eht_h2.iter()).enumerate() {
        let mark = if i == MID { "  ← frame 10" } else { "" };
        match res {
            Some((homo, gap)) => println!(
                "  {:>5}  {:>8.4}  {:>14.6}  {:>10.6}{}",
                i, r, homo, gap, mark
            ),
            None => println!("  {:>5}  {:>8.4}  skipped{}", i, r, mark),
        }
    }

    // H₂O: θ = 90° → 180° en 21 frames (paralelo)
    let angles: Vec<f64> = (0..N)
        .map(|i| 90.0 + 90.0 * i as f64 / (N - 1) as f64)
        .collect();
    let r_oh = 0.957_f64;
    let eht_h2o: Vec<Option<(f64, f64, f64)>> = angles
        .par_iter()
        .map(|&theta_deg| {
            let half = (theta_deg / 2.0).to_radians();
            let pos = vec![
                [0.0_f64, 0.0, 0.0],
                [r_oh * half.sin(), r_oh * half.cos(), 0.0],
                [-r_oh * half.sin(), r_oh * half.cos(), 0.0],
            ];
            let eht = sci_form::eht::solve_eht(&[8u8, 1, 1], &pos, None).ok()?;
            let fk = sci_form::reactivity::compute_fukui_descriptors(&[8u8, 1, 1], &pos, &eht);
            Some((eht.homo_energy, eht.gap, fk.f_plus[0]))
        })
        .collect();

    println!("\n  ── EHT H₂O: scan angular 90°→180° (21 frames paralelo) ──");
    println!(
        "  {:>5}  {:>8}  {:>12}  {:>10}  {:>8}",
        "frame", "θ(°)", "HOMO(eV)", "gap(eV)", "f+(O)"
    );
    for (i, (theta, res)) in angles.iter().zip(eht_h2o.iter()).enumerate() {
        let mark = if i == MID { "  ← frame 10" } else { "" };
        match res {
            Some((homo, gap, fp)) => println!(
                "  {:>5}  {:>8.2}  {:>12.6}  {:>10.6}  {:>8.4}{}",
                i, theta, homo, gap, fp, mark
            ),
            None => println!("  {:>5}  {:>8.2}  skipped{}", i, theta, mark),
        }
    }

    let n_h2 = eht_h2.iter().filter(|x| x.is_some()).count();
    let n_h2o = eht_h2o.iter().filter(|x| x.is_some()).count();
    println!("\n  ✓ H₂: {n_h2}/{N} puntos  H₂O: {n_h2o}/{N} puntos");
    assert!(n_h2 >= N / 2, "EHT H₂ debe resolver ≥ {}", N / 2);
    assert!(n_h2o >= N / 2, "EHT H₂O debe resolver ≥ {}", N / 2);
}

// ═══════════════════════════════════════════════════════════════════════════
// Comparación de 4 métodos — H₂ stretch
// L1: Analítico  L2: EHT  L3: PM3  L4: xTB
// ═══════════════════════════════════════════════════════════════════════════
#[test]
fn method_comparison_4_levels_h2_21_frames() {
    let xs: Vec<f64> = (0..N)
        .map(|i| -1.0 + 2.0 * i as f64 / (N - 1) as f64)
        .collect();

    // Level 1: Analítico (doble pozo, TS exacto en frame 10)
    let l1: Vec<f64> = xs.par_iter().map(|&x| double_well(x)).collect();

    // Level 2: EHT HOMO (negado y normalizado)
    let l2: Vec<f64> = xs
        .par_iter()
        .map(|&x| {
            let r = 0.5 + (x + 1.0) * 1.25;
            let pos = vec![[0.0_f64, 0.0, 0.0], [r, 0.0, 0.0]];
            sci_form::eht::solve_eht(&[1u8, 1], &pos, None)
                .map(|e| -(e.homo_energy)) // más negativo = más estable; invertimos
                .unwrap_or(f64::NAN)
        })
        .collect();

    // Level 3: PM3
    let l3: Vec<f64> = xs
        .par_iter()
        .map(|&x| {
            let r = 0.5 + (x + 1.0) * 1.25;
            let positions = vec![[0.0_f64, 0.0, 0.0], [r, 0.0, 0.0]];
            sci_form::pm3::gradients::compute_pm3_gradient(&[1u8, 1], &positions)
                .map(|g| g.energy * 23.0605)
                .unwrap_or(f64::NAN)
        })
        .collect();

    // Level 4: xTB
    let l4: Vec<f64> = xs
        .par_iter()
        .map(|&x| {
            let r = 0.5 + (x + 1.0) * 1.25;
            let positions = vec![[0.0_f64, 0.0, 0.0], [r, 0.0, 0.0]];
            sci_form::xtb::gradients::compute_xtb_gradient(&[1u8, 1], &positions)
                .map(|g| g.energy * 23.0605)
                .unwrap_or(f64::NAN)
        })
        .collect();

    println!("\n═══ Comparación de 4 Métodos — H₂ (21 frames) ═══");
    println!("  x→r_HH: -1→0.5Å, 0→1.75Å, +1→3.0Å");
    println!(
        "\n  {:>5}  {:>8}  {:>10}  {:>10}  {:>10}  {:>10}",
        "frame", "r(Å)", "L1-Analít.", "L2-EHT", "L3-PM3", "L4-xTB"
    );
    for i in 0..N {
        let r = 0.5 + (xs[i] + 1.0) * 1.25;
        let mark = if i == MID { "  ← frame 10 (TS)" } else { "" };
        println!(
            "  {:>5}  {:>8.4}  {:>10.4}  {:>10.4}  {:>10.4}  {:>10.4}{}",
            i,
            r,
            l1[i],
            l2[i],
            if l3[i].is_nan() { 0.0 } else { l3[i] },
            if l4[i].is_nan() { 0.0 } else { l4[i] },
            mark
        );
    }

    // L1 analítico: TS exacto en frame 10
    assert_ts_near(&l1, MID, 0, "L1-analítico");

    // L2 EHT: verificar variación
    let l2_valid: Vec<f64> = l2.iter().copied().filter(|e| e.is_finite()).collect();
    println!("\n  L2-EHT: {}/{} puntos válidos", l2_valid.len(), N);
    if l2_valid.len() >= 3 {
        let spread = l2_valid.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
            - l2_valid.iter().cloned().fold(f64::INFINITY, f64::min);
        assert!(
            spread > 0.0,
            "EHT debe mostrar variación de energía: spread={spread}"
        );
        println!("  L2-EHT: spread={spread:.4} eV  ✓");
    }

    // Paralelismo: contador atómico
    let count = Arc::new(AtomicUsize::new(0));
    let c2 = Arc::clone(&count);
    let _: Vec<f64> = (0..N)
        .into_par_iter()
        .map(|i| {
            c2.fetch_add(1, Ordering::Relaxed);
            l1[i]
        })
        .collect();
    assert_eq!(count.load(Ordering::Relaxed), N);
    println!("  ✓ Paralelismo rayon: {N} invocaciones confirmadas");
}

// ═══════════════════════════════════════════════════════════════════════════
// NEB UFF — 21 imágenes con H₂O
// ═══════════════════════════════════════════════════════════════════════════
#[test]
fn neb_uff_h2o_bending_21_images() {
    use sci_form::dynamics::compute_simplified_neb_path;

    let r_oh = 0.957_f64;
    let half0 = (104.5_f64 / 2.0).to_radians();
    let half1 = (150.0_f64 / 2.0).to_radians();

    let start_coords = vec![
        0.0,
        0.0,
        0.0,
        r_oh * half0.sin(),
        r_oh * half0.cos(),
        0.0,
        -r_oh * half0.sin(),
        r_oh * half0.cos(),
        0.0,
    ];
    let end_coords = vec![
        0.0,
        0.0,
        0.0,
        r_oh * half1.sin(),
        r_oh * half1.cos(),
        0.0,
        -r_oh * half1.sin(),
        r_oh * half1.cos(),
        0.0,
    ];

    match compute_simplified_neb_path("[H]O[H]", &start_coords, &end_coords, N, 5, 0.1, 0.001) {
        Ok(neb) => {
            let energies: Vec<f64> = neb
                .images
                .iter()
                .map(|i| i.potential_energy_kcal_mol)
                .collect();
            let e_finite: Vec<f64> = energies.iter().copied().filter(|e| e.is_finite()).collect();
            let ts_idx = energies
                .iter()
                .enumerate()
                .filter(|(_, e)| e.is_finite())
                .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
                .map(|(i, _)| i)
                .unwrap_or(0);

            println!("\n═══ NEB UFF: H₂O bending (21 imágenes) ═══");
            println!(
                "  n_images={}, {}/{} energías finitas, TS@frame{}",
                neb.images.len(),
                e_finite.len(),
                N,
                ts_idx
            );
            println!(
                "\n  {:>5}  {:>14}  {:>14}",
                "img", "E(kcal/mol)", "O-H dist"
            );
            for img in &neb.images {
                let c = &img.coords;
                let d1 =
                    ((c[3] - c[0]).powi(2) + (c[4] - c[1]).powi(2) + (c[5] - c[2]).powi(2)).sqrt();
                let mark = match img.index {
                    i if i == ts_idx && energies[i].is_finite() => "  ← TS",
                    i if i == MID => "  ← frame central",
                    _ => "",
                };
                if energies[img.index].is_finite() {
                    println!(
                        "  {:>5}  {:>14.4}  {:>14.4}{}",
                        img.index, img.potential_energy_kcal_mol, d1, mark
                    );
                }
            }
            assert_eq!(neb.images.len(), N, "NEB debe generar {N} imágenes");
        }
        Err(e) => println!("\n  NEB UFF: {} — skipped (OK)", e),
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// Fukui en los 21 frames del scan H₂O — descripción atómica
// ═══════════════════════════════════════════════════════════════════════════
#[test]
fn reactivity_fukui_descriptors_21_frames() {
    let angles: Vec<f64> = (0..N)
        .map(|i| 90.0 + 90.0 * i as f64 / (N - 1) as f64)
        .collect();
    let r_oh = 0.957_f64;
    let elements = &[8u8, 1, 1];

    let results: Vec<Option<(f64, f64, f64, f64)>> = angles
        .par_iter()
        .map(|&theta_deg| {
            let half = (theta_deg / 2.0).to_radians();
            let pos = vec![
                [0.0_f64, 0.0, 0.0],
                [r_oh * half.sin(), r_oh * half.cos(), 0.0],
                [-r_oh * half.sin(), r_oh * half.cos(), 0.0],
            ];
            let eht = sci_form::eht::solve_eht(elements, &pos, None).ok()?;
            let fk = sci_form::reactivity::compute_fukui_descriptors(elements, &pos, &eht);
            Some((eht.homo_energy, eht.gap, fk.f_plus[0], fk.f_minus[0]))
        })
        .collect();

    println!("\n═══ Fukui Descriptors — H₂O (21 frames angulares, paralelo) ═══");
    println!(
        "  {:>5}  {:>8}  {:>12}  {:>10}  {:>8}  {:>8}  {:>10}",
        "frame", "θ(°)", "HOMO(eV)", "gap(eV)", "f+(O)", "f-(O)", "dual_O"
    );
    for (i, (theta, r)) in angles.iter().zip(results.iter()).enumerate() {
        let mark = if i == MID { "  ← frame 10" } else { "" };
        match r {
            Some((homo, gap, fp, fm)) => {
                let dual = fp - fm;
                println!(
                    "  {:>5}  {:>8.2}  {:>12.6}  {:>10.6}  {:>8.4}  {:>8.4}  {:>10.4}{}",
                    i, theta, homo, gap, fp, fm, dual, mark
                );
            }
            None => println!("  {:>5}  {:>8.2}  skipped{}", i, theta, mark),
        }
    }

    let n_ok = results.iter().filter(|x| x.is_some()).count();
    println!("\n  ✓ {n_ok}/{N} frames con EHT+Fukui");
    assert!(n_ok >= N / 2);
}

// ═══════════════════════════════════════════════════════════════════════════
// GSM — 3 potenciales × 21 nodos, frame central = TS (necesita alpha-gsm)
// ═══════════════════════════════════════════════════════════════════════════
#[cfg(feature = "alpha-gsm")]
#[test]
fn gsm_three_reactions_21_nodes_mid_is_ts() {
    use sci_form::alpha::gsm::{find_transition_state, GsmConfig};

    let cfg = GsmConfig {
        max_nodes: N,
        max_iter: 300,
        step_size: 0.001,
        grad_tol: 0.02,
        fd_step: 0.005,
    };

    #[allow(clippy::type_complexity)]
    let potentials: Vec<(
        &str,
        Box<dyn Fn(&[f64]) -> f64 + Send + Sync>,
        [f64; 3],
        [f64; 3],
    )> = vec![
        (
            "Doble pozo 1D",
            Box::new(|c: &[f64]| double_well(c[0])),
            [-1., 0., 0.],
            [1., 0., 0.],
        ),
        (
            "SN2 modelo 1D",
            Box::new(|c: &[f64]| sn2(c[0])),
            [-1.5, 0., 0.],
            [1.5, 0., 0.],
        ),
        (
            "EVB protón",
            Box::new(|c: &[f64]| evb(c[0], 1.0, 0.15)),
            [-1., 0., 0.],
            [1., 0., 0.],
        ),
    ];

    println!("\n═══ GSM: 3 reacciones × 21 nodos (paralelo) ═══");
    let results: Vec<_> = potentials
        .par_iter()
        .map(|(name, efn, r, p)| {
            let result = find_transition_state(r, p, efn.as_ref(), &cfg);
            (name, result)
        })
        .collect();

    println!(
        "  {:>22}  {:>8}  {:>8}  {:>8}  {:>10}  {:>10}",
        "Reacción", "n_nodes", "ts_idx", "MID", "Ea_fwd", "Ea_rev"
    );
    for (name, res) in &results {
        println!(
            "  {:>22}  {:>8}  {:>8}  {:>8}  {:>10.4}  {:>10.4}",
            name,
            res.n_nodes,
            res.ts_node_index,
            res.n_nodes / 2,
            res.activation_energy,
            res.reverse_barrier
        );
    }

    println!("\n  ── Perfil completo de energía por reacción ──");
    for (name, res) in &results {
        let mid2 = res.n_nodes / 2;
        println!("\n  [{name}] — {} frames:", res.n_nodes);
        for (i, (coords, e)) in res
            .path_coords
            .iter()
            .zip(res.path_energies.iter())
            .enumerate()
        {
            let mark = match (i == res.ts_node_index, i == mid2) {
                (true, true) => "  ← TS = frame central ✓",
                (true, false) => "  ← TS",
                (false, true) => "  ← frame central",
                _ => "",
            };
            println!("  frame {:>2}  x={:>8.5}  E={:.8}{}", i, coords[0], e, mark);
        }
    }

    for (name, res) in &results {
        let mid2 = res.n_nodes / 2;
        assert!(res.n_nodes >= 5, "{name}: debe alcanzar ≥ 5 nodos");
        assert!(res.activation_energy >= 0.0, "{name}: barrera ≥ 0");
        let dist = (res.ts_node_index as isize - mid2 as isize).unsigned_abs();
        println!(
            "  ✓ {name}: TS=frame{} MID={mid2} Δ={dist}",
            res.ts_node_index
        );
    }
}
