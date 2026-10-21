//! Integration tests for E10: Constant Potential Method (core `sci_form::cpm`).

mod cpm_tests {
    use sci_form::cpm::*;

    fn water() -> (Vec<u8>, Vec<[f64; 3]>) {
        (
            vec![8, 1, 1],
            vec![[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]],
        )
    }

    fn ethanol() -> (Vec<u8>, Vec<[f64; 3]>) {
        (
            vec![6, 6, 8, 1, 1, 1, 1, 1, 1],
            vec![
                [0.0, 0.0, 0.0],
                [1.52, 0.0, 0.0],
                [2.14, 1.21, 0.0],
                [-0.36, 1.03, 0.0],
                [-0.36, -0.52, 0.89],
                [-0.36, -0.52, -0.89],
                [1.88, -0.52, 0.89],
                [1.88, -0.52, -0.89],
                [3.10, 1.21, 0.0],
            ],
        )
    }

    #[test]
    fn test_cpm_converges() {
        let (el, pos) = water();
        let result = compute_cpm_charges(&el, &pos, &CpmConfig::default());
        assert!(result.converged);
    }

    #[test]
    fn test_cpm_nonzero_charge() {
        let (el, pos) = water();
        // At non-neutral potential, total charge should deviate from zero
        let r1 = compute_cpm_charges(
            &el,
            &pos,
            &CpmConfig {
                mu_ev: -3.0,
                ..Default::default()
            },
        );
        let r2 = compute_cpm_charges(
            &el,
            &pos,
            &CpmConfig {
                mu_ev: -6.0,
                ..Default::default()
            },
        );
        assert!(
            (r1.total_charge - r2.total_charge).abs() > 0.01,
            "Different μ should give different Q"
        );
    }

    #[test]
    fn test_cpm_charge_increases_with_mu() {
        let (el, pos) = ethanol();
        let r_low = compute_cpm_charges(
            &el,
            &pos,
            &CpmConfig {
                mu_ev: -5.5,
                ..Default::default()
            },
        );
        let r_high = compute_cpm_charges(
            &el,
            &pos,
            &CpmConfig {
                mu_ev: -3.5,
                ..Default::default()
            },
        );
        assert!(r_high.total_charge > r_low.total_charge);
    }

    #[test]
    fn test_cpm_grand_potential_finite() {
        let (el, pos) = water();
        let result = compute_cpm_charges(&el, &pos, &CpmConfig::default());
        assert!(result.grand_potential.is_finite());
        assert!(result.electrostatic_energy.is_finite());
    }

    #[test]
    fn test_cpm_surface_scan() {
        let (el, pos) = water();
        let surface = compute_cpm_surface(&el, &pos, -5.5, -3.5, 10, 78.5);
        assert_eq!(surface.mu_values.len(), 10);
        assert_eq!(surface.total_charge.len(), 10);
        assert_eq!(surface.capacitance.len(), 10);
    }

    #[test]
    fn test_cpm_surface_charge_monotonic() {
        let (el, pos) = ethanol();
        let surface = compute_cpm_surface(&el, &pos, -5.5, -3.5, 20, 78.5);
        // Q should generally increase with μ
        let first = surface.total_charge[0];
        let last = *surface.total_charge.last().unwrap();
        assert!(
            last > first,
            "Q should increase: first={}, last={}",
            first,
            last
        );
    }

    #[test]
    fn test_cpm_surface_capacitance_positive() {
        let (el, pos) = water();
        let surface = compute_cpm_surface(&el, &pos, -5.0, -4.0, 15, 78.5);
        let positive_count = surface.capacitance.iter().filter(|&&c| c > 0.0).count();
        assert!(
            positive_count > surface.capacitance.len() / 2,
            "Most capacitance values should be positive"
        );
    }

    #[test]
    fn test_cpm_dielectric_effect() {
        let (el, pos) = water();
        let r_water = compute_cpm_charges(
            &el,
            &pos,
            &CpmConfig {
                dielectric: 78.5,
                mu_ev: -4.0,
                ..Default::default()
            },
        );
        let r_vacuum = compute_cpm_charges(
            &el,
            &pos,
            &CpmConfig {
                dielectric: 1.0,
                mu_ev: -4.0,
                ..Default::default()
            },
        );
        // In vacuum (ε=1), Coulomb coupling is stronger → different charge distribution
        assert!(
            (r_water.total_charge - r_vacuum.total_charge).abs() > 0.001,
            "Dielectric should affect charges"
        );
    }

    // ── F95: bilayer capacitor vs LAMMPS `fix electrode/conp` ──────────
    #[derive(serde::Deserialize)]
    struct CpmLammpsVoltage {
        voltage: f64,
        lammps_qbot: f64,
        lammps_qtop: f64,
        lammps_charges: Vec<f64>,
    }

    #[derive(serde::Deserialize)]
    struct CpmLammpsFixture {
        elements: Vec<u8>,
        coords: Vec<f64>,
        nplate: usize,
        voltages: Vec<CpmLammpsVoltage>,
    }

    fn load_cpm_lammps() -> CpmLammpsFixture {
        let path = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/tests/fixtures/validation/cpm_lammps.json"
        );
        let text = std::fs::read_to_string(path).expect("cpm_lammps.json");
        serde_json::from_str(&text).expect("parse cpm fixture")
    }

    fn cpm_run(elements: &[u8], pos: &[[f64; 3]], mus: Vec<f64>) -> Vec<f64> {
        compute_cpm_charges(
            elements,
            pos,
            &CpmConfig {
                mu_ev: 0.0,
                dielectric: 1.0,
                max_iter: 200,
                charge_tol: 1e-9,
                mu_per_atom: Some(mus),
            },
        )
        .charges
    }

    /// F95: our CPM (per-atom μ, direct LU solve) against LAMMPS conp on a
    /// graphene bilayer capacitor. Our model is point-charge + hardness
    /// while LAMMPS uses Gaussian-smeared Ewald electrostatics (no χ/η),
    /// so the gate compares the voltage-INDUCED charges
    /// (ΔQ = Q(±V) − Q(0), where χ cancels exactly by linearity) and locks
    /// in the quantified model gap: totals within 2.5x, pattern cosine
    /// similarity ≥ 0.3, exact antisymmetry/neutrality. Measured at V=±1:
    /// ours ΔQ = ±0.666e vs LAMMPS ±0.379e (1.76x), similarity 0.44.
    /// Path forward: Gaussian smearing in our J (ELECTRODE-style) and/or
    /// retuned onsite hardness (ours 5.0 eV vs Gaussian self ~22 eV).
    #[test]
    fn test_cpm_bilayer_matches_lammps() {
        let fix = load_cpm_lammps();
        assert!(fix.voltages.len() >= 3);
        let pos: Vec<[f64; 3]> = fix.coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();
        let n = pos.len();
        let np = fix.nplate;

        // Fixture quality: LAMMPS response is exactly linear; cell neutral.
        let caps: Vec<f64> = fix
            .voltages
            .iter()
            .map(|v| v.lammps_qtop / v.voltage)
            .collect();
        for c in &caps {
            assert!(
                (c - caps[0]).abs() / caps[0] < 0.01,
                "LAMMPS Q(V) not linear: {caps:?}"
            );
        }
        for v in &fix.voltages {
            assert!(
                (v.lammps_qbot + v.lammps_qtop).abs() < 1e-8,
                "LAMMPS cell not neutral at V={}",
                v.voltage
            );
        }

        let q0 = cpm_run(&fix.elements, &pos, vec![0.0; n]);
        let mut worst_ratio = 0.0f64;
        let mut worst_sim = 1.0f64;
        for v in &fix.voltages {
            let mut mus = vec![0.0; n];
            for m in mus.iter_mut().take(np) {
                *m = -v.voltage;
            }
            for m in mus.iter_mut().skip(np) {
                *m = v.voltage;
            }
            let qv = cpm_run(&fix.elements, &pos, mus);
            let ind: Vec<f64> = qv.iter().zip(q0.iter()).map(|(a, b)| a - b).collect();
            let qb: f64 = ind[..np].iter().sum();
            let qt: f64 = ind[np..].iter().sum();
            // Exact antisymmetry (symmetric plates, uniform χ).
            assert!(
                (qb + qt).abs() < 1e-8,
                "ours not antisymmetric at V={}: {qb} {qt}",
                v.voltage
            );
            // Scale: documented model gap (Gaussian self vs hardness).
            let ratio = (qb / v.lammps_qbot).abs();
            worst_ratio = worst_ratio.max(ratio.max(1.0 / ratio));
            // Pattern: cosine similarity of induced per-atom charges.
            let dot: f64 = ind
                .iter()
                .zip(v.lammps_charges.iter())
                .map(|(a, b)| a * b)
                .sum();
            let no: f64 = ind.iter().map(|a| a * a).sum::<f64>().sqrt();
            let nl: f64 = v.lammps_charges.iter().map(|a| a * a).sum::<f64>().sqrt();
            let sim = dot / (no * nl);
            worst_sim = worst_sim.min(sim);
            eprintln!(
                "F95 V=±{}: ours ΔQ={:+.4}/{:+.4} LAMMPS {:+.4}/{:+.4} ratio={:.3} sim={:.3}",
                v.voltage, qb, qt, v.lammps_qbot, v.lammps_qtop, ratio, sim
            );
        }
        assert!(
            worst_ratio <= 2.5,
            "plate-Q scale ratio {worst_ratio:.3} exceeds 2.5"
        );
        assert!(
            worst_sim >= 0.3,
            "pattern similarity {worst_sim:.3} below 0.3"
        );
    }
}
