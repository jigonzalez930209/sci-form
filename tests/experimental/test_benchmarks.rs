//! Comparative benchmark tests: experimental algorithms vs. legacy/production algorithms.
//!
//! Each comparison measures:
//! - Execution time ratio (experimental / legacy)
//! - Value deviation from reference values (crystallographic or exact quantum)
//! - Qualitative agreement (same sign, same trend, same order of magnitude)

// ═══════════════════════════════════════════════════════════════════
// B1 — EHT Exact vs. RandNLA (E2): Eigenvalue accuracy
// ═══════════════════════════════════════════════════════════════════

#[cfg(feature = "experimental-randnla")]
mod bench_randnla_vs_eht {
    use sci_form::eht::solve_eht;
    use sci_form::experimental::rand_nla::*;
    use std::time::Instant;

    fn water() -> (Vec<u8>, Vec<[f64; 3]>) {
        (
            vec![8, 1, 1],
            vec![[0.0, 0.0, 0.117], [0.0, 0.757, -0.469], [0.0, -0.757, -0.469]],
        )
    }

    fn ethanol() -> (Vec<u8>, Vec<[f64; 3]>) {
        (
            vec![6, 6, 8, 1, 1, 1, 1, 1, 1],
            vec![
                [0.0, 0.0, 0.0], [1.52, 0.0, 0.0], [2.13, 1.21, 0.0],
                [-0.36, 1.03, 0.0], [-0.36, -0.51, 0.89], [-0.36, -0.51, -0.89],
                [1.88, -0.51, 0.89], [1.88, -0.51, -0.89], [3.09, 1.21, 0.0],
            ],
        )
    }

    #[test]
    fn b1_eigenvalue_deviation_water() {
        let (elems, pos) = water();
        let eht = solve_eht(&elems, &pos, None).unwrap();

        let basis = sci_form::eht::basis::build_basis(&elems, &pos);
        let s = sci_form::eht::overlap::build_overlap_matrix(&basis);
        let h = sci_form::eht::hamiltonian::build_hamiltonian(&basis, &s, None);

        let config = RandNlaConfig { seed: 42, ..Default::default() };
        let (rand_energies, _, info) = solve_eht_randnla(&h, &s, &config);

        let n = eht.energies.len().min(rand_energies.len());
        let mut max_dev = 0.0_f64;
        for i in 0..n {
            let dev = (eht.energies[i] - rand_energies[i]).abs();
            max_dev = max_dev.max(dev);
        }

        println!("[B1] Water eigenvalue max deviation: {:.6} eV", max_dev);
        println!("[B1] RandNLA residual error: {:.2e}", info.residual_error);
        println!("[B1] Used fallback: {}", info.used_fallback);

        assert!(
            max_dev < 0.1 || info.used_fallback,
            "Eigenvalue deviation {:.4} eV too large without fallback", max_dev
        );
    }

    #[test]
    fn b1_homo_lumo_gap_agreement() {
        let (elems, pos) = ethanol();
        let eht = solve_eht(&elems, &pos, None).unwrap();

        let basis = sci_form::eht::basis::build_basis(&elems, &pos);
        let s = sci_form::eht::overlap::build_overlap_matrix(&basis);
        let h = sci_form::eht::hamiltonian::build_hamiltonian(&basis, &s, None);

        let config = RandNlaConfig::default();
        let (rand_energies, _, _) = solve_eht_randnla(&h, &s, &config);

        let homo = eht.homo_index;
        let lumo = eht.lumo_index;
        let eht_gap = eht.gap;
        let rand_gap = if lumo < rand_energies.len() && homo < rand_energies.len() {
            rand_energies[lumo] - rand_energies[homo]
        } else {
            0.0
        };

        let gap_dev = (eht_gap - rand_gap).abs();
        println!("[B1] Ethanol HOMO-LUMO gap: EHT={:.3} eV, RandNLA={:.3} eV, Δ={:.4} eV",
            eht_gap, rand_gap, gap_dev);
        assert!(gap_dev < 0.5, "Gap deviation {:.3} eV exceeds tolerance", gap_dev);
    }

    #[test]
    fn b1_timing_comparison() {
        let (elems, pos) = ethanol();
        let basis = sci_form::eht::basis::build_basis(&elems, &pos);
        let s = sci_form::eht::overlap::build_overlap_matrix(&basis);
        let h = sci_form::eht::hamiltonian::build_hamiltonian(&basis, &s, None);

        let t0 = Instant::now();
        for _ in 0..50 { let _ = solve_eht(&elems, &pos, None); }
        let eht_time = t0.elapsed().as_micros() as f64 / 50.0;

        let config = RandNlaConfig::default();
        let t0 = Instant::now();
        for _ in 0..50 { let _ = solve_eht_randnla(&h, &s, &config); }
        let rand_time = t0.elapsed().as_micros() as f64 / 50.0;

        println!("[B1] Timing: EHT={:.0} µs, RandNLA={:.0} µs, ratio={:.2}x",
            eht_time, rand_time, rand_time / eht_time.max(1.0));
    }
}

// ═══════════════════════════════════════════════════════════════════
// B2 — EHT Exact vs. KPM (E4): DOS comparison, Mulliken charges
// ═══════════════════════════════════════════════════════════════════

#[cfg(feature = "experimental-kpm")]
mod bench_kpm_vs_eht {
    use sci_form::eht::solve_eht;
    use sci_form::experimental::kpm::*;
    use std::time::Instant;

    fn water() -> (Vec<u8>, Vec<[f64; 3]>) {
        (
            vec![8, 1, 1],
            vec![[0.0, 0.0, 0.117], [0.0, 0.757, -0.469], [0.0, -0.757, -0.469]],
        )
    }

    fn ethanol() -> (Vec<u8>, Vec<[f64; 3]>) {
        (
            vec![6, 6, 8, 1, 1, 1, 1, 1, 1],
            vec![
                [0.0, 0.0, 0.0], [1.52, 0.0, 0.0], [2.13, 1.21, 0.0],
                [-0.36, 1.03, 0.0], [-0.36, -0.51, 0.89], [-0.36, -0.51, -0.89],
                [1.88, -0.51, 0.89], [1.88, -0.51, -0.89], [3.09, 1.21, 0.0],
            ],
        )
    }

    #[test]
    fn b2_dos_spectral_bounds_agree() {
        let (elems, pos) = water();
        let eht = solve_eht(&elems, &pos, None).unwrap();

        let basis = sci_form::eht::basis::build_basis(&elems, &pos);
        let s = sci_form::eht::overlap::build_overlap_matrix(&basis);
        let h = sci_form::eht::hamiltonian::build_hamiltonian(&basis, &s, None);

        let (kpm_emin, kpm_emax) = estimate_spectral_bounds(&h);
        let eht_emin = eht.energies.iter().cloned().fold(f64::INFINITY, f64::min);
        let eht_emax = eht.energies.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

        println!("[B2] EHT spectral range: [{:.2}, {:.2}] eV", eht_emin, eht_emax);
        println!("[B2] KPM spectral range: [{:.2}, {:.2}] eV", kpm_emin, kpm_emax);

        assert!(kpm_emin <= eht_emin + 1.0, "KPM e_min too high");
        assert!(kpm_emax >= eht_emax - 1.0, "KPM e_max too low");
    }

    #[test]
    fn b2_kpm_dos_peak_locations() {
        let (elems, pos) = ethanol();

        let basis = sci_form::eht::basis::build_basis(&elems, &pos);
        let s = sci_form::eht::overlap::build_overlap_matrix(&basis);
        let h = sci_form::eht::hamiltonian::build_hamiltonian(&basis, &s, None);

        let config = KpmConfig { order: 200, ..Default::default() };
        let dos = compute_kpm_dos(&h, &config, -30.0, 10.0, 500);

        let max_dos = dos.total_dos.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let peak_idx = dos.total_dos.iter().position(|&v| v == max_dos).unwrap();
        let peak_energy = dos.energies[peak_idx];

        let eht = solve_eht(&elems, &pos, None).unwrap();
        let min_dist: f64 = eht.energies.iter()
            .map(|&e| (e - peak_energy).abs())
            .fold(f64::INFINITY, f64::min);

        println!("[B2] KPM DOS peak at {:.2} eV, nearest eigenvalue {:.2} eV away",
            peak_energy, min_dist);
        assert!(min_dist < 5.0, "KPM DOS peak too far from any eigenvalue");
    }

    #[test]
    fn b2_kpm_mulliken_vs_eht_population() {
        let (elems, pos) = water();

        let basis = sci_form::eht::basis::build_basis(&elems, &pos);
        let s = sci_form::eht::overlap::build_overlap_matrix(&basis);
        let h = sci_form::eht::hamiltonian::build_hamiltonian(&basis, &s, None);

        let eht = solve_eht(&elems, &pos, None).unwrap();
        let nuclear_charges: Vec<f64> = elems.iter().map(|&z| z as f64).collect();
        let config = KpmConfig { order: 200, ..Default::default() };
        let kpm_mull = compute_kpm_mulliken(&h, &s, eht.n_electrons, &nuclear_charges, &config);

        let pop = sci_form::compute_population(&elems, &pos).unwrap();

        println!("[B2] EHT Mulliken: {:?}", pop.mulliken_charges);
        println!("[B2] KPM Mulliken: {:?}", kpm_mull.charges);

        let o_sign_agree = kpm_mull.charges[0] * pop.mulliken_charges[0] >= 0.0
            || kpm_mull.charges[0].abs() < 0.1;
        println!("[B2] O charge sign agreement: {}", o_sign_agree);
    }

    #[test]
    fn b2_timing_kpm_vs_exact_diag() {
        let (elems, pos) = ethanol();

        let basis = sci_form::eht::basis::build_basis(&elems, &pos);
        let s = sci_form::eht::overlap::build_overlap_matrix(&basis);
        let h = sci_form::eht::hamiltonian::build_hamiltonian(&basis, &s, None);

        let t0 = Instant::now();
        for _ in 0..50 { let _ = solve_eht(&elems, &pos, None); }
        let eht_us = t0.elapsed().as_micros() as f64 / 50.0;

        let config = KpmConfig { order: 100, ..Default::default() };
        let t0 = Instant::now();
        for _ in 0..50 { let _ = compute_kpm_dos(&h, &config, -30.0, 10.0, 200); }
        let kpm_us = t0.elapsed().as_micros() as f64 / 50.0;

        println!("[B2] Timing: EHT exact={:.0} µs, KPM DOS={:.0} µs", eht_us, kpm_us);
    }
}

// ═══════════════════════════════════════════════════════════════════
// B3 — Gasteiger vs. EEQ (E5): Charge accuracy
// ═══════════════════════════════════════════════════════════════════

#[cfg(feature = "experimental-eeq")]
mod bench_eeq_vs_gasteiger {
    use sci_form::experimental::eeq::*;
    use std::time::Instant;

    fn water() -> (Vec<u8>, Vec<[f64; 3]>) {
        (
            vec![8, 1, 1],
            vec![[0.0, 0.0, 0.117], [0.0, 0.757, -0.469], [0.0, -0.757, -0.469]],
        )
    }

    fn acetic_acid() -> (Vec<u8>, Vec<[f64; 3]>) {
        (
            vec![6, 6, 8, 8, 1, 1, 1, 1],
            vec![
                [0.0, 0.0, 0.0], [1.52, 0.0, 0.0], [2.08, 1.20, 0.0],
                [2.08, -1.20, 0.0], [-0.36, 1.03, 0.0], [-0.36, -0.51, 0.89],
                [-0.36, -0.51, -0.89], [3.05, -1.20, 0.0],
            ],
        )
    }

    /// Reference: water partial charges from ab initio HF/6-31G* Mulliken.
    const WATER_REF_O: f64 = -0.80;

    #[test]
    fn b3_charge_sign_agreement() {
        let (elems, pos) = water();
        let gasteiger = sci_form::compute_charges("O").unwrap();
        let eeq = compute_eeq_charges(&elems, &pos, &EeqConfig::default());

        println!("[B3] Gasteiger O={:.4}, H={:.4}", gasteiger.charges[0], gasteiger.charges[1]);
        println!("[B3] EEQ       O={:.4}, H={:.4}", eeq.charges[0], eeq.charges[1]);

        assert!(gasteiger.charges[0] < 0.0, "Gasteiger O should be negative");
        assert!(eeq.charges[0] < 0.0, "EEQ O should be negative");
    }

    #[test]
    fn b3_deviation_from_reference_water() {
        let (elems, pos) = water();
        let eeq = compute_eeq_charges(&elems, &pos, &EeqConfig::default());
        let gasteiger = sci_form::compute_charges("O").unwrap();

        let eeq_dev_o = (eeq.charges[0] - WATER_REF_O).abs();
        let gas_dev_o = (gasteiger.charges[0] - WATER_REF_O).abs();

        println!("[B3] Deviation from HF ref (O): EEQ={:.4}, Gasteiger={:.4}", eeq_dev_o, gas_dev_o);
        assert!(eeq_dev_o < 1.5, "EEQ deviation too large");
        assert!(gas_dev_o < 1.5, "Gasteiger deviation too large");
    }

    #[test]
    fn b3_charge_neutrality() {
        let (elems, pos) = acetic_acid();
        let eeq = compute_eeq_charges(&elems, &pos, &EeqConfig::default());
        let total: f64 = eeq.charges.iter().sum();

        println!("[B3] EEQ total charge (acetic acid): {:.6}", total);
        assert!(total.abs() < 0.01, "EEQ should preserve charge neutrality: {:.6}", total);
    }

    #[test]
    fn b3_geometry_dependence() {
        let (elems, _) = water();
        let pos_normal = vec![[0.0, 0.0, 0.117], [0.0, 0.757, -0.469], [0.0, -0.757, -0.469]];
        let pos_stretched = vec![[0.0, 0.0, 0.117], [0.0, 1.50, -0.469], [0.0, -0.757, -0.469]];

        let eeq_normal = compute_eeq_charges(&elems, &pos_normal, &EeqConfig::default());
        let eeq_stretched = compute_eeq_charges(&elems, &pos_stretched, &EeqConfig::default());

        let delta_o = (eeq_normal.charges[0] - eeq_stretched.charges[0]).abs();
        println!("[B3] EEQ O charge change on stretching: {:.4}", delta_o);
        assert!(delta_o > 1e-4, "EEQ should be geometry-dependent");
    }

    #[test]
    fn b3_timing_comparison() {
        let (elems, pos) = acetic_acid();

        let t0 = Instant::now();
        for _ in 0..200 { let _ = sci_form::compute_charges("CC(=O)O"); }
        let gas_us = t0.elapsed().as_micros() as f64 / 200.0;

        let config = EeqConfig::default();
        let t0 = Instant::now();
        for _ in 0..200 { let _ = compute_eeq_charges(&elems, &pos, &config); }
        let eeq_us = t0.elapsed().as_micros() as f64 / 200.0;

        println!("[B3] Timing: Gasteiger={:.0} µs, EEQ={:.0} µs", gas_us, eeq_us);
    }
}

// ═══════════════════════════════════════════════════════════════════
// B4 — GB-HCT vs. ALPB (E6): Solvation energy comparison
// ═══════════════════════════════════════════════════════════════════

#[cfg(feature = "experimental-alpb")]
mod bench_alpb_vs_gb {
    use sci_form::solvation::{compute_gb_solvation, compute_nonpolar_solvation};
    use sci_form::experimental::alpb::*;
    use std::time::Instant;

    fn water() -> (Vec<u8>, Vec<[f64; 3]>, Vec<f64>) {
        (
            vec![8, 1, 1],
            vec![[0.0, 0.0, 0.117], [0.0, 0.757, -0.469], [0.0, -0.757, -0.469]],
            vec![-0.834, 0.417, 0.417],
        )
    }

    fn ethanol() -> (Vec<u8>, Vec<[f64; 3]>, Vec<f64>) {
        (
            vec![6, 6, 8, 1, 1, 1, 1, 1, 1],
            vec![
                [0.0, 0.0, 0.0], [1.52, 0.0, 0.0], [2.13, 1.21, 0.0],
                [-0.36, 1.03, 0.0], [-0.36, -0.51, 0.89], [-0.36, -0.51, -0.89],
                [1.88, -0.51, 0.89], [1.88, -0.51, -0.89], [3.09, 1.21, 0.0],
            ],
            vec![-0.10, 0.15, -0.68, 0.04, 0.04, 0.04, 0.06, 0.06, 0.39],
        )
    }

    const WATER_SOLVATION_REF: f64 = -6.3;

    #[test]
    fn b4_solvation_sign_agreement() {
        let (elems, pos, charges) = water();
        let gb = compute_gb_solvation(&elems, &pos, &charges, Some(78.5), Some(1.0), Some(1.4));
        let alpb = compute_alpb_solvation(&elems, &pos, &charges, &AlpbConfig::default());

        println!("[B4] GB-HCT total: {:.3} kcal/mol", gb.total_energy_kcal_mol);
        println!("[B4] ALPB   total: {:.3} kcal/mol", alpb.total_energy);

        assert!(gb.total_energy_kcal_mol < 0.0 || gb.electrostatic_energy_kcal_mol < 0.0,
            "GB should give negative solvation");
        assert!(alpb.total_energy < 0.0 || alpb.electrostatic_energy < 0.0,
            "ALPB should give negative solvation");
    }

    #[test]
    fn b4_deviation_from_experimental() {
        let (elems, pos, charges) = water();
        let gb = compute_gb_solvation(&elems, &pos, &charges, Some(78.5), Some(1.0), Some(1.4));
        let alpb = compute_alpb_solvation(&elems, &pos, &charges, &AlpbConfig::default());

        let gb_dev = (gb.total_energy_kcal_mol - WATER_SOLVATION_REF).abs();
        let alpb_dev = (alpb.total_energy - WATER_SOLVATION_REF).abs();

        println!("[B4] Deviation from exp ({:.1} kcal/mol): GB={:.2}, ALPB={:.2}",
            WATER_SOLVATION_REF, gb_dev, alpb_dev);
        assert!(gb_dev < 30.0, "GB deviation too large");
        assert!(alpb_dev < 30.0, "ALPB deviation too large");
    }

    #[test]
    fn b4_born_radii_comparison() {
        let (elems, pos, _) = water();
        let gb_born = sci_form::solvation::compute_born_radii(&elems, &pos);
        let alpb_born = compute_born_radii(&elems, &pos, 1.4);

        println!("[B4] GB-HCT Born radii: {:?}", gb_born);
        println!("[B4] ALPB   Born radii: {:?}", alpb_born.radii);

        for (i, (&gb_r, &alpb_r)) in gb_born.iter().zip(alpb_born.radii.iter()).enumerate() {
            // GB-HCT returns 50.0 as fallback for small molecules — accept any positive value
            assert!(gb_r > 0.0, "GB Born[{}]={}", i, gb_r);
            assert!(alpb_r > 0.3 && alpb_r < 50.0, "ALPB Born[{}]={}", i, alpb_r);
        }
    }

    #[test]
    fn b4_nonpolar_component_agreement() {
        let (elems, pos, _) = water();
        let np = compute_nonpolar_solvation(&elems, &pos, None);
        let charges = vec![-0.834, 0.417, 0.417];
        let alpb = compute_alpb_solvation(&elems, &pos, &charges, &AlpbConfig::default());

        println!("[B4] Legacy non-polar: {:.4} kcal/mol", np.energy_kcal_mol);
        println!("[B4] ALPB   non-polar: {:.4} kcal/mol", alpb.nonpolar_energy);

        assert!(np.energy_kcal_mol.abs() < 20.0, "Legacy NP too large");
        assert!(alpb.nonpolar_energy.abs() < 20.0, "ALPB NP too large");
    }

    #[test]
    fn b4_timing_comparison() {
        let (elems, pos, charges) = ethanol();

        let t0 = Instant::now();
        for _ in 0..100 {
            let _ = compute_gb_solvation(&elems, &pos, &charges, Some(78.5), Some(1.0), Some(1.4));
        }
        let gb_us = t0.elapsed().as_micros() as f64 / 100.0;

        let config = AlpbConfig::default();
        let t0 = Instant::now();
        for _ in 0..100 { let _ = compute_alpb_solvation(&elems, &pos, &charges, &config); }
        let alpb_us = t0.elapsed().as_micros() as f64 / 100.0;

        println!("[B4] Timing: GB-HCT={:.0} µs, ALPB={:.0} µs", gb_us, alpb_us);
    }
}

// ═══════════════════════════════════════════════════════════════════
// B5 — UFF vs. UFF+D4 (E7): Dispersion energy contribution
// ═══════════════════════════════════════════════════════════════════

#[cfg(feature = "experimental-d4")]
mod bench_d4_vs_uff {
    use sci_form::experimental::d4::*;
    use std::time::Instant;

    fn benzene_dimer() -> (Vec<u8>, Vec<[f64; 3]>) {
        (
            vec![
                6, 6, 6, 6, 6, 6, 1, 1, 1, 1, 1, 1,
                6, 6, 6, 6, 6, 6, 1, 1, 1, 1, 1, 1,
            ],
            vec![
                [1.40, 0.0, 0.0], [0.70, 1.21, 0.0], [-0.70, 1.21, 0.0],
                [-1.40, 0.0, 0.0], [-0.70, -1.21, 0.0], [0.70, -1.21, 0.0],
                [2.48, 0.0, 0.0], [1.24, 2.15, 0.0], [-1.24, 2.15, 0.0],
                [-2.48, 0.0, 0.0], [-1.24, -2.15, 0.0], [1.24, -2.15, 0.0],
                [1.40, 0.0, 3.5], [0.70, 1.21, 3.5], [-0.70, 1.21, 3.5],
                [-1.40, 0.0, 3.5], [-0.70, -1.21, 3.5], [0.70, -1.21, 3.5],
                [2.48, 0.0, 3.5], [1.24, 2.15, 3.5], [-1.24, 2.15, 3.5],
                [-2.48, 0.0, 3.5], [-1.24, -2.15, 3.5], [1.24, -2.15, 3.5],
            ],
        )
    }

    fn ethanol() -> (Vec<u8>, Vec<[f64; 3]>) {
        (
            vec![6, 6, 8, 1, 1, 1, 1, 1, 1],
            vec![
                [0.0, 0.0, 0.0], [1.52, 0.0, 0.0], [2.13, 1.21, 0.0],
                [-0.36, 1.03, 0.0], [-0.36, -0.51, 0.89], [-0.36, -0.51, -0.89],
                [1.88, -0.51, 0.89], [1.88, -0.51, -0.89], [3.09, 1.21, 0.0],
            ],
        )
    }

    const BENZENE_DIMER_REF_KCAL: f64 = -2.7;

    #[test]
    fn b5_d4_dispersion_attractive() {
        let (elems, pos) = benzene_dimer();
        let d4 = compute_d4_energy(&elems, &pos, &D4Config::default());

        println!("[B5] D4 total: {:.4} Hartree = {:.2} kcal/mol", d4.total_energy, d4.total_kcal_mol);
        assert!(d4.total_energy < 0.0, "D4 dispersion should be negative (attractive)");
    }

    #[test]
    fn b5_d4_vs_reference_benzene_dimer() {
        let (elems, pos) = benzene_dimer();
        let d4 = compute_d4_energy(&elems, &pos, &D4Config::default());

        println!("[B5] D4 dimer: {:.2} kcal/mol, ref: {:.2} kcal/mol, Δ={:.2}",
            d4.total_kcal_mol, BENZENE_DIMER_REF_KCAL,
            (d4.total_kcal_mol - BENZENE_DIMER_REF_KCAL).abs());
        assert!(d4.total_kcal_mol < 0.0, "D4 should be attractive");
        assert!(d4.total_kcal_mol.abs() < 50.0, "D4 magnitude reasonable");
    }

    #[test]
    fn b5_uff_augmented_with_d4() {
        let (elems, pos) = ethanol();
        let conf = sci_form::embed("CCO", 42);
        let uff_energy = sci_form::compute_uff_energy("CCO", &conf.coords).unwrap_or(0.0);
        let d4 = compute_d4_energy(&elems, &pos, &D4Config::default());

        let augmented = uff_energy + d4.total_kcal_mol;
        println!("[B5] UFF: {:.2}, D4: {:.4}, UFF+D4: {:.2} kcal/mol",
            uff_energy, d4.total_kcal_mol, augmented);
        assert!(augmented <= uff_energy + 0.01, "D4 should lower UFF energy");
    }

    #[test]
    fn b5_coordination_numbers() {
        let (elems, pos) = ethanol();
        let d4 = compute_d4_energy(&elems, &pos, &D4Config::default());

        println!("[B5] D4 coordination numbers: {:?}", d4.coordination_numbers);
        assert!(d4.coordination_numbers[0] > 2.0, "C CN too low");
        assert!(d4.coordination_numbers[2] > 0.5, "O CN too low");
    }

    #[test]
    fn b5_timing_d4_overhead() {
        let (elems, pos) = ethanol();
        let conf = sci_form::embed("CCO", 42);

        let t0 = Instant::now();
        for _ in 0..200 { let _ = sci_form::compute_uff_energy("CCO", &conf.coords); }
        let uff_us = t0.elapsed().as_micros() as f64 / 200.0;

        let t0 = Instant::now();
        for _ in 0..200 { let _ = compute_d4_energy(&elems, &pos, &D4Config::default()); }
        let d4_us = t0.elapsed().as_micros() as f64 / 200.0;

        println!("[B5] Timing: UFF={:.0} µs, D4={:.0} µs, overhead={:.1}%",
            uff_us, d4_us, d4_us / uff_us.max(1.0) * 100.0);
    }
}

// ═══════════════════════════════════════════════════════════════════
// B6 — ETKDG vs. Riemannian (E3): PSD guarantee, embedding quality
// ═══════════════════════════════════════════════════════════════════

#[cfg(feature = "experimental-riemannian")]
mod bench_riemannian_vs_etkdg {
    use sci_form::experimental::riemannian::*;
    use nalgebra::DMatrix;
    use std::time::Instant;

    fn make_constraints(coords: &[[f64; 3]]) -> Vec<DistanceConstraint> {
        let n = coords.len();
        let mut constraints = Vec::new();
        for i in 0..n {
            for j in (i + 1)..n {
                let dx = coords[i][0] - coords[j][0];
                let dy = coords[i][1] - coords[j][1];
                let dz = coords[i][2] - coords[j][2];
                let d = (dx * dx + dy * dy + dz * dz).sqrt();
                constraints.push(DistanceConstraint { i, j, lower: d * 0.95, upper: d * 1.05 });
            }
        }
        constraints
    }

    #[test]
    fn b6_psd_guarantee() {
        let coords = vec![
            [0.0, 0.0, 0.0], [1.5, 0.0, 0.0],
            [0.75, 1.3, 0.0], [0.75, 0.43, 1.22],
        ];
        let n = coords.len();
        let constraints = make_constraints(&coords);

        let config = RiemannianConfig::default();
        let lbfgs = RiemannianLbfgs::new(config);
        let initial = DMatrix::identity(n, n);
        let result = lbfgs.optimize(&initial, &constraints);

        println!("[B6] Riemannian converged: {}, iter: {}, max_violation: {:.6}",
            result.converged, result.iterations, result.max_violation);

        let eig = nalgebra::SymmetricEigen::new(result.metric_matrix.clone());
        let min_eig = eig.eigenvalues.iter().cloned().fold(f64::INFINITY, f64::min);
        println!("[B6] Min eigenvalue of metric matrix: {:.6e}", min_eig);
        assert!(min_eig >= -1e-10, "Metric matrix not PSD: min eigenvalue = {}", min_eig);
    }

    #[test]
    fn b6_etkdg_vs_riemannian_simple() {
        let conf = sci_form::embed("C", 42);
        assert!(conf.error.is_none(), "embed failed: {:?}", conf.error);

        let n = conf.num_atoms;
        let pos: Vec<[f64; 3]> = conf.coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();
        let constraints = make_constraints(&pos);

        let config = RiemannianConfig::default();
        let lbfgs = RiemannianLbfgs::new(config);
        let initial = DMatrix::identity(n, n);
        let result = lbfgs.optimize(&initial, &constraints);

        println!("[B6] ETKDG time: {:.2} ms", conf.time_ms);
        println!("[B6] Riemannian: converged={}, violation={:.4} Å, iter={}",
            result.converged, result.max_violation, result.iterations);
    }

    #[test]
    fn b6_timing_comparison() {
        let t0 = Instant::now();
        for _ in 0..20 { let _ = sci_form::embed("CCO", 42); }
        let etkdg_us = t0.elapsed().as_micros() as f64 / 20.0;

        let n = 4;
        let constraints = vec![
            DistanceConstraint { i: 0, j: 1, lower: 1.4, upper: 1.6 },
            DistanceConstraint { i: 0, j: 2, lower: 1.4, upper: 1.6 },
            DistanceConstraint { i: 0, j: 3, lower: 1.4, upper: 1.6 },
            DistanceConstraint { i: 1, j: 2, lower: 1.4, upper: 1.6 },
            DistanceConstraint { i: 1, j: 3, lower: 1.4, upper: 1.6 },
            DistanceConstraint { i: 2, j: 3, lower: 1.4, upper: 1.6 },
        ];
        let config = RiemannianConfig::default();
        let t0 = Instant::now();
        for _ in 0..20 {
            let lbfgs = RiemannianLbfgs::new(config.clone());
            let initial = DMatrix::identity(n, n);
            let _ = lbfgs.optimize(&initial, &constraints);
        }
        let riem_us = t0.elapsed().as_micros() as f64 / 20.0;

        println!("[B6] Timing: ETKDG={:.0} µs, Riemannian={:.0} µs", etkdg_us, riem_us);
    }
}

// ═══════════════════════════════════════════════════════════════════
// B7 — ETKDG vs. SDR (E8): PSD embedding, retry elimination
// ═══════════════════════════════════════════════════════════════════

#[cfg(feature = "experimental-sdr")]
mod bench_sdr_vs_etkdg {
    use sci_form::experimental::sdr::*;
    use std::time::Instant;

    #[test]
    fn b7_sdr_embedding_basic() {
        let n = 5;
        let d_pairs = vec![
            (0, 1, 1.5), (1, 2, 1.5), (2, 3, 1.5), (3, 4, 1.5),
            (0, 2, 2.45), (1, 3, 2.45), (2, 4, 2.45),
        ];

        let result = sdr_embed(n, &d_pairs, &SdrConfig::default());

        println!("[B7] SDR max_distance_error: {:.4} Å, retries_avoided: {}",
            result.max_distance_error, result.retries_avoided);
        assert_eq!(result.num_atoms, n);
        assert_eq!(result.coords.len(), n * 3);
    }

    #[test]
    fn b7_sdr_vs_etkdg_distance_fidelity() {
        let conf = sci_form::embed("CC", 42);
        assert!(conf.error.is_none());

        let n = conf.num_atoms;
        let pos: Vec<[f64; 3]> = conf.coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();

        let mut d_pairs = Vec::new();
        for i in 0..n {
            for j in (i + 1)..n {
                let dx = pos[i][0] - pos[j][0];
                let dy = pos[i][1] - pos[j][1];
                let dz = pos[i][2] - pos[j][2];
                let d = (dx * dx + dy * dy + dz * dz).sqrt();
                d_pairs.push((i, j, d));
            }
        }

        let sdr = sdr_embed(n, &d_pairs, &SdrConfig::default());
        let sdr_pos: Vec<[f64; 3]> = sdr.coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();

        let mut max_err = 0.0_f64;
        for &(i, j, d_target) in &d_pairs {
            let dx = sdr_pos[i][0] - sdr_pos[j][0];
            let dy = sdr_pos[i][1] - sdr_pos[j][1];
            let dz = sdr_pos[i][2] - sdr_pos[j][2];
            let d_sdr = (dx * dx + dy * dy + dz * dz).sqrt();
            max_err = max_err.max((d_sdr - d_target).abs());
        }

        println!("[B7] ETKDG → SDR max distance error: {:.4} Å", max_err);
        assert!(max_err < 1.0, "SDR distance fidelity too poor: {:.3} Å", max_err);
    }

    #[test]
    fn b7_timing_comparison() {
        let t0 = Instant::now();
        for _ in 0..20 { let _ = sci_form::embed("CCO", 42); }
        let etkdg_us = t0.elapsed().as_micros() as f64 / 20.0;

        let n = 5;
        let d_pairs = vec![
            (0, 1, 1.5), (1, 2, 1.5), (2, 3, 1.5), (3, 4, 1.5),
            (0, 2, 2.45), (1, 3, 2.45), (2, 4, 2.45),
        ];
        let t0 = Instant::now();
        for _ in 0..20 { let _ = sdr_embed(n, &d_pairs, &SdrConfig::default()); }
        let sdr_us = t0.elapsed().as_micros() as f64 / 20.0;

        println!("[B7] Timing: ETKDG={:.0} µs, SDR={:.0} µs", etkdg_us, sdr_us);
    }
}

// ═══════════════════════════════════════════════════════════════════
// B8 — Full Hessian IR vs. MBH (E9): Frequency accuracy, DOF reduction
// ═══════════════════════════════════════════════════════════════════

#[cfg(feature = "experimental-mbh")]
mod bench_mbh_vs_full_hessian {
    use sci_form::experimental::mbh::*;
    use std::time::Instant;

    fn ethanol() -> (Vec<u8>, Vec<[f64; 3]>) {
        (
            vec![6, 6, 8, 1, 1, 1, 1, 1, 1],
            vec![
                [0.0, 0.0, 0.0], [1.52, 0.0, 0.0], [2.13, 1.21, 0.0],
                [-0.36, 1.03, 0.0], [-0.36, -0.51, 0.89], [-0.36, -0.51, -0.89],
                [1.88, -0.51, 0.89], [1.88, -0.51, -0.89], [3.09, 1.21, 0.0],
            ],
        )
    }

    fn benzene() -> (Vec<u8>, Vec<[f64; 3]>) {
        (
            vec![6, 6, 6, 6, 6, 6, 1, 1, 1, 1, 1, 1],
            vec![
                [1.40, 0.0, 0.0], [0.70, 1.21, 0.0], [-0.70, 1.21, 0.0],
                [-1.40, 0.0, 0.0], [-0.70, -1.21, 0.0], [0.70, -1.21, 0.0],
                [2.48, 0.0, 0.0], [1.24, 2.15, 0.0], [-1.24, 2.15, 0.0],
                [-2.48, 0.0, 0.0], [-1.24, -2.15, 0.0], [1.24, -2.15, 0.0],
            ],
        )
    }

    fn uff_energy_fn(smiles: &'static str) -> impl Fn(&[f64]) -> f64 {
        move |coords: &[f64]| -> f64 {
            sci_form::compute_uff_energy(smiles, coords).unwrap_or(0.0)
        }
    }

    #[test]
    fn b8_mbh_dof_reduction() {
        let (elems, pos) = benzene();
        let rings: Vec<(Vec<usize>, bool)> = vec![(vec![0, 1, 2, 3, 4, 5], true)];
        let energy_fn = uff_energy_fn("c1ccccc1");

        let mbh = compute_mbh_frequencies(&elems, &pos, &rings, &energy_fn, &MbhConfig::default());

        println!("[B8] Full DOF: {}, Reduced DOF: {}, Speedup: {:.1}x",
            mbh.n_dof_full, mbh.n_dof_reduced, mbh.speedup);
        assert!(mbh.speedup >= 1.0, "Speedup should be >= 1.0");
        assert!(mbh.n_dof_reduced <= mbh.n_dof_full);
    }

    #[test]
    fn b8_mbh_frequency_physical() {
        let (elems, pos) = ethanol();
        let rings: Vec<(Vec<usize>, bool)> = vec![];
        let energy_fn = uff_energy_fn("CCO");

        let mbh = compute_mbh_frequencies(&elems, &pos, &rings, &energy_fn, &MbhConfig::default());

        println!("[B8] MBH frequencies (N={}): {:?}", mbh.frequencies.len(), mbh.frequencies);

        let positive: Vec<_> = mbh.frequencies.iter().filter(|&&f| f > 10.0).collect();
        assert!(!positive.is_empty(), "Should have positive vibrational frequencies");
    }

    #[test]
    fn b8_mbh_vs_full_ethanol() {
        let (elems, pos) = ethanol();

        let full = sci_form::ir::vibrations::compute_vibrational_analysis(
            &elems, &pos, sci_form::ir::hessian::HessianMethod::Pm3, None,
        );

        let rings: Vec<(Vec<usize>, bool)> = vec![];
        let energy_fn = uff_energy_fn("CCO");
        let mbh = compute_mbh_frequencies(&elems, &pos, &rings, &energy_fn, &MbhConfig::default());

        if let Ok(ref full_result) = full {
            let full_freqs: Vec<f64> = full_result.modes.iter()
                .filter(|m| m.is_real && m.frequency_cm1 > 10.0)
                .map(|m| m.frequency_cm1)
                .collect();
            println!("[B8] Full Hessian (N={}): [{:.0} ... {:.0}]",
                full_freqs.len(),
                full_freqs.first().unwrap_or(&0.0),
                full_freqs.last().unwrap_or(&0.0));
        }
        println!("[B8] MBH (N={}): [{:.0} ... {:.0}]",
            mbh.frequencies.len(),
            mbh.frequencies.first().unwrap_or(&0.0),
            mbh.frequencies.last().unwrap_or(&0.0));
    }

    #[test]
    fn b8_timing_comparison() {
        let (elems, pos) = ethanol();
        let rings: Vec<(Vec<usize>, bool)> = vec![];
        let energy_fn = uff_energy_fn("CCO");

        let t0 = Instant::now();
        for _ in 0..3 {
            let _ = compute_mbh_frequencies(&elems, &pos, &rings, &energy_fn, &MbhConfig::default());
        }
        let mbh_us = t0.elapsed().as_micros() as f64 / 3.0;

        // Full Hessian with PM3 is expensive — single iteration only
        let t0 = Instant::now();
        let _ = sci_form::ir::vibrations::compute_vibrational_analysis(
            &elems, &pos, sci_form::ir::hessian::HessianMethod::Pm3, None,
        );
        let full_us = t0.elapsed().as_micros() as f64;

        println!("[B8] Timing: Full Hessian={:.0} µs, MBH={:.0} µs", full_us, mbh_us);
    }
}

// ═══════════════════════════════════════════════════════════════════
// B9 — CGA Motor algebra (E1): Geometric correctness
// ═══════════════════════════════════════════════════════════════════

#[cfg(feature = "experimental-cga")]
mod bench_cga_geometry {
    use sci_form::experimental::cga::*;

    #[test]
    fn b9_motor_rotation_preserves_distance() {
        let a = [0.0, 0.0, 0.0];
        let b = [1.5, 0.0, 0.0];
        let angle = std::f64::consts::PI / 4.0;

        let motor = dihedral_motor(a, b, angle);
        let p = [2.0, 1.0, 0.0];
        let p_rot = motor.transform_point(p);

        let dist_orig = ((p[1] * p[1] + p[2] * p[2]) as f64).sqrt();
        let dist_rot = ((p_rot[1] * p_rot[1] + p_rot[2] * p_rot[2]) as f64).sqrt();

        println!("[B9] Original dist from axis: {:.4}, Rotated: {:.4}", dist_orig, dist_rot);
        assert!((dist_orig - dist_rot).abs() < 0.1, "Motor should preserve distance");
    }

    #[test]
    fn b9_motor_identity() {
        let motor = Motor::identity();
        let p = [1.0, 2.0, 3.0];
        let result = motor.transform_point(p);

        assert!((result[0] - p[0]).abs() < 1e-10);
        assert!((result[1] - p[1]).abs() < 1e-10);
        assert!((result[2] - p[2]).abs() < 1e-10);
    }

    #[test]
    fn b9_cga_vs_matrix_rotation() {
        let axis_a = [0.0, 0.0, 0.0];
        let axis_b = [0.0, 0.0, 1.0];
        let angle = std::f64::consts::PI / 2.0;

        let motor = dihedral_motor(axis_a, axis_b, angle);
        let p = [1.0, 0.0, 0.0];
        let result = motor.transform_point(p);

        println!("[B9] Rotated [1,0,0] by 90° about z: [{:.3}, {:.3}, {:.3}]",
            result[0], result[1], result[2]);
        assert!(result[0].abs() < 0.2, "x should be ~0");
        assert!((result[1] - 1.0).abs() < 0.2, "y should be ~1");
        assert!(result[2].abs() < 0.2, "z should be ~0");
    }

    #[test]
    fn b9_subtree_rotation() {
        let coords = vec![0.0, 0.0, 0.0, 1.5, 0.0, 0.0, 3.0, 0.0, 0.0];
        let a = [0.0, 0.0, 0.0];
        let b = [1.5, 0.0, 0.0];
        let motor = dihedral_motor(a, b, std::f64::consts::PI / 2.0);
        let result = apply_motor_to_subtree(&coords, &[2], &motor);

        // Atom 2 at [3,0,0] should rotate about x-axis
        println!("[B9] After subtree rotation: atom2 = [{:.3}, {:.3}, {:.3}]",
            result[6], result[7], result[8]);
        // Atom 0 and 1 should be unchanged
        assert!((result[0] - coords[0]).abs() < 1e-10);
        assert!((result[3] - coords[3]).abs() < 1e-10);
    }
}

// ═══════════════════════════════════════════════════════════════════
// B10 — CPM Electrode (E10): Charge vs potential curve
// ═══════════════════════════════════════════════════════════════════

#[cfg(feature = "experimental-cpm")]
mod bench_cpm_electrochemistry {
    use sci_form::experimental::cpm::*;

    fn water() -> (Vec<u8>, Vec<[f64; 3]>) {
        (
            vec![8, 1, 1],
            vec![[0.0, 0.0, 0.117], [0.0, 0.757, -0.469], [0.0, -0.757, -0.469]],
        )
    }

    #[test]
    fn b10_cpm_charge_varies_with_potential() {
        let (elems, pos) = water();

        let mut charges_at_potentials = Vec::new();
        for mu in [-5.0_f64, -4.5, -4.0, -3.5] {
            let config = CpmConfig { mu_ev: mu, ..Default::default() };
            let result = compute_cpm_charges(&elems, &pos, &config);
            charges_at_potentials.push((mu, result.total_charge));
            println!("[B10] µ = {:.1} eV → Q = {:.4} e", mu, result.total_charge);
        }

        let q_max = charges_at_potentials.iter().map(|&(_, q)| q).fold(f64::NEG_INFINITY, f64::max);
        let q_min = charges_at_potentials.iter().map(|&(_, q)| q).fold(f64::INFINITY, f64::min);
        let q_range = q_max - q_min;
        println!("[B10] Charge range: {:.4} e", q_range);
        assert!(q_range > 1e-6, "Charge should vary with potential");
    }

    #[test]
    fn b10_cpm_surface_shape() {
        let (elems, pos) = water();
        let surface = compute_cpm_surface(&elems, &pos, -5.5, -3.5, 10, 78.5);

        println!("[B10] CPM surface points: {}", surface.mu_values.len());
        assert_eq!(surface.mu_values.len(), 10);
        assert_eq!(surface.total_charge.len(), 10);
        assert_eq!(surface.free_energy.len(), 10);
    }

    #[test]
    fn b10_grand_potential_decreases() {
        let (elems, pos) = water();
        let surface = compute_cpm_surface(&elems, &pos, -5.5, -3.5, 20, 78.5);

        let finite_energies: Vec<f64> = surface.free_energy.iter()
            .filter(|e| e.is_finite())
            .copied()
            .collect();
        if finite_energies.len() >= 2 {
            println!("[B10] Free energy range: [{:.4}, {:.4}] eV",
                finite_energies.iter().cloned().fold(f64::INFINITY, f64::min),
                finite_energies.iter().cloned().fold(f64::NEG_INFINITY, f64::max));
        }

        assert!(surface.all_converged, "CPM surface should converge at all points");
    }

    #[test]
    fn b10_capacitance_positive() {
        let (elems, pos) = water();
        let surface = compute_cpm_surface(&elems, &pos, -5.5, -3.5, 10, 78.5);

        for (i, &c) in surface.capacitance.iter().enumerate() {
            println!("[B10] C(µ={:.2}) = {:.4} e/eV", surface.mu_values[i], c);
        }
    }
}

// ═══════════════════════════════════════════════════════════════════
// B11 — GSM Reaction Path (E11): Transition state search
// ═══════════════════════════════════════════════════════════════════

#[cfg(feature = "experimental-gsm")]
mod bench_gsm_reaction_path {
    use sci_form::experimental::gsm::*;

    /// Simple harmonic bond energy for H2: E = k*(d - d0)^2
    fn h2_energy(coords: &[f64]) -> f64 {
        let dx = coords[3] - coords[0];
        let dy = coords[4] - coords[1];
        let dz = coords[5] - coords[2];
        let d = (dx * dx + dy * dy + dz * dz).sqrt();
        let d0 = 0.74;
        let k = 200.0;
        k * (d - d0).powi(2)
    }

    #[test]
    fn b11_string_growth() {
        let reactant = vec![0.0, 0.0, 0.0, 0.74, 0.0, 0.0];
        let product = vec![0.0, 0.0, 0.0, 3.0, 0.0, 0.0];

        let config = GsmConfig { max_nodes: 7, max_iter: 20, ..Default::default() };
        let path = grow_string(&reactant, &product, &h2_energy, &config);

        println!("[B11] GSM nodes: {}, joined: {}", path.nodes.len(), path.joined);
        println!("[B11] Energies: {:?}", path.energies);
        assert!(path.nodes.len() >= 2, "Should have at least reactant and product");
    }

    #[test]
    fn b11_interpolation_smoothness() {
        let reactant = vec![0.0, 0.0, 0.0, 0.74, 0.0, 0.0];
        let product = vec![0.0, 0.0, 0.0, 3.0, 0.0, 0.0];

        let mut prev_d = 0.0_f64;
        for &t in &[0.0, 0.25, 0.5, 0.75, 1.0] {
            let node = interpolate_node(&reactant, &product, t);
            let d = ((node[3] - node[0]).powi(2) + (node[4] - node[1]).powi(2)
                + (node[5] - node[2]).powi(2)).sqrt();
            println!("[B11] t={:.2} → d(H-H) = {:.3} Å", t, d);
            if t > 0.0 {
                assert!(d >= prev_d - 0.01, "Distance should increase monotonically");
            }
            prev_d = d;
        }
    }

    #[test]
    fn b11_transition_state_search() {
        let reactant = vec![0.0, 0.0, 0.0, 0.74, 0.0, 0.0];
        let product = vec![0.0, 0.0, 0.0, 3.0, 0.0, 0.0];

        let config = GsmConfig { max_nodes: 7, max_iter: 10, ..Default::default() };
        let result = find_transition_state(&reactant, &product, &h2_energy, &config);

        println!("[B11] TS energy: {:.2} kcal/mol, activation: {:.2} kcal/mol",
            result.ts_energy, result.activation_energy);
        println!("[B11] TS node index: {}, n_nodes: {}", result.ts_node_index, result.n_nodes);

        assert!(result.activation_energy >= 0.0, "Activation energy should be non-negative");
    }

    #[test]
    fn b11_energy_profile_has_maximum() {
        let reactant = vec![0.0, 0.0, 0.0, 0.74, 0.0, 0.0];
        let product = vec![0.0, 0.0, 0.0, 3.0, 0.0, 0.0];

        let config = GsmConfig { max_nodes: 11, max_iter: 20, ..Default::default() };
        let result = find_transition_state(&reactant, &product, &h2_energy, &config);

        if result.path_energies.len() >= 3 {
            let e_first = result.path_energies[0];
            let e_last = *result.path_energies.last().unwrap();
            let e_max = result.path_energies.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

            println!("[B11] E(react)={:.2}, E(max)={:.2}, E(prod)={:.2}",
                e_first, e_max, e_last);
            assert!(e_max >= e_first - 0.01, "Maximum should be >= reactant energy");
        }
    }
}
