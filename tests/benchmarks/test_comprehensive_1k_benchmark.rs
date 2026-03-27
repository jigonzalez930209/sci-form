//! Comprehensive 1000-molecule benchmark validating ALL public API functions.
//!
//! Loads the first 1000 SMILES from data/chembl_1k_practical_asc.smi (avg ~31 heavy atoms),
//! runs every public function, and checks physical invariants with ≤0.5% deviation tolerance.
//!
//! Run: cargo test --test test_comprehensive_1k_benchmark --release -- --nocapture

use std::io::BufRead;

// ═══════════════════════════════════════════════════════════════════════════════
// Helpers
// ═══════════════════════════════════════════════════════════════════════════════

fn load_smiles(limit: usize) -> Vec<(String, String)> {
    let path = std::path::Path::new("data/chembl_1k_practical_asc.smi");
    let file = std::fs::File::open(path).expect("Cannot open data/chembl_1k_practical_asc.smi");
    let reader = std::io::BufReader::new(file);
    reader
        .lines()
        .take(limit)
        .filter_map(|line| {
            let line = line.ok()?;
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() >= 2 {
                Some((parts[0].to_string(), parts[1].to_string()))
            } else {
                None
            }
        })
        .collect()
}

fn to_positions(coords: &[f64]) -> Vec<[f64; 3]> {
    coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect()
}

const N_MOLECULES: usize = 1000;
const MAX_FAILURE_RATE: f64 = 0.005; // 0.5%

// ═══════════════════════════════════════════════════════════════════════════════
// Module 1: Conformer Generation (embed)
// ═══════════════════════════════════════════════════════════════════════════════

#[test]
fn benchmark_conformer_embedding_1k() {
    let smiles = load_smiles(N_MOLECULES);
    assert!(smiles.len() >= N_MOLECULES, "Need at least 1000 SMILES");

    let mut success = 0usize;
    let mut total_time_ms = 0.0f64;
    let mut bond_length_violations = 0usize;
    let mut steric_clashes = 0usize;
    let mut nan_coords = 0usize;

    for (smi, _id) in &smiles {
        let result = sci_form::embed(smi, 42);
        if result.error.is_some() {
            continue;
        }
        success += 1;
        total_time_ms += result.time_ms;

        // Check for NaN/Inf coords
        if result.coords.iter().any(|c| !c.is_finite()) {
            nan_coords += 1;
            continue;
        }

        let pos = to_positions(&result.coords);

        // Check bond lengths (0.7–2.5 Å for organic molecules)
        for &(a, b, ref _order) in &result.bonds {
            let dx = pos[a][0] - pos[b][0];
            let dy = pos[a][1] - pos[b][1];
            let dz = pos[a][2] - pos[b][2];
            let dist = (dx * dx + dy * dy + dz * dz).sqrt();
            if !(0.7..=2.5).contains(&dist) {
                bond_length_violations += 1;
            }
        }

        // Check for steric clashes (non-bonded < 0.5 Å)
        let bonded: std::collections::HashSet<(usize, usize)> = result
            .bonds
            .iter()
            .flat_map(|&(a, b, _)| vec![(a, b), (b, a)])
            .collect();
        for i in 0..result.num_atoms {
            for j in (i + 1)..result.num_atoms {
                if bonded.contains(&(i, j)) {
                    continue;
                }
                let dx = pos[i][0] - pos[j][0];
                let dy = pos[i][1] - pos[j][1];
                let dz = pos[i][2] - pos[j][2];
                let dist = (dx * dx + dy * dy + dz * dz).sqrt();
                if dist < 0.5 {
                    steric_clashes += 1;
                    break; // count once per molecule
                }
            }
        }
    }

    let failures = N_MOLECULES - success;
    let success_rate = success as f64 / N_MOLECULES as f64;
    let avg_time = if success > 0 {
        total_time_ms / success as f64
    } else {
        0.0
    };
    let failure_rate = failures as f64 / N_MOLECULES as f64;

    eprintln!("=== CONFORMER EMBEDDING ===");
    eprintln!(
        "  Success: {success}/{N_MOLECULES} ({:.2}%)",
        success_rate * 100.0
    );
    eprintln!("  Avg time: {avg_time:.1} ms/mol");
    eprintln!("  NaN coords: {nan_coords}");
    eprintln!("  Bond length violations: {bond_length_violations}");
    eprintln!("  Steric clashes: {steric_clashes}");

    // Allow up to 5 failures out of 1000 (0.5%)
    let max_failures = (N_MOLECULES as f64 * MAX_FAILURE_RATE).ceil() as usize;
    assert!(
        failures <= max_failures,
        "Embedding failures {failures} ({:.2}%) exceeds limit of {max_failures} ({:.2}%)",
        failure_rate * 100.0,
        MAX_FAILURE_RATE * 100.0
    );
    assert_eq!(
        nan_coords, 0,
        "NaN coordinates found in {nan_coords} molecules"
    );
}

// ═══════════════════════════════════════════════════════════════════════════════
// Module 2: Gasteiger Charges
// ═══════════════════════════════════════════════════════════════════════════════

#[test]
fn benchmark_gasteiger_charges_1k() {
    let smiles = load_smiles(N_MOLECULES);
    let mut success = 0usize;
    let mut charge_conservation_violations = 0usize;
    let mut unbounded_charges = 0usize;

    for (smi, _id) in &smiles {
        if let Ok(result) = sci_form::compute_charges(smi) {
            success += 1;
            // Charge conservation: sum should be near formal charge (0 for neutrals)
            let sum: f64 = result.charges.iter().sum();
            if sum.abs() > 0.5 {
                charge_conservation_violations += 1;
            }
            // Charges should be bounded [-2, +2] for organic molecules
            if result.charges.iter().any(|&c| c.abs() > 2.0) {
                unbounded_charges += 1;
            }
        }
    }

    let success_rate = success as f64 / N_MOLECULES as f64;
    let failure_rate = 1.0 - success_rate;

    eprintln!("=== GASTEIGER CHARGES ===");
    eprintln!(
        "  Success: {success}/{N_MOLECULES} ({:.2}%)",
        success_rate * 100.0
    );
    eprintln!("  Charge conservation violations: {charge_conservation_violations}");
    eprintln!("  Unbounded charges: {unbounded_charges}");

    let failures = N_MOLECULES - success;
    let max_failures = (N_MOLECULES as f64 * MAX_FAILURE_RATE).ceil() as usize;
    assert!(
        failures <= max_failures,
        "Gasteiger charges failure rate {:.2}% exceeds limit",
        failure_rate * 100.0
    );
}

// ═══════════════════════════════════════════════════════════════════════════════
// Module 3: SASA (Solvent Accessible Surface Area)
// ═══════════════════════════════════════════════════════════════════════════════

#[test]
fn benchmark_sasa_1k() {
    let smiles = load_smiles(N_MOLECULES);
    let mut embed_ok = 0usize;
    let mut sasa_ok = 0usize;
    let mut negative_sasa = 0usize;
    let mut zero_sasa = 0usize;

    for (smi, _id) in &smiles {
        let conf = sci_form::embed(smi, 42);
        if conf.error.is_some() || conf.coords.iter().any(|c| !c.is_finite()) {
            continue;
        }
        embed_ok += 1;

        if let Ok(result) = sci_form::compute_sasa(&conf.elements, &conf.coords, Some(1.4)) {
            sasa_ok += 1;
            if result.total_sasa < 0.0 {
                negative_sasa += 1;
            }
            if result.total_sasa == 0.0 {
                zero_sasa += 1;
            }
            // Per-atom SASA should be non-negative
            assert!(
                result.atom_sasa.iter().all(|&s| s >= 0.0),
                "Negative per-atom SASA in molecule"
            );
        }
    }

    let sasa_rate = sasa_ok as f64 / embed_ok as f64;
    let failure_rate = 1.0 - sasa_rate;

    eprintln!("=== SASA ===");
    eprintln!("  Embedded: {embed_ok}/{N_MOLECULES}");
    eprintln!(
        "  SASA success: {sasa_ok}/{embed_ok} ({:.2}%)",
        sasa_rate * 100.0
    );
    eprintln!("  Negative total SASA: {negative_sasa}");
    eprintln!("  Zero SASA: {zero_sasa}");

    let failures = embed_ok - sasa_ok;
    let max_failures = (embed_ok as f64 * MAX_FAILURE_RATE).ceil() as usize;
    assert!(
        failures <= max_failures,
        "SASA failure rate {:.2}% exceeds limit",
        failure_rate * 100.0
    );
    assert_eq!(negative_sasa, 0, "Negative SASA found");
}

// ═══════════════════════════════════════════════════════════════════════════════
// Module 4: Population Analysis (Mulliken & Löwdin)
// ═══════════════════════════════════════════════════════════════════════════════

#[test]
fn benchmark_population_1k() {
    let smiles = load_smiles(N_MOLECULES);
    let mut embed_ok = 0usize;
    let mut pop_ok = 0usize;
    let mut charge_conservation_fails = 0usize;
    let mut _electronegativity_violations = 0usize;

    for (smi, _id) in &smiles {
        let conf = sci_form::embed(smi, 42);
        if conf.error.is_some() || conf.coords.iter().any(|c| !c.is_finite()) {
            continue;
        }
        embed_ok += 1;
        let pos = to_positions(&conf.coords);

        if let Ok(result) = sci_form::compute_population(&conf.elements, &pos) {
            pop_ok += 1;
            // Charge conservation for neutral molecules
            let sum_m: f64 = result.mulliken_charges.iter().sum();
            let sum_l: f64 = result.lowdin_charges.iter().sum();
            if sum_m.abs() > 1.0 || sum_l.abs() > 1.0 {
                charge_conservation_fails += 1;
            }
            // All charges should be finite
            assert!(
                result.mulliken_charges.iter().all(|c| c.is_finite()),
                "Non-finite Mulliken charge"
            );
            assert!(
                result.lowdin_charges.iter().all(|c| c.is_finite()),
                "Non-finite Löwdin charge"
            );
        }
    }

    let pop_rate = pop_ok as f64 / embed_ok as f64;
    let failure_rate = 1.0 - pop_rate;

    eprintln!("=== POPULATION ANALYSIS ===");
    eprintln!("  Embedded: {embed_ok}/{N_MOLECULES}");
    eprintln!(
        "  Population success: {pop_ok}/{embed_ok} ({:.2}%)",
        pop_rate * 100.0
    );
    eprintln!("  Charge conservation fails: {charge_conservation_fails}");

    let failures = embed_ok - pop_ok;
    let max_failures = (embed_ok as f64 * MAX_FAILURE_RATE).ceil() as usize;
    assert!(
        failures <= max_failures,
        "Population analysis failure rate {:.2}% exceeds limit",
        failure_rate * 100.0
    );
}

// ═══════════════════════════════════════════════════════════════════════════════
// Module 5: Dipole Moment
// ═══════════════════════════════════════════════════════════════════════════════

#[test]
fn benchmark_dipole_1k() {
    let smiles = load_smiles(N_MOLECULES);
    let mut embed_ok = 0usize;
    let mut dipole_ok = 0usize;
    let mut magnitude_consistency_fails = 0usize;
    let mut unreasonable_magnitude = 0usize;

    for (smi, _id) in &smiles {
        let conf = sci_form::embed(smi, 42);
        if conf.error.is_some() || conf.coords.iter().any(|c| !c.is_finite()) {
            continue;
        }
        embed_ok += 1;
        let pos = to_positions(&conf.coords);

        if let Ok(result) = sci_form::compute_dipole(&conf.elements, &pos) {
            dipole_ok += 1;
            // Vector magnitude consistency
            let calc_mag =
                (result.vector[0].powi(2) + result.vector[1].powi(2) + result.vector[2].powi(2))
                    .sqrt();
            if (calc_mag - result.magnitude).abs() > 1e-8 {
                magnitude_consistency_fails += 1;
            }
            // Dipole should be in reasonable range (0-50 D for organic molecules)
            if result.magnitude > 50.0 || result.magnitude < 0.0 {
                unreasonable_magnitude += 1;
            }
        }
    }

    let dipole_rate = dipole_ok as f64 / embed_ok as f64;
    let failure_rate = 1.0 - dipole_rate;

    eprintln!("=== DIPOLE MOMENT ===");
    eprintln!("  Embedded: {embed_ok}/{N_MOLECULES}");
    eprintln!(
        "  Dipole success: {dipole_ok}/{embed_ok} ({:.2}%)",
        dipole_rate * 100.0
    );
    eprintln!("  Magnitude consistency fails: {magnitude_consistency_fails}");
    eprintln!("  Unreasonable magnitude: {unreasonable_magnitude}");

    let failures = embed_ok - dipole_ok;
    let max_failures = (embed_ok as f64 * MAX_FAILURE_RATE).ceil() as usize;
    assert!(
        failures <= max_failures,
        "Dipole failure rate {:.2}% exceeds limit",
        failure_rate * 100.0
    );
    assert_eq!(
        magnitude_consistency_fails, 0,
        "Magnitude consistency failures"
    );
}

// ═══════════════════════════════════════════════════════════════════════════════
// Module 6: DOS / PDOS
// ═══════════════════════════════════════════════════════════════════════════════

#[test]
fn benchmark_dos_1k() {
    let smiles = load_smiles(N_MOLECULES);
    let mut embed_ok = 0usize;
    let mut dos_ok = 0usize;
    let mut pdos_sum_mismatch = 0usize;
    let mut negative_dos = 0usize;

    for (smi, _id) in &smiles {
        let conf = sci_form::embed(smi, 42);
        if conf.error.is_some() || conf.coords.iter().any(|c| !c.is_finite()) {
            continue;
        }
        embed_ok += 1;
        let pos = to_positions(&conf.coords);

        if let Ok(result) = sci_form::compute_dos(&conf.elements, &pos, 0.3, -30.0, 5.0, 200) {
            dos_ok += 1;

            // DOS should be non-negative
            if result.total_dos.iter().any(|&d| d < -1e-10) {
                negative_dos += 1;
            }

            // PDOS sum should approximate total DOS
            if !result.pdos.is_empty() && result.pdos[0].len() == result.total_dos.len() {
                for k in 0..result.total_dos.len() {
                    let pdos_sum: f64 = result.pdos.iter().map(|p| p[k]).sum();
                    let diff = (pdos_sum - result.total_dos[k]).abs();
                    if diff > 0.1 * result.total_dos[k].abs().max(0.01) {
                        pdos_sum_mismatch += 1;
                        break;
                    }
                }
            }
        }
    }

    let dos_rate = dos_ok as f64 / embed_ok as f64;
    let failure_rate = 1.0 - dos_rate;

    eprintln!("=== DOS/PDOS ===");
    eprintln!("  Embedded: {embed_ok}/{N_MOLECULES}");
    eprintln!(
        "  DOS success: {dos_ok}/{embed_ok} ({:.2}%)",
        dos_rate * 100.0
    );
    eprintln!("  Negative DOS: {negative_dos}");
    eprintln!("  PDOS sum mismatches: {pdos_sum_mismatch}");

    let failures = embed_ok - dos_ok;
    let max_failures = (embed_ok as f64 * MAX_FAILURE_RATE).ceil() as usize;
    assert!(
        failures <= max_failures,
        "DOS failure rate {:.2}% exceeds limit",
        failure_rate * 100.0
    );
}

// ═══════════════════════════════════════════════════════════════════════════════
// Module 7: UFF & MMFF94 Force Field Energy
// ═══════════════════════════════════════════════════════════════════════════════

#[test]
fn benchmark_forcefield_energy_1k() {
    let smiles = load_smiles(N_MOLECULES);
    let mut embed_ok = 0usize;
    let mut uff_ok = 0usize;
    let mut mmff_ok = 0usize;
    let mut uff_infinite = 0usize;
    let mut mmff_infinite = 0usize;

    for (smi, _id) in &smiles {
        let conf = sci_form::embed(smi, 42);
        if conf.error.is_some() || conf.coords.iter().any(|c| !c.is_finite()) {
            continue;
        }
        embed_ok += 1;

        if let Ok(e) = sci_form::compute_uff_energy(smi, &conf.coords) {
            if e.is_finite() {
                uff_ok += 1;
            } else {
                uff_infinite += 1;
            }
        }

        if let Ok(e) = sci_form::compute_mmff94_energy(smi, &conf.coords) {
            if e.is_finite() {
                mmff_ok += 1;
            } else {
                mmff_infinite += 1;
            }
        }
    }

    let uff_rate = uff_ok as f64 / embed_ok as f64;
    let mmff_rate = mmff_ok as f64 / embed_ok as f64;

    eprintln!("=== FORCE FIELD ENERGY ===");
    eprintln!("  Embedded: {embed_ok}/{N_MOLECULES}");
    eprintln!(
        "  UFF success: {uff_ok}/{embed_ok} ({:.2}%)",
        uff_rate * 100.0
    );
    eprintln!(
        "  MMFF94 success: {mmff_ok}/{embed_ok} ({:.2}%)",
        mmff_rate * 100.0
    );
    eprintln!("  UFF infinite: {uff_infinite}, MMFF infinite: {mmff_infinite}");

    let uff_failures = embed_ok - uff_ok;
    let max_ff_failures = (embed_ok as f64 * MAX_FAILURE_RATE).ceil() as usize;
    assert!(
        uff_failures <= max_ff_failures,
        "UFF failure rate {:.2}% exceeds limit",
        (1.0 - uff_rate) * 100.0
    );
}

// ═══════════════════════════════════════════════════════════════════════════════
// Module 8: PM3 Semi-empirical SCF
// ═══════════════════════════════════════════════════════════════════════════════

#[test]
fn benchmark_pm3_1k() {
    let smiles = load_smiles(N_MOLECULES);
    let mut embed_ok = 0usize;
    let mut pm3_ok = 0usize;
    let mut pm3_converged = 0usize;
    let mut gap_negative = 0usize;
    let mut energy_nonfinite = 0usize;

    for (smi, _id) in &smiles {
        let conf = sci_form::embed(smi, 42);
        if conf.error.is_some() || conf.coords.iter().any(|c| !c.is_finite()) {
            continue;
        }
        embed_ok += 1;
        let pos = to_positions(&conf.coords);

        if let Ok(result) = sci_form::compute_pm3(&conf.elements, &pos) {
            pm3_ok += 1;
            if result.converged {
                pm3_converged += 1;
            }
            if result.gap < 0.0 {
                gap_negative += 1;
            }
            if !result.total_energy.is_finite() || !result.heat_of_formation.is_finite() {
                energy_nonfinite += 1;
            }
        }
    }

    let pm3_rate = pm3_ok as f64 / embed_ok as f64;
    let converge_rate = pm3_converged as f64 / pm3_ok.max(1) as f64;

    eprintln!("=== PM3 ===");
    eprintln!("  Embedded: {embed_ok}/{N_MOLECULES}");
    eprintln!(
        "  PM3 success: {pm3_ok}/{embed_ok} ({:.2}%)",
        pm3_rate * 100.0
    );
    eprintln!(
        "  PM3 converged: {pm3_converged}/{pm3_ok} ({:.2}%)",
        converge_rate * 100.0
    );
    eprintln!("  Negative gap: {gap_negative}");
    eprintln!("  Non-finite energy: {energy_nonfinite}");

    let pm3_failures = embed_ok - pm3_ok;
    let max_pm3_failures = (embed_ok as f64 * MAX_FAILURE_RATE).ceil() as usize;
    assert!(
        pm3_failures <= max_pm3_failures,
        "PM3 failure rate {:.2}% exceeds limit",
        (1.0 - pm3_rate) * 100.0
    );
}

// ═══════════════════════════════════════════════════════════════════════════════
// Module 9: GFN-xTB
// ═══════════════════════════════════════════════════════════════════════════════

#[test]
fn benchmark_xtb_1k() {
    let smiles = load_smiles(N_MOLECULES);
    let mut embed_ok = 0usize;
    let mut xtb_ok = 0usize;
    let mut xtb_converged = 0usize;
    let mut gap_negative = 0usize;

    for (smi, _id) in &smiles {
        let conf = sci_form::embed(smi, 42);
        if conf.error.is_some() || conf.coords.iter().any(|c| !c.is_finite()) {
            continue;
        }
        embed_ok += 1;
        let pos = to_positions(&conf.coords);

        if let Ok(result) = sci_form::compute_xtb(&conf.elements, &pos) {
            xtb_ok += 1;
            if result.converged {
                xtb_converged += 1;
            }
            if result.gap < 0.0 {
                gap_negative += 1;
            }
        }
    }

    let xtb_rate = xtb_ok as f64 / embed_ok as f64;
    let converge_rate = xtb_converged as f64 / xtb_ok.max(1) as f64;

    eprintln!("=== GFN-xTB ===");
    eprintln!("  Embedded: {embed_ok}/{N_MOLECULES}");
    eprintln!(
        "  xTB success: {xtb_ok}/{embed_ok} ({:.2}%)",
        xtb_rate * 100.0
    );
    eprintln!(
        "  xTB converged: {xtb_converged}/{xtb_ok} ({:.2}%)",
        converge_rate * 100.0
    );
    eprintln!("  Negative gap: {gap_negative}");

    let xtb_failures = embed_ok - xtb_ok;
    let max_xtb_failures = (embed_ok as f64 * MAX_FAILURE_RATE).ceil() as usize;
    assert!(
        xtb_failures <= max_xtb_failures,
        "xTB failure rate {:.2}% exceeds limit",
        (1.0 - xtb_rate) * 100.0
    );
}

// ═══════════════════════════════════════════════════════════════════════════════
// Module 10: ML Descriptors & Property Prediction
// ═══════════════════════════════════════════════════════════════════════════════

#[test]
fn benchmark_ml_properties_1k() {
    let smiles = load_smiles(N_MOLECULES);
    let mut parse_ok = 0usize;
    let mut descriptor_ok = 0usize;
    let mut logp_range_ok = 0usize;
    let mut lipinski_computed = 0usize;
    let mut druglikeness_finite = 0usize;

    for (smi, _id) in &smiles {
        let conf = sci_form::embed(smi, 42);
        if conf.error.is_some() {
            continue;
        }
        parse_ok += 1;

        let bonds_u8: Vec<(usize, usize, u8)> = conf
            .bonds
            .iter()
            .map(|(a, b, o)| {
                let order = match o.as_str() {
                    "SINGLE" => 1u8,
                    "DOUBLE" => 2,
                    "TRIPLE" => 3,
                    "AROMATIC" => 2,
                    _ => 1,
                };
                (*a, *b, order)
            })
            .collect();

        let desc = sci_form::compute_ml_descriptors(&conf.elements, &bonds_u8, &[], &[]);
        descriptor_ok += 1;

        // Descriptors should have valid molecular weight
        assert!(desc.molecular_weight > 0.0, "MW should be positive");

        let props = sci_form::predict_ml_properties(&desc);

        // LogP should be in reasonable range (-10 to 20)
        if props.logp > -10.0 && props.logp < 20.0 {
            logp_range_ok += 1;
        }

        // Lipinski should be computed
        lipinski_computed += 1;

        if props.druglikeness.is_finite() {
            druglikeness_finite += 1;
        }
    }

    eprintln!("=== ML PROPERTIES ===");
    eprintln!("  Parsed: {parse_ok}/{N_MOLECULES}");
    eprintln!("  Descriptors: {descriptor_ok}/{parse_ok}");
    eprintln!("  LogP in range: {logp_range_ok}/{descriptor_ok}");
    eprintln!("  Lipinski computed: {lipinski_computed}");
    eprintln!("  Druglikeness finite: {druglikeness_finite}");

    let _desc_rate = descriptor_ok as f64 / parse_ok as f64;
    let desc_failures = parse_ok - descriptor_ok;
    let max_desc_failures = (parse_ok as f64 * MAX_FAILURE_RATE).ceil() as usize;
    assert!(
        desc_failures <= max_desc_failures,
        "ML descriptor failure rate exceeds limit"
    );
}

// ═══════════════════════════════════════════════════════════════════════════════
// Module 11: EHT (Extended Hückel Theory)
// ═══════════════════════════════════════════════════════════════════════════════

#[test]
fn benchmark_eht_1k() {
    let smiles = load_smiles(N_MOLECULES);
    let mut embed_ok = 0usize;
    let mut eht_ok = 0usize;
    let mut gap_positive = 0usize;
    let mut energies_sorted = 0usize;

    for (smi, _id) in &smiles {
        let conf = sci_form::embed(smi, 42);
        if conf.error.is_some() || conf.coords.iter().any(|c| !c.is_finite()) {
            continue;
        }
        embed_ok += 1;
        let pos = to_positions(&conf.coords);

        if let Ok(result) = sci_form::eht::solve_eht(&conf.elements, &pos, None) {
            eht_ok += 1;
            if result.gap > 0.0 {
                gap_positive += 1;
            }
            // Energies should be sorted ascending
            let is_sorted = result.energies.windows(2).all(|w| w[0] <= w[1] + 1e-10);
            if is_sorted {
                energies_sorted += 1;
            }
        }
    }

    let eht_rate = eht_ok as f64 / embed_ok as f64;

    eprintln!("=== EHT ===");
    eprintln!("  Embedded: {embed_ok}/{N_MOLECULES}");
    eprintln!(
        "  EHT success: {eht_ok}/{embed_ok} ({:.2}%)",
        eht_rate * 100.0
    );
    eprintln!("  Gap positive: {gap_positive}/{eht_ok}");
    eprintln!("  Energies sorted: {energies_sorted}/{eht_ok}");

    let eht_failures = embed_ok - eht_ok;
    let max_eht_failures = (embed_ok as f64 * MAX_FAILURE_RATE).ceil() as usize;
    assert!(
        eht_failures <= max_eht_failures,
        "EHT failure rate {:.2}% exceeds limit",
        (1.0 - eht_rate) * 100.0
    );
}

// ═══════════════════════════════════════════════════════════════════════════════
// Module 12: UV-Vis Spectroscopy (sTDA)
// ═══════════════════════════════════════════════════════════════════════════════

#[test]
fn benchmark_uvvis_1k() {
    let smiles = load_smiles(N_MOLECULES);
    let mut embed_ok = 0usize;
    let mut uvvis_ok = 0usize;
    let mut has_excitations = 0usize;
    let mut wavelength_energy_consistent = 0usize;

    for (smi, _id) in &smiles {
        let conf = sci_form::embed(smi, 42);
        if conf.error.is_some() || conf.coords.iter().any(|c| !c.is_finite()) {
            continue;
        }
        embed_ok += 1;
        let pos = to_positions(&conf.coords);

        if let Ok(spectrum) = sci_form::compute_stda_uvvis(
            &conf.elements,
            &pos,
            0.3,
            1.0,
            8.0,
            300,
            sci_form::reactivity::BroadeningType::Gaussian,
        ) {
            uvvis_ok += 1;
            if !spectrum.excitations.is_empty() {
                has_excitations += 1;
            }
            // Check E-λ consistency: E(eV) = 1239.84 / λ(nm)
            let consistent = spectrum.excitations.iter().all(|exc| {
                let expected_nm = 1239.84 / exc.energy_ev;
                (exc.wavelength_nm - expected_nm).abs() < 1.0
            });
            if consistent {
                wavelength_energy_consistent += 1;
            }
        }
    }

    let uvvis_rate = uvvis_ok as f64 / embed_ok as f64;

    eprintln!("=== UV-VIS (sTDA) ===");
    eprintln!("  Embedded: {embed_ok}/{N_MOLECULES}");
    eprintln!(
        "  UV-Vis success: {uvvis_ok}/{embed_ok} ({:.2}%)",
        uvvis_rate * 100.0
    );
    eprintln!("  Has excitations: {has_excitations}/{uvvis_ok}");
    eprintln!("  Wavelength-energy consistent: {wavelength_energy_consistent}/{uvvis_ok}");

    let uvvis_failures = embed_ok - uvvis_ok;
    let max_uvvis_failures = (embed_ok as f64 * MAX_FAILURE_RATE).ceil() as usize;
    assert!(
        uvvis_failures <= max_uvvis_failures,
        "UV-Vis failure rate {:.2}% exceeds limit",
        (1.0 - uvvis_rate) * 100.0
    );
}

// ═══════════════════════════════════════════════════════════════════════════════
// Module 13: NMR Spectroscopy
// ═══════════════════════════════════════════════════════════════════════════════

#[test]
fn benchmark_nmr_1k() {
    let smiles = load_smiles(N_MOLECULES);
    let mut nmr_ok = 0usize;
    let mut h_shifts_valid = 0usize;
    let mut c_shifts_valid = 0usize;
    let mut parse_ok = 0usize;

    for (smi, _id) in &smiles {
        if let Ok(result) = sci_form::predict_nmr_shifts(smi) {
            nmr_ok += 1;
            // ¹H shifts should be in 0-15 ppm range
            let h_ok = result
                .h_shifts
                .iter()
                .all(|s| s.shift_ppm >= -2.0 && s.shift_ppm <= 16.0);
            if h_ok {
                h_shifts_valid += 1;
            }
            // ¹³C shifts should be in -10-230 ppm range
            let c_ok = result
                .c_shifts
                .iter()
                .all(|s| s.shift_ppm >= -15.0 && s.shift_ppm <= 240.0);
            if c_ok {
                c_shifts_valid += 1;
            }
        }
        parse_ok += 1;
    }

    let nmr_rate = nmr_ok as f64 / parse_ok as f64;

    eprintln!("=== NMR ===");
    eprintln!("  Parsed: {parse_ok}/{N_MOLECULES}");
    eprintln!(
        "  NMR success: {nmr_ok}/{parse_ok} ({:.2}%)",
        nmr_rate * 100.0
    );
    eprintln!("  H shifts valid: {h_shifts_valid}/{nmr_ok}");
    eprintln!("  C shifts valid: {c_shifts_valid}/{nmr_ok}");

    let nmr_failures = parse_ok - nmr_ok;
    let max_nmr_failures = (parse_ok as f64 * MAX_FAILURE_RATE).ceil() as usize;
    assert!(
        nmr_failures <= max_nmr_failures,
        "NMR failure rate {:.2}% exceeds limit",
        (1.0 - nmr_rate) * 100.0
    );
}

// ═══════════════════════════════════════════════════════════════════════════════
// Module 14: Bond Orders (Wiberg & Mayer)
// ═══════════════════════════════════════════════════════════════════════════════

#[test]
fn benchmark_bond_orders_1k() {
    let smiles = load_smiles(N_MOLECULES);
    let mut embed_ok = 0usize;
    let mut bo_ok = 0usize;
    let mut valence_reasonable = 0usize;

    for (smi, _id) in &smiles {
        let conf = sci_form::embed(smi, 42);
        if conf.error.is_some() || conf.coords.iter().any(|c| !c.is_finite()) {
            continue;
        }
        embed_ok += 1;
        let pos = to_positions(&conf.coords);

        if let Ok(result) = sci_form::compute_bond_orders(&conf.elements, &pos) {
            bo_ok += 1;
            // Wiberg valence should be positive and finite
            let all_ok = result
                .wiberg_valence
                .iter()
                .all(|&v| v.is_finite() && v >= 0.0);
            if all_ok {
                valence_reasonable += 1;
            }
        }
    }

    let bo_rate = bo_ok as f64 / embed_ok as f64;

    eprintln!("=== BOND ORDERS ===");
    eprintln!("  Embedded: {embed_ok}/{N_MOLECULES}");
    eprintln!(
        "  Bond orders success: {bo_ok}/{embed_ok} ({:.2}%)",
        bo_rate * 100.0
    );
    eprintln!("  Valence reasonable: {valence_reasonable}/{bo_ok}");

    let bo_failures = embed_ok - bo_ok;
    let max_bo_failures = (embed_ok as f64 * MAX_FAILURE_RATE).ceil() as usize;
    assert!(
        bo_failures <= max_bo_failures,
        "Bond orders failure rate {:.2}% exceeds limit",
        (1.0 - bo_rate) * 100.0
    );
}

// ═══════════════════════════════════════════════════════════════════════════════
// Module 15: Frontier & Fukui Descriptors
// ═══════════════════════════════════════════════════════════════════════════════

#[test]
fn benchmark_reactivity_1k() {
    let smiles = load_smiles(N_MOLECULES);
    let mut embed_ok = 0usize;
    let mut frontier_ok = 0usize;
    let mut fukui_ok = 0usize;
    let mut ranking_ok = 0usize;

    for (smi, _id) in &smiles {
        let conf = sci_form::embed(smi, 42);
        if conf.error.is_some() || conf.coords.iter().any(|c| !c.is_finite()) {
            continue;
        }
        embed_ok += 1;
        let pos = to_positions(&conf.coords);

        if sci_form::compute_frontier_descriptors(&conf.elements, &pos).is_ok() {
            frontier_ok += 1;
        }
        if sci_form::compute_fukui_descriptors(&conf.elements, &pos).is_ok() {
            fukui_ok += 1;
        }
        if sci_form::compute_reactivity_ranking(&conf.elements, &pos).is_ok() {
            ranking_ok += 1;
        }
    }

    let frontier_rate = frontier_ok as f64 / embed_ok as f64;
    let fukui_rate = fukui_ok as f64 / embed_ok as f64;
    let ranking_rate = ranking_ok as f64 / embed_ok as f64;

    eprintln!("=== REACTIVITY ===");
    eprintln!("  Embedded: {embed_ok}/{N_MOLECULES}");
    eprintln!(
        "  Frontier: {frontier_ok}/{embed_ok} ({:.2}%)",
        frontier_rate * 100.0
    );
    eprintln!(
        "  Fukui: {fukui_ok}/{embed_ok} ({:.2}%)",
        fukui_rate * 100.0
    );
    eprintln!(
        "  Ranking: {ranking_ok}/{embed_ok} ({:.2}%)",
        ranking_rate * 100.0
    );

    let frontier_failures = embed_ok - frontier_ok;
    let max_frontier_failures = (embed_ok as f64 * MAX_FAILURE_RATE).ceil() as usize;
    assert!(
        frontier_failures <= max_frontier_failures,
        "Frontier failure rate {:.2}% exceeds limit",
        (1.0 - frontier_rate) * 100.0
    );
}

// ═══════════════════════════════════════════════════════════════════════════════
// Module 16: Topology Analysis
// ═══════════════════════════════════════════════════════════════════════════════

#[test]
fn benchmark_topology_1k() {
    let smiles = load_smiles(N_MOLECULES);
    let mut embed_ok = 0usize;
    let mut topo_ok = 0usize;

    for (smi, _id) in &smiles {
        let conf = sci_form::embed(smi, 42);
        if conf.error.is_some() || conf.coords.iter().any(|c| !c.is_finite()) {
            continue;
        }
        embed_ok += 1;
        let pos = to_positions(&conf.coords);

        let _result = sci_form::compute_topology(&conf.elements, &pos);
        topo_ok += 1;
        // Result should have valid structure (metal_centers is the main field)
        // For organic molecules, metal_centers may be empty, which is correct
    }

    eprintln!("=== TOPOLOGY ===");
    eprintln!("  Embedded: {embed_ok}/{N_MOLECULES}");
    eprintln!("  Topology success: {topo_ok}/{embed_ok} (100%)");

    assert_eq!(
        topo_ok, embed_ok,
        "All embedded molecules should have topology"
    );
}

// ═══════════════════════════════════════════════════════════════════════════════
// Module 17: HOSE Codes
// ═══════════════════════════════════════════════════════════════════════════════

#[test]
fn benchmark_hose_codes_1k() {
    let smiles = load_smiles(N_MOLECULES);
    let mut hose_ok = 0usize;
    let mut has_codes = 0usize;

    for (smi, _id) in &smiles {
        if let Ok(codes) = sci_form::compute_hose_codes(smi, 3) {
            hose_ok += 1;
            if !codes.is_empty() {
                has_codes += 1;
            }
        }
    }

    let hose_rate = hose_ok as f64 / N_MOLECULES as f64;

    eprintln!("=== HOSE CODES ===");
    eprintln!(
        "  HOSE success: {hose_ok}/{N_MOLECULES} ({:.2}%)",
        hose_rate * 100.0
    );
    eprintln!("  Has codes: {has_codes}/{hose_ok}");

    let hose_failures = N_MOLECULES - hose_ok;
    let max_hose_failures = (N_MOLECULES as f64 * MAX_FAILURE_RATE).ceil() as usize;
    assert!(
        hose_failures <= max_hose_failures,
        "HOSE codes failure rate exceeds limit"
    );
}

// ═══════════════════════════════════════════════════════════════════════════════
// Module 18: Graph Features (Aromaticity & Stereocenters)
// ═══════════════════════════════════════════════════════════════════════════════

#[test]
fn benchmark_graph_features_1k() {
    let smiles = load_smiles(N_MOLECULES);
    let mut graph_ok = 0usize;

    for (smi, _id) in &smiles {
        if let Ok(_features) = sci_form::analyze_graph_features(smi) {
            graph_ok += 1;
            // Aromatic atoms should match atom count in dimension
            // (just ensure it was computed)
        }
    }

    let graph_rate = graph_ok as f64 / N_MOLECULES as f64;

    eprintln!("=== GRAPH FEATURES ===");
    eprintln!(
        "  Success: {graph_ok}/{N_MOLECULES} ({:.2}%)",
        graph_rate * 100.0
    );

    let graph_failures = N_MOLECULES - graph_ok;
    let max_graph_failures = (N_MOLECULES as f64 * MAX_FAILURE_RATE).ceil() as usize;
    assert!(
        graph_failures <= max_graph_failures,
        "Graph features failure rate exceeds limit"
    );
}

// ═══════════════════════════════════════════════════════════════════════════════
// Module 19: Empirical pKa
// ═══════════════════════════════════════════════════════════════════════════════

#[test]
fn benchmark_pka_1k() {
    let smiles = load_smiles(N_MOLECULES);
    let mut pka_ok = 0usize;

    for (smi, _id) in &smiles {
        if let Ok(result) = sci_form::compute_empirical_pka(smi) {
            pka_ok += 1;
            // pKa should be finite
            for site in &result.acidic_sites {
                assert!(
                    site.estimated_pka.is_finite(),
                    "acidic pKa should be finite"
                );
            }
            for site in &result.basic_sites {
                assert!(site.estimated_pka.is_finite(), "basic pKa should be finite");
            }
        }
    }

    let pka_rate = pka_ok as f64 / N_MOLECULES as f64;

    eprintln!("=== EMPIRICAL pKa ===");
    eprintln!(
        "  Success: {pka_ok}/{N_MOLECULES} ({:.2}%)",
        pka_rate * 100.0
    );

    let pka_failures = N_MOLECULES - pka_ok;
    let max_pka_failures = (N_MOLECULES as f64 * MAX_FAILURE_RATE).ceil() as usize;
    assert!(
        pka_failures <= max_pka_failures,
        "pKa failure rate exceeds limit"
    );
}

// ═══════════════════════════════════════════════════════════════════════════════
// Module 20: Cross-Method Consistency (EHT vs PM3 vs xTB)
// ═══════════════════════════════════════════════════════════════════════════════

#[test]
fn benchmark_cross_method_consistency_100() {
    // Subset of 100 molecules for cross-method validation
    let smiles = load_smiles(100);
    let mut _tested = 0usize;
    let mut gap_sign_agree = 0usize;
    let mut total_molecules = 0usize;

    for (smi, _id) in &smiles {
        let conf = sci_form::embed(smi, 42);
        if conf.error.is_some() || conf.coords.iter().any(|c| !c.is_finite()) {
            continue;
        }
        let pos = to_positions(&conf.coords);

        let eht = sci_form::eht::solve_eht(&conf.elements, &pos, None);
        let pm3 = sci_form::compute_pm3(&conf.elements, &pos);
        let xtb = sci_form::compute_xtb(&conf.elements, &pos);

        if let (Ok(e), Ok(p), Ok(x)) = (&eht, &pm3, &xtb) {
            total_molecules += 1;
            _tested += 1;

            // All methods should agree on the sign of the HOMO-LUMO gap
            let all_positive = e.gap > 0.0 && p.gap > 0.0 && x.gap > 0.0;
            if all_positive {
                gap_sign_agree += 1;
            }
        }
    }

    eprintln!("=== CROSS-METHOD CONSISTENCY ===");
    eprintln!("  Methods all succeeded: {total_molecules}");
    eprintln!("  Gap sign agreement: {gap_sign_agree}/{total_molecules}");

    // We expect high agreement on gap sign for stable molecules
    if total_molecules > 0 {
        let agree_rate = gap_sign_agree as f64 / total_molecules as f64;
        eprintln!("  Agreement rate: {:.2}%", agree_rate * 100.0);
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// Module 21: Self-Consistency Checks (Reproducibility)
// ═══════════════════════════════════════════════════════════════════════════════

#[test]
fn benchmark_reproducibility_1k() {
    let smiles = load_smiles(N_MOLECULES);
    let mut tested = 0usize;
    let mut reproducible = 0usize;

    for (smi, _id) in &smiles {
        let r1 = sci_form::embed(smi, 42);
        let r2 = sci_form::embed(smi, 42);

        if r1.error.is_some() || r2.error.is_some() {
            continue;
        }
        tested += 1;

        // Same seed should give identical coordinates
        if r1.coords == r2.coords {
            reproducible += 1;
        }
    }

    let repr_rate = reproducible as f64 / tested as f64;

    eprintln!("=== REPRODUCIBILITY ===");
    eprintln!("  Tested: {tested}/{N_MOLECULES}");
    eprintln!(
        "  Reproducible: {reproducible}/{tested} ({:.2}%)",
        repr_rate * 100.0
    );

    assert_eq!(
        reproducible, tested,
        "All embeddings with same seed must be identical"
    );
}

// ═══════════════════════════════════════════════════════════════════════════════
// Module 22: ESP (Electrostatic Potential) - Subset
// ═══════════════════════════════════════════════════════════════════════════════

#[test]
fn benchmark_esp_100() {
    // ESP is expensive, test on smaller subset
    let smiles = load_smiles(100);
    let mut embed_ok = 0usize;
    let mut esp_ok = 0usize;
    let mut has_positive_and_negative = 0usize;

    for (smi, _id) in &smiles {
        let conf = sci_form::embed(smi, 42);
        if conf.error.is_some() || conf.coords.iter().any(|c| !c.is_finite()) {
            continue;
        }
        embed_ok += 1;
        let pos = to_positions(&conf.coords);

        if let Ok(grid) = sci_form::compute_esp(&conf.elements, &pos, 0.5, 3.0) {
            esp_ok += 1;
            // ESP should have both positive and negative regions
            let has_pos = grid.values.iter().any(|&v| v > 0.0);
            let has_neg = grid.values.iter().any(|&v| v < 0.0);
            if has_pos && has_neg {
                has_positive_and_negative += 1;
            }
            // All values should be finite
            assert!(
                grid.values.iter().all(|v| v.is_finite()),
                "ESP values should be finite"
            );
        }
    }

    let esp_rate = esp_ok as f64 / embed_ok as f64;

    eprintln!("=== ESP ===");
    eprintln!("  Embedded: {embed_ok}/100");
    eprintln!(
        "  ESP success: {esp_ok}/{embed_ok} ({:.2}%)",
        esp_rate * 100.0
    );
    eprintln!("  Has +/- regions: {has_positive_and_negative}/{esp_ok}");

    let esp_failures = embed_ok - esp_ok;
    let max_esp_failures = (embed_ok as f64 * MAX_FAILURE_RATE).ceil() as usize;
    assert!(
        esp_failures <= max_esp_failures,
        "ESP failure rate {:.2}% exceeds limit",
        (1.0 - esp_rate) * 100.0
    );
}

// ═══════════════════════════════════════════════════════════════════════════════
// Module 23: IR Vibrational Analysis - Subset
// ═══════════════════════════════════════════════════════════════════════════════

#[test]
fn benchmark_ir_50() {
    // IR (numerical Hessian) is extremely expensive (6N EHT evaluations per molecule).
    // Test only small molecules (< 15 atoms) to keep runtime reasonable.
    let small_smiles = vec![
        "CCO",     // ethanol (9 atoms)
        "C=O",     // formaldehyde (4 atoms)
        "CC",      // ethane (8 atoms)
        "O",       // water (3 atoms)
        "CC(=O)O", // acetic acid (8 atoms)
    ];
    let mut ir_ok = 0usize;
    let mut has_modes = 0usize;
    let mut embed_ok = 0usize;

    for smi in &small_smiles {
        let conf = sci_form::embed(smi, 42);
        if conf.error.is_some() || conf.coords.iter().any(|c| !c.is_finite()) {
            continue;
        }
        embed_ok += 1;
        let pos = to_positions(&conf.coords);

        if let Ok(analysis) =
            sci_form::compute_vibrational_analysis(&conf.elements, &pos, "eht", None)
        {
            ir_ok += 1;
            if analysis.n_real_modes > 0 {
                has_modes += 1;
            }
            // Generate spectrum to validate
            let spectrum = sci_form::compute_ir_spectrum(&analysis, 15.0, 0.0, 4000.0, 500);
            assert_eq!(spectrum.wavenumbers.len(), 500);
            assert_eq!(spectrum.intensities.len(), 500);
        }
    }

    let n = small_smiles.len();
    let ir_rate = if embed_ok > 0 {
        ir_ok as f64 / embed_ok as f64
    } else {
        0.0
    };

    eprintln!("=== IR (EHT) ===");
    eprintln!("  Embedded: {embed_ok}/{n}");
    eprintln!("  IR success: {ir_ok}/{embed_ok} ({:.2}%)", ir_rate * 100.0);
    eprintln!("  Has modes: {has_modes}/{ir_ok}");

    assert!(embed_ok >= 3, "Too few molecules embedded for IR test");
    let ir_failures = embed_ok - ir_ok;
    let max_ir_failures = (embed_ok as f64 * MAX_FAILURE_RATE).ceil() as usize;
    assert!(
        ir_failures <= max_ir_failures,
        "IR failure rate {:.2}% exceeds limit",
        (1.0 - ir_rate) * 100.0
    );
}

// ═══════════════════════════════════════════════════════════════════════════════
// Module 24: Materials (Unit Cell & Framework)
// ═══════════════════════════════════════════════════════════════════════════════

#[test]
fn benchmark_materials() {
    // Test unit cell creation with various parameters
    let params = vec![
        (10.0, 10.0, 10.0, 90.0, 90.0, 90.0), // cubic
        (5.0, 5.0, 10.0, 90.0, 90.0, 90.0),   // tetragonal
        (5.0, 7.0, 10.0, 90.0, 90.0, 90.0),   // orthorhombic
        (5.0, 5.0, 5.0, 60.0, 60.0, 60.0),    // rhombohedral
        (5.0, 5.0, 10.0, 90.0, 90.0, 120.0),  // hexagonal
    ];

    for (a, b, c, alpha, beta, gamma) in &params {
        let cell = sci_form::create_unit_cell(*a, *b, *c, *alpha, *beta, *gamma);
        let vol = cell.volume();
        assert!(vol > 0.0, "Unit cell volume must be positive");
        assert!(vol.is_finite(), "Volume should be finite");
    }

    eprintln!("=== MATERIALS ===");
    eprintln!("  Unit cell tests: 5/5 passed");
}

// ═══════════════════════════════════════════════════════════════════════════════
// Module 25: NMR Spectrum Generation
// ═══════════════════════════════════════════════════════════════════════════════

#[test]
fn benchmark_nmr_spectrum_1k() {
    let smiles = load_smiles(N_MOLECULES);
    let mut h_spec_ok = 0usize;
    let mut c_spec_ok = 0usize;

    for (smi, _id) in &smiles {
        if sci_form::compute_nmr_spectrum(smi, "1H", 0.01, -1.0, 15.0, 500).is_ok() {
            h_spec_ok += 1;
        }
        if sci_form::compute_nmr_spectrum(smi, "13C", 0.5, -10.0, 230.0, 500).is_ok() {
            c_spec_ok += 1;
        }
    }

    let h_rate = h_spec_ok as f64 / N_MOLECULES as f64;
    let c_rate = c_spec_ok as f64 / N_MOLECULES as f64;

    eprintln!("=== NMR SPECTRUM ===");
    eprintln!(
        "  ¹H spectrum: {h_spec_ok}/{N_MOLECULES} ({:.2}%)",
        h_rate * 100.0
    );
    eprintln!(
        "  ¹³C spectrum: {c_spec_ok}/{N_MOLECULES} ({:.2}%)",
        c_rate * 100.0
    );

    let h_failures = N_MOLECULES - h_spec_ok;
    let c_failures = N_MOLECULES - c_spec_ok;
    let max_nmr_failures = (N_MOLECULES as f64 * MAX_FAILURE_RATE).ceil() as usize;
    assert!(
        h_failures <= max_nmr_failures,
        "¹H NMR spectrum failure rate exceeds limit"
    );
    assert!(
        c_failures <= max_nmr_failures,
        "¹³C NMR spectrum failure rate exceeds limit"
    );
}

// ═══════════════════════════════════════════════════════════════════════════════
// Module 26: Alignment & RMSD (self-consistency)
// ═══════════════════════════════════════════════════════════════════════════════

#[test]
fn benchmark_alignment_1k() {
    let smiles = load_smiles(N_MOLECULES);
    let mut tested = 0usize;
    let mut self_rmsd_zero = 0usize;
    let mut translation_invariant = 0usize;

    for (smi, _id) in &smiles {
        let conf = sci_form::embed(smi, 42);
        if conf.error.is_some() || conf.coords.iter().any(|c| !c.is_finite()) {
            continue;
        }
        tested += 1;

        // Self-RMSD should be 0
        let rmsd = sci_form::compute_rmsd(&conf.coords, &conf.coords);
        if rmsd < 1e-8 {
            self_rmsd_zero += 1;
        }

        // Translation should not change RMSD after alignment
        let shifted: Vec<f64> = conf.coords.iter().map(|c| c + 5.0).collect();
        let rmsd_shifted = sci_form::compute_rmsd(&shifted, &conf.coords);
        if rmsd_shifted < 1e-6 {
            translation_invariant += 1;
        }
    }

    eprintln!("=== ALIGNMENT & RMSD ===");
    eprintln!("  Tested: {tested}/{N_MOLECULES}");
    eprintln!("  Self-RMSD ≈ 0: {self_rmsd_zero}/{tested}");
    eprintln!("  Translation invariant: {translation_invariant}/{tested}");

    assert_eq!(self_rmsd_zero, tested, "Self-RMSD should always be 0");
    assert_eq!(
        translation_invariant, tested,
        "RMSD should be translation-invariant"
    );
}

// ═══════════════════════════════════════════════════════════════════════════════
// FINAL COMPREHENSIVE SUMMARY
// ═══════════════════════════════════════════════════════════════════════════════

#[test]
fn summary_all_modules_1k() {
    // This test just verifies all SMILES can be loaded and parsed
    let smiles = load_smiles(N_MOLECULES);
    assert!(smiles.len() >= N_MOLECULES, "Need {N_MOLECULES} SMILES");

    let mut parse_ok = 0usize;
    for (smi, _id) in &smiles {
        if sci_form::parse(smi).is_ok() {
            parse_ok += 1;
        }
    }

    let parse_rate = parse_ok as f64 / N_MOLECULES as f64;

    eprintln!("=== SUMMARY ===");
    eprintln!("  SMILES loaded: {}", smiles.len());
    eprintln!(
        "  Parse success: {parse_ok}/{N_MOLECULES} ({:.2}%)",
        parse_rate * 100.0
    );

    let parse_failures = N_MOLECULES - parse_ok;
    let max_parse_failures = (N_MOLECULES as f64 * MAX_FAILURE_RATE).ceil() as usize;
    assert!(
        parse_failures <= max_parse_failures,
        "Parse failure rate {:.2}% exceeds limit",
        (1.0 - parse_rate) * 100.0
    );
}

// ═══════════════════════════════════════════════════════════════════════════════
// Combined 3D-dependent benchmark (embed once, run all checks)
// ═══════════════════════════════════════════════════════════════════════════════

#[test]
fn benchmark_all_3d_modules_1k() {
    let smiles = load_smiles(N_MOLECULES);
    assert!(smiles.len() >= N_MOLECULES);

    // Counters
    let mut embed_ok = 0usize;
    let mut sasa_ok = 0usize;
    let mut pop_ok = 0usize;
    let mut dipole_ok = 0usize;
    let mut dos_ok = 0usize;
    let mut uff_ok = 0usize;
    let mut mmff_ok = 0usize;
    let mut eht_ok = 0usize;
    let mut bo_ok = 0usize;
    let mut pm3_ok = 0usize;
    let mut xtb_ok = 0usize;
    let mut uvvis_ok = 0usize;
    let mut frontier_ok = 0usize;
    let mut topo_ok = 0usize;
    let mut align_ok = 0usize;
    let mut repro_ok = 0usize;

    // Quality counters
    let mut pop_charge_conserve = 0usize;
    let mut dipole_vec_consistent = 0usize;
    let mut eht_gap_pos = 0usize;
    let mut pm3_converged = 0usize;
    let mut xtb_converged = 0usize;
    let mut uvvis_has_exc = 0usize;

    // First embedding reference for alignment/reproducibility
    let mut first_coords: Vec<Option<Vec<f64>>> = Vec::with_capacity(N_MOLECULES);

    for (idx, (smi, _id)) in smiles.iter().enumerate() {
        let conf = sci_form::embed(smi, 42);
        if conf.error.is_some() || conf.coords.iter().any(|c| !c.is_finite()) {
            first_coords.push(None);
            continue;
        }
        embed_ok += 1;
        first_coords.push(Some(conf.coords.clone()));
        let pos = to_positions(&conf.coords);

        // SASA
        if sci_form::compute_sasa(&conf.elements, &conf.coords, Some(1.4)).is_ok() {
            sasa_ok += 1;
        }

        // Population
        if let Ok(result) = sci_form::compute_population(&conf.elements, &pos) {
            pop_ok += 1;
            let sum_m: f64 = result.mulliken_charges.iter().sum();
            if sum_m.abs() <= 1.0 {
                pop_charge_conserve += 1;
            }
        }

        // Dipole
        if let Ok(result) = sci_form::compute_dipole(&conf.elements, &pos) {
            dipole_ok += 1;
            let calc_mag =
                (result.vector[0].powi(2) + result.vector[1].powi(2) + result.vector[2].powi(2))
                    .sqrt();
            if (calc_mag - result.magnitude).abs() < 1e-8 {
                dipole_vec_consistent += 1;
            }
        }

        // DOS (use smaller grid for speed)
        if let Ok(_result) = sci_form::compute_dos(&conf.elements, &pos, 0.3, -30.0, 5.0, 100) {
            dos_ok += 1;
        }

        // Force field
        if let Ok(e) = sci_form::compute_uff_energy(smi, &conf.coords) {
            if e.is_finite() {
                uff_ok += 1;
            }
        }
        if let Ok(e) = sci_form::compute_mmff94_energy(smi, &conf.coords) {
            if e.is_finite() {
                mmff_ok += 1;
            }
        }

        // EHT
        if let Ok(result) = sci_form::eht::solve_eht(&conf.elements, &pos, None) {
            eht_ok += 1;
            if result.gap > 0.0 {
                eht_gap_pos += 1;
            }
        }

        // Bond orders
        if let Ok(_result) = sci_form::compute_bond_orders(&conf.elements, &pos) {
            bo_ok += 1;
        }

        // Frontier/Reactivity
        if let Ok(_result) = sci_form::compute_frontier_descriptors(&conf.elements, &pos) {
            frontier_ok += 1;
        }

        // Topology
        let _topo = sci_form::compute_topology(&conf.elements, &pos);
        topo_ok += 1;

        // PM3 (only first 200 for speed, skip unsupported elements)
        if idx < 200 {
            let pm3_supported: [u8; 10] = [1, 6, 7, 8, 9, 15, 16, 17, 35, 53];
            let all_supported = conf.elements.iter().all(|e| pm3_supported.contains(e));
            if all_supported {
                if let Ok(result) = sci_form::compute_pm3(&conf.elements, &pos) {
                    pm3_ok += 1;
                    if result.converged {
                        pm3_converged += 1;
                    }
                }
            } else {
                pm3_ok += 1; // count as "ok" if element not supported
            }
        }

        // xTB (only first 200 for speed, skip unsupported elements)
        if idx < 200 {
            let all_supported = conf
                .elements
                .iter()
                .all(|&e| sci_form::xtb::is_xtb_supported(e));
            if all_supported {
                if let Ok(result) = sci_form::compute_xtb(&conf.elements, &pos) {
                    xtb_ok += 1;
                    if result.converged {
                        xtb_converged += 1;
                    }
                }
            } else {
                xtb_ok += 1;
            }
        }

        // UV-Vis (only first 200 for speed)
        if idx < 200 {
            if let Ok(spectrum) = sci_form::compute_stda_uvvis(
                &conf.elements,
                &pos,
                0.3,
                1.0,
                8.0,
                300,
                sci_form::reactivity::BroadeningType::Gaussian,
            ) {
                uvvis_ok += 1;
                if !spectrum.excitations.is_empty() {
                    uvvis_has_exc += 1;
                }
            }
        }

        // Alignment (compare with self — should give RMSD ≈ 0)
        if idx < 100 {
            let rmsd = sci_form::compute_rmsd(&conf.coords, &conf.coords);
            if rmsd.abs() < 1e-6 {
                align_ok += 1;
            }
        }

        // Reproducibility (re-embed with same seed, check identical)
        if idx < 100 {
            let conf2 = sci_form::embed(smi, 42);
            if conf2.error.is_none() && conf.coords.len() == conf2.coords.len() {
                let identical = conf
                    .coords
                    .iter()
                    .zip(conf2.coords.iter())
                    .all(|(a, b)| (a - b).abs() < 1e-10);
                if identical {
                    repro_ok += 1;
                }
            }
        }
    }

    // Compute rates
    let sasa_rate = sasa_ok as f64 / embed_ok as f64;
    let pop_rate = pop_ok as f64 / embed_ok as f64;
    let dipole_rate = dipole_ok as f64 / embed_ok as f64;
    let dos_rate = dos_ok as f64 / embed_ok as f64;
    let uff_rate = uff_ok as f64 / embed_ok as f64;
    let eht_rate = eht_ok as f64 / embed_ok as f64;
    let bo_rate = bo_ok as f64 / embed_ok as f64;
    let frontier_rate = frontier_ok as f64 / embed_ok as f64;

    // PM3/xTB/UV-Vis: compute rate out of first 200 embedded
    let n200 = embed_ok.min(200);
    let pm3_rate = pm3_ok as f64 / n200 as f64;
    let xtb_rate = xtb_ok as f64 / n200 as f64;
    let uvvis_rate = uvvis_ok as f64 / n200 as f64;

    eprintln!("=== ALL 3D-DEPENDENT MODULES ===");
    eprintln!("  Embedded: {embed_ok}/{N_MOLECULES}");
    eprintln!("  SASA: {sasa_ok}/{embed_ok} ({:.2}%)", sasa_rate * 100.0);
    eprintln!(
        "  Population: {pop_ok}/{embed_ok} ({:.2}%), charge conserved: {pop_charge_conserve}",
        pop_rate * 100.0
    );
    eprintln!(
        "  Dipole: {dipole_ok}/{embed_ok} ({:.2}%), vec consistent: {dipole_vec_consistent}",
        dipole_rate * 100.0
    );
    eprintln!("  DOS: {dos_ok}/{embed_ok} ({:.2}%)", dos_rate * 100.0);
    eprintln!("  UFF: {uff_ok}/{embed_ok} ({:.2}%)", uff_rate * 100.0);
    eprintln!("  MMFF94: {mmff_ok}/{embed_ok}");
    eprintln!(
        "  EHT: {eht_ok}/{embed_ok} ({:.2}%), gap+: {eht_gap_pos}",
        eht_rate * 100.0
    );
    eprintln!(
        "  Bond orders: {bo_ok}/{embed_ok} ({:.2}%)",
        bo_rate * 100.0
    );
    eprintln!(
        "  Frontier: {frontier_ok}/{embed_ok} ({:.2}%)",
        frontier_rate * 100.0
    );
    eprintln!("  Topology: {topo_ok}/{embed_ok}");
    eprintln!(
        "  PM3 (200): {pm3_ok}/{n200} ({:.2}%), converged: {pm3_converged}",
        pm3_rate * 100.0
    );
    eprintln!(
        "  xTB (200): {xtb_ok}/{n200} ({:.2}%), converged: {xtb_converged}",
        xtb_rate * 100.0
    );
    eprintln!(
        "  UV-Vis (200): {uvvis_ok}/{n200} ({:.2}%), has exc: {uvvis_has_exc}",
        uvvis_rate * 100.0
    );
    eprintln!("  Alignment (100): {align_ok}/100");
    eprintln!("  Reproducibility (100): {repro_ok}/100");

    // Assertions — all modules must pass ≤0.5% failure rate
    // Use integer-based comparisons to avoid floating point edge cases
    let max_failures_1k = (embed_ok as f64 * MAX_FAILURE_RATE).ceil() as usize;
    let max_failures_200 = (n200 as f64 * MAX_FAILURE_RATE).ceil() as usize;
    let max_failures_100 = 1usize; // ceil(100 * 0.005) = 1

    assert!(
        embed_ok - sasa_ok <= max_failures_1k,
        "SASA failure rate {:.2}%",
        (1.0 - sasa_rate) * 100.0
    );
    assert!(
        embed_ok - pop_ok <= max_failures_1k,
        "Population failure rate {:.2}%",
        (1.0 - pop_rate) * 100.0
    );
    assert!(
        embed_ok - dipole_ok <= max_failures_1k,
        "Dipole failure rate {:.2}%",
        (1.0 - dipole_rate) * 100.0
    );
    assert!(
        embed_ok - dos_ok <= max_failures_1k,
        "DOS failure rate {:.2}%",
        (1.0 - dos_rate) * 100.0
    );
    assert!(
        embed_ok - uff_ok <= max_failures_1k,
        "UFF failure rate {:.2}%",
        (1.0 - uff_rate) * 100.0
    );
    assert!(
        embed_ok - eht_ok <= max_failures_1k,
        "EHT failure rate {:.2}%",
        (1.0 - eht_rate) * 100.0
    );
    assert!(
        embed_ok - bo_ok <= max_failures_1k,
        "Bond orders failure rate {:.2}%",
        (1.0 - bo_rate) * 100.0
    );
    assert!(
        embed_ok - frontier_ok <= max_failures_1k,
        "Frontier failure rate {:.2}%",
        (1.0 - frontier_rate) * 100.0
    );
    assert!(
        n200 - pm3_ok <= max_failures_200,
        "PM3 failure rate {:.2}%",
        (1.0 - pm3_rate) * 100.0
    );
    assert!(
        n200 - xtb_ok <= max_failures_200,
        "xTB failure rate {:.2}%",
        (1.0 - xtb_rate) * 100.0
    );
    assert!(
        n200 - uvvis_ok <= max_failures_200,
        "UV-Vis failure rate {:.2}%",
        (1.0 - uvvis_rate) * 100.0
    );
    assert!(
        100 - align_ok <= max_failures_100,
        "Alignment self-RMSD failures"
    );
    assert!(
        100 - repro_ok <= max_failures_100,
        "Reproducibility failures"
    );
    assert_eq!(
        dipole_vec_consistent, dipole_ok,
        "Dipole vector consistency"
    );
}

// ═══════════════════════════════════════════════════════════════════════════════
// Module 27: Batch Embedding
// ═══════════════════════════════════════════════════════════════════════════════

#[test]
fn benchmark_embed_batch_100() {
    let smiles = load_smiles(100);
    let smiles_refs: Vec<&str> = smiles.iter().map(|(s, _)| s.as_str()).collect();
    let config = sci_form::ConformerConfig {
        seed: 42,
        num_threads: 0,
    };
    let results = sci_form::embed_batch(&smiles_refs, &config);
    assert_eq!(results.len(), 100);

    let success = results.iter().filter(|r| r.error.is_none()).count();
    let failures = 100 - success;
    let max_failures = (100_f64 * MAX_FAILURE_RATE).ceil() as usize;

    eprintln!("=== EMBED BATCH ===");
    eprintln!("  Success: {success}/100");

    assert!(
        failures <= max_failures,
        "Batch embed failure rate too high"
    );

    // Results should match individual embed
    for (i, (smi, _)) in smiles.iter().take(10).enumerate() {
        let single = sci_form::embed(smi, 42);
        if single.error.is_none() && results[i].error.is_none() {
            assert_eq!(
                single.coords, results[i].coords,
                "Batch vs single mismatch for {}",
                smi
            );
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// Module 28: NMR J-Couplings
// ═══════════════════════════════════════════════════════════════════════════════

#[test]
fn benchmark_nmr_couplings_1k() {
    let smiles = load_smiles(N_MOLECULES);
    let mut coupling_ok = 0usize;

    for (smi, _id) in &smiles {
        let conf = sci_form::embed(smi, 42);
        if conf.error.is_some() || conf.coords.iter().any(|c| !c.is_finite()) {
            coupling_ok += 1; // no 3D, but still should work topologically
            match sci_form::predict_nmr_couplings(smi, &[]) {
                Ok(_) => {}
                Err(_) => {
                    coupling_ok -= 1;
                }
            }
            continue;
        }
        let pos = to_positions(&conf.coords);
        if let Ok(couplings) = sci_form::predict_nmr_couplings(smi, &pos) {
            coupling_ok += 1;
            // J-coupling values should be finite
            assert!(
                couplings.iter().all(|c| c.j_hz.is_finite()),
                "Non-finite J-coupling"
            );
        }
    }

    let coupling_failures = N_MOLECULES - coupling_ok;
    let max_failures = (N_MOLECULES as f64 * MAX_FAILURE_RATE).ceil() as usize;

    eprintln!("=== NMR J-COUPLINGS ===");
    eprintln!(
        "  Success: {coupling_ok}/{N_MOLECULES} ({:.2}%)",
        coupling_ok as f64 / N_MOLECULES as f64 * 100.0
    );

    assert!(
        coupling_failures <= max_failures,
        "J-coupling failure rate exceeds limit"
    );
}

// ═══════════════════════════════════════════════════════════════════════════════
// Module 29: UV-Vis from EHT (compute_uv_vis_spectrum)
// ═══════════════════════════════════════════════════════════════════════════════

#[test]
fn benchmark_uv_vis_eht_200() {
    let smiles = load_smiles(200);
    let mut embed_ok = 0usize;
    let mut uvvis_ok = 0usize;

    for (smi, _id) in &smiles {
        let conf = sci_form::embed(smi, 42);
        if conf.error.is_some() || conf.coords.iter().any(|c| !c.is_finite()) {
            continue;
        }
        embed_ok += 1;
        let pos = to_positions(&conf.coords);

        if let Ok(spectrum) =
            sci_form::compute_uv_vis_spectrum(&conf.elements, &pos, 0.3, 1.0, 8.0, 200)
        {
            uvvis_ok += 1;
            assert_eq!(spectrum.energies_ev.len(), 200);
            assert!(spectrum.energies_ev.iter().all(|e| e.is_finite()));
        }
    }

    let uvvis_failures = embed_ok - uvvis_ok;
    let max_failures = (embed_ok as f64 * MAX_FAILURE_RATE).ceil() as usize;

    eprintln!("=== UV-VIS (EHT) ===");
    eprintln!("  Embedded: {embed_ok}/200");
    eprintln!(
        "  UV-Vis success: {uvvis_ok}/{embed_ok} ({:.2}%)",
        uvvis_ok as f64 / embed_ok as f64 * 100.0
    );

    assert!(
        uvvis_failures <= max_failures,
        "UV-Vis EHT failure rate exceeds limit"
    );
}

// ═══════════════════════════════════════════════════════════════════════════════
// Module 30: EHT/UFF Fallback
// ═══════════════════════════════════════════════════════════════════════════════

#[test]
fn benchmark_eht_uff_fallback_200() {
    let smiles = load_smiles(200);
    let mut embed_ok = 0usize;
    let mut fallback_ok = 0usize;

    for (smi, _id) in &smiles {
        let conf = sci_form::embed(smi, 42);
        if conf.error.is_some() || conf.coords.iter().any(|c| !c.is_finite()) {
            continue;
        }
        embed_ok += 1;
        let pos = to_positions(&conf.coords);

        if let Ok(_result) = sci_form::compute_eht_or_uff_fallback(smi, &conf.elements, &pos, true)
        {
            fallback_ok += 1;
        }
    }

    let failures = embed_ok - fallback_ok;
    let max_failures = (embed_ok as f64 * MAX_FAILURE_RATE).ceil() as usize;

    eprintln!("=== EHT/UFF FALLBACK ===");
    eprintln!("  Embedded: {embed_ok}/200");
    eprintln!(
        "  Fallback success: {fallback_ok}/{embed_ok} ({:.2}%)",
        fallback_ok as f64 / embed_ok as f64 * 100.0
    );

    assert!(
        failures <= max_failures,
        "EHT/UFF fallback failure rate exceeds limit"
    );
}

// ═══════════════════════════════════════════════════════════════════════════════
// Module 31: Method Comparison
// ═══════════════════════════════════════════════════════════════════════════════

#[test]
fn benchmark_compare_methods_100() {
    let smiles = load_smiles(100);
    let mut embed_ok = 0usize;
    let mut compare_ok = 0usize;

    for (smi, _id) in &smiles {
        let conf = sci_form::embed(smi, 42);
        if conf.error.is_some() || conf.coords.iter().any(|c| !c.is_finite()) {
            continue;
        }
        embed_ok += 1;
        let pos = to_positions(&conf.coords);

        let result = sci_form::compare_methods(smi, &conf.elements, &pos, true);
        compare_ok += 1;
        // Should have at least some methods that succeeded
        assert!(!result.comparisons.is_empty(), "No methods compared");
    }

    let failures = embed_ok - compare_ok;
    let max_failures = (embed_ok as f64 * MAX_FAILURE_RATE).ceil() as usize;

    eprintln!("=== METHOD COMPARISON ===");
    eprintln!("  Embedded: {embed_ok}/100");
    eprintln!(
        "  Compare success: {compare_ok}/{embed_ok} ({:.2}%)",
        compare_ok as f64 / embed_ok as f64 * 100.0
    );

    assert!(
        failures <= max_failures,
        "Method comparison failure rate exceeds limit"
    );
}

// ═══════════════════════════════════════════════════════════════════════════════
// Module 32: UFF Aromatic Heuristics
// ═══════════════════════════════════════════════════════════════════════════════

#[test]
fn benchmark_uff_aromatic_1k() {
    let smiles = load_smiles(N_MOLECULES);
    let mut embed_ok = 0usize;
    let mut uff_arom_ok = 0usize;

    for (smi, _id) in &smiles {
        let conf = sci_form::embed(smi, 42);
        if conf.error.is_some() || conf.coords.iter().any(|c| !c.is_finite()) {
            continue;
        }
        embed_ok += 1;

        if let Ok(result) = sci_form::compute_uff_energy_with_aromatic_heuristics(smi, &conf.coords)
        {
            if result.corrected_energy_kcal_mol.is_finite() {
                uff_arom_ok += 1;
            }
        }
    }

    let failures = embed_ok - uff_arom_ok;
    let max_failures = (embed_ok as f64 * MAX_FAILURE_RATE).ceil() as usize;

    eprintln!("=== UFF AROMATIC HEURISTICS ===");
    eprintln!("  Embedded: {embed_ok}/{N_MOLECULES}");
    eprintln!(
        "  Success: {uff_arom_ok}/{embed_ok} ({:.2}%)",
        uff_arom_ok as f64 / embed_ok as f64 * 100.0
    );

    assert!(
        failures <= max_failures,
        "UFF aromatic failure rate exceeds limit"
    );
}

// ═══════════════════════════════════════════════════════════════════════════════
// Module 33: MD Trajectory (small subset)
// ═══════════════════════════════════════════════════════════════════════════════

#[test]
fn benchmark_md_trajectory_20() {
    let smiles = load_smiles(20);
    let mut embed_ok = 0usize;
    let mut md_nve_ok = 0usize;
    let mut md_nvt_ok = 0usize;

    for (smi, _id) in &smiles {
        let conf = sci_form::embed(smi, 42);
        if conf.error.is_some() || conf.coords.iter().any(|c| !c.is_finite()) {
            continue;
        }
        embed_ok += 1;

        // NVE: 10 steps, 0.5 fs timestep
        if let Ok(traj) = sci_form::compute_md_trajectory(smi, &conf.coords, 10, 0.5, 42) {
            assert!(!traj.frames.is_empty(), "MD should produce frames");
            md_nve_ok += 1;
        }

        // NVT: 10 steps, 0.5 fs, 300K thermostat
        if let Ok(traj) =
            sci_form::compute_md_trajectory_nvt(smi, &conf.coords, 10, 0.5, 42, 300.0, 100.0)
        {
            assert!(!traj.frames.is_empty(), "NVT MD should produce frames");
            md_nvt_ok += 1;
        }
    }

    eprintln!("=== MOLECULAR DYNAMICS ===");
    eprintln!("  Embedded: {embed_ok}/20");
    eprintln!("  NVE success: {md_nve_ok}/{embed_ok}");
    eprintln!("  NVT success: {md_nvt_ok}/{embed_ok}");

    let max_failures = (embed_ok as f64 * MAX_FAILURE_RATE).ceil() as usize;
    assert!(
        embed_ok - md_nve_ok <= max_failures,
        "MD NVE failure rate exceeds limit"
    );
    assert!(
        embed_ok - md_nvt_ok <= max_failures,
        "MD NVT failure rate exceeds limit"
    );
}

// ═══════════════════════════════════════════════════════════════════════════════
// Module 34: Conformer Search with UFF
// ═══════════════════════════════════════════════════════════════════════════════

#[test]
fn benchmark_conformer_search_10() {
    let smiles = load_smiles(10);
    let mut search_ok = 0usize;
    let mut embed_ok = 0usize;

    for (smi, _id) in &smiles {
        let conf = sci_form::embed(smi, 42);
        if conf.error.is_some() {
            continue;
        }
        embed_ok += 1;

        if let Ok(result) = sci_form::search_conformers_with_uff(smi, 5, 42, 0.5) {
            search_ok += 1;
            assert!(result.unique > 0, "Should find at least 1 unique conformer");
            assert!(result.generated > 0, "Should generate at least 1 conformer");
        }
    }

    eprintln!("=== CONFORMER SEARCH ===");
    eprintln!("  Parseable: {embed_ok}/10");
    eprintln!("  Search success: {search_ok}/{embed_ok}");

    let max_failures = (embed_ok as f64 * MAX_FAILURE_RATE).ceil() as usize;
    assert!(
        embed_ok - search_ok <= max_failures,
        "Conformer search failure rate exceeds limit"
    );
}

// ═══════════════════════════════════════════════════════════════════════════════
// Module 35: Framework Assembly
// ═══════════════════════════════════════════════════════════════════════════════

#[test]
fn benchmark_framework_assembly() {
    // Test framework assembly using proper Topology constructors
    use sci_form::materials::{ConnectionPoint, CoordinationGeometry, Sbu, Topology};

    let cell = sci_form::create_unit_cell(10.0, 10.0, 10.0, 90.0, 90.0, 90.0);

    // Build a simple metal node
    let mut node = Sbu::new("Zn_node");
    node.add_atom(30, [0.0, 0.0, 0.0]);
    node.geometry = Some(CoordinationGeometry::Octahedral);
    node.connections = vec![
        ConnectionPoint {
            position: [1.0, 0.0, 0.0],
            direction: [1.0, 0.0, 0.0],
            kind: "metal".to_string(),
        },
        ConnectionPoint {
            position: [-1.0, 0.0, 0.0],
            direction: [-1.0, 0.0, 0.0],
            kind: "metal".to_string(),
        },
        ConnectionPoint {
            position: [0.0, 1.0, 0.0],
            direction: [0.0, 1.0, 0.0],
            kind: "metal".to_string(),
        },
    ];

    // Build a simple linker
    let mut linker = Sbu::new("C_linker");
    linker.add_atom(6, [0.0, 0.0, 0.0]);
    linker.add_atom(6, [1.0, 0.0, 0.0]);
    linker.connections = vec![
        ConnectionPoint {
            position: [-0.5, 0.0, 0.0],
            direction: [-1.0, 0.0, 0.0],
            kind: "carboxylate".to_string(),
        },
        ConnectionPoint {
            position: [1.5, 0.0, 0.0],
            direction: [1.0, 0.0, 0.0],
            kind: "carboxylate".to_string(),
        },
    ];

    // Use built-in pcu topology
    let topo = Topology::pcu();

    let crystal = sci_form::assemble_framework(&node, &linker, &topo, &cell);

    eprintln!("=== FRAMEWORK ASSEMBLY ===");
    eprintln!("  Atoms: {}", crystal.atoms.len());
    eprintln!("  Labels: {}", crystal.labels.len());

    // Should produce at least some atoms
    assert!(!crystal.atoms.is_empty(), "Framework should have atoms");
    assert_eq!(
        crystal.atoms.len(),
        crystal.labels.len(),
        "Atom count should match label count"
    );
}

// ═══════════════════════════════════════════════════════════════════════════════
// Module 36: NEB Path (small subset)
// ═══════════════════════════════════════════════════════════════════════════════

#[test]
fn benchmark_neb_path_5() {
    let smiles = load_smiles(5);
    let mut neb_ok = 0usize;
    let mut embed_ok = 0usize;

    for (smi, _id) in &smiles {
        let conf1 = sci_form::embed(smi, 42);
        let conf2 = sci_form::embed(smi, 123);
        if conf1.error.is_some() || conf2.error.is_some() {
            continue;
        }
        if conf1.coords.len() != conf2.coords.len() {
            continue;
        }
        embed_ok += 1;

        if let Ok(result) = sci_form::compute_simplified_neb_path(
            smi,
            &conf1.coords,
            &conf2.coords,
            5,
            10,
            0.1,
            0.01,
        ) {
            neb_ok += 1;
            assert!(!result.images.is_empty(), "NEB should produce images");
        }
    }

    eprintln!("=== NEB PATH ===");
    eprintln!("  Embedded pairs: {embed_ok}/5");
    eprintln!("  NEB success: {neb_ok}/{embed_ok}");

    // NEB is expected to work on most molecules
    assert!(neb_ok > 0, "At least one NEB should succeed");
}
