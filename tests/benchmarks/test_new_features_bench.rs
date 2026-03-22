//! Benchmarks for newly implemented features.
//!
//! Measures throughput and validates correctness on a subset of molecules:
//! - SMARTS batch matching (sequential vs parallel)
//! - 3D descriptors computation
//! - Pharmacophore fingerprints
//! - XRD simulation
//! - Space group symmetry expansion
//! - Porosity calculation
//! - ML uncertainty estimation
//! - Semi-analytical Hessian (small molecules)
//!
//! Run: cargo test --test benchmarks test_bench_new_features --release -- --nocapture

use std::time::Instant;

fn embed_or_skip(smiles: &str) -> Option<sci_form::ConformerResult> {
    let r = sci_form::embed(smiles, 42);
    if r.error.is_none() {
        Some(r)
    } else {
        None
    }
}

fn _to_pos(coords: &[f64]) -> Vec<[f64; 3]> {
    coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect()
}

// ─── SMARTS batch matching benchmark ──────────────────────────────────────────

#[test]
fn test_bench_smarts_batch_matching() {
    let smiles_list = [
        "CCO",
        "c1ccccc1",
        "CC(=O)O",
        "CC(=O)Nc1ccc(O)cc1",
        "c1ccc(cc1)C(=O)O",
        "C1CCCCC1",
        "CC(C)O",
        "CCCCCC",
        "c1ccncc1",
        "c1ccc2ccccc2c1",
        "OC(=O)c1ccccc1O",
        "CC(=O)OC1=CC=CC=C1C(=O)O",
        "c1ccc(NC(=O)C)cc1",
        "CCN(CC)CC",
        "CC(C)(C)O",
        "OCCO",
        "c1ccc(F)cc1",
        "CC=O",
        "CC(=O)c1ccccc1",
        "OC(=O)CCCCC(=O)O",
    ];

    let mols: Vec<sci_form::graph::Molecule> = smiles_list
        .iter()
        .filter_map(|s| sci_form::graph::Molecule::from_smiles(s).ok())
        .collect();
    let mol_refs: Vec<&sci_form::graph::Molecule> = mols.iter().collect();
    let pattern = sci_form::smarts::parse_smarts("[OX2H]").unwrap();

    // Sequential
    let t0 = Instant::now();
    let n_iters = 100;
    for _ in 0..n_iters {
        let _ = sci_form::smarts::substruct_match_batch(&mol_refs, &pattern);
    }
    let seq_us = t0.elapsed().as_micros() as f64 / n_iters as f64;

    println!(
        "SMARTS batch ({}×{} mols, sequential): {:.1} µs/batch",
        n_iters,
        mol_refs.len(),
        seq_us,
    );

    // Validate results
    let results = sci_form::smarts::substruct_match_batch(&mol_refs, &pattern);
    let matches = results.iter().filter(|r| !r.is_empty()).count();
    println!("  Matches: {} / {} molecules", matches, mol_refs.len());
    assert!(matches > 0, "Should find some OH matches");
}

// ─── 3D descriptors benchmark ─────────────────────────────────────────────────

#[test]
fn test_bench_3d_descriptors() {
    let smiles_list = [
        "CCO",
        "c1ccccc1",
        "CC(=O)O",
        "CC(=O)Nc1ccc(O)cc1",
        "c1ccc(cc1)C(=O)O",
        "C1CCCCC1",
        "CCCCCCCC",
        "c1ccncc1",
        "c1ccc2ccccc2c1",
        "OC(=O)CCCCC(=O)O",
    ];

    let conformers: Vec<sci_form::ConformerResult> = smiles_list
        .iter()
        .filter_map(|s| embed_or_skip(s))
        .collect();

    let t0 = Instant::now();
    let n_iters = 1000;
    for _ in 0..n_iters {
        for conf in &conformers {
            let _ = sci_form::ml::compute_3d_descriptors(&conf.elements, &conf.coords);
        }
    }
    let total_us = t0.elapsed().as_micros();
    let per_mol_us = total_us as f64 / (n_iters as f64 * conformers.len() as f64);

    println!(
        "3D descriptors: {:.2} µs/molecule ({}×{} mols)",
        per_mol_us,
        n_iters,
        conformers.len(),
    );

    // Validate one
    let desc = sci_form::ml::compute_3d_descriptors(&conformers[0].elements, &conformers[0].coords);
    assert!(desc.radius_of_gyration > 0.0);
    assert!(desc.span > 0.0);
}

// ─── Pharmacophore fingerprint benchmark ──────────────────────────────────────

#[test]
fn test_bench_pharmacophore_fingerprints() {
    use sci_form::ml::pharmacophore::{
        compute_pharmacophore_fingerprint, detect_features, pharmacophore_tanimoto,
    };

    let test_mols = [
        (
            &[6u8, 6, 8, 1, 1, 1, 1, 1, 1] as &[u8],
            vec![
                (0, 1, 1u8),
                (1, 2, 1),
                (0, 3, 1),
                (0, 4, 1),
                (0, 5, 1),
                (1, 6, 1),
                (1, 7, 1),
                (2, 8, 1),
            ],
        ),
        (
            &[6, 6, 6, 8, 1, 1, 1, 1, 1, 1, 1, 1],
            vec![
                (0, 1, 1),
                (1, 2, 1),
                (2, 3, 1),
                (0, 4, 1),
                (0, 5, 1),
                (0, 6, 1),
                (1, 7, 1),
                (1, 8, 1),
                (2, 9, 1),
                (2, 10, 1),
                (3, 11, 1),
            ],
        ),
    ];

    let t0 = Instant::now();
    let n_iters = 500;
    let mut fps = Vec::new();
    for _ in 0..n_iters {
        fps.clear();
        for (elems, bonds) in &test_mols {
            let features = detect_features(elems, bonds, &[], &[]);
            let fp = compute_pharmacophore_fingerprint(&features, 2048);
            fps.push(fp);
        }
    }
    let total_us = t0.elapsed().as_micros();
    let per_mol_us = total_us as f64 / (n_iters as f64 * test_mols.len() as f64);

    println!(
        "Pharmacophore FP: {:.2} µs/molecule ({}×{} mols)",
        per_mol_us,
        n_iters,
        test_mols.len(),
    );

    // Benchmark Tanimoto comparison
    let features1 = detect_features(test_mols[0].0, &test_mols[0].1, &[], &[]);
    let features2 = detect_features(test_mols[1].0, &test_mols[1].1, &[], &[]);
    let fp1 = compute_pharmacophore_fingerprint(&features1, 2048);
    let fp2 = compute_pharmacophore_fingerprint(&features2, 2048);

    let t0 = Instant::now();
    let n_tani = 10000;
    for _ in 0..n_tani {
        let _ = pharmacophore_tanimoto(&fp1, &fp2);
    }
    let tani_ns = t0.elapsed().as_nanos() as f64 / n_tani as f64;
    println!("Pharmacophore Tanimoto: {:.1} ns/comparison", tani_ns);

    let sim = pharmacophore_tanimoto(&fp1, &fp2);
    assert!(sim > 0.0 && sim <= 1.0);
}

// ─── Powder XRD benchmark ─────────────────────────────────────────────────────

#[test]
fn test_bench_powder_xrd() {
    use sci_form::materials::{simulate_powder_xrd, UnitCell};

    let cell = UnitCell::cubic(5.64);
    let elements = [11u8, 17, 11, 17, 11, 17, 11, 17];
    let frac_coords = [
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.5],
        [0.5, 0.5, 0.0],
        [0.0, 0.0, 0.5],
        [0.5, 0.0, 0.5],
        [0.0, 0.5, 0.0],
        [0.0, 0.5, 0.5],
        [0.5, 0.0, 0.0],
    ];

    let t0 = Instant::now();
    let n_iters = 50;
    for _ in 0..n_iters {
        let _ = simulate_powder_xrd(&cell, &elements, &frac_coords, 90.0);
    }
    let per_ms = t0.elapsed().as_millis() as f64 / n_iters as f64;
    println!("Powder XRD (8 atoms, 90°): {:.2} ms/simulation", per_ms);

    let xrd = simulate_powder_xrd(&cell, &elements, &frac_coords, 90.0);
    println!("  {} reflections", xrd.two_theta.len());
    assert!(!xrd.two_theta.is_empty());
}

// ─── Space group symmetry benchmark ───────────────────────────────────────────

#[test]
fn test_bench_symmetry_expansion() {
    use sci_form::materials::{expand_by_symmetry, get_space_group};

    let sg = get_space_group(14).unwrap(); // P2_1/c — 4 operations

    // 10 atoms in asymmetric unit
    let n_asym = 10;
    let asym_frac: Vec<[f64; 3]> = (0..n_asym)
        .map(|i| {
            let f = i as f64 / n_asym as f64;
            [f * 0.3 + 0.1, f * 0.2 + 0.05, f * 0.4 + 0.1]
        })
        .collect();
    let asym_elements: Vec<u8> = (0..n_asym)
        .map(|i| if i % 2 == 0 { 6 } else { 8 })
        .collect();

    let t0 = Instant::now();
    let n_iters = 1000;
    for _ in 0..n_iters {
        let _ = expand_by_symmetry(&sg, &asym_frac, &asym_elements);
    }
    let per_us = t0.elapsed().as_micros() as f64 / n_iters as f64;
    println!(
        "Symmetry expansion (P2_1/c, {} asym atoms): {:.2} µs/expansion",
        n_asym, per_us
    );

    let (expanded, _) = expand_by_symmetry(&sg, &asym_frac, &asym_elements);
    println!("  {} → {} atoms", n_asym, expanded.len());
    assert!(expanded.len() >= n_asym);
}

// ─── Porosity benchmark ──────────────────────────────────────────────────────

#[test]
fn test_bench_porosity() {
    use sci_form::materials::{compute_porosity, UnitCell};

    let cell = UnitCell::cubic(15.0);
    let elements = [30u8, 8, 8, 6, 6, 30, 8, 8]; // Simple MOF-like
    let frac_coords = [
        [0.0, 0.0, 0.0],
        [0.1, 0.0, 0.0],
        [-0.1, 0.0, 0.0],
        [0.25, 0.25, 0.0],
        [0.5, 0.5, 0.0],
        [0.5, 0.5, 0.5],
        [0.6, 0.5, 0.5],
        [0.4, 0.5, 0.5],
    ];

    let t0 = Instant::now();
    let n_iters = 10;
    for _ in 0..n_iters {
        let _ = compute_porosity(&cell, &elements, &frac_coords, 1.4, 1.0);
    }
    let per_ms = t0.elapsed().as_millis() as f64 / n_iters as f64;
    println!("Porosity (15Å cubic, 8 atoms, 1.0Å grid): {:.2} ms", per_ms);

    let result = compute_porosity(&cell, &elements, &frac_coords, 1.4, 1.0);
    println!(
        "  Porosity: {:.1}%, pore vol: {:.1} ų, LCD: {:.1} Å",
        result.porosity * 100.0,
        result.pore_volume,
        result.largest_cavity_diameter
    );
    assert!(result.porosity > 0.0);
}

// ─── ML uncertainty benchmark ─────────────────────────────────────────────────

#[test]
fn test_bench_ml_prediction_with_uncertainty() {
    use petgraph::visit::EdgeRef;

    let smiles_list = [
        "CCO",
        "c1ccccc1",
        "CC(=O)O",
        "CCCCCC",
        "c1ccncc1",
        "CC(C)O",
        "OCCO",
        "c1ccc(F)cc1",
        "CC=O",
        "CCC(=O)O",
    ];

    let t0 = Instant::now();
    let n_iters = 200;
    for _ in 0..n_iters {
        for smi in &smiles_list {
            let mol = sci_form::parse(smi).unwrap();
            let bonds: Vec<(usize, usize, u8)> = mol
                .graph
                .edge_references()
                .map(|e| {
                    let order = match e.weight().order {
                        sci_form::graph::BondOrder::Single => 1,
                        sci_form::graph::BondOrder::Double => 2,
                        sci_form::graph::BondOrder::Triple => 3,
                        _ => 0,
                    };
                    (e.source().index(), e.target().index(), order)
                })
                .collect();
            let elements: Vec<u8> = (0..mol.graph.node_count())
                .map(|i| mol.graph[petgraph::graph::NodeIndex::new(i)].element)
                .collect();
            let desc = sci_form::compute_ml_descriptors(&elements, &bonds, &[], &[]);
            let _ = sci_form::predict_ml_properties(&desc);
        }
    }
    let total_us = t0.elapsed().as_micros();
    let per_mol_us = total_us as f64 / (n_iters as f64 * smiles_list.len() as f64);

    println!(
        "ML predict + uncertainty: {:.2} µs/molecule ({}×{} mols)",
        per_mol_us,
        n_iters,
        smiles_list.len(),
    );
}

// ─── Semi-analytical Hessian benchmark ────────────────────────────────────────

#[test]
fn test_bench_semianalytical_hessian_pm3() {
    let elements = [8u8, 1, 1]; // Water
    let positions = [[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]];

    // Semi-analytical (PM3 gradient-based) — should be faster
    let t0 = Instant::now();
    let n_iters = 3;
    for _ in 0..n_iters {
        let _ = sci_form::ir::hessian::compute_semianalytical_hessian(
            &elements,
            &positions,
            sci_form::ir::hessian::HessianMethod::Pm3,
            None,
        );
    }
    let semi_ms = t0.elapsed().as_millis() as f64 / n_iters as f64;

    // Numerical (PM3 energy-based) — for comparison
    let t0 = Instant::now();
    for _ in 0..n_iters {
        let _ = sci_form::ir::hessian::compute_numerical_hessian(
            &elements,
            &positions,
            sci_form::ir::hessian::HessianMethod::Pm3,
            None,
        );
    }
    let numer_ms = t0.elapsed().as_millis() as f64 / n_iters as f64;

    println!(
        "Hessian (water, PM3):\n  Semi-analytical: {:.1} ms\n  Numerical:       {:.1} ms\n  Speedup:         {:.1}x",
        semi_ms,
        numer_ms,
        numer_ms / semi_ms.max(0.1)
    );
}

// ─── MMFF94 vs MMFF94s benchmark ──────────────────────────────────────────────

#[test]
fn test_bench_mmff94_vs_mmff94s() {
    use sci_form::forcefield::mmff94::{Mmff94Builder, Mmff94Variant};

    let elements = [6u8, 6, 8, 1, 1, 1, 1, 1, 1]; // ethanol
    let bonds = vec![
        (0, 1, 1u8),
        (1, 2, 1),
        (0, 3, 1),
        (0, 4, 1),
        (0, 5, 1),
        (1, 6, 1),
        (1, 7, 1),
        (2, 8, 1),
    ];

    let t0 = Instant::now();
    let n_iters = 5000;
    for _ in 0..n_iters {
        let _ = Mmff94Builder::build(&elements, &bonds);
    }
    let mmff94_us = t0.elapsed().as_micros() as f64 / n_iters as f64;

    let t0 = Instant::now();
    for _ in 0..n_iters {
        let _ = Mmff94Builder::build_variant(&elements, &bonds, Mmff94Variant::Mmff94s);
    }
    let mmff94s_us = t0.elapsed().as_micros() as f64 / n_iters as f64;

    println!(
        "MMFF94 build (ethanol): {:.2} µs\nMMFF94s build:          {:.2} µs",
        mmff94_us, mmff94s_us
    );
}

// ─── Stereo analysis with new features benchmark ──────────────────────────────

#[test]
fn test_bench_stereo_with_atropisomerism() {
    let smiles_list = [
        "C[C@H](F)Cl",
        "CCO",
        "c1ccccc1",
        "C(F)(Cl)(Br)I",
        "CC(=O)Nc1ccc(O)cc1",
    ];

    let t0 = Instant::now();
    let n_iters = 100;
    for _ in 0..n_iters {
        for smi in &smiles_list {
            let _ = sci_form::analyze_stereo(smi, &[]);
        }
    }
    let per_mol_us = t0.elapsed().as_micros() as f64 / (n_iters as f64 * smiles_list.len() as f64);

    println!(
        "Stereo analysis (with atropis/prochiral): {:.2} µs/molecule",
        per_mol_us
    );
}
