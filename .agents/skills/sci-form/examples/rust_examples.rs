// sci-form Rust Examples
// All major API groups demonstrated with real usage patterns.

use sci_form::{
    embed, embed_batch, parse, ConformerConfig,
    eht, alignment,
    compute_charges, compute_charges_configured, compute_sasa, compute_dipole, compute_esp,
    compute_population, compute_bond_orders,
    compute_frontier_descriptors, compute_fukui_descriptors,
    compute_reactivity_ranking, compute_empirical_pka,
    compute_uv_vis_spectrum,
    compute_stda_uvvis, reactivity::BroadeningType,
    compute_vibrational_analysis, compute_ir_spectrum, compute_ir_spectrum_broadened,
    predict_nmr_shifts, predict_nmr_couplings, compute_nmr_spectrum, compute_hose_codes,
    compute_ensemble_j_couplings,
    compute_pm3, compute_xtb,
    compute_ml_descriptors, predict_ml_properties,
    compute_uff_energy, compute_uff_energy_with_aromatic_heuristics, compute_mmff94_energy,
    analyze_graph_features, compute_topology,
    create_unit_cell, assemble_framework, materials,
    transport,
    get_system_capabilities, get_system_method_plan, compare_methods,
    analyze_stereo,
    compute_nonpolar_solvation, compute_gb_solvation,
    compute_sssr, compute_ecfp, compute_tanimoto,
    butina_cluster, compute_rmsd_matrix,
};

fn water_system() -> (Vec<u8>, Vec<[f64; 3]>, Vec<f64>) {
    let elements = vec![8u8, 1, 1];
    let positions = vec![[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]];
    let flat: Vec<f64> = positions.iter().flat_map(|p| p.iter().copied()).collect();
    (elements, positions, flat)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {

    // ── 1. Geometry / Embedding ─────────────────────────────────────────────

    // Single conformer
    let result = embed("C1=CC=CC=C1", 42);
    assert!(result.error.is_none(), "embed failed: {:?}", result.error);
    println!("Benzene: {} atoms, coords.len={}", result.num_atoms, result.coords.len());

    // Parallel batch (num_threads=0 = all cores)
    let smiles = vec!["O", "CCO", "c1ccccc1", "CC(=O)O"];
    let config = ConformerConfig { seed: 42, num_threads: 0 };
    let batch = embed_batch(&smiles, &config);
    println!("Batch: {}/{} succeeded", batch.iter().filter(|r| r.error.is_none()).count(), batch.len());

    // Topology parse
    let mol = parse("CC(=O)N")?;
    println!("Acetamide: {} atoms, {} bonds", mol.graph.node_count(), mol.graph.edge_count());

    // RMSD / Kabsch alignment
    let r1 = embed("C", 1);
    let r2 = embed("C", 2);
    let aligned = alignment::align_coordinates(&r1.coords, &r2.coords);
    println!("Methane RMSD: {:.4} Å", aligned.rmsd);

    // ── 2. Force Fields ─────────────────────────────────────────────────────
    let benz = embed("c1ccccc1", 42);
    let e_uff = compute_uff_energy("c1ccccc1", &benz.coords)?;
    println!("UFF energy benzene: {:.3} kcal/mol", e_uff);

    let e_mmff = compute_mmff94_energy("c1ccccc1", &benz.coords)?;
    println!("MMFF94 energy: {:.3} kcal/mol", e_mmff);

    let e_heur = compute_uff_energy_with_aromatic_heuristics("c1ccccc1", &benz.coords)?;
    println!("UFF corrected: {:.3} kcal/mol (aromatic bonds: {})",
        e_heur.corrected_energy_kcal_mol, e_heur.aromatic_bond_count);

    // ── 3. EHT / Quantum ────────────────────────────────────────────────────
    let (elements, positions, _flat) = water_system();
    let eht = eht::solve_eht(&elements, &positions, None)?;
    println!("H2O EHT: HOMO={:.3} eV LUMO={:.3} eV gap={:.3} eV",
        eht.homo_energy, eht.lumo_energy, eht.gap);

    // Orbital isosurface mesh
    let mesh = eht::marching_cubes(
        &eht::evaluate_orbital_on_grid(&eht::basis::build_basis(&elements, &positions),
            &eht.coefficients, eht.homo_index, &positions, 0.2, 3.0),
        0.02
    );
    println!("HOMO mesh: {} triangles", mesh.num_triangles);

    // DOS
    let dos = sci_form::compute_dos(&elements, &positions, 0.3, -30.0, 5.0, 500)?;
    println!("DOS: {} points, sigma={}", dos.energies.len(), dos.sigma);

    // ── 4. Charges & Surface ────────────────────────────────────────────────
    let charges = compute_charges("CC(=O)O")?;
    println!("Charges: {:#.3?}, total={:.4}", &charges.charges[..3], charges.total_charge);

    let (elements, _ps, flat) = water_system();
    let sasa = compute_sasa(&elements, &flat, None)?;
    println!("SASA: {:.2} Ų", sasa.total_sasa);

    let (elements, positions, _) = water_system();
    let dipole = compute_dipole(&elements, &positions)?;
    println!("Dipole: {:.3} D", dipole.magnitude);

    let esp = compute_esp(&elements, &positions, 0.5, 3.0)?;
    println!("ESP grid: {}x{}x{} = {} values", esp.dims[0], esp.dims[1], esp.dims[2], esp.values.len());

    // ── 5. Population & Bond Orders ─────────────────────────────────────────
    let (elements, positions, _) = water_system();
    let pop = compute_population(&elements, &positions)?;
    println!("Mulliken total: {:.4}", pop.total_charge_mulliken);

    let bo = compute_bond_orders(&elements, &positions)?;
    for (i, (a, b)) in bo.bonds.iter().take(3).enumerate() {
        println!("Bond {}-{}: Wiberg={:.3} Mayer={:.3}", a.atom_i, a.atom_j, a.wiberg, a.mayer);
    }

    // ── 6. Reactivity ───────────────────────────────────────────────────────
    let phenol = embed("Oc1ccccc1", 42);
    let (el, pos) = extract_elpos(&phenol);
    let frontier = compute_frontier_descriptors(&el, &pos)?;
    println!("Frontier gap: {:.3} eV", frontier.gap);

    let fukui = compute_fukui_descriptors(&el, &pos)?;
    println!("Fukui f+ max at atom {}", fukui.condensed.iter()
        .enumerate().max_by(|a, b| a.1.f_plus.partial_cmp(&b.1.f_plus).unwrap())
        .map(|(i, _)| i).unwrap_or(0));

    let reactml = compute_reactivity_ranking(&el, &pos)?;
    println!("Top nucleophilic: atom {}", reactml.nucleophilic_attack_sites.first()
        .map(|s| s.atom_index).unwrap_or(0));

    let pka = compute_empirical_pka("CC(=O)O")?;
    println!("pKa acidic sites: {}", pka.acidic_sites.len());

    // ── 7. UV-Vis Spectroscopy ──────────────────────────────────────────────
    let (el, pos, _) = water_system();
    let uvvis = compute_uv_vis_spectrum(&el, &pos, 0.25, 0.5, 8.0, 600)?;
    println!("UV-Vis: {} points, {} peaks", uvvis.energies_ev.len(), uvvis.peaks.len());

    let stda = compute_stda_uvvis(&el, &pos, 0.3, 1.0, 8.0, 500, BroadeningType::Gaussian)?;
    println!("sTDA: {} excitations", stda.excitations.len());

    // ── 8. IR Spectroscopy ──────────────────────────────────────────────────
    let (el, pos, _) = water_system();
    let vib = compute_vibrational_analysis(&el, &pos, "eht", None)?;
    println!("Vibrational: {} modes, ZPVE={:.4} eV", vib.modes.len(), vib.zpve_ev);

    let ir = compute_ir_spectrum(&vib, 15.0, 400.0, 4000.0, 1000);
    println!("IR spectrum: {} points, {} peaks", ir.wavenumbers.len(), ir.peaks.len());

    // ── 9. NMR Spectroscopy ─────────────────────────────────────────────────
    let shifts = predict_nmr_shifts("CC(=O)O")?;
    for s in &shifts.h_shifts { println!("H atom {}: {:.2} ppm ({})", s.atom_index, s.shift_ppm, s.environment); }
    for s in &shifts.c_shifts { println!("C atom {}: {:.2} ppm", s.atom_index, s.shift_ppm); }

    let couplings = predict_nmr_couplings("CCC", &[])?;
    for c in &couplings { println!("J({}-{}): {:.1} Hz [{}]", c.h1_index, c.h2_index, c.j_hz, c.coupling_type); }

    let nmr = compute_nmr_spectrum("CC(=O)O", "1H", 0.02, 0.0, 12.0, 1000)?;
    println!("1H NMR: {} points, {} peaks", nmr.ppm_axis.len(), nmr.peaks.len());

    let hose = compute_hose_codes("c1ccccc1", 2)?;
    for (idx, el, code) in &hose { println!("Atom {}: {}", idx, code); }

    // ── 10. PM3 and xTB ─────────────────────────────────────────────────────
    let (el, pos, _) = water_system();
    let pm3 = compute_pm3(&el, &pos)?;
    println!("PM3: Hf={:.3} kcal/mol, gap={:.3} eV, converged={}", pm3.heat_of_formation, pm3.gap, pm3.converged);

    let xtb = compute_xtb(&el, &pos)?;
    println!("xTB: E={:.6} eV, gap={:.3} eV, converged={}", xtb.total_energy, xtb.gap, xtb.converged);

    // ── 11. ML Properties ───────────────────────────────────────────────────
    let desc = {
        let mol = parse("Cc1ccccc1")?;
        let n = mol.graph.node_count();
        let elems: Vec<u8> = (0..n).map(|i| mol.graph[sci_form::graph::NodeIndex::new(i)].element).collect();
        let bonds: Vec<(usize, usize, u8)> = mol.graph.edge_indices().map(|e| {
            let (a, b) = mol.graph.edge_endpoints(e).unwrap();
            (a.index(), b.index(), 1u8)
        }).collect();
        compute_ml_descriptors(&elems, &bonds, &[], &[])
    };
    println!("Toluene MW={:.2}, HBA={}, HBD={}, LogP≈{:.2}", desc.molecular_weight, desc.n_hba, desc.n_hbd, predict_ml_properties(&desc).logp);

    // ── 12. Topology & Graph ─────────────────────────────────────────────────
    let gf = analyze_graph_features("c1ccccc1")?;
    println!("Aromatic atoms: {}", gf.aromaticity.aromatic_atoms.iter().filter(|&&a| a).count());

    // ── 13. Materials ────────────────────────────────────────────────────────
    let cell = create_unit_cell(10.0, 10.0, 10.0, 90.0, 90.0, 90.0);
    println!("Cell volume: {:.2} ų", cell.volume());

    let node = materials::Sbu::metal_node(30, 0.0, materials::CoordinationGeometry::Octahedral);
    let linker = materials::Sbu::linear_linker(&[6, 6], 1.4, "carboxylate");
    let topo = materials::Topology::pcu();
    let framework = assemble_framework(&node, &linker, &topo, &cell);
    println!("MOF: {} atoms", framework.num_atoms());

    // ── 14. Transport ────────────────────────────────────────────────────────
    let all_smiles = vec!["O".to_string(), "CCO".to_string(), "c1ccccc1".to_string(), "CC(=O)O".to_string()];
    let n = transport::estimate_workers(all_smiles.len(), 8);
    let tasks = transport::split_batch(&all_smiles, n, 42);
    println!("Transport: {} workers, {} tasks", n, tasks.len());

    // ── 15. System Planning ──────────────────────────────────────────────────
    let (el, _, _) = water_system();
    let caps = get_system_capabilities(&el);
    println!("H2O EHT available: {}", caps.eht.available);

    let plan = get_system_method_plan(&el);
    println!("Recommended orbitals method: {:?}", plan.orbitals.recommended);

    let cmp = compare_methods("O", &el, &[[0.0,0.0,0.0],[0.96,0.0,0.0],[-0.24,0.93,0.0]], false);
    for entry in &cmp.comparisons { println!("Method {:?}: status={:?}", entry.method, entry.status); }

    // ── 16. Stereochemistry ────────────────────────────────────────────────────
    let chiral = embed("C(F)(Cl)(Br)I", 42);
    let stereo = analyze_stereo("C(F)(Cl)(Br)I", &chiral.coords)?;
    println!("Stereocenters: {}, double bonds: {}", stereo.n_stereocenters, stereo.n_double_bonds);
    for sc in &stereo.stereocenters {
        println!("  atom {}: config={:?}  priorities={:?}", sc.atom_index, sc.configuration, sc.priorities);
    }
    let alk = embed("CC=CC", 42);
    let ez = analyze_stereo("CC=CC", &alk.coords)?;
    for db in &ez.double_bonds {
        println!("  double bond {}-{}: config={:?}", db.atom1, db.atom2, db.configuration);
    }

    // ── 17. Solvation ───────────────────────────────────────────────────────────
    let (el, pos, _) = water_system();
    let nonpolar = compute_nonpolar_solvation(&el, &pos, None);
    println!("Non-polar solvation: {:.4} kcal/mol  SASA={:.2} Å²", nonpolar.energy_kcal_mol, nonpolar.total_sasa);

    let charges_w = compute_charges("O")?;
    let gb = compute_gb_solvation(&el, &pos, &charges_w.charges, None, None, None);
    println!("GB electrostatic: {:.4} kcal/mol  total: {:.4} kcal/mol",
        gb.electrostatic_energy_kcal_mol, gb.total_energy_kcal_mol);

    // ── 18. Rings & Fingerprints ─────────────────────────────────────────────────
    let sssr = compute_sssr("c1ccccc1")?;
    println!("Benzene ring histogram: {:?}", sssr.ring_size_histogram);

    let fp_benz = compute_ecfp("c1ccccc1", 2, 2048)?;
    let fp_tol  = compute_ecfp("Cc1ccccc1", 2, 2048)?;
    println!("Benzene/Toluene Tanimoto: {:.4}", compute_tanimoto(&fp_benz, &fp_tol));
    println!("Benzene on-bits: {}", fp_benz.on_bits.len());

    // Multi-ring
    let sssr_naph = compute_sssr("c1ccc2ccccc2c1")?; // naphthalene
    println!("Naphthalene rings: {}", sssr_naph.rings.len());

    // ── 19. Clustering ───────────────────────────────────────────────────────────
    let smiles_set = ["c1ccccc1", "Cc1ccccc1", "CC(=O)O", "CCO"];
    let conformers: Vec<Vec<f64>> = smiles_set.iter()
        .map(|s| embed(s, 42).coords)
        .collect();
    let clusters = butina_cluster(&conformers, 2.0);
    println!("Clusters (cutoff=2.0Å): {} total", clusters.n_clusters);
    println!("Assignments: {:?}", clusters.assignments);

    let rmsd_mat = compute_rmsd_matrix(&conformers);
    println!("RMSD matrix: {}x{}", rmsd_mat.len(), rmsd_mat[0].len());
    for (i, row) in rmsd_mat.iter().enumerate() {
        println!("  row {}: {:?}", i, row.iter().map(|v| format!("{:.3}", v)).collect::<Vec<_>>());
    }

    // ── NMR ensemble J-couplings ───────────────────────────────────────────────
    let conf_coords: Vec<Vec<f64>> = (0..3u64).map(|s| embed("CCC", s).coords).collect();
    let energies_kcal = vec![0.0, 0.5, 1.2];
    let ens_j = compute_ensemble_j_couplings("CCC", &conf_coords, &energies_kcal, 298.15)?;
    println!("Ensemble J-couplings: {} pairs", ens_j.len());
    for c in &ens_j { println!("  J({},{}) = {:.2} Hz  {}", c.h1_index, c.h2_index, c.j_hz, c.coupling_type); }

    // ── IR broadened (Gaussian) ──────────────────────────────────────────────
    let (el, pos, _) = water_system();
    let vib2 = compute_vibrational_analysis(&el, &pos, "eht", None)?;
    let ir_gauss = compute_ir_spectrum_broadened(&vib2, 20.0, 400.0, 4000.0, 1000, "gaussian");
    println!("IR Gaussian: {} points, {} peaks", ir_gauss.wavenumbers.len(), ir_gauss.peaks.len());

    // ── Charges configured ─────────────────────────────────────────────────
    let cfg = sci_form::charges::gasteiger::GasteigerConfig {
        max_iter: 10,
        initial_damping: 0.4,
        convergence_threshold: 1e-8,
    };
    let ch_cfg = compute_charges_configured("CC(=O)O", &cfg)?;
    println!("Configured charges: {:?}", ch_cfg.charges.iter().map(|c| format!("{:.3}", c)).collect::<Vec<_>>());

    Ok(())
}

fn extract_elpos(r: &sci_form::ConformerResult) -> (Vec<u8>, Vec<[f64; 3]>) {
    let positions = r.coords.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect();
    (r.elements.clone(), positions)
}
