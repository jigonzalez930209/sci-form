use sci_form::{
    analyze_stereo, compute_charges, compute_dipole, compute_ecfp,
    compute_nonpolar_solvation, compute_pm3, compute_population, compute_sssr,
    compute_tanimoto, compute_xtb, embed, embed_batch, parse, ConformerConfig,
};

fn positions_from_flat(coords: &[f64]) -> Vec<[f64; 3]> {
    coords.chunks_exact(3).map(|chunk| [chunk[0], chunk[1], chunk[2]]).collect()
}

fn main() -> Result<(), String> {
    let ethanol = embed("CCO", 42);
    if let Some(error) = ethanol.error.clone() {
        return Err(error);
    }

    let positions = positions_from_flat(&ethanol.coords);
    println!("ethanol atoms: {}", ethanol.num_atoms);

    let batch = embed_batch(
        &["O", "CCO", "c1ccccc1"],
        &ConformerConfig {
            seed: 42,
            num_threads: 0,
        },
    );
    let successes = batch.iter().filter(|entry| entry.error.is_none()).count();
    println!("batch success: {successes}/{}", batch.len());

    let parsed = parse("CC(=O)O")?;
    println!(
        "acetic acid topology: {} atoms, {} bonds",
        parsed.graph.node_count(),
        parsed.graph.edge_count()
    );

    let charge_result = compute_charges("CC(=O)O")?;
    println!("Gasteiger total charge: {:.4}", charge_result.total_charge);

    let population = compute_population(&ethanol.elements, &positions)?;
    println!(
        "Mulliken total charge: {:.4}",
        population.total_charge_mulliken
    );

    let dipole = compute_dipole(&ethanol.elements, &positions)?;
    println!("dipole magnitude: {:.3} D", dipole.magnitude);

    let pm3 = compute_pm3(&ethanol.elements, &positions)?;
    println!("PM3 gap: {:.3} eV, converged: {}", pm3.gap, pm3.converged);

    let xtb = compute_xtb(&ethanol.elements, &positions)?;
    println!("xTB gap: {:.3} eV, converged: {}", xtb.gap, xtb.converged);

    let chiral = embed("C(F)(Cl)(Br)I", 42);
    if let Some(error) = chiral.error.clone() {
        return Err(error);
    }
    let stereo = analyze_stereo("C(F)(Cl)(Br)I", &chiral.coords)?;
    println!("stereocenters: {}", stereo.n_stereocenters);

    let solvation = compute_nonpolar_solvation(&ethanol.elements, &positions, None);
    println!(
        "non-polar solvation: {:.3} kcal/mol",
        solvation.energy_kcal_mol
    );

    let rings = compute_sssr("c1ccccc1")?;
    println!("benzene rings: {}", rings.rings.len());

    let benzene = compute_ecfp("c1ccccc1", 2, 2048)?;
    let toluene = compute_ecfp("Cc1ccccc1", 2, 2048)?;
    println!(
        "benzene/toluene tanimoto: {:.3}",
        compute_tanimoto(&benzene, &toluene)
    );

    Ok(())
}

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
