use sci_form::{
    analyze_stereo, compute_charges, compute_dipole, compute_ecfp, compute_nonpolar_solvation,
    compute_pm3, compute_population, compute_sssr, compute_tanimoto, compute_xtb, embed,
    embed_batch, parse, ConformerConfig,
};

#[cfg(feature = "alpha-dft")]
use sci_form::dft::ks_fock::{solve_ks_dft, DftConfig};
#[cfg(feature = "alpha-mlff")]
use sci_form::mlff::{compute_aevs, SymmetryFunctionParams};
#[cfg(feature = "alpha-reaxff")]
use sci_form::forcefield::reaxff::{compute_reaxff_gradient, ReaxffParams};
#[cfg(feature = "beta-cpm")]
use sci_form::beta::cpm::{compute_cpm_charges, CpmConfig};
#[cfg(feature = "beta-kpm")]
use sci_form::beta::kpm::{compute_kpm_dos, KpmConfig};
#[cfg(feature = "beta-randnla")]
use sci_form::beta::rand_nla::{solve_eht_randnla, RandNlaConfig};
#[cfg(any(feature = "beta-kpm", feature = "beta-randnla"))]
use sci_form::eht::solver::solve_eht;

fn positions_from_flat(coords: &[f64]) -> Vec<[f64; 3]> {
    coords
        .chunks_exact(3)
        .map(|chunk| [chunk[0], chunk[1], chunk[2]])
        .collect()
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
    println!("Mulliken total charge: {:.4}", population.total_charge_mulliken);

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
    println!("non-polar solvation: {:.3} kcal/mol", solvation.energy_kcal_mol);

    let rings = compute_sssr("c1ccccc1")?;
    println!("benzene rings: {}", rings.rings.len());

    let benzene = compute_ecfp("c1ccccc1", 2, 2048)?;
    let toluene = compute_ecfp("Cc1ccccc1", 2, 2048)?;
    println!("benzene/toluene tanimoto: {:.3}", compute_tanimoto(&benzene, &toluene));

    #[cfg(feature = "alpha-dft")]
    {
        let dft = solve_ks_dft(&ethanol.elements, &positions, &DftConfig::default())?;
        println!("DFT gap: {:.3} eV, converged: {}", dft.gap, dft.converged);
    }

    #[cfg(feature = "alpha-reaxff")]
    {
        let params = ReaxffParams::default_chon();
        let (energy, gradient) = compute_reaxff_gradient(&ethanol.coords, &ethanol.elements, &params)?;
        println!("ReaxFF energy: {:.3} gradient length: {}", energy, gradient.len());
    }

    #[cfg(feature = "alpha-mlff")]
    {
        let aevs = compute_aevs(&ethanol.elements, &positions, &SymmetryFunctionParams::default());
        println!("AEV atoms: {}", aevs.len());
    }

    #[cfg(feature = "beta-kpm")]
    {
        let eht = solve_eht(&ethanol.elements, &positions, None)?;
        let config = KpmConfig::default();
        let kpm = compute_kpm_dos(&eht.hamiltonian, &config, -30.0, 5.0, 256);
        println!("KPM points: {} order: {}", kpm.energies.len(), kpm.order);
    }

    #[cfg(feature = "beta-randnla")]
    {
        let eht = solve_eht(&ethanol.elements, &positions, None)?;
        let config = RandNlaConfig::default();
        let (orbital_energies, _, info) = solve_eht_randnla(&eht.hamiltonian, &eht.overlap, &config);
        println!(
            "RandNLA orbitals: {} residual: {:.3e}",
            orbital_energies.len(),
            info.residual_error
        );
    }

    #[cfg(feature = "beta-cpm")]
    {
        let cpm = compute_cpm_charges(&ethanol.elements, &positions, &CpmConfig::default());
        println!("CPM total charge: {:.4}", cpm.total_charge);
    }

    Ok(())
}
