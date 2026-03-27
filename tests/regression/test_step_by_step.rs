#![allow(
    unused_imports,
    unused_variables,
    dead_code,
    clippy::unnecessary_cast,
    clippy::needless_range_loop,
    clippy::manual_repeat_n,
    clippy::manual_str_repeat,
    clippy::manual_is_multiple_of,
    clippy::redundant_field_names,
    clippy::useless_vec,
    clippy::single_range_in_vec_init
)]
use serde::Deserialize;
/// Trace the full pipeline step by step for a specific molecule.
/// Compares our intermediate values with expected behavior.
use std::fs;

#[derive(Deserialize)]
struct RefAtom {
    element: u8,
    x: f64,
    y: f64,
    z: f64,
    formal_charge: i8,
    hybridization: String,
}
#[derive(Deserialize)]
struct RefBond {
    start: usize,
    end: usize,
    order: String,
}
#[derive(Deserialize)]
struct RefTorsion {
    atoms: Vec<usize>,
    v: Vec<f64>,
    signs: Vec<f64>,
}
#[derive(Deserialize)]
struct RefMolecule {
    smiles: String,
    atoms: Vec<RefAtom>,
    bonds: Vec<RefBond>,
    #[serde(default)]
    torsions: Vec<RefTorsion>,
}

fn build_mol_from_ref(ref_mol: &RefMolecule) -> sci_form::graph::Molecule {
    let mut mol = sci_form::graph::Molecule::new(&ref_mol.smiles);
    let mut node_indices = Vec::with_capacity(ref_mol.atoms.len());
    for atom in &ref_mol.atoms {
        let hybridization = match atom.hybridization.as_str() {
            "SP" => sci_form::graph::Hybridization::SP,
            "SP2" => sci_form::graph::Hybridization::SP2,
            "SP3" => sci_form::graph::Hybridization::SP3,
            "SP3D" => sci_form::graph::Hybridization::SP3D,
            "SP3D2" => sci_form::graph::Hybridization::SP3D2,
            _ => sci_form::graph::Hybridization::Unknown,
        };
        let new_atom = sci_form::graph::Atom {
            element: atom.element,
            position: nalgebra::Vector3::zeros(),
            charge: 0.0,
            formal_charge: atom.formal_charge,
            hybridization,
            chiral_tag: sci_form::graph::ChiralType::Unspecified,
            explicit_h: if atom.element == 1 || atom.element == 0 {
                1
            } else {
                0
            },
        };
        node_indices.push(mol.add_atom(new_atom));
    }
    for bond in &ref_mol.bonds {
        let order = match bond.order.as_str() {
            "DOUBLE" => sci_form::graph::BondOrder::Double,
            "TRIPLE" => sci_form::graph::BondOrder::Triple,
            "AROMATIC" => sci_form::graph::BondOrder::Aromatic,
            _ => sci_form::graph::BondOrder::Single,
        };
        mol.add_bond(
            node_indices[bond.start],
            node_indices[bond.end],
            sci_form::graph::Bond {
                order,
                stereo: sci_form::graph::BondStereo::None,
            },
        );
    }
    mol
}

fn build_csd_torsions(
    ref_torsions: &[RefTorsion],
) -> Vec<sci_form::forcefield::etkdg_3d::M6TorsionContrib> {
    ref_torsions
        .iter()
        .filter_map(|t| {
            if t.atoms.len() < 4 || t.v.len() < 6 || t.signs.len() < 6 {
                return None;
            }
            let mut signs = [0.0f64; 6];
            let mut v = [0.0f64; 6];
            for k in 0..6 {
                signs[k] = t.signs[k] as f64;
                v[k] = t.v[k];
            }
            Some(sci_form::forcefield::etkdg_3d::M6TorsionContrib {
                i: t.atoms[0],
                j: t.atoms[1],
                k: t.atoms[2],
                l: t.atoms[3],
                signs,
                v,
            })
        })
        .collect()
}

#[test]
fn test_step_by_step_trace() {
    let ref_data = sci_form::fixture_io::read_text_fixture("tests/fixtures/gdb20_reference_1k.json")
        .unwrap();
    let ref_mols: Vec<RefMolecule> = serde_json::from_str(&ref_data).unwrap();

    // Target: mol[10] = N#Cc1ccc(Br)c(CC2CCCCC2)c1C=O (aromatic, RMSD=0.9280)
    let target_smi = "N#Cc1ccc(Br)c(CC2CCCCC2)c1C=O";
    let ref_mol = ref_mols
        .iter()
        .find(|m| m.smiles == target_smi)
        .unwrap_or_else(|| panic!("SMILES '{}' not found", target_smi));

    let mol = build_mol_from_ref(ref_mol);
    let n = mol.graph.node_count();
    eprintln!(
        "=== STEP-BY-STEP TRACE for '{}' (n={}) ===\n",
        target_smi, n
    );

    // Step 1: Bounds matrix
    let mut bounds = sci_form::distgeom::calculate_bounds_matrix_opts(&mol, true);
    let smooth_ok = sci_form::distgeom::triangle_smooth_tol(&mut bounds, 0.0);
    eprintln!("Step 1: Bounds matrix {}×{}, smoothing={}", n, n, smooth_ok);

    // Step 2: Check chiral sets
    let chiral_sets = sci_form::distgeom::identify_chiral_sets(&mol);
    let tet_centers = sci_form::distgeom::identify_tetrahedral_centers(&mol);
    let use_4d = !chiral_sets.is_empty();
    let embed_dim = if use_4d { 4 } else { 3 };
    eprintln!(
        "Step 2: chiral_sets={}, tet_centers={}, use_4d={}, embed_dim={}",
        chiral_sets.len(),
        tet_centers.len(),
        use_4d,
        embed_dim
    );

    // Step 3: Retry loop — trace each attempt
    let max_iterations = 10 * n;
    let mut rng = sci_form::distgeom::MinstdRand::new(42);

    for iter in 0..max_iterations {
        let rng_state_before = rng.get_state();
        let dists = sci_form::distgeom::pick_rdkit_distances(&mut rng, &bounds);
        let rng_state_after_dists = rng.get_state();

        let coords_opt =
            sci_form::distgeom::compute_initial_coords_rdkit(&mut rng, &dists, embed_dim);
        let rng_state_after_embed = rng.get_state();

        let rng_calls_dists = {
            // Count: n*(n-1)/2 calls for distance picking
            n * (n - 1) / 2
        };

        if coords_opt.is_none() {
            if iter < 5 || iter == max_iterations - 1 {
                eprintln!("  attempt {} → embedding FAILED (rng_state_before={}, after_dists={}, after_embed={})",
                    iter, rng_state_before, rng_state_after_dists, rng_state_after_embed);
            }
            continue;
        }

        let mut coords = coords_opt.unwrap();

        // Compute initial pairwise RMSD (before any optimization)
        let initial_rmsd = {
            let mut sq_sum = 0.0f64;
            let mut npairs = 0u64;
            for a in 0..n {
                for b in (a + 1)..n {
                    let dr = ((ref_mol.atoms[a].x - ref_mol.atoms[b].x).powi(2)
                        + (ref_mol.atoms[a].y - ref_mol.atoms[b].y).powi(2)
                        + (ref_mol.atoms[a].z - ref_mol.atoms[b].z).powi(2))
                    .sqrt();
                    let du = ((coords[(a, 0)] - coords[(b, 0)]).powi(2)
                        + (coords[(a, 1)] - coords[(b, 1)]).powi(2)
                        + (coords[(a, 2)] - coords[(b, 2)]).powi(2))
                    .sqrt();
                    sq_sum += (dr - du).powi(2);
                    npairs += 1;
                }
            }
            (sq_sum / npairs as f64).sqrt()
        };

        eprintln!(
            "  attempt {} → embedding OK (rng_before={}, initial_RMSD={:.4})",
            iter, rng_state_before, initial_rmsd
        );

        // Bounds FF #1
        let basin_thresh = 5.0;
        let force_tol = 1e-3;
        let error_tol = 1e-5;

        let initial_energy = sci_form::conformer::compute_total_bounds_energy_f64(
            &coords,
            &bounds,
            &chiral_sets,
            basin_thresh,
            0.1,
            1.0,
        );
        eprintln!("    bounds_ff_1: initial_energy={:.6}", initial_energy);

        if initial_energy > error_tol {
            let mut need_more = 1;
            let mut bfgs_calls = 0;
            while need_more != 0 {
                need_more = sci_form::forcefield::bounds_ff::minimize_bfgs_rdkit(
                    &mut coords,
                    &bounds,
                    &chiral_sets,
                    400,
                    force_tol as f64,
                    basin_thresh,
                    0.1,
                    1.0,
                );
                bfgs_calls += 1;
            }
            let post_energy = sci_form::conformer::compute_total_bounds_energy_f64(
                &coords,
                &bounds,
                &chiral_sets,
                basin_thresh,
                0.1,
                1.0,
            );
            eprintln!(
                "    bounds_ff_1: {} BFGS calls, final_energy={:.6}",
                bfgs_calls, post_energy
            );
        }

        // Energy check
        let energy = sci_form::conformer::compute_total_bounds_energy_f64(
            &coords,
            &bounds,
            &chiral_sets,
            basin_thresh,
            0.1,
            1.0,
        );
        let e_per_atom = energy / n as f64;
        if e_per_atom >= 0.05 {
            eprintln!("    REJECT: energy/atom={:.6} >= 0.05", e_per_atom);
            continue;
        }

        // Tet check
        if !sci_form::distgeom::check_tetrahedral_centers(&coords, &tet_centers) {
            eprintln!("    REJECT: tetrahedral check");
            continue;
        }

        // Chiral check
        if use_4d && !sci_form::distgeom::check_chiral_centers(&coords, &chiral_sets) {
            eprintln!("    REJECT: chiral check");
            continue;
        }

        // Drop to 3D
        let coords3d = coords.columns(0, 3).into_owned();

        // Post-bounds-FF RMSD
        let post_bounds_rmsd = {
            let mut sq_sum = 0.0f64;
            let mut npairs = 0u64;
            for a in 0..n {
                for b in (a + 1)..n {
                    let dr = ((ref_mol.atoms[a].x - ref_mol.atoms[b].x).powi(2)
                        + (ref_mol.atoms[a].y - ref_mol.atoms[b].y).powi(2)
                        + (ref_mol.atoms[a].z - ref_mol.atoms[b].z).powi(2))
                    .sqrt();
                    let du = ((coords3d[(a, 0)] - coords3d[(b, 0)]).powi(2)
                        + (coords3d[(a, 1)] - coords3d[(b, 1)]).powi(2)
                        + (coords3d[(a, 2)] - coords3d[(b, 2)]).powi(2))
                    .sqrt();
                    sq_sum += (dr - du).powi(2);
                    npairs += 1;
                }
            }
            (sq_sum / npairs as f64).sqrt()
        };
        eprintln!("    post_bounds_RMSD={:.4}", post_bounds_rmsd);

        // 3D FF
        let csd_torsions = build_csd_torsions(&ref_mol.torsions);
        let ff = sci_form::forcefield::etkdg_3d::build_etkdg_3d_ff_with_torsions(
            &mol,
            &coords3d,
            &bounds,
            &csd_torsions,
        );

        let flat_f64: Vec<f64> = {
            let mut flat = vec![0.0f64; n * 3];
            for a in 0..n {
                flat[a * 3] = coords3d[(a, 0)];
                flat[a * 3 + 1] = coords3d[(a, 1)];
                flat[a * 3 + 2] = coords3d[(a, 2)];
            }
            flat
        };
        let e3d = sci_form::forcefield::etkdg_3d::etkdg_3d_energy_f64(&flat_f64, n, &mol, &ff);
        eprintln!("    3d_ff: initial_energy={:.6}", e3d);

        let refined = if e3d > error_tol {
            sci_form::forcefield::etkdg_3d::minimize_etkdg_3d_bfgs(&mol, &coords3d, &ff, 300, 1e-3)
        } else {
            coords3d.clone()
        };

        // Final RMSD
        let final_rmsd = {
            let mut sq_sum = 0.0f64;
            let mut npairs = 0u64;
            for a in 0..n {
                for b in (a + 1)..n {
                    let dr = ((ref_mol.atoms[a].x - ref_mol.atoms[b].x).powi(2)
                        + (ref_mol.atoms[a].y - ref_mol.atoms[b].y).powi(2)
                        + (ref_mol.atoms[a].z - ref_mol.atoms[b].z).powi(2))
                    .sqrt();
                    let du = ((refined[(a, 0)] - refined[(b, 0)]).powi(2)
                        + (refined[(a, 1)] - refined[(b, 1)]).powi(2)
                        + (refined[(a, 2)] - refined[(b, 2)]).powi(2))
                    .sqrt();
                    sq_sum += (dr - du).powi(2);
                    npairs += 1;
                }
            }
            (sq_sum / npairs as f64).sqrt()
        };

        eprintln!(
            "    FINAL RMSD={:.4} (initial={:.4}, post_bounds={:.4})",
            final_rmsd, initial_rmsd, post_bounds_rmsd
        );

        // Planarity check
        let n_improper = ff.inversion_contribs.len() / 3;
        let flat_refined: Vec<f64> = {
            let nr = refined.nrows();
            let mut flat = vec![0.0f64; nr * 3];
            for a in 0..nr {
                flat[a * 3] = refined[(a, 0)];
                flat[a * 3 + 1] = refined[(a, 1)];
                flat[a * 3 + 2] = refined[(a, 2)];
            }
            flat
        };
        let plan_e =
            sci_form::forcefield::etkdg_3d::planarity_check_energy_f64(&flat_refined, n, &ff);
        if plan_e > n_improper as f64 * 0.7 {
            eprintln!(
                "    REJECT: planarity (energy={:.4} > {:.4})",
                plan_e,
                n_improper as f64 * 0.7
            );
            continue;
        }

        // Double bond check
        if !sci_form::distgeom::check_double_bond_geometry(&mol, &refined) {
            eprintln!("    REJECT: double bond geometry");
            continue;
        }

        eprintln!("    SUCCESS on attempt {}", iter);
        eprintln!("\n=== DONE ===");
        return;
    }

    eprintln!(
        "FAILED: no valid conformer after {} attempts",
        max_iterations
    );
}
