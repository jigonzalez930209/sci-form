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
//! Diagnostic test: compare 3D FF energy at our coordinates vs RDKit reference.
//!
//! Run: cargo test --release --test debug_ff_energy -- --nocapture

use serde::Deserialize;
use std::fs;

#[derive(Deserialize)]
struct RefAtom {
    element: u8,
    x: f32,
    y: f32,
    z: f32,
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
    signs: Vec<i32>,
}

#[derive(Deserialize)]
struct RefMolecule {
    smiles: String,
    atoms: Vec<RefAtom>,
    bonds: Vec<RefBond>,
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
            position: nalgebra::Vector3::new(0.0, 0.0, 0.0),
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
                v[k] = t.v[k] as f64;
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

fn pairwise_rmsd(coords: &nalgebra::DMatrix<f32>, ref_atoms: &[RefAtom]) -> f32 {
    let n = ref_atoms.len();
    let mut sq_sum = 0.0f64;
    let mut npairs = 0u64;
    for a in 0..n {
        for b in (a + 1)..n {
            let dr = ((ref_atoms[a].x - ref_atoms[b].x).powi(2)
                + (ref_atoms[a].y - ref_atoms[b].y).powi(2)
                + (ref_atoms[a].z - ref_atoms[b].z).powi(2))
            .sqrt() as f64;
            let du = ((coords[(a, 0)] - coords[(b, 0)]).powi(2)
                + (coords[(a, 1)] - coords[(b, 1)]).powi(2)
                + (coords[(a, 2)] - coords[(b, 2)]).powi(2))
            .sqrt() as f64;
            sq_sum += (dr - du).powi(2);
            npairs += 1;
        }
    }
    if npairs > 0 {
        (sq_sum / npairs as f64).sqrt() as f32
    } else {
        0.0
    }
}

#[test]
fn test_ff_energy_comparison() {
    use sci_form::distgeom::bounds::{calculate_bounds_matrix_opts, triangle_smooth_tol};
    use sci_form::forcefield::etkdg_3d::{build_etkdg_3d_ff_with_torsions, etkdg_3d_energy_f64};

    let ref_data = fs::read_to_string("tests/fixtures/gdb20_reference.json")
        .expect("Run scripts/generate_gdb20_reference.py first");
    let ref_mols: Vec<RefMolecule> =
        serde_json::from_str(&ref_data).expect("Invalid gdb20_reference.json");

    // Target: first 1000 molecules, find worst 10 failures
    let limit = std::env::var("GDB20_LIMIT")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(1000usize);
    let ref_mols = &ref_mols[..limit.min(ref_mols.len())];

    println!("\n=== FF ENERGY DIAGNOSTIC ===");

    // Find failing molecules
    let mut results: Vec<(usize, f32, String)> = Vec::new();
    for (idx, ref_mol) in ref_mols.iter().enumerate() {
        let mol = build_mol_from_ref(ref_mol);
        let csd_torsions = build_csd_torsions(&ref_mol.torsions);
        let result =
            sci_form::conformer::generate_3d_conformer_with_torsions(&mol, 42, &csd_torsions);
        if let Ok(coords) = result {
            let rmsd = pairwise_rmsd(&coords, &ref_mol.atoms);
            if rmsd > 0.5 {
                results.push((idx, rmsd, ref_mol.smiles.clone()));
            }
        }
    }
    results.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());

    println!("Found {} failing molecules", results.len());

    // Analyze top 10 failures
    for (rank, (idx, rmsd, smi)) in results.iter().take(10).enumerate() {
        let ref_mol = &ref_mols[*idx];
        let mol = build_mol_from_ref(ref_mol);
        let n = mol.graph.node_count();
        let csd_torsions = build_csd_torsions(&ref_mol.torsions);

        // Generate our conformer
        let our_coords =
            sci_form::conformer::generate_3d_conformer_with_torsions(&mol, 42, &csd_torsions)
                .unwrap();

        // Build reference coords (RDKit's)
        let ref_coords_f32 = nalgebra::DMatrix::from_fn(n, 3, |i, j| match j {
            0 => ref_mol.atoms[i].x,
            1 => ref_mol.atoms[i].y,
            2 => ref_mol.atoms[i].z,
            _ => unreachable!(),
        });

        // Build bounds matrix (same as in pipeline)
        let bounds = {
            let raw = calculate_bounds_matrix_opts(&mol, true);
            let mut b = raw;
            if triangle_smooth_tol(&mut b, 0.0) {
                b
            } else {
                let raw2 = calculate_bounds_matrix_opts(&mol, false);
                let mut b2 = raw2.clone();
                if triangle_smooth_tol(&mut b2, 0.0) {
                    b2
                } else {
                    let mut b3 = raw2;
                    triangle_smooth_tol(&mut b3, 0.05);
                    b3
                }
            }
        };

        // Build 3D FF using OUR coordinates (this is what the pipeline does)
        let ff_ours = build_etkdg_3d_ff_with_torsions(
            &mol,
            &our_coords.map(|v| v as f64),
            &bounds,
            &csd_torsions,
        );

        // Build 3D FF using REFERENCE coordinates (this is what RDKit would do)
        let ff_ref = build_etkdg_3d_ff_with_torsions(
            &mol,
            &ref_coords_f32.map(|v| v as f64),
            &bounds,
            &csd_torsions,
        );

        // Evaluate our FF at our coordinates
        let our_coords_f64: Vec<f64> = (0..n)
            .flat_map(|i| {
                vec![
                    our_coords[(i, 0)] as f64,
                    our_coords[(i, 1)] as f64,
                    our_coords[(i, 2)] as f64,
                ]
            })
            .collect();
        let e_our_at_our = etkdg_3d_energy_f64(&our_coords_f64, n, &mol, &ff_ours);

        // Evaluate our FF at reference coordinates
        let ref_coords_f64: Vec<f64> = (0..n)
            .flat_map(|i| {
                vec![
                    ref_mol.atoms[i].x as f64,
                    ref_mol.atoms[i].y as f64,
                    ref_mol.atoms[i].z as f64,
                ]
            })
            .collect();
        let e_our_at_ref = etkdg_3d_energy_f64(&ref_coords_f64, n, &mol, &ff_ours);

        // Evaluate reference FF at reference coordinates
        let e_ref_at_ref = etkdg_3d_energy_f64(&ref_coords_f64, n, &mol, &ff_ref);

        // Evaluate reference FF at our coordinates
        let e_ref_at_our = etkdg_3d_energy_f64(&our_coords_f64, n, &mol, &ff_ref);

        // Count constraints
        let n_torsions = ff_ours.torsion_contribs.len();
        let n_inversions = ff_ours.inversion_contribs.len();
        let n_dist = ff_ours.dist_12.len() + ff_ours.dist_13.len() + ff_ours.dist_long.len();
        let n_angles = ff_ours.angle_constraints.len();

        // Count constraint types for our FF
        let n_bond_constraints = mol.graph.edge_count();
        let n_long_range = n_dist - n_bond_constraints; // rough estimate

        println!("\n--- #{} SMILES={} RMSD={:.4} ---", rank + 1, smi, rmsd);
        println!("  Atoms: {}, Bonds: {}", n, mol.graph.edge_count());
        println!(
            "  CSD torsions: {}, Ring/chain torsions: {}",
            csd_torsions.len(),
            n_torsions - csd_torsions.len()
        );
        println!("  Inversion contribs: {}", n_inversions);
        println!(
            "  Dist constraints: {} (approx {} bond + {} other)",
            n_dist,
            n_bond_constraints,
            n_dist - n_bond_constraints
        );
        println!("  Angle constraints: {}", n_angles);
        println!("  FF(ours) @ our coords:   {:.6}", e_our_at_our);
        println!("  FF(ours) @ ref coords:   {:.6}", e_our_at_ref);
        println!("  FF(ref)  @ ref coords:   {:.6}", e_ref_at_ref);
        println!("  FF(ref)  @ our coords:   {:.6}", e_ref_at_our);

        // The key comparison:
        // If e_our_at_our < e_our_at_ref: our optimizer found a lower-energy state
        //   (for our FF). The FF landscape differs from RDKit's.
        // If e_our_at_our > e_our_at_ref: our optimizer missed a better minimum in
        //   our own FF landscape. Optimizer bug.
        if e_our_at_our < e_our_at_ref {
            println!("  → Our optimum is LOWER than ref in our FF (FF landscape differs)");
        } else {
            println!("  → Our optimum is HIGHER than ref in our FF (optimizer issue)");
        }
    }
}
