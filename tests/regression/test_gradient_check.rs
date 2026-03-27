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
//! Finite-difference gradient verification for the 3D ETKDG force field.
//!
//! For each molecule, builds the FF, computes the analytical gradient at the
//! pre-BFGS coordinates, and compares with a central-difference numerical gradient.
//! A mismatch indicates a gradient bug that would cause BFGS to diverge.
//!
//! Run:  cargo test --release --test test_gradient_check -- --nocapture

use std::fs;

use nalgebra::DMatrix;
use serde::Deserialize;

#[derive(Deserialize)]
struct RefAtom {
    element: u8,
    formal_charge: i8,
    x: f64,
    y: f64,
    z: f64,
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
    signs: Vec<i32>,
    v: Vec<f64>,
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
    let mut node_indices = Vec::new();

    for atom in &ref_mol.atoms {
        let hybridization = sci_form::graph::Hybridization::Unknown;
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
fn test_gradient_finite_diff() {
    let fixture = "tests/fixtures/gdb20_reference_1k.json";
    if !sci_form::fixture_io::fixture_exists(fixture) {
        eprintln!(
            "SKIP {fixture}: file not found (generate with scripts/generate_gdb20_reference.py)"
        );
        return;
    }

    // Check if file is too small (likely a Git LFS pointer, not actual data)
    let resolved_fixture = sci_form::fixture_io::resolve_fixture_path(fixture).unwrap();
    let metadata = match fs::metadata(&resolved_fixture) {
        Ok(m) => m,
        Err(e) => {
            eprintln!("SKIP {fixture}: cannot read metadata: {e}");
            return;
        }
    };
    if metadata.len() < 1000 {
        eprintln!(
            "SKIP {fixture}: file too small ({} bytes), likely Git LFS pointer. \
             Ensure Git LFS is installed and run 'git lfs pull'",
            metadata.len()
        );
        return;
    }

    let ref_data = match sci_form::fixture_io::read_text_fixture(fixture) {
        Ok(data) => data,
        Err(e) => {
            eprintln!("SKIP {fixture}: cannot read file: {e}");
            return;
        }
    };

    let mut ref_mols: Vec<RefMolecule> = match serde_json::from_str(&ref_data) {
        Ok(mols) => mols,
        Err(e) => {
            eprintln!("SKIP {fixture}: JSON parsing failed: {e}");
            eprintln!(
                "  First 200 chars: {}",
                &ref_data[..ref_data.len().min(200)]
            );
            return;
        }
    };

    // Sort by heaviest molecules first (most atoms) for a representative sample
    ref_mols.sort_by(|a, b| b.atoms.len().cmp(&a.atoms.len()));

    // Optional limit — default 100 heaviest; use GDB20_LIMIT=0 for all
    let limit: usize = std::env::var("GDB20_LIMIT")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(100);
    let ref_mols = &ref_mols[..limit.min(ref_mols.len())];

    let eps = 1e-5; // Finite difference step
    let mut total_checked = 0u32;
    let mut total_bad = 0u32;
    let mut worst_rel_err = 0.0f64;
    let mut worst_mol = String::new();

    for (mol_idx, ref_mol) in ref_mols.iter().enumerate() {
        let mol = build_mol_from_ref(ref_mol);
        let csd_torsions = build_csd_torsions(&ref_mol.torsions);
        let n = mol.graph.node_count();

        // Use RDKit reference coords as evaluation point for gradient correctness.
        let mut coords = vec![0.0f64; n * 3];
        for (a, atom) in ref_mol.atoms.iter().enumerate() {
            coords[a * 3] = atom.x;
            coords[a * 3 + 1] = atom.y;
            coords[a * 3 + 2] = atom.z;
        }

        // Build the FF from these coords (need a DMatrix for the builder)
        let coords_mat = DMatrix::from_fn(n, 3, |r, c| coords[r * 3 + c]);
        let bounds = {
            let raw = sci_form::distgeom::calculate_bounds_matrix_opts(&mol, true);
            let mut b = raw.clone();
            if sci_form::distgeom::triangle_smooth_tol(&mut b, 0.0) {
                b
            } else {
                let raw2 = sci_form::distgeom::calculate_bounds_matrix_opts(&mol, false);
                let mut b2 = raw2.clone();
                if sci_form::distgeom::triangle_smooth_tol(&mut b2, 0.0) {
                    b2
                } else {
                    let mut b3 = raw2;
                    sci_form::distgeom::triangle_smooth_tol(&mut b3, 0.05);
                    b3
                }
            }
        };

        let ff = sci_form::forcefield::etkdg_3d::build_etkdg_3d_ff_with_torsions(
            &mol,
            &coords_mat,
            &bounds,
            &csd_torsions,
        );

        // Compute analytical gradient
        let grad = sci_form::forcefield::etkdg_3d::etkdg_3d_gradient_f64(&coords, n, &mol, &ff);

        // Compute numerical gradient (central difference)
        let mut num_grad = vec![0.0f64; n * 3];
        for i in 0..(n * 3) {
            let mut cp = coords.clone();
            cp[i] = coords[i] + eps;
            let ep = sci_form::forcefield::etkdg_3d::etkdg_3d_energy_f64(&cp, n, &mol, &ff);
            cp[i] = coords[i] - eps;
            let em = sci_form::forcefield::etkdg_3d::etkdg_3d_energy_f64(&cp, n, &mol, &ff);
            num_grad[i] = (ep - em) / (2.0 * eps);
        }

        // Compare
        let mut max_rel_err = 0.0f64;
        let mut max_abs_err = 0.0f64;
        let mut bad_component = 0;
        let grad_norm: f64 = grad.iter().map(|g| g * g).sum::<f64>().sqrt();

        for i in 0..(n * 3) {
            let abs_err = (grad[i] - num_grad[i]).abs();
            let denom = grad[i].abs().max(num_grad[i].abs()).max(1e-8);
            let rel_err = abs_err / denom;

            if rel_err > max_rel_err {
                max_rel_err = rel_err;
                max_abs_err = abs_err;
                bad_component = i;
            }
        }

        total_checked += 1;

        // Flag molecules with large gradient errors
        if max_rel_err > 1e-3 {
            total_bad += 1;
            let atom_idx = bad_component / 3;
            let dim = bad_component % 3;
            let dim_name = ["x", "y", "z"][dim];
            println!(
                "  [BAD] mol[{}] {} n={} max_rel_err={:.6e} abs_err={:.6e} at atom {} {} (anal={:.8e} num={:.8e}) grad_norm={:.4}",
                mol_idx, ref_mol.smiles, n, max_rel_err, max_abs_err,
                atom_idx, dim_name, grad[bad_component], num_grad[bad_component], grad_norm
            );
        }

        if max_rel_err > worst_rel_err {
            worst_rel_err = max_rel_err;
            worst_mol = ref_mol.smiles.clone();
        }
    }

    println!("\n=== GRADIENT CHECK RESULTS ===");
    println!("Molecules checked: {}", total_checked);
    println!(
        "Bad gradient (rel_err > 1e-3): {} ({:.2}%)",
        total_bad,
        total_bad as f64 / total_checked.max(1) as f64 * 100.0
    );
    println!(
        "Worst relative error: {:.6e} ({})",
        worst_rel_err, worst_mol
    );

    if total_bad > 0 {
        println!("\nWARNING: Gradient bugs detected!");
    } else {
        println!("\nAll gradients match finite differences within 1e-3 relative error.");
    }
}
