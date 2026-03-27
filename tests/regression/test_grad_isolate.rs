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
//! Isolate which force field term has the gradient bug by testing each term type separately.
//! Run:  cargo test --release --test test_grad_isolate -- --nocapture

use nalgebra::DMatrix;
use serde::Deserialize;
use std::fs;

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
        let new_atom = sci_form::graph::Atom {
            element: atom.element,
            position: nalgebra::Vector3::new(0.0, 0.0, 0.0),
            charge: 0.0,
            formal_charge: atom.formal_charge,
            hybridization: sci_form::graph::Hybridization::Unknown,
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

fn build_csd_torsions(t: &[RefTorsion]) -> Vec<sci_form::forcefield::etkdg_3d::M6TorsionContrib> {
    t.iter()
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
fn test_grad_isolate() {
    let ref_data =
        sci_form::fixture_io::read_text_fixture("tests/fixtures/gdb20_reference_1k.json").unwrap();
    let ref_mols: Vec<RefMolecule> = serde_json::from_str(&ref_data).unwrap();

    // Test one of the worst molecules: mol[47]
    let test_indices = [4, 7, 47];
    let eps = 1e-5;

    for &mol_idx in &test_indices {
        if mol_idx >= ref_mols.len() {
            continue;
        }
        let ref_mol = &ref_mols[mol_idx];
        let mol = build_mol_from_ref(ref_mol);
        let csd_torsions = build_csd_torsions(&ref_mol.torsions);
        let n = mol.graph.node_count();

        let mut coords = vec![0.0f64; n * 3];
        for (a, atom) in ref_mol.atoms.iter().enumerate() {
            coords[a * 3] = atom.x;
            coords[a * 3 + 1] = atom.y;
            coords[a * 3 + 2] = atom.z;
        }

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

        let ff_full = sci_form::forcefield::etkdg_3d::build_etkdg_3d_ff_with_torsions(
            &mol,
            &coords_mat,
            &bounds,
            &csd_torsions,
        );

        println!("\n=== mol[{}] {} (n={}) ===", mol_idx, ref_mol.smiles, n);

        // Find worst component in full gradient
        let anal_full =
            sci_form::forcefield::etkdg_3d::etkdg_3d_gradient_f64(&coords, n, &mol, &ff_full);
        let num_full = numerical_gradient(
            &|c| sci_form::forcefield::etkdg_3d::etkdg_3d_energy_f64(c, n, &mol, &ff_full),
            &coords,
            eps,
        );
        let mut worst_i = 0;
        let mut worst_rel = 0.0f64;
        for i in 0..(n * 3) {
            let denom = anal_full[i].abs().max(num_full[i].abs()).max(1e-8);
            let rel = (anal_full[i] - num_full[i]).abs() / denom;
            if rel > worst_rel {
                worst_rel = rel;
                worst_i = i;
            }
        }
        let atom = worst_i / 3;
        let dim = ["x", "y", "z"][worst_i % 3];
        println!(
            "  FULL gradient worst: comp={} (atom {} {}) anal={:.8e} num={:.8e} rel={:.4e}",
            worst_i, atom, dim, anal_full[worst_i], num_full[worst_i], worst_rel
        );

        // Now test each term type in isolation by creating FFs with only one term
        let term_names = [
            "torsions",
            "inversions",
            "dist_12",
            "angles",
            "dist_13",
            "dist_long",
        ];

        for (ti, term_name) in term_names.iter().enumerate() {
            // Create a copy of FF with only one term type
            let ff_iso = sci_form::forcefield::etkdg_3d::Etkdg3DFF {
                torsion_contribs: if ti == 0 {
                    ff_full.torsion_contribs.clone()
                } else {
                    vec![]
                },
                inversion_contribs: if ti == 1 {
                    ff_full.inversion_contribs.clone()
                } else {
                    vec![]
                },
                dist_12: if ti == 2 {
                    ff_full.dist_12.clone()
                } else {
                    vec![]
                },
                angle_constraints: if ti == 3 {
                    ff_full.angle_constraints.clone()
                } else {
                    vec![]
                },
                dist_13: if ti == 4 {
                    ff_full.dist_13.clone()
                } else {
                    vec![]
                },
                dist_long: if ti == 5 {
                    ff_full.dist_long.clone()
                } else {
                    vec![]
                },
                oop_k: ff_full.oop_k,
                torsion_k_omega: ff_full.torsion_k_omega,
                bounds_force_scaling: ff_full.bounds_force_scaling,
                use_m6_torsions: ff_full.use_m6_torsions,
            };

            let anal =
                sci_form::forcefield::etkdg_3d::etkdg_3d_gradient_f64(&coords, n, &mol, &ff_iso);
            let num = numerical_gradient(
                &|c| sci_form::forcefield::etkdg_3d::etkdg_3d_energy_f64(c, n, &mol, &ff_iso),
                &coords,
                eps,
            );

            // Check worst component across ALL components (not just the global worst)
            let mut term_worst_i = 0;
            let mut term_worst_rel = 0.0f64;
            for i in 0..(n * 3) {
                let denom = anal[i].abs().max(num[i].abs()).max(1e-8);
                let rel = (anal[i] - num[i]).abs() / denom;
                if rel > term_worst_rel {
                    term_worst_rel = rel;
                    term_worst_i = i;
                }
            }

            let ta = term_worst_i / 3;
            let td = ["x", "y", "z"][term_worst_i % 3];

            // Also check the globally-worst component
            let denom_gw = anal[worst_i].abs().max(num[worst_i].abs()).max(1e-8);
            let rel_gw = (anal[worst_i] - num[worst_i]).abs() / denom_gw;

            if term_worst_rel > 1e-3 {
                println!("  [BAD]  {:12}: worst_rel={:.4e} at atom {} {} (A={:.4e} N={:.4e})  |  at global worst: A={:.4e} N={:.4e} rel={:.4e}",
                    term_name, term_worst_rel, ta, td,
                    anal[term_worst_i], num[term_worst_i],
                    anal[worst_i], num[worst_i], rel_gw);
            } else {
                println!(
                    "  [OK]   {:12}: worst_rel={:.4e}  |  at global worst: A={:.4e} N={:.4e}",
                    term_name, term_worst_rel, anal[worst_i], num[worst_i]
                );
            }
        }
    }
}

fn numerical_gradient(energy_fn: &dyn Fn(&[f64]) -> f64, coords: &[f64], eps: f64) -> Vec<f64> {
    let dim = coords.len();
    let mut grad = vec![0.0f64; dim];
    for i in 0..dim {
        let mut cp = coords.to_vec();
        cp[i] = coords[i] + eps;
        let ep = energy_fn(&cp);
        cp[i] = coords[i] - eps;
        let em = energy_fn(&cp);
        grad[i] = (ep - em) / (2.0 * eps);
    }
    grad
}
