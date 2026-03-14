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
//! Isolate WHICH torsion contribution has the gradient bug.
//! Run:  cargo test --release --test test_grad_torsion_detail -- --nocapture

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
fn test_torsion_gradient_detail() {
    let ref_data = fs::read_to_string("tests/fixtures/gdb20_reference.json").unwrap();
    let ref_mols: Vec<RefMolecule> = serde_json::from_str(&ref_data).unwrap();

    let mol_idx = 47; // C#CCNCC1(C)CC2C(C(C)O)CCCCC2(C)C1 — worst case
    let ref_mol = &ref_mols[mol_idx];
    let mol = build_mol_from_ref(ref_mol);
    let csd_torsions = build_csd_torsions(&ref_mol.torsions);
    let n = mol.graph.node_count();
    let eps = 1e-5;

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

    println!(
        "\n=== mol[{}] {} (n={}, {} torsions) ===",
        mol_idx,
        ref_mol.smiles,
        n,
        ff_full.torsion_contribs.len()
    );

    // Test each torsion individually
    for (ti, tc) in ff_full.torsion_contribs.iter().enumerate() {
        let ff_one = sci_form::forcefield::etkdg_3d::Etkdg3DFF {
            torsion_contribs: vec![tc.clone()],
            inversion_contribs: vec![],
            dist_12: vec![],
            angle_constraints: vec![],
            dist_13: vec![],
            dist_long: vec![],
            oop_k: ff_full.oop_k,
            torsion_k_omega: 0.0,
            bounds_force_scaling: ff_full.bounds_force_scaling,
            use_m6_torsions: ff_full.use_m6_torsions,
        };

        let anal = sci_form::forcefield::etkdg_3d::etkdg_3d_gradient_f64(&coords, n, &mol, &ff_one);
        let num = numerical_gradient(
            &|c| sci_form::forcefield::etkdg_3d::etkdg_3d_energy_f64(c, n, &mol, &ff_one),
            &coords,
            eps,
        );

        let mut worst_i = 0;
        let mut worst_rel = 0.0f64;
        for i in 0..(n * 3) {
            let denom = anal[i].abs().max(num[i].abs()).max(1e-8);
            let rel = (anal[i] - num[i]).abs() / denom;
            if rel > worst_rel {
                worst_rel = rel;
                worst_i = i;
            }
        }

        // Compute cos_phi for this torsion
        let c = |atom: usize, d: usize| -> f64 { coords[atom * 3 + d] };
        let r1 = [
            c(tc.i, 0) - c(tc.j, 0),
            c(tc.i, 1) - c(tc.j, 1),
            c(tc.i, 2) - c(tc.j, 2),
        ];
        let r2 = [
            c(tc.k, 0) - c(tc.j, 0),
            c(tc.k, 1) - c(tc.j, 1),
            c(tc.k, 2) - c(tc.j, 2),
        ];
        let r3 = [
            c(tc.j, 0) - c(tc.k, 0),
            c(tc.j, 1) - c(tc.k, 1),
            c(tc.j, 2) - c(tc.k, 2),
        ];
        let r4 = [
            c(tc.l, 0) - c(tc.k, 0),
            c(tc.l, 1) - c(tc.k, 1),
            c(tc.l, 2) - c(tc.k, 2),
        ];
        let t1 = [
            r1[1] * r2[2] - r1[2] * r2[1],
            r1[2] * r2[0] - r1[0] * r2[2],
            r1[0] * r2[1] - r1[1] * r2[0],
        ];
        let t2 = [
            r3[1] * r4[2] - r3[2] * r4[1],
            r3[2] * r4[0] - r3[0] * r4[2],
            r3[0] * r4[1] - r3[1] * r4[0],
        ];
        let d1 = (t1[0] * t1[0] + t1[1] * t1[1] + t1[2] * t1[2]).sqrt();
        let d2 = (t2[0] * t2[0] + t2[1] * t2[1] + t2[2] * t2[2]).sqrt();
        let cos_phi = if d1 > 1e-10 && d2 > 1e-10 {
            let n1 = [t1[0] / d1, t1[1] / d1, t1[2] / d1];
            let n2 = [t2[0] / d2, t2[1] / d2, t2[2] / d2];
            (n1[0] * n2[0] + n1[1] * n2[1] + n1[2] * n2[2]).clamp(-1.0, 1.0)
        } else {
            999.0
        };
        let sin_phi = (1.0 - cos_phi * cos_phi).max(0.0).sqrt();

        let atom = worst_i / 3;
        let dim = ["x", "y", "z"][worst_i % 3];
        if worst_rel > 1e-3 {
            println!("  [BAD] torsion[{}] ({}-{}-{}-{}) cos={:.6} sin={:.6} V=[{:.4},{:.4},{:.4},{:.4},{:.4},{:.4}] s=[{:.0},{:.0},{:.0},{:.0},{:.0},{:.0}]",
                ti, tc.i, tc.j, tc.k, tc.l,
                cos_phi, sin_phi,
                tc.v[0], tc.v[1], tc.v[2], tc.v[3], tc.v[4], tc.v[5],
                tc.signs[0], tc.signs[1], tc.signs[2], tc.signs[3], tc.signs[4], tc.signs[5]);
            println!(
                "        worst: atom {} {} A={:.6e} N={:.6e} rel={:.4e}",
                atom, dim, anal[worst_i], num[worst_i], worst_rel
            );
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
