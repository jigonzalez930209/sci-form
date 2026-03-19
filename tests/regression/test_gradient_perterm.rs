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
//! Per-term gradient verification to isolate which force field term has the bug.
//!
//! Run:  cargo test --release --test test_gradient_perterm -- --nocapture

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

/// Compute energy from ONLY torsion terms
fn torsion_only_energy(
    coords: &[f64],
    n: usize,
    ff: &sci_form::forcefield::etkdg_3d::Etkdg3DFF,
) -> f64 {
    let c = |atom: usize, d: usize| -> f64 { coords[atom * 3 + d] };
    let mut energy = 0.0f64;
    for tc in &ff.torsion_contribs {
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
        if d1 < 1e-10 || d2 < 1e-10 {
            continue;
        }
        let n1 = [t1[0] / d1, t1[1] / d1, t1[2] / d1];
        let n2 = [t2[0] / d2, t2[1] / d2, t2[2] / d2];
        let cos_phi = (n1[0] * n2[0] + n1[1] * n2[1] + n1[2] * n2[2]).clamp(-1.0, 1.0);
        let cp2 = cos_phi * cos_phi;
        let cp3 = cos_phi * cp2;
        let cp4 = cos_phi * cp3;
        let cp5 = cos_phi * cp4;
        let ss = &tc.signs;
        let vv = &tc.v;
        energy += vv[0] * (1.0 + ss[0] * cos_phi)
            + vv[1] * (1.0 + ss[1] * (2.0 * cp2 - 1.0))
            + vv[2] * (1.0 + ss[2] * (4.0 * cp3 - 3.0 * cos_phi))
            + vv[3] * (1.0 + ss[3] * (8.0 * cp4 - 8.0 * cp2 + 1.0))
            + vv[4] * (1.0 + ss[4] * (16.0 * cp5 - 20.0 * cp3 + 5.0 * cos_phi))
            + vv[5] * (1.0 + ss[5] * (32.0 * cp5 * cos_phi - 48.0 * cp4 + 18.0 * cp2 - 1.0));
    }
    energy
}

/// Compute energy from ONLY inversion terms
fn inversion_only_energy(
    coords: &[f64],
    n: usize,
    ff: &sci_form::forcefield::etkdg_3d::Etkdg3DFF,
) -> f64 {
    let c = |atom: usize, d: usize| -> f64 { coords[atom * 3 + d] };
    let mut energy = 0.0f64;
    for ic in &ff.inversion_contribs {
        let rji = [
            c(ic.at1, 0) - c(ic.at2, 0),
            c(ic.at1, 1) - c(ic.at2, 1),
            c(ic.at1, 2) - c(ic.at2, 2),
        ];
        let rjk = [
            c(ic.at3, 0) - c(ic.at2, 0),
            c(ic.at3, 1) - c(ic.at2, 1),
            c(ic.at3, 2) - c(ic.at2, 2),
        ];
        let rjl = [
            c(ic.at4, 0) - c(ic.at2, 0),
            c(ic.at4, 1) - c(ic.at2, 1),
            c(ic.at4, 2) - c(ic.at2, 2),
        ];
        let n_vec = [
            rji[1] * rjk[2] - rji[2] * rjk[1],
            rji[2] * rjk[0] - rji[0] * rjk[2],
            rji[0] * rjk[1] - rji[1] * rjk[0],
        ];
        let n_len = (n_vec[0] * n_vec[0] + n_vec[1] * n_vec[1] + n_vec[2] * n_vec[2]).sqrt();
        let rjl_len = (rjl[0] * rjl[0] + rjl[1] * rjl[1] + rjl[2] * rjl[2]).sqrt();
        if n_len < 1e-8 || rjl_len < 1e-8 {
            continue;
        }
        let sin_y = (n_vec[0] * rjl[0] + n_vec[1] * rjl[1] + n_vec[2] * rjl[2]) / (n_len * rjl_len);
        let sin_y = sin_y.clamp(-1.0, 1.0);
        let cos_y_sq = 1.0 - sin_y * sin_y;
        let cos2w = 2.0 * sin_y * sin_y - 1.0; // cos(2*asin(sinY)) = 1 - 2*sin²Y
        energy += ff.oop_k * ic.force_constant * (ic.c0 + ic.c1 * sin_y + ic.c2 * cos2w);
    }
    energy
}

/// Compute energy from ONLY distance constraint terms (all three groups)
fn dist_only_energy(
    coords: &[f64],
    n: usize,
    ff: &sci_form::forcefield::etkdg_3d::Etkdg3DFF,
) -> f64 {
    let c = |atom: usize, d: usize| -> f64 { coords[atom * 3 + d] };
    let mut energy = 0.0f64;
    for dc in ff
        .dist_12
        .iter()
        .chain(ff.dist_13.iter())
        .chain(ff.dist_long.iter())
    {
        let dx = c(dc.i, 0) - c(dc.j, 0);
        let dy = c(dc.i, 1) - c(dc.j, 1);
        let dz = c(dc.i, 2) - c(dc.j, 2);
        let d = (dx * dx + dy * dy + dz * dz).sqrt();
        let term = if d < dc.min_len {
            d - dc.min_len
        } else if d > dc.max_len {
            d - dc.max_len
        } else {
            0.0
        };
        energy += dc.k * term * term;
    }
    energy
}

/// Compute energy from ONLY angle constraint terms
fn angle_only_energy(
    coords: &[f64],
    n: usize,
    ff: &sci_form::forcefield::etkdg_3d::Etkdg3DFF,
) -> f64 {
    let c = |atom: usize, d: usize| -> f64 { coords[atom * 3 + d] };
    let mut energy = 0.0f64;
    for ac in &ff.angle_constraints {
        let p1 = [c(ac.i, 0), c(ac.i, 1), c(ac.i, 2)];
        let p2 = [c(ac.j, 0), c(ac.j, 1), c(ac.j, 2)];
        let p3 = [c(ac.k, 0), c(ac.k, 1), c(ac.k, 2)];
        let r1 = [p1[0] - p2[0], p1[1] - p2[1], p1[2] - p2[2]];
        let r2 = [p3[0] - p2[0], p3[1] - p2[1], p3[2] - p2[2]];
        let l1 = (r1[0] * r1[0] + r1[1] * r1[1] + r1[2] * r1[2]).sqrt();
        let l2 = (r2[0] * r2[0] + r2[1] * r2[1] + r2[2] * r2[2]).sqrt();
        if l1 < 1e-8 || l2 < 1e-8 {
            continue;
        }
        let cos_theta =
            ((r1[0] * r2[0] + r1[1] * r2[1] + r1[2] * r2[2]) / (l1 * l2)).clamp(-1.0, 1.0);
        let theta_deg = cos_theta.acos() * 180.0 / std::f64::consts::PI;
        let term = if theta_deg < ac.min_deg {
            theta_deg - ac.min_deg
        } else if theta_deg > ac.max_deg {
            theta_deg - ac.max_deg
        } else {
            0.0
        };
        energy += ac.force_k * term * term;
    }
    energy
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

#[test]
fn test_gradient_perterm() {
    let ref_data = fs::read_to_string("tests/fixtures/gdb20_reference.json").unwrap();
    let ref_mols: Vec<RefMolecule> = serde_json::from_str(&ref_data).unwrap();

    // Pick molecules that had BAD gradients: indices 10, 13, 47, 50, 89, 93
    let test_indices = [10, 13, 47, 50, 89, 93];

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

        let ff = sci_form::forcefield::etkdg_3d::build_etkdg_3d_ff_with_torsions(
            &mol,
            &coords_mat,
            &bounds,
            &csd_torsions,
        );

        println!("\n=== mol[{}] {} (n={}) ===", mol_idx, ref_mol.smiles, n);
        println!(
            "  FF: {} torsions, {} inversions, {} angles, {} dist_12, {} dist_13, {} dist_long",
            ff.torsion_contribs.len(),
            ff.inversion_contribs.len(),
            ff.angle_constraints.len(),
            ff.dist_12.len(),
            ff.dist_13.len(),
            ff.dist_long.len()
        );

        // Check total gradient first
        let anal_grad =
            sci_form::forcefield::etkdg_3d::etkdg_3d_gradient_f64(&coords, n, &mol, &ff);
        let num_grad = numerical_gradient(
            &|c| sci_form::forcefield::etkdg_3d::etkdg_3d_energy_f64(c, n, &mol, &ff),
            &coords,
            eps,
        );

        // Find worst component
        let mut worst_i = 0;
        let mut worst_rel = 0.0f64;
        for i in 0..(n * 3) {
            let denom = anal_grad[i].abs().max(num_grad[i].abs()).max(1e-8);
            let rel = (anal_grad[i] - num_grad[i]).abs() / denom;
            if rel > worst_rel {
                worst_rel = rel;
                worst_i = i;
            }
        }
        let atom = worst_i / 3;
        let dim = ["x", "y", "z"][worst_i % 3];
        println!(
            "  TOTAL worst: component={} (atom {} {}) anal={:.8e} num={:.8e} rel_err={:.4e}",
            worst_i, atom, dim, anal_grad[worst_i], num_grad[worst_i], worst_rel
        );

        // Now check each term type at the WORST component
        let check_term = |name: &str, energy_fn: &dyn Fn(&[f64]) -> f64| {
            let num_g = {
                let mut cp = coords.clone();
                cp[worst_i] = coords[worst_i] + eps;
                let ep = energy_fn(&cp);
                cp[worst_i] = coords[worst_i] - eps;
                let em = energy_fn(&cp);
                (ep - em) / (2.0 * eps)
            };
            println!("    {}: num_grad[{}]={:.8e}", name, worst_i, num_g);
        };

        let ff_ref = &ff;
        check_term("Torsions", &|c| torsion_only_energy(c, n, ff_ref));
        check_term("Inversions", &|c| inversion_only_energy(c, n, ff_ref));
        check_term("Distances", &|c| dist_only_energy(c, n, ff_ref));
        check_term("Angles", &|c| angle_only_energy(c, n, ff_ref));

        // Also check the FULL numerical gradient decomposed
        let torsion_num = {
            let mut cp = coords.clone();
            cp[worst_i] = coords[worst_i] + eps;
            let ep = torsion_only_energy(&cp, n, &ff);
            cp[worst_i] = coords[worst_i] - eps;
            let em = torsion_only_energy(&cp, n, &ff);
            (ep - em) / (2.0 * eps)
        };
        let inv_num = {
            let mut cp = coords.clone();
            cp[worst_i] = coords[worst_i] + eps;
            let ep = inversion_only_energy(&cp, n, &ff);
            cp[worst_i] = coords[worst_i] - eps;
            let em = inversion_only_energy(&cp, n, &ff);
            (ep - em) / (2.0 * eps)
        };
        let dist_num = {
            let mut cp = coords.clone();
            cp[worst_i] = coords[worst_i] + eps;
            let ep = dist_only_energy(&cp, n, &ff);
            cp[worst_i] = coords[worst_i] - eps;
            let em = dist_only_energy(&cp, n, &ff);
            (ep - em) / (2.0 * eps)
        };
        let angle_num = {
            let mut cp = coords.clone();
            cp[worst_i] = coords[worst_i] + eps;
            let ep = angle_only_energy(&cp, n, &ff);
            cp[worst_i] = coords[worst_i] - eps;
            let em = angle_only_energy(&cp, n, &ff);
            (ep - em) / (2.0 * eps)
        };
        let sum_parts = torsion_num + inv_num + dist_num + angle_num;
        println!(
            "  Sum-of-parts num: {:.8e} (torsion={:.4e} inv={:.4e} dist={:.4e} angle={:.4e})",
            sum_parts, torsion_num, inv_num, dist_num, angle_num
        );
        println!("  Full num:         {:.8e}", num_grad[worst_i]);
        println!("  Analytical:       {:.8e}", anal_grad[worst_i]);
    }
}
