/// Quick tool to dump which attempt our code uses for specific molecules.
/// Run: cargo run --release --bin find_attempt -- 0 1 10

use sci_form::distgeom::embedding::MinstdRand;
use sci_form::distgeom::bounds::{calculate_bounds_matrix_opts, triangle_smooth_tol};
use sci_form::distgeom::embedding::{pick_rdkit_distances, compute_initial_coords_rdkit};
use sci_form::forcefield::bounds_ff::{minimize_bfgs_rdkit, ChiralSet};
use sci_form::distgeom::chirality::identify_chiral_sets;
use sci_form::distgeom::validation::{identify_tetrahedral_centers, check_tetrahedral_centers, check_chiral_centers};

use serde::Deserialize;
use std::env;

#[derive(Deserialize)]
struct RefAtom {
    element: u8,
    x: f32, y: f32, z: f32,
    formal_charge: i8,
    hybridization: String,
}
#[derive(Deserialize)]
struct RefBond { start: usize, end: usize, order: String }
#[derive(Deserialize)]
struct RefTorsion { atoms: Vec<usize>, v: Vec<f64>, signs: Vec<i32> }
#[derive(Deserialize)]
struct RefMolecule { smiles: String, atoms: Vec<RefAtom>, bonds: Vec<RefBond>, torsions: Vec<RefTorsion> }

fn build_mol(r: &RefMolecule) -> sci_form::graph::Molecule {
    let mut mol = sci_form::graph::Molecule::new(&r.smiles);
    let mut nidx = Vec::new();
    for a in &r.atoms {
        let hyb = match a.hybridization.as_str() {
            "SP" => sci_form::graph::Hybridization::SP,
            "SP2" => sci_form::graph::Hybridization::SP2,
            "SP3" => sci_form::graph::Hybridization::SP3,
            "SP3D" => sci_form::graph::Hybridization::SP3D,
            "SP3D2" => sci_form::graph::Hybridization::SP3D2,
            _ => sci_form::graph::Hybridization::Unknown,
        };
        nidx.push(mol.add_atom(sci_form::graph::Atom {
            element: a.element, position: nalgebra::Vector3::zeros(),
            charge: 0.0, formal_charge: a.formal_charge, hybridization: hyb,
            chiral_tag: sci_form::graph::ChiralType::Unspecified,
            explicit_h: if a.element <= 1 { 1 } else { 0 },
        }));
    }
    for b in &r.bonds {
        let ord = match b.order.as_str() {
            "DOUBLE" => sci_form::graph::BondOrder::Double,
            "TRIPLE" => sci_form::graph::BondOrder::Triple,
            "AROMATIC" => sci_form::graph::BondOrder::Aromatic,
            _ => sci_form::graph::BondOrder::Single,
        };
        mol.add_bond(nidx[b.start], nidx[b.end], sci_form::graph::Bond { order: ord, stereo: sci_form::graph::BondStereo::None });
    }
    mol
}

fn build_torsions(t: &[RefTorsion]) -> Vec<sci_form::forcefield::etkdg_3d::M6TorsionContrib> {
    t.iter().filter_map(|t| {
        if t.atoms.len() < 4 || t.v.len() < 6 || t.signs.len() < 6 { return None; }
        let mut signs = [0.0f64; 6];
        let mut v = [0.0f64; 6];
        for k in 0..6 { signs[k] = t.signs[k] as f64; v[k] = t.v[k]; }
        Some(sci_form::forcefield::etkdg_3d::M6TorsionContrib { i: t.atoms[0], j: t.atoms[1], k: t.atoms[2], l: t.atoms[3], signs, v })
    }).collect()
}

fn main() {
    let data: Vec<RefMolecule> = serde_json::from_str(
        &std::fs::read_to_string("tests/fixtures/gdb20_reference.json").unwrap()
    ).unwrap();

    let indices: Vec<usize> = env::args().skip(1).map(|s| s.parse().unwrap()).collect();
    let indices = if indices.is_empty() { (0..20).collect() } else { indices };

    for idx in indices {
        let r = &data[idx];
        let mol = build_mol(r);
        let n = mol.graph.node_count();
        let torsions = build_torsions(&r.torsions);

        // Reproduce the conformer generation step by step
        let bounds = {
            let raw = calculate_bounds_matrix_opts(&mol, true);
            let mut b = raw;
            if triangle_smooth_tol(&mut b, 0.0) { b }
            else {
                let raw2 = calculate_bounds_matrix_opts(&mol, false);
                let mut b2 = raw2;
                triangle_smooth_tol(&mut b2, 0.0);
                b2
            }
        };
        let chiral_sets = identify_chiral_sets(&mol);
        let tet_centers = identify_tetrahedral_centers(&mol);
        let use_4d = !chiral_sets.is_empty();
        let embed_dim = if use_4d { 4 } else { 3 };
        let mut rng = MinstdRand::new(42);
        let max_iters = 10 * n;

        let mut attempt_found = None;
        for iter in 0..max_iters {
            let dists = pick_rdkit_distances(&mut rng, &bounds);
            let coords_opt = compute_initial_coords_rdkit(&mut rng, &dists, embed_dim);
            let mut coords = match coords_opt {
                Some(c) => c,
                None => { continue; }
            };

            // Bounds FF #1
            let energy_pre = sci_form::conformer::compute_total_bounds_energy_f64(&coords, &bounds, &chiral_sets, 5.0, 0.1, 1.0);
            if energy_pre > 1e-5 {
                let mut need_more = 1;
                while need_more != 0 {
                    need_more = minimize_bfgs_rdkit(&mut coords, &bounds, &chiral_sets, 400, 1e-3, 5.0, 0.1, 1.0);
                }
            }

            // Energy check
            let energy = sci_form::conformer::compute_total_bounds_energy_f64(&coords, &bounds, &chiral_sets, 5.0, 0.1, 1.0);
            if energy / n as f64 >= 0.05 { continue; }

            // Tetrahedral check
            if !check_tetrahedral_centers(&coords, &tet_centers) { continue; }

            // Chiral check  
            if !chiral_sets.is_empty() && !check_chiral_centers(&coords, &chiral_sets) { continue; }

            attempt_found = Some(iter);
            break;
        }

        match attempt_found {
            Some(a) => println!("mol[{}]: smiles={}, our_attempt={}", idx, r.smiles, a),
            None => println!("mol[{}]: smiles={}, FAILED (no valid attempt)", idx, r.smiles),
        }
    }
}
