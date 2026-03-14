#![allow(unused_imports, unused_variables, dead_code, clippy::unnecessary_cast, clippy::needless_range_loop, clippy::manual_repeat_n, clippy::manual_str_repeat, clippy::manual_is_multiple_of, clippy::redundant_field_names, clippy::useless_vec, clippy::single_range_in_vec_init)]
//! Energy comparison diagnostic: for failing molecules (RMSD > 0.5),
//! compare our FF energy at our final coords vs at RDKit's reference coords.
//!
//! Run:  GDB20_LIMIT=1000 cargo test --release --test test_energy_comparison -- --nocapture

use nalgebra::DMatrix;
use rayon::prelude::*;
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

fn pairwise_rmsd(our_coords: &DMatrix<f32>, ref_atoms: &[RefAtom]) -> f32 {
    let n = our_coords.nrows();
    let mut sum_sq = 0.0f64;
    let mut count = 0u32;
    for i in 0..n {
        for j in (i + 1)..n {
            let our_d = {
                let dx = (our_coords[(i, 0)] - our_coords[(j, 0)]) as f64;
                let dy = (our_coords[(i, 1)] - our_coords[(j, 1)]) as f64;
                let dz = (our_coords[(i, 2)] - our_coords[(j, 2)]) as f64;
                (dx * dx + dy * dy + dz * dz).sqrt()
            };
            let ref_d = {
                let dx = ref_atoms[i].x - ref_atoms[j].x;
                let dy = ref_atoms[i].y - ref_atoms[j].y;
                let dz = ref_atoms[i].z - ref_atoms[j].z;
                (dx * dx + dy * dy + dz * dz).sqrt()
            };
            sum_sq += (our_d - ref_d) * (our_d - ref_d);
            count += 1;
        }
    }
    (sum_sq / count.max(1) as f64).sqrt() as f32
}

#[test]
fn test_energy_comparison() {
    let ref_data = fs::read_to_string("tests/fixtures/gdb20_reference.json").unwrap();
    let ref_mols: Vec<RefMolecule> = serde_json::from_str(&ref_data).unwrap();

    let limit: usize = std::env::var("GDB20_LIMIT")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(1000);
    let ref_mols = &ref_mols[..limit.min(ref_mols.len())];

    struct EnergyResult {
        mol_idx: usize,
        smiles: String,
        rmsd: f32,
        our_energy: f64,
        rdkit_energy: f64,
        our_lower: bool,
    }

    let results: Vec<Option<EnergyResult>> = ref_mols
        .par_iter()
        .enumerate()
        .map(|(mol_idx, ref_mol)| {
            let mol = build_mol_from_ref(ref_mol);
            let csd_torsions = build_csd_torsions(&ref_mol.torsions);
            let n = mol.graph.node_count();

            // Generate our conformer
            let our_coords = match sci_form::conformer::generate_3d_conformer_with_torsions(
                &mol,
                42,
                &csd_torsions,
            ) {
                Ok(c) => c,
                Err(_) => return None,
            };

            // Compute RMSD
            let rmsd = pairwise_rmsd(&our_coords, &ref_mol.atoms);
            if rmsd <= 0.5 {
                return None;
            } // Only analyze failing molecules

            // Build FF from OUR pre-BFGS coords (as in the pipeline)
            let our_coords_f64 = our_coords.map(|v| v as f64);
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

            // Build FF using OUR final coords (this is NOT exactly what the pipeline does —
            // the pipeline builds FF from pre-BFGS coords, but this approximation is fine
            // for comparing energy levels)
            let ff = sci_form::forcefield::etkdg_3d::build_etkdg_3d_ff_with_torsions(
                &mol,
                &our_coords_f64,
                &bounds,
                &csd_torsions,
            );

            // Evaluate energy at our final coords
            let our_flat: Vec<f64> = (0..n)
                .flat_map(|a| {
                    vec![
                        our_coords_f64[(a, 0)],
                        our_coords_f64[(a, 1)],
                        our_coords_f64[(a, 2)],
                    ]
                })
                .collect();
            let our_energy =
                sci_form::forcefield::etkdg_3d::etkdg_3d_energy_f64(&our_flat, n, &mol, &ff);

            // Evaluate energy at RDKit's reference coords
            let rdkit_flat: Vec<f64> = (0..n)
                .flat_map(|a| vec![ref_mol.atoms[a].x, ref_mol.atoms[a].y, ref_mol.atoms[a].z])
                .collect();
            let rdkit_energy =
                sci_form::forcefield::etkdg_3d::etkdg_3d_energy_f64(&rdkit_flat, n, &mol, &ff);

            Some(EnergyResult {
                mol_idx,
                smiles: ref_mol.smiles.clone(),
                rmsd,
                our_energy,
                rdkit_energy,
                our_lower: our_energy <= rdkit_energy,
            })
        })
        .collect();

    let results: Vec<EnergyResult> = results.into_iter().flatten().collect();
    let total = results.len();
    let our_lower = results.iter().filter(|r| r.our_lower).count();
    let rdkit_lower = total - our_lower;

    println!("\n=== ENERGY COMPARISON FOR HIGH-RMSD MOLECULES ===");
    println!("Molecules with RMSD > 0.5: {}", total);
    println!(
        "Our energy <= RDKit energy (our FF): {} ({:.1}%)",
        our_lower,
        our_lower as f64 / total.max(1) as f64 * 100.0
    );
    println!(
        "RDKit energy < Our energy (our FF): {} ({:.1}%)",
        rdkit_lower,
        rdkit_lower as f64 / total.max(1) as f64 * 100.0
    );

    // Sort by energy difference
    let mut sorted = results;
    sorted.sort_by(|a, b| {
        let diff_a = a.our_energy - a.rdkit_energy;
        let diff_b = b.our_energy - b.rdkit_energy;
        diff_a.partial_cmp(&diff_b).unwrap()
    });

    println!("\n--- Molecules where OUR energy is HIGHER (optimizer stuck) ---");
    for r in sorted.iter().rev().take(20) {
        if !r.our_lower {
            println!(
                "  mol[{}] rmsd={:.4} our_E={:.4} rdkit_E={:.4} diff={:.4} {}",
                r.mol_idx,
                r.rmsd,
                r.our_energy,
                r.rdkit_energy,
                r.our_energy - r.rdkit_energy,
                r.smiles
            );
        }
    }

    println!("\n--- Molecules where OUR energy is LOWER (different basin) ---");
    for r in sorted.iter().take(20) {
        if r.our_lower {
            println!(
                "  mol[{}] rmsd={:.4} our_E={:.4} rdkit_E={:.4} diff={:.4} {}",
                r.mol_idx,
                r.rmsd,
                r.our_energy,
                r.rdkit_energy,
                r.our_energy - r.rdkit_energy,
                r.smiles
            );
        }
    }
}
