//! Debug: identify which atom pairs cause the largest distance errors in failing molecules.

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
            position: nalgebra::Vector3::zeros(),
            charge: 0.0,
            formal_charge: atom.formal_charge,
            hybridization,
            chiral_tag: sci_form::graph::ChiralType::Unspecified,
            explicit_h: if atom.element == 1 || atom.element == 0 { 1 } else { 0 },
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
            sci_form::graph::Bond { order, stereo: sci_form::graph::BondStereo::None },
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
            if t.atoms.len() < 4 || t.v.len() < 6 || t.signs.len() < 6 { return None; }
            let mut signs = [0.0f64; 6];
            let mut v = [0.0f64; 6];
            for k in 0..6 {
                signs[k] = t.signs[k] as f64;
                v[k] = t.v[k] as f64;
            }
            Some(sci_form::forcefield::etkdg_3d::M6TorsionContrib {
                i: t.atoms[0], j: t.atoms[1], k: t.atoms[2], l: t.atoms[3], signs, v,
            })
        })
        .collect()
}

#[test]
fn test_diagnose_failures() {
    let ref_data = fs::read_to_string("tests/fixtures/gdb20_reference.json")
        .expect("Run scripts/generate_gdb20_reference.py first");
    let ref_mols: Vec<RefMolecule> = serde_json::from_str(&ref_data).unwrap();

    let limit = std::env::var("GDB20_LIMIT")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(ref_mols.len());
    let ref_mols = &ref_mols[..limit.min(ref_mols.len())];

    println!("\n=== FAILURE DIAGNOSIS ===");

    let mut cases_analyzed = 0;
    for ref_mol in ref_mols {
        let mol = build_mol_from_ref(ref_mol);
        let csd_torsions = build_csd_torsions(&ref_mol.torsions);
        let result = sci_form::conformer::generate_3d_conformer_with_torsions(&mol, 42, &csd_torsions);

        let coords = match result {
            Ok(c) => c,
            Err(_) => continue,
        };

        let n = ref_mol.atoms.len();
        let mut sq_sum = 0.0f64;
        let mut npairs = 0u64;
        for a in 0..n {
            for b in (a + 1)..n {
                let dr = ((ref_mol.atoms[a].x - ref_mol.atoms[b].x).powi(2)
                    + (ref_mol.atoms[a].y - ref_mol.atoms[b].y).powi(2)
                    + (ref_mol.atoms[a].z - ref_mol.atoms[b].z).powi(2))
                .sqrt() as f64;
                let du = ((coords[(a, 0)] - coords[(b, 0)]).powi(2)
                    + (coords[(a, 1)] - coords[(b, 1)]).powi(2)
                    + (coords[(a, 2)] - coords[(b, 2)]).powi(2))
                .sqrt() as f64;
                sq_sum += (dr - du).powi(2);
                npairs += 1;
            }
        }
        let rmsd = if npairs > 0 {
            (sq_sum / npairs as f64).sqrt()
        } else {
            0.0
        };

        if rmsd < 0.5 { continue; }
        if cases_analyzed >= 5 { break; }
        cases_analyzed += 1;

        println!("\n--- {} (RMSD={:.3}) ---", ref_mol.smiles, rmsd);
        println!("  Heavy atoms: {}", ref_mol.atoms.iter().filter(|a| a.element != 1).count());
        println!("  CSD torsions: {}", ref_mol.torsions.len());
        
        // Find top 10 atom pairs with largest distance differences
        let mut pair_diffs: Vec<(usize, usize, f64, f64, f64)> = Vec::new();
        for a in 0..n {
            for b in (a + 1)..n {
                let dr = ((ref_mol.atoms[a].x - ref_mol.atoms[b].x).powi(2)
                    + (ref_mol.atoms[a].y - ref_mol.atoms[b].y).powi(2)
                    + (ref_mol.atoms[a].z - ref_mol.atoms[b].z).powi(2))
                .sqrt() as f64;
                let du = ((coords[(a, 0)] - coords[(b, 0)]).powi(2)
                    + (coords[(a, 1)] - coords[(b, 1)]).powi(2)
                    + (coords[(a, 2)] - coords[(b, 2)]).powi(2))
                .sqrt() as f64;
                pair_diffs.push((a, b, dr, du, (dr - du).abs()));
            }
        }
        pair_diffs.sort_by(|a, b| b.4.partial_cmp(&a.4).unwrap());

        println!("  Top 10 worst atom pairs:");
        for &(a, b, dr, du, diff) in pair_diffs.iter().take(10) {
            let ea = ref_mol.atoms[a].element;
            let eb = ref_mol.atoms[b].element;
            let ha = &ref_mol.atoms[a].hybridization;
            let hb = &ref_mol.atoms[b].hybridization;
            println!(
                "    ({:2},{:2}) e{}({})-e{}({}): ref={:.3} ours={:.3} Δ={:.3}",
                a, b, ea, ha, eb, hb, dr, du, diff
            );
        }

        // Count how many pairs involve heavy atoms only
        let heavy_only: Vec<_> = pair_diffs.iter()
            .filter(|&&(a, b, _, _, _)| ref_mol.atoms[a].element != 1 && ref_mol.atoms[b].element != 1)
            .collect();
        let heavy_sq_sum: f64 = heavy_only.iter().map(|&&(_, _, dr, du, _)| (dr - du).powi(2)).sum();
        let heavy_rmsd = (heavy_sq_sum / heavy_only.len() as f64).sqrt();
        println!("  Heavy-only RMSD: {:.3}", heavy_rmsd);

        // Identify which fragment is displaced
        // Check if large errors cluster around specific atoms
        let mut atom_error = vec![0.0f64; n];
        let mut atom_count = vec![0u32; n];
        for &(a, b, _, _, diff) in &pair_diffs {
            atom_error[a] += diff;
            atom_count[a] += 1;
            atom_error[b] += diff;
            atom_count[b] += 1;
        }
        let mut avg_err: Vec<(usize, f64)> = (0..n)
            .map(|i| (i, if atom_count[i] > 0 { atom_error[i] / atom_count[i] as f64 } else { 0.0 }))
            .collect();
        avg_err.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());

        println!("  Atoms with highest avg error:");
        for &(idx, err) in avg_err.iter().take(8) {
            let e = ref_mol.atoms[idx].element;
            let h = &ref_mol.atoms[idx].hybridization;
            println!("    atom {:2} e{}({}): avg_err={:.3}", idx, e, h, err);
        }
    }
}
