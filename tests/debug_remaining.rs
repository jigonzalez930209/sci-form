/// Debug test: compare coordinates for the two remaining failing molecules
/// Dumps pairwise distances to compare with RDKit reference

#[test]
fn debug_remaining_failures() {
    use std::fs;
    use serde::Deserialize;
    
    #[derive(Deserialize)]
    struct OracleAtom {
        element: u8,
        x: f32,
        y: f32,
        z: f32,
        formal_charge: i8,
        hybridization: String,
    }

    #[derive(Deserialize)]
    struct OracleBond {
        start: usize,
        end: usize,
        order: String,
    }

    #[derive(Deserialize)]
    struct OracleMolecule {
        smiles: String,
        atoms: Vec<OracleAtom>,
        bonds: Vec<OracleBond>,
    }
    
    let data = fs::read_to_string("tests/fixtures/reference_coords_no_mmff.json").unwrap();
    let mut molecules: Vec<OracleMolecule> = serde_json::from_str(&data).unwrap();
    
    use rand::seq::SliceRandom;
    use rand::SeedableRng;
    let mut rng = rand::rngs::StdRng::seed_from_u64(123);
    molecules.shuffle(&mut rng);
    molecules.truncate(100);
    
    let targets = vec!["CC1(O)CN2CC12", "CC1(O)C2OCC12C"];
    
    for (idx, mol) in molecules.iter().enumerate() {
        if !targets.iter().any(|t| mol.smiles == *t) { continue; }
        
        let n = mol.atoms.len();
        let mut our_mol = sci_form::graph::Molecule::new(&mol.smiles);
        let mut node_indices = Vec::new();
        
        for atom in &mol.atoms {
            let hybridization = match atom.hybridization.as_str() {
                "S" => sci_form::graph::Hybridization::SP,
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
            node_indices.push(our_mol.add_atom(new_atom));
        }

        for bond in &mol.bonds {
            let order = match bond.order.as_str() {
                "DOUBLE" => sci_form::graph::BondOrder::Double,
                "TRIPLE" => sci_form::graph::BondOrder::Triple,
                "AROMATIC" => sci_form::graph::BondOrder::Aromatic,
                _ => sci_form::graph::BondOrder::Single,
            };
            our_mol.add_bond(
                node_indices[bond.start],
                node_indices[bond.end],
                sci_form::graph::Bond { order, stereo: sci_form::graph::BondStereo::None },
            );
        }
        
        match sci_form::conformer::generate_3d_conformer(&our_mol, 42) {
            Ok(coords) => {
                println!("\n=== Mol {} ({}) ===", idx, mol.smiles);
                
                // Find the worst pairwise distance differences
                let mut diffs: Vec<(usize, usize, f32, f32, f32)> = Vec::new();
                for i in 0..n {
                    for j in (i+1)..n {
                        let dr = ((mol.atoms[i].x - mol.atoms[j].x).powi(2)
                            + (mol.atoms[i].y - mol.atoms[j].y).powi(2)
                            + (mol.atoms[i].z - mol.atoms[j].z).powi(2)).sqrt();
                        let du = ((coords[(i,0)] - coords[(j,0)]).powi(2)
                            + (coords[(i,1)] - coords[(j,1)]).powi(2)
                            + (coords[(i,2)] - coords[(j,2)]).powi(2)).sqrt();
                        diffs.push((i, j, dr, du, (dr - du).abs()));
                    }
                }
                diffs.sort_by(|a, b| b.4.partial_cmp(&a.4).unwrap());
                
                println!("Top 20 worst pairwise distance diffs:");
                for &(i, j, ref_d, our_d, diff) in diffs.iter().take(20) {
                    let ei = mol.atoms[i].element;
                    let ej = mol.atoms[j].element;
                    let elem_name = |e: u8| match e { 1 => "H", 6 => "C", 7 => "N", 8 => "O", _ => "?" };
                    println!("  ({:2}:{},{:2}:{}): ref={:.4} ours={:.4} diff={:.4}",
                        i, elem_name(ei), j, elem_name(ej), ref_d, our_d, diff);
                }
                
                // Print our coords
                println!("\nOur coords:");
                for i in 0..n {
                    println!("  atom {:2}: ({:8.4}, {:8.4}, {:8.4})",
                        i, coords[(i,0)], coords[(i,1)], coords[(i,2)]);
                }
                
                // Print ref coords
                println!("\nRef coords:");
                for i in 0..n {
                    println!("  atom {:2}: ({:8.4}, {:8.4}, {:8.4})",
                        i, mol.atoms[i].x, mol.atoms[i].y, mol.atoms[i].z);
                }
            }
            Err(e) => {
                println!("Mol {} ({}): FAILED: {}", idx, mol.smiles, e);
            }
        }
    }
}
