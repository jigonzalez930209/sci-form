use rand::SeedableRng;
use serde::Deserialize;
use std::fs;
use petgraph::visit::EdgeRef;

#[derive(Deserialize, Debug)]
struct OracleAtom {
    element: u8,
    x: f32,
    y: f32,
    z: f32,
    formal_charge: i8,
    hybridization: String,
}

#[derive(Deserialize, Debug)]
struct OracleBond {
    start: usize,
    end: usize,
    order: String,
}

#[derive(Deserialize, Debug)]
struct OracleMolecule {
    smiles: String,
    atoms: Vec<OracleAtom>,
    bonds: Vec<OracleBond>,
}

fn build_mol(mol: &OracleMolecule) -> sci_form::graph::Molecule {
    let mut our_mol = sci_form::graph::Molecule::new(&mol.smiles);
    let mut node_indices = Vec::new();
    for atom in &mol.atoms {
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
        node_indices.push(our_mol.add_atom(new_atom));
    }
    for bond in &mol.bonds {
        let order = match bond.order.as_str() {
            "DOUBLE" => sci_form::graph::BondOrder::Double,
            "TRIPLE" => sci_form::graph::BondOrder::Triple,
            "AROMATIC" => sci_form::graph::BondOrder::Aromatic,
            _ => sci_form::graph::BondOrder::Single,
        };
        our_mol.add_bond(node_indices[bond.start], node_indices[bond.end],
            sci_form::graph::Bond { order, stereo: sci_form::graph::BondStereo::None });
    }
    our_mol
}

#[test]
fn debug_embedding_failures() {
    let data = fs::read_to_string("tests/fixtures/reference_coords.json")
        .expect("Should be able to read reference_coords.json");
    let mut molecules: Vec<OracleMolecule> =
        serde_json::from_str(&data).expect("JSON was not well-formatted");

    use rand::seq::SliceRandom;
    let mut rng = rand::rngs::StdRng::seed_from_u64(42); // deterministic
    molecules.shuffle(&mut rng);
    molecules.truncate(30);

    let mut bond_errors: Vec<(f32, f32, f32, u8, u8, String)> = Vec::new(); // (our_bl, rdkit_bl, diff, e1, e2, bond_order)

    for (idx, mol) in molecules.iter().enumerate() {
        let our_mol = build_mol(mol);
        let n = our_mol.graph.node_count();

        // Compare bond lengths: ours vs oracle
        for edge in our_mol.graph.edge_references() {
            let i = edge.source().index();
            let j = edge.target().index();
            let our_bl = sci_form::distgeom::get_bond_length(&our_mol, edge.source(), edge.target()) as f32;
            let dx = mol.atoms[i].x - mol.atoms[j].x;
            let dy = mol.atoms[i].y - mol.atoms[j].y;
            let dz = mol.atoms[i].z - mol.atoms[j].z;
            let rdkit_bl = (dx * dx + dy * dy + dz * dz).sqrt();
            let diff = (our_bl - rdkit_bl).abs();
            let e1 = mol.atoms[i].element.min(mol.atoms[j].element);
            let e2 = mol.atoms[i].element.max(mol.atoms[j].element);
            let order = format!("{:?}", our_mol.graph[edge.id()].order);
            bond_errors.push((our_bl, rdkit_bl, diff, e1, e2, order));
        }
    }

    // Sort by absolute error descending
    bond_errors.sort_by(|a, b| b.2.partial_cmp(&a.2).unwrap());

    println!("=== TOP 20 BOND LENGTH ERRORS ===");
    for (our, rdkit, diff, e1, e2, order) in bond_errors.iter().take(20) {
        println!("  ({:2},{:2}) {} : ours={:.3} rdkit={:.3} diff={:.3}", e1, e2, order, our, rdkit, diff);
    }

    // Group by element pair and bond order, compute average error
    use std::collections::HashMap;
    let mut group_errors: HashMap<String, Vec<f32>> = HashMap::new();
    for (our, rdkit, diff, e1, e2, order) in &bond_errors {
        let key = format!("({},{}){}", e1, e2, order);
        group_errors.entry(key).or_default().push(*diff);
    }
    println!("\n=== AVERAGE BOND LENGTH ERRORS BY TYPE ===");
    let mut groups: Vec<_> = group_errors.iter().collect();
    groups.sort_by(|a, b| {
        let avg_a = a.1.iter().sum::<f32>() / a.1.len() as f32;
        let avg_b = b.1.iter().sum::<f32>() / b.1.len() as f32;
        avg_b.partial_cmp(&avg_a).unwrap()
    });
    for (key, errs) in &groups {
        let avg = errs.iter().sum::<f32>() / errs.len() as f32;
        let max = errs.iter().cloned().fold(0.0f32, f32::max);
        println!("  {:25} count={:3} avg_err={:.4} max_err={:.4}", key, errs.len(), avg, max);
    }
}
