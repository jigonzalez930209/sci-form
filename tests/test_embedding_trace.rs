use sci_form::distgeom::bounds::{calculate_bounds_matrix_opts, triangle_smooth_tol};
/// Diagnostic test: compare embedding intermediate values with Python reference
/// Run: cargo test --release --test test_embedding_trace -- --nocapture
use sci_form::distgeom::embedding::{
    compute_initial_coords_rdkit, pick_rdkit_distances, MinstdRand,
};
use serde::Deserialize;

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
struct RefMolecule {
    smiles: String,
    atoms: Vec<RefAtom>,
    bonds: Vec<RefBond>,
}

fn build_mol_from_ref(ref_mol: &RefMolecule) -> sci_form::graph::Molecule {
    let mut mol = sci_form::graph::Molecule::new(&ref_mol.smiles);
    let mut node_indices = Vec::new();
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

#[test]
fn test_embedding_trace_mol0() {
    let data: Vec<RefMolecule> = serde_json::from_str(
        &std::fs::read_to_string("tests/fixtures/gdb20_reference.json").unwrap(),
    )
    .unwrap();

    let ref_mol = &data[0];
    let n = ref_mol.atoms.len();

    eprintln!("=== Molecule 0: {} (n={}) ===", ref_mol.smiles, n);

    // Build bounds matrix
    let mol = build_mol_from_ref(ref_mol);
    let bounds = {
        let raw = calculate_bounds_matrix_opts(&mol, true);
        let mut b = raw;
        if triangle_smooth_tol(&mut b, 0.0) {
            b
        } else {
            let raw2 = calculate_bounds_matrix_opts(&mol, false);
            let mut b2 = raw2;
            triangle_smooth_tol(&mut b2, 0.0);
            b2
        }
    };

    // Initialize RNG same as conformer.rs
    let mut rng = MinstdRand::new(42);

    // Pick distances
    let dists = pick_rdkit_distances(&mut rng, &bounds);

    // Print first 10 distances
    eprintln!("\nFirst 10 picked distances:");
    let mut count = 0;
    for i in 1..n {
        for j in 0..i {
            if count < 10 {
                eprintln!("  ({},{}) -> {:.15}", i, j, dists[(i, j)]);
            }
            count += 1;
        }
    }

    // Compute initial coords (no chiral centers -> 3D)
    let embed_dim = 3;
    let coords_opt = compute_initial_coords_rdkit(&mut rng, &dists, embed_dim);

    match coords_opt {
        Some(coords) => {
            eprintln!("\nInitial coords (first 5 atoms):");
            for i in 0..5.min(n) {
                eprintln!(
                    "  atom[{}] = ({:.15}, {:.15}, {:.15})",
                    i,
                    coords[(i, 0)],
                    coords[(i, 1)],
                    coords[(i, 2)]
                );
            }

            // Python reference values for mol 0
            let py_coords: Vec<(f64, f64, f64)> = vec![
                (4.494861601277167, -0.365490219152904, -0.445524138813395),
                (4.058597298215237, 0.548756239671342, 0.104282875642444),
                (2.866596711452423, -1.298654014442495, -0.161334685200206),
                (2.765518713341003, -0.676153352164740, 0.213391849969879),
                (1.381417622688304, -1.516520310891142, 0.637962622230838),
            ];

            eprintln!("\nComparison with Python (first 5 atoms):");
            let mut max_diff = 0.0f64;
            for i in 0..5 {
                let dx = (coords[(i, 0)] - py_coords[i].0).abs();
                let dy = (coords[(i, 1)] - py_coords[i].1).abs();
                let dz = (coords[(i, 2)] - py_coords[i].2).abs();
                let diff = dx.max(dy).max(dz);
                if diff > max_diff {
                    max_diff = diff;
                }
                eprintln!("  atom[{}]: dx={:.2e}, dy={:.2e}, dz={:.2e}", i, dx, dy, dz);
            }
            eprintln!("Max coordinate diff (first 5 atoms): {:.2e}", max_diff);
        }
        None => {
            eprintln!("\nEmbedding FAILED (attempt 0) as expected");
            eprintln!("(Python also failed: eigenvalue[2]=-79.07 is negative,");
            eprintln!(" but randNegEig=true so should use random coords for dim 2)");
        }
    }
}
