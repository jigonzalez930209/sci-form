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
use serde::Deserialize;
/// Diagnostic test: dump bounds matrix and compare with RDKit
use std::fs;

#[derive(Deserialize)]
struct RefAtom {
    element: u8,
    x: f64,
    y: f64,
    z: f64,
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
    signs: Vec<f64>,
}

#[derive(Deserialize)]
struct RefMolecule {
    smiles: String,
    atoms: Vec<RefAtom>,
    bonds: Vec<RefBond>,
    #[serde(default)]
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

fn compare_bounds(ref_mol: &RefMolecule, rdkit_file: &str) {
    // Load RDKit bounds matrix
    let rdkit_data = fs::read_to_string(rdkit_file)
        .unwrap_or_else(|_| panic!("Missing RDKit bounds file: {}", rdkit_file));
    let rdkit_mols: Vec<serde_json::Value> = serde_json::from_str(&rdkit_data).unwrap();
    let rdkit_mol = &rdkit_mols[0];
    let rdkit_upper: Vec<Vec<f64>> = serde_json::from_value(rdkit_mol["upper"].clone()).unwrap();
    let rdkit_lower: Vec<Vec<f64>> = serde_json::from_value(rdkit_mol["lower"].clone()).unwrap();
    let n_rdkit = rdkit_mol["n_atoms"].as_u64().unwrap() as usize;

    assert_eq!(ref_mol.atoms.len(), n_rdkit, "Atom count mismatch");

    let mol = build_mol_from_ref(ref_mol);
    let mut bounds = sci_form::distgeom::calculate_bounds_matrix_opts(&mol, true);
    // Apply triangle smoothing (same as the pipeline)
    let smooth_ok = sci_form::distgeom::triangle_smooth_tol(&mut bounds, 0.0);
    eprintln!("Triangle smoothing succeeded: {}", smooth_ok);
    let n = mol.graph.node_count();

    // Compare bounds matrices
    // Our convention: upper triangle (i<j) = upper bound, lower triangle (j>i stored at bounds[(j,i)]) = lower bound
    let mut total_diffs = 0;
    let mut max_upper_diff = 0.0f64;
    let mut max_lower_diff = 0.0f64;
    let mut max_upper_pair = (0, 0);
    let mut max_lower_pair = (0, 0);

    eprintln!(
        "\n=== BOUNDS MATRIX COMPARISON (mol 0: {}) ===",
        ref_mol.smiles
    );
    eprintln!("n_atoms = {}", n);
    eprintln!();

    for i in 0..n {
        for j in (i + 1)..n {
            let our_upper = bounds[(i, j)];
            let our_lower = bounds[(j, i)];
            let rdk_upper = rdkit_upper[i][j];
            let rdk_lower = rdkit_lower[i][j];

            let du = (our_upper - rdk_upper).abs();
            let dl = (our_lower - rdk_lower).abs();

            if du > max_upper_diff {
                max_upper_diff = du;
                max_upper_pair = (i, j);
            }
            if dl > max_lower_diff {
                max_lower_diff = dl;
                max_lower_pair = (i, j);
            }

            if du > 0.001 || dl > 0.001 {
                total_diffs += 1;
                if total_diffs <= 40 {
                    let ai = &ref_mol.atoms[i];
                    let aj = &ref_mol.atoms[j];
                    eprintln!(
                        "DIFF [{:2}][{:2}] elem={}-{}  upper: ours={:.6} rdk={:.6} (d={:.6})  lower: ours={:.6} rdk={:.6} (d={:.6})",
                        i, j, ai.element, aj.element, our_upper, rdk_upper, du, our_lower, rdk_lower, dl
                    );
                }
            }
        }
    }

    eprintln!();
    eprintln!("Total pairs with diff > 0.001: {}", total_diffs);
    eprintln!(
        "Max upper diff: {:.6} at [{},{}]",
        max_upper_diff, max_upper_pair.0, max_upper_pair.1
    );
    eprintln!(
        "Max lower diff: {:.6} at [{},{}]",
        max_lower_diff, max_lower_pair.0, max_lower_pair.1
    );

    // Also dump specific 1-2 bounds for comparison
    eprintln!();
    eprintln!("=== 1-2 bonds (ours vs RDKit) ===");
    for bond in &ref_mol.bonds {
        let (i, j) = if bond.start < bond.end {
            (bond.start, bond.end)
        } else {
            (bond.end, bond.start)
        };
        let our_upper = bounds[(i, j)];
        let our_lower = bounds[(j, i)];
        let rdk_upper = rdkit_upper[i][j];
        let rdk_lower = rdkit_lower[i][j];
        let du = (our_upper - rdk_upper).abs();
        let dl = (our_lower - rdk_lower).abs();
        if du > 0.0001 || dl > 0.0001 {
            let ai = &ref_mol.atoms[i];
            let aj = &ref_mol.atoms[j];
            eprintln!(
                "  bond [{:2}]-[{:2}] ({}-{}, {}): upper ours={:.6} rdk={:.6}  lower ours={:.6} rdk={:.6}",
                i, j, ai.element, aj.element, bond.order,
                our_upper, rdk_upper, our_lower, rdk_lower
            );
        }
    }
}

#[test]
fn test_compare_bounds_mol0() {
    let ref_data = fs::read_to_string("tests/fixtures/gdb20_reference.json").unwrap();
    let ref_mols: Vec<RefMolecule> = serde_json::from_str(&ref_data).unwrap();
    compare_bounds(&ref_mols[0], "/tmp/rdkit_bounds_0.json");
}

/// Finds the first failing molecule by index in gdb20_reference.json
/// That matches the SMILES C#CC1C(C(O)(CC)CN)CCCC(O)C12CCCO2
#[test]
fn test_compare_bounds_fail1() {
    let ref_data = fs::read_to_string("tests/fixtures/gdb20_reference.json").unwrap();
    let ref_mols: Vec<RefMolecule> = serde_json::from_str(&ref_data).unwrap();
    let target = "C#CC1C(C(O)(CC)CN)CCCC(O)C12CCCO2";
    let ref_mol = ref_mols
        .iter()
        .find(|m| m.smiles == target)
        .unwrap_or_else(|| panic!("SMILES '{}' not found in gdb20_reference.json", target));
    compare_bounds(ref_mol, "/tmp/rdkit_bounds_fail1.json");
}
