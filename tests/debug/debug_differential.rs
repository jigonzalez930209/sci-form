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
use sci_form::distgeom;
use sci_form::graph::{Atom, Bond, BondOrder, BondStereo, ChiralType, Hybridization, Molecule};
use serde::Deserialize;
use std::fs;

#[derive(Deserialize, Debug)]
struct OracleAtom {
    element: u8,
}

#[derive(Deserialize, Debug)]
struct OracleBond {
    start: usize,
    end: usize,
}

#[derive(Deserialize, Debug)]
struct OracleMolecule {
    smiles: String,
    atoms: Vec<OracleAtom>,
    bonds: Vec<OracleBond>,
}

#[test]
fn test_bounds_matrix() {
    let data = fs::read_to_string("tests/fixtures/reference_coords_500.json").unwrap();
    let molecules: Vec<OracleMolecule> = serde_json::from_str(&data).unwrap();
    let mol = &molecules[47]; // C1COC1

    let mut our_mol = Molecule::new(&mol.smiles);
    let mut node_indices = Vec::new();
    for atom in &mol.atoms {
        let new_atom = Atom {
            element: atom.element,
            position: nalgebra::Vector3::zeros(),
            charge: 0.0,
            formal_charge: 0,
            hybridization: Hybridization::Unknown,
            chiral_tag: ChiralType::Unspecified,
            explicit_h: 0,
        };
        node_indices.push(our_mol.add_atom(new_atom));
    }
    for bond in &mol.bonds {
        our_mol.add_bond(
            node_indices[bond.start],
            node_indices[bond.end],
            Bond {
                order: BondOrder::Single,
                stereo: BondStereo::None,
            },
        );
    }

    let bounds = distgeom::calculate_bounds_matrix(&our_mol);

    println!("Before smooth: Upper bounds >= 1000.0:");
    for i in 0..bounds.nrows() {
        for j in (i + 1)..bounds.ncols() {
            if bounds[(i, j)] >= 1000.0 {
                println!("  {}-{}", i, j);
            }
        }
    }

    let smoothed = distgeom::smooth_bounds_matrix(bounds.clone());
    println!("After smooth: Upper bounds >= 1000.0:");
    for i in 0..smoothed.nrows() {
        for j in (i + 1)..smoothed.ncols() {
            if smoothed[(i, j)] >= 1000.0 {
                println!("  {}-{} = {}", i, j, smoothed[(i, j)]);
            }
        }
    }
}
