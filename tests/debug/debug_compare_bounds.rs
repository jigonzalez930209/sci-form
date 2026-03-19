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
//! Debug test to compare bounds matrix values with RDKit.
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
fn test_compare_bounds() {
    let ref_data = fs::read_to_string("tests/fixtures/debug_single_mol.json")
        .expect("Run the debug script first");
    let ref_mols: Vec<RefMolecule> = serde_json::from_str(&ref_data).unwrap();
    let ref_mol = &ref_mols[0];

    println!("\nMolecule: {}", ref_mol.smiles);
    println!("Atoms: {}", ref_mol.atoms.len());

    let mol = build_mol_from_ref(ref_mol);
    let bounds = sci_form::distgeom::calculate_bounds_matrix_opts(&mol, true);
    let n = ref_mol.atoms.len();

    // Print 1-2 bond lengths
    println!("\n=== Our 1-2 Bond Lengths ===");
    for bond in &ref_mol.bonds {
        let (i, j) = if bond.start < bond.end {
            (bond.start, bond.end)
        } else {
            (bond.end, bond.start)
        };
        let lower = bounds[(j, i)];
        let upper = bounds[(i, j)];
        let mid = (lower + upper) / 2.0;
        let ei = ref_mol.atoms[i].element;
        let ej = ref_mol.atoms[j].element;
        let hi = &ref_mol.atoms[i].hybridization;
        let hj = &ref_mol.atoms[j].hybridization;
        println!(
            "  ({:2},{:2}) e{}({})-e{}({}) [{}]: {:.6} ± {:.6}",
            i,
            j,
            ei,
            hi,
            ej,
            hj,
            bond.order,
            mid,
            (upper - lower) / 2.0
        );
    }

    // Load RDKit bounds
    let rdkit_data: serde_json::Value = serde_json::from_str(
        &fs::read_to_string("/tmp/rdkit_bounds_phenylacetylene.json").unwrap(),
    )
    .unwrap();
    let rdkit_bounds = rdkit_data["bounds"].as_array().unwrap();

    // Compare all pairs
    println!("\n=== Bounds Comparison (largest differences) ===");
    let mut diffs: Vec<(usize, usize, f64, f64, f64, f64)> = Vec::new();
    for entry in rdkit_bounds {
        let i = entry["i"].as_u64().unwrap() as usize;
        let j = entry["j"].as_u64().unwrap() as usize;
        let rdkit_lower = entry["lower"].as_f64().unwrap();
        let rdkit_upper = entry["upper"].as_f64().unwrap();

        let our_lower = bounds[(j, i)];
        let our_upper = bounds[(i, j)];

        let diff_lower = (our_lower - rdkit_lower).abs();
        let diff_upper = (our_upper - rdkit_upper).abs();
        let max_diff = diff_lower.max(diff_upper);

        diffs.push((
            i,
            j,
            diff_lower,
            diff_upper,
            our_lower - rdkit_lower,
            our_upper - rdkit_upper,
        ));
    }

    // Sort by max difference
    diffs.sort_by(|a, b| {
        let ma = a.2.max(a.3);
        let mb = b.2.max(b.3);
        mb.partial_cmp(&ma).unwrap()
    });

    println!(
        "  {:>5} | {:>12} {:>12} | {:>12} {:>12} | {:>10} {:>10}",
        "pair", "our_lower", "rdkit_lower", "our_upper", "rdkit_upper", "Δ_lower", "Δ_upper"
    );
    println!("  {}", "-".repeat(90));

    for (idx, &(i, j, dl, du, sl, su)) in diffs.iter().enumerate().take(30) {
        let ei = ref_mol.atoms[i].element;
        let ej = ref_mol.atoms[j].element;
        let rdkit_entry = rdkit_bounds
            .iter()
            .find(|e| e["i"].as_u64().unwrap() == i as u64 && e["j"].as_u64().unwrap() == j as u64)
            .unwrap();
        let rl = rdkit_entry["lower"].as_f64().unwrap();
        let ru = rdkit_entry["upper"].as_f64().unwrap();
        let ol = bounds[(j, i)];
        let ou = bounds[(i, j)];
        println!(
            "  ({:2},{:2}) | {:12.6} {:12.6} | {:12.6} {:12.6} | {:+10.6} {:+10.6}",
            i, j, ol, rl, ou, ru, sl, su
        );
    }

    // Also compute smoothed bounds
    let mut smoothed = bounds.clone();
    let ok = sci_form::distgeom::triangle_smooth_tol(&mut smoothed, 0.0);
    println!("\nTriangle smoothing: {}", if ok { "OK" } else { "FAILED" });

    // Compare smoothed bounds
    println!("\n=== Smoothed Bounds Comparison (largest differences) ===");
    let rdkit_smoothed: serde_json::Value = serde_json::from_str(
        &fs::read_to_string("/tmp/rdkit_bounds_phenylacetylene.json").unwrap(),
    )
    .unwrap();
    // Note: rdkit_bounds IS already smoothed (GetMoleculeBoundsMatrix returns smoothed)

    let mut sdiffs: Vec<(usize, usize, f64, f64)> = Vec::new();
    for entry in rdkit_bounds {
        let i = entry["i"].as_u64().unwrap() as usize;
        let j = entry["j"].as_u64().unwrap() as usize;
        let rdkit_lower = entry["lower"].as_f64().unwrap();
        let rdkit_upper = entry["upper"].as_f64().unwrap();

        let our_lower = smoothed[(j, i)];
        let our_upper = smoothed[(i, j)];

        let diff_lower = our_lower - rdkit_lower;
        let diff_upper = our_upper - rdkit_upper;

        sdiffs.push((i, j, diff_lower, diff_upper));
    }

    sdiffs.sort_by(|a, b| {
        let ma = a.2.abs().max(a.3.abs());
        let mb = b.2.abs().max(b.3.abs());
        mb.partial_cmp(&ma).unwrap()
    });

    println!(
        "  {:>5} | {:>12} {:>12} | {:>10} {:>10}",
        "pair", "Δ_lower", "Δ_upper", "our_lower", "our_upper"
    );
    println!("  {}", "-".repeat(65));
    for &(i, j, dl, du) in sdiffs.iter().take(30) {
        let ol = smoothed[(j, i)];
        let ou = smoothed[(i, j)];
        println!(
            "  ({:2},{:2}) | {:+12.6} {:+12.6} | {:10.6} {:10.6}",
            i, j, dl, du, ol, ou
        );
    }
}
