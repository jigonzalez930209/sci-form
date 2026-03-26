//! Debug: find the actual attempt number used by the torsion-aware
//! conformer pipeline.
//!
//! Usage: `LOG_ATTEMPTS=1 cargo run --release --bin find_real_attempt -- 1 17`
//!
//! Category: debug

/// Find the actual attempt number used by generate_3d_conformer_with_torsions
/// Run: LOG_ATTEMPTS=1 cargo run --release --bin find_real_attempt -- 1 17
use serde::Deserialize;
use std::env;

#[derive(Deserialize)]
#[allow(dead_code)]
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
            element: a.element,
            position: nalgebra::Vector3::zeros(),
            charge: 0.0,
            formal_charge: a.formal_charge,
            hybridization: hyb,
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
        mol.add_bond(
            nidx[b.start],
            nidx[b.end],
            sci_form::graph::Bond {
                order: ord,
                stereo: sci_form::graph::BondStereo::None,
            },
        );
    }
    mol
}

fn build_torsions(t: &[RefTorsion]) -> Vec<sci_form::forcefield::etkdg_3d::M6TorsionContrib> {
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

fn main() {
    let data: Vec<RefMolecule> = serde_json::from_str(
        &std::fs::read_to_string("tests/fixtures/gdb20_reference.json").unwrap(),
    )
    .unwrap();

    let indices: Vec<usize> = env::args().skip(1).map(|s| s.parse().unwrap()).collect();
    let indices = if indices.is_empty() {
        vec![1, 17]
    } else {
        indices
    };

    for idx in indices {
        let r = &data[idx];
        let mol = build_mol(r);
        let torsions = build_torsions(&r.torsions);

        eprintln!("=== mol[{}]: {} ===", idx, r.smiles);
        let result = sci_form::conformer::generate_3d_conformer_with_torsions(&mol, 42, &torsions);
        match result {
            Ok(_) => eprintln!("  Result: SUCCESS"),
            Err(e) => eprintln!("  Result: FAIL ({})", e),
        }
        eprintln!();
    }
}
