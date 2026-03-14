/// Output per-molecule RMSD for first N molecules
use sci_form::conformer::generate_3d_conformer_with_torsions;
use sci_form::graph::Molecule;
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
struct RefTorsion {
    atoms: Vec<usize>,
    v: Vec<f64>,
    signs: Vec<i32>,
}

#[derive(Deserialize)]
#[allow(dead_code)]
struct RefMolecule {
    smiles: String,
    atoms: Vec<RefAtom>,
    bonds: Vec<RefBond>,
    torsions: Vec<RefTorsion>,
}

fn main() {
    let limit: usize = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(1000);

    let data: Vec<RefMolecule> = serde_json::from_str(
        &std::fs::read_to_string("tests/fixtures/gdb20_reference.json").unwrap(),
    )
    .unwrap();

    let mut results = Vec::new();
    for (_i, entry) in data.iter().enumerate().take(limit.min(data.len())) {

        let mol = build_mol_from_ref(entry);
        let csd_torsions = build_csd_torsions(&entry.torsions);
        let result = generate_3d_conformer_with_torsions(&mol, 42, &csd_torsions);

        let rmsd = match result {
            Ok(coords) => pairwise_rmsd(&coords, &entry.atoms) as f64,
            Err(_) => -1.0,
        };
        results.push(rmsd);
    }

    let json = serde_json::to_string(&results).unwrap();
    println!("{}", json);
}

fn build_mol_from_ref(r: &RefMolecule) -> Molecule {
    use sci_form::graph::{Atom, Bond, BondOrder, BondStereo, Hybridization};
    let mut mol = Molecule::new("");
    for a in &r.atoms {
        let hyb = match a.hybridization.as_str() {
            "SP" => Hybridization::SP,
            "SP2" => Hybridization::SP2,
            "SP3" => Hybridization::SP3,
            "SP3D" => Hybridization::SP3D,
            "SP3D2" => Hybridization::SP3D2,
            _ => Hybridization::Unknown,
        };
        let mut atom = Atom::new(a.element, 0.0, 0.0, 0.0);
        atom.hybridization = hyb;
        atom.formal_charge = a.formal_charge;
        mol.add_atom(atom);
    }
    for b in &r.bonds {
        let order = match b.order.as_str() {
            "SINGLE" => BondOrder::Single,
            "DOUBLE" => BondOrder::Double,
            "TRIPLE" => BondOrder::Triple,
            "AROMATIC" => BondOrder::Aromatic,
            _ => BondOrder::Single,
        };
        mol.add_bond(
            petgraph::graph::NodeIndex::new(b.start),
            petgraph::graph::NodeIndex::new(b.end),
            Bond {
                order,
                stereo: BondStereo::None,
            },
        );
    }
    mol
}

fn build_csd_torsions(
    torsions: &[RefTorsion],
) -> Vec<sci_form::forcefield::etkdg_3d::M6TorsionContrib> {
    torsions
        .iter()
        .map(|t| sci_form::forcefield::etkdg_3d::M6TorsionContrib {
            i: t.atoms[0],
            j: t.atoms[1],
            k: t.atoms[2],
            l: t.atoms[3],
            signs: [
                t.signs[0] as f64,
                t.signs[1] as f64,
                t.signs[2] as f64,
                t.signs[3] as f64,
                t.signs[4] as f64,
                t.signs[5] as f64,
            ],
            v: [t.v[0], t.v[1], t.v[2], t.v[3], t.v[4], t.v[5]],
        })
        .collect()
}

fn pairwise_rmsd(coords: &nalgebra::DMatrix<f32>, ref_atoms: &[RefAtom]) -> f32 {
    let n = ref_atoms.len();
    let mut sq_sum = 0.0f64;
    let mut npairs = 0u64;
    for a in 0..n {
        for b in (a + 1)..n {
            let dr = ((ref_atoms[a].x - ref_atoms[b].x).powi(2)
                + (ref_atoms[a].y - ref_atoms[b].y).powi(2)
                + (ref_atoms[a].z - ref_atoms[b].z).powi(2))
            .sqrt() as f64;
            let du = ((coords[(a, 0)] - coords[(b, 0)]).powi(2)
                + (coords[(a, 1)] - coords[(b, 1)]).powi(2)
                + (coords[(a, 2)] - coords[(b, 2)]).powi(2))
            .sqrt() as f64;
            sq_sum += (dr - du).powi(2);
            npairs += 1;
        }
    }
    if npairs > 0 {
        (sq_sum / npairs as f64).sqrt() as f32
    } else {
        0.0
    }
}
