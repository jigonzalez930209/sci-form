use crate::graph::Molecule;
use nalgebra::DMatrix;

pub fn calc_chiral_volume(
    idx1: usize,
    idx2: usize,
    idx3: usize,
    idx4: usize,
    coords: &DMatrix<f32>,
) -> f32 {
    let dim = coords.ncols();
    assert!(dim >= 3);
    let v1 = nalgebra::Vector3::new(
        coords[(idx1, 0)] - coords[(idx4, 0)],
        coords[(idx1, 1)] - coords[(idx4, 1)],
        coords[(idx1, 2)] - coords[(idx4, 2)],
    );
    let v2 = nalgebra::Vector3::new(
        coords[(idx2, 0)] - coords[(idx4, 0)],
        coords[(idx2, 1)] - coords[(idx4, 1)],
        coords[(idx2, 2)] - coords[(idx4, 2)],
    );
    let v3 = nalgebra::Vector3::new(
        coords[(idx3, 0)] - coords[(idx4, 0)],
        coords[(idx3, 1)] - coords[(idx4, 1)],
        coords[(idx3, 2)] - coords[(idx4, 2)],
    );
    v1.dot(&v2.cross(&v3))
}

pub fn identify_chiral_sets(mol: &Molecule) -> Vec<crate::forcefield::bounds_ff::ChiralSet> {
    let mut sets = Vec::new();
    for i in 0..mol.graph.node_count() {
        let ni = petgraph::graph::NodeIndex::new(i);
        let atom = &mol.graph[ni];
        if atom.chiral_tag != crate::graph::ChiralType::Unspecified {
            let mut neighbors: Vec<_> = mol.graph.neighbors(ni).map(|n| n.index()).collect();
            if neighbors.len() == 4 {
                neighbors.sort();
                let (lower, upper) = match atom.chiral_tag {
                    crate::graph::ChiralType::TetrahedralCW => (2.0, 50.0),
                    crate::graph::ChiralType::TetrahedralCCW => (-50.0, -2.0),
                    _ => (0.0, 0.0),
                };
                if lower != 0.0 || upper != 0.0 {
                    sets.push(crate::forcefield::bounds_ff::ChiralSet {
                        center: i,
                        neighbors: [neighbors[0], neighbors[1], neighbors[2], neighbors[3]],
                        lower_vol: lower,
                        upper_vol: upper,
                    });
                }
            }
        }
    }
    sets
}
