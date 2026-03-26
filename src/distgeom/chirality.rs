use crate::graph::{BondOrder, Molecule};
use nalgebra::DMatrix;
use petgraph::graph::NodeIndex;

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

pub fn calc_chiral_volume_f64(
    idx1: usize,
    idx2: usize,
    idx3: usize,
    idx4: usize,
    coords: &DMatrix<f64>,
) -> f64 {
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
        let ni = NodeIndex::new(i);
        let atom = &mol.graph[ni];
        if atom.chiral_tag != crate::graph::ChiralType::Unspecified {
            let neighbors: Vec<_> = mol.graph.neighbors(ni).map(|n| n.index()).collect();
            if neighbors.len() == 4 {
                let canonical = canonicalize_chiral_neighbors(mol, ni, &neighbors);
                let parity_odd = permutation_is_odd(&neighbors, &canonical);
                let (mut lower, mut upper) = match atom.chiral_tag {
                    crate::graph::ChiralType::TetrahedralCW => (2.0, 50.0),
                    crate::graph::ChiralType::TetrahedralCCW => (-50.0, -2.0),
                    _ => (0.0, 0.0),
                };
                if parity_odd {
                    (lower, upper) = (-upper, -lower);
                }
                if lower != 0.0 || upper != 0.0 {
                    sets.push(crate::forcefield::bounds_ff::ChiralSet {
                        center: i,
                        neighbors: [canonical[0], canonical[1], canonical[2], canonical[3]],
                        lower_vol: lower,
                        upper_vol: upper,
                    });
                }
            }
        }
    }
    sets
}

fn canonicalize_chiral_neighbors(mol: &Molecule, center: NodeIndex, neighbors: &[usize]) -> Vec<usize> {
    let mut ordered = neighbors.to_vec();
    ordered.sort_by(|&left, &right| neighbor_priority(mol, center, left).cmp(&neighbor_priority(mol, center, right)));
    ordered
}

fn neighbor_priority(mol: &Molecule, center: NodeIndex, neighbor: usize) -> (u8, u8, usize, i8, u8, std::cmp::Reverse<usize>) {
    let node = NodeIndex::new(neighbor);
    let atom = &mol.graph[node];
    let degree = mol.graph.neighbors(node).count();
    let bond_rank = mol
        .graph
        .find_edge(center, node)
        .and_then(|edge| mol.graph.edge_weight(edge))
        .map(|bond| bond_order_rank(bond.order))
        .unwrap_or(0);

    (
        atom.element,
        bond_rank,
        degree,
        atom.formal_charge,
        atom.explicit_h,
        std::cmp::Reverse(neighbor),
    )
}

fn bond_order_rank(order: BondOrder) -> u8 {
    match order {
        BondOrder::Triple => 4,
        BondOrder::Double => 3,
        BondOrder::Aromatic => 2,
        BondOrder::Single => 1,
        BondOrder::Unknown => 0,
    }
}

fn permutation_is_odd(original: &[usize], reordered: &[usize]) -> bool {
    let mut inversions = 0usize;
    let positions: Vec<usize> = reordered
        .iter()
        .map(|idx| {
            original
                .iter()
                .position(|candidate| candidate == idx)
                .expect("reordered neighbor missing from original list")
        })
        .collect();

    for i in 0..positions.len() {
        for j in (i + 1)..positions.len() {
            if positions[i] > positions[j] {
                inversions += 1;
            }
        }
    }

    inversions % 2 == 1
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::{Atom, Bond, BondStereo, ChiralType, Hybridization};

    fn build_tetrahedral(order: [u8; 4], chiral_tag: ChiralType) -> Molecule {
        let mut mol = Molecule::new("tetrahedral");
        let center = mol.add_atom(Atom {
            element: 6,
            position: nalgebra::Vector3::zeros(),
            charge: 0.0,
            formal_charge: 0,
            hybridization: Hybridization::SP3,
            chiral_tag,
            explicit_h: 0,
        });

        let neighbors: Vec<_> = order
            .iter()
            .map(|&z| {
                mol.add_atom(Atom {
                    element: z,
                    position: nalgebra::Vector3::zeros(),
                    charge: 0.0,
                    formal_charge: 0,
                    hybridization: Hybridization::SP3,
                    chiral_tag: ChiralType::Unspecified,
                    explicit_h: 0,
                })
            })
            .collect();

        for neighbor in neighbors {
            mol.add_bond(
                center,
                neighbor,
                Bond {
                    order: BondOrder::Single,
                    stereo: BondStereo::None,
                },
            );
        }
        mol
    }

    #[test]
    fn test_permutation_parity() {
        assert!(!permutation_is_odd(&[1, 2, 3, 4], &[1, 2, 3, 4]));
        assert!(permutation_is_odd(&[1, 2, 3, 4], &[2, 1, 3, 4]));
        assert!(!permutation_is_odd(&[1, 2, 3, 4], &[4, 3, 2, 1]));
    }

    #[test]
    fn test_identify_chiral_sets_uses_chemical_priority_order() {
        let mol = build_tetrahedral([9, 53, 17, 35], ChiralType::TetrahedralCW);
        let sets = identify_chiral_sets(&mol);

        assert_eq!(sets.len(), 1);
        let elements: Vec<u8> = sets[0]
            .neighbors
            .iter()
            .map(|&idx| mol.graph[NodeIndex::new(idx)].element)
            .collect();
        assert_eq!(elements, vec![53, 35, 17, 9]);
    }

    #[test]
    fn test_odd_reordering_flips_volume_sign_bounds() {
        let mol = build_tetrahedral([53, 9, 17, 35], ChiralType::TetrahedralCW);
        let sets = identify_chiral_sets(&mol);

        assert_eq!(sets.len(), 1);
        assert!(sets[0].lower_vol < 0.0);
        assert!(sets[0].upper_vol < 0.0);
    }
}
