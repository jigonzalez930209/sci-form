/// Build a full UFF force field for a given molecule.
///
/// Walks all bonds, angles, torsions, out-of-plane inversions, and non-bonded
/// pairs and adds the appropriate [`ForceFieldContribution`] terms to a fresh
/// [`MolecularForceField`].
///
/// The assembly mirrors the canonical RDKit builder
/// (`ForceFieldHelpers/UFF/Builder.cpp`) so energies can be cross-validated:
/// - sp2 angles use the trigonal-planar functional (`order = 3`) with the
///   3-/4-membered ring `theta0` corrections,
/// - torsions are only built for non-terminal, non-triple central bonds and
///   their force constants are divided by the number of torsion paths,
/// - inversions add three contributions per sp2 centre (K already /3),
/// - van der Waals uses geometric-mean combining with no 1-4 scaling.
use crate::forcefield::atom_typer::{
    assign_uff_type, atom_total_valence, conjugated_bond_flags, is_atom_aromatic, order_value,
};
use crate::forcefield::params::{
    get_uff_angle_force_constant, get_uff_bond_force_constant, get_uff_bond_length,
    get_uff_inversion_params, get_uff_params, get_uff_torsion_params,
};
use crate::forcefield::traits::MolecularForceField;
use crate::forcefield::uff::{
    UffAngleBend, UffHarmonicBondStretch, UffInversion, UffLennardJones, UffTorsion,
};
use crate::graph::{BondOrder, Hybridization, Molecule};
use crate::rings::sssr::RingInfo;
use petgraph::graph::NodeIndex;
use petgraph::visit::EdgeRef;

fn atom_in_ring_of_size(rings: &[RingInfo], atom: usize, size: usize) -> bool {
    rings
        .iter()
        .any(|r| r.size == size && r.atoms.contains(&atom))
}

/// Build a complete UFF force field for `mol`, including:
/// - Bond stretching (harmonic)
/// - Angle bending (Fourier / order-specific)
/// - Torsions (scaled and special-cased like the canonical builder)
/// - Out-of-plane inversions (three contributions per sp2 centre)
/// - Non-bonded Lennard-Jones (12-6), skipping 1-2 / 1-3 pairs
pub fn build_uff_force_field(mol: &Molecule) -> MolecularForceField {
    let mut ff = MolecularForceField::new();
    let n = mol.graph.node_count();
    let topo = topological_distances(mol);
    let conj_flags = conjugated_bond_flags(mol);

    // ---- Atom typing (RDKit getAtomLabel semantics) ----
    let types: Vec<&'static str> = (0..n)
        .map(|i| {
            let node = NodeIndex::new(i);
            let atom = &mol.graph[node];
            let has_conj = mol.graph.edges(node).any(|e| conj_flags[e.id().index()]);
            assign_uff_type(
                atom.element,
                &atom.hybridization,
                is_atom_aromatic(mol, i),
                has_conj,
                atom_total_valence(mol, i),
            )
        })
        .collect();

    let elements: Vec<u8> = (0..n)
        .map(|i| mol.graph[NodeIndex::new(i)].element)
        .collect();

    // ---- Bond stretching ----
    for edge in mol.graph.edge_references() {
        let i = edge.source().index();
        let j = edge.target().index();
        let bo = &mol.graph[edge.id()].order;

        let (pi, pj) = match (get_uff_params(types[i]), get_uff_params(types[j])) {
            (Some(a), Some(b)) => (a, b),
            _ => continue,
        };

        let r0 = get_uff_bond_length(&pi, &pj, bo);
        let k_b = get_uff_bond_force_constant(&pi, &pj, r0);

        ff.insert_dynamic_term(Box::new(UffHarmonicBondStretch {
            atom_i_idx: i,
            atom_j_idx: j,
            force_constant_kb: k_b,
            equilibrium_r0: r0,
        }));
    }

    // ---- Build adjacency list for angle/torsion iteration ----
    let mut neighbors: Vec<Vec<usize>> = vec![Vec::new(); n];
    for edge in mol.graph.edge_references() {
        let i = edge.source().index();
        let j = edge.target().index();
        neighbors[i].push(j);
        neighbors[j].push(i);
    }

    let rings = crate::rings::compute_sssr(mol).rings;

    // ---- Angle bending (i–j–k, j = centre) ----
    for j in 0..n {
        let nb_j = &neighbors[j];
        if nb_j.len() < 2 {
            continue;
        }
        let pj = match get_uff_params(types[j]) {
            Some(p) => p,
            None => continue,
        };
        let hyb_j = mol.graph[NodeIndex::new(j)].hybridization;

        for ii in 0..nb_j.len() {
            for kk in (ii + 1)..nb_j.len() {
                let i = nb_j[ii];
                let k = nb_j[kk];

                let pi = match get_uff_params(types[i]) {
                    Some(p) => p,
                    None => continue,
                };
                let pk = match get_uff_params(types[k]) {
                    Some(p) => p,
                    None => continue,
                };

                let bo_ij = mol
                    .graph
                    .find_edge(NodeIndex::new(i), NodeIndex::new(j))
                    .map(|e| mol.graph[e].order)
                    .unwrap_or(BondOrder::Single);
                let bo_jk = mol
                    .graph
                    .find_edge(NodeIndex::new(j), NodeIndex::new(k))
                    .map(|e| mol.graph[e].order)
                    .unwrap_or(BondOrder::Single);

                // RDKit angle order + 3/4-membered ring theta0 corrections.
                let mut order = match hyb_j {
                    Hybridization::SP => 1usize,
                    Hybridization::SP2 => 3,
                    Hybridization::SP3D2 => 4,
                    _ => 0,
                };
                let mut theta0_deg = pj.theta0;
                if matches!(hyb_j, Hybridization::SP2) {
                    let j3 = atom_in_ring_of_size(&rings, j, 3);
                    let j4 = atom_in_ring_of_size(&rings, j, 4);
                    if j3 {
                        let i3 = atom_in_ring_of_size(&rings, i, 3);
                        let k3 = atom_in_ring_of_size(&rings, k, 3);
                        if i3 != k3 {
                            order = 30;
                            theta0_deg = 150.0;
                        } else if i3 && k3 {
                            order = 35;
                            theta0_deg = 60.0;
                        }
                    } else if j4 {
                        let i4 = atom_in_ring_of_size(&rings, i, 4);
                        let k4 = atom_in_ring_of_size(&rings, k, 4);
                        if i4 != k4 {
                            order = 40;
                            theta0_deg = 135.0;
                        } else if i4 && k4 {
                            order = 45;
                            theta0_deg = 90.0;
                        }
                    }
                    if order >= 30 {
                        order = 0; // ring hack uses the polynomial with corrected theta0
                    }
                }
                let theta0_rad = theta0_deg.to_radians();

                let r_ij = get_uff_bond_length(&pi, &pj, &bo_ij);
                let r_jk = get_uff_bond_length(&pj, &pk, &bo_jk);
                let k_a = get_uff_angle_force_constant(&pi, &pk, r_ij, r_jk, theta0_rad);

                ff.insert_dynamic_term(Box::new(UffAngleBend {
                    atom_i_idx: i,
                    atom_j_idx: j,
                    atom_k_idx: k,
                    force_constant_ka: k_a,
                    equilibrium_theta0: theta0_rad,
                    order,
                }));
            }
        }
    }

    // ---- Torsion (i–j–k–l) ----
    for edge_jk in mol.graph.edge_references() {
        let j = edge_jk.source().index();
        let k = edge_jk.target().index();
        let bo_jk = &mol.graph[edge_jk.id()].order;

        let (pj, pk) = match (get_uff_params(types[j]), get_uff_params(types[k])) {
            (Some(a), Some(b)) => (a, b),
            _ => continue,
        };
        let hyb_j = mol.graph[NodeIndex::new(j)].hybridization;
        let hyb_k = mol.graph[NodeIndex::new(k)].hybridization;

        // Central bond must be sp2/sp3 on both ends, non-terminal, not triple.
        let both_ok = matches!(hyb_j, Hybridization::SP2 | Hybridization::SP3)
            && matches!(hyb_k, Hybridization::SP2 | Hybridization::SP3)
            && neighbors[j].len() > 1
            && neighbors[k].len() > 1;
        let has_triple = |atom: usize| {
            mol.graph
                .edges(NodeIndex::new(atom))
                .any(|e| e.weight().order == BondOrder::Triple)
        };
        if !both_ok || has_triple(j) || has_triple(k) {
            continue;
        }

        let bo_value = order_value(bo_jk);
        let mut contribs: Vec<(usize, usize, usize, usize, f64, f64, f64)> = Vec::new();
        for &i in &neighbors[j] {
            if i == k {
                continue;
            }
            for &l in &neighbors[k] {
                if l == j || l == i {
                    continue; // includes the 3-membered ring exclusion
                }
                let hyb_i = mol.graph[NodeIndex::new(i)].hybridization;
                let hyb_l = mol.graph[NodeIndex::new(l)].hybridization;
                let end_sp2 =
                    matches!(hyb_i, Hybridization::SP2) || matches!(hyb_l, Hybridization::SP2);
                let (v, period, cos_term) = get_uff_torsion_params(
                    &hyb_j,
                    &hyb_k,
                    elements[j],
                    elements[k],
                    &pj,
                    &pk,
                    bo_value,
                    end_sp2,
                );
                if v.abs() < 1e-12 {
                    continue;
                }
                contribs.push((i, j, k, l, v, period, cos_term));
            }
        }

        // Canonical builder divides each torsion constant by the number of
        // paths about this central bond.
        let scale = contribs.len() as f64;
        for (i, jj, kk, l, v, period, cos_term) in contribs {
            ff.insert_dynamic_term(Box::new(UffTorsion {
                atom_i_idx: i,
                atom_j_idx: jj,
                atom_k_idx: kk,
                atom_l_idx: l,
                force_constant_v: v / scale,
                periodicity_n: period,
                cos_phi0: cos_term,
            }));
        }
    }

    // ---- Out-of-plane inversions ----
    for j in 0..n {
        let element = elements[j];
        if !matches!(element, 6 | 7 | 8 | 15 | 33 | 51 | 83) {
            continue;
        }
        let nb_j = &neighbors[j];
        if nb_j.len() != 3 {
            continue;
        }
        let hyb_j = mol.graph[NodeIndex::new(j)].hybridization;
        if matches!(element, 6 | 7 | 8) && !matches!(hyb_j, Hybridization::SP2) {
            continue;
        }
        let is_c_bound_to_o = element == 6
            && nb_j.iter().any(|&x| {
                elements[x] == 8
                    && matches!(
                        mol.graph[NodeIndex::new(x)].hybridization,
                        Hybridization::SP2
                    )
            });
        let (k_inv, c0, c1, c2) = get_uff_inversion_params(element, is_c_bound_to_o);
        if k_inv.abs() < 1e-12 {
            continue;
        }
        let (a, b, c) = (nb_j[0], nb_j[1], nb_j[2]);
        for (i, k, l) in [(a, b, c), (a, c, b), (b, c, a)] {
            ff.insert_dynamic_term(Box::new(UffInversion {
                idx_i: i,
                idx_j: j,
                idx_k: k,
                idx_l: l,
                k_inv,
                c0,
                c1,
                c2,
            }));
        }
    }

    // ---- Non-bonded Lennard-Jones ----
    // Skip 1-2 and 1-3 pairs; geometric-mean combining (x_ij = sqrt(x_i x_j),
    // D_ij = sqrt(D_i D_j)); no 1-4 scaling (canonical behavior).
    for i in 0..n {
        let pi = match get_uff_params(types[i]) {
            Some(p) => p,
            None => continue,
        };
        for j in (i + 1)..n {
            let d = topo[i][j];
            if d <= 2 {
                continue;
            }

            let pj = match get_uff_params(types[j]) {
                Some(p) => p,
                None => continue,
            };

            let r_star = (pi.x1 * pj.x1).sqrt();
            let epsilon = (pi.d1 * pj.d1).sqrt();

            ff.insert_dynamic_term(Box::new(UffLennardJones {
                atom_i_idx: i,
                atom_j_idx: j,
                r_star,
                epsilon,
            }));
        }
    }

    ff
}

/// Build a topological distance matrix using BFS from each node.
/// Returns `dist[i][j]` = number of bonds between atoms i and j (or 999 if not reachable).
fn topological_distances(mol: &Molecule) -> Vec<Vec<u32>> {
    let n = mol.graph.node_count();
    let mut dist = vec![vec![999u32; n]; n];
    for start in 0..n {
        dist[start][start] = 0;
        let mut queue = std::collections::VecDeque::new();
        queue.push_back(start);
        while let Some(u) = queue.pop_front() {
            let node_u = NodeIndex::new(u);
            for edge in mol.graph.edges(node_u) {
                let v = edge.target().index();
                if dist[start][v] == 999 {
                    dist[start][v] = dist[start][u] + 1;
                    queue.push_back(v);
                }
            }
        }
    }
    dist
}
