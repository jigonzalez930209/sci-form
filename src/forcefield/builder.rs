/// Build a full UFF force field for a given molecule.
///
/// Walks all bonds, angles, torsions, out-of-plane inversions, and non-bonded
/// pairs and adds the appropriate [`ForceFieldContribution`] terms to a fresh
/// [`MolecularForceField`].
use crate::forcefield::atom_typer::{assign_uff_type, is_atom_aromatic};
use crate::forcefield::params::{
    get_uff_params, get_uff_bond_force_constant, get_uff_bond_length,
    get_uff_angle_force_constant, get_uff_torsion_params,
};
use crate::forcefield::traits::MolecularForceField;
use crate::forcefield::uff::{
    UffHarmonicBondStretch, UffAngleBend, UffTorsion, UffInversion, UffLennardJones,
};
use crate::graph::{BondOrder, Hybridization, Molecule};
use petgraph::graph::NodeIndex;
use petgraph::visit::EdgeRef;

/// Assign atom types and look up UFF params for every atom in `mol`.
/// Returns a parallel `Vec<&'static str>` of atom type labels.
fn atom_types(mol: &Molecule) -> Vec<&'static str> {
    let n = mol.graph.node_count();
    (0..n)
        .map(|i| {
            let node = NodeIndex::new(i);
            let atom = &mol.graph[node];
            let aromatic = is_atom_aromatic(mol, i);
            assign_uff_type(atom.element, &atom.hybridization, aromatic)
        })
        .collect()
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

/// Build a complete UFF force field for `mol`, including:
/// - Bond stretching (harmonic)
/// - Angle bending (Fourier)
/// - Torsions (for all i-j-k-l paths)
/// - Out-of-plane inversions (sp2 atoms)
/// - Non-bonded Lennard-Jones (12-6), skipping 1-2 / 1-3, scaling 1-4 by 0.5
pub fn build_uff_force_field(mol: &Molecule) -> MolecularForceField {
    let mut ff = MolecularForceField::new();
    let types = atom_types(mol);
    let n = mol.graph.node_count();
    let topo = topological_distances(mol);

    // ---- Bond stretching ----
    for edge in mol.graph.edge_references() {
        let i = edge.source().index();
        let j = edge.target().index();
        let bo = &mol.graph[edge.id()].order;

        let (pi, pj) = match (get_uff_params(types[i]), get_uff_params(types[j])) {
            (Some(a), Some(b)) => (a, b),
            _ => continue,
        };

        let r0  = get_uff_bond_length(&pi, &pj, bo);
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

    // ---- Angle bending (i–j–k, j = centre) ----
    for j in 0..n {
        let nb_j = &neighbors[j];
        if nb_j.len() < 2 { continue; }
        let pj = match get_uff_params(types[j]) { Some(p) => p, None => continue };

        let theta0_deg = pj.theta0;
        let theta0_rad = theta0_deg.to_radians();

        // coordination_n: 0 = linear (sp), 3 = sp2, 4 = sp3
        let coord_n = match mol.graph[NodeIndex::new(j)].hybridization {
            Hybridization::SP  => 0,
            Hybridization::SP2 => 3,
            _                  => 4,
        };

        for ii in 0..nb_j.len() {
            for kk in (ii + 1)..nb_j.len() {
                let i = nb_j[ii];
                let k = nb_j[kk];

                let pi = match get_uff_params(types[i]) { Some(p) => p, None => continue };
                let pk = match get_uff_params(types[k]) { Some(p) => p, None => continue };

                let bo_ij = mol.graph.find_edge(NodeIndex::new(i), NodeIndex::new(j))
                    .map(|e| mol.graph[e].order.clone())
                    .unwrap_or(BondOrder::Single);
                let bo_jk = mol.graph.find_edge(NodeIndex::new(j), NodeIndex::new(k))
                    .map(|e| mol.graph[e].order.clone())
                    .unwrap_or(BondOrder::Single);

                let r_ij = get_uff_bond_length(&pi, &pj, &bo_ij);
                let r_jk = get_uff_bond_length(&pj, &pk, &bo_jk);
                let k_a  = get_uff_angle_force_constant(&pi, &pk, r_ij, r_jk, theta0_rad);

                ff.insert_dynamic_term(Box::new(UffAngleBend {
                    atom_i_idx: i,
                    atom_j_idx: j,
                    atom_k_idx: k,
                    force_constant_ka: k_a,
                    equilibrium_theta0: theta0_rad,
                    coordination_n: coord_n,
                }));
            }
        }
    }

    // ---- Torsion (i–j–k–l) ----
    for edge_jk in mol.graph.edge_references() {
        let j = edge_jk.source().index();
        let k = edge_jk.target().index();
        let bo_jk = &mol.graph[edge_jk.id()].order;

        let pj = match get_uff_params(types[j]) { Some(p) => p, None => continue };
        let pk = match get_uff_params(types[k]) { Some(p) => p, None => continue };

        let hyb_j = &mol.graph[NodeIndex::new(j)].hybridization;
        let hyb_k = &mol.graph[NodeIndex::new(k)].hybridization;

        let (v_barrier, period_n, cos_phi0) =
            get_uff_torsion_params(hyb_j, hyb_k, pj.v_tors, pk.v_tors, bo_jk);

        if v_barrier < 1e-10 {
            continue;
        }

        for &i in &neighbors[j] {
            if i == k { continue; }
            for &l in &neighbors[k] {
                if l == j || l == i { continue; }
                ff.insert_dynamic_term(Box::new(UffTorsion {
                    atom_i_idx: i,
                    atom_j_idx: j,
                    atom_k_idx: k,
                    atom_l_idx: l,
                    force_constant_v: v_barrier,
                    periodicity_n: period_n,
                    cos_phi0,
                }));
            }
        }
    }

    // ---- Out-of-plane inversions (sp2 centre, 3 substituents) ----
    for j in 0..n {
        if mol.graph[NodeIndex::new(j)].hybridization != Hybridization::SP2 { continue; }
        let nb_j = &neighbors[j];
        if nb_j.len() != 3 { continue; }

        // UFF inversion parameters for sp2: K = 6.0 kcal/mol for C_R/N_R, else 6.0 flat
        let k_inv = 6.0_f64;

        // c0, c1, c2 for ideal planar psi=0: c0=1, c1=0, c2=-1 → E = K(1 - cos(2*psi))
        let (c0, c1, c2) = (1.0, 0.0, -1.0);

        let (i, k, l) = (nb_j[0], nb_j[1], nb_j[2]);
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

    // ---- Non-bonded Lennard-Jones ----
    // Skip 1-2, 1-3 pairs; scale 1-4 by 0.5; full weight for 1-5+
    for i in 0..n {
        let pi = match get_uff_params(types[i]) { Some(p) => p, None => continue };
        for j in (i + 1)..n {
            let d = topo[i][j];
            if d <= 2 { continue; } // 1-2 and 1-3 excluded

            let pj = match get_uff_params(types[j]) { Some(p) => p, None => continue };

            let r_star = (pi.x1 + pj.x1) * 0.5;
            let eps_full = (pi.d1 * pj.d1).sqrt();
            let epsilon = if d == 3 { eps_full * 0.5 } else { eps_full };

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
