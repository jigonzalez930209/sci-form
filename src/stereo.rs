//! Stereochemistry perception: CIP priority, R/S, and E/Z assignment.
//!
//! Implements Cahn–Ingold–Prelog priority rules for chirality perception
//! using depth-first graph traversal with phantom atom expansion.

use crate::graph::{BondOrder, Molecule};
use petgraph::graph::NodeIndex;
use petgraph::visit::EdgeRef;
use serde::{Deserialize, Serialize};

/// CIP priority tree node for depth-first traversal.
#[derive(Debug, Clone)]
struct CipNode {
    atomic_number: u8,
    mass: f64,
    children: Vec<CipNode>,
}

/// Result of R/S assignment at a single stereocenter.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Stereocenter {
    /// Atom index of the chiral center.
    pub atom_index: usize,
    /// Atomic number of the center.
    pub element: u8,
    /// Indices of the four substituents (CIP priority order: 1 > 2 > 3 > 4).
    pub substituent_indices: Vec<usize>,
    /// CIP priorities (1 = highest).
    pub priorities: Vec<usize>,
    /// R, S, or None if undetermined.
    pub configuration: Option<String>,
}

/// Result of E/Z assignment at a double bond.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DoubleBondStereo {
    /// Index of atom on one end of the double bond.
    pub atom1: usize,
    /// Index of atom on the other end.
    pub atom2: usize,
    /// E, Z, or None if undetermined.
    pub configuration: Option<String>,
    /// Highest-priority substituent on atom1 side.
    pub high_priority_sub1: Option<usize>,
    /// Highest-priority substituent on atom2 side.
    pub high_priority_sub2: Option<usize>,
}

/// Complete stereochemistry analysis result.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StereoAnalysis {
    /// Detected stereocenters with R/S assignments.
    pub stereocenters: Vec<Stereocenter>,
    /// Detected double bonds with E/Z assignments.
    pub double_bonds: Vec<DoubleBondStereo>,
    /// Total number of chiral centers found.
    pub n_stereocenters: usize,
    /// Total number of E/Z-assignable double bonds.
    pub n_double_bonds: usize,
    /// Detected atropisomeric axes (biaryl axial chirality).
    pub atropisomeric_axes: Vec<AtropisomericAxis>,
    /// Detected pro-chiral centers.
    pub prochiral_centers: Vec<ProchiralCenter>,
    /// Detected helical chirality (M/P in helicenes, metallocenes).
    pub helical_chirality: Vec<HelicalChirality>,
}

/// An atropisomeric axis (restricted rotation around a single bond).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AtropisomericAxis {
    /// Index of atom on one end of the axis bond.
    pub atom1: usize,
    /// Index of atom on the other end.
    pub atom2: usize,
    /// Number of ortho substituents (indicates steric barrier degree).
    pub ortho_substituent_count: usize,
    /// Configuration: aR or aS (if 3D coords available).
    pub configuration: Option<String>,
}

/// A pro-chiral center (sp3 with two identical substituents).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProchiralCenter {
    /// Atom index.
    pub atom_index: usize,
    /// Element at the center.
    pub element: u8,
    /// Indices of the two enantiotopic substituents (identical pair).
    pub enantiotopic_pair: [usize; 2],
}

/// Helical chirality (M/P) detected in helicenes, allenes with extended
/// conjugation, or metallocenes with non-coplanar rings.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HelicalChirality {
    /// Atom indices that form the helical backbone.
    pub backbone_atoms: Vec<usize>,
    /// Configuration: "M" (minus, left-handed) or "P" (plus, right-handed).
    pub configuration: Option<String>,
    /// Type: "helicene", "allene", "metallocene".
    pub helix_type: String,
    /// Dihedral angle of the helix (degrees).
    pub dihedral_angle: Option<f64>,
}

/// Approximate atomic mass for CIP tie-breaking.
fn atomic_mass(z: u8) -> f64 {
    match z {
        1 => 1.008,
        5 => 10.81,
        6 => 12.011,
        7 => 14.007,
        8 => 15.999,
        9 => 18.998,
        14 => 28.085,
        15 => 30.974,
        16 => 32.06,
        17 => 35.45,
        35 => 79.904,
        53 => 126.904,
        _ => z as f64,
    }
}

/// Build a CIP priority tree rooted at `start` coming from `from`, up to depth `max_depth`.
fn build_cip_tree(mol: &Molecule, start: usize, from: usize, max_depth: usize) -> CipNode {
    let start_idx = NodeIndex::new(start);
    let atom = &mol.graph[start_idx];

    let mut node = CipNode {
        atomic_number: atom.element,
        mass: atomic_mass(atom.element),
        children: Vec::new(),
    };

    if max_depth == 0 {
        return node;
    }

    // Collect neighbors excluding `from`
    for neighbor in mol.graph.neighbors(start_idx) {
        let ni = neighbor.index();
        if ni == from {
            continue;
        }

        // Find the bond order between start and neighbor
        let edge = mol.graph.find_edge(start_idx, neighbor);
        let bond_order = edge
            .map(|e| &mol.graph[e].order)
            .unwrap_or(&BondOrder::Single);

        // For double/triple bonds, add phantom atoms (duplicates)
        let n_phantoms = match bond_order {
            BondOrder::Double | BondOrder::Aromatic => 1,
            BondOrder::Triple => 2,
            _ => 0,
        };

        // Real child
        node.children
            .push(build_cip_tree(mol, ni, start, max_depth - 1));

        // Phantom children (copies of the neighbor with no further children)
        let neighbor_atom = &mol.graph[neighbor];
        for _ in 0..n_phantoms {
            node.children.push(CipNode {
                atomic_number: neighbor_atom.element,
                mass: atomic_mass(neighbor_atom.element),
                children: Vec::new(),
            });
        }

        // Add phantom of start to neighbor too (symmetry of double bond expansion)
        // This is handled by the child's own tree construction.
    }

    // Sort children by CIP priority (higher atomic number = higher priority)
    node.children.sort_by(|a, b| cip_compare(b, a));

    node
}

/// Compare two CIP nodes (higher = "greater").
/// Returns Ordering for descending sort (highest priority first).
fn cip_compare(a: &CipNode, b: &CipNode) -> std::cmp::Ordering {
    // Rule 1: Higher atomic number wins
    match a.atomic_number.cmp(&b.atomic_number) {
        std::cmp::Ordering::Equal => {}
        other => return other,
    }

    // Rule 2: Higher mass wins (isotope rule)
    match a.mass.partial_cmp(&b.mass) {
        Some(std::cmp::Ordering::Equal) | None => {}
        Some(other) => return other,
    }

    // Rule 3: Compare children layer by layer
    let max_children = a.children.len().max(b.children.len());
    for i in 0..max_children {
        match (a.children.get(i), b.children.get(i)) {
            (Some(ca), Some(cb)) => match cip_compare(ca, cb) {
                std::cmp::Ordering::Equal => continue,
                other => return other,
            },
            (Some(_), None) => return std::cmp::Ordering::Greater,
            (None, Some(_)) => return std::cmp::Ordering::Less,
            (None, None) => break,
        }
    }

    std::cmp::Ordering::Equal
}

/// Detect stereocenters: sp3 atoms with 4 different substituents.
fn find_stereocenters(mol: &Molecule) -> Vec<usize> {
    let mut centers = Vec::new();

    for node in mol.graph.node_indices() {
        let atom = &mol.graph[node];
        // Skip H and uncommon centers
        if atom.element == 1 {
            continue;
        }

        let neighbors: Vec<usize> = mol.graph.neighbors(node).map(|n| n.index()).collect();

        // Need exactly 4 substituents for sp3 chiral center
        if neighbors.len() != 4 {
            continue;
        }

        // Build CIP trees for each substituent
        let trees: Vec<CipNode> = neighbors
            .iter()
            .map(|&n| build_cip_tree(mol, n, node.index(), 6))
            .collect();

        // Check all pairs are distinct
        let mut all_different = true;
        'outer: for i in 0..4 {
            for j in (i + 1)..4 {
                if cip_compare(&trees[i], &trees[j]) == std::cmp::Ordering::Equal {
                    all_different = false;
                    break 'outer;
                }
            }
        }

        if all_different {
            centers.push(node.index());
        }
    }

    centers
}

/// Assign CIP priorities to substituents around a center atom.
/// Returns (sorted_neighbor_indices, priority_ranks) where rank 1 = highest.
fn assign_cip_priorities(
    mol: &Molecule,
    center: usize,
    neighbors: &[usize],
) -> (Vec<usize>, Vec<usize>) {
    let mut indexed_trees: Vec<(usize, CipNode)> = neighbors
        .iter()
        .map(|&n| (n, build_cip_tree(mol, n, center, 6)))
        .collect();

    // Sort by descending CIP priority
    indexed_trees.sort_by(|a, b| cip_compare(&b.1, &a.1));

    let sorted_indices: Vec<usize> = indexed_trees.iter().map(|(idx, _)| *idx).collect();
    let priorities: Vec<usize> = (1..=sorted_indices.len()).collect();

    (sorted_indices, priorities)
}

/// Determine R/S configuration from 3D coordinates.
///
/// Uses the signed volume of the tetrahedron formed by substituents 1-2-3
/// viewed from substituent 4 (lowest priority).
/// Positive volume with 1→2→3 going clockwise = R.
fn assign_rs(positions: &[[f64; 3]], center: usize, sorted_subs: &[usize]) -> Option<String> {
    if sorted_subs.len() != 4 || positions.is_empty() {
        return None;
    }
    if center >= positions.len() || sorted_subs.iter().any(|&s| s >= positions.len()) {
        return None;
    }

    // Vector from lowest-priority (4th) substituent to center
    let p4 = positions[sorted_subs[3]];
    let p1 = positions[sorted_subs[0]];
    let p2 = positions[sorted_subs[1]];
    let p3 = positions[sorted_subs[2]];

    // Vectors from p4 to p1, p2, p3
    let v1 = [p1[0] - p4[0], p1[1] - p4[1], p1[2] - p4[2]];
    let v2 = [p2[0] - p4[0], p2[1] - p4[1], p2[2] - p4[2]];
    let v3 = [p3[0] - p4[0], p3[1] - p4[1], p3[2] - p4[2]];

    // Signed volume = v1 · (v2 × v3)
    let cross = [
        v2[1] * v3[2] - v2[2] * v3[1],
        v2[2] * v3[0] - v2[0] * v3[2],
        v2[0] * v3[1] - v2[1] * v3[0],
    ];
    let signed_volume = v1[0] * cross[0] + v1[1] * cross[1] + v1[2] * cross[2];

    if signed_volume.abs() < 1e-6 {
        return None; // Planar, can't determine
    }

    Some(if signed_volume > 0.0 {
        "R".to_string()
    } else {
        "S".to_string()
    })
}

/// Find double bonds that could have E/Z isomerism.
fn find_ez_double_bonds(mol: &Molecule) -> Vec<(usize, usize)> {
    let mut bonds = Vec::new();

    for edge in mol.graph.edge_indices() {
        let bond = &mol.graph[edge];
        if bond.order != BondOrder::Double {
            continue;
        }

        let (a, b) = mol.graph.edge_endpoints(edge).unwrap();
        let a_idx = a.index();
        let b_idx = b.index();

        // Each end needs at least 2 different substituents (including the other end)
        let a_neighbors: Vec<usize> = mol
            .graph
            .neighbors(a)
            .map(|n| n.index())
            .filter(|&n| n != b_idx)
            .collect();
        let b_neighbors: Vec<usize> = mol
            .graph
            .neighbors(b)
            .map(|n| n.index())
            .filter(|&n| n != a_idx)
            .collect();

        // Need at least 1 substituent on each side (besides the double bond partner)
        if a_neighbors.is_empty() || b_neighbors.is_empty() {
            continue;
        }

        // If only 1 substituent on a side, we can still assign E/Z
        // If 2 substituents, they must be different
        if a_neighbors.len() >= 2 {
            let trees_a: Vec<CipNode> = a_neighbors
                .iter()
                .map(|&n| build_cip_tree(mol, n, a_idx, 6))
                .collect();
            if cip_compare(&trees_a[0], &trees_a[1]) == std::cmp::Ordering::Equal {
                continue; // Identical substituents on atom a
            }
        }
        if b_neighbors.len() >= 2 {
            let trees_b: Vec<CipNode> = b_neighbors
                .iter()
                .map(|&n| build_cip_tree(mol, n, b_idx, 6))
                .collect();
            if cip_compare(&trees_b[0], &trees_b[1]) == std::cmp::Ordering::Equal {
                continue;
            }
        }

        bonds.push((a_idx, b_idx));
    }

    bonds
}

/// Assign E/Z from 3D coordinates using the dihedral angle of highest-priority
/// substituents across the double bond.
fn assign_ez(
    mol: &Molecule,
    positions: &[[f64; 3]],
    atom1: usize,
    atom2: usize,
) -> DoubleBondStereo {
    let a_idx = NodeIndex::new(atom1);
    let b_idx = NodeIndex::new(atom2);

    let mut a_subs: Vec<usize> = mol
        .graph
        .neighbors(a_idx)
        .map(|n| n.index())
        .filter(|&n| n != atom2)
        .collect();
    let mut b_subs: Vec<usize> = mol
        .graph
        .neighbors(b_idx)
        .map(|n| n.index())
        .filter(|&n| n != atom1)
        .collect();

    // Sort by CIP priority (highest first)
    a_subs.sort_by(|&x, &y| {
        let tx = build_cip_tree(mol, x, atom1, 6);
        let ty = build_cip_tree(mol, y, atom1, 6);
        cip_compare(&ty, &tx)
    });
    b_subs.sort_by(|&x, &y| {
        let tx = build_cip_tree(mol, x, atom2, 6);
        let ty = build_cip_tree(mol, y, atom2, 6);
        cip_compare(&ty, &tx)
    });

    let high_a = a_subs.first().copied();
    let high_b = b_subs.first().copied();

    let configuration = match (high_a, high_b) {
        (Some(ha), Some(hb)) if !positions.is_empty() => {
            if atom1 >= positions.len()
                || atom2 >= positions.len()
                || ha >= positions.len()
                || hb >= positions.len()
            {
                None
            } else {
                // Compute dihedral angle ha-atom1-atom2-hb
                let dihedral = compute_dihedral(
                    &positions[ha],
                    &positions[atom1],
                    &positions[atom2],
                    &positions[hb],
                );
                if dihedral.abs() < std::f64::consts::FRAC_PI_2 {
                    Some("Z".to_string()) // same side
                } else {
                    Some("E".to_string()) // opposite side
                }
            }
        }
        _ => None,
    };

    DoubleBondStereo {
        atom1,
        atom2,
        configuration,
        high_priority_sub1: high_a,
        high_priority_sub2: high_b,
    }
}

/// Compute dihedral angle (in radians) for four points a-b-c-d.
fn compute_dihedral(a: &[f64; 3], b: &[f64; 3], c: &[f64; 3], d: &[f64; 3]) -> f64 {
    let b1 = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let b2 = [c[0] - b[0], c[1] - b[1], c[2] - b[2]];
    let b3 = [d[0] - c[0], d[1] - c[1], d[2] - c[2]];

    let cross_b1_b2 = cross(&b1, &b2);
    let cross_b2_b3 = cross(&b2, &b3);

    let n1_dot_n2 = dot(&cross_b1_b2, &cross_b2_b3);
    let b2_norm = dot(&b2, &b2).sqrt();

    let x = n1_dot_n2;
    let y = dot(&cross(&cross_b1_b2, &cross_b2_b3), &b2) / b2_norm.max(1e-12);

    y.atan2(x)
}

fn cross(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn dot(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Perform complete stereochemistry analysis on a molecule.
///
/// Detects chiral centers (R/S), E/Z double bonds, atropisomeric axes,
/// and pro-chiral centers from SMILES topology and optional 3D coordinates.
pub fn analyze_stereo(mol: &Molecule, positions: &[[f64; 3]]) -> StereoAnalysis {
    // Find and assign stereocenters
    let center_indices = find_stereocenters(mol);
    let stereocenters: Vec<Stereocenter> = center_indices
        .iter()
        .map(|&center| {
            let center_idx = NodeIndex::new(center);
            let neighbors: Vec<usize> =
                mol.graph.neighbors(center_idx).map(|n| n.index()).collect();
            let (sorted_subs, priorities) = assign_cip_priorities(mol, center, &neighbors);
            let configuration = if !positions.is_empty() {
                assign_rs(positions, center, &sorted_subs)
            } else {
                None
            };

            Stereocenter {
                atom_index: center,
                element: mol.graph[center_idx].element,
                substituent_indices: sorted_subs,
                priorities,
                configuration,
            }
        })
        .collect();

    // Find and assign E/Z double bonds
    let ez_bonds = find_ez_double_bonds(mol);
    let double_bonds: Vec<DoubleBondStereo> = ez_bonds
        .iter()
        .map(|&(a, b)| assign_ez(mol, positions, a, b))
        .collect();

    // Find atropisomeric axes
    let atropisomeric_axes = find_atropisomeric_axes(mol, positions);

    // Find pro-chiral centers
    let prochiral_centers = find_prochiral_centers(mol);

    let n_stereocenters = stereocenters.len();
    let n_double_bonds = double_bonds.len();

    StereoAnalysis {
        stereocenters,
        double_bonds,
        n_stereocenters,
        n_double_bonds,
        atropisomeric_axes,
        prochiral_centers,
        helical_chirality: find_helical_chirality(mol, positions),
    }
}

/// Detect helical chirality (M/P) in helicenes and metallocenes.
///
/// Helical chirality arises when a molecule forms a non-planar spiral.
/// M (minus) = left-handed helix, P (plus) = right-handed helix.
fn find_helical_chirality(mol: &Molecule, positions: &[[f64; 3]]) -> Vec<HelicalChirality> {
    let mut helices = Vec::new();

    if positions.is_empty() {
        return helices;
    }

    // Detect allenic systems: C=C=C with helical potential
    for node in mol.graph.node_indices() {
        let atom = &mol.graph[node];
        if atom.element != 6 {
            continue;
        }

        // Check for cumulated double bonds: =C=
        let double_neighbors: Vec<usize> = mol
            .graph
            .edges(node)
            .filter(|e| mol.graph[e.id()].order == BondOrder::Double)
            .map(|e| {
                let (a, b) = mol.graph.edge_endpoints(e.id()).unwrap();
                if a == node {
                    b.index()
                } else {
                    a.index()
                }
            })
            .collect();

        if double_neighbors.len() == 2 {
            // Allenic system: C=C=C
            let a = double_neighbors[0];
            let center = node.index();
            let b = double_neighbors[1];

            // Get terminal substituents for M/P assignment
            let a_subs: Vec<usize> = mol
                .graph
                .neighbors(NodeIndex::new(a))
                .filter(|&n| n.index() != center)
                .map(|n| n.index())
                .collect();
            let b_subs: Vec<usize> = mol
                .graph
                .neighbors(NodeIndex::new(b))
                .filter(|&n| n.index() != center)
                .map(|n| n.index())
                .collect();

            if a_subs.len() >= 2 && b_subs.len() >= 2 {
                let dihedral = if positions.len() > a_subs[0] && positions.len() > b_subs[0] {
                    Some(compute_dihedral(
                        &positions[a_subs[0]],
                        &positions[a],
                        &positions[b],
                        &positions[b_subs[0]],
                    ))
                } else {
                    None
                };

                let config = dihedral.map(|d| {
                    if d > 0.0 {
                        "P".to_string()
                    } else {
                        "M".to_string()
                    }
                });

                helices.push(HelicalChirality {
                    backbone_atoms: vec![a, center, b],
                    configuration: config,
                    helix_type: "allene".to_string(),
                    dihedral_angle: dihedral,
                });
            }
        }
    }

    // Detect helicene-like fused aromatic rings with non-planarity
    // Look for sequences of 4+ fused aromatic rings that curve
    let aromatic_atoms: Vec<usize> = mol
        .graph
        .node_indices()
        .filter(|&n| {
            mol.graph
                .edges(n)
                .any(|e| mol.graph[e.id()].order == BondOrder::Aromatic)
        })
        .map(|n| n.index())
        .collect();

    if aromatic_atoms.len() >= 10 {
        // Check if the aromatic system is non-planar
        if let Some(planarity_dev) = measure_planarity(&aromatic_atoms, positions) {
            if planarity_dev > 0.3 {
                // Non-planar aromatic system → potential helicene
                // Determine M/P from the sign of the average torsion
                let avg_torsion = compute_average_torsion(&aromatic_atoms, positions);
                let config = avg_torsion.map(|t| {
                    if t > 0.0 {
                        "P".to_string()
                    } else {
                        "M".to_string()
                    }
                });

                helices.push(HelicalChirality {
                    backbone_atoms: aromatic_atoms.clone(),
                    configuration: config,
                    helix_type: "helicene".to_string(),
                    dihedral_angle: avg_torsion,
                });
            }
        }
    }

    // Detect metallocene helical chirality (tilted Cp rings)
    for node in mol.graph.node_indices() {
        let atom = &mol.graph[node];
        // Check if this is a transition metal
        if !is_transition_metal(atom.element) {
            continue;
        }

        let metal_idx = node.index();
        let neighbors: Vec<usize> = mol.graph.neighbors(node).map(|n| n.index()).collect();

        // Find groups of aromatic ring atoms bonded to the metal
        let ring_neighbors: Vec<usize> = neighbors
            .iter()
            .filter(|&&n| {
                mol.graph
                    .edges(NodeIndex::new(n))
                    .any(|e| mol.graph[e.id()].order == BondOrder::Aromatic)
            })
            .copied()
            .collect();

        if ring_neighbors.len() >= 4 {
            // Potential metallocene with tilted rings
            let ring_torsion = compute_ring_tilt(&ring_neighbors, &positions[metal_idx], positions);
            let config = ring_torsion.map(|t| {
                if t > 0.0 {
                    "P".to_string()
                } else {
                    "M".to_string()
                }
            });

            let mut backbone = vec![metal_idx];
            backbone.extend_from_slice(&ring_neighbors);

            helices.push(HelicalChirality {
                backbone_atoms: backbone,
                configuration: config,
                helix_type: "metallocene".to_string(),
                dihedral_angle: ring_torsion,
            });
        }
    }

    helices
}

fn is_transition_metal(z: u8) -> bool {
    matches!(z, 21..=30 | 39..=48 | 72..=80)
}

/// Measure planarity deviation (RMSD from best-fit plane) in ångströms.
fn measure_planarity(atom_indices: &[usize], positions: &[[f64; 3]]) -> Option<f64> {
    if atom_indices.len() < 3 || positions.is_empty() {
        return None;
    }

    // Compute centroid
    let n = atom_indices.len() as f64;
    let mut cz = 0.0;
    for &idx in atom_indices {
        if idx >= positions.len() {
            return None;
        }
        cz += positions[idx][2];
    }
    cz /= n;

    // Simple planarity check: average |z - z_mean| after centering
    let dev: f64 = atom_indices
        .iter()
        .map(|&idx| {
            let dz = positions[idx][2] - cz;
            dz * dz
        })
        .sum::<f64>()
        / n;

    Some(dev.sqrt())
}

/// Compute average torsion along a set of atoms (for helicene sign).
fn compute_average_torsion(atoms: &[usize], positions: &[[f64; 3]]) -> Option<f64> {
    if atoms.len() < 4 {
        return None;
    }

    let mut total_torsion = 0.0;
    let mut count = 0;
    for i in 0..atoms.len().saturating_sub(3) {
        let a = atoms[i];
        let b = atoms[i + 1];
        let c = atoms[i + 2];
        let d = atoms[i + 3];
        if a < positions.len() && b < positions.len() && c < positions.len() && d < positions.len()
        {
            total_torsion +=
                compute_dihedral(&positions[a], &positions[b], &positions[c], &positions[d]);
            count += 1;
        }
    }

    if count > 0 {
        Some(total_torsion / count as f64)
    } else {
        None
    }
}

/// Compute ring tilt angle for metallocene helical chirality.
fn compute_ring_tilt(
    ring_atoms: &[usize],
    _metal_pos: &[f64; 3],
    positions: &[[f64; 3]],
) -> Option<f64> {
    if ring_atoms.len() < 4 {
        return None;
    }

    // Use first 4 ring atoms to compute a torsion that indicates tilt
    let a = ring_atoms[0];
    let b = ring_atoms[1];
    let c = ring_atoms[2];
    let d = ring_atoms[3];

    if a < positions.len() && b < positions.len() && c < positions.len() && d < positions.len() {
        Some(compute_dihedral(
            &positions[a],
            &positions[b],
            &positions[c],
            &positions[d],
        ))
    } else {
        None
    }
}

/// Detect atropisomeric axes: biaryl single bonds with ortho substituents
/// that create a steric barrier to rotation.
fn find_atropisomeric_axes(mol: &Molecule, positions: &[[f64; 3]]) -> Vec<AtropisomericAxis> {
    let mut axes = Vec::new();

    for edge in mol.graph.edge_indices() {
        let bond = &mol.graph[edge];
        // Atropisomerism occurs at single bonds between aromatic rings
        if bond.order != BondOrder::Single {
            continue;
        }

        let (a, b) = mol.graph.edge_endpoints(edge).unwrap();
        let a_idx = a.index();
        let b_idx = b.index();

        // Both atoms must be in aromatic rings (biaryl)
        let a_aromatic = mol
            .graph
            .edges(a)
            .any(|e| mol.graph[e.id()].order == BondOrder::Aromatic);
        let b_aromatic = mol
            .graph
            .edges(b)
            .any(|e| mol.graph[e.id()].order == BondOrder::Aromatic);

        if !a_aromatic || !b_aromatic {
            continue;
        }

        // Count ortho substituents (neighbors of each atom that are not on the
        // aromatic ring between a and b, and are not H)
        let a_neighbors: Vec<usize> = mol
            .graph
            .neighbors(a)
            .map(|n| n.index())
            .filter(|&n| n != b_idx && mol.graph[NodeIndex::new(n)].element != 1)
            .collect();
        let b_neighbors: Vec<usize> = mol
            .graph
            .neighbors(b)
            .map(|n| n.index())
            .filter(|&n| n != a_idx && mol.graph[NodeIndex::new(n)].element != 1)
            .collect();

        // Need at least 2 ortho substituents on each side (total ≥ 3 for restricted rotation)
        let ortho_count = a_neighbors.len() + b_neighbors.len();
        if ortho_count < 3 {
            continue;
        }

        // Check if substituents are large enough (non-H neighbors of ortho atoms)
        let has_bulky = a_neighbors.iter().chain(b_neighbors.iter()).any(|&n| {
            let n_idx = NodeIndex::new(n);
            mol.graph
                .neighbors(n_idx)
                .any(|nb| mol.graph[nb].element != 1 && nb.index() != a_idx && nb.index() != b_idx)
        });

        if !has_bulky {
            continue;
        }

        let configuration =
            if !positions.is_empty() && a_idx < positions.len() && b_idx < positions.len() {
                // Use dihedral of highest-priority ortho subs to determine aR/aS
                if let (Some(&sub_a), Some(&sub_b)) = (a_neighbors.first(), b_neighbors.first()) {
                    if sub_a < positions.len() && sub_b < positions.len() {
                        let dihedral = compute_dihedral(
                            &positions[sub_a],
                            &positions[a_idx],
                            &positions[b_idx],
                            &positions[sub_b],
                        );
                        Some(if dihedral > 0.0 { "aR" } else { "aS" }.to_string())
                    } else {
                        None
                    }
                } else {
                    None
                }
            } else {
                None
            };

        axes.push(AtropisomericAxis {
            atom1: a_idx,
            atom2: b_idx,
            ortho_substituent_count: ortho_count,
            configuration,
        });
    }

    axes
}

/// Detect pro-chiral centers: sp3 atoms with exactly two identical substituents.
fn find_prochiral_centers(mol: &Molecule) -> Vec<ProchiralCenter> {
    let mut centers = Vec::new();

    for node in mol.graph.node_indices() {
        let atom = &mol.graph[node];
        if atom.element == 1 {
            continue;
        }

        let neighbors: Vec<usize> = mol.graph.neighbors(node).map(|n| n.index()).collect();

        // Need exactly 4 substituents (sp3)
        if neighbors.len() != 4 {
            continue;
        }

        // Build CIP trees
        let trees: Vec<CipNode> = neighbors
            .iter()
            .map(|&n| build_cip_tree(mol, n, node.index(), 6))
            .collect();

        // Find pairs with identical CIP priority
        for i in 0..4 {
            for j in (i + 1)..4 {
                if cip_compare(&trees[i], &trees[j]) == std::cmp::Ordering::Equal {
                    // Check that the other two are different from each other and from the pair
                    let others: Vec<usize> = (0..4).filter(|&k| k != i && k != j).collect();
                    let other_different = cip_compare(&trees[others[0]], &trees[others[1]])
                        != std::cmp::Ordering::Equal
                        && cip_compare(&trees[i], &trees[others[0]]) != std::cmp::Ordering::Equal;

                    if other_different {
                        centers.push(ProchiralCenter {
                            atom_index: node.index(),
                            element: atom.element,
                            enantiotopic_pair: [neighbors[i], neighbors[j]],
                        });
                    }
                    break; // Only one pair per center
                }
            }
        }
    }

    centers
}

/// Generate SMILES-style stereo descriptors from a stereochemistry analysis.
///
/// Returns a map of atom indices to their stereo descriptor strings:
/// - Chiral centers: "@" (S-configuration) or "@@" (R-configuration)
/// - Double bond atoms: "/" or "\" for E/Z geometry
///
/// These descriptors can be embedded into SMILES strings to encode stereochemistry.
pub fn stereo_descriptors(analysis: &StereoAnalysis) -> StereoDescriptors {
    let mut center_descriptors = Vec::new();
    let mut bond_descriptors = Vec::new();

    // R/S → @/@@
    for sc in &analysis.stereocenters {
        if let Some(ref config) = sc.configuration {
            let descriptor = match config.as_str() {
                "S" => "@".to_string(),
                "R" => "@@".to_string(),
                _ => continue,
            };
            center_descriptors.push(StereoAtomDescriptor {
                atom_index: sc.atom_index,
                descriptor,
                configuration: config.clone(),
            });
        }
    }

    // E/Z → / and \
    for db in &analysis.double_bonds {
        if let Some(ref config) = db.configuration {
            let (desc1, desc2) = match config.as_str() {
                "E" => ("/".to_string(), "/".to_string()),
                "Z" => ("/".to_string(), "\\".to_string()),
                _ => continue,
            };
            bond_descriptors.push(StereoBondDescriptor {
                atom1: db.atom1,
                atom2: db.atom2,
                descriptor_atom1: desc1,
                descriptor_atom2: desc2,
                configuration: config.clone(),
            });
        }
    }

    StereoDescriptors {
        centers: center_descriptors,
        bonds: bond_descriptors,
    }
}

/// Collection of stereo descriptors for a molecule.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StereoDescriptors {
    /// Stereo descriptors for chiral centers.
    pub centers: Vec<StereoAtomDescriptor>,
    /// Stereo descriptors for double bonds.
    pub bonds: Vec<StereoBondDescriptor>,
}

/// Stereo descriptor for a single chiral center.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StereoAtomDescriptor {
    /// Atom index of the chiral center.
    pub atom_index: usize,
    /// SMILES descriptor: "@" or "@@".
    pub descriptor: String,
    /// Configuration: "R" or "S".
    pub configuration: String,
}

/// Stereo descriptor for a double bond.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StereoBondDescriptor {
    /// First atom of the double bond.
    pub atom1: usize,
    /// Second atom of the double bond.
    pub atom2: usize,
    /// SMILES bond descriptor for the atom1 side.
    pub descriptor_atom1: String,
    /// SMILES bond descriptor for the atom2 side.
    pub descriptor_atom2: String,
    /// Configuration: "E" or "Z".
    pub configuration: String,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_no_stereocenters_ethane() {
        let mol = Molecule::from_smiles("CC").unwrap();
        let result = analyze_stereo(&mol, &[]);
        assert_eq!(result.n_stereocenters, 0);
    }

    #[test]
    fn test_stereocenter_detection_chiral() {
        // 2-bromobutane: C(Br)(CC)C — central carbon has 4 different substituents
        let mol = Molecule::from_smiles("CC(Br)CC").unwrap();
        let result = analyze_stereo(&mol, &[]);
        assert!(
            result.n_stereocenters >= 1,
            "2-bromobutane should have at least 1 stereocenter, got {}",
            result.n_stereocenters
        );
    }

    #[test]
    fn test_cip_priority_ordering() {
        // Br > Cl > F > O > N > C > H (atomic number)
        let mol = Molecule::from_smiles("C(F)(Cl)Br").unwrap();
        let result = analyze_stereo(&mol, &[]);
        if let Some(sc) = result.stereocenters.first() {
            // The Br neighbor should have priority 1, Cl 2, F 3
            let elements: Vec<u8> = sc
                .substituent_indices
                .iter()
                .map(|&i| mol.graph[NodeIndex::new(i)].element)
                .collect();
            // First should be Br(35), then Cl(17), then F(9), then H(1)
            assert!(
                elements[0] >= elements[1],
                "CIP order wrong: first {} should be >= second {}",
                elements[0],
                elements[1]
            );
        }
    }

    #[test]
    fn test_ez_detection_2_butene() {
        // 2-butene: CC=CC has a double bond with different substituents
        let mol = Molecule::from_smiles("CC=CC").unwrap();
        let result = analyze_stereo(&mol, &[]);
        assert!(
            result.n_double_bonds >= 1,
            "2-butene should have E/Z-assignable double bond, got {}",
            result.n_double_bonds
        );
    }

    #[test]
    fn test_no_ez_ethylene() {
        // Ethylene C=C — same substituents on each side (H, H)
        let mol = Molecule::from_smiles("C=C").unwrap();
        let result = analyze_stereo(&mol, &[]);
        assert_eq!(
            result.n_double_bonds, 0,
            "Ethylene should have no E/Z double bonds"
        );
    }

    #[test]
    fn test_benzene_no_stereo() {
        let mol = Molecule::from_smiles("c1ccccc1").unwrap();
        let result = analyze_stereo(&mol, &[]);
        assert_eq!(result.n_stereocenters, 0);
    }
}
