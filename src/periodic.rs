//! Periodic molecular graph support for 3D crystalline systems.
//!
//! Extends the molecular graph with periodic boundary conditions (PBC),
//! enabling representation of infinite crystalline structures and
//! metallocene/hapticity detection for organometallic complexes.

use crate::materials::UnitCell;
use serde::{Deserialize, Serialize};

/// A periodic molecular graph: molecule + unit cell + image bonds.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PeriodicMolecule {
    /// Atomic numbers in the asymmetric unit.
    pub elements: Vec<u8>,
    /// Fractional coordinates in the unit cell.
    pub frac_coords: Vec<[f64; 3]>,
    /// Cartesian coordinates (Å).
    pub cart_coords: Vec<[f64; 3]>,
    /// Bonds within the unit cell: (atom_i, atom_j, order, image_vector).
    pub bonds: Vec<PeriodicBond>,
    /// Unit cell definition.
    pub cell: UnitCell,
    /// Number of atoms in the asymmetric unit.
    pub n_atoms: usize,
}

/// A bond that may cross periodic boundaries.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PeriodicBond {
    /// Source atom index (in asymmetric unit).
    pub atom_i: usize,
    /// Target atom index (in asymmetric unit).
    pub atom_j: usize,
    /// Bond order.
    pub order: String,
    /// Image vector [na, nb, nc] — [0,0,0] for bonds within the cell.
    pub image: [i32; 3],
    /// Bond length in Å.
    pub distance: f64,
}

/// Hapticity descriptor for metal-ring interactions.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HapticInteraction {
    /// Metal atom index.
    pub metal_index: usize,
    /// Metal element (atomic number).
    pub metal_element: u8,
    /// Ring atom indices bound to the metal.
    pub ring_atoms: Vec<usize>,
    /// Hapticity (η number): number of contiguous ring atoms coordinated.
    pub hapticity: usize,
    /// Centroid of the ring atoms (Cartesian Å).
    pub centroid: [f64; 3],
    /// Metal-centroid distance (Å).
    pub metal_centroid_distance: f64,
}

/// Result of metallocene/hapticity detection.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HapticityAnalysis {
    /// Detected haptic interactions.
    pub interactions: Vec<HapticInteraction>,
    /// Whether the molecule is a metallocene (two η5-Cp rings).
    pub is_metallocene: bool,
    /// Number of haptic interactions found.
    pub n_interactions: usize,
}

/// Common transition metals for organometallic detection.
fn is_transition_metal(z: u8) -> bool {
    matches!(z, 21..=30 | 39..=48 | 57..=80)
}

/// Build a periodic molecular graph from elements, fractional coordinates, and a unit cell.
///
/// Automatically detects bonds using distance criteria with periodic images.
pub fn build_periodic_molecule(
    elements: &[u8],
    frac_coords: &[[f64; 3]],
    cell: &UnitCell,
    bond_tolerance: Option<f64>,
) -> PeriodicMolecule {
    let tol = bond_tolerance.unwrap_or(0.3);
    let n_atoms = elements.len();

    let cart_coords: Vec<[f64; 3]> = frac_coords.iter().map(|f| cell.frac_to_cart(*f)).collect();

    let mut bonds = Vec::new();

    // Check all atom pairs including periodic images (-1, 0, +1 in each direction)
    for i in 0..n_atoms {
        for j in i..n_atoms {
            for na in -1i32..=1 {
                for nb in -1i32..=1 {
                    for nc in -1i32..=1 {
                        if i == j && na == 0 && nb == 0 && nc == 0 {
                            continue;
                        }

                        let fj_image = [
                            frac_coords[j][0] + na as f64,
                            frac_coords[j][1] + nb as f64,
                            frac_coords[j][2] + nc as f64,
                        ];
                        let cj_image = cell.frac_to_cart(fj_image);
                        let dist = distance_3d(&cart_coords[i], &cj_image);

                        let r_sum = crate::graph::get_covalent_radius(elements[i])
                            + crate::graph::get_covalent_radius(elements[j]);

                        // Tighter tolerance for TM-TM pairs to avoid spurious bonds
                        let pair_tol = if is_transition_metal(elements[i])
                            && is_transition_metal(elements[j])
                        {
                            tol * 0.3
                        } else {
                            tol
                        };

                        if dist < r_sum + pair_tol && dist > 0.4 {
                            let order = if dist < r_sum * 0.85 {
                                "DOUBLE"
                            } else if dist < r_sum * 0.75 {
                                "TRIPLE"
                            } else {
                                "SINGLE"
                            };

                            bonds.push(PeriodicBond {
                                atom_i: i,
                                atom_j: j,
                                order: order.to_string(),
                                image: [na, nb, nc],
                                distance: dist,
                            });
                        }
                    }
                }
            }
        }
    }

    PeriodicMolecule {
        elements: elements.to_vec(),
        frac_coords: frac_coords.to_vec(),
        cart_coords,
        bonds,
        cell: cell.clone(),
        n_atoms,
    }
}

/// Detect metallocene and haptic (η) metal-ring interactions.
///
/// Scans for transition metals bonded to contiguous ring atoms
/// and determines the hapticity (η1–η8).
pub fn detect_hapticity(
    elements: &[u8],
    coords: &[[f64; 3]],
    bonds: &[(usize, usize)],
) -> HapticityAnalysis {
    let n_atoms = elements.len();

    // Find ring atoms using simple cycle detection
    let ring_atoms = find_ring_atoms(n_atoms, bonds);

    // Find transition metals
    let metals: Vec<usize> = (0..n_atoms)
        .filter(|&i| is_transition_metal(elements[i]))
        .collect();

    let mut interactions = Vec::new();

    for &metal_idx in &metals {
        // Find all ring atoms within bonding distance of this metal
        let metal_pos = coords[metal_idx];
        let mut coordinated_ring_atoms: Vec<usize> = Vec::new();

        for &ring_atom in &ring_atoms {
            // Skip other metals — they aren't part of haptic rings
            if ring_atom == metal_idx || is_transition_metal(elements[ring_atom]) {
                continue;
            }
            let dist = distance_3d(&metal_pos, &coords[ring_atom]);
            let max_dist = metal_ring_cutoff(elements[metal_idx], elements[ring_atom]);
            if dist < max_dist {
                coordinated_ring_atoms.push(ring_atom);
            }
        }

        if coordinated_ring_atoms.is_empty() {
            continue;
        }

        // Group contiguous ring atoms into haptic sets
        let haptic_groups = group_contiguous_ring_atoms(&coordinated_ring_atoms, bonds);

        for group in haptic_groups {
            if group.is_empty() {
                continue;
            }

            let centroid = compute_centroid(&group, coords);
            let metal_centroid_dist = distance_3d(&metal_pos, &centroid);

            interactions.push(HapticInteraction {
                metal_index: metal_idx,
                metal_element: elements[metal_idx],
                ring_atoms: group.clone(),
                hapticity: group.len(),
                centroid,
                metal_centroid_distance: metal_centroid_dist,
            });
        }
    }

    // A metallocene has exactly two η5-Cp rings on one metal
    let is_metallocene = metals.len() == 1
        && interactions.len() == 2
        && interactions.iter().all(|h| h.hapticity == 5);

    let n_interactions = interactions.len();

    HapticityAnalysis {
        interactions,
        is_metallocene,
        n_interactions,
    }
}

/// Maximum metal-ring atom distance for haptic coordination (Å).
fn metal_ring_cutoff(metal_z: u8, ring_z: u8) -> f64 {
    let r_metal = crate::graph::get_covalent_radius(metal_z);
    let r_ring = crate::graph::get_covalent_radius(ring_z);
    // Haptic bonds are typically 0.3–0.5 Å longer than σ bonds
    r_metal + r_ring + 0.8
}

/// Find atoms that belong to rings via simple DFS cycle detection.
fn find_ring_atoms(n_atoms: usize, bonds: &[(usize, usize)]) -> Vec<usize> {
    let mut adj = vec![vec![]; n_atoms];
    for &(a, b) in bonds {
        adj[a].push(b);
        adj[b].push(a);
    }

    let mut in_ring = vec![false; n_atoms];
    let mut visited = vec![false; n_atoms];
    let mut parent = vec![usize::MAX; n_atoms];

    for start in 0..n_atoms {
        if visited[start] {
            continue;
        }
        let mut stack = vec![(start, 0usize)];
        while let Some((node, nb_idx)) = stack.last_mut() {
            let node = *node;
            if !visited[node] {
                visited[node] = true;
            }
            if *nb_idx < adj[node].len() {
                let neighbor = adj[node][*nb_idx];
                *nb_idx += 1;
                if !visited[neighbor] {
                    parent[neighbor] = node;
                    stack.push((neighbor, 0));
                } else if neighbor != parent[node] {
                    // Found a cycle — mark atoms on the path
                    in_ring[neighbor] = true;
                    let mut cur = node;
                    while cur != neighbor && cur != usize::MAX {
                        in_ring[cur] = true;
                        cur = parent[cur];
                    }
                }
            } else {
                stack.pop();
            }
        }
    }

    (0..n_atoms).filter(|&i| in_ring[i]).collect()
}

/// Group ring atoms into contiguous sets based on bonding connectivity.
fn group_contiguous_ring_atoms(ring_atoms: &[usize], bonds: &[(usize, usize)]) -> Vec<Vec<usize>> {
    if ring_atoms.is_empty() {
        return vec![];
    }

    let atom_set: std::collections::HashSet<usize> = ring_atoms.iter().copied().collect();
    let mut visited = std::collections::HashSet::new();
    let mut groups = Vec::new();

    for &start in ring_atoms {
        if visited.contains(&start) {
            continue;
        }
        let mut group = Vec::new();
        let mut queue = std::collections::VecDeque::new();
        queue.push_back(start);
        visited.insert(start);

        while let Some(node) = queue.pop_front() {
            group.push(node);
            for &(a, b) in bonds {
                let neighbor = if a == node {
                    b
                } else if b == node {
                    a
                } else {
                    continue;
                };
                if atom_set.contains(&neighbor) && !visited.contains(&neighbor) {
                    visited.insert(neighbor);
                    queue.push_back(neighbor);
                }
            }
        }

        groups.push(group);
    }

    groups
}

fn compute_centroid(atoms: &[usize], coords: &[[f64; 3]]) -> [f64; 3] {
    let n = atoms.len() as f64;
    let mut c = [0.0, 0.0, 0.0];
    for &i in atoms {
        c[0] += coords[i][0];
        c[1] += coords[i][1];
        c[2] += coords[i][2];
    }
    [c[0] / n, c[1] / n, c[2] / n]
}

fn distance_3d(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::materials::UnitCell;

    #[test]
    fn test_periodic_molecule_cubic() {
        let cell = UnitCell::cubic(5.0);
        let elements = vec![26, 8]; // Fe, O
        let frac = vec![[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]];
        let pm = build_periodic_molecule(&elements, &frac, &cell, None);
        assert_eq!(pm.n_atoms, 2);
    }

    #[test]
    fn test_hapticity_detection() {
        // Ferrocene-like: Fe at origin, two Cp rings
        let mut elements = vec![26u8]; // Fe
        let mut coords = vec![[0.0, 0.0, 0.0f64]];
        let mut bonds = Vec::new();

        // Top Cp ring (5 carbons)
        for i in 0..5 {
            let angle = 2.0 * std::f64::consts::PI * i as f64 / 5.0;
            let x = 1.2 * angle.cos();
            let y = 1.2 * angle.sin();
            elements.push(6); // C
            coords.push([x, y, 1.65]);
            // Bond C to Fe
            bonds.push((0, i + 1));
        }
        // Ring bonds for top Cp
        for i in 0..5 {
            bonds.push((i + 1, (i + 1) % 5 + 1));
        }

        // Bottom Cp ring
        for i in 0..5 {
            let angle = 2.0 * std::f64::consts::PI * i as f64 / 5.0;
            let x = 1.2 * angle.cos();
            let y = 1.2 * angle.sin();
            elements.push(6);
            coords.push([x, y, -1.65]);
            bonds.push((0, i + 6));
        }
        for i in 0..5 {
            bonds.push((i + 6, (i + 1) % 5 + 6));
        }

        let result = detect_hapticity(&elements, &coords, &bonds);
        assert!(result.n_interactions >= 2);
        assert!(result.interactions.iter().any(|h| h.hapticity == 5));
    }
}
