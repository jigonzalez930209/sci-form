//! Pharmacophore fingerprints: encode molecular pharmacophoric features.
//!
//! Detects pharmacophoric feature points (HBD, HBA, aromatic, hydrophobic,
//! positive, negative) and encodes pairwise distance bins into a fingerprint.

use serde::{Deserialize, Serialize};
use std::collections::BTreeSet;

/// Type of pharmacophoric feature.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
pub enum PharmFeatureType {
    /// Hydrogen bond donor (N-H, O-H).
    HBondDonor,
    /// Hydrogen bond acceptor (N, O lone pairs).
    HBondAcceptor,
    /// Aromatic ring centroid.
    Aromatic,
    /// Hydrophobic group.
    Hydrophobic,
    /// Positively ionizable (amines).
    Positive,
    /// Negatively ionizable (carboxylates, phosphates, sulfonates).
    Negative,
}

/// A detected pharmacophoric feature.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PharmFeature {
    /// Feature type.
    pub feature_type: PharmFeatureType,
    /// Atom indices contributing to this feature.
    pub atom_indices: Vec<usize>,
    /// 3D centroid of the feature (if coordinates available).
    pub centroid: Option<[f64; 3]>,
}

/// Pharmacophore fingerprint result.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PharmacophoreFingerprint {
    /// Number of bits in the fingerprint.
    pub n_bits: usize,
    /// Set of on-bits.
    pub on_bits: BTreeSet<usize>,
    /// Detected features.
    pub features: Vec<PharmFeature>,
    /// Bit density (fraction of on bits).
    pub density: f64,
}

/// Detect pharmacophoric features from molecular topology.
///
/// `elements`: atomic numbers.
/// `bonds`: (atom_i, atom_j, bond_order) list.
/// `coords`: flat xyz coordinates (optional, for 3D fingerpoints). Empty = 2D only.
/// `aromatic_atoms`: boolean aromatic flags per atom (or empty).
pub fn detect_features(
    elements: &[u8],
    bonds: &[(usize, usize, u8)],
    coords: &[f64],
    aromatic_atoms: &[bool],
) -> Vec<PharmFeature> {
    let n = elements.len();
    let has_3d = coords.len() >= n * 3;

    // Build adjacency
    let mut adj: Vec<Vec<(usize, u8)>> = vec![vec![]; n];
    for &(i, j, ord) in bonds {
        if i < n && j < n {
            adj[i].push((j, ord));
            adj[j].push((i, ord));
        }
    }

    let mut features = Vec::new();

    for i in 0..n {
        let z = elements[i];
        let centroid = if has_3d {
            Some([coords[i * 3], coords[i * 3 + 1], coords[i * 3 + 2]])
        } else {
            None
        };

        // H-bond donor: N-H, O-H, S-H
        if z == 7 || z == 8 || z == 16 {
            let has_h = adj[i].iter().any(|&(nb, _)| elements[nb] == 1);
            if has_h {
                features.push(PharmFeature {
                    feature_type: PharmFeatureType::HBondDonor,
                    atom_indices: vec![i],
                    centroid,
                });
            }
        }

        // H-bond acceptor: N, O with lone pairs
        if z == 7 || z == 8 {
            features.push(PharmFeature {
                feature_type: PharmFeatureType::HBondAcceptor,
                atom_indices: vec![i],
                centroid,
            });
        }

        // Positive ionizable: primary, secondary, tertiary amines
        if z == 7 {
            let heavy_neighbors = adj[i].iter().filter(|&&(nb, _)| elements[nb] != 1).count();
            if heavy_neighbors <= 3 {
                let has_h = adj[i].iter().any(|&(nb, _)| elements[nb] == 1);
                if has_h {
                    features.push(PharmFeature {
                        feature_type: PharmFeatureType::Positive,
                        atom_indices: vec![i],
                        centroid,
                    });
                }
            }
        }

        // Negative ionizable: carboxylate O, phosphate O, sulfonate O
        if z == 8 {
            let is_terminal = adj[i].len() == 1;
            if is_terminal {
                let &(nb, ord) = &adj[i][0];
                // Check if neighbor C/P/S has double-bonded O
                if (elements[nb] == 6 || elements[nb] == 15 || elements[nb] == 16)
                    && (ord == 1 || ord == 2)
                {
                    let neighbor_has_other_o =
                        adj[nb].iter().any(|&(k, _)| k != i && elements[k] == 8);
                    if neighbor_has_other_o {
                        features.push(PharmFeature {
                            feature_type: PharmFeatureType::Negative,
                            atom_indices: vec![i],
                            centroid,
                        });
                    }
                }
            }
        }

        // Hydrophobic: C, S, halogen not bonded to polar atoms
        if z == 6 || z == 16 || z == 9 || z == 17 || z == 35 || z == 53 {
            let is_polar_neighbor = adj[i]
                .iter()
                .any(|&(nb, _)| elements[nb] == 7 || elements[nb] == 8);
            if !is_polar_neighbor && z == 6 {
                features.push(PharmFeature {
                    feature_type: PharmFeatureType::Hydrophobic,
                    atom_indices: vec![i],
                    centroid,
                });
            }
            if z == 17 || z == 35 || z == 53 {
                features.push(PharmFeature {
                    feature_type: PharmFeatureType::Hydrophobic,
                    atom_indices: vec![i],
                    centroid,
                });
            }
        }
    }

    // Aromatic ring centroids
    if !aromatic_atoms.is_empty() {
        // Detect aromatic ring groups via connected aromatic atoms
        let mut visited = vec![false; n];
        for start in 0..n {
            if !aromatic_atoms.get(start).copied().unwrap_or(false) || visited[start] {
                continue;
            }
            // BFS to find connected aromatic cluster
            let mut ring_atoms = Vec::new();
            let mut queue = std::collections::VecDeque::new();
            queue.push_back(start);
            visited[start] = true;
            while let Some(u) = queue.pop_front() {
                ring_atoms.push(u);
                for &(v, _) in &adj[u] {
                    if !visited[v] && aromatic_atoms.get(v).copied().unwrap_or(false) {
                        visited[v] = true;
                        queue.push_back(v);
                    }
                }
            }
            if ring_atoms.len() >= 5 {
                let centroid = if has_3d {
                    let mut c = [0.0; 3];
                    for &a in &ring_atoms {
                        for k in 0..3 {
                            c[k] += coords[a * 3 + k];
                        }
                    }
                    let n_r = ring_atoms.len() as f64;
                    Some([c[0] / n_r, c[1] / n_r, c[2] / n_r])
                } else {
                    None
                };
                features.push(PharmFeature {
                    feature_type: PharmFeatureType::Aromatic,
                    atom_indices: ring_atoms,
                    centroid,
                });
            }
        }
    }

    features
}

/// Compute pharmacophore fingerprint from detected features.
///
/// Encodes feature-type pairs and (optionally) distance bins into a bit vector.
/// Without 3D coords, encodes only feature-pair counts.
pub fn compute_pharmacophore_fingerprint(
    features: &[PharmFeature],
    n_bits: usize,
) -> PharmacophoreFingerprint {
    let mut on_bits = BTreeSet::new();

    // One-point features: hash feature type × count
    let type_counts: std::collections::HashMap<PharmFeatureType, usize> =
        features
            .iter()
            .fold(std::collections::HashMap::new(), |mut m, f| {
                *m.entry(f.feature_type).or_insert(0) += 1;
                m
            });

    for (&ft, &count) in &type_counts {
        let base = ft as usize * 31;
        on_bits.insert((base + count.min(10)) % n_bits);
    }

    // Two-point features: hash pairs of feature types
    for i in 0..features.len() {
        for j in (i + 1)..features.len() {
            let ft1 = features[i].feature_type as usize;
            let ft2 = features[j].feature_type as usize;
            let (lo, hi) = if ft1 <= ft2 { (ft1, ft2) } else { (ft2, ft1) };

            // If 3D, add distance bin encoding
            if let (Some(c1), Some(c2)) = (&features[i].centroid, &features[j].centroid) {
                let dist =
                    ((c1[0] - c2[0]).powi(2) + (c1[1] - c2[1]).powi(2) + (c1[2] - c2[2]).powi(2))
                        .sqrt();
                // Bin distance: 0-2Å, 2-4Å, 4-6Å, 6-8Å, 8-10Å, 10+Å
                let bin = ((dist / 2.0) as usize).min(5);
                let hash = (lo * 7 + hi * 13 + bin * 37 + 1000) % n_bits;
                on_bits.insert(hash);
            } else {
                let hash = (lo * 7 + hi * 13 + 500) % n_bits;
                on_bits.insert(hash);
            }
        }
    }

    let density = on_bits.len() as f64 / n_bits as f64;

    PharmacophoreFingerprint {
        n_bits,
        on_bits,
        features: features.to_vec(),
        density,
    }
}

/// Compute Tanimoto similarity between two pharmacophore fingerprints.
pub fn pharmacophore_tanimoto(
    fp1: &PharmacophoreFingerprint,
    fp2: &PharmacophoreFingerprint,
) -> f64 {
    let intersection = fp1.on_bits.intersection(&fp2.on_bits).count();
    let union = fp1.on_bits.union(&fp2.on_bits).count();
    if union == 0 {
        0.0
    } else {
        intersection as f64 / union as f64
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_detect_features_ethanol() {
        // Ethanol: CCO → has HBD (O-H), HBA (O), hydrophobic (C)
        let elements = [6u8, 6, 8, 1, 1, 1, 1, 1, 1];
        let bonds = [
            (0, 1, 1u8),
            (1, 2, 1),
            (0, 3, 1),
            (0, 4, 1),
            (0, 5, 1),
            (1, 6, 1),
            (1, 7, 1),
            (2, 8, 1),
        ];
        let features = detect_features(&elements, &bonds, &[], &[]);
        let has_hbd = features
            .iter()
            .any(|f| f.feature_type == PharmFeatureType::HBondDonor);
        let has_hba = features
            .iter()
            .any(|f| f.feature_type == PharmFeatureType::HBondAcceptor);
        assert!(has_hbd, "Ethanol should have an H-bond donor");
        assert!(has_hba, "Ethanol should have an H-bond acceptor");
    }

    #[test]
    fn test_pharmacophore_fingerprint() {
        let elements = [6u8, 6, 8, 1, 1, 1, 1, 1, 1];
        let bonds = [
            (0, 1, 1u8),
            (1, 2, 1),
            (0, 3, 1),
            (0, 4, 1),
            (0, 5, 1),
            (1, 6, 1),
            (1, 7, 1),
            (2, 8, 1),
        ];
        let features = detect_features(&elements, &bonds, &[], &[]);
        let fp = compute_pharmacophore_fingerprint(&features, 1024);
        assert!(!fp.on_bits.is_empty(), "Fingerprint should have on-bits");
        assert!(fp.density > 0.0 && fp.density < 1.0);
    }

    #[test]
    fn test_pharmacophore_tanimoto_self() {
        let elements = [6u8, 8, 1, 1, 1];
        let bonds = [(0, 1, 1u8), (0, 2, 1), (0, 3, 1), (1, 4, 1)];
        let features = detect_features(&elements, &bonds, &[], &[]);
        let fp = compute_pharmacophore_fingerprint(&features, 1024);
        let sim = pharmacophore_tanimoto(&fp, &fp);
        assert!((sim - 1.0).abs() < 1e-10, "Self-similarity should be 1.0");
    }
}
