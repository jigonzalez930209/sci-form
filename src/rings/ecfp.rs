//! Extended-Connectivity Fingerprints (ECFP / Morgan fingerprints).
//!
//! Implements Rogers & Hahn 2010 circular fingerprints:
//! 1. Initialize atom invariants (element, degree, H-count, charge, ring)
//! 2. Iteratively aggregate neighbor identifiers up to a given radius
//! 3. Fold into a fixed-length bit vector

use crate::graph::{BondOrder, Molecule};
use petgraph::graph::NodeIndex;
use petgraph::visit::EdgeRef;
use serde::{Deserialize, Serialize};
use std::collections::BTreeSet;

/// An ECFP fingerprint represented as a set of on-bits in a folded bit vector.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ECFPFingerprint {
    /// Number of bits in the fingerprint.
    pub n_bits: usize,
    /// Set of on-bit positions.
    pub on_bits: BTreeSet<usize>,
    /// Radius used for generation.
    pub radius: usize,
    /// Raw feature identifiers (before folding).
    pub raw_features: Vec<u64>,
}

impl ECFPFingerprint {
    /// Get the density (fraction of on-bits).
    pub fn density(&self) -> f64 {
        self.on_bits.len() as f64 / self.n_bits as f64
    }

    /// Convert to a dense bit vector.
    pub fn to_bitvec(&self) -> Vec<bool> {
        let mut bv = vec![false; self.n_bits];
        for &bit in &self.on_bits {
            bv[bit] = true;
        }
        bv
    }
}

/// FNV-1a hash for deterministic fingerprint generation.
fn fnv1a_hash(data: &[u64]) -> u64 {
    let mut hash: u64 = 0xcbf29ce484222325; // FNV offset basis
    for &val in data {
        let bytes = val.to_le_bytes();
        for &byte in &bytes {
            hash ^= byte as u64;
            hash = hash.wrapping_mul(0x100000001b3); // FNV prime
        }
    }
    hash
}

/// Compute initial atom invariant for ECFP generation.
///
/// Encodes: atomic number, degree, H-count, formal charge, ring membership.
fn atom_invariant(mol: &Molecule, idx: NodeIndex) -> u64 {
    let atom = &mol.graph[idx];
    let element = atom.element as u64;
    let degree = mol.graph.edges(idx).count() as u64;
    let h_count = mol
        .graph
        .neighbors(idx)
        .filter(|n| mol.graph[*n].element == 1)
        .count() as u64;
    let formal_charge = (atom.formal_charge as i64 + 10) as u64; // offset to avoid negative
    let is_aromatic = mol
        .graph
        .edges(idx)
        .any(|e| matches!(e.weight().order, BondOrder::Aromatic));
    let ring_flag = if crate::graph::atom_in_ring(mol, idx) {
        1u64
    } else {
        0u64
    };
    let aromatic_flag = if is_aromatic { 1u64 } else { 0u64 };

    // Combine into a single hash
    fnv1a_hash(&[
        element,
        degree,
        h_count,
        formal_charge,
        ring_flag,
        aromatic_flag,
    ])
}

/// Compute ECFP fingerprint for a molecule.
///
/// # Arguments
/// - `mol`: parsed molecular graph
/// - `radius`: maximum radius (ECFP4 = radius 2, ECFP6 = radius 3)
/// - `n_bits`: bit vector length (typically 1024 or 2048)
///
/// # Returns
/// `ECFPFingerprint` with on-bits and raw features.
pub fn compute_ecfp(mol: &Molecule, radius: usize, n_bits: usize) -> ECFPFingerprint {
    let n = mol.graph.node_count();
    let n_bits = n_bits.max(64);

    // 1. Initialize invariants
    let mut current_ids: Vec<u64> = (0..n)
        .map(|i| atom_invariant(mol, NodeIndex::new(i)))
        .collect();

    let mut all_features: Vec<u64> = current_ids.clone();

    // 2. Iterative neighborhood aggregation
    for _r in 0..radius {
        let mut next_ids = Vec::with_capacity(n);

        for i in 0..n {
            let node = NodeIndex::new(i);

            // Collect neighbor identifiers with bond info
            let mut neighbor_data: Vec<u64> = Vec::new();
            for edge in mol.graph.edges(node) {
                let neighbor = if edge.source() == node {
                    edge.target()
                } else {
                    edge.source()
                };
                let bond_type: u64 = match edge.weight().order {
                    BondOrder::Single => 1,
                    BondOrder::Double => 2,
                    BondOrder::Triple => 3,
                    BondOrder::Aromatic => 4,
                    BondOrder::Unknown => 0,
                };
                // Hash bond_type with neighbor's current ID
                neighbor_data.push(fnv1a_hash(&[bond_type, current_ids[neighbor.index()]]));
            }

            // Sort for order independence
            neighbor_data.sort();

            // Hash: [current_id, sorted_neighbor_hashes...]
            let mut hash_input = vec![current_ids[i]];
            hash_input.extend_from_slice(&neighbor_data);
            let new_id = fnv1a_hash(&hash_input);

            next_ids.push(new_id);
            all_features.push(new_id);
        }

        current_ids = next_ids;
    }

    // 3. Fold into bit vector
    let mut on_bits = BTreeSet::new();
    for &feature in &all_features {
        let bit = (feature % n_bits as u64) as usize;
        on_bits.insert(bit);
    }

    ECFPFingerprint {
        n_bits,
        on_bits,
        radius,
        raw_features: all_features,
    }
}

/// Compute the Tanimoto similarity between two ECFP fingerprints.
///
/// $$T(A,B) = \frac{|A \cap B|}{|A \cup B|}$$
///
/// Returns a value in [0.0, 1.0] where 1.0 = identical fingerprints.
pub fn compute_tanimoto(fp1: &ECFPFingerprint, fp2: &ECFPFingerprint) -> f64 {
    let intersection = fp1.on_bits.intersection(&fp2.on_bits).count();
    let union = fp1.on_bits.union(&fp2.on_bits).count();

    if union == 0 {
        return 1.0; // both empty
    }

    intersection as f64 / union as f64
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ecfp_benzene() {
        let mol = Molecule::from_smiles("c1ccccc1").unwrap();
        let fp = compute_ecfp(&mol, 2, 1024);

        assert_eq!(fp.n_bits, 1024);
        assert_eq!(fp.radius, 2);
        assert!(!fp.on_bits.is_empty(), "Benzene should have on-bits");
        assert!(fp.density() > 0.0 && fp.density() < 1.0);
    }

    #[test]
    fn test_self_similarity() {
        let mol = Molecule::from_smiles("CCO").unwrap();
        let fp = compute_ecfp(&mol, 2, 1024);
        let tanimoto = compute_tanimoto(&fp, &fp);
        assert!(
            (tanimoto - 1.0).abs() < 1e-10,
            "Self-similarity should be 1.0, got {}",
            tanimoto
        );
    }

    #[test]
    fn test_tanimoto_similar_molecules() {
        let benzene = Molecule::from_smiles("c1ccccc1").unwrap();
        let toluene = Molecule::from_smiles("Cc1ccccc1").unwrap();
        let fp1 = compute_ecfp(&benzene, 2, 2048);
        let fp2 = compute_ecfp(&toluene, 2, 2048);
        let tanimoto = compute_tanimoto(&fp1, &fp2);

        assert!(
            tanimoto > 0.3,
            "Benzene-toluene similarity should be >0.3, got {}",
            tanimoto
        );
    }

    #[test]
    fn test_tanimoto_dissimilar_molecules() {
        let benzene = Molecule::from_smiles("c1ccccc1").unwrap();
        let hexane = Molecule::from_smiles("CCCCCC").unwrap();
        let fp1 = compute_ecfp(&benzene, 2, 2048);
        let fp2 = compute_ecfp(&hexane, 2, 2048);
        let tanimoto = compute_tanimoto(&fp1, &fp2);

        // Dissimilar molecules should have lower Tanimoto
        assert!(
            tanimoto < 0.5,
            "Benzene-hexane similarity should be <0.5, got {}",
            tanimoto
        );
    }

    #[test]
    fn test_deterministic_fingerprint() {
        let mol = Molecule::from_smiles("CC(=O)O").unwrap();
        let fp1 = compute_ecfp(&mol, 2, 1024);
        let fp2 = compute_ecfp(&mol, 2, 1024);
        assert_eq!(fp1.on_bits, fp2.on_bits, "ECFP should be deterministic");
        assert_eq!(fp1.raw_features, fp2.raw_features);
    }

    #[test]
    fn test_different_radii() {
        let mol = Molecule::from_smiles("c1ccc(O)cc1").unwrap();
        let fp2 = compute_ecfp(&mol, 1, 1024);
        let fp4 = compute_ecfp(&mol, 2, 1024);
        let fp6 = compute_ecfp(&mol, 3, 1024);

        // Higher radius should produce more features
        assert!(
            fp4.raw_features.len() >= fp2.raw_features.len(),
            "ECFP4 should have >= features than ECFP2"
        );
        assert!(
            fp6.raw_features.len() >= fp4.raw_features.len(),
            "ECFP6 should have >= features than ECFP4"
        );
    }

    #[test]
    fn test_fnv1a_deterministic() {
        let a = fnv1a_hash(&[1, 2, 3]);
        let b = fnv1a_hash(&[1, 2, 3]);
        assert_eq!(a, b);

        let c = fnv1a_hash(&[3, 2, 1]);
        assert_ne!(a, c, "Different inputs should produce different hashes");
    }
}
