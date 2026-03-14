// Scientific/numerical code patterns that are idiomatic in this domain
#![allow(clippy::too_many_arguments)]
#![allow(clippy::needless_range_loop)]

pub mod conformer;
pub mod distgeom;
pub mod etkdg;
pub mod forcefield;
pub mod graph;
pub mod optimization;
pub mod smarts;
pub mod smiles;

use serde::{Deserialize, Serialize};

// ─── Public API Types ────────────────────────────────────────────────────────

/// Result of a 3D conformer generation for a single molecule.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConformerResult {
    /// Input SMILES string.
    pub smiles: String,
    /// Number of atoms in the molecule.
    pub num_atoms: usize,
    /// Flat xyz coordinates: [x0, y0, z0, x1, y1, z1, ...].
    /// Empty on failure.
    pub coords: Vec<f64>,
    /// Atom elements (atomic numbers) in the same order as coords.
    pub elements: Vec<u8>,
    /// Bond list as (start_atom, end_atom, order_string).
    pub bonds: Vec<(usize, usize, String)>,
    /// Error message if generation failed.
    pub error: Option<String>,
    /// Generation time in milliseconds.
    pub time_ms: f64,
}

/// Configuration for conformer generation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConformerConfig {
    /// RNG seed (same seed = reproducible output).
    pub seed: u64,
    /// Number of threads for batch processing (0 = auto-detect).
    pub num_threads: usize,
}

impl Default for ConformerConfig {
    fn default() -> Self {
        Self {
            seed: 42,
            num_threads: 0,
        }
    }
}

// ─── Public API Functions ────────────────────────────────────────────────────

/// Library version string.
pub fn version() -> String {
    format!("sci-form {}", env!("CARGO_PKG_VERSION"))
}

/// Generate a 3D conformer from a SMILES string.
pub fn embed(smiles: &str, seed: u64) -> ConformerResult {
    #[cfg(not(target_arch = "wasm32"))]
    let start = std::time::Instant::now();

    let mol = match graph::Molecule::from_smiles(smiles) {
        Ok(m) => m,
        Err(e) => {
            return ConformerResult {
                smiles: smiles.to_string(),
                num_atoms: 0,
                coords: vec![],
                elements: vec![],
                bonds: vec![],
                error: Some(e),
                #[cfg(not(target_arch = "wasm32"))]
                time_ms: start.elapsed().as_secs_f64() * 1000.0,
                #[cfg(target_arch = "wasm32")]
                time_ms: 0.0,
            };
        }
    };

    let n = mol.graph.node_count();
    let elements: Vec<u8> = (0..n)
        .map(|i| mol.graph[petgraph::graph::NodeIndex::new(i)].element)
        .collect();
    let bonds: Vec<(usize, usize, String)> = mol
        .graph
        .edge_indices()
        .map(|e| {
            let (a, b) = mol.graph.edge_endpoints(e).unwrap();
            let order = match mol.graph[e].order {
                graph::BondOrder::Single => "SINGLE",
                graph::BondOrder::Double => "DOUBLE",
                graph::BondOrder::Triple => "TRIPLE",
                graph::BondOrder::Aromatic => "AROMATIC",
                graph::BondOrder::Unknown => "UNKNOWN",
            };
            (a.index(), b.index(), order.to_string())
        })
        .collect();

    match conformer::generate_3d_conformer(&mol, seed) {
        Ok(coords) => {
            let mut flat = Vec::with_capacity(n * 3);
            for i in 0..n {
                flat.push(coords[(i, 0)] as f64);
                flat.push(coords[(i, 1)] as f64);
                flat.push(coords[(i, 2)] as f64);
            }
            ConformerResult {
                smiles: smiles.to_string(),
                num_atoms: n,
                coords: flat,
                elements,
                bonds,
                error: None,
                #[cfg(not(target_arch = "wasm32"))]
                time_ms: start.elapsed().as_secs_f64() * 1000.0,
                #[cfg(target_arch = "wasm32")]
                time_ms: 0.0,
            }
        }
        Err(e) => ConformerResult {
            smiles: smiles.to_string(),
            num_atoms: n,
            coords: vec![],
            elements,
            bonds,
            error: Some(e),
            #[cfg(not(target_arch = "wasm32"))]
            time_ms: start.elapsed().as_secs_f64() * 1000.0,
            #[cfg(target_arch = "wasm32")]
            time_ms: 0.0,
        },
    }
}

/// Batch-embed multiple SMILES in parallel.
///
/// Uses rayon thread pool for parallel processing.
/// `config.num_threads = 0` means auto-detect CPU count.
#[cfg(feature = "parallel")]
pub fn embed_batch(smiles_list: &[&str], config: &ConformerConfig) -> Vec<ConformerResult> {
    use rayon::prelude::*;

    if config.num_threads > 0 {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(config.num_threads)
            .build()
            .unwrap();
        pool.install(|| {
            smiles_list
                .par_iter()
                .map(|smi| embed(smi, config.seed))
                .collect()
        })
    } else {
        smiles_list
            .par_iter()
            .map(|smi| embed(smi, config.seed))
            .collect()
    }
}

/// Batch-embed multiple SMILES sequentially (no rayon dependency).
#[cfg(not(feature = "parallel"))]
pub fn embed_batch(smiles_list: &[&str], config: &ConformerConfig) -> Vec<ConformerResult> {
    smiles_list
        .iter()
        .map(|smi| embed(smi, config.seed))
        .collect()
}

/// Parse a SMILES string and return molecular structure (no 3D generation).
pub fn parse(smiles: &str) -> Result<graph::Molecule, String> {
    graph::Molecule::from_smiles(smiles)
}
