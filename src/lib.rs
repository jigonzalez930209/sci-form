// Scientific/numerical code patterns that are idiomatic in this domain
#![allow(clippy::too_many_arguments)]
#![allow(clippy::needless_range_loop)]

pub mod alignment;
pub mod charges;
pub mod conformer;
pub mod dipole;
pub mod distgeom;
pub mod dos;
pub mod eht;
pub mod esp;
pub mod etkdg;
pub mod forcefield;
pub mod graph;
pub mod materials;
pub mod optimization;
pub mod population;
pub mod smarts;
pub mod smiles;
pub mod surface;
pub mod transport;

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

/// Compute Gasteiger-Marsili partial charges from a SMILES string.
///
/// Parses the molecule, extracts bonds and elements, then runs 6 iterations
/// of electronegativity equalization.
pub fn compute_charges(smiles: &str) -> Result<charges::gasteiger::ChargeResult, String> {
    let mol = graph::Molecule::from_smiles(smiles)?;
    let n = mol.graph.node_count();
    let elements: Vec<u8> = (0..n)
        .map(|i| mol.graph[petgraph::graph::NodeIndex::new(i)].element)
        .collect();
    let formal_charges: Vec<i8> = (0..n)
        .map(|i| mol.graph[petgraph::graph::NodeIndex::new(i)].formal_charge)
        .collect();
    let bonds: Vec<(usize, usize)> = mol
        .graph
        .edge_indices()
        .map(|e| {
            let (a, b) = mol.graph.edge_endpoints(e).unwrap();
            (a.index(), b.index())
        })
        .collect();
    charges::gasteiger::gasteiger_marsili_charges(&elements, &bonds, &formal_charges, 6)
}

/// Compute solvent-accessible surface area from SMILES + 3D coordinates.
///
/// The `coords_flat` parameter is a flat [x0,y0,z0, x1,y1,z1,...] array.
pub fn compute_sasa(
    elements: &[u8],
    coords_flat: &[f64],
    probe_radius: Option<f64>,
) -> Result<surface::sasa::SasaResult, String> {
    if coords_flat.len() != elements.len() * 3 {
        return Err(format!(
            "coords length {} != 3 * elements {}",
            coords_flat.len(),
            elements.len()
        ));
    }
    let positions: Vec<[f64; 3]> = coords_flat
        .chunks_exact(3)
        .map(|c| [c[0], c[1], c[2]])
        .collect();

    #[cfg(feature = "parallel")]
    {
        Ok(surface::sasa::compute_sasa_parallel(
            elements,
            &positions,
            probe_radius,
            None,
        ))
    }

    #[cfg(not(feature = "parallel"))]
    {
        Ok(surface::sasa::compute_sasa(
            elements,
            &positions,
            probe_radius,
            None,
        ))
    }
}

/// Compute Mulliken & Löwdin population analysis from atomic elements and positions.
pub fn compute_population(
    elements: &[u8],
    positions: &[[f64; 3]],
) -> Result<population::PopulationResult, String> {
    let eht_result = eht::solve_eht(elements, positions, None)?;
    Ok(population::compute_population(
        elements,
        positions,
        &eht_result.coefficients,
        eht_result.n_electrons,
    ))
}

/// Compute molecular dipole moment (Debye) from atomic elements and positions.
pub fn compute_dipole(
    elements: &[u8],
    positions: &[[f64; 3]],
) -> Result<dipole::DipoleResult, String> {
    let eht_result = eht::solve_eht(elements, positions, None)?;
    Ok(dipole::compute_dipole_from_eht(
        elements,
        positions,
        &eht_result.coefficients,
        eht_result.n_electrons,
    ))
}

/// Compute ESP grid from atomic elements, positions and Mulliken charges.
pub fn compute_esp(
    elements: &[u8],
    positions: &[[f64; 3]],
    spacing: f64,
    padding: f64,
) -> Result<esp::EspGrid, String> {
    let pop = compute_population(elements, positions)?;
    #[cfg(feature = "parallel")]
    {
        Ok(esp::compute_esp_grid_parallel(
            elements,
            positions,
            &pop.mulliken_charges,
            spacing,
            padding,
        ))
    }

    #[cfg(not(feature = "parallel"))]
    {
        Ok(esp::compute_esp_grid(
            elements,
            positions,
            &pop.mulliken_charges,
            spacing,
            padding,
        ))
    }
}

/// Compute DOS/PDOS from EHT orbital energies.
pub fn compute_dos(
    elements: &[u8],
    positions: &[[f64; 3]],
    sigma: f64,
    e_min: f64,
    e_max: f64,
    n_points: usize,
) -> Result<dos::DosResult, String> {
    let eht_result = eht::solve_eht(elements, positions, None)?;
    let flat: Vec<f64> = positions.iter().flat_map(|p| p.iter().copied()).collect();

    #[cfg(feature = "parallel")]
    {
        Ok(dos::compute_pdos_parallel(
            elements,
            &flat,
            &eht_result.energies,
            &eht_result.coefficients,
            eht_result.n_electrons,
            sigma,
            e_min,
            e_max,
            n_points,
        ))
    }

    #[cfg(not(feature = "parallel"))]
    {
        Ok(dos::compute_pdos(
            elements,
            &flat,
            &eht_result.energies,
            &eht_result.coefficients,
            eht_result.n_electrons,
            sigma,
            e_min,
            e_max,
            n_points,
        ))
    }
}

/// Compute RMSD between two coordinate sets after Kabsch alignment.
pub fn compute_rmsd(coords: &[f64], reference: &[f64]) -> f64 {
    alignment::compute_rmsd(coords, reference)
}

/// Compute UFF force field energy for a molecule.
pub fn compute_uff_energy(smiles: &str, coords: &[f64]) -> Result<f64, String> {
    let mol = graph::Molecule::from_smiles(smiles)?;
    let n = mol.graph.node_count();
    if coords.len() != n * 3 {
        return Err(format!("coords length {} != 3 * atoms {}", coords.len(), n));
    }
    let ff = forcefield::builder::build_uff_force_field(&mol);
    let mut gradient = vec![0.0f64; n * 3];
    let energy = ff.compute_system_energy_and_gradients(coords, &mut gradient);
    Ok(energy)
}

/// Create a periodic unit cell from lattice parameters (a, b, c in Å; α, β, γ in degrees).
pub fn create_unit_cell(
    a: f64,
    b: f64,
    c: f64,
    alpha: f64,
    beta: f64,
    gamma: f64,
) -> materials::UnitCell {
    materials::UnitCell::from_parameters(&materials::CellParameters {
        a,
        b,
        c,
        alpha,
        beta,
        gamma,
    })
}

/// Assemble a framework crystal structure from node/linker SBUs on a topology.
///
/// Returns the assembled crystal structure as a JSON-serializable `CrystalStructure`.
pub fn assemble_framework(
    node: &materials::Sbu,
    linker: &materials::Sbu,
    topology: &materials::Topology,
    cell: &materials::UnitCell,
) -> materials::CrystalStructure {
    materials::assemble_framework(node, linker, topology, cell)
}
