use pyo3::prelude::*;
use pyo3::types::PyDict;

/// A single conformer result returned to Python.
#[pyclass]
#[derive(Clone)]
struct ConformerResult {
    #[pyo3(get)]
    smiles: String,
    #[pyo3(get)]
    num_atoms: usize,
    #[pyo3(get)]
    coords: Vec<f64>,
    #[pyo3(get)]
    elements: Vec<u8>,
    #[pyo3(get)]
    bonds: Vec<(usize, usize, String)>,
    #[pyo3(get)]
    error: Option<String>,
    #[pyo3(get)]
    time_ms: f64,
}

#[pymethods]
impl ConformerResult {
    /// Get coordinates as list of (x, y, z) tuples.
    fn get_positions(&self) -> Vec<(f64, f64, f64)> {
        self.coords
            .chunks_exact(3)
            .map(|c| (c[0], c[1], c[2]))
            .collect()
    }

    /// Check if generation succeeded.
    fn is_ok(&self) -> bool {
        self.error.is_none()
    }

    fn __repr__(&self) -> String {
        if let Some(ref e) = self.error {
            format!("ConformerResult(smiles='{}', error='{}')", self.smiles, e)
        } else {
            format!(
                "ConformerResult(smiles='{}', atoms={}, time={:.1}ms)",
                self.smiles, self.num_atoms, self.time_ms
            )
        }
    }
}

impl From<sci_form_core::ConformerResult> for ConformerResult {
    fn from(r: sci_form_core::ConformerResult) -> Self {
        ConformerResult {
            smiles: r.smiles,
            num_atoms: r.num_atoms,
            coords: r.coords,
            elements: r.elements,
            bonds: r.bonds,
            error: r.error,
            time_ms: r.time_ms,
        }
    }
}

/// Generate a 3D conformer from a SMILES string.
///
/// Args:
///     smiles: SMILES string of the molecule.
///     seed: RNG seed for reproducibility (default: 42).
///
/// Returns:
///     ConformerResult with 3D coordinates.
#[pyfunction]
#[pyo3(signature = (smiles, seed=42))]
fn embed(smiles: &str, seed: u64) -> ConformerResult {
    sci_form_core::embed(smiles, seed).into()
}

/// Batch-generate 3D conformers for multiple SMILES in parallel.
///
/// Args:
///     smiles_list: List of SMILES strings.
///     seed: RNG seed (default: 42).
///     num_threads: Number of threads, 0 = auto (default: 0).
///
/// Returns:
///     List of ConformerResult objects.
#[pyfunction]
#[pyo3(signature = (smiles_list, seed=42, num_threads=0))]
fn embed_batch(smiles_list: Vec<String>, seed: u64, num_threads: usize) -> Vec<ConformerResult> {
    let refs: Vec<&str> = smiles_list.iter().map(|s| s.as_str()).collect();
    let config = sci_form_core::ConformerConfig { seed, num_threads };
    sci_form_core::embed_batch(&refs, &config)
        .into_iter()
        .map(|r| r.into())
        .collect()
}

/// Parse a SMILES string without generating 3D coordinates.
///
/// Returns a dict with atom/bond info.
#[pyfunction]
fn parse(py: Python<'_>, smiles: &str) -> PyResult<PyObject> {
    match sci_form_core::parse(smiles) {
        Ok(mol) => {
            let dict = PyDict::new_bound(py);
            let n = mol.graph.node_count();
            dict.set_item("num_atoms", n)?;
            dict.set_item("num_bonds", mol.graph.edge_count())?;

            let mut atoms_list: Vec<PyObject> = Vec::new();
            for i in 0..n {
                let idx = sci_form_core::graph::NodeIndex::new(i);
                let atom = &mol.graph[idx];
                let adict = PyDict::new_bound(py);
                adict.set_item("element", atom.element)?;
                adict.set_item("hybridization", format!("{:?}", atom.hybridization))?;
                adict.set_item("formal_charge", atom.formal_charge)?;
                atoms_list.push(adict.into());
            }
            dict.set_item("atoms", atoms_list)?;
            Ok(dict.into())
        }
        Err(e) => Err(pyo3::exceptions::PyValueError::new_err(e)),
    }
}

/// Get the library version string.
#[pyfunction]
fn version() -> String {
    sci_form_core::version()
}

/// sci_form Python module
#[pymodule]
fn sci_form(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(embed, m)?)?;
    m.add_function(wrap_pyfunction!(embed_batch, m)?)?;
    m.add_function(wrap_pyfunction!(parse, m)?)?;
    m.add_function(wrap_pyfunction!(version, m)?)?;
    m.add_function(wrap_pyfunction!(eht_calculate, m)?)?;
    m.add_function(wrap_pyfunction!(eht_orbital_mesh, m)?)?;
    m.add_function(wrap_pyfunction!(charges, m)?)?;
    m.add_function(wrap_pyfunction!(sasa, m)?)?;
    m.add_function(wrap_pyfunction!(population, m)?)?;
    m.add_function(wrap_pyfunction!(dipole, m)?)?;
    m.add_function(wrap_pyfunction!(dos, m)?)?;
    m.add_function(wrap_pyfunction!(rmsd, m)?)?;
    m.add_function(wrap_pyfunction!(uff_energy, m)?)?;
    m.add_function(wrap_pyfunction!(unit_cell, m)?)?;
    m.add_function(wrap_pyfunction!(assemble_framework, m)?)?;
    m.add_function(wrap_pyfunction!(pack_conformers, m)?)?;
    m.add_function(wrap_pyfunction!(split_worker_tasks, m)?)?;
    m.add_function(wrap_pyfunction!(estimate_workers, m)?)?;
    m.add_class::<ConformerResult>()?;
    m.add_class::<EhtResultPy>()?;
    m.add_class::<ChargeResultPy>()?;
    m.add_class::<SasaResultPy>()?;
    m.add_class::<PopulationResultPy>()?;
    m.add_class::<DipoleResultPy>()?;
    m.add_class::<DosResultPy>()?;
    m.add_class::<AlignmentResultPy>()?;
    m.add_class::<UnitCellPy>()?;
    m.add_class::<CrystalStructurePy>()?;
    m.add_class::<RecordBatchPy>()?;
    Ok(())
}

// ─── EHT (Extended Hückel Theory) API ────────────────────────────────────────

/// EHT calculation result.
#[pyclass]
#[derive(Clone)]
struct EhtResultPy {
    #[pyo3(get)]
    energies: Vec<f64>,
    #[pyo3(get)]
    n_electrons: usize,
    #[pyo3(get)]
    homo_index: usize,
    #[pyo3(get)]
    lumo_index: usize,
    #[pyo3(get)]
    homo_energy: f64,
    #[pyo3(get)]
    lumo_energy: f64,
    #[pyo3(get)]
    gap: f64,
}

#[pymethods]
impl EhtResultPy {
    fn __repr__(&self) -> String {
        format!(
            "EhtResult(n_mo={}, gap={:.3} eV, HOMO={:.3} eV, LUMO={:.3} eV)",
            self.energies.len(),
            self.gap,
            self.homo_energy,
            self.lumo_energy,
        )
    }
}

/// Run an EHT calculation on a molecule.
///
/// Args:
///     elements: List of atomic numbers (e.g. [8, 1, 1] for water).
///     coords: Flat list of xyz coordinates [x0,y0,z0, x1,y1,z1, ...].
///     k: Wolfsberg-Helmholtz constant (default: 1.75).
///
/// Returns:
///     EhtResult with orbital energies, HOMO/LUMO info.
#[pyfunction]
#[pyo3(signature = (elements, coords, k=1.75))]
fn eht_calculate(elements: Vec<u8>, coords: Vec<f64>, k: f64) -> PyResult<EhtResultPy> {
    if coords.len() != elements.len() * 3 {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "coords length must be 3 * len(elements)",
        ));
    }
    let positions: Vec<[f64; 3]> = coords.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect();
    let k_opt = if (k - 1.75).abs() < 1e-10 {
        None
    } else {
        Some(k)
    };

    match sci_form_core::eht::solve_eht(&elements, &positions, k_opt) {
        Ok(r) => Ok(EhtResultPy {
            energies: r.energies,
            n_electrons: r.n_electrons,
            homo_index: r.homo_index,
            lumo_index: r.lumo_index,
            homo_energy: r.homo_energy,
            lumo_energy: r.lumo_energy,
            gap: r.gap,
        }),
        Err(e) => Err(pyo3::exceptions::PyRuntimeError::new_err(e)),
    }
}

/// Generate an orbital isosurface mesh as a dict.
///
/// Args:
///     elements: List of atomic numbers.
///     coords: Flat xyz coordinates.
///     mo_index: Molecular orbital index (0-based).
///     spacing: Grid spacing in Ångström (default: 0.2).
///     isovalue: Isosurface cutoff (default: 0.02).
///
/// Returns:
///     Dict with 'vertices', 'normals', 'indices', 'num_triangles'.
#[pyfunction]
#[pyo3(signature = (elements, coords, mo_index, spacing=0.2, isovalue=0.02))]
fn eht_orbital_mesh(
    py: Python<'_>,
    elements: Vec<u8>,
    coords: Vec<f64>,
    mo_index: usize,
    spacing: f64,
    isovalue: f32,
) -> PyResult<PyObject> {
    if coords.len() != elements.len() * 3 {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "coords length must be 3 * len(elements)",
        ));
    }
    let positions: Vec<[f64; 3]> = coords.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect();

    let result = sci_form_core::eht::solve_eht(&elements, &positions, None)
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e))?;

    let basis = sci_form_core::eht::basis::build_basis(&elements, &positions);
    let grid = sci_form_core::eht::evaluate_orbital_on_grid(
        &basis,
        &result.coefficients,
        mo_index,
        &positions,
        spacing,
        3.0,
    );
    let mesh = sci_form_core::eht::marching_cubes(&grid, isovalue);

    let dict = PyDict::new_bound(py);
    dict.set_item("vertices", mesh.vertices)?;
    dict.set_item("normals", mesh.normals)?;
    dict.set_item("indices", mesh.indices)?;
    dict.set_item("num_triangles", mesh.num_triangles)?;
    dict.set_item("isovalue", mesh.isovalue)?;
    Ok(dict.into())
}

// ─── Charges (Gasteiger-Marsili) API ─────────────────────────────────────────

/// Gasteiger-Marsili charge result.
#[pyclass]
#[derive(Clone)]
struct ChargeResultPy {
    #[pyo3(get)]
    charges: Vec<f64>,
    #[pyo3(get)]
    iterations: usize,
    #[pyo3(get)]
    total_charge: f64,
}

#[pymethods]
impl ChargeResultPy {
    fn __repr__(&self) -> String {
        format!(
            "ChargeResult(n_atoms={}, total_charge={:.4}, iterations={})",
            self.charges.len(),
            self.total_charge,
            self.iterations,
        )
    }
}

/// Compute Gasteiger-Marsili partial charges from a SMILES string.
///
/// Args:
///     smiles: SMILES string of the molecule.
///
/// Returns:
///     ChargeResult with partial charges per atom.
#[pyfunction]
fn charges(smiles: &str) -> PyResult<ChargeResultPy> {
    match sci_form_core::compute_charges(smiles) {
        Ok(r) => Ok(ChargeResultPy {
            charges: r.charges,
            iterations: r.iterations,
            total_charge: r.total_charge,
        }),
        Err(e) => Err(pyo3::exceptions::PyRuntimeError::new_err(e)),
    }
}

// ─── SASA (Shrake-Rupley) API ────────────────────────────────────────────────

/// SASA calculation result.
#[pyclass]
#[derive(Clone)]
struct SasaResultPy {
    #[pyo3(get)]
    total_sasa: f64,
    #[pyo3(get)]
    atom_sasa: Vec<f64>,
    #[pyo3(get)]
    probe_radius: f64,
    #[pyo3(get)]
    num_points: usize,
}

#[pymethods]
impl SasaResultPy {
    fn __repr__(&self) -> String {
        format!(
            "SasaResult(total={:.2} Ų, n_atoms={}, probe={:.2} Å)",
            self.total_sasa,
            self.atom_sasa.len(),
            self.probe_radius,
        )
    }
}

/// Compute solvent-accessible surface area.
///
/// Args:
///     elements: List of atomic numbers.
///     coords: Flat list of xyz coordinates [x0,y0,z0, x1,y1,z1, ...].
///     probe_radius: Probe radius in Å (default: 1.4).
///
/// Returns:
///     SasaResult with total and per-atom SASA values.
#[pyfunction]
#[pyo3(signature = (elements, coords, probe_radius=1.4))]
fn sasa(elements: Vec<u8>, coords: Vec<f64>, probe_radius: f64) -> PyResult<SasaResultPy> {
    let pr = if (probe_radius - 1.4).abs() < 1e-10 {
        None
    } else {
        Some(probe_radius)
    };
    match sci_form_core::compute_sasa(&elements, &coords, pr) {
        Ok(r) => Ok(SasaResultPy {
            total_sasa: r.total_sasa,
            atom_sasa: r.atom_sasa,
            probe_radius: r.probe_radius,
            num_points: r.num_points,
        }),
        Err(e) => Err(pyo3::exceptions::PyValueError::new_err(e)),
    }
}

// ─── Population Analysis (Mulliken & Löwdin) API ─────────────────────────────

#[pyclass]
#[derive(Clone)]
struct PopulationResultPy {
    #[pyo3(get)]
    mulliken_charges: Vec<f64>,
    #[pyo3(get)]
    lowdin_charges: Vec<f64>,
    #[pyo3(get)]
    num_atoms: usize,
    #[pyo3(get)]
    total_charge_mulliken: f64,
    #[pyo3(get)]
    total_charge_lowdin: f64,
}

#[pymethods]
impl PopulationResultPy {
    fn __repr__(&self) -> String {
        format!(
            "PopulationResult(n_atoms={}, total_mulliken={:.4})",
            self.num_atoms, self.total_charge_mulliken
        )
    }
}

/// Compute Mulliken & Löwdin population analysis.
#[pyfunction]
fn population(elements: Vec<u8>, coords: Vec<f64>) -> PyResult<PopulationResultPy> {
    let positions: Vec<[f64; 3]> = coords.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect();
    match sci_form_core::compute_population(&elements, &positions) {
        Ok(r) => Ok(PopulationResultPy {
            mulliken_charges: r.mulliken_charges,
            lowdin_charges: r.lowdin_charges,
            num_atoms: r.num_atoms,
            total_charge_mulliken: r.total_charge_mulliken,
            total_charge_lowdin: r.total_charge_lowdin,
        }),
        Err(e) => Err(pyo3::exceptions::PyRuntimeError::new_err(e)),
    }
}

// ─── Dipole Moment API ──────────────────────────────────────────────────────

#[pyclass]
#[derive(Clone)]
struct DipoleResultPy {
    #[pyo3(get)]
    vector: [f64; 3],
    #[pyo3(get)]
    magnitude: f64,
    #[pyo3(get)]
    unit: String,
}

#[pymethods]
impl DipoleResultPy {
    fn __repr__(&self) -> String {
        format!("DipoleResult({:.3} {})", self.magnitude, self.unit)
    }
}

/// Compute molecular dipole moment.
#[pyfunction]
fn dipole(elements: Vec<u8>, coords: Vec<f64>) -> PyResult<DipoleResultPy> {
    let positions: Vec<[f64; 3]> = coords.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect();
    match sci_form_core::compute_dipole(&elements, &positions) {
        Ok(r) => Ok(DipoleResultPy {
            vector: r.vector,
            magnitude: r.magnitude,
            unit: r.unit,
        }),
        Err(e) => Err(pyo3::exceptions::PyRuntimeError::new_err(e)),
    }
}

// ─── DOS/PDOS API ───────────────────────────────────────────────────────────

#[pyclass]
#[derive(Clone)]
struct DosResultPy {
    #[pyo3(get)]
    energies: Vec<f64>,
    #[pyo3(get)]
    total_dos: Vec<f64>,
    #[pyo3(get)]
    sigma: f64,
}

#[pymethods]
impl DosResultPy {
    fn __repr__(&self) -> String {
        format!(
            "DosResult(n_points={}, sigma={:.2} eV)",
            self.energies.len(),
            self.sigma
        )
    }
}

/// Compute density of states (DOS/PDOS).
#[pyfunction]
#[pyo3(signature = (elements, coords, sigma=0.3, e_min=-30.0, e_max=5.0, n_points=500))]
fn dos(
    elements: Vec<u8>,
    coords: Vec<f64>,
    sigma: f64,
    e_min: f64,
    e_max: f64,
    n_points: usize,
) -> PyResult<DosResultPy> {
    let positions: Vec<[f64; 3]> = coords.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect();
    match sci_form_core::compute_dos(&elements, &positions, sigma, e_min, e_max, n_points) {
        Ok(r) => Ok(DosResultPy {
            energies: r.energies,
            total_dos: r.total_dos,
            sigma: r.sigma,
        }),
        Err(e) => Err(pyo3::exceptions::PyRuntimeError::new_err(e)),
    }
}

// ─── RMSD / Kabsch Alignment API ────────────────────────────────────────────

#[pyclass]
#[derive(Clone)]
struct AlignmentResultPy {
    #[pyo3(get)]
    rmsd: f64,
    #[pyo3(get)]
    aligned_coords: Vec<f64>,
}

#[pymethods]
impl AlignmentResultPy {
    fn __repr__(&self) -> String {
        format!("AlignmentResult(rmsd={:.4} Å)", self.rmsd)
    }
}

/// Compute RMSD between two coordinate sets (Kabsch alignment).
#[pyfunction]
fn rmsd(coords: Vec<f64>, reference: Vec<f64>) -> PyResult<AlignmentResultPy> {
    let result = sci_form_core::alignment::align_coordinates(&coords, &reference);
    Ok(AlignmentResultPy {
        rmsd: result.rmsd,
        aligned_coords: result.aligned_coords,
    })
}

// ─── UFF Force Field API ────────────────────────────────────────────────────

/// Compute UFF force field energy for a molecule.
#[pyfunction]
fn uff_energy(smiles: &str, coords: Vec<f64>) -> PyResult<f64> {
    match sci_form_core::compute_uff_energy(smiles, &coords) {
        Ok(e) => Ok(e),
        Err(e) => Err(pyo3::exceptions::PyRuntimeError::new_err(e)),
    }
}

// ─── Materials Assembly API ──────────────────────────────────────────────────

/// Unit cell result with parameters and volume.
#[pyclass]
#[derive(Clone)]
struct UnitCellPy {
    #[pyo3(get)]
    a: f64,
    #[pyo3(get)]
    b: f64,
    #[pyo3(get)]
    c: f64,
    #[pyo3(get)]
    alpha: f64,
    #[pyo3(get)]
    beta: f64,
    #[pyo3(get)]
    gamma: f64,
    #[pyo3(get)]
    volume: f64,
    #[pyo3(get)]
    lattice: Vec<Vec<f64>>,
}

/// Crystal structure result.
#[pyclass]
#[derive(Clone)]
struct CrystalStructurePy {
    #[pyo3(get)]
    num_atoms: usize,
    #[pyo3(get)]
    elements: Vec<u8>,
    #[pyo3(get)]
    frac_coords: Vec<Vec<f64>>,
    #[pyo3(get)]
    cart_coords: Vec<Vec<f64>>,
    #[pyo3(get)]
    labels: Vec<String>,
    #[pyo3(get)]
    lattice: Vec<Vec<f64>>,
}

/// Create a unit cell from crystallographic parameters.
#[pyfunction]
#[pyo3(signature = (a, b, c, alpha=90.0, beta=90.0, gamma=90.0))]
fn unit_cell(a: f64, b: f64, c: f64, alpha: f64, beta: f64, gamma: f64) -> UnitCellPy {
    let cell = sci_form_core::create_unit_cell(a, b, c, alpha, beta, gamma);
    let p = cell.parameters();
    let vol = cell.volume();
    UnitCellPy {
        a: p.a,
        b: p.b,
        c: p.c,
        alpha: p.alpha,
        beta: p.beta,
        gamma: p.gamma,
        volume: vol,
        lattice: cell.lattice.iter().map(|r| r.to_vec()).collect(),
    }
}

/// Assemble a framework crystal structure.
#[pyfunction]
#[pyo3(signature = (topology="pcu", metal=30, geometry="octahedral", lattice_a=10.0, supercell=1))]
fn assemble_framework(
    topology: &str,
    metal: u8,
    geometry: &str,
    lattice_a: f64,
    supercell: usize,
) -> PyResult<CrystalStructurePy> {
    let geom = match geometry {
        "linear" => sci_form_core::materials::CoordinationGeometry::Linear,
        "trigonal" => sci_form_core::materials::CoordinationGeometry::Trigonal,
        "tetrahedral" => sci_form_core::materials::CoordinationGeometry::Tetrahedral,
        "square_planar" => sci_form_core::materials::CoordinationGeometry::SquarePlanar,
        "octahedral" => sci_form_core::materials::CoordinationGeometry::Octahedral,
        _ => {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "Unknown geometry: {}",
                geometry
            )))
        }
    };
    let topo = match topology {
        "pcu" => sci_form_core::materials::Topology::pcu(),
        "dia" => sci_form_core::materials::Topology::dia(),
        "sql" => sci_form_core::materials::Topology::sql(),
        _ => {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "Unknown topology: {}",
                topology
            )))
        }
    };
    let node = sci_form_core::materials::Sbu::metal_node(metal, 0.0, geom);
    let linker = sci_form_core::materials::Sbu::linear_linker(&[6, 6], 1.4, "carboxylate");
    let cell = sci_form_core::materials::UnitCell::cubic(lattice_a);
    let mut structure = sci_form_core::assemble_framework(&node, &linker, &topo, &cell);
    if supercell > 1 {
        structure = structure.make_supercell(supercell, supercell, supercell);
    }

    let cart = structure.cartesian_coords();
    Ok(CrystalStructurePy {
        num_atoms: structure.num_atoms(),
        elements: structure.elements(),
        frac_coords: structure
            .atoms
            .iter()
            .map(|a| a.frac_coords.to_vec())
            .collect(),
        cart_coords: cart.iter().map(|c| c.to_vec()).collect(),
        labels: structure.labels,
        lattice: structure.cell.lattice.iter().map(|r| r.to_vec()).collect(),
    })
}

// ─── C9: Transport / Streaming API ──────────────────────────────────────────

/// Arrow-compatible record batch result.
#[pyclass]
#[derive(Clone)]
struct RecordBatchPy {
    #[pyo3(get)]
    num_rows: usize,
    #[pyo3(get)]
    num_columns: usize,
    #[pyo3(get)]
    byte_size: usize,
    #[pyo3(get)]
    column_names: Vec<String>,
    #[pyo3(get)]
    float_data: std::collections::HashMap<String, Vec<f64>>,
    #[pyo3(get)]
    int_data: std::collections::HashMap<String, Vec<i32>>,
    #[pyo3(get)]
    uint8_data: std::collections::HashMap<String, Vec<u8>>,
}

/// Pack conformer results into Arrow-compatible columnar format.
#[pyfunction]
fn pack_conformers(results: Vec<ConformerResult>) -> RecordBatchPy {
    let core_results: Vec<sci_form_core::ConformerResult> = results
        .iter()
        .map(|r| sci_form_core::ConformerResult {
            smiles: r.smiles.clone(),
            num_atoms: r.num_atoms,
            coords: r.coords.clone(),
            elements: r.elements.clone(),
            bonds: vec![],
            error: r.error.clone(),
            time_ms: r.time_ms,
        })
        .collect();
    let batch = sci_form_core::transport::pack_conformers(&core_results);
    let mut float_data = std::collections::HashMap::new();
    let mut int_data = std::collections::HashMap::new();
    let mut uint8_data = std::collections::HashMap::new();
    let column_names: Vec<String> = batch.schema.iter().map(|s| s.name.clone()).collect();
    for c in &batch.float_columns {
        float_data.insert(c.name.clone(), c.values.clone());
    }
    for c in &batch.int_columns {
        int_data.insert(c.name.clone(), c.values.clone());
    }
    for c in &batch.uint8_columns {
        uint8_data.insert(c.name.clone(), c.values.clone());
    }
    RecordBatchPy {
        num_rows: batch.num_rows,
        num_columns: batch.num_columns(),
        byte_size: batch.byte_size(),
        column_names,
        float_data,
        int_data,
        uint8_data,
    }
}

/// Split SMILES batch into worker tasks.
#[pyfunction]
#[pyo3(signature = (smiles, n_workers=4, seed=42))]
fn split_worker_tasks(
    smiles: Vec<String>,
    n_workers: usize,
    seed: u64,
) -> Vec<std::collections::HashMap<String, String>> {
    let tasks = sci_form_core::transport::split_batch(&smiles, n_workers, seed);
    tasks
        .iter()
        .map(|t| {
            let mut m = std::collections::HashMap::new();
            m.insert("id".to_string(), t.id.to_string());
            // Format SMILES as JSON-like array string
            let smiles_str = format!(
                "[{}]",
                t.smiles
                    .iter()
                    .map(|s| format!("\"{}\"", s))
                    .collect::<Vec<_>>()
                    .join(",")
            );
            m.insert("smiles".to_string(), smiles_str);
            m.insert("kind".to_string(), format!("{:?}", t.kind));
            m
        })
        .collect()
}

/// Estimate optimal number of workers.
#[pyfunction]
#[pyo3(signature = (n_items, max_workers=8))]
fn estimate_workers(n_items: usize, max_workers: usize) -> usize {
    sci_form_core::transport::estimate_workers(n_items, max_workers)
}
