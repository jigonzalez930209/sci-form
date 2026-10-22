//! Unit cell and framework assembly bindings.

use pyo3::prelude::*;
use pyo3::types::PyDict;

#[pyclass(name = "UnitCell")]
#[derive(Clone)]
pub(crate) struct UnitCellPy {
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

#[pyclass(name = "CrystalStructure")]
#[derive(Clone)]
pub(crate) struct CrystalStructurePy {
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

#[pyclass(name = "SymmetryOp")]
#[derive(Clone)]
pub(crate) struct SymmetryOpPy {
    #[pyo3(get)]
    rotation: Vec<Vec<i8>>,
    #[pyo3(get)]
    translation: [f64; 3],
}

#[pyclass(name = "SpaceGroup")]
#[derive(Clone)]
pub(crate) struct SpaceGroupPy {
    #[pyo3(get)]
    number: u16,
    #[pyo3(get)]
    symbol: String,
    #[pyo3(get)]
    crystal_system: String,
    #[pyo3(get)]
    lattice_type: String,
    #[pyo3(get)]
    point_group: String,
    #[pyo3(get)]
    operations: Vec<SymmetryOpPy>,
}
#[pymethods]
impl SpaceGroupPy {
    fn __repr__(&self) -> String {
        format!(
            "SpaceGroup(number={}, symbol='{}', n_ops={})",
            self.number,
            self.symbol,
            self.operations.len()
        )
    }
}

#[pyclass(name = "PeriodicBond")]
#[derive(Clone)]
pub(crate) struct PeriodicBondPy {
    #[pyo3(get)]
    atom_i: usize,
    #[pyo3(get)]
    atom_j: usize,
    #[pyo3(get)]
    order: String,
    #[pyo3(get)]
    image: [i32; 3],
    #[pyo3(get)]
    distance: f64,
}

#[pyclass(name = "PeriodicMolecule")]
#[derive(Clone)]
pub(crate) struct PeriodicMoleculePy {
    #[pyo3(get)]
    n_atoms: usize,
    #[pyo3(get)]
    elements: Vec<u8>,
    #[pyo3(get)]
    frac_coords: Vec<[f64; 3]>,
    #[pyo3(get)]
    cart_coords: Vec<[f64; 3]>,
    #[pyo3(get)]
    bonds: Vec<PeriodicBondPy>,
    #[pyo3(get)]
    cell: UnitCellPy,
}
#[pymethods]
impl PeriodicMoleculePy {
    fn __repr__(&self) -> String {
        format!(
            "PeriodicMolecule(n_atoms={}, n_bonds={})",
            self.n_atoms,
            self.bonds.len()
        )
    }
}

#[pyclass(name = "HapticInteraction")]
#[derive(Clone)]
pub(crate) struct HapticInteractionPy {
    #[pyo3(get)]
    metal_index: usize,
    #[pyo3(get)]
    metal_element: u8,
    #[pyo3(get)]
    ring_atoms: Vec<usize>,
    #[pyo3(get)]
    hapticity: usize,
    #[pyo3(get)]
    centroid: [f64; 3],
    #[pyo3(get)]
    metal_centroid_distance: f64,
}

#[pyclass(name = "HapticityAnalysis")]
#[derive(Clone)]
pub(crate) struct HapticityAnalysisPy {
    #[pyo3(get)]
    interactions: Vec<HapticInteractionPy>,
    #[pyo3(get)]
    is_metallocene: bool,
    #[pyo3(get)]
    n_interactions: usize,
}
#[pymethods]
impl HapticityAnalysisPy {
    fn __repr__(&self) -> String {
        format!(
            "HapticityAnalysis(n_interactions={}, is_metallocene={})",
            self.n_interactions, self.is_metallocene
        )
    }
}

#[pyclass(name = "FrameworkOptResult")]
#[derive(Clone)]
pub(crate) struct FrameworkOptResultPy {
    #[pyo3(get)]
    positions: Vec<[f64; 3]>,
    #[pyo3(get)]
    energy: f64,
    #[pyo3(get)]
    forces: Vec<[f64; 3]>,
    #[pyo3(get)]
    max_force: f64,
    #[pyo3(get)]
    n_iterations: usize,
    #[pyo3(get)]
    converged: bool,
    #[pyo3(get)]
    energy_history: Vec<f64>,
}
#[pymethods]
impl FrameworkOptResultPy {
    fn __repr__(&self) -> String {
        format!(
            "FrameworkOptResult(n_atoms={}, iterations={}, converged={}, E={:.4} eV)",
            self.positions.len(),
            self.n_iterations,
            self.converged,
            self.energy
        )
    }
}

/// Convert a list of `[i, j]` atom-index pairs into bond tuples.
fn parse_bond_pairs(bonds: &[Vec<usize>]) -> PyResult<Vec<(usize, usize)>> {
    let mut out = Vec::with_capacity(bonds.len());
    for (k, b) in bonds.iter().enumerate() {
        if b.len() != 2 {
            return Err(crate::errors::value_error(format!(
                "bonds[{}] must have exactly 2 atom indices, got {}",
                k,
                b.len()
            )));
        }
        out.push((b[0], b[1]));
    }
    Ok(out)
}

fn map_space_group(sg: sci_form_core::materials::space_groups::SpaceGroup) -> SpaceGroupPy {
    SpaceGroupPy {
        number: sg.number,
        symbol: sg.symbol.clone(),
        crystal_system: format!("{:?}", sg.crystal_system),
        lattice_type: format!("{:?}", sg.lattice_type),
        point_group: sg.point_group.clone(),
        operations: sg
            .operations
            .iter()
            .map(|op| SymmetryOpPy {
                rotation: op.rotation.iter().map(|r| r.to_vec()).collect(),
                translation: op.translation,
            })
            .collect(),
    }
}

#[pyfunction]
#[pyo3(signature = (a, b, c, alpha=90.0, beta=90.0, gamma=90.0))]
fn create_unit_cell(a: f64, b: f64, c: f64, alpha: f64, beta: f64, gamma: f64) -> UnitCellPy {
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
            return Err(crate::errors::value_error(format!(
                "Invalid 'geometry' '{}'. Supported: 'linear', 'trigonal', 'tetrahedral', 'square_planar', 'octahedral'.",
                geometry
            )))
        }
    };
    let topo = match topology {
        "pcu" => sci_form_core::materials::Topology::pcu(),
        "dia" => sci_form_core::materials::Topology::dia(),
        "sql" => sci_form_core::materials::Topology::sql(),
        _ => {
            return Err(crate::errors::value_error(format!(
                "Invalid 'topology' '{}'. Supported: 'pcu', 'dia', 'sql'.",
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

/// Look up a crystallographic space group by ITC number (1–230).
#[pyfunction]
fn space_group_by_number(number: u16) -> Option<SpaceGroupPy> {
    sci_form_core::materials::space_groups::space_group_by_number(number).map(map_space_group)
}

/// Look up a crystallographic space group by Hermann-Mauguin symbol.
#[pyfunction]
fn space_group_by_symbol(symbol: &str) -> Option<SpaceGroupPy> {
    sci_form_core::materials::space_groups::space_group_by_symbol(symbol).map(map_space_group)
}

#[pyclass(name = "SpaceGroupEntry")]
#[derive(Clone)]
pub(crate) struct SpaceGroupEntryPy {
    #[pyo3(get)]
    number: u16,
    #[pyo3(get)]
    symbol: String,
}
#[pymethods]
impl SpaceGroupEntryPy {
    fn __repr__(&self) -> String {
        format!(
            "SpaceGroupEntry(number={}, symbol='{}')",
            self.number, self.symbol
        )
    }
}

/// List all 230 space groups with their ITC number and Hermann-Mauguin symbol.
#[pyfunction]
fn all_space_groups() -> Vec<SpaceGroupEntryPy> {
    sci_form_core::materials::space_groups::all_space_groups()
        .into_iter()
        .map(|(number, symbol)| SpaceGroupEntryPy {
            number,
            symbol: symbol.to_string(),
        })
        .collect()
}

/// Build a periodic molecular graph from fractional coordinates and a unit cell.
///
/// `frac_coords` is a list of `[x, y, z]` fractional triples; the cell uses
/// the same `(a, b, c, alpha, beta, gamma)` parameters as [`create_unit_cell`].
/// `bond_tolerance` is the extra bonding tolerance in Å (default 0.3).
#[pyfunction]
#[pyo3(signature = (elements, frac_coords, a, b, c, alpha=90.0, beta=90.0, gamma=90.0, bond_tolerance=None, options=None))]
fn build_periodic_molecule(
    elements: Vec<u8>,
    frac_coords: Vec<Vec<f64>>,
    a: f64,
    b: f64,
    c: f64,
    alpha: f64,
    beta: f64,
    gamma: f64,
    bond_tolerance: Option<f64>,
    options: Option<&Bound<'_, PyDict>>,
) -> PyResult<PeriodicMoleculePy> {
    // The optional `options` mapping overrides the positional arguments
    // (blocker C-M3).
    const KEYS: &[&str] = &["a", "b", "c", "alpha", "beta", "gamma", "bond_tolerance"];
    crate::options::reject_unknown(options, KEYS)?;
    let a = crate::options::f64_from(options, "a", a)?;
    let b = crate::options::f64_from(options, "b", b)?;
    let c = crate::options::f64_from(options, "c", c)?;
    let alpha = crate::options::f64_from(options, "alpha", alpha)?;
    let beta = crate::options::f64_from(options, "beta", beta)?;
    let gamma = crate::options::f64_from(options, "gamma", gamma)?;
    let bond_tolerance = crate::options::opt_f64_from(options, "bond_tolerance", bond_tolerance)?;
    let mut frac: Vec<[f64; 3]> = Vec::with_capacity(frac_coords.len());
    for (i, f) in frac_coords.iter().enumerate() {
        if f.len() != 3 {
            return Err(crate::errors::value_error(format!(
                "frac_coords[{}] must have 3 components, got {}",
                i,
                f.len()
            )));
        }
        frac.push([f[0], f[1], f[2]]);
    }
    let unit = sci_form_core::create_unit_cell(a, b, c, alpha, beta, gamma);
    let pm =
        sci_form_core::periodic::build_periodic_molecule(&elements, &frac, &unit, bond_tolerance);
    let p = unit.parameters();
    Ok(PeriodicMoleculePy {
        n_atoms: pm.n_atoms,
        elements: pm.elements,
        frac_coords: pm.frac_coords.clone(),
        cart_coords: pm.cart_coords.clone(),
        bonds: pm
            .bonds
            .iter()
            .map(|b| PeriodicBondPy {
                atom_i: b.atom_i,
                atom_j: b.atom_j,
                order: b.order.clone(),
                image: b.image,
                distance: b.distance,
            })
            .collect(),
        cell: UnitCellPy {
            a: p.a,
            b: p.b,
            c: p.c,
            alpha: p.alpha,
            beta: p.beta,
            gamma: p.gamma,
            volume: unit.volume(),
            lattice: unit.lattice.iter().map(|r| r.to_vec()).collect(),
        },
    })
}

/// Detect haptic (η) metal–ring interactions and metallocene motifs.
///
/// `bonds` is an optional list of `[atom_i, atom_j]` index pairs; when
/// omitted the bonds are perceived from the geometry.
#[pyfunction]
#[pyo3(signature = (elements, coords, bonds=None))]
fn detect_hapticity(
    elements: Vec<u8>,
    coords: Vec<f64>,
    bonds: Option<Vec<Vec<usize>>>,
) -> PyResult<HapticityAnalysisPy> {
    let positions = crate::system::coords_to_positions(&coords);
    let bond_list = match bonds {
        Some(b) => parse_bond_pairs(&b)?,
        None => crate::system::perceive_bonds(&elements, &positions, 0.3),
    };
    let analysis = sci_form_core::periodic::detect_hapticity(&elements, &positions, &bond_list);
    Ok(HapticityAnalysisPy {
        interactions: analysis
            .interactions
            .iter()
            .map(|h| HapticInteractionPy {
                metal_index: h.metal_index,
                metal_element: h.metal_element,
                ring_atoms: h.ring_atoms.clone(),
                hapticity: h.hapticity,
                centroid: h.centroid,
                metal_centroid_distance: h.metal_centroid_distance,
            })
            .collect(),
        is_metallocene: analysis.is_metallocene,
        n_interactions: analysis.n_interactions,
    })
}

/// Optimize a molecular/framework geometry with the shared BFGS engine.
///
/// The energy/force backend is the UFF force field built from `smiles`;
/// energies are reported in eV and forces in eV/Å (UFF kcal/mol converted
/// with 1 eV = 23.0605 kcal/mol). Supplying `cell_a`/`cell_b`/`cell_c`
/// enables periodic (minimum-image) optimization using the given
/// `(a, b, c, alpha, beta, gamma)` unit cell.
#[pyfunction]
#[pyo3(signature = (smiles, coords, method="bfgs", max_iter=200, force_tol=0.05, energy_tol=1e-6, max_step=0.2, fixed_atoms=None, cell_a=None, cell_b=None, cell_c=None, cell_alpha=90.0, cell_beta=90.0, cell_gamma=90.0, options=None))]
fn optimize_framework(
    smiles: &str,
    coords: Vec<f64>,
    method: &str,
    max_iter: usize,
    force_tol: f64,
    energy_tol: f64,
    max_step: f64,
    fixed_atoms: Option<Vec<usize>>,
    cell_a: Option<f64>,
    cell_b: Option<f64>,
    cell_c: Option<f64>,
    cell_alpha: f64,
    cell_beta: f64,
    cell_gamma: f64,
    options: Option<&Bound<'_, PyDict>>,
) -> PyResult<FrameworkOptResultPy> {
    // The optional `options` mapping overrides the positional arguments
    // (blocker C-M3); the positional form keeps working unchanged.
    const KEYS: &[&str] = &[
        "method",
        "max_iter",
        "force_tol",
        "energy_tol",
        "max_step",
        "fixed_atoms",
        "cell_a",
        "cell_b",
        "cell_c",
        "cell_alpha",
        "cell_beta",
        "cell_gamma",
    ];
    crate::options::reject_unknown(options, KEYS)?;
    let method = &crate::options::str_from(options, "method", method)?;
    let method = method.as_str();
    let max_iter = crate::options::usize_from(options, "max_iter", max_iter)?;
    let force_tol = crate::options::f64_from(options, "force_tol", force_tol)?;
    let energy_tol = crate::options::f64_from(options, "energy_tol", energy_tol)?;
    let max_step = crate::options::f64_from(options, "max_step", max_step)?;
    let fixed_atoms = crate::options::usize_list_from(options, "fixed_atoms", fixed_atoms)?;
    let cell_a = crate::options::opt_f64_from(options, "cell_a", cell_a)?;
    let cell_b = crate::options::opt_f64_from(options, "cell_b", cell_b)?;
    let cell_c = crate::options::opt_f64_from(options, "cell_c", cell_c)?;
    let cell_alpha = crate::options::f64_from(options, "cell_alpha", cell_alpha)?;
    let cell_beta = crate::options::f64_from(options, "cell_beta", cell_beta)?;
    let cell_gamma = crate::options::f64_from(options, "cell_gamma", cell_gamma)?;
    let mol = sci_form_core::parse(smiles)
        .map_err(|e| crate::errors::value_error(crate::errors::smiles_error(&e, smiles)))?;
    let n = mol.graph.node_count();
    if coords.len() != n * 3 {
        return Err(crate::errors::value_error(format!(
            "coords length {} != 3 * atoms {}",
            coords.len(),
            n
        )));
    }
    let elements: Vec<u8> = (0..n)
        .map(|i| mol.graph[sci_form_core::graph::NodeIndex::new(i)].element)
        .collect();
    let initial_positions = crate::system::coords_to_positions(&coords);

    let opt_method = match method.to_ascii_lowercase().as_str() {
        "bfgs" => sci_form_core::materials::geometry_opt::OptMethod::Bfgs,
        "steepest_descent" | "sd" => {
            sci_form_core::materials::geometry_opt::OptMethod::SteepestDescent
        }
        other => {
            return Err(crate::errors::value_error(format!(
                "Invalid 'method' '{}'. Supported: 'bfgs', 'steepest_descent', 'sd'.",
                other
            )))
        }
    };

    let lattice = match (cell_a, cell_b, cell_c) {
        (Some(a), Some(b), Some(c)) => Some(
            sci_form_core::create_unit_cell(a, b, c, cell_alpha, cell_beta, cell_gamma).lattice,
        ),
        _ => None,
    };

    let config = sci_form_core::materials::geometry_opt::FrameworkOptConfig {
        method: opt_method,
        max_iter,
        force_tol,
        energy_tol,
        max_step,
        fixed_atoms: fixed_atoms.unwrap_or_default(),
        lattice,
    };

    // UFF-backed energy/force closure (UFF works in kcal/mol and kcal/mol/Å).
    let ff = sci_form_core::forcefield::builder::build_uff_force_field(&mol);
    const KCAL_PER_EV: f64 = 23.0605;
    let energy_force = |_els: &[u8], pos: &[[f64; 3]]| -> Result<(f64, Vec<[f64; 3]>), String> {
        let flat: Vec<f64> = pos.iter().flat_map(|p| p.iter().copied()).collect();
        let mut grad = vec![0.0f64; flat.len()];
        let energy_kcal = ff.compute_system_energy_and_gradients(&flat, &mut grad);
        let mut forces = vec![[0.0f64; 3]; pos.len()];
        for (i, f) in forces.iter_mut().enumerate() {
            f[0] = -grad[3 * i] / KCAL_PER_EV;
            f[1] = -grad[3 * i + 1] / KCAL_PER_EV;
            f[2] = -grad[3 * i + 2] / KCAL_PER_EV;
        }
        Ok((energy_kcal / KCAL_PER_EV, forces))
    };

    let result = sci_form_core::materials::geometry_opt::optimize_framework(
        &elements,
        &initial_positions,
        &energy_force,
        &config,
    )
    .map_err(|e| crate::errors::smiles_or_engine_error(e, smiles))?;

    Ok(FrameworkOptResultPy {
        positions: result.positions.clone(),
        energy: result.energy,
        forces: result.forces.clone(),
        max_force: result.max_force,
        n_iterations: result.n_iterations,
        converged: result.converged,
        energy_history: result.energy_history,
    })
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(create_unit_cell, m)?)?;
    m.add_function(wrap_pyfunction!(assemble_framework, m)?)?;
    m.add_function(wrap_pyfunction!(space_group_by_number, m)?)?;
    m.add_function(wrap_pyfunction!(space_group_by_symbol, m)?)?;
    m.add_function(wrap_pyfunction!(all_space_groups, m)?)?;
    m.add_function(wrap_pyfunction!(build_periodic_molecule, m)?)?;
    m.add_function(wrap_pyfunction!(detect_hapticity, m)?)?;
    m.add_function(wrap_pyfunction!(optimize_framework, m)?)?;
    m.add_class::<UnitCellPy>()?;
    m.add_class::<CrystalStructurePy>()?;
    m.add_class::<SymmetryOpPy>()?;
    m.add_class::<SpaceGroupPy>()?;
    m.add_class::<SpaceGroupEntryPy>()?;
    m.add_class::<PeriodicBondPy>()?;
    m.add_class::<PeriodicMoleculePy>()?;
    m.add_class::<HapticInteractionPy>()?;
    m.add_class::<HapticityAnalysisPy>()?;
    m.add_class::<FrameworkOptResultPy>()?;
    Ok(())
}
