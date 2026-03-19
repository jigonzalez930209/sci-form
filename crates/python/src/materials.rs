//! Unit cell and framework assembly bindings.

use pyo3::prelude::*;

#[pyclass]
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

#[pyclass]
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

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(unit_cell, m)?)?;
    m.add_function(wrap_pyfunction!(assemble_framework, m)?)?;
    m.add_class::<UnitCellPy>()?;
    m.add_class::<CrystalStructurePy>()?;
    Ok(())
}
