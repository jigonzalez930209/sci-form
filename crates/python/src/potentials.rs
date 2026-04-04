use pyo3::prelude::*;
use sci_form_core::potentials;

#[pyfunction]
fn double_well_1d(x: f64) -> PyResult<f64> {
    Ok(potentials::double_well_1d(x))
}

#[pyfunction]
fn asymmetric_1d(x: f64) -> PyResult<f64> {
    Ok(potentials::asymmetric_1d(x))
}

#[pyfunction]
fn sn2_model_1d(x: f64) -> PyResult<f64> {
    Ok(potentials::sn2_model_1d(x))
}

#[pyfunction]
fn muller_brown_2d(x: f64, y: f64) -> PyResult<f64> {
    Ok(potentials::muller_brown_2d(x, y))
}

#[pyfunction]
fn proton_transfer_evb(s: f64, k: f64, coupling: f64) -> PyResult<f64> {
    Ok(potentials::proton_transfer_evb(s, k, coupling))
}

#[pyfunction]
fn proton_transfer_2d_morse(r_ah: f64, r_hb: f64) -> PyResult<f64> {
    Ok(potentials::proton_transfer_2d_morse(r_ah, r_hb))
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    let sub = PyModule::new_bound(m.py(), "potentials")?;

    sub.add_function(wrap_pyfunction!(double_well_1d, &sub)?)?;
    sub.add_function(wrap_pyfunction!(asymmetric_1d, &sub)?)?;
    sub.add_function(wrap_pyfunction!(sn2_model_1d, &sub)?)?;
    sub.add_function(wrap_pyfunction!(muller_brown_2d, &sub)?)?;
    sub.add_function(wrap_pyfunction!(proton_transfer_evb, &sub)?)?;
    sub.add_function(wrap_pyfunction!(proton_transfer_2d_morse, &sub)?)?;

    m.add_submodule(&sub)?;

    // Register it to sys.modules so `from sci_form.potentials import xyz` works
    let sys = m.py().import_bound("sys")?;
    let sys_modules = sys.getattr("modules")?;
    sys_modules.set_item("sci_form.potentials", sub)?;

    Ok(())
}
