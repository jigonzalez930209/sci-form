//! Transport: Arrow-compatible batch packing, worker task splitting.

use crate::embed::ConformerResult;
use pyo3::prelude::*;

#[pyclass]
#[derive(Clone)]
pub(crate) struct RecordBatchPy {
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

#[pyfunction]
#[pyo3(signature = (n_items, max_workers=8))]
fn estimate_workers(n_items: usize, max_workers: usize) -> usize {
    sci_form_core::transport::estimate_workers(n_items, max_workers)
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(pack_conformers, m)?)?;
    m.add_function(wrap_pyfunction!(split_worker_tasks, m)?)?;
    m.add_function(wrap_pyfunction!(estimate_workers, m)?)?;
    m.add_class::<RecordBatchPy>()?;
    Ok(())
}
