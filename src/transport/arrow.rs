//! Arrow-compatible columnar memory layouts for zero-copy data transfer.
//!
//! Provides flat, typed-array compatible buffers that can be directly
//! transferred to JavaScript TypedArrays (Float64Array, Int32Array, etc.)
//! without serialization overhead.
//!
//! Memory layout matches Apache Arrow IPC format for interoperability:
//! - Values buffer: contiguous typed array
//! - Offsets buffer: for variable-length data (strings, nested arrays)
//! - Null bitmap: optional validity buffer

use serde::{Deserialize, Serialize};

/// A columnar buffer of f64 values with shape metadata.
///
/// Designed for zero-copy transfer to JavaScript `Float64Array`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Float64Column {
    /// Column name/label.
    pub name: String,
    /// Flat f64 values.
    pub values: Vec<f64>,
    /// Shape of the data: e.g., \[n_rows\] for 1D, \[n_rows, n_cols\] for 2D.
    pub shape: Vec<usize>,
}

/// A columnar buffer of i32 values.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Int32Column {
    /// Column name/label.
    pub name: String,
    /// Flat i32 values.
    pub values: Vec<i32>,
    /// Shape of the data.
    pub shape: Vec<usize>,
}

/// A columnar buffer of u8 values (for element arrays, flags, etc.).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Uint8Column {
    /// Column name/label.
    pub name: String,
    /// Flat u8 values.
    pub values: Vec<u8>,
    /// Shape of the data.
    pub shape: Vec<usize>,
}

/// A record batch: a named collection of typed columns.
///
/// This is the primary unit for zero-copy data transfer.
/// Analogous to Arrow RecordBatch.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct RecordBatch {
    /// Number of rows in the batch.
    pub num_rows: usize,
    /// Schema: column names and types.
    pub schema: Vec<ColumnSchema>,
    /// Float64 columns.
    pub float_columns: Vec<Float64Column>,
    /// Int32 columns.
    pub int_columns: Vec<Int32Column>,
    /// Uint8 columns.
    pub uint8_columns: Vec<Uint8Column>,
}

/// Schema entry for a column.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ColumnSchema {
    pub name: String,
    pub dtype: DataType,
    pub shape: Vec<usize>,
}

/// Supported data types.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
pub enum DataType {
    Float64,
    Int32,
    Uint8,
}

impl RecordBatch {
    /// Create an empty batch.
    pub fn new() -> Self {
        Self {
            num_rows: 0,
            schema: Vec::new(),
            float_columns: Vec::new(),
            int_columns: Vec::new(),
            uint8_columns: Vec::new(),
        }
    }

    /// Add a Float64 column.
    pub fn add_float64(&mut self, name: &str, values: Vec<f64>, shape: Vec<usize>) {
        if !shape.is_empty() {
            self.num_rows = shape[0];
        }
        self.schema.push(ColumnSchema {
            name: name.to_string(),
            dtype: DataType::Float64,
            shape: shape.clone(),
        });
        self.float_columns.push(Float64Column {
            name: name.to_string(),
            values,
            shape,
        });
    }

    /// Add an Int32 column.
    pub fn add_int32(&mut self, name: &str, values: Vec<i32>, shape: Vec<usize>) {
        if !shape.is_empty() {
            self.num_rows = shape[0];
        }
        self.schema.push(ColumnSchema {
            name: name.to_string(),
            dtype: DataType::Int32,
            shape: shape.clone(),
        });
        self.int_columns.push(Int32Column {
            name: name.to_string(),
            values,
            shape,
        });
    }

    /// Add a Uint8 column.
    pub fn add_uint8(&mut self, name: &str, values: Vec<u8>, shape: Vec<usize>) {
        if !shape.is_empty() {
            self.num_rows = shape[0];
        }
        self.schema.push(ColumnSchema {
            name: name.to_string(),
            dtype: DataType::Uint8,
            shape: shape.clone(),
        });
        self.uint8_columns.push(Uint8Column {
            name: name.to_string(),
            values,
            shape,
        });
    }

    /// Total byte size of all buffers (for transfer cost estimation).
    pub fn byte_size(&self) -> usize {
        let f64_bytes: usize = self
            .float_columns
            .iter()
            .map(|c| c.values.len() * 8)
            .sum();
        let i32_bytes: usize = self.int_columns.iter().map(|c| c.values.len() * 4).sum();
        let u8_bytes: usize = self.uint8_columns.iter().map(|c| c.values.len()).sum();
        f64_bytes + i32_bytes + u8_bytes
    }

    /// Number of columns.
    pub fn num_columns(&self) -> usize {
        self.schema.len()
    }
}

/// Pack conformer results into an Arrow-compatible batch.
///
/// Columns: elements (u8), coords_x/y/z (f64), success (i32), time_ms (f64)
pub fn pack_conformers(results: &[crate::ConformerResult]) -> RecordBatch {
    let n = results.len();
    let mut batch = RecordBatch::new();

    // Success flag (1=ok, 0=failed)
    let success: Vec<i32> = results
        .iter()
        .map(|r| if r.error.is_none() { 1 } else { 0 })
        .collect();
    batch.add_int32("success", success, vec![n]);

    // Number of atoms per molecule
    let num_atoms: Vec<i32> = results.iter().map(|r| r.num_atoms as i32).collect();
    batch.add_int32("num_atoms", num_atoms, vec![n]);

    // Timing
    let times: Vec<f64> = results.iter().map(|r| r.time_ms).collect();
    batch.add_float64("time_ms", times, vec![n]);

    // All coordinates flattened (for successful molecules)
    let all_coords: Vec<f64> = results
        .iter()
        .flat_map(|r| r.coords.iter().copied())
        .collect();
    let total_atoms: usize = results.iter().map(|r| r.num_atoms).sum();
    batch.add_float64("coords", all_coords, vec![total_atoms, 3]);

    // All elements flattened
    let all_elements: Vec<u8> = results
        .iter()
        .flat_map(|r| r.elements.iter().copied())
        .collect();
    batch.add_uint8("elements", all_elements, vec![total_atoms]);

    batch
}

/// Pack an ESP grid into an Arrow-compatible batch.
pub fn pack_esp_grid(grid: &crate::esp::EspGrid) -> RecordBatch {
    let mut batch = RecordBatch::new();

    batch.add_float64(
        "values",
        grid.values.clone(),
        vec![grid.dims[0], grid.dims[1], grid.dims[2]],
    );
    batch.add_float64("origin", grid.origin.to_vec(), vec![3]);
    batch.add_int32(
        "dims",
        grid.dims.iter().map(|&d| d as i32).collect(),
        vec![3],
    );

    batch
}

/// Pack DOS data into an Arrow-compatible batch.
pub fn pack_dos(dos: &crate::dos::DosResult) -> RecordBatch {
    let n_points = dos.energies.len();
    let n_atoms = dos.pdos.len();
    let mut batch = RecordBatch::new();

    batch.add_float64("energies", dos.energies.clone(), vec![n_points]);
    batch.add_float64("total_dos", dos.total_dos.clone(), vec![n_points]);

    // PDOS flattened: [atom0_point0, atom0_point1, ..., atom1_point0, ...]
    let flat_pdos: Vec<f64> = dos.pdos.iter().flat_map(|p| p.iter().copied()).collect();
    batch.add_float64("pdos", flat_pdos, vec![n_atoms, n_points]);

    batch
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_record_batch_empty() {
        let batch = RecordBatch::new();
        assert_eq!(batch.num_rows, 0);
        assert_eq!(batch.num_columns(), 0);
        assert_eq!(batch.byte_size(), 0);
    }

    #[test]
    fn test_record_batch_add_columns() {
        let mut batch = RecordBatch::new();
        batch.add_float64("x", vec![1.0, 2.0, 3.0], vec![3]);
        batch.add_int32("id", vec![0, 1, 2], vec![3]);
        batch.add_uint8("flags", vec![1, 0, 1], vec![3]);

        assert_eq!(batch.num_rows, 3);
        assert_eq!(batch.num_columns(), 3);
        assert_eq!(batch.byte_size(), 3 * 8 + 3 * 4 + 3); // 24 + 12 + 3 = 39
    }

    #[test]
    fn test_pack_conformers() {
        let results = vec![crate::ConformerResult {
            smiles: "C".to_string(),
            num_atoms: 5,
            coords: vec![0.0; 15],
            elements: vec![6, 1, 1, 1, 1],
            bonds: vec![],
            error: None,
            time_ms: 1.5,
        }];
        let batch = pack_conformers(&results);
        assert_eq!(batch.num_columns(), 5);
        assert!(batch.byte_size() > 0);
    }

    #[test]
    fn test_pack_esp_grid() {
        let grid = crate::esp::EspGrid {
            origin: [0.0, 0.0, 0.0],
            spacing: 0.5,
            dims: [3, 3, 3],
            values: vec![0.1; 27],
        };
        let batch = pack_esp_grid(&grid);
        assert_eq!(batch.float_columns[0].values.len(), 27);
    }

    #[test]
    fn test_column_schema() {
        let mut batch = RecordBatch::new();
        batch.add_float64("coords", vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0], vec![2, 3]);
        assert_eq!(batch.schema[0].dtype, DataType::Float64);
        assert_eq!(batch.schema[0].shape, vec![2, 3]);
    }
}
