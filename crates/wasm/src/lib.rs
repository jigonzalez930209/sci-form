use wasm_bindgen::prelude::*;

// ─── Typed-array WASM APIs ───────────────────────────────────────────────────
// These avoid JSON serialization overhead for large numeric data, returning
// Float64Array / Float32Array directly for zero-copy transfer to JavaScript.

/// Generate 3D coordinates as a Float64Array (typed-array transfer).
///
/// Returns a Float64Array [x0,y0,z0, x1,y1,z1,...] or empty on failure.
#[wasm_bindgen]
pub fn embed_coords_typed(smiles: &str, seed: u32) -> js_sys::Float64Array {
    let result = sci_form::embed(smiles, seed as u64);
    if result.error.is_some() || result.coords.is_empty() {
        return js_sys::Float64Array::new_with_length(0);
    }
    let arr = js_sys::Float64Array::new_with_length(result.coords.len() as u32);
    arr.copy_from(&result.coords);
    arr
}

/// Compute ESP grid and return values as Float64Array (typed-array transfer).
///
/// Returns Float64Array of grid values. Grid metadata available via compute_esp_grid_info().
#[wasm_bindgen]
pub fn compute_esp_grid_typed(
    elements: &str,
    coords_flat: &str,
    spacing: f64,
    padding: f64,
) -> js_sys::Float64Array {
    let elems: Vec<u8> = match serde_json::from_str(elements) {
        Ok(v) => v,
        Err(_) => return js_sys::Float64Array::new_with_length(0),
    };
    let flat: Vec<f64> = match serde_json::from_str(coords_flat) {
        Ok(v) => v,
        Err(_) => return js_sys::Float64Array::new_with_length(0),
    };
    let positions: Vec<[f64; 3]> = flat.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect();
    match sci_form::compute_esp(&elems, &positions, spacing, padding) {
        Ok(grid) => {
            let arr = js_sys::Float64Array::new_with_length(grid.values.len() as u32);
            arr.copy_from(&grid.values);
            arr
        }
        Err(_) => js_sys::Float64Array::new_with_length(0),
    }
}

/// Get ESP grid metadata (origin, spacing, dimensions) as JSON.
#[wasm_bindgen]
pub fn compute_esp_grid_info(
    elements: &str,
    coords_flat: &str,
    spacing: f64,
    padding: f64,
) -> String {
    let elems: Vec<u8> = match serde_json::from_str(elements) {
        Ok(v) => v,
        Err(e) => return format!("{{\"error\":\"bad elements: {}\"}}", e),
    };
    let flat: Vec<f64> = match serde_json::from_str(coords_flat) {
        Ok(v) => v,
        Err(e) => return format!("{{\"error\":\"bad coords: {}\"}}", e),
    };
    let positions: Vec<[f64; 3]> = flat.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect();
    match sci_form::compute_esp(&elems, &positions, spacing, padding) {
        Ok(grid) => format!(
            "{{\"origin\":[{:.4},{:.4},{:.4}],\"spacing\":{:.4},\"dims\":[{},{},{}]}}",
            grid.origin[0],
            grid.origin[1],
            grid.origin[2],
            grid.spacing,
            grid.dims[0],
            grid.dims[1],
            grid.dims[2]
        ),
        Err(e) => format!("{{\"error\":\"{}\"}}", e),
    }
}

/// Compute orbital grid as Float32Array (typed-array transfer).
///
/// Returns the flat volumetric grid as Float32Array for GPU volume rendering.
#[wasm_bindgen]
pub fn eht_orbital_grid_typed(
    elements: &str,
    coords_flat: &str,
    mo_index: usize,
    spacing: f64,
) -> js_sys::Float32Array {
    let elems: Vec<u8> = match serde_json::from_str(elements) {
        Ok(v) => v,
        Err(_) => return js_sys::Float32Array::new_with_length(0),
    };
    let flat: Vec<f64> = match serde_json::from_str(coords_flat) {
        Ok(v) => v,
        Err(_) => return js_sys::Float32Array::new_with_length(0),
    };
    let positions: Vec<[f64; 3]> = flat.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect();

    let result = match sci_form::eht::solve_eht(&elems, &positions, None) {
        Ok(r) => r,
        Err(_) => return js_sys::Float32Array::new_with_length(0),
    };
    let basis = sci_form::eht::basis::build_basis(&elems, &positions);
    let grid = sci_form::eht::evaluate_orbital_on_grid(
        &basis,
        &result.coefficients,
        mo_index,
        &positions,
        spacing,
        3.0,
    );
    let arr = js_sys::Float32Array::new_with_length(grid.values.len() as u32);
    arr.copy_from(&grid.values);
    arr
}

/// Library version.
#[wasm_bindgen]
pub fn version() -> String {
    sci_form::version()
}

/// Generate a 3D conformer from a SMILES string.
///
/// Returns a JSON string with the full result including atoms, bonds, and coordinates.
#[wasm_bindgen]
pub fn embed(smiles: &str, seed: u32) -> String {
    let result = sci_form::embed(smiles, seed as u64);
    serde_json::to_string(&result).unwrap_or_else(|e| format!("{{\"error\":\"{}\"}}", e))
}

/// Generate 3D coordinates only (compact format).
///
/// Returns JSON: {"coords": [x0,y0,z0,...], "num_atoms": N} or {"error": "..."}
#[wasm_bindgen]
pub fn embed_coords(smiles: &str, seed: u32) -> String {
    let result = sci_form::embed(smiles, seed as u64);
    if let Some(ref e) = result.error {
        return format!("{{\"error\":\"{}\"}}", e.replace('"', "\\'"));
    }
    format!(
        "{{\"coords\":[{}],\"num_atoms\":{}}}",
        result
            .coords
            .iter()
            .map(|v| format!("{:.4}", v))
            .collect::<Vec<_>>()
            .join(","),
        result.num_atoms
    )
}

/// Batch-embed multiple molecules from newline-separated SMILES.
///
/// Returns a JSON array of results.
#[wasm_bindgen]
pub fn embed_batch(smiles_list: &str, seed: u32) -> String {
    let lines: Vec<&str> = smiles_list
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| l.trim())
        .collect();
    let config = sci_form::ConformerConfig {
        seed: seed as u64,
        num_threads: 0,
    };
    let results = sci_form::embed_batch(&lines, &config);
    serde_json::to_string(&results).unwrap_or_else(|e| format!("[{{\"error\":\"{}\"}}]", e))
}

/// Parse a SMILES string and return molecular info (no 3D).
#[wasm_bindgen]
pub fn parse_smiles(smiles: &str) -> String {
    match sci_form::parse(smiles) {
        Ok(mol) => {
            let n = mol.graph.node_count();
            let nb = mol.graph.edge_count();
            format!("{{\"num_atoms\":{},\"num_bonds\":{}}}", n, nb)
        }
        Err(e) => format!("{{\"error\":\"{}\"}}", e.replace('"', "\\'")),
    }
}

// ─── EHT (Extended Hückel Theory) API ────────────────────────────────────────

/// Run an EHT calculation on a molecule given its elements and flat XYZ coordinates.
///
/// - `elements`: JSON array of atomic numbers, e.g. `[8,1,1]`
/// - `coords_flat`: JSON array of flat xyz coordinates `[x0,y0,z0,x1,y1,z1,...]`
/// - `k`: Wolfsberg-Helmholtz constant (pass 0.0 for default 1.75)
///
/// Returns JSON with energies, coefficients, HOMO/LUMO info.
#[wasm_bindgen]
pub fn eht_calculate(elements: &str, coords_flat: &str, k: f64) -> String {
    let elems: Vec<u8> = match serde_json::from_str(elements) {
        Ok(v) => v,
        Err(e) => return format!("{{\"error\":\"bad elements: {}\"}}", e),
    };
    let flat: Vec<f64> = match serde_json::from_str(coords_flat) {
        Ok(v) => v,
        Err(e) => return format!("{{\"error\":\"bad coords: {}\"}}", e),
    };
    if flat.len() != elems.len() * 3 {
        return format!(
            "{{\"error\":\"coords length {} != 3 * elements {}\"}}",
            flat.len(),
            elems.len()
        );
    }
    let positions: Vec<[f64; 3]> = flat.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect();
    let k_opt = if k <= 0.0 { None } else { Some(k) };

    match sci_form::eht::solve_eht(&elems, &positions, k_opt) {
        Ok(result) => {
            serde_json::to_string(&result).unwrap_or_else(|e| format!("{{\"error\":\"{}\"}}", e))
        }
        Err(e) => format!("{{\"error\":\"{}\"}}", e),
    }
}

/// Generate an orbital volumetric grid and extract an isosurface mesh.
///
/// - `elements`, `coords_flat`: same as `eht_calculate()`
/// - `mo_index`: molecular orbital index (0-based)
/// - `spacing`: grid spacing in Ångström (e.g. 0.2)
/// - `isovalue`: isosurface cutoff (e.g. 0.02)
///
/// Returns JSON with mesh vertices, normals, indices.
#[wasm_bindgen]
pub fn eht_orbital_mesh(
    elements: &str,
    coords_flat: &str,
    mo_index: usize,
    spacing: f64,
    isovalue: f32,
) -> String {
    let elems: Vec<u8> = match serde_json::from_str(elements) {
        Ok(v) => v,
        Err(e) => return format!("{{\"error\":\"bad elements: {}\"}}", e),
    };
    let flat: Vec<f64> = match serde_json::from_str(coords_flat) {
        Ok(v) => v,
        Err(e) => return format!("{{\"error\":\"bad coords: {}\"}}", e),
    };
    if flat.len() != elems.len() * 3 {
        return format!("{{\"error\":\"coords length mismatch\"}}");
    }
    let positions: Vec<[f64; 3]> = flat.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect();

    let result = match sci_form::eht::solve_eht(&elems, &positions, None) {
        Ok(r) => r,
        Err(e) => return format!("{{\"error\":\"{}\"}}", e),
    };

    let basis = sci_form::eht::basis::build_basis(&elems, &positions);
    let grid = sci_form::eht::evaluate_orbital_on_grid(
        &basis,
        &result.coefficients,
        mo_index,
        &positions,
        spacing,
        3.0,
    );
    let mesh = sci_form::eht::marching_cubes(&grid, isovalue);
    serde_json::to_string(&mesh).unwrap_or_else(|e| format!("{{\"error\":\"{}\"}}", e))
}

// ─── Charges (Gasteiger-Marsili) API ─────────────────────────────────────────

/// Compute Gasteiger-Marsili partial charges from a SMILES string.
///
/// Returns JSON with charges, iterations, total_charge.
#[wasm_bindgen]
pub fn compute_charges(smiles: &str) -> String {
    match sci_form::compute_charges(smiles) {
        Ok(result) => {
            serde_json::to_string(&result).unwrap_or_else(|e| format!("{{\"error\":\"{}\"}}", e))
        }
        Err(e) => format!("{{\"error\":\"{}\"}}", e),
    }
}

// ─── SASA (Shrake-Rupley) API ────────────────────────────────────────────────

/// Compute solvent-accessible surface area from atomic numbers and coordinates.
///
/// - `elements`: JSON array of atomic numbers, e.g. `[8,1,1]`
/// - `coords_flat`: JSON array of flat xyz coordinates
/// - `probe_radius`: probe radius in Å (pass 0.0 for default 1.4)
///
/// Returns JSON with total_sasa, atom_sasa, probe_radius, num_points.
#[wasm_bindgen]
pub fn compute_sasa(elements: &str, coords_flat: &str, probe_radius: f64) -> String {
    let elems: Vec<u8> = match serde_json::from_str(elements) {
        Ok(v) => v,
        Err(e) => return format!("{{\"error\":\"bad elements: {}\"}}", e),
    };
    let flat: Vec<f64> = match serde_json::from_str(coords_flat) {
        Ok(v) => v,
        Err(e) => return format!("{{\"error\":\"bad coords: {}\"}}", e),
    };
    let pr = if probe_radius <= 0.0 {
        None
    } else {
        Some(probe_radius)
    };
    match sci_form::compute_sasa(&elems, &flat, pr) {
        Ok(result) => {
            serde_json::to_string(&result).unwrap_or_else(|e| format!("{{\"error\":\"{}\"}}", e))
        }
        Err(e) => format!("{{\"error\":\"{}\"}}", e),
    }
}

/// Compute Mulliken & Löwdin population analysis.
/// Elements: JSON array of atomic numbers. Coords: JSON array of flat xyz.
#[wasm_bindgen]
pub fn compute_population(elements: &str, coords_flat: &str) -> String {
    let elems: Vec<u8> = match serde_json::from_str(elements) {
        Ok(v) => v,
        Err(e) => return format!("{{\"error\":\"bad elements: {}\"}}", e),
    };
    let flat: Vec<f64> = match serde_json::from_str(coords_flat) {
        Ok(v) => v,
        Err(e) => return format!("{{\"error\":\"bad coords: {}\"}}", e),
    };
    let positions: Vec<[f64; 3]> = flat.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect();
    match sci_form::compute_population(&elems, &positions) {
        Ok(result) => {
            serde_json::to_string(&result).unwrap_or_else(|e| format!("{{\"error\":\"{}\"}}", e))
        }
        Err(e) => format!("{{\"error\":\"{}\"}}", e),
    }
}

/// Compute molecular dipole moment.
#[wasm_bindgen]
pub fn compute_dipole(elements: &str, coords_flat: &str) -> String {
    let elems: Vec<u8> = match serde_json::from_str(elements) {
        Ok(v) => v,
        Err(e) => return format!("{{\"error\":\"bad elements: {}\"}}", e),
    };
    let flat: Vec<f64> = match serde_json::from_str(coords_flat) {
        Ok(v) => v,
        Err(e) => return format!("{{\"error\":\"bad coords: {}\"}}", e),
    };
    let positions: Vec<[f64; 3]> = flat.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect();
    match sci_form::compute_dipole(&elems, &positions) {
        Ok(result) => {
            serde_json::to_string(&result).unwrap_or_else(|e| format!("{{\"error\":\"{}\"}}", e))
        }
        Err(e) => format!("{{\"error\":\"{}\"}}", e),
    }
}

/// Compute DOS/PDOS from EHT orbital energies.
#[wasm_bindgen]
pub fn compute_dos(
    elements: &str,
    coords_flat: &str,
    sigma: f64,
    e_min: f64,
    e_max: f64,
    n_points: usize,
) -> String {
    let elems: Vec<u8> = match serde_json::from_str(elements) {
        Ok(v) => v,
        Err(e) => return format!("{{\"error\":\"bad elements: {}\"}}", e),
    };
    let flat: Vec<f64> = match serde_json::from_str(coords_flat) {
        Ok(v) => v,
        Err(e) => return format!("{{\"error\":\"bad coords: {}\"}}", e),
    };
    let positions: Vec<[f64; 3]> = flat.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect();
    match sci_form::compute_dos(&elems, &positions, sigma, e_min, e_max, n_points) {
        Ok(result) => {
            let energies_json: Vec<String> = result
                .energies
                .iter()
                .map(|v| format!("{:.4}", v))
                .collect();
            let dos_json: Vec<String> = result
                .total_dos
                .iter()
                .map(|v| format!("{:.6}", v))
                .collect();
            format!(
                "{{\"sigma\":{},\"energies\":[{}],\"total_dos\":[{}],\"n_atoms_pdos\":{}}}",
                result.sigma,
                energies_json.join(","),
                dos_json.join(","),
                result.pdos.len()
            )
        }
        Err(e) => format!("{{\"error\":\"{}\"}}", e),
    }
}

/// Compute RMSD between two coordinate sets after Kabsch alignment.
#[wasm_bindgen]
pub fn compute_rmsd(coords: &str, reference: &str) -> String {
    let c: Vec<f64> = match serde_json::from_str(coords) {
        Ok(v) => v,
        Err(e) => return format!("{{\"error\":\"bad coords: {}\"}}", e),
    };
    let r: Vec<f64> = match serde_json::from_str(reference) {
        Ok(v) => v,
        Err(e) => return format!("{{\"error\":\"bad reference: {}\"}}", e),
    };
    let result = sci_form::alignment::align_coordinates(&c, &r);
    format!("{{\"rmsd\":{:.6}}}", result.rmsd)
}

/// Compute UFF force field energy.
#[wasm_bindgen]
pub fn compute_uff_energy(smiles: &str, coords: &str) -> String {
    let flat: Vec<f64> = match serde_json::from_str(coords) {
        Ok(v) => v,
        Err(e) => return format!("{{\"error\":\"bad coords: {}\"}}", e),
    };
    match sci_form::compute_uff_energy(smiles, &flat) {
        Ok(energy) => format!("{{\"energy\":{:.6},\"unit\":\"kcal/mol\"}}", energy),
        Err(e) => format!("{{\"error\":\"{}\"}}", e),
    }
}

/// Create a unit cell and return parameters + volume as JSON.
#[wasm_bindgen]
pub fn create_unit_cell(a: f64, b: f64, c: f64, alpha: f64, beta: f64, gamma: f64) -> String {
    let cell = sci_form::create_unit_cell(a, b, c, alpha, beta, gamma);
    let vol = cell.volume();
    let p = cell.parameters();
    format!(
        "{{\"a\":{:.4},\"b\":{:.4},\"c\":{:.4},\"alpha\":{:.2},\"beta\":{:.2},\"gamma\":{:.2},\"volume\":{:.4}}}",
        p.a, p.b, p.c, p.alpha, p.beta, p.gamma, vol
    )
}

/// Assemble a framework crystal structure.
///
/// `topology`: "pcu", "dia", or "sql"
/// `metal`: atomic number of the metal center
/// `geometry`: "linear", "tetrahedral", "octahedral", "square_planar", "trigonal"
/// `lattice_a`: cubic lattice parameter in Å
/// `supercell`: replication factor (1 = no replication)
#[wasm_bindgen]
pub fn assemble_framework(
    topology: &str,
    metal: u8,
    geometry: &str,
    lattice_a: f64,
    supercell: usize,
) -> String {
    let geom = match geometry {
        "linear" => sci_form::materials::CoordinationGeometry::Linear,
        "trigonal" => sci_form::materials::CoordinationGeometry::Trigonal,
        "tetrahedral" => sci_form::materials::CoordinationGeometry::Tetrahedral,
        "square_planar" => sci_form::materials::CoordinationGeometry::SquarePlanar,
        "octahedral" => sci_form::materials::CoordinationGeometry::Octahedral,
        _ => return format!("{{\"error\":\"unknown geometry: {}\"}}", geometry),
    };
    let topo = match topology {
        "pcu" => sci_form::materials::Topology::pcu(),
        "dia" => sci_form::materials::Topology::dia(),
        "sql" => sci_form::materials::Topology::sql(),
        _ => return format!("{{\"error\":\"unknown topology: {}\"}}", topology),
    };
    let node = sci_form::materials::Sbu::metal_node(metal, 0.0, geom);
    let linker = sci_form::materials::Sbu::linear_linker(&[6, 6], 1.4, "carboxylate");
    let cell = sci_form::materials::UnitCell::cubic(lattice_a);
    let mut structure = sci_form::assemble_framework(&node, &linker, &topo, &cell);
    if supercell > 1 {
        structure = structure.make_supercell(supercell, supercell, supercell);
    }
    serde_json::to_string(&structure).unwrap_or_else(|e| format!("{{\"error\":\"{}\"}}", e))
}

// ─── C9: Transport / Streaming ────────────────────────────────────────────

/// Pack a batch of conformer results into Arrow-compatible columnar format.
///
/// Input: JSON array of ConformerResult objects.
/// Output: JSON RecordBatch with typed columns for zero-copy transfer.
#[wasm_bindgen]
pub fn pack_batch_arrow(results_json: &str) -> String {
    let results: Vec<sci_form::ConformerResult> = match serde_json::from_str(results_json) {
        Ok(r) => r,
        Err(e) => return format!("{{\"error\":\"bad JSON: {}\"}}", e),
    };
    let batch = sci_form::transport::pack_conformers(&results);
    serde_json::to_string(&batch).unwrap_or_else(|e| format!("{{\"error\":\"{}\"}}", e))
}

/// Split a batch of SMILES into worker tasks for Web Worker dispatch.
///
/// Input: JSON array of SMILES strings.
/// Output: JSON array of WorkerTask objects.
#[wasm_bindgen]
pub fn split_worker_tasks(smiles_json: &str, n_workers: usize, seed: u32) -> String {
    let smiles: Vec<String> = match serde_json::from_str(smiles_json) {
        Ok(s) => s,
        Err(e) => return format!("{{\"error\":\"bad JSON: {}\"}}", e),
    };
    let tasks = sci_form::transport::split_batch(&smiles, n_workers, seed as u64);
    serde_json::to_string(&tasks).unwrap_or_else(|e| format!("{{\"error\":\"{}\"}}", e))
}

/// Estimate optimal number of Web Workers for a batch size.
#[wasm_bindgen]
pub fn estimate_workers(n_items: usize, max_workers: usize) -> usize {
    sci_form::transport::estimate_workers(n_items, max_workers)
}
