pub mod graph;
pub mod distgeom;
pub mod etkdg;
pub mod forcefield;
pub mod smiles;

use wasm_bindgen::prelude::*;

#[wasm_bindgen]
pub fn greet() -> String {
    "Hello from sci-form WASM!".into()
}
