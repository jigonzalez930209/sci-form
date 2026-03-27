//! **ALPHA** — Live interactive molecular dynamics WASM bindings.
//!
//! Provides a persistent `LiveSimulation` handle that the browser can step
//! frame-by-frame, query positions via a zero-copy pointer, and apply
//! interactive steering forces (IMD).

use crate::helpers::*;
use wasm_bindgen::prelude::*;

/// Persistent live simulation state accessible from JavaScript.
///
/// The browser creates one instance via `live_simulation_create`, then calls
/// `live_simulation_step` and reads the flat coordinate buffer each
/// animation frame.
#[wasm_bindgen]
pub struct LiveSimulation {
    system: sci_form::dynamics_live::state::LiveMolecularSystem,
}

#[wasm_bindgen]
impl LiveSimulation {
    /// Create a new live simulation from SMILES, flat coordinates, and backend name.
    ///
    /// `backend`: `"uff"` | `"pm3"` | `"xtb"`.
    #[wasm_bindgen(constructor)]
    pub fn new(smiles: &str, coords_flat: &str, backend: &str) -> Result<LiveSimulation, JsValue> {
        let flat: Vec<f64> = parse_flat_coords(coords_flat).map_err(|e| JsValue::from_str(&e))?;

        let md_backend = match backend {
            "pm3" => sci_form::dynamics::MdBackend::Pm3,
            "xtb" => sci_form::dynamics::MdBackend::Xtb,
            _ => sci_form::dynamics::MdBackend::Uff,
        };

        let conf = sci_form::embed(smiles, 42);
        if conf.error.is_some() && flat.is_empty() {
            return Err(JsValue::from_str(
                "Failed to parse SMILES and no coords provided",
            ));
        }

        let mut system = if flat.is_empty() {
            sci_form::dynamics_live::state::LiveMolecularSystem::from_conformer(
                &conf.smiles,
                &conf.elements,
                &conf.coords,
                &conf.bonds,
                md_backend,
            )
        } else {
            let elements: Vec<u8> = if conf.elements.is_empty() {
                match sci_form::parse(smiles) {
                    Ok(mol) => {
                        let n = mol.graph.node_count();
                        (0..n)
                            .map(|i| mol.graph[petgraph::graph::NodeIndex::new(i)].element)
                            .collect()
                    }
                    Err(e) => return Err(JsValue::from_str(&e)),
                }
            } else {
                conf.elements.clone()
            };

            let bonds: Vec<(usize, usize, String)> = if conf.bonds.is_empty() {
                match sci_form::parse(smiles) {
                    Ok(mol) => mol
                        .graph
                        .edge_indices()
                        .map(|e| {
                            let (a, b) = mol.graph.edge_endpoints(e).unwrap();
                            let order = match mol.graph[e].order {
                                sci_form::graph::BondOrder::Single => "SINGLE",
                                sci_form::graph::BondOrder::Double => "DOUBLE",
                                sci_form::graph::BondOrder::Triple => "TRIPLE",
                                sci_form::graph::BondOrder::Aromatic => "AROMATIC",
                                sci_form::graph::BondOrder::Unknown => "UNKNOWN",
                            };
                            (a.index(), b.index(), order.to_string())
                        })
                        .collect(),
                    Err(_) => Vec::new(),
                }
            } else {
                conf.bonds.clone()
            };

            sci_form::dynamics_live::state::LiveMolecularSystem::from_conformer(
                smiles, &elements, &flat, &bonds, md_backend,
            )
        };

        system.initialize_velocities(300.0, 42);

        Ok(LiveSimulation { system })
    }

    /// Advance the simulation by `n_steps` with time step `dt_fs` (femtoseconds).
    pub fn step(&mut self, n_steps: usize, dt_fs: f64) {
        for _ in 0..n_steps {
            self.system.integrate(dt_fs);
        }
    }

    /// Apply a Berendsen thermostat coupling towards `target_temp` K.
    pub fn apply_thermostat(&mut self, target_temp: f64, tau_fs: f64, dt_fs: f64) {
        self.system.berendsen_thermostat(target_temp, tau_fs, dt_fs);
    }

    /// Apply a harmonic steering force pulling atom `atom_index` toward `(tx, ty, tz)`.
    pub fn apply_steering_force(
        &mut self,
        atom_index: usize,
        tx: f64,
        ty: f64,
        tz: f64,
        spring_k: f64,
    ) {
        if atom_index < self.system.atoms.len() {
            let sf = sci_form::dynamics_live::steering::SteeringForce {
                atom_index,
                target_xyz: [tx, ty, tz],
                spring_k,
            };
            let mut forces = sci_form::dynamics_live::steering::SteeringForces::new();
            forces.add(sf);
            forces.apply(&mut self.system.atoms);
        }
    }

    /// Clear all accumulated steering forces (call before next IMD frame).
    pub fn clear_steering(&mut self) {
        // Forces are ephemeral — applied per-step in apply_steering_force.
        // This is a no-op but kept for API clarity.
    }

    /// Get the number of atoms.
    pub fn num_atoms(&self) -> usize {
        self.system.atoms.len()
    }

    /// Return current flat positions as a JSON array.
    pub fn positions_json(&mut self) -> String {
        self.system.sync_flat_from_atoms();
        serde_json::to_string(&self.system.positions_flat).unwrap_or_else(|_| "[]".to_string())
    }

    /// Sync internal buffers and return a pointer to the flat positions array.
    ///
    /// The caller (JS side) can construct a `Float64Array` view over WASM memory
    /// at this pointer, with length = `num_atoms * 3`.
    pub fn positions_ptr(&mut self) -> *const f64 {
        self.system.sync_flat_from_atoms();
        self.system.positions_flat.as_ptr()
    }

    /// Return current kinetic temperature in Kelvin.
    pub fn temperature(&self) -> f64 {
        let (_, temp) = self.system.kinetic_energy_and_temperature();
        temp
    }

    /// Return current potential energy (from last force computation).
    pub fn potential_energy(&self) -> f64 {
        self.system.potential_energy
    }

    /// Return current kinetic energy.
    pub fn kinetic_energy(&self) -> f64 {
        let (ke, _) = self.system.kinetic_energy_and_temperature();
        ke
    }

    /// Return simulation step count.
    pub fn step_count(&self) -> usize {
        self.system.step
    }

    /// Set velocities from a Maxwell-Boltzmann distribution at `temperature` K.
    pub fn set_temperature(&mut self, temperature: f64, seed: u64) {
        self.system.initialize_velocities(temperature, seed);
    }
}
