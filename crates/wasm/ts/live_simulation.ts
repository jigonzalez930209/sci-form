/**
 * **ALPHA** — Zero-copy live simulation wrapper for sci-form WASM.
 *
 * Reads the internal Float64Array positions buffer directly from WASM memory
 * instead of JSON-parsing each frame. Designed for requestAnimationFrame loops.
 */

import type { LiveSimulation } from '../pkg/sci_form_wasm';

/**
 * Wrap a WASM LiveSimulation and expose a fast Float64Array view
 * over WASM linear memory for the position buffer.
 */
export class LiveSimulationView {
  private readonly sim: LiveSimulation;
  private readonly wasmMemory: WebAssembly.Memory;
  private readonly nAtoms: number;

  constructor(sim: LiveSimulation, wasmMemory: WebAssembly.Memory) {
    this.sim = sim;
    this.wasmMemory = wasmMemory;
    this.nAtoms = sim.num_atoms();
  }

  /** Number of atoms. */
  get numAtoms(): number {
    return this.nAtoms;
  }

  /**
   * Step the simulation forward and return a Float64Array view of current positions.
   *
   * The returned view reads directly from WASM memory — no JSON parsing or copying.
   * View is invalidated on next WASM memory growth; this method rebuilds it if needed.
   */
  step(nSteps: number, dtFs: number): Float64Array {
    this.sim.step(nSteps, dtFs);
    return this.getPositions();
  }

  /** Get current positions as Float64Array (zero-copy). */
  getPositions(): Float64Array {
    const ptr = this.sim.positions_ptr();
    // WASM memory may grow; rebuild view if buffer changes
    return new Float64Array(this.wasmMemory.buffer, ptr, this.nAtoms * 3);
  }

  /** Apply Berendsen thermostat. */
  applyThermostat(targetTempK: number, tauFs: number, dtFs: number): void {
    this.sim.apply_thermostat(targetTempK, tauFs, dtFs);
  }

  /** Apply interactive steering force on one atom. */
  applySteering(atomIndex: number, target: [number, number, number], springK: number): void {
    this.sim.apply_steering_force(atomIndex, target[0], target[1], target[2], springK);
  }

  /** Current kinetic temperature (K). */
  get temperature(): number {
    return this.sim.temperature();
  }

  /** Current potential energy (eV). */
  get potentialEnergy(): number {
    return this.sim.potential_energy();
  }

  /** Current kinetic energy (eV). */
  get kineticEnergy(): number {
    return this.sim.kinetic_energy();
  }

  /** Current step count. */
  get stepCount(): number {
    return this.sim.step_count();
  }

  /** Set velocities from Maxwell-Boltzmann distribution. */
  setTemperature(tempK: number, seed: number): void {
    this.sim.set_temperature(tempK, seed);
  }
}
