/**
 * **ALPHA** — requestAnimationFrame simulation loop for live molecular dynamics.
 *
 * Drives a LiveSimulationView at a target physics rate while rendering
 * at display refresh rate. Supports pause/resume and per-frame callbacks.
 */

import { LiveSimulationView } from './live_simulation';

export interface SimulationLoopConfig {
  /** Physics steps per animation frame (default: 10). */
  stepsPerFrame?: number;
  /** Time step in femtoseconds (default: 1.0). */
  dtFs?: number;
  /** Target temperature for Berendsen thermostat (K). Omit to disable thermostat. */
  targetTempK?: number;
  /** Thermostat coupling time (fs, default: 100). */
  thermostatTauFs?: number;
  /** Called each frame with the current positions Float64Array. */
  onFrame?: (positions: Float64Array, frameInfo: FrameInfo) => void;
}

export interface FrameInfo {
  /** Simulation step count. */
  step: number;
  /** Current temperature (K). */
  temperature: number;
  /** Potential energy (eV). */
  potentialEnergy: number;
  /** Kinetic energy (eV). */
  kineticEnergy: number;
  /** Real wall-clock time since loop start (ms). */
  wallTimeMs: number;
}

/**
 * Manage a requestAnimationFrame loop for live molecular dynamics.
 */
export class SimulationLoop {
  private readonly view: LiveSimulationView;
  private readonly config: Required<Pick<SimulationLoopConfig, 'stepsPerFrame' | 'dtFs' | 'thermostatTauFs'>>;
  private readonly targetTempK: number | null;
  private readonly onFrame: SimulationLoopConfig['onFrame'];

  private rafId: number | null = null;
  private running = false;
  private startTime = 0;

  constructor(view: LiveSimulationView, config: SimulationLoopConfig = {}) {
    this.view = view;
    this.config = {
      stepsPerFrame: config.stepsPerFrame ?? 10,
      dtFs: config.dtFs ?? 1.0,
      thermostatTauFs: config.thermostatTauFs ?? 100.0,
    };
    this.targetTempK = config.targetTempK ?? null;
    this.onFrame = config.onFrame;
  }

  /** Start the animation loop. */
  start(): void {
    if (this.running) return;
    this.running = true;
    this.startTime = performance.now();
    this.tick();
  }

  /** Pause the animation loop. */
  pause(): void {
    this.running = false;
    if (this.rafId !== null) {
      cancelAnimationFrame(this.rafId);
      this.rafId = null;
    }
  }

  /** Resume after pause. */
  resume(): void {
    this.start();
  }

  /** Whether the loop is currently running. */
  get isRunning(): boolean {
    return this.running;
  }

  /** Stop the loop permanently. */
  stop(): void {
    this.pause();
  }

  private tick = (): void => {
    if (!this.running) return;

    // Physics step
    const positions = this.view.step(this.config.stepsPerFrame, this.config.dtFs);

    // Thermostat
    if (this.targetTempK !== null) {
      this.view.applyThermostat(this.targetTempK, this.config.thermostatTauFs, this.config.dtFs);
    }

    // Callback
    if (this.onFrame) {
      this.onFrame(positions, {
        step: this.view.stepCount,
        temperature: this.view.temperature,
        potentialEnergy: this.view.potentialEnergy,
        kineticEnergy: this.view.kineticEnergy,
        wallTimeMs: performance.now() - this.startTime,
      });
    }

    this.rafId = requestAnimationFrame(this.tick);
  };
}
