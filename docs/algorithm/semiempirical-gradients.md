# PM3 and xTB: energy surfaces, SCF and gradients

This page documents what the semi-empirical entry points of sci-form
evaluate, how the PM3 self-consistent field converges, and how the energy
gradients are obtained. It exists because three different things used to be
called "xTB" in the code base and two gradients did not belong to the energy
they were attached to.

> The reference, measured deviation, confidence region and limitations for this
> method are recorded in [Scientific validation](/validation/scientific-validation).

## 1. Which energy surface does a name evaluate?

<SvgDiagram src="/svg/energy-surface-map.svg" alt="Every public method name maps to one energy surface and one gradient path." />

| Public name | Where | Energy surface | Gradient |
|---|---|---|---|
| `pm3`, `compute_pm3`, `pm3_calculate` | all four surfaces | MOPAC-parity PM3 (`pm3::faithful`) when every element has `Z ≤ 17` and the SCF accepts the molecule; the legacy NDDO path for `Z > 17` and as a **recorded** fallback (`notes`) | `compute_pm3_gradient` — stationary-density (Pulay-free) gradient: fixed-density integral rebuild, no SCF per displacement |
| `xtb`, `compute_xtb`, `xtb_calculate` | all four surfaces | **GFN2-xTB** (`solve_gfn2`) repackaged as `XtbResult`; numerically identical to `gfn2` | `compute_xtb_gradient` — warm-started central differences (same energy, fewer SCC iterations) |
| `gfn1`, `gfn2` | all four surfaces | `solve_gfn1`, `solve_gfn2` | none public; NEB/MD use finite differences |
| `"xtb"` (legacy spellings `"gfn0"`, `"gfn0-xtb"`) in **NEB, MD and reaction dynamics** | `NebBackend::Gfn2`, `MdBackend::Xtb` | **GFN2** (`solve_gfn2`) — the same surface as `compute_xtb` | central differences of that energy |
| `"xtb"` / `"gfn0"` in **GSM** | `GsmEnergyBackend::Xtb` | `compute_xtb`, i.e. **GFN2** | — (GSM evaluates energies) |
| `xtb::solve_gfn0` | Rust only | faithful GFN0 (xtb 6.7.1, validation F43), takes a `charge` argument, closed shell only | none |
| IR/Hessian `"xtb"` | `HessianMethod::Xtb` | GFN2 | energy differences |

Two consequences worth knowing:

* `sci-form xtb` is **GFN2-xTB**, not GFN0. Its output equals `sci-form gfn2`
  to the last digit. Earlier text in the documentation, the CLI `--help` and
  the README called it GFN0; that was wrong.
* `"xtb"` now means **GFN2 on every surface** (NEB, MD, reaction dynamics, GSM,
  `compute_xtb`, IR, DOS, mesh). Before v1.0 NEB/MD used a different,
  simplified model (`xtb::solve_xtb`) under the same name. The `"gfn0"`
  spellings are kept only as legacy aliases of the same GFN2 surface; the real
  GFN0 is `xtb::solve_gfn0` (Rust only).

### Units

| Quantity | Unit |
|---|---|
| Coordinates | Å |
| PM3 total, electronic, nuclear-repulsion and orbital energies | eV |
| PM3 heat of formation | kcal/mol |
| xTB / GFN1 / GFN2 total, electronic, repulsive, dispersion and orbital energies | eV |
| PM3 and xTB gradients | eV/Å |
| NEB / MD energies and forces | kcal/mol and kcal/mol/Å |
| KS-DFT energies and gradients | Hartree and Hartree/Bohr |

## 2. The PM3 SCF

NDDO methods work in an orthonormal basis (`S = 1`): the overlap matrix only
enters the core Hamiltonian through the resonance integral
`H_μν = ½(β_μ + β_ν) S_μν`. The SCF therefore solves the ordinary eigenproblem
`F C = C ε` and builds the density `P = 2 Σ_occ C Cᵀ`.

The MOPAC-parity solver converges it in two stages:

1. **Pulay DIIS** (`diis_scf`). The error vector is the commutator
   `e = F P − P F` (the `S = 1` form of `F P S − S P F`), the history holds the
   last 8 Fock matrices, the extrapolation solves the usual
   `(m+1)×(m+1)` system, and the loop stops when `‖e‖∞ < 1·10⁻⁹` (at most 120
   iterations). A singular DIIS system restarts the history; a non-finite
   error abandons DIIS.
2. **Damped fixed point** (fallback). Density mixing with an adaptive factor
   `α ∈ [0.02, 0.9]` (halved when the electronic energy rises, grown by 5 %
   otherwise), at most 200 iterations. It converges when
   `|ΔE| < 1·10⁻⁹ eV` **and** the density residual `max|P_new − P| < 1·10⁻⁷`;
   if the budget runs out while the energy has stalled, the result is still
   accepted as converged (the behaviour before DIIS existed).

<SvgDiagram src="/svg/pm3-scf-diis.svg" alt="Pulay DIIS converges the PM3 self-consistent field in 27 iterations while damped fixed-point mixing exhausts 200." />

On the H₂O₂ gate geometry the damped loop used all 200 iterations, DIIS needs
27, and the two reach the same energy to 1·10⁻⁷ eV. This is not cosmetic. The
energy-only stopping rule left roughly 2·10⁻³ eV/Å of scatter between
finite-difference steps, because the energy had stalled while the density was
still moving; the finite-difference gradient amplified that by `1/(2h)`. With
DIIS the scatter is about 6·10⁻⁶ eV/Å for steps up to 5·10⁻⁴ Å.

`Pm3Result.converged`, `scf_iterations` and `notes` report what happened:
`notes` says whether the MOPAC-parity or the legacy solver produced the
numbers, and why the fallback was taken.

## 3. Gradients

### PM3: stationary-density (Pulay-free) gradient

At a converged SCF density the electronic energy is stationary with respect to
the density, so the Pulay terms vanish. The molecular gradient is the
derivative of

```
E(P, R) = ½ Tr[P (H(R) + G(P, R))] + E_nuc(R)
```

at fixed `P`. `compute_pm3_gradient` therefore:

1. runs `solve_pm3` once to obtain the converged density;
2. for each Cartesian component, rebuilds `H`, `G(P)` and `E_nuc` at
   `R ± h·e_k` **without running an SCF**;
3. forms `(E(P,R+h·e_k) − E(P,R−h·e_k))/2h`, `h = 5·10⁻⁴ Å`.

The remaining integral derivatives are evaluated by central differences of the
fixed-density energy. The expensive SCF is not repeated, and the result
matches the full central difference of `solve_pm3` to better than
`1·10⁻³ eV/Å`. The legacy NDDO path (`Z > 17`) keeps the full central
difference because it does not expose a stationary density.

### xTB: warm-started central differences

GFN2-xTB does not expose a fixed-state energy evaluator, so
`compute_xtb_gradient` still differentiates `solve_gfn2` by central
differences. Each displaced geometry is warm-started from the reference shell
charges and atomic multipoles (`Gfn2Result::warm_start`), so the SCC converges
in a few iterations instead of from a zero state. The energy is identical to a
cold SCC run; only the iteration count drops.

<SvgDiagram src="/svg/gradient-central-difference.svg" alt="A central difference displaces one coordinate by plus and minus h and divides the energy difference by 2h." />

* The `3N` coordinates are independent, so they run on the rayon pool when the
  `parallel` feature is on. Each component is collected in order, so the result
  does not depend on the number of threads.
* Truncation error is `O(h²·E‴) ≈ 5·10⁻⁶ eV/Å`; the SCC is smooth enough (see
  §2) that the noise term `ε/2h` is smaller than that.
* **Cost:** PM3 `6N` integral rebuilds + one SCF; xTB `6N` warm SCC solves +
  one reference solve.

### Why the remaining derivative is still numerical

Until 1.0 the repository carried hand-written Hellmann-Feynman/Pulay
expressions. They differentiated the *legacy* simplified models, not the
energies the API reports (the MOPAC-parity PM3 and GFN2), so they could not
agree with those energies however carefully they were written:

| Gate (`tests/gradient_parity.rs`, H₂O₂) | Before | Now |
|---|---|---|
| PM3 max \|Δg\| against the numerical derivative of `solve_pm3` | 18.4 eV/Å | 1.1·10⁻⁵ eV/Å (tolerance 1·10⁻⁴) |
| xTB max \|Δg\| against the numerical derivative of the xTB energy | 209–263 eV/Å (of a different, simplified model) | 4.8·10⁻⁵ eV/Å (tolerance 1·10⁻⁴) |

The symbolic NDDO derivative (rotated two-centre integrals) and the GFN2
derivative (multipole electrostatics, D4) remain research work. When they land
they will replace the fixed-density / warm-start steps behind the same
function signature.

### Behaviour you should know about

* **Solver switching.** If a displaced geometry makes the MOPAC-parity solver
  fail and the legacy path takes over (see `notes` of `compute_pm3`), the two
  energies belong to different surfaces and the difference is meaningless. In
  practice this only happens for geometries the closed-shell SCF cannot
  represent; the gradient does not expose `notes`, so check `compute_pm3`
  first if a gradient looks implausible.
* **Open-shell and odd-electron systems.** `compute_pm3_gradient` returns the
  SCF error (for example an odd electron count). Use `compute_uhf` / `compute_rohf` for
  open-shell energies; they have no gradient.
* **Errors.** An SCF failure at the reference geometry or at any displaced
  geometry is returned unchanged (`Err(SciFormError)` in Rust; in Python
  `ValueError` for a malformed geometry, `NotConvergedError` for a
  non-convergent SCF and `RuntimeError` for any other engine failure; an
  `{"error": ...}` JSON payload in WASM; a non-zero exit code in the CLI).
* **Translational invariance.** The net force over all atoms is below
  `1·10⁻³ eV/Å` in the unit tests.

## 4. Where the gradients are consumed

| Consumer | Backend `"pm3"` | Backend `"xtb"` |
|---|---|---|
| `neb-gradient`, simplified NEB, MD, reaction dynamics | `compute_pm3_gradient` | finite differences of `solve_gfn2` (§1) |
| `neb-energy` and the final image energies of the simplified NEB | `solve_pm3` — **energy only, one SCF** | `solve_gfn2` — energy only, one SCF |
| `compute_semianalytical_hessian` | second finite difference of the gradient (`36·N²` energies) | the same, on GFN2 |

The semi-analytical Hessian used to be the cheap alternative to
`compute_numerical_hessian` because it differentiated analytic gradients. That
is no longer true, so prefer `compute_numerical_hessian` (or the `ir`
functions, which use it).

NEB, MD and reaction dynamics call the gradient once per image per iteration.
With PM3 the call costs `6N` integral rebuilds plus one SCF; with `xtb` it
costs `6N` warm SCC solves plus one reference solve. Use `uff` or `mmff94`
for exploratory paths and `gfn2` or `pm3` for the final refinement.

## 5. API

### Rust

```rust
use sci_form::{compute_pm3_gradient, compute_xtb_gradient};

let pm3 = compute_pm3_gradient(&elements, &positions)?;   // Pm3GradientResult
let xtb = compute_xtb_gradient(&elements, &positions)?;   // XtbGradientResult
```

| Type | Field | Unit |
|---|---|---|
| `Pm3GradientResult` | `gradients: Vec<[f64; 3]>` | eV/Å, one triple per atom, same order as the input |
| | `energy: f64` | eV (the `solve_pm3` total energy at the input geometry) |
| | `heat_of_formation: f64` | kcal/mol |
| `XtbGradientResult` | `gradients: Vec<[f64; 3]>` | eV/Å |
| | `energy: f64` | eV (the GFN2 total energy) |

Both functions validate the geometry first (length mismatch, non-finite
coordinate, empty system).

### Python

```python
g = sci_form.compute_pm3_gradient(elements, coords)   # coords: flat [x0, y0, z0, x1, ...]
g.gradients, g.energy, g.heat_of_formation
x = sci_form.compute_xtb_gradient(elements, coords)
x.gradients, x.energy
```

### TypeScript / WASM

```typescript
const g = JSON.parse(compute_pm3_gradient(JSON.stringify(elements), JSON.stringify(coordsFlat)));
// {"gradients": [[gx, gy, gz], ...], "energy": ..., "heat_of_formation": ...}
const x = JSON.parse(compute_xtb_gradient(JSON.stringify(elements), JSON.stringify(coordsFlat)));
// {"gradients": [[gx, gy, gz], ...], "energy": ...}
```

### CLI

```
sci-form pm3-gradient <ELEMENTS_JSON> <COORDS_JSON>
sci-form xtb-gradient <ELEMENTS_JSON> <COORDS_JSON>
```

Both print the JSON object above.

### Try it

Both gradients run here, in your browser, on a water geometry:

<LiveExample id="compute_pm3_gradient" />

<LiveExample id="compute_xtb_gradient" />

## 6. Validation

| Gate | What it checks | Result |
|---|---|---|
| `tests/gradient_parity.rs::f80_pm3` | PM3 gradient against a 1·10⁻⁴ Å difference of `solve_pm3`, tolerance `max(1·10⁻⁴, 1·10⁻⁵·max\|g\|)` eV/Å | 1.1·10⁻⁵ eV/Å (tolerance 1·10⁻⁴) |
| `tests/gradient_parity.rs::f80_xtb` | xTB gradient against a 1·10⁻⁴ Å difference of `solve_gfn2` | 4.8·10⁻⁵ eV/Å (tolerance 1·10⁻⁴) |
| unit tests `pm3::gradients`, `xtb::gradients` | Newton's third law, zero net force, agreement with an independent step, finiteness | pass |
| `fd_gradient` unit tests | the helper recovers the analytic gradient of a quartic; errors propagate | pass |
| `cross_validation` F42 | PM3 energies against MOPAC 23.2.5 (unchanged by DIIS; same fixed point) | pass |

See also [Kohn–Sham DFT](/algorithm/ks-dft) for the gradient of the DFT
solver (which *is* analytic, and was never wrong), and the
[reaction backends guide](/guide/reaction-backends).
