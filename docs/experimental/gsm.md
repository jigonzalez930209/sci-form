# E11: Growing String Method (GSM)

**Module:** `sci_form::experimental::gsm`
**Feature flag:** `experimental-gsm`

---

## Overview

Implements a **reaction path and transition state finder** driven purely by force-field gradients. Unlike NEB (Nudged Elastic Band), GSM does not require an initial guess of the full path — it grows the string from both reactant and product endpoints, finding the minimum energy path and saddle point with fewer energy evaluations.

---

## Theory

### String Initialization

Given reactant $R$ and product $P$ coordinates, the string starts with 3 nodes: $R$, an interpolated midpoint, and $P$.

### Growth Step

At each iteration, new nodes are added along the string:

1. **Tangent:** $\hat{\tau} = \frac{x_{i+1} - x_{i-1}}{|x_{i+1} - x_{i-1}|}$ (central difference)
2. **Perpendicular gradient:** $g^\perp = g - (g \cdot \hat{\tau})\hat{\tau}$
3. **Optimization:** Interior nodes are optimized along the perpendicular gradient
4. **New node:** Added at distance $\Delta s$ along the interpolation direction

### Junction and Saddle Detection

- **Junction:** When tips from both ends meet (distance < $\Delta s / 2$)
- **Saddle point:** Node with maximum energy along the converged string
- **Refinement:** Steepest descent on the perpendicular gradient around the saddle candidate

---

## API

```rust
use sci_form::experimental::gsm::*;

// Interpolate between two geometries
let midpoint = interpolate_node(&reactant, &product, 0.5);

// Energy function (any force field)
let energy_fn = |coords: &[f64]| -> f64 { /* UFF, MMFF94, EHT, etc. */ };

// Configure GSM
let config = GsmConfig {
    max_nodes: 11,      // nodes along the string
    step_size: 0.01,    // optimization step size
    grad_tol: 0.5,      // gradient tolerance
    max_iter: 100,       // max optimization iterations per growth
    fd_step: 0.005,      // finite difference step for gradient
};

// Grow the string
let path = grow_string(&reactant, &product, &energy_fn, &config);
// path.nodes: Vec<Vec<f64>>
// path.energies: Vec<f64>
// path.joined: bool
// path.iterations: usize

// Find transition state
let result = find_transition_state(&reactant, &product, &energy_fn, &config);
// result.ts_coords: Vec<f64>
// result.ts_energy: f64
// result.activation_energy: f64     (forward barrier)
// result.reverse_barrier: f64       (reverse barrier)
// result.path_coords: Vec<Vec<f64>>
// result.n_nodes: usize
// result.energy_evaluations: usize

// Refine saddle point
let (refined_coords, refined_energy) = refine_saddle(
    &ts_guess, &energy_fn, (&prev_node, &next_node),
    max_steps, fd_step, step_size,
);
```

---

## Energy Backends

GSM works with any energy function that maps coordinates to a scalar:

| Backend | Speed | Accuracy | Best For |
|---------|-------|----------|----------|
| UFF | Fast | Low | Quick screening |
| MMFF94 | Fast | Medium | Organic reactions |
| EHT | Medium | Medium | π-conjugated systems |
| PM3 | Slow | High | Reaction barriers |
| xTB | Slow | High | Organometallic reactions |

---

## Output

The transition state result includes:

- **Transition state geometry** and energy
- **Activation energy** (forward barrier: $E_{TS} - E_R$)
- **Reverse barrier** ($E_{TS} - E_P$)
- **Full reaction path** as a sequence of geometries with energies
- **Energy evaluation count** for benchmarking efficiency

---

## Tests

```bash
cargo test --features experimental-gsm --test regression -- test_gsm
```

8 integration tests covering: linear interpolation, quarter-point interpolation, string growth, transition state finding on a double-well potential, path ordering, energy evaluation counting, symmetric barrier detection, and saddle point refinement.
