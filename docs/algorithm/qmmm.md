# QM/MM with mechanical embedding

The first v2.2 multiscale prototype combines a quantum region with a classical
point-charge environment. It is the smallest useful QM/MM model: the QM region
is computed in vacuum, its Mulliken charges are coupled to the MM point
charges, and the MM-MM point-charge term is added.

> **Status.** Implemented in `feat/wave-b`. This is **mechanical embedding**,
> not electrostatic embedding. Link atoms, bonded MM terms and polarizable MM
> are follow-ups. See the [v2.x roadmap](/v2-roadmap).

## Why mechanical embedding first

Electrostatic embedding adds the MM point charges to the QM Hamiltonian:

$$
H_{\mathrm{QM/MM}} =
H_{\mathrm{QM}}
+
\sum_m\sum_i \frac{q_m}{|\mathbf{r}_i-\mathbf{R}_m|}
+
\sum_m\sum_A \frac{q_mZ_A}{|\mathbf{R}_A-\mathbf{R}_m|}.
$$

That requires a QM core that accepts an external potential. The current SCF
path does not. Mechanical embedding is the honest first step: it keeps the QM
calculation unchanged and couples the two subsystems through the converged QM
Mulliken charges.

## Energy expression

The total energy is

$$
E_{\mathrm{QM/MM}} =
E_{\mathrm{QM}}
+
E_{\mathrm{QM-MM}}
+
E_{\mathrm{MM-MM}},
$$

with the point-charge coupling

$$
E_{\mathrm{QM-MM}} =
k_C
\sum_{i\in\mathrm{QM}}
\sum_{m\in\mathrm{MM}}
\frac{q_i q_m}{r_{im}},
$$

and the MM-MM term

$$
E_{\mathrm{MM-MM}} =
k_C
\sum_{m<n}
\frac{q_m q_n}{r_{mn}}.
$$

`k_C = 332.06371 kcal·Å/(mol·e²)`. Positions are Å, charges are elementary
charge and energies are kcal/mol.

## API

```rust
use sci_form::qmmm::{QmMmConfig, QmMmSystem};
use sci_form::compute_qmmm_energy;

let system = QmMmSystem {
    qm_elements: vec![8, 1, 1],
    qm_positions: vec![
        [0.0, 0.0, 0.0],
        [0.9572, 0.0, 0.0],
        [-0.239987, 0.927297, 0.0],
    ],
    mm_elements: vec![8],
    mm_positions: vec![[0.0, 0.0, 5.0]],
    mm_charges: vec![-0.5],
};

let result = compute_qmmm_energy(&system, &QmMmConfig::default())?;
println!("{:.3} kcal/mol", result.total_energy_kcal_mol);
```

### Backends

| Method | QM region | Charge support |
|---|---|---|
| `QmMethod::Uhf` | UHF/STO-3G | charged and open-shell |
| `QmMethod::Xtb` | GFN2-xTB | neutral only |

The xTB backend returns `Unsupported` for a charged QM region rather than
silently using the wrong electron count.

### Output

| Field | Meaning |
|---|---|
| `qm_energy_kcal_mol` | QM energy |
| `qm_mm_energy_kcal_mol` | QM-MM point-charge coupling |
| `mm_mm_energy_kcal_mol` | MM-MM point-charge coupling |
| `total_energy_kcal_mol` | Sum of the three terms |
| `qm_charges` | Mulliken charges of the QM region |
| `n_qm`, `n_mm` | Region sizes |

## Validation

| Gate | Observable | Result |
|---|---|---|
| Empty MM region | Total energy equals QM energy | < 1e-12 kcal/mol |
| QM-MM Coulomb | Energy vs manual charge sum | < 1e-8 kcal/mol |
| Length mismatch | MM arrays | rejected |
| Charged xTB | Neutral-only backend | `Unsupported` |

These are implementation gates. A QM/MM benchmark against ONIOM or an
electrostatic-embedding reference is still required.

## Bindings

| Surface | Entry point |
|---|---|
| CLI | `sci-form qmmm '<QmMmSystem-json>' '<QmMmConfig-json>'` |
| WASM | `compute_qmmm_energy(system_json, config_json)` |
| Python | `compute_qmmm_energy(system_json, config_json)` |

The JSON payloads are the serialized `QmMmSystem` and `QmMmConfig`.
## Limitations

- Mechanical embedding only; the MM charges do not polarize the QM density.
- No link atoms, no bonded MM terms and no LJ term. The MM region is a set of
  point charges.
- No periodic boundaries.
- The QM-MM coupling uses Mulliken charges, which are basis-set dependent.
- The QM and MM subsystems must be supplied explicitly; automatic partitioning
  is not implemented.
- No gradients or geometry optimization yet.
- No Python, WASM or CLI wrapper.

## Next steps

1. Add an external point-charge potential to the SCF Hamiltonian (electrostatic
   embedding).
2. Add link atoms and boundary validation.
3. Add LJ and bonded MM terms.
4. Add the QM/MM gradient and geometry optimization.
5. Add PBC for surface and enzyme models.
6. Benchmark against ONIOM/ORCA.
7. Add Python, WASM and CLI wrappers.
