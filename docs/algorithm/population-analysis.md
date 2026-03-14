# Population Analysis (Mulliken & Löwdin)

Population analysis extracts **per-atom partial charges** from the molecular orbital coefficients produced by the EHT solver. sci-form implements both Mulliken and Löwdin partitioning schemes.

## Overview

Given the density matrix $P$ and overlap matrix $S$ from the EHT solution, population analysis answers: *how many electrons does each atom "own"?*

## Mulliken Population Analysis

<img src="/svg/population-mulliken.svg" alt="Mulliken population analysis" class="svg-diagram" />

### Density Matrix

The density matrix is constructed from occupied molecular orbital coefficients:

$$
P_{\mu\nu} = \sum_{i}^{\text{occ}} n_i \, C_{\mu i} \, C_{\nu i}
$$

where $n_i$ is the occupation number (2 for closed-shell doubly occupied orbitals) and $C_{\mu i}$ is the coefficient of atomic orbital $\mu$ in molecular orbital $i$.

### Mulliken Charge Formula

The Mulliken charge on atom $A$ is:

$$
q_A = Z_A - \sum_{\mu \in A} (PS)_{\mu\mu}
$$

where:
- $Z_A$ is the nuclear charge (number of valence electrons)
- $(PS)_{\mu\mu} = \sum_\nu P_{\mu\nu} S_{\mu\nu}$ is the **gross orbital population** for orbital $\mu$
- The sum runs over all atomic orbitals $\mu$ centered on atom $A$

### Charge Conservation

The total Mulliken charges always sum to the net molecular charge:

$$
\sum_A q_A = Q_{\text{total}}
$$

This follows from $\text{Tr}(PS) = N_{\text{electrons}}$.

## Löwdin Population Analysis

<img src="/svg/population-lowdin.svg" alt="Löwdin population analysis" class="svg-diagram" />

### Symmetric Orthogonalization

Löwdin analysis first orthogonalizes the basis by diagonalizing the overlap matrix:

$$
S = U \Lambda U^T \quad \Rightarrow \quad S^{-1/2} = U \Lambda^{-1/2} U^T
$$

### Löwdin Charge Formula

The Löwdin density matrix in the orthogonalized basis is:

$$
\tilde{P} = S^{-1/2} P \, S^{-1/2}
$$

The Löwdin charge on atom $A$:

$$
q_A = Z_A - \sum_{\mu \in A} \tilde{P}_{\mu\mu}
$$

### Advantages Over Mulliken

| Property | Mulliken | Löwdin |
|----------|----------|--------|
| Basis-set dependence | Strong — charges shift with basis size | Weaker — more stable |
| Orbital populations | Can be negative | Always $0 \leq \tilde{P}_{\mu\mu} \leq 2$ |
| Overlap handling | Splits 50/50 between atoms | Orthogonalizes first |
| Rotational invariance | No | Yes |

## API

### Rust

```rust
use sci_form::compute_population;

let result = compute_population("O", None);
// result.mulliken_charges: Vec<f64>
// result.lowdin_charges: Vec<f64>
// result.gross_populations: Vec<f64>
// result.bond_orders: Vec<Vec<f64>>
```

### CLI

```bash
sci-form population "O"
# Output: JSON with mulliken_charges, lowdin_charges, etc.
```

### Python

```python
import sci_form
result = sci_form.population("O")
print(result.mulliken_charges)  # [-0.33, 0.17, 0.17]
```

### WASM

```typescript
import { compute_population } from 'sci-form';
const result = compute_population("O");
// result.mulliken_charges, result.lowdin_charges
```

## Validation

- **Charge conservation**: $\sum_i q_i = Q_{\text{total}}$ tested for all molecules
- **Electronegativity ordering**: O charges < C charges < H charges (expected from EN)
- **Equivalence**: Symmetry-equivalent atoms get equal charges (e.g., H atoms in CH₄)
- **Boundedness**: Löwdin populations always in $[0, 2]$
