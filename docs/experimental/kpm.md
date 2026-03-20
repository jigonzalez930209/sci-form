# E4: Kernel Polynomial Method (KPM)

**Module:** `sci_form::experimental::kpm`
**Feature flag:** `experimental-kpm`

---

## Overview

Computes the density of states (DOS), density matrix, and Mulliken populations via **Chebyshev polynomial expansion** of the Hamiltonian, avoiding diagonalization entirely. Achieves $O(N)$ scaling for sparse Hamiltonians, enabling electronic structure calculations on systems of thousands of atoms.

<SvgDiagram src="/svg/experimental-kpm.svg" alt="KPM Pipeline" />

---

## Theory

### Chebyshev Expansion

Any function of a Hermitian matrix can be expanded in Chebyshev polynomials:

$$f(H) \approx \sum_{k=0}^{M} c_k T_k(\tilde{H})$$

where $\tilde{H} = (H - bI)/a$ is rescaled to $[-1, 1]$.

### Recursion

The Chebyshev polynomials are computed iteratively:

$$T_0 = I, \quad T_1 = \tilde{H}, \quad T_{k+1} = 2\tilde{H}T_k - T_{k-1}$$

Each step costs $O(\text{nnz})$ (one sparse matrix-vector product).

### Jackson Kernel

Gibbs oscillations are suppressed by the Jackson damping kernel:

$$g_k^{(M)} = \frac{(M-k+1)\cos\frac{\pi k}{M+1} + \sin\frac{\pi k}{M+1}\cot\frac{\pi}{M+1}}{M+1}$$

### Stochastic Trace

The trace is estimated stochastically:

$$\text{Tr}[A] \approx \frac{1}{n_v} \sum_{i=1}^{n_v} r_i^T A\, r_i$$

using $n_v$ random probe vectors $r_i$ with entries $\pm 1$.

---

## API

```rust
use sci_form::experimental::kpm::*;

// Spectral bounds
let (e_min, e_max) = estimate_spectral_bounds(&hamiltonian);

// Rescale matrix to [-1, 1]
let h_scaled = rescale_matrix(&hamiltonian, e_min, e_max);

// Chebyshev expansion (exact trace for small systems)
let expansion = ChebyshevExpansion::from_matrix_exact(&h_scaled, 100);

// Jackson kernel damping
let kernel = jackson_kernel(100);

// DOS at a specific energy
let dos = expansion.dos_at_energy(0.0, &kernel);

// Full KPM DOS computation
let config = KpmConfig { order: 100, n_random: 20, ..Default::default() };
let dos_result = compute_kpm_dos(&hamiltonian, &config);
// dos_result.energies, dos_result.dos, dos_result.homo_lumo_gap

// KPM Mulliken charges
let mulliken = compute_kpm_mulliken(&hamiltonian, &overlap, n_electrons, &config);
// mulliken.charges, mulliken.total_electrons
```

---

## Configuration

| Parameter | Default | Description |
|-----------|---------|-------------|
| `order` | 100 | Number of Chebyshev moments |
| `n_random` | 20 | Random vectors for stochastic trace |
| `temperature` | 300.0 | Electronic temperature (K) for Fermi smearing |

---

## Scaling

| System Size | Exact Diag. | KPM |
|-------------|------------|-----|
| $N = 100$ | $O(10^6)$ | $O(10^4)$ |
| $N = 1000$ | $O(10^9)$ | $O(10^5)$ |
| $N = 10000$ | Infeasible | $O(10^6)$ |

---

## Tests

```bash
cargo test --features experimental-kpm --test regression -- test_kpm
```

9 integration tests covering: spectral bounds estimation, Chebyshev recursion stability, Jackson kernel positivity, DOS peak detection, band gap identification, and Mulliken charge conservation.
