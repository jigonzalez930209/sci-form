# Conventional MP2 correlation energy

The v2.3 wavefunction prototype adds a closed-shell MP2 correlation energy on
top of the existing RHF/STO-3G orbitals. It is a conventional four-index MP2,
not yet the RI/local variant.

> **Status.** Implemented in `feat/wave-c`. The RI approximation, local PNO
> correlation and a PySCF parity benchmark are pending. See the
> [v2.x roadmap](/v2-roadmap).

## Theory

For closed-shell spatial orbitals,

$$
E_{\mathrm{MP2}} =
\sum_{ij}^{\mathrm{occ}}
\sum_{ab}^{\mathrm{vir}}
\frac{
\left[2(ia|jb)-(ib|ja)\right](ia|jb)
}{
\epsilon_i+\epsilon_j-\epsilon_a-\epsilon_b
}.
$$

The implementation transforms the AO ERIs to the MO basis one index at a time:

$$
(ia|jb) =
\sum_{\mu\nu\lambda\sigma}
C_{\mu i}C_{\nu a}(\mu\nu|\lambda\sigma)C_{\lambda j}C_{\sigma b}.
$$

The quarter transforms keep the memory at `O(N² N_occ N_vir)` instead of
materialising the full `N⁴` MO tensor.

## API

```rust
use sci_form::compute_mp2_energy;

let result = compute_mp2_energy(&[8, 1, 1], &water)?;
println!("HF   = {:.6} Ha", result.hf_energy_hartree);
println!("MP2  = {:.6} Ha", result.correlation_energy_hartree);
println!("Total= {:.6} Ha", result.total_energy_hartree);
```

### Output

| Field | Meaning |
|---|---|
| `hf_energy_hartree` | RHF total energy |
| `correlation_energy_hartree` | MP2 correlation energy |
| `correlation_energy_kcal_mol` | Correlation energy in kcal/mol |
| `total_energy_hartree` | HF + MP2 |
| `n_basis`, `n_occupied`, `n_virtual` | Orbital dimensions |

## Scope

| Capability | Status |
|---|---|
| Closed-shell conventional MP2 | implemented |
| RHF/STO-3G reference | implemented |
| RI-MP2 / auxiliary basis | not implemented |
| SCS/SOS scaling | not implemented |
| Local PNO correlation | not implemented |
| Frozen core | not implemented |
| Open-shell MP2 | not implemented |

## Validation

| Gate | Observable | Result |
|---|---|---|
| Water correlation | Negative and finite | passes |
| Energy decomposition | Total = HF + correlation | < 1e-12 |
| Invalid geometry | Length mismatch | rejected |

A PySCF MP2 comparison and a basis-convergence study are still required before
promotion.

## Bindings

| Surface | Entry point |
|---|---|
| CLI | `sci-form mp2 '<elements>' '<coords>'` |
| WASM | `compute_mp2_energy(elements, coords)` |
| Python | `compute_mp2_energy(elements, coords)` |

The CLI prints the full JSON result.

## Limitations

- Conventional `O(N⁵)` transform; not suitable for large molecules.
- STO-3G only.
- All occupied orbitals are correlated; no frozen core.
- No RI, SCS/SOS or local correlation.
- No analytical MP2 gradient or response properties.
- No Python, WASM or CLI wrapper.
- No benchmark against PySCF/ORCA yet.

## Next steps

1. Add def2-SVP/TZVP and auxiliary basis sets.
2. Add RI-MP2 and SCS/SOS scaling.
3. Add frozen-core handling.
4. Add local PNO correlation.
5. Add PySCF parity and basis-convergence gates.
6. Add Python, WASM and CLI wrappers.
