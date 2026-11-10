# Hybrid DFT: PBE0 exact exchange

The v2.3 electronic-structure prototype adds **PBE0** to the gated KS-DFT
engine. PBE0 is PBE exchange-correlation plus 25 % exact Hartree–Fock exchange.
The exact-exchange path is currently restricted to closed-shell systems.

> **Status.** Implemented in `feat/wave-c` behind `alpha-dft`. Range-separated
> hybrids, unrestricted exact exchange and the ORCA/PySCF parity gate are still
> pending. See the [v2.x roadmap](/v2-roadmap).

## Theory

The hybrid exchange-correlation energy is

$$
E_{xc} = aE_x^{\mathrm{HF}} + (1-a)E_x^{\mathrm{PBE}} + E_c^{\mathrm{PBE}},
\qquad a = 0.25,
$$

with

$$
E_x^{\mathrm{HF}} = -\frac{1}{4}
\sum_{\mu\nu\lambda\sigma}
P_{\mu\nu}P_{\lambda\sigma}(\mu\lambda|\nu\sigma).
$$

The corresponding Fock matrix is

$$
F = H_{\mathrm{core}} + J + V_{xc} - \frac{a}{2}K,
$$

where

$$
K_{\mu\nu} =
\sum_{\lambda\sigma}P_{\lambda\sigma}(\mu\lambda|\nu\sigma).
$$

The implementation reuses the existing DFT ERI tensor and adds the exchange
matrix once per SCF iteration.

## API

```rust
use sci_form::dft::ks_fock::{solve_ks_dft, DftConfig, DftMethod};

let config = DftConfig {
    method: DftMethod::Pbe0,
    ..DftConfig::default()
};
let result = solve_ks_dft(&elements, &positions, &config)?;
println!("E = {} Hartree", result.total_energy);
```

PBE0 is available through `DftMethod::Pbe0`; `pbe0`, `pbe-0` and `hybrid` are
accepted aliases. The module is gated behind `alpha-dft`, so the caller must
enable that feature.

## Scope

| Capability | Status |
|---|---|
| Restricted closed-shell PBE0 | implemented |
| Unrestricted/open-shell PBE0 | returns `Unsupported` |
| Range-separated hybrids | not implemented |
| B3LYP | not implemented |
| Analytical hybrid gradient | not implemented |

The current gradient path is for the pure DFT functional; PBE0 gradients are a
follow-up.

## Validation

| Gate | Observable | Result |
|---|---|---|
| Water PBE0 SCF | Convergence | passes |
| Water PBE0 energy | Finite and negative | passes |
| Electron count | 10 electrons | passes |
| Unrestricted PBE0 | Explicit error | guarded |

The important missing gate is **method parity** against PySCF/ORCA for the same
PBE0/basis. Until that exists, the implementation is a prototype.

## Bindings

| Surface | Entry point |
|---|---|
| CLI | pending |
| WASM | pending |
| Python | `compute_pbe0_energy(elements, coords, max_iterations=100)` behind the `alpha-dft` feature |

The DFT engine remains feature-gated, so the CLI/WASM surfaces are not built by
default.
## Limitations

- The exact-exchange matrix is built with the existing `O(N^4)` ERI loop; no RI
  or Cholesky exchange yet.
- Restricted closed-shell only.
- PBE0 uses the existing PBE exchange-correlation evaluator; no B3LYP.
- No range separation, no exact-exchange gradient and no response properties.
- The DFT module remains feature-gated and is not exposed on Python/WASM/CLI.
- No benchmark against PySCF/ORCA yet.

## Next steps

1. Add an unrestricted exact-exchange path.
2. Add RI/Cholesky exact exchange for larger systems.
3. Add range-separated hybrids.
4. Add the analytical PBE0 gradient.
5. Add PySCF/ORCA parity and a 20-molecule reaction-energy set.
6. Expose the method on Python/WASM/CLI.
