# Kohn–Sham DFT

sci-form carries its own Kohn–Sham solver behind the `alpha-dft` feature. It is not a binding to REST, hartree, or Libxc. The Fock matrix is

\[
F = H_\text{core} + J + V_{xc}
\]

with SVWN or PBE on a Becke grid. There is no exact exchange.

## Basis this cycle

Hehre STO-3G is used unchanged for H, He, C, N, O, F, P, S, and Cl.

Every other element through Zn gets an all-electron minimal basis. The exponent of each Slater group is

\[
\zeta = (Z - \sigma) / n^*
\]

with Slater's shielding rules. The radial function is the STO-3G three-Gaussian fit, \(\alpha_i = \alpha_i(\zeta=1)\,\zeta^2\). A valence p shell is kept even when it is empty, so the calculation has virtual orbitals. 3d functions are present from Sc onward.

Z > 30 is rejected. An all-electron minimal basis is the wrong tool for 4d, 5d, and f-block metals; the next step is a small-core ECP.

The SCF is closed-shell. An odd electron count drops the unpaired electron the same way restricted HF does. UKS is still ahead.

## Entry point

```rust
use sci_form::dft::{solve_ks_dft, DftConfig, DftMethod};

let result = solve_ks_dft(&elements, &positions, &DftConfig {
    method: DftMethod::Pbe,
    ..DftConfig::default()
})?;
```

Compile with `--features alpha-dft`.
