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

Z > 30 uses a small-core ECP: n_core electrons are removed, the nuclear charge is Z_eff, and a local Gaussian replaces the core. f electrons stay in the core. This is not a fitted Stuttgart/LANL potential.

The Coulomb matrix uses Cartesian Obara–Saika ERIs through d Gaussians (memoized VRR). The PySCF gate (F77, 2026-10-14) covers H₂/H₂O/CH₄/NH₃/N₂/CO × SVWN/PBE vs PySCF RKS/STO-3G (`generate_dft_reference.py`, `tests/alpha_dft_pyscf.rs --features alpha-dft`): dE ≤ 7.2e-5 Ha, dgap ≤ 9.7e-3 eV. Sc–Zn all-electron bases include 3d; a Zn atom RKS smoke checks that the d Coulomb is finite and the SCF runs.

Closed-shell RKS is the default. UKS runs when `multiplicity > 1`, the valence electron count is odd, or `unrestricted` is set. Spin-polarized Slater + unpolarized VWN is used for SVWN.

`DftConfig.basis = SplitValence` selects 3-21G for H, C, and O.

`ks_nuclear_gradient` is the frozen-P KS force: analytic integrals are differentiated by central difference and the Pulay term −W:∂S/∂R is included.

## Entry point

```rust
use sci_form::dft::{solve_ks_dft, DftConfig, DftMethod};

let result = solve_ks_dft(&elements, &positions, &DftConfig {
    method: DftMethod::Pbe,
    ..DftConfig::default()
})?;
```

Compile with `--features alpha-dft`.
