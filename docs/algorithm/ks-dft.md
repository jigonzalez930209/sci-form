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

## Nuclear gradient

`ks_nuclear_gradient(elements, positions, &config)` returns a `GradientResult`
with the nuclear gradient \(dE/dR\) in **Hartree/Bohr** (one `[x, y, z]` per
atom, `positions` in Å), the RMS and maximum components, and the total energy
at the evaluated geometry (Hartree).

It is the frozen-density force: the AO integrals (core Hamiltonian, ERIs, the
exchange-correlation matrix and the overlap) are differentiated by central
differences of 10⁻⁴ Bohr at the converged density, and the Pulay term
\(-\mathrm{Tr}(W\,\partial S/\partial R)\) with the energy-weighted density
\(W\) is added. It costs one SCF plus the integral differences, not \(6N\) SCFs.

<SvgDiagram src="/svg/ks-gradient-units.svg" alt="The KS-DFT analytic gradient in Hartree per Bohr equals the finite-difference energy derivative per angstrom multiplied by 0.5292" />

**Validation (F80).** The gradient agrees with the central-difference
derivative of the converged `solve_ks_dft` energy, converted to Hartree/Bohr,
to \(3\cdot10^{-7}\) Ha/Bohr on H₂O₂ (largest component 0.227 Ha/Bohr,
tolerance \(2.3\cdot10^{-6}\)). The unit test
`dft::gradient::tests::h2_ks_gradient_is_the_energy_derivative_in_ha_per_bohr`
checks the same on H₂ with an independent step.

**History.** For a while the F80 gate was `#[ignore]`d and the module was
described as missing an "adiabatic XC response term", because the analytic
gradient differed from the numerical one by a factor of exactly 0.5292 on every
component. That factor is Å per Bohr: the numerical reference differentiated
the energy per Å and was compared with a gradient in Hartree/Bohr. The
gradient was correct; the reference is now converted, the gate is enforced and
no response term is needed at this level of theory. (The PM3 and xTB gradients
of sci-form were genuinely wrong and were rebuilt — see
[PM3 and xTB gradients](/algorithm/semiempirical-gradients).)

**Open-shell (UKS).** When the solution is unrestricted (`result.unrestricted`),
the gradient uses the two spin densities \(P^\alpha, P^\beta\): the Coulomb and
core terms use \(P^\alpha+P^\beta\), the exchange-correlation energy is the spin-polarised
one (`build_vxc_uks`), and the Pulay weight is
\(W=\sum_\sigma\sum_{i\in\sigma}\varepsilon_{i\sigma}\,c_{i\sigma}c_{i\sigma}^{T}\)
(occupation 1 per spin orbital, not 2). `DftResult` therefore carries
`orbital_energies_beta` and `coefficients_beta` (empty for a restricted
solution). A second F80 gate covers the OH radical (doublet): \(10^{-6}\) Ha/Bohr
with `energy_threshold = 1e-13` (largest component 0.28 Ha/Bohr). The UKS SCF stops on the
energy change only, so a looser threshold leaves ~7·10⁻⁶ Ha/Bohr of density noise in the
gradient; tighten `energy_threshold` when you need forces.

*History.* Before this fix the UKS gradient reused the restricted formulas (one density, occupation 2,
alpha orbital energies only) and was off by 1.4·10⁻² Ha/Bohr (5 %) on OH.

### Try it

<LiveExample id="compute_dft" />

## Entry point

```rust
use sci_form::dft::{solve_ks_dft, DftConfig, DftMethod};

let result = solve_ks_dft(&elements, &positions, &DftConfig {
    method: DftMethod::Pbe,
    ..DftConfig::default()
})?;
```

Compile with `--features alpha-dft`.
