# EPR hyperfine coupling

The first v2.2 EPR prototype evaluates the isotropic Fermi-contact coupling
from the UHF spin density. It is a light-element, single-reference model; the
g-tensor, anisotropic tensor, relativistic corrections and paramagnetic shifts
are documented follow-ups.

> **Status.** Implemented in `feat/wave-b` (Rust). Python, WASM and CLI wrappers
> are pending. See the [v2.x roadmap](/v2-roadmap).

## Fermi-contact coupling

For a nucleus `A`, the isotropic hyperfine coupling is proportional to the
spin density at the nucleus:

$$
A_{\mathrm{iso}} \propto g_N\,\rho_s(\mathbf{R}_A),
\qquad
\rho_s(\mathbf{r}) =
\sum_{\mu\nu}
\left(P^{\alpha}_{\mu\nu}-P^{\beta}_{\mu\nu}\right)
\phi_\mu(\mathbf{r})\phi_\nu(\mathbf{r}).
$$

The implementation uses the proton as the reference:

$$
A_{\mathrm{iso}} =
1420.405751\ \mathrm{MHz}
\cdot
\frac{g_N}{g_p}
\cdot
\frac{\rho_s(\mathbf{R}_A)}{\rho_{1s}^{\mathrm{H}}(0)}.
$$

`ρ_{1s}^{H}(0)` is evaluated from the same STO-3G basis, so the conversion is
consistent with the basis instead of assuming an analytic hydrogen density.

## Nuclear g-factors

The built-in table covers the common light nuclei and several heavier probes:

| Z | Nucleus | Z | Nucleus |
|---|---|---|---|
| 1 | ¹H | 19 | ³⁹K |
| 6 | ¹³C | 25 | ⁵⁵Mn |
| 7 | ¹⁴N / ¹⁵N | 26 | ⁵⁷Fe |
| 8 | ¹⁷O | 30 | ⁶⁷Zn |
| 9 | ¹⁹F | 31 | ⁶⁹Ga / ⁷¹Ga |
| 15 | ³¹P | 53 | ¹²⁷I |
| 16 | ³³S | 78 | ¹⁹⁵Pt |

The table is intentionally small and explicit. An element without a g-factor is
reported as skipped, not silently assigned a default.

## API

```rust
use sci_form::epr::EprConfig;
use sci_form::compute_epr_hyperfine;

let result = compute_epr_hyperfine(
    &[1],
    &[[0.0, 0.0, 0.0]],
    &EprConfig {
        nucleus: Some(1),
        multiplicity: 2,
        ..EprConfig::default()
    },
)?;

println!("A_iso = {:.3} MHz", result.entries[0].a_iso_mhz);
```

### Input

| Field | Meaning |
|---|---|
| `nucleus` | Target atomic number; `None` evaluates every known nucleus |
| `charge` | Total molecular charge |
| `multiplicity` | Spin multiplicity; default doublet |

### Output

| Type | Meaning |
|---|---|
| `HyperfineEntry` | Atom index, element, g-factor, spin density and `A_iso` |
| `EprHyperfineResult` | Entries, multiplicity and notes |

Spin density is in electrons/Bohr³; coupling is in MHz.

## Validation

| Gate | Observable | Result |
|---|---|---|
| Hydrogen atom | `A_iso` | 1420.405751 MHz within 5 MHz |
| Closed-shell water | H spin density | < 1 MHz |
| Unknown g-factor | Unsupported nucleus | explicit error |
| Multiplicity | Invalid value | rejected |

The hydrogen test is a real calibration: the STO-3G proton reference density is
computed from the basis, so the result is not hard-coded.

## Bindings

| Surface | Entry point |
|---|---|
| CLI | `sci-form epr '<elements>' '<coords>' --nucleus 1 --multiplicity 2` |
| WASM | `compute_epr_hyperfine(elements, coords, nucleus, charge, multiplicity)` |
| Python | `compute_epr_hyperfine(elements, coords, nucleus=None, charge=0, multiplicity=2)` |

The CLI prints the full JSON result. `nucleus = 0` in WASM means "all known
nuclei".

## Limitations

- Only the isotropic Fermi-contact term; the anisotropic hyperfine tensor and
  the g-tensor are not implemented.
- UHF/STO-3G spin polarization is approximate; spin contamination is not
  corrected.
- Relativistic effects are absent, so heavy-element couplings are screening
  values at best.
- The g-factor table is finite; an unknown nucleus is skipped.
- No paramagnetic contact shift, no zero-field splitting, no EPR spectrum
  simulation.
- No conformer averaging and no Python/WASM/CLI wrapper yet.

## Next steps

1. Add the anisotropic spin-dipole tensor.
2. Add the g-tensor response.
3. Add relativistic corrections and ECP support.
4. Add spin-contamination diagnostics.
5. Add conformer averaging.
6. Validate against ORCA for organic radicals.
7. Add Python, WASM and CLI wrappers.
