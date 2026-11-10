# Redox potentials and pKa

This page documents the first implementation of the v2.x electrochemistry
workflow: solution pKa values and redox potentials from an explicit pair of
molecular states. It is a **thermodynamic-cycle wrapper**, not a new electronic
structure engine. The quantum energy comes from an existing closed-shell SCF and
the solvent term comes from the existing ALPB model.

> **Status.** Implemented in the `feat/wave-a` development branch. The
> validation set and the Python/WASM/CLI wrappers are still in progress; see
> the [v2.x roadmap](/v2-roadmap).

## Why an explicit pair

Automatic protonation-site enumeration and automatic redox-state assignment are
separate cheminformatics problems. This API asks the caller for the two states:

| Workflow | State 1 | State 2 |
|---|---|---|
| pKa | acid `HA` | conjugate base `A-` |
| Redox | oxidised `Ox` | reduced `Red` |

That keeps the electronic-structure input unambiguous: the formal charge is read
from the SMILES (`[O-]`, `[NH4+]`, …), and the geometry is embedded independently
for each state.

## pKa cycle

For a deprotonation

$$
\mathrm{HA} \rightarrow \mathrm{A^-} + \mathrm{H^+},
$$

the absolute cycle is

$$
\Delta G_{\mathrm{dep}} =
G(\mathrm{A^-}) + G(\mathrm{H^+}) - G(\mathrm{HA}),
\qquad
\mathrm{p}K_a = \frac{\Delta G_{\mathrm{dep}}}{RT\ln 10}.
$$

`G(H+)` is the proton solvation free energy. Its absolute value depends on the
standard state and on the solvent, so the result always reports the value used.
The default is **-265.9 kcal/mol** for water.

### Relative pKa is the recommended path

Absolute pKa inherits the systematic error of the electronic method and the
proton convention. The relative cycle against a reference acid cancels most of
that error:

$$
\mathrm{p}K_a(\mathrm{target}) =
\mathrm{p}K_a(\mathrm{ref}) +
\frac{\Delta G_{\mathrm{target}}-\Delta G_{\mathrm{ref}}}{RT\ln 10}.
$$

`PkaReference::acetic_acid()` supplies the reference pair and `pKa = 4.76`.
Supplying the target as its own reference reproduces the reference value
exactly, which is the unit-test invariant of the implementation.

## Redox cycle

For the reduction

$$
\mathrm{Ox} + n e^- \rightarrow \mathrm{Red},
$$

the reaction free energy is

$$
\Delta G = G(\mathrm{Red}) - G(\mathrm{Ox}),
$$

with the electron free energy set to zero by convention. The absolute potential
and the potential against a reference electrode are

$$
E_{\mathrm{abs}} = -\frac{\Delta G}{nF},
\qquad
E_{\mathrm{ref}} = E_{\mathrm{abs}} - E_{\mathrm{ref,abs}}.
$$

The default reference is the standard hydrogen electrode at **4.44 V**. A
relative cycle against a couple with a known experimental potential is also
available:

$$
E(\mathrm{target}) =
E(\mathrm{ref}) -
\frac{\Delta G_{\mathrm{target}}-\Delta G_{\mathrm{ref}}}{nF}.
$$

## Method

For each state:

1. Parse the SMILES and read the total formal charge.
2. Embed a conformer ensemble with the existing ETKDG pipeline.
3. Evaluate each conformer with closed-shell UHF/STO-3G.
4. Evaluate ALPB solvation with the SCF Mulliken charges.
5. Add the electronic and solvation energies.
6. Select the minimum or the Boltzmann average at the configured temperature.

The default electronic method is UHF/STO-3G because the current GFN2-xTB path is
neutral-only. The `Xtb` backend is available for neutral states and is faster;
it returns `Unsupported` for a charged state instead of silently using the
wrong electron count.

## API

### Rust

```rust
use sci_form::electrochemistry::{
    ElectrochemistryConfig, PkaReference, RedoxReference,
};
use sci_form::{compute_pka, compute_redox_potential};

let reference = PkaReference::acetic_acid();
let config = ElectrochemistryConfig::default();

let result = compute_pka(
    "FC(F)(F)C(=O)O",
    "FC(F)(F)C(=O)[O-]",
    Some(&reference),
    &config,
)?;

println!("pKa = {:.2}", result.pka);
println!("ΔG = {:.2} kcal/mol", result.delta_g_kcal_mol);
println!("conformer spread = {:.2} kcal/mol", result.acid.spread_kcal_mol);
```

A redox call has the same shape:

```rust
let result = compute_redox_potential(
    oxidized_smiles,
    reduced_smiles,
    /* n_electrons = */ 1,
    Some(&reference_couple),
    &config,
)?;
```

### Python, WASM and CLI

The same workflow is exposed on every surface:

| Surface | Entry points |
|---|---|
| Python | `compute_pka`, `compute_redox_potential`, `PkaResult`, `RedoxResult` |
| WASM | `compute_pka`, `compute_pka_relative`, `compute_redox_potential`, `compute_redox_potential_relative` |
| CLI | `sci-form pka ...`, `sci-form redox ...` |

Python example:

```python
import sci_form

result = sci_form.compute_pka(
    "FC(F)(F)C(=O)O",
    "FC(F)(F)C(=O)[O-]",
    reference_acid="CC(=O)O",
    reference_base="CC(=O)[O-]",
    reference_pka=4.76,
)
print(result.pka, result.relative)
```

CLI example:

```bash
sci-form pka "CC(=O)O" "CC(=O)[O-]" \
  --reference-acid "CC(=O)O" --reference-base "CC(=O)[O-]" \
  --reference-pka 4.76 --n-conformers 1
```

The CLI prints the full JSON result, including every conformer and note.
The Python classes expose the top-level fields; the WASM functions return the
full JSON payload.

### Result fields

| Type | Key fields |
|---|---|
| `StateEnergy` | `charge`, `selected_kcal_mol`, `boltzmann_kcal_mol`, `spread_kcal_mol`, `conformers` |
| `PkaResult` | `pka`, `delta_g_kcal_mol`, `relative`, `reference`, `notes` |
| `RedoxResult` | `potential_vs_reference_v`, `absolute_potential_v`, `delta_g_kcal_mol`, `relative`, `reference`, `notes` |

Energies are in kcal/mol, potentials in volt, temperature in kelvin and pKa is
dimensionless.

## Units and conventions

| Quantity | Unit / convention |
|---|---|
| Electronic energy | kcal/mol after conversion from Hartree |
| Solvation energy | kcal/mol |
| Proton solvation | kcal/mol, default -265.9 |
| Reference electrode | V, default SHE at 4.44 V |
| Temperature | K, default 298.15 |
| pKa | dimensionless |
| Molar gas constant | 0.001987204258 kcal/(mol·K) |
| Faraday constant | 23.060548867 kcal/(mol·V) |

Every result carries the conventions it used. A change in the proton solvation
term or the reference electrode is visible in the output, not hidden in a
constant.

## Validation

The current gates are:

| Gate | Observable | Result |
|---|---|---|
| Formula | Acetic-acid-scale ΔG → pKa | 4.76 ± 0.05 |
| Self-reference | Target = reference acid | exact reference pKa |
| Charge parsing | `[O-]`, `[NH4+]`, neutral | -1, +1, 0 |
| Redox formula | ΔG = -23.06 kcal/mol, n = 1 | E_abs = 1 V |
| TFA trend | Relative to acetic acid | pKa 0.82 (experimental 0.23) |
| Heavy integration | TFA + reference, UHF/STO-3G | ignored by default; run with `--ignored` |

The TFA value shows the expected trend and the size of the systematic error.
It is **not** a claim of 0.5 pKa accuracy. A production validation set with
experimental acids, bases and redox couples is still required before promotion.

## Limitations

- UHF/STO-3G is a minimal-basis, uncorrelated method. Absolute values are
  screening estimates; use the relative cycle for comparisons.
- Thermal and vibrational corrections are not included. The result is an
  electronic-plus-solvation free energy, not a full RRHO free energy.
- The workflow requires a closed-shell, even-electron state.
- GFN2-xTB is neutral-only in the current core. Charged states use UHF.
- Conformer sampling is the existing ETKDG ensemble; a single conformer is the
  default for speed and can miss a lower state.
- The caller must supply the correct conjugate base or redox pair. The API does
  not enumerate protonation or oxidation states.
- The relative result is only as good as the reference value and the similarity
  between target and reference.

## Next steps

1. Add Python, WASM and CLI wrappers.
2. Add a curated acid/base/redox benchmark with experimental references.
3. Add thermal corrections from the existing vibrational analysis.
4. Investigate charged GFN2-xTB support so the fast backend can replace UHF for
   the production path.
5. Add the validation rows to the full matrix and the scientific report.
