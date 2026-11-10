# Core-level XPS binding energies

The first v2.3 core-spectroscopy prototype estimates core-level binding
energies from the UHF/STO-3G orbital energies. It is a Koopmans screening model:
it does **not** perform a ΔSCF core-hole calculation and does **not** implement
XANES.

> **Status.** Implemented in `feat/wave-c` (Rust). Python, WASM and CLI wrappers
> are pending. See the [v2.x roadmap](/v2-roadmap).

## Koopmans estimate

For a core-like orbital `c`,

$$
E_b \approx -\epsilon_c,
$$

where `ε_c` is the orbital energy in Hartree and `E_b` is reported in eV. The
implementation picks, for each atom, the occupied s orbital with the largest
weight on that atom's s basis functions. That is the most atom-localized core-
like state in the minimal basis.

This is the simplest useful XPS estimate. It omits core-hole relaxation, final-
state screening and relativistic corrections, so absolute values are screening
numbers.

## API

```rust
use sci_form::xps::XpsConfig;
use sci_form::compute_xps_binding_energies;

let result = compute_xps_binding_energies(
    &[6, 1, 1, 1, 1],
    &methane,
    &XpsConfig {
        elements: Some(vec![6]),
        ..XpsConfig::default()
    },
)?;

for binding in result.bindings {
    println!(
        "atom {} Z={} E_b={:.1} eV (localization {:.2})",
        binding.atom_index, binding.element, binding.binding_energy_ev, binding.localization
    );
}
```

### Output

| Field | Meaning |
|---|---|
| `atom_index` | Atom index |
| `element` | Atomic number |
| `orbital_energy_hartree` | Orbital energy used |
| `binding_energy_ev` | Koopmans binding energy |
| `localization` | Weight on the atom's s functions |

## Validation

| Gate | Observable | Result |
|---|---|---|
| Methane C 1s | Binding energy | 200–400 eV |
| Invalid geometry | Length mismatch | rejected |

The methane test checks the core-level scale, not chemical accuracy. A ΔSCF
reference and an experimental C/N/O 1s set are still required.

## Limitations

- Koopmans approximation, no core-hole relaxation or final-state screening.
- UHF/STO-3G only; absolute values are not benchmarked.
- The selected orbital is the most atom-localized s orbital, not a symmetry-
  labelled core state.
- No XANES, transition potential or core excitation.
- No relativistic corrections for heavy elements.
- No spectrum broadening or peak assignment.
- No Python, WASM or CLI wrapper.

## Next steps

1. Add ΔSCF core-hole energies.
2. Add a transition-potential or restricted-excitation XANES path.
3. Add relativistic corrections and ECPs for heavy elements.
4. Add C/N/O 1s experimental benchmarks.
5. Add spectrum broadening.
6. Add Python, WASM and CLI wrappers.
