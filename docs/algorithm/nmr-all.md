# All-nuclide NMR registry

The existing NMR module already knows 85 NMR-active isotopes and can run a GIAO
shielding calculation for a requested nucleus. This page documents the
orchestration layer that answers two practical questions:

1. which isotopes are active for an element?
2. can the current GIAO route evaluate every element present in this molecule?

> **Status.** Registry and multi-nucleus driver implemented in `feat/wave-c`.
> Quantitative all-element GIAO still depends on basis/ECP coverage; the driver
> reports unsupported elements explicitly. See the [v2.x roadmap](/v2-roadmap).

## Registry

`active_isotopes(z)` returns every isotope in `NmrNucleus::ALL` with that atomic
number:

```rust
use sci_form::nmr_all::active_isotopes;

let carbon = active_isotopes(6);
assert!(carbon.iter().any(|isotope| isotope.label == "13C"));
```

`active_elements()` returns the sorted, deduplicated list of elements that have
at least one registered isotope.

## Multi-nucleus driver

`compute_all_nmr` finds the unique elements in the molecule, selects the default
isotope of each, and calls the existing GIAO route once per element:

```rust
use sci_form::nmr_all::compute_all_nmr;

let result = compute_all_nmr(
    &[8, 1, 1],
    &water,
    &sci_form::GiaoNmrConfig::default(),
)?;

for nucleus in &result.nuclei {
    println!("{}: {} shift(s)", nucleus.nucleus, nucleus.chemical_shifts.len());
}
for failure in &result.failures {
    println!("skipped Z={}: {}", failure.element, failure.error);
}
```

### Result

| Type | Meaning |
|---|---|
| `NmrIsotopeInfo` | Element, mass number and canonical label |
| `NmrFailure` | Element, attempted label and error |
| `AllNmrResult` | Successful GIAO results, failures and notes |

A failure is explicit. An unsupported element is never reported as a zero shift.

## Bindings

| Surface | Entry point |
|---|---|
| CLI | not exposed |
| WASM | `active_nmr_isotopes(element)`, `compute_all_nmr(elements, coords, config_json)` |
| Python | `active_nmr_isotopes(element)`, `compute_all_nmr(elements, coords, config_json)` |

Both return JSON. `compute_all_nmr` reports unsupported elements in `failures`.
## Coverage and limits

- The registry is the existing 85-isotope table; it is not expanded by this
  module.
- The GIAO route is closed-shell and STO-3G. Elements without a supported basis
  appear in `failures`.
- One SCF is run per element. This is simple and transparent, not the fastest
  possible multi-nucleus driver.
- Quantitative heavy-element shifts need relativistic corrections, ECPs and a
  larger basis; that is follow-up work.
- The driver selects the default isotope only. A caller that needs a specific
  isotope can call `compute_giao_nmr` directly.
- No Python, WASM or CLI wrapper yet.

## Validation

| Gate | Observable | Result |
|---|---|---|
| Carbon registry | ¹³C present | yes |
| Element list | sorted and unique | yes |
| Water driver | H result or explicit failure | covered |
| Invalid geometry | length mismatch | rejected |

## Next steps

1. Reuse one SCF across every nucleus instead of one SCF per element.
2. Add basis/ECP coverage for the full registry.
3. Add isotope selection and reference standards per nuclide.
4. Add relativistic corrections for heavy elements.
5. Add Python, WASM and CLI wrappers.
6. Benchmark against ORCA for a multi-element set.
