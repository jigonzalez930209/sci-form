# Universal MLIP model contract

The v2.2 MLIP prototype defines the **contract** for a user-supplied universal
potential. The library does not bundle MACE, CHGNet, Orb or a private checkpoint:
model weights have their own licences and training provenance. The registry and
callback backend make the contract executable and testable without pretending
that a random or test-weight model is a production potential.

> **Status.** Contract and callback backend implemented in `feat/wave-b`. The
> ONNX checkpoint adapter is the remaining piece; it needs a real exported model
> and a licence review. See the [v2.x roadmap](/v2-roadmap).

## Contract

A model evaluates

```text
(elements, positions, cell?) -> (energy, forces, stress?)
```

with:

| Quantity | Unit |
|---|---|
| positions | Å |
| energy | kcal/mol |
| forces | kcal/mol/Å |
| stress | kcal/mol/Å³ |
| cutoff | Å |

The contract uses kcal/mol because that is the native unit of the existing
force fields and MD integrators. A model trained in eV must convert at the
adapter boundary; the registry does not hide the conversion.

## Manifest

Every model carries a manifest:

```json
{
  "schema_version": 1,
  "name": "mace-off23-small",
  "elements": [1, 6, 7, 8, 9, 15, 16, 17],
  "cutoff_angstrom": 5.0,
  "license": "MIT",
  "supports_stress": false,
  "source": "sha256:..."
}
```

Validation rejects an unknown schema version, an empty element list, an invalid
cutoff, an empty licence and an out-of-range atomic number.

## API

```rust
use sci_form::mlip::{CallbackMlip, MlipManifest, MlipRegistry};

let manifest = MlipManifest {
    schema_version: 1,
    name: "my-model".into(),
    elements: vec![1, 6, 8],
    cutoff_angstrom: 5.0,
    license: "user-supplied".into(),
    supports_stress: false,
    source: None,
};

let model = CallbackMlip::new(manifest, |elements, positions, cell| {
    // Call a user-supplied ONNX runtime or a native backend here.
    Ok((0.0, vec![[0.0, 0.0, 0.0]; elements.len()], None))
})?;

let result = model.evaluate(&[8, 1, 1], &water, None)?;
```

`MlipRegistry` stores named models and rejects duplicate names. The registry is
deliberately small; it is a contract boundary, not a model zoo.

## ONNX adapter contract

An ONNX checkpoint adapter must document:

| ONNX input | Shape | Meaning |
|---|---|---|
| `positions` | `[N, 3]` | Cartesian coordinates in Å |
| `atomic_numbers` | `[N]` | Atomic numbers |
| `cell` | `[3, 3]` optional | Lattice vectors for periodic models |
| `periodic` | `[1]` optional | PBC flag |

| ONNX output | Shape | Meaning |
|---|---|---|
| `energy` | `[1]` | Energy in the model's native unit |
| `forces` | `[N, 3]` | Forces in the model's native unit per Å |
| `stress` | `[3, 3]` optional | Stress in the model's native unit per Å³ |

The adapter must convert eV → kcal/mol at the boundary and reject a model whose
element list does not cover the requested system. It must not silently fall back
to a different model.

## Validation

| Gate | Observable | Result |
|---|---|---|
| Callback contract | Energy and force shape | finite, one force per atom |
| Element coverage | Missing element | `Unsupported` |
| Registry | Duplicate name | rejected |
| Manifest | Empty elements / licence | rejected |
| Stress | Model declares stress but returns none | rejected |

A production gate needs a real checkpoint and a reference trajectory
(MD17/SPICE or a materials benchmark). That is why the adapter is not marked
promoted.

## Limitations

- No ONNX runtime or model weights are included.
- The callback backend is a contract test, not a universal potential.
- No neighbour list, PBC wrapping or stress validation beyond finiteness.
- No training, fine-tuning or uncertainty calibration.
- The registry is in-memory and not serializable.
- No Python, WASM or CLI wrapper for the registry.
- Model licence and checkpoint provenance are the caller's responsibility.

## Next steps

1. Add an optional `tract-onnx` adapter behind a feature flag.
2. Define a fixed ONNX metadata block and a conversion test.
3. Add the neighbour list and PBC cell handling.
4. Add stress and periodic benchmark cases.
5. Add a model registry JSON manifest.
6. Add Python, WASM and CLI wrappers.
