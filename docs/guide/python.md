# Python API

## Installation

```bash
pip install sci-form
```

Or build from source:

```bash
cd crates/python
pip install maturin
maturin develop --release
```

## Functions

### `embed(smiles, seed=42) → ConformerResult`

Generate a single 3D conformer.

```python
import sci_form

result = sci_form.embed("CCO")
print(f"Atoms: {result.num_atoms}, Time: {result.time_ms:.1f}ms")

# With a specific seed for reproducibility
result = sci_form.embed("c1ccccc1", seed=123)
```

### `embed_batch(smiles_list, seed=42, num_threads=0) → list[ConformerResult]`

Generate conformers for multiple molecules in parallel.

```python
smiles = ["CCO", "c1ccccc1", "CC(=O)O", "CC(C)CC"]
results = sci_form.embed_batch(smiles, num_threads=4)

for r in results:
    if r.is_ok():
        print(f"{r.smiles}: {r.num_atoms} atoms in {r.time_ms:.1f}ms")
    else:
        print(f"{r.smiles}: FAILED - {r.error}")
```

### `parse(smiles) → dict`

Parse a SMILES string into a molecular structure (no 3D generation).

```python
info = sci_form.parse("c1ccccc1")
print(f"Atoms: {info['num_atoms']}, Bonds: {info['num_bonds']}")

for atom in info['atoms']:
    print(f"  Z={atom['element']}, hybrid={atom['hybridization']}, charge={atom['formal_charge']}")
```

Returns a dict:
```python
{
    "num_atoms": 12,
    "num_bonds": 12,
    "atoms": [
        {"element": 6, "hybridization": "SP2", "formal_charge": 0},
        ...
    ]
}
```

### `version() → str`

```python
print(sci_form.version())  # "sci-form 0.1.0"
```

## ConformerResult

| Attribute | Type | Description |
|-----------|------|-------------|
| `smiles` | `str` | Input SMILES string |
| `num_atoms` | `int` | Number of atoms |
| `coords` | `list[float]` | Flat coordinates [x₀, y₀, z₀, ...] |
| `elements` | `list[int]` | Atomic numbers |
| `bonds` | `list[tuple]` | Bond list as (a, b, order) |
| `error` | `str \| None` | Error message or None |
| `time_ms` | `float` | Generation time in ms |

### Methods

| Method | Returns | Description |
|--------|---------|-------------|
| `get_positions()` | `list[tuple[float, float, float]]` | Coordinates as (x, y, z) tuples |
| `is_ok()` | `bool` | True if no error |

## Integration with RDKit

Convert sci-form output to an RDKit molecule:

```python
from rdkit import Chem
from rdkit.Chem import AllChem
import sci_form

result = sci_form.embed("CCO")
positions = result.get_positions()

# Create RDKit molecule
mol = Chem.MolFromSmiles("CCO")
mol = Chem.AddHs(mol)
conf = Chem.Conformer(mol.GetNumAtoms())
for i, (x, y, z) in enumerate(positions):
    conf.SetAtomPosition(i, (x, y, z))
mol.AddConformer(conf, assignId=True)
```

## Integration with NumPy

```python
import numpy as np
import sci_form

result = sci_form.embed("c1ccccc1")
coords = np.array(result.coords).reshape(-1, 3)
print(coords.shape)  # (12, 3)
```
