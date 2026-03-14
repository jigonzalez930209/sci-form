# Python API Reference

## Installation

```bash
pip install sci-form
```

## Functions

### `embed`

```python
def embed(smiles: str, seed: int = 42) -> ConformerResult
```

Generate a 3D conformer for a single molecule.

```python
from sci_form import embed

result = embed("CCO", seed=42)
print(result.num_atoms)   # 9
print(result.coords[:3])  # [x₀, y₀, z₀]
```

### `embed_batch`

```python
def embed_batch(
    smiles_list: list[str],
    seed: int = 42,
    num_threads: int = 0
) -> list[ConformerResult]
```

Generate conformers for multiple molecules. Uses Rust-level parallelism when `num_threads > 1` (or `0` for auto-detect).

```python
from sci_form import embed_batch

results = embed_batch(["CCO", "c1ccccc1", "CC(=O)O"], seed=42)
for r in results:
    print(f"{r.smiles}: {r.num_atoms} atoms, {r.time_ms:.1f}ms")
```

### `parse`

```python
def parse(smiles: str) -> dict
```

Parse SMILES into a molecular graph representation.

```python
from sci_form import parse

mol = parse("CCO")
# Returns dict with atoms, bonds, ring info
```

### `version`

```python
def version() -> str
```

Returns version string, e.g. `"sci-form 0.1.0"`.

## `ConformerResult`

Returned by `embed()` and `embed_batch()`.

### Attributes

| Attribute | Type | Description |
|-----------|------|-------------|
| `smiles` | `str` | Input SMILES string |
| `num_atoms` | `int` | Total atom count including implicit H |
| `coords` | `list[float]` | Flat 3D coordinates `[x₀, y₀, z₀, ...]` |
| `elements` | `list[int]` | Atomic numbers |
| `bonds` | `list[tuple]` | `(atom_a, atom_b, order_string)` |
| `error` | `str \| None` | Error message, or `None` on success |
| `time_ms` | `float` | Generation time in milliseconds |

### Methods

#### `to_xyz()`

```python
def to_xyz(self) -> str
```

Returns XYZ-format string:
```
9
CCO
C    0.123   0.456   0.789
...
```

#### `to_dict()`

```python
def to_dict(self) -> dict
```

Serializes result to a Python dictionary.

## Integration Examples

### RDKit Conversion

```python
from rdkit import Chem
from rdkit.Chem import AllChem
from sci_form import embed

result = embed("c1ccccc1")
mol = Chem.MolFromSmiles("c1ccccc1")
mol = Chem.AddHs(mol)

conf = Chem.Conformer(result.num_atoms)
for i in range(result.num_atoms):
    x = result.coords[i * 3]
    y = result.coords[i * 3 + 1]
    z = result.coords[i * 3 + 2]
    conf.SetAtomPosition(i, (x, y, z))

mol.AddConformer(conf, assignId=True)
```

### NumPy Array

```python
import numpy as np
from sci_form import embed

result = embed("CCO")
coords = np.array(result.coords).reshape(-1, 3)
# coords.shape = (9, 3)
```

### Pandas DataFrame

```python
import pandas as pd
from sci_form import embed_batch

smiles = ["CCO", "c1ccccc1", "CC(=O)O"]
results = embed_batch(smiles)

df = pd.DataFrame([r.to_dict() for r in results])
print(df[["smiles", "num_atoms", "time_ms"]])
```
