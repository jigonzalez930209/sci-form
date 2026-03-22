# Python API

## Installation

The package is published as **`sciforma`** on PyPI. The import name inside Python is `sci_form`.

```bash
pip install sciforma
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
print(sci_form.version())  # "sci-form 0.9.1"
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

## Molecular Properties & Analysis

### Gasteiger-Marsili Partial Charges

```python
import sci_form

result = sci_form.compute_charges("CCO")
print(f"Charges: {result['charges']}")          # [q₀, q₁, ...]
print(f"Total charge: {result['total_charge']}")
print(f"Iterations: {result['iterations']}")
```

### Extended Hückel Theory (EHT)

Compute electronic structure and orbital energies:

```python
import sci_form

result = sci_form.embed("CCO")
elements = result.elements
coords = result.get_positions()

eht = sci_form.eht_calculate(elements, coords, k=1.75)
print(f"HOMO energy: {eht['homo_energy']:.3f} eV")
print(f"LUMO energy: {eht['lumo_energy']:.3f} eV")
print(f"Gap: {eht['gap']:.3f} eV")
print(f"All orbital energies: {eht['energies']}")
```

### Population Analysis (Mulliken & Löwdin)

```python
import sci_form

pop = sci_form.compute_population(elements, coords)
print(f"Mulliken charges: {pop['mulliken_charges']}")
print(f"Löwdin charges: {pop['lowdin_charges']}")
print(f"HOMO energy: {pop['homo_energy']:.3f} eV")
print(f"LUMO energy: {pop['lumo_energy']:.3f} eV")
print(f"HOMO–LUMO gap: {pop['homo_lumo_gap']:.3f} eV")
```

### Molecular Dipole Moment

```python
import sci_form

dipole = sci_form.compute_dipole(elements, coords)
print(f"Dipole magnitude: {dipole['magnitude']:.3f} Debye")
print(f"Dipole vector: {dipole['vector']} D")  # [x, y, z]
```

### Electrostatic Potential (ESP)

```python
import sci_form

esp = sci_form.compute_esp(elements, coords, spacing=0.3, padding=3.0)
print(f"Grid dimensions: {esp['dims']}")   # (nx, ny, nz)
print(f"Grid origin: {esp['origin']}")     # [x, y, z]
print(f"Grid spacing: {esp['spacing']} Ų")
print(f"Grid values shape: ({esp['nx']} × {esp['ny']} × {esp['nz']})")
```

### Solvent-Accessible Surface Area (SASA)

```python
import sci_form

sasa = sci_form.compute_sasa(elements, coords, probe_radius=1.4)
print(f"Total SASA: {sasa['total_sasa']:.2f} Ų")
print(f"Per-atom SASA: {sasa['per_atom_sasa']}")  # [a₀, a₁, ...]
```

### Density of States (DOS) & Projected DOS

```python
import sci_form

dos = sci_form.compute_dos(
    elements, coords,
    sigma=0.3,        # Gaussian smearing width
    e_min=-30.0,      # Energy range (eV)
    e_max=5.0,
    n_points=500
)

print(f"HOMO–LUMO gap: {dos['homo_lumo_gap']:.3f} eV")
print(f"Total DOS curve: {dos['total_dos']}")      # [dos₀, dos₁, ...]
print(f"Projected DOS (per-orbital): {dos['pdos']}")
print(f"Energy grid: {dos['energies']}")
print(f"Orbital energies: {dos['orbital_energies']}")
```

### Molecular Alignment & RMSD

```python
import sci_form

coords_flat = [x0, y0, z0, x1, y1, z1, ...]  # from result.coords
reference_flat = [...]                         # reference conformation

rmsd_result = sci_form.compute_rmsd(coords_flat, reference_flat)
print(f"RMSD: {rmsd_result['rmsd']:.3f} Ų")
print(f"Aligned coordinates: {rmsd_result['aligned_coords']}")
```

### Force Field Energy (UFF)

```python
import sci_form

# From SMILES directly
energy = sci_form.compute_uff_energy("CCO", coords_flat)
print(f"UFF Energy: {energy:.3f} kcal/mol")
```

## Parallel Processing

All compute functions automatically use **parallel rayon thread pool** when available (default build). This provides:

- **Automatic parallelization** — grid evaluation, per-atom loops, force field terms all use `rayon::par_iter()`
- **CPU-aware scheduling** — work-stealing queue distributes tasks across all cores
- **Zero configuration** — no thread pool setup required

Example:

```python
import sci_form
import time

# Large molecule with 100+ atoms
result = sci_form.embed("c1ccc(cc1)C(=O)Nc1ccc(O)cc1" * 2, seed=42)

# These run in parallel automatically
start = time.time()
dos = sci_form.compute_dos(result.elements, result.coords, sigma=0.3, n_points=500)
print(f"DOS computed in {time.time() - start:.3f}s (parallel)")

sasa = sci_form.compute_sasa(result.elements, result.coords)
print(f"SASA computed in {time.time() - start:.3f}s (parallel)")
```

To disable parallelization, build from source without the `parallel` feature:

```bash
cd crates/python
pip install maturin
maturin develop --release --no-default-features
```

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

## Quantum Chemistry

```python
from sci_form import embed, pm3_calculate, xtb_calculate

conf = embed("CCO", seed=42)

# PM3 semi-empirical (supports transition metals via PM3(tm))
pm3 = pm3_calculate(conf.elements, conf.coords)
print(f"PM3 HOF: {pm3.heat_of_formation:.2f} kcal/mol, gap: {pm3.gap:.3f} eV")

# GFN-xTB
xtb = xtb_calculate(conf.elements, conf.coords)
print(f"xTB gap: {xtb.gap:.3f} eV, converged: {xtb.converged}")

# HF-3c
from sci_form import hf3c_calculate
hf = hf3c_calculate(conf.elements, conf.coords)
```

## Spectroscopy

```python
from sci_form import ir_spectrum, stda_uvvis, nmr_spectrum

conf = embed("CCO", seed=42)

# IR spectrum
ir = ir_spectrum(conf.elements, conf.coords, method="pm3")
print(f"IR peaks: {len(ir.peaks)}")

# UV-Vis
uv = stda_uvvis(conf.elements, conf.coords)

# NMR
nmr = nmr_spectrum("CCO", nucleus="1H")
```

## Reactivity & Stereochemistry

```python
from sci_form import stereo_analysis, fukui_descriptors, empirical_pka

# Stereochemistry (R/S, E/Z, helical, atropisomeric)
stereo = stereo_analysis("C(F)(Cl)(Br)I")
print(f"Stereocenters: {stereo.n_stereocenters}")

# Fukui functions
fukui = fukui_descriptors(conf.elements, conf.coords)

# Empirical pKa
pka = empirical_pka("CC(=O)O")
```

→ See the [Python API reference](/api/python) for all functions.
