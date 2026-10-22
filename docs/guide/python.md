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

## Conventions

- The canonical name of every capability is the name of the Rust root function
  (`compute_charges`, `solve_eht`, `predict_nmr_shifts`, …). The old bare-noun names
  (`charges`, `sasa`, `population`, `dipole`, `dos`, `esp`, `rmsd`, …) and the `*_calculate`
  forms (`pm3_calculate`, `xtb_calculate`, …) are **deprecated aliases**: they still work but
  emit a `DeprecationWarning` and will be removed in sci-form 2.0.
- Coordinates are flat `[x0, y0, z0, ...]` arrays in **Å**.
- Partial charges are in **elementary charge (e)**, dipole moments in **Debye**, IR/NMR
  frequencies in **cm⁻¹**/**ppm**, J-couplings in **Hz** and UFF/MMFF94 energies in
  **kcal/mol**. PM3/xTB energies are in **eV**; HF-3c and ANI energies are in **Hartree**.
- `ValueError` means malformed input, `RuntimeError` means the engine failed, and
  `NotConvergedError` (a subclass of `RuntimeError`) means an iterative engine did not converge.

See the [Python API reference](/api/python) for the complete surface.

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

### `embed_coords_into(smiles, coords, seed=42) → int` — zero copy

Embed straight into a caller-owned NumPy buffer instead of copying the coordinates across the
Python boundary.

```python
import numpy as np, sci_form

buf = np.zeros(3 * sci_form.parse("CCO").num_atoms, dtype=np.float64)
n_atoms = sci_form.embed_coords_into("CCO", buf, seed=42)
print(n_atoms, buf.reshape(-1, 3)[:2])
```

### `parse(smiles) → ParseResult`

Parse a SMILES string into a molecular structure (no 3D generation). Returns a `ParseResult`
object with attribute access.

```python
info = sci_form.parse("c1ccccc1")
print(f"Atoms: {info.num_atoms}, Bonds: {info.num_bonds}")

for atom in info.atoms:
    print(f"  Z={atom.element}, hybrid={atom.hybridization}, charge={atom.formal_charge}")
```

### `version() → str`

```python
print(sci_form.version())         # "sci-form 1.0.0"
print(sci_form.__version__)       # "1.0.0"
```

## ConformerResult

| Attribute | Type | Description |
|-----------|------|-------------|
| `smiles` | `str` | Input SMILES string |
| `num_atoms` | `int` | Number of atoms |
| `coords` | `list[float]` | Flat coordinates [x₀, y₀, z₀, ...] in Å |
| `elements` | `list[int]` | Atomic numbers |
| `bonds` | `list[tuple]` | Bond list as (a, b, order) |
| `error` | `str \| None` | Error message or None |
| `time_ms` | `float` | Generation time in ms |

### Methods

| Method | Returns | Description |
|--------|---------|-------------|
| `get_positions()` | `list[tuple[float, float, float]]` | Coordinates as (x, y, z) tuples in Å |
| `is_ok()` | `bool` | True if no error |

## ParseResult

| Attribute | Type | Description |
|-----------|------|-------------|
| `num_atoms` | `int` | Number of atoms |
| `num_bonds` | `int` | Number of bonds |
| `atoms` | `list[ParsedAtom]` | Per-atom records |

`ParsedAtom` has `element` (int), `hybridization` (str), `formal_charge` (int) and
`explicit_h` (int).

## Molecular Properties & Analysis

### Gasteiger-Marsili partial charges

```python
import sci_form

result = sci_form.compute_charges("CCO")
print(result.charges)        # [-0.0418, 0.0414, -0.3953, 0.0252, ...]  (e)
print(result.total_charge)   # ~0.0
print(result.iterations)     # 6
```

### Extended Hückel Theory (EHT)

Compute electronic structure and orbital energies:

```python
import sci_form

result = sci_form.embed("CCO")
elements = result.elements
coords = result.coords

eht = sci_form.solve_eht(elements, coords, k=1.75)
print(f"HOMO energy: {eht.homo_energy:.3f} eV")
print(f"LUMO energy: {eht.lumo_energy:.3f} eV")
print(f"Gap: {eht.gap:.3f} eV")
print(f"All orbital energies: {eht.energies}")
```

### Population analysis (Mulliken & Löwdin)

```python
pop = sci_form.compute_population(elements, coords)
print(f"Mulliken charges: {pop.mulliken_charges}")   # e
print(f"Löwdin charges: {pop.lowdin_charges}")       # e
print(f"HOMO energy: {pop.homo_energy:.3f} eV")
print(f"LUMO energy: {pop.lumo_energy:.3f} eV")
print(f"HOMO–LUMO gap: {pop.gap:.3f} eV")
```

### Molecular dipole moment

```python
dipole = sci_form.compute_dipole(elements, coords)
print(f"Dipole magnitude: {dipole.magnitude:.3f} Debye")
print(f"Dipole vector: {dipole.vector} D")   # [x, y, z]
```

### Electrostatic potential (ESP)

```python
esp = sci_form.compute_esp(elements, coords, spacing=0.3, padding=3.0)
print(f"Grid dimensions: {esp.dims}")   # (nx, ny, nz)
print(f"Grid origin: {esp.origin}")     # [x, y, z] in Å
print(f"Grid spacing: {esp.spacing} Å")
print(f"Grid values: {len(esp.values)} (Hartree/e)")
```

### Solvent-accessible surface area (SASA)

```python
sasa = sci_form.compute_sasa(elements, coords, probe_radius=1.4)
print(f"Total SASA: {sasa.total_sasa:.2f} Å²")
print(f"Per-atom SASA: {sasa.atom_sasa}")   # [a₀, a₁, ...]
```

### Density of states (DOS)

```python
dos = sci_form.compute_dos(
    elements, coords,
    sigma=0.3,        # Gaussian smearing width (eV)
    e_min=-30.0,      # Energy range (eV)
    e_max=5.0,
    n_points=500
)

print(f"Total DOS curve: {dos.total_dos}")    # 1/eV
print(f"Energy grid: {dos.energies}")         # eV
print(f"Smearing width: {dos.sigma} eV")
```

For a DOS built from another method use `compute_dos_multimethod(elements, coords, method="xtb")`.

### Molecular alignment & RMSD

```python
coords_flat = result.coords
reference_flat = [...]                      # reference conformation

rmsd_result = sci_form.compute_rmsd(coords_flat, reference_flat)
print(f"RMSD: {rmsd_result.rmsd:.3f} Å")
print(f"Aligned coordinates: {rmsd_result.aligned_coords}")
```

### Force field energy (UFF / MMFF94)

```python
energy = sci_form.compute_uff_energy("CCO", coords_flat)
print(f"UFF energy: {energy:.3f} kcal/mol")

mmff = sci_form.compute_mmff94_energy("CCO", coords_flat)
print(f"MMFF94 energy: {mmff:.3f} kcal/mol")
```

## Parallel processing

All compute functions automatically use the **parallel rayon thread pool** when available
(default build). This provides:

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
print(coords.shape)   # (12, 3)
```

## Quantum chemistry

```python
from sci_form import embed, compute_pm3, compute_xtb, compute_gfn1, compute_gfn2, compute_hf3c

conf = embed("CCO", seed=42)

# PM3 semi-empirical (supports transition metals via PM3(tm))
pm3 = compute_pm3(conf.elements, conf.coords)
print(f"PM3 HOF: {pm3.heat_of_formation:.2f} kcal/mol, gap: {pm3.gap:.3f} eV")
print(f"PM3 total energy: {pm3.total_energy:.3f} eV")

# GFN-xTB
xtb = compute_xtb(conf.elements, conf.coords)
print(f"xTB gap: {xtb.gap:.3f} eV, converged: {xtb.converged}")

# HF-3c — energies are in Hartree
hf = compute_hf3c(conf.elements, conf.coords)
print(f"HF-3c energy: {hf.energy:.4f} Hartree")
```

## Spectroscopy

```python
from sci_form import embed, compute_ir_spectrum, compute_stda_uvvis, compute_nmr_spectrum

conf = embed("CCO", seed=42)

# IR spectrum — wavenumbers in cm⁻¹
ir = compute_ir_spectrum(conf.elements, conf.coords, method="pm3")
print(f"IR peaks: {len(ir.peaks)}")
print(f"Dominant band: {ir.wavenumbers[0]:.1f} cm⁻¹")

# UV-Vis — energies in eV, wavelengths in nm
uv = compute_stda_uvvis(conf.elements, conf.coords)

# NMR — shifts in ppm
nmr = compute_nmr_spectrum("CCO", nucleus="1H")
print(f"NMR peaks: {len(nmr.peaks)}")
```

## Reactivity & stereochemistry

```python
from sci_form import analyze_stereo, compute_fukui_descriptors, compute_empirical_pka

# Stereochemistry (R/S, E/Z, helical, atropisomeric)
stereo = analyze_stereo("C(F)(Cl)(Br)I")
print(f"Stereocenters: {stereo.n_stereocenters}")

# Fukui functions
fukui = compute_fukui_descriptors(conf.elements, conf.coords)

# Empirical pKa
pka = compute_empirical_pka("CC(=O)O")
```

## Error handling

```python
import sci_form

try:
    sci_form.parse("C(")
except ValueError as e:
    print("bad input:", e)
# bad input: Invalid 'smiles' 'C('. Expected ')' at position 2 (found 'end of input') in "C("

try:
    sci_form.compute_hf3c(sci_form.embed("CCO").elements, [0.0, 0.0, 0.0])
except sci_form.NotConvergedError as e:
    print("did not converge:", e)
except RuntimeError as e:
    print("engine failure:", e)
```

→ See the [Python API reference](/api/python) for all functions.
