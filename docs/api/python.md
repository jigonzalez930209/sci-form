# Python API Reference

## Installation

```bash
pip install sciforma        # package name on PyPI
import sci_form             # import name
```

---

## Conformer Generation

### `embed`

```python
def embed(smiles: str, seed: int = 42) -> ConformerResult
```

Generate a 3D conformer for a single molecule via ETKDGv2.

```python
from sci_form import embed

result = embed("CCO", seed=42)
print(result)  # ConformerResult(smiles='CCO', atoms=9, time=8.3ms)

if result.is_ok():
    positions = result.get_positions()  # list of (x, y, z) tuples
    print(f"O atom: {positions[2]}")
```

### `embed_batch`

```python
def embed_batch(
    smiles_list: list[str],
    seed: int = 42,
    num_threads: int = 0,
) -> list[ConformerResult]
```

Batch-generate conformers. `num_threads=0` auto-detects CPU count.

```python
from sci_form import embed_batch

results = embed_batch(["CCO", "c1ccccc1", "CC(=O)O"], seed=42)
for r in results:
    if r.is_ok():
        print(f"{r.smiles}: {r.num_atoms} atoms, {r.time_ms:.1f}ms")
```

### `parse`

```python
def parse(smiles: str) -> dict
```

Parse SMILES to a molecular graph without generating 3D coordinates.

```python
from sci_form import parse

mol = parse("CCO")
print(mol["num_atoms"])   # 9
print(mol["atoms"][0])    # {'element': 8, 'hybridization': 'SP3', 'formal_charge': 0}
```

---

## Gasteiger Charges

### `charges`

```python
def charges(smiles: str) -> ChargeResult
```

Compute Gasteiger-Marsili partial charges (6 iterations).

```python
from sci_form import charges

result = charges("CCO")
print(result)  # ChargeResult(n_atoms=9, total_charge=0.0000, iterations=6)
print(result.charges)      # [-0.387, -0.042, -0.228, 0.123, ...]
print(result.total_charge) # ~0.0
```

---

## SASA

### `sasa`

```python
def sasa(
    elements: list[int],
    coords: list[float],
    probe_radius: float = 1.4,
) -> SasaResult
```

Compute solvent-accessible surface area (Shrake-Rupley).

```python
from sci_form import embed, sasa

conf = embed("CCO")
result = sasa(conf.elements, conf.coords, probe_radius=1.4)
print(result)             # SasaResult(total=82.31 Ų, n_atoms=9, probe=1.40 Å)
print(result.total_sasa)  # Å²
print(result.atom_sasa)   # per-atom list
```

---

## EHT — Extended Hückel Theory

### `eht_calculate`

```python
def eht_calculate(
    elements: list[int],
    coords: list[float],
    k: float = 1.75,
) -> EhtResult
```

Run an EHT semi-empirical calculation. `k` is the Wolfsberg-Helmholtz constant (default 1.75).

```python
from sci_form import embed, eht_calculate

conf = embed("CCO")
eht = eht_calculate(conf.elements, conf.coords)
print(eht)              # EhtResult(n_mo=9, gap=7.831 eV, HOMO=-12.453 eV, LUMO=-4.622 eV)
print(eht.homo_energy)  # eV
print(eht.gap)          # HOMO-LUMO gap in eV
print(eht.energies)     # all orbital energies (eV)
```

### `eht_orbital_mesh`

```python
def eht_orbital_mesh(
    elements: list[int],
    coords: list[float],
    mo_index: int,
    spacing: float = 0.2,
    isovalue: float = 0.02,
) -> dict
```

Generate a 3D isosurface mesh for a molecular orbital at the given isovalue.

```python
from sci_form import embed, eht_orbital_mesh

conf = embed("CCO")
eht = eht_calculate(conf.elements, conf.coords)

# HOMO mesh
homo_idx = eht.homo_index
mesh = eht_orbital_mesh(conf.elements, conf.coords, homo_idx, isovalue=0.02)
print(mesh.keys())  # dict_keys(['vertices', 'normals', 'indices', 'num_triangles', 'isovalue'])

# Use in Three.js / Open3D / matplotlib
vertices = mesh["vertices"]   # flat [x,y,z, ...] array
triangles = mesh["indices"]   # flat [i0,i1,i2, ...] triangle index array
```

---

## Population Analysis

### `population`

```python
def population(
    elements: list[int],
    coords: list[float],
) -> PopulationResult
```

Mulliken and Löwdin population analysis from EHT.

```python
from sci_form import embed, population

conf = embed("CCO")
pop = population(conf.elements, conf.coords)
print(pop)                       # PopulationResult(n_atoms=9, total_mulliken=0.0000)
print(pop.mulliken_charges)      # per-atom Mulliken charges
print(pop.lowdin_charges)        # per-atom Löwdin charges
print(pop.total_charge_mulliken) # sum ≈ 0 for neutral
```

---

## Dipole Moment

### `dipole`

```python
def dipole(
    elements: list[int],
    coords: list[float],
) -> DipoleResult
```

Compute molecular dipole moment in Debye.

```python
from sci_form import embed, dipole

conf = embed("CCO")
d = dipole(conf.elements, conf.coords)
print(d)              # DipoleResult(1.847 Debye)
print(d.magnitude)    # 1.847 (Debye)
print(d.vector)       # [x, y, z] in Debye
print(d.unit)         # "Debye"
```

---

## Density of States

### `dos`

```python
def dos(
    elements: list[int],
    coords: list[float],
    sigma: float = 0.3,
    e_min: float = -30.0,
    e_max: float = 5.0,
    n_points: int = 500,
) -> DosResult
```

Compute total DOS from EHT orbital energies with Gaussian smearing $\sigma$.

```python
from sci_form import embed, dos
import matplotlib.pyplot as plt

conf = embed("c1ccccc1")
result = dos(conf.elements, conf.coords, sigma=0.2, e_min=-25.0, e_max=5.0, n_points=300)

plt.plot(result.energies, result.total_dos)
plt.xlabel("Energy (eV)")
plt.ylabel("DOS (states/eV)")
plt.show()
```

---

## Molecular Alignment

### `rmsd`

```python
def rmsd(
    coords: list[float],
    reference: list[float],
) -> AlignmentResult
```

Kabsch optimal alignment and RMSD computation.

```python
from sci_form import embed, rmsd

conf1 = embed("CCO", seed=42)
conf2 = embed("CCO", seed=123)
result = rmsd(conf1.coords, conf2.coords)
print(result)              # AlignmentResult(rmsd=0.0342 Å)
print(result.aligned_coords)  # conf1 coords after optimal rotation onto conf2
```

---

## Force Field Energy

### `uff_energy`

```python
def uff_energy(smiles: str, coords: list[float]) -> float
```

Evaluate UFF force field energy in kcal/mol.

```python
from sci_form import embed, uff_energy

conf = embed("CCO")
e = uff_energy("CCO", conf.coords)
print(f"UFF energy: {e:.3f} kcal/mol")
```

---

## Materials

### `unit_cell`

```python
def unit_cell(
    a: float, b: float, c: float,
    alpha: float = 90.0,
    beta: float = 90.0,
    gamma: float = 90.0,
) -> UnitCell
```

Create a periodic unit cell from crystallographic parameters (a, b, c in Å; angles in degrees).

```python
from sci_form import unit_cell

cell = unit_cell(10.0, 10.0, 10.0)           # cubic
cell = unit_cell(5.0, 5.0, 12.0, 90, 90, 120)  # hexagonal
print(cell.volume)   # Å³
print(cell.lattice)  # 3×3 matrix
```

### `assemble_framework`

```python
def assemble_framework(
    topology: str = "pcu",
    metal: int = 30,
    geometry: str = "octahedral",
    lattice_a: float = 10.0,
    supercell: int = 1,
) -> CrystalStructure
```

Assemble a MOF-type framework crystal structure. `topology` ∈ {"pcu", "dia", "sql"}. `geometry` ∈ {"linear", "trigonal", "tetrahedral", "square_planar", "octahedral"}.

```python
from sci_form import assemble_framework

mof = assemble_framework(topology="pcu", metal=30, geometry="octahedral", lattice_a=26.3)
print(mof.num_atoms)    # total atoms
print(mof.elements)     # atomic numbers
print(mof.cart_coords)  # Cartesian coordinates
print(mof.frac_coords)  # fractional coordinates
```

---

## Spectroscopy (Track D)

### `stda_uvvis`

```python
def stda_uvvis(
    elements: list[int],
    coords: list[float],
    sigma: float = 0.3,
    e_min: float = 1.0,
    e_max: float = 8.0,
    n_points: int = 500,
    broadening: str = "gaussian",  # "gaussian" or "lorentzian"
) -> StdaUvVisSpectrumPy
```

Compute UV-Vis spectrum via sTDA on EHT MO transitions. Returns a `StdaUvVisSpectrumPy` with `.wavelengths_nm`, `.intensities`, `.excitations` (list of `StdaExcitationPy`), `.gap`, `.homo_energy`, `.lumo_energy`, `.n_transitions`, `.notes`.

```python
from sci_form import embed, stda_uvvis

conf = embed("c1ccccc1", seed=42)
spec = stda_uvvis(conf.elements, conf.coords, sigma=0.3, e_min=1.0, e_max=8.0)
print(f"Gap: {spec.gap:.2f} eV")
max_i  = max(spec.intensities)
lam_max = spec.wavelengths_nm[spec.intensities.index(max_i)]
print(f"λ_max ≈ {lam_max:.1f} nm")
for ex in sorted(spec.excitations, key=lambda e: -e.oscillator_strength)[:3]:
    print(f"  {ex.wavelength_nm:.1f} nm  f={ex.oscillator_strength:.4f}  MO{ex.from_mo}→MO{ex.to_mo}")
```

### `vibrational_analysis`

```python
def vibrational_analysis(
    elements: list[int],
    coords: list[float],
    method: str = "eht",   # "eht", "pm3", or "xtb"
    step_size: float = 0.01,
) -> VibrationalAnalysisPy
```

Numerical Hessian ($6N$ energy evaluations) + diagonalization. Returns `VibrationalAnalysisPy` with `.n_atoms`, `.modes` (list of `VibrationalModePy`), `.n_real_modes`, `.zpve_ev`, `.method`, `.notes`.

Each `VibrationalModePy` has `.frequency_cm1`, `.ir_intensity` (km/mol), `.displacement` (list, length 3N), `.is_real`.

```python
from sci_form import embed, vibrational_analysis

conf = embed("CCO", seed=42)
vib  = vibrational_analysis(conf.elements, conf.coords, method="xtb")
print(f"ZPVE: {vib.zpve_ev:.4f} eV, {vib.n_real_modes} real modes")

real = [m for m in vib.modes if m.is_real]
strongest = max(real, key=lambda m: m.ir_intensity)
print(f"Strongest IR: {strongest.frequency_cm1:.1f} cm⁻¹  ({strongest.ir_intensity:.1f} km/mol)")
```

### `ir_spectrum`

```python
def ir_spectrum(
    analysis: VibrationalAnalysisPy,
    gamma: float = 15.0,
    wn_min: float = 600.0,
    wn_max: float = 4000.0,
    n_points: int = 1000,
) -> IrSpectrumPy
```

Lorentzian-broadened IR spectrum. Returns `IrSpectrumPy` with `.wavenumbers`, `.intensities`, `.peaks` (list of `IrPeakPy`), `.gamma`.

```python
from sci_form import ir_spectrum

spec = ir_spectrum(vib, gamma=15.0, wn_min=600.0, wn_max=4000.0, n_points=1000)
import matplotlib.pyplot as plt
plt.plot(spec.wavenumbers, spec.intensities)
plt.xlabel("Wavenumber (cm⁻¹)")
plt.gca().invert_xaxis()
plt.show()
```

### `nmr_shifts`

```python
def nmr_shifts(smiles: str) -> NmrShiftResultPy
```

Predict ¹H and ¹³C chemical shifts using HOSE code environment matching. Returns `NmrShiftResultPy` with `.h_shifts` and `.c_shifts` (lists of `ChemicalShiftPy`). No 3D geometry needed.

Each `ChemicalShiftPy` has `.atom_index`, `.element`, `.shift_ppm`, `.environment`, `.confidence`.

```python
from sci_form import nmr_shifts

shifts = nmr_shifts("CCO")
for h in shifts.h_shifts:
    print(f"H#{h.atom_index}: {h.shift_ppm:.2f} ppm  [{h.environment}]  conf={h.confidence:.2f}")
for c in shifts.c_shifts:
    print(f"C#{c.atom_index}: {c.shift_ppm:.1f} ppm  [{c.environment}]")
```

### `nmr_couplings`

```python
def nmr_couplings(
    smiles: str,
    coords: list[float] = [],  # flat [x0,y0,z0,...]; empty → free-rotation average
) -> list[JCouplingPy]
```

Predict H–H J-coupling constants. With 3D coords uses the Karplus equation for ³J; topology only gives free-rotation averages (~5.3 Hz). Each `JCouplingPy` has `.h1_index`, `.h2_index`, `.j_hz`, `.n_bonds`, `.coupling_type`.

```python
from sci_form import embed, nmr_couplings

conf = embed("CC", seed=42)
couplings = nmr_couplings("CC", conf.coords)
for j in couplings:
    print(f"{j.n_bonds}J(H{j.h1_index},H{j.h2_index}) = {j.j_hz:.2f} Hz  [{j.coupling_type}]")
```

### `nmr_spectrum`

```python
def nmr_spectrum(
    smiles: str,
    nucleus: str = "1H",   # "1H" or "13C"
    gamma: float = 0.01,
    ppm_min: float = -2.0,
    ppm_max: float = 14.0,
    n_points: int = 2000,
) -> NmrSpectrumPy
```

Full NMR spectrum pipeline. Returns `NmrSpectrumPy` with `.ppm_axis`, `.intensities`, `.peaks` (list of `NmrPeakPy`), `.nucleus`, `.gamma`.

```python
from sci_form import nmr_spectrum

spec = nmr_spectrum("CCO", nucleus="1H")
import matplotlib.pyplot as plt
plt.plot(spec.ppm_axis, spec.intensities)
plt.xlabel("δ (ppm)")
plt.gca().invert_xaxis()
plt.show()
```

### `hose_codes`

```python
def hose_codes(smiles: str, max_radius: int = 4) -> list[HoseCodePy]
```

Generate HOSE code fingerprints for all atoms. Each `HoseCodePy` has `.atom_index`, `.element`, `.code_string`, `.radius`.

---

## Transport / Batch Streaming

### `pack_conformers`

```python
def pack_conformers(results: list[ConformerResult]) -> RecordBatch
```

Pack conformer results into Arrow-compatible columnar format for efficient transfer.

```python
from sci_form import embed_batch, pack_conformers

results = embed_batch(["CCO", "c1ccccc1"])
batch = pack_conformers(results)
print(batch.num_rows)        # 2
print(batch.column_names)    # ['smiles', 'num_atoms', 'coords', ...]
print(batch.float_data)      # {'coords': [...], 'time_ms': [...]}
```

### `split_worker_tasks`

```python
def split_worker_tasks(
    smiles: list[str],
    n_workers: int = 4,
    seed: int = 42,
) -> list[dict]
```

Split a SMILES list into balanced worker tasks for multi-process/web-worker dispatch.

### `estimate_workers`

```python
def estimate_workers(n_items: int, max_workers: int = 8) -> int
```

Estimate the optimal number of workers for a given batch size.

---

## Return Types

### `ConformerResult`

| Attribute | Type | Description |
|-----------|------|-------------|
| `smiles` | `str` | Input SMILES |
| `num_atoms` | `int` | Atom count |
| `coords` | `list[float]` | Flat coords `[x₀,y₀,z₀,…]` |
| `elements` | `list[int]` | Atomic numbers |
| `bonds` | `list[tuple]` | `(a, b, order_str)` |
| `error` | `str \| None` | Error message or `None` |
| `time_ms` | `float` | Generation time ms |

**Methods:** `get_positions() → list[tuple[float,float,float]]`, `is_ok() → bool`

### `EhtResult`

| Attribute | Type | Description |
|-----------|------|-------------|
| `energies` | `list[float]` | All MO energies (eV) |
| `n_electrons` | `int` | Total electron count |
| `homo_index` | `int` | Index of HOMO orbital |
| `lumo_index` | `int` | Index of LUMO orbital |
| `homo_energy` | `float` | HOMO energy (eV) |
| `lumo_energy` | `float` | LUMO energy (eV) |
| `gap` | `float` | HOMO-LUMO gap (eV) |

### `ChargeResult`

| Attribute | Type | Description |
|-----------|------|-------------|
| `charges` | `list[float]` | Per-atom partial charges |
| `iterations` | `int` | Convergence iterations |
| `total_charge` | `float` | Sum of all charges |

### `SasaResult`

| Attribute | Type | Description |
|-----------|------|-------------|
| `total_sasa` | `float` | Total SASA in Å² |
| `atom_sasa` | `list[float]` | Per-atom SASA in Å² |
| `probe_radius` | `float` | Probe radius used (Å) |
| `num_points` | `int` | Fibonacci sphere points per atom |

### `PopulationResult`

| Attribute | Type | Description |
|-----------|------|-------------|
| `mulliken_charges` | `list[float]` | Per-atom Mulliken charges |
| `lowdin_charges` | `list[float]` | Per-atom Löwdin charges |
| `num_atoms` | `int` | Atom count |
| `total_charge_mulliken` | `float` | Total Mulliken charge |

### `DipoleResult`

| Attribute | Type | Description |
|-----------|------|-------------|
| `vector` | `list[float]` | Dipole vector [x, y, z] in Debye |
| `magnitude` | `float` | Dipole magnitude in Debye |
| `unit` | `str` | Always `"Debye"` |

### `DosResult`

| Attribute | Type | Description |
|-----------|------|-------------|
| `energies` | `list[float]` | Energy axis (eV) |
| `total_dos` | `list[float]` | Total DOS curve |
| `sigma` | `float` | Smearing width used (eV) |

### `AlignmentResult`

| Attribute | Type | Description |
|-----------|------|-------------|
| `rmsd` | `float` | RMSD after Kabsch alignment (Å) |
| `aligned_coords` | `list[float]` | Transformed coordinates |

---

## Integration Examples

### NumPy Array

```python
import numpy as np
from sci_form import embed

result = embed("CCO")
coords = np.array(result.coords).reshape(-1, 3)  # (9, 3)
elements = np.array(result.elements)              # (9,)
```

### Pandas DataFrame

```python
import pandas as pd
from sci_form import embed_batch, charges

smiles = ["CCO", "c1ccccc1", "CC(=O)O"]
results = embed_batch(smiles)

df = pd.DataFrame([
    {
        "smiles": r.smiles,
        "num_atoms": r.num_atoms,
        "time_ms": r.time_ms,
        "ok": r.is_ok(),
    }
    for r in results
])
```

### RDKit Interop

```python
from rdkit import Chem
from rdkit.Chem import AllChem
from sci_form import embed

result = embed("c1ccccc1")
mol = Chem.MolFromSmiles("c1ccccc1")
mol = Chem.AddHs(mol)

conf = Chem.Conformer(result.num_atoms)
for i in range(result.num_atoms):
    conf.SetAtomPosition(i, (
        result.coords[i*3],
        result.coords[i*3+1],
        result.coords[i*3+2],
    ))

mol.AddConformer(conf, assignId=True)
```

### DOS + matplotlib

```python
import matplotlib.pyplot as plt
from sci_form import embed, dos

conf = embed("c1ccccc1")   # benzene
result = dos(conf.elements, conf.coords, sigma=0.3)

plt.figure(figsize=(8, 4))
plt.plot(result.energies, result.total_dos, 'b-', lw=1.5)
plt.axvline(x=0, color='k', ls='--', label='Fermi level (HOMO)')
plt.xlabel("Energy (eV)")
plt.ylabel("DOS (states/eV)")
plt.title("Benzene DOS")
plt.legend()
plt.tight_layout()
plt.show()
```

