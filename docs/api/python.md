# Python API Reference

The Python package is published as **`sciforma`** on PyPI; the import name is `sci_form`.

```bash
pip install sciforma
```

```python
import sci_form

print(sci_form.__version__)   # "1.0.0"
print(sci_form.version())      # "sci-form 1.0.0"
```

---

## Naming convention

**The canonical name of every capability is the name of the Rust root function**, and every
action carries an explicit verb (`compute_*`, `analyze_*`, `predict_*`). The old Python names —
bare nouns such as `charges`, `sasa`, `population`, `dipole`, `dos`, and the `*_calculate` forms
such as `pm3_calculate`, `xtb_calculate`, `gfn1_calculate`, `gfn2_calculate`, `hf3c_calculate`,
`ani_calculate` — are **deprecated aliases**. They still work, but they emit a
`DeprecationWarning` and will be removed in sci-form 2.0:

```python
>>> import warnings, sci_form
>>> with warnings.catch_warnings(record=True) as w:
...     warnings.simplefilter("always")
...     sci_form.charges("CCO")
...     print(w[0].category.__name__, w[0].message)
DeprecationWarning sci_form.charges is deprecated and will be removed in sci-form 2.0; use sci_form.compute_charges instead.
```

Every deprecated alias is listed next to its canonical function below, and the complete table
is in [`crates/python/src/naming.rs`](https://github.com/jigonzalez930209/sci-form/blob/main/crates/python/src/naming.rs).

The module also exposes `sci_form.__all__` — the grouped public surface — so
`dir(sci_form)` and `help(sci_form)` are meaningful:

```python
>>> len(sci_form.__all__)
247
```

The list holds 145 callables (93 canonical functions plus 52 deprecated aliases), 98 result
classes, the three submodules `alpha`, `beta` and `potentials`, and the `__version__` string.

---

## Units

Units are **not** uniform across the library. Check the unit of every returned number:

| Quantity | Unit |
|---|---|
| Geometry (coordinates, bond lengths, cell parameters) | Å (cell angles in degrees) |
| UFF / MMFF94 energies | kcal/mol |
| UFF / MMFF94 gradients | kcal/mol/Å |
| EHT, PM3, xTB orbital energies and HOMO–LUMO gaps | eV |
| PM3 / xTB total, electronic, repulsive, dispersion energies | eV |
| PM3 heat of formation | kcal/mol |
| HF-3c energy, HF energy, nuclear repulsion, D3/gCP/SRB corrections | Hartree |
| HF-3c orbital energies, HOMO/LUMO, gap | Hartree |
| ANI energy and per-atom energies | Hartree |
| ANI forces | Hartree/Å |
| CISD excitation `energy`, `correlation_energy` | Hartree (`energy_ev` is eV) |
| EHT gradients, `optimize_geometry_eht` energy | eV, eV/Å |
| `optimize_framework` energy and forces | eV, eV/Å |
| Vibrational frequencies, IR wavenumbers | cm⁻¹ |
| NMR chemical shifts and shieldings | ppm |
| NMR J-couplings | Hz |
| Dipole moment | Debye |
| Partial charges | elementary charge (e) |
| Electrostatic potential | Hartree/e |
| DOS / PDOS | 1/eV |
| Surface area | Å² |
| MD energies | kcal/mol; times in fs; temperatures in K |
| Pressure | bar |

Energies are always **total system energies**, not per-atom normalised, unless the field name
says otherwise. Gradients are `dE/dR` with `R` in Å, so they carry the sign convention
"force = −gradient".

---

## Errors

| Situation | Exception |
|---|---|
| The caller's input is malformed (bad SMILES, wrong coordinate count, unknown method) | `ValueError` |
| The engine ran but failed (SCF blew up, optimisation failed) | `RuntimeError` |
| The engine ran but did not converge in the allowed number of iterations | `NotConvergedError` |

`NotConvergedError` is a **subclass of `RuntimeError`**, so `except RuntimeError` still catches
it while `except NotConvergedError` lets you separate "the physics did not settle" from "you
passed nonsense".

```python
import sci_form

try:
    sci_form.parse("C(")
except ValueError as e:
    print(e)
# Invalid 'smiles' 'C('. Expected ')' at position 2 (found 'end of input') in "C("
```

Messages follow one template everywhere:

```text
<What failed> '<actual value>'. Supported: 'a', 'b', 'c'.
```

Rust internal names (`positions`, `unwrap`, `None`, a `file:line`) never leak into a Python
message, and a NaN is never returned.

---

## Package metadata

| Name | Type | Description |
|---|---|---|
| `sci_form.__version__` | `str` | Crate version, e.g. `"1.0.0"` |
| `version` | `str` | `"sci-form 1.0.0"` (also available as `sci_form.version()`) |
| `sci_form.NotConvergedError` | exception | Subclass of `RuntimeError` |

---

## Conformer generation

### `embed`

```python
def embed(smiles: str, seed: int = 42) -> ConformerResult
```

Generate a 3D conformer for one molecule with ETKDGv2. Coordinates are flat and in **Å**.

```python
from sci_form import embed

result = embed("CCO", seed=42)
print(result)                # ConformerResult(smiles='CCO', atoms=9, time=1.4ms)
print(result.num_atoms)      # 9
print(result.get_positions() # [(x, y, z), ...] in Å
print(result.is_ok())        # True
```

### `embed_batch`

```python
def embed_batch(smiles_list: list[str], seed: int = 42, num_threads: int = 0) -> list[ConformerResult]
```

Batch-generate conformers. `num_threads=0` auto-detects the CPU count.

```python
results = embed_batch(["CCO", "c1ccccc1", "CC(=O)O"], seed=42)
for r in results:
    if r.is_ok():
        print(f"{r.smiles}: {r.num_atoms} atoms, {r.time_ms:.1f}ms")
```

### `embed_diverse`

```python
def embed_diverse(
    smiles: str,
    n_conformers: int = 10,
    rmsd_cutoff: float = 1.0,
    base_seed: int = 42,
) -> ConformerSearchResult
```

Embed `n_conformers` times with consecutive seeds, keep only the Butina RMSD cluster centroids
(cutoff in **Å**) and report the UFF energy of every member (**kcal/mol**).

```python
result = embed_diverse("CCCC", n_conformers=10, rmsd_cutoff=1.0)
print(result.generated, result.unique, len(result.clusters))
for member in result.conformers:
    print(member.seed, member.cluster_id, member.energy_kcal_mol)
```

### `embed_coords_into` — zero copy

```python
def embed_coords_into(smiles: str, coords: "numpy.ndarray", seed: int = 42) -> int
```

Embed a SMILES and write the coordinates **straight into a caller-owned buffer** instead of
copying them across the Python boundary. `coords` must be a writable, C-contiguous
`float64` array of length `3 * num_atoms`; otherwise a `ValueError` is raised and the buffer
is left untouched. Returns the number of atoms written.

```python
import numpy as np, sci_form

buf = np.zeros(3 * sci_form.parse("CCO").num_atoms, dtype=np.float64)
n_atoms = sci_form.embed_coords_into("CCO", buf, seed=42)
print(n_atoms, buf[:3])
# 9 [-0.89306372  0.06926389  0.06360101]
```

### `parse`

```python
def parse(smiles: str) -> ParseResult
```

Parse a SMILES string into a molecular graph **without** generating 3D coordinates. Returns a
`ParseResult` object (it used to return a plain `dict`).

```python
mol = sci_form.parse("c1ccccc1")
print(mol.num_atoms, mol.num_bonds)     # 12 12
for atom in mol.atoms:
    print(atom.element, atom.hybridization, atom.formal_charge, atom.explicit_h)
```

### `to_canonical_smiles`

```python
def to_canonical_smiles(smiles: str) -> str
```

```python
print(sci_form.to_canonical_smiles("C(F)(Cl)(Br)I"))  # BrC(Cl)(F)I
```

### Return types

#### `ConformerResult`

| Attribute | Type | Description |
|---|---|---|
| `smiles` | `str` | Input SMILES |
| `num_atoms` | `int` | Atom count (implicit H included) |
| `coords` | `list[float]` | Flat coordinates `[x₀,y₀,z₀,…]` in Å |
| `elements` | `list[int]` | Atomic numbers |
| `bonds` | `list[tuple]` | `(a, b, order_str)` with order ∈ `SINGLE/DOUBLE/TRIPLE/AROMATIC` |
| `error` | `str \| None` | Error message or `None` |
| `time_ms` | `float` | Generation time in ms |

Methods: `get_positions() -> list[tuple[float, float, float]]`, `is_ok() -> bool`.

#### `ParseResult`

| Attribute | Type | Description |
|---|---|---|
| `num_atoms` | `int` | Atom count |
| `num_bonds` | `int` | Bond count |
| `atoms` | `list[ParsedAtom]` | Per-atom records |

#### `ParsedAtom`

| Attribute | Type | Description |
|---|---|---|
| `element` | `int` | Atomic number |
| `hybridization` | `str` | `SP3`, `SP2`, `SP`, `Unknown` |
| `formal_charge` | `int` | Formal charge |
| `explicit_h` | `int` | Explicitly attached hydrogens |

---

## Charges, surface and electrostatics

### `compute_charges`

```python
def compute_charges(smiles: str) -> ChargeResult
```

Gasteiger-Marsili partial charges, 6 iterations. Charges are in **elementary charge (e)**.

```python
result = sci_form.compute_charges("CCO")
print(result.charges)       # [-0.0418, 0.0414, -0.3953, 0.0252, ...]
print(result.total_charge)  # ~0.0
print(result.iterations, result.converged)
```

> **Deprecated alias:** `charges` → `compute_charges`.

### `compute_charges_configured`

```python
def compute_charges_configured(
    smiles: str,
    max_iter: int = 6,
    initial_damping: float = 0.5,
    convergence_threshold: float = 1e-10,
) -> ChargeResult
```

The same solver with an explicit convergence configuration: `max_iter` electronegativity
equalisation sweeps, `initial_damping` the first damping factor (halved every iteration) and
`convergence_threshold` the maximum charge-transfer convergence criterion.

```python
result = sci_form.compute_charges_configured("CCO", max_iter=12, initial_damping=0.7)
print(result.iterations)  # 12
```

### `compute_sasa`

```python
def compute_sasa(elements: list[int], coords: list[float], probe_radius: float = 1.4) -> SasaResult
```

Solvent-accessible surface area (Shrake-Rupley) in **Å²**.

```python
conf = sci_form.embed("CCO")
sasa = sci_form.compute_sasa(conf.elements, conf.coords, probe_radius=1.4)
print(sasa.total_sasa)   # 198.45192083308308  Å²
print(sasa.num_points)    # 960
```

> **Deprecated alias:** `sasa` → `compute_sasa`.

### `compute_population`

```python
def compute_population(elements: list[int], coords: list[float]) -> PopulationResult
```

Mulliken and Löwdin population analysis from the EHT density matrix. Charges in **e**, orbital
energies in **eV**.

```python
pop = sci_form.compute_population(conf.elements, conf.coords)
print(pop.mulliken_charges, pop.lowdin_charges)
print(pop.total_charge_mulliken, pop.charge_conservation_error)
print(pop.homo_energy, pop.lumo_energy, pop.gap)   # eV
```

> **Deprecated alias:** `population` → `compute_population`.

### `compute_dipole`

```python
def compute_dipole(elements: list[int], coords: list[float]) -> DipoleResult
```

Molecular dipole moment in **Debye**.

```python
d = sci_form.compute_dipole(conf.elements, conf.coords)
print(d.magnitude)   # 4.8790738513386565
print(d.vector)      # [x, y, z] in Debye
print(d.unit)        # "Debye"
```

> **Deprecated alias:** `dipole` → `compute_dipole`.

### `compute_esp_grid`

```python
def compute_esp_grid(elements: list[int], coords: list[float], spacing: float = 0.5, padding: float = 3.0) -> EspGrid
```

Electrostatic potential on a 3D grid, in **Hartree/e**. `spacing` and `padding` are in Å.

```python
esp = sci_form.compute_esp_grid(conf.elements, conf.coords, spacing=1.0, padding=2.0)
print(esp.dims, esp.spacing, esp.origin)
print(len(esp.values))
```

> **Deprecated aliases:** `compute_esp` and `esp` → `compute_esp_grid`.

### `compute_esp_grid_into` — zero copy

```python
def compute_esp_grid_into(
    elements: list[int],
    coords: list[float],
    out: "numpy.ndarray",
    spacing: float = 0.5,
    padding: float = 3.0,
) -> int
```

Compute the ESP grid and write the values **straight into a caller-owned buffer**. `out` must
be a writable, C-contiguous `float64` array of length `dims[0] * dims[1] * dims[2]`; a
`ValueError` is raised if the length does not match.

```python
import numpy as np, sci_form

info = sci_form.compute_esp([6, 6, 8], [0.0, 0.0, 0.0, 1.5, 0, 0, 3.0, 0, 0], 1.0, 1.0)
n = info.dims[0] * info.dims[1] * info.dims[2]
buf = np.zeros(n, dtype=np.float64)
sci_form.compute_esp_grid_into([6, 6, 8], [0.0, 0.0, 0.0, 1.5, 0, 0, 3.0, 0, 0], buf, 1.0, 1.0)
print(buf[:4])
```

### `compute_dos`

```python
def compute_dos(
    elements: list[int],
    coords: list[float],
    sigma: float = 0.3,
    e_min: float = -30.0,
    e_max: float = 5.0,
    n_points: int = 500,
) -> DosResult
```

Total DOS from the EHT orbital energies with Gaussian smearing. Energy axis in **eV**, DOS in
**1/eV**.

```python
dos = sci_form.compute_dos(conf.elements, conf.coords, sigma=0.3, e_min=-30.0, e_max=5.0)
print(dos.energies[0], dos.energies[-1], dos.sigma)
```

> **Deprecated alias:** `dos` → `compute_dos`.

### `compute_dos_multimethod`

```python
def compute_dos_multimethod(
    elements: list[int],
    coords: list[float],
    method: str = "eht",
    sigma: float = 0.3,
    e_min: float = -30.0,
    e_max: float = 5.0,
    n_points: int = 500,
) -> MultiMethodDosResult
```

The same DOS, but the orbital energies come from the chosen electronic-structure method.
`method` is one of `eht`, `pm3`, `xtb`, `gfn1`, `gfn2`, `hf3c`.

```python
dos = sci_form.compute_dos_multimethod(conf.elements, conf.coords, method="eht")
print(dos.method, dos.gap)   # 'Eht', 3.1085595007296085  eV
```

### `compute_rmsd`

```python
def compute_rmsd(coords: list[float], reference: list[float]) -> AlignmentResult
```

Kabsch optimal alignment and RMSD. RMSD in **Å**.

```python
a = sci_form.embed("CCO", seed=42)
b = sci_form.embed("CCO", seed=123)
result = sci_form.compute_rmsd(a.coords, b.coords)
print(result.rmsd)            # Å
print(result.aligned_coords)  # `coords` after the optimal rotation
print(result.rotation, result.translation)
```

> **Deprecated alias:** `rmsd` → `compute_rmsd`.

### `compute_rmsd_matrix`

```python
def compute_rmsd_matrix(conformers: list[list[float]]) -> list[list[float]]
```

All-pairs RMSD matrix (Å) with Kabsch alignment.

> **Deprecated alias:** `rmsd_matrix` → `compute_rmsd_matrix`.

### Return types

| Class | Attributes |
|---|---|
| `ChargeResult` | `charges` (e), `iterations`, `total_charge` (e), `converged` |
| `SasaResult` | `total_sasa` (Å²), `atom_sasa` (Å²), `probe_radius` (Å), `num_points` |
| `DipoleResult` | `vector` (Debye), `magnitude` (Debye), `unit`, `nuclear_dipole`, `electronic_dipole` |
| `DosResult` | `energies` (eV), `total_dos` (1/eV), `sigma` (eV) |
| `MultiMethodDosResult` | `energies` (eV), `total_dos` (1/eV), `method`, `homo_energy`, `lumo_energy`, `gap`, `orbital_energies`, `sigma` |
| `AlignmentResult` | `rmsd` (Å), `aligned_coords` (Å), `rotation`, `translation` |
| `EspGrid` | `dims`, `spacing` (Å), `origin` (Å), `values` (Hartree/e) |
| `PopulationResult` | `mulliken_charges` (e), `lowdin_charges` (e), `mulliken_populations`, `num_atoms`, `total_charge_mulliken` (e), `total_charge_lowdin` (e), `charge_conservation_error`, `homo_energy` (eV), `lumo_energy` (eV), `gap` (eV), `warnings` |

---

## Electronic structure

### `solve_eht`

```python
def solve_eht(elements: list[int], coords: list[float], k: float = 1.75) -> EhtResult
```

Extended-Hückel calculation. `k` is the Wolfsberg-Helmholtz constant. Orbital energies and the
HOMO-LUMO gap are in **eV**.

```python
eht = sci_form.solve_eht(conf.elements, conf.coords)
print(eht.homo_energy, eht.lumo_energy, eht.gap)   # eV
print(eht.energies)                                  # all MO energies, eV
print(eht.homo_index, eht.lumo_index, eht.n_electrons)
```

> **Deprecated alias:** `eht_calculate` → `solve_eht`.

### `eht_support`

```python
def eht_support(elements: list[int]) -> EhtSupport
```

Report whether an element set is inside the EHT domain and with which confidence.

```python
support = sci_form.eht_support(conf.elements)
print(support.level, support.has_transition_metals)
print(support.supported_elements, support.provisional_elements, support.unsupported_elements)
```

### `compute_eht_or_uff_fallback`

```python
def compute_eht_or_uff_fallback(
    smiles: str,
    elements: list[int],
    coords: list[float],
    allow_experimental_eht: bool = False,
) -> dict
```

Run EHT when the element set is supported, otherwise fall back to UFF. The returned mapping has
a `mode` key (`"eht"` or `"uff"`) plus the method payload; for the EHT path the payload carries
the orbital energies in **eV**, `homo_energy`, `lumo_energy`, `gap`, `n_electrons`,
`homo_index`, `lumo_index`, `support_level` and `warnings`.

```python
out = sci_form.compute_eht_or_uff_fallback("CCO", conf.elements, conf.coords)
print(out["mode"], out.get("gap"))
```

> **Deprecated alias:** `eht_or_uff_fallback` → `compute_eht_or_uff_fallback`.

### `compute_eht_gradient`

```python
def compute_eht_gradient(elements: list[int], coords: list[float]) -> EhtGradient
```

Central-finite-difference gradients on the EHT energy surface, in **eV/Å**.

```python
g = sci_form.compute_eht_gradient(conf.elements, conf.coords)
print(g.rms_gradient, g.max_gradient)   # eV/Å
print(g.energy)                          # eV
```

### `optimize_geometry_eht`

```python
def optimize_geometry_eht(
    elements: list[int],
    coords: list[float],
    max_steps: int = 200,
    grad_threshold: float = 0.01,
    energy_threshold: float = 1e-6,
    step_size: float = 0.05,
) -> EhtOptResult
```

Optimise a geometry on the EHT surface with the shared BFGS engine. Energies in **eV**,
gradients in **eV/Å**.

```python
opt = sci_form.optimize_geometry_eht(conf.elements, conf.coords)
print(opt.converged, opt.n_steps, opt.energy, opt.rms_gradient)
```

### `eht_orbital_mesh`

```python
def eht_orbital_mesh(
    elements: list[int],
    coords: list[float],
    mo_index: int,
    spacing: float = 0.2,
    isovalue: float = 0.02,
) -> OrbitalMeshResult
```

Isosurface mesh for one EHT molecular orbital. Vertices and normals in **Å**.

### `compute_pm3`

```python
def compute_pm3(elements: list[int], coords: list[float]) -> Pm3Result
```

PM3 NDDO semi-empirical SCF. Supports transition metals via PM3(tm). **Total, electronic and
nuclear-repulsion energies are in eV**; `heat_of_formation` is in **kcal/mol**; orbital energies
and the gap are in **eV**.

```python
pm3 = sci_form.compute_pm3(conf.elements, conf.coords)
print(pm3.total_energy, pm3.electronic_energy, pm3.nuclear_repulsion)   # eV
print(pm3.heat_of_formation)                                           # kcal/mol
print(pm3.homo_energy, pm3.lumo_energy, pm3.gap)                       # eV
print(pm3.converged, pm3.scf_iterations)
```

> **Deprecated alias:** `pm3_calculate` → `compute_pm3`.

### `compute_xtb`

```python
def compute_xtb(elements: list[int], coords: list[float]) -> XtbResult
```

GFN0-xTB tight binding with self-consistent charges (25 elements, including transition
metals). All energies in **eV**, orbital energies and gap in **eV**.

```python
xtb = sci_form.compute_xtb(conf.elements, conf.coords)
print(xtb.total_energy, xtb.electronic_energy, xtb.repulsion_energy)   # eV
print(xtb.homo_energy, xtb.lumo_energy, xtb.gap)                       # eV
```

> **Deprecated alias:** `xtb_calculate` → `compute_xtb`.

### `compute_gfn1`

```python
def compute_gfn1(elements: list[int], coords: list[float]) -> Gfn1Result
```

GFN1-xTB with shell-resolved charges and D3-style dispersion. All energies in **eV**.

```python
gfn1 = sci_form.compute_gfn1(conf.elements, conf.coords)
print(gfn1.total_energy, gfn1.dispersion_energy)   # eV
print(gfn1.shell_charges)                            # per-atom s/p/d shell charges, e
```

> **Deprecated alias:** `gfn1_calculate` → `compute_gfn1`.

### `compute_gfn2`

```python
def compute_gfn2(elements: list[int], coords: list[float]) -> Gfn2Result
```

GFN2-xTB with multipole electrostatics, D4-style dispersion and halogen-bond corrections. All
energies in **eV**; `atomic_dipoles` in **Debye**.

```python
gfn2 = sci_form.compute_gfn2(conf.elements, conf.coords)
print(gfn2.total_energy, gfn2.dispersion_energy, gfn2.halogen_bond_energy)   # eV
print(gfn2.atomic_dipoles)                                                   # Debye
```

> **Deprecated alias:** `gfn2_calculate` → `compute_gfn2`.

### `compute_hf3c`

```python
def compute_hf3c(
    elements: list[int],
    coords: list[float],
    max_scf_iter: int = 300,
    n_cis_states: int = 5,
    corrections: bool = True,
) -> Hf3cResult
```

Minimal-basis Hartree-Fock with D3, gCP and SRB corrections. **All energies and orbital energies
are in Hartree.**

```python
hf = sci_form.compute_hf3c(conf.elements, conf.coords)
print(hf.energy, hf.hf_energy, hf.nuclear_repulsion)   # Hartree
print(hf.d3_energy, hf.gcp_energy, hf.srb_energy)      # Hartree
print(hf.homo_energy, hf.lumo_energy, hf.gap)          # Hartree
```

> **Deprecated alias:** `hf3c_calculate` → `compute_hf3c`.

### `compute_ani`

```python
def compute_ani(elements: list[int], coords: list[float]) -> AniResult
```

ANI neural-network potential for H, C, N, O, F, S, Cl. `energy` and `atomic_energies` are in
**Hartree**, `forces` in **Hartree/Å**. The bundled weights are *test* weights, so the energies
are not physically meaningful ANI-2x energies.

> **Deprecated alias:** `ani_calculate` → `compute_ani`.

### `compute_cisd`

```python
def compute_cisd(elements: list[int], coords: list[float], n_states: int = 5) -> CisdResult
```

CISD excited states from an internal STO-3G RHF calculation. `correlation_energy` and each
excitation's `energy` are in **Hartree**; `energy_ev` is the same excitation in **eV**.

```python
cisd = sci_form.compute_cisd(conf.elements, conf.coords, n_states=2)
print(cisd.approximation, cisd.n_csfs, cisd.correlation_energy)
for ex in cisd.excitations:
    print(ex.energy_ev, ex.wavelength_nm, ex.oscillator_strength)
```

### `compute_uhf`

```python
def compute_uhf(elements: list[int], coords: list[float], charge: int = 0, multiplicity: int = 1) -> UhfResult
```

Unrestricted Hartree-Fock for open-shell systems. The energies are in **Hartree**, the
HOMO-LUMO gap in **eV** and the Mulliken charges in **e**. `charge` is the total
molecular charge and `multiplicity` the spin multiplicity (1 = singlet, 2 = doublet,
3 = triplet); a multiplicity that does not match the electron count is an error.

```python
uhf = sci_form.compute_uhf([6, 1, 1, 1], coords, 0, 2)   # a methyl radical
print(uhf.converged, uhf.total_energy, uhf.s2_expectation)
```

### `compute_rohf`

```python
def compute_rohf(elements: list[int], coords: list[float], charge: int = 0, multiplicity: int = 1) -> UhfResult
```

Restricted open-shell Hartree-Fock: the same spatial orbitals for both spins, so the
spin contamination is zero by construction.

### `compute_pm3_gradient`

```python
def compute_pm3_gradient(elements: list[int], coords: list[float]) -> Pm3GradientResult
```

Analytical PM3 energy gradient. The gradient is in **eV/Å** (one `[x, y, z]` per
atom), the energy in **eV** and the heat of formation in **kcal/mol**.

### `compute_xtb_gradient`

```python
def compute_xtb_gradient(elements: list[int], coords: list[float]) -> XtbGradientResult
```

Analytical xTB energy gradient, in **eV/Å** per atom; the energy is in **eV**.

### `compute_aevs_tm`

```python
def compute_aevs_tm(elements: list[int], coords: list[float]) -> AniTmAevs
```

ANI-TM atomic environment vectors over the extended element set (24 elements). Dimensionless.

### Return types

| Class | Attributes (unit) |
|---|---|
| `EhtResult` | `energies` (eV), `n_electrons`, `homo_index`, `lumo_index`, `homo_energy` (eV), `lumo_energy` (eV), `gap` (eV), `support_level`, `has_transition_metals`, `supported_elements`, `provisional_elements`, `unsupported_elements`, `warnings` |
| `UhfResult` | `alpha_orbital_energies`, `beta_orbital_energies` (Hartree), `electronic_energy`, `nuclear_repulsion`, `total_energy` (Hartree), `homo_energy` (Hartree), `lumo_energy` (Hartree), `gap_ev` (eV), `mulliken_charges` (e), `scf_iterations`, `converged`, `n_basis`, `n_alpha`, `n_beta`, `s2_expectation`, `spin_contamination`, `alpha_coefficients`, `beta_coefficients` |
| `Pm3GradientResult` | `gradients` (eV/Å), `energy` (eV), `heat_of_formation` (kcal/mol) |
| `XtbGradientResult` | `gradients` (eV/Å), `energy` (eV) |
| `EhtSupport` | `level`, `has_transition_metals`, `supported_elements`, `provisional_elements`, `unsupported_elements`, `warnings` |
| `EhtGradient` | `gradients` (eV/Å), `rms_gradient` (eV/Å), `max_gradient` (eV/Å), `energy` (eV) |
| `EhtOptResult` | `positions` (Å), `energy` (eV), `n_steps`, `converged`, `rms_gradient` (eV/Å), `energies` (eV) |
| `Pm3Result` | `orbital_energies` (eV), `electronic_energy` (eV), `nuclear_repulsion` (eV), `total_energy` (eV), `heat_of_formation` (kcal/mol), `n_basis`, `n_electrons`, `homo_energy` (eV), `lumo_energy` (eV), `gap` (eV), `mulliken_charges` (e), `scf_iterations`, `converged` |
| `XtbResult` | `orbital_energies` (eV), `electronic_energy` (eV), `repulsion_energy` (eV), `total_energy` (eV), `n_basis`, `n_electrons`, `homo_energy` (eV), `lumo_energy` (eV), `gap` (eV), `mulliken_charges` (e), `scf_iterations`, `converged` |
| `Gfn1Result` | as `XtbResult` plus `dispersion_energy` (eV), `shell_charges` (e) |
| `Gfn2Result` | as `XtbResult` plus `dispersion_energy` (eV), `dispersion_2_body_energy` (eV), `halogen_bond_energy` (eV), `atomic_dipoles` (Debye), `atomic_quadrupoles` |
| `Hf3cResult` | `energy`, `hf_energy`, `nuclear_repulsion`, `d3_energy`, `gcp_energy`, `srb_energy`, `orbital_energies` (Hartree), `scf_iterations`, `converged`, `n_basis`, `n_electrons`, `homo_energy`, `lumo_energy`, `gap` (Hartree), `mulliken_charges` (e) |
| `AniResult` | `energy` (Hartree), `forces` (Hartree/Å), `species`, `atomic_energies` (Hartree), `aevs` |
| `CisdResult` | `excitations`, `correlation_energy` (Hartree), `n_csfs`, `approximation` |
| `CisdExcitation` | `energy` (Hartree), `energy_ev` (eV), `wavelength_nm`, `oscillator_strength`, `dominant_transition`, `character`, `singles_weight` |
| `AniTmAevs` | `aev_length`, `aevs`, `n_atoms` |

---

## Force fields

### `compute_uff_energy`

```python
def compute_uff_energy(smiles: str, coords: list[float]) -> float
```

UFF force-field energy in **kcal/mol**.

```python
print(sci_form.compute_uff_energy("CCO", conf.coords))   # 6.993048...
```

> **Deprecated alias:** `uff_energy` → `compute_uff_energy`.

### `compute_uff_energy_with_aromatic_heuristics`

```python
def compute_uff_energy_with_aromatic_heuristics(smiles: str, coords: list[float]) -> UffHeuristicEnergy
```

UFF energy with an aromaticity-informed correction. All energies in **kcal/mol**.

> **Deprecated alias:** `uff_energy_with_aromatic_heuristics` →
> `compute_uff_energy_with_aromatic_heuristics`.

### `compute_mmff94_energy`

```python
def compute_mmff94_energy(smiles: str, coords: list[float]) -> float
```

MMFF94 force-field energy in **kcal/mol**.

> **Deprecated alias:** `mmff94_energy` → `compute_mmff94_energy`.

### `compute_uff_energy_breakdown`

```python
def compute_uff_energy_breakdown(smiles: str, coords: list[float]) -> list[tuple[str, float]]
```

The UFF energy split by term type (`bond`, `angle`, `torsion`, `inversion`, `vdw`).
The values are in **kcal/mol** and sum to `compute_uff_energy`.

### `compute_mmff94_energy_breakdown`

```python
def compute_mmff94_energy_breakdown(smiles: str, coords: list[float]) -> list[tuple[str, float]]
```

The same split for MMFF94; the values sum to `compute_mmff94_energy`.

### `compute_empirical_pka`

```python
def compute_empirical_pka(smiles: str) -> EmpiricalPkaResult
```

Graph-heuristic pKa estimation. The pKa values are coarse screening hints, not
publication-grade predictions.

> **Deprecated alias:** `empirical_pka` → `compute_empirical_pka`.

### Return types

| Class | Attributes |
|---|---|
| `UffHeuristicEnergy` | `raw_energy_kcal_mol`, `aromatic_stabilization_kcal_mol`, `corrected_energy_kcal_mol` (kcal/mol), `aromatic_bond_count`, `notes` |
| `EmpiricalPkaResult` | `acidic_sites`, `basic_sites`, `notes` |
| `EmpiricalPkaSite` | `atom_index`, `estimated_pka`, `site_type`, `environment`, `confidence` |

---

## Molecular descriptors and machine learning

### `compute_ml_descriptors`

```python
def compute_ml_descriptors(smiles: str) -> MolecularDescriptors
```

Topology-derived descriptors. No 3D geometry needed.

> **Deprecated alias:** `ml_descriptors` → `compute_ml_descriptors`.

### `predict_ml_properties`

```python
def predict_ml_properties(smiles: str) -> MlPropertyResult
```

Lipinski / ADME proxy prediction.

> **Deprecated alias:** `ml_predict` → `predict_ml_properties`.

### `compute_tpsa`

```python
def compute_tpsa(smiles: str) -> float
```

Topological polar surface area in **Å²** (Ertl fragment contributions). It is a
topological property: no 3D geometry is needed.

### `predict_ensemble`

```python
def predict_ensemble(smiles: str) -> EnsembleResult
```

Ensemble property prediction: the consensus LogP with its uncertainty, the
solubility, the TPSA, the pKa estimates, the Veber rules, the BBB permeability
and an overall confidence score (0–1).

### `compute_whim`

```python
def compute_whim(elements: list[int], coords: list[float], weighting: str = "mass") -> WhimDescriptors
```

WHIM 3D molecular descriptors. `weighting` is one of `unit`, `mass`, `volume`,
`electronegativity`, `polarizability`.

### `compute_rdf`

```python
def compute_rdf(
    elements: list[int],
    coords: list[float],
    charges: list[float] | None = None,
    beta: float = 100.0,
    r_max: float = 12.0,
    n_points: int = 128,
) -> RdfDescriptors
```

Radial distribution function descriptors. `radii` are in **Å**, `beta` in **Å⁻²**. When
`charges` is `None` the charges are derived internally with Gasteiger.

### `compute_getaway`

```python
def compute_getaway(
    elements: list[int],
    coords: list[float],
    bonds: list[tuple[int, int]] | None = None,
    max_lag: int = 8,
) -> GetawayDescriptors
```

GETAWAY 3D molecular descriptors. When `bonds` is `None` the bonds are perceived from the
geometry.

### `train_random_forest`

```python
def train_random_forest(
    features: list[list[float]],
    targets: list[float],
    n_trees: int = 100,
    max_depth: int = 10,
    min_samples_split: int = 5,
    min_samples_leaf: int = 2,
    sample_fraction: float = 1.0,
    seed: int = 42,
    options: dict | None = None,
) -> RegressionModel
```

Train a Random Forest regressor and return its training-set predictions.

### `train_gradient_boosting`

```python
def train_gradient_boosting(
    features: list[list[float]],
    targets: list[float],
    n_estimators: int = 100,
    learning_rate: float = 0.1,
    max_depth: int = 5,
    min_samples_split: int = 5,
    min_samples_leaf: int = 2,
    subsample: float = 0.8,
    seed: int = 42,
    options: dict | None = None,
) -> RegressionModel
```

Train a Gradient Boosting regressor and return its training-set predictions.

### Return types

| Class | Attributes |
|---|---|
| `MolecularDescriptors` | `molecular_weight`, `n_heavy_atoms`, `n_hydrogens`, `n_bonds`, `n_rotatable_bonds`, `n_hbd`, `n_hba`, `fsp3`, `wiener_index`, `n_rings`, `n_aromatic`, `sum_electronegativity`, `sum_polarizability` |
| `MlPropertyResult` | `logp`, `molar_refractivity`, `log_solubility`, `druglikeness`, `lipinski_violations`, `lipinski_passes` |
| `WhimDescriptors` | `eigenvalues`, `lambda_`, `nu`, `gamma`, `eta`, `total_spread`, `anisotropy`, `directional`, `kurtosis` |
| `RdfDescriptors` | `radii` (Å), `rdf_unweighted`, `rdf_mass`, `rdf_electronegativity`, `rdf_charge`, `beta` (Å⁻²), `n_pairs` |
| `GetawayDescriptors` | `leverages`, `hat_autocorrelation`, `r_autocorrelation`, `hat_total`, `r_total`, `h_information`, `h_max`, `h_mean`, `max_lag` |
| `RegressionModel` | `predictions`, `train_rmse`, `model_type`, `n_estimators`, `n_samples` |
| `EnsembleResult` | `logp`, `logp_std`, `log_solubility`, `tpsa` (Å²), `pka_acidic`, `pka_basic`, `veber`, `bbb_permeable`, `bbb_score` (0–1), `confidence` (0–1) |
| `VeberResult` | `tpsa_ok`, `rotb_ok`, `passes` |

---

## Spectroscopy

### `compute_stda_uvvis`

```python
def compute_stda_uvvis(
    elements: list[int],
    coords: list[float],
    sigma: float = 0.3,
    e_min: float = 1.0,
    e_max: float = 8.0,
    n_points: int = 500,
    broadening: str = "gaussian",
) -> StdaUvVisSpectrum
```

sTDA UV-Vis spectrum on EHT MO transitions. `energies_ev` in **eV**, `wavelengths_nm` in **nm**,
`absorptivity` in **L mol⁻¹ cm⁻¹**. `broadening` is `gaussian` or `lorentzian`.

> **Deprecated alias:** `stda_uvvis` → `compute_stda_uvvis`.

### `compute_fluorescence`

```python
def compute_fluorescence(
    elements: list[int],
    coords: list[float],
    sigma: float = 0.25,
    e_min: float = 0.5,
    e_max: float = 8.0,
    n_points: int = 400,
    stokes_ev: float = 0.0,
    refractive_index: float = 1.0,
    dark_threshold: float = 1e-3,
    broadening: str = "gaussian",
    options: dict | None = None,
) -> FluorescenceSpectrum
```

Fluorescence from the same sTDA singlets as the absorption spectrum. `radiative_rate_hz` in
**Hz**, `radiative_lifetime_ns` in **ns**, energies in **eV**, wavelengths in **nm**.

> **Deprecated alias:** `fluorescence_spectrum` → `compute_fluorescence`.

### `compute_vibrational_analysis`

```python
def compute_vibrational_analysis(
    elements: list[int],
    coords: list[float],
    method: str = "eht",
    step_size: float | None = None,
    smiles: str | None = None,
) -> VibrationalAnalysis
```

Numerical Hessian (`6N` energy evaluations) plus diagonalisation; translation/rotation modes are
projected out. Frequencies in **cm⁻¹**, IR intensities in **km/mol**, `zpve_ev` in **eV**.

```python
vib = sci_form.compute_vibrational_analysis(conf.elements, conf.coords, method="xtb")
print(vib.n_real_modes, vib.zpve_ev)
for mode in vib.modes:
    print(mode.frequency_cm1, mode.ir_intensity, mode.is_real)
```

::: warning
The EHT Hessian is a one-electron model with no electron-electron repulsion: the surface has no
minimum, so imaginary (negative) frequencies are expected. Use `pm3` or `xtb` for vibrational
analysis.
:::

> **Deprecated alias:** `vibrational_analysis` → `compute_vibrational_analysis`.

### `compute_ir_spectrum`

```python
def compute_ir_spectrum(
    elements: list[int],
    coords: list[float],
    method: str = "eht",
    gamma: float = 15.0,
    wn_min: float = 400.0,
    wn_max: float = 4000.0,
    n_points: int = 1000,
    step_size: float | None = None,
    smiles: str | None = None,
    options: dict | None = None,
) -> IrSpectrum
```

Runs the vibrational analysis internally and applies the default **Lorentzian** broadening.
`wavenumbers` in **cm⁻¹**, `gamma` in **cm⁻¹**, `intensities` in **km/mol**.

> **Deprecated alias:** `ir_spectrum` → `compute_ir_spectrum`.

### `compute_ir_spectrum_broadened`

```python
def compute_ir_spectrum_broadened(
    elements: list[int],
    coords: list[float],
    method: str = "eht",
    gamma: float = 15.0,
    wn_min: float = 400.0,
    wn_max: float = 4000.0,
    n_points: int = 1000,
    broadening: str = "lorentzian",
    step_size: float | None = None,
    smiles: str | None = None,
    options: dict | None = None,
) -> IrSpectrum
```

An unknown `broadening` raises `ValueError`: the value used to fall back to
the Lorentzian shape, so the typo `gauassian` silently changed the spectrum.

Same as `compute_ir_spectrum` but with a selectable broadening type (`gaussian` or
`lorentzian`). `step_size` is the finite-difference step in **Å** (`None` uses the engine
default) and `smiles` is required when `method="uff"`.

### `assign_peaks`

```python
def assign_peaks(
    frequencies: list[float],
    intensities: list[float],
    elements: list[int],
    mode_displacements: list[list[list[float]]] | None = None,
) -> AssignmentResult
```

Assign IR peaks to functional groups using the standard group-frequency tables. `frequencies`
are in **cm⁻¹**.

```python
assignment = sci_form.assign_peaks([1000.0, 3000.0], [10.0, 20.0], conf.elements)
for peak in assignment.assignments:
    print(peak.frequency, peak.confidence)
    for match in peak.assignments:
        print("  ", match.group, match.vibration_type, match.match_quality)
print(assignment.functional_groups)
```

### Return types

| Class | Attributes (unit) |
|---|---|
| `StdaUvVisSpectrum` | `energies_ev` (eV), `wavelengths_nm` (nm), `absorptivity`, `excitations`, `sigma` (eV), `broadening`, `notes` |
| `StdaExcitation` | `energy_ev` (eV), `wavelength_nm` (nm), `oscillator_strength`, `from_mo`, `to_mo`, `transition_dipole` |
| `FluorescenceSpectrum` | `energies_ev` (eV), `wavelengths_nm` (nm), `intensity`, `absorption_energy_ev` (eV), `emission_energy_ev` (eV), `emission_wavelength_nm` (nm), `oscillator_strength`, `radiative_rate_hz` (Hz), `radiative_lifetime_ns` (ns), `kasha`, `from_mo`, `to_mo`, `notes` |
| `VibrationalAnalysis` | `n_atoms`, `elements`, `modes`, `n_real_modes`, `n_imaginary_modes`, `zpve_ev` (eV), `thermochemistry`, `method`, `notes` |
| `VibrationalMode` | `frequency_cm1` (cm⁻¹), `ir_intensity` (km/mol), `displacement`, `is_real` |
| `IrSpectrum` | `wavenumbers` (cm⁻¹), `intensities` (km/mol), `peaks`, `gamma` (cm⁻¹), `notes` |
| `IrPeak` | `frequency_cm1` (cm⁻¹), `ir_intensity` (km/mol), `mode_index` |
| `PeakAssignment` | `frequency` (cm⁻¹), `intensity`, `assignments`, `confidence` |
| `AssignmentResult` | `assignments`, `functional_groups` |
| `FunctionalGroupMatch` | `group`, `vibration_type`, `match_quality` |

---

## NMR

### `predict_nmr_shifts`

```python
def predict_nmr_shifts(smiles: str) -> NmrShiftResult
```

Representative chemical shifts for the default nucleus of each element. Shifts in **ppm**, no 3D
geometry needed. The legacy per-element fields (`h_shifts`, `c_shifts`, `f_shifts`, `p_shifts`,
`n_shifts`, `b_shifts`, `si_shifts`, `se_shifts`, `o_shifts`, `s_shifts`) plus `other_shifts`
cover the whole registry.

> **Deprecated alias:** `nmr_shifts` → `predict_nmr_shifts`.

### `predict_nmr_shifts_for_nucleus`

```python
def predict_nmr_shifts_for_nucleus(smiles: str, nucleus: str) -> NucleusShiftSeries
```

Shifts for one requested nucleus, including the expanded registry (`2H`, `35Cl`, `79Br`,
`195Pt`, `207Pb`, …).

> **Deprecated alias:** `nmr_shifts_for_nucleus` → `predict_nmr_shifts_for_nucleus`.

### `predict_nmr_couplings`

```python
def predict_nmr_couplings(smiles: str, coords: list[float] = []) -> list[JCoupling]
```

J-coupling constants in **Hz**. With 3D coordinates the Karplus relation is used for ³J;
topology only gives free-rotation averages. Parameterised pathways: H-C-C-H
(Altona-Sundaralingam), H-C-N-H (Bystrov), H-C-O-H, H-C-S-H.

```python
cc = sci_form.embed("CC")
for j in sci_form.predict_nmr_couplings("CC", cc.coords):
    print(j.n_bonds, j.h1_index, j.h2_index, j.j_hz, j.coupling_type)
```

> **Deprecated alias:** `nmr_couplings` → `predict_nmr_couplings`.

### `compute_ensemble_j_couplings`

```python
def compute_ensemble_j_couplings(
    smiles: str,
    conformer_coords: list[list[float]],
    energies_kcal: list[float],
    temperature_k: float = 298.15,
) -> list[JCoupling]
```

J-couplings in **Hz**, averaged over a conformer ensemble with Boltzmann weighting.
`energies_kcal` holds the matching conformer energies in **kcal/mol**.

### `compute_nmr_spectrum`

```python
def compute_nmr_spectrum(
    smiles: str,
    coords: list[float] = [],
    nucleus: str = "1H",
    gamma: float = 0.02,
    ppm_min: float = 0.0,
    ppm_max: float = 12.0,
    n_points: int = 1000,
) -> NmrSpectrum
```

Full NMR spectrum pipeline. `ppm_axis` in **ppm**, `gamma` in **ppm**.

> **Deprecated alias:** `nmr_spectrum` → `compute_nmr_spectrum`.

### `compute_giao_nmr`

```python
def compute_giao_nmr(
    elements: list[int],
    coords: list[float],
    nucleus: str = "1H",
    charge: int = 0,
    multiplicity: int = 1,
    max_scf_iter: int = 100,
    allow_basis_fallback: bool = False,
    use_parallel_eri: bool = False,
    options: dict | None = None,
) -> GiaoNmrResult
```

The public SCF-backed GIAO route for one nucleus. `total_energy_hartree` in **Hartree**,
`chemical_shifts` in **ppm**, `shieldings` in **ppm**. The current public path is restricted to
singlet closed-shell systems and by default rejects elements that only have the hydrogen-like
fallback basis in the SCF stack.

> **Deprecated alias:** `giao_nmr` → `compute_giao_nmr`.

### `compute_hose_codes`

```python
def compute_hose_codes(smiles: str, max_radius: int = 2) -> list[tuple[int, int, str]]
```

HOSE code fingerprints for all atoms: one `(atom_index, element, full_code)` tuple per atom.

> **Deprecated alias:** `hose_codes` → `compute_hose_codes`.

### Return types

| Class | Attributes (unit) |
|---|---|
| `NmrShiftResult` | `h_shifts`, `c_shifts`, `f_shifts`, `p_shifts`, `n_shifts`, `b_shifts`, `si_shifts`, `se_shifts`, `o_shifts`, `s_shifts`, `other_shifts`, `notes` |
| `ChemicalShift` | `atom_index`, `element`, `shift_ppm` (ppm), `environment`, `confidence` |
| `NucleusShiftSeries` | `nucleus`, `shifts` |
| `JCoupling` | `h1_index`, `h2_index`, `j_hz` (Hz), `n_bonds`, `coupling_type` |
| `NmrPeak` | `atom_index`, `environment`, `intensity`, `multiplicity`, `shift_ppm` (ppm) |
| `NmrSpectrum` | `ppm_axis` (ppm), `intensities`, `peaks`, `nucleus`, `gamma` (ppm), `notes` |
| `GiaoNmrResult` | `nucleus`, `target_atomic_number`, `method`, `basis_set`, `charge`, `multiplicity`, `scf_converged`, `scf_iterations`, `total_energy_hartree` (Hartree), `n_basis`, `n_target_atoms`, `chemical_shifts` (ppm), `shieldings`, `fallback_elements`, `notes` |
| `GiaoNmrEntry` | `atom_index`, `element`, `isotropic` (ppm), `anisotropy` (ppm), `tensor`, `chemical_shift` (ppm) |

---

## Solvation

### `compute_nonpolar_solvation`

```python
def compute_nonpolar_solvation(elements: list[int], coords: list[float], probe_radius: float = 1.4) -> NonPolarSolvation
```

Non-polar solvation `ΔG_np = Σ_i σ_i A_i` from per-element ASP and Shrake-Rupley SASA.
Energy in **kcal/mol**, areas in **Å²**.

> **Deprecated alias:** `nonpolar_solvation` → `compute_nonpolar_solvation`.

### `compute_gb_solvation`

```python
def compute_gb_solvation(
    elements: list[int],
    coords: list[float],
    charges: list[float],
    solvent_dielectric: float = 78.5,
    solute_dielectric: float = 1.0,
    probe_radius: float = 1.4,
) -> GbSolvation
```

GB/SA electrostatic solvation (HCT Born radii + Still GB equation). Energies in **kcal/mol**,
born radii in **Å**.

**Raises `ValueError`** when the input is malformed: `coords` does not have
`3 * len(elements)` values, `charges` does not have one value per atom, a
charge is `NaN`/`±Inf`, a dielectric constant is not finite or `< 1.0`
(`ε = 1.0` is the vacuum limit and is accepted), or the probe radius is
negative or not finite.

> **Deprecated alias:** `gb_solvation` → `compute_gb_solvation`.

### Return types

| Class | Attributes (unit) |
|---|---|
| `NonPolarSolvation` | `energy_kcal_mol` (kcal/mol), `atom_contributions` (kcal/mol), `atom_sasa` (Å²), `total_sasa` (Å²) |
| `GbSolvation` | `electrostatic_energy_kcal_mol`, `nonpolar_energy_kcal_mol`, `total_energy_kcal_mol` (kcal/mol), `born_radii` (Å), `charges` (e), `solvent_dielectric`, `solute_dielectric` |

---

## Rings, fingerprints and clustering

### `compute_sssr`

```python
def compute_sssr(smiles: str) -> SssrResult
```

Smallest Set of Smallest Rings (Horton) with Hückel aromaticity.

```python
r = sci_form.compute_sssr("c1ccc2ccccc2c1")
print(len(r.rings), r.ring_size_histogram)
```

> **Deprecated alias:** `sssr` → `compute_sssr`.

### `compute_ecfp`

```python
def compute_ecfp(smiles: str, radius: int = 2, n_bits: int = 2048) -> ECFPFingerprint
```

Extended-Connectivity Fingerprint (Morgan). `radius=2` is ECFP4, `radius=3` is ECFP6.

> **Deprecated alias:** `ecfp` → `compute_ecfp`.

### `compute_ecfp_batch`

```python
def compute_ecfp_batch(smiles_list: list[str], radius: int = 2, n_bits: int = 2048) -> list[ECFPFingerprint]
```

Batch ECFP generation, in input order.

### `compute_tanimoto`

```python
def compute_tanimoto(fp1: ECFPFingerprint, fp2: ECFPFingerprint) -> float
```

Jaccard-Tanimoto similarity `T = |A ∩ B| / |A ∪ B|`, dimensionless in `[0, 1]`.

> **Deprecated alias:** `tanimoto` → `compute_tanimoto`.

### `butina_cluster`

```python
def butina_cluster(conformers: list[list[float]], rmsd_cutoff: float = 1.0) -> ClusterResult
```

Taylor-Butina single-linkage RMSD clustering. Cutoff in **Å**.

### `filter_diverse`

```python
def filter_diverse(conformers: list[list[float]], rmsd_cutoff: float = 1.0) -> list[int]
```

Indices of the cluster centroid (diverse representative) conformers.

### Return types

| Class | Attributes |
|---|---|
| `SssrResult` | `rings`, `atom_ring_count`, `atom_ring_sizes`, `ring_size_histogram` |
| `RingInfo` | `atoms`, `size`, `is_aromatic` |
| `ECFPFingerprint` | `n_bits`, `on_bits`, `radius`, `density` |
| `ClusterResult` | `n_clusters`, `assignments`, `centroid_indices`, `cluster_sizes`, `rmsd_cutoff` (Å) |

---

## Stereochemistry, topology and reactivity

### `analyze_stereo`

```python
def analyze_stereo(smiles: str, coords: list[float] = []) -> StereoAnalysis
```

Assign CIP priorities and determine R/S at stereocenters and E/Z at double bonds. With an
empty `coords` the analysis is topology-only and `configuration` stays `None`.

> **Deprecated alias:** `stereo_analysis` → `analyze_stereo`.

### `analyze_graph_features`

```python
def analyze_graph_features(smiles: str) -> GraphFeatureAnalysis
```

Aromaticity and graph-level stereocentre analysis.

> **Deprecated alias:** `graph_features` → `analyze_graph_features`.

### `compute_topology`

```python
def compute_topology(elements: list[int], coords: list[float]) -> TopologyAnalysisResult
```

Metal coordination environments from 3D coordinates: ligand spheres, coordination numbers and a
geometric guess of the coordination polyhedron.

> **Deprecated alias:** `topology_analysis` → `compute_topology`.

### `compute_uv_vis_spectrum`

```python
def compute_uv_vis_spectrum(
    elements: list[int],
    coords: list[float],
    sigma: float = 0.25,
    e_min: float = 0.5,
    e_max: float = 8.0,
    n_points: int = 600,
) -> UvVisSpectrum
```

Exploratory UV-Vis-like spectrum from EHT MO energy differences. `energies_ev` in **eV**,
`wavelengths_nm` in **nm**.

> **Deprecated alias:** `uvvis_spectrum` → `compute_uv_vis_spectrum`.

### `compute_frontier_descriptors`

```python
def compute_frontier_descriptors(elements: list[int], coords: list[float]) -> FrontierDescriptors
```

HOMO/LUMO atom contributions and the dual descriptor. Energies in **eV**.

> **Deprecated alias:** `frontier_descriptors` → `compute_frontier_descriptors`.

### `compute_fukui_descriptors`

```python
def compute_fukui_descriptors(elements: list[int], coords: list[float]) -> FukuiDescriptors
```

Condensed Fukui functions f⁺, f⁻, f⁰. Energies in **eV**.

> **Deprecated alias:** `fukui_descriptors` → `compute_fukui_descriptors`.

### `compute_reactivity_ranking`

```python
def compute_reactivity_ranking(elements: list[int], coords: list[float]) -> ReactivityRanking
```

Ranked atomic sites for nucleophilic, electrophilic and radical attack. The scores are empirical
composites, not calibrated kinetic barriers.

> **Deprecated alias:** `reactivity_ranking` → `compute_reactivity_ranking`.

### `compute_bond_orders`

```python
def compute_bond_orders(elements: list[int], coords: list[float]) -> BondOrderResult
```

Wiberg and Mayer bond orders from the EHT density matrix. Distances in **Å**.

> **Deprecated alias:** `bond_orders` → `compute_bond_orders`.

### Return types

| Class | Attributes (unit) |
|---|---|
| `StereoAnalysis` | `stereocenters`, `double_bonds`, `n_stereocenters`, `n_double_bonds` |
| `Stereocenter` | `atom_index`, `element`, `substituent_indices`, `priorities`, `configuration` (`"R"`/`"S"`/`None`) |
| `DoubleBondStereo` | `atom1`, `atom2`, `configuration` (`"E"`/`"Z"`/`None`), `high_priority_sub1`, `high_priority_sub2` |
| `GraphFeatureAnalysis` | `aromatic_atoms`, `aromatic_bonds`, `tagged_tetrahedral_centers`, `inferred_tetrahedral_centers` |
| `MetalCoordinationCenter` | `atom_index`, `element`, `ligand_indices`, `coordination_number`, `geometry`, `geometry_score` |
| `TopologyAnalysisResult` | `metal_centers`, `warnings` |
| `UvVisSpectrum` | `energies_ev` (eV), `intensities`, `peaks`, `sigma` (eV), `notes` |
| `UvVisPeak` | `energy_ev` (eV), `wavelength_nm` (nm), `intensity`, `from_mo`, `to_mo` |
| `FrontierDescriptors` | `homo_atom_contributions`, `lumo_atom_contributions`, `dual_descriptor`, `homo_energy` (eV), `lumo_energy` (eV), `gap` (eV), `num_atoms` |
| `BondOrderResult` | `atom_pairs`, `distances` (Å), `wiberg`, `mayer`, `wiberg_valence`, `mayer_valence`, `num_atoms` |
| `FukuiDescriptors` | `f_plus`, `f_minus`, `f_radical`, `dual_descriptor`, `condensed_atom_indices`, `condensed_f_plus`, `condensed_f_minus`, `condensed_f_radical`, `condensed_dual_descriptor`, `homo_energy` (eV), `lumo_energy` (eV), `gap` (eV), `num_atoms`, `validity_notes` |
| `ReactivityRanking` | `nucleophilic_attack_sites`, `electrophilic_attack_sites`, `radical_attack_sites`, `notes` |
| `ReactivitySiteScore` | `atom_index`, `score` |

---

## Orbital meshes

### `compute_orbital_mesh`

```python
def compute_orbital_mesh(
    elements: list[int],
    coords: list[float],
    method: str = "eht",
    mo_index: int = 0,
    spacing: float = 0.2,
    padding: float = 3.0,
    isovalue: float = 0.02,
) -> OrbitalMeshResult
```

Isosurface mesh for any supported method (`eht`, `pm3`, `xtb`, `hf3c`). `spacing` and `padding`
in **Å**.

> **Deprecated alias:** `orbital_mesh` → `compute_orbital_mesh`.

### `OrbitalMeshResult`

| Attribute | Type | Description |
|---|---|---|
| `vertices` | `list[float]` | Flat mesh vertices in Å |
| `normals` | `list[float]` | Flat normals |
| `indices` | `list[int]` | Flat triangle indices |
| `num_triangles` | `int` | Triangle count |
| `method` | `str` | Method used |
| `mo_index` | `int` | Requested orbital index |
| `homo_index` | `int` | HOMO index |
| `orbital_energies` | `list[float]` | Orbital energies (eV) |
| `gap` | `float` | HOMO-LUMO gap (eV) |

---

## Materials and periodic systems

### `create_unit_cell`

```python
def create_unit_cell(a: float, b: float, c: float, alpha: float = 90.0, beta: float = 90.0, gamma: float = 90.0) -> UnitCell
```

Periodic unit cell from crystallographic parameters. Lengths in **Å**, angles in **degrees**,
`volume` in **Å³**.

> **Deprecated alias:** `unit_cell` → `create_unit_cell`.

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

Assemble a MOF-type framework. `topology` ∈ `{pcu, dia, sql}`; `geometry` ∈ `{linear, trigonal,
tetrahedral, square_planar, octahedral}`.

### `space_group_by_number`

```python
def space_group_by_number(number: int) -> SpaceGroup | None
```

Look up a crystallographic space group by its ITC number (1–230).

### `space_group_by_symbol`

```python
def space_group_by_symbol(symbol: str) -> SpaceGroup | None
```

Look up a space group by Hermann-Mauguin symbol, e.g. `"P 21/c"`.

### `all_space_groups`

```python
def all_space_groups() -> list[SpaceGroupEntry]
```

All 230 space groups as `(number, symbol)` records.

### `build_periodic_molecule`

```python
def build_periodic_molecule(
    elements: list[int],
    frac_coords: list[list[float]],
    a: float,
    b: float,
    c: float,
    alpha: float = 90.0,
    beta: float = 90.0,
    gamma: float = 90.0,
    bond_tolerance: float | None = None,
    options: dict | None = None,
) -> PeriodicMolecule
```

Build a periodic molecular graph from fractional coordinates and a unit cell. `frac_coords` is a
list of `[x, y, z]` triples; the lattice lengths are in **Å**.

### `detect_hapticity`

```python
def detect_hapticity(elements: list[int], coords: list[float], bonds: list[tuple[int, int]] | None = None) -> HapticityAnalysis
```

Detect haptic (η) metal–ring interactions and metallocene motifs. When `bonds` is `None` the
bonds are perceived from the geometry.

### `optimize_framework`

```python
def optimize_framework(
    smiles: str,
    coords: list[float],
    method: str = "bfgs",
    max_iter: int = 200,
    force_tol: float = 0.05,
    energy_tol: float = 1e-6,
    max_step: float = 0.2,
    fixed_atoms: list[int] | None = None,
    cell_a: float | None = None,
    cell_b: float | None = None,
    cell_c: float | None = None,
    cell_alpha: float = 90.0,
    cell_beta: float = 90.0,
    cell_gamma: float = 90.0,
    options: dict | None = None,
) -> FrameworkOptResult
```

Optimise a geometry with the shared BFGS engine (UFF backend). `energy` is reported in **eV**
and `forces` in **eV/Å**; `cell_a`/`cell_b`/`cell_c` enable periodic (minimum-image)
optimisation. `method` is `bfgs` or `steepest_descent`.

### `parse_cif` and `write_cif`

```python
def parse_cif(cif_text: str) -> CifStructure
def write_cif(structure: CifStructure) -> str
```

Read a CIF file into a crystal structure and write it back. The cell lengths are
in **Å**, the cell angles in **degrees** and the atom-site coordinates are
fractional (dimensionless).

```python
structure = sci_form.parse_cif(open("NaCl.cif").read())
print(structure.data_block, structure.a, len(structure.atom_sites))
text = sci_form.write_cif(structure)
```

A `CifStructure` can also be built directly:

```python
structure = sci_form.CifStructure(
    data_block="test",
    lattice=[[5.0, 0.0, 0.0], [0.0, 5.0, 0.0], [0.0, 0.0, 5.0]],
    a=5.0, b=5.0, c=5.0, alpha=90.0, beta=90.0, gamma=90.0,
)
```

### Return types

| Class | Attributes (unit) |
|---|---|
| `UnitCell` | `a`, `b`, `c` (Å), `alpha`, `beta`, `gamma` (degrees), `volume` (Å³), `lattice` |
| `CifStructure` | `data_block`, `lattice` (Å), `a`, `b`, `c` (Å), `alpha`, `beta`, `gamma` (degrees), `space_group_hm`, `space_group_number`, `atom_sites` |
| `CifAtomSite` | `label`, `type_symbol`, `atomic_number`, `frac_x`, `frac_y`, `frac_z`, `occupancy` |
| `CrystalStructure` | `elements`, `cart_coords` (Å), `frac_coords`, `labels`, `lattice`, `num_atoms` |
| `SpaceGroup` | `number`, `symbol`, `crystal_system`, `lattice_type`, `point_group`, `operations` |
| `SymmetryOp` | `rotation`, `translation` |
| `SpaceGroupEntry` | `number`, `symbol` |
| `PeriodicMolecule` | `elements`, `frac_coords`, `cart_coords` (Å), `bonds`, `cell`, `n_atoms` |
| `PeriodicBond` | `atom_i`, `atom_j`, `order`, `image`, `distance` (Å) |
| `HapticityAnalysis` | `interactions`, `is_metallocene`, `n_interactions` |
| `HapticInteraction` | `metal_index`, `metal_element`, `ring_atoms`, `hapticity`, `centroid` (Å), `metal_centroid_distance` (Å) |
| `FrameworkOptResult` | `positions` (Å), `energy` (eV), `forces` (eV/Å), `max_force` (eV/Å), `n_iterations`, `converged`, `energy_history` (eV) |

---

## Dynamics, reactions and transport

### `compute_md_trajectory`

```python
def compute_md_trajectory(smiles: str, coords: list[float], n_steps: int = 100, dt_fs: float = 1.0, seed: int = 42) -> MdTrajectory
```

Short exploratory NVE molecular dynamics (velocity Verlet + UFF). Times in **fs**, energies in
**kcal/mol**, temperatures in **K**.

> **Deprecated alias:** `md_trajectory` → `compute_md_trajectory`.

### `compute_md_trajectory_nvt`

```python
def compute_md_trajectory_nvt(
    smiles: str,
    coords: list[float],
    n_steps: int = 100,
    dt_fs: float = 1.0,
    seed: int = 42,
    target_temp_k: float = 300.0,
    thermostat_tau_fs: float = 100.0,
) -> MdTrajectory
```

NVT molecular dynamics with the Berendsen thermostat. `target_temp_k` in **K**,
`thermostat_tau_fs` in **fs**.

> **Deprecated alias:** `md_trajectory_nvt` → `compute_md_trajectory_nvt`.

### `compute_simplified_neb_path`

```python
def compute_simplified_neb_path(
    smiles: str,
    start_coords: list[float],
    end_coords: list[float],
    n_images: int = 8,
    n_iter: int = 50,
    spring_k: float = 0.1,
    step_size: float = 0.01,
) -> NebPathResult
```

Simplified NEB-like path between two geometries: linear interpolation plus spring-coupled UFF
gradient relaxation. This is a low-cost exploratory path tool, not a full climbing-image /
tangent-projected NEB.

> **Deprecated alias:** `simplified_neb_path` → `compute_simplified_neb_path`.

### `search_conformers_with_uff`

```python
def search_conformers_with_uff(smiles: str, n_samples: int = 10, seed: int = 42, rmsd_threshold: float = 0.5) -> ConformerSearchResult
```

Multi-seed conformer exploration ranked by UFF energy and deduplicated by Kabsch-aligned RMSD
(threshold in **Å**). Energies in **kcal/mol**.

> **Deprecated alias:** `search_conformers` → `search_conformers_with_uff`.

### `parse_smirks` and `apply_smirks`

```python
def parse_smirks(smirks: str) -> SmirksTransform
def apply_smirks(smirks: str, smiles: str) -> SmirksResult
```

Parse a SMIRKS reaction string and apply it to a molecule.

### `pack_conformers`

```python
def pack_conformers(results: list[ConformerResult]) -> RecordBatch
```

Pack conformer results into sci-form's columnar (Arrow-like) format for efficient transfer. This
is **not** the Apache Arrow IPC format.

### `split_worker_tasks`

```python
def split_worker_tasks(smiles: list[str], n_workers: int = 4, seed: int = 42) -> list[dict]
```

Split a SMILES list into balanced worker tasks for multi-process / web-worker dispatch.

### `estimate_workers`

```python
def estimate_workers(n_items: int, max_workers: int = 8) -> int
```

Estimate the optimal number of workers for a given batch size.

### Return types

| Class | Attributes (unit) |
|---|---|
| `MdTrajectory` | `frames`, `dt_fs` (fs), `notes` |
| `MdFrame` | `step`, `time_fs` (fs), `coords` (Å), `potential_energy_kcal_mol`, `kinetic_energy_kcal_mol` (kcal/mol), `temperature_k` (K) |
| `NebPathResult` | `images`, `notes` |
| `NebImage` | `index`, `coords` (Å), `potential_energy_kcal_mol` (kcal/mol) |
| `ConformerSearchResult` | `generated`, `unique`, `rotatable_bonds`, `conformers`, `clusters`, `notes` |
| `ConformerEnsembleMember` | `seed`, `cluster_id`, `coords` (Å), `energy_kcal_mol` (kcal/mol) |
| `ConformerClusterSummary` | `cluster_id`, `size`, `representative_seed`, `member_seeds` |
| `SmirksTransform` | `reactant_smarts`, `product_smarts`, `smirks` |
| `SmirksResult` | `products`, `atom_mapping`, `n_transforms`, `success`, `messages`, `is_ok` |
| `RecordBatch` | `num_rows`, `num_columns`, `byte_size`, `column_names`, `float_data`, `int_data`, `uint8_data` |

---

## System planning

### `system_capabilities`

```python
def system_capabilities(elements: list[int]) -> SystemCapabilities
```

Report which core workflows are usable for an element set.

### `system_method_plan`

```python
def system_method_plan(elements: list[int]) -> SystemMethodPlan
```

Recommend a method (and a fallback) per property domain.

### `compare_methods`

```python
def compare_methods(smiles: str, elements: list[int], coords: list[float], allow_experimental_eht: bool = False) -> MethodComparisonResult
```

Compare the available methods on the same system and geometry.

### Return types

| Class | Attributes |
|---|---|
| `SystemCapabilities` | `embed`, `uff`, `eht`, `population`, `orbital_grid` (each a `MethodCapability`) |
| `MethodCapability` | `available`, `confidence`, `unsupported_elements`, `warnings` |
| `SystemMethodPlan` | `capabilities`, `geometry`, `force_field_energy`, `orbitals`, `population`, `orbital_grid` (each a `PropertyMethodPlan`) |
| `PropertyMethodPlan` | `property`, `recommended`, `fallback`, `rationale`, `methods` |
| `MethodMetadata` | `method`, `available`, `confidence`, `confidence_score`, `warnings`, `limitations` |
| `MethodComparisonResult` | `plan`, `comparisons` |
| `MethodComparisonEntry` | `method`, `status`, `available`, `confidence`, `confidence_score`, `warnings`, `limitations`, `error`, `support_level`, `homo_energy` (eV), `lumo_energy` (eV), `gap` (eV), `energy_kcal_mol` (kcal/mol) |

---

## Options mappings

Eight over-argumented functions accept an optional trailing `options` mapping. The positional
form keeps working exactly as before; a key present in the mapping **overrides** the positional
value. An unknown key is an error, not a silent no-op.

| Function | Accepted option keys |
|---|---|
| `compute_fluorescence` | `sigma`, `e_min`, `e_max`, `n_points`, `stokes_ev`, `refractive_index`, `dark_threshold`, `broadening` |
| `compute_ir_spectrum` | `method`, `gamma`, `wn_min`, `wn_max`, `n_points`, `step_size`, `smiles` |
| `compute_ir_spectrum_broadened` | `method`, `gamma`, `wn_min`, `wn_max`, `n_points`, `broadening`, `step_size`, `smiles` |
| `compute_giao_nmr` | `nucleus`, `charge`, `multiplicity`, `max_scf_iter`, `allow_basis_fallback`, `use_parallel_eri` |
| `build_periodic_molecule` | `a`, `b`, `c`, `alpha`, `beta`, `gamma`, `bond_tolerance` |
| `optimize_framework` | `method`, `max_iter`, `force_tol`, `energy_tol`, `max_step`, `fixed_atoms`, `cell_a`, `cell_b`, `cell_c`, `cell_alpha`, `cell_beta`, `cell_gamma` |
| `train_random_forest` | `n_trees`, `max_depth`, `min_samples_split`, `min_samples_leaf`, `sample_fraction`, `seed` |
| `train_gradient_boosting` | `n_estimators`, `learning_rate`, `max_depth`, `min_samples_split`, `min_samples_leaf`, `subsample`, `seed` |

```python
result = sci_form.optimize_framework("CCO", coords, options={"max_iter": 5})
print(result.n_iterations)   # 5

sci_form.optimize_framework("CCO", coords, options={"nope": 1})
# ValueError: Invalid 'options' key 'nope'. Supported: 'method', 'max_iter', 'force_tol',
#             'energy_tol', 'max_step', 'fixed_atoms', 'cell_a', 'cell_b', 'cell_c',
#             'cell_alpha', 'cell_beta', 'cell_gamma'.
```

---

## Experimental, feature-gated surface

These functions only exist when the matching Cargo feature is enabled on the `sci-form-python`
package. They are still listed in `sci_form.__all__` when present.

| Capability | Feature | Canonical functions | Deprecated aliases |
|---|---|---|---|
| EEQ charges | `experimental-eeq` | `compute_eeq_charges`, `compute_eeq_energy` | `eeq_charges`, `eeq_energy` |
| ALPB solvation | `experimental-alpb` | `compute_alpb_solvation`, `compute_alpb_born_radii` | `alpb_solvation`, `alpb_born_radii` |
| D4 dispersion | `experimental-d4` | `compute_d4_energy` | `d4_energy` |
| CPM electrochemistry | `experimental-cpm` | `compute_cpm_charges`, `compute_cpm_surface` | `cpm_charges`, `cpm_surface` |
| Electrical double layer | `alpha-edl` | `compute_edl_profile`, `compute_edl_capacitance_scan`, `compute_helmholtz_profile`, `compute_gouy_chapman_profile`, `compute_gcs_profile` | `edl_profile`, `edl_capacitance_scan`, `edl_helmholtz`, `edl_gouy_chapman`, `edl_gcs` |
| Periodic linear algebra | `alpha-periodic-linear` | `compute_kmesh`, `compute_periodic_dos`, `solve_periodic_randnla`, `bloch_phase`, `build_bloch_hamiltonian`, `assemble_periodic_operators`, `bz_integrate_scalar`, `bz_integrate_vector`, `validate_electron_count`, `validate_kmesh_weights` | `kmesh` |
| Kinetics | `alpha-kinetics` | `compute_htst_rate`, `compute_htst_sweep`, `extract_kinetics_diagnostics`, `solve_microkinetic_network`, `solve_microkinetic_steady_state`, `analyze_gsm_mbh_htst_step` | `htst_rate`, `htst_temperature_sweep` |
| Chart utilities | `alpha-render-bridge` | `chart_to_json`, `chart_from_json`, `validate_chart_schema`, `pack_chart_payload`, `compute_edl_chart`, `compute_capacitance_chart`, `compute_arrhenius_chart`, `compute_thermal_prop_chart`, `compute_phase_portrait_chart`, `compute_band_structure_chart`, `compute_dos_chart`, `compute_kpoint_path_chart`, `compute_fermi_surface_chart`, `compute_trajectory_chart`, `compute_reaction_coordinate_chart` | `edl_chart`, `capacitance_chart`, `arrhenius_chart`, `thermal_prop_chart`, `phase_portrait_chart`, `band_structure_chart`, `density_of_states_chart`, `kpoint_path_chart`, `fermi_surface_chart`, `trajectory_chart`, `reaction_coordinate_chart` |

The deprecated aliases follow exactly the same rule as the stable ones: they still work, they emit
a `DeprecationWarning`, and they will be removed in sci-form 2.0.


## Submodules

### `sci_form.alpha` and `sci_form.beta`

Empty in a default build; populated when the matching `alpha-*` / `beta-*` features are enabled.
The `alpha` module covers DFT, ReaxFF, MLFF, CGA, GSM, SDR and live dynamics; `beta` covers KPM,
MBH, RandNLA, Riemannian methods and CPM.

`alpha` — canonical name and deprecated alias:

| Canonical | Deprecated alias |
|---|---|
| `alpha_compute_dft` | `dft_calculate` |
| `alpha_compute_reaxff_gradient` | `reaxff_gradient` |
| `alpha_compute_reaxff_energy` | `reaxff_energy` |
| `alpha_compute_eem_charges` | `eem_charges` |
| `alpha_rotate_dihedral_cga` | `rotate_dihedral_cga` |
| `alpha_refine_torsion_cga` | `refine_torsion_cga` |
| `alpha_gsm_interpolate` | `gsm_interpolate` |
| `alpha_gsm_backend_plan` | `gsm_backend_plan` |
| `alpha_gsm_compare_backends` | `gsm_compare_backends` |
| `alpha_gsm_find_ts` | `gsm_find_ts` |
| `alpha_sdr_embed` | `sdr_embed_py` |

`beta` — canonical name and deprecated alias:

| Canonical | Deprecated alias |
|---|---|
| `beta_compute_kpm_dos` | `kpm_dos` |
| `beta_compute_kpm_mulliken` | `kpm_mulliken` |
| `beta_compute_mbh_frequencies` | `mbh_frequencies` |
| `beta_solve_eht_randnla` | `eht_randnla` |
| `beta_compute_cpm_charges` | `cpm_charges` |

### `sci_form.potentials`

Analytic model potentials, useful as test targets for the reaction-dynamics code.

| Function | Signature |
|---|---|
| `double_well_1d` | `(x) -> float` |
| `asymmetric_1d` | `(x) -> float` |
| `sn2_model_1d` | `(x) -> float` |
| `muller_brown_2d` | `(x, y) -> float` |
| `proton_transfer_evb` | `(s, k, coupling) -> float` |
| `proton_transfer_2d_morse` | `(r_ah, r_hb) -> float` |

---

## Integration examples

### NumPy

```python
import numpy as np
from sci_form import embed

result = embed("CCO")
coords = np.array(result.coords).reshape(-1, 3)   # (9, 3), Å
elements = np.array(result.elements)               # (9,)
```

### Zero-copy embedding

```python
import numpy as np, sci_form

buf = np.zeros(3 * 9, dtype=np.float64)
sci_form.embed_coords_into("CCO", buf, 42)
print(buf.reshape(-1, 3))
```

### Pandas

```python
import pandas as pd
from sci_form import embed_batch

results = embed_batch(["CCO", "c1ccccc1", "CC(=O)O"])
df = pd.DataFrame([
    {"smiles": r.smiles, "num_atoms": r.num_atoms, "time_ms": r.time_ms, "ok": r.is_ok()}
    for r in results
])
```

### RDKit interop

```python
from rdkit import Chem
from sci_form import embed

result = embed("c1ccccc1")
mol = Chem.MolFromSmiles("c1ccccc1")
mol = Chem.AddHs(mol)
conf = Chem.Conformer(mol.GetNumAtoms())
for i in range(result.num_atoms):
    conf.SetAtomPosition(i, (result.coords[i*3], result.coords[i*3+1], result.coords[i*3+2]))
mol.AddConformer(conf, assignId=True)
```

### DOS + matplotlib

```python
import matplotlib.pyplot as plt
from sci_form import embed, compute_dos

conf = embed("c1ccccc1")
result = compute_dos(conf.elements, conf.coords, sigma=0.3)

plt.plot(result.energies, result.total_dos)
plt.xlabel("Energy (eV)")
plt.ylabel("DOS (states/eV)")
plt.show()
```

---

## See also

- [CLI reference](./cli.md)
- [TypeScript / JavaScript reference](./typescript.md)
- [Rust reference](./rust.md)
- [Complete supported reference](./complete-supported-reference.md)
- [Python guide](../guide/python.md)
