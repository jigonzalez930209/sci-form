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

### `charges_configured`

```python
def charges_configured(
    smiles: str,
    max_iter: int = 6,
    initial_damping: float = 1.0,
    convergence_threshold: float = 1e-8,
) -> ChargeResult
```

Gasteiger-Marsili charges with full configuration control. Covers all main-group elements up to period 4 (Li, Be, Na, Mg, Al, Si, …) and handles formal charges on charged species.

```python
from sci_form import charges_configured

# Charged species
result = charges_configured("[NH4+]", max_iter=12, convergence_threshold=1e-10)
print(f"Total charge: {result.total_charge:.2f}")  # ~+1.0
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
) -> VibrationalAnalysisPy
```

Numerical Hessian ($6N$ energy evaluations with mass-aware step size) + diagonalization. Removes translation/rotation modes automatically. Returns `VibrationalAnalysisPy` with `.n_atoms`, `.modes` (list of `VibrationalModePy`), `.n_real_modes`, `.zpve_ev`, `.method`, `.notes`, `.thermochemistry`.

Each `VibrationalModePy` has `.frequency_cm1`, `.ir_intensity` (km/mol), `.displacement` (list, length 3N), `.is_real`, `.label` (functional group annotation).

`ThermochemistryPy` (RRHO, 298.15 K): `.zpve_kcal`, `.thermal_energy_kcal`, `.entropy_vib_cal`, `.gibbs_correction_kcal`.

```python
from sci_form import embed, vibrational_analysis

conf = embed("CCO", seed=42)
vib  = vibrational_analysis(conf.elements, conf.coords, method="xtb")
print(f"ZPVE: {vib.zpve_ev:.4f} eV, {vib.n_real_modes} real modes")
print(f"Gibbs correction: {vib.thermochemistry.gibbs_correction_kcal:.3f} kcal/mol")

real = [m for m in vib.modes if m.is_real]
strongest = max(real, key=lambda m: m.ir_intensity)
print(f"Strongest IR: {strongest.frequency_cm1:.1f} cm⁻¹  "
      f"({strongest.ir_intensity:.1f} km/mol)  [{strongest.label or '?'}]")
```

### `ir_spectrum`

```python
def ir_spectrum(
    elements: list[int],
    coords: list[float],
    method: str = "xtb",
) -> IrSpectrumPy
```

Convenience wrapper: runs `vibrational_analysis` and applies default Lorentzian broadening. Returns `IrSpectrumPy` with `.wavenumbers`, `.intensities`, `.transmittance`, `.peaks` (list of `IrPeakPy`), `.gamma`.

```python
from sci_form import embed, ir_spectrum

conf = embed("CCO", seed=42)
spec = ir_spectrum(conf.elements, conf.coords, method="xtb")
import matplotlib.pyplot as plt
plt.plot(spec.wavenumbers, spec.intensities)
plt.xlabel("Wavenumber (cm⁻¹)")
plt.gca().invert_xaxis()
plt.show()
```

### `ir_spectrum_broadened`

```python
def ir_spectrum_broadened(
    analysis: VibrationalAnalysisPy,
    gamma: float = 15.0,
    wn_min: float = 600.0,
    wn_max: float = 4000.0,
    n_points: int = 1000,
    broadening: str = "lorentzian",  # "lorentzian" or "gaussian"
) -> IrSpectrumPy
```

Broadened IR spectrum from an existing `VibrationalAnalysisPy` with full parameter control. Includes both absorbance (`.intensities`) and transmittance (`.transmittance`) axes.

```python
from sci_form import ir_spectrum_broadened

# Gaussian broadening for publication-quality output
spec = ir_spectrum_broadened(vib, gamma=10.0, wn_min=400.0, wn_max=4000.0,
                              n_points=2000, broadening="gaussian")
print(f"{len(spec.peaks)} labelled peaks")
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

Predict H–H J-coupling constants. With 3D coords uses the Karplus equation for ³J; topology only gives free-rotation averages (~5.3 Hz). Parameterized pathways: H-C-C-H (Altona-Sundaralingam), H-C-N-H (Bystrov), H-C-O-H, H-C-S-H.

Each `JCouplingPy` has `.h1_index`, `.h2_index`, `.j_hz`, `.n_bonds`, `.coupling_type`.

```python
from sci_form import embed, nmr_couplings

conf = embed("CC", seed=42)
couplings = nmr_couplings("CC", conf.coords)
for j in couplings:
    print(f"{j.n_bonds}J(H{j.h1_index},H{j.h2_index}) = {j.j_hz:.2f} Hz  [{j.coupling_type}]")
```

### `ensemble_j_couplings`

```python
def ensemble_j_couplings(
    smiles: str,
    conformer_coords: list[list[float]],
    energies_kcal: list[float],
    temperature_k: float = 298.15,
) -> list[JCouplingPy]
```

Boltzmann-average ³J couplings over a conformer ensemble.

```python
from sci_form import embed_batch, ensemble_j_couplings

results = embed_batch(["CCCC"] * 20, seed=42)
coords = [r.coords for r in results if r.is_ok()]
energies = [0.0] * len(coords)  # or compute from uff_energy
couplings = ensemble_j_couplings("CCCC", coords, energies, temperature_k=298.15)
for jc in couplings:
    print(f"H{jc.h1_index}–H{jc.h2_index}: {jc.j_hz:.2f} Hz")
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

## Semi-Empirical QM

### `pm3_calculate`

```python
def pm3_calculate(
    elements: list[int],  # atomic numbers
    coords: list[float],  # flat [x0,y0,z0,...]
) -> Pm3Result
```

PM3 NDDO semi-empirical SCF. Supports transition metals via PM3(tm): Ti,Cr,Mn,Fe,Co,Ni,Cu,Zn,Al,Si,Ge.

Returns `Pm3Result` with `.orbital_energies`, `.electronic_energy`, `.nuclear_repulsion`, `.total_energy`, `.heat_of_formation`, `.n_basis`, `.n_electrons`, `.homo_energy`, `.lumo_energy`, `.gap`, `.mulliken_charges`, `.scf_iterations`, `.converged`.

```python
from sci_form import embed, pm3_calculate

conf = embed("CCO", seed=42)
pm3 = pm3_calculate(conf.elements, conf.coords)
print(f"HOF: {pm3.heat_of_formation:.2f} kcal/mol, gap: {pm3.gap:.3f} eV")
```

### `xtb_calculate`

```python
def xtb_calculate(
    elements: list[int],
    coords: list[float],
) -> XtbResult
```

GFN0-xTB tight-binding with self-consistent charges. Supports 25 elements including transition metals (Ti,Cr,Mn,Fe,Co,Ni,Cu,Zn,Ru,Pd,Ag,Pt,Au).

Returns `XtbResult` with `.orbital_energies`, `.electronic_energy`, `.repulsive_energy`, `.total_energy`, `.n_basis`, `.n_electrons`, `.homo_energy`, `.lumo_energy`, `.gap`, `.mulliken_charges`, `.scc_iterations`, `.converged`.

```python
from sci_form import embed, xtb_calculate

conf = embed("CCO", seed=42)
xtb = xtb_calculate(conf.elements, conf.coords)
print(f"xTB gap: {xtb.gap:.3f} eV, converged: {xtb.converged}")
```

---

## Ab-initio / Neural Potentials

### `hf3c_calculate`

```python
def hf3c_calculate(
    elements: list[int],
    coords: list[float],
    max_scf_iter: int = 100,
    n_cis_states: int = 5,
    corrections: bool = True,
) -> Hf3cResultPy
```

Minimal-basis Hartree-Fock with D3, gCP, and SRB corrections.

### `ani_calculate`

```python
def ani_calculate(
    elements: list[int],
    coords: list[float],
) -> AniResultPy
```

ANI-2x neural network potential for H, C, N, O, F, S, Cl.

---

## Bond Orders & Frontier Analysis

### `bond_orders`

```python
def bond_orders(
    elements: list[int],
    coords: list[float],
) -> BondOrderResultPy
```

Wiberg and Mayer bond orders from EHT density matrix. Returns `.atom_pairs`, `.distances`, `.wiberg`, `.mayer`, `.wiberg_valence`, `.mayer_valence`.

### `frontier_descriptors`

```python
def frontier_descriptors(
    elements: list[int],
    coords: list[float],
) -> FrontierDescriptorsPy
```

HOMO/LUMO atom contributions and dual descriptor. Returns `.homo_atom_contributions`, `.lumo_atom_contributions`, `.dual_descriptor`, `.homo_energy`, `.lumo_energy`, `.gap`.

---

## Reactivity

### `fukui_descriptors`

```python
def fukui_descriptors(
    elements: list[int],
    coords: list[float],
) -> FukuiDescriptorsPy
```

Condensed Fukui functions f⁺, f⁻, f⁰. Returns `.condensed_f_plus`, `.condensed_f_minus`, `.condensed_f_radical`, `.condensed_dual_descriptor`, `.gap`, `.validity_notes`.

### `reactivity_ranking`

```python
def reactivity_ranking(
    elements: list[int],
    coords: list[float],
) -> ReactivityRankingPy
```

Ranked atomic sites for nucleophilic, electrophilic, and radical attack.

### `empirical_pka`

```python
def empirical_pka(smiles: str) -> EmpiricalPkaResultPy
```

Graph-heuristic pKa estimation. Returns `.acidic_sites`, `.basic_sites`.

---

## ML Properties

### `ml_descriptors`

```python
def ml_descriptors(smiles: str) -> MolecularDescriptors
```

Topology-derived descriptors. No 3D required. Returns `.molecular_weight`, `.n_heavy_atoms`, `.n_hbd`, `.n_hba`, `.fsp3`, `.wiener_index`, etc.

### `ml_predict`

```python
def ml_predict(smiles: str) -> MlPropertyResult
```

Lipinski/ADME prediction. Returns `.logp`, `.molar_refractivity`, `.log_solubility`, `.druglikeness`, `.lipinski_violations`, `.lipinski_passes`.

---

## Stereochemistry

### `stereo_analysis`

```python
def stereo_analysis(
    smiles: str,
    coords: list[float] = [],  # flat [x0,y0,z0,...]; empty → topology-only
) -> StereoAnalysisPy
```

Assign CIP priorities and determine R/S at stereocenters and E/Z at double bonds.

```python
class StereoAnalysisPy:
    stereocenters: list[StereocenterPy]
    double_bonds:  list[DoubleBondStereoPy]
    n_stereocenters: int
    n_double_bonds: int

class StereocenterPy:
    atom_index: int
    element: int
    substituent_indices: list[int]
    priorities: list[int]       # CIP rank (1 = highest)
    configuration: str | None   # "R" or "S"

class DoubleBondStereoPy:
    atom1: int
    atom2: int
    configuration: str | None   # "E" or "Z"
    high_priority_sub1: int | None
    high_priority_sub2: int | None
```

**Example:**

```python
from sci_form import embed, stereo_analysis

# Chiral center
conf = embed("C(F)(Cl)(Br)I", seed=42)
stereo = stereo_analysis("C(F)(Cl)(Br)I", conf.coords)
print(f"{stereo.n_stereocenters} stereocenters")
for sc in stereo.stereocenters:
    print(f"  Atom {sc.atom_index}: {sc.configuration}")

# Double bond geometry
stereo2 = stereo_analysis("C/C=C/C")  # trans-2-butene
print(f"E/Z: {stereo2.double_bonds[0].configuration}")
```

---

## Implicit Solvation

### `nonpolar_solvation`

```python
def nonpolar_solvation(
    elements: list[int],
    coords: list[float],
    probe_radius: float = 1.4,
) -> NonPolarSolvationPy
```

Non-polar solvation $\Delta G_{\mathrm{np}} = \sum_i \sigma_i A_i$ using per-element ASP and Shrake-Rupley SASA.

```python
class NonPolarSolvationPy:
    energy_kcal_mol: float
    atom_contributions: list[float]
    atom_sasa: list[float]
    total_sasa: float
```

### `gb_solvation`

```python
def gb_solvation(
    elements: list[int],
    coords: list[float],
    charges: list[float],
    solvent_dielectric: float = 78.5,
    solute_dielectric: float = 1.0,
    probe_radius: float = 1.4,
) -> GbSolvationPy
```

GB/SA electrostatic solvation (HCT Born radii + Still GB equation).

```python
class GbSolvationPy:
    electrostatic_energy_kcal_mol: float
    nonpolar_energy_kcal_mol: float
    total_energy_kcal_mol: float
    born_radii: list[float]
    charges: list[float]
    solvent_dielectric: float
    solute_dielectric: float
```

**Example:**

```python
from sci_form import embed, charges as gasteiger_charges, nonpolar_solvation, gb_solvation

conf = embed("CCO", seed=42)
q = gasteiger_charges("CCO").charges

np_solv = nonpolar_solvation(conf.elements, conf.coords)
print(f"Non-polar ΔG: {np_solv.energy_kcal_mol:.2f} kcal/mol")

gb = gb_solvation(conf.elements, conf.coords, q)
print(f"Total solvation: {gb.total_energy_kcal_mol:.2f} kcal/mol")
```

---

## Ring Perception

### `sssr`

```python
def sssr(smiles: str) -> SssrResultPy
```

Smallest Set of Smallest Rings (Horton's algorithm) with Hückel aromaticity.

```python
class SssrResultPy:
    rings: list[RingInfoPy]
    atom_ring_count: list[int]          # rings containing each atom
    atom_ring_sizes: list[list[int]]    # ring sizes per atom
    ring_size_histogram: dict[int, int] # size → count

class RingInfoPy:
    atoms: list[int]    # atom indices
    size: int
    is_aromatic: bool
```

**Example:**

```python
from sci_form import sssr

r = sssr("c1ccc2ccccc2c1")  # naphthalene
print(f"{r.rings[0].size}-membered ring, aromatic={r.rings[0].is_aromatic}")
print(f"Histogram: {r.ring_size_histogram}")  # {6: 2}
```

---

## Fingerprints (ECFP)

### `ecfp`

```python
def ecfp(
    smiles: str,
    radius: int = 2,   # ECFP4 = radius 2
    n_bits: int = 2048,
) -> ECFPFingerprintPy
```

Extended-Connectivity Fingerprint (Morgan algorithm). Atom invariants: element, degree, valence, ring membership, aromaticity, formal charge.

```python
class ECFPFingerprintPy:
    n_bits: int
    on_bits: list[int]   # set bit indices
    radius: int
    density: float       # fraction of set bits
```

### `tanimoto`

```python
def tanimoto(
    fp1: ECFPFingerprintPy,
    fp2: ECFPFingerprintPy,
) -> float
```

Jaccard-Tanimoto similarity $T = |A \cap B| / |A \cup B|$.

**Example:**

```python
from sci_form import ecfp, tanimoto

fp1 = ecfp("c1ccccc1", radius=2, n_bits=2048)
fp2 = ecfp("Cc1ccccc1", radius=2, n_bits=2048)
print(f"Tanimoto: {tanimoto(fp1, fp2):.3f}")  # ~0.5–0.7
print(f"Self-sim: {tanimoto(fp1, fp1):.3f}")  # 1.0
```

---

## Conformer Clustering

### `butina_cluster`

```python
def butina_cluster(
    conformers: list[list[float]],  # list of flat coord arrays
    rmsd_cutoff: float = 1.0,
) -> ClusterResultPy
```

Taylor-Butina single-linkage RMSD clustering.

```python
class ClusterResultPy:
    n_clusters: int
    assignments: list[int]        # cluster index per conformer
    centroid_indices: list[int]   # representative conformer per cluster
    cluster_sizes: list[int]
    rmsd_cutoff: float
```

### `rmsd_matrix`

```python
def rmsd_matrix(conformers: list[list[float]]) -> list[list[float]]
```

O(N²) pairwise RMSD matrix with Kabsch alignment.

### `filter_diverse`

```python
def filter_diverse(
    conformers: list[list[float]],
    rmsd_cutoff: float = 1.0,
) -> list[int]
```

Return indices of cluster centroid (diverse representative) conformers.

**Example:**

```python
from sci_form import embed_batch, butina_cluster, filter_diverse

results = embed_batch(["CCCC"] * 20, seed=42)
coords = [r.coords for r in results if r.is_ok()]

clusters = butina_cluster(coords, rmsd_cutoff=0.5)
print(f"{clusters.n_clusters} clusters, sizes: {clusters.cluster_sizes}")

diverse_idx = filter_diverse(coords, rmsd_cutoff=1.0)
print(f"Diverse set: {len(diverse_idx)} conformers")
```

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

