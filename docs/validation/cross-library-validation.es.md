# Cross-library validation plan

> **Historical route, English content.** This URL is kept for existing links.
> The original Spanish working notebook is archived outside the documentation
> tree; the maintained quantitative report is
> [Scientific validation](./scientific-validation.md), and the raw rows are in
> the [full matrix](./full-matrix.md).

<h2 id="indice">Index</h2>

1. [Purpose](#purpose)
2. [Validation levels](#validation-levels)
3. [Comparison rules](#comparison-rules)
4. [Tolerance guide](#tolerance-guide)
5. [Method families and gates](#method-families-and-gates)
6. [Traceability](#traceability)
7. [Execution phases](#execution-phases)
8. [Risks and anti-false-positive rules](#risks-and-anti-false-positive-rules)
9. [Current status](#current-status)
10. [Where the numbers live](#where-the-numbers-live)

## Purpose

For every capability in sci-form, the plan records:

1. which tested program or library performs the same task;
2. which test compares the two implementations;
3. which fixture or dataset stores the reference values;
4. how large the measured separation is, and where that measurement holds.

The plan is deliberately separate from the release narrative. A green row means
the observable agrees with the named reference inside the stated confidence
region. It does not mean the two methods are the same theory.

## Validation levels

| Level | What it proves | Example |
|---|---|---|
| **L0** | An invariant that must always hold | charge conservation; DOS integrates to the orbital count; analytic gradient equals a finite difference |
| **L1** | The same result on another surface | Rust vs Python vs WASM vs CLI with the same seed |
| **L2** | Same method, canonical implementation | UFF vs RDKit UFF; GFN2 vs tblite; PM3 vs MOPAC |
| **L3** | Same observable, independent theory | sTDA vs TD-SCF; EHT trends vs HF |
| **L4** | Experiment or a frozen dataset | IR frequencies vs experiment; NMR shifts vs a database |

Release-blocking rows are L0–L2. L3 and L4 rows document accuracy and method
gaps; they are not used to claim parity with a different theory.

## Comparison rules

1. **Same geometry first.** Energies, charges and gradients are compared on the
   same frozen coordinates. If both sides optimize, compare after re-optimizing
   with the same method.
2. **Explicit units.** EHT in eV, PM3 heat of formation in kcal/mol, xTB/HF/DFT
   in eV or Hartree, UFF/MMFF in kcal/mol, frequencies in cm⁻¹, shifts in ppm.
   Convert before comparing.
3. **Eigenstates.** Compare magnitudes or project degenerate subspaces. Never
   compare element-by-element signatures without aligning phase and order.
4. **Atom permutation.** Align by graph isomorphism, not by file index.
5. **Orientation and centre.** Use Kabsch/RMSD or compare internal distance
   matrices when global alignment is ambiguous.
6. **Randomness.** Fix the seed on every side; repeat 20–50 seeds and compare
   distributions, not a single draw.
7. **Non-uniqueness.** SSSR is not unique (compare the size multiset). Canonical
   SMILES is not bit-identical across toolkits (compare isomorphism and
   idempotence).
8. **Topology first.** Before comparing a property, verify that elements, bonds,
   bond orders and formal charges agree.

## Tolerance guide

| Observable | Level | Target tolerance | Note |
|---|---|---|---|
| Coordinates vs another surface, same seed | L1 | bit-identical when the path is identical; ≤1e-12 Å otherwise | compare internal coordinates |
| Conformers vs ETKDG | L4 | median heavy-atom RMSD ≤0.25 Å; P90 ≤0.5 Å | independent random streams |
| Bond lengths vs experiment/CSD | L4 | ≤0.02–0.03 Å after optimisation | not raw embedding |
| Overlap matrix, same STO-3G basis | L2 | 1e-10 abs vs PySCF | Gaussian integrals |
| EHT Hamiltonian, same VSIP table | L2 | 1e-8 abs | K = 1.75 |
| Electron-repulsion integrals | L2 | 1e-10 rel vs PySCF/libcint | s, p, d shells |
| HF eigenvalues | L2 | 1e-6 Ha vs PySCF | same basis |
| HF total energy | L2 | 1e-7 Ha vs PySCF | |
| DFT total energy | L2 | 1e-5 Ha vs PySCF | same grid |
| EHT orbital energy / gap | L2 | 1e-6 eV | |
| PM3 heat of formation | L2 | ≤0.5 kcal/mol vs MOPAC | same geometry |
| GFN1/GFN2 energy | L2 | ≤0.5 kcal/mol vs tblite | same geometry |
| GFN0 energy | L2 | ≤1 kcal/mol vs xtb | |
| D3/D4 dispersion | L2 | 1e-3 kcal/mol vs dftd3/dftd4 | same parameters |
| UFF/MMFF energy | L2 | ≤0.01 kcal/mol vs RDKit | same topology |
| Analytic gradient vs finite difference | L0 | relative error ≤1e-4 | h = 1e-3 Å |
| Gasteiger charges | L2 | ≤0.02 e vs RDKit | same iteration count |
| Mulliken charges | L2 | ≤0.005 e | same definition |
| Dipole | L4 | ≤0.1 D or 5% vs NIST/CCBDB | |
| SASA | L2 | ≤1% vs FreeSASA | same radii and probe |
| ESP grid | L2 | RMSE ≤1e-3 a.u. vs same density | |
| DOS/PDOS | L2 | MSE ≤0.01 vs exact diagonalization | KPM: ≤2–5% |
| Harmonic frequencies, same method | L2 | ≤2 cm⁻¹ | |
| Frequencies vs experiment | L4 | ±20 cm⁻¹ harmonic; ±10 cm⁻¹ scaled | document the scale factor |
| UV-Vis sTDA vs TD-DFT | L3 | ≤0.3 eV per dominant transition | |
| NMR GIAO vs ORCA | L2 | ≤0.5 ppm | same basis |
| NMR HOSE vs experiment | L4 | MAE ¹H ≤0.3 ppm; ¹³C ≤3 ppm | |
| J coupling vs experiment | L4 | ≤1 Hz | ensemble-averaged conformers |
| ANI energy/forces vs torchani | L2 | ≤1e-4 eV/atom; ≤1e-3 eV/Å | same weights |
| Descriptors vs RDKit/mordred | L2 | exact where definitions coincide; ±0.1 logP Crippen | |
| ECFP Tanimoto | L2 | same invariants; ±1e-9 similarity | ranking is the gate |
| KPM DOS vs diagonalization | L2 | L2 relative ≤5% | order 100–500 |

## Method families and gates

### Parsing, topology and conformers

| Gate | Capability | Reference | Fixture or test | Status |
|---|---|---|---|---|
| F01 | SMILES graph invariants | RDKit 2026.09.1 | `tests/cross_validation.rs` | PASS |
| F05 | Ring-size multiset | RDKit | same | PASS |
| F75 | SMARTS match sets | RDKit `GetSubstructMatches` | 837 CSD torsion patterns × 10 molecules | PASS |
| F81 | CIP stereo labels | RDKit CIP labeller | same 3D geometry | PASS |
| F83 | Canonical SMILES | RDKit `MolToSmiles` | isomorphism and round-trip | PASS with definition caveat |
| F71 | Distance bounds | RDKit `GetMoleculeBoundsMatrix` | per topological class | PASS with class tolerances |
| F72 | Embedding | RDKit ETKDGv3 | 29 molecules | PASS in distribution |
| F73 | ETKDG force field | RDKit force field | 19 shared geometries | PASS |
| F74 | BFGS minimisation | RDKit `ForceField::minimize` | 28 molecules, same start | PASS |

### Force fields and gradients

| Gate | Capability | Reference | Fixture or test | Status |
|---|---|---|---|---|
| F16 | UFF energy | RDKit UFF | 400 molecules | PASS |
| F17 | UFF atom types | RDKit AtomTyper | 33641 labels | PASS |
| F70 | UFF minimisation | RDKit UFF optimiser | 10 molecules | PASS |
| F100 | MMFF94 energy and charges | RDKit MMFF94 | 7 molecules | PASS |
| F80 | Analytic gradients | central differences | H₂O₂, OH, term level | PASS |
| F78 | ReaxFF bond orders | LAMMPS `reaxff` | published points | PASS |
| F42 | PM3 energy | MOPAC 23.2.5 | 24 molecules | PASS |
| F43 | GFN0 energy | xtb 6.7.1 | 18 cases | PASS |
| F44 | GFN1 energy | tblite 0.7.0 | 18 closed-shell + 5 open-shell | PASS |
| F45 | GFN2 energy | tblite 0.7.0 | 18 closed-shell + 5 open-shell | PASS |
| F46 | D3/D4 dispersion | dftd4 4.3.0 | 25 molecules | PASS with standalone D4 gap |

### Charges, electrostatics and solvation

| Gate | Capability | Reference | Fixture or test | Status |
|---|---|---|---|---|
| F21 | Gasteiger charges | RDKit | 251 molecules | PASS |
| F22 | EEQ charges and CN | multicharge 0.6.0 | 23 molecules | PASS |
| F87-EEQ2019 | EEQ2019 charges | independent NumPy oracle | 14 cases | PASS |
| F23 | Mulliken charges | PySCF 2.14.0 | 11 molecules | PASS |
| F26 | SCF dipole | PySCF | 11 molecules | PASS |
| F28 | SASA | FreeSASA | 60 molecules | PASS |
| F84 | ALPB solvation | tblite 0.7.0 | 16 solutes in water | PASS with model caveat |
| F78 | ESP grid | PySCF density | far-field points | PASS |
| F59 | ECFP similarity ranking | RDKit Morgan/ECFP4 | 780 pairs | PASS |

### Electronic structure

| Gate | Capability | Reference | Fixture or test | Status |
|---|---|---|---|---|
| F86-EHT | EHT overlap, Hamiltonian, orbitals | independent NumPy EHT | 9 molecules | PASS |
| F31 | EHT metal parameters | YAeHMOP table | 76 parameters | PASS |
| F36 | HF/SCF matrices and energy | PySCF | 11 molecules, STO-3G | PASS |
| F37 | Electron-repulsion integrals | PySCF/libcint | 13074 elements | PASS |
| F38 | UHF/ROHF energy and S² | PySCF | 8 cases | PASS |
| F39 | AO→MO transform | PySCF `ao2mo.full` | 13058 elements | PASS |
| F41 | CIS/TDA excitations | PySCF | 8 molecules | PASS |
| F77 | sTDA UV-Vis | PySCF TDA/TDHF | 7 molecules | PASS with method gap |
| F77 | Kohn-Sham DFT | PySCF | 6 molecules, SVWN/PBE | PASS with target gap |
| F79 | ML descriptors and forces | independent NumPy BP oracle | 8 molecules | PASS with weight gap |

### Spectroscopy

| Gate | Capability | Reference | Fixture or test | Status |
|---|---|---|---|---|
| F48 | IR frequencies | xtb GFN2 `vibspectrum` | 39 modes | PASS |
| F79b | IR intensities | PySCF Hessian/dipoles | 2 molecules | PASS with reference caveat |
| F54 | NMR J couplings | reference HLA calculator | 7 + 49 cases | PASS |
| F77 | UV-Vis | PySCF TDA/TDHF | 7 molecules | PASS with method gap |
| F47 | EAS site ranking | textbook/DFT rules | 8 substrates | PASS |

### Materials and periodic systems

| Gate | Capability | Reference | Fixture or test | Status |
|---|---|---|---|---|
| F64 | Space groups | spglib 2.8.0 | 230 cells | PASS |
| F65 | CIF round-trip | pymatgen 2026.9.24 | 4 structures | PASS |
| F66 | Minimum-image distance | ASE 3.29.0 | 4 cells | PASS |
| F67 | MOF-5 assembly and topology | pormake / RCSR | one framework, 16 topologies | PASS |
| F68 | Porosity | pyzeo / Zeo++ | MOF-5 | PASS with one-framework limit |
| F69 | Powder XRD | pymatgen | 4 structures | PASS |
| F85 | Variable-cell stress | ASE Lennard-Jones + UnitCellFilter | 3 cells | PASS with convention note |

### Dynamics, reaction paths and electrochemistry

| Gate | Capability | Reference | Fixture or test | Status |
|---|---|---|---|---|
| F82 | Growing String Method | ASE CI-NEB and analytic surfaces | Müller-Brown + 5 wells | PASS |
| F87 | Kinetics / transition-state theory | Eyring and SciPy `solve_ivp` | small networks | PASS |
| F84 | Live molecular dynamics | RDKit UFF forces | UFF systems | PASS |
| F86 | Electrical double layer | analytic Gouy-Chapman/Stern | 3T × 3I × 6 potentials | PASS |
| F95 | Constant Potential Method | LAMMPS `electrode/conp` | graphene bilayer | PASS with hardness gap |
| F92 | Mobile Block Hessian | analytic diatomic + full Hessian | diatomic/benzene | PASS with experimental gap |

### Experimental modules promoted to core

| Module | Core path | Promotion evidence | Main limitation |
|---|---|---|---|
| Conformal geometric algebra | `cga` | motors 3.4e-15; dihedrals 6.6e-15 | refinement and MOF assembly not covered |
| Semidefinite relaxation | `sdr` | projections 4.4e-16; embeds 7e-15 | SDP scaling not benchmarked |
| Kernel Polynomial Method | `kpm` | moments 3.5e-15; DOS 1.7e-14 | T=0 gapless Gibbs |
| Randomized eigensolver | `rand_nla` | rank-N/4 3.8e-4; N=600 9.6e-4 | full-spectrum fallback |
| Riemannian geometry | `riemannian` | projection 8.2e-12 | Frobenius proxy, first-order retraction |
| Growing String Method | `gsm` | barrier 1.6e-3 vs ASE | no pyGSM parity |
| Kinetics | `kinetics` | Eyring 1.9e-16 | no Cantera I/O |
| Periodic adapters | `periodic_linear` | phases 3.6e-16; bands 4.4e-16 | complex H(k) not implemented |
| Live MD | `dynamics_live` | forces 3.6e-13 vs RDKit | short trajectories |

## Traceability

The maintained traceability record is the
[full validation matrix](./full-matrix.md): 129 rows, 122 PASS, 7 method or
provenance gaps, 0 FAIL. Each row names the gate, the observable, the reference,
the metric, the measured value, the tolerance, the level and the owning test.

The curated interpretation is in
[Scientific validation](./scientific-validation.md). It groups the rows by
method family and states the confidence region and limitation of each method.

## Execution phases

| Phase | Work | Output |
|---|---|---|
| 1 | Pin reference versions and build the oracle environment | `validation/refs.lock.json` |
| 2 | Write independent oracles for the observables | `scripts/validation/oracle_*.py` |
| 3 | Regenerate fixtures and run the Rust gates | `tests/*_validation.rs`, `tests/cross_validation.rs` |
| 4 | Collect the rows and review failures | `reports/validation/summary.json` |
| 5 | Publish the report and the full matrix | this documentation site |

## Risks and anti-false-positive rules

- **Same-method first.** Compare against the canonical implementation of the
  same method before comparing against experiment.
- **Frozen geometry.** A different optimiser path is not a method failure.
- **Degenerate subspaces.** Compare projected subspaces, not arbitrary
  eigenvector order.
- **Random seeds.** A single draw is not a distribution.
- **Unit conversion.** Most historical false alarms were unit mismatches, not
  missing physics.
- **Definition drift.** Canonical SMILES, HBA counts and IR intensities are
  definition-dependent; record the definition next to the number.

## Current status

| Item | Value |
|---|---|
| Release audited | `sci-form 2.0.0` |
| Recorded rows | 129 |
| Passing gates | 122 |
| Method or provenance gaps | 7 |
| Failing gates | 0 |
| Oracle scripts | 36/36 |
| Rust validation gates | 117/117 |
| Oracle runtime | Python 3.14.4 |

The seven non-passing rows are provenance or method-gap markers, not silent
failures. Their reasons are listed in
[Scientific validation](./scientific-validation.md#recorded-gaps).

## Where the numbers live

- [Scientific validation](./scientific-validation.md) — curated report with
  reference, deviation, confidence region and limitation.
- [Full matrix](./full-matrix.md) — 129 recorded rows with tolerances.
- `reports/validation/summary.json` — machine-readable source of truth.
- `validation/refs.lock.json` — pinned reference versions.
