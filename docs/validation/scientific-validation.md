# Scientific validation

This page is the quantitative companion to the
[cross-library plan](./cross-library-validation.es.md). It answers four
questions for each method:

1. **What was it compared with?** A canonical implementation, an independent
   method, or an analytic closed form.
2. **How far is it from the reference?** The measured separation, in the units
   of the observable.
3. **Where does that measurement hold?** The confidence region: molecule set,
   element coverage, system size, basis or parameter set, and thermodynamic
   state.
4. **What is still missing?** The limitations that the gate does not cover.

The raw rows, including every tolerance and the owning test, are in the
[full validation matrix](./full-matrix.md).

## Campaign snapshot

| Item | Value |
|---|---|
| Release audited | `sci-form 2.0.0` |
| Campaign date | 2026-11-07, with the MMFF94 addendum on 2026-11-08 |
| Recorded rows | 129 |
| Passing gates | 122 |
| Method or provenance gaps | 7 |
| Failing gates | 0 |
| Oracle scripts | 36/36 |
| Rust validation gates | 117/117 |
| Oracle runtime | Python 3.14.4 |

Reference versions are pinned in `validation/refs.lock.json`:

| Tool | Version | Used for |
|---|---|---|
| RDKit | 2026.09.1 | parsing, rings, SMARTS, UFF, MMFF94, ETKDG, CIP, descriptors |
| PySCF | 2.14.0 | SCF, ERI, AO→MO, CIS/TDA, ESP, Hessian, IR intensities |
| tblite | 0.7.0 | GFN1/GFN2 and GB/ALPB |
| dftd4 | 4.3.0 | D3/D4 dispersion |
| xtb | 6.7.1 | GFN0 and GFN2 vibration reference |
| multicharge | 0.6.0 | EEQ charges |
| MOPAC | 23.2.5 | PM3 |
| LAMMPS | 29Aug2024_update2 | ReaxFF, CPM electrode |
| ASE | 3.29.0 | CI-NEB, k-meshes, stress |
| pymatgen | 2026.9.24 | CIF and XRD |
| pormake | 0.2.3 | MOF assembly |
| pyzeo / Zeo++ | 0.2.0 | porosity |
| NumPy / SciPy | 2.5.4 / 1.18.1 | independent oracles |
| FreeSASA | python-bindings | solvent-accessible surface |

## How to read the numbers

### Validation levels

| Level | What it proves | Example |
|---|---|---|
| **L0** | An invariant that must always hold | charge conservation, DOS integrates to the orbital count, analytic gradient equals a finite difference |
| **L1** | The same result on another surface | Rust vs Python vs WASM vs CLI with the same seed |
| **L2** | Same method, canonical implementation | UFF vs RDKit UFF; GFN2 vs tblite; PM3 vs MOPAC |
| **L3** | Same observable, independent theory | sTDA vs TD-SCF; EHT trends vs HF |
| **L4** | Experiment or a frozen dataset | IR frequencies vs experiment; NMR shifts vs a database |

Most release-blocking rows are L0–L2. L3/L4 rows document accuracy and method
gaps; they are not used to claim parity with a different theory.

### Distance

"Measured separation" is the metric named in the row: maximum absolute error,
maximum relative error, mismatched cases, Spearman correlation, or a fraction.
The units are explicit. A value of `0` for a mismatch metric means every case
agreed; it is not a floating-point zero.

### Confidence region

The confidence region is the **tested input envelope**, not a statistical
confidence interval. It records the molecules, elements, sizes, basis sets,
parameter sets and thermodynamic states over which the deviation was measured.
Inside that region the number in the table is the number to expect. Outside it,
the deviation is not characterized and should not be inferred.

### Limitations

A limitation is a known reason the gate cannot be read as a broader claim. It
can be a method gap (sTDA is not TD-DFT), a coverage gap (seven MMFF94
molecules), a boundary condition (no periodic EEQ), or a reference gap (no
external HOSE database). A limitation is not a hidden failure; it is the edge
of the confidence region.

## Parsing, topology and conformers

| Method | Compared with | Measured separation | Confidence region | Main limitation |
|---|---|---|---|---|
| SMILES graph parsing | RDKit 2026.09.1 | 0 mismatched cases / 159 | neutral organic molecules; heavy-atom graph, bond multiset, formal charges, hybridization | one multi-fragment counterion case skipped; no arbitrary valence or organometallic coverage |
| Ring perception | RDKit 2026.09.1 | 0 mismatched ring-size multisets / 159 | same corpus | SSSR is not unique; the gate compares the multiset, not one canonical ring set |
| SMARTS matching | RDKit `GetSubstructMatches` | 100% canonical match sets; 99.94% raw counts / 8370 pairs | 837 CSD torsion patterns × 10 molecules | arbitrary user SMARTS and recursive patterns are not covered |
| CIP stereo labels | RDKit CIP labeller | 0 mismatched labels / 19 | R/S and E/Z derived from the same 3D geometry | SMILES-tag-only labels differ by design; no axial or planar chirality |
| Canonical SMILES | RDKit `MolToSmiles(canonical=True)` | 25% identical strings by design; 0 round-trip/reordering mismatches / 36 | graph isomorphism and idempotence | string equality across toolkits is not a goal; the canonical form is internal |
| Distance bounds | RDKit `GetMoleculeBoundsMatrix` | topo-1 0.0259 Å; topo-2 0.337 Å; topo-3 0.373 Å; topo-5 0.537 Å; topo-6 0.95 Å | 29 molecules, per topological class | tolerances are class-dependent; 4 divergent pairs out of 388 at topo-4 |
| Embedding | RDKit ETKDGv3 | heavy-atom RMSD mean 0.0286 Å; bond length 0.505 of the bounds window | 29 molecules; independent random streams | coordinates are not expected to match; only connectivity and distribution are gated |
| ETKDG force field | RDKit force field | energy 9.1e-13 kcal/mol; inversion 4.3e-14 kcal/mol | 19 shared perturbed geometries | 1-2/1-3 constraints and inversions; long-range constraints are reported but not gated |
| BFGS minimization | RDKit `ForceField::minimize` | energy 0.0057 kcal/mol; RMSD 0.1355 Å | 28 molecules, same start and constraints | two local minima can be different even when both are valid |
| Butina clustering | RDKit Butina | 0 mismatched partitions / 12 | Euclidean point clouds | cluster labels are compared up to relabeling |

## Force fields and gradients

| Method | Compared with | Measured separation | Confidence region | Main limitation |
|---|---|---|---|---|
| UFF energy | RDKit UFF | 2.27e-13 kcal/mol max abs | 400 molecules | organic element types; no periodic UFF benchmark |
| UFF atom typing | RDKit UFF AtomTyper | 0 mismatched labels / 33641 | per-atom hybridization, conjugation and valence flags | exotic valences are not covered |
| UFF minimization | RDKit UFF optimiser | 7.84e-10 kcal/mol/atom | 10 molecules, same start | local minima and convergence paths can differ |
| MMFF94 energy | RDKit MMFF94 | 2.806e-11 kcal/mol; charges exactly 0 | 7 molecules | small molecule set; only MMFF94 types present in the fixture |
| Analytic gradients | central differences of the same energy | UFF 4.0e-8; MMFF94 term 8.3e-8; EHT 3.1e-8; ANI 8.3e-7; SCF nuclear 3.9e-10 | H₂O₂, OH and term-level checks | not every element or every energy term has a gradient gate |
| ReaxFF | LAMMPS `reaxff` | bond-order parity on the published points; analytic gradients self-consistent | CHON parameter set | no full LAMMPS energy parity; parameter coverage is limited to the default CHON set |
| PM3 | MOPAC 23.2.5 | 0.0083 kcal/mol max abs | 24 molecules, same seed geometry | the fixture is geometry-dependent; `Z > 17` uses the legacy NDDO fallback |
| GFN1-xTB | tblite 0.7.0 | 8.06e-5 eV; UKS 0.0336 eV | 18 closed-shell / 5 open-shell cases | no periodic GFN1; no element-wide coverage |
| GFN2-xTB | tblite 0.7.0 | 2.25e-5 eV; UKS 0.0346 eV | 18 closed-shell / 5 open-shell cases | no periodic GFN2; no element-wide coverage |
| GFN0-xTB | xtb 6.7.1 | 1.65e-4 eV | 18 cases | closed shell only |
| D3 / D4 dispersion | dftd4 4.3.0 | D3(HF-3c) 6.8e-9 Ha; D4(PBE) 6.2e-10 Ha | 25 molecules | the standalone D4 path uses element-level C6 coefficients without coordination or charge dependence |

The PM3 and xTB gradients are not symbolic derivatives of every NDDO/GFN2
integral. PM3 differentiates a **stationary-density** energy (no SCF per
displacement); xTB uses **warm-started central differences**. Both are
documented as method choices, not as exact analytic gradients.

## Charges, electrostatics and solvation

| Method | Compared with | Measured separation | Confidence region | Main limitation |
|---|---|---|---|---|
| Gasteiger charges | RDKit | 8.33e-17 e | 251 neutral organic molecules | no metals, no periodic systems |
| EEQ charges | multicharge 0.6.0 | 6.90e-12 e; coordination numbers 1.91e-11 | 23 molecules | main-group plus selected metals; no periodic boundary conditions |
| EEQ2019 charges | independent NumPy oracle | 6.08e-7 e; total charge 8.33e-17 e | 14 cases | same equations, not an experimental reference |
| Mulliken charges | PySCF 2.14.0 | 8.56e-8 e | 11 molecules, STO-3G | basis-set dependent; no NPA/NBO |
| SCF dipole | PySCF 2.14.0 | 1.83e-7 D | 11 molecules, STO-3G | no correlation; minimal basis |
| SASA | FreeSASA | 1.53% maximum relative error | same Bondi radii, probe 1.4 Å, 960 points; 60 molecules | residual is a point-set difference, not a model difference |
| ALPB solvation | tblite 0.7.0 | electrostatic 3.14 kcal/mol (tolerance 5); total 11.66 (tolerance 11.66) | 16 solutes in water | surface, CDS and empirical-shift terms differ; no solvent-specific fitting |
| ESP grid | PySCF density | far field 4.39e-3 Ha/e (0.6% tolerance on the max); relative 3.2% | 2 molecules, far-field points | sci-form uses Mulliken point charges; near field is not a density parity test |
| ECFP similarity | RDKit Morgan/ECFP4 | Spearman 0.756 | 780 pairs | bit vectors are not identical by design; ranking is the gate |

## Electronic structure

| Method | Compared with | Measured separation | Confidence region | Main limitation |
|---|---|---|---|---|
| EHT overlap | independent NumPy EHT | 7.53e-13 | 9 molecules, STO-3G basis | Wolfsberg-Helmholtz K = 1.75 is empirical |
| EHT Hamiltonian | independent NumPy EHT | 1.42e-11 eV | 9 molecules | no self-consistency, no correlation |
| EHT orbital energies | independent NumPy EHT | 1.25e-11 eV | 9 molecules | qualitative frontier orbitals; not quantitative spectroscopy |
| EHT metal parameters | YAeHMOP table | 0.0 max abs (tolerance 0.02) | 76 parameters | table transcription, not a physics validation |
| HF / SCF | PySCF 2.14.0 | total energy 1.22e-8 Ha; ERI 2.14e-11; AO→MO 1.83e-14 | 11 molecules, STO-3G, RHF/UHF/ROHF | minimal basis; no correlation; small molecule set |
| CIS / TDA | PySCF 2.14.0 | excitation energy 1.72e-6 eV; oscillator 0.0041 | 8 molecules | STO-3G, no solvent, no triplets |
| sTDA UV-Vis | PySCF TDA/TDHF | implementation parity 3.3e-6 eV; vs full TD-SCF 16.0 eV (tolerance 20), oscillator 0.93 | 7 molecules | monopole Coulomb is a method gap, not a bug; systematic blue shift |
| Kohn-Sham DFT | PySCF 2.14.0 | total energy 7.2e-5 Ha; gap 9.7e-3 eV | 6 small molecules, STO-3G, SVWN/PBE | target is 1e-5 Ha; `coulomb_energy` is zero; ECP is not fitted; no exact exchange |
| Machine-learned force field | independent NumPy Behler-Parrinello oracle | AEV exact; energy 5.1e-14; force 6.4e-6 | 8 molecules, two parameter sets, tiny networks | no trained production weights; no MD17/SPICE; the ANI engine has no external gate |

## Spectroscopy

| Method | Compared with | Measured separation | Confidence region | Main limitation |
|---|---|---|---|---|
| IR frequencies | xtb GFN2 `vibspectrum` | 0.0089 cm⁻¹ | 39 modes, same Hessian on both sides | harmonic approximation; no anharmonicity |
| IR intensities | PySCF analytic Hessian + SCF dipoles | 0.58 km/mol vs the Mulliken reference; 37.15 vs the exact-dipole reference | 2 molecules | intensity definition is reference-dependent; small set |
| NMR J couplings (Karplus/HLA) | reference HLA calculator | baseline 0.004 Hz; substituted 0.42 Hz | 7 baseline + 49 substituted cases | no ensemble averaging against experiment; solvent and temperature are not modeled |
| NMR chemical shifts | HOSE empirical table | not externally gated in this campaign | common ¹H/¹³C environments | NMRShiftDB2/SDBS data are not bundled because of licensing; rare environments fall back to defaults |
| UV-Vis | sTDA vs PySCF TDA/TDHF | 16.04 eV max abs (tolerance 20); oscillator 0.93 (tolerance 5) | 7 molecules | method gap; use for trends, not absolute λmax |
| EAS site ranking | textbook/DFT-confirmed orientation rules | 87.5% top-1 (7/8) | 8 substrates | qualitative ranking only |

## Materials and periodic systems

| Method | Compared with | Measured separation | Confidence region | Main limitation |
|---|---|---|---|---|
| Minimum-image distance | ASE 3.29.0 | 8.88e-16 Å | 4 cells, orthorhombic and triclinic | no long-range electrostatics |
| Space groups | spglib 2.8.0 | 0 translation mismatches / 230 | 230 cells | no magnetic space groups |
| CIF round-trip | pymatgen 2026.9.24 | 0 Å / degree | 4 structures | no disorder or partial occupancy |
| MOF-5 assembly | pormake Builder | 0 atom/equation difference | one pcu topology, same building blocks | one framework, not a MOF benchmark |
| Topology degrees | RCSR `.cgd` | 0 mismatches / 16 | 16 topologies | no nets outside the table |
| Porosity | pyzeo / Zeo++ | LCD 0.0037; PLD 0.2904; void fraction 0.00078 | MOF-5 | one framework; PLD tolerance is 0.5 Å |
| Powder XRD | pymatgen | 4.26e-14° in 2θ | 4 structures | peak positions only; intensities are a correlation check |
| Periodic stress | ASE Lennard-Jones + UnitCellFilter | 4.29e-7 eV/Å³ | 3 cells | numerical stress, not a general stress functional |
| Variable-cell relaxation | ASE BFGS | 3.78e-4 Å lattice parameter | 3 cells | one optimizer path |

## Dynamics, reaction paths and electrochemistry

| Method | Compared with | Measured separation | Confidence region | Main limitation |
|---|---|---|---|---|
| Growing String Method | ASE CI-NEB and analytic surfaces | Müller-Brown barrier 1.6e-3, TS 0.028; five curved wells ≤0.8% | Müller-Brown plus five curved wells | no pyGSM/Transition1x parity; backend accuracy limits the path |
| Kinetics / transition-state theory | Eyring expression and SciPy `solve_ivp` | rate 1.9e-16; equilibrium 8e-16 | small elementary networks | no Cantera mechanism I/O; the MBH adapter stays behind `beta-mbh` |
| Live molecular dynamics | RDKit UFF forces | forces 3.6e-13; NVE drift 4.4e-4; reversibility 5.7e-17 | UFF systems, short trajectories | no long-trajectory thermostat benchmark; force field limited |
| Electrical double layer | analytic Gouy-Chapman/Stern | worst 2.0e-14 | 3 temperatures × 3 ionic strengths × 6 potentials | 1:1 electrolyte, mean-field, no finite ion size beyond Stern |
| Constant Potential Method | LAMMPS `fix electrode/conp` | plate charge scale 1.76×; per-atom cosine 0.44 | graphene bilayer, 3×3 cells, V = ±0.5/1/2 V | point-charge hardness 5 eV vs ~22 eV Gaussian self-energy; over-polarizes dense conductors |
| Mobile Block Hessian | analytic diatomic + full Hessian | stretch 1.3e-5 SI; zero modes ≤5 cm⁻¹; 24/36 degrees of freedom | diatomic and benzene | experimental reference ±2 cm⁻¹ (Ghysels); global minima not integrated |

## Experimental modules: promotion gates

These are the modules that moved out of the alpha/beta tier. The table keeps
the promotion evidence in one place; the per-module pages keep the API and the
limitations.

| Module | Core path | Promotion evidence | Confidence region | Main limitation |
|---|---|---|---|---|
| Conformal geometric algebra | `cga` | motors 3.4e-15; dihedral 6.6e-15 | 16 rigid motions, seed 42 | refinement vs matrix method and MOF assembly not covered |
| Semidefinite relaxation | `sdr` | projections 4.4e-16; exact embeds 7e-15 | small distance sets | SDP scaling not benchmarked |
| Electrical double layer | `edl` | worst 2.0e-14 | 3T × 3I × 6 potentials | mean-field, 1:1 electrolyte |
| Kernel Polynomial Method | `kpm` | moments 3.5e-15; DOS 1.7e-14 | gapped spectra, up to 5000 sites | T=0 gapless Gibbs; stochastic noise |
| Constant Potential Method | `cpm` | LAMMPS bilayer comparison | graphene bilayer | hardness gap |
| Randomized eigensolver | `rand_nla` | rank-N/4 3.8e-4; N=600 9.6e-4 | near-identity/gapped EHT | full-spectrum fallback |
| Riemannian geometry | `riemannian` | projection 8.2e-12 | small PSD matrices | Frobenius proxy, first-order retraction |
| Growing String Method | `gsm` | barrier 1.6e-3 vs ASE | Müller-Brown + 5 wells | no pyGSM |
| Kinetics | `kinetics` | Eyring 1.9e-16 | small networks | no Cantera I/O |
| Periodic adapters | `periodic_linear` | phases 3.6e-16; bands 4.4e-16 | 1D chain, Monkhorst-Pack | complex H(k) not implemented |
| Live MD | `dynamics_live` | forces 3.6e-13 vs RDKit | UFF systems | short trajectories |
| Render bridge | `render_bridge` | JSON Schema + pack round-trip | payload contracts | schema, not visual correctness |

## What PASS does not mean

A passing row means the measured observable agrees with the named reference in
the stated region. It does **not** mean:

- the method is the same theory as the reference (sTDA is not TD-DFT);
- the element or molecule set is complete;
- the result is independent of the basis, functional or force-field parameters;
- the implementation is fast enough for production;
- a failed or missing gate would be hidden — the seven non-passing rows are
  listed below.

## Recorded gaps

The seven rows with `pass = false` are provenance or method-gap markers, not
silent failures. The full text is in the matrix.

| Row | What it marks | Why it is not a failure |
|---|---|---|
| F44 provenance | GFN1/GFN2 fixture provenance | The measured rows F44/F45 are green; this row records where the fixture came from |
| F46 provenance | D3/D4 fixture provenance | The measured rows are green; standalone D4 remains a documented method gap |
| F71 bounds | per-class bounds tolerance | Tolerances are read per topological class; class 4 has 4 divergent pairs out of 388 |
| F82 descriptors | descriptor-by-descriptor comparison | Tolerances live with each measured row; not every descriptor has the same definition |
| F83 canonical SMILES | string equality across toolkits | Canonical strings are intentionally not bit-identical; isomorphism and round-trip are the gates |
| F84 ALPB | total vs tblite ALPB | The electrostatic term is gated at 3.14 kcal/mol; surface terms are a different model |
| F85 stress | stress convention | The convention matches ASE; the relaxed cell is a separate optimizer path |

## Methods outside this campaign

Some methods have unit or integration tests but are not part of the 129-row
external campaign. They are listed here so the absence of a row is not read as
a passing gate.

| Method | What is covered | What is missing |
|---|---|---|
| HF-3c full energy | the D3(BJ) component is gated against dftd4 (6.8e-9 Ha); the SCF and integral layers share the HF/PySCF gates | no independent full-HF-3c energy oracle in this campaign; the historical "0.01%" claim is not a substitute |
| GIAO NMR shieldings | unit and integration tests of the RHF/GIAO path | no ORCA or experimental shielding parity in this campaign |
| Fluorescence | implementation tests | no external emission benchmark |
| ANI production weights | architecture and finite-difference force tests | no official ANI-2x weight comparison; no MD17/SPICE benchmark |
| MMFF94 parameter coverage | full-energy parity on 7 molecules | no broad parameter-table sweep across all MMFF94 atom types |
| Large-system scaling | KPM at 5000 sites; RandNLA at N=600 | no 10k+ atom or periodic production benchmark in this campaign |

## Reproduce

```bash
# The pinned reference environment
cat validation/refs.lock.json

# The individual oracle scripts
python3 scripts/validation/oracle_rdkit.py
python3 scripts/validation/oracle_pyscf.py
python3 scripts/validation/oracle_tblite.py
python3 scripts/validation/oracle_kpm.py
python3 scripts/validation/oracle_randnla.py
python3 scripts/validation/oracle_gsm.py

# Collect the rows into the machine-readable summary
python3 scripts/validation/collect.py

# Re-render the full matrix
python3 scripts/docs/generate_validation_matrix.py

# Re-render the animated diagrams (theme-aware, transparent)
python3 scripts/docs/generate_animated_svgs.py
```

The Rust gates named in the matrix are run by the regular test targets, for
example:

```bash
cargo test --release --test cross_validation
cargo test --release --test alpha_gsm_validation
cargo test --release --test beta_randnla_validation
cargo test --release --test alpha_edl_validation
```
