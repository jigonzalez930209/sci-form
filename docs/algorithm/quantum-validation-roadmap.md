# Quantum Validation Roadmap

This note collects the current validation surface for PM3, GFN0/xTB, GFN1, GFN2, and HF-3c, and defines a practical reference strategy for turning the current implementations into a defensible benchmarked stack.

## Current state in the repo

- PM3 lives in [src/pm3/solver.rs](../../src/pm3/solver.rs) and already exposes SCF results, gradients, and heat-of-formation estimates.
- GFN0/xTB lives in [src/xtb/solver.rs](../../src/xtb/solver.rs) and is the current tight-binding base.
- GFN1 and GFN2 are currently layered on top of GFN0 in [src/xtb/gfn1.rs](../../src/xtb/gfn1.rs) and [src/xtb/gfn2.rs](../../src/xtb/gfn2.rs); they are not yet independent, literature-faithful parameterizations.
- HF-3c lives in [src/hf/api.rs](../../src/hf/api.rs) and is already integrated into the public API.

## What the repo already validates

- [tests/regression/test_experimental_comparison.rs](../../tests/regression/test_experimental_comparison.rs) already contains NIST/CCCBDB reference energies for H2, H2O, CH4, NH3, HF, CO, C2H4, and CH2O.
- The same test file also stores NIST ionization potentials for H2, H2O, CH4, NH3, and HF.
- [tests/fixtures/experimental_reference.json](../../tests/fixtures/experimental_reference.json) contains experimental bond lengths, dipoles, ionization energies, IR frequencies, and NMR shifts for a compact set of light molecules.
- [tests/fixtures/gdb20_reference.json.gz](../../tests/fixtures/gdb20_reference.json.gz) is the current geometry reference source for conformer validation.

## Observed data snapshot

The following values were recovered from the existing regression suite on 2026-03-24, without adding new benchmark code.

### Current gap table from the regression suite

Source: existing comprehensive table in [tests/regression/test_experimental_comparison.rs](../../tests/regression/test_experimental_comparison.rs).

| Molecule | PM3 gap (eV) | xTB gap (eV, current `solve_xtb` = GFN0 SCC) | STO-3G HF reference energy (Ha) |
|----------|--------------|----------------------------------------------|----------------------------------|
| H2 | 44.054 | 18.980 | -1.1175 |
| H2O | 1.839 | 25.579 | -74.9659 |
| CH4 | 10.814 | 26.026 | -39.7269 |
| NH3 | 20.383 | 5.288 | -55.4544 |
| HF | 3.318 | 24.285 | -98.5708 |
| CO | 12.273 | 7.758 | -111.2255 |
| C2H4 | 7.175 | 4.972 | -77.0739 |
| CH2O | 10.928 | 25.976 | -112.3522 |

What this table tells us:

- The repo already exposes a reproducible compact comparison surface for PM3 and the current xTB path.
- The current `xTB` regression output is GFN0-style SCC from `solve_xtb`, not GFN1 or GFN2.
- These values are useful for trend monitoring, but they are not directly comparable to STO-3G HF total energies because they are different observables.

### Current HF-3c energy comparison snapshot

Source: existing focused tests in [tests/regression/test_experimental_comparison.rs](../../tests/regression/test_experimental_comparison.rs).

| Molecule | HF-3c pure HF energy (Ha) | Experimental SCF energy reported by current test (Ha) | Difference (mHa) |
|----------|----------------------------|--------------------------------------------------------|------------------|
| H2 | -1.116684 | -4.148070 | 3031.39 |
| H2O | -74.966720 | -267.610738 | 192644.02 |

This is an important warning signal:

- The current focused HF-3c vs experimental SCF printout is not yet a valid acceptance gate.
- Either the compared quantities are not normalized to the same definition, or the current experimental SCF path is reporting a value that is not directly comparable to `hf_energy`.
- Until that mismatch is resolved, CCCBDB/PySCF-style references remain the safer validation anchor for HF-3c.

## Validation rule of thumb

Do not force a single 1% rule onto every observable.

- 1% is reasonable only for same-method or same-basis scalar quantities, such as STO-3G HF energies against a CCCBDB/PySCF reference, or a fully specified benchmark dataset from the same model.
- For semi-empirical and tight-binding methods, geometry and relative ordering are usually more meaningful than raw experimental agreement.
- For charges, gaps, and vibrational spectra, use property-specific tolerances rather than a single percentage threshold.

## Recommended reference targets

### HF-3c

- Reference source: CCCBDB / PySCF-style STO-3G HF values.
- Primary observables: total energy, orbital energies, bond lengths, dipole moment, CIS excitation order.
- Current data note: the repo already has direct HF-3c comparison prints, but the present `HF-3c vs experimental SCF` numbers are inconsistent enough that they should be treated as diagnostics, not release criteria.
- Suggested thresholds:
  - total energy: <= 0.5% relative error for small closed-shell systems
  - bond lengths: <= 0.01-0.02 Å absolute error
  - dipole moments: <= 5% relative error or <= 0.1 D absolute error
  - orbital ordering: exact HOMO/LUMO ordering, not percentage-based

### PM3

- Reference source: original Stewart PM3 parameterization targets, plus a small public molecule set with published heats of formation and geometries.
- Primary observables: heat of formation, optimized geometry, Mulliken charges, HOMO/LUMO gap, convergence behavior.
- Current data note: the compact regression suite already yields PM3 gaps for H2, H2O, CH4, NH3, HF, CO, C2H4, and CH2O, so PM3 trend regression can start immediately on those molecules.
- Suggested thresholds:
  - heat of formation: <= 1% only when the reference is a like-for-like PM3 value; otherwise use kcal/mol absolute error
  - geometry: <= 0.05 Å RMSD for light molecules; looser for flexible species
  - Mulliken charges: <= 0.05 e per atom or <= 5%
  - gap: compare in eV with an absolute tolerance, not a percent threshold

### GFN0 / xTB

- Reference source: official xTB/CREST outputs or published benchmark sets, not experiment alone.
- Primary observables: conformer energetics, charges, dipoles, gradients, and structural preferences.
- Current data note: the existing regression table is already producing GFN0-style gaps through `solve_xtb`; this is enough to establish a first compact regression baseline before introducing external benchmark corpora.
- Suggested benchmark families: S66/S66x8, WATER27, HAL59, ICONF, MOR41, and a small metal-organic subset already supported by the parameter table.
- Suggested thresholds:
  - conformer energy ranking: correct ordering on benchmark subsets
  - geometries: RMSD or pairwise-distance RMSD against reference conformers
  - charges and dipoles: absolute tolerances, with per-atom and per-molecule checks

### GFN1 / GFN2

- Reference source: official GFN1-xTB and GFN2-xTB implementations, plus literature benchmark sets.
- Primary observables: shell-resolved charges, multipole response, noncovalent interaction energies, halogen-bond trends, and geometry optimization quality.
- Suggested thresholds:
  - same-method energy comparison against reference implementation: <= 1%
  - geometry and conformer ranking: benchmark-specific absolute thresholds
  - charges/dipoles/multipoles: component-wise absolute thresholds

## Roadmap to build a defensible GFN1/GFN2 base

### Phase 1: Freeze the reference corpus

- Create a compact benchmark suite of representative molecules:
  - small organics: H2, H2O, CH4, NH3, HF, CO, C2H4, CH2O, ethanol, acetic acid
  - noncovalent interactions: dimers from S66/S66x8 and water clusters
  - halogen-bond cases: chlorinated / brominated / iodinated probes
  - transition-metal cases: the supported metal subset already present in the xTB parameter table
- Store the reference observable for each molecule: energy, geometry, charge distribution, dipole, and gradient norm when available.
- Immediate action: promote the current compact gap table above into a checked-in baseline artifact so the project stops depending on console output as the only source of truth.

### Phase 2: Separate the models cleanly

- Keep GFN0 as the minimal base.
- Move GFN1 to its own parameter table and shell-resolved self-consistency path.
- Move GFN2 to its own parameter table, multipole response, and halogen-bond correction path.
- Avoid reusing GFN0 defaults as a hidden proxy for missing GFN1/GFN2 parameters.

### Phase 3: Benchmark each observable independently

- Energies: compare against reference implementation outputs and benchmark datasets.
- Geometries: compare against reference conformers using pairwise-distance RMSD and heavy-atom RMSD.
- Charges and dipoles: compare component-wise against reference outputs.
- Gradients: validate finite-difference consistency and force direction.
- Convergence: enforce deterministic iteration counts on the small reference set where possible.

### Phase 4: Add CI gates

- Add a fast smoke suite for the compact corpus.
- Add a slower regression suite for the benchmark families.
- Fail the build when a method regresses on energy ranking, convergence, or geometry quality.
- Report thresholds per observable instead of one global percentage.
- Add a specific guard that rejects comparisons when units or energy definitions differ, so Hartree-vs-eV or corrected-vs-uncorrected values cannot silently pass through the benchmark pipeline.

### Phase 5: Document supported accuracy

- Publish a per-method accuracy table in the docs.
- Mark which observables are validated against experiment, which ones are validated against a reference implementation, and which ones are only invariant-checked.

## Practical conclusion

The current repo is already strong on geometry validation and on NIST/CCCBDB comparisons for the HF/STO-3G-like layer. PM3, xTB, GFN1, GFN2, and HF-3c should not all be judged with the same 1% experimental rule. The correct approach is a hybrid validation plan:

- use experiment where it is physically meaningful and well tabulated,
- use same-method references where the method is semi-empirical,
- and use absolute or structural tolerances for charges, geometries, gaps, and spectra.

The current measured data already supports three concrete decisions:

- PM3 and the current GFN0/xTB path can start using the compact eight-molecule gap panel as an immediate regression baseline.
- HF-3c should stay anchored to CCCBDB/PySCF-like references until the current `HF-3c vs experimental SCF` mismatch is resolved.
- GFN1 and GFN2 still need an external same-method reference corpus before any 1% target becomes technically defensible.