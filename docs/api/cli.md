# CLI Reference

The `sci-form` binary exposes **90 subcommands**: 56 are always available and 34 are behind Cargo feature flags. This page is the exhaustive reference; the [CLI guide](../guide/cli.md) is the task-oriented walkthrough.

## Installation

### Pre-built Binaries

Download from [GitHub Releases](https://github.com/jigonzalez930209/sci-form/releases):
- `sci-form-x86_64-unknown-linux-gnu` — Linux x86_64
- `sci-form-x86_64-apple-darwin` — macOS Intel
- `sci-form-aarch64-apple-darwin` — macOS Apple Silicon
- `sci-form-x86_64-pc-windows-msvc.exe` — Windows x86_64
- `sci-form-aarch64-unknown-linux-gnu` — Linux ARM64

### Build from Source

```bash
cargo install sci-form-cli
```

### Feature-gated subcommands

The optional subcommands are compiled in only when the corresponding Cargo feature is enabled on the `sci-form-cli` package:

```bash
cargo install sci-form-cli --features alpha-gsm,alpha-edl,alpha-kinetics
```

| Feature | Unlocks |
|---------|---------|
| `experimental-alpb` | `alpb` |
| `experimental-d4` | `d4` |
| `experimental-cpm` | `cpm` |
| `experimental-eeq` | (none — `eeq` is always built) |
| `alpha-gsm` | `gsm-backend-plan`, `gsm-compare-backends` |
| `alpha-edl` | `edl-profile`, `edl-capacitance-scan`, `edl-helmholtz`, `edl-gouy-chapman`, `edl-gcs` |
| `alpha-kinetics` | `htst-rate`, `extract-kinetics-diagnostics`, `solve-microkinetic-network`, `solve-microkinetic-steady-state`, `analyze-gsm-mbh-htst-step` |
| `alpha-periodic-linear` | `kmesh`, `compute-periodic-dos`, `solve-periodic-randnla`, `bloch-phase`, `build-bloch-hamiltonian`, `assemble-periodic-operators`, `validate-electron-count` |
| `alpha-render-bridge` | `edl-chart`, `capacitance-chart`, `arrhenius-chart`, `thermal-prop-chart`, `phase-portrait-chart`, `band-structure-chart`, `dos-chart`, `kpoint-path-chart`, `fermi-surface-chart`, `trajectory-chart`, `reaction-coordinate-chart` |
| `alpha-reaction-dynamics` | `reaction-dynamics3d` |

The `alpha-render-bridge` chart commands need a second feature as well: `edl-chart` and `capacitance-chart` also require `alpha-edl`, `arrhenius-chart`, `thermal-prop-chart` and `phase-portrait-chart` require `alpha-kinetics`, `band-structure-chart`, `dos-chart`, `kpoint-path-chart` and `fermi-surface-chart` require `alpha-periodic-linear`, and `trajectory-chart` / `reaction-coordinate-chart` require `alpha-gsm`.

### Negative numeric option values

Several options take negative numbers (`--e-min`, `--e-max`, `--mu-ev`, `--v-min`, `--v-max`, `--reaction-free-energy`, …). The argument parser treats a bare `-4.44` as a flag, so negative values must be attached with `=`:

```bash
sci-form cpm "[6,8,1,1,1,1,1,1,1]" "[...]" --mu-ev=-4.44     # correct
sci-form cpm "[6,8,1,1,1,1,1,1,1]" "[...]" --mu-ev -4.44      # usage error
```

---

## Conformer generation and structure

### `embed`

Generate a 3D conformer for a single SMILES.

```
sci-form embed <SMILES> [OPTIONS]
```

**Options:**

| Option | Default | Description |
|--------|---------|-------------|
| `-s, --seed <N>` | `42` | RNG seed |
| `-f, --format <FMT>` | `json` | Output format: `json`, `xyz`, `sdf` |

**Examples:**

```bash
sci-form embed "CCO"                # JSON
sci-form embed "c1ccccc1" -f xyz    # XYZ
sci-form embed "CC(=O)O" -f sdf -s 123
```

**JSON output:**

```json
{
  "smiles": "CCO",
  "num_atoms": 9,
  "coords": [0.123, 0.456, 0.789, "..."],
  "elements": [6, 6, 8, 1, 1, 1, 1, 1, 1],
  "bonds": [[0, 1, "SINGLE"], [1, 2, "SINGLE"]],
  "error": null,
  "time_ms": 2.34
}
```

---

### `batch`

Process multiple SMILES from a file or stdin. The input path is an **option**, not a positional argument.

```
sci-form batch [OPTIONS]
```

**Options:**

| Option | Default | Description |
|--------|---------|-------------|
| `-i, --input <FILE>` | stdin | Input file (omit to read stdin) |
| `-o, --output <FILE>` | stdout | Output file (omit to write stdout) |
| `-s, --seed <N>` | `42` | RNG seed |
| `-t, --threads <N>` | `0` | Number of threads (`0` = auto) |
| `-f, --format <FMT>` | `json` | Output format: `json`, `xyz`, `sdf` |

**Input file format:** One SMILES per line (name optional):
```
CCO ethanol
c1ccccc1 benzene
CC(=O)O acetic_acid
```

**Examples:**

```bash
sci-form batch -i molecules.smi
sci-form batch -i molecules.smi -f xyz -t 8
cat molecules.smi | sci-form batch
sci-form batch -i molecules.smi > conformers.json
```

---

### `parse`

Parse SMILES and print the molecular graph. This command writes **plain text**, not JSON.

```
sci-form parse <SMILES>
```

```bash
$ sci-form parse "c1ccccc1"
Atoms: 12, Bonds: 12
    0: elem= 6 hyb=SP2 charge=0 h=0
    1: elem= 6 hyb=SP2 charge=0 h=0
    ...
```

---

### `info`

Print the build information: version, feature summary, supported output formats and the available bindings.

```
sci-form info
```

```bash
$ sci-form info
sci-form 0.15.2
Features: ETKDG distance geometry, CSD torsion patterns, SMARTS matching
  EHT, PM3, xTB, GFN1, GFN2, HF-3c, ANI, Gasteiger charges, SASA, ESP,
   DOS, MMFF94, UFF, stereochemistry, solvation, SSSR, ECFP, Tanimoto, cell, MOF assembly
Formats: JSON, XYZ, SDF
Bindings: CLI, Python (PyO3), TypeScript (WASM)
```

---

### `embed-diverse`

Generate a diverse conformer ensemble and cluster it by Butina RMSD.

```
sci-form embed-diverse [OPTIONS] <SMILES>
```

| Option | Default | Description |
|--------|---------|-------------|
| `--n-conformers <N>` | `10` | Number of embedding attempts (different seeds) |
| `--rmsd-cutoff <F>` | `1.0` | RMSD clustering cutoff in Å |
| `--base-seed <N>` | `42` | Base RNG seed (seeds are `base_seed..base_seed+n_conformers`) |

```bash
sci-form embed-diverse "CCO" --n-conformers 20 --rmsd-cutoff 0.8
```

---

### `filter-diverse`

Keep only the diverse cluster centroids of a conformer ensemble.

```
sci-form filter-diverse [OPTIONS] <CONFORMERS>
```

| Option | Default | Description |
|--------|---------|-------------|
| `--rmsd-cutoff <F>` | `1.0` | RMSD clustering cutoff in Å |

`CONFORMERS` is a JSON array of conformers, each one a flat xyz coordinate array. The command prints the JSON array of the retained indices.

```bash
sci-form filter-diverse "[[0,0,0, 1,0,0], [0,0,0, 1.05,0,0], [0,0,0, 5,0,0]]"
```

---

### `canonical-smiles`

Print the canonical SMILES representation of a molecule.

```
sci-form canonical-smiles <SMILES>
```

```bash
$ sci-form canonical-smiles "C(F)(Cl)(Br)I"
BrC(Cl)(F)I
```

---

## Charges, surfaces and stereochemistry

### `charges`

Compute Gasteiger-Marsili partial charges.

```
sci-form charges <SMILES>
```

```bash
$ sci-form charges "CCO"
{
  "charges": [-0.387, -0.042, -0.228, 0.123],
  "iterations": 6,
  "total_charge": 0.0,
  "converged": true
}
```

---

### `charges-configured`

Gasteiger-Marsili charges with an explicit convergence configuration.

```
sci-form charges-configured [OPTIONS] <SMILES>
```

| Option | Default | Description |
|--------|---------|-------------|
| `--max-iter <N>` | `6` | Electronegativity-equalization iterations |
| `--initial-damping <F>` | `0.5` | Initial damping factor (halved every iteration) |
| `--convergence-threshold <F>` | `1e-10` | Convergence threshold on the maximum charge transfer |

```bash
sci-form charges-configured "CCO" --max-iter 12 --initial-damping 0.7
```

---

### `sasa`

Compute solvent-accessible surface area.

```
sci-form sasa [OPTIONS] <ELEMENTS_JSON> <COORDS_JSON>
```

| Option | Default | Description |
|--------|---------|-------------|
| `-p, --probe-radius <F>` | `1.4` | Probe radius in Å |

---

### `solvation`

Compute non-polar (SASA-based) solvation energy; with charges it also returns the Generalized Born electrostatic term.

```
sci-form solvation [OPTIONS] <ELEMENTS_JSON> <COORDS_JSON>
```

| Option | Default | Description |
|--------|---------|-------------|
| `-c, --charges <JSON>` | *(empty)* | Per-atom charges as JSON array. When omitted only the non-polar term is computed |
| `-p, --probe-radius <F>` | `1.4` | Probe radius in Å |

The GB term uses a solvent dielectric of 78.5 and an internal dielectric of 1.0; these are fixed in the CLI.

```bash
# Using the SMILES pipeline
RESULT=$(sci-form embed "CCO")
ELEM=$(echo "$RESULT" | jq '.elements')
COORDS=$(echo "$RESULT" | jq '.coords')
sci-form solvation "$ELEM" "$COORDS"

# With GB electrostatics
sci-form solvation "$ELEM" "$COORDS" --charges "[-0.39,-0.04,-0.23,0.13,0.13,0.13,0.13,0.13,0.13]"
```

Output (non-polar only):

```json
{
  "energy_kcal_mol": -0.72,
  "atom_contributions": [-0.31, 0.02, 0.02],
  "atom_sasa": [53.68, 33.80, 33.71],
  "total_sasa": 121.20
}
```

---

### `eeq`

Compute EEQ (electronegativity-equalization) geometry-dependent charges.

```
sci-form eeq [OPTIONS] <ELEMENTS_JSON> <COORDS_JSON>
```

| Option | Default | Description |
|--------|---------|-------------|
| `-t, --total-charge <F>` | `0.0` | Total molecular charge |

---

### `stereo`

Analyze stereochemistry: R/S stereocenters and E/Z double bonds.

```
sci-form stereo [OPTIONS] <SMILES>
```

| Option | Default | Description |
|--------|---------|-------------|
| `-c, --coords <JSON>` | `[]` | Flat coordinate array as JSON string (enables 3D R/S and E/Z assignment) |

```bash
sci-form stereo "C(F)(Cl)(Br)I"
sci-form stereo "C/C=C/C" --coords "[$(sci-form embed 'C/C=C/C' | jq '.coords')]"
```

Output (the `configuration` field is only filled in when 3D coordinates are supplied):

```json
{
  "stereocenters": [
    {
      "atom_index": 0,
      "element": 6,
      "substituent_indices": [4, 3, 2, 1],
      "priorities": [1, 2, 3, 4],
      "configuration": "R"
    }
  ],
  "double_bonds": [],
  "n_stereocenters": 1,
  "n_double_bonds": 0,
  "atropisomeric_axes": [],
  "prochiral_centers": [],
  "helical_chirality": []
}
```

---

### `sssr`

Perceive the Smallest Set of Smallest Rings.

```
sci-form sssr <SMILES>
```

```bash
sci-form sssr "c1ccc2ccccc2c1"  # naphthalene
```

Output:

```json
{
  "rings": [
    { "atoms": [0, 9, 8, 3, 2, 1], "size": 6, "is_aromatic": true },
    { "atoms": [3, 8, 7, 6, 5, 4], "size": 6, "is_aromatic": true }
  ],
  "ring_size_histogram": [0, 0, 0, 0, 0, 0, 2]
}
```

---

### `ecfp`

Compute an ECFP (Morgan) fingerprint.

```
sci-form ecfp [OPTIONS] <SMILES>
```

| Option | Default | Description |
|--------|---------|-------------|
| `-r, --radius <N>` | `2` | Morgan radius (ECFP4 = 2, ECFP6 = 3) |
| `-n, --n-bits <N>` | `2048` | Bit-vector size |

```bash
sci-form ecfp "c1ccccc1"
sci-form ecfp "c1ccccc1" --radius 3 --n-bits 1024
```

---

### `ecfp-batch`

Compute ECFP fingerprints for a batch of SMILES given as a JSON array.

```
sci-form ecfp-batch [OPTIONS] <SMILES_LIST>
```

| Option | Default | Description |
|--------|---------|-------------|
| `-r, --radius <N>` | `2` | Morgan radius |
| `-n, --n-bits <N>` | `2048` | Bit-vector size |

```bash
sci-form ecfp-batch '["CCO", "c1ccccc1"]'
```

---

### `tanimoto`

Compute the Tanimoto similarity between two molecules.

```
sci-form tanimoto [OPTIONS] <SMILES1> <SMILES2>
```

| Option | Default | Description |
|--------|---------|-------------|
| `-r, --radius <N>` | `2` | ECFP radius |
| `-n, --n-bits <N>` | `2048` | Bit-vector size |

```bash
$ sci-form tanimoto "c1ccccc1" "Cc1ccccc1"
{"tanimoto": 0.333333}
```

---

### `rmsd`

Compute the RMSD between two coordinate sets (Kabsch alignment).

```
sci-form rmsd <COORDS_JSON> <REFERENCE_JSON>
```

Output:

```json
{ "rmsd": 0.034, "rotation": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]] }
```

---

## Electronic structure

All of these take the same two positional arguments: `ELEMENTS_JSON` (JSON array of atomic numbers) and `COORDS_JSON` (flat `[x, y, z, ...]` array in Å).

### `eht`

Run Extended Hückel Theory.

```
sci-form eht [OPTIONS] <ELEMENTS_JSON> <COORDS_JSON>
```

| Option | Default | Description |
|--------|---------|-------------|
| `-k, --k <F>` | `1.75` | Wolfsberg-Helmholtz constant (`0` selects the default) |

```bash
sci-form eht "[6,1,1,1,1]" "[0,0,0,1,0,0,0,1,0,0,0,1]" --k 1.75
```

---

### `eht-gradient`

Compute EHT energy-surface gradients by central finite differences.

```
sci-form eht-gradient <ELEMENTS_JSON> <COORDS_JSON>
```

---

### `optimize-geometry-eht`

Optimize a geometry on the EHT energy surface with the shared BFGS engine.

```
sci-form optimize-geometry-eht [OPTIONS] <ELEMENTS_JSON> <COORDS_JSON>
```

| Option | Default | Description |
|--------|---------|-------------|
| `--max-steps <N>` | `200` | Maximum number of optimization steps |
| `--grad-threshold <F>` | `0.01` | RMS gradient convergence threshold in eV/Å |
| `--energy-threshold <F>` | `1e-6` | Energy convergence threshold in eV |
| `--step-size <F>` | `0.05` | Initial steepest-descent step size in Å |

---

### `pm3`, `xtb`, `gfn1`, `gfn2`, `hf3c`

Semi-empirical and tight-binding single points.

```
sci-form pm3  <ELEMENTS_JSON> <COORDS_JSON>   # PM3
sci-form xtb  <ELEMENTS_JSON> <COORDS_JSON>   # GFN0-xTB
sci-form gfn1 <ELEMENTS_JSON> <COORDS_JSON>   # GFN1-xTB
sci-form gfn2 <ELEMENTS_JSON> <COORDS_JSON>   # GFN2-xTB
sci-form hf3c <ELEMENTS_JSON> <COORDS_JSON>   # HF-3c
```

---

### `cisd`

Compute CISD excited states from an internal STO-3G RHF calculation.

```
sci-form cisd [OPTIONS] <ELEMENTS_JSON> <COORDS_JSON>
```

| Option | Default | Description |
|--------|---------|-------------|
| `--n-states <N>` | `5` | Number of excited states to report |

---

### `population`

Mulliken & Löwdin population analysis (requires the EHT engine).

```
sci-form population <ELEMENTS_JSON> <COORDS_JSON>
```

---

### `dipole`

Compute the molecular dipole moment (Debye).

```
sci-form dipole <ELEMENTS_JSON> <COORDS_JSON>
```

---

### `dos`

Compute the density of states from the EHT spectrum.

```
sci-form dos [OPTIONS] <ELEMENTS_JSON> <COORDS_JSON>
```

| Option | Default | Description |
|--------|---------|-------------|
| `-s, --sigma <F>` | `0.3` | Broadening width in eV |
| `--e-min <F>` | `-30.0` | Energy window minimum in eV |
| `--e-max <F>` | `5.0` | Energy window maximum in eV |
| `-n, --n-points <N>` | `500` | Grid points |

Output: `{"sigma": ..., "data": [[energy, dos], ...]}`.

---

### `dos-multimethod`

Same as `dos`, but the orbital energies come from the chosen electronic-structure method.

```
sci-form dos-multimethod [OPTIONS] <ELEMENTS_JSON> <COORDS_JSON>
```

| Option | Default | Description |
|--------|---------|-------------|
| `--method <M>` | `eht` | Method: `eht`, `pm3`, `xtb`, `gfn1`, `gfn2` or `hf3c` |
| `-s, --sigma <F>` | `0.3` | Broadening width in eV |
| `--e-min <F>` | `-30.0` | Energy window minimum in eV |
| `--e-max <F>` | `5.0` | Energy window maximum in eV |
| `-n, --n-points <N>` | `500` | Grid points |

---

### `esp`

Compute the electrostatic potential on a 3D grid.

```
sci-form esp [OPTIONS] <ELEMENTS_JSON> <COORDS_JSON>
```

| Option | Default | Description |
|--------|---------|-------------|
| `-s, --spacing <F>` | `0.5` | Grid spacing in Å |
| `-p, --padding <F>` | `3.0` | Padding beyond the molecule in Å |

Output: `{"origin": [...], "spacing": ..., "dims": [...], "n_values": ..., "values": [...]}`.

---

### `ani`

Run an ANI ML-potential calculation.

```
sci-form ani <ELEMENTS_JSON> <COORDS_JSON>
```

The CLI builds a small internally generated network whose input dimension matches the real AEV length (4032 with the default parameters), so the reported energies are **not** physically meaningful ANI-2x energies — see [Machine Learning Potentials](../algorithm/machine-learning-potentials.md).

---

### `aevs-tm`

Compute ANI-TM atomic environment vectors (extended element set).

```
sci-form aevs-tm <ELEMENTS_JSON> <COORDS_JSON>
```

---

## Force fields and geometry optimization

### `uff`

Compute the UFF force-field energy.

```
sci-form uff <SMILES> <COORDS_JSON>
```

Output: `{"energy": 12.3, "unit": "kcal/mol"}`.

---

### `optimize-framework`

Optimize a geometry with the shared BFGS engine (UFF backend).

```
sci-form optimize-framework [OPTIONS] <SMILES> <COORDS_JSON>
```

| Option | Default | Description |
|--------|---------|-------------|
| `--method <M>` | `bfgs` | Optimization method: `bfgs` or `steepest_descent` |
| `--max-iter <N>` | `200` | Maximum number of iterations |
| `--force-tol <F>` | `0.05` | Force convergence threshold in eV/Å |
| `--energy-tol <F>` | `1e-6` | Energy convergence threshold in eV |
| `--max-step <F>` | `0.2` | Maximum step size in Å |
| `--fixed-atoms <JSON>` | *(empty)* | JSON array of fixed atom indices |
| `--a`, `--b`, `--c` | — | Lattice lengths in Å (all three enable periodic optimization) |
| `--alpha`, `--beta`, `--gamma` | `90.0` | Lattice angles in degrees |

---

## Spectroscopy

### `ir-spectrum-broadened`

Compute a broadened IR spectrum; the vibrational analysis is run internally.

```
sci-form ir-spectrum-broadened [OPTIONS] <ELEMENTS_JSON> <COORDS_JSON>
```

| Option | Default | Description |
|--------|---------|-------------|
| `--method <M>` | `eht` | Hessian method: `eht`, `pm3`, `xtb`, `uff`, ... |
| `--gamma <F>` | `15.0` | Broadening width in cm⁻¹ |
| `--wn-min <F>` | `400.0` | Minimum wavenumber in cm⁻¹ |
| `--wn-max <F>` | `4000.0` | Maximum wavenumber in cm⁻¹ |
| `--n-points <N>` | `1000` | Number of grid points |
| `--broadening <T>` | `lorentzian` | Broadening type: `gaussian` or `lorentzian` |
| `--step-size <F>` | `0.0` | Finite-difference step in Å (`0` = engine default) |
| `--smiles <S>` | *(empty)* | SMILES string (required for `method=uff`) |

---

### `ensemble-j-couplings`

Predict J-couplings averaged over a conformer ensemble with Boltzmann weighting.

```
sci-form ensemble-j-couplings [OPTIONS] <SMILES> <CONFORMER_COORDS> <ENERGIES_KCAL>
```

| Option | Default | Description |
|--------|---------|-------------|
| `--temperature-k <F>` | `298.15` | Temperature in K |

`CONFORMER_COORDS` is a JSON array of conformers (each a flat xyz array) and `ENERGIES_KCAL` a JSON array of the matching energies in kcal/mol.

---

### `assign-peaks`

Assign IR peaks to functional groups using the group-frequency tables.

```
sci-form assign-peaks [OPTIONS] <FREQUENCIES> <INTENSITIES> <ELEMENTS>
```

| Option | Default | Description |
|--------|---------|-------------|
| `--mode-displacements <JSON>` | *(empty)* | JSON array of modes, each an array of `[x, y, z]` triples |

`FREQUENCIES` are vibrational wavenumbers in cm⁻¹ and `INTENSITIES` the matching IR intensities.

---

## 3D descriptors and machine learning

### `whim`

Compute WHIM 3D molecular descriptors.

```
sci-form whim [OPTIONS] <ELEMENTS_JSON> <COORDS_JSON>
```

| Option | Default | Description |
|--------|---------|-------------|
| `--weighting <W>` | `mass` | Weighting: `unit`, `mass`, `volume`, `electronegativity` or `polarizability` |

---

### `rdf`

Compute radial distribution function (RDF) descriptors.

```
sci-form rdf [OPTIONS] <ELEMENTS_JSON> <COORDS_JSON>
```

| Option | Default | Description |
|--------|---------|-------------|
| `--charges <JSON>` | *(empty)* | Per-atom charges; when omitted Gasteiger charges are derived internally |
| `--beta <F>` | `100.0` | Gaussian smoothing factor in Å⁻² |
| `--r-max <F>` | `12.0` | Maximum radius in Å |
| `--n-points <N>` | `128` | Number of sampling points |

---

### `getaway`

Compute GETAWAY 3D molecular descriptors.

```
sci-form getaway [OPTIONS] <ELEMENTS_JSON> <COORDS_JSON>
```

| Option | Default | Description |
|--------|---------|-------------|
| `--bonds <JSON>` | *(empty)* | JSON array of `[i, j]` bond pairs; when omitted bonds are perceived from the geometry |
| `--max-lag <N>` | `8` | Maximum topological distance for the autocorrelations |

---

### `train-random-forest`

Train a Random Forest regressor on a feature matrix.

```
sci-form train-random-forest [OPTIONS] <FEATURES> <TARGETS>
```

| Option | Default | Description |
|--------|---------|-------------|
| `--n-trees <N>` | `100` | Number of trees |
| `--max-depth <N>` | `10` | Maximum tree depth |
| `--min-samples-split <N>` | `5` | Minimum samples required to split a node |
| `--min-samples-leaf <N>` | `2` | Minimum samples required for a leaf |
| `--sample-fraction <F>` | `1.0` | Bootstrap sample fraction (0.0–1.0) |
| `--seed <N>` | `42` | RNG seed |

`FEATURES` is a JSON array of feature vectors and `TARGETS` a JSON array of target values.

---

### `train-gradient-boosting`

Train a Gradient Boosting regressor on a feature matrix.

```
sci-form train-gradient-boosting [OPTIONS] <FEATURES> <TARGETS>
```

| Option | Default | Description |
|--------|---------|-------------|
| `--n-estimators <N>` | `100` | Number of boosting stages |
| `--learning-rate <F>` | `0.1` | Learning rate |
| `--max-depth <N>` | `5` | Maximum tree depth |
| `--min-samples-split <N>` | `5` | Minimum samples required to split a node |
| `--min-samples-leaf <N>` | `2` | Minimum samples required for a leaf |
| `--subsample <F>` | `0.8` | Row subsample fraction (0.0–1.0) |
| `--seed <N>` | `42` | RNG seed |

---

## Materials and periodicity

### `cell`

Create a unit cell and print its parameters.

```
sci-form cell [OPTIONS] --a <A> --b <B> --c <C>
```

| Option | Default | Description |
|--------|---------|-------------|
| `--a`, `--b`, `--c` | (required) | Lattice parameters in Å |
| `--alpha`, `--beta`, `--gamma` | `90.0` | Lattice angles in degrees |

```bash
sci-form cell --a 5.43 --b 5.43 --c 5.43
```

Output: `{"a": ..., "b": ..., "c": ..., "alpha": ..., "beta": ..., "gamma": ..., "volume": ..., "lattice": [...]}`.

---

### `assemble`

Assemble a framework crystal structure from a topology.

```
sci-form assemble [OPTIONS]
```

| Option | Default | Description |
|--------|---------|-------------|
| `-t, --topology <T>` | `pcu` | Framework topology: `pcu`, `dia` or `sql` |
| `--a <F>` | `10.0` | Lattice parameter a in Å |
| `--metal <N>` | `30` | Metal atomic number |
| `--geometry <G>` | `octahedral` | Coordination geometry: `linear`, `trigonal`, `tetrahedral`, `square_planar` or `octahedral` |
| `--supercell <N>` | `1` | Replication factor |

```bash
sci-form assemble --topology pcu --metal 30 --geometry octahedral --a 10.0
```

---

### `space-group`

Look up a crystallographic space group by number or symbol.

```
sci-form space-group [OPTIONS]
```

| Option | Description |
|--------|-------------|
| `--number <N>` | ITC number (1–230) |
| `--symbol <S>` | Hermann-Mauguin symbol, e.g. `"P 21/c"` |

---

### `all-space-groups`

List all 230 space groups.

```
sci-form all-space-groups
```

---

### `periodic-molecule`

Build a periodic molecular graph from fractional coordinates and a unit cell.

```
sci-form periodic-molecule [OPTIONS] --a <A> --b <B> --c <C> <ELEMENTS> <FRAC_COORDS>
```

| Option | Default | Description |
|--------|---------|-------------|
| `--a`, `--b`, `--c` | (required) | Lattice lengths in Å |
| `--alpha`, `--beta`, `--gamma` | `90.0` | Lattice angles in degrees |
| `--bond-tolerance <F>` | `0.3` | Extra bonding tolerance in Å (`0` selects the default) |

---

### `hapticity`

Detect haptic (η) metal–ring interactions and metallocene motifs.

```
sci-form hapticity [OPTIONS] <ELEMENTS_JSON> <COORDS_JSON>
```

| Option | Default | Description |
|--------|---------|-------------|
| `--bonds <JSON>` | *(empty)* | JSON array of `[i, j]` bond pairs; when omitted bonds are perceived from the geometry |

---

## Reaction paths

### `simplified-neb-path`

Compute a simplified NEB-like path between two geometries.

```
sci-form simplified-neb-path [OPTIONS] <SMILES> <START_COORDS> <END_COORDS>
```

| Option | Default | Description |
|--------|---------|-------------|
| `--n-images <N>` | `8` | Number of images along the path |
| `--n-iter <N>` | `50` | Number of relaxation iterations |
| `--spring-k <F>` | `0.1` | Harmonic spring constant between adjacent images |
| `--step-size <F>` | `0.01` | Cartesian update step size |
| `--method <M>` | `uff` | Energy backend: `uff`, `mmff94`, `pm3`, `xtb`, `gfn1`, `gfn2`, `hf3c` |

---

### `neb-energy` and `neb-gradient`

Single-point energy and energy + gradient with any NEB-capable backend.

```
sci-form neb-energy   [OPTIONS] <SMILES> <COORDS_JSON>
sci-form neb-gradient [OPTIONS] <SMILES> <COORDS_JSON>
```

| Option | Default | Description |
|--------|---------|-------------|
| `--method <M>` | `uff` | Energy backend: `uff`, `mmff94`, `pm3`, `xtb`, `gfn1`, `gfn2`, `hf3c` |

---

### `gsm-backend-plan` *(feature `alpha-gsm`)*

Show which GSM reaction-path backends are available for a molecule.

```
sci-form gsm-backend-plan <SMILES>
```

---

### `gsm-compare-backends` *(feature `alpha-gsm`)*

Compare GSM backends on the same geometry.

```
sci-form gsm-compare-backends [OPTIONS] <SMILES> <COORDS_JSON>
```

| Option | Default | Description |
|--------|---------|-------------|
| `--methods <JSON>` | `[]` | JSON array of backend names; `[]` means all |

---

### `analyze-gsm-mbh-htst-step` *(feature `alpha-kinetics`)*

Analyze a complete GSM + MBH + HTST elementary step.

```
sci-form analyze-gsm-mbh-htst-step [OPTIONS] <SMILES> <REACTANT_COORDS> <PRODUCT_COORDS>
```

| Option | Default | Description |
|--------|---------|-------------|
| `--method <M>` | `uff` | GSM backend: `uff`, `mmff94`, `pm3`, `xtb`, `gfn1`, `gfn2`, `hf3c`, `reaxff` |
| `--step-id <S>` | `elementary-step` | Stable identifier for the elementary step |
| `--temperature-k <F>` | `298.15` | Temperature in K |
| `--pressure-bar <F>` | `1.0` | Pressure in bar |
| `--n-nodes <N>` | `11` | Maximum number of GSM nodes |
| `--mbh-fd-step <F>` | `0.005` | MBH finite-difference step in Å |

---

### `trajectory-chart` and `reaction-coordinate-chart` *(feature `alpha-render-bridge` + `alpha-gsm`)*

Chart utilities for a GSM reaction path.

```
sci-form trajectory-chart          [OPTIONS] <SMILES> <REACTANT_COORDS> <PRODUCT_COORDS>
sci-form reaction-coordinate-chart [OPTIONS] <SMILES> <REACTANT_COORDS> <PRODUCT_COORDS>
```

| Option | Default | Description |
|--------|---------|-------------|
| `-m, --method <M>` | `uff` | GSM backend: `uff`, `mmff94`, `pm3`, `xtb`, `gfn1`, `gfn2`, `hf3c` |

---

### `reaction-dynamics3d` *(feature `alpha-reaction-dynamics`)*

Compute a full 3D reaction-dynamics path.

```
sci-form reaction-dynamics3d [OPTIONS] <REACTANT_SMILES> <PRODUCT_SMILES>
```

| Option | Default | Description |
|--------|---------|-------------|
| `--method <M>` | `gfn2` | Energy backend: `uff`, `xtb`, `gfn1`, `gfn2` (default), `pm3`, `hf3c` |
| `--n-images <N>` | `30` | Number of NEB images |
| `--n-approach <N>` | `15` | Number of approach frames |
| `--n-departure <N>` | `15` | Number of departure frames |
| `--smirks <S>` | — | Optional SMIRKS reaction string |
| `--seed <N>` | `42` | Random seed |

Both positional arguments are JSON arrays of SMILES, e.g. `'["[Cl-]", "CBr"]'` and `'["ClC", "[Br-]"]'`.

---

## Solvation, dispersion and electrochemistry

### `alpb` *(feature `experimental-alpb`)*

Compute ALPB implicit-solvation energy.

```
sci-form alpb [OPTIONS] <ELEMENTS_JSON> <COORDS_JSON>
```

| Option | Default | Description |
|--------|---------|-------------|
| `-c, --charges <JSON>` | *(empty)* | Per-atom charges |
| `-d, --dielectric <F>` | `78.5` | Solvent dielectric constant |

---

### `d4` *(feature `experimental-d4`)*

Compute the DFT-D4 dispersion correction.

```
sci-form d4 [OPTIONS] <ELEMENTS_JSON> <COORDS_JSON>
```

| Option | Default | Description |
|--------|---------|-------------|
| `--three-body` | off | Include the three-body ATM term |

---

### `cpm` *(feature `experimental-cpm`)*

Compute CPM charges at a fixed electrochemical potential.

```
sci-form cpm [OPTIONS] <ELEMENTS_JSON> <COORDS_JSON>
```

| Option | Default | Description |
|--------|---------|-------------|
| `-m, --mu-ev <F>` | `-4.44` | Electrochemical potential in eV |
| `-d, --dielectric <F>` | `78.5` | Solvent dielectric constant |

---

## Electrical double layer *(features `alpha-edl`)*

### `edl-profile`

Compute an EDL profile.

```
sci-form edl-profile [OPTIONS]
```

| Option | Default | Description |
|--------|---------|-------------|
| `-p, --potential <F>` | `0.05` | Surface potential (or surface charge) in V |
| `-m, --model <M>` | `gouy-chapman-stern` | EDL model: `helmholtz`, `gouy-chapman`, `gouy-chapman-stern` |
| `-i, --ionic-strength <F>` | `0.1` | Bulk ionic strength in mol/L |
| `-t, --temperature <F>` | `298.15` | Temperature in K |
| `-n, --n-points <N>` | `128` | Number of grid points |
| `-e, --extent <F>` | `12.0` | Maximum profile extent in Å |

### `edl-capacitance-scan`

Scan the EDL capacitance over a potential range.

```
sci-form edl-capacitance-scan [OPTIONS]
```

| Option | Default | Description |
|--------|---------|-------------|
| `--v-min <F>` | `-0.5` | Minimum potential in V |
| `--v-max <F>` | `0.5` | Maximum potential in V |
| `-n, --n-points <N>` | `50` | Number of grid points |
| `-i, --ionic-strength <F>` | `0.1` | Bulk ionic strength in mol/L |
| `-m, --model <M>` | `gouy-chapman-stern` | EDL model |

### `edl-helmholtz`, `edl-gouy-chapman`

Single-model EDL profiles.

```
sci-form edl-helmholtz    [OPTIONS]
sci-form edl-gouy-chapman [OPTIONS]
```

| Option | Default | Description |
|--------|---------|-------------|
| `-p, --potential <F>` | `0.05` | Surface potential in V |
| `-i, --ionic-strength <F>` | `0.1` | Bulk ionic strength in mol/L |
| `-t, --temperature <F>` | `298.15` | Temperature in K |

### `edl-gcs`

Gouy-Chapman-Stern EDL profile.

```
sci-form edl-gcs [OPTIONS]
```

| Option | Default | Description |
|--------|---------|-------------|
| `-p, --potential <F>` | `0.05` | Surface potential in V |
| `-i, --ionic-strength <F>` | `0.1` | Bulk ionic strength in mol/L |
| `-t, --temperature <F>` | `298.15` | Temperature in K |
| `-s, --stern-thickness <F>` | `3.0` | Stern-layer thickness in Å |

---

## Kinetics *(feature `alpha-kinetics`)*

### `htst-rate`

Evaluate an HTST transition rate.

```
sci-form htst-rate --activation-free-energy <EV> --reaction-free-energy <EV>
```

| Option | Default | Description |
|--------|---------|-------------|
| `--activation-free-energy <F>` | (required) | Activation free energy in eV |
| `--reaction-free-energy <F>` | (required) | Reaction free energy in eV |
| `-t, --temperature <F>` | `298.15` | Temperature in K |

### `solve-microkinetic-network`

Integrate a microkinetic network.

```
sci-form solve-microkinetic-network [OPTIONS]
```

| Option | Default | Description |
|--------|---------|-------------|
| `--network <JSON>` | built-in demo network | `MicrokineticNetworkConfig` as JSON (file path or inline) |
| `--n-steps <N>` | `5` | Number of integration frames |

### `solve-microkinetic-steady-state`

Solve the microkinetic steady state.

```
sci-form solve-microkinetic-steady-state [OPTIONS]
```

| Option | Default | Description |
|--------|---------|-------------|
| `--network <JSON>` | built-in demo network | `MicrokineticNetworkConfig` as JSON (file path or inline) |
| `--n-steps <N>` | `5` | Number of integration frames |
| `--tol <F>` | `1e-10` | Convergence tolerance on the population change |

### `extract-kinetics-diagnostics`

Extract kinetics diagnostics from a microkinetic trace.

```
sci-form extract-kinetics-diagnostics [OPTIONS]
```

| Option | Default | Description |
|--------|---------|-------------|
| `--network <JSON>` | built-in demo network | `MicrokineticNetworkConfig` as JSON (file path or inline) |

---

## Periodic linear algebra *(feature `alpha-periodic-linear`)*

### `kmesh`

Generate a k-mesh.

```
sci-form kmesh [OPTIONS]
```

| Option | Default | Description |
|--------|---------|-------------|
| `--nx <N>` | `4` | Mesh divisions along the first reciprocal axis |
| `--ny <N>` | `4` | Mesh divisions along the second reciprocal axis |
| `--nz <N>` | `1` | Mesh divisions along the third reciprocal axis |
| `-c, --centering <C>` | `monkhorst-pack` | Centering: `monkhorst-pack` or `gamma` |

### `compute-periodic-dos`

Compute the k-averaged periodic DOS with the KPM solver.

```
sci-form compute-periodic-dos [OPTIONS] <HAMILTONIAN_BLOCKS>
```

| Option | Default | Description |
|--------|---------|-------------|
| `--nx`, `--ny`, `--nz` | `2` | Mesh divisions along the reciprocal axes |
| `-c, --centering <C>` | `monkhorst-pack` | Centering |
| `--order <N>` | `100` | Chebyshev expansion order |
| `--e-min <F>` | `-20.0` | Output minimum energy in eV |
| `--e-max <F>` | `10.0` | Output maximum energy in eV |
| `--n-points <N>` | `512` | Number of output grid points |
| `--fermi-energy-ev <F>` | — | Fermi energy in eV |

### `solve-periodic-randnla`

Solve the periodic RandNLA problem.

```
sci-form solve-periodic-randnla [OPTIONS] <HAMILTONIAN_BLOCKS> <OVERLAP_BLOCKS>
```

| Option | Default | Description |
|--------|---------|-------------|
| `--n-electrons <N>` | `2` | Number of valence electrons |
| `--nx`, `--ny`, `--nz` | `2` | Mesh divisions along the reciprocal axes |
| `-c, --centering <C>` | `monkhorst-pack` | Centering |
| `--sketch-size <N>` | engine default | Nyström sketch rank |

### `bloch-phase`

Compute a Bloch phase.

```
sci-form bloch-phase [OPTIONS]
```

| Option | Default | Description |
|--------|---------|-------------|
| `--k-x`, `--k-y`, `--k-z` | `0.0` | Fractional k coordinates |
| `--t-x`, `--t-y`, `--t-z` | `1`, `0`, `0` | Lattice translation indices |

### `build-bloch-hamiltonian`

Build the Bloch Hamiltonian.

```
sci-form build-bloch-hamiltonian [OPTIONS] <BLOCKS>
```

| Option | Default | Description |
|--------|---------|-------------|
| `--k-x`, `--k-y`, `--k-z` | `0.0` | Fractional k coordinates |

### `assemble-periodic-operators`

Assemble the periodic operators for every k-point.

```
sci-form assemble-periodic-operators [OPTIONS] <HAMILTONIAN_BLOCKS> [OVERLAP_BLOCKS]
```

| Option | Default | Description |
|--------|---------|-------------|
| `--nx`, `--ny`, `--nz` | `2` | Mesh divisions along the reciprocal axes |
| `-c, --centering <C>` | `monkhorst-pack` | Centering |

### `validate-electron-count`

Validate the electron count of a periodic calculation.

```
sci-form validate-electron-count --n-electrons <N> --n-bands <N> --n-kpoints <N>
```

All three options are required.

---

## Chart utilities *(feature `alpha-render-bridge`)*

### `chart-to-json`, `chart-from-json`, `validate-chart-schema`

Serialize, deserialize and validate chart payloads.

```
sci-form chart-to-json        [JSON]
sci-form chart-from-json      [JSON]
sci-form validate-chart-schema [JSON]
```

The `JSON` positional argument is optional; it accepts an inline JSON string.

### `edl-chart` *(also needs `alpha-edl`)*

```
sci-form edl-chart [OPTIONS]
```

| Option | Default | Description |
|--------|---------|-------------|
| `-p, --potential <F>` | `0.1` | Surface potential in V |
| `-i, --ionic-strength <F>` | `0.1` | Bulk ionic strength in mol/L |

### `capacitance-chart` *(also needs `alpha-edl`)*

```
sci-form capacitance-chart [OPTIONS]
```

| Option | Default | Description |
|--------|---------|-------------|
| `--v-min <F>` | `-0.5` | Minimum potential in V |
| `--v-max <F>` | `0.5` | Maximum potential in V |
| `-n, --n-points <N>` | `50` | Number of grid points |
| `-i, --ionic-strength <F>` | `0.1` | Bulk ionic strength in mol/L |
| `-m, --model <M>` | `gouy-chapman` | EDL model |

### `arrhenius-chart`, `thermal-prop-chart` *(also need `alpha-kinetics`)*

```
sci-form arrhenius-chart --temperatures <JSON> --rates <JSON>
sci-form thermal-prop-chart --activation-free-energy <EV> --reaction-free-energy <EV> --temperatures <JSON>
```

### `phase-portrait-chart` *(also needs `alpha-kinetics`)*

```
sci-form phase-portrait-chart [OPTIONS]
```

| Option | Default | Description |
|--------|---------|-------------|
| `--network <JSON>` | built-in demo network | `MicrokineticNetworkConfig` JSON (file path or inline) |
| `--n-steps <N>` | `50` | Number of integration frames |

### `band-structure-chart`, `dos-chart`, `kpoint-path-chart`, `fermi-surface-chart` *(also need `alpha-periodic-linear`)*

```
sci-form band-structure-chart [OPTIONS] <ELEMENTS> <COORDS> <LATTICE>
sci-form dos-chart --energies <JSON> --dos-values <JSON>
sci-form kpoint-path-chart [OPTIONS]
sci-form fermi-surface-chart [OPTIONS] <BLOCKS>
```

| Command | Options |
|---------|---------|
| `band-structure-chart` | `--n-kpoints-per-segment <N>` (default `20`) |
| `dos-chart` | `-e, --energies <JSON>`, `-d, --dos-values <JSON>` (both required) |
| `kpoint-path-chart` | `--n-points-per-segment <N>` (default `20`) |
| `fermi-surface-chart` | `--fermi-energy-ev <F>`, `--n-kx <N>` (default `24`), `--n-ky <N>` (default `24`) |

---

## Output Formats

### JSON

Default format. Returns the full `ConformerResult` as compact JSON.

### XYZ

Standard XYZ format for molecular viewers (Avogadro, VESTA, VMD):

```
9
CCO
C     -0.893064     0.069264     0.063601
C      0.525542    -0.422006    -0.147153
O      1.389964     0.450482     0.467233
...
```

### SDF

MDL Structure Data File format. Compatible with RDKit, OpenBabel, MarvinSuite.

---

## Exit Codes

| Code | Meaning |
|------|---------|
| `0` | Success |
| `1` | Runtime failure: invalid SMILES, embedding failure, missing input file, I/O error, or a failed calculation |
| `2` | Usage error reported by the argument parser (unknown subcommand, missing or invalid flag) |

---

## Pipeline Examples

### Filter successful results

```bash
sci-form batch -i molecules.smi | jq -c 'select(.error == null)'
```

### Get batch statistics

```bash
sci-form batch -i molecules.smi | jq -s '{
  total: length,
  success: [.[] | select(.error == null)] | length,
  avg_time_ms: ([.[] | .time_ms] | add / length)
}'
```

### Convert batch to XYZ

```bash
sci-form batch -i molecules.smi -f xyz > all_conformers.xyz
```

### Compute charges for all molecules in a file

```bash
while read smiles; do
  sci-form charges "$smiles"
done < molecules.smi | jq -s '.'
```

### Parallel processing with GNU Parallel

```bash
cat huge_set.smi | parallel -j8 --pipe -N100 sci-form batch /dev/stdin > results.json
```
