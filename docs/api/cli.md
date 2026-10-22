# CLI Reference

The `sci-form` binary exposes **100 subcommands**: 66 are always built and 34 are behind Cargo
feature flags. This page is the exhaustive reference; the [CLI guide](../guide/cli.md) is the
task-oriented walkthrough.

---

## Installation

### Pre-built binaries

Download from [GitHub Releases](https://github.com/jigonzalez930209/sci-form/releases):

- `sci-form-x86_64-unknown-linux-gnu` — Linux x86_64
- `sci-form-x86_64-apple-darwin` — macOS Intel
- `sci-form-aarch64-apple-darwin` — macOS Apple Silicon
- `sci-form-x86_64-pc-windows-msvc.exe` — Windows x86_64
- `sci-form-aarch64-unknown-linux-gnu` — Linux ARM64

### Build from source

```bash
cargo install sci-form-cli
```

### Feature-gated subcommands

The optional subcommands are compiled in only when the corresponding Cargo feature is enabled on
the `sci-form-cli` package:

```bash
cargo install sci-form-cli --features alpha-gsm,alpha-edl,alpha-kinetics
```

| Feature | Unlocks |
|---|---|
| `experimental-alpb` | `alpb` |
| `experimental-d4` | `d4` |
| `experimental-cpm` | `cpm` |
| `alpha-gsm` | `gsm-backend-plan`, `gsm-compare-backends` |
| `alpha-edl` | `edl-profile`, `edl-capacitance-scan`, `edl-helmholtz`, `edl-gouy-chapman`, `edl-gcs` |
| `alpha-kinetics` | `htst-rate`, `extract-kinetics-diagnostics`, `solve-microkinetic-network`, `solve-microkinetic-steady-state`, `analyze-gsm-mbh-htst-step` |
| `alpha-periodic-linear` | `kmesh`, `compute-periodic-dos`, `solve-periodic-randnla`, `bloch-phase`, `build-bloch-hamiltonian`, `assemble-periodic-operators`, `validate-electron-count` |
| `alpha-render-bridge` | `edl-chart`, `capacitance-chart`, `arrhenius-chart`, `thermal-prop-chart`, `phase-portrait-chart`, `band-structure-chart`, `dos-chart`, `kpoint-path-chart`, `fermi-surface-chart`, `trajectory-chart`, `reaction-coordinate-chart` |
| `alpha-reaction-dynamics` | `reaction-dynamics3d` |

The `alpha-render-bridge` chart commands need a second feature as well: `edl-chart` and
`capacitance-chart` also require `alpha-edl`; `arrhenius-chart`, `thermal-prop-chart` and
`phase-portrait-chart` require `alpha-kinetics`; `band-structure-chart`, `dos-chart`,
`kpoint-path-chart` and `fermi-surface-chart` require `alpha-periodic-linear`; and
`trajectory-chart` / `reaction-coordinate-chart` require `alpha-gsm`.

---

## Conventions

### Optional arguments are absent

An optional argument is simply **absent**. There are no magic sentinels: you never pass `""`,
`"[]"` or `0.0` to mean "not given". Omit the flag and the engine default applies.

```bash
sci-form stereo "C(F)(Cl)(Br)I"                          # topology-only
sci-form stereo "C(F)(Cl)(Br)I" --coords "[1.0, 0, 0]"   # with 3D coordinates
```

### Error messages

Every user-facing error follows one template:

```text
<What failed> '<actual value>'. Supported: 'a', 'b', 'c'.
```

The offending input is always included and the valid values are always listed when there is a
fixed set.

```bash
$ sci-form dos-multimethod "[6,6,8]" "[0,0,0,1,0,0,2,0,0]" --method bogus
Invalid 'method' 'bogus'. Supported: 'eht', 'pm3', 'xtb', 'gfn1', 'gfn2', 'hf3c'.

$ sci-form sasa "notjson" "[]"
Invalid 'elements' JSON 'notjson'. (expected ident at line 1 column 2)

$ sci-form charges "C("
Charges error: Expected ')' at position 2 (found 'end of input') in "C("
```

### Negative numeric option values

Several options take negative numbers (`--e-min`, `--e-max`, `--mu-ev`, `--v-min`, `--v-max`,
`--reaction-free-energy`, …). The argument parser treats a bare `-4.44` as a flag, so negative
values must be attached with `=`:

```bash
sci-form cpm "$EL" "$COORDS" --mu-ev=-4.44     # correct
sci-form cpm "$EL" "$COORDS" --mu-ev -4.44      # usage error
```

### A reusable geometry

Most commands take a geometry as two positional JSON arguments. Build it once and reuse it:

```bash
RESULT=$(sci-form embed "CCO")
EL=$(echo "$RESULT" | jq -c '.elements')      # [6, 6, 8, 1, 1, 1, 1, 1, 1]
COORDS=$(echo "$RESULT" | jq -c '.coords')    # flat [x, y, z, ...] array in Å
```

The examples below use `$EL` and `$COORDS` for those two arguments.

### Units

| Quantity | Unit |
|---|---|
| Geometry (coordinates, bond lengths, cell parameters) | Å (cell angles in degrees) |
| UFF / MMFF94 energies (`uff`, `neb-energy`, `neb-gradient`) | kcal/mol |
| UFF / MMFF94 gradients (`neb-gradient`) | kcal/mol/Å |
| EHT, PM3, xTB orbital energies and HOMO–LUMO gaps | eV |
| PM3 / xTB total, electronic, repulsive, dispersion energies | eV |
| PM3 heat of formation | kcal/mol |
| HF-3c energy, HF energy, nuclear repulsion, D3/gCP/SRB | Hartree |
| HF-3c orbital energies, HOMO/LUMO, gap | Hartree |
| ANI energy and per-atom energies | Hartree (forces Hartree/Å) |
| CISD excitation `energy`, `correlation_energy` | Hartree (`energy_ev` is eV) |
| EHT gradients, `optimize-geometry-eht` energy | eV, eV/Å |
| `optimize-framework` energy and forces | eV, eV/Å |
| IR wavenumbers, vibrational frequencies | cm⁻¹ |
| J-couplings (`ensemble-j-couplings`) | Hz |
| Dipole moment | Debye |
| Partial charges | elementary charge (e) |
| Electrostatic potential | Hartree/e |
| DOS / PDOS | 1/eV |
| Surface area | Å² |
| MD energies | kcal/mol; times in fs; temperatures in K |
| Pressure | bar |

---

## Command list

### Always built (66)

`embed`, `batch`, `parse`, `info`, `eht`, `eht-gradient`, `optimize-geometry-eht`, `charges`,
`charges-configured`, `sasa`, `solvation`, `eeq`, `population`, `dipole`, `esp`, `dos`,
`dos-multimethod`, `pm3`, `xtb`, `gfn1`, `gfn2`, `hf3c`, `cisd`, `ani`, `aevs-tm`, `uff`,
`optimize-framework`, `stereo`, `sssr`, `ecfp`, `ecfp-batch`, `tanimoto`, `rmsd`,
`canonical-smiles`, `embed-diverse`, `filter-diverse`, `ir-spectrum-broadened`, `assign-peaks`,
`ensemble-j-couplings`, `whim`, `rdf`, `getaway`, `train-random-forest`,
`train-gradient-boosting`, `cell`, `assemble`, `space-group`, `all-space-groups`,
`periodic-molecule`, `hapticity`, `simplified-neb-path`, `neb-energy`, `neb-gradient`,
`chart-to-json`, `chart-from-json`, `validate-chart-schema`, `uhf`, `rohf`,
`pm3-gradient`, `xtb-gradient`, `uff-breakdown`, `mmff94-breakdown`, `tpsa`,
`ensemble`, `parse-cif`, `write-cif`.

### Feature-gated (34)

`alpb` (`experimental-alpb`), `d4` (`experimental-d4`), `cpm` (`experimental-cpm`),
`edl-profile`, `edl-capacitance-scan`, `edl-helmholtz`, `edl-gouy-chapman`, `edl-gcs`
(`alpha-edl`), `htst-rate`, `extract-kinetics-diagnostics`, `solve-microkinetic-network`,
`solve-microkinetic-steady-state`, `analyze-gsm-mbh-htst-step` (`alpha-kinetics`), `kmesh`,
`compute-periodic-dos`, `solve-periodic-randnla`, `bloch-phase`, `build-bloch-hamiltonian`,
`assemble-periodic-operators`, `validate-electron-count` (`alpha-periodic-linear`),
`gsm-backend-plan`, `gsm-compare-backends` (`alpha-gsm`), `edl-chart`, `capacitance-chart`,
`arrhenius-chart`, `thermal-prop-chart`, `phase-portrait-chart`, `band-structure-chart`,
`dos-chart`, `kpoint-path-chart`, `fermi-surface-chart`, `trajectory-chart`,
`reaction-coordinate-chart` (`alpha-render-bridge` + the matching domain feature),
`reaction-dynamics3d` (`alpha-reaction-dynamics`).

---

## Conformer generation and structure

### `embed`

Generate a 3D conformer for a single SMILES.

```
sci-form embed <SMILES> [OPTIONS]
```

| Option | Default | Description |
|--------|---------|-------------|
| `-s, --seed <N>` | `42` | RNG seed |
| `-f, --format <FMT>` | `json` | Output format: `json`, `xyz`, `sdf` |

```bash
$ sci-form embed "CCO"
{
  "smiles": "CCO",
  "num_atoms": 9,
  "coords": [-0.8930637240409851, 0.06926389038562775, 0.06360100954771042, "..."],
  "elements": [6, 6, 8, 1, 1, 1, 1, 1, 1],
  "bonds": [[0, 1, "SINGLE"], [1, 2, "SINGLE"], [0, 3, "SINGLE"], "..."],
  "error": null,
  "time_ms": 1.530118
}
```

(`time_ms` is the wall-clock time of that call and varies from run to run.)

```bash
$ sci-form embed "c1ccccc1" -f xyz
12
c1ccccc1
C      0.801870    -1.121113    -0.008347
C      1.386838     0.149298    -0.028242
C      0.584813     1.253631    -0.019782
...
```

---

### `batch`

Process multiple SMILES from a file or stdin. The input path is an **option**, not a positional
argument.

```
sci-form batch [OPTIONS]
```

| Option | Default | Description |
|--------|---------|-------------|
| `-i, --input <FILE>` | stdin | Input file (omit to read stdin) |
| `-o, --output <FILE>` | stdout | Output file (omit to write stdout) |
| `-s, --seed <N>` | `42` | RNG seed |
| `-t, --threads <N>` | `0` | Number of threads (`0` = auto) |
| `-f, --format <FMT>` | `json` | Output format: `json`, `xyz`, `sdf` |

**Input file format:** exactly one SMILES per line. A trailing name is **not** supported — the
whole line is parsed as SMILES:

```text
CCO
c1ccccc1
CC(=O)O
```

```bash
$ sci-form batch -i molecules.smi
Processing molecules from molecules.smi (auto threads)...
[{"smiles":"CCO","num_atoms":9,"coords":[-0.8930637240409851,...],"elements":[6,6,8,1,1,1,1,1,1],"bonds":[[0,1,"SINGLE"],...],"error":null,"time_ms":7.475605}, ...]
Done: 3/3 OK, 0 failed, 0.0s total (2.9 ms/mol)
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
    2: elem= 6 hyb=SP2 charge=0 h=0
    3: elem= 6 hyb=SP2 charge=0 h=0
    4: elem= 6 hyb=SP2 charge=0 h=0
    5: elem= 6 hyb=SP2 charge=0 h=0
    6: elem= 1 hyb=Unknown charge=0 h=0
    ...
```

---

### `info`

Print the build information: version, feature summary, supported output formats and the
available bindings.

```
sci-form info
```

```bash
$ sci-form info
sci-form 1.0.0
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
$ sci-form embed-diverse "CCO" --n-conformers 2
[
  {
    "smiles": "CCO",
    "num_atoms": 9,
    "coords": [-0.8822969794273376, -0.19749876856803894, -0.05994144827127457, "..."],
    "elements": [6, 6, 8, 1, 1, 1, 1, 1, 1],
    "bonds": [[0, 1, "SINGLE"], [1, 2, "SINGLE"], [0, 3, "SINGLE"], "..."],
    "error": null,
    "time_ms": 1.09874
  },
  ...
]
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

`CONFORMERS` is a JSON array of conformers, each one a flat xyz coordinate array. The command
prints the JSON array of the retained indices.

```bash
$ sci-form filter-diverse "[$COORDS]"
[
 0
]
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

Compute Gasteiger-Marsili partial charges. Charges are in **elementary charge (e)**.

```
sci-form charges <SMILES>
```

```bash
$ sci-form charges "CCO"
{
  "charges": [
    -0.04182525705714308,
    0.041384322456729204,
    -0.39527742524721593,
    0.025164166823372726,
    0.025164166823372726,
    0.025164166823372726,
    0.05542120163292069,
    0.05542120163292069,
    0.20938345611167028
  ],
  "iterations": 6,
  "total_charge": 8.326672684688674e-17,
  "converged": false
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
$ sci-form charges-configured "CCO" --max-iter 12
{
  "charges": [-0.04183848615477239, 0.04022058174908642, -0.3966637147812598, "..."],
  "iterations": 12,
  "total_charge": -2.7755575615628914e-17,
  "converged": false
}
```

---

### `sasa`

Compute solvent-accessible surface area. Areas are in **Å²**.

```
sci-form sasa [OPTIONS] <ELEMENTS_JSON> <COORDS_JSON>
```

| Option | Default | Description |
|--------|---------|-------------|
| `-p, --probe-radius <F>` | `1.4` | Probe radius in Å |

```bash
$ sci-form sasa "$EL" "$COORDS"
{
  "total_sasa": 198.45192083308308,
  "atom_sasa": [7.296087138574494, 3.3964543576122646, 30.692836626796677, "..."],
  "probe_radius": 1.4,
  "num_points": 960
}
```

---

### `solvation`

Compute the non-polar (SASA-based) solvation energy; with `--charges` it also returns the
Generalized Born electrostatic term. Energies in **kcal/mol**.

```
sci-form solvation [OPTIONS] <ELEMENTS_JSON> <COORDS_JSON>
```

| Option | Default | Description |
|--------|---------|-------------|
| `-c, --charges <JSON>` | *absent* | Per-atom charges as a JSON array. When absent only the non-polar term is computed |
| `-p, --probe-radius <F>` | `1.4` | Probe radius in Å |

The GB term uses a solvent dielectric of 78.5 and an internal dielectric of 1.0; these are fixed
in the CLI.

```bash
# Using the SMILES pipeline
RESULT=$(sci-form embed "CCO")
ELEM=$(echo "$RESULT" | jq '.elements')
COORDS=$(echo "$RESULT" | jq '.coords')
sci-form solvation "$ELEM" "$COORDS"

# With GB electrostatics
sci-form solvation "$ELEM" "$COORDS" --charges "$(sci-form charges CCO | jq -c '.charges')"
```

Output (non-polar only):

```json
{
  "energy_kcal_mol": -3.86723458312331,
  "atom_contributions": [0.08755304566289393, 0.040757452291347175, -5.095010880048249, "..."],
  "atom_sasa": [7.296087138574494, 3.3964543576122646, 30.692836626796677, "..."],
  "total_sasa": 198.45192083308308
}
```

Output (with charges):

```json
{
  "electrostatic_energy_kcal_mol": -8.72532248578084e-05,
  "nonpolar_energy_kcal_mol": -3.86723458312331,
  "total_energy_kcal_mol": -3.8673218363481676,
  "born_radii": [50.0, 50.0, 50.0, "..."],
  "charges": [-0.04182525705714308, 0.0413843224567292, -0.39527742524721593, "..."],
  "solvent_dielectric": 78.5,
  "solute_dielectric": 1.0
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

```bash
$ sci-form eeq "$EL" "$COORDS"
{
  "charges": [-0.24266719355854416, -0.025089103301529116, -0.4990266160632522, "..."],
  "coordination_numbers": [1.8598524075112612, 2.2522901147968173, 1.6024019714740372, "..."],
  "total_charge": -1.5265566588595902e-15
}
```

---

### `stereo`

Analyze stereochemistry: R/S stereocenters and E/Z double bonds.

```
sci-form stereo [OPTIONS] <SMILES>
```

| Option | Default | Description |
|--------|---------|-------------|
| `-c, --coords <JSON>` | *absent* | Flat coordinate array as a JSON string. When absent the analysis is topology-only |

```bash
$ sci-form stereo "C(F)(Cl)(Br)I"
{
  "stereocenters": [
    {
      "atom_index": 0,
      "element": 6,
      "substituent_indices": [4, 3, 2, 1],
      "priorities": [1, 2, 3, 4],
      "configuration": null
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

With `--coords` the `configuration` field is filled in (`"R"` here).

---

### `sssr`

Perceive the Smallest Set of Smallest Rings.

```
sci-form sssr <SMILES>
```

```bash
$ sci-form sssr "c1ccc2ccccc2c1"
{
  "rings": [
    { "atoms": [0, 9, 8, 3, 2, 1], "size": 6, "is_aromatic": true },
    { "atoms": [3, 8, 7, 6, 5, 4], "size": 6, "is_aromatic": true }
  ],
  "atom_ring_count": [1, 1, 1, "..."],
  "atom_ring_sizes": [[6], [6], [6], "..."],
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
$ sci-form ecfp "c1ccccc1"
{
  "n_bits": 2048,
  "on_bits": [78, 373, 587, 1798, 1963, 2040],
  "radius": 2,
  "raw_features": [9995780334538829739, 9995780334538829739, "..."]
}
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
$ sci-form ecfp-batch '["CCO", "c1ccccc1"]'
[
  { "n_bits": 2048, "on_bits": [78, 236, 298, "..."], "radius": 2, "raw_features": [...] },
  { "n_bits": 2048, "on_bits": [78, 373, 587, "..."], "radius": 2, "raw_features": [...] }
]
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
{
  "tanimoto": 0.333333
}
```

---

### `rmsd`

Compute the RMSD between two coordinate sets (Kabsch alignment). RMSD in **Å**.

```
sci-form rmsd <COORDS_JSON> <REFERENCE_JSON>
```

```bash
$ sci-form rmsd "$COORDS" "$COORDS"
{
  "rmsd": 0.0,
  "rotation": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
}
```

---

## Electronic structure

All of these take the same two positional arguments: `ELEMENTS_JSON` (JSON array of atomic
numbers) and `COORDS_JSON` (flat `[x, y, z, ...]` array in Å).

### `eht`

Run Extended Hückel Theory. Orbital energies and the gap are in **eV**.

```
sci-form eht [OPTIONS] <ELEMENTS_JSON> <COORDS_JSON>
```

| Option | Default | Description |
|--------|---------|-------------|
| `-k, --k <F>` | `1.75` | Wolfsberg-Helmholtz constant (omit the flag for the default) |

```bash
$ sci-form eht "$EL" "$COORDS"
{
  "energies": [-32.67566514175308, -24.111381812665172, -21.774536672315957, "..."],
  "coefficients": [[0.006099259985427686, -0.5943455200430089, "..."], "..."],
  "n_electrons": 20,
  "homo_index": 9,
  "lumo_index": 10,
  "homo_energy": -13.41902554862098,
  "lumo_energy": -10.31046604789137,
  "gap": 3.1085595007296085,
  "support": {
    "level": "supported",
    "has_transition_metals": false,
    "supported_elements": [1, 6, 8],
    "provisional_elements": [],
    "unsupported_elements": [],
    "warnings": []
  }
}
```

---

### `eht-gradient`

Compute EHT energy-surface gradients by central finite differences. Gradients in **eV/Å**.

```
sci-form eht-gradient <ELEMENTS_JSON> <COORDS_JSON>
```

```bash
$ sci-form eht-gradient "$EL" "$COORDS"
{
  "gradients": [[303.8476992855976, -41.757622227578395, -9.075754780951684], "..."],
  "rms_gradient": 141.26950682295026,
  "max_gradient": 407.27235699478114,
  "energy": 1886.0591362006721
}
```

---

### `optimize-geometry-eht`

Optimize a geometry on the EHT energy surface with the shared BFGS engine. Energy in **eV**.

```
sci-form optimize-geometry-eht [OPTIONS] <ELEMENTS_JSON> <COORDS_JSON>
```

| Option | Default | Description |
|--------|---------|-------------|
| `--max-steps <N>` | `200` | Maximum number of optimization steps |
| `--grad-threshold <F>` | `0.01` | RMS gradient convergence threshold in eV/Å |
| `--energy-threshold <F>` | `1e-6` | Energy convergence threshold in eV |
| `--step-size <F>` | `0.05` | Initial step size for steepest descent in Å |

```bash
$ sci-form optimize-geometry-eht "$EL" "$COORDS"
{
  "positions": [[-2.228681350704961, 3.742950765522979, -8.305902124169542], "..."],
  "energy": -56.40417039919288,
  "n_steps": 200,
  "converged": false,
  "rms_gradient": 2.7638612117764696,
  "energies": [1886.0591362006721, 1751.3898905746394, "..."]
}
```

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

PM3 (all energies in **eV**, `heat_of_formation` in **kcal/mol**):

```bash
$ sci-form pm3 "$EL" "$COORDS"
{
  "orbital_energies": [-38.82711576123044, -30.929875185107743, -22.646733347120808, "..."],
  "electronic_energy": -1761.7957179078226,
  "nuclear_repulsion": 1138.2523453052277,
  "total_energy": -623.5433726025949,
  "heat_of_formation": -53.997668255763756,
  "n_basis": 18,
  "n_electrons": 20,
  "homo_energy": -11.091818031628284,
  "lumo_energy": 3.266643650908672,
  "gap": 14.358461682536955,
  "mulliken_charges": [-0.1724405535934732, 0.08956201679898346, "..."],
  "scf_iterations": 32,
  "converged": true,
  "notes": ["Produced by the faithful MOPAC-parity PM3 solver."]
}
```

GFN0-xTB (all energies in **eV**):

```bash
$ sci-form xtb "$EL" "$COORDS"
{
  "orbital_energies": [-18.709580590018355, -17.072433084117943, "..."],
  "electronic_energy": -312.96057303129896,
  "repulsive_energy": 3.0694987989800326,
  "total_energy": -309.8910594776795,
  "n_basis": 18,
  "n_electrons": 20,
  "homo_energy": -11.113088168797308,
  "lumo_energy": 1.592666900269751,
  "gap": 12.70575506906706,
  "mulliken_charges": [-0.11257146455801315, 0.09593650567551615, "..."],
  "scc_iterations": 9,
  "converged": true,
  "notes": []
}
```

GFN1-xTB adds `dispersion_energy`, `dispersion_e6`, `dispersion_e8`, `halogen_bond_energy` and
`shell_charges`; GFN2-xTB adds `dispersion_2_body_energy`, `atomic_dipoles` (Debye) and
`atomic_quadrupoles`.

HF-3c (all energies in **Hartree**):

```bash
$ sci-form hf3c "$EL" "$COORDS"
{
  "energy": -151.23033482038505,
  "hf_energy": -152.1174873874554,
  "nuclear_repulsion": 82.58742066152247,
  "d3_energy": -0.023907644680519688,
  "gcp_energy": 0.009399471038504277,
  "srb_energy": 0.9016607407123806,
  "orbital_energies": [-20.2397700816775, -11.097490693337386, "..."],
  "scf_iterations": 16,
  "converged": true,
  "cis": { "excitations": [...], "error": null },
  "n_basis": 21,
  "n_electrons": 26,
  "homo_energy": -0.33991491972198973,
  "lumo_energy": 0.5640256248572378,
  "gap": 0.9039405445792276,
  "mulliken_charges": [-0.1726611892372265, 0.035806063383814646, "..."]
}
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

```bash
$ sci-form cisd "$EL" "$COORDS" --n-states 1
{
  "excitations": [
    {
      "energy": 0.4585478005172503,
      "energy_ev": 12.477721208496853,
      "wavelength_nm": 99.3644519927016,
      "oscillator_strength": 0.30569853367816624,
      "dominant_transition": [[12, 13]],
      "character": "singles",
      "singles_weight": 1.0
    }
  ],
  "correlation_energy": -0.13837407597225154,
  "n_csfs": 15393,
  "approximation": "perturbative"
}
```

---

### `population`

Mulliken & Löwdin population analysis (requires the EHT engine).

```
sci-form population <ELEMENTS_JSON> <COORDS_JSON>
```

```bash
$ sci-form population "$EL" "$COORDS"
{
  "mulliken_charges": [-0.09449102702381129, 0.4563487412194829, -1.2894376098578126, "..."],
  "lowdin_charges": [-0.05092419844160023, 0.4799987561029999, -1.2376280130457857, "..."],
  "mulliken_populations": [1.4166064141980024, 0.8184169414172641, "..."],
  "num_atoms": 9,
  "total_charge_mulliken": -1.4988010832439613e-14,
  "total_charge_lowdin": -2.5757174171303632e-14,
  "charge_conservation_error": 1.4988010832439613e-14,
  "homo_energy": -13.41902554862098,
  "lumo_energy": -10.31046604789137,
  "gap": 3.1085595007296085,
  "warnings": []
}
```

---

### `dipole`

Compute the molecular dipole moment (Debye).

```
sci-form dipole <ELEMENTS_JSON> <COORDS_JSON>
```

```bash
$ sci-form dipole "$EL" "$COORDS"
{
  "vector": [-3.409698969668939, -0.7030268699182814, -3.418342844600006],
  "magnitude": 4.879073851338651,
  "unit": "Debye",
  "nuclear_dipole": [28.085597083396912, 5.735898877275735, 10.017142454130425],
  "electronic_dipole": [-31.495296053065843, -6.438925747194016, -13.435485298730434]
}
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

Output: `{"sigma": ..., "data": [[energy, dos], ...]}` with the energy in **eV** and the DOS in
**1/eV**.

```bash
$ sci-form dos "$EL" "$COORDS" -n 3
{
  "sigma": 0.3,
  "data": [[-30.0, 0.0], [-12.5, 0.0122], [5.0, 0.0]]
}
```

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

```bash
$ sci-form dos-multimethod "$EL" "$COORDS" -n 3
{
  "energies": [-30.0, -12.5, 5.0],
  "total_dos": [7.087422781333044e-18, 0.012199607342232369, 5.750197910390419e-187],
  "method": "Eht",
  "homo_energy": -13.41902554862098,
  "lumo_energy": -10.31046604789137,
  "gap": 3.1085595007296085,
  "orbital_energies": [-32.67566514175308, "..."],
  "sigma": 0.3
}
```

---

### `esp`

Compute the electrostatic potential on a 3D grid. Values in **Hartree/e**.

```
sci-form esp [OPTIONS] <ELEMENTS_JSON> <COORDS_JSON>
```

| Option | Default | Description |
|--------|---------|-------------|
| `-s, --spacing <F>` | `0.5` | Grid spacing in Å |
| `-p, --padding <F>` | `3.0` | Padding beyond the molecule in Å |

Output: `{"origin": [...], "spacing": ..., "dims": [...], "n_values": ..., "values": [...]}`.

```bash
$ sci-form esp "$EL" "$COORDS" --spacing 1.0 --padding 2.0
{
  "origin": [-3.6034250259399414, -3.4132779836654663, -3.206077814102173],
  "spacing": 1,
  "dims": [9, 8, 8],
  "n_values": 576,
  "values": [0.011519911150848217, 0.012910413639238533, 0.013471356599164, "..."]
}
```

---

### `ani`

Run an ANI ML-potential calculation. Energy in **Hartree**, forces in **Hartree/Å**.

```
sci-form ani <ELEMENTS_JSON> <COORDS_JSON>
```

The CLI builds a small internally generated network whose input dimension matches the real AEV
length (4032 with the default parameters), so the reported energies are **not** physically
meaningful ANI-2x energies — see
[Machine Learning Potentials](../algorithm/machine-learning-potentials.md).

```bash
$ sci-form ani "$EL" "$COORDS"
{
  "energy": 4.403698521286101,
  "forces": [[-1.5078480065994087, 0.29679825434186385, 0.29679355927487694], "..."],
  "species": [6, 6, 8, "..."],
  "atomic_energies": [0.5741402307662766, 0.4447433910543191, 0.698298154587427, "..."],
  "aevs": [[0.40703985420692895, 0.7576920324602711, "..."], "..."]
}
```

---

### `aevs-tm`

Compute ANI-TM atomic environment vectors (extended element set).

```
sci-form aevs-tm <ELEMENTS_JSON> <COORDS_JSON>
```

```bash
$ sci-form aevs-tm "$EL" "$COORDS"
{
  "aev_length": 176256,
  "aevs": [[0.6041681147195205, 1.5593699558625387, "..."], "..."],
  "n_atoms": 9
}
```

---

## Force fields and geometry optimization

### `uff`

Compute the UFF force-field energy in **kcal/mol**.

```
sci-form uff <SMILES> <COORDS_JSON>
```

```bash
$ sci-form uff "CCO" "$COORDS"
{
  "energy": 6.993048,
  "unit": "kcal/mol"
}
```

---

### `optimize-framework`

Optimize a geometry with the shared BFGS engine (UFF backend). Energy in **eV**, forces in
**eV/Å**.

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
| `--fixed-atoms <JSON>` | *absent* | JSON array of fixed atom indices |
| `--a`, `--b`, `--c` | *absent* | Lattice lengths in Å (all three enable periodic optimization) |
| `--alpha`, `--beta`, `--gamma` | `90.0` | Lattice angles in degrees |

```bash
$ sci-form optimize-framework "CCO" "$COORDS"
{
  "positions": [[-0.9213744617887727, 0.06510163058925493, 0.07289555794502865], "..."],
  "energy": 0.06718227655173616,
  "forces": [[0.009298812066820062, 0.0030568574379609397, -0.009864256301979856], "..."],
  "max_force": 0.04905923800286457,
  "n_iterations": 21,
  "converged": true,
  "energy_history": [0.3032479101885856, 0.2589600738340934, "..."]
}
```

---

### `uhf`, `rohf`

Open-shell Hartree-Fock. The energies are in **Hartree**, the gap in **eV**.

```
sci-form uhf  <ELEMENTS_JSON> <COORDS_JSON> --charge 0 --multiplicity 2
sci-form rohf <ELEMENTS_JSON> <COORDS_JSON> --charge 0 --multiplicity 2
```

### `pm3-gradient`, `xtb-gradient`

Analytical energy gradients, in **eV/Å** per atom.

```
sci-form pm3-gradient <ELEMENTS_JSON> <COORDS_JSON>
sci-form xtb-gradient <ELEMENTS_JSON> <COORDS_JSON>
```

### `uff-breakdown`, `mmff94-breakdown`

The force-field energy split by term type, in **kcal/mol**. The values sum to
the corresponding total energy.

```
sci-form uff-breakdown    <SMILES> <COORDS_JSON>
sci-form mmff94-breakdown <SMILES> <COORDS_JSON>
```

### `tpsa`

Topological polar surface area, in **Å²**.

```
sci-form tpsa <SMILES>
```

### `ensemble`

Ensemble property prediction: the consensus LogP with its uncertainty, the
solubility, the TPSA, the pKa estimates, the Veber rules, the BBB permeability
and an overall confidence score.

```
sci-form ensemble <SMILES>
```

### `parse-cif`, `write-cif`

Read a CIF file into a JSON crystal structure and write it back. The path `-`
reads standard input.

```
sci-form parse-cif structure.cif      # CIF  -> JSON
sci-form write-cif structure.json     # JSON -> CIF
```

## Spectroscopy

### `ir-spectrum-broadened`

Compute a broadened IR spectrum; the vibrational analysis is run internally. Wavenumbers in
**cm⁻¹**.

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
| `--step-size <F>` | *absent* | Finite-difference step in Å (omit for the engine default) |
| `--smiles <S>` | *absent* | SMILES string (required for `method=uff`) |

```bash
$ sci-form ir-spectrum-broadened "$EL" "$COORDS" --n-points 3
{
  "wavenumbers": [400.0, 2200.0, 4000.0],
  "intensities": [0.3158458937328895, 0.007591114452204938, 0.0011914799455048945],
  "peaks": [
    { "frequency_cm1": 567.1107379279891, "ir_intensity": 476.9858310295087, "mode_index": 14, "assignment": "skeletal mode" },
    { "frequency_cm1": 1009.3246264224712, "ir_intensity": 476.4248130923807, "mode_index": 17, "assignment": "C-O stretch" },
    "..."
  ],
  "gamma": 15.0,
  "notes": [
    "IR spectrum generated with Lorentzian broadening (γ = 15 cm⁻¹).",
    "Vibrational analysis method: EHT."
  ]
}
```

::: warning
The EHT Hessian is a one-electron model with no electron-electron repulsion, so the surface has
no minimum and imaginary frequencies are expected. Use `--method pm3` or `--method xtb`.
:::

---

### `ensemble-j-couplings`

Predict J-couplings averaged over a conformer ensemble with Boltzmann weighting. Couplings in
**Hz**.

```
sci-form ensemble-j-couplings [OPTIONS] <SMILES> <CONFORMER_COORDS> <ENERGIES_KCAL>
```

| Option | Default | Description |
|--------|---------|-------------|
| `--temperature-k <F>` | `298.15` | Temperature in K |

`CONFORMER_COORDS` is a JSON array of conformers (each a flat xyz array) and `ENERGIES_KCAL` a
JSON array of the matching energies in kcal/mol.

```bash
$ sci-form ensemble-j-couplings "CC" "[$COORDS]" "[0.0]"
[
  { "h1_index": 2, "h2_index": 3, "j_hz": 12.0, "n_bonds": 2, "coupling_type": "geminal_2J" },
  { "h1_index": 2, "h2_index": 4, "j_hz": 12.0, "n_bonds": 2, "coupling_type": "geminal_2J" },
  { "h1_index": 2, "h2_index": 5, "j_hz": 7.0, "n_bonds": 3, "coupling_type": "vicinal_3J_H-C-C-H" },
  "..."
]
```

---

### `assign-peaks`

Assign IR peaks to functional groups using the group-frequency tables.

```
sci-form assign-peaks [OPTIONS] <FREQUENCIES> <INTENSITIES> <ELEMENTS>
```

| Option | Default | Description |
|--------|---------|-------------|
| `--mode-displacements <JSON>` | *absent* | JSON array of modes, each an array of `[x, y, z]` triples |

`FREQUENCIES` are vibrational wavenumbers in cm⁻¹ and `INTENSITIES` the matching IR
intensities.

```bash
$ sci-form assign-peaks "[1000.0,3000.0]" "[10.0,20.0]" "[6,6,8,1,1,1,1,1,1]"
{
  "assignments": [
    {
      "frequency": 1000.0,
      "intensity": 10.0,
      "assignments": [{ "group": "C-O (alcohol)", "vibration_type": "stretch", "match_quality": 0.019607843137254943 }],
      "confidence": 0.019607843137254943
    },
    {
      "frequency": 3000.0,
      "intensity": 20.0,
      "assignments": [{ "group": "O-H (carboxylic)", "vibration_type": "stretch", "match_quality": 0.7506234413965087 }],
      "confidence": 0.7506234413965087
    }
  ],
  "functional_groups": ["C-O (alcohol)", "O-H (carboxylic)"]
}
```

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

```bash
$ sci-form whim "$EL" "$COORDS"
{
  "eigenvalues": [1.033237341984902, 0.25228372095377516, 0.08062139400291679],
  "lambda": [1.033237341984902, 0.25228372095377516, 0.08062139400291679],
  "nu": [0.7563174226339673, 0.18466867761255804, 0.05901389975347466],
  "gamma": [-0.24577459077723623, -0.6156797284172242, -0.12340579807615487],
  "eta": [1.6438974110884448, 2.849790381993459, 9.148746251418517],
  "total_spread": 1.3661424569415939,
  "anisotropy": 0.6767367754460833,
  "directional": 0.008242361511986561,
  "kurtosis": 13.64243404450042
}
```

---

### `rdf`

Compute radial distribution function (RDF) descriptors.

```
sci-form rdf [OPTIONS] <ELEMENTS_JSON> <COORDS_JSON>
```

| Option | Default | Description |
|--------|---------|-------------|
| `--charges <JSON>` | *absent* | Per-atom charges; when absent Gasteiger charges are derived internally |
| `--beta <F>` | `100.0` | Gaussian smoothing factor in Å⁻² |
| `--r-max <F>` | `12.0` | Maximum radius in Å |
| `--n-points <N>` | `128` | Number of sampling points |

```bash
$ sci-form rdf "$EL" "$COORDS" --n-points 3
{
  "radii": [2.0, 6.0, 10.0],
  "rdf_unweighted": [0.057364339808192256, 1.0595680824025278e-238, 0.0],
  "rdf_mass": [0.7803947670122672, 1.0765889840782421e-238, 0.0],
  "rdf_electronegativity": [0.37220077878129904, 5.128309518828235e-238, 0.0],
  "rdf_charge": [-0.00041614178062904564, 5.5828220761880704e-241, 0.0],
  "beta": 100.0,
  "n_pairs": 36
}
```

---

### `getaway`

Compute GETAWAY 3D molecular descriptors.

```
sci-form getaway [OPTIONS] <ELEMENTS_JSON> <COORDS_JSON>
```

| Option | Default | Description |
|--------|---------|-------------|
| `--bonds <JSON>` | *absent* | JSON array of `[i, j]` bond pairs; when absent bonds are perceived from the geometry |
| `--max-lag <N>` | `8` | Maximum topological distance for the autocorrelations |

```bash
$ sci-form getaway "$EL" "$COORDS"
{
  "leverages": [0.07676409836167276, 0.06192187725391481, "..."],
  "hat_autocorrelation": [1.402739873882024, 3.5290159966929386, "..."],
  "r_autocorrelation": [1.269920321983992, 1.8085469206448157, "..."],
  "hat_total": 4.121156515642049,
  "r_total": 2.030294077854183,
  "h_information": 2.054885008317817,
  "h_max": 0.48705220350105116,
  "h_mean": 0.3333333333333333,
  "max_lag": 8
}
```

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

```bash
$ sci-form train-random-forest "[[1,2],[3,4],[5,6],[7,8]]" "[1.0,2.0,3.0,4.0]"
{
  "model_type": "random_forest",
  "n_estimators": 100,
  "n_samples": 4,
  "predictions": [2.505, 2.505, 2.505, 2.505],
  "train_rmse": 1.118045169033881
}
```

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

```bash
$ sci-form train-gradient-boosting "[[1,2],[3,4],[5,6],[7,8]]" "[1.0,2.0,3.0,4.0]"
{
  "model_type": "gradient_boosting",
  "n_estimators": 100,
  "n_samples": 4,
  "predictions": [2.510925865628069, 2.510925865628069, 2.510925865628069, 2.510925865628069],
  "train_rmse": 1.1180873733924923
}
```

---

## Materials and periodicity

### `cell`

Create a unit cell and print its parameters. Lengths in **Å**, angles in **degrees**, volume in
**Å³**.

```
sci-form cell [OPTIONS] --a <A> --b <B> --c <C>
```

| Option | Default | Description |
|--------|---------|-------------|
| `--a`, `--b`, `--c` | (required) | Lattice parameters in Å |
| `--alpha`, `--beta`, `--gamma` | `90.0` | Lattice angles in degrees |

```bash
$ sci-form cell --a 5.43 --b 5.43 --c 5.43
{
  "a": 5.43,
  "b": 5.43,
  "c": 5.43,
  "alpha": 90.0,
  "beta": 90.0,
  "gamma": 90.0,
  "volume": 160.103,
  "lattice": [[5.43, 0.0, 0.0], [0.0, 5.43, 0.0], [0.0, 0.0, 5.43]]
}
```

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
$ sci-form assemble --topology pcu --metal 30 --geometry octahedral --a 10.0
{
  "cell": { "lattice": [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]] },
  "atoms": [
    { "element": 30, "frac_coords": [0.0, 0.0, 0.0] },
    { "element": 6, "frac_coords": [0.43, 0.0, 0.0] },
    "..."
  ],
  "labels": ["node_0", "linker_0", "linker_0", "..."]
}
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

```bash
$ sci-form space-group --number 14
{
  "number": 14,
  "symbol": "P 21/c",
  "crystal_system": "Monoclinic",
  "lattice_type": "P",
  "point_group": "2/m",
  "operations": [{ "rotation": [[1, 0, 0], [0, -1, 0], [0, 0, 1]], "translation": [0.0, -0.5, -0.5] }, "..."]
}
```

---

### `all-space-groups`

List all 230 space groups.

```
sci-form all-space-groups
```

```bash
$ sci-form all-space-groups
[[1, "P 1"], [2, "P -1"], [3, "P 2"], "..."]
```

---

### `periodic-molecule`

Build a periodic molecular graph from fractional coordinates and a unit cell.
`FRAC_COORDS` is a JSON array of `[x, y, z]` triples.

```
sci-form periodic-molecule [OPTIONS] --a <A> --b <B> --c <C> <ELEMENTS> <FRAC_COORDS>
```

| Option | Default | Description |
|--------|---------|-------------|
| `--a`, `--b`, `--c` | (required) | Lattice lengths in Å |
| `--alpha`, `--beta`, `--gamma` | `90.0` | Lattice angles in degrees |
| `--bond-tolerance <F>` | *absent* | Extra bonding tolerance in Å (omit for the default 0.3) |

```bash
$ sci-form periodic-molecule "$EL" "[[0.1,0.2,0.3],[0.1,0.2,0.3],[0.1,0.2,0.3],[0.1,0.2,0.3],[0.1,0.2,0.3],[0.1,0.2,0.3],[0.1,0.2,0.3],[0.1,0.2,0.3],[0.1,0.2,0.3]]" --a 10 --b 10 --c 10
{
  "elements": [6, 6, 8, "..."],
  "frac_coords": [[0.1, 0.2, 0.3], "..."],
  "cart_coords": [[1.0000000000000004, 2.0, 3.0], "..."],
  "bonds": [],
  "cell": { "lattice": [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]] },
  "n_atoms": 9
}
```

---

### `hapticity`

Detect haptic (η) metal–ring interactions and metallocene motifs.

```
sci-form hapticity [OPTIONS] <ELEMENTS_JSON> <COORDS_JSON>
```

| Option | Default | Description |
|--------|---------|-------------|
| `--bonds <JSON>` | *absent* | JSON array of `[i, j]` bond pairs; when absent bonds are perceived from the geometry |

```bash
$ sci-form hapticity "$EL" "$COORDS"
{
  "interactions": [],
  "is_metallocene": false,
  "n_interactions": 0
}
```

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

```bash
$ sci-form simplified-neb-path "CCO" "$COORDS" "$COORDS" --n-images 2 --n-iter 2
{
  "images": [
    { "index": 0, "coords": [-0.8930637240409851, "..."], "potential_energy_kcal_mol": 6.993048432903878 },
    { "index": 1, "coords": [-0.8930637240409851, "..."], "potential_energy_kcal_mol": 6.993048432903878 }
  ],
  "notes": [
    "Simplified NEB (uff) with spring-coupled gradient relaxation on 0 internal images.",
    "Low-cost exploratory path; not a full climbing-image / tangent-projected NEB."
  ]
}
```

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

```bash
$ sci-form neb-energy "CCO" "$COORDS"
{
  "energy_kcal_mol": 6.993048432903878,
  "method": "uff"
}

$ sci-form neb-gradient "CCO" "$COORDS"
{
  "energy_kcal_mol": 6.993048432903878,
  "gradient_kcal_mol_ang": [11.829034662841224, 5.8505849421547325, "..."],
  "method": "uff"
}
```

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
| `--methods <JSON>` | *absent* | JSON array of backend names; when absent every available backend is used |

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
| `--smirks <S>` | *absent* | Optional SMIRKS reaction string |
| `--seed <N>` | `42` | Random seed |

Both positional arguments are JSON arrays of SMILES, e.g. `'["[Cl-]", "CBr"]'` and
`'["ClC", "[Br-]"]'`.

---

## Solvation, dispersion and electrochemistry

### `alpb` *(feature `experimental-alpb`)*

Compute ALPB implicit-solvation energy.

```
sci-form alpb [OPTIONS] <ELEMENTS_JSON> <COORDS_JSON>
```

| Option | Default | Description |
|--------|---------|-------------|
| `-c, --charges <JSON>` | *absent* | Per-atom charges |
| `-d, --dielectric <F>` | `78.5` | Solvent dielectric constant |

```bash
$ sci-form alpb "$EL" "$COORDS" --charges "$CH"
{"alpb_factor": 0.9871227020284441, "born_radii": [2.742389270922675, "..."], "electrostatic_energy": -0.6413647814752553, "nonpolar_energy": 0.9922596041654154, "total_energy": 0.35089482269016004}
```

---

### `d4` *(feature `experimental-d4`)*

Compute the DFT-D4 dispersion correction.

```
sci-form d4 [OPTIONS] <ELEMENTS_JSON> <COORDS_JSON>
```

| Option | Default | Description |
|--------|---------|-------------|
| `--three-body` | off | Include the three-body ATM term |

```bash
$ sci-form d4 "$EL" "$COORDS"
{"coordination_numbers": [3.7422894434484193, "..."], "e2_body": -0.006217973865799453, "e3_body": 0.0, "total_energy": -0.006217973865799453, "total_kcal_mol": -3.901834562553949}
```

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

```bash
$ sci-form cpm "$EL" "$COORDS" --mu-ev=-4.44
{"charges": [-1.8305450096916451, "..."], "converged": true, "electrostatic_energy": -20.21137656783163, "grand_potential": -87.85108248313057, "iterations": 1, "mu_ev": -4.44, "total_charge": -15.2341679989412}
```

---

## Electrical double layer *(feature `alpha-edl`)*

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

```bash
$ sci-form htst-rate --activation-free-energy 0.5 --reaction-free-energy 0.1
{"equilibrium_constant": 0.020400936963103865, "forward_rate_s_inv": 21953.948275021547, "reverse_rate_s_inv": 1076124.5091206536, "step_id": "cli", "temperature_k": 298.15}
```

### `solve-microkinetic-network`

Integrate a microkinetic network.

```
sci-form solve-microkinetic-network [OPTIONS]
```

| Option | Default | Description |
|--------|---------|-------------|
| `--network <JSON>` | *absent* (built-in demo network) | `MicrokineticNetworkConfig` as JSON (file path or inline) |
| `--n-steps <N>` | `5` | Number of integration frames |

### `solve-microkinetic-steady-state`

Solve the microkinetic steady state.

```
sci-form solve-microkinetic-steady-state [OPTIONS]
```

| Option | Default | Description |
|--------|---------|-------------|
| `--network <JSON>` | *absent* (built-in demo network) | `MicrokineticNetworkConfig` as JSON (file path or inline) |
| `--n-steps <N>` | `5` | Number of integration frames |
| `--tol <F>` | `1e-10` | Convergence tolerance on the population change |

### `extract-kinetics-diagnostics`

Extract kinetics diagnostics from a microkinetic trace.

```
sci-form extract-kinetics-diagnostics [OPTIONS]
```

| Option | Default | Description |
|--------|---------|-------------|
| `--network <JSON>` | *absent* (built-in demo network) | `MicrokineticNetworkConfig` as JSON (file path or inline) |

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

```bash
$ sci-form kmesh
{"grid": [4, 4, 1], "n_points": 16, "weights_sum": 1.0}
```

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
| `--fermi-energy-ev <F>` | *absent* | Fermi energy in eV |

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
| `--sketch-size <N>` | *absent* (engine default) | Nyström sketch rank |

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

`OVERLAP_BLOCKS` is optional; when absent a non-overlapping operator is assembled.

### `validate-electron-count`

Validate the electron count of a periodic calculation.

```
sci-form validate-electron-count --n-electrons <N> --n-bands <N> --n-kpoints <N>
```

All three options are required.

```bash
$ sci-form validate-electron-count --n-electrons 2 --n-bands 4 --n-kpoints 8
{"n_bands": 4, "n_electrons": 2, "n_kpoints": 8, "valid": true}
```

---

## Chart utilities

### `chart-to-json`, `chart-from-json`, `validate-chart-schema`

Serialize, deserialize and validate chart payloads.

```
sci-form chart-to-json        [JSON]
sci-form chart-from-json      [JSON]
sci-form validate-chart-schema [JSON]
```

The `JSON` positional argument is optional; it accepts an inline JSON string. The payload
schema is `{title, series: [{series_id, label, x, y, x_unit, y_unit}]}`.

```bash
$ sci-form chart-to-json '{"title":"t","series":[{"series_id":"s","label":"S","x":[0,1],"y":[0,1],"x_unit":"V","y_unit":"A"}]}'
{
  "data": {
    "series": [{ "label": "S", "series_id": "s", "x": [0.0, 1.0], "x_unit": "V", "y": [0.0, 1.0], "y_unit": "A" }],
    "title": "t"
  },
  "format": "json",
  "version": "1.0"
}

$ sci-form chart-from-json '{"title":"t","series":[{"series_id":"s","label":"S","x":[0,1],"y":[0,1],"x_unit":"V","y_unit":"A"}]}'
{
  "title": "t",
  "series": [{ "series_id": "s", "label": "S", "x": [0.0, 1.0], "x_unit": "V", "y": [0.0, 1.0], "y_unit": "A" }]
}

$ sci-form validate-chart-schema '{"title":"t","series":[{"series_id":"s","label":"S","x":[0,1],"y":[0,1],"x_unit":"V","y_unit":"A"}]}'
{
  "valid": true
}
```

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

```bash
$ sci-form dos-chart --energies '[-1,0,1]' --dos-values '[0.1,0.5,0.1]'
{"series": [{"label": "Density of states", "series_id": "dos", "x": [-1.0, 0.0, 1.0], "x_unit": "eV", "y": [0.1, 0.5, 0.1], "y_unit": "a.u."}], "series_count": 1, "title": "Periodic DOS"}
```

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
| `--network <JSON>` | *absent* (built-in demo network) | `MicrokineticNetworkConfig` JSON (file path or inline) |
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
| `fermi-surface-chart` | `--fermi-energy-ev <F>` (*absent*), `--n-kx <N>` (default `24`), `--n-ky <N>` (default `24`) |

---

## Output formats

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

## Exit codes

| Code | Meaning |
|------|---------|
| `0` | Success |
| `1` | Runtime failure: invalid SMILES, embedding failure, missing input file, I/O error, or a failed calculation |
| `2` | Usage error reported by the argument parser (unknown subcommand, missing or invalid flag) |

---

## Pipeline examples

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

---

## See also

- [CLI guide](../guide/cli.md)
- [Python reference](./python.md)
- [TypeScript reference](./typescript.md)
- [Rust reference](./rust.md)
- [Complete supported reference](./complete-supported-reference.md)
