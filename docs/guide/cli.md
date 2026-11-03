# CLI

This is the **guide**: a task-oriented tour of the `sci-form` binary. For the exhaustive list of every subcommand, flag and default, see the [CLI reference](../api/cli.md).

## Installation

Every [release](https://github.com/jigonzalez930209/sci-form/releases/latest) carries the CLI as an installable package for each system:

| System | Package | Install |
|---|---|---|
| Debian / Ubuntu | `.deb` (`amd64`, `arm64`) | `sudo apt install ./sci-form_*.deb` |
| Fedora / RHEL | `.rpm` (`x86_64`, `aarch64`) | `sudo dnf install ./sci-form-*.rpm` |
| Linux (any) | `.tar.gz` (x86_64, aarch64) | unpack and copy `sci-form` to a folder on your `PATH` |
| macOS | `.tar.gz` (Apple Silicon, Intel) | unpack and copy `sci-form` to a folder on your `PATH` |
| Windows | `.zip` (x86_64) | unpack and add the folder to `PATH` |

The `.deb` and `.rpm` also install the bash, zsh and fish completions. A `SHA256SUMS` file lists the checksums. Or build from source:

```bash
cargo install sci-form-cli
```

The binary exposes 102 subcommands. 69 are always available; the other 33 are compiled in only when the corresponding Cargo feature is enabled:

```bash
cargo install sci-form-cli --features alpha-gsm,alpha-edl,alpha-kinetics,alpha-periodic-linear
```

Use `sci-form --help` to list what a given build actually contains.

## Finding your way around

The binary documents itself — you never need to leave the terminal:

```bash
sci-form                        # a short guide: three commands to try, then every command by topic
sci-form commands               # the same list, grouped by topic
sci-form commands quantum       # only one topic (any word of its name)
sci-form examples pm3           # a ready-to-run example for one command
sci-form examples pm3 --run     # …and its output
sci-form examples               # the example of every command
sci-form help pm3               # every option of one command (its last lines repeat the example)
sci-form completions bash       # tab completion: bash, zsh, fish, powershell, elvish
```

Every example in that help is executed by the test-suite, so it works exactly as printed. For
example outputs you can also run in the browser see the [live examples](/guide/live-examples).

## The two calling conventions

Almost every command takes one of two argument shapes:

| Shape | Meaning | Example |
|-------|---------|---------|
| `<SMILES>` | a molecule, parsed from SMILES | `sci-form charges "CCO"` |
| `<ELEMENTS_JSON> <COORDS_JSON>` | a geometry: atomic numbers plus a flat `[x, y, z, ...]` array in Å | `sci-form eht "[8,1,1]" "[0,0,0,0.757,0.586,0,-0.757,0.586,0]"` |

Most commands print JSON on stdout; `parse`, `info` and the chart commands print plain text.

**Negative option values** must be attached with `=`, because the argument parser would otherwise read `-4.44` as a flag:

```bash
sci-form cpm "[6,8,1,1,1,1,1,1,1]" "[...]" --mu-ev=-4.44     # correct
sci-form cpm "[6,8,1,1,1,1,1,1,1]" "[...]" --mu-ev -4.44      # usage error
```

---

## Embedding

### `embed` — one molecule

```bash
sci-form embed "CCO" --format xyz
```

| Option | Default | Description |
|--------|---------|-------------|
| `--seed, -s` | 42 | RNG seed for reproducibility |
| `--format, -f` | json | Output format: `json`, `xyz`, `sdf` |

**JSON output** (default):
```bash
$ sci-form embed "CCO"
{
  "smiles": "CCO",
  "num_atoms": 9,
  "coords": [-0.893, 0.069, 0.064, "..."],
  "elements": [6, 6, 8, 1, 1, 1, 1, 1, 1],
  "bonds": [[0, 1, "SINGLE"], "..."],
  "error": null,
  "time_ms": 9.27
}
```

**XYZ output**:
```bash
$ sci-form embed "CCO" -f xyz
9
CCO
C     -0.893064     0.069264     0.063601
C      0.525542    -0.422006    -0.147153
O      1.389964     0.450482     0.467233
...
```

**SDF output**:
```bash
$ sci-form embed "CCO" -f sdf > ethanol.sdf
```

### `batch` — many molecules

The input path is an option, not a positional argument.

```bash
# From file
sci-form batch -i molecules.smi -o conformers.sdf -f sdf -t 8

# From stdin
cat molecules.smi | sci-form batch -f xyz > output.xyz

# No output file → stdout
sci-form batch -i molecules.smi -f json
```

| Option | Default | Description |
|--------|---------|-------------|
| `--input, -i` | stdin | Input file (one SMILES per line) |
| `--output, -o` | stdout | Output file |
| `--seed, -s` | 42 | RNG seed |
| `--threads, -t` | 0 (auto) | Number of threads |
| `--format, -f` | json | Output format: `json`, `xyz`, `sdf` |

### `embed-diverse` and `filter-diverse`

Generate an ensemble and cluster it by Butina RMSD, then keep only the cluster centroids.

```bash
sci-form embed-diverse "CCO" --n-conformers 20 --rmsd-cutoff 0.8
sci-form filter-diverse '[[0,0,0, 1,0,0], [0,0,0, 1.05,0,0]]' --rmsd-cutoff 1.0
```

### `parse` — SMILES analysis

Plain text, not JSON:

```bash
$ sci-form parse "c1ccccc1"
Atoms: 12, Bonds: 12
    0: elem= 6 hyb=SP2 charge=0 h=0
    ...
```

---

## Electronic structure

```bash
# EHT (Extended Hückel)
sci-form eht "[6,8,1,1,1,1,1,1,1]" "[0.76,0.01,-0.08, ...]" --k 1.75

# Semi-empirical and tight binding
sci-form pm3  "[6,8,1,1,1,1,1,1,1]" "[0.76,0.01,-0.08, ...]"
sci-form xtb  "[6,8,1,1,1,1,1,1,1]" "[0.76,0.01,-0.08, ...]"   # GFN2-xTB (same as gfn2)
sci-form gfn1 "[6,8,1,1,1,1,1,1,1]" "[0.76,0.01,-0.08, ...]"
sci-form gfn2 "[6,8,1,1,1,1,1,1,1]" "[0.76,0.01,-0.08, ...]"
sci-form hf3c "[6,8,1,1,1,1,1,1,1]" "[0.76,0.01,-0.08, ...]"

# CISD excited states from an internal STO-3G RHF calculation
sci-form cisd "[6,8,1,1,1,1,1,1,1]" "[0.76,0.01,-0.08, ...]" --n-states 5

# Gradients and geometry optimization on the EHT surface
sci-form eht-gradient        "[6,8,1,1,1,1,1,1,1]" "[0.76,0.01,-0.08, ...]"
sci-form optimize-geometry-eht "[6,8,1,1,1,1,1,1,1]" "[0.76,0.01,-0.08, ...]"
```

Orbitals and derived quantities:

```bash
sci-form population "[6,8,1,1,1,1,1,1,1]" "[0.76,0.01,-0.08, ...]"   # Mulliken + Löwdin
sci-form dipole     "[6,8,1,1,1,1,1,1,1]" "[0.76,0.01,-0.08, ...]"   # Debye
sci-form esp        "[6,8,1,1,1,1,1,1,1]" "[0.76,0.01,-0.08, ...]" --spacing 0.5 --padding 3.0
sci-form dos        "[6,8,1,1,1,1,1,1,1]" "[0.76,0.01,-0.08, ...]" --sigma 0.3 -n 500
sci-form dos-multimethod "[6,8,1,1,1,1,1,1,1]" "[0.76,0.01,-0.08, ...]" --method gfn2
```

`dos-multimethod` accepts `--method eht|pm3|xtb|gfn1|gfn2|hf3c`; `dos` is EHT-only.

> The `ani` and `aevs-tm` commands exercise the ANI engine, but the CLI builds an internally generated test network, so **the reported energies are not physically meaningful ANI-2x values**. See [Machine Learning Potentials](../algorithm/machine-learning-potentials.md).

---

## Force fields and optimization

```bash
# UFF energy (kcal/mol)
sci-form uff "CCO" "[0.76,0.01,-0.08, ...]"

# BFGS minimization with the UFF backend
sci-form optimize-framework "CCO" "[0.76,0.01,-0.08, ...]" \
  --method bfgs --max-iter 200 --force-tol 0.05

# Constrained / periodic optimization
sci-form optimize-framework "CCO" "[...]" --fixed-atoms "[0,1]" --a 10 --b 10 --c 10
```

---

## Molecular properties

```bash
# Charges
sci-form charges             "CCO"
sci-form charges-configured  "CCO" --max-iter 12 --initial-damping 0.7
sci-form eeq                 "[6,8,1,1,1,1,1,1,1]" "[0.76,0.01,-0.08, ...]" --total-charge 0

# Surface and solvation
sci-form sasa      "[6,8,1,1,1,1,1,1,1]" "[...]" --probe-radius 1.4
sci-form solvation "[6,8,1,1,1,1,1,1,1]" "[...]" --charges "[...]"

# Topology
sci-form sssr "c1ccc2ccccc2c1"
sci-form stereo "C(F)(Cl)(Br)I"
sci-form stereo "C/C=C/C" --coords "[$(sci-form embed 'C/C=C/C' | jq '.coords')]"
sci-form canonical-smiles "C(F)(Cl)(Br)I"

# Fingerprints
sci-form ecfp        "c1ccccc1" --radius 2 --n-bits 2048
sci-form ecfp-batch  '["CCO", "c1ccccc1"]'
sci-form tanimoto    "c1ccccc1" "Cc1ccccc1"

# Alignment
sci-form rmsd "[0.1,0.2,0.3,...]" "[0.15,0.22,0.31,...]"
```

3D descriptors and the bundled regressors:

```bash
sci-form whim    "[6,8,1,1,1,1,1,1,1]" "[...]" --weighting mass
sci-form rdf     "[6,8,1,1,1,1,1,1,1]" "[...]" --r-max 12 --n-points 128
sci-form getaway "[6,8,1,1,1,1,1,1,1]" "[...]" --max-lag 8

sci-form train-random-forest      '[[0.1,0.2],[0.3,0.4]]' '[1.0,1.2]' --n-trees 50
sci-form train-gradient-boosting  '[[0.1,0.2],[0.3,0.4]]' '[1.0,1.2]' --n-estimators 50
```

---

## Spectroscopy

```bash
# Broadened IR spectrum (the vibrational analysis runs internally)
sci-form ir-spectrum-broadened "[6,8,1,1,1,1,1,1,1]" "[...]" \
  --method eht --gamma 15 --broadening lorentzian --n-points 1000

# Peak assignment from frequencies + intensities
sci-form assign-peaks '[3000, 1600, 1400]' '[12.0, 30.0, 8.0]' "[6,8,1,1,1,1,1,1,1]"

# J-couplings averaged over a Boltzmann-weighted ensemble
sci-form ensemble-j-couplings "CCO" '[[0.0,0.0,0.0, ...]]' '[-12.3]' --temperature-k 298.15
```

---

## Materials

```bash
# Unit cell
sci-form cell --a 5.43 --b 5.43 --c 5.43 --alpha 90 --beta 90 --gamma 90

# Framework assembly (any Topology::from_name net)
sci-form assemble --topology pcu --metal 30 --geometry octahedral --a 10.0

# Crystallography
sci-form space-group --number 14
sci-form space-group --symbol "P 21/c"
sci-form all-space-groups

# Periodic molecular graph
sci-form periodic-molecule "[6,6]" "[0.0,0.0,0.0, 0.5,0.5,0.0]" --a 2.5 --b 2.5 --c 2.5

# Hapticity (η-counting of metal-ring contacts)
sci-form hapticity "[6,6,6,6,6,26]" "[...]" --bonds "[[0,1],[1,2],[2,3],[3,4],[4,0],[0,5]]"
```

---

## Reactions

```bash
# Simplified NEB-like path between two geometries
sci-form simplified-neb-path "CCO" "[start coords]" "[end coords]" \
  --n-images 8 --n-iter 50 --method uff

# Single-point energy / gradient with any NEB-capable backend
sci-form neb-energy   "CCO" "[...]" --method mmff94
sci-form neb-gradient "CCO" "[...]" --method gfn2
```

The GSM (growing string method) commands need the `alpha-gsm` feature:

```bash
sci-form gsm-backend-plan  "CCO"
sci-form gsm-compare-backends "CCO" "[...]" --methods '["uff","gfn2"]'
```

The full 3D reaction-dynamics path (`alpha-reaction-dynamics`) is available only from the Rust core; it has no CLI command.

---

## Experimental and alpha surfaces

These commands are feature-gated. Build the binary with the matching feature:

```bash
cargo install sci-form-cli --features experimental-alpb,experimental-d4,experimental-cpm,alpha-edl,alpha-kinetics
```

```bash
# Solvation, dispersion and electrochemistry
sci-form alpb "[6,8,1,1,1,1,1,1,1]" "[...]" --dielectric 78.5
sci-form d4   "[6,8,1,1,1,1,1,1,1]" "[...]" --three-body
sci-form cpm  "[6,8,1,1,1,1,1,1,1]" "[...]" --mu-ev=-4.44

# Electrical double layer
sci-form edl-profile --model gouy-chapman-stern --potential 0.05 --ionic-strength 0.1
sci-form edl-capacitance-scan --v-min=-0.5 --v-max=0.5

# Kinetics
sci-form htst-rate --activation-free-energy 0.35 --reaction-free-energy=-0.12
sci-form solve-microkinetic-network --n-steps 50
```

---

## Exit Codes

| Code | Meaning |
|------|---------|
| 0 | Success |
| 1 | Runtime failure: invalid SMILES, embedding failure, missing input file, I/O error, or a failed calculation |
| 2 | Usage error reported by the argument parser (unknown subcommand, missing or invalid flag) |

## Pipeline Examples

```bash
# Generate SDF file from SMILES list
sci-form batch -i ligands.smi -o conformers.sdf -f sdf -t 8

# Filter only successful embeddings
sci-form batch -i molecules.smi | jq -c 'select(.error == null)'

# Batch statistics
sci-form batch -i molecules.smi | jq -s '{total:length, success:[.[]|select(.error==null)]|length}'

# Embed → charges → EHT pipeline
RESULT=$(sci-form embed "CCO")
COORDS=$(echo "$RESULT" | jq -c '.coords')
ELEMS=$(echo "$RESULT" | jq -c '.elements')
sci-form charges "CCO"
sci-form eht "$ELEMS" "$COORDS" --k 1.75

# Parallel batch with GNU Parallel
cat huge.smi | parallel -j8 --pipe -N100 sci-form batch /dev/stdin > results.json

# Stereochemistry analysis of all molecules
sci-form batch -i molecules.smi | jq -r '.smiles' | while read smi; do
  sci-form stereo "$smi"
done

# Fingerprint similarity matrix
sci-form tanimoto "c1ccccc1" "CCO"
sci-form tanimoto "c1ccccc1" "CC(=O)O"
```
