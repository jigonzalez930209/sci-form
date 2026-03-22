# CLI

## Installation

Download prebuilt binaries from [GitHub Releases](https://github.com/jigonzalez930209/sci-form/releases), or build from source:

```bash
cargo install sci-form-cli
```

## Commands

### `embed` — Single molecule

Generate a 3D conformer for one SMILES string.

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
  "coords": [0.123, -0.456, ...],
  "elements": [6, 6, 8, 1, 1, 1, 1, 1, 1],
  "bonds": [[0, 1, "SINGLE"], ...],
  "error": null,
  "time_ms": 12.3
}
```

**XYZ output**:
```bash
$ sci-form embed "CCO" -f xyz
9
CCO
 C     0.764735     0.013411    -0.082862
 C    -0.718153     0.048620     0.089507
 O    -1.168543    -0.490543     1.323891
 ...
```

**SDF output**:
```bash
$ sci-form embed "CCO" -f sdf > ethanol.sdf
```

### `batch` — Multiple molecules

Process many molecules from a file or stdin.

```bash
# From file
sci-form batch -i molecules.smi -o output.sdf -f sdf -t 8

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

### `parse` — SMILES analysis

Parse a SMILES string and show molecular structure (no 3D generation).

```bash
$ sci-form parse "c1ccccc1" | jq .
{
  "num_atoms": 12,
  "num_bonds": 12,
  "smiles": "c1ccccc1"
}
```

### `charges` — Gasteiger-Marsili partial charges

```bash
$ sci-form charges "CCO"
{
  "charges": [-0.032, 0.081, -0.329, ...],
  "total_charge": 0.0,
  "iterations": 8
}
```

### `eht` — Extended Hückel Theory

Compute EHT electronic structure from atomic numbers and coordinates.

```bash
sci-form eht "[6,8,1,1,1,1,1,1,1]" "[0.76,0.01,-0.08,-0.72,0.05,0.09,...]" --k 1.75
```

| Option | Default | Description |
|--------|---------|-------------|
| `--k` | 1.75 | Wolfsberg-Helmholz parameter |

### `population` — Mulliken & Löwdin analysis

```bash
sci-form population "[6,8,1,1,1,1,1,1,1]" "[0.76,0.01,-0.08,...]"
# Output: mulliken_charges, lowdin_charges, homo_energy, lumo_energy, gap
```

### `dipole` — Molecular dipole moment

```bash
sci-form dipole "[6,8,1,1,1,1,1,1,1]" "[0.76,0.01,-0.08,...]"
# Output: magnitude (Debye), vector [x, y, z]
```

### `esp` — Electrostatic potential grid

```bash
sci-form esp "[6,8,1,1,1,1,1,1,1]" "[...]" --spacing 0.5 --padding 3.0
```

| Option | Default | Description |
|--------|---------|-------------|
| `--spacing` | 0.5 | Grid spacing (Å) |
| `--padding` | 3.0 | Padding around molecule (Å) |

### `dos` — Density of states

```bash
sci-form dos "[6,8,1,1,1,1,1,1,1]" "[...]" --sigma 0.3 --e-min -30 --e-max 5 -n 500
```

| Option | Default | Description |
|--------|---------|-------------|
| `--sigma` | 0.3 | Gaussian broadening (eV) |
| `--e-min` | -30.0 | Energy range minimum (eV) |
| `--e-max` | 5.0 | Energy range maximum (eV) |
| `-n` | 500 | Number of grid points |

### `rmsd` — RMSD alignment

Compute RMSD between two coordinate sets (flat arrays).

```bash
sci-form rmsd "[0.1,0.2,0.3,...]" "[0.15,0.22,0.31,...]"
# Output: {"rmsd": 0.034, "aligned_coords": [...]}
```

### `uff` — UFF force field energy

```bash
sci-form uff "CCO" "[0.76,0.01,-0.08,...]"
# Output: {"energy": 12.3}  (kcal/mol)
```

### `sasa` — Solvent-accessible surface area

```bash
sci-form sasa "[6,8,1,1,1,1,1,1,1]" "[...]" --probe-radius 1.4
# Output: total_sasa (ų), per_atom_sasa
```

### `ani` — ANI neural network potential

ANI-2x supports H, C, N, O, F, S, Cl.

```bash
sci-form ani "[6,8,1,1,1,1,1,1,1]" "[0.76,0.01,-0.08,...]"
# Output: energy (eV), forces per atom
```

### `hf3c` — Hartree-Fock 3c calculation

Minimal-basis HF with D3 dispersion, gCP counterpoise, and SRB short-range corrections.

```bash
sci-form hf3c "[6,8,1,1,1,1,1,1,1]" "[0.76,0.01,-0.08,...]"
# Output: energy, hf_energy, d3_energy, gcp_energy, srb_energy, orbital_energies, converged
```

### `stereo` — Stereochemistry analysis

Detect R/S stereocenters, E/Z double bonds, atropisomeric axes, and helical chirality.

```bash
$ sci-form stereo "C(F)(Cl)(Br)I"
{
  "stereocenters": [{"atom_index": 0, "configuration": "R", ...}],
  "double_bonds": [],
  "n_stereocenters": 1,
  "n_double_bonds": 0
}

# With 3D coordinates for geometric assignment
sci-form stereo "C/C=C/C" --coords "[...]"
```

### `solvation` — Solvation energy

Non-polar (SASA-based) solvation, optionally with Generalized Born electrostatics.

```bash
# Non-polar only
sci-form solvation "[8,1,1]" "[0,0,0,0.757,0.586,0,-0.757,0.586,0]"

# With GB electrostatic solvation (provide charges)
sci-form solvation "[8,1,1]" "[...]" --charges "[-0.8,0.4,0.4]"
```

### `sssr` — Ring perception

Smallest Set of Smallest Rings (Horton's algorithm).

```bash
$ sci-form sssr "c1ccc2ccccc2c1"   # naphthalene
{
  "rings": [[0,1,2,3,4,9], [4,5,6,7,8,9]],
  "ring_size_histogram": {"6": 2}
}
```

### `ecfp` — ECFP fingerprints

Extended-Connectivity Fingerprints (Morgan algorithm).

```bash
sci-form ecfp "CCO" --radius 2 --n-bits 2048
# Output: {n_bits, on_bits, radius, raw_features}
```

### `tanimoto` — Fingerprint similarity

Compute Tanimoto similarity between two SMILES.

```bash
$ sci-form tanimoto "c1ccccc1" "Cc1ccccc1"
{"tanimoto": 0.714}
```

### `cell` — Create unit cell

```bash
sci-form cell --a 5.43 --b 5.43 --c 5.43 --alpha 90 --beta 90 --gamma 90
# Output: volume, lattice vectors, fractional ↔ Cartesian matrices
```

### `assemble` — Framework assembly

Build MOF/COF crystal structures from topology and building blocks.

```bash
sci-form assemble --topology pcu --metal 30 --geometry octahedral --a 10.0
# Output: num_atoms, elements, coordinates, unit cell, space group
```

| Option | Default | Description |
|--------|---------|-------------|
| `--topology` | required | `pcu`, `dia`, `sql` |
| `--metal` | required | Atomic number of metal node |
| `--geometry` | required | `linear`, `trigonal`, `tetrahedral`, `square_planar`, `octahedral` |
| `--a` | required | Lattice parameter a (Å) |
| `--supercell` | 1 | Supercell expansion factor |

### `info` — Version

```bash
$ sci-form info
sci-form 0.9.1
```

## Exit Codes

| Code | Meaning |
|------|---------|
| 0 | Success |
| 1 | Invalid SMILES or embedding failure |
| 2 | I/O error |

## Pipeline Examples

```bash
# Generate SDF file from SMILES list
sci-form batch -i ligands.smi -o conformers.sdf -f sdf -t 8

# Filter only successful embeddings
sci-form batch -i molecules.smi | jq -c 'select(.error == null)'

# Batch statistics
sci-form batch -i molecules.smi | jq -s '{total:length, success:[.[]|select(.error==null)]|length}'

# Embed → PM3 pipeline (using coordinate extraction with jq)
COORDS=$(sci-form embed "CCO" | jq -c '.coords')
ELEMS=$(sci-form embed "CCO" | jq -c '.elements')

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
