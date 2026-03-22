# CLI Reference

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

---

## Commands

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
  "coords": [0.123, 0.456, 0.789, ...],
  "elements": [6, 6, 8, 1, 1, 1, 1, 1, 1],
  "bonds": [[0, 1, "SINGLE"], [1, 2, "SINGLE"], ...],
  "error": null,
  "time_ms": 2.34
}
```

---

### `batch`

Process multiple SMILES from a file.

```
sci-form batch <FILE> [OPTIONS]
```

**Options:**

| Option | Default | Description |
|--------|---------|-------------|
| `-s, --seed <N>` | `42` | RNG seed |
| `-f, --format <FMT>` | `json` | Output format: `json`, `xyz`, `sdf` |
| `-t, --threads <N>` | `0` | Number of threads (`0` = auto) |

**Input file format:** One SMILES per line (name optional):
```
CCO ethanol
c1ccccc1 benzene
CC(=O)O acetic_acid
```

**Examples:**

```bash
sci-form batch molecules.smi
sci-form batch molecules.smi -f xyz -t 8
cat molecules.smi | sci-form batch /dev/stdin
sci-form batch molecules.smi > conformers.json
```

---

### `parse`

Parse SMILES and return molecular graph info (no 3D generation).

```
sci-form parse <SMILES>
```

```bash
sci-form parse "c1ccccc1"
```

Output:
```json
{
  "num_atoms": 12,
  "num_bonds": 12,
  "elements": [6, 6, 6, 6, 6, 6, 1, 1, 1, 1, 1, 1],
  "bonds": [[0, 1, "AROMATIC"], ...]
}
```

---

### `charges`

Compute Gasteiger-Marsili partial charges.

```
sci-form charges <SMILES>
```

```bash
sci-form charges "CCO"
```

Output:
```json
{
  "charges": [-0.387, -0.042, -0.228, 0.123, ...],
  "iterations": 6,
  "total_charge": 0.0
}
```

---

### `eht`

Run Extended Hückel Theory calculation.

```
sci-form eht <ELEMENTS_JSON> <COORDS_JSON> [OPTIONS]
```

**Options:**

| Option | Default | Description |
|--------|---------|-------------|
| `-k` | `1.75` | Wolfsberg-Helmholtz constant (0 = default 1.75) |

```bash
sci-form eht "[6,1,1,1,1]" "[0,0,0,1,0,0,0,1,0,0,0,1]" --k 1.75
```

---

### `population`

Mulliken & Löwdin population analysis.

```
sci-form population <ELEMENTS_JSON> <COORDS_JSON>
```

---

### `dipole`

Compute molecular dipole moment (Debye).

```
sci-form dipole <ELEMENTS_JSON> <COORDS_JSON>
```

---

### `esp`

Compute electrostatic potential on a 3D grid.

```
sci-form esp <ELEMENTS_JSON> <COORDS_JSON> [OPTIONS]
```

| Option | Default | Description |
|--------|---------|-------------|
| `--spacing` | `0.5` | Grid spacing in Å |
| `--padding` | `3.0` | Padding beyond molecule in Å |

---

### `dos`

Compute density of states.

```
sci-form dos <ELEMENTS_JSON> <COORDS_JSON> [OPTIONS]
```

| Option | Default | Description |
|--------|---------|-------------|
| `--sigma` | `0.3` | Broadening width in eV |
| `--e-min` | `-30.0` | Energy window min in eV |
| `--e-max` | `5.0` | Energy window max in eV |
| `-n` | `500` | Grid points |

---

### `rmsd`

Compute RMSD between two coordinate sets (Kabsch alignment).

```
sci-form rmsd <COORDS_JSON> <REFERENCE_JSON>
```

---

### `uff`

Compute UFF force field energy.

```
sci-form uff <SMILES> <COORDS_JSON>
```

```bash
sci-form uff "CCO" "[...]"
```

---

### `sasa`

Compute solvent-accessible surface area.

```
sci-form sasa <ELEMENTS_JSON> <COORDS_JSON> [OPTIONS]
```

| Option | Default | Description |
|--------|---------|-------------|
| `--probe-radius` | `1.4` | Probe radius in Å |

---

### `ani`

Run ANI ML potential calculation.

```
sci-form ani <ELEMENTS_JSON> <COORDS_JSON>
```

---

### `hf3c`

Run HF-3c quantum calculation.

```
sci-form hf3c <ELEMENTS_JSON> <COORDS_JSON>
```

---

### `cell`

Create a unit cell and show parameters.

```
sci-form cell [OPTIONS]
```

| Option | Default | Description |
|--------|---------|-------------|
| `--a`, `--b`, `--c` | (required) | Lattice parameters in Å |
| `--alpha`, `--beta`, `--gamma` | `90.0` | Lattice angles in degrees |

```bash
sci-form cell --a 5.43 --b 5.43 --c 5.43
```

---

### `assemble`

Assemble a framework crystal structure.

```
sci-form assemble [OPTIONS]
```

| Option | Default | Description |
|--------|---------|-------------|
| `--topology` | `pcu` | Framework topology |
| `--metal` | `30` | Metal atomic number |
| `--geometry` | `octahedral` | Coordination geometry |
| `--a` | `10.0` | Lattice parameter Å |
| `--supercell` | `1` | Replication factor |

```bash
sci-form assemble --topology pcu --metal 30 --geometry octahedral --a 10.0
```

---

### `info`

Display build and feature information.

```bash
sci-form info
```

Output includes: version, compiled features (`parallel`, etc.), supported platforms.

---

### `stereo`

Analyze stereochemistry: CIP priorities, R/S, E/Z.

```
sci-form stereo <SMILES> [OPTIONS]
```

**Options:**

| Option | Description |
|--------|-------------|
| `--coords <JSON>` | Flat coordinate array as JSON string (enables 3D R/S and E/Z assignment) |

```bash
sci-form stereo "C(F)(Cl)(Br)I"
sci-form stereo "C/C=C/C" --coords "[$(sci-form embed 'C/C=C/C' | jq '.coords')]"
```

Output:
```json
{
  "stereocenters": [
    { "atom_index": 0, "element": 6, "priorities": [1,2,3,4], "configuration": "R" }
  ],
  "double_bonds": [],
  "n_stereocenters": 1,
  "n_double_bonds": 0
}
```

---

### `solvation`

Compute GB/SA solvation energy.

```
sci-form solvation <ELEMENTS_JSON> <COORDS_JSON> [OPTIONS]
```

**Arguments:**
- `ELEMENTS_JSON` — Atomic numbers as JSON array, e.g. `"[8,1,1]"`
- `COORDS_JSON` — Flat coordinates as JSON array (Å)

**Options:**

| Option | Default | Description |
|--------|---------|-------------|
| `--charges <JSON>` | Gasteiger auto | Per-atom charges as JSON array |
| `--probe-radius <F>` | `1.4` | Probe radius in Å |
| `--solvent-dielectric <F>` | `78.5` | Solvent dielectric constant |
| `--nonpolar-only` | — | Skip GB, return only non-polar SASA term |

```bash
# Using SMILES pipeline
RESULT=$(sci-form embed "CCO")
ELEM=$(echo "$RESULT" | jq '.elements')
COORDS=$(echo "$RESULT" | jq '.coords')
sci-form solvation "$ELEM" "$COORDS"
```

Output:
```json
{
  "nonpolar_energy_kcal_mol": 0.72,
  "electrostatic_energy_kcal_mol": -3.45,
  "total_energy_kcal_mol": -2.73,
  "total_sasa": 187.4
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
    { "atoms": [0,1,2,3,4,5], "size": 6, "is_aromatic": true },
    { "atoms": [4,5,6,7,8,9], "size": 6, "is_aromatic": true }
  ],
  "ring_size_histogram": { "6": 2 }
}
```

---

### `ecfp`

Compute an ECFP fingerprint.

```
sci-form ecfp <SMILES> [OPTIONS]
```

**Options:**

| Option | Default | Description |
|--------|---------|-------------|
| `--radius <N>` | `2` | Morgan radius (ECFP4 = 2, ECFP6 = 3) |
| `--n-bits <N>` | `2048` | Bit-vector size |

```bash
sci-form ecfp "c1ccccc1"
sci-form ecfp "c1ccccc1" --radius 3 --n-bits 1024
```

Output:
```json
{
  "n_bits": 2048,
  "on_bits": [42, 317, 891, ...],
  "radius": 2,
  "density": 0.0137
}
```

---

### `tanimoto`

Compute Tanimoto similarity between two molecules.

```
sci-form tanimoto <SMILES1> <SMILES2> [OPTIONS]
```

**Options:**

| Option | Default | Description |
|--------|---------|-------------|
| `--radius <N>` | `2` | ECFP radius |
| `--n-bits <N>` | `2048` | Bit-vector size |

```bash
sci-form tanimoto "c1ccccc1" "Cc1ccccc1"
```

Output:
```json
{ "tanimoto": 0.583, "smiles1": "c1ccccc1", "smiles2": "Cc1ccccc1" }
```

---

## Output Formats

### JSON

Default format. Returns the full `ConformerResult` as compact JSON.

### XYZ

Standard XYZ format for molecular viewers (Avogadro, VESTA, VMD):

```
9
CCO seed=42
C    0.123   0.456   0.789
O    1.234   0.567   0.890
...
```

### SDF

MDL Structure Data File format. Compatible with RDKit, OpenBabel, MarvinSuite.

---

## Exit Codes

| Code | Meaning |
|------|---------|
| `0` | Success |
| `1` | Invalid SMILES or embedding failure |
| `2` | File not found or I/O error |

---

## Pipeline Examples

### Filter successful results

```bash
sci-form batch molecules.smi | jq -c 'select(.error == null)'
```

### Get batch statistics

```bash
sci-form batch molecules.smi | jq -s '{
  total: length,
  success: [.[] | select(.error == null)] | length,
  avg_time_ms: ([.[] | .time_ms] | add / length)
}'
```

### Convert batch to XYZ

```bash
sci-form batch molecules.smi -f xyz > all_conformers.xyz
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

