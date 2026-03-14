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

### `energy`

Evaluate UFF force field energy.

```
sci-form energy <SMILES> [OPTIONS]
```

**Options:**

| Option | Description |
|--------|-------------|
| `--coords <FILE>` | JSON file with flat coordinate array |
| `--from-smiles` | Generate conformer first, then evaluate energy |

```bash
# Auto-generate + evaluate
sci-form energy "CCO" --from-smiles

# From existing coords
sci-form energy "CCO" --coords coords.json
```

Output:
```json
{"energy": 12.345, "unit": "kcal/mol"}
```

---

### `info`

Display build and feature information.

```bash
sci-form info
```

Output includes: version, compiled features (`parallel`, etc.), supported platforms.

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

