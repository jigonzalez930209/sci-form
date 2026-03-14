# CLI Reference

## Installation

### Pre-built Binaries

Download from [GitHub Releases](https://github.com/lestard/sci-form/releases):
- `sci-form-x86_64-unknown-linux-gnu` — Linux x86_64
- `sci-form-x86_64-apple-darwin` — macOS Intel
- `sci-form-aarch64-apple-darwin` — macOS Apple Silicon
- `sci-form-x86_64-pc-windows-msvc.exe` — Windows x86_64
- `sci-form-aarch64-unknown-linux-gnu` — Linux ARM64

### Build from Source

```bash
cargo install sci-form-cli
```

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
# JSON output (default)
sci-form embed "CCO"

# XYZ format
sci-form embed "c1ccccc1" -f xyz

# SDF format with custom seed
sci-form embed "CC(=O)O" -f sdf -s 123
```

**JSON output structure:**

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

**Input file format:** One SMILES per line, optionally followed by a name:

```
CCO ethanol
c1ccccc1 benzene
CC(=O)O acetic_acid
```

**Examples:**

```bash
# Process file
sci-form batch molecules.smi

# XYZ with 4 threads
sci-form batch molecules.smi -f xyz -t 4

# Pipe from stdin
cat molecules.smi | sci-form batch /dev/stdin

# Save JSON output
sci-form batch molecules.smi > conformers.json
```

### `parse`

Parse SMILES without generating coordinates.

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
  "bonds": [[0, 1, "AROMATIC"], [1, 2, "AROMATIC"], ...]
}
```

### `info`

Display build and version information.

```bash
sci-form info
```

## Exit Codes

| Code | Meaning |
|------|---------|
| `0` | Success |
| `1` | Invalid SMILES or embedding failure |
| `2` | File not found or I/O error |

## Pipeline Examples

### Extract coordinates with jq

```bash
sci-form embed "CCO" | jq '.coords | [range(0; length; 3)] | map(. as $i | [input.coords[$i], input.coords[$i+1], input.coords[$i+2]])'
```

### Filter successful results

```bash
sci-form batch molecules.smi | jq -c 'select(.error == null)'
```

### Convert batch to XYZ file

```bash
sci-form batch molecules.smi -f xyz > all_conformers.xyz
```

### Get statistics

```bash
sci-form batch molecules.smi | jq -s '{
  total: length,
  success: [.[] | select(.error == null)] | length,
  avg_time_ms: ([.[] | .time_ms] | add / length)
}'
```

### Parallel processing with GNU Parallel

```bash
cat huge_set.smi | parallel -j8 --pipe -N100 sci-form batch /dev/stdin > results.json
```
