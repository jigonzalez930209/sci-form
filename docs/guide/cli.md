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

### `info` — Version

```bash
$ sci-form info
sci-form 0.1.0
```

## Exit Codes

| Code | Meaning |
|------|---------|
| 0 | Success |
| 1 | Conformer generation failed |

## Pipeline Examples

```bash
# Generate SDF file from SMILES list
sci-form batch -i ligands.smi -o conformers.sdf -f sdf -t 8

# Pipe through jq to extract coordinates
sci-form embed "CCO" | jq '.coords | [range(0; length; 3)] | map(. as $i | [input[$i], input[$i+1], input[$i+2]])'

# Count embed failures in a batch
sci-form batch -i molecules.smi | jq '[.[] | select(.error != null)] | length'
```
