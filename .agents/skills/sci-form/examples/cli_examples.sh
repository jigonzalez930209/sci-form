#!/usr/bin/env bash
# sci-form CLI Examples — Complete Reference
# Binary: sci-form (install with: cargo install sci-form-cli)
# All commands print JSON to stdout, stats to stderr.

set -euo pipefail
SF="sci-form"

# ─── Data ──────────────────────────────────────────────────────────────────────
ELEMENTS='[8,1,1]'
COORDS='[0,0,0, 0.96,0,0, -0.24,0.93,0]'
BENZ='[6,6,6,6,6,6,1,1,1,1,1,1]'

# ─── 1. Geometry / Embedding ───────────────────────────────────────────────────

# Single conformer — JSON output
$SF embed "C1=CC=CC=C1"
# JSON output
$SF embed "CC(=O)O" --seed 42

# XYZ format
$SF embed "c1ccccc1" --format xyz > benzene.xyz

# SDF format
$SF embed "CC(N)C(=O)O" --format sdf > alanine.sdf

# Parse topology only (no 3D)
$SF parse "CC(=O)N"

# Info / version
$SF info

# ─── 2. Batch Processing (Parallel) ───────────────────────────────────────────

# From file — auto-detects cores with --threads 0
echo -e "O\nCCO\nc1ccccc1\nCC(=O)O" > smiles.txt
$SF batch --input smiles.txt --threads 0 --format json

# From stdin
cat smiles.txt | $SF batch --threads 0 --format xyz

# Save to file
$SF batch --input smiles.txt --output conformers.sdf --format sdf --threads 4

# With custom seed
$SF batch --input smiles.txt --seed 2024 --threads 0

# ─── 3. EHT Calculation ────────────────────────────────────────────────────────

# Default k=1.75
$SF eht "$ELEMENTS" "$COORDS"

# Custom Wolfsberg-Helmholtz k constant
$SF eht "$ELEMENTS" "$COORDS" --k 1.8

# ─── 4. Charges (Gasteiger-Marsili) ───────────────────────────────────────────

$SF charges "CC(=O)O"
$SF charges "c1ccccc1"
$SF charges "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"   # Caffeine

# ─── 5. SASA (Shrake-Rupley) ──────────────────────────────────────────────────

# Default probe radius 1.4 Å
$SF sasa "$ELEMENTS" "$COORDS"

# Custom probe radius
$SF sasa "$ELEMENTS" "$COORDS" --probe-radius 1.6

# ─── 6. Population Analysis ────────────────────────────────────────────────────

# Mulliken + Löwdin analysis from EHT
$SF population "$ELEMENTS" "$COORDS"

# ─── 7. Dipole Moment ─────────────────────────────────────────────────────────

$SF dipole "$ELEMENTS" "$COORDS"

# ─── 8. ESP Grid ──────────────────────────────────────────────────────────────

# Default: spacing=0.5 Å, padding=3.0 Å (prints metadata only, too large for stdout)
$SF esp "$ELEMENTS" "$COORDS"

# Finer grid
$SF esp "$ELEMENTS" "$COORDS" --spacing 0.3 --padding 4.0

# ─── 9. Density of States ─────────────────────────────────────────────────────

# Default: sigma=0.3 eV, e_min=-30 eV, e_max=5 eV, 500 points
$SF dos "$ELEMENTS" "$COORDS"

# Custom window and resolution
$SF dos "$ELEMENTS" "$COORDS" --sigma 0.2 --e-min -20 --e-max 2 --n-points 1000

# ─── 10. RMSD / Kabsch Alignment ──────────────────────────────────────────────

C1_COORDS='[0,0,0, 0,0,1.09, 1.02,0,-0.36, -0.51,0.89,-0.36, -0.51,-0.89,-0.36]'
C2_COORDS='[0.1,0.1,0, 0.1,0.1,1.09, 1.12,0.1,-0.36, -0.41,0.99,-0.36, -0.41,-0.79,-0.36]'
$SF rmsd "$C1_COORDS" "$C2_COORDS"

# ─── 11. Force Fields ─────────────────────────────────────────────────────────

BENZ_COORDS='[-1.21,0.7,0, -1.21,-0.7,0, 0,-1.4,0, 1.21,-0.7,0, 1.21,0.7,0, 0,1.4,0, -2.16,1.25,0, -2.16,-1.25,0, 0,-2.5,0, 2.16,-1.25,0, 2.16,1.25,0, 0,2.5,0]'
$SF uff "c1ccccc1" "$BENZ_COORDS"

# ─── 12. Unit Cell ────────────────────────────────────────────────────────────

# Cubic cell
$SF cell --a 10 --b 10 --c 10

# Monoclinic cell
$SF cell --a 8.5 --b 12.3 --c 7.1 --alpha 90 --beta 102.5 --gamma 90

# ─── 13. Framework Assembly ────────────────────────────────────────────────────

# Default: PCU, Zn(II) octahedral, 10 Å lattice
$SF assemble

# MOF-5 style: Zn(40) octahedral, PCU topology, 25 Å unit cell
$SF assemble --topology pcu --metal 30 --geometry octahedral --a 25.0

# Diamond topology (dia)
$SF assemble --topology dia --metal 30 --geometry tetrahedral --a 20.0

# 2×2×2 supercell
$SF assemble --topology pcu --metal 30 --geometry octahedral --a 10.0 --supercell 2

# ─── 14. ANI Neural Network Potential ─────────────────────────────────────────

# Uses test model weights (real weights require separate integration)
$SF ani "$ELEMENTS" "$COORDS"

# ─── 15. HF-3c Quantum Chemistry ──────────────────────────────────────────────

$SF hf3c "$ELEMENTS" "$COORDS"

# ─── 16. Stereochemistry ──────────────────────────────────────────────────────

# R/S chirality (topology-only, no 3D)
$SF stereo "C(F)(Cl)(Br)I"

# R/S chirality with 3D coords (more accurate)
CHIRAL_COORDS=$($SF embed "C(F)(Cl)(Br)I" --seed 42 | python3 -c "import sys,json; print(json.dumps(json.load(sys.stdin)['coords']))")
$SF stereo "C(F)(Cl)(Br)I" --coords "$CHIRAL_COORDS"

# E/Z double bond configuration with 3D
ALK_COORDS=$($SF embed "CC=CC" --seed 42 | python3 -c "import sys,json; print(json.dumps(json.load(sys.stdin)['coords']))")
$SF stereo "CC=CC" --coords "$ALK_COORDS"

# ─── 17. Solvation ────────────────────────────────────────────────────────────

# Non-polar solvation (SASA ASP model)
$SF solvation "$ELEMENTS" "$COORDS"

# With explicit charges for GB
WATER_CHARGES='[-0.834, 0.417, 0.417]'
$SF solvation "$ELEMENTS" "$COORDS" --charges "$WATER_CHARGES"

# Custom probe radius
$SF solvation "$ELEMENTS" "$COORDS" --probe-radius 1.6

# ─── 18. SSSR — Ring Perception ───────────────────────────────────────────────

# Benzene (1 ring, size 6)
$SF sssr "c1ccccc1"

# Naphthalene (2 fused 6-membered rings)
$SF sssr "c1ccc2ccccc2c1"

# Cyclopentadiene (5-membered ring)
$SF sssr "C1=CC=CC1"

# Spiro compound
$SF sssr "C1CCC2(CC1)CCCC2"

# ─── 19. ECFP Fingerprints ────────────────────────────────────────────────────

# ECFP4 (radius=2, 2048 bits)
$SF ecfp "c1ccccc1"
$SF ecfp "c1ccccc1" --radius 2 --n-bits 2048

# ECFP2 (radius=1, 1024 bits)
$SF ecfp "CC(=O)O" --radius 1 --n-bits 1024

# ─── 20. Tanimoto Similarity ──────────────────────────────────────────────────

# Benzene vs Toluene (should be > 0.5)
$SF tanimoto "c1ccccc1" "Cc1ccccc1"

# Benzene vs Acetic acid (should be low)
$SF tanimoto "c1ccccc1" "CC(=O)O"

# Identical molecules (should be 1.0)
$SF tanimoto "c1ccccc1" "c1ccccc1"

# Custom radius / bit length
$SF tanimoto "c1ccccc1" "Cc1ccccc1" --radius 3 --n-bits 1024

# ─── Combined Pipeline Examples ───────────────────────────────────────────────

# Pipeline: embed → get coords → run EHT
COORDS_JSON=$($SF embed "O" --seed 42 | python3 -c "
import sys, json
r = json.load(sys.stdin)
print(json.dumps(r['coords']))
")
ELEMS_JSON='[8,1,1]'
$SF eht "$ELEMS_JSON" "$COORDS_JSON"

# Pipeline: batch embed → count successes
$SF batch --input smiles.txt | python3 -c "
import sys, json
results = json.load(sys.stdin)
ok = sum(1 for r in results if not r.get('error'))
print(f'{ok}/{len(results)} succeeded')
"

# Mass spectroscopy pipeline
for SMI in "c1ccccc1" "CC(=O)O" "CCN"; do
    RESULT=$($SF embed "$SMI" --seed 42)
    COORDS=$(echo "$RESULT" | python3 -c "import sys,json; print(json.dumps(json.load(sys.stdin)['coords']))")
    ELEMS=$(echo "$RESULT" | python3 -c "import sys,json; print(json.dumps(json.load(sys.stdin)['elements']))")
    echo "=== $SMI ==="
    $SF charges "$SMI"
    $SF sasa "$ELEMS" "$COORDS"
done
