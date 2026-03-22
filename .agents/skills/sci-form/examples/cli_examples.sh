#!/usr/bin/env bash

set -euo pipefail

SF="sci-form"
ELEMENTS='[8,1,1]'
COORDS='[0,0,0,0.96,0,0,-0.24,0.93,0]'

# Geometry
$SF embed "CCO"
$SF parse "CC(=O)O"
$SF info

# Batch embedding
printf 'O
CCO
c1ccccc1
' | $SF batch --threads 0 --format json

# Properties and electronic structure
$SF charges "CC(=O)O"
$SF population "$ELEMENTS" "$COORDS"
$SF dipole "$ELEMENTS" "$COORDS"
$SF dos "$ELEMENTS" "$COORDS" --sigma 0.3 --e-min -30 --e-max 5 --n-points 500
$SF ani "$ELEMENTS" "$COORDS"
$SF hf3c "$ELEMENTS" "$COORDS"

# Materials
$SF cell --a 10 --b 10 --c 10
$SF assemble --topology pcu --metal 30 --geometry octahedral --a 10.0

# Stereo and solvation
$SF stereo "C(F)(Cl)(Br)I"
$SF solvation "$ELEMENTS" "$COORDS"

# Rings and similarity
$SF sssr "c1ccccc1"
$SF ecfp "c1ccccc1" --radius 2 --n-bits 2048
$SF tanimoto "c1ccccc1" "Cc1ccccc1" --radius 2 --n-bits 2048

# Simple pipeline: embed once, reuse coords/elements for downstream commands.
EMBED_JSON=$($SF embed "CCO" --seed 42)
MOLECULE_COORDS=$(printf '%s' "$EMBED_JSON" | python3 -c 'import json,sys; print(json.dumps(json.load(sys.stdin)["coords"]))')
MOLECULE_ELEMENTS=$(printf '%s' "$EMBED_JSON" | python3 -c 'import json,sys; print(json.dumps(json.load(sys.stdin)["elements"]))')

$SF population "$MOLECULE_ELEMENTS" "$MOLECULE_COORDS"
$SF dipole "$MOLECULE_ELEMENTS" "$MOLECULE_COORDS"
