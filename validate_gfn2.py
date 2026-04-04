#!/usr/bin/env python3
"""Compare sci-form GFN2 vs tblite GFN2 for several molecules."""
import subprocess, json, sys, os
import numpy as np

os.environ.setdefault("LD_PRELOAD", "/usr/lib/x86_64-linux-gnu/libgfortran.so.5")

from tblite.interface import Calculator

EV_PER_HARTREE = 27.21138505
KCAL_PER_HARTREE = 627.509474
SCIFORM_BIN = "./target/release/sci-form"

# Phase 1: Fixed-coordinate test cases
fixed_cases = {
    "H2O": {
        "smiles": "O",
        "elements": [8, 1, 1],
        "positions": [[0.0, 0.0, 0.11779], [0.0, 0.75545, -0.47116], [0.0, -0.75545, -0.47116]],
    },
    "CH4": {
        "smiles": "C",
        "elements": [6, 1, 1, 1, 1],
        "positions": [
            [0.0, 0.0, 0.0],
            [0.6276, 0.6276, 0.6276],
            [-0.6276, -0.6276, 0.6276],
            [-0.6276, 0.6276, -0.6276],
            [0.6276, -0.6276, -0.6276],
        ],
    },
    "H2": {
        "smiles": "[H][H]",
        "elements": [1, 1],
        "positions": [[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]],
    },
    "NH3": {
        "smiles": "N",
        "elements": [7, 1, 1, 1],
        "positions": [
            [0.0000, 0.0000, 0.1173],
            [0.0000, 0.9377, -0.2737],
            [0.8122, -0.4689, -0.2737],
            [-0.8122, -0.4689, -0.2737],
        ],
    },
    "C2H6": {
        "smiles": "CC",
        "elements": [6, 6, 1, 1, 1, 1, 1, 1],
        "positions": [
            [0.0, 0.0, 0.7680],
            [0.0, 0.0, -0.7680],
            [0.0, 1.0186, 1.1572],
            [0.8821, -0.5093, 1.1572],
            [-0.8821, -0.5093, 1.1572],
            [0.0, -1.0186, -1.1572],
            [-0.8821, 0.5093, -1.1572],
            [0.8821, 0.5093, -1.1572],
        ],
    },
}

# Phase 2: SMILES to embed (generate coordinates then compare)
embed_smiles = [
    "CCO",          # ethanol
    "CC=O",         # acetaldehyde
    "CC(=O)O",      # acetic acid
    "c1ccccc1",     # benzene
    "C1CCCCC1",     # cyclohexane
    "CC(C)C",       # isobutane
    "C=C",          # ethylene
    "C#C",          # acetylene
    "CCN",          # ethylamine
    "OC=O",         # formic acid
]

ELEMENT_SYMBOLS = {1: "H", 5: "B", 6: "C", 7: "N", 8: "O", 9: "F", 14: "Si", 15: "P", 16: "S", 17: "Cl", 35: "Br"}


def run_tblite_gfn2(elements, positions):
    numbers = np.array(elements)
    pos_bohr = np.array(positions) / 0.529177
    calc = Calculator("GFN2-xTB", numbers, pos_bohr)
    old_stdout_fd = os.dup(1)
    devnull_fd = os.open(os.devnull, os.O_WRONLY)
    os.dup2(devnull_fd, 1)
    try:
        res = calc.singlepoint()
    finally:
        os.dup2(old_stdout_fd, 1)
        os.close(old_stdout_fd)
        os.close(devnull_fd)
    return res.get("energy")


def run_sciform_gfn2_coords(smiles, positions):
    """Run sci-form GFN2 with explicit coordinates, return energy in Hartree."""
    flat = []
    for p in positions:
        flat.extend(p)
    coords_json = json.dumps(flat)
    cmd = [SCIFORM_BIN, "neb-energy", smiles, coords_json, "--method", "gfn2"]
    result = subprocess.run(cmd, capture_output=True, text=True, cwd="/home/lestad/github/sci-form")
    try:
        data = json.loads(result.stdout)
        if isinstance(data, list):
            data = data[0]
        e_kcal = data.get("energy_kcal_mol")
        if e_kcal is not None:
            return e_kcal / KCAL_PER_HARTREE
        return None
    except:
        return None


def embed_sciform(smiles):
    """Embed a SMILES and return (elements, positions, coords_flat)."""
    cmd = [SCIFORM_BIN, "embed", smiles]
    result = subprocess.run(cmd, capture_output=True, text=True, cwd="/home/lestad/github/sci-form")
    try:
        data = json.loads(result.stdout)
        elements = data["elements"]
        coords = data["coords"]
        positions = []
        for i in range(0, len(coords), 3):
            positions.append([coords[i], coords[i+1], coords[i+2]])
        return elements, positions, coords
    except:
        sys.stderr.write(f"Embed failed for {smiles}: {result.stderr[:200]}\n")
        return None, None, None


print("=" * 80)
print("GFN2-xTB Validation: sci-form vs tblite")
print("=" * 80)

results = []

# Phase 1: Fixed coordinates
print("\n--- Phase 1: Fixed-coordinate molecules ---")
for name, case in fixed_cases.items():
    elements = case["elements"]
    positions = case["positions"]
    smiles = case["smiles"]

    e_tblite = run_tblite_gfn2(elements, positions)
    e_sciform = run_sciform_gfn2_coords(smiles, positions)

    if e_sciform is None:
        print(f"{name:12s}  FAILED")
        results.append(("FAIL", name, None))
        continue

    diff = abs(e_sciform - e_tblite)
    pct = abs(diff / e_tblite) * 100
    status = "OK" if pct < 0.1 else ("CLOSE" if pct < 1.0 else "FAIL")
    results.append((status, name, pct))
    print(f"{name:12s}  tblite={e_tblite:12.6f}  sci-form={e_sciform:12.6f}  diff={diff:.6f} Ha  ({pct:.4f}%)  [{status}]")

# Phase 2: Embed + compare
print("\n--- Phase 2: Embedded molecules ---")
for smi in embed_smiles:
    elements, positions, coords_flat = embed_sciform(smi)
    if elements is None:
        print(f"{smi:15s}  EMBED FAILED")
        results.append(("FAIL", smi, None))
        continue

    e_tblite = run_tblite_gfn2(elements, positions)
    e_sciform = run_sciform_gfn2_coords(smi, positions)

    if e_sciform is None:
        print(f"{smi:15s}  SCIFORM FAILED")
        results.append(("FAIL", smi, None))
        continue

    diff = abs(e_sciform - e_tblite)
    pct = abs(diff / e_tblite) * 100
    status = "OK" if pct < 0.1 else ("CLOSE" if pct < 1.0 else "FAIL")
    results.append((status, smi, pct))
    n_atoms = len(elements)
    print(f"{smi:15s}  ({n_atoms:2d} atoms)  tblite={e_tblite:12.6f}  sci-form={e_sciform:12.6f}  diff={diff:.6f} Ha  ({pct:.4f}%)  [{status}]")

# Summary
print("\n" + "=" * 80)
ok_count = sum(1 for s, _, _ in results if s == "OK")
close_count = sum(1 for s, _, _ in results if s == "CLOSE")
fail_count = sum(1 for s, _, _ in results if s == "FAIL")
total = len(results)
pcts = [p for _, _, p in results if p is not None]
avg_pct = np.mean(pcts) if pcts else 0
max_pct = max(pcts) if pcts else 0
print(f"OK: {ok_count}/{total}  CLOSE: {close_count}/{total}  FAIL: {fail_count}/{total}")
print(f"Average error: {avg_pct:.4f}%  Max error: {max_pct:.4f}%")
print("=" * 80)
